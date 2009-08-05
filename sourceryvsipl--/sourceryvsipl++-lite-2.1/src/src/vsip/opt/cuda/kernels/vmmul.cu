/* Copyright (c) 2009 by CodeSourcery.  All rights reserved. */

/** @file    vmmul.cu
    @author  Don McCoy
    @date    2009-03-10
    @brief   VSIPL++ Library: CUDA Kernel for Vector-Matrix Multiplication
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <cuComplex.h>



/***********************************************************************
  Device functions (callable only via kernels)
***********************************************************************/

// complex multiply with scale
//   d = (a * b) * c   where a, b and d are complex and c is real
__device__ void cmuls(cuComplex& d, cuComplex a, cuComplex b, float c)
{
  d.x = (a.x * b.x - a.y * b.y) * c;
  d.y = (a.y * b.x + a.x * b.y) * c;
}

// complex multiply
//   c = a * b   where a, b and c are complex
__device__ void cmul(cuComplex& c, cuComplex a, cuComplex b)
{
  c.x = a.x * b.x - a.y * b.y;
  c.y = a.y * b.x + a.x * b.y;
}

// scalar-complex multiply
//   c = a * b   where b and c are complex and a is real
__device__ void scmul(cuComplex& c, float a, cuComplex b)
{
  c.x = a * b.x;
  c.y = a * b.y;
}


/***********************************************************************
  Device Kernels -- Each thread computes one element
***********************************************************************/

// Support provided for some (but not all) combinations of real and complex, 
// single and double precision.  The BLAS notation convention is used to 
// indicate the types of the two inputs and the result, i.e.:
//
//   S = single precision real      D = double precision real
//   C = single precision complex   Z = double precision complex


// Vector-Matrix multiply         S * S --> S
__global__ void 
vmmul_ss(float const* kernel, float const* input, float* output, 
          size_t size, bool by_row)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;
  int const vec_idx = by_row ? col : row;

  int const idx = __mul24(row, size) + col;
  output[idx] = kernel[vec_idx] * input[idx];
}

// Vector-Matrix multiply         S * C --> C
__global__ void 
vmmul_sc(float const* kernel, cuComplex const* input, cuComplex* output, 
         size_t size, bool by_row)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;
  int const vec_idx = by_row ? col : row;
  
  int const idx = __mul24(row, size) + col;
  scmul(output[idx], kernel[vec_idx], input[idx]);
}

// Vector-Matrix multiply         C * S --> C
__global__ void 
vmmul_cs(cuComplex const* kernel, float const* input, cuComplex* output, 
         size_t size, bool by_row)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;
  int const vec_idx = by_row ? col : row;
  
  int const idx = __mul24(row, size) + col;
  scmul(output[idx], input[idx], kernel[vec_idx]);
}

// Vector-Matrix multiply         C * C --> C
__global__ void 
vmmul_cc(cuComplex const* kernel, cuComplex const* input, cuComplex* output, 
         size_t size, bool by_row)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;
  int const vec_idx = by_row ? col : row;
  
  int const idx = __mul24(row, size) + col;
  cmul(output[idx], kernel[vec_idx], input[idx]);
}


// Vector-Matrix multiply with scale		C * C * s --> C
//   Computes A = A * B * c  where A is a complex matrix, B is a complex 
//   vector and c is real.  This is used to combine the scaling step from 
//   an inverse FFT with the vector-multiplication step when doing fast 
//   convolution.
__global__ void 
vmmuls_cc(cuComplex const* kernel, cuComplex const* input, cuComplex* output, 
          float scale, size_t size, bool by_row)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;
  int const vec_idx = by_row ? col : row;

  int const idx = __mul24(row, size) + col;
  cmuls(output[idx], kernel[vec_idx], input[idx], scale);
}


// Matrix-Matrix multiply with scale		C * C * s --> C
//   Computes A = A * B * c  where A is a complex matrix, B is a complex 
//   vector and c is real.  This is used to combine the scaling step from 
//   an inverse FFT with the vector-multiplication step when doing fast 
//   convolution.
__global__ void 
matmuls_cc(cuComplex const* kernel, cuComplex const* input, cuComplex* output, 
           float scale, size_t size)
{
  int const tx = threadIdx.x;
  int const ty = threadIdx.y;
  int const bx = blockIdx.x;
  int const by = blockIdx.y;
  int const row = __mul24(blockDim.y, by) + ty;
  int const col = __mul24(blockDim.x, bx) + tx;

  int const idx = __mul24(row, size) + col;
  cmuls(output[idx], kernel[idx], input[idx], scale);
}



/***********************************************************************
  Helper functions
***********************************************************************/

dim3
calculate_thread_block_size(size_t rows, size_t cols)
{
  // Calculate the optimal number of sub-columns and sub-rows per thread
  // block.  A sub-block size of 16 x 32 yields 512 threads, the maximum 
  // allowed.
  dim3 threads(32, 16);

  while (threads.y > 1)
  {
    if (rows % threads.y) threads.y /= 2;
    else break;
  }

  while (threads.x > 1)
  {
    if (cols % threads.x) threads.x /= 2;
    else break;
  }

  return threads;
}



/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cuda
{

void 
vmmul_row_ss(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols)
{
  bool by_row = true;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_ss<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_row_sc(
  float const*     kernel,
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = true;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_sc<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_row_cs(
  cuComplex const* kernel,  
  float const*     input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = true;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_cs<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_row_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = true;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_cc<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}


void 
vmmul_col_ss(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols)
{
  bool by_row = false;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_ss<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_col_sc(
  float const*     kernel,
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = false;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_sc<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_col_cs(
  cuComplex const* kernel,  
  float const*     input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = false;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_cs<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}

void 
vmmul_col_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols)
{
  bool by_row = false;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmul_cc<<<blocks, threads>>>(kernel, input, output, cols, by_row);
}


void 
vmmuls_row_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols)
{
  bool by_row = true;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmuls_cc<<<blocks, threads>>>(kernel, input, output, scale, cols, by_row);
}

void 
vmmuls_col_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols)
{
  bool by_row = false;
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  vmmuls_cc<<<blocks, threads>>>(kernel, input, output, scale, cols, by_row);
}


void 
mmmuls_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols)
{
  dim3 threads = calculate_thread_block_size(rows, cols);
  dim3 blocks(cols/threads.x, rows/threads.y);
  matmuls_cc<<<blocks, threads>>>(kernel, input, output, scale, cols);
}

} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip
