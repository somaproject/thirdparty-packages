/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/fastconv.cpp
    @author  Don McCoy
    @date    2009-03-22
    @brief   VSIPL++ Library: Wrapper for fast convolution using CUDA
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/fns_scalar.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/static_assert.hpp>
#include <vsip/math.hpp>
#include <vsip/opt/cuda/fastconv.hpp>
#include <vsip/opt/cuda/gpu_block.hpp>

#include <cufft.h>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cuda
{

// Fast convolution binding for interleaved complex data.


template <dimension_type D,
          typename       T,
	  typename       ComplexFmt>
void
Fastconv_base<D, T, ComplexFmt>::fconv
  (T const* in, T const* kernel, T* out, length_type rows, length_type columns, bool transform_kernel)
{
  size_t kernel_size = (D == 1) ? columns : rows * columns;

  // allocate device memory and copy input and kernel over from host
  Device_storage<T> dev_out(rows * columns);
  Device_storage<T> dev_kernel(kernel_size);
  Device_storage<T> dev_in(rows * columns);

  // If the kernel is a matrix, it is assumed to be row-major and dense.
  // As a result, it can be copied as one contiguous chunk.
  cudaMemcpy(
    dev_kernel.data(),
    kernel,
    kernel_size * sizeof(T), 
    cudaMemcpyHostToDevice);
  ASSERT_CUDA_OK();

  // Transfer the input (row major, dense)
  cudaMemcpy(
    dev_in.data(),
    in,
    rows * columns * sizeof(T), 
    cudaMemcpyHostToDevice);
  ASSERT_CUDA_OK();
 

  // convert pointers to types the CUFFT library accepts
  typedef cufftComplex ctype;
  ctype* d_out = reinterpret_cast<ctype*>(dev_out.data());
  ctype* d_kernel = reinterpret_cast<ctype*>(dev_kernel.data());
  ctype* d_in = reinterpret_cast<ctype*>(dev_in.data());

  cufftHandle plan;
  if (transform_kernel)
  {
    // Create a 1D FFT plan and transform the kernel
    cufftPlan1d(&plan, columns, CUFFT_C2C, 1);
    cufftExecC2C(plan, d_kernel, d_kernel, CUFFT_FORWARD);
    cufftDestroy(plan);
  }

  // Create a FFTM plan
  cufftPlan1d(&plan, columns, CUFFT_C2C, rows);

  // transform the data
  cufftExecC2C(plan, d_in, d_in, CUFFT_FORWARD);

  // convolve with kernel, combine with scaling needed for inverse FFT
  typedef typename impl::Scalar_of<T>::type scalar_type;
  scalar_type scale = 1 / static_cast<scalar_type>(columns);
  if (D == 1)
    vmmuls_row_cc(d_kernel, d_in, d_out, scale, rows, columns);
  else
    mmmuls_cc(d_kernel, d_in, d_out, scale, rows, columns);

  // inverse transform the signal
  cufftExecC2C(plan, d_out, d_out, CUFFT_INVERSE);
  cufftDestroy(plan);

  // Move data back to the host from the output buffer
  cudaMemcpy(
    out,
    dev_out.data(),
    rows * columns * sizeof(T), 
    cudaMemcpyDeviceToHost);
  ASSERT_CUDA_OK();
}



typedef std::complex<float> ctype;

template void
Fastconv_base<1, ctype, Cmplx_inter_fmt>::fconv(
  ctype const* in, ctype const* kernel, ctype* out, 
  length_type rows, length_type columns, bool transform_kernel);

template void
Fastconv_base<2, ctype, Cmplx_inter_fmt>::fconv(
  ctype const* in, ctype const* kernel, ctype* out, 
  length_type rows, length_type columns, bool transform_kernel);


          
} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip
