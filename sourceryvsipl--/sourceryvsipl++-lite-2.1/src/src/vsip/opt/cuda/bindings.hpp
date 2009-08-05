/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/bindings.hpp
    @author  Don McCoy
    @date    2009-02-05
    @brief   VSIPL++ Library: Bindings for CUDA's BLAS functions and
               for custom CUDA kernels.
*/

#ifndef VSIP_OPT_CUDA_BINDINGS_HPP
#define VSIP_OPT_CUDA_BINDINGS_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <complex>

#include <vsip/support.hpp>
#include <vsip/opt/cuda/kernels.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>


#ifndef NDEBUG
#include <iostream>
#define ASSERT_CUDA_OK()					\
{								\
  cudaError_t error = cudaGetLastError();			\
  if (error != cudaSuccess)					\
  {								\
    std::cerr << "CUDA problem encountered (error "		\
	      << error << ")" << std::endl;			\
    std::cerr << cudaGetErrorString(error) << std::endl;	\
  }								\
  assert(error == cudaSuccess);					\
}

#define ASSERT_CUBLAS_OK()					\
{								\
  cuda::cublasStatus status = cuda::cublasGetError();		\
  if (status != 0)						\
    std::cerr << "CUBLAS problem encountered (error "		\
	      << status << ")" << std::endl;			\
  assert(status == 0);						\
}

#define ASSERT_CUFFT_OK(result)					\
{								\
  if (result != 0)						\
    std::cerr << "CUFFT problem encountered (error "		\
	      << result << ")" << std::endl;			\
  assert(result == 0);						\
}

#else
#define ASSERT_CUDA_OK()
#define ASSERT_CUBLAS_OK()
#define ASSERT_CUFFT_OK(r)

#endif // CUDA_DEBUG


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{
// Functions to interface with CUDA
//
void initialize(int& argc, char**&argv);
void finalize();
  

extern "C"
{
// Prototypes for CUBLAS functions called directly (see cublas.h)
// 
typedef unsigned int cublasStatus;

float cublasSdot(int n, 
		 const float *x, int incx, 
		 const float *y, int incy);
float _Complex cublasCdotu(int n, 
			   const float _Complex *x, int incx, 
			   const float _Complex *y, int incy);
float _Complex cublasCdotc(int n, 
			   const float _Complex *x, int incx, 
			   const float _Complex *y, int incy);
cublasStatus cublasGetError();


// From cuda_runtime_api.h
//
enum cudaError
{ 
  cudaSuccess = 0
};
typedef cudaError cudaError_t;

enum cudaMemcpyKind
{
  cudaMemcpyHostToHost = 0,
  cudaMemcpyHostToDevice,
  cudaMemcpyDeviceToHost,
  cudaMemcpyDeviceToDevice
};

cudaError_t cudaMemcpy  (void *dst, const void *src, size_t count, enum cudaMemcpyKind kind);
cudaError_t cudaMemcpy2D(void *dst, size_t dpitch, const void *src, size_t spitch, size_t width, size_t height, enum cudaMemcpyKind kind);
cudaError_t cudaMalloc(void **devPtr, size_t size);
cudaError_t cudaFree(void *devPtr);
cudaError_t cudaGetLastError(void);
const char* cudaGetErrorString(cudaError_t error);
}


//
// C++ --> C interface functions
//

// cuda::dot()
#define VSIP_IMPL_CUBLAS_DOT(T, CUDA_T, VPPFCN, CUBLASFCN)	\
inline T							\
VPPFCN(int n,							\
    const T* x, int incx,					\
    const T* y, int incy)					\
{								\
  if (incx < 0) x += incx * (n-1);				\
  if (incy < 0) y += incy * (n-1);				\
  return CUBLASFCN(n, (const CUDA_T*)x, incx,                   \
                      (const CUDA_T*)y, incy);                  \
}

// Note: CUDA functions return the C99 Complex type.  The way the
// return value is handled when converting back to the C++ type relies
// on a GNU extension and may not work with all compilers.
VSIP_IMPL_CUBLAS_DOT(float,               float,          dot,  cublasSdot)
VSIP_IMPL_CUBLAS_DOT(std::complex<float>, float _Complex, dot,  cublasCdotu)
VSIP_IMPL_CUBLAS_DOT(std::complex<float>, float _Complex, dotc, cublasCdotc)
#undef VSIP_IMPL_CUBLAS_DOT


// Wrapper functions used for vmmul serial expression evaluator.
// These functions convert parameters to the proper (standard C) types
// and call the appropriate (non-overloaded) kernel entry point.

inline
void 
vmmul_row(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols)
{
  vmmul_row_ss(kernel, input, output, rows, cols);
}

inline
void 
vmmul_row(
  float const*               kernel,  
  std::complex<float> const* input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_row_sc(
    kernel,
    reinterpret_cast<cuComplex const*>(input),
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}

inline
void 
vmmul_row(
  std::complex<float> const* kernel,  
  float const*               input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_row_cs(
    reinterpret_cast<cuComplex const*>(kernel),
    input,
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}

inline
void 
vmmul_row(
  std::complex<float> const* kernel,  
  std::complex<float> const* input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_row_cc(
    reinterpret_cast<cuComplex const*>(kernel),
    reinterpret_cast<cuComplex const*>(input),
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}


inline
void 
vmmul_col(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols)
{
  vmmul_col_ss(kernel, input, output, rows, cols);
}

inline
void 
vmmul_col(
  float const*               kernel,  
  std::complex<float> const* input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_col_sc(
    kernel,
    reinterpret_cast<cuComplex const*>(input),
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}

inline
void 
vmmul_col(
  std::complex<float> const* kernel,  
  float const*               input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_col_cs(
    reinterpret_cast<cuComplex const*>(kernel),
    input,
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}

inline
void 
vmmul_col(
  std::complex<float> const* kernel,  
  std::complex<float> const* input,
  std::complex<float>*       output,
  size_t                     rows,
  size_t                     cols)
{
  vmmul_col_cc(
    reinterpret_cast<cuComplex const*>(kernel),
    reinterpret_cast<cuComplex const*>(input),
    reinterpret_cast<cuComplex*>(output),
    rows, cols);
}


inline
void
copy_device_memory(
  float const* src, 
  float* dest, 
  size_t size)
{
  copy_device_to_device(src, dest, size);
}

inline
void
copy_device_memory(
  std::complex<float> const* src, 
  std::complex<float>* dest, 
  size_t size)
{
  copy_device_to_device(
    reinterpret_cast<float const*>(src),
    reinterpret_cast<float*>(dest),
    size * 2);
}


inline
void
zero_device_memory(
  float* dest, 
  size_t size)
{
  copy_zeroes_to_device(dest, size);
}

inline
void
zero_device_memory(
  std::complex<float>* dest, 
  size_t size)
{
  copy_zeroes_to_device(
    reinterpret_cast<float*>(dest),
    size * 2);
}



/// Lightweight class used to copy data between blocks.  This class is
/// used to overload functions needed for differently dimensioned blocks.
/// The functions copy_host_to_dev() and copy_dev_to_host() should be
/// used to invoke these member functions rather than calling them directly.
template <dimension_type Dim,
          typename Order,
          typename Block>
struct copy_block;

template <typename Order, typename Block>
struct copy_block<1, Order, Block>
{
  typedef typename Block::value_type T;
  static inline void host_to_dev(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy(dest, src, block.size() * sizeof(T), cudaMemcpyHostToDevice);
    ASSERT_CUDA_OK();
  }
  static inline void dev_to_host(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy(dest, src, block.size() * sizeof(T), cudaMemcpyDeviceToHost);
    ASSERT_CUDA_OK();
  }
};

template <typename Block>
struct copy_block<2, row2_type, Block>
{
  typedef typename Block::value_type T;
  typedef typename Block_layout<Block>::order_type order_type;
  static void host_to_dev(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy2D( 
      dest, block.size(2, 1) * sizeof(T),
      src, block.impl_stride(2, 0) * sizeof(T),
      block.size(2, 1) * sizeof(T), block.size(2, 0),
      cudaMemcpyHostToDevice);
    ASSERT_CUDA_OK();
  }
  static void dev_to_host(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy2D( 
      dest, block.size(2, 1) * sizeof(T),
      src, block.impl_stride(2, 0) * sizeof(T),
      block.size(2, 1) * sizeof(T), block.size(2, 0),
      cudaMemcpyDeviceToHost);
    ASSERT_CUDA_OK();
  }
};


template <typename Block>
struct copy_block<2, col2_type, Block>
{
  typedef typename Block::value_type T;
  static inline void host_to_dev(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy2D( 
      dest, block.size(2, 0) * sizeof(T),
      src, block.impl_stride(2, 1) * sizeof(T),
      block.size(2, 0) * sizeof(T), block.size(2, 1),
      cudaMemcpyHostToDevice);
    ASSERT_CUDA_OK();
  }
  static inline void dev_to_host(Block const& block, T const* src, T* dest) 
  {
    cudaMemcpy2D( 
      dest, block.size(2, 0) * sizeof(T),
      src, block.impl_stride(2, 1) * sizeof(T),
      block.size(2, 0) * sizeof(T), block.size(2, 1),
      cudaMemcpyDeviceToHost);
    ASSERT_CUDA_OK();
  }
};


/// External interface for copy_block<> template member function
/// that copies from CPU host memory to GPU device memory.
template <typename Block>
inline void
copy_host_to_dev(Block const& block,
		 typename Block::value_type const* src,
		 typename Block::value_type*       dest)
{
  copy_block<Block::dim,
             typename Block_layout<Block>::order_type,
             Block>().
    host_to_dev(block, src, dest);
}

/// External interface for copy_block<> template member function
/// that copies from GPU device memory to CPU host memory.
template <typename Block>
inline void
copy_dev_to_host(Block const& block,
		 typename Block::value_type const* src,
		 typename Block::value_type*       dest)
{
  copy_block<Block::dim,
             typename Block_layout<Block>::order_type,
             Block>().
    dev_to_host(block, src, dest);
}
 


 
/// CUDA capabilities known at compile time are expressed as traits.
///
/// Note that support for double precision is included in CUDA, but
/// not present on all hardware.  The correct way to fill out these
/// traits would be with an inquiry run at configure time, or an
/// option for cross compiling that forces it to assume one way or
/// the other.
///
template <typename T>
struct Cuda_traits
{
  static bool const valid = false;
};

template <>
struct Cuda_traits<float>
{
  static bool const valid = true;
  static char const trans = 't';
};

template <>
struct Cuda_traits<double>
{
  static bool const valid = false;
  static char const trans = 't';
};

template <>
struct Cuda_traits<std::complex<float> >
{
  static bool const valid = true;
  static char const trans = 'c';
};

template <>
struct Cuda_traits<std::complex<double> >
{
  static bool const valid = false;
  static char const trans = 'c';
};

} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_BINDINGS_HPP
