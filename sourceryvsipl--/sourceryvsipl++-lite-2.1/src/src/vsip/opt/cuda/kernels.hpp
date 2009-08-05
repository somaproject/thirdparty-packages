/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/kernels.hpp
    @author  Don McCoy
    @date    2009-03-11
    @brief   VSIPL++ Library: Custom CUDA kernels
*/

#ifndef VSIP_OPT_CUDA_KERNELS_HPP
#define VSIP_OPT_CUDA_KERNELS_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <cuComplex.h>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{

// Host-side functions
//
// These are simple wrappers around the actual CUDA kernels.  Note that 
// input data resides in device memory allocated with cudaMalloc() and 
// filled from the host using cudaMemcpy().  Upon return, data may be 
// copied back to host or operated on by other kernels.



/// Vector-Matrix multiply, row-wise, S * S --> S
void 
vmmul_row_ss(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols);

/// Vector-Matrix multiply, row-wise, S * C --> C
void 
vmmul_row_sc(
  float const*     kernel,
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);

/// Vector-Matrix multiply, row-wise, C * S --> C
void 
vmmul_row_cs(
  cuComplex const* kernel,  
  float const*     input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);

/// Vector-Matrix multiply, row-wise, C * C --> C
void 
vmmul_row_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);


/// Vector-Matrix multiply, column-wise, S * S --> S
void 
vmmul_col_ss(
  float const* kernel,  
  float const* input,
  float*       output,
  size_t       rows,
  size_t       cols);

/// Vector-Matrix multiply, column-wise, S * C --> C
void 
vmmul_col_sc(
  float const*     kernel,
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);

/// Vector-Matrix multiply, column-wise, C * S --> C
void 
vmmul_col_cs(
  cuComplex const* kernel,  
  float const*     input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);

/// Vector-Matrix multiply, column-wise, C * C --> C
void 
vmmul_col_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  size_t           rows,
  size_t           cols);


/// Vector-Matrix multiply, with scale, row-wise, C * C --> C
void 
vmmuls_row_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols);

/// Vector-Matrix multiply, with scale, column-wise, C * C --> C
void 
vmmuls_col_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols);


/// Matrix-Matrix multiply, with scale, row-wise, C * C --> C
void 
mmmuls_cc(
  cuComplex const* kernel,  
  cuComplex const* input,
  cuComplex*       output,
  float            scale,
  size_t           rows,
  size_t           cols);


/// Memory copy, from device (global) to shared (fast, on-core) memory
void
copy_device_to_shared(
  float* src, 
  size_t size);

/// Memory copy, from device to device memory
void
copy_device_to_device(
  float const* src, 
  float* dest, 
  size_t size);

/// Memory fill device with zeroes
void
copy_zeroes_to_device(
  float* dest, 
  size_t size);

} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_KERNELS_HPP
