/* Copyright (c) 2009 by CodeSourcery.  All rights reserved. 

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/device_storage.hpp
    @author  Don McCoy
    @date    2009-04-20
    @brief   VSIPL++ Library: Class for accessing on-device memory with CUDA
*/

#ifndef VSIP_OPT_CUDA_DEVICE_STORAGE_HPP
#define VSIP_OPT_CUDA_DEVICE_STORAGE_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/cuda/bindings.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{

template <typename T>
class Device_storage
{
public:
  Device_storage(length_type size) VSIP_NOTHROW
    : dev_(NULL)
  {
    cudaError_t error = cudaMalloc((void**)&dev_, size * sizeof(T));
    ASSERT_CUDA_OK();
    if (error != cudaSuccess)
      VSIP_IMPL_THROW(std::bad_alloc());
  }

  ~Device_storage()
  {
    cudaFree(dev_);
    ASSERT_CUDA_OK();
  }

  T*       data()       { return dev_; }
  T const* data() const { return dev_; }

private:
  T* dev_;
};

} // namespace cuda
} // namespace impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_DEVICE_STORAGE_HPP
