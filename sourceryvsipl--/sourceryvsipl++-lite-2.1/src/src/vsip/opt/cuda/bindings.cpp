/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/bindings.cpp
    @author  Don McCoy
    @date    2009-02-05
    @brief   VSIPL++ Library: CUDA BLAS interface
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <ostream>

#include <vsip/core/extdata_common.hpp>
#include <vsip/opt/cuda/bindings.hpp>
#include <vsip/opt/cuda/gpu_block.hpp>

#include <cublas.h>


/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{

/// This function must be called prior to any other CUDA-related function.
/// Presently it only initializes the CUDA BLAS library.
void
initialize(int& /* argc */, char**& /* argv */)
{
  cublasStatus status = cublasInit();
  ASSERT_CUBLAS_OK();

  if (status != CUBLAS_STATUS_SUCCESS)
  {
    std::ostringstream message;
    message << "CUDA Library failed to initialize (error "
	    << status << ")" << std::endl;
    VSIP_IMPL_THROW(std::runtime_error(message.str()));
  }
}

/// Called to shut down and free all resources
void
finalize()
{
  cublasShutdown();
  ASSERT_CUBLAS_OK();
}

} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip
