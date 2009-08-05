/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/params/fused_param.hpp
    @author  Jules Bergmann
    @date    2008-06-13
    @brief   VSIPL++ Library: Parameters for Fused kernels.
*/

#ifndef VSIP_OPT_UKERNEL_KERNELS_PARAMS_FUSED_PARAM_HPP
#define VSIP_OPT_UKERNEL_KERNELS_PARAMS_FUSED_PARAM_HPP

template <typename Param1,
	  typename Param2,
	  typename Param3>
struct Uk_fused_params
{
  Param1 k1_params;
  Param2 k2_params;
  Param3 k3_params;
};

#endif // VSIP_OPT_UKERNEL_KERNELS_PARAMS_FUSED_PARAM_HPP
