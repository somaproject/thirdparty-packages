/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/params/box2_param.hpp
    @author  Jules Bergmann
    @date    2008-08-01
    @brief   VSIPL++ Library: Parameters for Box2 kernels.
*/

#ifndef VSIP_OPT_UKERNEL_KERNELS_PARAMS_BOX2_PARAM_HPP
#define VSIP_OPT_UKERNEL_KERNELS_PARAMS_BOX2_PARAM_HPP

typedef struct
{
  unsigned int  overlap0;
  unsigned int  overlap1;
} Uk_box2_params;

#endif // VSIP_OPT_UKERNEL_KERNELS_PARAMS_BOX2_PARAM_HPP
