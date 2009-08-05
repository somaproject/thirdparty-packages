/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/params/clobk_param.hpp
    @author  Jules Bergmann
    @date    2008-12-12
    @brief   VSIPL++ Library: Parameters for cblock user defined kernel
                              example.
*/

#ifndef VSIP_OPT_UKERNEL_KERNELS_PARAMS_CBLOCK_PARAM_HPP
#define VSIP_OPT_UKERNEL_KERNELS_PARAMS_CBLOCK_PARAM_HPP

typedef struct
{
  uint64_t      in0;
  int32_t       in0_stride;
  uint64_t      in1;
  int32_t       in1_stride;
  uint64_t      in2;
  int32_t       in2_stride;
  uint64_t      out;
  int32_t       out_stride;
  int32_t       rows, cols;
} Uk_cblock_params;

#endif // VSIP_OPT_UKERNEL_KERNELS_PARAMS_CBLOCK_PARAM_HPP
