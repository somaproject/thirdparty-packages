/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/id2_f.hpp
    @author  Jules Bergmann
    @date    2008-07-29
    @brief   VSIPL++ Library: Kernel to ID2.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/id2_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

struct Id2_kernel : Spu_kernel
{
  typedef float* in0_type;
  typedef float* out0_type;

  typedef Uk_id2_params param_type;

  void init(param_type& params)
  {
    rows = params.rows;
    cols = params.cols;
  }

  void compute(
    in0_type      in,
    out0_type     out,
    Pinfo const&  p_in,
    Pinfo const&  p_out)
  {
    int size0   = p_out.l_size[0];
    int size1   = p_out.l_size[1];
    int stride0 = p_out.l_stride[0];
    int stride1 = p_out.l_stride[1];

    for (int i = 0; i < size0; ++i)
      for (int j = 0; j < size1; ++j)
	out[i * stride0 + j] = in[i * stride0 + j]
	  + (p_out.g_offset[0] + i) * cols
	  + (p_out.g_offset[1] + j);
  }

  unsigned int rows, cols;
};
