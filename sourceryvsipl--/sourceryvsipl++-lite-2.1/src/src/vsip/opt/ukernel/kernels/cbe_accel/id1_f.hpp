/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/id1_f.hpp
    @author  Jules Bergmann
    @date    2008-01-23
    @brief   VSIPL++ Library: Kernel to perform vector copy.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

struct Id1_kernel : Spu_kernel
{
  typedef float* in0_type;
  typedef float* out0_type;

  void compute(
    in0_type      in,
    out0_type     out,
    Pinfo const&  p_in,
    Pinfo const&  p_out)
  {
    int length = p_out.l_total_size;

    for (int i = 0; i < length; ++i)
      out[i] = in[i] + p_out.g_offset[0] + i;
  }
};
