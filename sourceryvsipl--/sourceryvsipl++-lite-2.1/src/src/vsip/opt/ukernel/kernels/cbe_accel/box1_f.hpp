/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/box1_f.hpp
    @author  Jules Bergmann
    @date    2008-01-23
    @brief   VSIPL++ Library: 1D Box filter UKernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/box1_param.hpp>

#define min(A, B) ((A)<(B)?(A):(B))



/***********************************************************************
  Definitions
***********************************************************************/

struct Box1_kernel : Spu_kernel
{
  typedef float* in0_type;
  typedef float* out0_type;

  typedef Box1_params param_type;

  void init(param_type& params)
  {
    overlap = params.overlap;
  }

  void compute(
    in0_type      in,
    out0_type     out,
    Pinfo const&  p_in,
    Pinfo const&  p_out)
  {
    int size   = p_in.l_total_size;
    int offset = p_in.l_offset[0];

#if 1
    for (int r = 0; r < size; ++r)
    {
      float sum = 0;

      int rstart = -min(overlap, p_in.o_leading[0] + r);
      int rend   =  min(overlap, size-r-1+p_in.o_trailing[0]) + 1;

      for (int rr = rstart; rr < rend; ++rr)
	sum += in[r+offset+rr];
      out[r] = sum;
    }
#else
    if (p_in.o_leading[0] > 0)
      out[0] = in[offset - 1] 
	     + in[offset] 
	     + in[offset + 1];
    else
      out[0] = in[offset] 
	     + in[offset + 1];

    for (int i = 1; i < length-1; ++i)
      out[i] = in[i + offset - 1] 
	     + in[i + offset] 
	     + in[i + offset + 1];

    if (p_in.o_trailing[0] > 0)
      out[length-1] = in[length-1 + offset - 1] 
	            + in[length-1 + offset] 
	            + in[length-1 + offset + 1];
    else
      out[length-1] = in[length-1 + offset] 
	            + in[length-1 + offset - 1];
#endif
  }

  int       overlap;
};
