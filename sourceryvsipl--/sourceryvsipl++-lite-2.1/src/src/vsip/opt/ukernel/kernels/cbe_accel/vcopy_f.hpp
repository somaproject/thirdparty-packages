/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/vcopy_f.hpp
    @author  Jules Bergmann
    @date    2008-06-10
    @brief   VSIPL++ Library: Ukernel to perform vector copy.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <cstdio>
#include <stdint.h>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>

// Cell specific
#define DMA_ALIGNMENT 16
#define DMA_ALIGNMENT_OF(A) ((uintptr_t)(A) & (DMA_ALIGNMENT - 1))
#define IS_DMA_ALIGNED(A) (DMA_ALIGNMENT_OF(A) == 0)


/***********************************************************************
  Definitions
***********************************************************************/

struct Vcopy_kernel : Spu_kernel
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
    int offset = p_in.l_offset[0];

    if (IS_DMA_ALIGNED(in + offset))
    {
      vector float* a = (vector float*)(in + offset);
      vector float* z = (vector float*)out;

      for (int i = 0; i < (length / 4); ++i)
	*z++ = *a++;
    }
    else
    {
      for (int i = 0; i < length; ++i)
	out[i] = in[offset + i];
    }
  }
};
