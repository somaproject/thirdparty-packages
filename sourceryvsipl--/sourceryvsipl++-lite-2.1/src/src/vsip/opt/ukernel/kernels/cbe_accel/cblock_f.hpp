/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/cblock_f.hpp
    @author  Jules Bergmann
    @date    2008-12-12
    @brief   VSIPL++ Library: User-defined kernel example illustrating
                              a control-block approach.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <complex>
#include <spu_mfcio.h>

#include <cml.h>
#include <cml_core.h>

#include <vsip/opt/ukernel/cbe_accel/debug.hpp>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/cblock_param.hpp>

#define DEBUG 0


/***********************************************************************
  Definitions
***********************************************************************/

struct Cblock_kernel : Spu_kernel
{
  static unsigned int const in_argc  = 0;
  static unsigned int const out_argc = 0;
  typedef Uk_cblock_params param_type;

  void init_rank(int r, int n)
  {
    rank = r;
    nspe = n;
  }

  void init(param_type& params)
  {
    in0        = params.in0;
    in0_stride = params.in0_stride;
    in1        = params.in1;
    in1_stride = params.in1_stride;
    in2        = params.in2;
    in2_stride = params.in2_stride;
    out        = params.out;
    out_stride = params.out_stride;
    rows       = params.rows;
    cols       = params.cols;
  }

  void compute()
  {
    // Compute my portion to compute
    int my_rows = rows / nspe + (rank < rows % nspe);
    int offset  = rank * (rows / nspe) + std::min(rank, rows % nspe);

#if DEBUG
    printf("Compute (%d/%d %d, %d) %d/%d\n", my_rows, rows, offset,
	   cols, rank, nspe);
#endif
    int tag = 23;

    uint64_t pin0 = in0 + offset * cols * sizeof(float);
    uint64_t pin1 = in1 + offset * cols * sizeof(float);
    uint64_t pin2 = in2 + offset * cols * sizeof(float);
    uint64_t pout = out + offset * cols * sizeof(float);

    float  buf[4*cols];
    float* buf0 = buf + 0*cols;
    float* buf1 = buf + 1*cols;
    float* buf2 = buf + 2*cols;
    float* buf3 = buf + 3*cols;

    for (int r=0; r<my_rows; ++r)
    {
      mfc_get(buf0, pin0, cols*sizeof(float), tag, 0, 0);
      mfc_get(buf1, pin1, cols*sizeof(float), tag, 0, 0);
      mfc_get(buf2, pin2, cols*sizeof(float), tag, 0, 0);

      pin0 += cols * sizeof(float);
      pin1 += cols * sizeof(float);
      pin2 += cols * sizeof(float);

      // Wait for DMAs to complete
      mfc_write_tag_mask(1<<tag);
      mfc_read_tag_status_all();

      for (int c=0; c<cols; ++c)
	buf3[c] = buf0[c] * buf1[c] + buf2[c];

      mfc_put(buf3, pout, cols*sizeof(float), tag, 0, 0);
      pout += cols * sizeof(float);
    }

    mfc_write_tag_mask(1<<tag);
    mfc_read_tag_status_all();
  }

  int32_t  rank;
  int32_t  nspe;
  uint64_t in0;
  int32_t  in0_stride;
  uint64_t in1;
  int32_t  in1_stride;
  uint64_t in2;
  int32_t  in2_stride;
  uint64_t out;
  int32_t  out_stride;
  int32_t  rows;
  int32_t  cols;
};
