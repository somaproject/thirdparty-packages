/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/zvmul_f.hpp
    @author  Jules Bergmann
    @date    2008-06-12
    @brief   VSIPL++ Library: UKernel to compute split-complex vmul.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <cml.h>
#include <cml_core.h>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>


#define DEBUG 0

/***********************************************************************
  Definitions
***********************************************************************/

struct Zvmul_kernel : Spu_kernel
{
  typedef std::pair<float*, float*> in0_type;
  typedef std::pair<float*, float*> in1_type;
  typedef std::pair<float*, float*> out0_type;

  static unsigned int const in_argc  = 2;
  static unsigned int const out_argc = 1;

  static bool const in_place = true;

  void compute(
    in0_type  const& in0,
    in1_type  const& in1,
    out0_type const& out,
    Pinfo const&     p_in0,
    Pinfo const&     p_in1,
    Pinfo const&     p_out)
  {
    cml_zvmul1_f(in0.first, in0.second,
		 in1.first, in1.second,
		 out.first, out.second,
		 p_out.l_total_size);
  }
};

typedef Zvmul_kernel kernel_type;
