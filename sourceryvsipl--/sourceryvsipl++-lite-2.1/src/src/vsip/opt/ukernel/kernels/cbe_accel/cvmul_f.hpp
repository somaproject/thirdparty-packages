/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/cvmul_f.hpp
    @author  Jules Bergmann
    @date    2008-06-12
    @brief   VSIPL++ Library: UKernel to compute inter-complex vmul.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <complex>
#include <cml.h>
#include <cml_core.h>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

struct Cvmul_kernel : Spu_kernel
{
  typedef std::complex<float>* in0_type;
  typedef std::complex<float>* in1_type;
  typedef std::complex<float>* out0_type;

  static unsigned int const in_argc  = 2;
  static unsigned int const out_argc = 1;

  static bool const in_place = true;

  void compute(
    in0_type     in0,
    in1_type     in1,
    out0_type    out,
    Pinfo const& p_in0,
    Pinfo const& p_in1,
    Pinfo const& p_out)
  {
    cml_cvmul1_f((float const*)in0, (float const*)in1,
		 (float*)out, p_out.l_total_size);
  }
};
