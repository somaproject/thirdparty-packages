/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/cvmmul_f.hpp
    @author  Jules Bergmann
    @date    2008-06-12
    @brief   VSIPL++ Library: Split-complex vector-matrix multiply UKernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <complex>

#include <cml.h>
#include <cml_core.h>

#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/vmmul_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
struct Vmmul_kernel : Spu_kernel
{
  typedef T* in0_type;
  typedef T* in1_type;
  typedef T* out0_type;

  static unsigned int const pre_argc = 1;
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;
  typedef Uk_vmmul_params param_type;

  static bool const in_place = true;

  void pre_compute(
    in0_type     in0,
    Pinfo const& p_in0)
  {
    cml_core_vcopy1_f((float const*)in0, (float*)buf, 2*p_in0.l_total_size);
  }

  void compute(
    in1_type     in1,
    out0_type    out,
    Pinfo const& p_in1,
    Pinfo const& p_out)
  {
    cml_cvmul1_f((float const*)buf, (float const*)in1,
		 (float*)out, p_out.l_total_size);
  }

// Member data
  vector float buf[2*4096/4];	// Use vector float to force alignment
};
