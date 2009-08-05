/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/fused_kernel.hpp
    @author  Jules Bergmann
    @date    2008-06-24
    @brief   VSIPL++ Library: Fused ukernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <cstdio>
#include <cml.h>
#include <cml_core.h>

#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

template <typename K1,
	  typename K2,
	  typename K3>
struct Fused_kernel : Spu_kernel
{
  // typedef std::complex<float>* in0_type;
  // typedef std::complex<float>* in1_type;
  // typedef std::complex<float>* out0_type;
  typedef typename K1::in0_type in0_type;
  typedef typename K2::in0_type in1_type;
  typedef typename K3::out0_type out0_type;

  static unsigned int const pre_argc  = 1;
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;

  typedef Uk_fused_params<typename K1::param_type,
			  typename K2::param_type,
			  typename K3::param_type>
		param_type;

  static bool const in_place = true;

  void init(param_type& params)
  {
    kernel1.init(params.k1_params);
    kernel2.init(params.k2_params);
    kernel3.init(params.k3_params);
  }

  void pre_compute(
    in0_type     in0,
    Pinfo const& p_in0)
  {
    kernel2.pre_compute(in0, p_in0);
  }

  void compute(
    in0_type     in,
    out0_type    out,
    Pinfo const& p_in,
    Pinfo const& p_out)
  {
    in0_type arg;
    kernel1.compute(in, out, p_in, p_in);
    kernel2.compute(out, in, p_in, p_in);
    kernel3.compute(in, out, p_in, p_out);
  }

  void fini()
  {
    kernel1.fini();
    kernel2.fini();
    kernel3.fini();
  }

  K1 kernel1;
  K2 kernel2;
  K3 kernel3;
};
