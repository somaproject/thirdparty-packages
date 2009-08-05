/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/madd.hpp
    @author  Don McCoy
    @date    2008-08-26
    @brief   VSIPL++ Library: User-defined kernel for multiply-add.
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_MADD_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_MADD_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise multiply-add ukernel.

class Madd_kernel : public vsip::impl::ukernel::Host_kernel_base
{
  // Parameters.
  //  - 'tag_type' is used to select the appropriate kernel (via
  //    Task_manager's Task_map)
  //  - in_argc and out_argc describe the number of input and output
  //    streams.
  //  - param_type (inherited) defaults to 'Empty_params'.
public:
  static unsigned int const in_argc  = 3;
  static unsigned int const out_argc = 1;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides matrix into whole, single rows.

  Madd_kernel()
    : sp (vsip::impl::ukernel::Blocksize_sdist(1), 
          vsip::impl::ukernel::Whole_sdist())
  {}



  // Host-side compute kernel.  Used if accelerator is not available.

  template <typename View0,
	    typename View1,
	    typename View2,
	    typename View3>
  void compute(
    View0 in0,
    View1 in1,
    View2 in2,
    View3 out)
  {
    out = in0 * in1 + in2;
  }


  // Queury API:
  // - in_spatt()/out_spatt() allow VSIPL++ to determine streaming
  //   pattern for user-kernel.  Since both input and output have same
  //   streaming pattern, simply return 'sp'

  vsip::impl::ukernel::Stream_pattern const& in_spatt(vsip::index_type i) const
  { return sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return sp; }


  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  vsip::impl::ukernel::Stream_pattern sp;	
};

DEFINE_UKERNEL_TASK(
  Madd_kernel,
  void(float*, float*, float*, float*),
  "uk_plugin", madd_f)

DEFINE_UKERNEL_TASK(
  Madd_kernel,
  void(std::complex<float>*, std::complex<float>*, std::complex<float>*, 
    std::complex<float>*),
  "uk_plugin", cmadd_f)

DEFINE_UKERNEL_TASK(
  Madd_kernel,
  void(float*, std::complex<float>*, std::complex<float>*, 
    std::complex<float>*),
  "uk_plugin", scmadd_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_MADD_HPP
