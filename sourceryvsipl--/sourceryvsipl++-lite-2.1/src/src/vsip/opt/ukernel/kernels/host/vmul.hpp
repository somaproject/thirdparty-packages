/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/vmul.hpp
    @author  Jules Bergmann
    @date    2008-06-24
    @brief   VSIPL++ Library: Vector element-wise multiply Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMUL_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMUL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Vmul_kernel : public vsip::impl::ukernel::Host_kernel_base
{
  // Parameters.
  //  - 'tag_type' is used to select the appropriate kernel (via
  //    Task_manager's Task_map)
  //  - in_argc and out_argc describe the number of input and output
  //    streams.
  //  - param_type (inherited) defaults to 'Empty_params'.
public:
  static unsigned int const in_argc  = 2;
  static unsigned int const out_argc = 1;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Vmul_kernel()
    : sp    (vsip::impl::ukernel::Blocksize_sdist(1024, 256))
  {}



  // Host-side compute kernel.  Used if accelerator is not available.

  template <typename View0,
	    typename View1,
	    typename View2>
  void compute(
    View0 in0,
    View1 in1,
    View2 out)
  {
    out = in0 * in1;
  }


  // Queury API:  in_spatt()/out_spatt() allow VSIPL++ to determine
  // streaming pattern for user-kernel.  Since both input and output
  // have same streaming pattern, simply return 'sp'

  vsip::impl::ukernel::Stream_pattern const& in_spatt(vsip::index_type) const
  { return sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return sp; }

  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  vsip::impl::ukernel::Stream_pattern sp;	
};

DEFINE_UKERNEL_TASK(Vmul_kernel,
		    void(float*, float*, float*), "uk_plugin", vmul_f)
DEFINE_UKERNEL_TASK(Vmul_kernel,
		    void(std::complex<float>*, std::complex<float>*, std::complex<float>*),
		    "uk_plugin", cvmul_f)
DEFINE_UKERNEL_TASK(Vmul_kernel,
                    void(std::pair<float*, float*>, std::pair<float*, float*>, std::pair<float*, float*>),
		    "uk_plugin", zvmul_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMUL_HPP
