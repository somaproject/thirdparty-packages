/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/id1.hpp
    @author  Jules Bergmann
    @date    2008-06-18
    @brief   VSIPL++ Library: Vector copy Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID1_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID1_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Id1_ukernel : public vsip::impl::ukernel::Host_kernel_base
{
  // Parameters.
  //  - 'tag_type' is used to select the appropriate kernel (via
  //    Task_manager's Task_map)
  //  - in_argc and out_argc describe the number of input and output
  //    streams.
  //  - param_type (inherited) defaults to 'Empty_params'.
public:
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Id1_ukernel()
    : sp    (vsip::impl::ukernel::Blocksize_sdist(1024, 256))
  {}



  // Host-side compute kernel.  Used if accelerator is not available.

  template <typename View1,
	    typename View2>
  void compute(
    View1 in,
    View2 out)
  {
    out = in;
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

DEFINE_UKERNEL_TASK(Id1_ukernel, void(float*, float*), "uk_plugin", id1_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID1_HPP
