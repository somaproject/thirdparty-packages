/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/box1.hpp
    @author  Jules Bergmann
    @date    2008-08-01
    @brief   VSIPL++ Library: Box1 Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_BOX1_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_BOX1_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/box1_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Box1_ukernel : public vsip::impl::ukernel::Host_kernel_base
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
  typedef Box1_params param_type;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Box1_ukernel(int overlap, int same_support)
    : in_sp  (same_support ?
	      vsip::impl::ukernel::Blockoverlap_sdist(32, 32, overlap, overlap, 0, 0) :
	      vsip::impl::ukernel::Blockoverlap_sdist(32, overlap))
    , out_sp (vsip::impl::ukernel::Blocksize_sdist(32))
    , overlap_(overlap)
  {}

  void fill_params(param_type& param) const
  {
    param.overlap = overlap_;
  }



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
  { return in_sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return out_sp; }

  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  vsip::impl::ukernel::Stream_pattern in_sp;	
  vsip::impl::ukernel::Stream_pattern out_sp;	
  int                                 overlap_;
};

DEFINE_UKERNEL_TASK(Box1_ukernel, void(float*, float*), "uk_plugin", box1_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_BOX1_HPP
