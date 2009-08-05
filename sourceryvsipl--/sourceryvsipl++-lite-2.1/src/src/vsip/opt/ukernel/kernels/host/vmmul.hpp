/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/vmmul.hpp
    @author  Jules Bergmann
    @date    2008-06-27
    @brief   VSIPL++ Library: Vector-matrix element-wise multiply Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMMUL_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMMUL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/vmmul_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Vmmul_kernel : public vsip::impl::ukernel::Host_kernel_base
{
  // Parameters.
  //  - 'tag_type' is used to select the appropriate kernel (via
  //    Task_manager's Task_map)
  //  - in_argc and out_argc describe the number of input and output
  //    streams.
  //  - param_type (inherited) defaults to 'Empty_params'.
public:
  static unsigned int const pre_argc = 1;
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;

  typedef Uk_vmmul_params param_type;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Vmmul_kernel(unsigned int size)
    : size_ (size)
    , pre_sp(vsip::impl::ukernel::Whole_sdist())
    , io_sp (vsip::impl::ukernel::Blocksize_sdist(1), vsip::impl::ukernel::Whole_sdist())
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


  // Queury API:
  // - fill_params() fills the parameter block to be passed to the
  //   accelerators.
  // - in_spatt()/out_spatt() allow VSIPL++ to determine streaming
  //   pattern for user-kernel.  Since both input and output have same
  //   streaming pattern, simply return 'sp'

  void fill_params(param_type& param) const
  {
    param.size  = size_;
  }

  vsip::impl::ukernel::Stream_pattern const& in_spatt(vsip::index_type i) const
  { return (i == 0) ? pre_sp : io_sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return io_sp; }

  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  unsigned int               size_;
  vsip::impl::ukernel::Stream_pattern pre_sp;	
  vsip::impl::ukernel::Stream_pattern io_sp;	
};

DEFINE_UKERNEL_TASK(
  Vmmul_kernel,
  void(std::complex<float>*, std::complex<float>*, std::complex<float>*),
  "uk_plugin", cvmmul_f)

DEFINE_UKERNEL_TASK(
  Vmmul_kernel,
  void(std::pair<float*,float*>, std::pair<float*,float*>,
       std::pair<float*,float*>),
  "uk_plugin", zvmmul_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_VMMUL_HPP
