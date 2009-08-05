/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/interp.hpp
    @author  Don McCoy
    @date    2008-08-26
    @brief   VSIPL++ Library: User-defined polar to rectangular
               interpolation kernel for SSAR images.
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_INTERP_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_INTERP_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Interp_kernel : public vsip::impl::ukernel::Host_kernel_base
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
  // Streaming pattern divides matrix into whole, single columns (for now).
  //

  Interp_kernel()
    : sp (vsip::impl::ukernel::Blocksize_sdist(1), 
          vsip::impl::ukernel::Whole_sdist())
  {}



  // Host-side compute kernel.  Used if accelerator is not available.
  //
  // View sizes:
  //   in0 is N x M
  //   in1 is N x M x I
  //   in2 is N x M
  //   out is NX x M    
  // 
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
    out = View3::value_type();

    for (vsip::index_type j = 0; j < in0.size(1); ++j)
      for (vsip::index_type i = 0; i < in0.size(0); ++i)
        for (vsip::index_type h = 0; h < in1.size(2); ++h)
        {
          vsip::index_type ikxrows = in0.get(i, j) + h;

          out.put(ikxrows, j, out.get(ikxrows, j) + 
            (in2.get(i, j) * in1.get(i, j, h)));
        }
  }


  // Queury API:
  // - in_spatt()/out_spatt() allow VSIPL++ to determine streaming
  //   pattern for user-kernel.  Since both input and output have same
  //   streaming pattern, simply return 'sp'

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


DEFINE_UKERNEL_TASK(
  Interp_kernel,
  void(uint32_t*, float*, std::complex<float>*, std::complex<float>*),
  "uk_plugin", interp_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_INTERP_HPP
