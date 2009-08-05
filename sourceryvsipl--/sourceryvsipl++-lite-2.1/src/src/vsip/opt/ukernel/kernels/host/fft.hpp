/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/fft.hpp
    @author  Jules Bergmann
    @date    2008-06-18
    @brief   VSIPL++ Library: FFT Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_FFT_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_FFT_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/fft_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side FFT ukernel.

class Fft_ukernel : public vsip::impl::ukernel::Host_kernel_base
{
  // Parameters.
  //  - 'tag_type' is used to select the appropriate kernel (via
  //    Task_manager's Task_map)
  //  - in_argc and out_argc describe the number of input and output
  //    streams.
  //  - param_type describes the parameters sent to the accelerator.
public:
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;
  typedef Uk_fft_params param_type;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides matrix into blocks of 1 row, keeping
  // all elements in a row together.
public:
  Fft_ukernel(unsigned int size, int dir, float scale)
    : size_ (size)
    , dir_  (dir)
    , scale_(scale)
    , sp    (vsip::impl::ukernel::Blocksize_sdist(1), vsip::impl::ukernel::Whole_sdist())
  {}

  // Host-side compute kernel.  Used if accelerator is not available.

  template <typename View1,
	    typename View2>
  void compute(
    View1 in,
    View2 out)
  {
    assert(0); // TODO
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
    param.dir   = dir_;
    param.scale = scale_;
  }

  vsip::impl::ukernel::Stream_pattern const& in_spatt(vsip::index_type) const
  { return sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return sp; }

  vsip::length_type stack_size() const { return 4096; }

  // Member data.
private:
  unsigned int size_;
  int          dir_;
  float        scale_;

  vsip::impl::ukernel::Stream_pattern sp;	
};

DEFINE_UKERNEL_TASK(Fft_ukernel,
		    void(std::complex<float>*, std::complex<float>*),
		    "uk_plugin", cfft_f)
DEFINE_UKERNEL_TASK(Fft_ukernel,
		    void(split_float_type, split_float_type),
		    "uk_plugin", zfft_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_FFT_HPP
