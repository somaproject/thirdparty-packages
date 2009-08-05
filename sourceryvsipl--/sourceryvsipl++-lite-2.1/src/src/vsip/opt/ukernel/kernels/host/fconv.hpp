/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/host_uk_fconv.hpp
    @author  Jules Bergmann
    @date    2008-06-27
    @brief   VSIPL++ Library: Vector-matrix element-wise multiply Ukernel
*/

#ifndef VSIP_TESTS_HOST_UK_FCONV_HPP
#define VSIP_TESTS_HOST_UK_FCONV_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/fused_param.hpp>
#include <vsip/opt/ukernel/kernels/params/fft_param.hpp>
#include <vsip/opt/ukernel/kernels/params/vmmul_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

// Host-side vector elementwise copy ukernel.

class Fconv_kernel : public vsip::impl::ukernel::Host_kernel_base
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

  typedef Uk_fused_params<Uk_fft_params,
			  Uk_vmmul_params,
			  Uk_fft_params>
		param_type;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Fconv_kernel(unsigned int size)
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
    param.k1_params.size  = size_;
    param.k1_params.dir   = -1;
    param.k1_params.scale = 1.f;
    param.k2_params.size  = size_;
    param.k3_params.size  = size_;
    param.k3_params.dir   = +1;
    param.k3_params.scale = 1.f / size_;
  }

  vsip::impl::ukernel::Stream_pattern const& in_spatt(vsip::index_type i) const
  { return (i == 0) ? pre_sp : io_sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return io_sp; }

  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  unsigned int                        size_;
  vsip::impl::ukernel::Stream_pattern pre_sp;	
  vsip::impl::ukernel::Stream_pattern io_sp;	
};

DEFINE_UKERNEL_TASK(Fconv_kernel,
		    void(std::complex<float>*, std::complex<float>*,
			 std::complex<float>*),
		    "uk_plugin", cfconv_f)

DEFINE_UKERNEL_TASK(Fconv_kernel,
		    void(std::pair<float*,float*>, std::pair<float*,float*>,
			 std::pair<float*,float*>),
		    "uk_plugin", zfconv_f)

#endif // VSIP_TESTS_HOST_UK_FCONV_HPP
