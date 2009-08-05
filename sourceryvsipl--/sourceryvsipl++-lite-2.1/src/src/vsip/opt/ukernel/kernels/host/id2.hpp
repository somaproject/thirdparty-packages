/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    src/vsip/opt/ukernel/kernels/host/id2.hpp
    @author  Jules Bergmann
    @date    2008-07-29
    @brief   VSIPL++ Library: ID2 Ukernel
*/

#ifndef VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID2_HPP
#define VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID2_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/id2_param.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

vsip::impl::ukernel::Stream_pattern
spatt_shape(int shape)
{
  using vsip::impl::ukernel::Stream_pattern;
  using vsip::impl::ukernel::Blocksize_sdist;
  switch (shape)
  {
  case 1:
    return Stream_pattern(Blocksize_sdist(32),
		          Blocksize_sdist(32));
  case 2:
    return Stream_pattern(Blocksize_sdist(2),
			  Blocksize_sdist(32));
  case 3:
    return Stream_pattern(Blocksize_sdist(2),
			  Blocksize_sdist(16));
  case 4:
    return Stream_pattern(Blocksize_sdist(4),
			  Blocksize_sdist(32));
  default:
    return Stream_pattern(Blocksize_sdist(1),
		          Blocksize_sdist(1024));
  }
}


// Host-side vector elementwise copy ukernel.

class Id2_ukernel : public vsip::impl::ukernel::Host_kernel_base
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
  typedef Uk_id2_params param_type;

  // Host-side ukernel object initialization.
  //
  // Streaming pattern divides linear vector into blocks, with minimum
  // size 256, maximum size 1024.  (Blocksize_sdist has an implicit
  // quantum of 16 elements).

  Id2_ukernel(int shape, unsigned int rows, unsigned int cols)
    : sp    (spatt_shape(shape)),
      rows_ (rows),
      cols_ (cols)
  {}

  void fill_params(param_type& param) const
  {
    param.rows  = rows_;
    param.cols  = cols_;
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
  { return sp; }

  vsip::impl::ukernel::Stream_pattern const& out_spatt(vsip::index_type) const
  { return sp; }

  // Member data.
  //
  // 'sp' is the stream pattern.
private:
  vsip::impl::ukernel::Stream_pattern sp;	

  unsigned int rows_;
  unsigned int cols_;
};

DEFINE_UKERNEL_TASK(Id2_ukernel, void(float*, float*), "uk_plugin", id2_f)

#endif // VSIP_SRC_OPT_UKERNEL_KERNELS_HOST_ID2_HPP
