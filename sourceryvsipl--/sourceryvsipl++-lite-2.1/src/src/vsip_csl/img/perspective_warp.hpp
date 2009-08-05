/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip_csl/img/perspective_warp.hpp
    @author  Jules Bergmann
    @date    2007-11-01
    @brief   VSIPL++ Library: Image-processing perspective warp.

*/

#ifndef VSIP_CSL_IMG_PERSPECTIVE_WARP_HPP
#define VSIP_CSL_IMG_PERSPECTIVE_WARP_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/matrix.hpp>
#include <vsip/core/impl_tags.hpp>
#include <vsip_csl/img/impl/pwarp_common.hpp>
#include <vsip_csl/img/impl/pwarp_gen.hpp>
#include <vsip/opt/simd/simd.hpp>
#include <vsip_csl/img/impl/pwarp_simd.hpp>
#ifdef VSIP_IMPL_CBE_SDK
#  include <vsip_csl/img/impl/pwarp_cbe.hpp>
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip_csl
{

namespace img
{

namespace impl
{

// Implemtation chooser for perspective warp processing object.

struct Choose_pwarp_impl
{
  
  template <typename         CoeffT,
	    typename         T,
	    interpolate_type InterpT,
	    transform_dir    T_dir>
  struct choose_impl
  {
    typedef vsip::impl::Intel_ipp_tag    Intel_ipp_tag;
    typedef vsip::impl::Mercury_sal_tag  Mercury_sal_tag;
    typedef vsip::impl::Simd_builtin_tag Simd_builtin_tag;
    typedef vsip::impl::Cbe_sdk_tag      Cbe_sdk_tag;
    typedef vsip::impl::Generic_tag      Generic_tag;

    typedef typename
    vsip::impl::ITE_Type<
      Is_pwarp_impl_avail<Intel_ipp_tag, CoeffT, T, InterpT, T_dir>::value,
		         vsip::impl::As_type<Intel_ipp_tag>, 
    vsip::impl::ITE_Type<
        Is_pwarp_impl_avail<Mercury_sal_tag, CoeffT, T, InterpT, T_dir>::value,
        vsip::impl::As_type<Mercury_sal_tag>, 
    vsip::impl::ITE_Type<
        Is_pwarp_impl_avail<Cbe_sdk_tag, CoeffT, T, InterpT, T_dir>::value,
        vsip::impl::As_type<Cbe_sdk_tag>, 
    vsip::impl::ITE_Type<
        Is_pwarp_impl_avail<Simd_builtin_tag, CoeffT, T, InterpT, T_dir>::value,
        vsip::impl::As_type<Simd_builtin_tag>, 
        vsip::impl::As_type<Generic_tag> > > > >::type type;
  };
};

} // namespace vsip_csl::img::impl



// Perspective warp image processing object.

template <typename            CoeffT,
	  typename            T,
	  interpolate_type    InterpT,
	  transform_dir       T_dir,
	  unsigned            N_times = 0,
	  vsip::alg_hint_type A_hint = vsip::alg_time,
	  typename            ChooserT = impl::Choose_pwarp_impl>
class Perspective_warp
  : public impl::Pwarp_impl<CoeffT, T, InterpT, T_dir, N_times, A_hint,
           typename ChooserT::template
                    choose_impl<CoeffT, T, InterpT, T_dir>::type>
{
// Compile-time values and types.
public:
  typedef typename ChooserT::template
                    choose_impl<CoeffT, T, InterpT, T_dir>::type
		impl_tag;
  typedef impl::Pwarp_impl<CoeffT, T, InterpT, T_dir, N_times, A_hint,
			   impl_tag>
		base_type;
  static vsip::dimension_type const dim = 2;

// Constructors, copies, assignments, and destructor.
public:
  template <typename Block1>
  Perspective_warp(
    vsip::const_Matrix<CoeffT, Block1> coeff,
    vsip::Domain<2> const&             size)
  VSIP_THROW((std::bad_alloc))
    : base_type(coeff, size)
    {}

  Perspective_warp(Perspective_warp const&) VSIP_NOTHROW;
  Perspective_warp& operator=(Perspective_warp const&) VSIP_NOTHROW;
  ~Perspective_warp() VSIP_NOTHROW {}

// Operator
public:
  template <typename Block1,
            typename Block2>
  vsip::Matrix<T, Block2>
  operator()(
    vsip::const_Matrix<T, Block1> in,
    vsip::Matrix<T, Block2>       out)
    VSIP_NOTHROW
  {
    filter(in, out);
    return out;
  }

// Accessors
public:
  vsip::Domain<dim> const& input_size()  const VSIP_NOTHROW;
  vsip::Domain<dim> const& output_size() const VSIP_NOTHROW;
};



// Perspective warp image processing utility function.

template <typename CoeffT,
	  typename T,
	  typename Block1,
	  typename Block2,
	  typename Block3>
void
perspective_warp(
  vsip::const_Matrix<CoeffT, Block1> P,
  vsip::const_Matrix<T, Block2>      in,
  vsip::Matrix<T, Block3>            out)
{
  typedef Perspective_warp<CoeffT, T, interp_linear, forward, 1,
                           vsip::alg_time>
    pwarp_type;

  pwarp_type pwarp(P, vsip::Domain<2>(in.size(0), in.size(1)));
  pwarp(in, out);
}

} // namespace vsip_csl::img

} // namespace vsip_csl

#endif // VSIP_CSL_IMG_PERSPECTIVE_WARP_HPP
