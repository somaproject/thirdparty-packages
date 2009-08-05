/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip_csl/img/separable_filter.hpp
    @author  Jules Bergmann
    @date    2007-10-04
    @brief   VSIPL++ Library: Image-processing separable filter.

*/

#ifndef VSIP_CSL_IMG_SEPARABLE_FILTER_HPP
#define VSIP_CSL_IMG_SEPARABLE_FILTER_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/impl_tags.hpp>

#include <vsip_csl/img/impl/sfilt_gen.hpp>
#if VSIP_IMPL_HAVE_IPP
#  include <vsip_csl/img/impl/sfilt_ipp.hpp>
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

template <vsip::dimension_type Dim,
	  typename             T>
struct Choose_sfilt_impl
{
  typedef vsip::impl::Intel_ipp_tag   Intel_ipp_tag;
  typedef vsip::impl::Mercury_sal_tag Mercury_sal_tag;
  typedef vsip::impl::Generic_tag     Generic_tag;

  typedef typename
  vsip::impl::ITE_Type<Is_sfilt_impl_avail<Intel_ipp_tag, Dim, T>::value,
		       vsip::impl::As_type<Intel_ipp_tag>, 
  vsip::impl::ITE_Type<Is_sfilt_impl_avail<Mercury_sal_tag, Dim, T>::value,
		       vsip::impl::As_type<Mercury_sal_tag>, 
		       vsip::impl::As_type<Generic_tag> > >::type type;
};

} // namespace vsip_csl::img::impl

template <typename                  T,
	  vsip::support_region_type SuppT,
	  edge_handling_type        EdgeT,
	  unsigned                  N_times = 0,
	  vsip::alg_hint_type       A_hint = vsip::alg_time>
class Separable_filter
  : public impl::Sfilt_impl<T, SuppT, EdgeT, N_times, A_hint,
			    typename impl::Choose_sfilt_impl<2, T>::type>
{
  // Compile-time values and typedefs.
public:
  typedef typename impl::Choose_sfilt_impl<2, T>::type impl_tag;

private:
  typedef impl::Sfilt_impl<T, SuppT, EdgeT, N_times, A_hint, impl_tag>
		base_type;
  static vsip::dimension_type const dim = 2;

// Constructors, copies, assignments, and destructor.
public:
  template <typename Block1,
            typename Block2>
  Separable_filter(
    vsip::Vector<T, Block1> coeff0,	// coeffs for dimension 0
    vsip::Vector<T, Block2> coeff1,	// coeffs for dimension 1
    vsip::Domain<2> const&  input_size)
  VSIP_THROW((std::bad_alloc))
    : base_type(coeff0, coeff1, input_size)
    {}

  Separable_filter(Separable_filter const&) VSIP_NOTHROW;
  Separable_filter& operator=(Separable_filter const&) VSIP_NOTHROW;
  ~Separable_filter() VSIP_NOTHROW {}

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
  vsip::Domain<dim> const& kernel_size() const VSIP_NOTHROW;
  vsip::Domain<dim> const& filter_order() const VSIP_NOTHROW;
  vsip::Domain<dim> const& input_size() const VSIP_NOTHROW;
  vsip::Domain<dim> const& output_size() const VSIP_NOTHROW;
};

} // namespace vsip_csl::img

} // namespace vsip_csl

#endif // VSIP_CSL_IMG_SEPARABLE_FILTER_HPP
