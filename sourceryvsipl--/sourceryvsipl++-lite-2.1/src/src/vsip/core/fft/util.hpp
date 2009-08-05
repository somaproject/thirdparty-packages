/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fft/util.cpp
    @author  Stefan Seefeld
    @date    2006-02-21
    @brief   VSIPL++ Library: FFT common infrastructure used by all 
    implementations.
*/

#ifndef VSIP_CORE_FFT_UTIL_HPP
#define VSIP_CORE_FFT_UTIL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/fft/backend.hpp>
#include <vsip/core/fast_block.hpp>
#include <vsip/core/view_traits.hpp>
#ifndef VSIP_IMPL_REF_IMPL
#  include <vsip/opt/fft/return_functor.hpp>
#  include <vsip/opt/expr/return_block.hpp>
#endif

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace fft
{

/// Determine whether the FFT size is a power of two.
inline bool 
is_power_of_two(unsigned size)
{
  return (size & (size - 1)) == 0;
}
template <dimension_type D>
inline bool 
is_power_of_two(Domain<D> const &dom)
{
  for (dimension_type d = 0; d != D; ++d)
    if (!is_power_of_two(dom[d].size())) return false;
  return true;
}

/// Determine the exponent (forward or inverse) of a given Fft
/// from its parameters.
template <typename I, typename O, int sD> struct exponent;
template <typename T, int sD>
struct exponent<T, std::complex<T>, sD> { static int const value = -1;};
template <typename T, int sD>
struct exponent<std::complex<T>, T, sD> { static int const value = 1;};
template <typename T>
struct exponent<T, T, -2> { static int const value = -1;};
template <typename T>
struct exponent<T, T, -1> { static int const value = 1;};

/// Determine the 'special dimension' of a real fft
/// from its parameters.
template <typename I, typename O, int S> 
struct axis { static int const value = S;};
template <typename T, int S>
struct axis<T, T, S> { static int const value = 0;};

/// Device to calculate the size of the input and output blocks
/// for complex and real ffts.
template <dimension_type D, typename I, typename O, int A>
struct io_size
{
  static Domain<D> size(Domain<D> const &dom) { return dom;}
};
template <dimension_type D, typename T, int A>
struct io_size<D, std::complex<T>, T, A>
{
  static Domain<D> size(Domain<D> const &dom) 
  {
    Domain<D> retn(dom);
    Domain<1> &mod = retn.impl_at(A);
    mod = Domain<1>(0, 1, mod.size() / 2 + 1); 
    return retn;
  }
};


template <typename View>
View
new_view(vsip::Domain<1> const& dom) { return View(dom.size());} 

template <typename View>
View 
new_view(vsip::Domain<2> const& dom)
{ return View(dom[0].size(), dom[1].size());}

template <typename View>
View  
new_view(vsip::Domain<3> const& dom)
{ return View(dom[0].size(), dom[1].size(), dom[2].size());}



template <typename View>
View
new_view(
  vsip::Domain<1> const&                     dom,
  typename View::block_type::map_type const& map)
{ return View(dom.size(), map);} 

template <typename View>
View 
new_view(
  vsip::Domain<2> const&                     dom,
  typename View::block_type::map_type const& map)
{ return View(dom[0].size(), dom[1].size(), map);}

template <typename View>
View  
new_view(
  vsip::Domain<3> const&                     dom,
  typename View::block_type::map_type const& map)
{ return View(dom[0].size(), dom[1].size(), dom[2].size(), map);}



/// Traits class to determine block type returned by Fft and Fftm
/// by_value operators.

/// General case: when result is distributed, we use a Dense block
/// because Fast_blocks do support non-local maps (060512).
template<typename T,
	 typename BlockT,
	 typename MapT = typename BlockT::map_type>
struct result
{
  static dimension_type const dim = BlockT::dim;
  typedef Dense<dim, T, tuple<0,1,2>, MapT> block_type;

  typedef typename View_of_dim<dim, T, block_type>::type view_type;

  static view_type create(Domain<dim> const &dom, MapT const& map)
  { return new_view<view_type>(dom, map);}
};

/// Specialization: when result is local, we use a Fast_block,
/// because it lets us match complex format of the input block.
template<typename T, typename BlockT>
struct result<T, BlockT, Local_map>
{
  static dimension_type const dim = BlockT::dim;
  typedef typename
  impl::Fast_block<dim, T, 
		   Layout<dim, tuple<0,1,2>,
			  impl::Stride_unit_dense,
			  typename Block_layout<BlockT>::complex_type>,
		   typename BlockT::map_type> block_type;

  typedef typename View_of_dim<dim, T, block_type>::type view_type;

  static view_type create(Domain<dim> const &dom, Local_map const&)
  { return new_view<view_type>(dom);}
};



#ifndef VSIP_IMPL_REF_IMPL

/// Traits class to determine view type returned by Fft for
/// by_value operators with return-block optimization.

template <dimension_type Dim,
	  typename       InT,
	  typename       OutT,
	  typename       ViewT,
	  typename       WorkspaceT,
	  int            AxisV,
	  int            ExponentV>
struct Result_rbo
{
  typedef Fft_return_functor<Dim, OutT, typename ViewT::block_type,
                      fft::backend<Dim, InT, OutT, AxisV, ExponentV>,
		      WorkspaceT>
		functor_type;

  typedef Return_expr_block<Dim, OutT, functor_type>
		block_type;

  typedef typename View_of_dim<Dim, OutT, block_type const>::type
		view_type;
};



template <typename       InT,
	  typename       OutT,
	  typename       BlockT,
	  typename       WorkspaceT,
	  int            AxisV,
	  int            ExponentV>
struct Result_fftm_rbo
{
  static dimension_type const dim = 2;
  typedef const_Matrix<InT, BlockT> in_view_type;
  typedef Fft_return_functor<dim, OutT, BlockT,
                      fft::fftm<InT, OutT, AxisV, ExponentV>,
		      WorkspaceT>
		functor_type;

  typedef Return_expr_block<dim, OutT, functor_type>
		block_type;

  typedef typename View_of_dim<dim, OutT, block_type const>::type
		view_type;
};

#endif

} // namespace vsip::impl::fft
} // namespace vsip::impl
} // namespace vsip

#endif
