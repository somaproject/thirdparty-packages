/* Copyright (c) 2006, 2007, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/expr_evaluator.hpp
    @author  Stefan Seefeld
    @date    2006-07-25
    @brief   VSIPL++ Library: SIMD expression evaluator proxy factory.

*/

#ifndef VSIP_IMPL_SIMD_PROXY_FACTORY_HPP
#define VSIP_IMPL_SIMD_PROXY_FACTORY_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/opt/simd/simd.hpp>
#include <vsip/opt/simd/expr_iterator.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>

/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace simd
{

template <typename BlockT, bool A>
struct Proxy_factory
{
  typedef Direct_access_traits<typename BlockT::value_type> access_traits;
  typedef Proxy<access_traits, A> proxy_type;
  typedef typename Adjust_layout_dim<
                     1, typename Block_layout<BlockT>::layout_type>::type
		layout_type;

  static bool const ct_valid = Ext_data_cost<BlockT>::value == 0 &&
    !Is_split_block<BlockT>::value;

  static bool 
  rt_valid(BlockT const &b, int alignment)
  {
    Ext_data<BlockT, layout_type> dda(b, SYNC_IN);
    return dda.stride(0) == 1 && 
      (!A ||
       Simd_traits<typename BlockT::value_type>::alignment_of(dda.data()) ==
       alignment);
  }

  static int
  alignment(BlockT const &b)
  {
    Ext_data<BlockT, layout_type> dda(b, SYNC_IN);
    return Simd_traits<typename BlockT::value_type>::alignment_of(dda.data());
  }

  static proxy_type
  create(BlockT const &b) 
  {
    Ext_data<BlockT, layout_type> dda(b, SYNC_IN);
    return proxy_type(dda.data());
  }
};

template <typename T, bool A>
struct Proxy_factory<Scalar_block<1, T>, A>
{
  typedef Scalar_access_traits<T> access_traits;
  typedef Proxy<access_traits, A> proxy_type;
  static bool const ct_valid = true;

  static bool 
  rt_valid(Scalar_block<1, T> const &, int) {return true;}

  static proxy_type
  create(Scalar_block<1, T> const &b) 
  {
    return proxy_type(b.value());
  }
};

template <dimension_type D,
	  template <typename> class O,
	  typename B,
	  typename T,
	  bool A>
struct Proxy_factory<Unary_expr_block<D, O, B, T> const, A>
{
  typedef 
    Unary_access_traits<typename Proxy_factory<B,A>::proxy_type, O>
    access_traits;
  typedef Proxy<access_traits,A> proxy_type;

  static bool const ct_valid =
    Unary_operator_map<T, O>::is_supported &&
    Type_equal<typename B::value_type, T>::value &&
    Proxy_factory<B, A>::ct_valid;

  static bool 
  rt_valid(Unary_expr_block<D, O, B, T> const &b, int alignment)
  {
    return Proxy_factory<B, A>::rt_valid(b.op(), alignment);
  }

  static proxy_type
  create(Unary_expr_block<D, O, B, T> const &b)
  {
    return proxy_type(Proxy_factory<B, A>::create(b.op()));
  }
};

// This proxy is specialized for unaligned blocks. If the user specifies
// ualigned(block), this is a hint to switch to an unaligned proxy.
template <dimension_type D,
	  typename B,
	  typename T,
	  bool A>
struct Proxy_factory<Unary_expr_block<D, unaligned_functor, B, T> const, A>
{
  typedef typename Proxy_factory<B, false>::access_traits access_traits;
  typedef Proxy<access_traits,false> proxy_type;
  static bool const ct_valid = Proxy_factory<B,false>::ct_valid;


  static bool 
  rt_valid(Unary_expr_block<D, unaligned_functor, B, T> const &b, int alignment)
  {
    return Proxy_factory<B, false>::rt_valid(b.op(), alignment);
  }

  static proxy_type
  create(Unary_expr_block<D, unaligned_functor, B, T> const &b)
  {
    return proxy_type(Proxy_factory<B, false>::create(b.op()));
  }
};

template <dimension_type                D,
	  template <typename, typename> class O,
	  typename                      LB,
	  typename                      LT,
	  typename                      RB,
	  typename                      RT,
	  bool A>
struct Proxy_factory<Binary_expr_block<D, O, LB, LT, RB, RT> const, A>
{
  typedef
    Binary_access_traits<typename Proxy_factory<LB, A>::proxy_type,
			 typename Proxy_factory<RB, A>::proxy_type, O> 
    access_traits;
  typedef Proxy<access_traits, A> proxy_type;
  static bool const ct_valid = 
    Type_equal<typename LB::value_type, LT>::value &&
    Type_equal<typename RB::value_type, RT>::value &&
    Type_equal<LT, RT>::value &&
    Binary_operator_map<LT, O>::is_supported &&
    Proxy_factory<LB, A>::ct_valid &&
    Proxy_factory<RB, A>::ct_valid;

  static bool 
  rt_valid(Binary_expr_block<D, O, LB, LT, RB, RT> const &b, int alignment)
  {
    return Proxy_factory<LB, A>::rt_valid(b.left(), alignment) &&
           Proxy_factory<RB, A>::rt_valid(b.right(), alignment);
  }

  static proxy_type
  create(Binary_expr_block<D, O, LB, LT, RB, RT> const &b)
  {
    typename Proxy_factory<LB, A>::proxy_type lp =
      Proxy_factory<LB, A>::create(b.left());
    typename Proxy_factory<RB, A>::proxy_type rp =
      Proxy_factory<RB, A>::create(b.right());

    return proxy_type(lp, rp);
  }
};

template <dimension_type                         D,
	  template <typename, typename,typename> class O,
	  typename                               Block1, typename Type1,
	  typename                               Block2, typename Type2,
	  typename                               Block3, typename Type3,
	  bool A>
struct Proxy_factory<Ternary_expr_block<D, O,
  Block1,Type1,Block2,Type2,Block3,Type3> const, A>
{
  typedef Ternary_access_traits<typename Proxy_factory<Block1, A>::proxy_type,
                                typename Proxy_factory<Block2, A>::proxy_type,
                                typename Proxy_factory<Block3, A>::proxy_type,
		 	        O> 
    access_traits;

  typedef Ternary_expr_block<D, O, Block1,Type1,Block2,Type2,Block3,Type3>
    SrcBlock;

  typedef Proxy<access_traits, A> proxy_type;
  static bool const ct_valid = 
    Ternary_operator_map<Type1, O>::is_supported &&
    Proxy_factory<Block1, A>::ct_valid &&
    Proxy_factory<Block2, A>::ct_valid &&
    Proxy_factory<Block3, A>::ct_valid;

  static bool 
  rt_valid(SrcBlock const &b, int alignment)
  {
    return Proxy_factory<Block1, A>::rt_valid(b.first(), alignment) &&
           Proxy_factory<Block2, A>::rt_valid(b.second(), alignment) &&
           Proxy_factory<Block3, A>::rt_valid(b.third(), alignment);
  }

  static proxy_type
  create(SrcBlock const &b)
  {
    typename Proxy_factory<Block1, A>::proxy_type
      b1p = Proxy_factory<Block1, A>::create(b.first());
    typename Proxy_factory<Block2, A>::proxy_type
      b2p = Proxy_factory<Block2, A>::create(b.second());
    typename Proxy_factory<Block3, A>::proxy_type
      b3p = Proxy_factory<Block3, A>::create(b.third());

    return proxy_type(b1p,b2p,b3p);
  }
};


} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_SIMD_PROXY_FACTORY_HPP
