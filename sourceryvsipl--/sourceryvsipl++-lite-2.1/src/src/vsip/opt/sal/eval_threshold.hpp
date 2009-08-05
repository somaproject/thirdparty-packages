/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/eval_threshold.hpp
    @author  Jules Bergmann
    @date    2006-10-18
    @brief   VSIPL++ Library: Dispatch for Mercury SAL -- threshold.
*/

#ifndef VSIP_OPT_SAL_EVAL_THRESHOLD_HPP
#define VSIP_OPT_SAL_EVAL_THRESHOLD_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/expr/ternary_block.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/opt/sal/eval_util.hpp>
#include <vsip/core/adjust_layout.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

/***********************************************************************
  Threshold Expressions
***********************************************************************/

// Optimize threshold expression: ite(A > b, A, c)
#define VSIP_IMPL_SAL_GTE_THRESH_EXPR(T0, FUN, FUN0)			\
template <typename DstBlock,						\
	  typename T,							\
	  typename Block1>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
									\
         Ternary_expr_block<1, ite_functor,				\
           Binary_expr_block<1u, ge_functor,				\
			     Block1, T,					\
			     Scalar_block<1, T>, T> const, bool,	\
	   Block1, T,							\
	   Scalar_block<1, T>, T> const,				\
									\
         Mercury_sal_tag>						\
{									\
  static char const* name() { return "Expr_SAL_thresh"; }		\
									\
  typedef Ternary_expr_block<1, ite_functor,				\
            Binary_expr_block<1u, ge_functor,				\
			      Block1, T,				\
			      Scalar_block<1, T>, T> const, bool,	\
	   Block1, T,							\
	   Scalar_block<1, T>, T>					\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, T>::type eff_1_t;	\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block1>::layout_type>::type		\
    block1_lp;								\
  									\
  static bool const ct_valid =						\
     Type_equal<T, T0>::value &&					\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<Block1>::value == 0;					\
									\
  static bool rt_valid(DstBlock&, SrcBlock const& src)			\
  {									\
    return &(src.first().left()) == &(src.second()) &&			\
           (src.first().right().value() == src.third().value() ||	\
            src.third().value() == T(0));				\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    typedef Scalar_block<1, T> sb_type;					\
									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,        SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp> ext_A(src.second(), SYNC_IN);	\
    sal::Ext_wrapper<sb_type>           ext_b(src.first().right(), SYNC_IN);\
									\
    if (src.third().value() == T(0))					\
      FUN0(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_A),\
           typename sal::Ext_wrapper<sb_type>::sal_type(ext_b),		\
	   typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),\
	   dst.size());							\
    else								\
      FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_A),\
           typename sal::Ext_wrapper<sb_type>::sal_type(ext_b),		\
	   typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),\
	   dst.size());							\
  }									\
};

// This definition (and the above ifdef) will go away once the code
// has benchmarked.
// VSIP_IMPL_SAL_GTE_THRESH_EXPR(float, sal::vthresh, sal::vthresh0)



// Common evaluator for a threshold expressions.

template <typename DstBlock,
	  typename T,
	  typename Block1>
struct Thresh_expr_evaluator
{
  static char const* name() { return "Expr_SAL_thresh"; }

  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;
  typedef typename sal::Effective_value_type<Block1, T>::type eff_1_t;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;
  
  typedef typename Adjust_layout_dim<
      1, typename Block_layout<Block1>::layout_type>::type
    block1_lp;

  static bool const ct_valid =
     Type_equal<T, float>::value &&
     /* check that direct access is supported */
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<Block1>::value == 0;

  static bool rt_valid(DstBlock&,
		       Block1 const&             a1,
		       Scalar_block<1, T> const& b,
		       Block1 const&             a2,
		       Scalar_block<1, T> const& c)
  {
    return &a1 == &a2 && (b.value() == c.value() || c.value() == T(0));
  }

  static void exec(
    DstBlock&                 dst,
    Block1 const&             a,
    Scalar_block<1, T> const& b,
    Block1 const&,
    Scalar_block<1, T> const& c)
  {
    typedef Scalar_block<1, T> sb_type;

    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst, SYNC_OUT);
    sal::Ext_wrapper<Block1, block1_lp> ext_A(a,     SYNC_IN);
    sal::Ext_wrapper<sb_type>           ext_b(b,     SYNC_IN);

    if (c.value() == T(0))
      sal::vthresh0(
	typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_A),
        typename sal::Ext_wrapper<sb_type>::sal_type(ext_b),
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),
	dst.size());
    else
      sal::vthresh(
	typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_A),
	typename sal::Ext_wrapper<sb_type>::sal_type(ext_b),
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),
	dst.size());
  }
};



/// Frontend for threshold expressions like:
///
///   ite(A >= b, A, c)

template <typename DstBlock,
	  typename T,
	  typename Block1>
struct Serial_expr_evaluator<
         1, DstBlock, 

         Ternary_expr_block<1, ite_functor,
           Binary_expr_block<1u, ge_functor,
			     Block1, T,
			     Scalar_block<1, T>, T> const, bool,
	   Block1, T,
	   Scalar_block<1, T>, T> const,

         Mercury_sal_tag>
  : Thresh_expr_evaluator<DstBlock, T, Block1>
{
  typedef Thresh_expr_evaluator<DstBlock, T, Block1> base_type;

  typedef Ternary_expr_block<1, ite_functor,
            Binary_expr_block<1u, ge_functor,
			      Block1, T,
			      Scalar_block<1, T>, T> const, bool,
	   Block1, T,
	   Scalar_block<1, T>, T>
	SrcBlock;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    return base_type::rt_valid(dst,
			src.first().left(),
			src.first().right(),
			src.second(),
			src.third());
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    base_type::exec(dst,
		    src.first().left(),
		    src.first().right(),
		    src.second(),
		    src.third());
  }
};



/// Frontend for threshold expressions like:
///
///   ite(A < b, c, A)

template <typename DstBlock,
	  typename T,
	  typename Block1>
struct Serial_expr_evaluator<
         1, DstBlock, 

         Ternary_expr_block<1, ite_functor,
           Binary_expr_block<1u, lt_functor,
			     Block1, T,
			     Scalar_block<1, T>, T> const, bool,
	   Scalar_block<1, T>, T,
	   Block1, T> const,

         Mercury_sal_tag>
  : Thresh_expr_evaluator<DstBlock, T, Block1>
{
  typedef Thresh_expr_evaluator<DstBlock, T, Block1> base_type;

  typedef Ternary_expr_block<1, ite_functor,
            Binary_expr_block<1u, lt_functor,
			      Block1, T,
			      Scalar_block<1, T>, T> const, bool,
	   Scalar_block<1, T>, T,
	   Block1, T>
	SrcBlock;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    return base_type::rt_valid(dst,
			src.first().left(),
			src.first().right(),
			src.third(),
			src.second());
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    base_type::exec(dst,
		    src.first().left(),
		    src.first().right(),
		    src.third(),
		    src.second());
  }
};



/// Frontend for threshold expressions like:
///
///   ite(b <= A, A, c)

template <typename DstBlock,
	  typename T,
	  typename Block1>
struct Serial_expr_evaluator<
         1, DstBlock, 

         Ternary_expr_block<1, ite_functor,
           Binary_expr_block<1u, le_functor,
			     Scalar_block<1, T>, T,
			     Block1, T> const, bool,
	   Block1, T,
	   Scalar_block<1, T>, T> const,

         Mercury_sal_tag>
  : Thresh_expr_evaluator<DstBlock, T, Block1>
{
  typedef Thresh_expr_evaluator<DstBlock, T, Block1> base_type;

  typedef Ternary_expr_block<1, ite_functor,
            Binary_expr_block<1u, le_functor,
			      Scalar_block<1, T>, T,
			      Block1, T> const, bool,
	   Block1, T,
	   Scalar_block<1, T>, T>
	SrcBlock;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    return base_type::rt_valid(dst,
			src.first().right(),
			src.first().left(),
			src.second(),
			src.third());
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    base_type::exec(dst,
		    src.first().right(),
		    src.first().left(),
		    src.second(),
		    src.third());
  }
};



/// Frontend for threshold expressions like:
///
///   ite(b > A, c, A)

template <typename DstBlock,
	  typename T,
	  typename Block1>
struct Serial_expr_evaluator<
         1, DstBlock, 

         Ternary_expr_block<1, ite_functor,
           Binary_expr_block<1u, gt_functor,
			     Scalar_block<1, T>, T,
			     Block1, T> const, bool,
	   Scalar_block<1, T>, T,
	   Block1, T> const,

         Mercury_sal_tag>
  : Thresh_expr_evaluator<DstBlock, T, Block1>
{
  typedef Thresh_expr_evaluator<DstBlock, T, Block1> base_type;

  typedef Ternary_expr_block<1, ite_functor,
            Binary_expr_block<1u, gt_functor,
			      Scalar_block<1, T>, T,
			      Block1, T> const, bool,
	   Scalar_block<1, T>, T,
	   Block1, T>
	SrcBlock;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    return base_type::rt_valid(dst,
			src.first().right(),
			src.first().left(),
			src.third(),
			src.second());
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    base_type::exec(dst,
		    src.first().right(),
		    src.first().left(),
		    src.third(),
		    src.second());
  }
};

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_SAL_EVAL_THRESHOLD_HPP
