/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/eval_vcmp.hpp
    @author  Jules Bergmann
    @date    2006-10-26
    @brief   VSIPL++ Library: Dispatch for Mercury SAL -- vector comparisons.
*/

#ifndef VSIP_OPT_SAL_EVAL_VCMP_HPP
#define VSIP_OPT_SAL_EVAL_VCMP_HPP

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
#include <vsip/opt/sal/is_op_supported.hpp>



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

// Optimize threshold expression: ite(A > B, VAL1, VAL0)
#define VSIP_IMPL_SAL_VCMP_EXPR(FUNCTOR, OPTOKEN, FUN, VAL1, VAL0)	\
template <typename DstBlock,						\
	  typename T,							\
	  typename Block2,						\
	  typename Block1>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
									\
         Ternary_expr_block<1, ite_functor,				\
           Binary_expr_block<1u, FUNCTOR,				\
			     Block1, T,					\
			     Block2, T> const, bool,			\
	   Scalar_block<1, T>, T,					\
	   Scalar_block<1, T>, T> const,				\
									\
         Mercury_sal_tag>						\
{									\
  static char const* name() { return "Expr_SAL_vsmp-" # FUNCTOR; }	\
									\
  typedef Ternary_expr_block<1, ite_functor,				\
            Binary_expr_block<1u, FUNCTOR,				\
			      Block1, T,				\
			      Block2, T> const, bool,			\
	   Scalar_block<1, T>, T,					\
	   Scalar_block<1, T>, T>					\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type  eff_d_t;	\
  typedef typename sal::Effective_value_type<Block1, T>::type eff_1_t;	\
  typedef typename sal::Effective_value_type<Block2, T>::type eff_2_t;	\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block1>::layout_type>::type		\
    block1_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block2>::layout_type>::type		\
    block2_lp;								\
  									\
  static bool const ct_valid =						\
     sal::Is_op2_supported<OPTOKEN, eff_1_t, eff_2_t, eff_d_t>::value &&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<Block1>::value == 0 &&				\
     Ext_data_cost<Block2>::value == 0;					\
									\
  static bool rt_valid(DstBlock&, SrcBlock const& src)			\
  {									\
    return src.second().value() == T(VAL1) &&				\
            src.third().value() == T(VAL0);				\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    typedef Scalar_block<1, T> sb_type;					\
									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,        SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp> ext_A(src.first().left(), SYNC_IN);\
    sal::Ext_wrapper<Block2, block2_lp> ext_B(src.first().right(), SYNC_IN);\
									\
    FUN(								\
      typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_A),	\
      typename sal::Ext_wrapper<Block2, block2_lp>::sal_type(ext_B),	\
      typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
      dst.size());							\
  }									\
};

VSIP_IMPL_SAL_VCMP_EXPR(eq_functor, sal::veq_token, sal::lveq, 1, 0)
VSIP_IMPL_SAL_VCMP_EXPR(ne_functor, sal::vne_token, sal::lvne, 1, 0)
VSIP_IMPL_SAL_VCMP_EXPR(gt_functor, sal::vgt_token, sal::lvgt, 1, 0)
VSIP_IMPL_SAL_VCMP_EXPR(ge_functor, sal::vge_token, sal::lvge, 1, 0)
VSIP_IMPL_SAL_VCMP_EXPR(lt_functor, sal::vlt_token, sal::lvlt, 1, 0)
VSIP_IMPL_SAL_VCMP_EXPR(le_functor, sal::vle_token, sal::lvle, 1, 0)

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_SAL_EVAL_VCMP_HPP
