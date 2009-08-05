/* Copyright (c) 2006, 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/eval_elementwise.hpp
    @author  Jules Bergmann
    @date    2006-05-26
    @brief   VSIPL++ Library: Dispatch for Mercury SAL.
*/

#ifndef VSIP_OPT_SAL_EVAL_ELEMENTWISE_HPP
#define VSIP_OPT_SAL_EVAL_ELEMENTWISE_HPP

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
#include <vsip/core/view_cast.hpp>
#include <vsip/opt/sal/is_op_supported.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace sal
{

/***********************************************************************
  Serial expression evaluator base classes
***********************************************************************/

/// Evaluator base class for view-view expressions with mixed value types
/// I.e. complex<float> = float * complex<float>

template <template <typename, typename> class Operator,
	  typename DstBlock,
	  typename LBlock,
	  typename RBlock,
	  typename LType,
	  typename RType>
struct Serial_expr_evaluator_base_mixed
{
  typedef Binary_expr_block<1, Operator, LBlock, LType, RBlock, RType>
    SrcBlock;

  typedef typename DstBlock::value_type dst_type;

  typedef typename sal::Effective_value_type<DstBlock>::type d_eff_t;
  typedef typename sal::Effective_value_type<LBlock>::type l_eff_t;
  typedef typename sal::Effective_value_type<RBlock>::type r_eff_t;

  static bool const ct_valid = 
    Type_equal<typename LBlock::value_type, LType>::value &&
    Type_equal<typename RBlock::value_type, RType>::value &&
    (!Is_expr_block<LBlock>::value || Is_scalar_block<LBlock>::value) &&
    (!Is_expr_block<RBlock>::value || Is_scalar_block<RBlock>::value) &&
     Is_op2_supported<Operator, l_eff_t, r_eff_t, d_eff_t>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     (Ext_data_cost<LBlock>::value == 0 || Is_scalar_block<LBlock>::value) &&
     (Ext_data_cost<RBlock>::value == 0 || Is_scalar_block<RBlock>::value);
  
  static bool rt_valid(DstBlock&, SrcBlock const&)
  {
    // SAL supports all strides
    return true;
  }
};

} // namespace vsip::impl::sal



/***********************************************************************
  Copy expression evaluators
***********************************************************************/

#define VSIP_IMPL_SAL_COPY_EXPR(OP, FUN)				\
template <typename DstBlock,						\
	  typename SrcBlock>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         SrcBlock,							\
         typename Type_if<Mercury_sal_tag,				\
                          Is_leaf_block<SrcBlock>::value>::type>	\
{									\
  static char const* name() { return "Expr_SAL_COPY"; }			\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<SrcBlock>::type eff_src_t;	\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<SrcBlock>::layout_type>::type		\
    src_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<SrcBlock>::value || Is_scalar_block<SrcBlock>::value) &&\
     sal::Is_op1_supported<OP, eff_src_t, eff_dst_t>::value&&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<SrcBlock>::value == 0 ||				\
      Is_scalar_block<SrcBlock>::value);				\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);		\
    sal::Ext_wrapper<SrcBlock, src_lp> ext_src(src, SYNC_IN);		\
    FUN(typename sal::Ext_wrapper<SrcBlock, dst_lp>::sal_type(ext_src),	\
	typename sal::Ext_wrapper<DstBlock, src_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_COPY_EXPR(sal::copy_token, vcopy)



/***********************************************************************
  Unary expression evaluators
***********************************************************************/

#define VSIP_IMPL_SAL_V_EXPR(OP, FUN)					\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         Unary_expr_block<1, OP, Block1, Type1> const,			\
         typename Type_if<Mercury_sal_tag,				\
                          Is_leaf_block<Block1>::value>::type>		\
{									\
  static char const* name() { return "Expr_SAL_V"; }			\
									\
  typedef Unary_expr_block<1, OP, Block1, Type1> const			\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block1>::layout_type>::type		\
    block1_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
     sal::Is_op1_supported<OP, eff_1_t, eff_dst_t>::value&&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,     SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp> ext_1(src.op(),  SYNC_IN);	\
    FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_V_EXPR(magsq_functor, sal::vmagsq)
VSIP_IMPL_SAL_V_EXPR(mag_functor,   sal::vmag)
VSIP_IMPL_SAL_V_EXPR(op::Minus,     sal::vneg)
VSIP_IMPL_SAL_V_EXPR(cos_functor,   sal::vcos)
VSIP_IMPL_SAL_V_EXPR(sin_functor,   sal::vsin)
VSIP_IMPL_SAL_V_EXPR(tan_functor,   sal::vtan)
VSIP_IMPL_SAL_V_EXPR(atan_functor,  sal::vatan)
VSIP_IMPL_SAL_V_EXPR(log_functor,   sal::vlog)
VSIP_IMPL_SAL_V_EXPR(log10_functor, sal::vlog10)
VSIP_IMPL_SAL_V_EXPR(exp_functor,   sal::vexp)
VSIP_IMPL_SAL_V_EXPR(exp10_functor, sal::vexp10)
VSIP_IMPL_SAL_V_EXPR(sqrt_functor,  sal::vsqrt)
VSIP_IMPL_SAL_V_EXPR(rsqrt_functor, sal::vrsqrt)
VSIP_IMPL_SAL_V_EXPR(sq_functor,    sal::vsq)
VSIP_IMPL_SAL_V_EXPR(recip_functor, sal::vrecip)

VSIP_IMPL_SAL_V_EXPR(Cast_closure<long          >::Cast, sal::vconv)
VSIP_IMPL_SAL_V_EXPR(Cast_closure<short         >::Cast, sal::vconv)
VSIP_IMPL_SAL_V_EXPR(Cast_closure<char          >::Cast, sal::vconv)
VSIP_IMPL_SAL_V_EXPR(Cast_closure<unsigned long >::Cast, sal::vconv)
VSIP_IMPL_SAL_V_EXPR(Cast_closure<unsigned short>::Cast, sal::vconv)
VSIP_IMPL_SAL_V_EXPR(Cast_closure<unsigned char >::Cast, sal::vconv)

VSIP_IMPL_SAL_V_EXPR(Cast_closure<float>::Cast, sal::vconv)



/***********************************************************************
  Binary expression evaluators
***********************************************************************/

#define VSIP_IMPL_SAL_VV_EXPR(OP, FUN)					\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename RBlock,						\
	  typename LType,						\
	  typename RType>						\
struct Serial_expr_evaluator<						\
    1, DstBlock, 							\
    const Binary_expr_block<1, OP, LBlock, LType, RBlock, RType>,	\
    typename Type_if<Mercury_sal_tag,					\
                     Is_leaf_block<LBlock>::value &&			\
                     Is_leaf_block<RBlock>::value>::type>		\
  : sal::Serial_expr_evaluator_base_mixed<OP, DstBlock, LBlock, RBlock, LType, RType>		\
{									\
  static char const* name() { return "Expr_SAL_VV"; }			\
									\
  typedef Binary_expr_block<1, OP, LBlock, LType, RBlock, RType>	\
    SrcBlock;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<LBlock>::layout_type>::type		\
    lblock_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<RBlock>::layout_type>::type		\
    rblock_lp;								\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);	\
    sal::Ext_wrapper<LBlock, lblock_lp> ext_l(src.left(),  SYNC_IN);	\
    sal::Ext_wrapper<RBlock, rblock_lp> ext_r(src.right(), SYNC_IN);	\
									\
    assert(dst.size() <= src.left().size()  || src.left().size() == 0);	\
    assert(dst.size() <= src.right().size() || src.right().size() == 0); \
									\
    VSIP_IMPL_COVER_BLK("SAL_VV", SrcBlock);				\
    FUN(								\
      typename sal::Ext_wrapper<LBlock, lblock_lp>::sal_type(ext_l),	\
      typename sal::Ext_wrapper<RBlock, lblock_lp>::sal_type(ext_r),	\
      typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
      dst.size()							\
    );									\
  }									\
};

VSIP_IMPL_SAL_VV_EXPR(op::Add,  sal::vadd)
VSIP_IMPL_SAL_VV_EXPR(op::Sub,  sal::vsub)
VSIP_IMPL_SAL_VV_EXPR(op::Mult, sal::vmul)
VSIP_IMPL_SAL_VV_EXPR(op::Div,  sal::vdiv)

VSIP_IMPL_SAL_VV_EXPR(max_functor, sal::vmax)
VSIP_IMPL_SAL_VV_EXPR(min_functor, sal::vmin)

VSIP_IMPL_SAL_VV_EXPR(band_functor, sal::vband)
VSIP_IMPL_SAL_VV_EXPR(bor_functor,  sal::vbor)


/***********************************************************************
  Ternary expression evaluators
***********************************************************************/

#define VSIP_IMPL_SAL_VVV_EXPR(OP, FUN)					\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1,						\
	  typename Block2,						\
	  typename Type2,						\
	  typename Block3,						\
	  typename Type3>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         const Ternary_expr_block<1, OP,				\
                                 Block1, Type1,				\
                                 Block2, Type2,				\
                                 Block3, Type3>,			\
         Mercury_sal_tag>						\
{									\
  static char const* name() { return "Expr_SAL_VVV"; }			\
									\
  typedef Ternary_expr_block<1, OP,					\
                                 Block1, Type1,				\
                                 Block2, Type2,				\
                                 Block3, Type3>				\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
  typedef typename sal::Effective_value_type<Block2, Type2>::type eff_2_t;\
  typedef typename sal::Effective_value_type<Block3, Type3>::type eff_3_t;\
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
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block3>::layout_type>::type		\
    block3_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
    (!Is_expr_block<Block2>::value || Is_scalar_block<Block2>::value) &&\
    (!Is_expr_block<Block3>::value || Is_scalar_block<Block3>::value) &&\
     sal::Is_op3_supported<OP, eff_1_t, eff_2_t, eff_3_t, eff_dst_t>::value&&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value) &&\
     (Ext_data_cost<Block2>::value == 0 || Is_scalar_block<Block2>::value) &&\
     (Ext_data_cost<Block3>::value == 0 || Is_scalar_block<Block3>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,        SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp> ext_1(src.first(),  SYNC_IN);	\
    sal::Ext_wrapper<Block2, block2_lp> ext_2(src.second(), SYNC_IN);	\
    sal::Ext_wrapper<Block3, block3_lp> ext_3(src.third(),  SYNC_IN);	\
    FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
        typename sal::Ext_wrapper<Block2, block1_lp>::sal_type(ext_2),	\
        typename sal::Ext_wrapper<Block3, block1_lp>::sal_type(ext_3),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_VVV_EXPR(ma_functor,  sal::vma)
VSIP_IMPL_SAL_VVV_EXPR(msb_functor, sal::vmsb)
VSIP_IMPL_SAL_VVV_EXPR(am_functor,  sal::vam)
VSIP_IMPL_SAL_VVV_EXPR(sbm_functor, sal::vsbm)



// Ternary expressions, F(V).V.V

#define VSIP_IMPL_SAL_fVVV_EXPR(OP, UOP, FUN, LOOKUP_OP)		\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1,						\
	  typename Block2,						\
	  typename Type2,						\
	  typename Block3,						\
	  typename Type3>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         Ternary_expr_block<1, OP,					\
           Unary_expr_block<1, UOP, Block1, Type1> const, Type1,	\
           Block2, Type2,						\
           Block3, Type3> const,					\
         Mercury_sal_tag>						\
{									\
  static char const* name() { return "Expr_SAL_fVVV"; }			\
									\
  typedef Ternary_expr_block<1, OP,					\
           Unary_expr_block<1, UOP, Block1, Type1> const, Type1,	\
           Block2, Type2,						\
           Block3, Type3>						\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
  typedef typename sal::Effective_value_type<Block2, Type2>::type eff_2_t;\
  typedef typename sal::Effective_value_type<Block3, Type3>::type eff_3_t;\
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
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block3>::layout_type>::type		\
    block3_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
    (!Is_expr_block<Block2>::value || Is_scalar_block<Block2>::value) &&\
    (!Is_expr_block<Block3>::value || Is_scalar_block<Block3>::value) &&\
     sal::Is_op3_supported<LOOKUP_OP, eff_1_t, eff_2_t, eff_3_t,	\
                           eff_dst_t>::value&&				\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value) &&\
     (Ext_data_cost<Block2>::value == 0 || Is_scalar_block<Block2>::value) &&\
     (Ext_data_cost<Block3>::value == 0 || Is_scalar_block<Block3>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,            SYNC_OUT);\
    sal::Ext_wrapper<Block1, block1_lp> ext_1(src.first().op(), SYNC_IN);\
    sal::Ext_wrapper<Block2, block2_lp> ext_2(src.second(),     SYNC_IN);\
    sal::Ext_wrapper<Block3, block3_lp> ext_3(src.third(),      SYNC_IN);\
    FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
        typename sal::Ext_wrapper<Block2, block2_lp>::sal_type(ext_2),	\
        typename sal::Ext_wrapper<Block3, block3_lp>::sal_type(ext_3),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_fVVV_EXPR(ma_functor, conj_functor, sal::vcma, sal::cma_token)



// Nested binary expressions, VV_V

#define VSIP_IMPL_SAL_VV_V_EXPR(OP, OP1, OP2, FUN)			\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1,						\
	  typename Block2,						\
	  typename Type2,						\
	  typename Block3,						\
	  typename Type3,						\
	  typename TypeB>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         Binary_expr_block<						\
                 1, OP2,						\
                 Binary_expr_block<					\
                         1, OP1,					\
                         Block1, Type1,					\
                         Block2, Type2> const, TypeB,			\
                 Block3, Type3> const,					\
         typename Type_if<Mercury_sal_tag,				\
                     Is_leaf_block<Block1>::value &&			\
                     Is_leaf_block<Block2>::value &&			\
                     Is_leaf_block<Block3>::value>::type>		\
{									\
  static char const* name() { return "Expr_SAL_VV_V"; }			\
									\
  typedef Binary_expr_block<						\
                 1, OP2,						\
                 Binary_expr_block<					\
                         1, OP1,					\
                         Block1, Type1,					\
                         Block2, Type2> const, TypeB,			\
                 Block3, Type3>						\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
  typedef typename sal::Effective_value_type<Block2, Type2>::type eff_2_t;\
  typedef typename sal::Effective_value_type<Block3, Type3>::type eff_3_t;\
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
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block3>::layout_type>::type		\
    block3_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
    (!Is_expr_block<Block2>::value || Is_scalar_block<Block2>::value) &&\
    (!Is_expr_block<Block3>::value || Is_scalar_block<Block3>::value) &&\
     sal::Is_op3_supported<OP, eff_1_t, eff_2_t, eff_3_t, eff_dst_t>::value&&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value) &&\
     (Ext_data_cost<Block2>::value == 0 || Is_scalar_block<Block2>::value) &&\
     (Ext_data_cost<Block3>::value == 0 || Is_scalar_block<Block3>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp> ext_dst(dst,        SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp>   ext_1(src.left().left(),  SYNC_IN);\
    sal::Ext_wrapper<Block2, block2_lp>   ext_2(src.left().right(), SYNC_IN);\
    sal::Ext_wrapper<Block3, block3_lp>   ext_3(src.right(),  SYNC_IN);	\
    FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
        typename sal::Ext_wrapper<Block2, block2_lp>::sal_type(ext_2),	\
        typename sal::Ext_wrapper<Block3, block3_lp>::sal_type(ext_3),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

// Nested binary expressions, V_VV

#define VSIP_IMPL_SAL_V_VV_EXPR(OP, OP1, OP2, FUN)			\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1,						\
	  typename Block2,						\
	  typename Type2,						\
	  typename Block3,						\
	  typename Type3,						\
	  typename TypeB>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         Binary_expr_block<						\
                 1, OP2,						\
                 Block1, Type1,						\
                 Binary_expr_block<					\
                         1, OP1,					\
                         Block2, Type2,					\
                         Block3, Type3> const, TypeB> const,		\
         typename Type_if<Mercury_sal_tag,				\
                     Is_leaf_block<Block1>::value &&			\
                     Is_leaf_block<Block2>::value &&			\
                     Is_leaf_block<Block3>::value>::type>		\
{									\
  static char const* name() { return "Expr_SAL_V_VV"; }			\
									\
  typedef Binary_expr_block<						\
                 1, OP2,						\
                 Block1, Type1,						\
                 Binary_expr_block<					\
                         1, OP1,					\
                         Block2, Type2,					\
                         Block3, Type3> const, TypeB>			\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
  typedef typename sal::Effective_value_type<Block2, Type2>::type eff_2_t;\
  typedef typename sal::Effective_value_type<Block3, Type3>::type eff_3_t;\
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
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block3>::layout_type>::type		\
    block3_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
    (!Is_expr_block<Block2>::value || Is_scalar_block<Block2>::value) &&\
    (!Is_expr_block<Block3>::value || Is_scalar_block<Block3>::value) &&\
     sal::Is_op3_supported<OP, eff_2_t, eff_3_t, eff_1_t, eff_dst_t>::value&&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value) &&\
     (Ext_data_cost<Block2>::value == 0 || Is_scalar_block<Block2>::value) &&\
     (Ext_data_cost<Block3>::value == 0 || Is_scalar_block<Block3>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,               SYNC_OUT);\
    sal::Ext_wrapper<Block1, block1_lp> ext_1(src.left(),          SYNC_IN);\
    sal::Ext_wrapper<Block2, block2_lp> ext_2(src.right().left(),  SYNC_IN);\
    sal::Ext_wrapper<Block3, block3_lp> ext_3(src.right().right(), SYNC_IN);\
    FUN(typename sal::Ext_wrapper<Block2, block2_lp>::sal_type(ext_2),	\
        typename sal::Ext_wrapper<Block3, block3_lp>::sal_type(ext_3),	\
        typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_VV_V_EXPR(ma_functor,  op::Mult, op::Add,  sal::vma)
VSIP_IMPL_SAL_VV_V_EXPR(msb_functor, op::Mult, op::Sub,  sal::vmsb)
VSIP_IMPL_SAL_VV_V_EXPR(am_functor,  op::Add,  op::Mult, sal::vam)
VSIP_IMPL_SAL_VV_V_EXPR(sbm_functor, op::Sub,  op::Mult, sal::vsbm)

// V OP2 (V OP1 V)
VSIP_IMPL_SAL_V_VV_EXPR(ma_functor,  op::Mult, op::Add,  sal::vma)
VSIP_IMPL_SAL_V_VV_EXPR(am_functor,  op::Add,  op::Mult, sal::vam)
VSIP_IMPL_SAL_V_VV_EXPR(sbm_functor, op::Sub,  op::Mult, sal::vsbm)



// Nested binary expressions, f(V)V_V

#define VSIP_IMPL_SAL_fVV_V_EXPR(OP, OP1, OP2, UOP, FUN)		\
template <typename DstBlock,						\
	  typename Block1,						\
	  typename Type1,						\
	  typename Block2,						\
	  typename Type2,						\
	  typename Block3,						\
	  typename Type3,						\
	  typename TypeB>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         Binary_expr_block<						\
           1, OP2,							\
           Binary_expr_block<						\
             1, OP1,							\
             Unary_expr_block<1, UOP, Block1, Type1> const, Type1,	\
             Block2, Type2> const, TypeB,				\
           Block3, Type3> const,					\
         Mercury_sal_tag>						\
{									\
  static char const* name() { return "Expr_SAL_fVV_V"; }		\
									\
  typedef Binary_expr_block<						\
            1, OP2,							\
            Binary_expr_block<						\
              1, OP1,							\
              Unary_expr_block<1, UOP, Block1, Type1> const, Type1,	\
              Block2, Type2> const, TypeB,				\
            Block3, Type3>						\
	SrcBlock;							\
									\
  typedef typename DstBlock::value_type dst_type;			\
									\
  typedef typename sal::Effective_value_type<DstBlock>::type eff_dst_t;	\
  typedef typename sal::Effective_value_type<Block1, Type1>::type eff_1_t;\
  typedef typename sal::Effective_value_type<Block2, Type2>::type eff_2_t;\
  typedef typename sal::Effective_value_type<Block3, Type3>::type eff_3_t;\
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
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block3>::layout_type>::type		\
    block3_lp;								\
  									\
  static bool const ct_valid = 						\
    (!Is_expr_block<Block1>::value || Is_scalar_block<Block1>::value) &&\
    (!Is_expr_block<Block2>::value || Is_scalar_block<Block2>::value) &&\
    (!Is_expr_block<Block3>::value || Is_scalar_block<Block3>::value) &&\
     sal::Is_op3_supported<OP, eff_1_t, eff_2_t, eff_3_t, eff_dst_t>::value&&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     (Ext_data_cost<Block1>::value == 0 || Is_scalar_block<Block1>::value) &&\
     (Ext_data_cost<Block2>::value == 0 || Is_scalar_block<Block2>::value) &&\
     (Ext_data_cost<Block3>::value == 0 || Is_scalar_block<Block3>::value);\
									\
  static bool rt_valid(DstBlock&, SrcBlock const&)			\
  {									\
    /* SAL supports all strides */					\
    return true;							\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    sal::Ext_wrapper<DstBlock, dst_lp>  ext_dst(dst,        SYNC_OUT);	\
    sal::Ext_wrapper<Block1, block1_lp> ext_1(src.left().left().op(),SYNC_IN);\
    sal::Ext_wrapper<Block2, block2_lp> ext_2(src.left().right(),    SYNC_IN);\
    sal::Ext_wrapper<Block3, block3_lp> ext_3(src.right(),           SYNC_IN);\
    FUN(typename sal::Ext_wrapper<Block1, block1_lp>::sal_type(ext_1),	\
        typename sal::Ext_wrapper<Block2, block2_lp>::sal_type(ext_2),	\
        typename sal::Ext_wrapper<Block3, block3_lp>::sal_type(ext_3),	\
	typename sal::Ext_wrapper<DstBlock, dst_lp>::sal_type(ext_dst),	\
	dst.size());							\
  }									\
};

VSIP_IMPL_SAL_fVV_V_EXPR(sal::cma_token, op::Mult, op::Add, conj_functor,
			 sal::vcma)

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_SAL_EVAL_ELEMENTWISE_HPP
