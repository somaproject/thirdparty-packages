/* Copyright (c) 2006, 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/eval-generic.hpp
    @author  Jules Bergmann
    @date    2006-01-25
    @brief   VSIPL++ Library: Wrappers and traits to bridge with generic SIMD.
*/

#ifndef VSIP_IMPL_SIMD_EVAL_GENERIC_HPP
#define VSIP_IMPL_SIMD_EVAL_GENERIC_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/coverage.hpp>

#include <vsip/opt/simd/simd.hpp>
#include <vsip/opt/simd/vadd.hpp>
#include <vsip/opt/simd/vmul.hpp>
#include <vsip/opt/simd/rscvmul.hpp>
#include <vsip/opt/simd/vgt.hpp>
#include <vsip/opt/simd/vlogic.hpp>
#include <vsip/opt/simd/threshold.hpp>
#include <vsip/opt/simd/expr_iterator.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace simd
{

template <template <typename, typename> class Operator>
struct Map_operator_to_algorithm
{
  typedef Alg_none type;
};

template <>
struct Map_operator_to_algorithm<op::Add>  { typedef Alg_vadd type; };
template <>
struct Map_operator_to_algorithm<op::Mult> { typedef Alg_vmul type; };
template <>
struct Map_operator_to_algorithm<band_functor> { typedef Alg_vband type; };
template <>
struct Map_operator_to_algorithm<bor_functor> { typedef Alg_vbor type; };
template <>
struct Map_operator_to_algorithm<bxor_functor> { typedef Alg_vbxor type; };



template <template <typename, typename> class Operator,
	  typename DstBlock,
	  typename LBlock,
	  typename RBlock,
	  typename LType,
	  typename RType>
struct Serial_expr_evaluator_base
{
  typedef Binary_expr_block<1, Operator, LBlock, LType, RBlock, RType>
    SrcBlock;

  static bool const ct_valid = 
    !Is_expr_block<LBlock>::value &&
    !Is_expr_block<RBlock>::value &&
    simd::Is_algorithm_supported<
        typename DstBlock::value_type,
        Is_split_block<DstBlock>::value,
	typename Map_operator_to_algorithm<Operator>::type>::value &&
     Type_equal<typename DstBlock::value_type, LType>::value &&
     Type_equal<typename DstBlock::value_type, RType>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<LBlock>::value == 0 &&
     Ext_data_cost<RBlock>::value == 0 &&
     // Must have same complex interleaved/split format
     Type_equal<typename Block_layout<DstBlock>::complex_type,
		typename Block_layout<LBlock>::complex_type>::value &&
     Type_equal<typename Block_layout<DstBlock>::complex_type,
		typename Block_layout<RBlock>::complex_type>::value;

  
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock> ext_dst(dst, SYNC_OUT);
    Ext_data<LBlock> ext_l(src.left(), SYNC_IN);
    Ext_data<RBlock> ext_r(src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_l.stride(0) == 1 &&
	    ext_r.stride(0) == 1);
  }
};
} // namespace vsip::impl::simd



#define VSIP_IMPL_SIMD_V_EXPR(OP, ALG, FCN)				\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename LType>						\
struct Serial_expr_evaluator<						\
  1, DstBlock, 								\
  const Unary_expr_block<1, OP, LBlock, LType>,				\
  Simd_builtin_tag>							\
{									\
  typedef Unary_expr_block<1, OP, LBlock, LType>			\
    SrcBlock;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<LBlock>::layout_type>::type		\
    lblock_lp;								\
  									\
  static char const* name() { return "Expr_SIMD_V-" #FCN; }		\
  									\
  static bool const ct_valid = 						\
    !Is_expr_block<LBlock>::value &&					\
    simd::Is_algorithm_supported<					\
        typename DstBlock::value_type,					\
        Is_split_block<DstBlock>::value,				\
	ALG>::value &&							\
     Type_equal<typename DstBlock::value_type, LType>::value &&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<LBlock>::value == 0 &&				\
     /* Must have same complex interleaved/split format */		\
     Type_equal<typename Block_layout<DstBlock>::complex_type,		\
		typename Block_layout<LBlock>::complex_type>::value;	\
  									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);			\
    Ext_data<LBlock, lblock_lp> ext_l(src.op(), SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_l.stride(0) == 1);					\
  }									\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);			\
    Ext_data<LBlock, lblock_lp> ext_l(src.op(), SYNC_IN);		\
    FCN(ext_l.data(), ext_dst.data(), dst.size());			\
  }									\
};

#define VSIP_IMPL_SIMD_VV_EXPR(OP, FCN)					\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename RBlock,						\
	  typename LType,						\
	  typename RType>						\
struct Serial_expr_evaluator<						\
  1, DstBlock, 								\
  const Binary_expr_block<1, OP, LBlock, LType, RBlock, RType>,		\
  Simd_builtin_tag>							\
{									\
  typedef Binary_expr_block<1, OP, LBlock, LType, RBlock, RType>	\
    SrcBlock;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<LBlock>::layout_type>::type		\
    lblock_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<RBlock>::layout_type>::type		\
    rblock_lp;								\
  									\
  static char const* name() { return "Expr_SIMD_VV-" #FCN; }		\
  									\
  static bool const ct_valid = 						\
    !Is_expr_block<LBlock>::value &&					\
    !Is_expr_block<RBlock>::value &&					\
    simd::Is_algorithm_supported<					\
        typename DstBlock::value_type,					\
        Is_split_block<DstBlock>::value,				\
	typename simd::Map_operator_to_algorithm<OP>::type>::value &&	\
     Type_equal<typename DstBlock::value_type, LType>::value &&		\
     Type_equal<typename DstBlock::value_type, RType>::value &&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<LBlock>::value == 0 &&				\
     Ext_data_cost<RBlock>::value == 0 &&				\
     /* Must have same complex interleaved/split format */		\
     Type_equal<typename Block_layout<DstBlock>::complex_type,		\
		typename Block_layout<LBlock>::complex_type>::value &&	\
     Type_equal<typename Block_layout<DstBlock>::complex_type,		\
		typename Block_layout<RBlock>::complex_type>::value;	\
  									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst, SYNC_OUT);			\
    Ext_data<LBlock, lblock_lp> ext_l(src.left(), SYNC_IN);		\
    Ext_data<RBlock, rblock_lp> ext_r(src.right(), SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_l.stride(0) == 1 &&					\
	    ext_r.stride(0) == 1);					\
  }									\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst, SYNC_OUT);			\
    Ext_data<LBlock, lblock_lp> ext_l(src.left(), SYNC_IN);		\
    Ext_data<RBlock, rblock_lp> ext_r(src.right(), SYNC_IN);		\
    FCN(ext_l.data(), ext_r.data(), ext_dst.data(), dst.size());	\
  }									\
};


VSIP_IMPL_SIMD_V_EXPR (bnot_functor, simd::Alg_vbnot, simd::vbnot)

VSIP_IMPL_SIMD_VV_EXPR(op::Mult,     simd::vmul)
VSIP_IMPL_SIMD_VV_EXPR(op::Add,      simd::vadd)
VSIP_IMPL_SIMD_VV_EXPR(band_functor, simd::vband)
VSIP_IMPL_SIMD_VV_EXPR(bor_functor,  simd::vbor)
VSIP_IMPL_SIMD_VV_EXPR(bxor_functor, simd::vbxor)

#undef VSIP_IMPL_SIMD_V_EXPR
#undef VSIP_IMPL_SIMD_VV_EXPR



/***********************************************************************
  vgt: vector greater-than operator
***********************************************************************/

template <typename DstBlock,
	  typename LBlock,
	  typename RBlock,
	  typename LType,
	  typename RType>
struct Serial_expr_evaluator<
  1, DstBlock,
  const Binary_expr_block<1, gt_functor, LBlock, LType, RBlock, RType>,
  Simd_builtin_tag>
{
  typedef Binary_expr_block<1, gt_functor, LBlock, LType, RBlock, RType>
    SrcBlock;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;
  typedef typename Adjust_layout_dim<
      1, typename Block_layout<LBlock>::layout_type>::type
    lblock_lp;
  typedef typename Adjust_layout_dim<
      1, typename Block_layout<RBlock>::layout_type>::type
    rblock_lp;

  static char const* name() { return "Expr_SIMD_VV-simd::vgt"; }

  static bool const ct_valid = 
    !Is_expr_block<LBlock>::value &&
    !Is_expr_block<RBlock>::value &&
     Type_equal<typename DstBlock::value_type, bool>::value &&
     Type_equal<LType, RType>::value &&
     simd::Is_algorithm_supported<LType, false, simd::Alg_vgt>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<LBlock>::value == 0 &&
     Ext_data_cost<RBlock>::value == 0;

  
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_l.stride(0) == 1 &&
	    ext_r.stride(0) == 1);
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);
    simd::vgt(ext_l.data(), ext_r.data(), ext_dst.data(), dst.size());
  }
};



template <typename DstBlock,
	  typename LBlock,
	  typename RBlock,
	  typename LType,
	  typename RType>
struct Serial_expr_evaluator<
  1, DstBlock,
  const Binary_expr_block<1, lt_functor, LBlock, LType, RBlock, RType>,
  Simd_builtin_tag>
{
  typedef Binary_expr_block<1, lt_functor, LBlock, LType, RBlock, RType>
    SrcBlock;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;
  typedef typename Adjust_layout_dim<
      1, typename Block_layout<LBlock>::layout_type>::type
    lblock_lp;
  typedef typename Adjust_layout_dim<
      1, typename Block_layout<RBlock>::layout_type>::type
    rblock_lp;

  static char const* name() { return "Expr_SIMD_VV-simd::vlt_as_gt"; }

  static bool const ct_valid = 
    !Is_expr_block<LBlock>::value &&
    !Is_expr_block<RBlock>::value &&
     Type_equal<typename DstBlock::value_type, bool>::value &&
     Type_equal<LType, RType>::value &&
     simd::Is_algorithm_supported<LType, false, simd::Alg_vgt>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<LBlock>::value == 0 &&
     Ext_data_cost<RBlock>::value == 0;

  
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_l.stride(0) == 1 &&
	    ext_r.stride(0) == 1);
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);
    // Swap left and right arguments to vgt
    simd::vgt(ext_r.data(), ext_l.data(), ext_dst.data(), dst.size());
  }
};



/***********************************************************************
  vector logical operators
***********************************************************************/

#define VSIP_IMPL_SIMD_LOGIC_V_EXPR(OP, ALG, FCN)			\
template <typename DstBlock,						\
	  typename BlockT>						\
struct Serial_expr_evaluator<						\
  1, DstBlock,								\
  const Unary_expr_block<1, OP, BlockT, bool>,				\
  Simd_builtin_tag>							\
{									\
  typedef Unary_expr_block<1, OP, BlockT, bool>				\
    SrcBlock;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<BlockT>::layout_type>::type		\
    block_lp;								\
									\
  static char const* name() { return "Expr_SIMD_V-" #FCN; }		\
  									\
  static bool const ct_valid = 						\
    !Is_expr_block<BlockT>::value &&					\
     Type_equal<typename DstBlock::value_type, bool>::value &&		\
     simd::Is_algorithm_supported<bool, false, ALG>::value &&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<BlockT>::value == 0;					\
									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);		\
    Ext_data<BlockT, block_lp>  ext_l  (src.op(),  SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_l.stride(0) == 1);					\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);		\
    Ext_data<BlockT, block_lp>  ext_l  (src.op(),  SYNC_IN);		\
    FCN(ext_l.data(), ext_dst.data(), dst.size());			\
  }									\
};

#define VSIP_IMPL_SIMD_LOGIC_VV_EXPR(OP, ALG, FCN)			\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename RBlock>						\
struct Serial_expr_evaluator<						\
  1, DstBlock,								\
  const Binary_expr_block<1, OP, LBlock, bool, RBlock, bool>,		\
  Simd_builtin_tag>							\
{									\
  typedef Binary_expr_block<1, OP, LBlock, bool, RBlock, bool>		\
    SrcBlock;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<LBlock>::layout_type>::type		\
    lblock_lp;								\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<RBlock>::layout_type>::type		\
    rblock_lp;								\
									\
  static char const* name() { return "Expr_SIMD_VV-" #FCN; }		\
  									\
  static bool const ct_valid = 						\
    !Is_expr_block<LBlock>::value &&					\
    !Is_expr_block<RBlock>::value &&					\
     Type_equal<typename DstBlock::value_type, bool>::value &&		\
     simd::Is_algorithm_supported<bool, false, ALG>::value &&\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<LBlock>::value == 0 &&				\
     Ext_data_cost<RBlock>::value == 0;					\
									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);		\
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);		\
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_l.stride(0) == 1 &&					\
	    ext_r.stride(0) == 1);					\
  }									\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,         SYNC_OUT);		\
    Ext_data<LBlock, lblock_lp> ext_l  (src.left(),  SYNC_IN);		\
    Ext_data<RBlock, rblock_lp> ext_r  (src.right(), SYNC_IN);		\
    FCN(ext_l.data(), ext_r.data(), ext_dst.data(), dst.size());	\
  }									\
};

VSIP_IMPL_SIMD_LOGIC_V_EXPR (lnot_functor, simd::Alg_vlnot, simd::vlnot)
VSIP_IMPL_SIMD_LOGIC_VV_EXPR(land_functor, simd::Alg_vland, simd::vland)
VSIP_IMPL_SIMD_LOGIC_VV_EXPR(lor_functor,  simd::Alg_vlor,  simd::vlor)
VSIP_IMPL_SIMD_LOGIC_VV_EXPR(lxor_functor, simd::Alg_vlxor, simd::vlxor)

#undef VSIP_IMPL_SIMD_LOGIC_V_EXPR
#undef VSIP_IMPL_SIMD_LOGIC_VV_EXPR


/***********************************************************************
  Scalar-view element-wise operations
***********************************************************************/

// Evaluate real-scalar * complex-view

template <typename DstBlock,
	  typename T,
	  typename VBlock>
struct Serial_expr_evaluator<
         1, DstBlock, 
         const Binary_expr_block<1, op::Mult,
                                 Scalar_block<1, T>, T,
                                 VBlock, std::complex<T> >,
         Simd_builtin_tag>
{
  typedef Binary_expr_block<1, op::Mult,
			    Scalar_block<1, T>, T,
			    VBlock, complex<T> >
	SrcBlock;

  typedef typename Adjust_layout_dim<
    1, typename Block_layout<DstBlock>::layout_type>::type
  dst_lp;
  typedef typename Adjust_layout_dim<
    1, typename Block_layout<VBlock>::layout_type>::type
  vblock_lp;

  static char const* name() { return "Expr_SIMD_V-simd::rscvmul"; }

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
    simd::Is_algorithm_supported<
        T,
        Is_split_block<DstBlock>::value,
	typename simd::Map_operator_to_algorithm<op::Mult>::type>::value &&

    Type_equal<typename DstBlock::value_type, std::complex<T> >::value &&
    // check that direct access is supported
    Ext_data_cost<DstBlock>::value == 0 &&
    Ext_data_cost<VBlock>::value == 0 &&
    // Must have same complex interleaved/split format
    Type_equal<typename Block_layout<DstBlock>::complex_type,
	       typename Block_layout<VBlock>::complex_type>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_r(src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 && ext_r.stride(0) == 1);
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_r(src.right(), SYNC_IN);
    simd::rscvmul(src.left().value(), ext_r.data(), ext_dst.data(),
		  dst.size());
  }
};



// Evaluate complex-view * real-scalar

template <typename DstBlock,
	  typename T,
	  typename VBlock>
struct Serial_expr_evaluator<
         1, DstBlock, 
         const Binary_expr_block<1, op::Mult,
                                 VBlock, std::complex<T>,
                                 Scalar_block<1, T>, T>,
         Simd_builtin_tag>
{
  typedef Binary_expr_block<1, op::Mult,
			    VBlock, std::complex<T>,
			    Scalar_block<1, T>, T>
	SrcBlock;

  typedef typename Adjust_layout_dim<
    1, typename Block_layout<DstBlock>::layout_type>::type
  dst_lp;
  typedef typename Adjust_layout_dim<
    1, typename Block_layout<VBlock>::layout_type>::type
  vblock_lp;

  static char const* name() { return "Expr_SIMD_V-simd::rscvmul"; }

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
    simd::Is_algorithm_supported<
        T,
        Is_split_block<DstBlock>::value,
	typename simd::Map_operator_to_algorithm<op::Mult>::type>::value &&

    Type_equal<typename DstBlock::value_type, std::complex<T> >::value &&
    // check that direct access is supported
    Ext_data_cost<DstBlock>::value == 0 &&
    Ext_data_cost<VBlock>::value == 0 &&
    // Must have same complex interleaved/split format
    Type_equal<typename Block_layout<DstBlock>::complex_type,
	       typename Block_layout<VBlock>::complex_type>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_l(src.left(), SYNC_IN);
    return (ext_dst.stride(0) == 1 && ext_l.stride(0) == 1);
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp> ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_l(src.left(), SYNC_IN);
    simd::rscvmul(src.right().value(), ext_l.data(), ext_dst.data(),
		  dst.size());
  }
};

/***********************************************************************
  threshold: vector threshold operator
  ite(A > B, A, k)
  ite(A < B, A, k)
  ite(A <= B, A, k)
  ite(A >= B, A, k)
***********************************************************************/

template <typename DstBlock,
          typename T,
	  typename Block1,
	  typename Block2,
	  template <typename,typename> class O>
struct Serial_expr_evaluator<
  1, DstBlock,

  Ternary_expr_block<1, ite_functor,
    Binary_expr_block<1u, O,
                      Block1, T,
		      Block2, T> const, bool,
    Block1, T,
    Scalar_block<1,T>, T> const,

    Simd_builtin_tag>
                      
{

  typedef Ternary_expr_block<1, ite_functor,
    Binary_expr_block<1u, O,
                      Block1, T,
		      Block2, T> const, bool,
    Block1, T,
    Scalar_block<1,T>, T> SrcBlock;


  static char const* name() { return "Expr_SIMD_threshold"; }

  typedef typename Adjust_layout_dim<
    1, typename Block_layout<DstBlock>::layout_type>::type
  dst_lp;
  typedef typename Adjust_layout_dim<
    1, typename Block_layout<Block1>::layout_type>::type
  a_lp;
  typedef typename Adjust_layout_dim<
    1, typename Block_layout<Block2>::layout_type>::type
  b_lp;

  static bool const ct_valid = 
    // Check that LHS & RHS have same type.
    Type_equal<typename DstBlock::value_type, T>::value &&
    // Make sure algorithm/op is supported.
    simd::Is_algorithm_supported<T, false, simd::Alg_threshold>::value &&
    simd::Binary_operator_map<T,O>::is_supported &&
    // Check that direct access is supported.
    Ext_data_cost<DstBlock>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;
  
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    typedef simd::Simd_traits<typename SrcBlock::value_type> simd;

    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>     ext_dst(dst,              SYNC_OUT);
    Ext_data<Block1,   a_lp>       ext_a(src.first().left(), SYNC_IN);
    Ext_data<Block2,   b_lp>       ext_b(src.first().right(),SYNC_IN);
    return(ext_dst.stride(0) == 1 &&
           ext_a.stride(0) == 1 &&
	   ext_b.stride(0) == 1 &&
	   // make sure (A op B, A, k)
	   (&(src.first().left()) == &(src.second())) &&
	   // make sure everyting has same alignment
           (simd::alignment_of(ext_dst.data()) ==
	    simd::alignment_of(ext_a.data())) &&
           (simd::alignment_of(ext_dst.data()) ==
	    simd::alignment_of(ext_b.data())));
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp>     ext_dst(dst,              SYNC_OUT);
    Ext_data<Block1,   a_lp>       ext_a(src.first().left(), SYNC_IN);
    Ext_data<Block2,   b_lp>       ext_b(src.first().right(),SYNC_IN);

    simd::threshold<O>(ext_dst.data(), ext_a.data(), ext_b.data(),
                       src.third().value(), ext_dst.size());
  }
};




} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_SIMD_EVAL_GENERIC_HPP
