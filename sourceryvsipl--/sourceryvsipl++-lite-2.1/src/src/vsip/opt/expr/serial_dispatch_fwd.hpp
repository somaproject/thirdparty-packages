/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/expr/serial_dispatch_fwd.hpp
    @author  Stefan Seefeld
    @date    2005-08-05
    @brief   VSIPL++ Library: Forward Decl of Serial_dispatch_helper.
*/

#ifndef VSIP_OPT_EXPR_SERIAL_DISPATCH_FWD_HPP
#define VSIP_OPT_EXPR_SERIAL_DISPATCH_FWD_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/type_list.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

/// The list of evaluators to be tried, in that specific order.

/// Note that the VSIP_IMPL_TAG_LIST macro will include its own comma.

typedef Make_type_list<Intel_ipp_tag,
                       Cml_tag,
		       Transpose_tag,
                       Mercury_sal_tag,
                       Cbe_sdk_tag,
                       Cuda_tag,
                       VSIP_IMPL_SIMD_TAG_LIST
#if VSIP_IMPL_ENABLE_EVAL_DENSE_EXPR
		       Dense_expr_tag,
#endif
		       Copy_tag,
		       Op_expr_tag,
		       Simd_loop_fusion_tag,
		       Simd_unaligned_loop_fusion_tag,
		       Fc_expr_tag,
		       Rbo_expr_tag,
		       Mdim_expr_tag,
		       Loop_fusion_tag>::type LibraryTagList;



/// Serial_dispatch_helper dispatches the evaluation of an expression along
/// a type list of potential backends.
/// Whether a given backend is actually used depends on its compile-time
/// and run-time validity checks.

/// Requires:
///   PROFILEP to be a profiling policy.  Used to optionaly add
///      profiling or coverage to dispatch.

template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  template <typename> class ProfileP,
	  typename TagList,
	  typename Tag = typename TagList::first,
	  typename Rest = typename TagList::rest,
	  typename EvalExpr = Serial_expr_evaluator<
  Dim, DstBlock, SrcBlock, Tag>,
	  bool CtValid = EvalExpr::ct_valid>
struct Serial_dispatch_helper;



/// Front-end to Serial_dispatch_helper.  This should be used
/// instead of S_d_h directly.

template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  typename TagList,
	  bool     Internal = false>
struct Serial_dispatch;

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_EXPR_SERIAL_DISPATCH_FWD_HPP
