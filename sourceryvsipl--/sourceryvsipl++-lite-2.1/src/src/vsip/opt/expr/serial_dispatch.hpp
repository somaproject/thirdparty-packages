/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/expr/serial_dispatch.hpp
    @author  Stefan Seefeld
    @date    2005-08-05
    @brief   VSIPL++ Library: Serial dispatch harness to allow to bind
    various backends to be bound to particular expressions.
*/

#ifndef VSIP_OPT_EXPR_SERIAL_DISPATCH_HPP
#define VSIP_OPT_EXPR_SERIAL_DISPATCH_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/type_list.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/opt/expr/serial_dispatch_fwd.hpp>
#include <vsip/opt/expr/eval_mcopy.hpp>
#include <vsip/opt/expr/eval_mdim.hpp>
#include <vsip/opt/expr/eval_dense.hpp>
#include <vsip/opt/expr/eval_return_block.hpp>
#include <vsip/opt/expr/eval_fastconv.hpp>
#include <vsip/opt/expr/ops_info.hpp>
#include <vsip/core/profile.hpp>

#ifdef VSIP_IMPL_HAVE_IPP
#  include <vsip/opt/ipp/bindings.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/bindings.hpp>
#endif
#ifdef VSIP_IMPL_CBE_SDK
#  include <vsip/opt/cbe/cml/transpose.hpp>
#  include <vsip/opt/cbe/ppu/bindings.hpp>
#  include <vsip/opt/cbe/ppu/eval_fastconv.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_CUDA
#  include <vsip/opt/cuda/eval_vmmul.hpp>
#  include <vsip/opt/cuda/eval_fastconv.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SIMD_LOOP_FUSION
#  include <vsip/opt/simd/expr_evaluator.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SIMD_UNALIGNED_LOOP_FUSION
#  include <vsip/opt/simd/eval_unaligned.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SIMD_GENERIC
#  include <vsip/opt/simd/eval_generic.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SIMD_3DNOWEXT
#  include <vsip/opt/simd/eval_simd_3dnowext.hpp>
#endif

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

// Policy for profiling an expression evaluator.
template <typename EvalExpr>
struct Eval_profile_policy
{
  typedef profile::Scope<profile::fns> scope_type;

  template <typename DstBlock,
	    typename SrcBlock>
  Eval_profile_policy(DstBlock const&, SrcBlock const& src)
    : scope_(Expr_op_name<EvalExpr, SrcBlock>::tag(src), 
             Expr_ops_per_point<SrcBlock>::value == 0
             // If ops_per_point is 0, then assume that operations
             // is a copy and record the number of bytes written.
             ? sizeof(typename DstBlock::value_type) *
             Expr_ops_per_point<SrcBlock>::size(src) 
             // Otherwise, record the number of flops.
             : Expr_ops_per_point<SrcBlock>::value * 
             Expr_ops_per_point<SrcBlock>::size(src))
  {}

private:
  scope_type scope_;
};



// Lightweight policy for profiling an expression evaluator.

template <typename EvalExpr>
struct Eval_lw_profile_policy
{
  typedef profile::Scope<profile::fns> scope_type;

  template <typename DstBlock,
	    typename SrcBlock>
  Eval_lw_profile_policy(DstBlock const&, SrcBlock const& src)
    : scope_(EvalExpr::name(),
             Expr_ops_per_point<SrcBlock>::value == 0
             // If ops_per_point is 0, then assume that operations
             // is a copy and record the number of bytes written.
             ? sizeof(typename DstBlock::value_type) *
	       Expr_ops_per_point<SrcBlock>::size(src) 
             // Otherwise, record the number of flops.
             : Expr_ops_per_point<SrcBlock>::value * 
               Expr_ops_per_point<SrcBlock>::size(src))
  {}

private:
  scope_type scope_;
};



// Policy for recording coverage for an expression evaluator.
template <typename EvalExpr>
struct Eval_coverage_policy
{
  template <typename DstBlock,
	    typename SrcBlock>
  Eval_coverage_policy(DstBlock const&, SrcBlock const&)
  {
    char const* evaluator_name = EvalExpr::name();
    VSIP_IMPL_COVER_BLK(evaluator_name, SrcBlock);
  }
};



// Policy for doing nothing special for an expression evaluator.
template <typename EvalExpr>
struct Eval_nop_policy
{
  template <typename DstBlock,
	    typename SrcBlock>
  Eval_nop_policy(DstBlock const&, SrcBlock const&)
  {}
};



// Choose appropriate profiling policy for user generated
// Serial_dispatch, based on settings of VSIP_IMPL_DO_COVERAGE and
// profile::mask.

template <typename E>
struct Eval_sd_policy
{
  typedef typename ITE_Type<VSIP_IMPL_DO_COVERAGE,
                            As_type<Eval_coverage_policy<E> >,
                            ITE_Type<profile::mask & profile::fns,
                                     As_type<Eval_profile_policy<E> >,
                                     As_type<Eval_nop_policy<E> > > >::type 
    base_type;

  template <typename DstBlock, typename SrcBlock>
  Eval_sd_policy(DstBlock const &dst, SrcBlock const &src) : base(dst, src) {}

  base_type base;
};



// Choose appropriate profiling policy for user generated
// Serial_dispatch, based on settings of VSIP_IMPL_DO_COVERAGE and
// profile::mask.

template <typename E>
struct Eval_int_sd_policy
{
  typedef typename ITE_Type<VSIP_IMPL_DO_COVERAGE,
                            As_type<Eval_coverage_policy<E> >,
                            ITE_Type<profile::mask & profile::fns_int,
                                     As_type<Eval_lw_profile_policy<E> >,
                                     As_type<Eval_nop_policy<E> > > >::type 
    base_type;

  template <typename DstBlock, typename SrcBlock>
  Eval_int_sd_policy(DstBlock const &dst, SrcBlock const &src)
    : base(dst, src)
  {}

  base_type base;
};



/// In case the compile-time check passes, we decide at run-time whether
/// or not to use this backend.
template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  template <typename> class ProfileP,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename EvalExpr>
struct Serial_dispatch_helper<Dim, DstBlock, SrcBlock, ProfileP, TagList,
		       Tag, Rest, EvalExpr, true>
{
  static void exec(DstBlock& dst, SrcBlock const& src)
    VSIP_NOTHROW
  {
    if (EvalExpr::rt_valid(dst, src))
    {
      ProfileP<EvalExpr> profile(dst, src);
      EvalExpr::exec(dst, src);
    }
    else
      Serial_dispatch_helper<Dim, DstBlock, SrcBlock, ProfileP, Rest>
	::exec(dst, src);
  }
};

/// In case the compile-time check fails, we continue the search
/// directly at the next entry in the type list.
template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  template <typename> class ProfileP,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename EvalExpr>
struct Serial_dispatch_helper<Dim, DstBlock, SrcBlock, ProfileP, TagList,
			      Tag, Rest, EvalExpr,
			      false>
  : Serial_dispatch_helper<Dim, DstBlock, SrcBlock, ProfileP, Rest>
{};

/// Terminator. Instead of passing on to the next element
/// it aborts the program. It is a program error to define
/// callback lists that can't handle a given expression.
template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  template <typename> class ProfileP,
	  typename TagList,
	  typename Tag,
	  typename EvalExpr>
struct Serial_dispatch_helper<Dim, DstBlock, SrcBlock, ProfileP, TagList,
			      Tag, None_type, EvalExpr, true>
{
  static void exec(DstBlock& dst, SrcBlock const& src)
    VSIP_NOTHROW
  {
    if (EvalExpr::rt_valid(dst, src))
    {
      ProfileP<EvalExpr> profile(dst, src);
      EvalExpr::exec(dst, src);
    }
    else assert(0);
  }
};



/// Front-end to Serial_dispatch_helper.  Uses S_d_h's ProfileP to
/// attach bits like profiling, coverage etc.

template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  typename TagList,
	  bool     Internal>
struct Serial_dispatch
{
  static void exec(DstBlock& dst, SrcBlock const& src)
    VSIP_NOTHROW
  {
    Serial_dispatch_helper<Dim, DstBlock, SrcBlock, Eval_sd_policy, TagList>::
      exec(dst, src);
  }
};

template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  typename TagList>
struct Serial_dispatch<Dim, DstBlock, SrcBlock, TagList, true>
{
  static void exec(DstBlock& dst, SrcBlock const& src)
    VSIP_NOTHROW
  {
    Serial_dispatch_helper<Dim, DstBlock, SrcBlock, Eval_int_sd_policy,
                           TagList>::
      exec(dst, src);
  }
};

} // namespace vsip::impl
} // namespace vsip

#endif
