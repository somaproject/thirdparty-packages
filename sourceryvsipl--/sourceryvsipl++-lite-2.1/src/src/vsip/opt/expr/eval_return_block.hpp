/* Copyright (c) 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/expr/eval_return_block.hpp
    @author  Jules Bergmann
    @date    2006-09-02
    @brief   VSIPL++ Library: Evaluator for return-block optimization.

*/

#ifndef VSIP_OPT_EXPR_EVAL_RETURN_BLOCK_HPP
#define VSIP_OPT_EXPR_EVAL_RETURN_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/expr/return_block.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Evaluator for return expression block.

template <dimension_type Dim,
	  typename       DstBlock,
	  typename       T,
	  typename       ReturnFunctor>
struct Serial_expr_evaluator<Dim, DstBlock,
	const Return_expr_block<Dim, T, ReturnFunctor>,
	Rbo_expr_tag>
{
  static char const* name() { return "Rbo_expr_tag"; }

  typedef Return_expr_block<Dim, T, ReturnFunctor> SrcBlock;

  typedef typename DstBlock::value_type dst_type;
  typedef typename SrcBlock::value_type src_type;

  static bool const ct_valid = true;

  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  {
    return true;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    src.apply(dst);
  }
};

} // namespace vsip::impl

} // namespace vsip

#endif // VSIP_OPT_EXPR_EVAL_RETURN_BLOCK_HPP
