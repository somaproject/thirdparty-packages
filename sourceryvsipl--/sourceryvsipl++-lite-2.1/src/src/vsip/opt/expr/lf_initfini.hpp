/* Copyright (c) 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/expr/lf_initfini.hpp
    @author  Jules Bergmann
    @date    2006-08-04
    @brief   VSIPL++ Library: Loop-fusion init/fini
*/

#ifndef VSIP_OPT_EXPR_LF_INITFINI_HPP
#define VSIP_OPT_EXPR_LF_INITFINI_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/expr/ternary_block.hpp>
#include <vsip/core/expr/vmmul_block.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/coverage.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

// fwd decl
template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
class Return_expr_block;



/// Reduction to count the number operations per point of an expression.

template <typename LeafT>
struct Apply_leaf
{
public:
  template <typename BlockT>
  struct transform
  {
    static void apply(BlockT const& block)
    {
      LeafT::template leaf_node<BlockT>::apply(block);
    }
  };

  template <typename BlockT>
  struct transform<BlockT const> : public transform<BlockT> {};

  template <dimension_type            Dim0,
	    template <typename> class Op,
	    typename                  BlockT,
	    typename                  T>
  struct transform<Unary_expr_block<Dim0, Op, BlockT, T> >
  {
    typedef Unary_expr_block<Dim0, Op, BlockT, T>
		block_type;

    static void apply(block_type const& block)
    {
      transform<BlockT>::apply(block.op());
    }
  };

  template <dimension_type                Dim0,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  struct transform<Binary_expr_block<Dim0, Op, LBlock, LType,
				     RBlock, RType> >
  {
    typedef Binary_expr_block<Dim0, Op, LBlock, LType, RBlock, RType> 
		block_type;

    static void apply(block_type const& block)
    {
      transform<LBlock>::apply(block.left());
      transform<RBlock>::apply(block.right());
    }
  };

  template <dimension_type                Dim0,
	    template <typename, typename, typename> class Op,
	    typename                      Block1,
	    typename                      Type1,
	    typename                      Block2,
	    typename                      Type2,
	    typename                      Block3,
	    typename                      Type3>
  struct transform<Ternary_expr_block<Dim0, Op, Block1, Type1,
				     Block2, Type2, Block3, Type3> >
  {
    typedef Ternary_expr_block<Dim0, Op, Block1, Type1,
			       Block2, Type2, Block3, Type3> const
		block_type;

    static void apply(block_type const& block)
    {
      transform<Block1>::apply(block.first());
      transform<Block2>::apply(block.second());
      transform<Block3>::apply(block.third());
    }
  };

  template <dimension_type VecDim,
	    typename       Block0,
	    typename       Block1>
  struct transform<Vmmul_expr_block<VecDim, Block0, Block1> >
  {
    typedef Vmmul_expr_block<VecDim, Block0, Block1> const block_type;

    static void apply(block_type const& block)
    {
      transform<Block0>::apply(block.get_vblk());
      transform<Block1>::apply(block.get_mblk());
    }
  };
};



struct Do_loop_fusion_init : Apply_leaf<Do_loop_fusion_init>
{
  template <typename BlockT>
  struct leaf_node
  {
    static void apply(BlockT const&) {} // no-op
  };

  template <dimension_type Dim,
	    typename       T,
	    typename       ReturnFunctor>
  struct leaf_node<Return_expr_block<Dim, T, ReturnFunctor> >
  {
    typedef Return_expr_block<Dim, T, ReturnFunctor>
		block_type;

    static void apply(block_type const& block) 
    {
      const_cast<block_type&>(block).loop_fusion_init();
    }
  };
};



struct Do_loop_fusion_fini : Apply_leaf<Do_loop_fusion_fini>
{
  template <typename BlockT>
  struct leaf_node
  {
    static void apply(BlockT const&) {} // no-op
  };

  template <dimension_type Dim,
	    typename       T,
	    typename       ReturnFunctor>
  struct leaf_node<Return_expr_block<Dim, T, ReturnFunctor> >
  {
    typedef Return_expr_block<Dim, T, ReturnFunctor>
		block_type;

    static void apply(block_type const& block) 
    {
      const_cast<block_type&>(block).loop_fusion_fini();
    }
  };
};


template <typename BlockT>
inline void
do_loop_fusion_init(BlockT const& block)
{
  Do_loop_fusion_init::transform<BlockT>::apply(block);
}



template <typename BlockT>
inline void
do_loop_fusion_fini(BlockT const& block)
{
  Do_loop_fusion_fini::transform<BlockT>::apply(block);
}



} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_EXPR_LF_INITFINI_HPP
