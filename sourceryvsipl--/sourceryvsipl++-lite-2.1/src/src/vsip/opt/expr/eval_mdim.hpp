/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/expr/eval_mdim.hpp
    @author  Jules Bergmann
    @date    2007-08-13
    @brief   VSIPL++ Library: Evaluate a multi-dimensional expression
                              as multiply vector expression.
*/

#ifndef VSIP_OPT_EXPR_EVAL_MDIM_HPP
#define VSIP_OPT_EXPR_EVAL_MDIM_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
// #include <vsip/core/extdata.hpp>
// #include <vsip/core/expr/scalar_block.hpp>
// #include <vsip/core/expr/unary_block.hpp>
// #include <vsip/core/expr/binary_block.hpp>
// #include <vsip/core/expr/ternary_block.hpp>
#include <vsip/core/coverage.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{



/***********************************************************************
  Expression Reductions
***********************************************************************/

// Expression reduction to determine if an expression contains
// difficult blocks to handle, such as Rbo blocks.

struct Reduce_is_expr_difficult
{
public:
  template <typename BlockT>
  struct leaf_node
  {
    typedef Bool_type<false> type;
  };

  template <dimension_type Dim,
	    typename       T,
	    typename       ReturnFunctor>
  struct leaf_node<Return_expr_block<Dim, T, ReturnFunctor> const>
  {
    typedef Bool_type<true> type;
  };

  template <typename BlockT>
  struct transform
  {
    typedef typename leaf_node<BlockT>::type type;
  };

  template <dimension_type            Dim,
	    template <typename> class Op,
	    typename                  BlockT,
	    typename                  T>
  struct transform<Unary_expr_block<Dim, Op, BlockT, T> const>
  {
    typedef typename transform<BlockT>::type type;
  };

  template <dimension_type                Dim0,
	    typename                      LBlock,
	    typename                      RBlock>
  struct transform<Vmmul_expr_block<Dim0, LBlock, RBlock> const>
  {
//    typedef Bool_type<true> type;
    typedef Bool_type<transform<LBlock>::type::value ||
                      transform<RBlock>::type::value> type;
  };

  template <dimension_type                Dim,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  struct transform<Binary_expr_block<Dim, Op, LBlock, LType,
				     RBlock, RType> const>
  {
    typedef Bool_type<transform<LBlock>::type::value ||
                      transform<RBlock>::type::value> type;
  };

  template <dimension_type                Dim,
	    template <typename, typename, typename> class Op,
	    typename                      Block1,
	    typename                      Type1,
	    typename                      Block2,
	    typename                      Type2,
	    typename                      Block3,
	    typename                      Type3>
  struct transform<Ternary_expr_block<Dim, Op, Block1, Type1,
				      Block2, Type2, Block3, Type3> const>
  {
    typedef Bool_type<transform<Block1>::type::value ||
                      transform<Block2>::type::value ||
                      transform<Block2>::type::value> type;
  };
};



// Reduction to extract a 1-dimensional subview from an expression
// with n-dimensional (where n > 1), taking care to push subview
// as "deep" as possible so that library evaluators can still be
// applied.

template <dimension_type FixedDim>
class Subdim_expr
{
public:
  template <typename BlockT>
  struct leaf_node
  {
    typedef Sliced_block<BlockT, FixedDim> type;
  };

  template <dimension_type Dim,
	    typename       T>
  struct leaf_node<Scalar_block<Dim, T> >
  {
    typedef Scalar_block<1, T> type;
  };

  template <dimension_type            Dim,
	    template <typename> class Op,
	    typename                  NewBlockT,
	    typename                  NewT>
  struct unary_node
  {
    typedef Unary_expr_block<1, Op, NewBlockT, NewT> const type;
  };

  template <dimension_type                Dim,
	    template <typename, typename> class Op,
	    typename                      NewLBlock,
	    typename                      NewLType,
	    typename                      NewRBlock,
	    typename                      NewRType>
  struct binary_node
  {
    typedef Binary_expr_block<1, Op,
			      NewLBlock, NewLType,
			      NewRBlock, NewRType> const type;
  };

  template <dimension_type                          Dim,
	    template <typename, typename, typename> class Op,
	    typename                                NewBlock1,
	    typename                                NewType1,
	    typename                                NewBlock2,
	    typename                                NewType2,
	    typename                                NewBlock3,
	    typename                                NewType3>
  struct ternary_node
  {
    typedef Ternary_expr_block<1, Op,
			       NewBlock1, NewType1,
			       NewBlock2, NewType2,
			       NewBlock3, NewType3> const type;
  };

  template <typename BlockT>
  struct transform
  {
    typedef typename leaf_node<BlockT>::type type;
  };

  template <dimension_type            Dim,
	    template <typename> class Op,
	    typename                  BlockT,
	    typename                  T>
  struct transform<Unary_expr_block<Dim, Op, BlockT, T> const>
  {
    typedef typename unary_node<Dim, Op,
				typename transform<BlockT>::type,
				T>::type type;
  };

  template <dimension_type                Dim,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  struct transform<Binary_expr_block<Dim, Op, LBlock, LType,
				     RBlock, RType> const>
  {
    typedef typename binary_node<Dim, Op,
				typename transform<LBlock>::type, LType,
				typename transform<RBlock>::type, RType>
				::type type;
  };

  template <dimension_type                Dim,
	    template <typename, typename, typename> class Op,
	    typename                      Block1,
	    typename                      Type1,
	    typename                      Block2,
	    typename                      Type2,
	    typename                      Block3,
	    typename                      Type3>
  struct transform<Ternary_expr_block<Dim, Op, Block1, Type1,
				      Block2, Type2, Block3, Type3> const>
  {
    typedef typename ternary_node<Dim, Op,
				typename transform<Block1>::type, Type1,
				typename transform<Block2>::type, Type2,
				typename transform<Block3>::type, Type3>
				::type type;
  };


  template <dimension_type            Dim,
	    template <typename> class Op,
	    typename                  BlockT,
	    typename                  T>
  typename transform<Unary_expr_block<Dim, Op, BlockT, T> const>::type
  apply(Unary_expr_block<Dim, Op, BlockT, T> const& blk)
  {
    typedef typename
      transform<Unary_expr_block<Dim, Op, BlockT, T> const>::type
        block_type;
    return block_type(apply(const_cast<BlockT&>(blk.op())), blk);
  }

  template <dimension_type                Dim,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  typename transform<Binary_expr_block<Dim, Op, LBlock, LType,
				       RBlock, RType> const>::type
  apply(Binary_expr_block<Dim, Op, LBlock, LType, RBlock, RType> const& blk)
  {
    typedef typename
      transform<Binary_expr_block<Dim, Op, LBlock, LType,
                                  RBlock, RType> const>::type
        block_type;
    return block_type(apply(const_cast<LBlock&>(blk.left())),
		      apply(const_cast<RBlock&>(blk.right())));
  }

  template <dimension_type                          Dim,
	    template <typename, typename, typename> class Op,
	    typename                                Block1,
	    typename                                Type1,
	    typename                                Block2,
	    typename                                Type2,
	    typename                                Block3,
	    typename                                Type3>
  typename transform<Ternary_expr_block<Dim, Op, Block1, Type1, Block2, Type2,
					Block3, Type3> const>::type
  apply(Ternary_expr_block<Dim, Op, Block1, Type1, Block2, Type2,
	Block3, Type3> const& blk)
  {
    typedef typename
      transform<Ternary_expr_block<Dim, Op, Block1, Type1,
                                   Block2, Type2, Block3, Type3> const>::type
        block_type;
    return block_type(apply(const_cast<Block1&>(blk.first ())),
		      apply(const_cast<Block2&>(blk.second())),
		      apply(const_cast<Block3&>(blk.third ())));
  }

  // Leaf combine function for Scalar_block.
  template <dimension_type Dim,
	    typename       T>
  typename transform<Scalar_block<Dim, T> >::type
  apply(Scalar_block<Dim, T>& block) const
  {
    typedef typename transform<Scalar_block<Dim, T> >::type type;
    return type(block.value());
  }


  // Leaf combine function.
  template <typename BlockT>
  typename transform<BlockT>::type
  apply(BlockT& block) const
  {
    typedef typename transform<BlockT>::type block_type;
    return block_type(block, idx_);
  }

  // Constructors.
public:
  Subdim_expr(index_type idx) : idx_(idx) {}

  index_type idx_;
};



// Evaluator to convert multi-dimensional expressions into multiple 1
// dimensional expressions.  Intended to handle non-dense multi-dim
// expressions that cannot be handled by Dense_expr_tag.

template <typename       DstBlock,
	  typename       SrcBlock>
struct Serial_expr_evaluator<2, DstBlock, SrcBlock, Mdim_expr_tag>
{
  static char const* name() { return "Expr_Mdim<2>"; }

  static bool const ct_valid =
    !Reduce_is_expr_difficult::template transform<SrcBlock>::type::value;

  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  { return true; }

  typedef typename Block_layout<DstBlock>::order_type order_type;
  static dimension_type const fixed_dim = 0;

  typedef typename Proper_type_of<SrcBlock>::type proper_src_block_type;

  typedef typename Subdim_expr<fixed_dim>::template
                   transform<proper_src_block_type>::type
	  new_src_type;
  typedef typename Subdim_expr<fixed_dim>::template transform<DstBlock>::type
	  new_dst_type;

  static new_dst_type diag_helper_dst(DstBlock& dst)
  {
    Subdim_expr<fixed_dim> reduce(0);
    typename View_block_storage<new_dst_type>::plain_type
      new_dst = reduce.apply(dst);
    return new_dst;
  }

  static new_src_type diag_helper_src(SrcBlock const& src)
  {
    Subdim_expr<fixed_dim> reduce(0);
    typename View_block_storage<new_src_type>::plain_type
      new_src = reduce.apply(const_cast<SrcBlock&>(src));
    return new_src;
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    for (index_type r=0; r<dst.size(2, fixed_dim); ++r)
    {
      Subdim_expr<fixed_dim> reduce(r);

      // Serial_dispatch_helper::exec takes the 'dst' block as
      // non-const reference.  We cannot pass redim.apply(...)
      // directly to exec() because if the result is a by-value block
      // (such as a Redim_block of a Sliced_block), it returns a
      // temporary object.

      typename View_block_storage<new_dst_type>::plain_type
	  new_dst = reduce.apply(const_cast<DstBlock&>(dst));

      Serial_dispatch<1, new_dst_type, new_src_type, LibraryTagList, true>
	::exec(new_dst, reduce.apply(const_cast<proper_src_block_type&>(src)));
    }
  }
};



} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_EXPR_EVAL_MDIM_HPP
