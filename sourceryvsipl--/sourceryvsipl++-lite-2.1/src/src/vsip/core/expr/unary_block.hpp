/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/expr/unary_block.hpp
    @author  Stefan Seefeld
    @date    2005-01-20
    @brief   VSIPL++ Library: Unary expression block class templates.

    This file defines the Unary_expr_block class templates.
*/

#ifndef VSIP_CORE_EXPR_UNARY_BLOCK_HPP
#define VSIP_CORE_EXPR_UNARY_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/noncopyable.hpp>

namespace vsip
{
namespace impl
{

/***********************************************************************
  Declarations
***********************************************************************/

/// Expression template block for unary expressions.
///
/// Template parameters:
///   :D: to be a dimension with range 0 < D <= VSIP_MAX_DIMENSION
///   :Operator: to be a class template representing a unary
///              operation on operand of type Operand::value_type.
///   :Block: to be a Block.
template <dimension_type D,
	  template <typename> class Operator,
	  typename Block,
	  typename Type>
class Unary_expr_block : public Operator<Type>,
			 public Non_assignable
{
public:
  static dimension_type const dim = D;
  typedef typename Operator<Type>::result_type value_type;

  typedef value_type& reference_type;
  typedef value_type const& const_reference_type;
  typedef typename Block::map_type map_type;

  Unary_expr_block(Block const& block) : block_(block) {}
  Unary_expr_block(Block const& block, Operator<Type> const& op)
    : Operator<Type>(op), block_(block) {}

  length_type size() const VSIP_NOTHROW { return block_.size();}
  length_type size(dimension_type block_dim, dimension_type d) const VSIP_NOTHROW;

  void increment_count() const VSIP_NOTHROW {}
  void decrement_count() const VSIP_NOTHROW {}
  map_type const& map() const VSIP_NOTHROW { return block_.map();}

  Block const& op() const VSIP_NOTHROW {return block_;}

  value_type get(index_type i) const;
  value_type get(index_type i, index_type j) const;
  value_type get(index_type i, index_type j, index_type k) const;

  // copy-constructor: default is OK.

private:
  typename View_block_storage<Block>::expr_type block_;
};



/// Specialize Is_expr_block for unary expr blocks.
template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Is_expr_block<Unary_expr_block<Dim, Op, Block, Type> >
{ static bool const value = true; };



/// Specialize View_block_storage to control how views store binary
/// expression template blocks.
template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct View_block_storage<const Unary_expr_block<Dim, Op, Block, Type> >
  : By_value_block_storage<const Unary_expr_block<Dim, Op, Block, Type> >
{};



template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct View_block_storage<Unary_expr_block<Dim, Op, Block, Type> >
{
  // No typedef provided, non-const Unar_expr_block should be an error.
};



/***********************************************************************
  Parallel traits and functions
***********************************************************************/

template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Distributed_local_block<
  Unary_expr_block<Dim, Op, Block, Type> const>
{
  typedef Unary_expr_block<Dim, Op,
			   typename Distributed_local_block<Block>::type,
			   Type> const
		type;
  typedef Unary_expr_block<Dim, Op,
			   typename Distributed_local_block<Block>::proxy_type,
			   Type> const
		proxy_type;
};

template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Distributed_local_block<
  Unary_expr_block<Dim, Op, Block, Type> >
{
  typedef Unary_expr_block<Dim, Op,
			   typename Distributed_local_block<Block>::type,
			   Type> 
		type;
  typedef Unary_expr_block<Dim, Op,
			   typename Distributed_local_block<Block>::proxy_type,
			   Type> 
		proxy_type;
};



template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
Unary_expr_block<Dim, Op,
		 typename Distributed_local_block<Block>::type, Type>
get_local_block(
  Unary_expr_block<Dim, Op, Block, Type> const& block)
{
  typedef Unary_expr_block<Dim, Op,
		  typename Distributed_local_block<Block>::type, Type>
		  block_type;

  return block_type(get_local_block(block.op()));
}



template <typename                  CombineT,
	  dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Combine_return_type<CombineT,
			   Unary_expr_block<Dim, Op, Block, Type> const>
{
  typedef Unary_expr_block<Dim, Op,
		typename Combine_return_type<CombineT, Block>::tree_type,
		Type> const tree_type;
  typedef tree_type type;
};



template <typename                  CombineT,
	  dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Combine_return_type<CombineT,
			   Unary_expr_block<Dim, Op, Block, Type> >
  : Combine_return_type<CombineT,
			Unary_expr_block<Dim, Op, Block, Type> const>
{};



template <typename                  CombineT,
	  dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
typename Combine_return_type<CombineT,
			     Unary_expr_block<Dim, Op,
					      Block, Type> const>::type
apply_combine(
  CombineT const&                               combine,
  Unary_expr_block<Dim, Op, Block, Type> const& block)
{
  typedef typename Combine_return_type<
    CombineT,
    Unary_expr_block<Dim, Op, Block, Type> const>::type
		block_type;

  return block_type(apply_combine(combine, block.op()));
}



template <typename                  VisitorT,
	  dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
void
apply_leaf(
  VisitorT const&                               visitor,
  Unary_expr_block<Dim, Op, Block, Type> const& block)
{
  apply_leaf(visitor, block.op());
}



template <dimension_type            MapDim,
	  typename                  MapT,
	  dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Is_par_same_map<MapDim, MapT,
		       const Unary_expr_block<Dim, Op,
					      Block, Type> >
{
  typedef const Unary_expr_block<Dim, Op,
				 Block, Type> block_type;

  static bool value(MapT const& map, block_type& block)
  {
    return Is_par_same_map<MapDim, MapT, Block>::value(map, block.op());
  }
};



template <dimension_type            Dim,
	  template <typename> class Op,
	  typename                  Block,
	  typename                  Type>
struct Is_par_reorg_ok<Unary_expr_block<Dim, Op, Block, Type> const>
{
  static bool const value = Is_par_reorg_ok<Block>::value;
};



/***********************************************************************
  Definitions
***********************************************************************/

template <dimension_type D,
	  template <typename> class Operator,
	  typename Block,
	  typename Type>
inline length_type
Unary_expr_block<D, Operator, Block, Type>::size(dimension_type block_dim,
					       dimension_type d) const
  VSIP_NOTHROW
{
  return block_.size(block_dim, d);
}

template <dimension_type D,
	  template <typename> class Operator,
	  typename Block,
	  typename Type>
inline typename Unary_expr_block<D, Operator, Block, Type>::value_type
Unary_expr_block<D, Operator, Block, Type>::get(index_type i) const
{
  return (*this)(static_cast<Type>(this->op().get(i)));
}

template <dimension_type D,
	  template <typename> class Operator,
	  typename Block,
	  typename Type>
inline typename Unary_expr_block<D, Operator, Block, Type>::value_type
Unary_expr_block<D, Operator, Block, Type>::get(index_type i,
					      index_type j) const
{
  return (*this)(static_cast<Type>(this->op().get(i, j)));
}

template <dimension_type D,
	  template <typename> class Operator,
	  typename Block,
	  typename Type>
inline typename Unary_expr_block<D, Operator, Block, Type>::value_type
Unary_expr_block<D, Operator, Block, Type>::get(index_type i,
					      index_type j,
					      index_type k) const
{
  return (*this)(static_cast<Type>(this->op().get(i, j, k)));
}

} // namespace vsip::impl
} // namespace vsip

#endif
