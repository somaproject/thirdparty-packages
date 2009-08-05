/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/expr/binary_block.hpp
    @author  Stefan Seefeld
    @date    2005-01-20
    @brief   VSIPL++ Library: Binary expression block class templates.

    This file defines the Binary_expr_block class templates.
*/

#ifndef VSIP_CORE_EXPR_BINARY_BLOCK_HPP
#define VSIP_CORE_EXPR_BINARY_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/noncopyable.hpp>
#include <vsip/core/metaprogramming.hpp>

namespace vsip
{
namespace impl
{

/***********************************************************************
  Declarations
***********************************************************************/

/// Expression template block for binary expressions.
///
/// Template parameters:
///   :D: to be a dimension with range 0 < D <= VSIP_MAX_DIMENSION
///   :Operator: to be a class template representing a binary
///              operation on operands of type LBlock::value_type and RBlock::value_type.
///   :LBlock: to be a Block.
///   :RBlock: to be a Block.
template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
class Binary_expr_block : private Operator<LType, RType>,
			  public Non_assignable
{
public:
  static dimension_type const dim = D;
  typedef typename Operator<LType, RType>::result_type value_type;

  typedef value_type& reference_type;
  typedef value_type const& const_reference_type;
  typedef typename LBlock::map_type map_type;

  Binary_expr_block(LBlock const& lhs, RBlock const& rhs);
  Binary_expr_block(LBlock const& lhs, RBlock const& rhs,
		    Operator<LType, RType> const& op);

  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type Dim, dimension_type d) const VSIP_NOTHROW;

  void increment_count() const VSIP_NOTHROW {}
  void decrement_count() const VSIP_NOTHROW {}
  map_type const& map() const VSIP_NOTHROW { return lhs_.map();}

  LBlock const& left() const VSIP_NOTHROW {return lhs_;}
  RBlock const& right() const VSIP_NOTHROW {return rhs_;}

  value_type get(index_type i) const;
  value_type get(index_type i, index_type j) const;
  value_type get(index_type i, index_type j, index_type k) const;

  // copy-constructor: default is OK.

private:
  typename View_block_storage<LBlock>::expr_type lhs_;
  typename View_block_storage<RBlock>::expr_type rhs_;
};



/// Specialize Is_expr_block for binary expr blocks.
template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
struct Is_expr_block<Binary_expr_block<D, Operator, LBlock, LType,
				       RBlock, RType> >
{ static bool const value = true; };



/// Specialize View_block_storage to control how views store binary
/// expression template blocks.
template <dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct View_block_storage<const Binary_expr_block<D, Operator,
					      LBlock, LType,
					      RBlock, RType> >
  : By_value_block_storage<const Binary_expr_block<D, Operator,
					      LBlock, LType,
					      RBlock, RType> >
{
};



template <dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct View_block_storage<Binary_expr_block<D, Operator,
					LBlock, LType,
					RBlock, RType> >
{
  // No typedef provided.  A non-const expresstion template block is
  // an error.
};



template <dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Distributed_local_block<
  Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> const>
{
  typedef Binary_expr_block<D, Operator,
			    typename Distributed_local_block<LBlock>::type,
			    LType,
			    typename Distributed_local_block<RBlock>::type,
			    RType> const
		type;
  typedef Binary_expr_block<D, Operator,
			typename Distributed_local_block<LBlock>::proxy_type,
			LType,
			typename Distributed_local_block<RBlock>::proxy_type,
			RType> const
		proxy_type;
};



template <dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Distributed_local_block<
  Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> >
{
  typedef Binary_expr_block<D, Operator,
			    typename Distributed_local_block<LBlock>::type,
			    LType,
			    typename Distributed_local_block<RBlock>::type,
			    RType>
		type;
  typedef Binary_expr_block<D, Operator,
			typename Distributed_local_block<LBlock>::proxy_type,
			LType,
			typename Distributed_local_block<RBlock>::proxy_type,
			RType>
		proxy_type;
};




template <dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
Binary_expr_block<D, Operator,
		  typename Distributed_local_block<LBlock>::type, LType,
		  typename Distributed_local_block<RBlock>::type, RType>
get_local_block(
  Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> const& block)
{
  typedef Binary_expr_block<D, Operator,
		  typename Distributed_local_block<LBlock>::type, LType,
		  typename Distributed_local_block<RBlock>::type, RType>
		  block_type;

  return block_type(get_local_block(block.left()),
		    get_local_block(block.right()));
}



template <typename                            CombineT,
	  dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Combine_return_type<CombineT,
			   Binary_expr_block<D, Operator, LBlock, LType,
					     RBlock, RType> const>
{
  typedef Binary_expr_block<D, Operator,
		typename Combine_return_type<CombineT, LBlock>::tree_type,
		LType,
		typename Combine_return_type<CombineT, RBlock>::tree_type,
		RType> const tree_type;
  typedef tree_type type;
};



template <typename                            CombineT,
	  dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Combine_return_type<CombineT,
			   Binary_expr_block<D, Operator, LBlock, LType,
					     RBlock, RType> >
: Combine_return_type<CombineT,
		      Binary_expr_block<D, Operator, LBlock, LType,
					RBlock, RType> const>
{};



template <typename                            CombineT,
	  dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
typename Combine_return_type<CombineT,
			     Binary_expr_block<D, Operator, LBlock, LType,
					       RBlock, RType> const>::type
apply_combine(
  CombineT const&                                                     combine,
  Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> const& block)
{
  typedef typename Combine_return_type<
    CombineT,
    Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> const>::type
		block_type;

  return block_type(apply_combine(combine, block.left()),
		    apply_combine(combine, block.right()));
}



template <typename                            VisitorT,
	  dimension_type                      D,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
void
apply_leaf(
  VisitorT const&                                                     visitor,
  Binary_expr_block<D, Operator, LBlock, LType, RBlock, RType> const& block)
{
  apply_leaf(visitor, block.left());
  apply_leaf(visitor, block.right());
}



template <dimension_type                      MapDim,
	  typename                            MapT,
	  dimension_type                      Dim,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Is_par_same_map<MapDim, MapT,
		       const Binary_expr_block<Dim, Operator,
					       LBlock, LType,
					       RBlock, RType> >
{
  typedef const Binary_expr_block<Dim, Operator,
				  LBlock, LType,
				  RBlock, RType> block_type;

  static bool value(MapT const& map, block_type& block)
  {
    return Is_par_same_map<MapDim, MapT, LBlock>::value(map, block.left()) &&
           Is_par_same_map<MapDim, MapT, RBlock>::value(map, block.right());
  }
};



template <dimension_type                      Dim,
	  template <typename, typename> class Operator,
	  typename                            LBlock,
	  typename                            LType,
	  typename                            RBlock,
	  typename                            RType>
struct Is_par_reorg_ok<Binary_expr_block<Dim, Operator,
				      LBlock, LType,
				      RBlock, RType> const>
{
  static bool const value = Is_par_reorg_ok<LBlock>::value &&
                            Is_par_reorg_ok<RBlock>::value;
};



/***********************************************************************
  Definitions
***********************************************************************/

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline
Binary_expr_block<D, Operator, 
		LBlock, LType,
		RBlock, RType>::Binary_expr_block
(LBlock const& lhs, RBlock const& rhs)
  : lhs_(lhs), rhs_(rhs) 
{
  assert(!Is_sized_block<LBlock>::value ||
	 !Is_sized_block<RBlock>::value ||
	 extent<D>(lhs_) == extent<D>(rhs_));
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline
Binary_expr_block<D, Operator, 
		LBlock, LType,
		RBlock, RType>::Binary_expr_block
(LBlock const& lhs, RBlock const& rhs, Operator<LType, RType> const& op)
  : Operator<LType, RType>(op), lhs_(lhs), rhs_(rhs) 
{
  assert(!Is_sized_block<LBlock>::value ||
	 !Is_sized_block<RBlock>::value ||
	 extent<D>(lhs_) == extent<D>(rhs_));
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline length_type
Binary_expr_block<D, Operator,
		LBlock, LType,
		RBlock, RType>::size() const VSIP_NOTHROW
{
  if (!Is_sized_block<LBlock>::value)
    return rhs_.size();
  else if (!Is_sized_block<RBlock>::value)
    return lhs_.size();
  else
  {
    assert(lhs_.size() == rhs_.size());
    return lhs_.size(); 
  }
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline length_type
Binary_expr_block<D, Operator,
		LBlock, LType,
		RBlock, RType>::size(dimension_type Dim,
				     dimension_type d) const VSIP_NOTHROW
{
  if (!Is_sized_block<LBlock>::value)
    return rhs_.size(Dim, d); 
  else if (!Is_sized_block<RBlock>::value)
    return lhs_.size(Dim, d); 
  else
  {
    assert(lhs_.size(Dim, d) == rhs_.size(Dim, d));
    return lhs_.size(Dim, d); 
  }
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline typename Binary_expr_block<D, Operator,
				LBlock, LType,
				RBlock, RType>::value_type
Binary_expr_block<D, Operator,
		LBlock, LType,
		RBlock, RType>::get(index_type i) const
{
  return (*this)(static_cast<LType>(this->left().get(i)),
		 static_cast<RType>(this->right().get(i)));
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline typename Binary_expr_block<D, Operator,
				LBlock, LType,
				RBlock, RType>::value_type
Binary_expr_block<D, Operator,
		LBlock, LType,
		RBlock, RType>::get(index_type i,
				    index_type j) const
{
  return (*this)(static_cast<LType>(this->left().get(i, j)),
		 static_cast<RType>(this->right().get(i, j)));
}

template <dimension_type D,
	  template <typename, typename> class Operator,
	  typename LBlock, typename LType,
	  typename RBlock, typename RType>
inline typename Binary_expr_block<D, Operator, 
				LBlock, LType,
				RBlock, RType>::value_type
Binary_expr_block<D, Operator, 
		LBlock, LType,
		RBlock, RType>::get(index_type i,
				    index_type j,
				    index_type k) const
{
  return (*this)(static_cast<LType>(this->left().get(i, j, k)),
		 static_cast<RType>(this->right().get(i, j, k)));
}

} // namespace vsip::impl
} // namespace vsip

#endif
