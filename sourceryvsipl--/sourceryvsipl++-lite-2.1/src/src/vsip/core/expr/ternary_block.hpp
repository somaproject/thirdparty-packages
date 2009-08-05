/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/expr/ternary_block.hpp
    @author  Stefan Seefeld
    @date    2005-04-25
    @brief   VSIPL++ Library: Ternary expression block class templates.

    This file defines the Ternary_expr_block class templates.
*/

#ifndef VSIP_CORE_EXPR_TERNARY_BLOCK_HPP
#define VSIP_CORE_EXPR_TERNARY_BLOCK_HPP

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

/// Expression template block for ternary expressions.
///
/// Template parameters:
///   :D: to be a dimension with range 0 < D <= VSIP_MAX_DIMENSION
///   :Functor: to be a class template representing a ternary
///             operation on operands of type Block1::value_type, Block2::value_type,
///             and Block3::value_type.
///   :Block1: to be a Block.
///   :Block2: to be a Block.
///   :Block3: to be a Block.
template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
class Ternary_expr_block : private Functor<Type1, Type2, Type3>,
			   public Non_assignable
{
public:
  static dimension_type const dim = D;
  typedef typename Functor<Type1, Type2, Type3>::result_type value_type;

  typedef value_type& reference_type;
  typedef value_type const& const_reference_type;
  typedef typename Block1::map_type map_type;

  Ternary_expr_block(Block1 const& first,
		     Block2 const& second,
		     Block3 const& third);
  Ternary_expr_block(Block1 const& first,
		     Block2 const& second,
		     Block3 const& third,
		     Functor<Type1, Type2, Type3> const&);

  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type Dim, dimension_type d) const VSIP_NOTHROW;

  void increment_count() const VSIP_NOTHROW {}
  void decrement_count() const VSIP_NOTHROW {}
  map_type const& map() const VSIP_NOTHROW { return first_.map();}

  Block1 const& first() const VSIP_NOTHROW {return first_;}
  Block2 const& second() const VSIP_NOTHROW {return second_;}
  Block3 const& third() const VSIP_NOTHROW {return third_;}

  value_type get(index_type i) const;
  value_type get(index_type i, index_type j) const;
  value_type get(index_type i, index_type j, index_type k) const;

  // copy-constructor: default is OK.

private:
  typename View_block_storage<Block1>::expr_type first_;
  typename View_block_storage<Block2>::expr_type second_;
  typename View_block_storage<Block3>::expr_type third_;
};



/// Specialize Is_expr_block for ternary expr blocks.
template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Is_expr_block<Ternary_expr_block<D, Functor,
					Block1, Type1,
					Block2, Type2,
					Block3, Type3> >
{ static bool const value = true; };



/// Specialize View_block_storage to control how views store binary
/// expression template blocks.
template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct View_block_storage<const Ternary_expr_block<D, Functor,
					       Block1, Type1,
					       Block2, Type2,
					       Block3, Type3> >
  : By_value_block_storage<const Ternary_expr_block<D, Functor,
					       Block1, Type1,
					       Block2, Type2,
					       Block3, Type3> >
{};

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct View_block_storage<Ternary_expr_block<D, Functor,
					 Block1, Type1,
					 Block2, Type2,
					 Block3, Type3> >
{
  // No typedef provided.  A non-const expresstion template block is
  // an error.
};



/***********************************************************************
  Parallel traits and functions
***********************************************************************/

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Distributed_local_block<
  Ternary_expr_block<D, Functor,
		     Block1, Type1,
		     Block2, Type2,
		     Block3, Type3> const>
{
  typedef Ternary_expr_block<D, Functor,
		typename Distributed_local_block<Block1>::type, Type1,
		typename Distributed_local_block<Block2>::type, Type2,
		typename Distributed_local_block<Block3>::type, Type3> const
		type;
  typedef Ternary_expr_block<D, Functor,
		typename Distributed_local_block<Block1>::proxy_type, Type1,
		typename Distributed_local_block<Block2>::proxy_type, Type2,
		typename Distributed_local_block<Block3>::proxy_type, Type3>
          const proxy_type;
};

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Distributed_local_block<
  Ternary_expr_block<D, Functor,
		     Block1, Type1,
		     Block2, Type2,
		     Block3, Type3> >
{
  typedef Ternary_expr_block<D, Functor,
		typename Distributed_local_block<Block1>::type, Type1,
		typename Distributed_local_block<Block2>::type, Type2,
		typename Distributed_local_block<Block3>::type, Type3>
		type;
};



template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
Ternary_expr_block<D, Functor,
		typename Distributed_local_block<Block1>::type, Type1,
		typename Distributed_local_block<Block2>::type, Type2,
		typename Distributed_local_block<Block3>::type, Type3>
get_local_block(
  Ternary_expr_block<D, Functor,
		     Block1, Type1,
		     Block2, Type2,
		     Block3, Type3> const& block)
{
  typedef Ternary_expr_block<D, Functor,
		typename Distributed_local_block<Block1>::type, Type1,
		typename Distributed_local_block<Block2>::type, Type2,
		typename Distributed_local_block<Block3>::type, Type3>
		  block_type;

  return block_type(get_local_block(block.first()),
		    get_local_block(block.second()),
		    get_local_block(block.third()));
}



template <typename                  CombineT,
	  dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Combine_return_type<CombineT,
			   Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> const>
{
  typedef Ternary_expr_block<D, Functor,
	typename Combine_return_type<CombineT, Block1>::tree_type, Type1,
	typename Combine_return_type<CombineT, Block2>::tree_type, Type2,
	typename Combine_return_type<CombineT, Block3>::tree_type, Type3>
	const tree_type;
  typedef tree_type type;
};



template <typename                  CombineT,
	  dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Combine_return_type<CombineT,
			   Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> >
  : Combine_return_type<CombineT,
			Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> const>
{};



template <typename                  CombineT,
	  dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
typename Combine_return_type<CombineT,
			   Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> const>::type
apply_combine(
  CombineT const&                               combine,
  Ternary_expr_block<D, Functor,
                     Block1, Type1,
                     Block2, Type2,
                     Block3, Type3> const& block)
{
  typedef typename Combine_return_type<
    CombineT,
    Ternary_expr_block<D, Functor,
      Block1, Type1,
      Block2, Type2,
      Block3, Type3> const>::type
		block_type;

  return block_type(apply_combine(combine, block.first()),
		    apply_combine(combine, block.second()),
		    apply_combine(combine, block.third()));
}



template <typename                  VisitorT,
	  dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
void
apply_leaf(
  VisitorT const&                               visitor,
  Ternary_expr_block<D, Functor,
                     Block1, Type1,
                     Block2, Type2,
                     Block3, Type3> const& block)
{
  apply_leaf(visitor, block.first());
  apply_leaf(visitor, block.second());
  apply_leaf(visitor, block.third());
}



template <dimension_type MapDim,
	  typename MapT,
	  dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Is_par_same_map<MapDim, MapT,
		       Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> const>
{
  typedef Ternary_expr_block<D, Functor,
			     Block1, Type1,
			     Block2, Type2,
			     Block3, Type3> const
		block_type;

  static bool value(MapT const& map, block_type& block)
  {
    return Is_par_same_map<MapDim, MapT, Block1>::value(map, block.first()) &&
           Is_par_same_map<MapDim, MapT, Block2>::value(map, block.second()) &&
           Is_par_same_map<MapDim, MapT, Block3>::value(map, block.third());
  }
};



template <dimension_type            D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
struct Is_par_reorg_ok<Ternary_expr_block<D, Functor,
					      Block1, Type1,
					      Block2, Type2,
					      Block3, Type3> const>
{
  static bool const value = Is_par_reorg_ok<Block1>::value &&
                            Is_par_reorg_ok<Block2>::value &&
                            Is_par_reorg_ok<Block3>::value;
};



/***********************************************************************
  Definitions
***********************************************************************/

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline
Ternary_expr_block<D, Functor, 
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::Ternary_expr_block
(Block1 const& first, Block2 const& second, Block3 const& third)
  : first_(first), second_(second), third_(third)
{
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline
Ternary_expr_block<D, Functor, 
		   Block1, Type1,
		   Block2, Type2,
		   Block3, Type3>::Ternary_expr_block
(Block1 const& first, Block2 const& second, Block3 const& third,
 Functor<Type1, Type2, Type3> const& f)
  : Functor<Type1, Type2, Type3>(f), 
    first_(first), second_(second), third_(third)
{
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline length_type
Ternary_expr_block<D, Functor,
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::size() const VSIP_NOTHROW
{
  if (Is_sized_block<Block1>::value)
  {
    assert(!Is_sized_block<Block2>::value || first_.size() == second_.size());
    assert(!Is_sized_block<Block3>::value || first_.size() == third_.size());
    return first_.size();
  }
  else if (Is_sized_block<Block2>::value)
  {
    assert(!Is_sized_block<Block3>::value || second_.size() == third_.size());
    return second_.size();
  }
  else
  {
    return third_.size();
  }
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline length_type
Ternary_expr_block<D, Functor,
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::size(dimension_type Dim,
				      dimension_type d) const VSIP_NOTHROW
{
  if (Is_sized_block<Block1>::value)
  {
    assert(!Is_sized_block<Block2>::value || first_.size(Dim, d) ==
	                                     second_.size(Dim, d));
    assert(!Is_sized_block<Block3>::value || first_.size(Dim, d) ==
		                             third_.size(Dim, d));
    return first_.size(Dim, d);
  }
  else if (Is_sized_block<Block2>::value)
  {
    assert(!Is_sized_block<Block3>::value || second_.size(Dim, d) ==
                                             third_.size(Dim, d));
    return second_.size(Dim, d);
  }
  else
  {
    return third_.size(Dim, d);
  }
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline typename Ternary_expr_block<D, Functor,
				 Block1, Type1,
				 Block2, Type2,
				 Block3, Type3>::value_type
Ternary_expr_block<D, Functor,
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::get(index_type i) const
{
  return (*this)(static_cast<Type1>(this->first().get(i)),
		 static_cast<Type2>(this->second().get(i)),
		 static_cast<Type3>(this->third().get(i)));
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline typename Ternary_expr_block<D, Functor,
				 Block1, Type1,
				 Block2, Type2,
				 Block3, Type3>::value_type
Ternary_expr_block<D, Functor,
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::get(index_type i,
				     index_type j) const
{
  return (*this)(static_cast<Type1>(this->first().get(i, j)),
		 static_cast<Type2>(this->second().get(i, j)),
		 static_cast<Type3>(this->third().get(i, j)));
}

template <dimension_type D,
	  template <typename, typename, typename> class Functor,
	  typename Block1, typename Type1,
	  typename Block2, typename Type2,
	  typename Block3, typename Type3>
inline typename Ternary_expr_block<D, Functor,
				 Block1, Type1,
				 Block2, Type2,
				 Block3, Type3>::value_type
Ternary_expr_block<D, Functor,
		 Block1, Type1,
		 Block2, Type2,
		 Block3, Type3>::get(index_type i,
				     index_type j,
				     index_type k) const
{
  return (*this)(static_cast<Type1>(this->first().get(i, j, k)),
		 static_cast<Type2>(this->second().get(i, j, k)),
		 static_cast<Type3>(this->third().get(i, j, k)));
}


} // namespace vsip::impl
} // namespace vsip

#endif
