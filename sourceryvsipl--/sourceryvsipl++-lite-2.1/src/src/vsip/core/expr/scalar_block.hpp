/* Copyright (c) 2005, 2006, 2007 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/expr/scalar_block.hpp
    @author  Stefan Seefeld
    @date    2005-01-20
    @brief   VSIPL++ Library: Scalar block class template.

    This file defines the Scalar_block class templates.
*/

#ifndef VSIP_CORE_EXPR_SCALAR_BLOCK_HPP
#define VSIP_CORE_EXPR_SCALAR_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/noncopyable.hpp>
#include <vsip/core/length.hpp>
#include <vsip/core/map_fwd.hpp>

namespace vsip
{
namespace impl
{

/***********************************************************************
  Declarations
***********************************************************************/

template <dimension_type D>
struct Scalar_block_shared_map
{
  // typedef Local_or_global_map<D> type;
  typedef Scalar_block_map<D> type;
  static type map;
};

// These are explicitly instantiated in scalar_block.cpp
#if defined (__ghs__)
#pragma do_not_instantiate Scalar_block_shared_map<1>::map
#pragma do_not_instantiate Scalar_block_shared_map<2>::map
#pragma do_not_instantiate Scalar_block_shared_map<3>::map
#endif


/// An adapter presenting a scalar as a block. This is useful when constructing
/// Binary_expr_block objects (which expect two block operands) taking a block and
/// a scalar.
///
/// Template parameters:
///   :D: to be a dimension with range 0 < D <= VSIP_MAX_DIMENSION
///   :Scalar: to be a builtin scalar type.
template <dimension_type D, typename Scalar>
class Scalar_block_base : public Non_assignable
{
public:
  typedef Scalar value_type;
  typedef value_type& reference_type;
  typedef value_type const& const_reference_type;
  typedef typename Scalar_block_shared_map<D>::type map_type;

  static dimension_type const dim = D;

  Scalar_block_base(Scalar s) : value_(s) {}
#if (defined(__GNUC__) && __GNUC__ < 4)
  // GCC 3.4.4 appears to over-optimize multiple scalar values on
  // stack when optimization & strong inlining are enabled, causing
  // threshold.cpp and other tests to fail.  (070618)
  Scalar_block_base(Scalar_block_base const& b) : value_(b.value_) {}
#endif

  void increment_count() const VSIP_NOTHROW {}
  void decrement_count() const VSIP_NOTHROW {}
  map_type const& map() const VSIP_NOTHROW
  { return Scalar_block_shared_map<D>::map; }

  Scalar value() const VSIP_NOTHROW {return value_;}

private:
  Scalar const         value_;
};

template <dimension_type D, typename Scalar>
class Scalar_block;

/// Scalar_block specialization for 1-dimension.
template <typename Scalar>
class Scalar_block<1, Scalar> : public Scalar_block_base<1, Scalar>
{
public:
  Scalar_block(Scalar s)
    : Scalar_block_base<1, Scalar>(s) {}
  length_type size() const VSIP_NOTHROW { return 0; }
  length_type size(dimension_type block_dim, dimension_type d) const VSIP_NOTHROW;

  Scalar get(index_type idx) const VSIP_NOTHROW;

  // No member data.
};

/// Scalar_block specialization for 2-dimension.
template <typename Scalar>
class Scalar_block<2, Scalar> : public Scalar_block_base<2, Scalar>
{
public:
  Scalar_block(Scalar s)
    : Scalar_block_base<2, Scalar>(s) {}

  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type block_dim, dimension_type d) const VSIP_NOTHROW;

  Scalar get(index_type idx) const VSIP_NOTHROW;
  Scalar get(index_type x, index_type y) const VSIP_NOTHROW;

  // No member data.
};

/// Scalar_block specialization for 3-dimension.
template <typename Scalar>
class Scalar_block<3, Scalar> : public Scalar_block_base<3, Scalar>
{
public:
  Scalar_block(Scalar s)
    : Scalar_block_base<3, Scalar>(s) {}

  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type block_dim, dimension_type d) const VSIP_NOTHROW;

  Scalar get(index_type idx) const VSIP_NOTHROW;
  Scalar get(index_type x, index_type y, index_type z) const VSIP_NOTHROW;

  // No member data.
};



/// Specialize Is_expr_block for scalar expr blocks.
template <dimension_type D,
	  typename       Scalar>
struct Is_expr_block<Scalar_block<D, Scalar> >
{ static bool const value = true; };



template <dimension_type D,
	  typename       Scalar>
struct Is_sized_block<Scalar_block<D, Scalar> >
{ static bool const value = false; };



template <dimension_type D,
	  typename       Scalar>
struct Is_scalar_block<Scalar_block<D, Scalar> >
{ static bool const value = true; };



template <dimension_type D,
	  typename       Scalar>
struct Distributed_local_block<Scalar_block<D, Scalar> const>
{
  typedef Scalar_block<D, Scalar> const type;
  typedef Scalar_block<D, Scalar> const proxy_type;
};

template <dimension_type D,
	  typename       Scalar>
struct Distributed_local_block<Scalar_block<D, Scalar> >
{
  typedef Scalar_block<D, Scalar> type;
  typedef Scalar_block<D, Scalar> proxy_type;
};



template <dimension_type D,
	  typename       Scalar>
Scalar_block<D, Scalar>
get_local_block(
  Scalar_block<D, Scalar> const& block)
{
  return block;
}



template <typename       CombineT,
	  dimension_type D,
	  typename       Scalar>
struct Combine_return_type<CombineT, Scalar_block<D, Scalar> const>
{
  typedef Scalar_block<D, Scalar> block_type;
  typedef typename CombineT::template return_type<block_type>::type
		type;
  typedef typename CombineT::template tree_type<block_type>::type
		tree_type;
};



template <typename       CombineT,
	  dimension_type D,
	  typename       Scalar>
struct Combine_return_type<CombineT, Scalar_block<D, Scalar> >
  : Combine_return_type<CombineT, Scalar_block<D, Scalar> const>
{};



template <typename       CombineT,
	  dimension_type D,
	  typename       Scalar>
typename Combine_return_type<CombineT,
			     Scalar_block<D, Scalar> const>::type
apply_combine(
  CombineT const&                combine,
  Scalar_block<D, Scalar> const& block)
{
  return combine.apply(block);
}



template <typename       VisitorT,
	  dimension_type D,
	  typename       Scalar>
void
apply_leaf(
  VisitorT const&                visitor,
  Scalar_block<D, Scalar> const& block)
{
  visitor.apply(block);
}



template <dimension_type MapDim,
	  typename       MapT,
	  dimension_type D,
	  typename       Scalar>
struct Is_par_same_map<MapDim, MapT,
		       Scalar_block<D, Scalar> >
{
  typedef Scalar_block<D, Scalar> block_type;

  static bool value(MapT const&, block_type const&)
  {
    return true;
  }
};


// Default Is_par_reorg_ok is OK.



/// Assert that subblock is local to block (overload).
template <dimension_type D,
	  typename       Scalar>
void
assert_local(
  Scalar_block<D, Scalar> const& block,
  index_type                     sb)
{
  // Scalar_block is always valid locally.
}



template <dimension_type D, typename ScalarT>
struct Choose_peb<Scalar_block<D, ScalarT> >
{ typedef Peb_reuse_tag type; };



/***********************************************************************
  Definitions
***********************************************************************/

template <typename Scalar>
inline length_type
Scalar_block<1, Scalar>::size(dimension_type block_dim, dimension_type d) const
  VSIP_NOTHROW
{
  assert(block_dim == 1);
  assert(d == 0);
  return 0;
}

template <typename Scalar>
inline Scalar 
Scalar_block<1, Scalar>::get(index_type) const VSIP_NOTHROW
{
  return this->value();
}

template <typename Scalar>
inline length_type
Scalar_block<2, Scalar>::size() const VSIP_NOTHROW
{
  return 0;
}

template <typename Scalar>
inline length_type
Scalar_block<2, Scalar>::size(dimension_type block_dim, dimension_type d) const
  VSIP_NOTHROW
{
  assert((block_dim == 1 || block_dim == 2) && d < block_dim);
  return 0;
}

template <typename Scalar>
inline Scalar 
Scalar_block<2, Scalar>::get(index_type idx) const VSIP_NOTHROW
{
  assert(idx < size());
  return this->value();
}

template <typename Scalar>
inline Scalar 
Scalar_block<2, Scalar>::get(index_type, index_type) const VSIP_NOTHROW
{
  return this->value();
}



template <typename Scalar>
inline length_type
Scalar_block<3, Scalar>::size() const VSIP_NOTHROW
{
  return 0;
}

template <typename Scalar>
inline length_type
Scalar_block<3, Scalar>::size(
  dimension_type block_dim,
  dimension_type d) const
  VSIP_NOTHROW
{
  assert((block_dim == 1 || block_dim == 3) && d < block_dim);
  return 0;
}

template <typename Scalar>
inline Scalar
Scalar_block<3, Scalar>::get(index_type idx) const VSIP_NOTHROW
{
  assert(idx < size());
  return this->value();
}

template <typename Scalar>
inline Scalar
Scalar_block<3, Scalar>::get(index_type, index_type, index_type)
  const VSIP_NOTHROW
{
  return this->value();
}



/// Store Scalar_blocks by-value.
template <dimension_type D, typename Scalar>
struct View_block_storage<Scalar_block<D, Scalar> >
  : By_value_block_storage<Scalar_block<D, Scalar> >
{};

template <dimension_type D, typename Scalar>
struct Expr_block_storage<Scalar_block<D, Scalar> >
{
  typedef Scalar_block<D, Scalar> type;
};

} // namespace vsip::impl
} // namespace vsip

#endif
