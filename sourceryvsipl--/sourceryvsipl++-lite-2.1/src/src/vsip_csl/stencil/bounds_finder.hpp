/* Copyright (c) 2007 by CodeSourcery, Inc.  All rights reserved. */

/** @file    vsip_csl/stencil/bounds_finder.hpp
    @author  Stefan Seefeld
    @date    2007-05-07
    @brief   VSIPL++ Library: Find kernel bounds.

*/

#ifndef VSIP_CSL_STENCIL_BOUNDS_FINDER_HPP
#define VSIP_CSL_STENCIL_BOUNDS_FINDER_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip_csl/stencil/expr.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip_csl
{
namespace stencil
{

//
// Bound finder: Determine bounds for the stencil kernel.
//
struct Bounds
{
  Bounds() : x_prev(0), x_next(0), y_prev(0), y_next(0) {}
  vsip::length_type x_prev, x_next, y_prev, y_next;
};

template <typename E>
struct Bounds_finder
{
  // Traverse the expression E and extract the bounds of the kernel
  // to be constructed from it.
  static void apply(E e, Bounds& b);
};

template <typename S>
void find_bounds(S, Bounds &) {}

template <typename L, typename R, template <typename, typename> class O>
void find_bounds(Binary_expr<L, R, O> e, Bounds &b)
{
  Bounds_finder<L>::apply(e.left(), b);
  Bounds_finder<R>::apply(e.right(), b);
}

template <typename V, typename I, typename J>
void find_bounds(Call<V, I, J> c, Bounds& b)
{
  Bounds_finder<I>::apply(c.i(), b);
  Bounds_finder<J>::apply(c.j(), b);
}

template <vsip::dimension_type D>
void find_bounds(Iterator<D>, Bounds&) {}

void find_bounds(Offset<0> o, Bounds& b)
{
  if (o.i > 0 && o.i > b.y_next) b.y_next = o.i;
  else if (o.i < 0 && -o.i > b.y_prev) b.y_prev = -o.i;
}

void find_bounds(Offset<1> o, Bounds& b)
{
  if (o.i > 0 && o.i > b.x_next) b.x_next = o.i;
  else if (o.i < 0 && -o.i > b.x_prev) b.x_prev = -o.i;
}

template <typename E>
void 
Bounds_finder<E>::apply(E e, Bounds &b)
{ find_bounds(e, b);}

} // namespace vsip_csl::stencil
} // namespace vsip_csl

#endif
