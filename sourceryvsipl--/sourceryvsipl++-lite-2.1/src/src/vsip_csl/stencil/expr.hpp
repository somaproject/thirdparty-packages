/* Copyright (c) 2007 by CodeSourcery, Inc.  All rights reserved. */

/** @file    vsip_csl/stencil/expr.hpp
    @author  Stefan Seefeld
    @date    2007-05-07
    @brief   VSIPL++ Library: Harness to process stencil expressions.

*/

#ifndef VSIP_CSL_STENCIL_EXPR_HPP
#define VSIP_CSL_STENCIL_EXPR_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/expr/operations.hpp>
#include <cassert>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip_csl
{
namespace stencil
{

// Different sub-expressions we need to recognize:
//
//  o offset (I + 1)
//  o call (image(I))
//  o binary (2 * image(I))

template <typename T> struct Is_expr { static bool const value = false;};

template <vsip::dimension_type D>
struct Iterator {};

template <vsip::dimension_type D>
struct Offset
{
  Offset(vsip::stride_type ii) : i(ii) {}
  vsip::stride_type i;
};

template <typename V, typename I, typename J> struct Call;

template <typename V, typename I, typename J>
struct Is_expr<Call<V, I, J> > { static bool const value = true;};

template <typename V, typename I, typename J, typename RHS>
void assign(Call<V, I, J> lhs, RHS rhs);

template <typename V, typename I, typename J>
class Call
{
public:
  typedef V view_type;

  Call(V v, I i, J j) : view_(v), i_(i), j_(j) {}
  view_type view() { return view_;}
  I i() { return i_;}
  J j() { return j_;}

private:
  view_type view_;
  I i_;
  J j_;
};

template <typename V>
class Call<V, Iterator<0>, Iterator<1> >
{
public:
  typedef V view_type;

  Call(V v, Iterator<0>, Iterator<1>) : view_(v) {}
  template <typename RHS>
  void operator= (RHS rhs) { assign(*this, rhs);}
  view_type view() { return view_;}
  Iterator<0> i() { return Iterator<0>();}
  Iterator<1> j() { return Iterator<1>();}

private:
  view_type view_;
};

template <vsip::dimension_type D>
Offset<D> operator+ (Iterator<D>, vsip::stride_type i)
{
  return Offset<D>(i);
}

template <vsip::dimension_type D>
Offset<D> operator- (Iterator<D>, vsip::stride_type i)
{
  return Offset<D>(-i);
}

// Generic binary expression
template <typename L, typename R, template <typename, typename> class O>
class Binary_expr
{
public:
  typedef typename L::view_type view_type;

  Binary_expr(L l, R r) : left_(l), right_(r) {}

  view_type view() { return left_.view();}
  L left() { return left_;}
  R right() { return right_;}

private:
  L left_;
  R right_;
};

template <typename L, typename R, template <typename, typename> class O>
struct Is_expr<Binary_expr<L, R, O> > { static bool const value = true;};

// Binary expression: the left operand is an expression, the right a scalar.
// This is valid only for O in (Mult, Div)
template <typename L, typename R, template <typename, typename> class O,
          typename S>
class Binary_expr<Binary_expr<L, R, O>, S, vsip::impl::op::Mult>
{
public:
  typedef typename Binary_expr<L, R, O>::view_type view_type;

  Binary_expr(Binary_expr<L, R, O> e, S s) : left_(e), right_(s) {}
  view_type view() { return left_.view();}
  Binary_expr<L, R, O> left() { return left_;}
  S right() { return right_;}

private:
  Binary_expr<L, R, O> left_;
  S right_;
};

template <typename L, typename R, template <typename, typename> class O,
          typename S>
class Binary_expr<Binary_expr<L, R, O>, S, vsip::impl::op::Div>
{
public:
  typedef typename Binary_expr<L, R, O>::view_type view_type;

  Binary_expr(Binary_expr<L, R, O> e, S s) : left_(e), right_(s) {}
  view_type view() { return right_.view();}
  Binary_expr<L, R, O> left() { return left_;}
  S right() { return right_;}

private:
  Binary_expr<L, R, O> left_;
  S right_;
};

// Binary expression: the left operand is a scalar, the right an expression.
// This is valid only for O == Mult.
template <typename S,
          typename L, typename R, template <typename, typename> class O>
class Binary_expr<S, Binary_expr<L, R, O>, vsip::impl::op::Mult>
{
public:
  typedef typename Binary_expr<L, R, O>::view_type view_type;

  Binary_expr(S s, Binary_expr<L, R, O> e) : left_(s), right_(e) {}
  view_type view() { return right_.view();}
  S left() { return left_;}
  Binary_expr<L, R, O> right() { return right_;}

private:
  S left_;
  Binary_expr<L, R, O> right_;
};

// Binary expression: the left operand is a call, the right a scalar.
// This is valid only for O in (Mult, Div)
template <typename V, typename I, typename J, typename S>
class Binary_expr<Call<V, I, J>, S, vsip::impl::op::Mult>
{
public:
  typedef V view_type;

  Binary_expr(Call<V, I, J> c, S s) : left_(c), right_(s) {}
  view_type view() { return left_.view();}
  Call<V, I, J> left() { return left_;}
  S right() { return right_;}

private:
  Call<V, I, J> left_;
  S right_;
};

template <typename V, typename I, typename J, typename S>
class Binary_expr<Call<V, I, J>, S, vsip::impl::op::Div>
{
public:
  typedef V view_type;

  Binary_expr(Call<V, I, J> c, S s) : left_(c), right_(s) {}
  view_type view() { return left_.view();}
  Call<V, I, J> left() { return left_;}
  S right() { return right_;}

private:
  Call<V, I, J> left_;
  S right_;
};

// Binary expression: the left operand is a scalar, the right a call.
// This is valid only for O == Mult.
template <typename S, typename V, typename I, typename J>
class Binary_expr<S, Call<V, I, J>, vsip::impl::op::Mult>
{
public:
  typedef V view_type;

  Binary_expr(S s, Call<V, I, J> c) : left_(s), right_(c) {}
  view_type view() { return right_.view();}
  S left() { return left_;}
  Call<V, I, J> right() { return right_;}

private:
  S left_;
  Call<V, I, J> right_;
};

template <typename L, typename R>
typename vsip::impl::Type_if<Binary_expr<L, R, vsip::impl::op::Add>, 
                             Is_expr<L>::value || Is_expr<R>::value>::type
operator+ (L l, R r)
{
  return Binary_expr<L, R, vsip::impl::op::Add>(l, r);
}

template <typename L, typename R>
typename vsip::impl::Type_if<Binary_expr<L, R, vsip::impl::op::Sub>,
                             Is_expr<L>::value || Is_expr<R>::value>::type
operator- (L l, R r)
{
  return Binary_expr<L, R, vsip::impl::op::Sub>(l, r);
}

template <typename L, typename R>
typename vsip::impl::Type_if<Binary_expr<L, R, vsip::impl::op::Mult>,
                             Is_expr<L>::value || Is_expr<R>::value>::type
operator* (L l, R r)
{
  return Binary_expr<L, R, vsip::impl::op::Mult>(l, r);
}

template <typename L, typename R>
typename vsip::impl::Type_if<Binary_expr<L, R, vsip::impl::op::Div>,
                             Is_expr<L>::value || Is_expr<R>::value>::type
operator/ (L l, R r)
{
  return Binary_expr<L, R, vsip::impl::op::Div>(l, r);
}

} // namespace vsip_csl::stencil
} // namespace vsip_csl

#endif
