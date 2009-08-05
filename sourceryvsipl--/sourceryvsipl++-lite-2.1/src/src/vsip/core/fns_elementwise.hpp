/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fns_elementwise.hpp
    @author  Stefan Seefeld
    @date    2005-04-21
    @brief   VSIPL++ Library: [math.fns.elementwise].

    This file declares functions to be used with expression templates.
*/

#ifndef VSIP_CORE_FNS_ELEMENTWISE_HPP
#define VSIP_CORE_FNS_ELEMENTWISE_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/expr/functor.hpp>
#include <vsip/core/promote.hpp>
#include <vsip/core/fns_scalar.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/lvalue_proxy.hpp>

namespace vsip
{
namespace impl
{

/// Macro to define a unary function on views in terms of
/// its homologe on scalars.
#define VSIP_IMPL_UNARY_FUNCTOR(fname)                                    \
template <typename T>                                                     \
struct fname##_functor                                                    \
{                                                                         \
  typedef T result_type;                                                  \
  static char const* name() { return #fname; }                            \
  static result_type apply(T t) { return fn::fname(t);}                   \
  result_type operator()(T t) const { return apply(t);}                   \
};

#define VSIP_IMPL_UNARY_FUNCTOR_RETN(fname, retn)                         \
template <typename T>                                                     \
struct fname##_functor                                                    \
{                                                                         \
  typedef retn result_type;                                               \
  static char const* name() { return #fname; }                            \
  static result_type apply(T t) { return fn::fname(t);}                   \
  result_type operator()(T t) const { return apply(t);}                   \
};

#define VSIP_IMPL_UNARY_DISPATCH(fname)                                   \
template <typename T>                                                     \
struct Dispatch_##fname :                                                 \
  ITE_Type<Is_view_type<T>::value,                                        \
           As_type<Unary_func_view<fname##_functor, T> >,                 \
  ITE_Type<Is_lvalue_proxy_type<T>::value,			          \
           As_type<fname##_functor<typename Is_lvalue_proxy_type<T>::value_type> >,	\
           As_type<fname##_functor<T> > > >::type                         \
{                                                                         \
};

#define VSIP_IMPL_UNARY_FUNCTION(fname)                                   \
template <typename T>                                                     \
inline                                                                    \
typename Dispatch_##fname<T>::result_type                                 \
fname(T const& t) { return Dispatch_##fname<T>::apply(t);}

/// This function gateway is roughly specialized to VSIPL++ view types.
/// This prevents it from competing with more general function overloads.
/// For example, cmath defines template <typename T> bool isnan(T& t)
/// which is ambiguous with the normal VSIP_IMPL_UNARY_FUNCTION.
#define VSIP_IMPL_UNARY_VIEW_FUNCTION(fname)				  \
template <template <typename, typename> class V,                          \
          typename T, typename B>                                         \
inline                                                                    \
typename Dispatch_##fname<V<T,B> >::result_type				  \
fname(V<T,B> const& t)							  \
{ return Dispatch_##fname<V<T,B> >::apply(t);}

#define VSIP_IMPL_UNARY_FUNC(fname)                                       \
VSIP_IMPL_UNARY_FUNCTOR(fname)                                            \
VSIP_IMPL_UNARY_DISPATCH(fname)                                           \
VSIP_IMPL_UNARY_FUNCTION(fname)

#define VSIP_IMPL_UNARY_FUNC_RETN(fname, retn)                            \
VSIP_IMPL_UNARY_FUNCTOR_RETN(fname, retn)                                 \
VSIP_IMPL_UNARY_DISPATCH(fname)                                           \
VSIP_IMPL_UNARY_FUNCTION(fname)

#define VSIP_IMPL_UNARY_VIEW_FUNC_RETN(fname, retn)                       \
VSIP_IMPL_UNARY_FUNCTOR_RETN(fname, retn)                                 \
VSIP_IMPL_UNARY_DISPATCH(fname)                                           \
VSIP_IMPL_UNARY_VIEW_FUNCTION(fname)

/// Define a unary operator. Assume the associated Dispatch 
/// is already defined.
#define VSIP_IMPL_UNARY_OP(op, fname)                                     \
template <typename T>                                                     \
typename Dispatch_##fname<typename Is_view_type<T>::type>::result_type    \
operator op(T t)                                                          \
{ return Dispatch_##fname<T>::apply(t);}

/// Macro to define a binary function on views in terms of
/// its homologe on scalars.
#define VSIP_IMPL_BINARY_FUNCTOR(fname)                                   \
template <typename T1, typename T2>                                       \
struct fname##_functor                                                    \
{                                                                         \
  typedef typename Promotion<T1, T2>::type result_type;                   \
  static char const* name() { return #fname; }                            \
  static result_type apply(T1 t1, T2 t2) { return fn::fname(t1, t2);}     \
  result_type operator()(T1 t1, T2 t2) const { return apply(t1, t2);}     \
};

#define VSIP_IMPL_BINARY_FUNCTOR_RETN(fname, retn)                        \
template <typename T1, typename T2>                                       \
struct fname##_functor                                                    \
{                                                                         \
  typedef retn result_type;                                               \
  static char const* name() { return #fname; }                            \
  static result_type apply(T1 t1, T2 t2) { return fn::fname(t1, t2);}     \
  result_type operator()(T1 t1, T2 t2) const { return apply(t1, t2);}     \
};

#define VSIP_IMPL_BINARY_FUNCTOR_SCALAR_RETN(fname)                       \
template <typename T1, typename T2>                                       \
struct fname##_functor                                                    \
{                                                                         \
  typedef typename Scalar_of<typename Promotion<T1, T2>::type>::type      \
                result_type;                                              \
  static char const* name() { return #fname; }                            \
  static result_type apply(T1 t1, T2 t2) { return fn::fname(t1, t2);}     \
  result_type operator()(T1 t1, T2 t2) const { return apply(t1, t2);}     \
};

#define VSIP_IMPL_BINARY_DISPATCH(fname)                                  \
template <typename T1, typename T2>                                       \
struct Dispatch_##fname :                                                 \
  ITE_Type<Is_view_type<T1>::value || Is_view_type<T2>::value,            \
           As_type<Binary_func_view<fname##_functor, T1, T2> >,           \
           As_type<fname##_functor<T1, T2> > >::type                      \
{                                                                         \
};

#define VSIP_IMPL_BINARY_DISPATCH_USEOP(fname, opname)			  \
template <typename T1, typename T2>                                       \
struct Dispatch_##fname :                                                 \
  ITE_Type<Is_view_type<T1>::value || Is_view_type<T2>::value,            \
           As_type<Binary_func_view<opname, T1, T2> >,			  \
           As_type<opname<T1, T2> > >::type				  \
{                                                                         \
};

/// Define a dispatcher that only matches if at least one of the arguments
/// is a view type.
#define VSIP_IMPL_BINARY_OP_DISPATCH(fname)                               \
template <typename T1, typename T2,                                       \
          bool P = Is_view_type<T1>::value || Is_view_type<T2>::value>    \
struct Dispatch_op_##fname                                                \
  : As_type<Binary_func_view<fname##_functor, T1, T2> >::type {};         \
template <typename T1, typename T2>                                       \
struct Dispatch_op_##fname<T1, T2, false> {};                             \

#define VSIP_IMPL_BINARY_FUNCTION(fname)                                  \
template <typename T1, typename T2>                                       \
inline                                                                    \
typename Dispatch_##fname<T1, T2>::result_type                            \
fname(T1 const& t1, T2 const& t2)					  \
{ return Dispatch_##fname<T1, T2>::apply(t1, t2); }

#define VSIP_IMPL_BINARY_OPERATOR_ONE(op, fname)                          \
template <typename T1, typename T2>                                       \
inline                                                                    \
typename Dispatch_op_##fname<T1, T2>::result_type                         \
operator op(T1 const& t1, T2 const& t2)					  \
{ return Dispatch_op_##fname<T1, T2>::apply(t1, t2);}

#define VSIP_IMPL_BINARY_OPERATOR_TWO(op, fname)                          \
template <template <typename, typename> class View,                       \
          typename T1, typename Block1, typename T2>                      \
inline                                                                    \
typename Dispatch_op_##fname<View<T1, Block1>, T2>::result_type           \
operator op(View<T1, Block1> const& t1, T2 t2)				  \
{ return Dispatch_op_##fname<View<T1, Block1>, T2>::apply(t1, t2);}       \
                                                                          \
template <template <typename, typename> class View,                       \
          typename T1, typename T2, typename Block2>                      \
inline                                                                    \
typename Dispatch_op_##fname<T1, View<T2, Block2> >::result_type          \
operator op(T1 t1, View<T2, Block2> const& t2)				  \
{ return Dispatch_op_##fname<T1, View<T2, Block2> >::apply(t1, t2);}      \
                                                                          \
template <template <typename, typename> class LView,                      \
          template <typename, typename> class RView,                      \
          typename T1, typename Block1,                                   \
          typename T2, typename Block2>                                   \
inline                                                                    \
typename Dispatch_op_##fname<LView<T1, Block1>,                           \
                             RView<T2, Block2> >::result_type             \
operator op(LView<T1, Block1> const& t1, RView<T2, Block2> const& t2)	  \
{ return Dispatch_op_##fname<LView<T1, Block1>,                           \
                             RView<T2, Block2> >::apply(t1, t2);}

#if (defined(__GNUC__) && __GNUC__ < 4) || defined(__ghs__) || defined(__ICL)
# define VSIP_IMPL_BINARY_OPERATOR(op, fname)                             \
VSIP_IMPL_BINARY_OPERATOR_ONE(op, fname)
#else
# define VSIP_IMPL_BINARY_OPERATOR(op, fname)                             \
VSIP_IMPL_BINARY_OPERATOR_TWO(op, fname)
#endif


#define VSIP_IMPL_BINARY_VIEW_FUNCTION(fname)                             \
template <template <typename, typename> class V,                          \
          typename T, typename B>                                         \
inline                                                                    \
typename Dispatch_##fname<V<T,B>, V<T,B> >::result_type                   \
fname(V<T,B> const& t1, V<T,B> const& t2)				  \
{ return Dispatch_##fname<V<T,B>, V<T,B> >::apply(t1, t2);}

#define VSIP_IMPL_BINARY_FUNC(fname)                                      \
VSIP_IMPL_BINARY_FUNCTOR(fname)                                           \
VSIP_IMPL_BINARY_DISPATCH(fname)                                          \
VSIP_IMPL_BINARY_FUNCTION(fname)                                          \
VSIP_IMPL_BINARY_VIEW_FUNCTION(fname)

/// Binary function that can use an existing op instead of a functor
/// For example, arithmetic operations like mul() can us op::Mult.
#define VSIP_IMPL_BINARY_FUNC_USEOP(fname, opname)			  \
VSIP_IMPL_BINARY_DISPATCH_USEOP(fname, opname) 				  \
VSIP_IMPL_BINARY_FUNCTION(fname)                                          \
VSIP_IMPL_BINARY_VIEW_FUNCTION(fname)

#define VSIP_IMPL_BINARY_FUNC_RETN(fname, retn)                           \
VSIP_IMPL_BINARY_FUNCTOR_RETN(fname, retn)                                \
VSIP_IMPL_BINARY_DISPATCH(fname)                                          \
VSIP_IMPL_BINARY_FUNCTION(fname)

#define VSIP_IMPL_BINARY_FUNC_SCALAR_RETN(fname)                          \
VSIP_IMPL_BINARY_FUNCTOR_SCALAR_RETN(fname)                               \
VSIP_IMPL_BINARY_DISPATCH(fname)                                          \
VSIP_IMPL_BINARY_FUNCTION(fname)

#define VSIP_IMPL_BINARY_OP(op, fname)                                    \
VSIP_IMPL_BINARY_OP_DISPATCH(fname)                                       \
VSIP_IMPL_BINARY_OPERATOR(op, fname)

/// Macro to define a ternary function on views in terms of
/// its homologe on scalars.
#define VSIP_IMPL_TERNARY_FUNC(fname)                                     \
template <typename T1, typename T2, typename T3>                          \
struct fname##_functor                                                    \
{                                                                         \
  typedef typename Promotion<typename Promotion<T1, T2>::type,            \
                             T3>::type result_type;                       \
  static char const* name() { return #fname; }                            \
  static result_type apply(T1 t1, T2 t2, T3 t3)                           \
  { return fn::fname(t1, t2, t3);}                                        \
  result_type operator()(T1 t1, T2 t2, T3 t3) const                       \
  { return apply(t1, t2, t3);}                                            \
};                                                                        \
                                                                          \
template <typename T1, typename T2, typename T3>                          \
struct Dispatch_##fname :                                                 \
  ITE_Type<Is_view_type<T1>::value ||                                     \
           Is_view_type<T2>::value ||                                     \
           Is_view_type<T3>::value,                                       \
           As_type<Ternary_func_view<fname##_functor, T1, T2, T3> >,      \
           As_type<fname##_functor<T1, T2, T3> > >::type                  \
{                                                                         \
};                                                                        \
                                                                          \
template <typename T1, typename T2, typename T3>                          \
inline									  \
typename Dispatch_##fname<T1, T2, T3>::result_type                        \
fname(T1 const& t1, T2 const& t2, T3 const& t3)				  \
{ return Dispatch_##fname<T1, T2, T3>::apply(t1, t2, t3);}


/***********************************************************************
  Unary Functions
***********************************************************************/

VSIP_IMPL_UNARY_FUNC(acos)

template <typename T> struct arg_functor {};

template <typename T>
struct arg_functor<std::complex<T> >
{
  typedef T result_type;
  static char const* name() { return "arg"; }                
  static result_type apply(std::complex<T> t) { return fn::arg(t);}
  result_type operator()(std::complex<T> t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(arg)
VSIP_IMPL_UNARY_FUNCTION(arg)

VSIP_IMPL_UNARY_FUNC(asin)
VSIP_IMPL_UNARY_FUNC(atan)
VSIP_IMPL_UNARY_FUNC(bnot)
VSIP_IMPL_UNARY_FUNC(ceil)
VSIP_IMPL_UNARY_FUNC(conj)
VSIP_IMPL_UNARY_FUNC(cos)
VSIP_IMPL_UNARY_FUNC(cosh)

template <typename T>
struct euler_functor
{
  typedef std::complex<T> result_type;
  static char const* name() { return "euler"; }                
  static result_type apply(T t) { return fn::euler(t);}
  result_type operator()(T t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(euler)
VSIP_IMPL_UNARY_FUNCTION(euler)

VSIP_IMPL_UNARY_FUNC(exp)
VSIP_IMPL_UNARY_FUNC(exp10)
VSIP_IMPL_UNARY_FUNC(floor)

template <typename T> struct imag_functor {};

template <typename T>
struct imag_functor<std::complex<T> >
{
  typedef T result_type;
  static char const* name() { return "imag"; }                
  static result_type apply(std::complex<T> t) { return fn::imag(t);}
  result_type operator()(std::complex<T> t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(imag)
VSIP_IMPL_UNARY_FUNCTION(imag)

VSIP_IMPL_UNARY_FUNC_RETN(is_finite, bool)
VSIP_IMPL_UNARY_FUNC_RETN(is_nan, bool)
VSIP_IMPL_UNARY_FUNC_RETN(is_normal, bool)

VSIP_IMPL_UNARY_FUNC_RETN(lnot, bool)
VSIP_IMPL_UNARY_FUNC(log)
VSIP_IMPL_UNARY_FUNC(log10)
VSIP_IMPL_UNARY_FUNC_RETN(mag, typename impl::Scalar_of<T>::type)
VSIP_IMPL_UNARY_FUNC_RETN(magsq, typename impl::Scalar_of<T>::type)
VSIP_IMPL_UNARY_FUNC(neg)

template <typename T> struct real_functor {};

template <typename T>
struct real_functor<std::complex<T> >
{
  typedef T result_type;
  static char const* name() { return "real"; }                
  static result_type apply(std::complex<T> t) { return fn::real(t);}
  result_type operator()(std::complex<T> t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(real)
VSIP_IMPL_UNARY_FUNCTION(real)

VSIP_IMPL_UNARY_FUNC(recip)
VSIP_IMPL_UNARY_FUNC(rsqrt)
VSIP_IMPL_UNARY_FUNC(sin)
VSIP_IMPL_UNARY_FUNC(sinh)
VSIP_IMPL_UNARY_FUNC(sq)
VSIP_IMPL_UNARY_FUNC(sqrt)
VSIP_IMPL_UNARY_FUNC(tan)
VSIP_IMPL_UNARY_FUNC(tanh)

VSIP_IMPL_UNARY_FUNC(impl_conj)

template <typename T>
struct impl_real_functor
{
  typedef typename Scalar_of<T>::type result_type;
  static char const* name() { return "impl_real"; }                
  static result_type apply(T t) { return fn::impl_real(t);}
  result_type operator()(T t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(impl_real)
VSIP_IMPL_UNARY_FUNCTION(impl_real)

template <typename T>
struct impl_imag_functor
{
  typedef typename Scalar_of<T>::type result_type;
  static char const* name() { return "impl_imag"; }                
  static result_type apply(T t) { return fn::impl_imag(t);}
  result_type operator()(T t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(impl_imag)
VSIP_IMPL_UNARY_FUNCTION(impl_imag)

// This unary operator gives a hint to the compiler that this block is
// unaligned
template <typename T>
struct unaligned_functor
{
  typedef T result_type;
  static char const* name() { return "unaligned"; }                
  static result_type apply(T t) { return t;}
  result_type operator()(T t) const { return apply(t);}
};

VSIP_IMPL_UNARY_DISPATCH(unaligned)
VSIP_IMPL_UNARY_FUNCTION(unaligned)

/***********************************************************************
  Binary Functions
***********************************************************************/

VSIP_IMPL_BINARY_FUNC(add)
VSIP_IMPL_BINARY_FUNC(atan2)
VSIP_IMPL_BINARY_FUNC(band)
VSIP_IMPL_BINARY_FUNC(bor)
VSIP_IMPL_BINARY_FUNC(bxor)
VSIP_IMPL_BINARY_FUNC(div)
VSIP_IMPL_BINARY_FUNC_RETN(eq, bool)
VSIP_IMPL_BINARY_FUNC(fmod)
VSIP_IMPL_BINARY_FUNC_RETN(ge, bool)
VSIP_IMPL_BINARY_FUNC_RETN(gt, bool)
VSIP_IMPL_BINARY_FUNC(hypot)
VSIP_IMPL_BINARY_FUNC(jmul)
VSIP_IMPL_BINARY_FUNC_RETN(land, bool)
VSIP_IMPL_BINARY_FUNC_RETN(le, bool)
VSIP_IMPL_BINARY_FUNC_RETN(lt, bool)
VSIP_IMPL_BINARY_FUNC_RETN(lor, bool)
VSIP_IMPL_BINARY_FUNC_RETN(lxor, bool)
VSIP_IMPL_BINARY_FUNC_USEOP(mul, op::Mult)
VSIP_IMPL_BINARY_FUNC(max)
VSIP_IMPL_BINARY_FUNC(maxmg)
VSIP_IMPL_BINARY_FUNC_SCALAR_RETN(maxmgsq)
VSIP_IMPL_BINARY_FUNC(min)
VSIP_IMPL_BINARY_FUNC(minmg)
VSIP_IMPL_BINARY_FUNC_SCALAR_RETN(minmgsq)
VSIP_IMPL_BINARY_FUNC_RETN(ne, bool)
VSIP_IMPL_BINARY_FUNC(pow)
VSIP_IMPL_BINARY_FUNC(sub)

/***********************************************************************
  Ternary Functions
***********************************************************************/

VSIP_IMPL_TERNARY_FUNC(am)
VSIP_IMPL_TERNARY_FUNC(expoavg)
VSIP_IMPL_TERNARY_FUNC(ma)
VSIP_IMPL_TERNARY_FUNC(msb)
VSIP_IMPL_TERNARY_FUNC(sbm)
VSIP_IMPL_TERNARY_FUNC(ite)

/***********************************************************************
  Unary Operators
***********************************************************************/

VSIP_IMPL_UNARY_OP(!, lnot)
VSIP_IMPL_UNARY_OP(~, bnot)

/***********************************************************************
  Binary Operators
***********************************************************************/

VSIP_IMPL_BINARY_OP(==, eq)
VSIP_IMPL_BINARY_OP(>=, ge)
VSIP_IMPL_BINARY_OP(>, gt)
VSIP_IMPL_BINARY_OP(<=, le)
VSIP_IMPL_BINARY_OP(<, lt)
VSIP_IMPL_BINARY_OP(!=, ne)
VSIP_IMPL_BINARY_OP(&&, land)
VSIP_IMPL_BINARY_OP(&, band)
VSIP_IMPL_BINARY_OP(||, lor)
VSIP_IMPL_BINARY_OP(|, bor)

template <typename T1, typename T2>
struct bxor_or_lxor_functor
{
  typedef typename Promotion<T1, T2>::type result_type;
  static char const* name() { return "bxor"; }                
  static result_type apply(T1 t1, T2 t2) { return fn::bxor(t1, t2);}
  result_type operator()(T1 t1, T2 t2) const { return apply(t1, t2);}
};

template <>
struct bxor_or_lxor_functor<bool, bool>
{
  typedef bool result_type;
  static char const* name() { return "lxor"; }                
  static result_type apply(bool t1, bool t2) { return fn::lxor(t1, t2);}
  result_type operator()(bool t1, bool t2) const { return apply(t1, t2);}
};

VSIP_IMPL_BINARY_OP(^, bxor_or_lxor)

} // namespace vsip::impl

using impl::acos;
using impl::arg;
using impl::asin;
using impl::atan;
using impl::bnot;
using impl::ceil;
using impl::conj;
using impl::cos;
using impl::cosh;
using impl::euler;
using impl::exp;
using impl::exp10;
using impl::floor;
using impl::imag;
using impl::lnot;
using impl::log;
using impl::log10;
using impl::mag;
using impl::magsq;
using impl::neg;
using impl::real;
using impl::recip;
using impl::rsqrt;
using impl::sin;
using impl::sinh;
using impl::sq;
using impl::sqrt;
using impl::tan;
using impl::tanh;

using impl::add;
using impl::atan2;
using impl::band;
using impl::bor;
using impl::bxor;
using impl::div;
using impl::eq;
using impl::fmod;
using impl::ge;
using impl::gt;
using impl::hypot;
using impl::jmul;
using impl::land;
using impl::le;
using impl::lt;
using impl::lor;
using impl::lxor;
using impl::max;
using impl::maxmg;
using impl::maxmgsq;
using impl::min;
using impl::minmg;
using impl::minmgsq;
using impl::mul;
using impl::ne;
using impl::pow;
using impl::sub;

using impl::am;
using impl::expoavg;
using impl::ma;
using impl::msb;
using impl::sbm;

using impl::operator!;
using impl::operator~;
using impl::operator==;
using impl::operator>=;
using impl::operator>;
using impl::operator<=;
using impl::operator<;
using impl::operator!=;
using impl::operator&&;
using impl::operator&;
using impl::operator||;
using impl::operator|;

} // namespace vsip

#undef VSIP_IMPL_TERNARY_FUNC
#undef VSIP_IMPL_BINARY_FUNC_RETN
#undef VSIP_IMPL_BINARY_FUNC
#undef VSIP_IMPL_UNARY_FUNC_RETN
#undef VSIP_IMPL_UNARY_FUNC
#undef VSIP_IMPL_BINARY_FUNCTOR_RETN
#undef VSIP_IMPL_BINARY_FUNCTOR
#undef VSIP_IMPL_BINARY_DISPATCH
#undef VSIP_IMPL_BINARY_FUNCTION
#undef VSIP_IMPL_UNARY_FUNCTOR_RETN
#undef VSIP_IMPL_UNARY_FUNCTOR
#undef VSIP_IMPL_UNARY_DISPATCH
#undef VSIP_IMPL_UNARY_FUNCTION

#endif
