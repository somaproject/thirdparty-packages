/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fft/factory.hpp
    @author  Stefan Seefeld
    @date    2006-02-20
    @brief   VSIPL++ Library: Fft backend dispatch harness.
*/

#ifndef VSIP_CORE_FFT_FACTORY_HPP
#define VSIP_CORE_FFT_FACTORY_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#define VSIP_IMPL_VERBOSE_FFT_EXCEPTION 1

#if VSIP_IMPL_VERBOSE_FFT_EXCEPTION
#  include <sstream>
#endif

#include <vsip/core/fft/backend.hpp>
#include <vsip/core/fft/util.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/type_list.hpp>
#include <memory>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace fft
{

/// evaluator template.
/// This needs to be provided for each tag in the LibraryTagList.
template <dimension_type D,
	  typename I,
	  typename O,
	  int sD,
	  vsip::return_mechanism_type rmT,
	  unsigned nT,
	  typename Tag>
struct evaluator
{
  static bool const ct_valid = false;
};

template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag = typename TagList::first,
	  typename Rest = typename TagList::rest,
	  typename Eval = evaluator<D, I, O, S, R, N, Tag>,
	  bool CtValid = Eval::ct_valid>
struct factory;

/// In case the compile-time check passes, we decide at run-time whether
/// or not to use this backend.
template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename Eval>
struct factory<D, I, O, S, R, N, TagList, Tag, Rest, Eval, true>
{
  typedef backend<D, I, O,
		  axis<I, O, S>::value,
		  exponent<I, O, S>::value> interface;
  static std::auto_ptr<interface> 
  create(vsip::Domain<D> const& dom, typename interface::scalar_type scale)
  {
    if (Eval::rt_valid(dom)) return Eval::create(dom, scale);
    else return factory<D, I, O, S, R, N, Rest>::create(dom, scale);
  }
};

/// In case the compile-time check fails, we continue the search
/// directly at the next entry in the type list.
template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename Eval>
struct factory<D, I, O, S, R, N, TagList, Tag, Rest, Eval, false>
  : factory<D, I, O, S, R, N, Rest>
{};

/// Terminator. Instead of passing on to the next element
/// it aborts the program. It is a program error to define
/// callback lists that can't handle a given expression.
template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Eval>
struct factory<D, I, O, S, R, N, TagList, Tag, None_type, Eval, true>
{
  typedef backend<D, I, O,
		  axis<I, O, S>::value,
		  exponent<I, O, S>::value> interface;
  static std::auto_ptr<interface> 
  create(Domain<D> const &dom, typename interface::scalar_type scale)
  {
    if (Eval::rt_valid(dom)) return Eval::create(dom, scale);
#if VSIP_IMPL_VERBOSE_FFT_EXCEPTION
    std::ostringstream msg;
    msg << "Requested Fft "
	<< "(dim: " << D << "  size: " << dom[0].size();
    for (index_type i=1; i<D; ++i)
      msg << "," << dom[i].size();
    msg << ") not supported by available backends";
    VSIP_IMPL_THROW(std::runtime_error(msg.str()));
#else
    VSIP_IMPL_THROW(std::runtime_error(
		      "Requested Fft not supported by available backends"));
#endif
    return std::auto_ptr<interface>();
  }
};

} // namespace vsip::impl::fft

namespace fftm
{

/// evaluator template.
/// This needs to be provided for each tag in the LibraryTagList.
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N,
	  typename Tag>
struct evaluator
{
  static bool const ct_valid = false;
};

template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag = typename TagList::first,
	  typename Rest = typename TagList::rest,
	  typename Eval = evaluator<I, O, A, E, R, N, Tag>,
	  bool CtValid = Eval::ct_valid>
struct factory;

/// In case the compile-time check passes, we decide at run-time whether
/// or not to use this backend.
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename Eval>
struct factory<I, O, A, E, R, N, TagList, Tag, Rest, Eval, true>
{
  typedef fft::fftm<I, O, A, E> interface;
  static std::auto_ptr<interface> 
  create(vsip::Domain<2> const& dom, typename interface::scalar_type scale)
  {
    if (Eval::rt_valid(dom)) return Eval::create(dom, scale);
    else return factory<I, O, A, E, R, N, Rest>::create(dom, scale);
  }
};

/// In case the compile-time check fails, we continue the search
/// directly at the next entry in the type list.
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Rest,
	  typename Eval>
struct factory<I, O, A, E, R, N, TagList, Tag, Rest, Eval, false>
  : factory<I, O, A, E, R, N, Rest>
{};

/// Terminator. Instead of passing on to the next element
/// it aborts the program. It is a program error to define
/// callback lists that can't handle a given expression.
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N,
	  typename TagList,
	  typename Tag,
	  typename Eval>
struct factory<I, O, A, E, R, N, TagList, Tag, None_type, Eval, true>
{
  typedef fft::fftm<I, O, A, E> interface;
  static std::auto_ptr<interface> 
  create(Domain<2> const &dom, typename interface::scalar_type scale)
  {
    if (Eval::rt_valid(dom)) return Eval::create(dom, scale);
#if VSIP_IMPL_VERBOSE_FFT_EXCEPTION
    std::ostringstream msg;
    msg << "Requested Fftm "
	<< "(size: " << dom[0].size() << "," << dom[1].size() 
	<< ") not supported by available backends";
    VSIP_IMPL_THROW(std::runtime_error(msg.str()));
#else
    VSIP_IMPL_THROW(std::runtime_error(
		      "Requested Fftm not supported by available backends"));
#endif
    return std::auto_ptr<interface>();
  }
};

} // namespace vsip::impl::fftm

} // namespace vsip::impl
} // namespace vsip

#endif
