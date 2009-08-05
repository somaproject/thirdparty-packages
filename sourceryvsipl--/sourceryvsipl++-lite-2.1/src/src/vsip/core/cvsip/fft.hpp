/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/cvsip/fft.hpp
    @author  Stefan Seefeld
    @date    2006-10-16
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             C-VSIPL.
*/

#ifndef VSIP_CORE_CVSIP_FFT_HPP
#define VSIP_CORE_CVSIP_FFT_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/support.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/fft/factory.hpp>
#include <vsip/core/fft/util.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cvsip
{

template <typename I, dimension_type D, typename S>
std::auto_ptr<I>
create(Domain<D> const &dom, S scale, unsigned int);

template <>
std::auto_ptr<fft::backend<1, float, std::complex<float>, 0, -1> >
create(Domain<1> const &, float, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, float, 0, 1> >
create(Domain<1> const &, float, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, std::complex<float>, 0, -1> >
create(Domain<1> const &, float, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, std::complex<float>, 0, 1> >
create(Domain<1> const &, float, unsigned int);
template <>
std::auto_ptr<fft::backend<1, double, std::complex<double>, 0, -1> >
create(Domain<1> const &, double, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, double, 0, 1> >
create(Domain<1> const &, double, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, std::complex<double>, 0, -1> >
create(Domain<1> const &, double, unsigned int);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, std::complex<double>, 0, 1> >
create(Domain<1> const &, double, unsigned int);
} // namespace vsip::impl::cvsip


namespace fft
{
struct Cvsip_tag;

template <typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<1, I, O, S, R, N, Cvsip_tag>
{
  static bool const has_float =
#if VSIP_IMPL_CVSIP_HAVE_FLOAT
    true
#else
    false
#endif
    ;
  static bool const has_double =
#if VSIP_IMPL_CVSIP_HAVE_DOUBLE
    true
#else
    false
#endif
    ;
  static bool const ct_valid = (has_float && 
                                Type_equal<typename Scalar_of<I>::type,
                                           float>::value) ||
                               (has_double && 
                                Type_equal<typename Scalar_of<I>::type,
                                           double>::value);
  static bool rt_valid(Domain<1> const &/*dom*/) { return true;}
  static std::auto_ptr<backend<1, I, O,
			       axis<I, O, S>::value,
			       exponent<I, O, S>::value> >
  create(Domain<1> const &dom, typename Scalar_of<I>::type scale)
  {
    return cvsip::create<backend<1, I, O,
                                 axis<I, O, S>::value,
                                 exponent<I, O, S>::value> >
      (dom, scale, N);
  }
};

} // namespace vsip::impl::fft

namespace fftm
{
template <typename I,
	  typename O,
	  int A,
	  int E,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<I, O, A, E, R, N, fft::Cvsip_tag>
{
  static bool const ct_valid = !Type_equal<typename Scalar_of<I>::type,
                                           long double>::value;
  static bool rt_valid(Domain<2> const& /*dom*/) { return true;}
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename Scalar_of<I>::type scale)
  {
    return cvsip::create<fft::fftm<I, O, A, E> >(dom, scale, N);
  }
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif
