/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ipp/fft.hpp
    @author  Stefan Seefeld
    @date    2006-05-05
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             Intel's IPP.
*/

#ifndef VSIP_IMPL_IPP_FFT_HPP
#define VSIP_IMPL_IPP_FFT_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

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
namespace ipp
{

/// These are the entry points into the IPP FFT bridge.
template <typename I, dimension_type D, typename S>
std::auto_ptr<I>
create(Domain<D> const &dom, S scale, bool fast);

template <dimension_type D, typename I, typename O, int S>
struct evaluator
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<D> const &) { return true;}
  static std::auto_ptr<fft::backend<D, I, O,
				    fft::axis<I, O, S>::value,
				    fft::exponent<I, O, S>::value> >
  create(Domain<D> const &dom, typename Scalar_of<I>::type scale)
  {
    bool fast = fft::is_power_of_two(dom);
    return ipp::create<fft::backend<D, I, O,
                       fft::axis<I, O, S>::value,
                       fft::exponent<I, O, S>::value> >
      (dom, scale, fast);
  }
};

} // namespace vsip::impl::ipp

namespace fft
{
struct Intel_ipp_tag;

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<float>, std::complex<float>, S, R, N,
		 Intel_ipp_tag>
  : ipp::evaluator<1, std::complex<float>, std::complex<float>, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, float, std::complex<float>, S, R, N, Intel_ipp_tag>
  : ipp::evaluator<1, float, std::complex<float>, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<float>, float, S, R, N, Intel_ipp_tag>
  : ipp::evaluator<1, std::complex<float>, float, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<double>, std::complex<double>, S, R, N,
		 Intel_ipp_tag>
  : ipp::evaluator<1, std::complex<double>, std::complex<double>, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, double, std::complex<double>, S, R, N, Intel_ipp_tag>
  : ipp::evaluator<1, double, std::complex<double>, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<double>, double, S, R, N, Intel_ipp_tag>
  : ipp::evaluator<1, std::complex<double>, double, S> {};

// Complex 2D FFTs

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<2, std::complex<float>, std::complex<float>, S, R, N,
		 Intel_ipp_tag>
  : ipp::evaluator<2, std::complex<float>, std::complex<float>, S> {};

} // namespace vsip::impl::fft

namespace fftm
{
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<I, O, A, E, R, N, fft::Intel_ipp_tag>
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<2> const &) { return true;}
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename Scalar_of<I>::type scale)
  {
    bool fast = fft::is_power_of_two(dom);
    return ipp::create<fft::fftm<I, O, A, E> >(dom, scale, fast);
  }
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<std::complex<long double>, std::complex<long double>,
		 A, E, R, N, fft::Intel_ipp_tag>
{
  static bool const ct_valid = false;
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<long double, std::complex<long double>,
		 A, E, R, N, fft::Intel_ipp_tag>
{
  static bool const ct_valid = false;
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<std::complex<long double>, long double,
		 A, E, R, N, fft::Intel_ipp_tag>
{
  static bool const ct_valid = false;
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif

