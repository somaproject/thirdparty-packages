/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/fft.hpp
    @author  Stefan Seefeld
    @date    2006-02-02
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             Mercury's SAL.
*/

#ifndef VSIP_IMPL_SAL_FFT_HPP
#define VSIP_IMPL_SAL_FFT_HPP

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
namespace sal
{

template <typename I, dimension_type D, typename S>
std::auto_ptr<I>
create(Domain<D> const &dom, S scale);

#if VSIP_IMPL_HAVE_SAL_FLOAT
template <>
std::auto_ptr<fft::backend<1, float, std::complex<float>, 0, -1> >
create(Domain<1> const &, float);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, float, 0, 1> >
create(Domain<1> const &, float);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, std::complex<float>, 0, -1> >
create(Domain<1> const &, float);
template <>
std::auto_ptr<fft::backend<1, std::complex<float>, std::complex<float>, 0, 1> >
create(Domain<1> const &, float);
#endif

#if VSIP_IMPL_HAVE_SAL_DOUBLE
template <>
std::auto_ptr<fft::backend<1, double, std::complex<double>, 0, -1> >
create(Domain<1> const &, double);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, double, 0, 1> >
create(Domain<1> const &, double);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, std::complex<double>, 0, -1> >
create(Domain<1> const &, double);
template <>
std::auto_ptr<fft::backend<1, std::complex<double>, std::complex<double>, 0, 1> >
create(Domain<1> const &, double);
#endif

#if VSIP_IMPL_HAVE_SAL_FLOAT
template <>
std::auto_ptr<fft::backend<2, float, std::complex<float>, 0, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, float, std::complex<float>, 1, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, float, 0, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, float, 1, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, std::complex<float>, 0, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, std::complex<float>, 1, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, std::complex<float>, 0, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::backend<2, std::complex<float>, std::complex<float>, 1, 1> >
create(Domain<2> const &, float);
#endif

#if VSIP_IMPL_HAVE_SAL_DOUBLE
template <>
std::auto_ptr<fft::backend<2, double, std::complex<double>, 0, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, double, std::complex<double>, 1, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, double, 0, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, double, 1, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, std::complex<double>, 0, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, std::complex<double>, 1, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, std::complex<double>, 0, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::backend<2, std::complex<double>, std::complex<double>, 1, 1> >
create(Domain<2> const &, double);
#endif

#if VSIP_IMPL_HAVE_SAL_FLOAT
template <>
std::auto_ptr<fft::fftm<float, std::complex<float>, 0, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<float, std::complex<float>, 1, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, float, 0, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, float, 1, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, std::complex<float>, 0, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, std::complex<float>, 1, -1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, std::complex<float>, 0, 1> >
create(Domain<2> const &, float);
template <>
std::auto_ptr<fft::fftm<std::complex<float>, std::complex<float>, 1, 1> >
create(Domain<2> const &, float);
#endif

#if VSIP_IMPL_HAVE_SAL_DOUBLE
template <>
std::auto_ptr<fft::fftm<double, std::complex<double>, 0, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<double, std::complex<double>, 1, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, double, 0, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, double, 1, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, std::complex<double>, 0, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, std::complex<double>, 1, -1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, std::complex<double>, 0, 1> >
create(Domain<2> const &, double);
template <>
std::auto_ptr<fft::fftm<std::complex<double>, std::complex<double>, 1, 1> >
create(Domain<2> const &, double);
#endif



// Traits class to indicate whether SAL FFTM supports a given type for
// input and output.

template <typename T> struct Is_fft_avail
{ static bool const value = false; };

#if VSIP_IMPL_HAVE_SAL_FLOAT
template <> struct Is_fft_avail<float>
{ static bool const value = true; };

template <> struct Is_fft_avail<complex<float> >
{ static bool const value = true; };
#endif

#if VSIP_IMPL_HAVE_SAL_DOUBLE
template <> struct Is_fft_avail<double>
{ static bool const value = true; };

template <> struct Is_fft_avail<complex<double> >
{ static bool const value = true; };
#endif

} // namespace vsip::impl::sal


namespace fft
{
struct Mercury_sal_tag;

template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<D, I, O, S, R, N, Mercury_sal_tag>
{
  static bool const ct_valid = sal::Is_fft_avail<I>::value &&
                               sal::Is_fft_avail<O>::value;
  static bool rt_valid(Domain<D> const &dom)
  {
    for (dimension_type d = 0; d != D; ++d)
      // SAL can only deal with powers of 2.
      if (dom[d].size() & (dom[d].size() - 1) ||
	  // SAL requires a minimum block size.
	  dom[d].size() < 8) return false;

    return true;
  }
  static std::auto_ptr<backend<D, I, O,
			       axis<I, O, S>::value,
			       exponent<I, O, S>::value> >
  create(Domain<D> const &dom, typename Scalar_of<I>::type scale)
  {
    return sal::create<backend<D, I, O,
                               axis<I, O, S>::value,
                               exponent<I, O, S>::value> >
      (dom, scale);
  }
};

template <dimension_type D,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<D, std::complex<long double>, std::complex<long double>,
		 S, R, N, Mercury_sal_tag>
{
  // No long double support.
  static bool const ct_valid = false;
};

template <dimension_type D,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<D, long double, std::complex<long double>,
		 S, R, N, Mercury_sal_tag>
{
  // No long double support.
  static bool const ct_valid = false;
};

template <dimension_type D,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<D, std::complex<long double>, long double,
		 S, R, N, Mercury_sal_tag>
{
  // No long double support.
  static bool const ct_valid = false;
};

template <typename I,
	  typename O,
	  int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<3, I, O, S, R, N, Mercury_sal_tag>
{
  // No FFT 3D yet.
  static bool const ct_valid = false;
};

template <int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<3, std::complex<long double>, long double,
		 S, R, N, Mercury_sal_tag>
{
  // No FFT 3D yet.
  static bool const ct_valid = false;
};

template <int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<3, long double, std::complex<long double>,
		 S, R, N, Mercury_sal_tag>
{
  // No FFT 3D yet.
  static bool const ct_valid = false;
};

template <int S,
	  return_mechanism_type R,
	  unsigned N>
struct evaluator<3, std::complex<long double>, std::complex<long double>,
		 S, R, N, Mercury_sal_tag>
{
  // No FFT 3D yet.
  static bool const ct_valid = false;
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
struct evaluator<I, O, A, E, R, N, fft::Mercury_sal_tag>
{
  static bool const ct_valid = sal::Is_fft_avail<I>::value &&
                               sal::Is_fft_avail<O>::value;
  static bool rt_valid(Domain<2> const &dom)
  {
    // SAL can only deal with powers of 2.
    if (dom[A].size() & (dom[A].size() - 1)) return false;
    // SAL requires a minimum block size.
    if (dom[A].size() < 8) return false;
    else return true;
  }
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename Scalar_of<I>::type scale)
  {
    return sal::create<fft::fftm<I, O, A, E> >(dom, scale);
  }
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<std::complex<long double>, long double,
		 A, E, R, N, fft::Mercury_sal_tag>
{
  static bool const ct_valid = false;
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<long double, std::complex<long double>,
		 A, E, R, N, fft::Mercury_sal_tag>
{
  static bool const ct_valid = false;
};

template <int A, int E, return_mechanism_type R, unsigned N>
struct evaluator<std::complex<long double>, std::complex<long double>,
		 A, E, R, N, fft::Mercury_sal_tag>
{
  static bool const ct_valid = false;
};
} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif
