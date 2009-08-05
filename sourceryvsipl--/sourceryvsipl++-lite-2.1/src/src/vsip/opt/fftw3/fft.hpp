/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/fftw3/fft.hpp
    @author  Stefan Seefeld
    @date    2006-03-06
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with FFTW3.
*/

#ifndef VSIP_OPT_FFTW3_FFT_HPP
#define VSIP_OPT_FFTW3_FFT_HPP

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
namespace fftw3
{

template <typename I, dimension_type D>
std::auto_ptr<I>
create(Domain<D> const &dom, unsigned);

#define VSIP_IMPL_FFT_DECL(D,I,O,A,E)                          \
template <>                                                    \
std::auto_ptr<fft::backend<D,I,O,A,E> >                        \
create(Domain<D> const &, unsigned);

#define VSIP_IMPL_FFT_DECL_T(T)				       \
VSIP_IMPL_FFT_DECL(1, T, std::complex<T>, 0, -1)               \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, T, 0, 1)                \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, std::complex<T>, 0, -1) \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, std::complex<T>, 0, 1)  \
VSIP_IMPL_FFT_DECL(2, T, std::complex<T>, 0, -1)               \
VSIP_IMPL_FFT_DECL(2, T, std::complex<T>, 1, -1)               \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, T, 0, 1)                \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, T, 1, 1)                \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, std::complex<T>, 0, -1) \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, std::complex<T>, 1, -1) \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, std::complex<T>, 0, 1)  \
VSIP_IMPL_FFT_DECL(2, std::complex<T>, std::complex<T>, 1, 1)  \
VSIP_IMPL_FFT_DECL(3, T, std::complex<T>, 0, -1)               \
VSIP_IMPL_FFT_DECL(3, T, std::complex<T>, 1, -1)               \
VSIP_IMPL_FFT_DECL(3, T, std::complex<T>, 2, -1)               \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, T, 0, 1)                \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, T, 1, 1)                \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, T, 2, 1)                \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 0, -1) \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 1, -1) \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 2, -1) \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 0, 1)  \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 1, 1)  \
VSIP_IMPL_FFT_DECL(3, std::complex<T>, std::complex<T>, 2, 1)

#ifdef VSIP_IMPL_FFTW3_HAVE_FLOAT
VSIP_IMPL_FFT_DECL_T(float)
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_DOUBLE
VSIP_IMPL_FFT_DECL_T(double)
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_LONG_DOUBLE
VSIP_IMPL_FFT_DECL_T(long double)
#endif

#undef VSIP_IMPL_FFT_DECL_T
#undef VSIP_IMPL_FFT_DECL

#define VSIP_IMPL_FFT_DECL(I,O,A,E)                            \
template <>                                                    \
std::auto_ptr<fft::fftm<I,O,A,E> >                             \
create(Domain<2> const &, unsigned);

#define VSIP_IMPL_FFT_DECL_T(T)				       \
VSIP_IMPL_FFT_DECL(T, std::complex<T>, 0, -1)                  \
VSIP_IMPL_FFT_DECL(T, std::complex<T>, 1, -1)                  \
VSIP_IMPL_FFT_DECL(std::complex<T>, T, 0, 1)                   \
VSIP_IMPL_FFT_DECL(std::complex<T>, T, 1, 1)                   \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 0, -1)    \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 1, -1)    \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 0, 1)     \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 1, 1)

#ifdef VSIP_IMPL_FFTW3_HAVE_FLOAT
VSIP_IMPL_FFT_DECL_T(float)
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_DOUBLE
VSIP_IMPL_FFT_DECL_T(double)
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_LONG_DOUBLE
VSIP_IMPL_FFT_DECL_T(long double)
#endif

#undef VSIP_IMPL_FFT_DECL_T
#undef VSIP_IMPL_FFT_DECL

} // namespace vsip::impl::fftw3

namespace fft
{
struct Fftw3_tag;

template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<D, I, O, S, R, N, Fftw3_tag>
{
  static bool const ct_valid =
#ifdef VSIP_IMPL_FFTW3_HAVE_FLOAT
    Type_equal<typename Scalar_of<I>::type, float>::value ||
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_DOUBLE
    Type_equal<typename Scalar_of<I>::type, double>::value ||
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_LONG_DOUBLE
    Type_equal<typename Scalar_of<I>::type, long double>::value ||
#endif
    false;

  static bool rt_valid(Domain<D> const &) { return true;}
  static std::auto_ptr<backend<D, I, O,
 			       axis<I, O, S>::value,
 			       exponent<I, O, S>::value> >
  create(Domain<D> const &dom, typename Scalar_of<I>::type)
  {
    return fftw3::create<backend<D, I, O,
      axis<I, O, S>::value,
      exponent<I, O, S>::value> >
      (dom, N);
  }
};

} // namespace vsip::impl::fft

namespace fftm
{
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<I, O, A, E, R, N, fft::Fftw3_tag>
{
  static bool const ct_valid =
#ifdef VSIP_IMPL_FFTW3_HAVE_FLOAT
    Type_equal<typename Scalar_of<I>::type, float>::value ||
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_DOUBLE
    Type_equal<typename Scalar_of<I>::type, double>::value ||
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_LONG_DOUBLE
    Type_equal<typename Scalar_of<I>::type, long double>::value ||
#endif
    false;

  static bool rt_valid(Domain<2> const &/*dom*/) { return true;}
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename impl::Scalar_of<I>::type /*scale*/)
  {
    return fftw3::create<fft::fftm<I, O, A, E> >(dom, N);
  }
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif

