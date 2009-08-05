/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/fftw3/fft.cpp
    @author  Stefan Seefeld
    @date    2006-04-10
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             FFTW3.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/support.hpp>
#include <fftw3.h>

// We need to include this create_plan.hpp header file because fft_impl.cpp
// uses this file. We cannot include this file in fft_impl.cpp because
// fft_impl.cpp gets included multiple times here.
#include <vsip/opt/fftw3/create_plan.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace fftw3
{

inline int
convert_NoT(unsigned int number)
{
  // a number value of '0' means 'infinity', and so is captured
  // by a wrap-around.
  if (number - 1 > 30) return FFTW_PATIENT;
  if (number - 1 > 10) return FFTW_MEASURE;
  return FFTW_ESTIMATE;
}

template <dimension_type D, typename I, typename O> struct Fft_base;
template <dimension_type D, typename I, typename O, int A, int E> class Fft_impl;
template <typename I, typename O, int A, int E> class Fftm_impl;

} // namespace vsip::impl::fftw3
} // namespace vsip::impl
} // namespace vsip

#ifdef VSIP_IMPL_FFTW3_HAVE_FLOAT
#  define FFTW(fun) fftwf_##fun
#  define SCALAR_TYPE float
#  include "fft_impl.cpp"
#  undef SCALAR_TYPE
#  undef FFTW
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_DOUBLE
#  define FFTW(fun) fftw_##fun
#  define SCALAR_TYPE double
#  include "fft_impl.cpp"
#  undef SCALAR_TYPE
#  undef FFTW
#endif
#ifdef VSIP_IMPL_FFTW3_HAVE_LONG_DOUBLE
#  define FFTW(fun) fftwl_##fun
#  define SCALAR_TYPE long double
#  include "fft_impl.cpp"
#  undef SCALAR_TYPE
#  undef FFTW
#endif
