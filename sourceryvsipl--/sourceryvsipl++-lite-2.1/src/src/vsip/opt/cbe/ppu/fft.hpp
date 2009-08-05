/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/fft.hpp
    @author  Stefan Seefeld
    @date    2007-01-31
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with the CBE SDK.
*/

#ifndef VSIP_OPT_CBE_PPU_FFT_HPP
#define VSIP_OPT_CBE_PPU_FFT_HPP

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
#include <vsip/opt/cbe/fft_params.h>
#include <vsip/opt/cbe/ppu/task_manager.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cbe
{

template <typename I, dimension_type D, typename S>
std::auto_ptr<I>
create(Domain<D> const &dom, S scale);

#define VSIP_IMPL_FFT_DECL(D,I,O,A,E)                          \
template <>                                                    \
std::auto_ptr<fft::backend<D,I,O,A,E> >                        \
create(Domain<D> const &, fft::backend<D, I, O, A, E>::scalar_type);

#define VSIP_IMPL_FFT_DECL_T(T)				       \
VSIP_IMPL_FFT_DECL(1, T, std::complex<T>, 0, -1)               \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, T, 0, 1)                \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, std::complex<T>, 0, -1) \
VSIP_IMPL_FFT_DECL(1, std::complex<T>, std::complex<T>, 0, 1)

VSIP_IMPL_FFT_DECL_T(float)

#undef VSIP_IMPL_FFT_DECL_T
#undef VSIP_IMPL_FFT_DECL

#define VSIP_IMPL_FFT_DECL(I,O,A,E)                            \
template <>                                                    \
std::auto_ptr<fft::fftm<I,O,A,E> >                             \
create(Domain<2> const &, fft::backend<2, I, O, A, E>::scalar_type);

#define VSIP_IMPL_FFT_DECL_T(T)				       \
VSIP_IMPL_FFT_DECL(T, std::complex<T>, 0, -1)                  \
VSIP_IMPL_FFT_DECL(T, std::complex<T>, 1, -1)                  \
VSIP_IMPL_FFT_DECL(std::complex<T>, T, 0, 1)                   \
VSIP_IMPL_FFT_DECL(std::complex<T>, T, 1, 1)                   \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 0, -1)    \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 1, -1)    \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 0, 1)     \
VSIP_IMPL_FFT_DECL(std::complex<T>, std::complex<T>, 1, 1)

VSIP_IMPL_FFT_DECL_T(float)

#undef VSIP_IMPL_FFT_DECL_T
#undef VSIP_IMPL_FFT_DECL

} // namespace vsip::impl::cbe

struct Cbe_sdk_tag;

namespace fft
{
template <typename I,
	  typename O,
	  int S,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<1, I, O, S, R, N, Cbe_sdk_tag>
{
  static bool const ct_valid = 
    Type_equal<I, complex<float> >::value &&
    Type_equal<I, O>::value;
  static bool rt_valid(Domain<1> const &dom) 
  { 
    return
      (dom.size() >= MIN_FFT_1D_SIZE) &&
      (dom.size() <= MAX_FFT_1D_SIZE) &&
      (fft::is_power_of_two(dom)) &&
      vsip::impl::cbe::Task_manager::instance()->num_spes() > 0;
  }
  static std::auto_ptr<backend<1, I, O,
 			       axis<I, O, S>::value,
 			       exponent<I, O, S>::value> >
  create(Domain<1> const &dom, typename Scalar_of<I>::type scale)
  {
    return cbe::create<backend<1, I, O, 
      axis<I, O, S>::value,
      exponent<I, O, S>::value> >
      (dom, scale);
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
struct evaluator<I, O, A, E, R, N, Cbe_sdk_tag>
{
  static bool const ct_valid = 
    Type_equal<I, complex<float> >::value &&
    Type_equal<I, O>::value;
  static bool rt_valid(Domain<2> const &dom) 
  { 
    length_type size = A ? dom[1].size() : dom[0].size();
    return
      (size >= MIN_FFT_1D_SIZE) &&
      (size <= MAX_FFT_1D_SIZE) &&
      (fft::is_power_of_two(size)) &&
      vsip::impl::cbe::Task_manager::instance()->num_spes() > 0;
  }
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename impl::Scalar_of<I>::type scale)
  {
    return cbe::create<fft::fftm<I, O, A, E> >(dom, scale);
  }
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif

