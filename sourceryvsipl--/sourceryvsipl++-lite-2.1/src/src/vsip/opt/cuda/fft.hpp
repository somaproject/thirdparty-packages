/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/fft.hpp
    @author  Don McCoy
    @date    2009-02-26
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             NVidia's CUDA FFT library.
*/

#ifndef VSIP_IMPL_CUDA_FFT_HPP
#define VSIP_IMPL_CUDA_FFT_HPP

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

#define CUFFT_1D_TRANSFORM_LIMIT (8e6)


namespace vsip
{
namespace impl
{
namespace cuda
{

/// These are the entry points into the CUDA FFT bridge.
template <typename I, dimension_type D, typename S>
std::auto_ptr<I>
create(Domain<D> const &dom, S scale);

}



namespace fft
{

/// Tag used for CUDA FFT backend
struct Cuda_tag;

template <dimension_type D, typename I, typename O, int S>
struct base_evaluator
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<D> const &dom)
  {
    return (dom.size() <= CUFFT_1D_TRANSFORM_LIMIT);
  }

  static std::auto_ptr<fft::backend<D, I, O,
				    fft::axis<I, O, S>::value,
				    fft::exponent<I, O, S>::value> >
  create(Domain<D> const &dom, typename Scalar_of<I>::type scale)
  {
    return cuda::create<fft::backend<D, I, O,
                       fft::axis<I, O, S>::value,
                       fft::exponent<I, O, S>::value> >
      (dom, scale);
  }
};


// Define only those evaluators that are currently supported.

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, float, std::complex<float>, S, R, N, Cuda_tag>
  : base_evaluator<1, float, std::complex<float>, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<float>, float, S, R, N, Cuda_tag>
  : base_evaluator<1, std::complex<float>, float, S> {};

template <int S, return_mechanism_type R, unsigned N>
struct evaluator<1, std::complex<float>, std::complex<float>, S, R, N,
		 Cuda_tag>
  : base_evaluator<1, std::complex<float>, std::complex<float>, S> {};

} // namespace vsip::impl::fft




namespace fftm
{

template <typename I,
	  typename O,
	  int A,
	  int E>
struct base_evaluator
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<2> const &dom)
  {
    return (dom.size() <= CUFFT_1D_TRANSFORM_LIMIT);
  }

  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &dom, typename Scalar_of<I>::type scale)
  {
    return cuda::create<fft::fftm<I, O, A, E> >(dom, scale);
  }
};


// Define only those evaluators that are currently supported.

template <int A, int E, vsip::return_mechanism_type R, unsigned N>
struct evaluator<float, std::complex<float>, A, E, R, N, 
                 fft::Cuda_tag>
  : base_evaluator<float, std::complex<float>, A, E> {};

template <int A, int E, vsip::return_mechanism_type R, unsigned N>
struct evaluator<std::complex<float>, float, A, E, R, N, 
                 fft::Cuda_tag>
  : base_evaluator<std::complex<float>, float, A, E> {};

template <int A, int E, vsip::return_mechanism_type R, unsigned N>
struct evaluator<std::complex<float>, std::complex<float>, A, E, R, N, 
                 fft::Cuda_tag>
  : base_evaluator<std::complex<float>, std::complex<float>, A, E> {};

} // namespace vsip::impl::fftm

} // namespace vsip::impl
} // namespace vsip

#endif

