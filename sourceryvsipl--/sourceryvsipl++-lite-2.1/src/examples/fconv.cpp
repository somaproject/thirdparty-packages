/* Copyright (c) 2007 by CodeSourcery.  All rights reserved. */

/** @file    examples/fconv.cpp
    @author  Don McCoy
    @date    2007-02-23
    @brief   VSIPL++ Library: Simple fast convolution example
                using the Cell Broadband Engine.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/signal.hpp>
#include <vsip/math.hpp>
#include <vsip/core/profile.hpp>
#ifdef VSIP_IMPL_CBE_SDK
#include <vsip/opt/cbe/ppu/fastconv.hpp>
#endif

#include <vsip_csl/error_db.hpp>
#include <vsip_csl/ref_dft.hpp>


/***********************************************************************
  Definitions
***********************************************************************/

using namespace vsip;
using namespace vsip_csl;
using namespace vsip::impl::profile;


/***********************************************************************
  Functions
***********************************************************************/

void
fconv_example()
{
  typedef std::complex<float> T;
  typedef Fft<const_Vector, T, T, fft_fwd, by_reference> for_fft_type;
  typedef Fft<const_Vector, T, T, fft_inv, by_reference> inv_fft_type;

  const length_type M = 32;
  const length_type N = 64;
  Matrix<T> in(M, N, T(1));
  Matrix<T> out(M, N, T());
  Vector<T> kernel(N, T());  Vector<T> coeffs(4, T());
  Vector<T> replica(N, T());
  Vector<T> ref(N, T());
  Vector<T> tmp(N, T());

  // Create the FFT objects.
  for_fft_type for_fft(Domain<1>(N), 1.0);
  inv_fft_type inv_fft(Domain<1>(N), 1.0/(N));

  // Initialize
  //   Note: the size of coeffs is less than kernel, which
  //   is zero-padded to the length of the input/output vectors.
  kernel(0) = T(1);  coeffs(0) = T(1);
  kernel(1) = T(1);  coeffs(1) = T(1);
  kernel(2) = T(1);  coeffs(2) = T(1);
  kernel(3) = T(1);  coeffs(3) = T(1);

  // Compute reference
  for_fft(kernel, replica);

  for_fft(in.row(0), tmp);
  tmp *= replica;
  inv_fft(tmp, ref);

#ifdef VSIP_IMPL_CBE_SDK
  // Create Fast Convolution object
  typedef vsip::impl::cbe::Fastconv<1, T, vsip::impl::dense_complex_type> fconv_type;
  fconv_type fconv(coeffs, N);


  // Compute convolution on a vector
  fconv(in.row(0), out.row(0));
  test_assert(error_db(ref, out.row(0)) < -100);


  // And now run the convolution over the rows of a matrix
  fconv(in, out);
  for (index_type i = 0; i < M; ++i)
    test_assert(error_db(ref, out.row(i)) < -100);
#endif
}


int
main(int argc, char **argv)
{
  vsipl init(argc, argv);
  
  fconv_example();

  return 0;
}
