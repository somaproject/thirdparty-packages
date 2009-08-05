/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved. */

/** @file    examples/fft.cpp
    @author  Jules Bergmann, Don McCoy
    @date    2006-07-02
    @brief   VSIPL++ Library: Simple FFT example
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

#include <vsip_csl/error_db.hpp>


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
fft_example()
{
  typedef Fft<const_Vector, cscalar_f, cscalar_f, fft_fwd> f_fft_type;
  typedef Fft<const_Vector, cscalar_f, cscalar_f, fft_inv> i_fft_type;

  // Create FFT objects
  vsip::length_type N = 2048;
  f_fft_type f_fft(Domain<1>(N), 1.0);
  i_fft_type i_fft(Domain<1>(N), 1.0 / N);

  // Allocate input and output buffers
  Vector<cscalar_f> in(N, cscalar_f(1.f));
  Vector<cscalar_f> inv(N);
  Vector<cscalar_f> out(N);
  Vector<cscalar_f> ref(N, cscalar_f());
  ref(0) = cscalar_f(N);
  
  // Compute forward and inverse FFT's
  out = f_fft(in);
  inv = i_fft(out);
  
  // Validate the results (allowing for small numerical errors)
  test_assert(error_db(ref, out) < -100);
  test_assert(error_db(inv, in) < -100);

  std::cout << "forward fft mops: " << f_fft.impl_performance("mops") << std::endl;
  std::cout << "inverse fft mops: " << i_fft.impl_performance("mops") << std::endl;
}


int
main(int argc, char **argv)
{
  vsipl init(argc, argv);
  
  Profile profile("fft_profile.txt", pm_accum);

  fft_example();

  return 0;
}
