/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/tutorial/profile_example.cpp
    @author  Jules Bergmann
    @date    2008-01-29
    @brief   VSIPL++ Library: Test tutorial example for profile example.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/signal.hpp>
#include <vsip/math.hpp>

using namespace vsip;



/***********************************************************************
  Main Program
***********************************************************************/

int
main(int argc, char** argv)
{
  // Initialize the library.
  vsipl vpp(argc, argv);

  typedef complex<float> value_type;

  // Parameters.
  length_type npulse = 64;	// number of pulses
  length_type nrange = 256;	// number of range cells

  // Views.
  Vector<value_type> replica(nrange);
  Matrix<value_type> data(npulse, nrange);
  Matrix<value_type> tmp(npulse, nrange);

  // A forward Fft for computing the frequency-domain version of
  // the replica.
  typedef Fft<const_Vector, value_type, value_type, fft_fwd, by_reference>
		for_fft_type;
  for_fft_type  for_fft (Domain<1>(nrange), 1.0);

  // A forward Fftm for converting the time-domain data matrix to the
  // frequency domain.
  typedef Fftm<value_type, value_type, row, fft_fwd, by_reference>
	  	for_fftm_type;
  for_fftm_type for_fftm(Domain<2>(npulse, nrange), 1.0);

  // An inverse Fftm for converting the frequency-domain data back to
  // the time-domain.
  typedef Fftm<value_type, value_type, row, fft_inv, by_reference>
	  	inv_fftm_type;
  inv_fftm_type inv_fftm(Domain<2>(npulse, nrange), 1.0/(nrange));

  // Initialize data to zero.
  data    = value_type();
  replica = value_type();

  // Before fast convolution, convert the replica to the the
  // frequency domain
  for_fft(replica);

  #include <../doc/users-guide/src/profile_example.cpp>
}
