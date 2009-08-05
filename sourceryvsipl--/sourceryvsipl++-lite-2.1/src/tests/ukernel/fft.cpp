/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernel/fft.cpp
    @author  Jules Bergmann
    @date    2008-06-10
    @brief   VSIPL++ Library: Test Fft Ukernel
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/fft.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;



/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
test_ukernel(length_type rows, length_type cols)
{
  Fft_ukernel cuk(cols, -1, 1.f);

  vsip_csl::ukernel::Ukernel<Fft_ukernel> uk(cuk);

  Matrix<T> in(rows, cols);
  Matrix<T> out(rows, cols);

  in = T(1);

  uk(in, out);

  for (index_type r=0; r<rows; ++r)
  {
    if (out.get(r, 0) != T(cols))
      std::cout << "error row: " << r
		<< " got: " << out.get(r, 0)
		<< " exp: " << T(cols)
		<< std::endl;
    test_assert(out.get(r, 0) == T(cols));
  }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  bool do_large = !VSIP_IMPL_PREFER_SPLIT_COMPLEX;

  if (do_large)
    test_ukernel<complex<float> >(32, 4096);
  test_ukernel<complex<float> >(32, 2048);

  for (index_type i=0; i<100; ++i)
  {
    if (do_large)
      test_ukernel<complex<float> >(32, 4096);
    test_ukernel<complex<float> >(32, 2048);
    test_ukernel<complex<float> >(32, 1024);
  }

  return 0;
}
