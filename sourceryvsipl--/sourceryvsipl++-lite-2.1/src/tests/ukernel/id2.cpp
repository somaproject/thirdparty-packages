/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/x-uk-id1.cpp
    @author  Jules Bergmann
    @date    2008-07-29
    @brief   VSIPL++ Library: Test ID2 Ukernel
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip/random.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/id2.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;
using vsip_csl::equal;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
test_ukernel(int shape, length_type rows, length_type cols)
{
  Id2_ukernel cuk(shape, rows, cols);

  vsip_csl::ukernel::Ukernel<Id2_ukernel> uk(cuk);

  Matrix<T> in(rows, cols);
  Matrix<T> out(rows, cols);

  Rand<T> rnd(0);

  in = T(0);
  // in = rnd.randu(rows, cols);

  uk(in, out);

  int misco = 0;
  for (index_type i=0; i<rows; ++i)
    for (index_type j=0; j<cols; ++j)
    {
      if (!equal(out.get(i, j), in.get(i, j) + T(i * cols + j)))
      {
	misco++;
#if VERBOSE >= 2
	std::cout << i << ": " << out.get(i, j) << " != "
		  << in.get(i, j) << " + " << T(i*cols + j)
		  << std::endl;
#endif
      }
  }
  std::cout << "id2: size " << rows << " x " << cols 
	    << "  shape " << shape
	    << "  misco " << misco << std::endl;
  test_assert(misco == 0);
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_ukernel<float>(1, 32, 32);
  test_ukernel<float>(2, 32, 32);
  test_ukernel<float>(3, 32, 32);
  test_ukernel<float>(4, 32, 32);

  test_ukernel<float>(1, 64, 64);
  test_ukernel<float>(2, 64, 64);
  test_ukernel<float>(3, 64, 64);
  test_ukernel<float>(4, 64, 64);

  test_ukernel<float>(0, 1024, 1024);
  test_ukernel<float>(1, 1024, 1024);
  test_ukernel<float>(2, 1024, 1024);
  test_ukernel<float>(3, 1024, 1024);
  test_ukernel<float>(4, 1024, 1024);

  return 0;
}
