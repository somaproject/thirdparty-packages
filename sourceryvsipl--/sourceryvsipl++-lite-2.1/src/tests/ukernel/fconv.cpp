/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernel/fconv.cpp
    @author  Jules Bergmann
    @date    2008-07-23
    @brief   VSIPL++ Library: Test Fastconv Ukernel
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
#include <vsip/opt/ukernel/kernels/host/fconv.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/error_db.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using vsip_csl::equal;
using vsip_csl::error_db;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
test_ukernel(length_type rows, length_type cols, float scale, int tc)
{
  Fconv_kernel obj(cols);

  vsip_csl::ukernel::Ukernel<Fconv_kernel> uk(obj);

  Vector<T> in0(cols);
  Matrix<T> in1(rows, cols);
  Matrix<T> out(rows, cols, T(-100));

  Rand<T> gen(0, 0);

  in0 = T(scale);

  switch(tc)
  {
  case 0:
    in1 = gen.randu(rows, cols);
    in1.row(0)(Domain<1>(16)) = ramp<T>(0, 0.1, 16);
    break;
  case 1:
    for (index_type r=0; r<rows; ++r)
      in1.row(r) = ramp<T>(r, 0.1, cols);
    break;
  }

  uk(in0, in1, out);

  for (index_type r=0; r<rows; ++r)
  {
    float e = error_db(scale * in1.row(r), out.row(r));
    if (e > -100)
      std::cout << "row " << r << ":  error " << e << std::endl;
    test_assert(e <= -100);
  }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  int tc = 0;

  test_ukernel<complex<float> >(4, 64, 0.5, tc);
  for (index_type i=0; i<100; ++i)
  {
    test_ukernel<complex<float> >(4, 2048, 0.5, tc);
    test_ukernel<complex<float> >(4, 1024, 1.0, tc);
    test_ukernel<complex<float> >(4,  512, 1.5, tc);
    test_ukernel<complex<float> >(4,  512, 3.5, tc);
    test_ukernel<complex<float> >(4, 1024, 2.5, tc);
    test_ukernel<complex<float> >(4,  256, 4.5, tc);
  }

  return 0;
}
