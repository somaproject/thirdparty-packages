/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernel/vmmul.cpp
    @author  Jules Bergmann
    @date    2008-06-24
    @brief   VSIPL++ Library: Test Vmul Ukernel
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
#include <vsip/opt/ukernel/kernels/host/vmmul.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using vsip_csl::equal;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
test_ukernel(length_type rows, length_type cols)
{
  Vmmul_kernel obj(cols);

  vsip_csl::ukernel::Ukernel<Vmmul_kernel> uk(obj);

  Vector<T> in0(cols);
  Matrix<T> in1(rows, cols);
  Matrix<T> out(rows, cols);

  Rand<T> gen(0, 0);

  in0 = gen.randu(cols);
  in1 = gen.randu(rows, cols);

  in0(Domain<1>(16)) = ramp<T>(0, 1, 16);
  in1.row(0)(Domain<1>(16)) = ramp<T>(0, 0.1, 16);

  uk(in0, in1, out);

  for (index_type r=0; r<rows; ++r)
    for (index_type c=0; c<cols; ++c)
    {
      if (!equal(in0.get(c) * in1.get(r, c), out.get(r, c)))
      {
	std::cout << "index " << r << ", " << c << ": "
		  << in0.get(c) << " * "
		  << in1.get(r, c) << " = "
		  << in0.get(c) * in1.get(r, c) << "  vs  "
		  << out.get(r, c)
		  << std::endl;
      }
      test_assert(equal(in0.get(c) * in1.get(r, c), out.get(r, c)));
    }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  // test_ukernel<float>();
  test_ukernel<complex<float> >(128, 1024);

  return 0;
}
