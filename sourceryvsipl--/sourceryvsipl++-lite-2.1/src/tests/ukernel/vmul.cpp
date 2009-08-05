/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/x-uk-vmul.cpp
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
#include <vsip/opt/ukernel/kernels/host/vmul.hpp>

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
test_ukernel(length_type size)
{
  Vmul_kernel obj;

  vsip_csl::ukernel::Ukernel<Vmul_kernel> uk(obj);

  Vector<T> in0(size);
  Vector<T> in1(size);
  Vector<T> out(size);

  Rand<T> gen(0, 0);

  in0 = gen.randu(size);
  in1 = gen.randu(size);

  in0(Domain<1>(16)) = ramp<T>(0, 1, 16);
  in1(Domain<1>(16)) = ramp<T>(0, 0.1, 16);

  uk(in0, in1, out);

  for (index_type i=0; i<size; ++i)
  {
    if (!equal(in0.get(i) * in1.get(i), out.get(i)))
    {
      std::cout << "index " << i << ": "
		<< in0.get(i) << " * "
		<< in1.get(i) << " = "
		<< in0.get(i) * in1.get(i) << "  vs  "
		<< out.get(i)
		<< std::endl;
    }
    test_assert(equal(in0.get(i) * in1.get(i), out.get(i)));
  }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_ukernel<float>(1024+32);
  test_ukernel<float>(16384);
  test_ukernel<complex<float> >(16384);

  return 0;
}
