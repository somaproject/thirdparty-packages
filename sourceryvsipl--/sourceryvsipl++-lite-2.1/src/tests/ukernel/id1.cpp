/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/x-uk-id1.cpp
    @author  Jules Bergmann
    @date    2008-07-10
    @brief   VSIPL++ Library: Test ID1 Ukernel
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/id1.hpp>

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
  Id1_ukernel cuk;

  vsip_csl::ukernel::Ukernel<Id1_ukernel> uk(cuk);

  Vector<T> in(size);
  Vector<T> out(size);

  in = ramp<T>(0, 0.1, size);

  uk(in, out);

  for (index_type i=0; i<size; ++i)
  {
    if (!equal(out.get(i), in.get(i) + T(i)))
    {
      std::cout << i << ": " << out.get(i) << " != "
		<< in.get(i) << " + " << T(i)
		<< std::endl;
    }
    test_assert(equal(out.get(i), in.get(i) + T(i)));
  }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_ukernel<float>(16384);

  return 0;
}
