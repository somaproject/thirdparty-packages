/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/x-ukernel.cpp
    @author  Jules Bergmann
    @date    2008-06-10
    @brief   VSIPL++ Library: Test Ukernel
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>

#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/vcopy.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
test_ukernel(length_type size, length_type offset)
{
  Copy_ukernel cuk;

  vsip_csl::ukernel::Ukernel<Copy_ukernel> uk(cuk);

  Vector<T> big_in(size + offset);
  Vector<T> out(size);

  typename Vector<T>::subview_type in(big_in(Domain<1>(offset, 1, size)));

  in = ramp<T>(0, 1, size);

  uk(in, out);

  for (index_type i=0; i<size; ++i)
  {
    if (in.get(i) != out.get(i))
    {
      std::cout << "Miscompare: i = " << i
		<< "  in: " << in.get(i)
		<< "  out: " << out.get(i)
		<< std::endl;
    }
    test_assert(in.get(i) == out.get(i));
  }
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_ukernel<float>(1024, 0);
  test_ukernel<float>(1024, 1);
  test_ukernel<float>(1024, 2);
  test_ukernel<float>(1024, 3);
  test_ukernel<float>(1024, 4);

  test_ukernel<float>(16384, 0);
  test_ukernel<float>(16384, 1);
  test_ukernel<float>(16384, 2);
  test_ukernel<float>(16384, 3);
  test_ukernel<float>(16384, 4);

  return 0;
}
