/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernel/box1.cpp
    @author  Jules Bergmann
    @date    2008-07-10
    @brief   VSIPL++ Library: Test Box1 Ukernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/box1.hpp>

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
test_ukernel(length_type size, length_type overlap)
{
  Box1_ukernel cuk(overlap, 0); // min-support

  vsip_csl::ukernel::Ukernel<Box1_ukernel> uk(cuk);

  Vector<T> in(size);
  Vector<T> out(size);

  in = ramp<T>(0, 0.1, size);

  uk(in, out);

  int misco = 0;
  for (index_type i=0; i<size; ++i)
  {
    T exp = in.get(i);
    if (i > 0) exp += in.get(i-1);
    if (i < size-1) exp += in.get(i+1);

    if (!equal(out.get(i), exp))
    {
      std::cout << i << ": " << out.get(i) << " != "
		<< exp
		<< std::endl;
      misco++;
    }
  }
  std::cout << "box1: size " << size
	    << "  overlap " << overlap
	    << "  misco " << misco << std::endl;
  test_assert(misco == 0);
}



template <typename T>
void
test_ukernel_subset(length_type size, length_type overlap, length_type offset)
{
  Box1_ukernel cuk(overlap, 1); // same-support

  vsip_csl::ukernel::Ukernel<Box1_ukernel> uk(cuk);

  Vector<T> in(size+offset+overlap);
  Vector<T> out(size);

  Domain<1> dom(offset, 1, size);

  in = ramp<T>(0, 0.1, size + offset + 1);

  uk(in(dom), out);

  int misco = 0;
  for (index_type i=0; i<size; ++i)
  {
    T exp = in.get(offset + i-1) + in.get(offset + i) + in.get(offset + i+1);

    if (!equal(out.get(i), exp))
    {
      std::cout << i << ": " << out.get(i) << " != "
		<< exp
		<< std::endl;
      misco++;
    }
  }
  std::cout << "box1-subset: size " << size
	    << "  overlap " << overlap
	    << "  offset " << offset
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

  test_ukernel<float>(256, 1);

  test_ukernel_subset<float>(256, 1, 1);
  test_ukernel_subset<float>(256, 1, 2);
  test_ukernel_subset<float>(256, 1, 3);
  test_ukernel_subset<float>(256, 1, 4);
  test_ukernel_subset<float>(256, 1, 16);

  return 0;
}
