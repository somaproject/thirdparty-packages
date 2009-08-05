/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/x-uk-box2.cpp
    @author  Jules Bergmann
    @date    2008-08-05
    @brief   VSIPL++ Library: Test Box2 Ukernel ID1
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/box2.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>


using namespace std;
using namespace vsip;
using namespace vsip_csl;
using vsip_csl::equal;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename ViewT>
typename ViewT::value_type
get(ViewT view, stride_type r, stride_type c, 
    stride_type rr, stride_type cc)
{
  if (r + rr >= 0 && r + rr < (stride_type)view.size(0) &&
      c + cc >= 0 && c + cc < (stride_type)view.size(1))
    return view.get(r + rr, c + cc);
  else
    return typename ViewT::value_type(0);
}

template <typename T>
void
test_ukernel(length_type rows, length_type cols, length_type overlap)
{
  Box2_ukernel cuk(overlap);

  vsip_csl::ukernel::Ukernel<Box2_ukernel> uk(cuk);

  Matrix<T> in(rows, cols);
  Matrix<T> out(rows, cols);

  for (index_type r=0; r<rows; ++r)
    in.row(r) = ramp<T>(T(r), 0.1, cols);

  in = T(1);
  for (index_type r=0; r<rows; ++r)
    in(r, 0) = 100 + r;

  for (index_type r=0; r<rows; ++r)
    in(r, r) = 200 + r;

  uk(in, out);

  int misco = 0;
  for (index_type r=0; r<rows; ++r)
    for (index_type c=0; c<cols; ++c)
    {
      T exp = T(0);
      for (stride_type rr=-overlap; rr<=+(stride_type)overlap; ++rr)
	for (stride_type cc=-overlap; cc<=+(stride_type)overlap; ++cc)
	  exp += get(in, r, c, rr, cc);

      if (!equal(out.get(r, c), exp))
      {
#if VERBOSE >= 2
	std::cout << r << ", " << c << ": " << out.get(r, c) << " != "
		  << exp
		  << std::endl;
#endif
	misco++;
      }
    }
  std::cout << "box2: size " << rows << " x " << cols 
	    << "  overlap " << overlap
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

  test_ukernel<float>(16, 16, 1);
  test_ukernel<float>(16, 16, 2);
  test_ukernel<float>(16, 16, 5);

  test_ukernel<float>(32, 32, 1);
  test_ukernel<float>(32, 32, 2);
  test_ukernel<float>(32, 32, 5);

  test_ukernel<float>(48, 48, 1);
  test_ukernel<float>(48, 48, 2);
  test_ukernel<float>(48, 48, 5);

  test_ukernel<float>(48, 64, 1);
  test_ukernel<float>(48, 64, 2);
  test_ukernel<float>(48, 64, 5);

  test_ukernel<float>(256, 256, 1);
  test_ukernel<float>(256, 256, 2);
  test_ukernel<float>(256, 256, 5);

  return 0;
}
