/* Copyright (c) 2005, 2006, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/extdata_subviews_v.cpp
    @author  Jules Bergmann
    @date    2005-07-22
    @brief   VSIPL++ Library: Unit tests for DDI to vector subviews.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/initfin.hpp>
#include <vsip/dense.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/tensor.hpp>

#include <vsip_csl/plainblock.hpp>

#include "extdata-subviews.hpp"

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
void
vector_test(Domain<1> const& dom)
{
  // Regular Dense
  typedef Dense<1, T>           block_type;
  typedef Vector<T, block_type> view_type;
  view_type view(dom[0].size());

  test_vector(view);
}



template <typename T,
	  typename ComplexFmt>
void
vector_fastblock_test(Domain<1> const& dom)
{
  typedef vsip::impl::Layout<1, row1_type,
    vsip::impl::Stride_unit_dense, ComplexFmt>
				              layout_type;
  typedef vsip::impl::Fast_block<1, T, layout_type> block_type;
  typedef Vector<T, block_type>               view_type;

  view_type view(dom[0].size());
  
  test_vector(view);
}



template <typename T>
void
test_for_type()
{
  vector_test<T>(Domain<1>(7));

  vector_fastblock_test<T, vsip::impl::Cmplx_inter_fmt>(Domain<1>(7));
  vector_fastblock_test<T, vsip::impl::Cmplx_split_fmt>(Domain<1>(7));
}




int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

#if VSIP_IMPL_TEST_LEVEL == 0
  vector_test<float>(Domain<1>(7));
#else
  test_for_type<float>();
  test_for_type<complex<float> >();
#endif

  return 0;
}
