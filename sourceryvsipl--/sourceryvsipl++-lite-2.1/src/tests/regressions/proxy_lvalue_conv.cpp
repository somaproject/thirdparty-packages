/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/reg_proxy_lvalue_conv.cpp
    @author  Jules Bergmann
    @date    2005-08-08
    @brief   VSIPL++ Library: Regression tests for missing proxy lvalue
	     conversions.
*/

/***********************************************************************
  Included Files
***********************************************************************/

// Force Plain_block to use proxy lvalues
#define PLAINBLOCK_ENABLE_IMPL_REF 0

#include <iostream>
#include <cassert>
#include <vsip/support.hpp>
#include <vsip/initfin.hpp>
#include <vsip/vector.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/plainblock.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions
***********************************************************************/

/// Vector subview of a vector.

template <typename       T,
	  typename       BlockT>
void
test_equal_op(
  Vector<T, BlockT> view)
{
  for (index_type i=0; i<view.size(); ++i)
    view(i) = T(i);

  for (index_type i=0; i<view.size(); ++i)
  {
    test_assert(static_cast<T>(i) == view(i));
    test_assert(view(i) == static_cast<T>(i));
    test_assert(view(i) == T(i));
  }
}


template <typename       T,
	  typename       BlockT>
void
test_equal_fn(
  Vector<T, BlockT> view)
{
  for (index_type i=0; i<view.size(); ++i)
    view(i) = T(i);

  for (index_type i=0; i<view.size(); ++i)
  {
    test_assert(equal(view(i), static_cast<T>(i)));
    test_assert(equal(static_cast<T>(i), view(i)));
    test_assert(equal(view(i), T(i)));
  }
}


int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  // Dense uses real lvalues.
  Vector<float> v_float(5);
  Vector<complex<float> > v_complex(5);

  // Plainblock configured to use proxy lvalues.
  Vector<float, Plain_block<1, float> > v_float_pb(5);
  Vector<complex<float>, Plain_block<1, complex<float> > > v_complex_pb(5);

  test_equal_op(v_float);	// works
  test_equal_fn(v_float);	// works

  test_equal_op(v_complex);	// works
  test_equal_fn(v_complex);	// works

  test_equal_op(v_float_pb);	// works
  test_equal_fn(v_float_pb);	// FAILS
  /* (gcc-4.0) line 64: error: no matching function for call to '
   equal(
      vsip::impl::Lvalue_proxy<
         vsip::Plain_block<1u,float,vsip::row1_type,vsip::Local_map>,
         1
      >,
      float
   )
  */

  test_equal_op(v_complex_pb);	// FAILS
  /* (gcc-4.0) line 47: error: no match for 'operator==' in '

   view. vsip::Vector<T,B>
   ::operator()
    [with T = std::complex<float>,
   Block = vsip::Plain_block<
      1u,
      std::complex<float>,
      vsip::row1_type,
      vsip::Local_map
   >
   ](i)
    == std::complex<float>
   (#'float_expr' not supported by dump_expr#<expression error>, 0.0f)
  */

  test_equal_fn(v_complex_pb);	// FAILS
  /* (gcc-4.0): line 64: error: no matching function for call to '
   equal(
      vsip::impl::Lvalue_proxy<
         vsip::Plain_block<
            1u,
            std::complex<float>,
            vsip::row1_type,
            vsip::Local_map
         >,
         1
      >,
      std::complex<float>
   )
  */


}
