/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/coverage_unary_impl.hpp
    @author  Jules Bergmann
    @date    2005-09-13
    @brief   VSIPL++ Library: Coverage tests for Sourcery VSIPL++
             implementation specific unary expressions.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>

#include <vsip/support.hpp>
#include <vsip/initfin.hpp>
#include <vsip/vector.hpp>
#include <vsip/math.hpp>
#include <vsip/random.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/test-storage.hpp>
#include "coverage_common.hpp"

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions
***********************************************************************/

// These C99 functions are unavailable on Windows.
#if !defined(_MSC_VER)
TEST_UNARY(is_nan,    vsip::impl::is_nan,    isnan,    anyval)
TEST_UNARY(is_finite, vsip::impl::is_finite, isfinite, anyval)
TEST_UNARY(is_normal, vsip::impl::is_normal, isnormal, anyval)
#else
TEST_UNARY(is_nan,    vsip::impl::is_nan,   _isnan,    anyval)
#endif


/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  // Unary operators
  vector_cases2_rt<Test_is_nan,    float,  bool>();
#if !defined(_MSC_VER)
  vector_cases2_rt<Test_is_finite, float,  bool>();
  vector_cases2_rt<Test_is_normal, float,  bool>();
#endif
}
