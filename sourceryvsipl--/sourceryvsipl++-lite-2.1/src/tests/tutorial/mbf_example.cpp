/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/tutorial/matlab_bin_formatter.cpp
    @author  Jules Bergmann
    @date    2007-07-13
    @brief   VSIPL++ Library: Test tutorial example for using
             matlab_bin_formatter.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <../doc/users-guide/src/matlab_bin_formatter_prelude.cpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/error_db.hpp>


void
test_write()
{
#include <../doc/users-guide/src/matlab_bin_formatter_write.cpp>
}

void
test_read()
{
#include <../doc/users-guide/src/matlab_bin_formatter_read.cpp>

  // Initialize matrix 'm'.
  Matrix<float> chk_m(3, 3);
  for(index_type i=0;i<3;i++)
    chk_m.row(i) = ramp<float>(3*i, 1, 3);

  // Initialize vector 'v'.
  Vector<float> chk_v(3);
  chk_v = ramp<float>(0, 1, 3);

  test_assert(vsip_csl::error_db(m, chk_m) < -100);
  test_assert(vsip_csl::error_db(v, chk_v) < -100);
}


int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_write();
  test_read();
}


