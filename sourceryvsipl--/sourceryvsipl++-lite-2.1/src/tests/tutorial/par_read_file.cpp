/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/tutorial/par_read_file.cpp
    @author  Jules Bergmann
    @date    2008-01-29
    @brief   VSIPL++ Library: Test tutorial example for read_file example.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <fstream>

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/selgen.hpp>
#include <vsip/math.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/error_db.hpp>

using namespace vsip;



/***********************************************************************
  Main Program
***********************************************************************/

#include <../doc/users-guide/src/par/read_file.hpp>
#include <../doc/users-guide/src/par/write_file.hpp>

int
main(int argc, char** argv)
{
  // Initialize the library.
  vsipl vpp(argc, argv);

  typedef complex<float> value_type;

  // Parameters.
  length_type size = 256;

  // Views.
  Vector<value_type> ref(size);
  Vector<value_type> chk(size, value_type(-100));

  ref = ramp<value_type>(0., 1., size);

  write_file(ref, "read_file-view.raw");
  read_file(chk, "read_file-view.raw");

  float error = vsip_csl::error_db(ref, chk);

  test_assert(error < -150);
}
