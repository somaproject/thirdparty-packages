/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/tutorial/mtf_example.cpp
    @author  Jules Bergmann
    @date    2007-07-12
    @brief   VSIPL++ Library: Test tutorial example for using
             matlab_test_formatter.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <../doc/users-guide/src/matlab_text_formatter_prelude.cpp>


void
test()
{
#include <../doc/users-guide/src/matlab_text_formatter.cpp>
}


int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test();
}


