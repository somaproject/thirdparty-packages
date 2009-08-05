/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/diag_eval.cpp
    @author  Jules Bergmann
    @date    2006-11-13
    @brief   VSIPL++ Library: Test Serial_expr_evaluator diagnostics.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/parallel.hpp>
#include <vsip/opt/diag/eval.hpp>

#include <vsip_csl/test.hpp>

using namespace std;
using namespace vsip;


/***********************************************************************
  Definitions
***********************************************************************/


template <typename T>
void
test_expr()
{
  length_type size = 64;

  Vector<T> A(size, T(1));
  Vector<T> B(size, T(2));
  Vector<T> Z(size, T(0));

  // Diagnose how dispatch will handle Z = A + B.
  vsip::impl::diagnose_eval_dispatch(Z, A + B);

  // Diagnose how tags in standard list will handle Z = A + B.
  // (assumes that Z = A + B is a local expression)
  vsip::impl::diagnose_eval_list_std(Z, A + B);

  // Zoom in on the Intel_ipp_tag evaluator.
  vsip::impl::diagnose_eval_tag<vsip::impl::Loop_fusion_tag>(Z, A + B);
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_expr<complex<float> >();

  return 0;
}
