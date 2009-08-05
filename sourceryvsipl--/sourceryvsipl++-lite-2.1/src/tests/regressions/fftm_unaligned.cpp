/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/fftm_unaligned.cpp
    @author  Jules Bergmann
    @date    2007-03-18
    @brief   VSIPL++ Library: Test Fftm on unaligned views.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>

#include <vsip_csl/test.hpp>

using namespace std;
using namespace vsip;


/***********************************************************************
  Definitions
***********************************************************************/

// Test FFTM by-reference out-of-place with given alignment.

// Requires
//   ROWS to be the number of rows
//   SIZE to be the number of elements per row.
//   GAP to be the extra space between rows.
//   ALIGN to be the alignment of the subview (ALIGN <= GAP)

template <typename T,
	  typename ComplexFmt>
void
test_fftm_op_align(
  length_type rows,
  length_type size,
  length_type gap,
  length_type align)
{
  typedef impl::Stride_unit_dense sud_type;
  typedef impl::Layout<2, row2_type, sud_type, ComplexFmt> lp_type;

  typedef impl::Fast_block<2, T, lp_type> block_type;

  typedef Fftm<T, T, row, fft_fwd, by_reference, 1, alg_space> fftm_type;

  fftm_type fftm(Domain<2>(rows, size), 1.f);

  Matrix<T, block_type> in (rows, size + gap);
  Matrix<T, block_type> out(rows, size + gap);

  in = T(1);

  fftm(in (Domain<2>(rows, Domain<1>(align, 1, size))),
       out(Domain<2>(rows, Domain<1>(align, 1, size))));

  for (index_type i=0; i<rows; ++i)
    test_assert(out.get(i, align) == T(size));
}



// Test FFTM by-reference in-place with given alignment.

// Requires
//   ROWS to be the number of rows
//   SIZE to be the number of elements per row.
//   GAP to be the extra space between rows.
//   ALIGN to be the alignment of the subview (ALIGN <= GAP)

template <typename T,
	  typename ComplexFmt>
void
test_fftm_ip_align(
  length_type rows,
  length_type size,
  length_type gap,
  length_type align)
{
  typedef impl::Stride_unit_dense sud_type;
  typedef impl::Layout<2, row2_type, sud_type, ComplexFmt> lp_type;

  typedef impl::Fast_block<2, T, lp_type> block_type;

  typedef Fftm<T, T, row, fft_fwd, by_reference, 1, alg_space> fftm_type;

  fftm_type fftm(Domain<2>(rows, size), 1.f);

  Matrix<T, block_type> inout(rows, size + gap);

  inout = T(1);

  fftm(inout(Domain<2>(rows, Domain<1>(align, 1, size))));

  for (index_type i=0; i<rows; ++i)
    test_assert(inout.get(i, align) == T(size));
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_fftm_op_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 0, 0);
  test_fftm_op_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 1, 1);
  test_fftm_op_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 1, 0);

  test_fftm_ip_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 0, 0);
  test_fftm_ip_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 1, 1);
  test_fftm_ip_align<complex<float>, impl::Cmplx_inter_fmt>(64, 256, 1, 0);

  return 0;
}
