/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/

/** @file    tests/fft_inter_split.cpp
    @author  Jules Bergmann
    @date    2007-03-16
    @brief   VSIPL++ Library: Test Fft between split and interleaved views.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

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


// Test FFT by-reference

template <typename T,
	  typename SrcComplexFmt,
	  typename DstComplexFmt>
void
test_fft_br(length_type size)
{
  typedef impl::Stride_unit_dense sud_type;
  typedef impl::Layout<1, row1_type, sud_type, SrcComplexFmt> src_lp_type;
  typedef impl::Layout<1, row1_type, sud_type, DstComplexFmt> dst_lp_type;

  typedef impl::Fast_block<1, T, src_lp_type> src_block_type;
  typedef impl::Fast_block<1, T, dst_lp_type> dst_block_type;

  typedef Fft<const_Vector, T, T, fft_fwd, by_reference, 1, alg_space>
	fft_type;

  fft_type fft(Domain<1>(size), 1.f);

  Vector<T, src_block_type> in(size);
  Vector<T, dst_block_type> out(size);

  in = T(1);

  fft(in, out);

  test_assert(out.get(0) == T(size));
}



// Test FFT by-value

template <typename T,
	  typename SrcComplexFmt,
	  typename DstComplexFmt>
void
test_fft_bv(length_type size)
{
  typedef impl::Stride_unit_dense sud_type;
  typedef impl::Layout<1, row1_type, sud_type, SrcComplexFmt> src_lp_type;
  typedef impl::Layout<1, row1_type, sud_type, DstComplexFmt> dst_lp_type;

  typedef impl::Fast_block<1, T, src_lp_type> src_block_type;
  typedef impl::Fast_block<1, T, dst_lp_type> dst_block_type;

  typedef Fft<const_Vector, T, T, fft_fwd, by_value, 1, alg_space>
	fft_type;

  fft_type fft(Domain<1>(size), 1.f);

  Vector<T, src_block_type> in(size);
  Vector<T, dst_block_type> out(size);

  in = T(1);

  out = fft(in);

  test_assert(out.get(0) == T(size));
}



// Test FFT by-value in an expression

template <typename T,
	  typename SrcComplexFmt,
	  typename DstComplexFmt>
void
test_fft_bv_expr(length_type size)
{
  typedef impl::Stride_unit_dense sud_type;
  typedef impl::Layout<1, row1_type, sud_type, SrcComplexFmt> src_lp_type;
  typedef impl::Layout<1, row1_type, sud_type, DstComplexFmt> dst_lp_type;

  typedef impl::Fast_block<1, T, src_lp_type> src_block_type;
  typedef impl::Fast_block<1, T, dst_lp_type> dst_block_type;

  typedef Fft<const_Vector, T, T, fft_fwd, by_value, 1, alg_space>
	fft_type;

  fft_type fft(Domain<1>(size), 1.f);

  Vector<T, src_block_type> in(size);
  Vector<T, dst_block_type> out(size);

  in = T(1);
  out = T(0);

  out = out + fft(in);

  test_assert(out.get(0) == T(size));
}


template <typename T>
void
test_set(length_type size)
{
  typedef impl::Cmplx_inter_fmt Cif;
  typedef impl::Cmplx_split_fmt Csf;

  test_fft_br<complex<float>, Cif, Cif>(size);
  test_fft_br<complex<float>, Csf, Csf>(size);
  test_fft_br<complex<float>, Cif, Csf>(size);
  test_fft_br<complex<float>, Csf, Cif>(size);

  test_fft_bv<complex<float>, Cif, Cif>(size);
  test_fft_bv<complex<float>, Csf, Csf>(size);
  test_fft_bv<complex<float>, Cif, Csf>(size);
  test_fft_bv<complex<float>, Csf, Cif>(size);

  test_fft_bv_expr<complex<float>, Cif, Cif>(size);
  test_fft_bv_expr<complex<float>, Csf, Csf>(size);
  test_fft_bv_expr<complex<float>, Cif, Csf>(size);
  test_fft_bv_expr<complex<float>, Csf, Cif>(size);
}



/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  // test_set<complex<float> >(16);
  test_set<complex<float> >(256);

  return 0;
}
