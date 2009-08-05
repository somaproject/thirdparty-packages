/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/fast-block.cpp
    @author  Jules Bergmann
    @date    2005-04-12
    @brief   VSIPL++ Library: Unit tests for Fast_blocks.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vsip/support.hpp>
#include <vsip/core/fast_block.hpp>
#include <vsip/core/length.hpp>
#include <vsip/core/domain_utils.hpp>
#include <vsip_csl/test.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;

using vsip::impl::Length;
using vsip::impl::extent;
using vsip::impl::valid;
using vsip::impl::next;



/***********************************************************************
  Definitions
***********************************************************************/


template <typename T>
inline T
identity(
  Length<1> /*extent*/,
  Index<1> idx,
  int      k)
{
  return static_cast<T>(k*idx[0] + 1);
}


template <typename T>
inline T
identity(
  Length<2> extent,
  Index<2>  idx,
  int       k)
{
  Index<2> offset;
  index_type i = (idx[0]+offset[0])*extent[1] + (idx[1]+offset[1]);
  return static_cast<T>(k*i+1);
}


template <dimension_type Dim,
	  typename       Block>
void
fill_block(Block& blk, int k)
{
  typedef typename Block::value_type value_type;

  Length<Dim> ex = extent<Dim>(blk);
  for (Index<Dim> idx; valid(ex,idx); next(ex, idx))
  {
    put(blk, idx, identity<value_type>(ex, idx, k));
  }
}


template <dimension_type Dim,
	  typename       Block>
void
check_block(Block& blk, int k)
{
  typedef typename Block::value_type value_type;

  Length<Dim> ex = extent<Dim>(blk);
  for (Index<Dim> idx; valid(ex,idx); next(ex, idx))
  {
    test_assert(equal( get(blk, idx),
		  identity<value_type>(ex, idx, k)));
  }
}


template <dimension_type Dim,
	  typename       Block>
void
test(Domain<Dim> const& dom)
{
  Block block(dom);

  fill_block<Dim>(block, 5);
  check_block<Dim>(block, 5);
}

template <typename T,
	  typename PackType,
	  typename CmplxType>
void
test_wrap()
{
  using vsip::impl::Layout;

  Domain<1> dom1(15);
  test<1, vsip::impl::Fast_block<1, T, Layout<1, row1_type, PackType, CmplxType> > >(dom1);

  Domain<2> dom2(15, 17);
  test<2, vsip::impl::Fast_block<2, T, Layout<2, row2_type, PackType, CmplxType> > >(dom2);
  test<2, vsip::impl::Fast_block<2, T, Layout<2, col2_type, PackType, CmplxType> > >(dom2);
}


int
main()
{
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Stride_unit_align;
  using vsip::impl::Cmplx_inter_fmt;
  using vsip::impl::Cmplx_split_fmt;

  test_wrap<float, Stride_unit_dense,     Cmplx_inter_fmt>();
  test_wrap<float, Stride_unit_align<16>, Cmplx_inter_fmt>();
  test_wrap<float, Stride_unit_dense,     Cmplx_split_fmt>();
  test_wrap<float, Stride_unit_align<16>, Cmplx_split_fmt>();

  test_wrap<complex<float>, Stride_unit_dense,     Cmplx_inter_fmt>();
  test_wrap<complex<float>, Stride_unit_align<16>, Cmplx_inter_fmt>();
  test_wrap<complex<float>, Stride_unit_dense,     Cmplx_split_fmt>();
  test_wrap<complex<float>, Stride_unit_align<16>, Cmplx_split_fmt>();
}
