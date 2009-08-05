/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/regr_ext_subview_split.cpp
    @author  Jules Bergmann
    @date    2005-09-01
    @brief   VSIPL++ Library: Regression test for Ext_data access to
	     real/imaginary subviews of split complex data.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/tensor.hpp>
#include <vsip/core/fast_block.hpp>

#include <vsip_csl/test.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions
***********************************************************************/

template <typename       ComplexFmt,
	  dimension_type Dim>
void
test(stride_type        component_stride,
     Domain<Dim> const& dom)
{
  typedef complex<float> T;

  typedef impl::Layout<Dim, row1_type, impl::Stride_unit_dense, ComplexFmt> 
		                                            layout_t;
  typedef impl::Fast_block<Dim, T, layout_t>                block_t;
  typedef typename impl::View_of_dim<Dim, T, block_t>::type view_t;

  block_t block(dom);
  view_t  view(block);


  impl::Ext_data<block_t> ext(view.block());
  test_assert(ext.cost()        == 0);
  test_assert(ext.stride(Dim-1) == 1);

  stride_type str = 1;
  for (dimension_type d=0; d<Dim; ++d)
  {
    test_assert(ext.stride(Dim-d-1) == str);
    str *= dom[Dim-d-1].size();
  }


  typename view_t::realview_type real = view.real();

  impl::Ext_data<typename view_t::realview_type::block_type>
    extr(real.block());

  test_assert(extr.cost()    == 0);
  test_assert(extr.stride(Dim-1) == component_stride);

  str = component_stride;
  for (dimension_type d=0; d<Dim; ++d)
  {
    test_assert(extr.stride(Dim-d-1) == str);
    str *= dom[Dim-d-1].size();
  }



  typename view_t::imagview_type imag = view.imag();

  impl::Ext_data<typename view_t::imagview_type::block_type>
    exti(imag.block());

  test_assert(exti.cost()    == 0);
  test_assert(exti.stride(Dim-1) == component_stride);

  str = component_stride;
  for (dimension_type d=0; d<Dim; ++d)
  {
    test_assert(exti.stride(Dim-d-1) == str);
    str *= dom[Dim-d-1].size();
  }
}



int
main()
{
  test<impl::Cmplx_inter_fmt>(2, Domain<1>(5));
  test<impl::Cmplx_split_fmt>(1, Domain<1>(5));

  test<impl::Cmplx_inter_fmt>(2, Domain<2>(5, 7));
  test<impl::Cmplx_split_fmt>(1, Domain<2>(5, 7));

  test<impl::Cmplx_inter_fmt>(2, Domain<3>(5, 7, 9));
  test<impl::Cmplx_split_fmt>(1, Domain<3>(5, 7, 9));
}
