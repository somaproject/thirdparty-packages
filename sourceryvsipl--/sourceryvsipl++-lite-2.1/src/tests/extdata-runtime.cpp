/* Copyright (c) 2005, 2006, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/extdata-runtime.cpp
    @author  Jules Bergmann
    @date    2005-07-26
    @brief   VSIPL++ Library: Unit tests for runtime determination of
             direct versus copy access.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vsip/support.hpp>
#include <vsip/initfin.hpp>
#include <vsip/dense.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/tensor.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/plainblock.hpp>
#include "extdata-output.hpp"

using namespace std;
using namespace vsip;
using namespace vsip_csl;

using vsip::impl::ITE_Type;
using vsip::impl::Type_equal;
using vsip::impl::As_type;

using vsip::impl::Cmplx_inter_fmt;
using vsip::impl::Cmplx_split_fmt;

length_type g_size = 10;

length_type g_rows = 10;
length_type g_cols = 15;

length_type g_dim0 = 8;
length_type g_dim1 = 12;
length_type g_dim2 = 14;



/***********************************************************************
  Definitions
***********************************************************************/

template <typename BlockT>
void
dump_access_details()
{
  cout << "Access details (block_type = " << Type_name<BlockT>::name() << endl;

  typedef vsip::impl::Desired_block_layout<BlockT> dbl_type;
  typedef typename dbl_type::access_type access_type;
  typedef typename dbl_type::order_type  order_type;
  typedef typename dbl_type::pack_type   pack_type;
  typedef typename dbl_type::layout_type layout_type;

  cout << "  dbl access_type = " << Type_name<access_type>::name() << endl;
  cout << "  dbl order_type  = " << Type_name<order_type>::name() << endl;
  cout << "  dbl pack_type   = " << Type_name<pack_type>::name() << endl;

  typedef vsip::impl::Block_layout<BlockT> bl_type;
  typedef typename bl_type::access_type bl_access_type;
  typedef typename bl_type::order_type  bl_order_type;
  typedef typename bl_type::pack_type   bl_pack_type;
  typedef typename bl_type::layout_type bl_layout_type;

  cout << "  bl  access_type = " << Type_name<bl_access_type>::name() << endl;
  cout << "  bl  order_type  = " << Type_name<bl_order_type>::name() << endl;
  cout << "  bl  pack_type   = " << Type_name<bl_pack_type>::name() << endl;

  typedef typename vsip::impl::Choose_access<BlockT, layout_type>::type
    use_access_type;

  cout << "  use_access_type = " << Type_name<use_access_type>::name() << endl;

  cout << "  cost            = " << vsip::impl::Ext_data_cost<BlockT>::value
       << endl;
}



struct Same {};
struct Diff {};

struct Full {};		// full dimension
struct Cont {};		// continuous - first half of dimension
struct Off1 {};		// ^ ^ offset by 1
struct Spar {};		// sparse - every other element
struct Spar4 {};		// sparse - every other element
struct Sing {};		// single - one element



/***********************************************************************
  Vector test harnass
***********************************************************************/

template <typename PackT,
	  typename CmplxT,
	  typename T,
	  typename BlockT>
void
test_vector(
  int               ct_cost,
  int               rt_cost,
  Vector<T, BlockT> view)
{
  // length_type size = view.size(0);

  dimension_type const dim = BlockT::dim;
  typedef typename vsip::impl::Block_layout<BlockT>::access_type  access_type;
  typedef typename vsip::impl::Block_layout<BlockT>::order_type   order_type;
  typedef typename vsip::impl::Block_layout<BlockT>::pack_type    blk_pack_type;
  typedef PackT                                             pack_type;
  typedef CmplxT                                            complex_type;

  typedef vsip::impl::Layout<dim, order_type, pack_type, complex_type> LP;

  typedef typename vsip::impl::Choose_access<BlockT, LP>::type real_access_type;

  for (index_type i=0; i<view.size(0); ++i)
    view.put(i, T(i));

  {
    vsip::impl::Ext_data<BlockT, LP> ext(view.block());

#if 0
    cout << "Block (" << Type_name<BlockT>::name() << ")" << endl
	 << "  Blk AT   = " << Type_name<access_type>::name() << endl
	 << "  Blk pack = " << Type_name<blk_pack_type>::name() << endl

	 << "  Req AT   = " << Type_name<real_access_type>::name() << endl
	 << "  Req pack = " << Type_name<pack_type>::name() << endl
	 << "  ct_cost  = " << vsip::impl::Ext_data_cost<BlockT, LP>::value << endl
	 << "  rt_cost  = " << ext.cost() << endl
      ;
#endif
    
    test_assert((ct_cost == vsip::impl::Ext_data_cost<BlockT, LP>::value));
    test_assert(rt_cost == ext.cost());
    test_assert(rt_cost == vsip::impl::cost<LP>(view.block()));

    // Check that rt_cost == 0 implies mem_required == 0
    test_assert((rt_cost == 0 && vsip::impl::mem_required<LP>(view.block()) == 0) ||
		(rt_cost != 0 && vsip::impl::mem_required<LP>(view.block()) > 0));

    // Check that rt_cost == 0 implies xfer_required == false
    test_assert((rt_cost == 0 && !vsip::impl::xfer_required<LP>(view.block())) ||
		(rt_cost != 0 &&  vsip::impl::xfer_required<LP>(view.block())) );

    test_assert(ext.size(0) == view.size(0));

    T* ptr              = ext.data();
    stride_type stride0 = ext.stride(0);

    for (index_type i=0; i<ext.size(0); ++i)
    {
      test_assert(equal(ptr[i*stride0], T(i)));
      ptr[i*stride0] = T(i+100);
      }
  }

  for (index_type i=0; i<view.size(0); ++i)
    test_assert(equal(view.get(i), T(i+100)));
}



template <typename BlockT,
	  typename Dim0,
	  typename PackT,
	  typename CmplxT>
struct Test_vector
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Vector<typename BlockT::value_type, BlockT> view(g_size);

    typedef Domain<1> D1;

    D1 dom0(view.size(0));

    if      (Type_equal<Dim0, Cont>::value) dom0 = D1(0, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Off1>::value) dom0 = D1(1, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Spar>::value) dom0 = D1(0, 2, view.size(0)/2);

    test_vector<PackT, CmplxT>(ct_cost, rt_cost, view(dom0));
  }
};


template <typename BlockT,
	  typename PackT,
	  typename CmplxT>
struct Test_vector<BlockT, Full, PackT, CmplxT>
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Vector<typename BlockT::value_type, BlockT> view(g_size);

    test_vector<PackT, CmplxT>(ct_cost, rt_cost, view);
  }
};



template <typename BlockT,
	  typename Dim0,
	  typename ReqPackT,
	  typename ReqCmplxT,
	  
	  typename BlkPackT =
	  typename vsip::impl::Block_layout<BlockT>::pack_type,

	  typename BlkCmplxT =
	  typename vsip::impl::Block_layout<BlockT>::complex_type,

	  typename ActPackT =
	  typename ITE_Type<Type_equal<ReqPackT, Same>::value,
			    As_type<BlkPackT>, As_type<ReqPackT> >::type,

	  typename ActCmplxT = typename
	  ITE_Type<Type_equal<ReqCmplxT, Same>::value, As_type<BlkCmplxT>,
	  ITE_Type<Type_equal<ReqCmplxT, Diff>::value,
		   ITE_Type<Type_equal<BlkCmplxT, Cmplx_split_fmt>::value,
			    As_type<Cmplx_inter_fmt>, As_type<Cmplx_split_fmt> >,
	           As_type<ReqCmplxT> > >::type>
struct Tv : Test_vector<BlockT, Dim0, ActPackT, ActCmplxT>
{};



/***********************************************************************
  Matrix test harnass
***********************************************************************/



template <typename PackT,
	  typename OrderT,
	  typename CmplxT,
	  typename T,
	  typename BlockT>
void
test_matrix(
  int               ct_cost,
  int               rt_cost,
  Matrix<T, BlockT> view)
{
  // length_type size = view.size(0);

  dimension_type const dim = BlockT::dim;
  typedef typename vsip::impl::Block_layout<BlockT>::access_type  access_type;
  typedef OrderT  order_type;
  typedef PackT   pack_type;
  typedef CmplxT  complex_type;

  typedef vsip::impl::Layout<dim, order_type, pack_type, complex_type> LP;

  typedef typename vsip::impl::Choose_access<BlockT, LP>::type real_access_type;

  for (index_type i=0; i<view.size(0); ++i)
    for (index_type j=0; j<view.size(1); ++j)
      view.put(i, j, T(i*view.size(1)+j));

  {
    vsip::impl::Ext_data<BlockT, LP> ext(view.block());

#if 0
    cout << "Block (" << Type_name<BlockT>::name() << ")" << endl
	 << "  AT      = " << Type_name<access_type>::name() << endl
	 << "  RAT     = " << Type_name<real_access_type>::name() << endl
	 << "  ct_cost = " << vsip::impl::Ext_data_cost<BlockT, LP>::value << endl
	 << "  rt_cost = " << ext.cost() << endl
      ;
#endif
    
    test_assert((ct_cost == vsip::impl::Ext_data_cost<BlockT, LP>::value));
    test_assert(rt_cost == ext.cost());

    // Check that rt_cost == 0 implies mem_required == 0
    test_assert((rt_cost == 0 && vsip::impl::mem_required<LP>(view.block()) == 0) ||
		(rt_cost != 0 && vsip::impl::mem_required<LP>(view.block()) > 0));

    // Check that rt_cost == 0 implies xfer_required == false
    test_assert((rt_cost == 0 && !vsip::impl::xfer_required<LP>(view.block())) ||
		(rt_cost != 0 &&  vsip::impl::xfer_required<LP>(view.block())) );

    test_assert(ext.size(0) == view.size(0));
    test_assert(ext.size(1) == view.size(1));

    T* ptr              = ext.data();
    stride_type stride0 = ext.stride(0);
    stride_type stride1 = ext.stride(1);

    for (index_type i=0; i<ext.size(0); ++i)
      for (index_type j=0; j<ext.size(1); ++j)
      {
	test_assert(equal(ptr[i*stride0 + j*stride1], T(i*view.size(1)+j)));
	ptr[i*stride0 + j*stride1] = T(i+j*view.size(0));
      }
  }

  for (index_type i=0; i<view.size(0); ++i)
    for (index_type j=0; j<view.size(1); ++j)
      test_assert(equal(view.get(i, j), T(i+j*view.size(0))));
}



template <typename BlockT,
	  typename Dim0,
	  typename Dim1,
	  typename PackT,
	  typename OrderT,
	  typename CmplxT>
struct Test_matrix
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Matrix<typename BlockT::value_type, BlockT> view(g_rows, g_cols);

    typedef Domain<1> D1;

    Domain<1> dom0(view.size(0));
    Domain<1> dom1(view.size(1));

    if      (Type_equal<Dim0, Cont>::value)  dom0 = D1(0, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Off1>::value)  dom0 = D1(1, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Spar>::value)  dom0 = D1(0, 2, view.size(0)/2);
    else if (Type_equal<Dim0, Spar4>::value) dom0 = D1(0, 4, view.size(0)/4);
    else if (Type_equal<Dim0, Sing>::value)  dom0 = D1(0, 1, 1);

    if      (Type_equal<Dim1, Cont>::value)  dom1 = D1(0, 1, view.size(1)/2);
    else if (Type_equal<Dim1, Off1>::value)  dom1 = D1(1, 1, view.size(1)/2);
    else if (Type_equal<Dim1, Spar>::value)  dom1 = D1(0, 2, view.size(1)/2);
    else if (Type_equal<Dim1, Spar4>::value) dom1 = D1(0, 4, view.size(1)/4);
    else if (Type_equal<Dim1, Sing>::value)  dom1 = D1(0, 1, 1);

    Domain<2> dom(dom0, dom1);

    test_matrix<PackT, OrderT, CmplxT>(ct_cost, rt_cost, view(dom));
  }
};


template <typename BlockT,
	  typename PackT,
	  typename OrderT,
	  typename CmplxT>
struct Test_matrix<BlockT, Full, Full, PackT, OrderT, CmplxT>
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Matrix<typename BlockT::value_type, BlockT> view(g_rows, g_cols);

    test_matrix<PackT, OrderT, CmplxT>(ct_cost, rt_cost, view);
  }
};


template <typename BlockT,
	  typename Dim0,
	  typename Dim1,
	  typename ReqPackT,
	  typename ReqOrderT,
	  typename ReqCmplxT,
	  
	  typename BlkPackT =
	  typename vsip::impl::Block_layout<BlockT>::pack_type,

	  typename BlkOrderT =
	  typename vsip::impl::Block_layout<BlockT>::order_type,

	  typename BlkCmplxT =
	  typename vsip::impl::Block_layout<BlockT>::complex_type,

	  typename ActPackT =
	  typename ITE_Type<Type_equal<ReqPackT, Same>::value,
			    As_type<BlkPackT>, As_type<ReqPackT> >::type,

	  typename ActOrderT = typename
	  ITE_Type<Type_equal<ReqOrderT, Same>::value, As_type<BlkOrderT>,
	  ITE_Type<Type_equal<ReqOrderT, Diff>::value,
		   ITE_Type<Type_equal<BlkOrderT, row2_type>::value,
			    As_type<col2_type>, As_type<row2_type> >,
	           As_type<ReqOrderT> > >::type,

	  typename ActCmplxT = typename
	  ITE_Type<Type_equal<ReqCmplxT, Same>::value, As_type<BlkCmplxT>,
	  ITE_Type<Type_equal<ReqCmplxT, Diff>::value,
		   ITE_Type<Type_equal<BlkCmplxT, Cmplx_split_fmt>::value,
			    As_type<Cmplx_inter_fmt>, As_type<Cmplx_split_fmt> >,
	           As_type<ReqOrderT> > >::type>
struct Tm : Test_matrix<BlockT, Dim0, Dim1, ActPackT, ActOrderT, ActCmplxT>
{};

	  





/***********************************************************************
  Tensor test harnass
***********************************************************************/

template <typename PackT,
	  typename OrderT,
	  typename ComplexT,
	  typename T,
	  typename BlockT>
void
test_tensor(
  int               ct_cost,
  int               rt_cost,
  Tensor<T, BlockT> view)
{
  dimension_type const dim = BlockT::dim;
  typedef typename vsip::impl::Block_layout<BlockT>::access_type  access_type;
  typedef OrderT                                            order_type;
  typedef PackT                                             pack_type;
  typedef ComplexT                                          complex_type;

  typedef vsip::impl::Layout<dim, order_type, pack_type, complex_type> LP;

  typedef typename vsip::impl::Choose_access<BlockT, LP>::type real_access_type;

  for (index_type i=0; i<view.size(0); ++i)
    for (index_type j=0; j<view.size(1); ++j)
      for (index_type k=0; k<view.size(2); ++k)
	view.put(i, j, k, T(i*view.size(1)*view.size(2) + j*view.size(2) + k));

  {
    vsip::impl::Ext_data<BlockT, LP> ext(view.block());

#if 0
    cout << "Block (" << Type_name<BlockT>::name() << ")" << endl
	 << "  AT      = " << Type_name<access_type>::name() << endl
	 << "  RAT     = " << Type_name<real_access_type>::name() << endl
	 << "  ct_cost = " << vsip::impl::Ext_data_cost<BlockT, LP>::value << endl
	 << "  rt_cost = " << ext.cost() << endl
      ;
#endif
    
    test_assert((ct_cost == vsip::impl::Ext_data_cost<BlockT, LP>::value));
    test_assert(rt_cost == ext.cost());

    // Check that rt_cost == 0 implies mem_required == 0
    test_assert((rt_cost == 0 && vsip::impl::mem_required<LP>(view.block()) == 0) ||
		(rt_cost != 0 && vsip::impl::mem_required<LP>(view.block()) > 0));

    // Check that rt_cost == 0 implies xfer_required == false
    test_assert((rt_cost == 0 && !vsip::impl::xfer_required<LP>(view.block())) ||
		(rt_cost != 0 &&  vsip::impl::xfer_required<LP>(view.block())) );

    test_assert(ext.size(0) == view.size(0));
    test_assert(ext.size(1) == view.size(1));
    test_assert(ext.size(2) == view.size(2));

    T* ptr              = ext.data();
    stride_type stride0 = ext.stride(0);
    stride_type stride1 = ext.stride(1);
    stride_type stride2 = ext.stride(2);

    for (index_type i=0; i<ext.size(0); ++i)
      for (index_type j=0; j<ext.size(1); ++j)
	for (index_type k=0; k<ext.size(2); ++k)
	{
	  test_assert(equal(ptr[i*stride0 + j*stride1 + k*stride2],
		       T(i*view.size(1)*view.size(2) + j*view.size(2) + k)));
	  ptr[i*stride0 + j*stride1 + k*stride2] =
	    T(i+j*view.size(0)+k*view.size(0)*view.size(1));
      }
  }

  for (index_type i=0; i<view.size(0); ++i)
    for (index_type j=0; j<view.size(1); ++j)
      for (index_type k=0; k<view.size(2); ++k)
	test_assert(equal(view.get(i, j, k),
		     T(i+j*view.size(0)+k*view.size(0)*view.size(1))));
}



template <typename BlockT,
	  typename Dim0,
	  typename Dim1,
	  typename Dim2,
	  typename PackT,
	  typename OrderT,
	  typename CmplxT>
struct Test_tensor
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Tensor<typename BlockT::value_type, BlockT> view(g_dim0, g_dim1, g_dim2);

    typedef Domain<1> D1;

    Domain<1> dom0(view.size(0));
    Domain<1> dom1(view.size(1));
    Domain<1> dom2(view.size(2));

    if      (Type_equal<Dim0, Cont>::value) dom0 = D1(0, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Off1>::value) dom0 = D1(1, 1, view.size(0)/2);
    else if (Type_equal<Dim0, Spar>::value) dom0 = D1(0, 2, view.size(0)/2);
    else if (Type_equal<Dim0, Sing>::value) dom0 = D1(0, 1, 1);

    if      (Type_equal<Dim1, Cont>::value) dom1 = D1(0, 1, view.size(1)/2);
    else if (Type_equal<Dim1, Off1>::value) dom1 = D1(1, 1, view.size(1)/2);
    else if (Type_equal<Dim1, Spar>::value) dom1 = D1(0, 2, view.size(1)/2);
    else if (Type_equal<Dim1, Sing>::value) dom1 = D1(1, 1, 1);

    if      (Type_equal<Dim2, Cont>::value) dom2 = D1(0, 1, view.size(2)/2);
    else if (Type_equal<Dim2, Off1>::value) dom2 = D1(1, 1, view.size(2)/2);
    else if (Type_equal<Dim2, Spar>::value) dom2 = D1(0, 2, view.size(2)/2);
    else if (Type_equal<Dim2, Sing>::value) dom2 = D1(1, 1, 1);

    Domain<3> dom(dom0, dom1, dom2);

    test_tensor<PackT, OrderT, CmplxT>(ct_cost, rt_cost, view(dom));
  }
};


template <typename BlockT,
	  typename PackT,
	  typename OrderT,
	  typename CmplxT>
struct Test_tensor<BlockT, Full, Full, Full, PackT, OrderT, CmplxT>
{
  static void
  test(int ct_cost, int rt_cost)
  {
    Tensor<typename BlockT::value_type, BlockT> view(g_dim0, g_dim1, g_dim2);

    test_tensor<PackT, OrderT, CmplxT>(ct_cost, rt_cost, view);
  }
};



template <typename BlockT,
	  typename Dim0,
	  typename Dim1,
	  typename Dim2,
	  typename ReqPackT,
	  typename ReqOrderT,
	  typename ReqCmplxT,
	  
	  typename BlkPackT =
	  typename vsip::impl::Block_layout<BlockT>::pack_type,

	  typename BlkOrderT =
	  typename vsip::impl::Block_layout<BlockT>::order_type,

	  typename BlkCmplxT =
	  typename vsip::impl::Block_layout<BlockT>::complex_type,

	  typename ActPackT =
	  typename ITE_Type<Type_equal<ReqPackT, Same>::value,
			    As_type<BlkPackT>, As_type<ReqPackT> >::type,

	  typename ActOrderT = typename
	  ITE_Type<Type_equal<ReqOrderT, Same>::value, As_type<BlkOrderT>,
	  ITE_Type<Type_equal<ReqOrderT, Diff>::value,
		   ITE_Type<Type_equal<BlkOrderT, row2_type>::value,
			    As_type<col2_type>, As_type<row2_type> >,
	           As_type<ReqOrderT> > >::type,

	  typename ActCmplxT = typename
	  ITE_Type<Type_equal<ReqCmplxT, Same>::value, As_type<BlkCmplxT>,
	  ITE_Type<Type_equal<ReqCmplxT, Diff>::value,
		   ITE_Type<Type_equal<BlkCmplxT, Cmplx_split_fmt>::value,
			    As_type<Cmplx_inter_fmt>, As_type<Cmplx_split_fmt> >,
	           As_type<ReqOrderT> > >::type>
struct Tt : Test_tensor<BlockT, Dim0, Dim1, Dim2,
			ActPackT, ActOrderT, ActCmplxT>
{};



/***********************************************************************
  Tests
***********************************************************************/

void
vector_tests()
{
  using vsip::impl::Stride_unknown;
  using vsip::impl::Stride_unit;
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Stride_unit_align;

  typedef Dense<1, float> blk_t;

  // Blk  , Dim0, Pack          , Cplx
  // -----,-----,---------------,-----

  // Asking for Stride_unknown packing indicates we don't care.
  //   If we ask for same complex format => get direct access
  //   If we ask for different complex format => get copy access
  Tv<blk_t, Full, Stride_unknown, Same>::test(0, 0);
  Tv<blk_t, Cont, Stride_unknown, Same>::test(0, 0);
  Tv<blk_t, Spar, Stride_unknown, Same>::test(0, 0);
  Tv<blk_t, Full, Stride_unknown, Diff>::test(2, 2);
  Tv<blk_t, Cont, Stride_unknown, Diff>::test(2, 2);
  Tv<blk_t, Spar, Stride_unknown, Diff>::test(2, 2);


  // Asking for Stride_unit packing indicates we want unit stride,
  //   (but don't care about overall dense-ness or alignment)
  // - Subviews will have a compile-time cost of 2, but a runtime cost
  //   of 0 iff they have stride-1.
  // - Different complex format => copy access
  Tv<blk_t, Full, Stride_unit, Same>::test(0, 0);
  Tv<blk_t, Cont, Stride_unit, Same>::test(2, 0);	// runtime stride-1
  Tv<blk_t, Spar, Stride_unit, Same>::test(2, 2);
  Tv<blk_t, Full, Stride_unit, Diff>::test(2, 2);
  Tv<blk_t, Cont, Stride_unit, Diff>::test(2, 2);
  Tv<blk_t, Spar, Stride_unit, Diff>::test(2, 2);

  Tv<blk_t, Full, Stride_unit_dense, Same>::test(0, 0);
  Tv<blk_t, Cont, Stride_unit_dense, Same>::test(2, 0);
  Tv<blk_t, Off1, Stride_unit_dense, Same>::test(2, 0);
  Tv<blk_t, Spar, Stride_unit_dense, Same>::test(2, 2);
  Tv<blk_t, Full, Stride_unit_dense, Diff>::test(2, 2);
  Tv<blk_t, Cont, Stride_unit_dense, Diff>::test(2, 2);
  Tv<blk_t, Spar, Stride_unit_dense, Diff>::test(2, 2);

  Tv<blk_t, Full, Stride_unit_align<16>, Same>::test(2, 0);
  Tv<blk_t, Cont, Stride_unit_align<16>, Same>::test(2, 0);
  Tv<blk_t, Off1, Stride_unit_align<16>, Same>::test(2, 2);
  Tv<blk_t, Spar, Stride_unit_align<16>, Same>::test(2, 2);
}



void
matrix_tests()
{
  using vsip::impl::Stride_unknown;
  using vsip::impl::Stride_unit;
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Stride_unit_align;

  typedef Dense<2, float, row2_type> row_t;
  typedef Dense<2, float, col2_type> col_t;

  Tm<row_t, Full, Full, Same, Same, Same>::test(0, 0);
  Tm<col_t, Full, Full, Same, Same, Same>::test(0, 0);


  // Rationale: By asking for something with Stride_unknown, we should
  //   get the lowest cost at both compile-time and run-time.  Since
  //   Dense is direct, cost should be 0.

  Tm<row_t, Full, Full, Stride_unknown, Same, Same>::test(0, 0);
  Tm<row_t, Cont, Cont, Stride_unknown, Same, Same>::test(0, 0);
  Tm<row_t, Spar, Spar, Stride_unknown, Same, Same>::test(0, 0);


  // Rationale: Asking for stride_unknown but with a particular
  //   dimension-ordering indicates that dimension-ordering is
  //   most important, only copy to rearrange.

  Tm<row_t, Full, Full, Stride_unknown, row2_type, Same>::test(0, 0);
  Tm<row_t, Full, Full, Stride_unknown, col2_type, Same>::test(2, 2);


  // Rationale: Asking for stride_unit, we should get compile-time
  //   direct access for whole views and run-time direct access
  //   if the subview has unit-stride in the lowest order dimension.

  Tm<row_t, Full, Full, Stride_unit, Same, Same>::test(0, 0);
  Tm<row_t, Cont, Cont, Stride_unit, Same, Same>::test(2, 0);
  Tm<row_t, Cont, Spar, Stride_unit, Same, Same>::test(2, 2);
  Tm<row_t, Spar, Cont, Stride_unit, Same, Same>::test(2, 0);
  Tm<row_t, Spar, Spar, Stride_unit, Same, Same>::test(2, 2);

  Tm<col_t, Full, Full, Stride_unit, Same, Same>::test(0, 0);
  Tm<col_t, Cont, Cont, Stride_unit, Same, Same>::test(2, 0);
  Tm<col_t, Cont, Spar, Stride_unit, Same, Same>::test(2, 0);
  Tm<col_t, Spar, Cont, Stride_unit, Same, Same>::test(2, 2);
  Tm<col_t, Spar, Spar, Stride_unit, Same, Same>::test(2, 2);


  // Rationale: Asking for stride_unit but with a particular
  //   dimension-ordering indicates that we should copy if either
  //   the dimension-ordering is wrong, or the lowest-dimension
  //   is not unit stride.

  Tm<row_t, Full, Full, Stride_unit, row2_type, Same>::test(0, 0);
  Tm<row_t, Full, Cont, Stride_unit, row2_type, Same>::test(2, 0);
  Tm<row_t, Full, Spar, Stride_unit, row2_type, Same>::test(2, 2);
  Tm<row_t, Full, Full, Stride_unit, col2_type, Same>::test(2, 2);
  Tm<row_t, Cont, Full, Stride_unit, col2_type, Same>::test(2, 2);
  Tm<row_t, Spar, Full, Stride_unit, col2_type, Same>::test(2, 2);

  Tm<col_t, Full, Full, Stride_unit, row2_type, Same>::test(2, 2);
  Tm<col_t, Full, Cont, Stride_unit, row2_type, Same>::test(2, 2);
  Tm<col_t, Full, Spar, Stride_unit, row2_type, Same>::test(2, 2);
  Tm<col_t, Full, Full, Stride_unit, col2_type, Same>::test(0, 0);
  Tm<col_t, Cont, Full, Stride_unit, col2_type, Same>::test(2, 0);
  Tm<col_t, Spar, Full, Stride_unit, col2_type, Same>::test(2, 2);


  // Rationale: Asking for stride_unit_dense
  //   compile-time cost of 0 iff whole-view
  //   run-time cost of 0 if:
  //    - dimension order subdomains: Single* Cont? Full*, AND
  //    - order is same, complex is same

  // Blk  , Dim0, Dim1, Pack             , Ord , Cplx
  // -----,-----,-----,------------------,-----,-------------------
  Tm<row_t, Full, Full, Stride_unit_dense, Same, Same>::test(0, 0);
  Tm<row_t, Cont, Cont, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<row_t, Cont, Full, Stride_unit_dense, Same, Same>::test(2, 0);
  Tm<row_t, Full, Cont, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<row_t, Sing, Full, Stride_unit_dense, Same, Same>::test(2, 0);
  Tm<row_t, Sing, Cont, Stride_unit_dense, Same, Same>::test(2, 0);
  Tm<row_t, Full, Sing, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<row_t, Full, Sing, Stride_unit_dense, Same, Same>::test(2, 2);

  Tm<row_t, Full, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Cont, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Cont, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Full, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Sing, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Sing, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Full, Sing, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<row_t, Full, Sing, Stride_unit_dense, Diff, Same>::test(2, 2);

  Tm<col_t, Full, Full, Stride_unit_dense, Same, Same>::test(0, 0);
  Tm<col_t, Cont, Cont, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<col_t, Cont, Full, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<col_t, Full, Cont, Stride_unit_dense, Same, Same>::test(2, 0);
  Tm<col_t, Sing, Full, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<col_t, Sing, Cont, Stride_unit_dense, Same, Same>::test(2, 2);
  Tm<col_t, Full, Sing, Stride_unit_dense, Same, Same>::test(2, 0);
  Tm<col_t, Full, Sing, Stride_unit_dense, Same, Same>::test(2, 0);

  Tm<col_t, Full, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Cont, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Cont, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Full, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Sing, Full, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Sing, Cont, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Full, Sing, Stride_unit_dense, Diff, Same>::test(2, 2);
  Tm<col_t, Full, Sing, Stride_unit_dense, Diff, Same>::test(2, 2);


  // Rationale: Dense blocks are not currently Stride_unit_aligned at
  // compile-time.

  // Row-major is not aligned.
  test_assert(g_cols * sizeof(float) % 16 != 0);
  Tm<row_t, Full, Full, Stride_unit_align<16>, Same, Same>::test(2, 2);

  // Every fourth column is aligned, for row-major:
  test_assert(4 * g_cols * sizeof(float) % 16 == 0);
  Tm<row_t, Spar4, Full, Stride_unit_align<16>, Same, Same>::test(2, 0);

  // Column-major is not aligned.
  test_assert(g_rows * sizeof(float) % 16 != 0);
  Tm<col_t, Full, Full, Stride_unit_align<16>, Same, Same>::test(2, 2);

  // Every other row is aligned, for column-major:
  test_assert(2 * g_rows * sizeof(float) % 16 == 0);
  Tm<col_t, Full, Spar, Stride_unit_align<16>, Same, Same>::test(2, 0);
}



void
tensor_tests()
{
  using vsip::impl::Stride_unknown;
  using vsip::impl::Stride_unit;
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Stride_unit_align;

  typedef Dense<3, float, row3_type> row_t;
  typedef Dense<3, float, col3_type> col_t;

  Tt<row_t, Full, Full, Full, Stride_unknown, Same, Same>::test(0, 0);
  Tt<row_t, Full, Spar, Cont, Stride_unknown, Same, Same>::test(0, 0);
  Tt<row_t, Spar, Spar, Spar, Stride_unknown, Same, Same>::test(0, 0);
  

  Tt<row_t, Full, Full, Full, Stride_unit, Same, Same>::test(0, 0);
  Tt<row_t, Spar, Spar, Full, Stride_unit, Same, Same>::test(2, 0);
  Tt<row_t, Spar, Spar, Full, Stride_unit, Same, Same>::test(2, 0);
  Tt<row_t, Spar, Spar, Cont, Stride_unit, Same, Same>::test(2, 0);
  Tt<row_t, Spar, Spar, Spar, Stride_unit, Same, Same>::test(2, 2);

  Tt<col_t, Full, Full, Full, Stride_unit, Same, Same>::test(0, 0);
  Tt<col_t, Full, Spar, Spar, Stride_unit, Same, Same>::test(2, 0);
  Tt<col_t, Full, Spar, Spar, Stride_unit, Same, Same>::test(2, 0);
  Tt<col_t, Cont, Spar, Spar, Stride_unit, Same, Same>::test(2, 0);
  Tt<col_t, Spar, Spar, Spar, Stride_unit, Same, Same>::test(2, 2);


  Tt<row_t, Full, Full, Full, Stride_unit_dense, Same, Same>::test(0, 0);
  Tt<row_t, Cont, Full, Full, Stride_unit_dense, Same, Same>::test(2, 0);
  Tt<row_t, Spar, Full, Full, Stride_unit_dense, Same, Same>::test(2, 2);
  Tt<row_t, Cont, Cont, Full, Stride_unit_dense, Same, Same>::test(2, 2);
  Tt<row_t, Sing, Cont, Full, Stride_unit_dense, Same, Same>::test(2, 0);

  Tt<col_t, Full, Full, Full, Stride_unit_dense, Same, Same>::test(0, 0);
  Tt<col_t, Full, Full, Cont, Stride_unit_dense, Same, Same>::test(2, 0);
  Tt<col_t, Full, Full, Spar, Stride_unit_dense, Same, Same>::test(2, 2);
  Tt<col_t, Full, Cont, Cont, Stride_unit_dense, Same, Same>::test(2, 2);
  Tt<col_t, Full, Cont, Sing, Stride_unit_dense, Same, Same>::test(2, 0);

  test_assert(g_dim2 * sizeof(float) % 16 != 0); // row_t is not aligned
  test_assert(g_dim0 * sizeof(float) % 16 == 0); // col_t is     aligned

  Tt<row_t, Full, Full, Full, Stride_unit_align<16>, Same, Same>::test(2, 2);
  Tt<col_t, Full, Full, Full, Stride_unit_align<16>, Same, Same>::test(2, 0);
  Tt<col_t, Full, Cont, Cont, Stride_unit_align<16>, Same, Same>::test(2, 0);
  Tt<col_t, Full, Spar, Spar, Stride_unit_align<16>, Same, Same>::test(2, 0);
}



int
main(int argc, char** argv)
{
  vsip::vsipl init(argc, argv);

  vector_tests();
  matrix_tests();
  tensor_tests();
}

