/* Copyright (c) 2005, 2006, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/extdata-matadd.cpp
    @author  Jules Bergmann
    @date    2005-02-14
    @brief   VSIPL++ Library: Element-wise matrix add exmples/tests for DDI.

    This file illustrates how data access may be used to perform
    element-wise matrix add, in three different ways:

     - A simple function for addition of matrices row by row, using
       row/col stride information, with no attempt to re-arrange data
       or select appropriate traversal based on layout. (matrix_add_1)

     - A simple function for addition of matrices using row/col stride
       information (relying on compiler strength reduction), with no
       attempt to re-arrange data or select appropriate traversal
       based on layout. (matrix_add_2)

     - A more complex class template and gateway function design to
       add matrices with an "appropriate" algorithm.  Algorithm is
       selected based on dimension ordering (selection of row-major
       traversal vs col-major traversal), and data layout (a single
       loop can be used for dense/continguous data). (matrix_add).
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/dense.hpp>
#include <vsip/matrix.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/plainblock.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions
***********************************************************************/


/// Element-wise matrix addition, example 1.

/// Matrices are traversed row-by-row, with preference for row-major
/// data, using row/column stride information.

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
void
matrix_add_1(
  const_Matrix<TR, BlockR> res,
  const_Matrix<T1, Block1> op1,
  const_Matrix<T2, Block2> op2)
{
  vsip::impl::Ext_data<BlockR> raw_res(res.block(), vsip::impl::SYNC_OUT);
  float*   p_raw       = raw_res.data();
  stride_type row_str_raw = raw_res.stride(0);
  stride_type col_str_raw = raw_res.stride(1);

  vsip::impl::Ext_data<Block1> raw1(op1.block(), vsip::impl::SYNC_IN);
  float*   p1       = raw1.data();
  stride_type row_str1 = raw1.stride(0);
  stride_type col_str1 = raw1.stride(1);

  vsip::impl::Ext_data<Block2> raw2(op2.block(), vsip::impl::SYNC_IN);
  float*   p2       = raw2.data();
  stride_type row_str2 = raw2.stride(0);
  stride_type col_str2 = raw2.stride(1);

  for (index_type r=0; r<res.size(0); ++r)
  {
    float* row_raw = p_raw;
    float* row_1   = p1;
    float* row_2   = p2;

    for (index_type c=0; c<res.size(1); ++c)
    {
      *row_raw = *row_1 + *row_2;

      row_1   += col_str1;
      row_2   += col_str2;
      row_raw += col_str_raw;
    }
    p_raw += row_str_raw;
    p1    += row_str1;
    p2    += row_str2;
  }
}



/// Element-wise matrix addition, example 2.

/// Matrices are traversed by computing index for each element
/// directly from row/col strides.  Preference for row-major data.

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
void
matrix_add_2(
  const_Matrix<TR, BlockR> res,
  const_Matrix<T1, Block1> op1,
  const_Matrix<T2, Block2> op2)
{
  vsip::impl::Ext_data<BlockR> raw_res(res.block(), vsip::impl::SYNC_OUT);
  vsip::impl::Ext_data<Block1> raw1(op1.block(), vsip::impl::SYNC_IN);
  vsip::impl::Ext_data<Block2> raw2(op2.block(), vsip::impl::SYNC_IN);

  float*   pR = raw_res.data();
  float*   p1 = raw1.data();
  float*   p2 = raw2.data();

  for (index_type r=0; r<res.size(0); ++r)
  {
    for (index_type c=0; c<res.size(1); ++c)
    {
      pR[r*raw_res.stride(0) + c*raw_res.stride(1)] =
	p1[r*raw1.stride(0) + c*raw1.stride(1)] +
	p2[r*raw2.stride(0) + c*raw2.stride(1)];
    }
  }
}



/***********************************************************************
  Definitions for extended example.
***********************************************************************/

/// Tag class to indicate preferred implementation.
///
/// Tag_plain<Order> indicates a plain two-loop implementation, with
///   traversal order matching ORDER.
///
/// Tag_contig indicates a one-loop implementation, possible when
///   matrices are dense/contiguous and has matching dimension ordering.

template <typename Order> struct Tag_plain  {};
struct Tag_contig {};



/// MASelect is a helper class to select the appropriate matrix-add
/// implementation tag based on the layout policies of the result and
/// operand matrices.

/// Requires:
///   LAYOUTPOLICYR to be the layout policy of the result matrix.
///   LAYOUTPOLICY1 to be the layout policy of the first operand.
///   LAYOUTPOLICY2 to be the layout policy of the second operand.
///
/// Provides:
///   TYPE, a tag selecting the preferred implementation.



/// General Case: select Tag_plain, using the dimension ordering of the
/// result to select order of traversal.

template <typename OrderR,
	  typename Order1,
	  typename Order2,
	  typename PackR,
	  typename Pack1,
	  typename Pack2>
struct MASelect
{
  typedef Tag_plain<OrderR> type;
};



/// Specialization: select Tag_contig when all matrices are dense and
/// have same dimension ordering.

template <typename Order>
struct MASelect<vsip::impl::Stride_unit_dense,
                vsip::impl::Stride_unit_dense,
		vsip::impl::Stride_unit_dense,
		Order, Order, Order>
{
  typedef Tag_contig type;
};



/// Specializations: select Tag_contig when all matrices are dense and
/// one operand has same dimension ordering as result.  This will require
/// one of the operands to be rearranged.

template <typename Order1,
	  typename Order2>
struct MASelect<vsip::impl::Stride_unit_dense,
                vsip::impl::Stride_unit_dense,
		vsip::impl::Stride_unit_dense,
		Order1, Order1, Order2>
{
  typedef Tag_contig type;
};

template <typename Order1,
	  typename Order2>
struct MASelect<vsip::impl::Stride_unit_dense,
                vsip::impl::Stride_unit_dense,
		vsip::impl::Stride_unit_dense,
		Order1, Order2, Order1>
{
  typedef Tag_contig type;
};



// Matrix_add is a class template used to dispatch and implement
// element-wise matrix-addition.
//
// Requires:
//   ... usual template parameters
//   TAG is a tag class that indicates which implementation should be
//      used.  Valid choices include Tag_plain and Tag_contig.
//      Default is to use the tag suggested by MASelect.

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2,
	  typename Tag = Tag_plain<typename BlockR::order_type> >
struct Matrix_add;



// Matrix_add specialization for Tag_plain<row2_type>.
//
// Element-wise matrix add using row/column strides to compute index,
// with traversal in row-major order (no attempt to rearrange data).

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
struct Matrix_add<TR, T1, T2, BlockR, Block1, Block2, Tag_plain<row2_type> >
{
  static void add(
    const_Matrix<TR, BlockR> res,
    const_Matrix<T1, Block1> op1,
    const_Matrix<T2, Block2> op2)
  {
    vsip::impl::Ext_data<BlockR> raw_res(res.block(), vsip::impl::SYNC_OUT);
    vsip::impl::Ext_data<Block1> raw1(op1.block(), vsip::impl::SYNC_IN);
    vsip::impl::Ext_data<Block2> raw2(op2.block(), vsip::impl::SYNC_IN);

    // int cost = raw_res.cost + raw1.cost + raw2.cost;
    // cout << "Tag_plain " << cost << endl;

    float*   pR = raw_res.data();
    float*   p1 = raw1.data();
    float*   p2 = raw2.data();

    for (index_type r=0; r<res.size(0); ++r)
    {
      for (index_type c=0; c<res.size(1); ++c)
      {
	pR[r*raw_res.stride(0) + c*raw_res.stride(1)] =
	  p1[r*raw1.stride(0) + c*raw1.stride(1)] +
	  p2[r*raw2.stride(0) + c*raw2.stride(1)];
      }
    }
  }
};



// Matrix_add specialization for Tag_plain<col2_type>.
//
// Element-wise matrix add using row/column strides to compute index,
// with traversal in col-major order (no attempt to rearrange data).

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
struct Matrix_add<TR, T1, T2, BlockR, Block1, Block2, Tag_plain<col2_type> >
{
  static void add(
    const_Matrix<TR, BlockR> res,
    const_Matrix<T1, Block1> op1,
    const_Matrix<T2, Block2> op2)
  {
    vsip::impl::Ext_data<BlockR> raw_res(res.block(), vsip::impl::SYNC_OUT);
    vsip::impl::Ext_data<Block1> raw1(op1.block(), vsip::impl::SYNC_IN);
    vsip::impl::Ext_data<Block2> raw2(op2.block(), vsip::impl::SYNC_IN);

    // int cost = raw_res.cost + raw1.cost + raw2.cost;
    // cout << "Tag_plain " << cost << endl;

    float*   pR = raw_res.data();
    float*   p1 = raw1.data();
    float*   p2 = raw2.data();

    for (index_type c=0; c<res.size(1); ++c)
    {
      for (index_type r=0; r<res.size(0); ++r)
      {
	pR[r*raw_res.stride(0) + c*raw_res.stride(1)] =
	  p1[r*raw1.stride(0) + c*raw1.stride(1)] +
	  p2[r*raw2.stride(0) + c*raw2.stride(1)];
      }
    }
  }
};



// Matrix_add specialization for Tag_contig.
//
// Element-wise matrix add for contiguous data with consistent
// dimension layout.  Dimension layout is choosen by result matrix.
// (Data may be rearranged).

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
struct Matrix_add<TR, T1, T2, BlockR, Block1, Block2, Tag_contig>
{
  static void add(
    const_Matrix<TR, BlockR> res,
    const_Matrix<T1, Block1> op1,
    const_Matrix<T2, Block2> op2)
  {
    typedef typename BlockR::layout_type layout_type;
    typedef vsip::impl::No_count_policy RP;

    // Check that no memory is required.
    // test_assert((impl::Ext_data<BlockR, layout_type>::CT_Mem_not_req));
    // test_assert((impl::Ext_data<Block1, layout_type>::CT_Mem_not_req));
    // test_assert((impl::Ext_data<Block2, layout_type>::CT_Mem_not_req));

    vsip::impl::Ext_data<BlockR, layout_type, RP> raw_res(res.block(), vsip::impl::SYNC_OUT);
    vsip::impl::Ext_data<Block1, layout_type, RP> raw1(op1.block(), vsip::impl::SYNC_IN);
    vsip::impl::Ext_data<Block2, layout_type, RP> raw2(op2.block(), vsip::impl::SYNC_IN);

    // int cost = raw_res.cost + raw1.cost + raw2.cost;
    // cout << "Tag_contig " << cost << endl;

    float*   pR = raw_res.data();
    float*   p1 = raw1.data();
    float*   p2 = raw2.data();

    for (index_type i=0; i<res.size(); ++i)
    {
      *pR = *p1 + *p2;
      ++pR;
      ++p1;
      ++p2;
    }
  }
};



/// Gateway function for matrix_add.

/// Uses Matrix_add class to dispatches to appropriate implementation.

template <typename TR,
	  typename T1,
	  typename T2,
	  typename BlockR,
	  typename Block1,
	  typename Block2>
void
matrix_add(
  const_Matrix<TR, BlockR> res,
  const_Matrix<T1, Block1> op1,
  const_Matrix<T2, Block2> op2)
{
  Matrix_add<TR, T1, T2, BlockR, Block1, Block2>::add(res, op1, op2);
}



template <typename OrderR,
	  typename Order1,
	  typename Order2>
void
test_matrix_add()
{
  length_type num_rows = 10;
  length_type num_cols = 15;

  Matrix<float, Dense<2, float, Order1> > view1(num_rows, num_cols);
  Matrix<float, Dense<2, float, Order2> > view2(num_rows, num_cols);
  Matrix<float, Dense<2, float, OrderR> > res(num_rows, num_cols);
  
  // Place data into operands.
  for (index_type r=0; r<num_rows; ++r)
    for (index_type c=0; c<num_cols; ++c)
    {
      view1.put(r, c, float(1*(r*num_cols+c)));
      view2.put(r, c, float(2*(r*num_cols+c)));
    }

  matrix_add_1(res, view1, view2);

  for (index_type r=0; r<num_rows; ++r)
    for (index_type c=0; c<num_cols; ++c)
      test_assert(equal(res.get(r, c), view1.get(r, c) + view2.get(r, c)));


  matrix_add_2(res, view1, view2);

  for (index_type r=0; r<num_rows; ++r)
    for (index_type c=0; c<num_cols; ++c)
      test_assert(equal(res.get(r, c), view1.get(r, c) + view2.get(r, c)));


  matrix_add(res, view1, view2);

  for (index_type r=0; r<num_rows; ++r)
    for (index_type c=0; c<num_cols; ++c)
      test_assert(equal(res.get(r, c), view1.get(r, c) + view2.get(r, c)));
}



int
main(int argc, char** argv)
{
  vsip::vsipl init(argc, argv);
  test_matrix_add<row2_type, row2_type, row2_type>();
  test_matrix_add<row2_type, col2_type, row2_type>();
  test_matrix_add<row2_type, row2_type, col2_type>();
  test_matrix_add<col2_type, row2_type, col2_type>();
  test_matrix_add<col2_type, row2_type, row2_type>();
}
