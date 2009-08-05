/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    extdata_dist.cpp
    @author  Jules Bergmann
    @date    2006-01-31
    @brief   VSIPL++ Library: Unit tests for Ext_data_dist.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/vector.hpp>
#include <vsip/core/extdata_dist.hpp>
#include <vsip/initfin.hpp>
#include <vsip/map.hpp>
#include <vsip/selgen.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using namespace vsip_csl;

using vsip::impl::Layout;
using vsip::impl::Stride_unit_dense;
using vsip::impl::Cmplx_inter_fmt;
using vsip::impl::Cmplx_split_fmt;
using vsip::impl::Any_type;
using vsip::impl::Fast_block;
using vsip::impl::Block_layout;
using vsip::impl::Adjust_layout;
using vsip::impl::Type_equal;
using vsip::impl::SYNC_IN;
using vsip::impl::SYNC_INOUT;
using vsip::impl::ITE_Type;
using vsip::impl::As_type;
using vsip::impl::Scalar_of;
using vsip::impl::Aligned_allocator;
using vsip::impl::Is_local_map;
using vsip::impl::Is_complex;



/***********************************************************************
  Definitions
***********************************************************************/

/// Test high-level data interface to 1-dimensional block.

template <typename T,
	  typename MapT,
	  typename GivenLP,
	  typename ReqLP>
void
test_1_ext(int expect_cost)
{
  length_type const size = 10;

  typedef Fast_block<1, T, GivenLP, MapT> block_type;

  typedef T value_type;
  typedef typename Adjust_layout<value_type, ReqLP, GivenLP>::type use_LP;
  typedef vsip::impl::Storage< typename use_LP::complex_type, T> storage_type;

  block_type block(size, 0.0);

  value_type val0 =  1.0f;
  value_type val1 =  2.78f;
  value_type val2 =  3.14f;
  value_type val3 = -1.5f;

  // Place values in block.
  block.put(0, val0);
  block.put(1, val1);

  typedef typename
      ITE_Type<Type_equal<typename use_LP::complex_type,
                          Cmplx_inter_fmt>::value,
               As_type<Aligned_allocator<T> >,
               As_type<Aligned_allocator<typename Scalar_of<T>::type> > >
      ::type alloc_type;
  alloc_type alloc;
  typename storage_type::type buffer = storage_type::allocate(alloc, size);

  {
    vsip::impl::Ext_data_dist<block_type, SYNC_INOUT, use_LP> raw(block, buffer);

    assert(raw.cost() == expect_cost);

    // Check properties of DDI.
    test_assert(raw.stride(0) == 1);
    test_assert(raw.size(0) == size);

    typename storage_type::type data = raw.data();

    // Check that block values are reflected.
    test_assert(equal(storage_type::get(data, 0), val0));
    test_assert(equal(storage_type::get(data, 1), val1));

    // Place values in raw data.
    storage_type::put(data, 1, val2);
    storage_type::put(data, 2, val3);
  }

  storage_type::deallocate(alloc, buffer, size);

  // Check that raw data values are reflected.
  test_assert(equal(block.get(1), val2));
  test_assert(equal(block.get(2), val3));
}



// Determine expected cost for Ext_data_dist.

template <typename T,
	  typename MapT,
	  typename OrderT,
	  typename use_LP,
	  typename GivenLP>
int expected_cost()
{
  bool is_local_equiv = Is_local_map<MapT>::value ||
                        Type_equal<MapT, Global_map<1> >::value;

  bool same_order = Type_equal<OrderT, typename use_LP::order_type>::value;

  bool same_complex_fmt = 
	(!Is_complex<T>::value ||
	 Type_equal<typename GivenLP::complex_type,
	            typename use_LP::complex_type>::value);

  return (is_local_equiv && same_order && same_complex_fmt) ? 0 : 2;
}



template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP>
void
test_1_dense(MapT const& map)
{
  length_type const size = 10;

  typedef Dense<1, T, OrderT, MapT> block_type;

  typedef typename Block_layout<block_type>::layout_type GivenLP;

  typedef T value_type;
  typedef typename Adjust_layout<value_type, ReqLP, GivenLP>::type use_LP;
  typedef vsip::impl::Storage< typename use_LP::complex_type, T> storage_type;

  block_type block(size, T(), map);

  value_type val0 =  1.0f;
  value_type val1 =  2.78f;
  value_type val2 =  3.14f;
  value_type val3 = -1.5f;

  // Place values in block.
  block.put(0, val0);
  block.put(1, val1);

  typedef typename
      ITE_Type<Type_equal<typename use_LP::complex_type,
                          Cmplx_inter_fmt>::value,
               As_type<Aligned_allocator<T> >,
               As_type<Aligned_allocator<typename Scalar_of<T>::type> > >
      ::type alloc_type;

  alloc_type alloc;
  typename storage_type::type buffer = storage_type::allocate(alloc, size);

  {
    vsip::impl::Ext_data_dist<block_type, SYNC_INOUT, use_LP> raw(block, buffer);

    int exp_cost = expected_cost<T, MapT, OrderT, use_LP, GivenLP>();
    assert(raw.cost() == exp_cost);

    // Check properties of DDI.
    test_assert(raw.stride(0) == 1);
    test_assert(raw.size(0) == size);

    typename storage_type::type data = raw.data();

    // Check that block values are reflected.
    test_assert(equal(storage_type::get(data, 0), val0));
    test_assert(equal(storage_type::get(data, 1), val1));

    // Place values in raw data.
    storage_type::put(data, 1, val2);
    storage_type::put(data, 2, val3);
  }

  storage_type::deallocate(alloc, buffer, size);

  // Check that raw data values are reflected.
  test_assert(equal(block.get(1), val2));
  test_assert(equal(block.get(2), val3));
}



template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP>
void
test_1_dense_const(MapT const& map)
{
  length_type const size = 10;

  typedef Dense<1, T, OrderT, MapT> block_type;

  typedef typename Block_layout<block_type>::layout_type GivenLP;

  typedef T value_type;
  typedef typename Adjust_layout<value_type, ReqLP, GivenLP>::type use_LP;
  typedef vsip::impl::Storage< typename use_LP::complex_type, T> storage_type;

  block_type block(size, T(), map);

  value_type val0 =  1.0f;
  value_type val1 =  2.78f;

  // Place values in block.
  block.put(0, val0);
  block.put(1, val1);

  typedef typename
      ITE_Type<Type_equal<typename use_LP::complex_type,
                          Cmplx_inter_fmt>::value,
               As_type<Aligned_allocator<T> >,
               As_type<Aligned_allocator<typename Scalar_of<T>::type> > >
      ::type alloc_type;

  alloc_type alloc;
  typename storage_type::type buffer = storage_type::allocate(alloc, size);

  {
    block_type const& ref = block;
    vsip::impl::Ext_data_dist<block_type, SYNC_IN, use_LP> raw(ref, buffer);

    int exp_cost = expected_cost<T, MapT, OrderT, use_LP, GivenLP>();
    assert(raw.cost() == exp_cost);

    // Check properties of DDI.
    test_assert(raw.stride(0) == 1);
    test_assert(raw.size(0) == size);

    typename storage_type::type data = raw.data();

    // Check that block values are reflected.
    test_assert(equal(storage_type::get(data, 0), val0));
    test_assert(equal(storage_type::get(data, 1), val1));
  }

  storage_type::deallocate(alloc, buffer, size);
}



// Helper function to test Ext_data_dist access to an expression block.
// Deduces type of expression block.
//

template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP,
	  typename ViewT>
void
test_1_expr_helper(ViewT view, T val0, T val1, bool alloc_buffer)
{
  typedef typename ViewT::block_type block_type;
  typedef typename ViewT::value_type value_type;

  typedef typename Block_layout<block_type>::layout_type GivenLP;

  typedef typename Adjust_layout<value_type, ReqLP, GivenLP>::type use_LP;
  typedef vsip::impl::Storage< typename use_LP::complex_type, T> storage_type;

  typedef typename
      ITE_Type<Type_equal<typename use_LP::complex_type,
                          Cmplx_inter_fmt>::value,
               As_type<Aligned_allocator<T> >,
               As_type<Aligned_allocator<typename Scalar_of<T>::type> > >
      ::type alloc_type;

  length_type size = view.size();

  alloc_type alloc;
  typedef typename storage_type::type ptr_type;
  ptr_type buffer = alloc_buffer ? storage_type::allocate(alloc, size)
                                 : ptr_type();

  {
    vsip::impl::Ext_data_dist<block_type, SYNC_IN, use_LP> raw(view.block(), buffer);

    // Because block is an expression block, access requires a copy.
    assert(raw.cost() == 2);

    // Check properties of DDI.
    test_assert(raw.stride(0) == 1);
    test_assert(raw.size(0) == size);

    typename storage_type::type data = raw.data();

    // Check that block values are reflected.
    test_assert(equal(storage_type::get(data, 0), val0));
    test_assert(equal(storage_type::get(data, 1), val1));
  }

  if (alloc_buffer) storage_type::deallocate(alloc, buffer, size);
}



// Test Ext_data_dist access to a simple expression block.
// (vector + vector).

template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP>
void
test_1_expr_1(MapT const& map)
{
  typedef Dense<1, T, OrderT, MapT> block_type;

  length_type const size = 10;

  T val0 =  1.0f;
  T val1 =  2.78f;

  Vector<T, block_type> view1(size, T(), map);
  Vector<T, block_type> view2(size, T(), map);

  // Place values in block.
  view1.put(0, val0);
  view2.put(1, val1);

  test_1_expr_helper<T, OrderT, MapT, ReqLP>(view1 + view2, val0, val1,
					     true);
  test_1_expr_helper<T, OrderT, MapT, ReqLP>(view1 + view2, val0, val1,
					     false);
}



// Test Ext_data_dist access to a more complex expression block.
// (vector + vector)(subset).

template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP>
void
test_1_expr_2(MapT const& map)
{
  typedef Dense<1, T, OrderT, MapT> block_type;

  length_type const size = 10;

  T val0 =  1.0f;
  T val1 =  2.78f;

  Vector<T, block_type> view1(size, T(), map);
  Vector<T, block_type> view2(size, T(), map);

  // Place values in block.
  view1.put(0, val0);
  view2.put(1, val1);

  test_1_expr_helper<T, OrderT, MapT, ReqLP>((view1 + view2)(Domain<1>(8)),
					     val0, val1, true);
  test_1_expr_helper<T, OrderT, MapT, ReqLP>((view1 + view2)(Domain<1>(8)),
					     val0, val1, false);
}



// Test Ext_data_dist access to an 'magsq(v)' expression block.

template <typename T,
	  typename OrderT,
	  typename MapT,
	  typename ReqLP>
void
test_1_expr_3(MapT const& map)
{
  typedef Dense<1, T, OrderT, MapT> block_type;

  length_type const size = 10;

  T val0 =  1.0f;
  T val1 =  2.78f;

  Vector<T, block_type> view(size, T(), map);

  // Place values in block.
  view.put(0, val0);
  view.put(1, val1);

  typedef typename vsip::impl::Scalar_of<T>::type scalar_type;

  test_1_expr_helper<scalar_type, OrderT, MapT, ReqLP>(
		magsq(view), magsq(val0), magsq(val1), true);
  test_1_expr_helper<scalar_type, OrderT, MapT, ReqLP>(
		magsq(view), magsq(val0), magsq(val1), false);
}



void
test()
{
  typedef Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP_1rdi;
  typedef Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP_1rds;

  typedef Layout<1, Any_type, Any_type, Cmplx_inter_fmt> LP_1xxi;
  typedef Layout<1, Any_type, Any_type, Cmplx_split_fmt> LP_1xxs;

  test_1_ext<float,          Local_map, LP_1rdi, LP_1xxi >(0);
  test_1_ext<float,          Local_map, LP_1rds, LP_1xxi >(0);
  test_1_ext<float,          Local_map, LP_1rdi, LP_1xxs >(0);
  test_1_ext<float,          Local_map, LP_1rds, LP_1xxs >(0);

  test_1_ext<complex<float>, Local_map, LP_1rdi, LP_1xxi >(0);
  test_1_ext<complex<float>, Local_map, LP_1rds, LP_1xxi >(2);
  test_1_ext<complex<float>, Local_map, LP_1rdi, LP_1xxs >(2);
  test_1_ext<complex<float>, Local_map, LP_1rds, LP_1xxs >(0);

  Local_map lmap;

  test_1_expr_1<float,          row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_1<complex<float>, row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_1<complex<float>, row1_type, Local_map, LP_1xxs >(lmap);

  test_1_expr_2<float,          row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_2<complex<float>, row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_2<complex<float>, row1_type, Local_map, LP_1xxs >(lmap);

  test_1_expr_3<float,          row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_3<complex<float>, row1_type, Local_map, LP_1xxi >(lmap);
  test_1_expr_3<complex<float>, row1_type, Local_map, LP_1xxs >(lmap);


  Global_map<1> gmap;

  test_1_dense<float,          row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_dense<float,          row1_type, Global_map<1>, LP_1xxs >(gmap);
  test_1_dense<complex<float>, row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_dense<complex<float>, row1_type, Global_map<1>, LP_1xxs >(gmap);

  test_1_dense_const<float,          row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_dense_const<float,          row1_type, Global_map<1>, LP_1xxs >(gmap);
  test_1_dense_const<complex<float>, row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_dense_const<complex<float>, row1_type, Global_map<1>, LP_1xxs >(gmap);

  test_1_expr_1<float,          row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_1<complex<float>, row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_1<complex<float>, row1_type, Global_map<1>, LP_1xxs >(gmap);

  test_1_expr_2<float,          row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_2<complex<float>, row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_2<complex<float>, row1_type, Global_map<1>, LP_1xxs >(gmap);

  test_1_expr_3<float,          row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_3<complex<float>, row1_type, Global_map<1>, LP_1xxi >(gmap);
  test_1_expr_3<complex<float>, row1_type, Global_map<1>, LP_1xxs >(gmap);


  Map<Block_dist> map(num_processors());

  test_1_dense<float,          row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_dense<float,          row1_type, Map<Block_dist>, LP_1xxs >(map);
  test_1_dense<complex<float>, row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_dense<complex<float>, row1_type, Map<Block_dist>, LP_1xxs >(map);

  test_1_dense_const<float,          row1_type, Map<Block_dist>, LP_1xxi>(map);
  test_1_dense_const<float,          row1_type, Map<Block_dist>, LP_1xxs>(map);
  test_1_dense_const<complex<float>, row1_type, Map<Block_dist>, LP_1xxi>(map);
  test_1_dense_const<complex<float>, row1_type, Map<Block_dist>, LP_1xxs>(map);

  test_1_expr_1<float,          row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_1<complex<float>, row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_1<complex<float>, row1_type, Map<Block_dist>, LP_1xxs >(map);

  test_1_expr_2<float,          row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_2<complex<float>, row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_2<complex<float>, row1_type, Map<Block_dist>, LP_1xxs >(map);

  test_1_expr_3<float,          row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_3<complex<float>, row1_type, Map<Block_dist>, LP_1xxi >(map);
  test_1_expr_3<complex<float>, row1_type, Map<Block_dist>, LP_1xxs >(map);
}



int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test();
}
