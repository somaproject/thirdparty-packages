/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/rt_extdata.cpp
    @author  Jules Bergmann
    @date    2006-05-03
    @brief   VSIPL++ Library: Unit-test for run-time external data access.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/map.hpp>
#include <vsip/opt/rt_extdata.hpp>

#include <vsip_csl/test.hpp>
#include "util.hpp"

using namespace vsip;
using namespace vsip_csl;

using vsip::impl::Rt_layout;
using vsip::impl::Rt_tuple;
using vsip::impl::rt_pack_type;
using vsip::impl::rt_complex_type;
using vsip::impl::Storage;
using vsip::impl::Cmplx_inter_fmt;
using vsip::impl::Cmplx_split_fmt;
using vsip::impl::Rt_ext_data;
using vsip::impl::Applied_layout;
using vsip::impl::Length;
using vsip::impl::sync_action_type;

using vsip::impl::stride_unit_dense;
using vsip::impl::stride_unit_align;
using vsip::impl::cmplx_inter_fmt;
using vsip::impl::cmplx_split_fmt;
using vsip::impl::SYNC_INOUT;
using vsip::impl::SYNC_IN_NOPRESERVE;

using vsip::impl::extent;



/***********************************************************************
  Definitions
***********************************************************************/

// Utility functions to return a unique value for each index in
// a view.  Overloaded so that tests can work for any-dimension view.

inline index_type value1(Index<1> const& idx) { return idx[0]; }
inline index_type value1(Index<2> const& idx) { return 100*idx[0] + idx[1]; }
inline index_type value1(Index<3> const& idx)
{ return 10000*idx[0] + 100*idx[1] + idx[2]; }

inline index_type value2(Index<1> const& idx) { return 2*idx[0]; }
inline index_type value2(Index<2> const& idx) { return 100*idx[1] + idx[0]; }
inline index_type value2(Index<3> const& idx)
{ return 10000*idx[2] + 100*idx[0] + idx[1]; }


template <dimension_type Dim,
          typename       ExtDataT>
inline stride_type
offset(Index<Dim> const& idx,
       ExtDataT const&   ext)
{
  stride_type off = stride_type();

  for (dimension_type d=0; d<Dim; ++d)
    off += idx[d] * ext.stride(d);

  return off;
}



// Test that Rt_layout matches Layout.
template <typename       LayoutT,
	  dimension_type Dim>
void
test_layout(Rt_layout<Dim> rtl)
{
  using vsip::impl::pack_format;
  using vsip::impl::complex_format;
  using vsip::impl::Is_stride_unit_align;

  typedef typename LayoutT::pack_type pack_type;
  typedef typename LayoutT::complex_type complex_type;
  typedef typename LayoutT::order_type order_type;
  typedef typename LayoutT::pack_type pack_type;

  test_assert(rtl.dim     == LayoutT::dim);
  test_assert(rtl.pack    == pack_format<pack_type>());
  test_assert(rtl.complex == complex_format<complex_type>());
  test_assert(rtl.order.impl_dim0 == LayoutT::order_type::impl_dim0);
  test_assert(rtl.order.impl_dim1 == LayoutT::order_type::impl_dim1);
  test_assert(rtl.order.impl_dim2 == LayoutT::order_type::impl_dim2);
  test_assert(rtl.align == Is_stride_unit_align<pack_type>::align);
}



// Test run-time external data access (assuming that data is either
// not complex or is interleaved-complex).

template <typename       T,
	  typename       LayoutT,
	  dimension_type Dim>
void
t_rtex(
  Domain<Dim> const& dom,
  Rt_tuple           order,
  rt_pack_type       pack,
  rt_complex_type    cformat,
  bool               alloc,
  sync_action_type   sync,
  int                cost)
{
  typedef impl::Fast_block<Dim, T, LayoutT> block_type;
  typedef typename impl::View_of_dim<Dim, T, block_type>::type view_type;

  view_type view = create_view<view_type>(dom);

  Rt_layout<Dim> blk_rtl = vsip::impl::block_layout<Dim>(view.block());
  test_layout<LayoutT>(blk_rtl);

  Length<Dim> len = impl::extent(dom);
  for (Index<Dim> idx; valid(len, idx); next(len, idx))
    put(view, idx, T(value1(idx)));

  Rt_layout<Dim> rt_layout;

  rt_layout.pack    = pack;
  rt_layout.order   = order; 
  rt_layout.complex = cformat;
  rt_layout.align   = (pack == stride_unit_align) ? 16 : 0;

  // Pre-allocate temporary buffer.
  T* buffer = 0;
  if (alloc)
  {
    Applied_layout<Rt_layout<Dim> > app_layout(rt_layout, len, sizeof(T));
    length_type total_size = app_layout.total_size();
    buffer = new T[total_size];
  }

  {
    Rt_ext_data<block_type> ext(view.block(), rt_layout, SYNC_INOUT, buffer);

    T* ptr = ext.data().as_inter();

    test_assert(cost == ext.cost());
    if (alloc && cost != 0)
      test_assert(ptr == buffer);

#if VERBOSE
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "cost: " << ext.cost() << std::endl;

    for (index_type i=0; i<rows*cols; ++i)
      std::cout << i << ": " << ptr[i] << std::endl;
#endif


    for (Index<Dim> idx; valid(len,idx); next(len, idx))
    {
      test_assert(equal(ptr[offset(idx, ext)], get(view, idx)));
      ptr[offset(idx, ext)] = T(value2(idx));
    }
  }

  if (alloc)
    delete[] buffer;

  if (sync == SYNC_INOUT)
  {
    for (Index<Dim> idx; valid(len,idx); next(len, idx))
      test_assert(equal(get(view, idx), T(value2(idx))));
  }
  else if (sync == SYNC_IN_NOPRESERVE)
  {
    for (Index<Dim> idx; valid(len,idx); next(len, idx))
      test_assert(equal(get(view, idx), T(value1(idx))));
  }
}



// Test run-time external data access (assuming that data is complex,
// either interleaved or split).

template <typename       T,
	  typename       LayoutT,
	  dimension_type Dim>
void
t_rtex_c(
  Domain<Dim> const& dom,
  Rt_tuple           order,
  rt_pack_type       pack,
  rt_complex_type    cformat,
  int                cost,
  bool               alloc,
  sync_action_type   sync)
{
  typedef impl::Fast_block<Dim, T, LayoutT>              block_type;
  typedef typename impl::View_of_dim<Dim, T, block_type>::type view_type;

  view_type view = create_view<view_type>(dom);

  Rt_layout<Dim> blk_rtl = vsip::impl::block_layout<Dim>(view.block());
  test_layout<LayoutT>(blk_rtl);

  Length<Dim> len = impl::extent(dom);
  for (Index<Dim> idx; valid(len, idx); next(len, idx))
    put(view, idx, T(value1(idx)));

  Rt_layout<Dim> rt_layout;

  rt_layout.pack    = pack;
  rt_layout.order   = order; 
  rt_layout.complex = cformat;
  rt_layout.align   = (pack == stride_unit_align) ? 16 : 0;

  // Pre-allocate temporary buffer.
  T* buffer = 0;
  if (alloc)
  {
    Applied_layout<Rt_layout<Dim> > app_layout(rt_layout, len, sizeof(T));
    length_type total_size = app_layout.total_size();
    buffer = new T[total_size];
  }

  {
    Rt_ext_data<block_type> ext(view.block(), rt_layout, sync, buffer);

#if VERBOSE
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "cost: " << ext.cost() << std::endl;
#endif

    test_assert(cost == ext.cost());

    if (rt_layout.complex == cmplx_inter_fmt)
    {
      T* ptr = ext.data().as_inter();
#if VERBOSE
      for (index_type i=0; i<rows*cols; ++i)
	std::cout << i << ": " << ptr[i] << std::endl;
#endif
      if (alloc && cost != 0)
	test_assert(ptr == buffer);

      for (Index<Dim> idx; valid(len,idx); next(len, idx))
      {
	test_assert(equal(ptr[offset(idx, ext)], get(view, idx)));
	ptr[offset(idx, ext)] = T(value2(idx));
      }
    }
    else /* rt_layout.complex == cmplx_split_fmt */
    {
      typedef Storage<Cmplx_split_fmt, T> storage_type;
      typedef typename vsip::impl::Scalar_of<T>::type scalar_type;
      std::pair<scalar_type*,scalar_type*> ptr = ext.data().as_split();
#if VERBOSE
      for (index_type i=0; i<rows*cols; ++i)
	std::cout << i << ": " << ptr.first[i] << "," << ptr.second[i]
		  << std::endl;
#endif
      if (alloc && cost != 0) 
	test_assert(reinterpret_cast<T*>(ptr.first) == buffer);

      for (Index<Dim> idx; valid(len,idx); next(len, idx))
      {
	test_assert(
	  equal(storage_type::get(ptr, offset(idx, ext)), get(view, idx)));
	storage_type::put(ptr, offset(idx, ext), T(value2(idx)));
	}
    }
  }

  if (alloc)
    delete[] buffer;

  if (sync == SYNC_INOUT)
  {
    for (Index<Dim> idx; valid(len,idx); next(len, idx))
      test_assert(equal(get(view, idx), T(value2(idx))));
  }
  else if (sync == SYNC_IN_NOPRESERVE)
  {
    for (Index<Dim> idx; valid(len,idx); next(len, idx))
      test_assert(equal(get(view, idx), T(value1(idx))));
  }
}
			


template <typename T>
void
test_noncomplex(
  Domain<2> const& d,		// size of matrix
  bool             a)		// pre-allocate buffer or not.
{
  using vsip::impl::Layout;
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Cmplx_inter_fmt;

  Rt_tuple r1_v = Rt_tuple(row1_type());
  Rt_tuple r2_v = Rt_tuple(row2_type());
  Rt_tuple c2_v = Rt_tuple(col2_type());

  typedef row1_type r1_t;
  typedef row2_type r2_t;
  typedef col2_type c2_t;

  typedef Stride_unit_dense sud_t;

  typedef Cmplx_inter_fmt cif_t;
  typedef Cmplx_split_fmt csf_t;

  rt_complex_type cif_v = cmplx_inter_fmt;
  rt_complex_type csf_v = cmplx_split_fmt;

  rt_pack_type sud_v = stride_unit_dense;

  sync_action_type sio = SYNC_INOUT;
  // sync_action_type sin = SYNC_IN_NOPRESERVE;

  Domain<1> d1(d[0]);

  // Ask for cmplx_inter_fmt
  t_rtex<T, Layout<1,r1_t,sud_t,cif_t> >(d1,r1_v,sud_v,cif_v,a,sio,0);
  t_rtex<T, Layout<1,r1_t,sud_t,csf_t> >(d1,r1_v,sud_v,cif_v,a,sio,0);

  // Check that cmplx_split_fmt is ignored since type is non-complex.
  t_rtex<T, Layout<1,r1_t,sud_t,cif_t> >(d1,r1_v,sud_v,csf_v,a,sio,0);
  t_rtex<T, Layout<1,r1_t,sud_t,csf_t> >(d1,r1_v,sud_v,csf_v,a,sio,0);

  t_rtex<T, Layout<2,r2_t,sud_t,cif_t> >(d,r2_v,sud_v,cif_v,a,sio,0);
  t_rtex<T, Layout<2,r2_t,sud_t,cif_t> >(d,c2_v,sud_v,cif_v,a,sio,2);
  t_rtex<T, Layout<2,c2_t,sud_t,cif_t> >(d,r2_v,sud_v,cif_v,a,sio,2);
  t_rtex<T, Layout<2,c2_t,sud_t,cif_t> >(d,c2_v,sud_v,cif_v,a,sio,0);
}
  

template <typename       T,
	  dimension_type D>
void
test(
  Domain<D> const& d,		// size of matrix
  bool             a)		// pre-allocate buffer or not.
{
  typedef complex<T> CT;

  using vsip::impl::Layout;
  using vsip::impl::Stride_unit_dense;
  using vsip::impl::Cmplx_inter_fmt;

  typedef typename impl::Row_major<D>::type r_t;
  typedef typename impl::Col_major<D>::type c_t;

  Rt_tuple r_v = Rt_tuple(r_t());
  Rt_tuple c_v = Rt_tuple(c_t());

  typedef Stride_unit_dense sud_t;

  typedef Cmplx_inter_fmt cif_t;
  typedef Cmplx_split_fmt csf_t;

  rt_pack_type sud_v = stride_unit_dense;

  rt_complex_type cif_v = cmplx_inter_fmt;
  rt_complex_type csf_v = cmplx_split_fmt;

  sync_action_type sio = SYNC_INOUT;
  sync_action_type sin = SYNC_IN_NOPRESERVE;

  // SYNC_IN_OUT cases --------------------------------------------------

  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,r_v,sud_v,cif_v,0,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,r_v,sud_v,csf_v,2,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,r_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,r_v,sud_v,csf_v,0,a,sio);

  if (D > 1)
  {
  // These tests only make sense if r_t and c_t are different.
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,c_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,c_v,sud_v,csf_v,2,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,c_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,c_v,sud_v,csf_v,2,a,sio);

  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,c_v,sud_v,cif_v,0,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,c_v,sud_v,csf_v,2,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,c_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,c_v,sud_v,csf_v,0,a,sio);

  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,r_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,r_v,sud_v,csf_v,2,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,r_v,sud_v,cif_v,2,a,sio);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,r_v,sud_v,csf_v,2,a,sio);
  }


  // SYNC_IN_NOPRESERVE cases -------------------------------------------

  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,r_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,r_v,sud_v,csf_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,r_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,r_v,sud_v,csf_v,2,a,sin);

  if (D > 1)
  {
  // These tests only make sense if r_t and c_t are different.
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,c_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,cif_t> > (d,c_v,sud_v,csf_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,c_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,r_t,sud_t,csf_t> > (d,c_v,sud_v,csf_v,2,a,sin);

  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,c_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,c_v,sud_v,csf_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,c_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,c_v,sud_v,csf_v,2,a,sin);

  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,r_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,cif_t> > (d,r_v,sud_v,csf_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,r_v,sud_v,cif_v,2,a,sin);
  t_rtex_c<CT, Layout<D,c_t,sud_t,csf_t> > (d,r_v,sud_v,csf_v,2,a,sin);
  }
}



template <dimension_type D>
void
test_types(Domain<D> const& dom, bool alloc)
{
  test<float> (dom, alloc);
#if VSIP_IMPL_TEST_LEVEL > 0
  test<short> (dom, alloc);
  test<int>   (dom, alloc);
  test<double>(dom, alloc);
#endif
}



int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  test_noncomplex<float>(Domain<2>(4, 8), true);
  test_noncomplex<float>(Domain<2>(4, 8), false);

  test_types(Domain<1>(4), true);
  test_types(Domain<1>(4), false);

  test_types(Domain<2>(4, 8), true);
  test_types(Domain<2>(4, 8), false);

  test_types(Domain<3>(6, 8, 12), true);
  test_types(Domain<3>(6, 8, 12), false);
}

