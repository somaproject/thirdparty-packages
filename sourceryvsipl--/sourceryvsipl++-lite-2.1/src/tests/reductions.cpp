/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/reductions.cpp
    @author  Jules Bergmann
    @date    2005-07-11
    @brief   VSIPL++ Library: Tests for math reductions.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/map.hpp>
#include <vsip/parallel.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/test-storage.hpp>

using namespace vsip;
using namespace vsip_csl;


void
simple_tests()
{
  Vector<float> vec(4);

  vec(0) = 0.;
  vec(1) = 1.;
  vec(2) = 2.;
  vec(3) = 3.;

  test_assert(equal(sumval(vec),    6.0f));
  test_assert(equal(meanval(vec),   1.5f));
  test_assert(equal(sumsqval(vec), 14.0f));
  test_assert(equal(meansqval(vec), 3.5f));

  test_assert(equal(sumval(vec+vec), 12.0f));

  Matrix<double> mat(2, 2);

  mat(0, 0) = 1.;
  mat(0, 1) = 2.;
  mat(1, 0) = 3.;
  mat(1, 1) = 4.;

  test_assert(equal(sumval(mat),   10.0));
  test_assert(equal(meanval(mat),   2.5));
  test_assert(equal(sumsqval(mat), 30.0));
  test_assert(equal(meansqval(mat), 7.5));

  Tensor<float> ten(2, 1, 2);

  ten(0, 0, 0) = 2.;
  ten(0, 0, 1) = 3.;
  ten(1, 0, 0) = 4.;
  ten(1, 0, 1) = 5.;

  test_assert(equal(sumval(ten),    14.0f));
  test_assert(equal(meanval(ten),    3.5f));
  test_assert(equal(sumsqval(ten),  54.0f));
  test_assert(equal(meansqval(ten), 13.5f));

  Vector<complex<float> > cvec(2);

  cvec(0) = complex<float>(3.f,  4.f); // -7 + 24i
  cvec(1) = complex<float>(3.f, -4.f); // -7 - 24i

  test_assert(equal(sumval(cvec),    complex<float>(6.0f, 0.0f)));
  // test_assert(equal(meanval(cvec), complex<float>(3.f, 0.f)));
  test_assert(equal(sumsqval(cvec),  complex<float>(-14.0f, 0.0f)));
  test_assert(equal(meansqval(cvec), 25.0f));


  Vector<bool> bvec(4);

  bvec(0) = true;
  bvec(1) = true;
  bvec(2) = false;
  bvec(3) = true;

  test_assert(equal(sumval(bvec), static_cast<length_type>(3)));
}



/***********************************************************************
  sumval tests.
***********************************************************************/

template <typename ViewT>
void
view_sumval(
  ViewT       view,
  length_type count)
{
  typedef typename ViewT::value_type T;

  view = T();

  index_type  i        = 0;
  T           expected = T();
  length_type size     = view.size();
  
  for (index_type c=0; c<count; ++c)
  {
    i      = (2*i+3) % size;
    T nval = T(i) - T(5);

    expected -= get_nth(view, i);
    expected += nval;

    put_nth(view, i, nval);
    
    T val = sumval(view);
    test_assert(equal(val, expected));
  }
}



template <typename       StoreT,
	  dimension_type Dim>
void
test_sumval(Domain<Dim> const& dom, length_type count)
{
  typedef typename StoreT::value_type T;

  StoreT      store(dom, T());

  view_sumval(store.view, count);
}



template <typename T>
void
cover_sumval()
{
  test_sumval<Storage<1, T> >(Domain<1>(15), 8);
  
  test_sumval<Storage<2, T, row2_type> >(Domain<2>(15, 17), 8);
  test_sumval<Storage<2, T, col2_type> >(Domain<2>(15, 17), 8);
  
  test_sumval<Storage<3, T, tuple<0, 1, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumval<Storage<3, T, tuple<0, 2, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumval<Storage<3, T, tuple<1, 0, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumval<Storage<3, T, tuple<1, 2, 0> > >(Domain<3>(15, 17, 7), 8);
  test_sumval<Storage<3, T, tuple<2, 0, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumval<Storage<3, T, tuple<2, 1, 0> > >(Domain<3>(15, 17, 7), 8);

  test_sumval<Storage<1, T, row1_type, Map<Block_dist> > >(Domain<1>(15), 8);
  test_sumval<Storage<1, T, row1_type, Global_map<1> > >(Domain<1>(15), 8);
}



template <typename T,
	  typename MapT>
void
par_cover_sumval()
{
  typedef Dense<1, T, row1_type, MapT> block_type;
  typedef Vector<T, block_type>        view_type;

  length_type size = 8;

  MapT      map = create_map<1, MapT>();
  view_type view(size, map);

  view_sumval(view, 8);
}



/***********************************************************************
  sumval bool tests.
***********************************************************************/

template <typename       StoreT,
	  dimension_type Dim>
void
test_sumval_bool(Domain<Dim> const& dom, length_type count)
{
  StoreT      store(dom, false);
  length_type size = store.view.size();

  index_type  i        = 0;
  length_type expected = 0;
  
  for (index_type c=0; c<count; ++c)
  {
    i         = (2*i+3) % size;
    bool nval = (3*i+1) % 2 == 0;

    if (get_nth(store.view, i))
      expected -= 1;

    if (nval)
      expected += 1;

    put_nth(store.view, i, nval);
    
    length_type val = sumval(store.view);
    test_assert(equal(val, expected));
  }
}



void
cover_sumval_bool()
{
  typedef bool T;

  test_sumval_bool<Storage<1, T> >(Domain<1>(15), 8);
  
  test_sumval_bool<Storage<2, T, row2_type> >(Domain<2>(15, 17), 8);
  test_sumval_bool<Storage<2, T, col2_type> >(Domain<2>(15, 17), 8);
  
  test_sumval_bool<Storage<3, T, tuple<0, 1, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumval_bool<Storage<3, T, tuple<0, 2, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumval_bool<Storage<3, T, tuple<1, 0, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumval_bool<Storage<3, T, tuple<1, 2, 0> > >(Domain<3>(15, 17, 7), 8);
  test_sumval_bool<Storage<3, T, tuple<2, 0, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumval_bool<Storage<3, T, tuple<2, 1, 0> > >(Domain<3>(15, 17, 7), 8);

  test_sumval_bool<Storage<1, T, row1_type, Map<Block_dist> > >(Domain<1>(15), 8);
  test_sumval_bool<Storage<1, T, row1_type, Global_map<1> > >(Domain<1>(15), 8);
}



/***********************************************************************
  sumsqval tests.

  Note: sumsqval returns the sum of squares of elements of a view.
***********************************************************************/

template <typename       StoreT,
	  dimension_type Dim>
void
test_sumsqval(Domain<Dim> const& dom, length_type count)
{
  typedef typename StoreT::value_type T;

  StoreT      store(dom, T());
  length_type size = store.view.size();

  index_type  i        = 0;
  T           expected = T();
  
  for (index_type c=0; c<count; ++c)
  {
    i      = (2*i+3) % size;
    T nval = T(i) - T(5);

    T nth  = get_nth(store.view, i);
    expected -= (nth  * nth);
    expected += (nval * nval);

    put_nth(store.view, i, nval);
    
    T val = sumsqval(store.view);
    test_assert(equal(val, expected));
  }
}



template <typename T>
void
cover_sumsqval()
{
  test_sumsqval<Storage<1, T> >(Domain<1>(15), 8);
  
  test_sumsqval<Storage<2, T, row2_type> >(Domain<2>(15, 17), 8);
  test_sumsqval<Storage<2, T, col2_type> >(Domain<2>(15, 17), 8);
  
  test_sumsqval<Storage<3, T, tuple<0, 1, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumsqval<Storage<3, T, tuple<0, 2, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumsqval<Storage<3, T, tuple<1, 0, 2> > >(Domain<3>(15, 17, 7), 8);
  test_sumsqval<Storage<3, T, tuple<1, 2, 0> > >(Domain<3>(15, 17, 7), 8);
  test_sumsqval<Storage<3, T, tuple<2, 0, 1> > >(Domain<3>(15, 17, 7), 8);
  test_sumsqval<Storage<3, T, tuple<2, 1, 0> > >(Domain<3>(15, 17, 7), 8);

  test_sumsqval<Storage<1, T, row1_type, Map<Block_dist> > >(Domain<1>(15), 8);
  test_sumsqval<Storage<1, T, row1_type, Global_map<1> > >(Domain<1>(15), 8);
}



/***********************************************************************
  meanval tests.
***********************************************************************/

template <typename       StoreT,
	  dimension_type Dim>
void
test_meanval(Domain<Dim> const& dom, length_type count)
{
  typedef typename StoreT::value_type T;

  StoreT      store(dom, T());
  length_type size = store.view.size();

  index_type  i        = 0;
  T           expected = T();
  
  for (index_type c=0; c<count; ++c)
  {
    i      = (2*i+3) % size;
    T nval = T(i) - T(5);

    expected -= get_nth(store.view, i);
    expected += nval;

    put_nth(store.view, i, nval);
    
    T sval = sumval(store.view);
    T mval = meanval(store.view);
    test_assert(equal(sval, expected));
    test_assert(equal(mval, expected/static_cast<T>(store.view.size())));
  }
}



template <typename T>
void
cover_meanval()
{
  test_meanval<Storage<1, T> >(Domain<1>(15), 8);
  
  test_meanval<Storage<2, T, row2_type> >(Domain<2>(15, 17), 8);
  test_meanval<Storage<2, T, col2_type> >(Domain<2>(15, 17), 8);
  
  test_meanval<Storage<3, T, tuple<0, 1, 2> > >(Domain<3>(15, 17, 7), 8);
  test_meanval<Storage<3, T, tuple<0, 2, 1> > >(Domain<3>(15, 17, 7), 8);
  test_meanval<Storage<3, T, tuple<1, 0, 2> > >(Domain<3>(15, 17, 7), 8);
  test_meanval<Storage<3, T, tuple<1, 2, 0> > >(Domain<3>(15, 17, 7), 8);
  test_meanval<Storage<3, T, tuple<2, 0, 1> > >(Domain<3>(15, 17, 7), 8);
  test_meanval<Storage<3, T, tuple<2, 1, 0> > >(Domain<3>(15, 17, 7), 8);

  test_meanval<Storage<1, T, row1_type, Map<Block_dist> > >(Domain<1>(15), 8);
  test_meanval<Storage<1, T, row1_type, Global_map<1> > >(Domain<1>(15), 8);
}



/***********************************************************************
  meansqval tests.
***********************************************************************/

template <typename       StoreT,
	  dimension_type Dim>
void
test_meansqval(Domain<Dim> const& dom, length_type count)
{
  typedef typename StoreT::value_type T;
  typedef typename impl::Scalar_of<T>::type result_type;

  StoreT      store(dom, T());
  length_type size = store.view.size();

  index_type  i        = 0;
  result_type expected = result_type();
  
  for (index_type c=0; c<count; ++c)
  {
    i      = (2*i+3) % size;
    T nval = T(i) - T(5);

    T nth  = get_nth(store.view, i);

    expected -= magsq(nth);
    expected += magsq(nval);

    put_nth(store.view, i, nval);
    
    result_type mval = meansqval(store.view);
    test_assert(equal(mval, expected/static_cast<result_type>(store.view.size())));
  }
}



template <typename T>
void
cover_meansqval()
{
  test_meansqval<Storage<1, T> >(Domain<1>(15), 8);
  
  test_meansqval<Storage<2, T, row2_type> >(Domain<2>(15, 17), 8);
  test_meansqval<Storage<2, T, col2_type> >(Domain<2>(15, 17), 8);
  
  test_meansqval<Storage<3, T, tuple<0, 1, 2> > >(Domain<3>(15, 17, 7), 8);
  test_meansqval<Storage<3, T, tuple<0, 2, 1> > >(Domain<3>(15, 17, 7), 8);
  test_meansqval<Storage<3, T, tuple<1, 0, 2> > >(Domain<3>(15, 17, 7), 8);
  test_meansqval<Storage<3, T, tuple<1, 2, 0> > >(Domain<3>(15, 17, 7), 8);
  test_meansqval<Storage<3, T, tuple<2, 0, 1> > >(Domain<3>(15, 17, 7), 8);
  test_meansqval<Storage<3, T, tuple<2, 1, 0> > >(Domain<3>(15, 17, 7), 8);

  test_meansqval<Storage<1, T, row1_type, Map<Block_dist> > >(Domain<1>(15), 8);
  test_meansqval<Storage<1, T, row1_type, Global_map<1> > >(Domain<1>(15), 8);
}



int
main(int argc, char** argv)
{
  vsipl init(argc, argv);
   
  simple_tests();

  par_cover_sumval<float, Global_map<1> >();
  par_cover_sumval<float, Map<Block_dist> >();

  cover_sumval<int>();
  cover_sumval<float>();
  cover_sumval<complex<float> >();

  cover_sumval_bool();

  cover_sumsqval<int>();
  cover_sumsqval<float>();
  cover_sumsqval<complex<float> >();

  cover_meanval<int>();
  cover_meanval<float>();
  cover_meanval<complex<float> >();

  cover_meansqval<int>();
  cover_meansqval<float>();
  cover_meansqval<complex<float> >();

#if VSIP_IMPL_TEST_DOUBLE
  cover_sumval<double>();
  cover_sumval<complex<double> >();

  cover_sumsqval<double>();
  cover_sumsqval<complex<double> >();

  cover_meanval<double>();

  cover_meansqval<double>();
  cover_meansqval<complex<double> >();
#endif
}
