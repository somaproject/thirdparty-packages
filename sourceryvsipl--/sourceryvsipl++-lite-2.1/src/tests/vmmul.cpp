/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/vmmul.cpp
    @author  Jules Bergmann
    @date    2005-08-15
    @brief   VSIPL++ Library: Unit tests for vector-matrix multiply.
*/

/***********************************************************************
  Included Files
***********************************************************************/


#include <iostream>

#include <vsip/support.hpp>
#include <vsip/initfin.hpp>
#include <vsip/vector.hpp>
#include <vsip/map.hpp>
#include <vsip/parallel.hpp>
#include <vsip/math.hpp>
#include <vsip/domain.hpp>
#include <vsip/signal.hpp>
#include <vsip/core/domain_utils.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/error_db.hpp>
#include "util-par.hpp"
#include "solver-common.hpp"

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Serial tests
***********************************************************************/

template <dimension_type Dim,
	  typename       OrderT2,       // input Matrix dim order
	  typename       OrderT3,       // output Matrix dim order
	  typename       T1,		// input Vector value type
	  typename       T2>		// input Matrix value type
void
test_vmmul(
  length_type rows,
  length_type cols)
{
  Vector<T1> v(Dim == 0 ? cols : rows);
  Matrix<T2, Dense<2, T2, OrderT2> > m(rows, cols);

  typedef typename Promotion<T1, T2>::type result_type;

  for (index_type r=0; r<rows; ++r)
    for (index_type c=0; c<cols; ++c)
      m(r, c) = T2(r*cols+c);

  v = test_ramp(T1(), T1(1), v.size());

  typedef Dense<2, T2, OrderT3> output_block_type;
  Matrix<result_type, output_block_type> res1 =  vmmul<Dim>(      v,       m);
  Matrix<result_type, output_block_type> res2 =  vmmul<Dim>(T1(2)*v,       m);
  Matrix<result_type, output_block_type> res3 =  vmmul<Dim>(      v, T2(3)*m);
  Matrix<result_type, output_block_type> res4 = -vmmul<Dim>(      v,       m);

  for (index_type r=0; r<rows; ++r)
    for (index_type c=0; c<cols; ++c)
      if (Dim == 0)
      {
	test_assert(equal(res1(r, c),  result_type(c   * (r*cols+c))));
	test_assert(equal(res2(r, c),  result_type(2*c * (r*cols+c))));
	test_assert(equal(res3(r, c),  result_type(3*c * (r*cols+c))));
	test_assert(equal(res4(r, c), -result_type(c   * (r*cols+c))));
      }
      else
      {
	test_assert(equal(res1(r, c),  result_type(r   * (r*cols+c))));
	test_assert(equal(res2(r, c),  result_type(2*r * (r*cols+c))));
	test_assert(equal(res3(r, c),  result_type(3*r * (r*cols+c))));
	test_assert(equal(res4(r, c), -result_type(r   * (r*cols+c))));
      }
}



template <typename T,
          typename OrderT,
	  int      SD>
struct test_vmmul_subview
{
  int ops(length_type rows, length_type cols)
    { return rows * cols * vsip::impl::Ops_info<T>::mul; }

  void operator()(length_type rows, length_type cols)
  {
    // Test with input and output dimension ordering the same
    typedef Dense<2, T, OrderT>   block_type;
    typedef Matrix<T, block_type> matrix_type;

    length_type const vector_length = (SD == row ? dom_[1].size() : dom_[0].size());

    Vector<T>   W(vector_length);
    matrix_type A(rows, cols, T());
    matrix_type Z(rows, cols, T());
    matrix_type ref(rows, cols, T());
    typename matrix_type::subview_type sub_A = A(dom_);
    typename matrix_type::subview_type sub_Z = Z(dom_);
    typename matrix_type::subview_type sub_ref = ref(dom_);

    // Compute reference
    W = ramp(T(1), T(1), W.size());
    A = T(1);

    if (SD == row)
      for (index_type r=0; r<rows; r += dom_[0].stride())
        ref.row(r)(dom_[1]) = W * A.row(r)(dom_[1]);
    else
      for (index_type c=0; c<cols; c += dom_[1].stride())
        ref.col(c)(dom_[0]) = W * A.col(c)(dom_[0]);

    // Compute result
    sub_Z = vmmul<SD>(W, sub_A);

    test_assert(error_db(ref, Z) < -100);
  }

  test_vmmul_subview(Domain<2> const& dom) { dom_ = dom; }

private:
  Domain<2> dom_;
};


template <typename T,
          typename OrderT,
	  int      SD>
void
test_subview(Domain<2> const& dom)
{
  const length_type rows = dom[0].stride() * dom[0].size();
  const length_type cols = dom[1].stride() * dom[1].size();

  test_vmmul_subview<T, OrderT, SD> test(dom);
  test(rows, cols);
}


template <typename T,
          typename OrderT,
	  int      SD>
void
test_subview_unaligned(index_type align, stride_type gap)
{
  test_assert(align <= static_cast<index_type>(gap));
  length_type const rows = 8;
  length_type const cols = 32 + gap;
  Domain<2> dom(rows, Domain<1>(align, 1, 32));

  test_vmmul_subview<T, OrderT, SD> test(dom);
  test(rows, cols);
}


template <typename T1,
	  typename T2>
void
vmmul_cases()
{
  // Tests SD, input order, output order, vector type, matrix type.
  test_vmmul<0, row2_type, row2_type, T1, T2>(5, 7);
  test_vmmul<0, col2_type, row2_type, T1, T2>(5, 7);
  test_vmmul<1, row2_type, row2_type, T1, T2>(5, 7);
  test_vmmul<1, col2_type, row2_type, T1, T2>(5, 7);

  test_vmmul<0, row2_type, col2_type, T1, T2>(5, 7);
  test_vmmul<0, col2_type, col2_type, T1, T2>(5, 7);
  test_vmmul<1, row2_type, col2_type, T1, T2>(5, 7);
  test_vmmul<1, col2_type, col2_type, T1, T2>(5, 7);

  // This tests the maximum vector length for the Cell BE dispatch.
  test_vmmul<0, row2_type, row2_type, T1, T2>(8, 8192);

  // Tests various subviews with well-behaved domains (these should
  // still be dispatched to backends that require unit-stride
  // in the major dimension, but not in the minor.
  test_subview<T1, row2_type, row>(Domain<2>(8, 32));
  test_subview<T1, row2_type, row>(Domain<2>(Domain<1>(0, 2, 4), 32));
  test_subview<T1, col2_type, col>(Domain<2>(32, 8));
  test_subview<T1, col2_type, col>(Domain<2>(32, Domain<1>(0, 2, 4)));

  // Tests to ensure alignment is checked by the backend.  These
  // are specifically for cases like Cell BE that cannot handle
  // unaligned transfers.
  index_type align;
  stride_type gap;

  align = 0;  gap = 1;
  test_subview_unaligned<T1, row2_type, row>(align, gap);

  align = 1;  gap = 1;
  test_subview_unaligned<T1, row2_type, row>(align, gap);
}



/***********************************************************************
  Parallel tests
***********************************************************************/

template <dimension_type VecDim,
	  dimension_type Dim>
class Check_vmmul
{
public:
  Check_vmmul(vsip::Domain<Dim> const& dom) : dom_(dom) {}

  template <typename T>
  T operator()(T value,
	       vsip::Index<2> const& /*local*/,
	       vsip::Index<2> const& global)
  {
    vsip::index_type i = global[0]*dom_[1].length()+global[1];
    T expected = (VecDim == 0) ? T(global[1] * i) : T(global[0] * i);

    if (value != expected)
    {
      std::cout << "Check_vmmul: MISCOMPARE" << std::endl
		<< "  global  = " << global[0] << ", " << global[1] 
		<< std::endl
		<< "  expected = " << expected << std::endl
		<< "  actual   = " << value << std::endl;
    }
    return value;
  }

private:
  vsip::Domain<Dim> dom_;
};

template <dimension_type Dim,
	  typename       OrderT,
	  typename       T,
	  typename       VecMapT,
	  typename       MatMapT>
void
test_par_vmmul(
  VecMapT const& vec_map,
  MatMapT const& mat_map,
  length_type    rows,
  length_type    cols)
{
  Matrix<T, Dense<2, T, OrderT, MatMapT> >    m  (rows, cols, mat_map);
  Matrix<T, Dense<2, T, OrderT, MatMapT> >    res(rows, cols, mat_map);
  Vector<T, Dense<1, T, row1_type, VecMapT> > v(Dim == 0 ? cols : rows,
						vec_map);

  foreach_point(m, Set_identity<2>(Domain<2>(rows, cols)));
  foreach_point(v, Set_identity<1>(Domain<1>(v.size())));

  res = vmmul<Dim>(v, m);

  foreach_point(res, Check_vmmul<Dim, 2>(Domain<2>(rows, cols)));
}



template <typename OrderT,
	  typename T>
void
par_vmmul_cases()
{
  length_type np, nr, nc;

  get_np_square(np, nr, nc);


  // -------------------------------------------------------------------
  // If vector is global (replicated on all processors),
  // The matrix must not be distributed along the vector

  Global_map<1> gmap;
  Map<Block_dist, Block_dist> row_map(np, 1);
  Map<Block_dist, Block_dist> col_map(1,  np);
  Map<Block_dist, Block_dist> chk_map(nr, nc);

  test_par_vmmul<0, OrderT, T>(gmap, row_map, 5, 7);
  // test_par_vmmul<1, OrderT, T>(gmap, row_map, 4, 3); // dist along vector

  // test_par_vmmul<0, OrderT, T>(gmap, col_map, 5, 7); // dist along vector
  test_par_vmmul<1, OrderT, T>(gmap, col_map, 5, 7);

  // test_par_vmmul<0, OrderT, T>(gmap, chk_map, 5, 7); // dist along vector
  // test_par_vmmul<1, OrderT, T>(gmap, chk_map, 5, 7); // dist along vector

  // Likewise for replicated_map
  Replicated_map<1> rmap;
  test_par_vmmul<0, OrderT, T>(rmap, row_map, 5, 7);
  test_par_vmmul<1, OrderT, T>(rmap, col_map, 5, 7);


  // -------------------------------------------------------------------
  // If vector is distributed (not replicated),
  // The matrix must
  //    have the same distribution along the vector
  //    not be distributed in the perpendicular to the vector

  Map<Block_dist> vmap(np);

  test_par_vmmul<0, OrderT, T>(vmap, col_map, 5, 7);
  test_par_vmmul<1, OrderT, T>(vmap, row_map, 5, 7);

  // -------------------------------------------------------------------
  // If vector and matrix are both on single processor
  for (processor_type p=0; p<np; ++p)
  {
    Vector<processor_type> pvec(1);
    pvec(0) = p;
    Map<Block_dist, Block_dist> p1_map(pvec, 1, 1); 

    test_par_vmmul<0, OrderT, T>(p1_map, p1_map, 5, 7);
    test_par_vmmul<1, OrderT, T>(p1_map, p1_map, 5, 7);
  }
}



int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  vmmul_cases<         float,          float >();
  vmmul_cases<complex<float>, complex<float> >();
  vmmul_cases<         float, complex<float> >();

  par_vmmul_cases<row2_type, float>();
  par_vmmul_cases<col2_type, float>();
  par_vmmul_cases<row2_type, complex<float> >();
  par_vmmul_cases<col2_type, complex<float> >();
}
