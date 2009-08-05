/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/vma_sal.cpp
    @author  Jules Bergmann
    @date    2006-06-01
    @brief   VSIPL++ Library: Benchmark for SAL vector multiply-add.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <complex>

#include <vsip/random.hpp>
#include <vsip/opt/profile.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/ops_info.hpp>
#include <sal.h>

#include "loop.hpp"
#include "benchmarks.hpp"

using namespace std;
using namespace vsip;

using impl::Stride_unit_dense;
using impl::Cmplx_inter_fmt;
using impl::Cmplx_split_fmt;


template <typename T,
	  typename ComplexFmt = Cmplx_inter_fmt>
struct t_vma_sal;

template <typename ComplexFmt>
struct t_vma_sal<float, ComplexFmt> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vma_sal"; }
  int ops_per_point(size_t)  
    { return vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add; }
  int riob_per_point(size_t) { return 3*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 4*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());
    Vector<T, block_type>   C(size, T());
    Vector<T, block_type>   Z(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B = gen.randu(size); B(0) = T(4);
    C = gen.randu(size); C(0) = T(5);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_z(Z.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();
    T* pZ = ext_z.data();

    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      vma_x( pA, 1, pB, 1, pC, 1, pZ, 1, size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(Z.get(i), A.get(i) * B.get(i) + C.get(i)));
    
    time = t1.delta();
  }
};



// vsma: (scalar * vector) + vector
template <typename T,
	  typename ComplexFmt = Cmplx_inter_fmt>
struct t_vsma_sal;

template <typename ComplexFmt>
struct t_vsma_sal<float, ComplexFmt> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vsma_sal"; }
  int ops_per_point(size_t)  
    { return vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    T                       b;
    Vector<T, block_type>   C(size, T());
    Vector<T, block_type>   Z(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    b = T(4);
    C = gen.randu(size); C(0) = T(5);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_z(Z.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = &b;
    T* pC = ext_c.data();
    T* pZ = ext_z.data();

    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      vsmax( pA, 1, pB, pC, 1, pZ, 1, size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(Z.get(i), A.get(i) * b + C.get(i)));
    
    time = t1.delta();
  }
};



template <>
struct t_vsma_sal<complex<float>, Cmplx_inter_fmt> : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vsma_sal complex<float> inter"; }
  int ops_per_point(size_t)  
    { return vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(1, T());
    Vector<T, block_type>   C(size, T());
    Vector<T, block_type>   Z(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B(0) = T(4);
    C = gen.randu(size); C(0) = T(5);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_z(Z.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();
    T* pZ = ext_z.data();
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      cvsmax((COMPLEX*)pA, 2,
	     (COMPLEX*)pB, 
	     (COMPLEX*)pC, 2,
	     (COMPLEX*)pZ, 2,
	     size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(Z.get(i), A.get(i) * B.get(0) + C.get(i)));
    
    time = t1.delta();
  }
};



template <>
struct t_vsma_sal<complex<float>, Cmplx_split_fmt> : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vsma_sal complex<float> split"; }
  int ops_per_point(size_t)  
    { return vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;
    typedef impl::Ext_data<block_type>::raw_ptr_type ptr_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());
    Vector<T, block_type>   C(size, T());
    Vector<T, block_type>   Z(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
                         B(0) = T(4);
    C = gen.randu(size); C(0) = T(5);

    vsip::impl::profile::Timer t1;
    
    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_z(Z.block(), impl::SYNC_OUT);
    
    ptr_type pA = ext_a.data();
    ptr_type pB = ext_b.data();
    ptr_type pC = ext_c.data();
    ptr_type pZ = ext_z.data();

    t1.start();
    for (size_t l=0; l<loop; ++l)
      zvsmax((COMPLEX_SPLIT*)&pA, 1,
	     (COMPLEX_SPLIT*)&pB,
	     (COMPLEX_SPLIT*)&pC, 1, 
	     (COMPLEX_SPLIT*)&pZ, 1, size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(Z.get(i), A.get(i) * B.get(0) + C.get(i)));
    
    time = t1.delta();
  }
};



void
defaults(Loop1P&)
{
}



void
test(Loop1P& loop, int what)
{
  switch (what)
  {
  case  1: loop(t_vma_sal<float>()); break;

  case  11: loop(t_vsma_sal<float>()); break;
  case  12: loop(t_vsma_sal<complex<float>, Cmplx_inter_fmt>()); break;
  case  13: loop(t_vsma_sal<complex<float>, Cmplx_split_fmt>()); break;
#if 0
  case  2: loop(t_vmul_sal<complex<float> >()); break;

  case  12: loop(t_vmul_sal_ip<1, complex<float> >()); break;
  case  22: loop(t_vmul_sal_ip<2, complex<float> >()); break;
#endif
  }
}
