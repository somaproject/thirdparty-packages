/* Copyright (c) 2006, 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/sal/vmul.cpp
    @author  Don McCoy
    @date    2005-01-23
    @brief   VSIPL++ Library: Benchmark for SAL vector multiply.

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



/***********************************************************************
  Definitions - vector element-wise multiply
***********************************************************************/

template <typename T,
	  typename ComplexFmt = Cmplx_inter_fmt>
struct t_vmul_sal;

template <int X, typename T>
struct t_vmul_sal_ip;

template <typename ComplexFmt>
struct t_vmul_sal<float, ComplexFmt> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vmul_sal"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<float>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(float); }
  int wiob_per_point(size_t) { return 1*sizeof(float); }
  int mem_per_point(size_t)  { return 3*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());
    Vector<T, block_type>   C(size);

    A(0) = T(3);
    B(0) = T(4);

    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      vmulx( pA, 1, pB, 1, pC, 1, size, 0 );
    t1.stop();
    
    if (pC[0] != 12.f)
    {
      std::cout << "t_vmul_sal: ERROR" << std::endl;
      abort();
    }
    
    time = t1.delta();
  }
};



template <>
struct t_vmul_sal<complex<float>, Cmplx_inter_fmt> : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vmul_sal complex<float> inter"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());
    Vector<T, block_type>   C(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B = gen.randu(size); B(0) = T(4);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();
    
    int conj_flag = 1;  // don't conjugate
    t1.start();
    for (size_t l=0; l<loop; ++l)
      cvmulx( (COMPLEX *)pA, 2, (COMPLEX *)pB, 2, 
                                (COMPLEX *)pC, 2, size, conj_flag, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(C.get(i), A.get(i) * B.get(i)));
    
    time = t1.delta();
  }
};



template <>
struct t_vmul_sal<complex<float>, Cmplx_split_fmt> : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vmul_sal complex<float> split"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;
    typedef impl::Ext_data<block_type>::raw_ptr_type ptr_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());
    Vector<T, block_type>   C(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B = gen.randu(size); B(0) = T(4);

    vsip::impl::profile::Timer t1;
    
    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    ptr_type pA = ext_a.data();
    ptr_type pB = ext_b.data();
    ptr_type pC = ext_c.data();

    int conj_flag = 1;  // don't conjugate
    t1.start();
    for (size_t l=0; l<loop; ++l)
      zvmulx((COMPLEX_SPLIT*)&pA, 1,
	     (COMPLEX_SPLIT*)&pB, 1, 
	     (COMPLEX_SPLIT*)&pC, 1, size, conj_flag, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(C.get(i), A.get(i) * B.get(i)));
    
    time = t1.delta();
  }
};



template <>
struct t_vmul_sal_ip<1, float> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vmul_sal_ip float"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 2*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T(1));

    A(0) = T(3);

    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      vmulx( pA, 1, pB, 1, pA, 1, size, 0 );
    t1.stop();
    
    if (pA[0] != 3.f)
    {
      std::cout << "t_vmul_sal: ERROR" << std::endl;
      abort();
    }
    
    time = t1.delta();
  }
};



template <>
struct t_vmul_sal_ip<1, complex<float> > : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vmul_sal complex<float>"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 2*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());

    A(0) = T(1, 0);
    B(0) = T(1, 0);

    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();

    vsip::impl::profile::Timer t1;
    
    int conj_flag = 1;  // don't conjugate
    t1.start();
    for (size_t l=0; l<loop; ++l)
      cvmulx( (COMPLEX *)pA, 2, (COMPLEX *)pB, 2, 
                                (COMPLEX *)pB, 2, size, conj_flag, 0 );
    t1.stop();
    
    if (pB[0].real() != 1.f)
    {
      std::cout << "t_vmul_sal: ERROR" << std::endl;
      abort();
    }
    
    time = t1.delta();
  }
};



template <>
struct t_vmul_sal_ip<2, complex<float> > : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_vmul_sal complex<float>"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 2*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(size, T());

    A(0) = T(1, 0);
    B(0) = T(1, 0);

    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();

    vsip::impl::profile::Timer t1;
    
    int conj_flag = 1;  // don't conjugate
    t1.start();
    for (size_t l=0; l<loop; ++l)
      cvmulx( (COMPLEX *)pA, 2, (COMPLEX *)pB, 2, 
                                (COMPLEX *)pA, 2, size, conj_flag, 0 );
    t1.stop();
    
    if (pB[0].real() != 1.f)
    {
      std::cout << "t_vmul_sal: ERROR" << std::endl;
      abort();
    }
    
    time = t1.delta();
  }
};



/***********************************************************************
  Definitions - scalar-vector element-wise multiply
***********************************************************************/

template <typename ScalarT,
	  typename T,
	  typename ComplexFmt = Cmplx_inter_fmt>
struct t_svmul_sal;

template <typename ComplexFmt>
struct t_svmul_sal<float, float, ComplexFmt> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_svmul_sal"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<float>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(float); }
  int wiob_per_point(size_t) { return 1*sizeof(float); }
  int mem_per_point(size_t)  { return 3*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(1, T());
    Vector<T, block_type>   C(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B(0) = T(4);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      vsmulx( pA, 1, pB, pC, 1, size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(C.get(i), A.get(i) * B.get(0)));
    
    time = t1.delta();
  }
};



template <>
struct t_svmul_sal<complex<float>, complex<float>, Cmplx_inter_fmt>
  : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_svmul_sal"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(1, T());
    Vector<T, block_type>   C(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B(0) = T(4);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    T* pA = ext_a.data();
    T* pB = ext_b.data();
    T* pC = ext_c.data();
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      cvcsmlx((COMPLEX*)pA, 2,
	      (COMPLEX*)pB,
	      (COMPLEX*)pC, 2,
	      size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(C.get(i), A.get(i) * B.get(0)));
    
    time = t1.delta();
  }
};



template <>
struct t_svmul_sal<complex<float>, complex<float>, Cmplx_split_fmt>
  : Benchmark_base
{
  typedef complex<float> T;

  char* what() { return "t_svmul_sal"; }
  int ops_per_point(size_t)  { return vsip::impl::Ops_info<T>::mul; }
  int riob_per_point(size_t) { return 2*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 3*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;
    typedef impl::Ext_data<block_type>::raw_ptr_type ptr_type;

    Vector<T, block_type>   A(size, T());
    Vector<T, block_type>   B(1, T());
    Vector<T, block_type>   C(size);

    Rand<T> gen(0, 0);
    A = gen.randu(size); A(0) = T(3);
    B(0) = T(4);

    vsip::impl::profile::Timer t1;

    {
    impl::Ext_data<block_type> ext_a(A.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_b(B.block(), impl::SYNC_IN);
    impl::Ext_data<block_type> ext_c(C.block(), impl::SYNC_OUT);
    
    ptr_type pA = ext_a.data();
    ptr_type pB = ext_b.data();
    ptr_type pC = ext_c.data();
    
    t1.start();
    for (size_t l=0; l<loop; ++l)
      zvzsmlx((COMPLEX_SPLIT*)&pA, 1,
	      (COMPLEX_SPLIT*)&pB,
	      (COMPLEX_SPLIT*)&pC, 1,
	      size, 0 );
    t1.stop();
    }

    for (index_type i=0; i<size; ++i)
      test_assert(equal(C.get(i), A.get(i) * B.get(0)));
    
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
  typedef complex<float> cf_type;
  switch (what)
  {
  case  1: loop(t_vmul_sal<float>()); break;
  case  2: loop(t_vmul_sal<complex<float>, Cmplx_inter_fmt>()); break;
  case  3: loop(t_vmul_sal<complex<float>, Cmplx_split_fmt>()); break;

  case 11: loop(t_svmul_sal<float, float>()); break;
  case 13: loop(t_svmul_sal<cf_type, cf_type, Cmplx_inter_fmt>()); break;
  case 14: loop(t_svmul_sal<cf_type, cf_type, Cmplx_split_fmt>()); break;

  case 31: loop(t_vmul_sal_ip<1, float>()); break;
  case 32: loop(t_vmul_sal_ip<1, complex<float> >()); break;
  case 33: loop(t_vmul_sal_ip<2, complex<float> >()); break;
  }
}
