/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/sal/vthresh.cpp
    @author  Jules Bergmann
    @date    2007-06-05
    @brief   VSIPL++ Library: Benchmark for vthresh

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <sal.h>

#include <vsip/random.hpp>
#include <vsip/opt/profile.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/ops_info.hpp>

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

template <typename T>
struct t_vthres_sal;

template <typename T>
struct t_vthr_sal;

template <>
struct t_vthres_sal<float> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vthres_sal"; }
  int ops_per_point(size_t)  { return 1; }

  int riob_per_point(size_t) { return 1*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 2*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A     (size, T());
    Vector<T, block_type>   result(size, T());
    T                       b = T(0.5);

    Rand<T> gen(0, 0);
    A = gen.randu(size);

    vsip::impl::profile::Timer t1;

    { // Control Ext_data scope.
      impl::Ext_data<block_type> ext_A     (A.block(),      impl::SYNC_IN);
      impl::Ext_data<block_type> ext_result(result.block(), impl::SYNC_OUT);
    
      T* p_A      = ext_A.data();
      T* p_result = ext_result.data();
      
      t1.start();
      marker1_start();
      for (size_t l=0; l<loop; ++l)
      {
	vthresx(p_A,1, &b, p_result,1, size, 0);
      }
      marker1_stop();
      t1.stop();
    }
    
    for (index_type i=0; i<size; ++i)
    {
      test_assert(equal<T>(result.get(i), (A(i) >= b) ? A(i) : 0.f));
    }
    
    time = t1.delta();
  }
};



#if VSIP_IMPL_SAL_HAVE_VTHRX
template <>
struct t_vthr_sal<float> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_vthr_sal"; }
  int ops_per_point(size_t)  { return 1; }

  int riob_per_point(size_t) { return 1*sizeof(T); }
  int wiob_per_point(size_t) { return 1*sizeof(T); }
  int mem_per_point(size_t)  { return 2*sizeof(T); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A     (size, T());
    Vector<T, block_type>   result(size, T());
    T                       b = T(0.5);

    Rand<T> gen(0, 0);
    A = gen.randu(size);

    vsip::impl::profile::Timer t1;

    { // Control Ext_data scope.
      impl::Ext_data<block_type> ext_A     (A.block(),      impl::SYNC_IN);
      impl::Ext_data<block_type> ext_result(result.block(), impl::SYNC_OUT);
    
      T* p_A      = ext_A.data();
      T* p_result = ext_result.data();
      
      t1.start();
      marker1_start();
      for (size_t l=0; l<loop; ++l)
      {
	vthrx(p_A,1, &b, p_result,1, size, 0);
      }
      marker1_stop();
      t1.stop();
    }
    
    for (index_type i=0; i<size; ++i)
    {
      test_assert(equal<T>(result.get(i), (A(i) >= b) ? A(i) : b));
    }
    
    time = t1.delta();
  }
};
#endif



template <typename T>
struct t_vthres_c : Benchmark_base
{
  char* what() { return "t_vthres_c"; }
  int ops_per_point(size_t)  { return 1; }

  int riob_per_point(size_t) { return 2*sizeof(float); }
  int wiob_per_point(size_t) { return 1*sizeof(float); }
  int mem_per_point(size_t)  { return 3*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A     (size, T());
    T                       b = T(0.5);
    Vector<T, block_type>   result(size, T());

    Rand<T> gen(0, 0);
    A = gen.randu(size);

    vsip::impl::profile::Timer t1;

    { // Control Ext_data scope.
      impl::Ext_data<block_type> ext_A     (A.block(),      impl::SYNC_IN);
      impl::Ext_data<block_type> ext_result(result.block(), impl::SYNC_OUT);
    
      T* p_A      = ext_A.data();
      T* p_result = ext_result.data();
      
      t1.start();
      for (size_t l=0; l<loop; ++l)
      {
	for (index_type i=0; i<size; ++i)
	  p_result[i] = (p_A[i] >= b) ? p_A[i] : 0.f;
      }
      t1.stop();
    }
    
    for (index_type i=0; i<size; ++i)
    {
      test_assert(equal<T>(result.get(i), (A(i) >= b) ? A(i) : 0.f));
    }
    
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
  case  1: loop(t_vthres_sal<float>()); break;
#if VSIP_IMPL_SAL_HAVE_VTHRX
  case  2: loop(t_vthr_sal<float>()); break;
#endif
  case 11: loop(t_vthres_c<float>()); break;
  case  0:
    std::cout
      << "SAL vthres\n"
      << "  -1 -- SAL vthresx (float) Z(i) = A(i) > b ? A(i) : 0\n"
#if VSIP_IMPL_SAL_HAVE_VTHRX
      << "  -2 -- SAL vthrx   (float) Z(i) = A(i) > b ? A(i) : b\n"
#else
      << " (-2 -- SAL vthrx function not available)\n"
#endif
      << " -11 -- C           (float) Z(i) = A(i) > b ? A(i) : 0\n"
      ;
  }
}
