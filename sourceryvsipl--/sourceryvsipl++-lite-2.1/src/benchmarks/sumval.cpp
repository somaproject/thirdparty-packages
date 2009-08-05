/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/sumval.cpp
    @author  Jules Bergmann
    @date    2005-07-11
    @brief   VSIPL++ Library: Benchmark for sumval reductions.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/random.hpp>
#include <vsip/opt/profile.hpp>

#include <vsip_csl/test.hpp>
#include "loop.hpp"

using namespace vsip;
using namespace vsip_csl;



/***********************************************************************
  VSIPL++ sumval
***********************************************************************/

template <typename T>
struct t_sumval1 : Benchmark_base
{
  char const* what() { return "t_sumval_vector"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return sizeof(T); }
  int wiob_per_point(length_type) { return 0; }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T>   view(size, T());
    T           val = T();
    
    Rand<T>     gen(0, 0);

    if (init_ == 0)
      view = gen.randu(size);
    else if (init_ == 1)
      view(0) = T(2);
    else if (init_ == 2)
      view(size-1) = T(2);
    
    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
      val = vsip::sumval(view);
    t1.stop();

    if (init_ == 1 || init_ == 2)
      test_assert(equal(val, T(2)));
    
    time = t1.delta();
  }

  t_sumval1(int init) : init_(init) {}

  // member data.
  int init_;
};



template <typename T>
struct t_sumval2 : Benchmark_base
{
  char const* what() { return "t_sumval_matrix32"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return sizeof(T); }
  int wiob_per_point(length_type) { return 0; }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Matrix<T>   view(32, size, T());
    T           val = T();
    
    view(0, 1) = T(2);
    
    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
      val = vsip::sumval(view);
    t1.stop();
    
    assert(equal(val, T(2)));
    if (!equal(val, T(2))) printf("t_sumval: ERROR\n");
    
    time = t1.delta();
  }
};



/***********************************************************************
  get/put sumval
***********************************************************************/

template <typename T>
struct t_sumval_gp : Benchmark_base
{
  char const* what() { return "t_sumval_gp"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return sizeof(T); }
  int wiob_per_point(length_type) { return 0; }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T>   view(size, T());
    T           val = T();
    
    Rand<T>     gen(0, 0);

    if (init_ == 0)
      view = gen.randu(size);
    else if (init_ == 1)
      view(0) = T(2);
    else if (init_ == 2)
      view(size-1) = T(2);
    
    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
    {
      val = T();
      for (index_type i=0; i<size; ++i)
	val += view.get(i);
    }
    t1.stop();

    if (init_ == 1 || init_ == 2)
      test_assert(equal(val, T(2)));
    
    time = t1.delta();
  }

  t_sumval_gp(int init) : init_(init) {}

  // member data.
  int init_;
};



void
defaults(Loop1P&)
{
}



int
test(Loop1P& loop, int what)
{
  switch (what)
  {
  case   1: loop(t_sumval1<float>(0)); break;
  case   2: loop(t_sumval1<float>(1)); break;
  case   3: loop(t_sumval1<float>(2)); break;

  case  11: loop(t_sumval1<int>(0)); break;

  case  21: loop(t_sumval_gp<float>(0)); break;

  case 101: loop(t_sumval2<float>()); break;
  default: return 0;
  }
  return 1;
}
