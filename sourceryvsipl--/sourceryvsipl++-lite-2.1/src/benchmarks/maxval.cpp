/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/maxval.cpp
    @author  Jules Bergmann
    @date    2006-06-01
    @brief   VSIPL++ Library: Benchmark for maxval reductions.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/random.hpp>
#include <vsip/selgen.hpp>
#include <vsip/opt/profile.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/selgen.hpp>

#include <vsip_csl/test.hpp>
#include "loop.hpp"

#include "maxval.hpp"

using namespace vsip;
using namespace vsip_csl;

template <typename T,
          typename MapT = Local_map>
struct t_maxval1
{
  char const* what() { return "t_maxval_vector"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return sizeof(T); }
  int wiob_per_point(length_type) { return 0; }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    using namespace vsip::impl;

    typedef Dense<1,T,row1_type,MapT>                               block_type;

    create_test_vector_helper<MapT,Vector,float,Dense<1,T,row1_type,MapT> >
      ctvh(size);

    T                      val = T();
    Index<1>               idx;

    Rand<T>     gen(0, 0);

    if (init_ == 0)
      ctvh.assign_view(gen.randu(size));
    else if (init_ == 1)
      ctvh.assign_view(ramp(T(0), T(1), size));
    else if (init_ == 2)
      ctvh.assign_view(ramp(T(size-1), T(-1), size));
   
    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
      val = maxval(ctvh.view, idx);
    t1.stop();

    if (init_ == 1)
    {
      test_assert(equal(val, T(size-1)));
      test_assert(idx == size-1);
    }
    else if (init_ == 2)
    {
      test_assert(equal(val, T(size-1)));
      test_assert(idx == 0);
    }
    
    time = t1.delta();
  }

  void diag()
  {
  }

  t_maxval1(int init) : init_(init) {}

  int init_;
};

template <typename T,
          typename MapT = Local_map,
	  typename Tag = impl::Cvsip_tag>
struct t_maxval2
{
  char const* what() { return "t_maxval_vector"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return sizeof(T); }
  int wiob_per_point(length_type) { return 0; }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    using namespace vsip::impl;

    typedef Dense<1,T,row1_type,MapT>                               block_type;
    typedef reduction_op_eval<Max_value,T,block_type,1,Tag> eval;

    create_test_vector_helper<MapT,Vector,float,Dense<1,T,row1_type,MapT> >
      ctvh(size);

    T                      val = T();
    Index<1>               idx;

    Rand<T>     gen(0, 0);

    if (init_ == 0)
      ctvh.assign_view(gen.randu(size));
    else if (init_ == 1)
      ctvh.assign_view(ramp(T(0), T(1), size));
    else if (init_ == 2)
      ctvh.assign_view(ramp(T(size-1), T(-1), size));
   
    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
      eval::exec(val, ctvh.view.block(), idx);
    t1.stop();

    if (init_ == 1)
    {
      test_assert(equal(val, T(size-1)));
      test_assert(idx == size-1);
    }
    else if (init_ == 2)
    {
      test_assert(equal(val, T(size-1)));
      test_assert(idx == 0);
    }
    
    time = t1.delta();
  }

  void diag()
  {
  }

  t_maxval2(int init) : init_(init) {}

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
  case  1: loop(t_maxval1<float>(0)); break;
  case  2: loop(t_maxval1<float>(1)); break;
  case  3: loop(t_maxval1<float>(2)); break;
  case  4: loop(t_maxval2<float,Map<>,impl::Parallel_tag>(0)); break;
  case  5: loop(t_maxval2<float,Map<>,impl::Parallel_tag>(1)); break;
  case  6: loop(t_maxval2<float,Map<>,impl::Parallel_tag>(2)); break;
  case  7: loop(t_maxval2<float,Map<>,impl::Generic_tag>(0)); break;
  case  8: loop(t_maxval2<float,Map<>,impl::Generic_tag>(1)); break;
  case  9: loop(t_maxval2<float,Map<>,impl::Generic_tag>(2)); break;
  default: return 0;
  }
  return 1;
}
