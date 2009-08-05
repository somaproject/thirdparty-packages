/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/corr.cpp
    @author  Jules Bergmann
    @date    2005-10-06
    @brief   VSIPL++ Library: Benchmark for Correlation.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/signal.hpp>

#include <vsip/core/profile.hpp>

#include <vsip_csl/test.hpp>
#include "loop.hpp"

using namespace vsip;


/***********************************************************************
  Definitions
***********************************************************************/

template <support_region_type Supp,
	  typename            T>
struct t_corr1 : Benchmark_base
{
  char const* what() { return "t_corr1"; }
  float ops_per_point(length_type size)
  {
    length_type output_size = this->my_output_size(size);
    float ops = ref_size_ * output_size *
      (vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add);

    return ops / size;
  }

  int riob_per_point(length_type) { return -1*(int)sizeof(T); }
  int wiob_per_point(length_type) { return -1*(int)sizeof(T); }
  int mem_per_point(length_type) { return 2*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    length_type output_size = this->my_output_size(size);

    Vector<T>   in (size, T());
    Vector<T>   out(output_size);
    Vector<T>   ref(ref_size_, T());

    ref(0) = T(1);
    ref(1) = T(2);

    typedef Correlation<const_Vector, Supp, T> corr_type;

    corr_type corr((Domain<1>(ref_size_)), Domain<1>(size));

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
     corr(bias_, ref, in, out);
    t1.stop();
    
    time = t1.delta();
  }

  t_corr1(length_type ref_size, bias_type bias)
    : ref_size_(ref_size),
      bias_    (bias)
  {}

  
  length_type my_output_size(length_type size)
  {
    if      (Supp == support_full)
      return size + ref_size_ - 1;
    else if (Supp == support_same)
      return size;
    else /* (Supp == support_min) */
      return size - ref_size_ + 1;
  }
  

  length_type ref_size_;
  bias_type   bias_;
};



// Benchmark performance of a Correlation_impl object
// (requires ImplTag to select implementation)

template <typename            ImplTag,
	  support_region_type Supp,
	  typename            T>
struct t_corr2 : Benchmark_base
{
  char const* what() { return "t_corr2"; }
  float ops_per_point(length_type size)
  {
    length_type output_size = this->my_output_size(size);
    float ops = ref_size_ * output_size *
      (vsip::impl::Ops_info<T>::mul + vsip::impl::Ops_info<T>::add);

    return ops / size;
  }

  int riob_per_point(length_type) { return -1*(int)sizeof(T); }
  int wiob_per_point(length_type) { return -1*(int)sizeof(T); }
  int mem_per_point(length_type) { return 2*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    length_type output_size = this->my_output_size(size);

    Vector<T>   in (size, T());
    Vector<T>   out(output_size);
    Vector<T>   ref(ref_size_, T());

    ref(0) = T(1);
    ref(1) = T(2);

    typedef typename vsip::impl::dispatcher::Evaluator<
      vsip::impl::dispatcher::Corr_tag<1, Supp, T, 0, alg_time>,
      ImplTag>::backend_type
      corr_type;

    corr_type corr((Domain<1>(ref_size_)), Domain<1>(size));

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
     corr.impl_correlate(bias_, ref, in, out);
    t1.stop();
    
    time = t1.delta();
  }

  t_corr2(length_type ref_size, bias_type bias)
    : ref_size_(ref_size),
      bias_    (bias)
  {}

  
  length_type my_output_size(length_type size)
  {
    if      (Supp == support_full)
      return size + ref_size_ - 1;
    else if (Supp == support_same)
      return size;
    else /* (Supp == support_min) */
      return size - ref_size_ + 1;
  }
  

  length_type ref_size_;
  bias_type   bias_;
};



void
defaults(Loop1P& loop)
{
  loop.loop_start_ = 5000;
  loop.start_ = 4;
  loop.user_param_ = 16;
}



int
test(Loop1P& loop, int what)
{
  length_type M = loop.user_param_;
  using vsip::impl::Opt_tag;
  using vsip::impl::Generic_tag;

  typedef float T;
  typedef complex<float> CT;

  switch (what)
  {
  case  1: loop(t_corr1<support_full, T>(M, biased)); break;
  case  2: loop(t_corr1<support_same, T>(M, biased)); break;
  case  3: loop(t_corr1<support_min,  T>(M, biased)); break;

  case  4: loop(t_corr1<support_full, T>(M, unbiased)); break;
  case  5: loop(t_corr1<support_same, T>(M, unbiased)); break;
  case  6: loop(t_corr1<support_min,  T>(M, unbiased)); break;

  case  7: loop(t_corr1<support_full, CT >(M, biased)); break;
  case  8: loop(t_corr1<support_same, CT >(M, biased)); break;
  case  9: loop(t_corr1<support_min,  CT >(M, biased)); break;

  case  10: loop(t_corr1<support_full, CT >(M, unbiased)); break;
  case  11: loop(t_corr1<support_same, CT >(M, unbiased)); break;
  case  12: loop(t_corr1<support_min,  CT >(M, unbiased)); break;

  case  13: loop(t_corr2<Opt_tag, support_full, T>(M, biased)); break;
  case  14: loop(t_corr2<Opt_tag, support_same, T>(M, biased)); break;
  case  15: loop(t_corr2<Opt_tag, support_min,  T>(M, biased)); break;

  case  16: loop(t_corr2<Opt_tag, support_full, T>(M, unbiased)); break;
  case  17: loop(t_corr2<Opt_tag, support_same, T>(M, unbiased)); break;
  case  18: loop(t_corr2<Opt_tag, support_min,  T>(M, unbiased)); break;

  case  19: loop(t_corr2<Opt_tag, support_full, CT >(M, biased)); break;
  case  20: loop(t_corr2<Opt_tag, support_same, CT >(M, biased)); break;
  case  21: loop(t_corr2<Opt_tag, support_min,  CT >(M, biased)); break;

  case  22: loop(t_corr2<Opt_tag, support_full, CT >(M, unbiased)); break;
  case  23: loop(t_corr2<Opt_tag, support_same, CT >(M, unbiased)); break;
  case  24: loop(t_corr2<Opt_tag, support_min,  CT >(M, unbiased)); break;

  case  25: loop(t_corr2<Generic_tag, support_full, T>(M, biased)); break;
  case  26: loop(t_corr2<Generic_tag, support_same, T>(M, biased)); break;
  case  27: loop(t_corr2<Generic_tag, support_min,  T>(M, biased)); break;

  case  28: loop(t_corr2<Generic_tag, support_full, T>(M, unbiased)); break;
  case  29: loop(t_corr2<Generic_tag, support_same, T>(M, unbiased)); break;
  case  30: loop(t_corr2<Generic_tag, support_min,  T>(M, unbiased)); break;

  case  31: loop(t_corr2<Generic_tag, support_full, CT >(M, biased)); break;
  case  32: loop(t_corr2<Generic_tag, support_same, CT >(M, biased)); break;
  case  33: loop(t_corr2<Generic_tag, support_min,  CT >(M, biased)); break;

  case  34: loop(t_corr2<Generic_tag, support_full, CT >(M, unbiased)); break;
  case  35: loop(t_corr2<Generic_tag, support_same, CT >(M, unbiased)); break;
  case  36: loop(t_corr2<Generic_tag, support_min,  CT >(M, unbiased)); break;


  default: return 0;
  }
  return 1;
}
