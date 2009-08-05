/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/memwrite_sal.cpp
    @author  Jules Bergmann
    @date    2006-10-12
    @brief   VSIPL++ Library: Benchmark for SAL memory write.

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
struct t_memwrite_sal;

template <typename ComplexFmt>
struct t_memwrite_sal<float, ComplexFmt> : Benchmark_base
{
  typedef float T;

  char* what() { return "t_memwrite_sal"; }
  int ops_per_point(size_t)  { return 1; }
  int riob_per_point(size_t) { return 1*sizeof(float); }
  int wiob_per_point(size_t) { return 1*sizeof(float); }
  int mem_per_point(size_t)  { return 2*sizeof(float); }

  void operator()(size_t size, size_t loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   Z(size, T(2));
    T val = T(3);

    vsip::impl::profile::Timer t1;

    {
      impl::Ext_data<block_type> ext_z(Z.block(), impl::SYNC_OUT);
    
      T* pZ = ext_z.data();
    
      t1.start();
      for (size_t l=0; l<loop; ++l)
	vfillx(&val, pZ, 1, size, 0);
      t1.stop();
    }
    
    for (index_type i=0; i<size; ++i)
      test_assert(Z.get(i) == val);
    
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
  case  1: loop(t_memwrite_sal<float>()); break;
  }
}
