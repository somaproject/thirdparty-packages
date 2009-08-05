/* Copyright (c) 2007 by CodeSourcery, LLC.  All rights reserved. */

/** @file    tests/vmul_sizes.cpp
    @author  Jules Bergmann
    @date    2007-03-16
    @brief   VSIPL++ Library: Check that range of vmul sizes are handled.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#define VERBOSE 0

#if VERBOSE
#  include <iostream>
#endif

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/domain.hpp>
#include <vsip/random.hpp>

#include <vsip_csl/test.hpp>

using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Definitions - Utility Functions
***********************************************************************/

template <typename T,
	  typename ComplexFmt>
void
test_vmul(length_type len)
{
  typedef impl::Layout<1, row1_type, impl::Stride_unit_dense, ComplexFmt>
		LP;
  typedef impl::Fast_block<1, T, LP> block_type;

  Rand<T> gen(0, 0);

  Vector<T, block_type> A(len);
  Vector<T, block_type> B(len);
  Vector<T, block_type> Z(len);

  A = gen.randu(len);
  B = gen.randu(len);

  Z = A * B;

  for (index_type i=0; i<len; ++i)
  {
#if VERBOSE
    if (!equal(Z.get(i), A.get(i) * B.get(i)))
    {
      std::cout << "Z(" << i << ")        = " << Z(i) << std::endl;
      std::cout << "A(" << i << ") * B(" << i << ") = "
		<< A(i) * B(i) << std::endl;
    }
#endif
    test_assert(Almost_equal<T>::eq(Z.get(i), A.get(i) * B.get(i)));
  }
}




template <typename T,
	  typename ComplexFmt>
void
test_sweep()
{
  for (index_type i=1; i<=128; ++i)
    test_vmul<T, ComplexFmt>(i);
}




int
main(int argc, char** argv)
{
  typedef impl::Cmplx_inter_fmt cif;
  typedef impl::Cmplx_split_fmt csf;

  vsipl init(argc, argv);

  test_sweep<float,          impl::Cmplx_inter_fmt>();
  test_sweep<complex<float>, impl::Cmplx_inter_fmt>();
  test_sweep<complex<float>, impl::Cmplx_split_fmt>();
}
