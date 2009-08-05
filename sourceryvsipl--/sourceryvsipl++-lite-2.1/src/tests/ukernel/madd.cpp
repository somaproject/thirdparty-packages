/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernels/madd.cpp
    @author  Don McCoy
    @date    2008-08-26
    @brief   VSIPL++ Library: User-defined kernel for multiply-add.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/signal.hpp>
#include <vsip/random.hpp>
#include <vsip_csl/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/host/madd.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using vsip_csl::equal;


/***********************************************************************
  Definitions
***********************************************************************/

// Performs a elementwise multiply-add for the expression
//   Z = A * B + C
// where T1 is the type of A and T2 is the type for B, C and D
//
template <typename T1,
          typename T2>
void
test_ukernel(length_type rows, length_type cols)
{
  Madd_kernel obj;

  vsip_csl::ukernel::Ukernel<Madd_kernel> madd_uk(obj);

  Matrix<T1> in0(rows, cols);
  Matrix<T2> in1(rows, cols);
  Matrix<T2> in2(rows, cols);
  Matrix<T2> out(rows, cols);

  Rand<T1> gen1(0, 0);
  in0 = gen1.randu(rows, cols);

  Rand<T2> gen2(1, 0);
  in1 = gen2.randu(rows, cols);
  in2 = gen2.randu(rows, cols);

  madd_uk(in0, in1, in2, out);


  for (index_type i=0; i < rows; ++i)
    for (index_type j=0; j < cols; ++j)
    {
      T2 madd = in0.get(i, j) * in1.get(i, j) + in2.get(i, j);
      if (!equal(madd, out.get(i, j)))
      {
        std::cout << "index " << i << ", " << j << " : "
                  << in0.get(i, j) << " * "
                  << in1.get(i, j) << " + "
                  << in2.get(i, j) << " = "
                  << in0.get(i, j) * in1.get(i, j) + in2.get(i, j) << "  vs  "
                  << out.get(i, j)
                  << std::endl;
      }
      test_assert(equal(
          in0.get(i, j) * in1.get(i, j) + in2.get(i, j), 
          out.get(i, j)));
    }
}


/***********************************************************************
  Main
***********************************************************************/

int
main(int argc, char** argv)
{
  vsipl init(argc, argv);

  // Parameters are rows then cols
  test_ukernel<float, float>(64, 1024);

// This kernel is presently only implemented for interleaved complex
#if !VSIP_IMPL_PREFER_SPLIT_COMPLEX
  test_ukernel<float, complex<float> >(64, 1024);
  test_ukernel<complex<float>, complex<float> >(64, 1024);
#endif

  return 0;
}
