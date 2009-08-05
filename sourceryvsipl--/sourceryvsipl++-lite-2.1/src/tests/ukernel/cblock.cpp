/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    tests/ukernels/cblock.cpp
    @author  Jules Bergmann
    @date    2008-12-12
    @brief   VSIPL++ Library: User-defined kernel example illustrating
                              a control-block approach.
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
#include <vsip/opt/ukernel/kernels/host/cblock.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/output.hpp>

using namespace std;
using namespace vsip;
using vsip_csl::equal;


/***********************************************************************
  Definitions
***********************************************************************/

// Applies cblock ukernel to perform multiply add:
//   out = in0 * in1 + in2;

template <typename View0,
	  typename View1,
	  typename View2,
	  typename View3>
void
apply_ukernel(
  View0       in0,
  View1       in1,
  View2       in2,
  View3       out,
  length_type req_accel)
{
  impl::Ext_data<typename View0::block_type> ext0(in0.block());
  impl::Ext_data<typename View1::block_type> ext1(in1.block());
  impl::Ext_data<typename View2::block_type> ext2(in2.block());
  impl::Ext_data<typename View3::block_type> ext3(out.block());

  test_assert(ext0.stride(1) == 1);
  test_assert(ext1.stride(1) == 1);
  test_assert(ext2.stride(1) == 1);
  test_assert(ext3.stride(1) == 1);

  Cblock_kernel obj(
    (uintptr_t)ext0.data(), ext0.stride(0),
    (uintptr_t)ext1.data(), ext1.stride(0),
    (uintptr_t)ext2.data(), ext2.stride(0),
    (uintptr_t)ext3.data(), ext3.stride(0),
    out.size(0), out.size(1),
    req_accel);

  vsip_csl::ukernel::Ukernel<Cblock_kernel> cblock_uk(obj);


  cblock_uk();
}



template <typename T>
void
test_ukernel(length_type rows, length_type cols, length_type req_accel)
{
  Matrix<T> in0(rows, cols);
  Matrix<T> in1(rows, cols);
  Matrix<T> in2(rows, cols);
  Matrix<T> out(rows, cols);

  Rand<T> gen1(0, 0);
  in0 = gen1.randu(rows, cols);

  Rand<T> gen2(1, 0);
  in1 = gen2.randu(rows, cols);
  in2 = gen2.randu(rows, cols);

  apply_ukernel(in0, in1, in2, out, req_accel);

  for (index_type i=0; i < rows; ++i)
    for (index_type j=0; j < cols; ++j)
    {
      T madd = in0.get(i, j) * in1.get(i, j) + in2.get(i, j);
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
  test_ukernel<float>(63, 1024, 0);
  test_ukernel<float>(64, 1024, 1);
  test_ukernel<float>(64, 1024, 2);
  test_ukernel<float>(64, 1024, 0);

  return 0;
}
