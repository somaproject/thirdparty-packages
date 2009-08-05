/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/rscvmul.cpp
    @author  Jules Bergmann
    @date    2006-03-28
    @brief   VSIPL++ Library: SIMD element-wise vector multiplication.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/simd/rscvmul.hpp>



/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace simd
{

#if !VSIP_IMPL_INLINE_LIBSIMD

template <typename T>
void
rscvmul(
  T                op1,
  std::complex<T>* op2,
  std::complex<T>* res,
  int size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_rscvmul>::value;
  Simd_rscvmul<std::complex<T>, Is_vectorized>::exec(op1, op2, res, size);
}

template void rscvmul(float, std::complex<float>*, std::complex<float>*, int);
template void rscvmul(double, std::complex<double>*,
		      std::complex<double>*, int);



template <typename T>
void
rscvmul(
  T                op1,
  std::pair<T*,T*> op2,
  std::pair<T*,T*> res,
  int              size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, true, Alg_rscvmul>::value;
  Simd_rscvmul<std::pair<T,T>, Is_vectorized>::exec(op1, op2, res, size);
}

template void rscvmul(float,
		   std::pair<float*,float*>,
		   std::pair<float*,float*>, int);
template void rscvmul(double,
		   std::pair<double*,double*>,
		   std::pair<double*,double*>, int);

#endif

} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip
