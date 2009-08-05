/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/vadd.cpp
    @author  Jules Bergmann
    @date    2006-06-08
    @brief   VSIPL++ Library: SIMD element-wise vector multiplication.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/simd/vadd.hpp>



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
vadd(
  T*  op1,
  T*  op2,
  T*  res,
  int size)
{
  static bool const Is_vectorized = Is_algorithm_supported<T, false, Alg_vadd>
                                      ::value;
  Simd_vadd<T, Is_vectorized>::exec(op1, op2, res, size);
}

template void vadd(short*, short*, short*, int);
template void vadd(float*, float*, float*, int);
template void vadd(double*, double*, double*, int);
template void vadd(std::complex<float>*, std::complex<float>*,
		   std::complex<float>*, int);
template void vadd(std::complex<double>*, std::complex<double>*,
		   std::complex<double>*, int);



template <typename T>
void
vadd(
  std::pair<T*,T*>  op1,
  std::pair<T*,T*>  op2,
  std::pair<T*,T*>  res,
  int size)
{
  static bool const Is_vectorized = Is_algorithm_supported<T, true, Alg_vadd>
                                      ::value;
  Simd_vadd<std::pair<T,T>, Is_vectorized>::exec(op1, op2, res, size);
}

template void vadd(std::pair<float*,float*>,
		   std::pair<float*,float*>,
		   std::pair<float*,float*>, int);
template void vadd(std::pair<double*,double*>,
		   std::pair<double*,double*>,
		   std::pair<double*,double*>, int);

#endif

} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip
