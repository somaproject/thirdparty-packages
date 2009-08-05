/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/vgt.cpp
    @author  Jules Bergmann
    @date    2006-07-26
    @brief   VSIPL++ Library: SIMD element-wise vector greater-than.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/simd/vgt.hpp>



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
vgt(
  T const* op1,
  T const* op2,
  bool*    res,
  int      size)
{
  static bool const Is_vectorized = Is_algorithm_supported<T, false, Alg_vgt>
                                      ::value;
  Simd_vgt<T, Is_vectorized>::exec(op1, op2, res, size);
}

template void vgt(float const*, float const*, bool*, int);



#endif

} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip
