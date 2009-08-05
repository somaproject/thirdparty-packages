/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/threshold.cpp
    @author  Assem Salama
    @date    2007-05-02
    @brief   VSIPL++ Library: SIMD threshold

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/simd/threshold.hpp>



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

template <template <typename,typename> class O,
          typename T>
void 
threshold(T* Z, T* A, T* B, T k, int n)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_threshold>::value &&
    Binary_operator_map<T,O>::is_supported;
  Simd_threshold<T, O, Is_vectorized>::exec(Z,A,B,k,n);
}

#define VSIP_OPT_DECL_THRESH(O) \
template void threshold<O>(float* Z, float* A, float* B, float k, int n); \
template void threshold<O>(double* Z, double* A, double* B, double k, int n);

VSIP_OPT_DECL_THRESH(gt_functor)
VSIP_OPT_DECL_THRESH(lt_functor)
VSIP_OPT_DECL_THRESH(ge_functor)
VSIP_OPT_DECL_THRESH(le_functor)

#undef VSIP_OPT_DECL_THRESH

#endif

}
}
}
