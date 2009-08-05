/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/vlogic.cpp
    @author  Jules Bergmann
    @date    2006-07-28
    @brief   VSIPL++ Library: SIMD element-wise vector operations.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/simd/vlogic.hpp>



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
vband(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vband>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vband>::
    exec(op1, op2, res, size);
}

template void vband(int const*, int const*, int*, int);

template <typename T>
void
vbor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbor>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vbor>::
    exec(op1, op2, res, size);
}

template void vbor(int const*, int const*, int*, int);

template <typename T>
void
vbxor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbxor>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vbxor>::
    exec(op1, op2, res, size);
}

template void vbxor(int const*, int const*, int*, int);

template <typename T>
void
vbnot(
  T const* op1,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbnot>::value;
  Simd_vlogic_unary<T, T, Is_vectorized, Fun_vbnot>::exec(op1, res, size);
}

template void vbnot(int const*, int*, int);

void
vland(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int         size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<bool, false, Alg_vland>::value;
  Simd_vlogic_binary<bool, signed char, Is_vectorized, Fun_vland>::
    exec(op1, op2, res, size);
}

void
vlor(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<bool, false, Alg_vlor>::value;
  Simd_vlogic_binary<bool, signed char, Is_vectorized, Fun_vlor>::
    exec(op1, op2, res, size);
}

void
vlxor(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int         size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<bool, false, Alg_vlxor>::value;
  Simd_vlogic_binary<bool, signed char, Is_vectorized, Fun_vlxor>::
    exec(op1, op2, res, size);
}

void
vlnot(
  bool const* op1,
  bool*       res,
  int         size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<bool, false, Alg_vlnot>::value;
  Simd_vlogic_unary<bool, signed char, Is_vectorized, Fun_vlnot>::
    exec(op1, res, size);
}

#endif

} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip
