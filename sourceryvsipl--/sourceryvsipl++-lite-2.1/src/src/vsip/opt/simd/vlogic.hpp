/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/vlogic.hpp
    @author  Jules Bergmann
    @date    2006-07-28
    @brief   VSIPL++ Library: SIMD element-wise vector logic operations.

*/

#ifndef VSIP_OPT_SIMD_VLOGIC_HPP
#define VSIP_OPT_SIMD_VLOGIC_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <complex>

#include <vsip/opt/simd/simd.hpp>
#include <vsip/core/metaprogramming.hpp>

#define VSIP_IMPL_INLINE_LIBSIMD 0



/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace simd
{

// Define value_types for which vband is optimized.
//  - float

template <typename T>
struct Is_algorithm_supported<T, false, Alg_vband>
{
  static bool const value = Simd_traits<T>::is_accel;
};

template <typename T>
struct Is_algorithm_supported<T, false, Alg_vbor>
{
  static bool const value = Simd_traits<T>::is_accel;
};

template <typename T>
struct Is_algorithm_supported<T, false, Alg_vbxor>
{
  static bool const value = Simd_traits<T>::is_accel;
};

template <typename T>
struct Is_algorithm_supported<T, false, Alg_vbnot>
{
  static bool const value = Simd_traits<T>::is_accel;
};

template <>
struct Is_algorithm_supported<bool, false, Alg_vland>
{
  static bool const value = Simd_traits<signed char>::is_accel;
};

template <>
struct Is_algorithm_supported<bool, false, Alg_vlor>
{
  static bool const value = Simd_traits<signed char>::is_accel;
};

template <>
struct Is_algorithm_supported<bool, false, Alg_vlxor>
{
  static bool const value = Simd_traits<signed char>::is_accel;
};

template <>
struct Is_algorithm_supported<bool, false, Alg_vlnot>
{
  static bool const value = Simd_traits<signed char>::is_accel;
};





// bitwise-and operation

struct Fun_vband
{
  template <typename T>
  static T exec(T const& A, T const& B)
  { return A & B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::band(A, B); }
};



// bitwise-or operation

struct Fun_vbor
{
  template <typename T>
  static T exec(T const& A, T const& B)
  { return A | B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::bor(A, B); }
};



// bitwise-xor operation

struct Fun_vbxor
{
  template <typename T>
  static T exec(T const& A, T const& B)
  { return A ^ B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::bxor(A, B); }
};



// bitwise-not operation

struct Fun_vbnot
{
  template <typename T>
  static T exec(T const& A)
  { return ~A; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A)
  { return SimdTraits::bnot(A); }
};



// logical-and operation

struct Fun_vland
{
  static bool exec(bool A, bool B)
  { return A && B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::band(A, B); }
};



// logical-or operation

struct Fun_vlor
{
  static bool exec(bool A, bool B)
  { return A || B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::bor(A, B); }
};



// logical-xor operation

struct Fun_vlxor
{
  static bool exec(bool A, bool B)
  { return A ^ B; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A, SimdValueT const& B)
  { return SimdTraits::bxor(A, B); }
};



// logical-not operation

struct Fun_vlnot
{
  static bool exec(bool A)
  { return !A; }

  template <typename SimdTraits, typename SimdValueT>
  static SimdValueT exec_simd(SimdValueT const& A)
  { return SimdTraits::bxor(A, SimdTraits::load_scalar_all(0x01)); }
};



/***********************************************************************
  Definitions -- Generic unary logical operations (boolean and bitwise)
***********************************************************************/

// General class template for unary logical operations

template <typename T,
	  typename SimdValueT,
	  bool     Is_vectorized,
	  typename FunctionT>
struct Simd_vlogic_unary;



// Generic, non-vectorized implementation of logical operations.

template <typename T,
	  typename SimdValueT,
          typename FunctionT>
struct Simd_vlogic_unary<T, SimdValueT, false, FunctionT>
{
  static void exec(T const* A, T* R, int n)
  {
    while (n)
    {
      *R = FunctionT::exec(*A);
      R++; A++;
      n--;
    }
  }
};



// Vectorized implementation of logical operations.

// Works under the following combinations:
//  - Fun_bnot: linux ia32   sse     GCC 3.4, T=int  (060728)
//  - Fun_lnot: linux ia32   sse     GCC 3.4, T=bool (060730)

template <typename T,
	  typename SimdValueT,
          typename FunctionT>
struct Simd_vlogic_unary<T, SimdValueT, true, FunctionT>
{
  static void exec(T const* A, T* R, int n)
  {
    typedef vsip::impl::simd::Simd_traits<SimdValueT> traits;
    typedef typename traits::simd_type                simd_type;

    // handle mis-aligned vectors
    if (   traits::alignment_of((SimdValueT*)R) !=
	   traits::alignment_of((SimdValueT*)A))
    {
      // PROFILE
      while (n)
      {
	*R = FunctionT::exec(*A);
	R++; A++;
	n--;
      }
      return;
    }

    // clean up initial unaligned values
    while (n && traits::alignment_of((SimdValueT*)A) != 0)
    {
      *R = FunctionT::exec(*A);
      R++; A++;
      n--;
    }
  
    if (n == 0) return;

    traits::enter();

    int const unroll = 1;
    while (n >= unroll*traits::vec_size)
    {
      n -= unroll*traits::vec_size;

      simd_type regA0 = traits::load((SimdValueT*)A);
      simd_type res   = FunctionT::template exec_simd<traits>(regA0);
      traits::store((SimdValueT*)R, res);
      
      A += unroll*traits::vec_size;
      R += unroll*traits::vec_size;
    }
    
    traits::exit();

    while (n)
    {
      *R = FunctionT::exec(*A);
      R++; A++;
      n--;
    }
  }
};



/***********************************************************************
  Definitions -- Generic binary logical operations (boolean and bitwise)
***********************************************************************/

// General class template for binary logical operations

template <typename T,
	  typename SimdValueT,
	  bool     Is_vectorized,
	  typename FunctionT>
struct Simd_vlogic_binary;



// Generic, non-vectorized implementation of logical operations.

template <typename T,
	  typename SimdValueT,
          typename FunctionT>
struct Simd_vlogic_binary<T, SimdValueT, false, FunctionT>
{
  static void exec(T const* A, T const* B, T* R, int n)
  {
    while (n)
    {
      *R = FunctionT::exec(*A, *B);
      R++; A++; B++;
      n--;
    }
  }
};



// Vectorized implementation of logical operations.

// Works under the following combinations:
//  - Fun_band: linux ia32   sse     GCC 3.4, T=int  (060728)
//  - Fun_bor : linux ia32   sse     GCC 3.4, T=int  (060728)
//  - Fun_bxor: linux ia32   sse     GCC 3.4, T=int  (060728) ?????
//  - Fun_land: linux ia32   sse     GCC 3.4, T=bool (060730)
//  - Fun_lor : linux ia32   sse     GCC 3.4, T=bool (060730)
//  - Fun_lxor: linux ia32   sse     GCC 3.4, T=bool (060730) ?????

template <typename T,
	  typename SimdValueT,
          typename FunctionT>
struct Simd_vlogic_binary<T, SimdValueT, true, FunctionT>
{
  static void exec(T const* A, T const* B, T* R, int n)
  {
    typedef vsip::impl::simd::Simd_traits<SimdValueT> traits;
    typedef typename traits::simd_type                simd_type;

    // handle mis-aligned vectors
    if (   traits::alignment_of((SimdValueT*)R) !=
	   traits::alignment_of((SimdValueT*)A)
	|| traits::alignment_of((SimdValueT*)R) !=
	   traits::alignment_of((SimdValueT*)B))
    {
      // PROFILE
      while (n)
      {
	*R = FunctionT::exec(*A, *B);
	R++; A++; B++;
	n--;
      }
      return;
    }

    // clean up initial unaligned values
    while (n && traits::alignment_of((SimdValueT*)A) != 0)
    {
      *R = FunctionT::exec(*A, *B);
      R++; A++; B++;
      n--;
    }
  
    if (n == 0) return;

    traits::enter();

    int const unroll = 1;
    while (n >= unroll*traits::vec_size)
    {
      n -= unroll*traits::vec_size;

      simd_type regA0 = traits::load((SimdValueT*)A);
      simd_type regB0 = traits::load((SimdValueT*)B);
      simd_type res   = FunctionT::template exec_simd<traits>(regA0, regB0);
      traits::store((SimdValueT*)R, res);
      
      A += unroll*traits::vec_size;
      B += unroll*traits::vec_size;
      R += unroll*traits::vec_size;
    }
    
    traits::exit();

    while (n)
    {
      *R = FunctionT::exec(*A, *B);
      R++; A++; B++;
      n--;
    }
  }
};



/***********************************************************************
  Definitions -- Specific gateway functions.
***********************************************************************/

// Depending on VSIP_IMPL_LIBSIMD_INLINE macro, either provide these
// functions inline, or provide non-inline functions in the libvsip.a.

#if VSIP_IMPL_INLINE_LIBSIMD
template <typename T>
inline void
vband(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vband>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vband>::exec(op1, op2, res, size);
}

template <typename T>
inline void
vbor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbor>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vbor>::exec(op1, op2, res, size);
}

template <typename T>
inline void
vbxor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbxor>::value;
  Simd_vlogic_binary<T, T, Is_vectorized, Fun_vbxor>::exec(op1, op2, res, size);
}

template <typename T>
inline void
vbnot(
  T const* op1,
  T*       res,
  int      size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_vbnot>::value;
  Simd_vlogic_unary<T, T, Is_vectorized, Fun_vbnot>::exec(op1, res, size);
}

inline void
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

inline void
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

inline void
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

inline void
vlnot(
  bool const* op1,
  bool*       res,
  int         size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<bool, false, Alg_vlnot>::value;
  Simd_vlogic_binary<bool, signed char, Is_vectorized, Fun_vlnot>::
    exec(op1, res, size);
}
#else
template <typename T>
void
vband(
  T const* op1,
  T const* op2,
  T*       res,
  int      size);

template <typename T>
void
vbor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size);

template <typename T>
void
vbxor(
  T const* op1,
  T const* op2,
  T*       res,
  int      size);

template <typename T>
void
vbnot(
  T const* op1,
  T*       res,
  int      size);

void
vland(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int         size);

void
vlor(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int         size);

void
vlxor(
  bool const* op1,
  bool const* op2,
  bool*       res,
  int         size);

void
vlnot(
  bool const* op1,
  bool*       res,
  int         size);
#endif // VSIP_IMPL_INLINE_LIBSIMD


} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_SIMD_VLOGIC_HPP
