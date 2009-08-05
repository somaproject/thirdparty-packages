/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/rscvmul.hpp
    @author  Jules Bergmann
    @date    2006-03-28
    @brief   VSIPL++ Library: SIMD element-wise vector multiplication.

*/

#ifndef VSIP_OPT_SIMD_RSCVMUL_HPP
#define VSIP_OPT_SIMD_RSCVMUL_HPP

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

// Define value_types for which vmul is optimized.
//  - complex<float>
//  - complex<double>

template <typename T,
	  bool     IsSplit>
struct Is_algorithm_supported<T, IsSplit, Alg_rscvmul>
{
  static bool const value =
    Simd_traits<T>::is_accel &&
    (Type_equal<T, float>::value ||
     Type_equal<T, double>::value);
};



// Class for vmul - vector element-wise multiplication.

template <typename T,
	  bool     Is_vectorized>
struct Simd_rscvmul;



// Generic, non-vectorized implementation of vector element-wise multiply.

template <typename T>
struct Simd_rscvmul<std::complex<T>, false>
{
  static void exec(T alpha, std::complex<T>* B, std::complex<T>* R, int n)
  {
    while (n)
    {
      *R = alpha * *B;
      R++; B++;
      n--;
    }
  }
};



// Vectorized implementation of vector element-wise multiply for
// interleaved complex (complex<float>, complex<double>, etc).

template <typename T>
struct Simd_rscvmul<std::complex<T>, true>
{
  static void exec(
    T                alpha,
    std::complex<T>* B,
    std::complex<T>* R,
    int n)
  {
    typedef Simd_traits<T> simd;
    typedef typename simd::simd_type simd_type;
    
    // handle mis-aligned vectors
    if (simd::alignment_of((T*)R) != simd::alignment_of((T*)B))
    {
      // PROFILE
      while (n)
      {
	*R = alpha * *B;
	R++; B++;
	n--;
      }
      return;
    }

    // clean up initial unaligned values
    while (n && simd::alignment_of((T*)B) != 0)
    {
      *R = alpha * *B;
      R++; B++;
      n--;
    }
  
    if (n == 0) return;

    simd::enter();

    simd_type regA = simd::load_scalar_all(alpha);

    while (n >= simd::vec_size)
    {
      n -= simd::vec_size;

      simd_type regB1 = simd::load((T*)B);
      simd_type regB2 = simd::load((T*)B + simd::vec_size);

      simd_type regR1 = simd::mul(regA, regB1);
      simd_type regR2 = simd::mul(regA, regB2);
      
      simd::store((T*)R,                  regR1);
      simd::store((T*)R + simd::vec_size, regR2);
      
      B+=simd::vec_size; R+=simd::vec_size;
    }

    simd::exit();

    while (n)
    {
      *R = alpha * *B;
      R++; B++;
      n--;
    }
  }
};



// Generic, non-vectorized implementation of vector element-wise multiply for
// split complex (as represented by pair<float*, float*>, etc).

template <typename T>
struct Simd_rscvmul<std::pair<T, T>, false>
{
  static void exec(
    T                        alpha,
    std::pair<T*, T*> const& B,
    std::pair<T*, T*> const& R,
    int n)
  {
    T const* pBr = B.first;
    T const* pBi = B.second;

    T* pRr = R.first;
    T* pRi = R.second;

    while (n)
    {
      *pRr = alpha * *pBr;
      *pRi = alpha * *pBi;
      pRr++; pRi++;
      pBr++; pBi++;
      n--;
    }
  }
};



// Vectorized implementation of vector element-wise multiply for
// split complex (as represented by pair<float*, float*>, etc).

template <typename T>
struct Simd_rscvmul<std::pair<T, T>, true>
{
  static void exec(
    T                        alpha,
    std::pair<T*, T*> const& B,
    std::pair<T*, T*> const& R,
    int                      n)
  {
    assert(n >= 0);

    typedef Simd_traits<T> simd;
    typedef typename simd::simd_type simd_type;
    
    T const* pBr = B.first;
    T const* pBi = B.second;

    T* pRr = R.first;
    T* pRi = R.second;

    // handle mis-aligned vectors
    if (simd::alignment_of(pRr) != simd::alignment_of(pRi) ||
	simd::alignment_of(pRr) != simd::alignment_of(pBr) ||
	simd::alignment_of(pRr) != simd::alignment_of(pBi))
    {
      // PROFILE
      while (n)
      {
	*pRr = alpha * *pBr;
	*pRi = alpha * *pBi;
	pRr++; pRi++;
	pBr++; pBi++;
	n--;
      }
      return;
    }

    // clean up initial unaligned values
    while (n && simd::alignment_of(pRr) != 0)
    {
      *pRr = alpha * *pBr;
      *pRi = alpha * *pBi;
      pRr++; pRi++;
      pBr++; pBi++;
      n--;
    }
  
    if (n == 0) return;

    simd::enter();

    simd_type regA = simd::load_scalar_all(alpha);

    while (n >= simd::vec_size)
    {
      n -= simd::vec_size;

      simd_type Br = simd::load((T*)pBr);
      simd_type Bi = simd::load((T*)pBi);
      
      simd_type Rr   = simd::mul(regA, Br);
      simd_type Ri   = simd::mul(regA, Bi);

      simd::store_stream(pRr, Rr);
      simd::store_stream(pRi, Ri);
      
      pRr += simd::vec_size; pRi += simd::vec_size;
      pBr += simd::vec_size; pBi += simd::vec_size;
    }

    simd::exit();

    while (n)
    {
      *pRr = alpha * *pBr;
      *pRi = alpha * *pBi;
      pRr++; pRi++;
      pBr++; pBi++;
      n--;
    }
  }
};



// Depending on VSIP_IMPL_LIBSIMD_INLINE macro, either provide these
// functions inline, or provide non-inline functions in the libvsip.a.

#if VSIP_IMPL_INLINE_LIBSIMD

template <typename T>
inline void
rscvmul(
  T                op1,
  std::complex<T>* op2,
  std::complex<T>* res,
  int size)
{
  static bool const Is_vectorized =
    Is_algorithm_supported<T, false, Alg_rscvmul>::value;
  Simd_rscvmul<T, Is_vectorized>::exec(op1, op2, res, size);
}

template <typename T>
inline void
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

#else

template <typename T>
void
rscvmul(
  T                op1,
  std::complex<T>* op2,
  std::complex<T>* res,
  int              size);

template <typename T>
void
rscvmul(
  T                op1,
  std::pair<T*,T*> op2,
  std::pair<T*,T*> res,
  int              size);

#endif // VSIP_IMPL_INLINE_LIBSIMD


} // namespace vsip::impl::simd
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_SIMD_VMUL_HPP
