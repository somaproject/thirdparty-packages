/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/cholesky.hpp
    @author  Assem Salama
    @date    2006-04-13
    @brief   VSIPL++ Library: Cholesky linear system solver using SAL.

*/

#ifndef VSIP_IMPL_SAL_SOLVER_CHOLESKY_HPP
#define VSIP_IMPL_SAL_SOLVER_CHOLESKY_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/working_view.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/solver/common.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

// A structure that tells us if sal chold impl is available
// for certain types
template <>
struct Is_chold_impl_avail<Mercury_sal_tag, float>
{
  static bool const value = true;
};

template <>
struct Is_chold_impl_avail<Mercury_sal_tag, complex<float> >
{
  static bool const value = true;
};

// SAL Cholesky decomposition
#define VSIP_IMPL_SAL_CHOL_DEC( T, D_T, SAL_T, SALFCN ) \
inline void                          \
sal_mat_chol_dec(                    \
  T *a,                              \
  D_T *d, int n)                     \
{                                    \
  SALFCN(                            \
   (SAL_T*) a, n,                    \
   (SAL_T*) a, n,                    \
   (D_T*) d, n);                     \
}

#define VSIP_IMPL_SAL_CHOL_DEC_SPLIT( T, D_T, SAL_T, SALFCN ) \
inline void                          \
sal_mat_chol_dec(                    \
  std::pair<T*,T*> a,                \
  D_T *d, int n)                     \
{                                    \
  SALFCN(                            \
   (SAL_T*) &a, n,                   \
   (SAL_T*) &a, n,                   \
   (D_T*) d, n);                     \
}


VSIP_IMPL_SAL_CHOL_DEC(float,                   float,float,        matchold)
VSIP_IMPL_SAL_CHOL_DEC(complex<float>,          float,COMPLEX,      cmatchold)
VSIP_IMPL_SAL_CHOL_DEC_SPLIT(float,             float,COMPLEX_SPLIT,zmatchold)

// SAL Cholesky solver
#define VSIP_IMPL_SAL_CHOL_SOL( T, D_T, SAL_T, SALFCN ) \
inline void                          \
sal_mat_chol_sol(                    \
  T *a, int atcols,                  \
  D_T *d, T *b, T *x, int n)         \
{                                    \
  SALFCN(                            \
   (SAL_T*) a, atcols,               \
   (D_T*) d, (SAL_T*)b, (SAL_T*)x,   \
   n);                               \
}

#define VSIP_IMPL_SAL_CHOL_SOL_SPLIT( T, D_T, SAL_T, SALFCN ) \
inline void                                        \
sal_mat_chol_sol(                                  \
  std::pair<T*,T*> a, int atcols,                  \
  D_T *d, std::pair<T*,T*> b,                      \
  std::pair<T*,T*> x, int n)                       \
{                                                  \
  SALFCN(                                          \
   (SAL_T*) &a, atcols,                            \
   (D_T*) d, (SAL_T*)&b, (SAL_T*)&x,               \
   n);                                             \
}

VSIP_IMPL_SAL_CHOL_SOL(float,         float,float,        matchols)
VSIP_IMPL_SAL_CHOL_SOL(complex<float>,float,COMPLEX,      cmatchols)
VSIP_IMPL_SAL_CHOL_SOL_SPLIT(float,   float,COMPLEX_SPLIT,zmatchols)

/// Cholesky factorization implementation class.  Common functionality
/// for chold by-value and by-reference classes.

template <typename T>
class Chold_impl<T,Mercury_sal_tag>
{
  // The matrix to be decomposed using SAL must be in ROW major format. The
  // other matrix B will be in COL major format so that we can pass each
  // column to the solver. SAL supports both split and interleaved format.
  typedef vsip::impl::dense_complex_type   complex_type;
  typedef Storage<complex_type, T>         storage_type;
  typedef typename storage_type::type      ptr_type;

  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  typedef Layout<2, col2_type, Stride_unit_dense, complex_type> b_data_LP;
  typedef Fast_block<2, T, b_data_LP> b_data_block_type;

  // Constructors, copies, assignments, and destructors.
public:
  Chold_impl(mat_uplo, length_type)
    VSIP_THROW((std::bad_alloc));
  Chold_impl(Chold_impl const&)
    VSIP_THROW((std::bad_alloc));

  Chold_impl& operator=(Chold_impl const&) VSIP_NOTHROW;
  ~Chold_impl() VSIP_NOTHROW;

  // Accessors.
public:
  mat_uplo    uplo()  const VSIP_NOTHROW { return uplo_; }
  length_type length()const VSIP_NOTHROW { return length_; }

  // Solve systems.
public:
  template <typename Block>
  bool decompose(Matrix<T, Block>) VSIP_NOTHROW;

protected:
  template <typename Block0,
	    typename Block1>
  bool impl_solve(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  // Member data.
private:
  typedef std::vector<float, Aligned_allocator<float> > vector_type;

  mat_uplo     uplo_;			// A upper/lower triangular
  length_type  length_;			// Order of A.
  vector_type  idv_;			// Daignal vector from decompose

  Matrix<T, data_block_type> data_;	// Factorized Cholesky matrix (A)
};

/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
Chold_impl<T,Mercury_sal_tag>::Chold_impl(
  mat_uplo    uplo,
  length_type length
  )
VSIP_THROW((std::bad_alloc))
  : uplo_   (uplo),
    length_ (length),
    idv_    (length_),
    data_   (length_, length_)
{
  assert(length_ > 0);
}



template <typename T>
Chold_impl<T,Mercury_sal_tag>::Chold_impl(Chold_impl const& qr)
VSIP_THROW((std::bad_alloc))
  : uplo_       (qr.uplo_),
    length_     (qr.length_),
    idv_        (length_),
    data_       (length_, length_)
{
  data_ = qr.data_;
}



template <typename T>
Chold_impl<T,Mercury_sal_tag>::~Chold_impl()
  VSIP_NOTHROW
{
}



/// Form Cholesky factorization of matrix A
///
/// Requires
///   A to be a square matrix, either
///     symmetric positive definite (T real), or
///     hermitian positive definite (T complex).
///
/// FLOPS:
///   real   : (1/3) n^3
///   complex: (4/3) n^3

template <typename T>
template <typename Block>
bool
Chold_impl<T,Mercury_sal_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == length_ && m.size(1) == length_);

  data_ = m;
  Ext_data<data_block_type> ext(data_.block());

  if(length_ > 1) 
    sal_mat_chol_dec(
                  ext.data(),               // matrix A, will also store output
		  &idv_[0],                // diagnal vector
		  length_);		    // order of matrix A
  return true;
}


/// Solve A x = b (where A previously given to decompose)

template <typename T>
template <typename Block0,
	  typename Block1>
bool
Chold_impl<T,Mercury_sal_tag>::impl_solve(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  assert(b.size(0) == length_);
  assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

  Matrix<T, b_data_block_type> b_int(b.size(0), b.size(1));
  Matrix<T, b_data_block_type> x_int(b.size(0), b.size(1));
  b_int = b;

  if (length_ > 1) 
  {
    Ext_data<b_data_block_type> b_ext(b_int.block());
    Ext_data<b_data_block_type> x_ext(x_int.block());
    Ext_data<data_block_type> a_ext(data_.block());

    ptr_type b_ptr = b_ext.data();
    ptr_type x_ptr = x_ext.data();

    for(index_type i=0;i<b.size(1);i++) {
      sal_mat_chol_sol(
                 a_ext.data(), a_ext.stride(0),
		 &idv_[0],
		 storage_type::offset(b_ptr,i*length_),
		 storage_type::offset(x_ptr,i*length_),
		 length_);
    }
  }
  else 
  {
    for(index_type i=0;i<b.size(1);i++)
      x_int.put(0,i,b.get(0,i)/data_.get(0,0));
  }
  x = x_int;
  return true;
}

} // namespace vsip::impl
} // namespace vsip

#endif
