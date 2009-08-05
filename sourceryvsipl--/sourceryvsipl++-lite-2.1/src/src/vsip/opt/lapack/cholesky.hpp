/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/lapack/cholesky.hpp
    @author  Assem Salama
    @date    2006-04-13
    @brief   VSIPL++ Library: Cholesky Linear system solver using LAPACK.

*/

#ifndef VSIP_OPT_LAPACK_CHOLESKY_HPP
#define VSIP_OPT_LAPACK_CHOLESKY_HPP

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
#include <vsip/opt/lapack/bindings.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/solver/common.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// The Lapack Cholesky solver supports all BLAS types.
template <typename T>
struct Is_chold_impl_avail<Lapack_tag, T>
{
  static bool const value = blas::Blas_traits<T>::valid;
};



/// Cholesky factorization implementation class.  Common functionality
/// for chold by-value and by-reference classes.

template <typename T>
class Chold_impl<T, Lapack_tag>
  : Compile_time_assert<blas::Blas_traits<T>::valid>
{
  // BLAS/LAPACK require complex data to be in interleaved format.
  typedef Layout<2, col2_type, Stride_unit_dense, Cmplx_inter_fmt> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

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
  length_type length()const VSIP_NOTHROW { return length_; }
  mat_uplo    uplo()  const VSIP_NOTHROW { return uplo_; }

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
  typedef std::vector<T, Aligned_allocator<T> > vector_type;

  mat_uplo     uplo_;			// A upper/lower triangular
  length_type  length_;			// Order of A.

  Matrix<T, data_block_type> data_;	// Factorized Cholesky matrix (A)
};



template <typename T>
Chold_impl<T,Lapack_tag>::Chold_impl(
  mat_uplo    uplo,
  length_type length
  )
VSIP_THROW((std::bad_alloc))
  : uplo_   (uplo),
    length_ (length),
    data_   (length_, length_)
{
  assert(length_ > 0);
  assert(uplo_ == upper || uplo_ == lower);
}



template <typename T>
Chold_impl<T,Lapack_tag>::Chold_impl(Chold_impl const& qr)
VSIP_THROW((std::bad_alloc))
  : uplo_       (qr.uplo_),
    length_     (qr.length_),
    data_       (length_, length_)
{
  data_ = qr.data_;
}



template <typename T>
Chold_impl<T,Lapack_tag>::~Chold_impl()
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
Chold_impl<T,Lapack_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == length_ && m.size(1) == length_);

  data_ = m;

  Ext_data<data_block_type> ext(data_.block());

  bool success = lapack::potrf(
		uplo_ == upper ? 'U' : 'L', // A upper/lower lower triangular
		length_,		    // order of matrix A
		ext.data(),		    // matrix A
		ext.stride(1));		    // lda - first dim of A

  return success;
}



/// Solve A x = b (where A previously given to decompose)

template <typename T>
template <typename Block0,
	  typename Block1>
bool
Chold_impl<T,Lapack_tag>::impl_solve(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  assert(b.size(0) == length_);
  assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

  Matrix<T, data_block_type> b_int(b.size(0), b.size(1));
  b_int = b;

  {
    Ext_data<data_block_type> b_ext(b_int.block());
    Ext_data<data_block_type> a_ext(data_.block());

    lapack::potrs(uplo_ == upper ? 'U' : 'L',
		  length_,
		  b.size(1),		    // number of RHS systems
		  a_ext.data(), a_ext.stride(1), // A, lda
		  b_ext.data(), b_ext.stride(1));  // B, ldb
  }
  x = b_int;

  return true;
}

} // namespace vsip::impl
} // namespace vsip


#endif // VSIP_OPT_LAPACK_CHOLESKY_HPP
