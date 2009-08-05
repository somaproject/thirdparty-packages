/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/lapack/lu.hpp
    @author  Assem Salama
    @date    2006-04-13
    @brief   VSIPL++ Library: LU linear system solver using lapack.

*/

#ifndef VSIP_CORE_LAPACK_LU_HPP
#define VSIP_CORE_LAPACK_LU_HPP

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

// The Lapack LU solver supports all BLAS types.
template <typename T>
struct Is_lud_impl_avail<Lapack_tag, T>
{
  static bool const value = blas::Blas_traits<T>::valid;
};



/// LU factorization implementation class.  Common functionality
/// for lud by-value and by-reference classes.

template <typename T>
class Lud_impl<T, Lapack_tag>
  : Compile_time_assert<blas::Blas_traits<T>::valid>
{
  // BLAS/LAPACK require complex data to be in interleaved format.
  typedef Layout<2, col2_type, Stride_unit_dense, Cmplx_inter_fmt> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  // Constructors, copies, assignments, and destructors.
public:
  Lud_impl(length_type)
    VSIP_THROW((std::bad_alloc));
  Lud_impl(Lud_impl const&)
    VSIP_THROW((std::bad_alloc));

  Lud_impl& operator=(Lud_impl const&) VSIP_NOTHROW;
  ~Lud_impl() VSIP_NOTHROW;

  // Accessors.
public:
  length_type length()const VSIP_NOTHROW { return length_; }

  // Solve systems.
public:
  template <typename Block>
  bool decompose(Matrix<T, Block>) VSIP_NOTHROW;

protected:
  template <mat_op_type tr,
	    typename    Block0,
	    typename    Block1>
  bool impl_solve(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  // Member data.
private:
  typedef std::vector<int, Aligned_allocator<int> > vector_type;

  length_type  length_;			// Order of A.
  vector_type  ipiv_;			// Additional info on Q

  Matrix<T, data_block_type> data_;	// Factorized Cholesky matrix (A)
};

} // namespace vsip::impl


/***********************************************************************
  Definitions
***********************************************************************/

namespace impl
{

template <typename T>
Lud_impl<T, Lapack_tag>::Lud_impl(
  length_type length
  )
VSIP_THROW((std::bad_alloc))
  : length_ (length),
    ipiv_   (length_),
    data_   (length_, length_)
{
  assert(length_ > 0);
}



template <typename T>
Lud_impl<T, Lapack_tag>::Lud_impl(Lud_impl const& lu)
VSIP_THROW((std::bad_alloc))
  : length_ (lu.length_),
    ipiv_   (length_),
    data_   (length_, length_)
{
  data_ = lu.data_;
  for (index_type i=0; i<length_; ++i)
    ipiv_[i] = lu.ipiv_[i];
}



template <typename T>
Lud_impl<T, Lapack_tag>::~Lud_impl()
  VSIP_NOTHROW
{
}



/// Form LU factorization of matrix A
///
/// Requires
///   A to be a square matrix, either
///
/// FLOPS:
///   real   : UPDATE
///   complex: UPDATE

template <typename T>
template <typename Block>
bool
Lud_impl<T, Lapack_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == length_ && m.size(1) == length_);

  assign_local(data_, m);

  Ext_data<data_block_type> ext(data_.block());

  bool success = lapack::getrf(
		length_, length_,
		ext.data(), ext.stride(1),	// matrix A, ldA
		&ipiv_[0]);			// pivots

  return success;
}



/// Solve Op(A) x = b (where A previously given to decompose)
///
/// Op(A) is
///   A   if tr == mat_ntrans
///   A^T if tr == mat_trans
///   A'  if tr == mat_herm (valid for T complex only)
///
/// Requires
///   B to be a (length, P) matrix
///   X to be a (length, P) matrix
///
/// Effects:
///   X contains solution to Op(A) X = B

template <typename T>
template <mat_op_type tr,
	  typename    Block0,
	  typename    Block1>
bool
Lud_impl<T, Lapack_tag>::impl_solve(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  assert(b.size(0) == length_);
  assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

  char trans;

  Matrix<T, data_block_type> b_int(b.size(0), b.size(1));
  assign_local(b_int, b);

  if (tr == mat_ntrans)
    trans = 'N';
  else if (tr == mat_trans)
    trans = 'T';
  else if (tr == mat_herm)
  {
    assert(Is_complex<T>::value);
    trans = 'C';
  }

  {
    Ext_data<data_block_type> b_ext(b_int.block());
    Ext_data<data_block_type> a_ext(data_.block());

    lapack::getrs(trans,
		  length_,			  // order of A
		  b.size(1),			  // nrhs: number of RH sides
		  a_ext.data(), a_ext.stride(1),  // A, lda
		  &ipiv_[0],			  // pivots
		  b_ext.data(), b_ext.stride(1)); // B, ldb
  }
  assign_local(x, b_int);

  return true;
}

} // namespace vsip::impl

} // namespace vsip


#endif // VSIP_CORE_LAPACK_LU_HPP
