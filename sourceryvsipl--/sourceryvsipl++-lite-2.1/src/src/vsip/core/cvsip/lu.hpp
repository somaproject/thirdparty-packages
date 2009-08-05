/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/cvsip/lu.hpp
    @author  Assem Salama
    @date    2006-04-13
    @brief   VSIPL++ Library: LU linear system solver using cvsip.

*/

#ifndef VSIP_CORE_CVSIP_LU_HPP
#define VSIP_CORE_CVSIP_LU_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/solver/common.hpp>
#include <vsip/core/cvsip/solver.hpp>
#include <vsip/core/cvsip/view.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// LU factorization implementation class.  Common functionality
/// for lud by-value and by-reference classes.

template <typename T>
class Lud_impl<T, Cvsip_tag>
  : Compile_time_assert<cvsip::Solver_traits<T>::valid>
{
  typedef cvsip::Solver_traits<T> traits;
  typedef Layout<2, row2_type, Stride_unit_dense, Cmplx_inter_fmt> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

public:
  Lud_impl(length_type length)
    : length_(length),
      data_(length_, length_),
      cvsip_data_(data_.block().impl_data(), length_, length_, true),
      lu_(traits::lu_create(length_))
  { assert(length_ > 0);}
  Lud_impl(Lud_impl const& lu)
    : length_(lu.length_),
      data_(length_, length_),
      cvsip_data_(data_.block().impl_data(), length_, length_, true),
      lu_(traits::lu_create(length_))
  { data_ = lu.data_;}
  ~Lud_impl() VSIP_NOTHROW { traits::lu_destroy(lu_);}

  length_type length()const VSIP_NOTHROW { return length_;}

  /// Form LU factorization of matrix A
  ///
  /// Requires
  ///   A to be a square matrix, either
  ///
  /// FLOPS:
  ///   real   : UPDATE
  ///   complex: UPDATE
  //
  template <typename Block>
  bool decompose(Matrix<T, Block> m) VSIP_NOTHROW
  {
    assert(m.size(0) == length_ && m.size(1) == length_);

    cvsip_data_.block().release(false);
    assign_local(data_, m);
    cvsip_data_.block().admit(true);
    bool success = !traits::lu_decompose(lu_, cvsip_data_.ptr());
    return success;
  }

protected:
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
  template <mat_op_type tr, typename Block0, typename Block1>
  bool impl_solve(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
  {
    typedef typename Block_layout<Block0>::order_type order_type;
    typedef typename Block_layout<Block0>::complex_type complex_type;
    typedef Layout<2, order_type, Stride_unit_dense, complex_type> data_LP;
    typedef Fast_block<2, T, data_LP, Local_map> block_type;

    assert(b.size(0) == length_);
    assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

    Matrix<T, block_type> b_int(b.size(0), b.size(1));
    assign_local(b_int, b);

    if (tr == mat_conj || 
        (tr == mat_trans && Is_complex<T>::value) ||
        (tr == mat_herm && !Is_complex<T>::value))
      VSIP_IMPL_THROW(unimplemented(
        "LU solver (CVSIP backend) does not implement this transformation"));
    {
      Ext_data<block_type> b_ext(b_int.block());

      cvsip::View<2,T,true>
        cvsip_b_int(b_ext.data(),0,b_ext.stride(0),b_ext.size(0),
                    b_ext.stride(1),b_ext.size(1));

      cvsip_b_int.block().admit(true);
      traits::lu_solve(lu_, tr, cvsip_b_int.ptr());
      cvsip_b_int.block().release(true);
    }
    assign_local(x, b_int);
    return true;
  }

private:
  Lud_impl& operator=(Lud_impl const&) VSIP_NOTHROW;

  length_type  length_;			// Order of A.
  Matrix<T, data_block_type> data_;	// Factorized Cholesky matrix (A)
  cvsip::View<2,T,true>      cvsip_data_;
  typename traits::lud_type *lu_;
};

// The CVSIP LU solver supports all CVSIP types.
template <typename T>
struct Is_lud_impl_avail<Cvsip_tag, T>
{
  static bool const value = cvsip::Solver_traits<T>::valid;
};

} // namespace vsip::impl
} // namespace vsip


#endif
