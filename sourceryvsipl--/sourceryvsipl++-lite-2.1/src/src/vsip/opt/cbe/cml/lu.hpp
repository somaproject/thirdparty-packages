/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/cml/lu.hpp
    @author  Brooks Moses
    @date    2008-07-07
    @brief   VSIPL++ Library: LU linear system solver using CML.
*/

#ifndef VSIP_OPT_CBE_CML_LU_HPP
#define VSIP_OPT_CBE_CML_LU_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/signal/fir_backend.hpp>
#include <vsip/opt/dispatch.hpp>

#include <cml.h>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{


// A structure that tells us if CML LU impl is available
// for certain types
template <>
struct Is_lud_impl_avail<Cml_tag, float>
{
  static bool const value = true;
};


/// LUD implementation using CML library.

/// Requires:
///   T to be a value type supported by CML's LU routines (i.e., float)

template <typename T>
class Lud_impl<T, Cml_tag>
{
  typedef vsip::impl::dense_complex_type   complex_type;

  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  typedef Layout<2, col2_type, Stride_unit_dense, complex_type> t_data_LP;
  typedef Fast_block<2, T, t_data_LP> t_data_block_type;

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
  length_type length() const VSIP_NOTHROW { return length_; }

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

  length_type max_decompose_size();

  // Member data.
private:
  typedef std::vector<T, Aligned_allocator<T> > vector_type;

  length_type  length_;                   // Order of A.

  Matrix<T, data_block_type> data_;	// Factorized matrix A

  cml_lud_f lud_obj_handle;
};


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T>
Lud_impl<T, Cml_tag>::Lud_impl(
  length_type  rows
  )
VSIP_THROW((std::bad_alloc))
  : length_     (rows),
    data_       (length_, length_)
{
  assert(length_ > 0);

  if (!cml_lud_create_f(&lud_obj_handle, rows))
    VSIP_IMPL_THROW(std::bad_alloc());
}



template <typename T>
Lud_impl<T, Cml_tag>::Lud_impl(Lud_impl const& lu)
VSIP_THROW((std::bad_alloc))
  : length_     (lu.length_),
    data_       (length_, length_)
{
  // Copy the LU product data.
  data_ = lu.data_;

  // Create a new LU object of the appropriate size.
  if (!cml_lud_create_f(&lud_obj_handle, length_))
    VSIP_IMPL_THROW(std::bad_alloc());

  // Set an appropriate stride value for the LU matrix.  (Note that
  // the actual LU matrix pointer is always reset before use.)
  lud_obj_handle.LU_stride_0 = length_;

  // Copy permutation matrix.
  memcpy(lu.lud_obj_handle.P, lud_obj_handle.P, length_ * sizeof(size_t));
}



template <typename T>
Lud_impl<T, Cml_tag>::~Lud_impl()
  VSIP_NOTHROW
{
  cml_lud_destroy_f(&lud_obj_handle);
}



/// Form LU factorization of matrix M
///
/// Requires
///   M to be a square matrix

template <typename T>
template <typename Block>
bool
Lud_impl<T, Cml_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == length_ && m.size(1) == length_);

  assign_local(data_, m);

  Ext_data<data_block_type> ext(data_.block());

  return cml_lud_decompose_f(&lud_obj_handle, ext.data(), length_); 
}



/// Solve Op(A) x = b (where A previously given to decompose)
///
/// Op(A) is
///   A   if tr == mat_ntrans
///   A'T if tr == mat_trans
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
Lud_impl<T,Cml_tag>::impl_solve(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  // CML only supports float LU decompositions, so calling this with
  // a hermitian is invalid.
  assert(tr != mat_herm);

  // Check that input and output sizes match.
  assert(x.size(0) == length_);
  assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

  // We need to put the LU storage array back into ext_data form,
  // and make sure the object pointer points to the correct place
  // still.
  Ext_data<data_block_type> ext(data_.block());
  lud_obj_handle.LU = ext.data();

  // CML does an in-place solve, so we must first do a copy.  We
  // use a local temporary to ensure it has the appropriate block
  // type.
  Matrix<T, data_block_type> x_local(x.size(0), x.size(1));
  assign_local(x_local, b);

  {
    Ext_data<data_block_type> x_ext(x_local.block());

    if(tr == mat_ntrans) 
    {
      cml_lud_sol_f(&lud_obj_handle,
        x_ext.data(),
        x.size(1),
        x.size(0),
        x.size(1));
    }
    else
    {
      cml_lud_solt_f(&lud_obj_handle,
        x_ext.data(),
        x.size(1),
        x.size(0),
        x.size(1));
    }
  }

  assign_local(x, x_local);
  
  return true;
}

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CBE_CML_LU_HPP
