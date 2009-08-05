/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/cml/qr.hpp
    @author  Brooks Moses
    @date    2008-06-12
    @brief   VSIPL++ Library: QR linear system solver using CML.
*/

#ifndef VSIP_OPT_CBE_CML_QR_HPP
#define VSIP_OPT_CBE_CML_QR_HPP

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


// A structure that tells us if CML qr impl is available
// for certain types
template <>
struct Is_qrd_impl_avail<Cml_tag, float>
{
  static bool const value = true;
};


/// Qrd implementation using CML library.

/// Requires:
///   T to be a value type supported by CML's QR routines (i.e., float)
///   BLOCKED is not used (it is used by the Lapack QR implementation
///      class).

template <typename T,
	  bool     Blocked>
class Qrd_impl<T, Blocked, Cml_tag>
{
  typedef vsip::impl::dense_complex_type   complex_type;

  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  typedef Layout<2, col2_type, Stride_unit_dense, complex_type> t_data_LP;
  typedef Fast_block<2, T, t_data_LP> t_data_block_type;

  // Qrd types supported.
protected:
  static bool const supports_qrd_saveq1  = false;
  static bool const supports_qrd_saveq   = true;
  static bool const supports_qrd_nosaveq = true;

  // Constructors, copies, assignments, and destructors.
public:
  Qrd_impl(length_type, length_type, storage_type)
    VSIP_THROW((std::bad_alloc));
  Qrd_impl(Qrd_impl const&)
    VSIP_THROW((std::bad_alloc));

  Qrd_impl& operator=(Qrd_impl const&) VSIP_NOTHROW;
  ~Qrd_impl() VSIP_NOTHROW;

  // Accessors.
public:
  length_type  rows()     const VSIP_NOTHROW { return m_; }
  length_type  columns()  const VSIP_NOTHROW { return n_; }
  storage_type qstorage() const VSIP_NOTHROW { return st_; }

  // Solve systems.
public:
  template <typename Block>
  bool decompose(Matrix<T, Block>) VSIP_NOTHROW;

protected:
  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block0,
	    typename          Block1>
  bool impl_prodq(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  template <mat_op_type       tr,
	    typename          Block0,
	    typename          Block1>
  bool impl_rsol(const_Matrix<T, Block0>, T const, Matrix<T, Block1>)
    VSIP_NOTHROW;

  template <typename          Block0,
	    typename          Block1>
  bool impl_covsol(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  template <typename          Block0,
	    typename          Block1>
  bool impl_lsqsol(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  // Member data.
private:
  typedef std::vector<T, Aligned_allocator<T> > vector_type;

  length_type  m_;			// Number of rows.
  length_type  n_;			// Number of cols.
  storage_type st_;			// Q storage type

  Matrix<T, data_block_type> data_;	// Factorized QR(mxn) matrix
  Matrix<T, t_data_block_type> t_data_;	// R (nxn) matrix, transposed

  cml_qrd_f qrd_obj_handle;
};


/***********************************************************************
  Definitions
***********************************************************************/

template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cml_tag>::Qrd_impl(
  length_type  rows,
  length_type  cols,
  storage_type st
  )
VSIP_THROW((std::bad_alloc))
  : m_          (rows),
    n_          (cols),
    st_         (st),
    data_       (m_, n_),
    t_data_     (n_, n_)
{
  assert(m_ > 0 && n_ > 0 && m_ >= n_);
  assert(st_ == qrd_nosaveq || st_ == qrd_saveq || st_ == qrd_saveq1);

  // CML stores Q as a set of reflectors, not as a constructed matrix,
  // and the result is equivalent to an m*m matrix for product purposes.
  if (st_ == qrd_saveq1)
    VSIP_IMPL_THROW(impl::unimplemented(
	      "Qrd does not support thin storage of Q (qrd_saveq1) when using CML"));

  if (!cml_qrd_create_f(&qrd_obj_handle, rows, cols))
    VSIP_IMPL_THROW(std::bad_alloc());
}



template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cml_tag>::Qrd_impl(Qrd_impl const& qr)
VSIP_THROW((std::bad_alloc))
  : m_          (qr.m_),
    n_          (qr.n_),
    st_         (qr.st_),
    data_       (m_, n_),
    t_data_     (n_, n_)
{
  // Copy the QR product data.
  data_ = qr.data_;

  // Create a new QR object of the appropriate size.
  if (!cml_qrd_create_f(&qrd_obj_handle, m_, n_))
    VSIP_IMPL_THROW(std::bad_alloc());

  // Set an appropriate stride value for the QR matrix.  (Note that
  // the actual QR matrix pointer is always reset before use.)
  qrd_obj_handle.QR_stride_0 = n_;

  // Copy tau matrix.
  cml_vcopy_f(qr.qrd_obj_handle.tau, 1, qrd_obj_handle.tau, 1, n_);
}



template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cml_tag>::~Qrd_impl()
  VSIP_NOTHROW
{
  cml_qrd_destroy_f(&qrd_obj_handle);
}



/// Decompose matrix M into QR form.
///
/// Requires
///   M to be a full rank, modifiable matrix of ROWS x COLS.

template <typename T,
	  bool     Blocked>
template <typename Block>
bool
Qrd_impl<T, Blocked, Cml_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == m_ && m.size(1) == n_);

  assign_local(data_, m);

  Ext_data<data_block_type> ext(data_.block());

  return cml_qrd_decompose_f(&qrd_obj_handle, ext.data(), n_); 
}



/// Compute product of Q and b

template <typename T,
	  bool     Blocked>
template <mat_op_type       tr,
	  product_side_type ps,
	  typename          Block0,
	  typename          Block1>
bool
Qrd_impl<T, Blocked, Cml_tag>::impl_prodq(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  // The impl_prodq function shouldn't be called for nosaveq, and we shouldn't
  // ever have a storage type of qrd_saveq1 for CML.
  assert(this->qstorage() == qrd_saveq);

  // Similarly, CML only supports float QR decompositions, so calling this with
  // a hermitian is invalid.
  assert(tr != mat_herm);

  // Check input and output sizes.  
  assert(b.size(0) == x.size(0));
  assert(b.size(1) == x.size(1));

  // We need to put the QR storage array back into ext_data form,
  // and make sure the object pointer points to the correct place
  // still.
  Ext_data<data_block_type> ext(data_.block());
  qrd_obj_handle.QR = ext.data();

  if (ps == mat_lside)
  {
    // CML does an in-place product, so we must first do a copy.  We
    // use a local temporary to ensure it has the appropriate block
    // type.
    assert(x.size(0) == m_);

    Matrix<T, data_block_type> x_local(x.size(0), x.size(1));
    assign_local(x_local, b);

    {
      Ext_data<data_block_type> x_ext(x_local.block());

      if(tr == mat_trans) 
      {
        cml_qrd_prodqt_f(&qrd_obj_handle,
          x_ext.data(),
          x.size(1),
          x.size(0),
          x.size(1));
      }
      else 
      {
        cml_qrd_prodq_f(&qrd_obj_handle,
          x_ext.data(),
          x.size(1),
          x.size(0),
          x.size(1));
      }
    }

    assign_local(x, x_local);
  }
  else  // ps == mat_rside
  {
    // We do a right-side product (X = BQ or X = BQ') by transposing
    // both sides, to get X' = Q'B' or X' = QB'.  The transposes are
    // handled implicitly by putting the X/B matrix into a column-major
    // block.  Note that the logic for whether to use prodq or prodqt
    // is inverted compared to the left-side case.

    assert(x.size(1) == m_);
    Matrix<T, t_data_block_type> x_local(x.size(0), x.size(1));
    assign_local(x_local, b);

    {
      Ext_data<t_data_block_type> x_ext(x_local.block());

      if(tr == mat_trans) 
      {
        cml_qrd_prodq_f(&qrd_obj_handle,
          x_ext.data(),
          x.size(0),
          x.size(1),
          x.size(0));
      }
      else 
      {
        cml_qrd_prodqt_f(&qrd_obj_handle,
          x_ext.data(),
          x.size(0),
          x.size(1),
          x.size(0));
      }
    }

    assign_local(x, x_local);
  }
  
  return true;
}



/// Solve op(R) x = alpha b

template <typename T,
	  bool     Blocked>
template <mat_op_type tr,
	  typename    Block0,
	  typename    Block1>
bool
Qrd_impl<T, Blocked, Cml_tag>::impl_rsol(
  const_Matrix<T, Block0> b,
  T const                 alpha,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  assert(b.size(0) == n_);
  assert(b.size(0) == x.size(0));
  assert(b.size(1) == x.size(1));

  // CML does an in-place solve, so we must first do a copy.  We
  // use a local temporary to ensure it has the appropriate block
  // type.
  Matrix<T, data_block_type> x_local(x.size(0), x.size(1));

  if(tr == mat_ntrans)
  {
    // We need to put the QR storage array back into ext_data form,
    // and make sure the object pointer points to the correct place
    // still.
    Ext_data<data_block_type> ext(data_.block());
    qrd_obj_handle.QR = ext.data();

    assign_local(x_local, b);

    // multiply b by alpha
    x_local *= alpha;

    {
      Ext_data<data_block_type> x_ext(x_local.block());

      cml_qrd_solr_f(&qrd_obj_handle,
        x_ext.data(),
        x.size(1),
        x.size(0),
        x.size(1));
    }

    assign_local(x, x_local);
  }
  else
  {
    // CML currently doesn't directly support transposed Rsolve. 
    // However, we know that R is stored in the upper triangle of the
    // QR matrix, so we can create a transpose of this and then use
    // the trisolve routine to compute the solution.

    t_data_ = data_(Domain<2>(Domain<1>(0, 1, n_),
			      Domain<1>(0, 1, n_)));
    Ext_data<t_data_block_type> ext(t_data_.block());

    assign_local(x_local, b);

    // multiply b by alpha
    x_local *= alpha;

    {
      Ext_data<data_block_type> x_ext(x_local.block());

      cml_trisolve_lower_f(ext.data(), n_,
        x_ext.data(), x.size(1),
        x_ext.data(), x.size(1),
        x.size(0),
	x.size(1),
	true, false);
    }

    assign_local(x,x_local);
  }
 
  return true;
}


/// Solve covariance system for x:
///   A' A X = B

template <typename T,
	  bool     Blocked>
template <typename Block0,
	  typename Block1>
bool
Qrd_impl<T, Blocked, Cml_tag>::
impl_covsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
  VSIP_NOTHROW
{
  assert(b.size(0) == n_);
  assert(b.size(0) == x.size(0));
  assert(b.size(1) == x.size(1));

  // We need to put the QR storage array back into ext_data form,
  // and make sure the object pointer points to the correct place
  // still.
  Ext_data<data_block_type> ext(data_.block());
  qrd_obj_handle.QR = ext.data();

  // CML currently doesn't directly support transposed Rsolve. 
  // However, we know that R is stored in the upper triangle of the
  // QR matrix, so we can create a transpose of this and then use
  // the trisolve routine to compute the solution.

  t_data_ = data_(Domain<2>(Domain<1>(0, 1, n_),
			      Domain<1>(0, 1, n_)));
  Ext_data<t_data_block_type> t_ext(t_data_.block());

  // CML does an in-place solve, so we must first do a copy.  We
  // use a local temporary to ensure it has the appropriate block
  // type.
  Matrix<T, data_block_type> x_local(x.size(0), x.size(1));
  assign_local(x_local, b);

  {
    Ext_data<data_block_type> x_ext(x_local.block());

    // First, we do a transposed rsolve, to find AX.
    
    cml_trisolve_lower_f(t_ext.data(), n_,
      x_ext.data(), x.size(1),
      x_ext.data(), x.size(1),
      x.size(0),
      x.size(1),
      true, false);

    // Then, a non-transposed rsolve, to find X.

    cml_qrd_solr_f(&qrd_obj_handle,
      x_ext.data(),
      x.size(1),
      x.size(0),
      x.size(1));
  }

  assign_local(x, x_local);
  return true;
}



/// Solve linear least squares problem for x:
///   min_x norm-2( A x - b )

template <typename T,
	  bool     Blocked>
template <typename Block0,
	  typename Block1>
bool
Qrd_impl<T, Blocked, Cml_tag>::
impl_lsqsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
  VSIP_NOTHROW
{

  assert(b.size(0) == m_);
  assert(x.size(0) == n_);
  assert(x.size(1) == b.size(1));
 
  // Solve  A X = B  for X
  //
  // 0. factor:             QR X = B
  //    mult by Q'        Q'QR X = Q'B
  //    simplify             R X = Q'B
  //
  // 1. compute C = Q'B:     R X = C
  // 2. solve for X:         R X = C

  // We need to put the QR storage array back into ext_data form,
  // and make sure the object pointer points to the correct place
  // still.
  Ext_data<data_block_type> ext(data_.block());
  qrd_obj_handle.QR = ext.data();

  // CML does an in-place solve, so we must first do a copy.  We
  // use a local temporary to ensure it has the appropriate block
  // type.
  Matrix<T, data_block_type> x_local(b.size(0), b.size(1));
  assign_local(x_local, b);

  {
    Ext_data<data_block_type> x_ext(x_local.block());

    cml_qrd_prodqt_f(&qrd_obj_handle,
      x_ext.data(),
      b.size(1),
      b.size(0),
      b.size(1));

   // Here, we do some sleight of hand; we really want the "thin Q"
   // product, which would have been only the first x.size(0) columns
   // of Q, and thus would have produced a result with only that many
   // rows.  However, what CML gives us is a "full Q" product with
   // b.size(0) rows.  But, since x_ext.data() is stored in row-major
   // format, all we have to do to get what we want is just ignore
   // the trailing rows of x_ext.data().

    cml_qrd_solr_f(&qrd_obj_handle,
      x_ext.data(),
      x.size(1),
      x.size(0),
      x.size(1));
  }

  // Of course, we then must assign only the relevant part of x_local
  // to x.
  assign_local(x, x_local(Domain<2>(Domain<1>(0, 1, n_),
			            Domain<1>(0, 1, x.size(1)))));

  return true;
}

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CBE_CML_QR_HPP
