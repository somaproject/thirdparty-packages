/* Copyright (c) 2006, 2008 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/cvsip/qr.hpp
    @author  Assem Salama
    @date    2006-10-26
    @brief   VSIPL++ Library: QR linear system solver using CVSIP.

*/

#ifndef VSIP_CORE_CVSIP_SOLVER_QR_HPP
#define VSIP_CORE_CVSIP_SOLVER_QR_HPP

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
#include <vsip/core/cvsip/solver.hpp>
#include <vsip/core/cvsip/view.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cvsip
{

template <typename T>
class Qrd : Non_copyable
{
  typedef Solver_traits<T> traits;

public:
  Qrd(length_type m, length_type n, storage_type op) : qr_(traits::qr_create(m, n, op)) {}
  ~Qrd() { traits::qr_destroy(qr_);}
  int decompose(View<2,T,true> &a) 
  { return !traits::qr_decompose(qr_, a.ptr());}
  int solve(View<2,T,true> &xb, vsip_qrd_prob p) 
  { return traits::qr_solve(qr_, p, xb.ptr());}
  int solve_r(mat_op_type op, T alpha, View<2,T,true> &m)
  { return traits::qr_solve_r(qr_, op, alpha, m.ptr());}
  int prodq(mat_op_type op, product_side_type side, View<2,T,true> &m) 
  { return traits::qr_prodq(qr_, op, side, m.ptr());}

private:
  typename traits::qr_type *qr_;
};

}

/// Qrd implementation using CVSIP

/// Requires:
///   T to be a value type supported by SAL's QR routines

template <typename T,
	  bool     Blocked>
class Qrd_impl<T, Blocked, Cvsip_tag>
{
  typedef vsip::impl::dense_complex_type   complex_type;
  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  // Qrd types supported.
protected:
  static bool const supports_qrd_saveq1  = true;
  static bool const supports_qrd_saveq   = false;
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
  cvsip::View<2,T,true>      cvsip_data_;
  cvsip::Qrd<T>              cvsip_qr_;
};



// The CVSIP QR solver supports all CVSIP types.

template <typename T>
struct Is_qrd_impl_avail<Cvsip_tag, T>
{
  static bool const value = cvsip::Solver_traits<T>::valid;
};



/***********************************************************************
  Definitions
***********************************************************************/

template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cvsip_tag>::Qrd_impl(
  length_type  rows,
  length_type  cols,
  storage_type st
  )
VSIP_THROW((std::bad_alloc))
  : m_          (rows),
    n_          (cols),
    st_         (st),
    data_       (m_, n_),
    cvsip_data_ (data_.block().impl_data(), m_, n_,true),
    cvsip_qr_   (m_, n_, st_)
{
  assert(m_ > 0 && n_ > 0 && m_ >= n_);
  assert(st_ == qrd_nosaveq || st_ == qrd_saveq || st_ == qrd_saveq1);
}



template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cvsip_tag>::Qrd_impl(Qrd_impl const& qr)
VSIP_THROW((std::bad_alloc))
  : m_          (qr.m_),
    n_          (qr.n_),
    st_         (qr.st_),
    data_       (m_, n_),
    cvsip_data_ (data_.block().impl_data(), m_, n_,true),
    cvsip_qr_   (m_, n_, st_)
{
  data_ = qr.data_;
}



template <typename T,
	  bool     Blocked>
Qrd_impl<T, Blocked, Cvsip_tag>::~Qrd_impl()
  VSIP_NOTHROW
{
}



/// Decompose matrix M into QR form.
///
/// Requires
///   M to be a full rank, modifiable matrix of ROWS x COLS.

template <typename T,
	  bool     Blocked>
template <typename Block>
bool
Qrd_impl<T, Blocked, Cvsip_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == m_ && m.size(1) == n_);

  cvsip_data_.block().release(false);
  assign_local(data_, m);
  cvsip_data_.block().admit(true);


  cvsip_qr_.decompose(cvsip_data_);

  return true;
}



/// Compute product of Q and b
/// 
/// If qstorage == qrd_saveq1, Q is MxN.
/// If qstorage == qrd_saveq,  Q is MxM.
///
/// qstoarge   | ps        | tr         | product | b (in) | x (out)
/// qrd_saveq1 | mat_lside | mat_ntrans | Q b     | (n, p) | (m, p)
/// qrd_saveq1 | mat_lside | mat_trans  | Q' b    | (m, p) | (n, p)
/// qrd_saveq1 | mat_lside | mat_herm   | Q* b    | (m, p) | (n, p)
///
/// qrd_saveq1 | mat_rside | mat_ntrans | b Q     | (p, m) | (p, n)
/// qrd_saveq1 | mat_rside | mat_trans  | b Q'    | (p, n) | (p, m)
/// qrd_saveq1 | mat_rside | mat_herm   | b Q*    | (p, n) | (p, m)
///
/// qrd_saveq  | mat_lside | mat_ntrans | Q b     | (m, p) | (m, p)
/// qrd_saveq  | mat_lside | mat_trans  | Q' b    | (m, p) | (m, p)
/// qrd_saveq  | mat_lside | mat_herm   | Q* b    | (m, p) | (m, p)
///
/// qrd_saveq  | mat_rside | mat_ntrans | b Q     | (p, m) | (p, m)
/// qrd_saveq  | mat_rside | mat_trans  | b Q'    | (p, m) | (p, m)
/// qrd_saveq  | mat_rside | mat_herm   | b Q*    | (p, m) | (p, m)

template <typename T,
	  bool     Blocked>
template <mat_op_type       tr,
	  product_side_type ps,
	  typename          Block0,
	  typename          Block1>
bool
Qrd_impl<T, Blocked, Cvsip_tag>::impl_prodq(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  typedef typename Block_layout<Block0>::order_type order_type;
  typedef typename Block_layout<Block0>::complex_type complex_type;
  typedef Layout<2, order_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP, Local_map> block_type;

  assert(this->qstorage() == qrd_saveq1 || this->qstorage() == qrd_saveq);
  length_type q_rows;
  length_type q_cols;
  if (qstorage() == qrd_saveq1)
  {
    q_rows = m_;
    q_cols = n_;
  }
  else // (qstorage() == qrd_saveq1)
  {
    q_rows = m_;
    q_cols = m_;
  }

  // do we need a transpose?
  if(tr == mat_trans || tr == mat_herm) 
  {
    std::swap(q_rows, q_cols);
  }
  if(tr == mat_herm) 
  {
    std::swap(q_rows, q_cols);
  }

  // left or right?
  if(ps == mat_lside) 
  {
    assert(b.size(0) == q_cols);
    assert(x.size(0) == q_rows);
    assert(b.size(1) == x.size(1));
  }
  else
  {
    assert(b.size(1) == q_rows);
    assert(x.size(1) == q_cols);
    assert(b.size(0) == x.size(0));
  }

  Matrix<T,block_type> b_int(b.size(0), b.size(1));
  Ext_data<block_type> b_ext(b_int.block());
  cvsip::View<2,T,true>
      cvsip_b_int(b_ext.data(),0,b_ext.stride(0),b_ext.size(0),
		                 b_ext.stride(1),b_ext.size(1));


  cvsip_b_int.block().release(false);
  assign_local(b_int, b);
  cvsip_b_int.block().admit(true);
  int ret = cvsip_qr_.prodq(tr,ps,cvsip_b_int);

  // now, copy into x
  cvsip_b_int.block().release(true);
  assign_local(x, b_int(Domain<2>(x.size(0),x.size(1))));

  

  return ret;
}



/// Solve op(R) x = alpha b

template <typename T,
	  bool     Blocked>
template <mat_op_type tr,
	  typename    Block0,
	  typename    Block1>
bool
Qrd_impl<T, Blocked, Cvsip_tag>::impl_rsol(
  const_Matrix<T, Block0> b,
  T const                 alpha,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  typedef typename Block_layout<Block0>::order_type order_type;
  typedef typename Block_layout<Block0>::complex_type complex_type;
  typedef Layout<2, order_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP, Local_map> block_type;

  assert(b.size(0) == n_);
  assert(b.size(0) == x.size(0));
  assert(b.size(1) == x.size(1));

  Matrix<T, block_type> b_int(b.size(0), b.size(1));
  Ext_data<block_type>  b_ext(b_int.block());
  cvsip::View<2,T,true>
      cvsip_b_int(b_ext.data(),0,b_ext.stride(0),b_ext.size(0),
		                 b_ext.stride(1),b_ext.size(1));

  cvsip_b_int.block().release(false);
  assign_local(b_int, b);
  cvsip_b_int.block().admit(true);

  cvsip_qr_.solve_r(tr, alpha, cvsip_b_int);

  // copy b_int back into x
  cvsip_b_int.block().release(true);
  assign_local(x, b_int);


  return true;
}



/// Solve covariance system for x:
///   A' A X = B

template <typename T,
	  bool     Blocked>
template <typename Block0,
	  typename Block1>
bool
Qrd_impl<T, Blocked, Cvsip_tag>::
impl_covsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
  VSIP_NOTHROW
{
  typedef typename Block_layout<Block0>::order_type order_type;
  typedef typename Block_layout<Block0>::complex_type complex_type;
  typedef Layout<2, order_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP, Local_map> block_type;

  Matrix<T, block_type> b_int(b.size(0), b.size(1));
  Ext_data<block_type>  b_ext(b_int.block());
  cvsip::View<2,T,true>
      cvsip_b_int(b_ext.data(),0,b_ext.stride(0),b_ext.size(0),
		                 b_ext.stride(1),b_ext.size(1));


  cvsip_b_int.block().release(false);
  assign_local(b_int, b);
  cvsip_b_int.block().admit(true);

  cvsip_qr_.solve(cvsip_b_int, VSIP_COV);

  // copy b_int back into x
  cvsip_b_int.block().release(true);
  assign_local(x, b_int);


  return true;
}



/// Solve linear least squares problem for x:
///   min_x norm-2( A x - b )

template <typename T,
	  bool     Blocked>
template <typename Block0,
	  typename Block1>
bool
Qrd_impl<T, Blocked, Cvsip_tag>::
impl_lsqsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
  VSIP_NOTHROW
{
  typedef typename Block_layout<Block0>::order_type order_type;
  typedef typename Block_layout<Block0>::complex_type complex_type;
  typedef Layout<2, order_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP, Local_map> block_type;

  Matrix<T, block_type> b_int(b.size(0), b.size(1));
  Ext_data<block_type> b_ext(b_int.block());
  cvsip::View<2,T,true>
      cvsip_b_int(b_ext.data(),0,b_ext.stride(0),b_ext.size(0),
		                 b_ext.stride(1),b_ext.size(1));

  cvsip_b_int.block().release(false);
  assign_local(b_int, b);
  cvsip_b_int.block().admit(true);

  cvsip_qr_.solve(cvsip_b_int, VSIP_LLS);
  // copy b_int back into x
  cvsip_b_int.block().release(true);
  assign_local(x, b_int(Domain<2>(x.size(0),x.size(1))));

  return true;
}

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_CORE_CVSIP_SOLVER_QR_HPP
