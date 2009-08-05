/* Copyright (c) 2005, 2006, 2008 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/solver/qr.hpp
    @author  Jules Bergmann
    @date    2005-08-19
    @brief   VSIPL++ Library: QR Linear system solver.

*/

#ifndef VSIP_CORE_SOLVER_QR_HPP
#define VSIP_CORE_SOLVER_QR_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/working_view.hpp>
#include <vsip/core/solver/common.hpp>
#ifdef VSIP_IMPL_HAVE_LAPACK
#  include <vsip/opt/lapack/qr.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/qr.hpp>
#endif
#ifdef VSIP_IMPL_CBE_SDK
#  include <vsip/opt/cbe/cml/qr.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_CVSIP
#  include <vsip/core/cvsip/qr.hpp>
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// List of implementation tags to consider for QR.

typedef Make_type_list<
#ifdef VSIP_IMPL_HAVE_CVSIP
  Cvsip_tag,
#endif
#ifdef VSIP_IMPL_HAVE_SAL
  Mercury_sal_tag,
#endif
#ifdef VSIP_IMPL_CBE_SDK
  Cml_tag,
#endif
#ifdef VSIP_IMPL_HAVE_LAPACK
  Lapack_tag,
#endif
  None_type // None_type is treated specially by Make_type_list; it is
            // not put into the list.  Putting an explicit None_type
            // at the end of the list lets us put a ',' after each impl
            // tag.
  >::type Qrd_type_list;



// a structure to chose implementation type
template <typename T>
struct Choose_qrd_impl
{
#ifdef VSIP_IMPL_REF_IMPL
  typedef Cvsip_tag type;
  typedef Cvsip_tag use_type;
#else
  typedef typename Choose_solver_impl<
    Is_qrd_impl_avail,
    T,
    Qrd_type_list>::type type;

  typedef typename ITE_Type<
    Type_equal<type, None_type>::value,
    As_type<Error_no_solver_for_this_type>,
    As_type<type> >::type use_type;
#endif
};



// Qrd traits.  Determine which QR types (no-save, skinny, full) are
// supported.

template <typename QrT>
struct Qrd_traits
{
  static bool const supports_qrd_saveq1  = QrT::supports_qrd_saveq1;
  static bool const supports_qrd_saveq   = QrT::supports_qrd_saveq;
  static bool const supports_qrd_nosaveq = QrT::supports_qrd_nosaveq;
};

} // namespace vsip::impl

/// QR solver object.

template <typename              T               = VSIP_DEFAULT_VALUE_TYPE,
	  return_mechanism_type ReturnMechanism = by_value>
class qrd;



// QR solver object (by-reference).

template <typename T>
class qrd<T, by_reference>
  : public impl::Qrd_impl<T,true,typename impl::Choose_qrd_impl<T>::use_type>
{
  typedef typename impl::Choose_qrd_impl<T>::use_type use_type;
  typedef impl::Qrd_impl<T,true, use_type> base_type;

  // template <typename>
  friend struct impl::Qrd_traits<qrd<T, by_reference> >;

  // Constructors, copies, assignments, and destructors.
public:
  qrd(length_type rows, length_type cols, storage_type st)
    VSIP_THROW((std::bad_alloc))
      : base_type(rows, cols, st)
    { VSIP_IMPL_COVER_TAG("qrd", use_type); }

  ~qrd() VSIP_NOTHROW {}

  // By-reference solvers.
public:
  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block0,
	    typename          Block1>
  bool prodq(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->template impl_prodq<tr, ps>(b, x); }

  template <mat_op_type       tr,
	    typename          Block0,
	    typename          Block1>
  bool rsol(const_Matrix<T, Block0> b, T const alpha, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->template impl_rsol<tr>(b, alpha, x); }

  template <typename          Block0,
	    typename          Block1>
  bool covsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->impl_covsol(b, x); }

  template <typename          Block0,
	    typename          Block1>
  bool lsqsol(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->impl_lsqsol(b, x); }
};



// QR solver object (by-value).

template <typename T>
class qrd<T, by_value>
  : public impl::Qrd_impl<T,true,typename impl::Choose_qrd_impl<T>::use_type >
{
  typedef impl::Qrd_impl<T,true,typename impl::Choose_qrd_impl<T>::use_type >
    base_type;

  // Constructors, copies, assignments, and destructors.
public:
  qrd(length_type rows, length_type cols, storage_type st)
    VSIP_THROW((std::bad_alloc))
      : base_type(rows, cols, st)
    {}

  ~qrd() VSIP_NOTHROW {}

  // By-value solvers.
public:
  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block0>
  Matrix<T>
  prodq(const_Matrix<T, Block0> b)
    VSIP_NOTHROW
  {
    length_type x_rows, x_cols;
    if (ps == mat_lside)
    {
      x_rows = (this->qstorage() == qrd_saveq)     ? this->rows()    : 
	       (tr == mat_trans || tr == mat_herm) ? this->columns() :
	                                             this->rows();
      x_cols = b.size(1);
    }
    else
    {
      x_rows = b.size(0);
      x_cols = (this->qstorage() == qrd_saveq)     ? this->rows() : 
	       (tr == mat_trans || tr == mat_herm) ? this->rows() :
	                                             this->columns();
    }

    Matrix<T> x(x_rows, x_cols);
    this->template impl_prodq<tr, ps>(b, x);
    return x;
  }

  template <mat_op_type       tr,
	    typename          Block0>
  Matrix<T>
  rsol(const_Matrix<T, Block0> b, T const alpha)
    VSIP_NOTHROW
  {
    Matrix<T> x(b.size(0), b.size(1));
    this->template impl_rsol<tr>(b, alpha, x); 
    return x;
  }

  template <typename          Block0>
  Matrix<T>
  covsol(const_Matrix<T, Block0> b)
    VSIP_NOTHROW
  {
    Matrix<T> x(b.size(0), b.size(1));
    this->impl_covsol(b, x);
    return x;
  }

  template <typename          Block0>
  Matrix<T>
  lsqsol(const_Matrix<T, Block0> b)
    VSIP_NOTHROW
  {
    Matrix<T> x(this->columns(), b.size(1));
    this->impl_lsqsol(b, x);
    return x;
  }
};



} // namespace vsip


#endif // VSIP_CORE_SOLVER_QR_HPP
