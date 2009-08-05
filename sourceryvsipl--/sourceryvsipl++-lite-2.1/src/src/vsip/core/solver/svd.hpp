/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/solver/svd.hpp
    @author  Jules Bergmann
    @date    2005-09-11
    @brief   VSIPL++ Library: SVD Linear system solver.

*/

#ifndef VSIP_CORE_SOLVER_SVD_HPP
#define VSIP_CORE_SOLVER_SVD_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/math.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#ifdef VSIP_IMPL_HAVE_LAPACK
#  include <vsip/opt/lapack/bindings.hpp>
#  include <vsip/opt/lapack/svd.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/svd.hpp>
#endif

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// List of implementation tags to consider for LU.

typedef Make_type_list<
#ifdef VSIP_IMPL_HAVE_SAL
  Mercury_sal_tag,
#endif
#ifdef VSIP_IMPL_HAVE_LAPACK
  Lapack_tag,
#endif
  None_type // None_type is treated specially by Make_type_list, it is
            // not be put into the list.  Putting an explicit None_type
            // at the end of the list lets us put a ',' after each impl
            // tag.
  >::type Svd_type_list;



// a structure to chose implementation type
template <typename T>
struct Choose_svd_impl
{
#ifdef VSIP_IMPL_REF_IMPL
  typedef Cvsip_tag use_type;
  typedef Cvsip_tag type;
#else
  typedef typename Choose_solver_impl<
    Is_svd_impl_avail,
    T,
    Svd_type_list>::type type;

  typedef typename ITE_Type<
    Type_equal<type, None_type>::value,
    As_type<Error_no_solver_for_this_type>,
    As_type<type> >::type use_type;
#endif
};

} // namespace vsip::impl


/// SVD solver object.

template <typename              T               = VSIP_DEFAULT_VALUE_TYPE,
	  return_mechanism_type ReturnMechanism = by_value>
class svd;

template <typename T>
class svd<T, by_reference>
  : public impl::Svd_impl<T,true,typename impl::Choose_svd_impl<T>::use_type>
{
  typedef impl::Svd_impl<T,true, typename impl::Choose_svd_impl<T>::use_type>
    base_type;

  // Constructors, copies, assignments, and destructors.
public:
  svd(length_type rows, length_type cols, storage_type ust, storage_type vst)
    VSIP_THROW((std::bad_alloc))
    : base_type(rows, cols, ust, vst)
  {}

  ~svd() VSIP_NOTHROW {}

  // By-reference solvers.
public:
  template <typename Block0,
	    typename Block1>
  bool decompose(Matrix<T, Block0> m, Vector<scalar_f, Block1> dest)
    VSIP_NOTHROW
  { return this->impl_decompose(m, dest); }

  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block0,
	    typename          Block1>
  bool produ(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  // { return true; }
  { return this->template impl_produ<tr, ps>(b, x); }

  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block0,
	    typename          Block1>
  bool prodv(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->template impl_prodv<tr, ps>(b, x); }

  template <typename Block>
  bool u(index_type low, index_type high, Matrix<T, Block> dest)
    VSIP_NOTHROW
  { return this->template impl_u(low, high, dest); }

  template <typename Block>
  bool v(index_type low, index_type high, Matrix<T, Block> dest)
    VSIP_NOTHROW
  { return this->template impl_v(low, high, dest); }
};

template <typename T>
class svd<T, by_value>
  : public impl::Svd_impl<T,true,typename impl::Choose_svd_impl<T>::use_type>
{
  typedef impl::Svd_impl<T,true,typename impl::Choose_svd_impl<T>::use_type>
    base_type;

  // Constructors, copies, assignments, and destructors.
public:
  svd(length_type rows, length_type cols, storage_type ust, storage_type vst)
    VSIP_THROW((std::bad_alloc))
    : base_type(rows, cols, ust, vst)
  {}

  ~svd() VSIP_NOTHROW {}

  // By-value solvers.
public:
  template <typename Block0>
  Vector<scalar_f>
  decompose(Matrix<T, Block0> m)
    VSIP_THROW((std::bad_alloc, computation_error))
  {
    Vector<scalar_f> dest(this->impl_order());
    if (!this->impl_decompose(m, dest))
      VSIP_IMPL_THROW(computation_error("svd::decompose"));
    return dest;
  }

  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block>
  Matrix<T>
  produ(const_Matrix<T, Block> b)
    VSIP_THROW((std::bad_alloc, computation_error))
  {
    length_type q_rows = this->rows();
    length_type q_cols = this->ustorage() == svd_uvfull ? this->rows() :
                                                          this->impl_order();

    length_type x_rows, x_cols;
    if (ps == mat_lside)
    {
      x_rows = (tr == mat_ntrans) ? q_rows : q_cols;
      x_cols = b.size(1);
    }
    else /* (ps == mat_rside) */
    {
      x_rows = b.size(0);
      x_cols = (tr == mat_ntrans) ? q_cols : q_rows;
    }
    Matrix<T> x(x_rows, x_cols);
    this->template impl_produ<tr, ps>(b, x);
    return x;
  }

  template <mat_op_type       tr,
	    product_side_type ps,
	    typename          Block>
  Matrix<T>
  prodv(const_Matrix<T, Block> b)
    VSIP_THROW((std::bad_alloc, computation_error))
  { 
    length_type vt_rows = this->vstorage() == svd_uvfull ? this->columns() :
                                                           this->impl_order();
    length_type vt_cols = this->columns();

    length_type x_rows, x_cols;
    if (ps == mat_lside)
    {
      x_rows = (tr == mat_ntrans) ? vt_cols : vt_rows;
      x_cols = b.size(1);
    }
    else /* (ps == mat_rside) */
    {
      x_rows = b.size(0);
      x_cols = (tr == mat_ntrans) ? vt_rows : vt_cols;
    }
    Matrix<T> x(x_rows, x_cols);
    this->template impl_prodv<tr, ps>(b, x);
    return x;
  }

  Matrix<T>
  u(index_type low, index_type high)
    VSIP_THROW((std::bad_alloc, computation_error))
  {
    assert((this->ustorage() == svd_uvpart && high <= this->impl_order()) ||
	   (this->ustorage() == svd_uvfull && high <= this->rows()));

    Matrix<T> dest(this->rows(), high - low + 1);
    if (!this->template impl_u(low, high, dest))
      VSIP_IMPL_THROW(computation_error("svd::u"));
    return dest;
  }

  Matrix<T>
  v(index_type low, index_type high)
    VSIP_THROW((std::bad_alloc, computation_error))
  {
    assert((this->vstorage() == svd_uvpart && high <= this->impl_order()) ||
	   (this->vstorage() == svd_uvfull && high <= this->columns()));

    Matrix<T> dest(this->columns(), high - low + 1);
    if (!this->template impl_v(low, high, dest))
      VSIP_IMPL_THROW(computation_error("svd::v"));
    return dest;
  }
};

} // namespace vsip

#endif // VSIP_CORE_SOLVER_SVD_HPP
