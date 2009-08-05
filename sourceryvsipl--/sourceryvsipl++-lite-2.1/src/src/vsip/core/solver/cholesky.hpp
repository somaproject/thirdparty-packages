/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/solver/cholesky.hpp
    @author  Jules Bergmann
    @date    2005-09-09
    @brief   VSIPL++ Library: Cholesky Linear system solver.

*/

#ifndef VSIP_CORE_SOLVER_CHOLESKY_HPP
#define VSIP_CORE_SOLVER_CHOLESKY_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/solver/common.hpp>
#ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/cholesky.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_LAPACK
#  include <vsip/opt/lapack/cholesky.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_CVSIP
#  include <vsip/core/cvsip/cholesky.hpp>
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// List of implementation tags to consider for Cholesky.

typedef Make_type_list<
#ifdef VSIP_IMPL_HAVE_CVSIP
  Cvsip_tag,
#endif
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
  >::type Chold_type_list;
  


// a structure to chose implementation type
template <typename T>
struct Choose_chold_impl
{
#ifdef VSIP_IMPL_REF_IMPL
  typedef Cvsip_tag type;
  typedef Cvsip_tag use_type;
#else
  typedef typename Choose_solver_impl<
    Is_chold_impl_avail,
    T,
    Chold_type_list>::type type;

  typedef typename ITE_Type<
    Type_equal<type, None_type>::value,
    As_type<Error_no_solver_for_this_type>,
    As_type<type> >::type use_type;
#endif
};

} // namespace vsip::impl

/// CHOLESKY solver object.

template <typename              T               = VSIP_DEFAULT_VALUE_TYPE,
	  return_mechanism_type ReturnMechanism = by_value>
class chold;

template <typename T>
class chold<T, by_reference>
  : public impl::Chold_impl<T, typename impl::Choose_chold_impl<T>::use_type>
{
  typedef typename impl::Choose_chold_impl<T>::use_type use_type;
  typedef impl::Chold_impl<T, use_type> base_type;

  // Constructors, copies, assignments, and destructors.
public:
  chold(mat_uplo uplo, length_type length)
    VSIP_THROW((std::bad_alloc))
      : base_type(uplo, length)
    { VSIP_IMPL_COVER_TAG("chold", use_type); }

  ~chold() VSIP_NOTHROW {}

  // By-reference solvers.
public:
  template <typename          Block0,
	    typename          Block1>
  bool solve(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->template impl_solve<Block0, Block1>(b, x); }
};

template <typename T>
class chold<T, by_value>
  : public impl::Chold_impl<T, typename impl::Choose_chold_impl<T>::use_type>
{
  typedef impl::Chold_impl<T, typename impl::Choose_chold_impl<T>::use_type>
    base_type;

  // Constructors, copies, assignments, and destructors.
public:
  chold(mat_uplo uplo, length_type length)
    VSIP_THROW((std::bad_alloc))
      : base_type(uplo, length)
    {}

  ~chold() VSIP_NOTHROW {}

  // By-value solvers.
public:
  template <typename Block0>
  Matrix<T>
  solve(const_Matrix<T, Block0> b)
    VSIP_NOTHROW
  {
    Matrix<T> x(b.size(0), b.size(1));
    this->template impl_solve(b, x); 
    return x;
  }
};


} // namespace vsip


#endif // VSIP_OPT_SOLVER_CHOLESKY_HPP
