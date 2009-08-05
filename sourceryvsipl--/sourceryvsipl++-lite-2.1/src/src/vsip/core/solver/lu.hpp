/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/solver/lu.hpp
    @author  Jules Bergmann
    @date    2005-09-29
    @brief   VSIPL++ Library: LU linear system solver.

*/

#ifndef VSIP_CORE_SOLVER_LU_HPP
#define VSIP_CORE_SOLVER_LU_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/metaprogramming.hpp>
#  include <vsip/core/solver/common.hpp>
#ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/lu.hpp>
#endif
#ifdef VSIP_IMPL_CBE_SDK
#  include <vsip/opt/cbe/cml/lu.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_LAPACK
#  include <vsip/opt/lapack/lu.hpp>
#endif
#ifdef VSIP_IMPL_HAVE_CVSIP
#  include <vsip/core/cvsip/lu.hpp>
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
  None_type // None_type is treated specially by Make_type_list, it is
            // not put into the list.  Putting an explicit None_type
            // at the end of the list lets us put a ',' after each impl
            // tag.
  >::type Lud_type_list;



// a structure to chose implementation type
template <typename T>
struct Choose_lud_impl
{
#ifdef VSIP_IMPL_REF_IMPL
  typedef Cvsip_tag use_type;
  typedef Cvsip_tag type;
#else
  typedef typename Choose_solver_impl<
    Is_lud_impl_avail,
    T,
    Lud_type_list>::type type;

  typedef typename ITE_Type<
    Type_equal<type, None_type>::value,
    As_type<Error_no_solver_for_this_type>,
    As_type<type> >::type use_type;
#endif
};

} // namespace impl

/// LU solver object.

template <typename              T               = VSIP_DEFAULT_VALUE_TYPE,
	  return_mechanism_type ReturnMechanism = by_value>
class lud;



/// LU solver object (by-reference).

template <typename T>
class lud<T, by_reference>
  : public impl::Lud_impl<T,typename impl::Choose_lud_impl<T>::use_type>
{
  typedef typename impl::Choose_lud_impl<T>::use_type use_type;
  typedef impl::Lud_impl<T, use_type> base_type;

  // Constructors, copies, assignments, and destructors.
public:
  lud(length_type length)
    VSIP_THROW((std::bad_alloc))
      : base_type(length)
    { VSIP_IMPL_COVER_TAG("lud", use_type); }

  ~lud() VSIP_NOTHROW {}

  // By-reference solvers.
public:
  template <mat_op_type tr,
	    typename    Block0,
	    typename    Block1>
  bool solve(const_Matrix<T, Block0> b, Matrix<T, Block1> x)
    VSIP_NOTHROW
  { return this->template impl_solve<tr, Block0, Block1>(b, x); }
};



/// LU solver object (by-value).

template <typename T>
class lud<T, by_value>
  : public impl::Lud_impl<T,typename impl::Choose_lud_impl<T>::use_type>
{
  typedef impl::Lud_impl<T,typename impl::Choose_lud_impl<T>::use_type>
	  base_type;

  // Constructors, copies, assignments, and destructors.
public:
  lud(length_type length)
    VSIP_THROW((std::bad_alloc))
      : base_type(length)
    {}

  ~lud() VSIP_NOTHROW {}

  // By-value solvers.
public:
  template <mat_op_type tr,
	    typename    Block0>
  Matrix<T>
  solve(const_Matrix<T, Block0> b)
    VSIP_NOTHROW
  {
    Matrix<T> x(b.size(0), b.size(1));
    this->template impl_solve<tr>(b, x); 
    return x;
  }
};



} // namespace vsip


#endif // VSIP_CORE_SOLVER_LU_HPP
