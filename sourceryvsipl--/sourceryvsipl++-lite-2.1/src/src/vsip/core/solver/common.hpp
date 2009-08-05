/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/solver/common.hpp
    @author  Assem Salama
    @date    2005-04-13
    @brief   VSIPL++ Library: Common stuff for linear system solvers.

*/

#ifndef VSIP_CORE_SOLVER_COMMON_HPP
#define VSIP_CORE_SOLVER_COMMON_HPP

#include <vsip/core/impl_tags.hpp>
#include <vsip/core/type_list.hpp>

namespace vsip
{
namespace impl
{

// Structures for availability
template <typename   ImplTag,
          typename   T>
struct Is_lud_impl_avail
{
  static bool const value = false;
};

template <typename   ImplTag,
          typename   T>
struct Is_chold_impl_avail
{
  static bool const value = false;
};

template <typename   ImplTag,
          typename   T>
struct Is_qrd_impl_avail
{
  static bool const value = false;
};

template <typename   ImplTag,
          typename   T>
struct Is_svd_impl_avail
{
  static bool const value = false;
};



// LUD solver impl class
template <typename T,
          typename ImplTag>
class Lud_impl;

// CHOLESKY solver impl class
template <typename T,
          typename ImplTag>
class Chold_impl;

// QR solver impl class
template <typename T,
	  bool     Blocked,
	  typename ImplTag>
class Qrd_impl;

template <typename T,
	  bool     Blocked,
	  typename ImplTag>
class Svd_impl;



// Error tags
struct Error_no_solver_for_this_type;


} // namespace vsip::impl

// Common enums
enum mat_uplo
{
  lower,
  upper
};



namespace impl
{

template <typename List>
struct List_get
{
  typedef typename List::first first;
  typedef typename List::rest  rest;
};

template <>
struct List_get<None_type>
{
  typedef None_type first;
  typedef None_type rest;
};

/// Template class to determine which tag implements a solver.

/// Requires:
///   ISTYPEAVAIL to be a template class that, given an ImplTag and
///      value type, defines VALUE to be true if ImplTag can solve
///      value type, else false.
///   T is a value type.
///   TAGLIST is a list of implementation tags.
///
/// Provides:
///   TYPE to be the first implementation tag from TAGLIST that supports
///      value type T, else None_type.

template <template <typename, typename> class IsTypeAvail,
	  typename T,
	  typename TagList,
	  typename Tag  = typename List_get<TagList>::first,
	  typename Rest = typename List_get<TagList>::rest,
	  bool     Valid = IsTypeAvail<Tag, T>::value>
struct Choose_solver_impl;

/// Specialization for case where impl tag TAG supports type T.
template <template <typename, typename> class IsTypeAvail,
	  typename T,
	  typename TagList,
	  typename Tag,
	  typename Rest>
struct Choose_solver_impl<IsTypeAvail, T, TagList, Tag, Rest, true>
{
  typedef Tag type;
};

/// Specialization for case where impl tag TAG does not support type T.
/// Fall through to next entry in TAGLIST.
template <template <typename, typename> class IsTypeAvail,
	  typename T,
	  typename TagList,
	  typename Tag,
	  typename Rest>
struct Choose_solver_impl<IsTypeAvail, T, TagList, Tag, Rest, false>
  : Choose_solver_impl<IsTypeAvail, T, Rest>
{};

/// Terminator.  If REST is empty, define type to be None_type.
template <template <typename, typename> class IsTypeAvail,
	  typename T,
	  typename TagList,
	  typename Tag>
struct Choose_solver_impl<IsTypeAvail, T, TagList, Tag, None_type, false>
{
  typedef None_type type;
};


/// Special terminator.  If original list is empty, define type to
/// be None_type.

template <template <typename, typename> class IsTypeAvail,
	  typename T>
struct Choose_solver_impl<IsTypeAvail, T, None_type, None_type, None_type, false>
{
  typedef None_type type;
};

} // namespace vsip::impl
} // namespace vsip

#endif
