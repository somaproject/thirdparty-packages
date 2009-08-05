/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/dispatch_assign.hpp
    @author  Jules Bergmann
    @date    2007-10-26
    @brief   VSIPL++ Library: Assignment dispatch declarations.

    Separating declarations here allows eval/diag.hpp to be used
    from headers that are included by dispatch_assign.hpp, for
    example, in eval_fastconv.hpp.

*/

#ifndef VSIP_CORE_DISPATCH_ASSIGN_DECL_HPP
#define VSIP_CORE_DISPATCH_ASSIGN_DECL_HPP

/***********************************************************************
  Included Files
***********************************************************************/



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// Tags used by Dispatch_assign to select assignment implementation.

struct Tag_illegal_mix_of_local_and_global_in_assign;

struct Tag_serial_expr {};
template <typename ParAssignImpl> struct Tag_par_assign {};
struct Tag_par_expr_noreorg {};
struct Tag_par_expr {};



template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2,
	  bool           EarlyBinding>
struct Dispatch_assign_helper;


template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2,
	  typename       Tag
	  = typename Dispatch_assign_helper<Dim, Block1, Block2, false>::type>
struct Dispatch_assign;

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_CORE_DISPATCH_ASSIGN_DECL_HPP
