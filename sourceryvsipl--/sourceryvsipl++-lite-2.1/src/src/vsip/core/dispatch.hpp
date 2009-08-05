/* Copyright (c) 2008 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/dispatch.hpp
    @author  Don McCoy
    @date    2008-11-17
    @brief   VSIPL++ Library: Dispatcher harness basic definitions (see
               vsip/opt/dispatch.hpp for the actual dispatcher).
*/

#ifndef VSIP_CORE_DISPATCH_HPP
#define VSIP_CORE_DISPATCH_HPP


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace dispatcher
{


/// Define the operation-specific Evaluator signature.
template <typename O, typename R = void> 
struct Signature
{
  // The default signature is useful for a compile-time check only,
  // as there are no arguments to inspect at runtime.
  typedef R(type)();
};


/// An Evaluator determines whether an Operation can be performed
/// with a particular backend.
///
///   O: Operation tag
///   B: Backend tag
///   S: Signature
template <typename O,
          typename B,
          typename S = typename Signature<O>::type>
struct Evaluator
{
  static bool const ct_valid = false;
};


} // namespace vsip::impl::dispatcher
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_CORE_DISPATCH_HPP
