/* Copyright (c) 2006, 2007, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/simd/expr_evaluator.hpp
    @author  Stefan Seefeld
    @date    2006-07-25
    @brief   VSIPL++ Library: SIMD expression evaluator logic.

*/

#ifndef VSIP_IMPL_SIMD_EXPR_EVALUATOR_HPP
#define VSIP_IMPL_SIMD_EXPR_EVALUATOR_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/opt/simd/simd.hpp>
#include <vsip/opt/simd/expr_iterator.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/opt/simd/proxy_factory.hpp>

/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{

// SIMD Loop Fusion evaluator for aligned expressions.

template <typename LB,
	  typename RB>
struct Serial_expr_evaluator<1, LB, RB, Simd_loop_fusion_tag>
{
  typedef typename Adjust_layout_dim<
                     1, typename Block_layout<LB>::layout_type>::type
		layout_type;

  static char const* name() { return "Expr_SIMD_Loop"; }
  
  static bool const ct_valid =
    // Is SIMD supported at all ?
    simd::Simd_traits<typename LB::value_type>::is_accel &&
    // Check that direct access is possible.
    Ext_data_cost<LB>::value == 0 &&
    simd::Proxy_factory<RB, true>::ct_valid &&
    // Only allow float, double, complex<float>, and complex<double> at this time.
    (Type_equal<typename Scalar_of<typename LB::value_type>::type, float>::value ||
     Type_equal<typename Scalar_of<typename LB::value_type>::type, double>::value) &&
    // Make sure both sides have the same type.
    Type_equal<typename LB::value_type, typename RB::value_type>::value &&
    // Make sure the left side is not a complex split block.
    !Is_split_block<LB>::value;


  static bool rt_valid(LB& lhs, RB const& rhs)
  {
    Ext_data<LB, layout_type> dda(lhs, SYNC_OUT);
    return (dda.stride(0) == 1 &&
      simd::Proxy_factory<RB, true>::rt_valid(rhs, 
      simd::Simd_traits<typename LB::value_type>::alignment_of(dda.data())));
  }

  static void exec(LB& lhs, RB const& rhs)
  {
    typedef typename simd::LValue_access_traits<typename LB::value_type> WAT;
    typedef typename simd::Proxy_factory<RB, true>::access_traits EAT;
    typedef typename simd::Proxy_factory<RB, true>::proxy_type proxy_type;

    length_type const vec_size =
      simd::Simd_traits<typename LB::value_type>::vec_size;
    Ext_data<LB, layout_type> dda(lhs, SYNC_OUT);

    simd::Proxy<WAT,true> lp(dda.data());
    proxy_type rp(simd::Proxy_factory<RB,true>::create(rhs));

    length_type const size = dda.size(0);
    length_type n = size;

    // loop using proxy interface. This generates the best code
    // with gcc 3.4 (with gcc 4.1 the difference to the first case
    // above is negligible).

    // First, deal with unaligned pointers
    typename Ext_data<LB, layout_type>::raw_ptr_type  raw_ptr = dda.data();
    while(simd::Simd_traits<typename LB::value_type>::alignment_of(raw_ptr) &&
          n > 0)
    {
      lhs.put(size-n, rhs.get(size-n));
      n--;
      raw_ptr++;
      lp.increment_by_element(1);
      rp.increment_by_element(1);
    }

    while (n >= vec_size)
    {
      lp.store(rp.load());
      n -= vec_size;
      lp.increment();
      rp.increment();
    }

    // Process the remainder, using simple loop fusion.
    for (index_type i = size - n; i != size; ++i)
      lhs.put(i, rhs.get(i));
  }
};

} // namespace vsip::impl
} // namespace vsip

#endif
