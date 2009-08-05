/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/eval_reductions.hpp
    @author  Jules Bergmann
    @date    2006-05-30
    @brief   VSIPL++ Library: Reduction functions returning indices.
	     [math.fns.reductidx].

*/

#ifndef VSIP_IMPL_SAL_EVAL_REDUCTIONS_HPP
#define VSIP_IMPL_SAL_EVAL_REDUCTIONS_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/dispatch.hpp>
#include <vsip/core/reductions/functors.hpp>
#include <vsip/opt/sal/eval_util.hpp>
#include <vsip/opt/sal/reductions.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace sal
{

template <template <typename> class ReduceT,
	  typename                  T>
struct Is_reduct_supported
{
  static bool const value = false;
};

#define VSIP_IMPL_REDUCT_SUP(OP, T)					\
template <> struct Is_reduct_supported<OP, T>				\
{ static bool const value = true; };

VSIP_IMPL_REDUCT_SUP(Sum_value, float*)
VSIP_IMPL_REDUCT_SUP(Sum_value, double*)

VSIP_IMPL_REDUCT_SUP(Sum_sq_value, float*)
VSIP_IMPL_REDUCT_SUP(Sum_sq_value, double*)

VSIP_IMPL_REDUCT_SUP(Mean_value, float*)
VSIP_IMPL_REDUCT_SUP(Mean_value, double*)

VSIP_IMPL_REDUCT_SUP(Mean_magsq_value, float*)
VSIP_IMPL_REDUCT_SUP(Mean_magsq_value, double*)

VSIP_IMPL_REDUCT_SUP(Max_value, float*)
VSIP_IMPL_REDUCT_SUP(Max_value, double*)

VSIP_IMPL_REDUCT_SUP(Min_value, float*)
VSIP_IMPL_REDUCT_SUP(Min_value, double*)

#undef VSIP_IMPL_REDUCT_SUP

} // namespace vsip::impl::sal



/***********************************************************************
  Evaluators for reduction functions
***********************************************************************/

namespace dispatcher
{

#define VSIP_IMPL_SAL_REDUCT(OP, SALFCN)                                \
template <typename T,                                                   \
          typename Block>                                               \
struct Evaluator<Op_reduce<OP>, Mercury_sal_tag,                        \
                 void(T&, Block const&, row1_type, Int_type<1>)>        \
{                                                                       \
  typedef typename Block::map_type                        map_type;     \
  typedef typename sal::Effective_value_type<Block>::type eff_t;        \
                                                                        \
  static bool const ct_valid =                                          \
    !Is_expr_block<Block>::value &&                                     \
     Is_local_map<map_type>::value &&                                   \
     sal::Is_reduct_supported<OP, eff_t>::value &&                      \
     Ext_data_cost<Block>::value == 0;                                  \
                                                                        \
  static bool rt_valid(T&, Block const&, row1_type, Int_type<1>)        \
  { return true; }                                                      \
                                                                        \
  static void exec(T& r, Block const& blk, row1_type, Int_type<1>)      \
  {                                                                     \
    sal::Ext_wrapper<Block> ext(blk, SYNC_IN);                          \
                                                                        \
    SALFCN(                                                             \
      typename sal::Ext_wrapper<Block>::sal_type(ext),                  \
      r,                                                                \
      blk.size());                                                      \
  }                                                                     \
};

VSIP_IMPL_SAL_REDUCT(Sum_value,        sal::sumval)
VSIP_IMPL_SAL_REDUCT(Sum_sq_value,     sal::sumsqval)
VSIP_IMPL_SAL_REDUCT(Mean_value,       sal::meanval)
VSIP_IMPL_SAL_REDUCT(Mean_magsq_value, sal::meanmagsqval)



/***********************************************************************
  Evaluators for reduction-idx functions
***********************************************************************/

#define VSIP_IMPL_SAL_REDUCT_IDX(OP, SALFCN)                            \
template <typename T,                                                   \
          typename Block>                                               \
struct Evaluator<Op_reduce_idx<OP>, Mercury_sal_tag,                    \
                 void(T&, Block const&, Index<1>&, row1_type)>          \
{                                                                       \
  typedef typename Block::map_type                        map_type;     \
  typedef typename sal::Effective_value_type<Block>::type eff_t;        \
                                                                        \
  static bool const ct_valid =                                          \
    !Is_expr_block<Block>::value &&                                     \
     Is_local_map<map_type>::value &&                                   \
     sal::Is_reduct_supported<OP, eff_t>::value &&                      \
     Ext_data_cost<Block>::value == 0;                                  \
                                                                        \
  static bool rt_valid(T&, Block const&, Index<1>&, row1_type)          \
  { return true; }                                                      \
                                                                        \
  static void exec(T& r, Block const& blk, Index<1>& idx, row1_type)    \
  {                                                                     \
    sal::Ext_wrapper<Block> ext(blk, SYNC_IN);                          \
                                                                        \
    int i;                                                              \
    SALFCN(                                                             \
      typename sal::Ext_wrapper<Block>::sal_type(ext),                  \
      r,                                                                \
      i,                                                                \
      blk.size());                                                      \
    idx = Index<1>(i);                                                  \
  }                                                                     \
};

VSIP_IMPL_SAL_REDUCT_IDX(Max_value, sal::maxval)
VSIP_IMPL_SAL_REDUCT_IDX(Min_value, sal::minval)

VSIP_IMPL_SAL_REDUCT_IDX(Max_mag_value, sal::maxmgval)
VSIP_IMPL_SAL_REDUCT_IDX(Min_mag_value, sal::minmgval)

#undef VSIP_IMPL_SAL_REDUCT_IDX

} // namespace vsip::impl::dispatcher

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_SAL_REDUCTIONS_HPP
