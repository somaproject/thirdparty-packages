/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/expr/serial_evaluator.hpp
    @author  Stefan Seefeld
    @date    2005-08-10
    @brief   VSIPL++ Library: Helper templates for expression evaluation.
*/

#ifndef VSIP_OPT_EXPR_SERIAL_EVALUATOR_HPP
#define VSIP_OPT_EXPR_SERIAL_EVALUATOR_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/adjust_layout.hpp>
#include <vsip/core/coverage.hpp>
#include <vsip/core/impl_tags.hpp>
#include <vsip/opt/expr/lf_initfini.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

/// Serial_expr_evaluator template.
/// This needs to be provided for each tag in the LibraryTagList.
template <dimension_type Dim,
	  typename DstBlock,
	  typename SrcBlock,
	  typename Tag>
struct Serial_expr_evaluator
{
  static char const* name() { return "Expr"; }
  static bool const ct_valid = false;
};

template <typename DstBlock,
	  typename SrcBlock>
struct Serial_expr_evaluator<1, DstBlock, SrcBlock, Loop_fusion_tag>
{
  static char const* name() { return "Expr_Loop"; }
  static bool const ct_valid = true;
  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  { return true; }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    do_loop_fusion_init(src);
    length_type const size = dst.size(1, 0);
    for (index_type i=0; i<size; ++i)
      dst.put(i, src.get(i));
    do_loop_fusion_fini(src);
  }
};



template <typename DstBlock,
	  typename SrcBlock>
struct Serial_expr_evaluator<1, DstBlock, SrcBlock, Copy_tag>
{
  static char const* name() { return "Expr_Copy"; }

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<SrcBlock>::layout_type>::type
    src_lp;

  static bool const ct_valid = 
    Ext_data_cost<DstBlock>::value == 0 &&
    Ext_data_cost<SrcBlock>::value == 0 &&
    !Is_split_block<SrcBlock>::value &&
    !Is_split_block<DstBlock>::value;

  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  { return true; }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<DstBlock, dst_lp> ext_dst(dst, impl::SYNC_OUT);
    Ext_data<SrcBlock, src_lp> ext_src(src, impl::SYNC_IN);

    typename Ext_data<DstBlock, dst_lp>::raw_ptr_type ptr1 = ext_dst.data();
    typename Ext_data<SrcBlock, src_lp>::raw_ptr_type ptr2 = ext_src.data();

    stride_type stride1 = ext_dst.stride(0);
    stride_type stride2 = ext_src.stride(0);
    length_type size    = ext_dst.size(0);
    assert(size <= ext_src.size(0));

    if (Type_equal<typename DstBlock::value_type,
	           typename SrcBlock::value_type>::value &&
	stride1 == 1 && stride2 == 1)
    {
      memcpy(ptr1, ptr2, size*sizeof(typename DstBlock::value_type));
    }
    else
    {
      while (size--)
      {
	*ptr1 = *ptr2;
	ptr1 += stride1;
	ptr2 += stride2;
      }
    }
  }
};






template <typename DstBlock,
	  typename SrcBlock>
struct Serial_expr_evaluator<2, DstBlock, SrcBlock, Loop_fusion_tag>
{
  static char const* name() { return "Expr_Loop"; }
  static bool const ct_valid = true;
  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  { return true; }

  static void exec(DstBlock& dst, SrcBlock const& src, row2_type)
  {
    length_type const rows = dst.size(2, 0);
    length_type const cols = dst.size(2, 1);
    for (index_type i=0; i<rows; ++i)
      for (index_type j=0; j<cols; ++j)
	dst.put(i, j, src.get(i, j));
  }

  static void exec(DstBlock& dst, SrcBlock const& src, col2_type)
  {
    length_type const rows = dst.size(2, 0);
    length_type const cols = dst.size(2, 1);
    for (index_type j=0; j<cols; ++j)
      for (index_type i=0; i<rows; ++i)
	dst.put(i, j, src.get(i, j));
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    typedef typename Block_layout<DstBlock>::order_type dst_order_type;
    do_loop_fusion_init(src);
    exec(dst, src, dst_order_type());
    do_loop_fusion_fini(src);
  }
};



template <typename DstBlock,
	  typename SrcBlock>
struct Serial_expr_evaluator<3, DstBlock, SrcBlock, Loop_fusion_tag>
{
  static char const* name() { return "Expr_Loop"; }
  static bool const ct_valid = true;
  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  { return true; }

  static void exec(DstBlock& dst, SrcBlock const& src, tuple<0,1,2>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type i=0; i<size0; ++i)
    for (index_type j=0; j<size1; ++j)
    for (index_type k=0; k<size2; ++k)
      dst.put(i, j, k, src.get(i, j, k));
  }
  static void exec(DstBlock& dst, SrcBlock const& src, tuple<0,2,1>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type i=0; i<size0; ++i)
    for (index_type k=0; k<size2; ++k)
    for (index_type j=0; j<size1; ++j)
      dst.put(i, j, k, src.get(i, j, k));
  }
  static void exec(DstBlock& dst, SrcBlock const& src, tuple<1,0,2>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type j=0; j<size1; ++j)
    for (index_type i=0; i<size0; ++i)
    for (index_type k=0; k<size2; ++k)
      dst.put(i, j, k, src.get(i, j, k));
  }
  static void exec(DstBlock& dst, SrcBlock const& src, tuple<1,2,0>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type j=0; j<size1; ++j)
    for (index_type k=0; k<size2; ++k)
    for (index_type i=0; i<size0; ++i)
      dst.put(i, j, k, src.get(i, j, k));
  }
  static void exec(DstBlock& dst, SrcBlock const& src, tuple<2,0,1>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type k=0; k<size2; ++k)
    for (index_type i=0; i<size0; ++i)
    for (index_type j=0; j<size1; ++j)
      dst.put(i, j, k, src.get(i, j, k));
  }
  static void exec(DstBlock& dst, SrcBlock const& src, tuple<2,1,0>)
  {
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type k=0; k<size2; ++k)
    for (index_type j=0; j<size1; ++j)
    for (index_type i=0; i<size0; ++i)
      dst.put(i, j, k, src.get(i, j, k));
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    typedef typename Block_layout<DstBlock>::order_type dst_order_type;
    do_loop_fusion_init(src);
    exec(dst, src, dst_order_type());
    do_loop_fusion_fini(src);
  }
};



/// A general expression evaluator for IPP that doesn't match
/// anything and thus should be skipped by the dispatcher.
template <typename DstBlock,
	  typename SrcBlock>
struct Serial_expr_evaluator<1, DstBlock, SrcBlock, Intel_ipp_tag>
{
  static bool const ct_valid = false;
  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/) 
  { return false;}
  static void exec(DstBlock& /*dst*/, SrcBlock const& /*src*/) {}
};

} // namespace vsip::impl
} // namespace vsip

#endif
