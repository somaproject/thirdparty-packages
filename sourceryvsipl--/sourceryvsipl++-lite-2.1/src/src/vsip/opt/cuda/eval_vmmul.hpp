/* Copyright (c) 2009 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/cuda/eval_vmmul.hpp
    @author  Don McCoy
    @date    2009-04-07
    @brief   VSIPL++ Library: CUDA evaluator for vector-matrix multiply.

*/

#ifndef VSIP_OPT_CUDA_EVAL_VMMUL_HPP
#define VSIP_OPT_CUDA_EVAL_VMMUL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/opt/cuda/device_memory.hpp>
#include <vsip/opt/cuda/kernels.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Evaluator for vector-matrix multiply.
///
/// Dispatches cases where the dimension ordering matches the 
/// requested orientation to the SPU's (row-major/by-row and 
/// col-major/by-col).  The other cases are re-dispatched.
template <typename DstBlock,
	  typename VBlock,
	  typename MBlock,
	  dimension_type SD>
struct Serial_expr_evaluator<2, DstBlock,
			     const Vmmul_expr_block<SD, VBlock, MBlock>,
			     Cuda_tag>
{
  static char const* name() { return "CUDA_vmmul"; }

  typedef Vmmul_expr_block<SD, VBlock, MBlock> SrcBlock;
  typedef typename SrcBlock::value_type src_type;
  typedef typename DstBlock::value_type dst_type;
  typedef typename VBlock::value_type v_type;
  typedef typename MBlock::value_type m_type;
  typedef typename Block_layout<DstBlock>::layout_type dst_lp;
  typedef typename Block_layout<VBlock>::layout_type vblock_lp;
  typedef typename Block_layout<MBlock>::layout_type mblock_lp;
  typedef typename Block_layout<DstBlock>::order_type order_type;
  typedef typename Block_layout<MBlock>::order_type src_order_type;

  static bool const ct_valid = 
    // inputs must not be an expression blocks
    Is_expr_block<VBlock>::value == 0 &&
    Is_expr_block<MBlock>::value == 0 &&
    // split complex not supported
    Is_split_block<DstBlock>::value == 0 &&
    Is_split_block<VBlock>::value == 0 &&
    Is_split_block<MBlock>::value == 0 &&
    // ensure value types are supported
    cuda::Cuda_traits<dst_type>::valid &&
    cuda::Cuda_traits<v_type>::valid &&
    cuda::Cuda_traits<m_type>::valid &&
    // result type must match type expected (determined by promotion)
    Type_equal<dst_type, src_type>::value &&
    // check that direct access is supported
    Ext_data_cost<DstBlock>::value == 0 &&
    Ext_data_cost<VBlock>::value == 0 &&
    Ext_data_cost<MBlock>::value == 0 &&
    // dimension ordering must be the same
    Type_equal<order_type, src_order_type>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    VBlock const& vblock = src.get_vblk();
    MBlock const& mblock = src.get_mblk();

    Ext_data<DstBlock, dst_lp>  ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_v(vblock, SYNC_IN);
    Ext_data<MBlock, mblock_lp> ext_m(mblock, SYNC_IN);

    if (SD == row && Type_equal<order_type, row2_type>::value ||
	SD == col && Type_equal<order_type, col2_type>::value)
    {
      dimension_type const axis = SD == row ? 1 : 0;
      length_type dst_stride = static_cast<length_type>(abs(ext_dst.stride(axis == 0)));
      length_type m_stride = static_cast<length_type>(abs(ext_m.stride(axis == 0)));
      return 
        // make sure blocks are dense (major stride == minor size)
        (ext_dst.size(axis) == dst_stride) &&
        (ext_m.size(axis) == m_stride) &&
        // ensure unit stride along the dimension opposite the chosen one
	(ext_dst.stride(axis) == 1) &&
	(ext_m.stride(axis) == 1) &&
	(ext_v.stride(0) == 1);
    }
    else
    {
      dimension_type const axis = SD == row ? 0 : 1;
      length_type dst_stride = static_cast<length_type>(abs(ext_dst.stride(axis == 0)));
      length_type m_stride = static_cast<length_type>(abs(ext_m.stride(axis == 0)));
      return 
        // make sure blocks are dense (major stride == minor size)
        (ext_dst.size(axis) == dst_stride) &&
        (ext_m.size(axis) == m_stride) &&
        // ensure unit stride along the same dimension as the chosen one
	(ext_dst.stride(axis) == 1) &&
	(ext_m.stride(axis) == 1) &&
	(ext_v.stride(0) == 1);
    }
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    VBlock const& vblock = src.get_vblk();
    MBlock const& mblock = src.get_mblk();

    cuda::Device_memory<DstBlock> dev_dst(dst, SYNC_OUT);
    cuda::Device_memory<VBlock const> dev_v(vblock);
    cuda::Device_memory<MBlock const> dev_m(mblock);

    // The ct_valid check above ensures that the order taken 
    // matches the storage order if reaches this point.
    if (SD == row && Type_equal<order_type, row2_type>::value)
    {
      cuda::vmmul_row(
        dev_v.data(),
        dev_m.data(),
        dev_dst.data(),
        dst.size(2, 0),    // number of rows
        dst.size(2, 1)     // length of each row
      );
    }
    else if (SD == col && Type_equal<order_type, row2_type>::value)
    {
      cuda::vmmul_col(
        dev_v.data(),
        dev_m.data(),
        dev_dst.data(),
        dst.size(2, 0),    // number of rows
        dst.size(2, 1)     // length of each row
      );
    }
    else if (SD == col && Type_equal<order_type, col2_type>::value)
    {
      cuda::vmmul_row(
        dev_v.data(),
        dev_m.data(),
        dev_dst.data(),
        dst.size(2, 1),    // number of cols
        dst.size(2, 0)     // length of each col
      );
    }
    else // if (SD == row && Type_equal<order_type, col2_type>::value)
    {
      cuda::vmmul_col(
        dev_v.data(),
        dev_m.data(),
        dev_dst.data(),
        dst.size(2, 1),    // number of cols
        dst.size(2, 0)     // length of each col
      );
    }
  }
};


} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_EVAL_VMMUL_HPP
