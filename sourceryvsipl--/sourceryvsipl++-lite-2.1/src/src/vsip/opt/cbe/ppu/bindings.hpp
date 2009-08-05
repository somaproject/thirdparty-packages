/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/bindings.hpp
    @author  Stefan Seefeld
    @date    2006-12-29
    @brief   VSIPL++ Library: Wrappers and traits to bridge with IBMs CBE SDK.
*/

#ifndef VSIP_OPT_CBE_PPU_BINDINGS_HPP
#define VSIP_OPT_CBE_PPU_BINDINGS_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/adjust_layout.hpp>
#include <vsip/opt/cbe/vmmul_params.h>
#include <vsip/opt/cbe/ppu/task_manager.hpp>
#include <vsip/opt/cbe/ppu/util.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cbe
{

template <typename T> void vmul(T const* A, T const* B, T* R, length_type len);
template <typename T> void vmul(std::pair<T*, T*> const& A,
				std::pair<T*, T*> const& B,
				std::pair<T*, T*> const& R,
				length_type              len);
template <typename T> void vadd(T const* A, T const* B, T* R, length_type len);
template <typename T> void vadd(std::pair<T*, T*> const& A,
				std::pair<T*, T*> const& B,
				std::pair<T*, T*> const& R,
				length_type              len);
template <typename T> void
vmmul_row(T const* V, T const* M, T* R, 
      stride_type m_stride, stride_type r_stride,
      length_type length, length_type lines);

template <typename T>
void vmmul_row(
  T const* V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const&             R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);

template <typename T>
void vmmul_row(
  std::pair<T*, T*> const& V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const&             R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);

template <typename T> void
vmmul_col(T const* V, T const* M, T* R, 
      stride_type m_stride, stride_type r_stride,
      length_type length, length_type lines);

template <typename T>
void vmmul_col(
  T const* V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const&             R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);

template <typename T>
void vmmul_col(
  std::pair<T*, T*> const& V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const&             R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);


template <template <typename, typename> class Operator,
	  typename DstBlock,
	  typename LBlock,
	  typename RBlock,
	  typename LType,
	  typename RType>
struct Serial_expr_evaluator_base
{
  typedef Binary_expr_block<1, Operator, LBlock, LType, RBlock, RType>
    SrcBlock;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<LBlock>::layout_type>::type
    lblock_lp;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<RBlock>::layout_type>::type
    rblock_lp;

  static bool const ct_valid = 
    !Is_expr_block<LBlock>::value &&
    !Is_expr_block<RBlock>::value &&
     Is_split_block<DstBlock>::value == Is_split_block<LBlock>::value &&
     Is_split_block<DstBlock>::value == Is_split_block<RBlock>::value &&
     Type_equal<typename DstBlock::value_type, LType>::value &&
     Type_equal<typename DstBlock::value_type, RType>::value &&
     // check that type is supported.
     (Type_equal<typename DstBlock::value_type, float>::value ||
      Type_equal<typename DstBlock::value_type, complex<float> >::value) &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<LBlock>::value == 0 &&
     Ext_data_cost<RBlock>::value == 0;

  static length_type tunable_threshold()
  {
    typedef typename DstBlock::value_type T;

    if (VSIP_IMPL_TUNE_MODE)
      return 0;
    // Compare interleaved vmul -2 --svpp-num-spes {0,8}.
    else if (Type_equal<Operator<T, T>,
	     op::Mult<complex<float>, complex<float> > >::value)
      return 16384;
    // Compare vmul -1 --svpp-num-spes {0,8}.
    else if (Type_equal<Operator<T, T>, op::Mult<float, float> >::value)
      return 65536;
    else
      return 0;
  }

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>    ext_dst(dst,       SYNC_OUT);
    Ext_data<LBlock,   lblock_lp> ext_l(src.left(),  SYNC_IN);
    Ext_data<RBlock,   rblock_lp> ext_r(src.right(), SYNC_IN);
    return ext_dst.size(0) >= tunable_threshold() &&
           ext_dst.stride(0) == 1 &&
	   ext_l.stride(0) == 1   &&
	   ext_r.stride(0) == 1   &&
	   is_dma_addr_ok(ext_dst.data()) &&
	   is_dma_addr_ok(ext_l.data())   &&
	   is_dma_addr_ok(ext_r.data())   &&
           Task_manager::instance()->num_spes() > 0;
  }
};

} // namespace vsip::impl::cbe

/***********************************************************************
  Binary expression evaluators
***********************************************************************/

#define VSIP_IMPL_CBE_PPU_VV_EXPR(OP, FUN)				\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename RBlock,						\
	  typename LType,						\
	  typename RType>						\
struct Serial_expr_evaluator<						\
    1, DstBlock, 							\
    Binary_expr_block<1, OP, LBlock, LType, RBlock, RType> const,	\
    Cbe_sdk_tag>							\
  : cbe::Serial_expr_evaluator_base<OP, DstBlock,			\
				    LBlock, RBlock, LType, RType>	\
{									\
  static char const* name() { return "Expr_CBE_SDK_VV-" #FUN; }		\
									\
  typedef Binary_expr_block<1, OP, LBlock, LType, RBlock, RType>	\
    SrcBlock;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<LBlock>::layout_type>::type		\
    lblock_lp;								\
  									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<RBlock>::layout_type>::type		\
    rblock_lp;								\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);		\
    Ext_data<LBlock, lblock_lp> ext_l(src.left(),  SYNC_IN);		\
    Ext_data<RBlock, rblock_lp> ext_r(src.right(), SYNC_IN);		\
									\
    FUN(								\
      ext_l.data(),							\
      ext_r.data(),							\
      ext_dst.data(),							\
      dst.size()							\
    );									\
  } 									\
};

VSIP_IMPL_CBE_PPU_VV_EXPR(op::Add,  cbe::vadd)
// VSIP_IMPL_CBE_PPU_VV_EXPR(op::Sub,  cbe::vsub)
VSIP_IMPL_CBE_PPU_VV_EXPR(op::Mult,  cbe::vmul)
// VSIP_IMPL_CBE_PPU_VV_EXPR(op::Div,  cbe::vdiv)

#undef VSIP_IMPL_CBE_PPU_VV_EXPR




/// Evaluator for vector-matrix multiply.

/// Dispatches cases where the dimension ordering matches the 
/// requested orientation to the SPU's (row-major/by-row and 
/// col-major/by-col).  The other cases are re-dispatched.

template <typename DstBlock,
	  typename VBlock,
	  typename MBlock,
	  dimension_type SD>
struct Serial_expr_evaluator<2, DstBlock,
			     const Vmmul_expr_block<SD, VBlock, MBlock>,
			     Cbe_sdk_tag>
{
  static char const* name() { return "Cbe_Sdk_Vmmul"; }

  typedef Vmmul_expr_block<SD, VBlock, MBlock> SrcBlock;
  typedef typename DstBlock::value_type dst_type;
  typedef typename VBlock::value_type v_type;
  typedef typename MBlock::value_type m_type;
  typedef typename Block_layout<DstBlock>::layout_type dst_lp;
  typedef typename Block_layout<VBlock>::layout_type vblock_lp;
  typedef typename Block_layout<MBlock>::layout_type mblock_lp;
  typedef typename Block_layout<DstBlock>::order_type order_type;
  typedef typename Block_layout<MBlock>::order_type src_order_type;

  static bool const is_row_vmmul =
    (SD == row && Type_equal<order_type, row2_type>::value ||
     SD == col && Type_equal<order_type, col2_type>::value);

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
    !Is_expr_block<MBlock>::value &&
    (Is_split_block<DstBlock>::value == Is_split_block<VBlock>::value ||
     Type_equal<v_type,   float>::value) &&
    (Is_split_block<DstBlock>::value == Is_split_block<MBlock>::value) &&
    (( is_row_vmmul && Type_equal<v_type, std::complex<float> >::value) ||
     (!is_row_vmmul && (Type_equal<v_type,   float>::value ||
			Type_equal<v_type,   std::complex<float> >::value))) &&
     // At this time, only complex-complex vmmul operations are supported.
     Type_equal<v_type,   std::complex<float> >::value &&
     Type_equal<m_type,   std::complex<float> >::value &&
     Type_equal<dst_type, std::complex<float> >::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<VBlock>::value == 0 &&
     Ext_data_cost<MBlock>::value == 0 &&
     Type_equal<order_type, src_order_type>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    VBlock const& vblock = src.get_vblk();
    MBlock const& mblock = src.get_mblk();

    Ext_data<DstBlock, dst_lp>  ext_dst(dst,    SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_v(vblock, SYNC_IN);
    Ext_data<MBlock, mblock_lp> ext_m(mblock, SYNC_IN);

    if (is_row_vmmul)
    {
      dimension_type const axis = SD == row ? 1 : 0;
      length_type size = dst.size(2, axis);

      return 
	// (large sizes are broken down)
	(size >= VSIP_IMPL_MIN_VMMUL_SIZE) && 
	(ext_dst.stride(axis) == 1) &&
	(ext_m.stride(axis)   == 1) &&
	(ext_v.stride(0) == 1) &&
	cbe::is_dma_addr_ok(ext_dst.data()) &&
	cbe::is_dma_addr_ok(ext_v.data()) &&
	cbe::is_dma_addr_ok(ext_m.data()) &&
	cbe::is_dma_stride_ok<dst_type>(ext_dst.stride(SD == row ? 0 : 1)) &&
	cbe::is_dma_stride_ok<m_type>(ext_m.stride(SD == row ? 0 : 1)) &&
	// (non-granular sizes handled)
	cbe::Task_manager::instance()->num_spes() > 0;
    }
    else
    {
      dimension_type const axis = SD == row ? 0 : 1;
      length_type size = dst.size(2, axis);

      return 
	// (large sizes are broken down)
	(size >= VSIP_IMPL_MIN_VMMUL_SIZE) && 
	(ext_dst.stride(axis) == 1) &&
	(ext_m.stride(axis)   == 1) &&
	(ext_v.stride(0) == 1) &&
	cbe::is_dma_addr_ok(ext_dst.data()) &&
	// (V doesn't need to be DMA aligned)
	cbe::is_dma_addr_ok(ext_m.data()) &&
	cbe::is_dma_stride_ok<dst_type>(ext_dst.stride(SD == row ? 1 : 0)) &&
	cbe::is_dma_stride_ok<m_type>(ext_m.stride(SD == row ? 1 : 0)) &&
	cbe::is_dma_size_ok(size * sizeof(v_type)) &&
	cbe::Task_manager::instance()->num_spes() > 0;
    }
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    VBlock const& vblock = src.get_vblk();
    MBlock const& mblock = src.get_mblk();

    Matrix<dst_type, DstBlock> m_dst(dst);
    const_Vector<dst_type, VBlock>  v(const_cast<VBlock&>(vblock));
    const_Matrix<dst_type, MBlock>  m(const_cast<MBlock&>(mblock));

    Ext_data<DstBlock, dst_lp>  ext_dst(dst,    SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_v(vblock, SYNC_IN);
    Ext_data<MBlock, mblock_lp> ext_m(mblock, SYNC_IN);

    // The ct_valid check above ensures that the order taken 
    // matches the storage order if reaches this point.
    if (SD == row && Type_equal<order_type, row2_type>::value)
    {
      cbe::vmmul_row(
        ext_v.data(),
        ext_m.data(),
        ext_dst.data(),
        ext_m.stride(0),   // elements between rows of source matrix
        ext_dst.stride(0), // elements between rows of destination matrix
        dst.size(2, 0),    // number of rows
        dst.size(2, 1)     // length of each row
      );
    }
    else if (SD == col && Type_equal<order_type, row2_type>::value)
    {
      cbe::vmmul_col(
        ext_v.data(),
        ext_m.data(),
        ext_dst.data(),
        ext_m.stride(0),   // elements between rows of source matrix
        ext_dst.stride(0), // elements between rows of destination matrix
        dst.size(2, 0),    // number of rows
        dst.size(2, 1)     // length of each row
      );
    }
    else if (SD == col && Type_equal<order_type, col2_type>::value)
    {
      cbe::vmmul_row(
        ext_v.data(),
        ext_m.data(),
        ext_dst.data(),
        ext_m.stride(1),   // elements between cols of source matrix
        ext_dst.stride(1), // elements between cols of destination matrix
        dst.size(2, 1),    // number of cols
        dst.size(2, 0)     // length of each col
      );
    }
    else // if (SD == row && Type_equal<order_type, col2_type>::value)
    {
      cbe::vmmul_col(
        ext_v.data(),
        ext_m.data(),
        ext_dst.data(),
        ext_m.stride(1),   // elements between cols of source matrix
        ext_dst.stride(1), // elements between cols of destination matrix
        dst.size(2, 1),    // number of cols
        dst.size(2, 0)     // length of each col
      );
    }
  }
};



} // namespace vsip::impl
} // namespace vsip

#endif
