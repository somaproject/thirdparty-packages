/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/vmmul.hpp
    @author  Jules Bergmann
    @date    2005-08-15
    @brief   VSIPL++ Library: vector-matrix multiply

*/

#ifndef VSIP_CORE_VMMUL_HPP
#define VSIP_CORE_VMMUL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/block_traits.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/promote.hpp>
#include <vsip/core/expr/vmmul_block.hpp>
#if !VSIP_IMPL_REF_IMPL
#  include <vsip/opt/expr/serial_evaluator.hpp>
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Traits class to determines return type for vmmul.

template <dimension_type Dim,
	  typename       T0,
	  typename       T1,
	  typename       Block0,
	  typename       Block1>
struct Vmmul_traits
{
  typedef typename vsip::Promotion<T0, T1>::type      value_type;
  typedef const Vmmul_expr_block<Dim, Block0, Block1> block_type;
  typedef Matrix<value_type, block_type>              view_type;
};



#if !VSIP_IMPL_REF_IMPL
/// Evaluator for vector-matrix multiply.

/// Reduces vmmul into either vector element-wise multiply, or
/// scalar-vector multiply, depending on the dimension-ordering and
/// requested orientation.  These reduced cases are then
/// re-dispatched, allowing them to be handled by a vendor library,

template <typename DstBlock,
	  typename VBlock,
	  typename MBlock,
	  dimension_type SD>
struct Serial_expr_evaluator<2, DstBlock,
			     const Vmmul_expr_block<SD, VBlock, MBlock>,
			     Op_expr_tag>
{
  static char const* name() { return "Expr_Loop_Vmmul"; }

  typedef Vmmul_expr_block<SD, VBlock, MBlock> SrcBlock;

  typedef typename DstBlock::value_type dst_type;
  typedef typename VBlock::value_type v_type;
  typedef typename MBlock::value_type m_type;

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
    !Is_expr_block<MBlock>::value;

  static bool rt_valid(DstBlock& /*dst*/, SrcBlock const& /*src*/)
  {
    return true;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    VBlock const& vblock = src.get_vblk();
    MBlock const& mblock = src.get_mblk();

    typedef typename Block_layout<DstBlock>::order_type order_type;

    Matrix<dst_type, DstBlock>   m_dst(dst);
    const_Vector<v_type, VBlock> v(const_cast<VBlock&>(vblock));
    const_Matrix<m_type, MBlock> m(const_cast<MBlock&>(mblock));

    if (Type_equal<order_type, row2_type>::value)
    {
      if (SD == row)
      {
	for (index_type r=0; r<dst.size(2, 0); ++r)
	  m_dst.row(r) = v * m.row(r);
      }
      else
      {
	for (index_type r=0; r<dst.size(2, 0); ++r)
	  m_dst.row(r) = v.get(r) * m.row(r);
      }
    }
    else // col2_type
    {
      if (SD == row)
      {
	for (index_type c=0; c<dst.size(2, 1); ++c)
	  m_dst.col(c) = v.get(c) * m.col(c);
      }
      else
      {
	for (index_type c=0; c<dst.size(2, 1); ++c)
	  m_dst.col(c) = v * m.col(c);
      }
    }
  }
};
#endif

} // namespace vsip::impl



/// Vector-matrix element-wise multiplication

template <dimension_type Dim,
	  typename       T0,
	  typename       T1,
	  typename       Block0,
	  typename       Block1>
typename vsip::impl::Vmmul_traits<Dim, T0, T1, Block0, Block1>::view_type
vmmul(
  const_Vector<T0, Block0> v,
  const_Matrix<T1, Block1> m)
VSIP_NOTHROW
{
  typedef impl::Vmmul_traits<Dim, T0, T1, Block0, Block1> traits;
  typedef typename traits::block_type block_type;
  typedef typename traits::view_type  view_type;

  return view_type(block_type(v.block(), m.block()));
}

} // namespace vsip

#endif // VSIP_IMPL_VMMUL_HPP
