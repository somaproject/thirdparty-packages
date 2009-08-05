/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/cml/transpose.hpp
    @author  Don McCoy
    @date    2008-06-04
    @brief   VSIPL++ Library: Bindings for CML matrix transpose.
*/

#ifndef VSIP_OPT_CBE_CML_TRANSPOSE_HPP
#define VSIP_OPT_CBE_CML_TRANSPOSE_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/impl_tags.hpp>
#include <vsip/opt/cbe/cml/traits.hpp>
#include <vsip/opt/cbe/ppu/task_manager.hpp>

#include <cml.h>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cml
{

// These macros support scalar and interleaved complex types

#define VSIP_IMPL_CML_TRANS(T, FCN, CML_FCN)    \
  inline void                                   \
  FCN(                                          \
    T* a, ptrdiff_t rsa, ptrdiff_t csa,         \
    T* z, ptrdiff_t rsz, ptrdiff_t csz,         \
    size_t m, size_t n)                         \
  {                                             \
    typedef Scalar_of<T>::type CML_T;           \
    CML_FCN(                                    \
      reinterpret_cast<CML_T*>(a), rsa, csa,    \
      reinterpret_cast<CML_T*>(z), rsz, csz,    \
      m, n );                                   \
  }

VSIP_IMPL_CML_TRANS(float,               transpose, cml_mtrans_f)
VSIP_IMPL_CML_TRANS(std::complex<float>, transpose, cml_cmtrans_f)
#undef VSIP_IMPL_CML_TRANS


#define VSIP_IMPL_CML_TRANS_UNIT(T, FCN, CML_FCN)  \
  inline void                                      \
  FCN(                                             \
    T* a, ptrdiff_t rsa,                           \
    T* z, ptrdiff_t rsz,                           \
    size_t m, size_t n)                            \
  {                                                \
    typedef Scalar_of<T>::type CML_T;              \
    CML_FCN(                                       \
      reinterpret_cast<CML_T*>(a), rsa,            \
      reinterpret_cast<CML_T*>(z), rsz,            \
      m, n );                                      \
  }

VSIP_IMPL_CML_TRANS_UNIT(float,               transpose_unit, cml_mtrans1_f)
VSIP_IMPL_CML_TRANS_UNIT(std::complex<float>, transpose_unit, cml_cmtrans1_f)
#undef VSIP_IMPL_CML_TRANS_UNIT


#define VSIP_IMPL_CML_VCOPY(T, FCN, CML_FCN)   \
  inline void                                      \
  FCN(                                             \
    T* a, ptrdiff_t rsa,                           \
    T* z, ptrdiff_t rsz,                           \
    size_t n)                                      \
  {                                                \
    typedef Scalar_of<T>::type CML_T;              \
    CML_FCN(                                       \
      reinterpret_cast<CML_T*>(a), rsa,            \
      reinterpret_cast<CML_T*>(z), rsz,            \
      n * (Is_complex<T>::value ? 2 : 1));         \
  }

VSIP_IMPL_CML_VCOPY(float,          vcopy, cml_vcopy_f)
VSIP_IMPL_CML_VCOPY(complex<float>, vcopy, cml_vcopy_f)
#undef VSIP_IMPL_CML_VCOPY


// These macros support split complex types only

#define VSIP_IMPL_CML_TRANS_SPLIT(T, FCN, CML_FCN)     \
  inline void                                          \
  FCN(                                                 \
    std::pair<T*, T*> a, ptrdiff_t rsa, ptrdiff_t csa, \
    std::pair<T*, T*> z, ptrdiff_t rsz, ptrdiff_t csz, \
    size_t m, size_t n)                                \
  {                                                    \
    CML_FCN(                                           \
      a.first, a.second, rsa, csa,                     \
      z.first, z.second, rsz, csz,                     \
      m, n );                                          \
  }

VSIP_IMPL_CML_TRANS_SPLIT(float, transpose, cml_zmtrans_f)
#undef VSIP_IMPL_CML_TRANS_SPLIT


#define VSIP_IMPL_CML_TRANS_UNIT_SPLIT(T, FCN, CML_FCN) \
  inline void                                           \
  FCN(                                                  \
    std::pair<T*, T*> a, ptrdiff_t rsa,                 \
    std::pair<T*, T*> z, ptrdiff_t rsz,                 \
    size_t m, size_t n)                                 \
  {                                                     \
    CML_FCN(                                            \
      a.first, a.second, rsa,                           \
      z.first, z.second, rsz,                           \
      m, n );                                           \
  }

VSIP_IMPL_CML_TRANS_UNIT_SPLIT(float, transpose_unit, cml_zmtrans1_f)
#undef VSIP_IMPL_CML_TRANS_UNIT_SPLIT


#define VSIP_IMPL_CML_VCOPY_SPLIT(T, FCN, CML_FCN)  \
  inline void                                           \
  FCN(                                                  \
    std::pair<T*, T*> a, ptrdiff_t rsa,                 \
    std::pair<T*, T*> z, ptrdiff_t rsz,                 \
    size_t n)                                           \
  {                                                     \
    CML_FCN(                                            \
      a.first, rsa,                                     \
      z.first, rsz,                                     \
      n);                                               \
    CML_FCN(                                            \
      a.second, rsa,                                    \
      z.second, rsz,                                    \
      n);                                               \
  }

VSIP_IMPL_CML_VCOPY_SPLIT(float, vcopy, cml_vcopy_f)
#undef VSIP_IMPL_CML_VCOPY_SPLIT


} // namespace vsip::impl::cml



template <typename DstBlock,
          typename SrcBlock>
struct Serial_expr_evaluator<2, DstBlock, SrcBlock, Cml_tag>
{
  typedef typename Block_layout<SrcBlock>::order_type src_order_type;
  typedef typename Block_layout<DstBlock>::order_type dst_order_type;

  typedef typename DstBlock::value_type dst_value_type;
  typedef typename SrcBlock::value_type src_value_type;

  static char const* name()
  {
    char s = Type_equal<src_order_type, row2_type>::value ? 'r' : 'c';
    char d = Type_equal<dst_order_type, row2_type>::value ? 'r' : 'c';
    if      (s == 'r' && d == 'r')    return "Expr_CML_Trans (rr copy)";
    else if (s == 'r' && d == 'c')    return "Expr_CML_Trans (rc trans)";
    else if (s == 'c' && d == 'r')    return "Expr_CML_Trans (cr trans)";
    else /* (s == 'c' && d == 'c') */ return "Expr_CML_Trans (cc copy)";
  }

  static bool const is_rhs_expr   = Is_expr_block<SrcBlock>::value;

  static bool const is_lhs_split  = Is_split_block<DstBlock>::value;
  static bool const is_rhs_split  = Is_split_block<SrcBlock>::value;

  static int const  lhs_cost      = Ext_data_cost<DstBlock>::value;
  static int const  rhs_cost      = Ext_data_cost<SrcBlock>::value;

  static bool const ct_valid =
    // check that CML supports this data type and/or layout
    cml::Cml_supports_block<SrcBlock>::valid &&
    cml::Cml_supports_block<DstBlock>::valid &&
    // check that types are equal
    Type_equal<src_value_type, dst_value_type>::value &&
    // check that the source block is not an expression
    !is_rhs_expr &&
    // check that direct access is supported
    lhs_cost == 0 && rhs_cost == 0 &&
    // check complex layout is consistent
    is_lhs_split == is_rhs_split;

  static length_type tunable_threshold()
  {
    if (VSIP_IMPL_TUNE_MODE)
      return 0;
    // Copy is always faster with SPU
    else if (Type_equal<src_order_type, dst_order_type>::value)
      return 0;
    // Transpose not always faster with SPU
    // mcopy -6
    else if (Type_equal<dst_value_type, complex<float> >::value)
      return 128*128;
    // mcopy -2
    else if (Type_equal<dst_value_type, float>::value)
      return 256*256;

    return 0;
  }

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  { 
    bool rt = true;

    // If performing a copy, both source and destination blocks
    // must be unit stride and dense.
    if (Type_equal<src_order_type, dst_order_type>::value)
    {
      Ext_data<DstBlock> dst_ext(dst, SYNC_OUT);
      Ext_data<SrcBlock> src_ext(src, SYNC_IN);

      dimension_type const s_dim0 = src_order_type::impl_dim0;
      dimension_type const s_dim1 = src_order_type::impl_dim1;
      dimension_type const d_dim0 = dst_order_type::impl_dim0;
      dimension_type const d_dim1 = dst_order_type::impl_dim1;

      if (dst_ext.stride(d_dim1) != 1 ||
	  dst_ext.stride(d_dim0) != static_cast<stride_type>(dst.size(2, d_dim1)) ||
	  src_ext.stride(s_dim1) != 1 ||
	  src_ext.stride(s_dim0) != static_cast<stride_type>(src.size(2, s_dim1)))
        rt = false;
    }

    rt &= dst.size(2, 0) * dst.size(2, 1) > tunable_threshold();
    rt &= cbe::Task_manager::instance()->num_spes() > 0;

    return rt; 
  }

  static void exec(DstBlock& dst, SrcBlock const& src, row2_type, row2_type)
  {
    Ext_data<DstBlock> dst_ext(dst, SYNC_OUT);
    Ext_data<SrcBlock> src_ext(src, SYNC_IN);

    if (dst_ext.stride(1) == 1 && src_ext.stride(1) == 1)
    {
      assert(dst_ext.stride(0) == static_cast<stride_type>(dst.size(2, 1)));
      assert(src_ext.stride(0) == static_cast<stride_type>(src.size(2, 1)));

      cml::vcopy(
        src_ext.data(), 1,
        dst_ext.data(), 1,
        dst.size(2, 0) * dst.size(2, 1) );
    }
    else
      assert(0);
  }

  static void exec(DstBlock& dst, SrcBlock const& src, col2_type, col2_type)
  {
    Ext_data<DstBlock> dst_ext(dst, SYNC_OUT);
    Ext_data<SrcBlock> src_ext(src, SYNC_IN);

    if (dst_ext.stride(0) == 1 && src_ext.stride(0) == 1)
    {
      assert(dst_ext.stride(1) == static_cast<stride_type>(dst.size(2, 0)));
      assert(src_ext.stride(1) == static_cast<stride_type>(src.size(2, 0)));

      cml::vcopy(
        src_ext.data(), 1,
        dst_ext.data(), 1,
        dst.size(2, 0) * dst.size(2, 1) );
    }
    else
      assert(0);
  }

  static void exec(DstBlock& dst, SrcBlock const& src, col2_type, row2_type)
  {
    Ext_data<DstBlock> dst_ext(dst, SYNC_OUT);
    Ext_data<SrcBlock> src_ext(src, SYNC_IN);

    if (dst_ext.stride(0) == 1 && src_ext.stride(1) == 1)
    {
      cml::transpose_unit(
        src_ext.data(), src_ext.stride(0),
        dst_ext.data(), dst_ext.stride(1),
        dst.size(2, 1), dst.size(2, 0));
    }
    else
    {
      cml::transpose(
        src_ext.data(), src_ext.stride(0), src_ext.stride(1),
        dst_ext.data(), dst_ext.stride(1), dst_ext.stride(0),
        dst.size(2, 1), dst.size(2, 0));
    }
  }

  static void exec(DstBlock& dst, SrcBlock const& src, row2_type, col2_type)
  {
    Ext_data<DstBlock> dst_ext(dst, SYNC_OUT);
    Ext_data<SrcBlock> src_ext(src, SYNC_IN);

    if (dst_ext.stride(1) == 1 && src_ext.stride(0) == 1)
    {
      cml::transpose_unit(
        src_ext.data(), src_ext.stride(1),
        dst_ext.data(), dst_ext.stride(0),
        dst.size(2, 0), dst.size(2, 1));
    }
    else
    {
      cml::transpose(
        src_ext.data(), src_ext.stride(1), src_ext.stride(0),
        dst_ext.data(), dst_ext.stride(0), dst_ext.stride(1),
        dst.size(2, 0), dst.size(2, 1));
    }
  }

  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    exec(dst, src, dst_order_type(), src_order_type());
  }
  
};


} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CBE_CML_TRANSPOSE_HPP
