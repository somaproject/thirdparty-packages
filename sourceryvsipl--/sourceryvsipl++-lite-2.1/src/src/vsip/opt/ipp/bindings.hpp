/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ipp/bindings.hpp
    @author  Stefan Seefeld
    @date    2005-08-10
    @brief   VSIPL++ Library: Wrappers and traits to bridge with Intel's IPP.
*/

#ifndef VSIP_IMPL_IPP_HPP
#define VSIP_IMPL_IPP_HPP

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

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace ipp
{

template <typename Type>
struct Is_type_supported
{
  static bool const value = false;
};

template <>
struct Is_type_supported<float>
{
  static bool const value = true;
};

template <>
struct Is_type_supported<double>
{
  static bool const value = true;
};

template <>
struct Is_type_supported<std::complex<float> >
{
  static bool const value = true;
};

template <>
struct Is_type_supported<std::complex<double> >
{
  static bool const value = true;
};

#define VSIP_IMPL_IPP_DECL_V(FCN, T, IPPFCN, IPPT)			\
void									\
FCN(									\
  T const* A,								\
  T*       Z,								\
  length_type len);

#define VSIP_IMPL_IPP_DECL_V_CR(FCN, T, IPPFCN, IPPCT, IPPT)		\
void									\
FCN(									\
  complex<T> const* A,							\
  T*       Z,								\
  length_type len);

// Square
VSIP_IMPL_IPP_DECL_V(vsq,  float,           ippsSqr_32f,  Ipp32f)
VSIP_IMPL_IPP_DECL_V(vsq,  double,          ippsSqr_64f,  Ipp64f)
VSIP_IMPL_IPP_DECL_V(vsq,  complex<float>,  ippsSqr_32fc, Ipp32fc)
VSIP_IMPL_IPP_DECL_V(vsq,  complex<double>, ippsSqr_64fc, Ipp64fc)

// Square-root
VSIP_IMPL_IPP_DECL_V(vsqrt, float,           ippsSqrt_32f,  Ipp32f)
VSIP_IMPL_IPP_DECL_V(vsqrt, double,          ippsSqrt_64f,  Ipp64f)
VSIP_IMPL_IPP_DECL_V(vsqrt, complex<float>,  ippsSqrt_32fc, Ipp32fc)
VSIP_IMPL_IPP_DECL_V(vsqrt, complex<double>, ippsSqrt_64fc, Ipp64fc)

// Mag 
VSIP_IMPL_IPP_DECL_V(vmag, float,           ippsAbs_32f,  Ipp32f)
VSIP_IMPL_IPP_DECL_V(vmag, double,          ippsAbs_64f,  Ipp64f)
VSIP_IMPL_IPP_DECL_V_CR(vmag, float,  ippsMagnitude_32f, Ipp32fc, Ipp32f)
VSIP_IMPL_IPP_DECL_V_CR(vmag, double, ippsMagnitude_64f, Ipp64fc, Ipp64f)

// Mag-sq
VSIP_IMPL_IPP_DECL_V(vmagsq,    float,  ippsSqr_32f,  Ipp32f)
VSIP_IMPL_IPP_DECL_V(vmagsq,    double, ippsSqr_64f,  Ipp64f)

// functions for vector addition
void vadd(float const* A, float const* B, float* Z, length_type len);
void vadd(double const* A, double const* B, double* Z, length_type len);
void vadd(std::complex<float> const* A, std::complex<float> const* B,
          std::complex<float>* Z, length_type len);
void vadd(std::complex<double> const* A, std::complex<double> const* B,
          std::complex<double>* Z, length_type len);

// functions for vector copy
void vcopy(float const* A, float* Z, length_type len);
void vcopy(double const* A, double* Z, length_type len);
void vcopy(complex<float> const* A, complex<float>* Z, length_type len);
void vcopy(complex<double> const* A, complex<double>* Z, length_type len);

// functions for vector subtraction
void vsub(float const* A, float const* B, float* Z, length_type len);
void vsub(double const* A, double const* B, double* Z, length_type len);
void vsub(std::complex<float> const* A, std::complex<float> const* B,
          std::complex<float>* Z, length_type len);
void vsub(std::complex<double> const* A, std::complex<double> const* B,
          std::complex<double>* Z, length_type len);

// functions for vector multiply
void vmul(float const* A, float const* B, float* Z, length_type len);
void vmul(double const* A, double const* B, double* Z, length_type len);
void vmul(std::complex<float> const* A, std::complex<float> const* B,
          std::complex<float>* Z, length_type len);
void vmul(std::complex<double> const* A, std::complex<double> const* B,
          std::complex<double>* Z, length_type len);

void svadd(float A, float const* B, float* Z, length_type len);
void svadd(double A, double const* B, double* Z, length_type len);
void svadd(complex<float> A, complex<float> const* B,
	   complex<float>* Z, length_type len);
void svadd(complex<double> A, complex<double> const* B,
	   complex<double>* Z, length_type len);

// sub: scalar - vector
void svsub(float A, float const* B, float* Z, length_type len);
void svsub(double A, double const* B, double* Z, length_type len);
void svsub(complex<float> A, complex<float> const* B,
	   complex<float>* Z, length_type len);
void svsub(complex<double> A, complex<double> const* B,
	   complex<double>* Z, length_type len);

// sub: vector - scalar
void svsub(float const* A, float B, float* Z, length_type len);
void svsub(double const* A, double B, double* Z, length_type len);
void svsub(complex<float> const* A, complex<float> B,
	   complex<float>* Z, length_type len);
void svsub(complex<double> const* A, complex<double> B,
	   complex<double>* Z, length_type len);

void svmul(float A, float const* B, float* Z, length_type len);
void svmul(double A, double const* B, double* Z, length_type len);
void svmul(complex<float> A, complex<float> const* B,
	   complex<float>* Z, length_type len);
void svmul(complex<double> A, complex<double> const* B,
	   complex<double>* Z, length_type len);

// functions for scalar-view division: scalar / vector
void svdiv(float A, float const* B, float* Z, length_type len);

// functions for scalar-view division: vector / scalar
void svdiv(float const* A, float B, float* Z, length_type len);
void svdiv(double const* A, double B, double* Z, length_type len);
void svdiv(complex<float> const* A, complex<float> B,
	   complex<float>* Z, length_type len);
void svdiv(complex<double> const* A, complex<double> B,
	   complex<double>* Z, length_type len);

// functions for vector division
void vdiv(float const* A, float const* B, float* Z, length_type len);
void vdiv(double const* A, double const* B, double* Z, length_type len);
void vdiv(std::complex<float> const* A, std::complex<float> const* B,
          std::complex<float>* Z, length_type len);
void vdiv(std::complex<double> const* A, std::complex<double> const* B,
          std::complex<double>* Z, length_type len);

// Vector zero
void vzero(float*           Z, length_type len);
void vzero(double*          Z, length_type len);
void vzero(complex<float>*  Z, length_type len);
void vzero(complex<double>* Z, length_type len);

// functions for convolution
void conv(float* coeff, length_type coeff_size,
	  float* in,    length_type in_size,
	  float* out);
void conv(double* coeff, length_type coeff_size,
	  double* in,    length_type in_size,
	  double* out);

void
conv_full_2d(
  float*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  float*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  float*      out,
  length_type out_row_stride);

void
conv_full_2d(
  short*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  short*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  short*      out,
  length_type out_row_stride);

void
conv_valid_2d(
  float*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  float*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  float*      out,
  length_type out_row_stride);

void
conv_valid_2d(
  short*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  short*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  short*      out,
  length_type out_row_stride);



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
     ipp::Is_type_supported<LType>::value &&
     ipp::Is_type_supported<RType>::value &&
     ipp::Is_type_supported<typename DstBlock::value_type>::value &&
     Type_equal<typename DstBlock::value_type, LType>::value &&
     Type_equal<typename DstBlock::value_type, RType>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<LBlock>::value == 0 &&
     Ext_data_cost<RBlock>::value == 0 &&
     /* IPP does not support complex split */
     !Is_split_block<DstBlock>::value &&
     !Is_split_block<LBlock>::value &&
     !Is_split_block<RBlock>::value;

  
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>    ext_dst(dst,       SYNC_OUT);
    Ext_data<LBlock,   lblock_lp> ext_l(src.left(),  SYNC_IN);
    Ext_data<RBlock,   rblock_lp> ext_r(src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_l.stride(0) == 1 &&
	    ext_r.stride(0) == 1);
  }
};
} // namespace vsip::impl::ipp



/***********************************************************************
  Unary expression evaluators
***********************************************************************/

#define VSIP_IMPL_IPP_V_EXPR(OP, FUN)					\
template <typename DstBlock,						\
	  typename Block,						\
	  typename Type>						\
struct Serial_expr_evaluator<						\
    1, DstBlock, 							\
    Unary_expr_block<1, OP, Block, Type> const,				\
    Intel_ipp_tag>							\
{									\
  static char const* name() { return "Expr_IPP_V-" #FUN; }		\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block>::layout_type>::type		\
    blk_lp;								\
									\
  typedef Unary_expr_block<1, OP, Block, Type>				\
    SrcBlock;								\
									\
  static bool const ct_valid =						\
    !Is_expr_block<Block>::value &&					\
     ipp::Is_type_supported<typename DstBlock::value_type>::value &&	\
     Type_equal<typename DstBlock::value_type, Type>::value &&		\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<Block>::value == 0 &&				\
     /* IPP does not support complex split */				\
     !Is_split_block<DstBlock>::value &&				\
     !Is_split_block<Block>::value;					\
  									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp> ext_dst(dst,      SYNC_OUT);		\
    Ext_data<Block,    blk_lp> ext_src(src.op(), SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_src.stride(0) == 1);					\
  }									\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp> ext_dst(dst,      SYNC_OUT);		\
    Ext_data<Block,    blk_lp> ext_src(src.op(), SYNC_IN);		\
									\
    FUN(								\
      ext_src.data(),							\
      ext_dst.data(),							\
      dst.size()							\
    );									\
  }									\
};



#define VSIP_IMPL_IPP_V_CR_EXPR(OP, FUN)				\
template <typename DstBlock,						\
	  typename Block,						\
	  typename Type>						\
struct Serial_expr_evaluator<						\
    1, DstBlock, 							\
    Unary_expr_block<1, OP, Block, Type> const,				\
    Intel_ipp_tag>							\
{									\
  static char const* name() { return "Expr_IPP_V_CR-" #FUN; }		\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<Block>::layout_type>::type		\
    blk_lp;								\
									\
  typedef Unary_expr_block<1, OP, Block, Type>				\
    SrcBlock;								\
									\
  static bool const ct_valid =						\
    !Is_expr_block<Block>::value &&					\
     ipp::Is_type_supported<typename DstBlock::value_type>::value &&	\
     ipp::Is_type_supported<Type>::value &&				\
     /* Type_equal<typename DstBlock::value_type, Type>::value && */	\
     /* check that direct access is supported */			\
     Ext_data_cost<DstBlock>::value == 0 &&				\
     Ext_data_cost<Block>::value == 0 &&				\
     /* IPP does not support complex split */				\
     !Is_split_block<DstBlock>::value &&				\
     !Is_split_block<Block>::value;					\
  									\
  static bool rt_valid(DstBlock& dst, SrcBlock const& src)		\
  {									\
    /* check if all data is unit stride */				\
    Ext_data<DstBlock, dst_lp> ext_dst(dst,      SYNC_OUT);		\
    Ext_data<Block,    blk_lp> ext_src(src.op(), SYNC_IN);		\
    return (ext_dst.stride(0) == 1 &&					\
	    ext_src.stride(0) == 1);					\
  }									\
  									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp> ext_dst(dst,      SYNC_OUT);		\
    Ext_data<Block,    blk_lp> ext_src(src.op(), SYNC_IN);		\
									\
    FUN(								\
      ext_src.data(),							\
      ext_dst.data(),							\
      dst.size()							\
    );									\
  }									\
};

VSIP_IMPL_IPP_V_EXPR(sq_functor,   ipp::vsq)
VSIP_IMPL_IPP_V_EXPR(sqrt_functor, ipp::vsqrt)
VSIP_IMPL_IPP_V_CR_EXPR(mag_functor,  ipp::vmag)
// Don't dispatch for now since only real magsq is supported.
// VSIP_IMPL_IPP_V_CR_EXPR(magsq_functor,  ipp::vmagsq)



/***********************************************************************
  Binary expression evaluators
***********************************************************************/

#define VSIP_IMPL_IPP_VV_EXPR(OP, FUN)					\
template <typename DstBlock,						\
	  typename LBlock,						\
	  typename RBlock,						\
	  typename LType,						\
	  typename RType>						\
struct Serial_expr_evaluator<						\
    1, DstBlock, 							\
    Binary_expr_block<1, OP, LBlock, LType, RBlock, RType> const,	\
    Intel_ipp_tag>							\
  : ipp::Serial_expr_evaluator_base<OP, DstBlock,			\
				    LBlock, RBlock, LType, RType>	\
{									\
  static char const* name() { return "Expr_IPP_VV-" #FUN; }		\
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


VSIP_IMPL_IPP_VV_EXPR(op::Add,  ipp::vadd)
VSIP_IMPL_IPP_VV_EXPR(op::Sub,  ipp::vsub)
VSIP_IMPL_IPP_VV_EXPR(op::Mult, ipp::vmul)
VSIP_IMPL_IPP_VV_EXPR(op::Div,  ipp::vdiv)



/***********************************************************************
  Scalar-view element-wise operations
***********************************************************************/

namespace ipp
{

template <template <typename, typename> class Operator,
	  typename DstBlock,
	  typename SType,
	  typename VBlock,
	  typename VType,
	  bool     Right>
struct Scalar_view_evaluator_base
{
  typedef Binary_expr_block<1, Operator,
			    Scalar_block<1, SType>, SType,
			    VBlock, VType>
	SrcBlock;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<VBlock>::layout_type>::type
    vblock_lp;

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
     ipp::Is_type_supported<typename DstBlock::value_type>::value &&
     Type_equal<typename DstBlock::value_type, SType>::value &&
     Type_equal<typename DstBlock::value_type, VType>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<VBlock>::value == 0 &&
     // Complex split format is not supported.
     Type_equal<typename Block_layout<DstBlock>::complex_type,
		Cmplx_inter_fmt>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_r(src.right(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_r.stride(0) == 1);
  }

};



template <template <typename, typename> class Operator,
	  typename DstBlock,
	  typename SType,
	  typename VBlock,
	  typename VType>
struct Scalar_view_evaluator_base<Operator, DstBlock, SType, VBlock, VType,
				  false>
{
  typedef Binary_expr_block<1, Operator,
			    VBlock, VType,
			    Scalar_block<1, SType>, SType>
	SrcBlock;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<DstBlock>::layout_type>::type
    dst_lp;

  typedef typename Adjust_layout_dim<
      1, typename Block_layout<VBlock>::layout_type>::type
    vblock_lp;

  static bool const ct_valid = 
    !Is_expr_block<VBlock>::value &&
     ipp::Is_type_supported<typename DstBlock::value_type>::value &&
     Type_equal<typename DstBlock::value_type, SType>::value &&
     Type_equal<typename DstBlock::value_type, VType>::value &&
     // check that direct access is supported
     Ext_data_cost<DstBlock>::value == 0 &&
     Ext_data_cost<VBlock>::value == 0 &&
     // Complex split format is not supported.
     Type_equal<typename Block_layout<DstBlock>::complex_type,
		Cmplx_inter_fmt>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    // check if all data is unit stride
    Ext_data<DstBlock, dst_lp>  ext_dst(dst, SYNC_OUT);
    Ext_data<VBlock, vblock_lp> ext_l(src.left(), SYNC_IN);
    return (ext_dst.stride(0) == 1 &&
	    ext_l.stride(0) == 1);
  }

};

} // namespace vsip::impl::ipp


#define VSIP_IMPL_IPP_SV_EXPR(OP, FCN)					\
template <typename DstBlock,						\
	  typename SType,						\
	  typename VBlock,						\
	  typename VType>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         const Binary_expr_block<1, OP,					\
                                 Scalar_block<1, SType>, SType,		\
                                 VBlock, VType>,			\
         Intel_ipp_tag>							\
  : ipp::Scalar_view_evaluator_base<OP, DstBlock, SType,		\
                                    VBlock, VType, true>		\
{									\
  static char const* name() { return "Expr_IPP_SV-" #FCN; }		\
									\
  typedef Binary_expr_block<1, OP,					\
			    Scalar_block<1, SType>, SType,		\
			    VBlock, VType>				\
	SrcBlock;							\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<VBlock>::layout_type>::type		\
    vblock_lp;								\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);		\
    Ext_data<VBlock, vblock_lp> ext_r(src.right(), SYNC_IN);		\
    FCN(src.left().value(), ext_r.data(), ext_dst.data(), dst.size());	\
  }									\
};

#define VSIP_IMPL_IPP_SV_EXPR_FO(OP, FCN)				\
template <typename DstBlock,						\
	  typename VBlock>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         const Binary_expr_block<1, OP,					\
                                 Scalar_block<1, float>, float,		\
                                 VBlock, float>,			\
         Intel_ipp_tag>							\
  : ipp::Scalar_view_evaluator_base<OP, DstBlock, float,		\
                                    VBlock, float, true>		\
{									\
  static char const* name() { return "Expr_IPP_SV_FO-" #FCN; }		\
									\
  typedef float SType;							\
  typedef float VType;							\
									\
  typedef Binary_expr_block<1, OP,					\
			    Scalar_block<1, SType>, SType,		\
			    VBlock, VType>				\
	SrcBlock;							\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<VBlock>::layout_type>::type		\
    vblock_lp;								\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);		\
    Ext_data<VBlock, vblock_lp> ext_r(src.right(), SYNC_IN);		\
    FCN(src.left().value(), ext_r.data(), ext_dst.data(), dst.size());	\
  }									\
};



#define VSIP_IMPL_IPP_VS_EXPR(OP, FCN)					\
template <typename DstBlock,						\
	  typename SType,						\
	  typename VBlock,						\
	  typename VType>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         const Binary_expr_block<1, OP,					\
                                 VBlock, VType,				\
                                 Scalar_block<1, SType>, SType>,	\
         Intel_ipp_tag>							\
  : ipp::Scalar_view_evaluator_base<OP, DstBlock, SType,		\
                                    VBlock, VType, false>		\
{									\
  static char const* name() { return "Expr_IPP_VS-" #FCN; }		\
									\
  typedef Binary_expr_block<1, OP,					\
			    VBlock, VType,				\
			    Scalar_block<1, SType>, SType>		\
	SrcBlock;							\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<VBlock>::layout_type>::type		\
    vblock_lp;								\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);		\
    Ext_data<VBlock, vblock_lp> ext_l(src.left(), SYNC_IN);		\
    FCN(ext_l.data(), src.right().value(), ext_dst.data(), dst.size());	\
  }									\
};

#define VSIP_IMPL_IPP_VS_AS_SV_EXPR(OP, FCN)				\
template <typename DstBlock,						\
	  typename SType,						\
	  typename VBlock,						\
	  typename VType>						\
struct Serial_expr_evaluator<						\
         1, DstBlock, 							\
         const Binary_expr_block<1, OP,					\
                                 VBlock, VType,				\
                                 Scalar_block<1, SType>, SType>,	\
         Intel_ipp_tag>							\
  : ipp::Scalar_view_evaluator_base<OP, DstBlock, SType,		\
                                    VBlock, VType, false>		\
{									\
  static char const* name() { return "Expr_IPP_VS_AS_SV-" #FCN; }	\
									\
  typedef Binary_expr_block<1, OP,					\
			    VBlock, VType,				\
			    Scalar_block<1, SType>, SType>		\
	SrcBlock;							\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<DstBlock>::layout_type>::type		\
    dst_lp;								\
									\
  typedef typename Adjust_layout_dim<					\
      1, typename Block_layout<VBlock>::layout_type>::type		\
    vblock_lp;								\
									\
  static void exec(DstBlock& dst, SrcBlock const& src)			\
  {									\
    Ext_data<DstBlock, dst_lp>  ext_dst(dst,       SYNC_OUT);		\
    Ext_data<VBlock, vblock_lp> ext_l(src.left(), SYNC_IN);		\
    FCN(src.right().value(), ext_l.data(), ext_dst.data(), dst.size());	\
  }									\
};

VSIP_IMPL_IPP_SV_EXPR      (op::Add,  ipp::svadd)
VSIP_IMPL_IPP_VS_AS_SV_EXPR(op::Add,  ipp::svadd)
VSIP_IMPL_IPP_SV_EXPR      (op::Sub,  ipp::svsub)
VSIP_IMPL_IPP_VS_EXPR      (op::Sub,  ipp::svsub)
VSIP_IMPL_IPP_SV_EXPR      (op::Mult, ipp::svmul)
VSIP_IMPL_IPP_VS_AS_SV_EXPR(op::Mult, ipp::svmul)
VSIP_IMPL_IPP_SV_EXPR_FO   (op::Div,  ipp::svdiv)
VSIP_IMPL_IPP_VS_EXPR      (op::Div,  ipp::svdiv)

#undef VSIP_IMPL_IPP_SV_EXPR
#undef VSIP_IMPL_IPP_SV_EXPR_FO
#undef VSIP_IMPL_IPP_VS_EXPR
#undef VSIP_IMPL_IPP_VS_AS_SV_EXPR


} // namespace vsip::impl
} // namespace vsip

#endif
