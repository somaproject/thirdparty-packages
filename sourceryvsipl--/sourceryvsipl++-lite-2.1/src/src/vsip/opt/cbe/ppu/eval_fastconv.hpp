/* Copyright (c) 2007 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/cbe/ppu/eval_fastconv.hpp
    @author  Jules Bergmann
    @date    2007-03-05
    @brief   VSIPL++ Library: General evaluator for fast convolution

*/

#ifndef VSIP_OPT_CBE_PPU_EVAL_FASTCONV_HPP
#define VSIP_OPT_CBE_PPU_EVAL_FASTCONV_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/fft.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/opt/cbe/ppu/fastconv.hpp>
#include <vsip/opt/expr/return_block.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Evaluator for return expression block.

template <typename       DstBlock,
	  typename       T,
	  typename       VecBlockT,
	  typename       MatBlockT,
	  typename       Backend1T,
	  typename       Workspace1T,
	  typename       Backend2T,
	  typename       Workspace2T>
struct Serial_expr_evaluator<2, DstBlock,
  const Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Vmmul_expr_block<0,
        VecBlockT,
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
          >
        >,
      Backend1T, Workspace1T>
    >,
  Cbe_sdk_tag
  >
{
  static char const* name() { return "Cbe_sdk_tag"; }

  typedef
  Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Vmmul_expr_block<0,
        VecBlockT,
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
          >
        >,
      Backend1T, Workspace1T>
    >
    SrcBlock;

  typedef typename DstBlock::value_type dst_type;
  typedef typename SrcBlock::value_type src_type;
  typedef typename Block_layout<DstBlock>::complex_type complex_type;
  typedef impl::cbe::Fastconv<1, T, complex_type> fconv_type;

  static bool const ct_valid = 
    Type_equal<T, std::complex<float> >::value &&
    Type_equal<T, typename VecBlockT::value_type>::value &&
    Type_equal<T, typename MatBlockT::value_type>::value;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<VecBlockT> ext_kernel(src.functor().block().get_vblk());
    Ext_data<MatBlockT> ext_in    (src.functor().block().get_mblk().functor().block());
    Ext_data<DstBlock>  ext_out   (dst);

    return 
      fconv_type::rt_valid_size(dst.size(2, 1)) &&
      ext_kernel.stride(0) == 1 &&
      ext_in.stride(1) == 1 &&
      ext_out.stride(1) == 1;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type cols = dst.size(2, 1);
    Matrix<T> tmp(1, cols);

    Vector<T, VecBlockT> w 
      (const_cast<VecBlockT&>(src.functor().block().get_vblk()));
    Matrix<T, MatBlockT> in 
      (const_cast<MatBlockT&>(src.functor().block().get_mblk().functor().block()));
    Matrix<T, DstBlock> out(dst);

    fconv_type fconv(w, cols, false);

    fconv(in, out);
  }
};


template <typename       DstBlock,
	  typename       T,
	  typename       CoeffsMatBlockT,
	  typename       MatBlockT,
	  typename       Backend1T,
	  typename       Workspace1T,
	  typename       Backend2T,
	  typename       Workspace2T>
struct Serial_expr_evaluator<2, DstBlock,
  const Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Binary_expr_block<2, op::Mult,
        CoeffsMatBlockT, T, 
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
        >, T
      >,
      Backend1T, Workspace1T>
    >,
    Cbe_sdk_tag
  >
{
  static char const* name() { return "Cbe_sdk_tag"; }

  typedef
  Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Binary_expr_block<2, op::Mult,
        CoeffsMatBlockT, T,
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
        >, T
      >,
      Backend1T, Workspace1T>
    >
    SrcBlock;

  typedef typename DstBlock::value_type dst_type;
  typedef typename SrcBlock::value_type src_type;
  typedef typename Block_layout<DstBlock>::complex_type complex_type;
  typedef impl::cbe::Fastconv_base<2, T, complex_type> fconv_type;

  static bool const ct_valid = 
    Type_equal<T, std::complex<float> >::value &&
    Type_equal<T, typename CoeffsMatBlockT::value_type>::value &&
    Type_equal<T, typename MatBlockT::value_type>::value &&
    Ext_data_cost<CoeffsMatBlockT>::value == 0 &&
    Ext_data_cost<MatBlockT>::value == 0 &&
    Ext_data_cost<DstBlock>::value == 0;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<CoeffsMatBlockT> ext_kernel(src.functor().block().left());
    Ext_data<MatBlockT>       ext_in    (src.functor().block().right().functor().block());
    Ext_data<DstBlock>        ext_out   (dst);

    return 
      fconv_type::rt_valid_size(dst.size(2, 1)) &&
      ext_kernel.stride(1) == 1 &&
      ext_in.stride(1) == 1 &&
      ext_out.stride(1) == 1;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type cols = dst.size(2, 1);
    Matrix<T, CoeffsMatBlockT> w 
      (const_cast<CoeffsMatBlockT&>(src.functor().block().left()));
    Matrix<T, MatBlockT> in 
      (const_cast<MatBlockT&>(src.functor().block().right().functor().block()));
    Matrix<T, DstBlock> out(dst);

    fconv_type fconv(cols, false);

    fconv.convolve(in, w, out);
  }
};


template <typename       DstBlock,
	  typename       T,
	  typename       CoeffsMatBlockT,
	  typename       MatBlockT,
	  typename       Backend1T,
	  typename       Workspace1T,
	  typename       Backend2T,
	  typename       Workspace2T>
struct Serial_expr_evaluator<2, DstBlock,
  const Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Binary_expr_block<2, op::Mult,
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
        >, T, 
        CoeffsMatBlockT, T
      >,
      Backend1T, Workspace1T>
    >,
    Cbe_sdk_tag
  >
{
  static char const* name() { return "Cbe_sdk_tag"; }

  typedef
  Return_expr_block<2, T,
    fft::Fft_return_functor<2, T,
      const Binary_expr_block<2, op::Mult,
        const Return_expr_block<2, T,
          fft::Fft_return_functor<2, T,
            MatBlockT,
            Backend2T, Workspace2T>
        >, T,
        CoeffsMatBlockT, T
      >,
      Backend1T, Workspace1T>
    >
    SrcBlock;

  typedef typename DstBlock::value_type dst_type;
  typedef typename SrcBlock::value_type src_type;
  typedef typename Block_layout<DstBlock>::complex_type complex_type;
  typedef impl::cbe::Fastconv_base<2, T, complex_type> fconv_type;

  static bool const ct_valid = 
    Type_equal<T, std::complex<float> >::value &&
    Type_equal<T, typename CoeffsMatBlockT::value_type>::value &&
    Type_equal<T, typename MatBlockT::value_type>::value &&
    Ext_data_cost<CoeffsMatBlockT>::value == 0 &&
    Ext_data_cost<MatBlockT>::value == 0 &&
    Ext_data_cost<DstBlock>::value == 0;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    Ext_data<CoeffsMatBlockT> ext_kernel(src.functor().block().right());
    Ext_data<MatBlockT>       ext_in    (src.functor().block().left().functor().block());
    Ext_data<DstBlock>        ext_out   (dst);

    return 
      fconv_type::rt_valid_size(dst.size(2, 1)) &&
      ext_kernel.stride(1) == 1 &&
      ext_in.stride(1) == 1 &&
      ext_out.stride(1) == 1;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type cols = dst.size(2, 1);
    Matrix<T, CoeffsMatBlockT> w 
      (const_cast<CoeffsMatBlockT&>(src.functor().block().right()));
    Matrix<T, MatBlockT> in 
      (const_cast<MatBlockT&>(src.functor().block().left().functor().block()));
    Matrix<T, DstBlock> out(dst);

    fconv_type fconv(cols, false);

    fconv.convolve(in, w, out);
  }
};

} // namespace vsip::impl

} // namespace vsip

#endif // VSIP_OPT_CBE_PPU_EVAL_FASTCONV_HPP
