/* Copyright (c) 2007 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/expr/eval_fastconv.hpp
    @author  Jules Bergmann
    @date    2007-02-02
    @brief   VSIPL++ Library: General evaluator for fast convolution

*/

#ifndef VSIP_OPT_EXPR_EVAL_FASTCONV_HPP
#define VSIP_OPT_EXPR_EVAL_FASTCONV_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/fft.hpp>
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
  Fc_expr_tag
  >
{
  static char const* name() { return "Fc_expr_tag-vw"; }

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

  static bool const ct_valid = true;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    (void)dst;
#if 0
    // Check if the evaluator supports the scaling implied by
    // the FFTs.
    //
    // This evaluator uses the FFTs directly, so it will implicitly
    // follow the requested scaling.  However, if this evaluator is
    // adapted to use other fast convolution implementations that
    // have limited scaling (such as only unit), this check will
    // be necessary.
    
    typedef typename Scalar_of<T>::type scalar_type;

    Workspace2T const& fwd_workspace(
      src.functor().block().get_mblk().functor().workspace());
    Workspace1T const& inv_workspace(src.functor().workspace());

    // Check FFT scaling totals 1
    return almost_equal(fwd_workspace.scale() * inv_workspace.scale(),
			scalar_type(1));
#else
    (void)src;
    return true;
#endif
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type rows = dst.size(2, 0);
    length_type cols = dst.size(2, 1);
    Matrix<T> tmp(1, cols);

    Vector<T, VecBlockT> w  (
      const_cast<VecBlockT&>(src.functor().block().get_vblk()));
    Matrix<T, MatBlockT> in (
      const_cast<MatBlockT&>(src.functor().block().get_mblk().functor().block()));
    Matrix<T, DstBlock>        out(dst);

    Workspace2T const& fwd_workspace(
      src.functor().block().get_mblk().functor().workspace());

    Backend2T&         fwd_backend  (const_cast<Backend2T&>(
      src.functor().block().get_mblk().functor().backend()) );

    Workspace1T const& inv_workspace(src.functor().workspace());
    Backend1T&         inv_backend  (const_cast<Backend1T&>(src.functor().backend()));

    for (index_type r=0; r<rows; ++r)
    {
      fwd_workspace.by_reference(&fwd_backend,
				 in (Domain<2>(Domain<1>(r, 1, 1), cols)),
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)) );
      tmp.row(0) *= w;
      inv_workspace.by_reference(&inv_backend,
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)),
				 out(Domain<2>(Domain<1>(r, 1, 1), cols)) );
    }
  }
};



// Evaluator for matrix of coefficients (weights * fftm)

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
    Fc_expr_tag
  >
{
  static char const* name() { return "Fc_expr_tag-mwl"; }

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

  static bool const ct_valid = true;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    (void)dst;
    (void)src;
    return true;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type rows = dst.size(2, 0);
    length_type cols = dst.size(2, 1);
    Matrix<T> tmp(1, cols);

    Matrix<T, CoeffsMatBlockT> w 
      (const_cast<CoeffsMatBlockT&>(src.functor().block().left()));
    Matrix<T, MatBlockT> in 
      (const_cast<MatBlockT&>(src.functor().block().right().functor().block()));
    Matrix<T, DstBlock> out(dst);

    Workspace2T const& fwd_workspace(
      src.functor().block().right().functor().workspace());

    Backend2T&         fwd_backend  (const_cast<Backend2T&>(
      src.functor().block().right().functor().backend()) );

    Workspace1T const& inv_workspace(src.functor().workspace());
    Backend1T&         inv_backend  (const_cast<Backend1T&>(src.functor().backend()));

    for (index_type r=0; r<rows; ++r)
    {
      fwd_workspace.by_reference(&fwd_backend,
				 in (Domain<2>(Domain<1>(r, 1, 1), cols)),
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)) );
      tmp.row(0) *= w.row(r);
      inv_workspace.by_reference(&inv_backend,
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)),
				 out(Domain<2>(Domain<1>(r, 1, 1), cols)) );
    }
  }
};



// Evaluator for matrix of coefficients (fftm * weights)

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
    Fc_expr_tag
  >
{
  static char const* name() { return "Fc_expr_tag-mwr"; }

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

  static bool const ct_valid = true;

  static bool rt_valid(DstBlock& dst, SrcBlock const& src)
  {
    (void)dst;
    (void)src;
    return true;
  }
  
  static void exec(DstBlock& dst, SrcBlock const& src)
  {
    length_type rows = dst.size(2, 0);
    length_type cols = dst.size(2, 1);
    Matrix<T> tmp(1, cols);

    Matrix<T, CoeffsMatBlockT> w 
      (const_cast<CoeffsMatBlockT&>(src.functor().block().right()));
    Matrix<T, MatBlockT> in 
      (const_cast<MatBlockT&>(src.functor().block().left().functor().block()));
    Matrix<T, DstBlock> out(dst);

    Workspace2T const& fwd_workspace(
      src.functor().block().left().functor().workspace());

    Backend2T&         fwd_backend  (const_cast<Backend2T&>(
      src.functor().block().left().functor().backend()) );

    Workspace1T const& inv_workspace(src.functor().workspace());
    Backend1T&         inv_backend  (const_cast<Backend1T&>(src.functor().backend()));

    for (index_type r=0; r<rows; ++r)
    {
      fwd_workspace.by_reference(&fwd_backend,
				 in (Domain<2>(Domain<1>(r, 1, 1), cols)),
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)) );
      tmp.row(0) *= w.row(r);
      inv_workspace.by_reference(&inv_backend,
				 tmp(Domain<2>(Domain<1>(0, 1, 1), cols)),
				 out(Domain<2>(Domain<1>(r, 1, 1), cols)) );
    }
  }
};

} // namespace vsip::impl

} // namespace vsip

#endif // VSIP_OPT_EXPR_EVAL_FASTCONV_HPP
