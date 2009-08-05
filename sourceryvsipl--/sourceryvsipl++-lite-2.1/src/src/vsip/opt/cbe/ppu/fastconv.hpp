/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/fastconv.hpp
    @author  Don McCoy
    @date    2007-02-23
    @brief   VSIPL++ Library: Wrapper for fast convolution on the SPEs.
*/

#ifndef VSIP_OPT_CBE_PPU_FASTCONV_H
#define VSIP_OPT_CBE_PPU_FASTCONV_H

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/allocation.hpp>
#include <vsip/core/config.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/opt/cbe/fconv_params.h>
#include <vsip/opt/cbe/ppu/bindings.hpp>
#include <vsip/opt/cbe/ppu/task_manager.hpp>
extern "C"
{
#include <libspe2.h>
}
/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cbe
{

template <dimension_type D,
          typename T,
	  typename ComplexFmt> 
struct Fastconv_traits;
template <>
struct Fastconv_traits<1, std::complex<float>, Cmplx_inter_fmt>
{
  typedef Fastconv_tag tag_type;
  static length_type const min_size = VSIP_IMPL_MIN_FCONV_SIZE;
  static length_type const max_size = VSIP_IMPL_MAX_FCONV_SIZE;
};
template <>
struct Fastconv_traits<1, std::complex<float>, Cmplx_split_fmt>
{
  typedef Fastconv_tag tag_type;
  static length_type const min_size = VSIP_IMPL_MIN_FCONV_SPLIT_SIZE;
  static length_type const max_size = VSIP_IMPL_MAX_FCONV_SPLIT_SIZE;
};
template <>
struct Fastconv_traits<2, std::complex<float>, Cmplx_inter_fmt>
{
  typedef Fastconvm_tag tag_type;
  static length_type const min_size = VSIP_IMPL_MIN_FCONV_SIZE;
  static length_type const max_size = VSIP_IMPL_MAX_FCONV_SIZE;
};
template <>
struct Fastconv_traits<2, std::complex<float>, Cmplx_split_fmt>
{
  typedef Fastconvm_tag tag_type;
  static length_type const min_size = VSIP_IMPL_MIN_FCONV_SPLIT_SIZE;
  static length_type const max_size = VSIP_IMPL_MAX_FCONV_SPLIT_SIZE;
};

// Fast convolution object using SPEs to perform computation.

// Template parameters:
//   T to be the value type of data that will be processed.
//   COMPLEXFMT to be the complex format (either Cmplx_inter_fmt or
//     Cmplx_split_fmt) to be processed.

template <dimension_type D,
          typename       T,
	  typename       ComplexFmt>
class Fastconv_base
{
  static dimension_type const dim = D;

  typedef ComplexFmt complex_type;
  typedef Layout<1, row1_type, Stride_unit_dense, complex_type> layout1_type;
  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> layout2_type;

public:
  Fastconv_base(length_type const input_size, bool transform_kernel)
    : size_            (input_size),
      transform_kernel_(transform_kernel),
      instance_id_     (++instance_id_counter_)
  {
    assert(rt_valid_size(this->size_));
  }

  static bool rt_valid_size(length_type size)
  {
    return (size >= cbe::Fastconv_traits<dim, T, complex_type>::min_size &&
	    size <= cbe::Fastconv_traits<dim, T, complex_type>::max_size &&
	    fft::is_power_of_two(size));
  }


  template <typename Block0, typename Block1, typename Block2>
  void convolve(const_Vector<T, Block0> in, const_Vector<T, Block1> kernel, Vector<T, Block2> out)
  {
    Ext_data<Block0, layout1_type> ext_in    (in.block(),     SYNC_IN);
    Ext_data<Block1, layout1_type> ext_kernel(kernel.block(), SYNC_IN);
    Ext_data<Block2, layout1_type> ext_out   (out.block(),    SYNC_OUT);
    assert(dim == 1);
    assert(ext_in.stride(0) == 1);
    assert(ext_kernel.stride(0) == 1);
    assert(ext_out.stride(0) == 1);

    length_type rows = 1;
    fconv(ext_in.data(), ext_kernel.data(), ext_out.data(), rows, out.size(0));
  }

  template <typename Block0, typename Block1, typename Block2>
  void convolve(const_Matrix<T, Block0> in, const_Vector<T, Block1> kernel, Matrix<T, Block2> out)
  {
    Ext_data<Block0, layout2_type> ext_in    (in.block(),     SYNC_IN);
    Ext_data<Block1, layout1_type> ext_kernel(kernel.block(), SYNC_IN);
    Ext_data<Block2, layout2_type> ext_out   (out.block(),    SYNC_OUT);
    assert(dim == 1);
    assert(ext_in.stride(1) == 1);
    assert(ext_kernel.stride(0) == 1);
    assert(ext_out.stride(1) == 1);

    length_type rows = in.size(0);
    fconv(ext_in.data(), ext_kernel.data(), ext_out.data(), rows, out.size(1));
  }

  template <typename Block0, typename Block1, typename Block2>
  void convolve(const_Matrix<T, Block0> in, const_Matrix<T, Block1> kernel, Matrix<T, Block2> out)
  {
    Ext_data<Block0, layout2_type> ext_in    (in.block(),     SYNC_IN);
    Ext_data<Block1, layout2_type> ext_kernel(kernel.block(), SYNC_IN);
    Ext_data<Block2, layout2_type> ext_out   (out.block(),    SYNC_OUT);
    assert(dim == 2);
    assert(ext_in.stride(1) == 1);
    assert(ext_kernel.stride(1) == 1);
    assert(ext_out.stride(1) == 1);

    length_type rows = in.size(0);
    fconv(ext_in.data(), ext_kernel.data(), ext_out.data(), rows, out.size(1));
  }

  length_type size() { return size_; }

private:
  typedef typename Scalar_of<T>::type uT;
  void fconv(T* in, T* kernel, T* out, length_type rows, length_type length);
  void fconv(std::pair<uT*,uT*> in, std::pair<uT*,uT*> kernel,
	     std::pair<uT*,uT*> out, length_type rows, length_type length);

  // Member data.
  length_type size_;
  bool transform_kernel_;
  unsigned int instance_id_;

  // This counter is used to give each instance of this type
  // a unique ID that is passed to the SPEs as part of the 
  // parameter block.  When an SPE sees that the ID has changed,
  // it knows to reload the weights.
  static unsigned int instance_id_counter_;
};



template <dimension_type D,
          typename T,
	  typename ComplexFmt = Cmplx_inter_fmt>
class Fastconv;

template <typename T, typename ComplexFmt>
class Fastconv<1, T, ComplexFmt> : public Fastconv_base<1, T, ComplexFmt>
{
  // Constructors, copies, assignments, and destructors.
public:

  template <typename Block>
  Fastconv(Vector<T, Block> coeffs,
           length_type input_size,
	   bool transform_kernel = true)
    VSIP_THROW((std::bad_alloc))
    : Fastconv_base<1, T, ComplexFmt>(input_size, transform_kernel),
      kernel_(input_size)
  {
    assert(coeffs.size(0) <= this->size());
    if (transform_kernel)
    {
      kernel_ = T();
      kernel_(view_domain(coeffs.local())) = coeffs.local();
    }
    else
      kernel_ = coeffs.local();
  }
  ~Fastconv() VSIP_NOTHROW {}

  // Fastconv operators.
  template <typename Block1,
	    typename Block2>
  Vector<T, Block2>
  operator()(
    const_Vector<T, Block1> in, 
    Vector<T, Block2>       out)
      VSIP_NOTHROW
  {
    assert(in.size() == this->size());
    assert(out.size() == this->size());
    
    this->convolve(in.local(), this->kernel_, out.local());
    
    return out;
  }

  template <typename Block1,
	    typename Block2>
  Matrix<T, Block2>
  operator()(
    const_Matrix<T, Block1> in, 
    Matrix<T, Block2>       out)
    VSIP_NOTHROW
  {
    assert(in.size(1) == this->size());
    assert(out.size(1) == this->size());

    this->convolve(in.local(), this->kernel_, out.local());
    
    return out;
  }

private:
  typedef ComplexFmt complex_type;
  typedef Layout<1, row1_type, 
                 Stride_unit_dense, complex_type>   kernel_layout_type;
  typedef Fast_block<1, T, 
                     kernel_layout_type, Local_map> kernel_block_type;
  typedef Vector<T, kernel_block_type>              kernel_view_type;

  // Member data.
  kernel_view_type kernel_;
};



template <typename T, typename ComplexFmt>
class Fastconv<2, T, ComplexFmt> : public Fastconv_base<2, T, ComplexFmt>
{
  // Constructors, copies, assignments, and destructors.
public:

  template <typename Block>
  Fastconv(Matrix<T, Block> coeffs,
           length_type input_size,
	   bool transform_kernel = true)
    VSIP_THROW((std::bad_alloc))
    : Fastconv_base<2, T, ComplexFmt>(input_size, transform_kernel),
      kernel_(coeffs.local().size(0), input_size)
  {
    assert(coeffs.size(1) <= this->size());
    if (transform_kernel)
    {
      kernel_ = T();
      kernel_(view_domain(coeffs.local())) = coeffs.local();
    }
    else
      kernel_ = coeffs.local();
  }
  ~Fastconv() VSIP_NOTHROW {}

  // Fastconv operators.
  template <typename Block1,
	    typename Block2>
  Vector<T, Block2>
  operator()(
    const_Vector<T, Block1> in, 
    Vector<T, Block2>       out)
      VSIP_NOTHROW
  {
    assert(in.size() == this->size());
    assert(out.size() == this->size());
    
    this->convolve(in.local(), this->kernel_, out.local());
    
    return out;
  }

  template <typename Block1,
	    typename Block2>
  Matrix<T, Block2>
  operator()(
    const_Matrix<T, Block1> in, 
    Matrix<T, Block2>       out)
    VSIP_NOTHROW
  {
    assert(in.size(1) == this->size());
    assert(out.size(1) == this->size());
    
    this->convolve(in.local(), this->kernel_, out.local());
    
    return out;
  }

private:
  // Member data.
  typedef ComplexFmt complex_type;
  typedef Layout<2, row2_type, 
                 Stride_unit_dense, complex_type>   kernel_layout_type;
  typedef Fast_block<2, T, 
                     kernel_layout_type, Local_map> kernel_block_type;
  typedef Matrix<T, kernel_block_type>              kernel_view_type;

  kernel_view_type kernel_;
};


} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CBE_PPU_FASTCONV_H
