/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fft/backend.hpp
    @author  Stefan Seefeld
    @date    2006-02-17
    @brief   VSIPL++ Library: fft backend interfaces.

    This file defines fft backend interfaces to be implemented by
    third-party library bridges.
*/

#ifndef VSIP_CORE_FFT_BACKEND_HPP
#define VSIP_CORE_FFT_BACKEND_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/layout.hpp>
#include <vsip/core/metaprogramming.hpp>

namespace vsip
{
namespace impl
{
namespace fft
{

bool const forward = true;
bool const inverse = false;
  
template <dimension_type D, typename I, typename O, int A, int E> 
class backend;

template <dimension_type D, typename T>
struct backend_base
{
  static dimension_type const dim = D;
  typedef T scalar_type;
};

/// 1D real forward FFT
template <typename T, int A>
class backend<1, T, std::complex<T>, A, -1>
  : public backend_base<1, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-1D-real-forward"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<1> &) { return false;}
  /// real -> complex (interleaved)
  virtual void by_reference(T *in, stride_type in_stride,
			    std::complex<T> *out, stride_type out_stride,
			    length_type length) = 0;
  /// real -> complex (split)
  virtual void by_reference(T *in, stride_type in_stride,
			    std::pair<T *, T *> out, stride_type out_stride,
			    length_type length) = 0;
};

/// 1D real inverse FFT
template <typename T, int A>
class backend<1, std::complex<T>, T, A, 1>
  : public backend_base<1, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-1D-real-inverse"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<1> &) { return false;}
  /// complex (interleaved) -> real
  virtual void by_reference(std::complex<T> *in, stride_type in_stride,
			    T *out, stride_type out_stride,
			    length_type length) = 0;
  /// real -> complex (split)
  virtual void by_reference(std::pair<T *, T *> in, stride_type in_stride,
			    T *out, stride_type out_stride,
			    length_type length) = 0;
};

/// 1D complex FFT
template <typename T, int A, int E>
class backend<1, std::complex<T>, std::complex<T>, A, E>
  : public backend_base<1, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-1D-complex"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<1> &rtl_inout)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_inout.pack = stride_unit_dense;
    rtl_inout.order = tuple<0, 1, 2>();
    rtl_inout.complex = cmplx_inter_fmt;
  }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<1> &) { return false;}
  /// complex (interleaved) in-place
  virtual void in_place(std::complex<T> *, stride_type, length_type) = 0;
  /// complex (split) in-place
  virtual void in_place(std::pair<T *, T *>, stride_type, length_type) = 0;
  /// complex (interleaved) by-reference
  virtual void by_reference(std::complex<T> *in, stride_type in_stride,
			    std::complex<T> *out, stride_type out_stride,
			    length_type length) = 0;
  /// complex (split) by-reference
  virtual void by_reference(std::pair<T *, T *>, stride_type in_stride,
			    std::pair<T *, T *>, stride_type out_stride,
			    length_type length) = 0;
};

/// 2D real forward FFT
template <typename T, int A>
class backend<2, T, std::complex<T>, A, -1>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-real-forward"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// real -> complex (interleaved) by-reference
  virtual void by_reference(T *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::complex<T> *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// real -> complex (split) by-reference
  virtual void by_reference(T *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::pair<T *, T *>,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

/// 2D real inverse FFT
template <typename T, int A>
class backend<2, std::complex<T>, T, A, 1>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-real-inverse"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// complex (interleaved) -> real by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    T *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// complex (split) -> real by-reference
  virtual void by_reference(std::pair<T *, T *>,
			    stride_type in_r_stride, stride_type in_c_stride,
			    T *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

/// 2D complex FFT
template <typename T, int A, int E>
class backend<2, std::complex<T>, std::complex<T>, A, E>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-complex"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_inout)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_inout.pack = stride_unit_dense;
    rtl_inout.order = tuple<0, 1, 2>();
    rtl_inout.complex = cmplx_inter_fmt;
  }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// complex (interleaved) in-place
  virtual void in_place(std::complex<T> *inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols) = 0;
  /// complex (split) in-place
  virtual void in_place(std::pair<T *, T *>,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols) = 0;
  /// complex (interleaved) by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::complex<T> *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// complex (split) by-reference
  virtual void by_reference(std::pair<T *, T *> in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::pair<T *, T *> out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

/// 3D real forward FFT
template <typename T, int A>
class backend<3, T, std::complex<T>, A, -1>
  : public backend_base<3, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-real-forward"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<3> &) { return false;}
  /// real -> complex (interleaved) by-reference
  virtual void by_reference(T *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    std::complex<T> *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
  /// real -> complex (split) by-reference
  virtual void by_reference(T *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    std::pair<T *, T *> out,
			    stride_type out_x_stride,
			    stride_type out_y_stridey,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
};

/// 3D real inverse FFT
template <typename T, int A>
class backend<3, std::complex<T>, T, A, 1>
  : public backend_base<3, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-real-inverse"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<3> &) { return false;}
  /// complex (interleaved) -> real by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    T *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
  /// complex (split) -> real by-reference
  virtual void by_reference(std::pair<T *, T *> in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    T *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
};

/// 3D complex FFT
template <typename T, int A, int E>
class backend<3, std::complex<T>, std::complex<T>, A, E>
  : public backend_base<3, T>
{
public:
  static int const axis = A;

  virtual ~backend() {}
  virtual char const* name() { return "fft-backend-2D-complex"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<3> &rtl_inout)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_inout.pack = stride_unit_dense;
    rtl_inout.order = tuple<0, 1, 2>();
    rtl_inout.complex = cmplx_inter_fmt;
  }
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    // By default use unit_stride, tuple<0, 1, 2>, cmplx_inter_fmt
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    rtl_in.order = rtl_out.order = tuple<0, 1, 2>();
    rtl_in.complex = rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<3> &) { return false;}
  /// complex (interleaved) in-place
  virtual void in_place(std::complex<T> *inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length) = 0;
  /// complex (split) in-place
  virtual void in_place(std::pair<T *, T *> inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length) = 0;
  /// complex (interleaved) by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    std::complex<T> *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
  /// complex (split) by-reference
  virtual void by_reference(std::pair<T *, T *> in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    std::pair<T *, T *> out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length) = 0;
};

/// FFTM
template <typename I, typename O, int A, int E> class fftm;

/// real forward FFTM
template <typename T, int A>
class fftm<T, std::complex<T>, A, -1>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~fftm() {}
  virtual char const* name() { return "fftm-backend-real-forward"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride,
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    // an ordering that gives unit strides on the axis perpendicular to A,
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_out.order = rtl_in.order;
    // and interleaved complex.
    rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// real -> complex (interleaved) by-reference
  virtual void by_reference(T *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::complex<T> *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// real -> complex (split) by-reference
  virtual void by_reference(T *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::pair<T *, T *> out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

/// real inverse FFTM
template <typename T, int A>
class fftm<std::complex<T>, T, A, 1>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~fftm() {}
  virtual char const* name() { return "fftm-backend-real-inverse"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride,
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    // an ordering that gives unit strides on the axis perpendicular to A,
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_out.order = rtl_in.order;
    // and interleaved complex.
    rtl_in.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// complex (interleaved) -> real by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    T *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// complex (split) -> real by-reference
  virtual void by_reference(std::pair<T *, T *> in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    T *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

/// complex FFTM
template <typename T, int A, int E>
class fftm<std::complex<T>, std::complex<T>, A, E>
  : public backend_base<2, T>
{
public:
  static int const axis = A;

  virtual ~fftm() {}
  virtual char const* name() { return "fftm-backend-complex"; }
  virtual bool supports_scale() { return false;}
  virtual void query_layout(Rt_layout<2> &rtl_inout)
  {
    // By default use unit_stride,
    rtl_inout.pack = stride_unit_dense;
    // an ordering that gives unit strides on the axis perpendicular to A,
    if (A == 0) rtl_inout.order = tuple<1, 0, 2>();
    else rtl_inout.order = tuple<0, 1, 2>();
    // and interleaved complex.
    rtl_inout.complex = cmplx_inter_fmt;
  }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    // By default use unit_stride,
    rtl_in.pack = rtl_out.pack = stride_unit_dense;
    // an ordering that gives unit strides on the axis perpendicular to A,
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_out.order = rtl_in.order;
    // and interleaved complex.
    rtl_in.complex = rtl_out.complex = cmplx_inter_fmt;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return false;}
  /// complex (interleaved) in-place
  virtual void in_place(std::complex<T> *inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols) = 0;
  /// complex (split) in-place
  virtual void in_place(std::pair<T *, T *> inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols) = 0;
  /// complex (interleaved) by-reference
  virtual void by_reference(std::complex<T> *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::complex<T> *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
  /// complex (split) by-reference
  virtual void by_reference(std::pair<T *, T *> in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    std::pair<T *, T *> out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols) = 0;
};

} // namespace vsip::impl::fft
} // namespace vsip::impl
} // namespace vsip

#endif
