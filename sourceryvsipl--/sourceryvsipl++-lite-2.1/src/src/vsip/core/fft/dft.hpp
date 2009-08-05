/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fft/dft.hpp
    @author  Stefan Seefeld
    @date    2006-05-01
    @brief   VSIPL++ Library: DFT backend.
*/

#ifndef VSIP_CORE_FFT_DFT_HPP
#define VSIP_CORE_FFT_DFT_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/fft/factory.hpp>
#include <vsip/core/fft/util.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace fft
{
namespace
{
template <typename T>
inline vsip::complex<T>
sin_cos(double phi)
{
  return vsip::complex<T>(cos(phi), sin(phi));
}

template <typename T>
std::pair<T*,T*> offset(std::pair<T*,T*> data, int o)
{
  return std::make_pair(data.first + o, data.second + o);
}

}
template <dimension_type D, typename I, typename O, int A, int E> class dft;

// 1D complex -> complex DFT
template <typename T, int A, int E>
class dft<1, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<1, std::complex<T>, std::complex<T>, A, E>

{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-1D-complex"; }
  virtual void query_layout(Rt_layout<1> &) {}
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void in_place(ctype *inout, stride_type s, length_type l)
  {
    typedef double AccT;
    aligned_array<std::complex<T> > tmp(l);
    AccT const phi = E * 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l; ++w)
    {
      complex<AccT> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<AccT>(inout[k * s]) * sin_cos<AccT>(phi * k * w);
      tmp[w] = sum;
    }
    for (index_type w = 0; w < l; ++w) inout[w * s] = tmp[w];
  }
  virtual void in_place(ztype inout, stride_type s, length_type l)
  {
    typedef double AccT;
    aligned_array<std::complex<T> > tmp(l);
    AccT const phi = E * 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l; ++w)
    {
      complex<T> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<AccT>(inout.first[k * s], inout.second[k * s])
	  * sin_cos<AccT>(phi * k * w);
      tmp[w] = sum;
    }
    for (index_type w = 0; w < l; ++w)
    {
      inout.first[w * s] = tmp[w].real();
      inout.second[w * s] = tmp[w].imag();
    }
  }
  virtual void by_reference(ctype *in, stride_type in_s,
			    ctype *out, stride_type out_s,
			    length_type l)
  {
    typedef double AccT;
    AccT const phi = E * 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l; ++w)
    {
      complex<AccT> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<AccT>(in[k * in_s]) * sin_cos<AccT>(phi * k * w);
      out[w * out_s] = ctype(sum);
    }
  }
  virtual void by_reference(ztype in, stride_type in_s,
			    ztype out, stride_type out_s,
			    length_type l)
  {
    typedef double AccT;
    AccT const phi = E * 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l; ++w)
    {
      complex<AccT> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<AccT>(in.first[k * in_s], in.second[k * in_s])
	  * sin_cos<AccT>(phi * k * w);
      out.first[w * out_s] = sum.real();
      out.second[w * out_s] = sum.imag();
    }
  }
};

// 1D real -> complex DFT
template <typename T, int A>
class dft<1, T, std::complex<T>, A, -1>
  : public fft::backend<1, T, std::complex<T>, A, -1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-1D-real-forward"; }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(rtype *in, stride_type in_s,
			    ctype *out, stride_type out_s,
			    length_type l)
  {
    T const phi = - 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l/2 + 1; ++w)
    {
      complex<T> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<T>(in[k * in_s]) * sin_cos<T>(phi * k * w);
      out[w * out_s] = sum;
    }
  }
  virtual void by_reference(rtype *in, stride_type in_s,
			    ztype out, stride_type out_s,
			    length_type l)
  {
    T const phi = - 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l/2 + 1; ++w)
    {
      complex<T> sum;
      for (index_type k = 0; k < l; ++k)
	sum += vsip::complex<T>(in[k * in_s]) * sin_cos<T>(phi * k * w);
      out.first[w * out_s] = sum.real();
      out.second[w * out_s] = sum.imag();
    }
  }
};

// 1D complex -> real DFT
template <typename T, int A>
class dft<1, std::complex<T>, T, A, 1>
  : public fft::backend<1, std::complex<T>, T, A, 1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-1D-real-inverse"; }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(ctype *in, stride_type in_s,
			    rtype *out, stride_type out_s,
			    length_type l)
  {
    T const phi = 2.0 * VSIP_IMPL_PI/l;

    for (index_type w = 0; w < l; ++w)
    {
      complex<T> sum;
      for (index_type k = 0; k < l/2 + 1; ++k)
	sum += in[k * in_s] * sin_cos<T>(phi * k * w);
      for (index_type k = l/2 + 1; k < l; ++k)
	sum += conj(in[(l - k) * in_s]) * sin_cos<T>(phi * k * w);
      out[w * out_s] = sum.real();
    }
  }
  virtual void by_reference(ztype in, stride_type in_s,
			    rtype *out, stride_type out_s,
			    length_type l)
  {
    T const phi = 2.0 * VSIP_IMPL_PI/l;
    
    for (index_type w = 0; w < l; ++w)
    {
      complex<T> sum;
      for (index_type k = 0; k < l/2 + 1; ++k)
	sum += complex<T>(in.first[k * in_s], in.second[k * in_s])
	  * sin_cos<T>(phi * k * w);
      for (index_type k = l/2 + 1; k < l; ++k)
	sum += complex<T>(in.first[(l - k) * in_s], -in.second[(l - k) * in_s])
	  * sin_cos<T>(phi * k * w);
      out[w * out_s] = sum.real();
    }
  }
};

// 2D complex -> complex DFT
template <typename T, int A, int E>
class dft<2, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<2, std::complex<T>, std::complex<T>, A, E>

{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-2D-complex"; }
  virtual void query_layout(Rt_layout<2> &) {}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void in_place(ctype *inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (length_type r = 0; r != rows; ++r)
      dft_1d.in_place(inout + r * r_stride, c_stride, cols);
    for (length_type c = 0; c != cols; ++c)
      dft_1d.in_place(inout + c * c_stride, r_stride, rows);
  }
  virtual void in_place(ztype inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (length_type r = 0; r != rows; ++r)
    {
      ztype line = std::make_pair(inout.first + r * r_stride,
				  inout.second + r * r_stride);
      dft_1d.in_place(line, c_stride, cols);
    }
    for (length_type c = 0; c != cols; ++c)
    {
      ztype line = std::make_pair(inout.first + c * c_stride,
				  inout.second + c * c_stride);
      dft_1d.in_place(line, r_stride, rows);
    }
  }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (length_type r = 0; r != rows; ++r)
      dft_1d.by_reference(in + r * in_r_stride, in_c_stride,
			  out + r * out_r_stride, out_c_stride, cols);
    for (length_type c = 0; c != cols; ++c)
      dft_1d.in_place(out + c * out_c_stride, out_r_stride, rows);
  }
  virtual void by_reference(ztype in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (length_type r = 0; r != rows; ++r)
    {
      ztype in_line = std::make_pair(in.first + r * in_r_stride,
				     in.second + r * in_r_stride);
      ztype out_line = std::make_pair(out.first + r * out_r_stride,
				      out.second + r * out_r_stride);
      dft_1d.by_reference(in_line, in_c_stride,
			  out_line, out_c_stride, cols);
    }
    for (length_type c = 0; c != cols; ++c)
    {
      ztype line = std::make_pair(out.first + c * out_c_stride,
				  out.second + c * out_c_stride);
      dft_1d.in_place(line, out_r_stride, rows);
    }
  }
};

// 2D real -> complex DFT
template <typename T, int A>
class dft<2, T, std::complex<T>, A, -1>
  : public fft::backend<2, T, std::complex<T>, A, -1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-2D-real-forward"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, rtype, ctype, 0, -1> rdft_1d;
    dft<1, ctype, ctype, 0, -1> dft_1d;
    if (A == 0)
    {
      for (length_type c = 0; c != cols; ++c)
	rdft_1d.by_reference(in + c * in_c_stride, in_r_stride,
			     out + c * out_c_stride, out_r_stride, rows);
      for (length_type r = 0; r != rows/2 + 1; ++r)
	dft_1d.in_place(out + r * out_r_stride, out_c_stride, cols);
    }
    else
    {
      for (length_type r = 0; r != rows; ++r)
	rdft_1d.by_reference(in + r * in_r_stride, in_c_stride,
			     out + r * out_r_stride, out_c_stride, cols);
      for (length_type c = 0; c != cols/2 + 1; ++c)
	dft_1d.in_place(out + c * out_c_stride, out_r_stride, rows);
    }
  }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, rtype, ctype, 0, -1> rdft_1d;
    dft<1, ctype, ctype, 0, -1> dft_1d;
    if (A == 0)
    {
      for (length_type c = 0; c != cols; ++c)
      {
	ztype line = std::make_pair(out.first + c * out_c_stride,
				    out.second + c * out_c_stride);
	rdft_1d.by_reference(in + c * in_c_stride, in_r_stride,
			     line, out_r_stride, rows);
      }
      for (length_type r = 0; r != rows/2 + 1; ++r)
      {
	ztype line = std::make_pair(out.first + r * out_r_stride,
				    out.second + r * out_r_stride);
	dft_1d.in_place(line, out_c_stride, cols);
      }
    }
    else
    {
      for (length_type r = 0; r != rows; ++r)
      {
	ztype line = std::make_pair(out.first + r * out_r_stride,
				    out.second + r * out_r_stride);
	rdft_1d.by_reference(in + r * in_r_stride, in_c_stride,
			     line, out_c_stride, cols);
      }
      for (length_type c = 0; c != cols/2 + 1; ++c)
      {
	ztype line = std::make_pair(out.first + c * out_c_stride,
				    out.second + c * out_c_stride);
	dft_1d.in_place(line, out_r_stride, rows);
      }
    }
  }

};

// 2D complex -> real DFT
template <typename T, int A>
class dft<2, std::complex<T>, T, A, 1>
  : public fft::backend<2, std::complex<T>, T, A, 1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-2D-real-inverse"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, 1> dft_1d;
    dft<1, ctype, rtype, 0, 1> rdft_1d;
    if (A == 0)
    {
      length_type rows2 = rows/2 + 1;
      aligned_array<ctype> tmp(rows2 * cols); // row-major temp matrix.
      for (length_type r = 0; r != rows2; ++r)
	dft_1d.by_reference(in + r * in_r_stride, in_c_stride,
			    tmp.get() + r * cols, 1, cols);
      for (length_type c = 0; c != cols; ++c)
	rdft_1d.by_reference(tmp.get() + c, cols,
			     out + c * out_c_stride, out_r_stride, rows);
    }
    else
    {
      length_type cols2 = cols/2 + 1;
      aligned_array<ctype> tmp(rows * cols2); // col-major temp matrix.
      for (length_type c = 0; c != cols2; ++c)
	dft_1d.by_reference(in + c * in_c_stride, in_r_stride,
			    tmp.get() + c * rows, 1, rows);
      for (length_type r = 0; r != rows; ++r)
	rdft_1d.by_reference(tmp.get() + r, rows,
			     out + r * out_r_stride, out_c_stride, cols);
    }
  }
  virtual void by_reference(ztype in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, 1> dft_1d;
    dft<1, ctype, rtype, 0, 1> rdft_1d;
    if (A == 0)
    {
      length_type rows2 = rows/2 + 1;
      aligned_array<rtype> tmp_r(rows2 * cols); // col-major temp real matrix.
      aligned_array<rtype> tmp_i(rows2 * cols); // col-major temp imag matrix.
      for (length_type r = 0; r != rows2; ++r)
      {
	ztype line = std::make_pair(tmp_r.get() + r * cols,
				    tmp_i.get() + r * cols);
	dft_1d.by_reference(offset(in, r * in_r_stride), in_c_stride,
			    line, 1, cols);
      }
      for (length_type c = 0; c != cols; ++c)
      {
	ztype line = std::make_pair(tmp_r.get() + c,
				    tmp_i.get() + c);
	rdft_1d.by_reference(line, cols,
			     out + c * out_c_stride, out_r_stride, rows);
      }
    }
    else
    {
      length_type cols2 = cols/2 + 1;
      aligned_array<rtype> tmp_r(rows * cols2); // col-major temp real matrix.
      aligned_array<rtype> tmp_i(rows * cols2); // col-major temp imag matrix.
      for (length_type c = 0; c != cols2; ++c)
      {
	ztype line = std::make_pair(tmp_r.get() + c * rows,
				    tmp_i.get() + c * rows);
	dft_1d.by_reference(offset(in, c * in_c_stride), in_r_stride,
			    line, 1, rows);
      }
      for (length_type r = 0; r != rows; ++r)
      {
	ztype line = std::make_pair(tmp_r.get() + r,
				    tmp_i.get() + r);
	rdft_1d.by_reference(line, rows,
			     out + r * out_r_stride, out_c_stride, cols);
      }
    }
  }

};

// 3D complex -> complex DFT
template <typename T, int A, int E>
class dft<3, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<3, std::complex<T>, std::complex<T>, A, E>

{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-3D-complex"; }
  virtual void query_layout(Rt_layout<3> &) {}
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void in_place(ctype *inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length)
  {
    dft<2, ctype, ctype, 0, E> dft_2d;
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (index_type x = 0; x != x_length; ++x)
      dft_2d.in_place(inout + x * x_stride,
		      y_stride, z_stride, y_length, z_length);
    for (index_type y = 0; y != y_length; ++y)
      for (index_type z = 0; z != z_length; ++z)
	dft_1d.in_place(inout + y * y_stride + z * z_stride,
			x_stride, x_length);
  }
  virtual void in_place(ztype inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length)
  {
    dft<2, ctype, ctype, 0, E> dft_2d;
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (index_type x = 0; x != x_length; ++x)
      dft_2d.in_place(offset(inout, x * x_stride),
		      y_stride, z_stride, y_length, z_length);
    for (index_type y = 0; y != y_length; ++y)
      for (index_type z = 0; z != z_length; ++z)
	dft_1d.in_place(offset(inout, y * y_stride + z * z_stride),
			x_stride, x_length);
  }
  virtual void by_reference(ctype *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    ctype *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<2, ctype, ctype, 0, E> dft_2d;
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (index_type x = 0; x != x_length; ++x)
      dft_2d.by_reference(in + x * in_x_stride,
			  in_y_stride, in_z_stride,
			  out + x * out_x_stride,
			  out_y_stride, out_z_stride,
			  y_length, z_length);
    for (index_type y = 0; y != y_length; ++y)
      for (index_type z = 0; z != z_length; ++z)
	dft_1d.in_place(out + y * out_y_stride + z * out_z_stride,
			out_x_stride, x_length);
  }
  virtual void by_reference(ztype in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    ztype out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<2, ctype, ctype, 0, E> dft_2d;
    dft<1, ctype, ctype, 0, E> dft_1d;
    for (index_type x = 0; x != x_length; ++x)
      dft_2d.by_reference(offset(in, x * in_x_stride),
			  in_y_stride, in_z_stride,
			  offset(out, x * out_x_stride),
			  out_y_stride, out_z_stride,
			  y_length, z_length);
    for (index_type y = 0; y != y_length; ++y)
      for (index_type z = 0; z != z_length; ++z)
	dft_1d.in_place(offset(out, y * out_y_stride + z * out_z_stride),
			out_x_stride, x_length);
  }
};

// 3D real -> complex DFT
template <typename T, int A>
class dft<3, T, std::complex<T>, A, -1>
  : public fft::backend<3, T, std::complex<T>, A, -1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-3D-real-forward"; }
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(rtype *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    ctype *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<1, rtype, ctype, 0, -1> rdft_1d;
    dft<2, ctype, ctype, 0, -1> dft_2d;
    if (A == 0)
    {
      for (length_type y = 0; y != y_length; ++y)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(in + y * in_y_stride + z * in_z_stride,
			       in_x_stride,
			       out + y * out_y_stride + z * out_z_stride,
			       out_x_stride, x_length);
      for (length_type x = 0; x != x_length/2 + 1; ++x)
	dft_2d.in_place(out + x * out_x_stride,
			out_y_stride, out_z_stride,
			y_length, z_length);
    }
    else if (A == 1)
    {
      for (length_type x = 0; x != x_length; ++x)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(in + x * in_x_stride + z * in_z_stride,
			       in_y_stride,
			       out + x * out_x_stride + z * out_z_stride,
			       out_y_stride, y_length);
      for (length_type y = 0; y != y_length/2 + 1; ++y)
	dft_2d.in_place(out + y * out_y_stride,
			out_x_stride, out_z_stride,
			x_length, z_length);
    }
    else
    {
      for (length_type x = 0; x != x_length; ++x)
	for (length_type y = 0; y != y_length; ++y)
	  rdft_1d.by_reference(in + x * in_x_stride + y * in_y_stride,
			       in_z_stride,
			       out + x * out_x_stride + y * out_y_stride,
			       out_z_stride, z_length);
      for (length_type z = 0; z != z_length/2 + 1; ++z)
	dft_2d.in_place(out + z * out_z_stride,
			out_x_stride, out_y_stride,
			x_length, y_length);
    }
  }
  virtual void by_reference(rtype *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    ztype out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<1, rtype, ctype, 0, -1> rdft_1d;
    dft<2, ctype, ctype, 0, -1> dft_2d;
    if (A == 0)
    {
      for (length_type y = 0; y != y_length; ++y)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(in + y * in_y_stride + z * in_z_stride,
			       in_x_stride,
			       offset(out, y * out_y_stride + z * out_z_stride),
			       out_x_stride, x_length);
      for (length_type x = 0; x != x_length/2 + 1; ++x)
	dft_2d.in_place(offset(out, x * out_x_stride),
			out_y_stride, out_z_stride,
			y_length, z_length);
    }
    else if (A == 1)
    {
      for (length_type x = 0; x != x_length; ++x)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(in + x * in_x_stride + z * in_z_stride,
			       in_y_stride,
			       offset(out, x * out_x_stride + z * out_z_stride),
			       out_y_stride, y_length);
      for (length_type y = 0; y != y_length/2 + 1; ++y)
	dft_2d.in_place(offset(out, y * out_y_stride),
			out_x_stride, out_z_stride,
			x_length, z_length);
    }
    else
    {
      for (length_type x = 0; x != x_length; ++x)
	for (length_type y = 0; y != y_length; ++y)
	  rdft_1d.by_reference(in + x * in_x_stride + y * in_y_stride,
			       in_z_stride,
			       offset(out, x * out_x_stride + y * out_y_stride),
			       out_z_stride, z_length);
      for (length_type z = 0; z != z_length/2 + 1; ++z)
	dft_2d.in_place(offset(out, z * out_z_stride),
			out_x_stride, out_y_stride,
			x_length, y_length);
    }
  }

};

// 3D complex -> real DFT
template <typename T, int A>
class dft<3, std::complex<T>, T, A, 1>
  : public fft::backend<3, std::complex<T>, T, A, 1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fft-dft-3D-real-inverse"; }
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(ctype *in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    rtype *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<2, ctype, ctype, 0, 1> dft_2d;
    dft<1, ctype, rtype, 0, 1> rdft_1d;
    if (A == 0)
    {
      length_type x2 = x_length/2 + 1;
      aligned_array<ctype> tmp(x2 * y_length * z_length);
      for (length_type x = 0; x != x2; ++x)
	dft_2d.by_reference(in + x * in_x_stride,
			    in_y_stride, in_z_stride,
			    tmp.get() + x * y_length * z_length,
			    1, y_length,
			    y_length, z_length);
      for (length_type y = 0; y != y_length; ++y)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(tmp.get() + y + z * y_length,
			       y_length * z_length,
			       out + y * out_y_stride + z * out_z_stride,
			       out_x_stride, x_length);
    }
    else if (A == 1)
    {
      length_type y2 = y_length/2 + 1;
      aligned_array<ctype> tmp(y2 * x_length * z_length);
      for (length_type y = 0; y != y2; ++y)
	dft_2d.by_reference(in + y * in_y_stride,
			    in_x_stride, in_z_stride,
			    tmp.get() + y * x_length * z_length,
			    1, x_length,
			    x_length, z_length);
      for (length_type x = 0; x != x_length; ++x)
	for (length_type z = 0; z != z_length; ++z)
	  rdft_1d.by_reference(tmp.get() + x + z * x_length,
			       x_length * z_length,
			       out + x * out_x_stride + z * out_z_stride,
			       out_y_stride, y_length);
    }
    else
    {
      length_type z2 = z_length/2 + 1;
      aligned_array<ctype> tmp(z2 * y_length * x_length);
      for (length_type z = 0; z != z2; ++z)
	dft_2d.by_reference(in + z * in_z_stride,
			    in_y_stride, in_x_stride,
			    tmp.get() + z * y_length * x_length,
			    1, y_length,
			    y_length, x_length);
      for (length_type y = 0; y != y_length; ++y)
	for (length_type x = 0; x != x_length; ++x)
	  rdft_1d.by_reference(tmp.get() + y + x * y_length,
			       y_length * x_length,
			       out + y * out_y_stride + x * out_x_stride,
			       out_z_stride, z_length);
    }
  }
  virtual void by_reference(ztype in,
			    stride_type in_x_stride,
			    stride_type in_y_stride,
			    stride_type in_z_stride,
			    rtype *out,
			    stride_type out_x_stride,
			    stride_type out_y_stride,
			    stride_type out_z_stride,
			    length_type x_length,
			    length_type y_length,
			    length_type z_length)
  {
    dft<2, ctype, ctype, 0, 1> dft_2d;
    dft<1, ctype, rtype, 0, 1> rdft_1d;
    if (A == 0)
    {
      length_type x2 = x_length/2 + 1;
      aligned_array<rtype> tmp_r(x2 * y_length * z_length);
      aligned_array<rtype> tmp_i(x2 * y_length * z_length);
      for (length_type x = 0; x != x2; ++x)
      {
	ztype line = std::make_pair(tmp_r.get() + x * y_length * z_length,
				    tmp_i.get() + x * y_length * z_length);
	dft_2d.by_reference(offset(in, x * in_x_stride),
			    in_y_stride, in_z_stride,
			    line,
			    1, y_length,
			    y_length, z_length);
      }
      for (length_type y = 0; y != y_length; ++y)
	for (length_type z = 0; z != z_length; ++z)
	{
	  ztype line = std::make_pair(tmp_r.get() + y + z * y_length,
				      tmp_i.get() + y + z * y_length);
	  rdft_1d.by_reference(line,
			       y_length * z_length,
			       out + y * out_y_stride + z * out_z_stride,
			       out_x_stride, x_length);
	}
    }
    else if (A == 1)
    {
      length_type y2 = y_length/2 + 1;
      aligned_array<rtype> tmp_r(y2 * x_length * z_length);
      aligned_array<rtype> tmp_i(y2 * x_length * z_length);
      for (length_type y = 0; y != y2; ++y)
      {
	ztype line = std::make_pair(tmp_r.get() + y * x_length * z_length,
				    tmp_i.get() + y * x_length * z_length);
	dft_2d.by_reference(offset(in, y * in_y_stride),
			    in_x_stride, in_z_stride,
			    line,
			    1, x_length,
			    x_length, z_length);
      }
      for (length_type x = 0; x != x_length; ++x)
	for (length_type z = 0; z != z_length; ++z)
	{
	  ztype line = std::make_pair(tmp_r.get() + x + z * x_length,
				      tmp_i.get() + x + z * x_length);
	  rdft_1d.by_reference(line,
			       x_length * z_length,
			       out + x * out_x_stride + z * out_z_stride,
			       out_y_stride, y_length);
	}
    }
    else
    {
      length_type z2 = z_length/2 + 1;
      aligned_array<rtype> tmp_r(z2 * y_length * x_length);
      aligned_array<rtype> tmp_i(z2 * y_length * x_length);
      for (length_type z = 0; z != z2; ++z)
      {
	ztype line = std::make_pair(tmp_r.get() + z * y_length * x_length,
				    tmp_i.get() + z * y_length * x_length);
	dft_2d.by_reference(offset(in, z * in_z_stride),
			    in_y_stride, in_x_stride,
			    line,
			    1, y_length,
			    y_length, x_length);
      }
      for (length_type y = 0; y != y_length; ++y)
	for (length_type x = 0; x != x_length; ++x)
	{
	  ztype line = std::make_pair(tmp_r.get() + y + x * y_length,
				      tmp_i.get() + y + x * y_length);
	  rdft_1d.by_reference(line,
			       y_length * x_length,
			       out + y * out_y_stride + x * out_x_stride,
			       out_z_stride, z_length);
	}
    }
  }

};

template <typename I, typename O, int A, int E> class dftm;

// real -> complex DFTM
template <typename T, int A>
class dftm<T, std::complex<T>, A, -1>
  : public fft::fftm<T, std::complex<T>, A, -1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fftm-dft-real-forward"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, rtype, ctype, 0, -1> rdft;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
	rdft.by_reference(in + c * in_c_stride, in_r_stride,
			  out + c * out_c_stride, out_r_stride, rows);
    else
      for (length_type r = 0; r != rows; ++r)
	rdft.by_reference(in + r * in_r_stride, in_c_stride,
			  out + r * out_r_stride, out_c_stride, cols);
  }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, rtype, ctype, 0, -1> rdft;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
	rdft.by_reference(in + c * in_c_stride, in_r_stride,
			  offset(out, c * out_c_stride), out_r_stride, rows);
    else
      for (length_type r = 0; r != rows; ++r)
	rdft.by_reference(in + r * in_r_stride, in_c_stride,
			  offset(out, r * out_r_stride), out_c_stride, cols);
  }
};

// complex -> real DFTM
template <typename T, int A>
class dftm<std::complex<T>, T, A, 1>
  : public fft::fftm<std::complex<T>, T, A, 1>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fftm-dft-real-inverse"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, rtype, 0, 1> rdft;
    if (A == 0)
    {
      for (length_type c = 0; c != cols; ++c)
	rdft.by_reference(in + c * in_c_stride, in_r_stride,
			  out + c * out_c_stride, out_r_stride, rows);
    }
    else
    {
      for (length_type r = 0; r != rows; ++r)
	rdft.by_reference(in + r * in_r_stride, in_c_stride,
			     out + r * out_r_stride, out_c_stride, cols);
    }
  }
  virtual void by_reference(ztype in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, rtype, 0, 1> rdft;
    if (A == 0)
    {
      for (length_type c = 0; c != cols; ++c)
      {
	ztype line = std::make_pair(in.first + c * in_c_stride,
				    in.second + c * in_c_stride);
	rdft.by_reference(line, in_r_stride,
			  out + c * out_c_stride, out_r_stride, rows);
      }
    }
    else
    {
      for (length_type r = 0; r != rows; ++r)
      {
	ztype line = std::make_pair(in.first + r * in_r_stride,
				    in.second + r * in_r_stride);
	rdft.by_reference(line, in_c_stride,
			  out + r * out_r_stride, out_c_stride, cols);
      }
    }
  }
};

// complex -> complex DFTM
template <typename T, int A, int E>
class dftm<std::complex<T>, std::complex<T>, A, E>
  : public fft::fftm<std::complex<T>, std::complex<T>, A, E>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual char const* name() { return "fftm-dft-complex"; }
  virtual void query_layout(Rt_layout<2> &) {}
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  { rtl_in.complex = rtl_out.complex; }
  virtual void in_place(ctype *inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
	dft_1d.in_place(inout + c * c_stride, r_stride, rows);
    else
      for (length_type r = 0; r != rows; ++r)
	dft_1d.in_place(inout + r * r_stride, c_stride, cols);
  }

  virtual void in_place(ztype inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
      {
	ztype line = std::make_pair(inout.first + c * c_stride,
				    inout.second + c * c_stride);
	dft_1d.in_place(line, r_stride, rows);
      }
    else
      for (length_type r = 0; r != rows; ++r)
      {
	ztype line = std::make_pair(inout.first + r * r_stride,
				    inout.second + r * r_stride);
	dft_1d.in_place(line, c_stride, cols);
      }
  }

  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
	dft_1d.by_reference(in + c * in_c_stride, in_r_stride,
			    out + c * out_c_stride, out_r_stride, rows);
    else
      for (length_type r = 0; r != rows; ++r)
	dft_1d.by_reference(in + r * in_r_stride, in_c_stride,
			    out + r * out_r_stride, out_c_stride, cols);
  }
  virtual void by_reference(ztype in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    dft<1, ctype, ctype, 0, E> dft_1d;
    if (A == 0)
      for (length_type c = 0; c != cols; ++c)
      {
	ztype in_line = std::make_pair(in.first + c * in_c_stride,
				       in.second + c * in_c_stride);
	ztype out_line = std::make_pair(out.first + c * out_c_stride,
					out.second + c * out_c_stride);
	dft_1d.by_reference(in_line, in_r_stride,
			    out_line, out_r_stride, rows);
      }
    else
      for (length_type r = 0; r != rows; ++r)
      {
	ztype in_line = std::make_pair(in.first + r * in_r_stride,
				       in.second + r * in_r_stride);
	ztype out_line = std::make_pair(out.first + r * out_r_stride,
					out.second + r * out_r_stride);
	dft_1d.by_reference(in_line, in_c_stride,
			    out_line, out_c_stride, cols);
      }
  }
};

struct DFT_tag;

template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<D, I, O, S, R, N, DFT_tag>
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<D> const &/*dom*/) { return true;}
  static std::auto_ptr<backend<D, I, O,
 			       axis<I, O, S>::value,
 			       exponent<I, O, S>::value> >
  create(Domain<D> const &/*dom*/, typename Scalar_of<I>::type /*scale*/)
  {
    static int const A = axis<I, O, S>::value;
    static int const E = exponent<I, O, S>::value;
    return std::auto_ptr<backend<D, I, O, A, E> >(new dft<D, I, O, A, E>());
  }
};

} // namespace vsip::impl::fft

namespace fftm
{
template <typename I,
	  typename O,
	  int A,
	  int E,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<I, O, A, E, R, N, fft::DFT_tag>
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<2> const &/*dom*/) { return true;}
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &/*dom*/, typename Scalar_of<I>::type /*scale*/)
  {
    return std::auto_ptr<fft::fftm<I, O, A, E> > (new fft::dftm<I, O, A, E>());
  }
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif

