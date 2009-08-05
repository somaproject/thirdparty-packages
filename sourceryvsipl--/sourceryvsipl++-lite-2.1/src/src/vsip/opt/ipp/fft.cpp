/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ipp/fft.cpp
    @author  Stefan Seefeld, Nathan Myers
    @date    2006-05-05
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             Intel's IPP.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/support.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/fft/backend.hpp>
#include <vsip/core/fft/util.hpp>
#include <vsip/opt/ipp/fft.hpp>
#include <vsip/core/fns_scalar.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <ipps.h>
#include <ippi.h>
#include <complex>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace ipp
{
namespace
{
inline int
int_log2(unsigned size)    // assume size = 2^n, != 0, return n.
{
  int n = 0;
  while (size >>= 1) ++n;
  return n;
}
}

// Intel's type naming convention:
// suffix:  C type:         typedef:
//
// 32f      float           Ipp32f
// 32fc     complex float   Ipp32fc
// 64f      double          Ipp64f
// 64fc     complex double  Ipp64fc

// flag values
//  IPP_FFT_NODIV_BY_ANY, IPP_FFT_DIV_INV_BY_N, IPP_FFT_DIV_FWD_BY_N

template <typename T, //< Type used in the API.
	  typename I, //< IPP's corresponding type.
	  typename P, //< Plan type.
	  IppStatus (VSIP_IMPL_IPP_CALL *Plan)(P**, int, int, IppHintAlgorithm),
	  IppStatus (VSIP_IMPL_IPP_CALL *Dispose)(P*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Bufsize)(P const*, int*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Forward)(I const*, I*, P const*, Ipp8u*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Inverse)(I const*, I*, P const*, Ipp8u*)>
struct Driver_base_1d
{
  Driver_base_1d() : plan_(0), buffer_(0) {}
  void init(int x, int flags) VSIP_THROW((std::bad_alloc))
  {
    IppStatus result = (*Plan)(&plan_, x, flags, ippAlgHintFast);
    if (result != ippStsNoErr) VSIP_THROW(std::bad_alloc());
    buffer_ = alloc_align<Ipp8u>(VSIP_IMPL_ALLOC_ALIGNMENT, bufsize());
    if (!buffer_)
    {
      IppStatus result = (*Dispose)(plan_);
      assert(result == ippStsNoErr);
      VSIP_THROW(std::bad_alloc());
    }
  }
  void fini() VSIP_NOTHROW
  {
    free_align(buffer_);
    IppStatus result = (*Dispose)(plan_);
    assert(result == ippStsNoErr);
  }
  int bufsize() VSIP_NOTHROW
  {
    int size;
    IppStatus result = (*Bufsize)(plan_, &size);
    assert(result == ippStsNoErr);
    return size;
  }
  void forward(T const* in, T* out)
      VSIP_NOTHROW
  {
    IppStatus result = (*Forward)(reinterpret_cast<I const*>(in),
				  reinterpret_cast<I*>(out),
				  plan_, buffer_);
    assert(result == ippStsNoErr);
  }
  void inverse(T const* in, T* out)
    VSIP_NOTHROW
  {
    IppStatus result = (*Inverse)(reinterpret_cast<I const*>(in),
				  reinterpret_cast<I*>(out),
				  plan_, buffer_);
    assert(result == ippStsNoErr);
  }

  P *plan_;
  Ipp8u *buffer_;
};

template <typename T, //< Type used in the API.
	  typename I, //< IPP's corresponding type.
	  typename P, //< Plan type.
	  IppStatus (VSIP_IMPL_IPP_CALL *Plan)(P**, int, int, int, IppHintAlgorithm),
	  IppStatus (VSIP_IMPL_IPP_CALL *Dispose)(P*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Bufsize)(P const*, int*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Forward)(I const*, int, I*, int, P const*, Ipp8u*),
	  IppStatus (VSIP_IMPL_IPP_CALL *Inverse)(I const*, int, I*, int, P const*, Ipp8u*)>
struct Driver_base_2d
{
  Driver_base_2d() : plan_(0), buffer_(0) {}
  void init(int x, int y, int flags) VSIP_THROW((std::bad_alloc))
  {
    // Attention: IPP uses the opposite axis order, compared to VSIPL++:
    IppStatus result = (*Plan)(&plan_, y, x, flags, ippAlgHintFast);
    if (result != ippStsNoErr) VSIP_THROW(std::bad_alloc());
    buffer_ = alloc_align<Ipp8u>(VSIP_IMPL_ALLOC_ALIGNMENT, bufsize());
    if (!buffer_)
    {
      IppStatus result = (*Dispose)(plan_);
      assert(result == ippStsNoErr);
      VSIP_THROW(std::bad_alloc());
    }
  }
  void fini() VSIP_NOTHROW
  {
    free_align(buffer_);
    IppStatus result = (*Dispose)(plan_);
    assert(result == ippStsNoErr);
  }
  int bufsize() VSIP_NOTHROW
  {
    int size;
    IppStatus result = (*Bufsize)(plan_, &size);
    assert(result == ippStsNoErr);
    return size;
  }
  void forward(T const* in, unsigned in_stride,
	       T* out, unsigned out_stride)
    VSIP_NOTHROW
  {
    IppStatus result = (*Forward)(reinterpret_cast<I const*>(in), 
				  sizeof(I) * in_stride,
				  reinterpret_cast<I*>(out),
				  sizeof(I) * out_stride,
				  plan_, buffer_);
    assert(result == ippStsNoErr);
  }
  void inverse(T const* in, unsigned in_stride,
	       T* out, unsigned out_stride)
    VSIP_NOTHROW
  {
    IppStatus result = (*Inverse)(reinterpret_cast<I const*>(in),
				  sizeof(T) * in_stride,
				  reinterpret_cast<I*>(out),
				  sizeof(T) * out_stride,
				  plan_, buffer_);
    assert(result == ippStsNoErr);
  }

  P *plan_;
  Ipp8u *buffer_;
};

template <dimension_type D, //< Dimension
	  typename T,       //< Type
	  bool F>           //< Fast (Use FFT if true, DFT otherwise).
struct Driver;

// 1D, complex -> complex, float
template <bool F>
struct Driver<1, std::complex<float>, F>
  : ITE_Type<F,
	     As_type<Driver_base_1d<std::complex<float>,
				    Ipp32fc,
				    IppsFFTSpec_C_32fc,
				    ippsFFTInitAlloc_C_32fc,
				    ippsFFTFree_C_32fc,
				    ippsFFTGetBufSize_C_32fc,
				    ippsFFTFwd_CToC_32fc,
				    ippsFFTInv_CToC_32fc> >,
	     As_type<Driver_base_1d<std::complex<float>,
				    Ipp32fc,
				    IppsDFTSpec_C_32fc,
				    ippsDFTInitAlloc_C_32fc,
				    ippsDFTFree_C_32fc,
				    ippsDFTGetBufSize_C_32fc,
				    ippsDFTFwd_CToC_32fc,
				    ippsDFTInv_CToC_32fc> > >::type
{
  Driver(Domain<1> const &dom) 
  {
    int size = dom.size();
    // For FFTs we actually pass the 2's exponent of the size.
    if (F) size = int_log2(size);
    this->init(size, IPP_FFT_NODIV_BY_ANY);
  }
  ~Driver() { this->fini();}
};

// Provide function wrapper to adjust to desired uniform signature.
IppStatus 
VSIP_IMPL_IPP_CALL ippiDFTInitAlloc_C_32fc(IppiDFTSpec_C_32fc**plan, int x, int y,
                                           int flag, IppHintAlgorithm hint)
{
  IppiSize roi = {x, y};
  return ippiDFTInitAlloc_C_32fc(plan, roi, flag, hint);
}


// 2D, complex -> complex, float
template <bool F>
struct Driver<2, std::complex<float>, F>
  : ITE_Type<F,
	     As_type<Driver_base_2d<std::complex<float>,
				    Ipp32fc,
				    IppiFFTSpec_C_32fc,
				    ippiFFTInitAlloc_C_32fc,
				    ippiFFTFree_C_32fc,
				    ippiFFTGetBufSize_C_32fc,
				    ippiFFTFwd_CToC_32fc_C1R,
				    ippiFFTInv_CToC_32fc_C1R> >,
	     As_type<Driver_base_2d<std::complex<float>,
				    Ipp32fc,
				    IppiDFTSpec_C_32fc,
				    ippiDFTInitAlloc_C_32fc,
				    ippiDFTFree_C_32fc,
				    ippiDFTGetBufSize_C_32fc,
				    ippiDFTFwd_CToC_32fc_C1R,
				    ippiDFTInv_CToC_32fc_C1R> > >::type
{
  Driver(Domain<2> const &dom) 
  {
    int x = dom[0].size();
    int y = dom[1].size();
    // For FFTs we actually pass the 2's exponent of the size.
    if (F)
    {
      x = int_log2(x);
      y = int_log2(y);
    }
    this->init(x, y, IPP_FFT_NODIV_BY_ANY);
  }
  ~Driver() { this->fini();}
};

// 1D, complex -> complex, double
template <bool F>
struct Driver<1, std::complex<double>, F>
  : ITE_Type<F,
	     As_type<Driver_base_1d<std::complex<double>,
				    Ipp64fc,
				    IppsFFTSpec_C_64fc,
				    ippsFFTInitAlloc_C_64fc,
				    ippsFFTFree_C_64fc,
				    ippsFFTGetBufSize_C_64fc,
				    ippsFFTFwd_CToC_64fc,
				    ippsFFTInv_CToC_64fc> >,
	     As_type<Driver_base_1d<std::complex<double>,
				    Ipp64fc,
				    IppsDFTSpec_C_64fc,
				    ippsDFTInitAlloc_C_64fc,
				    ippsDFTFree_C_64fc,
				    ippsDFTGetBufSize_C_64fc,
				    ippsDFTFwd_CToC_64fc,
				    ippsDFTInv_CToC_64fc> > >::type
{
  Driver(Domain<1> const &dom) 
  {
    int size = dom.size();
    // For FFTs we actually pass the 2's exponent of the size.
    if (F) size = int_log2(size);
    this->init(size, IPP_FFT_NODIV_BY_ANY);
  }
  ~Driver() { this->fini();}
};

// 1D, complex -> real, float
// 1D, real -> complex, float
template <bool F>
struct Driver<1, float, F>
  : ITE_Type<F,
	     As_type<Driver_base_1d<float,
				    Ipp32f,
				    IppsFFTSpec_R_32f,
				    ippsFFTInitAlloc_R_32f,
				    ippsFFTFree_R_32f,
				    ippsFFTGetBufSize_R_32f,
				    ippsFFTFwd_RToCCS_32f,
				    ippsFFTInv_CCSToR_32f> >,
	     As_type<Driver_base_1d<float,
				    Ipp32f,
				    IppsDFTSpec_R_32f,
				    ippsDFTInitAlloc_R_32f,
				    ippsDFTFree_R_32f,
				    ippsDFTGetBufSize_R_32f,
				    ippsDFTFwd_RToCCS_32f,
				    ippsDFTInv_CCSToR_32f> > >::type
{
  Driver(Domain<1> const &dom) 
  {
    int size = dom.size();
    // For FFTs we actually pass the 2's exponent of the size.
    if (F) size = int_log2(size);
    this->init(size, IPP_FFT_NODIV_BY_ANY);
  }
  ~Driver() { this->fini();}
};

// 2D, complex -> real, float
// 2D, real -> complex, float
// 
// Not implemented yet. Right now we will typically fall back to another backend,
// such as fftw.

// 1D, complex -> real, double
// 1D, real -> complex, double
template <bool F>
struct Driver<1, double, F>
  : ITE_Type<F,
	     As_type<Driver_base_1d<double,
				    Ipp64f,
				    IppsFFTSpec_R_64f,
				    ippsFFTInitAlloc_R_64f,
				    ippsFFTFree_R_64f,
				    ippsFFTGetBufSize_R_64f,
				    ippsFFTFwd_RToCCS_64f,
				    ippsFFTInv_CCSToR_64f> >,
	     As_type<Driver_base_1d<double,
				    Ipp64f,
				    IppsDFTSpec_R_64f,
				    ippsDFTInitAlloc_R_64f,
				    ippsDFTFree_R_64f,
				    ippsDFTGetBufSize_R_64f,
				    ippsDFTFwd_RToCCS_64f,
				    ippsDFTInv_CCSToR_64f> > >::type
{
  Driver(Domain<1> const &dom) 
  {
    int size = dom.size();
    // For FFTs we actually pass the 2's exponent of the size.
    if (F) size = int_log2(size);
    this->init(size, IPP_FFT_NODIV_BY_ANY);
  }
  ~Driver() { this->fini();}
};

template <dimension_type D, //< Dimension
	  typename I,       //< Input type
	  typename O,       //< Output type
	  int A,            //< Axis
	  int E,            //< Exponent
	  bool F>           //< Fast (FFT, as opposed to DFT)
class impl;

// 1D complex -> complex FFT
template <typename T, int A, int E, bool F>
class impl<1, std::complex<T>, std::complex<T>, A, E, F>
  : public fft::backend<1, std::complex<T>, std::complex<T>, A, E>,
    private Driver<1, std::complex<T>, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<1> const &dom, rtype /*scale*/)
    : Driver<1, std::complex<T>, F>(dom)
  {
  }
  virtual void in_place(ctype *inout, stride_type s, length_type /*l*/)
  {
    assert(s == 1);
    if (E == -1) this->forward(inout, inout);
    else this->inverse(inout, inout);
  }
  virtual void in_place(ztype, stride_type, length_type)
  {
  }
  virtual void by_reference(ctype *in, stride_type in_s,
			    ctype *out, stride_type out_s,
			    length_type /*l*/)
  {
    assert(in_s == 1 && out_s == 1);
    if (E == -1) this->forward(in, out);
    else this->inverse(in, out);
  }
  virtual void by_reference(ztype, stride_type,
			    ztype, stride_type,
			    length_type)
  {
  }
};

// 1D real -> complex FFT
template <typename T, int A, int E, bool F>
class impl<1, T, std::complex<T>, A, E, F>
  : public fft::backend<1, T, std::complex<T>, A, E>,
    private Driver<1, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<1> const &dom, rtype /*scale*/)
    : Driver<1, T, F>(dom)
  {
  }
  virtual void by_reference(rtype *in, stride_type in_s,
			    ctype *out, stride_type out_s,
			    length_type /*l*/)
  {
    assert(in_s == 1 && out_s == 1);
    if (E == -1) this->forward(in, reinterpret_cast<rtype*>(out));
    else this->inverse(in, reinterpret_cast<rtype*>(out));
  }
  virtual void by_reference(rtype */*in*/, stride_type /*in_s*/,
			    ztype /*out*/, stride_type /*out_s*/,
			    length_type /*l*/)
  {
  }
};

// 1D complex -> real FFT
template <typename T, int A, int E, bool F>
class impl<1, std::complex<T>, T, A, E, F>
  : public fft::backend<1, std::complex<T>, T, A, E>,
    private Driver<1, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<1> const &dom, rtype /*scale*/)
    : Driver<1, T, F>(dom)
  {
  }
  virtual void by_reference(ctype *in, stride_type in_s,
			    rtype *out, stride_type out_s,
			    length_type /*l*/)
  {
    assert(in_s == 1 && out_s == 1);
    if (E == -1) this->forward(reinterpret_cast<rtype*>(in), out);
    else this->inverse(reinterpret_cast<rtype*>(in), out);
  }
  virtual void by_reference(ztype, stride_type,
			    rtype*, stride_type,
			    length_type)
  {
  }
};

// 2D complex -> complex FFT
template <typename T, int A, int E, bool F>
class impl<2, std::complex<T>, std::complex<T>, A, E, F>
  : public fft::backend<2, std::complex<T>, std::complex<T>, A, E>,
    private Driver<2, std::complex<T>, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<2> const &dom, rtype /*scale*/)
    : Driver<2, std::complex<T>, F>(dom)
  {
  }
  virtual void in_place(ctype *inout,
			stride_type r_stride, stride_type c_stride,
			length_type /*rows*/, length_type /*cols*/)
  {
    if (A == 0)
    {
      assert(c_stride == 1);
      if (E == -1) this->forward(inout, r_stride, inout, r_stride);
      else this->inverse(inout, r_stride, inout, r_stride);
    }
    else
    {
      assert(r_stride == 1);
      if (E == -1) this->forward(inout, c_stride, inout, c_stride);
      else this->inverse(inout, c_stride, inout, c_stride);
    }
  }
  virtual void in_place(ztype,
			stride_type, stride_type,
			length_type, length_type)
  {
  }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type /*rows*/, length_type /*cols*/)
  {
    if (A == 0)
    {
      assert(in_c_stride == 1 && out_c_stride);
      if (E == -1) this->forward(in, in_r_stride, out, out_r_stride);
      else this->inverse(in, in_r_stride, out, out_r_stride);
    }
    else
    {
      assert(in_r_stride == 1 && out_r_stride == 1);
      if (E == -1) this->forward(in, in_c_stride, out, out_c_stride);
      else this->inverse(in, in_c_stride, out, out_c_stride);
    }
  }
  virtual void by_reference(ztype,
			    stride_type, stride_type,
			    ztype,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
};

// 2D real -> complex FFT
template <typename T, int A, int E, bool F>
class impl<2, T, std::complex<T>, A, E, F>
  : public fft::backend<2, T, std::complex<T>, A, E>,
    private Driver<2, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<2> const &dom, rtype scale)
    : Driver<2, T, F>(dom)
  {
    VSIP_IMPL_THROW(unimplemented("IPP FFT backend does not implement 2D real->complex FFT"));
  }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // IPP assumes A is the final dimension.
    if (A == 0) rtl_in.order = tuple<0, 1, 2>();
    else rtl_in.order = tuple<1, 0, 2>();
    rtl_in.complex = cmplx_inter_fmt;
    rtl_out = rtl_in;
  }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    assert(0);
  }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
  }

};

// 2D complex -> real FFT
template <typename T, int A, int E, bool F>
class impl<2, std::complex<T>, T, A, E, F>
  : public fft::backend<2, std::complex<T>, T, A, E>,
    private Driver<2, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  impl(Domain<2> const &dom, rtype scale)
    : Driver<2, T, F>(dom)
  {
    VSIP_IMPL_THROW(unimplemented("IPP FFT backend does not implement 2D complex->real FFT"));
  }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // IPP assumes A is the final dimension.
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = cmplx_inter_fmt;
    rtl_out = rtl_in;
  }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    assert(0);
  }
  virtual void by_reference(ztype,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
  }

};

template <typename I, //< Input type
	  typename O, //< Output type
	  int A,      //< Axis
	  int E,      //< Exponent
	  bool F>     //< Fast (FFT as opposed to DFT)
class fftm;

// real -> complex FFTM
template <typename T, int A, bool F>
class fftm<T, std::complex<T>, A, -1, F>
  : public fft::fftm<T, std::complex<T>, A, -1>,
    private Driver<1, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  fftm(Domain<2> const &dom, rtype /*scalar*/)
    : Driver<1, T, F>(dom[A]),
      mult_(dom[1 - A].size())
  {
  }
  virtual void by_reference(rtype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    if (A == 0)
    {
      assert(cols == mult_);
      for (length_type m = 0; m != mult_;
	   ++m, in += in_c_stride, out += out_c_stride)
	this->forward(in, reinterpret_cast<rtype*>(out));
    }
    else
    {
      assert(rows == mult_);
      for (length_type m = 0; m != mult_;
	   ++m, in += in_r_stride, out += out_r_stride)
	this->forward(in, reinterpret_cast<rtype*>(out));
    }
  }
  virtual void by_reference(rtype*, stride_type, stride_type,
			    ztype, stride_type, stride_type,
			    length_type, length_type)
  {
  }

  length_type mult_;
};

// complex -> real FFTM
template <typename T, int A, bool F>
class fftm<std::complex<T>, T, A, 1, F>
  : public fft::fftm<std::complex<T>, T, A, 1>,
    private Driver<1, T, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  fftm(Domain<2> const &dom, rtype /*scalar*/)
    : Driver<1, T, F>(dom[A]),
      mult_(dom[1 - A].size())
  {
  }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    rtype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    if (A == 0)
    {
      assert(cols == mult_);
      for (length_type m = 0; m != mult_;
	   ++m, in += in_c_stride, out += out_c_stride)
	this->inverse(reinterpret_cast<rtype*>(in), out);
    }
    else
    {
      assert(rows == mult_);
      for (length_type m = 0; m != mult_;
	   ++m, in += in_r_stride, out += out_r_stride)
	this->inverse(reinterpret_cast<rtype*>(in), out);
    }
  }
  virtual void by_reference(ztype, stride_type, stride_type,
			    rtype *, stride_type, stride_type,
			    length_type, length_type)
  {
  }

  length_type mult_;
};

// complex -> complex FFTM
template <typename T, int A, int E, bool F>
class fftm<std::complex<T>, std::complex<T>, A, E, F>
  : public fft::fftm<std::complex<T>, std::complex<T>, A, E>,
    private Driver<1, std::complex<T>, F>
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  fftm(Domain<2> const &dom, rtype /*scale*/)
    : Driver<1, std::complex<T>, F>(dom[A]),
      mult_(dom[1 - A].size())
  {
  }

  virtual void in_place(ctype *inout,
			stride_type r_stride, stride_type c_stride,
			length_type rows, length_type cols)
  {
    if (A == 0)
    {
      assert(cols <= mult_);
      for (length_type m = 0; m != cols; ++m, inout += c_stride)
	if (E == -1) this->forward(inout, inout);
	else this->inverse(inout, inout);
    }
    else
    {
      assert(rows <= mult_);
      for (length_type m = 0; m != rows; ++m, inout += r_stride)
	if (E == -1) this->forward(inout, inout);
	else this->inverse(inout, inout);
    }
  }

  virtual void in_place(ztype, stride_type, stride_type,
			length_type, length_type)
  {
  }

  virtual void by_reference(ctype *in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type rows, length_type cols)
  {
    if (A == 0)
    {
      assert(cols <= mult_);
      for (length_type m = 0; m != cols;
	   ++m, in += in_c_stride, out += out_c_stride)
	if (E == -1) this->forward(in, out);
	else this->inverse(in, out);
    }
    else
    {
      assert(rows <= mult_);
      for (length_type m = 0; m != rows;
	   ++m, in += in_r_stride, out += out_r_stride)
	if (E == -1) this->forward(in, out);
	else this->inverse(in, out);
    }
  }
  virtual void by_reference(ztype, stride_type, stride_type,
			    ztype, stride_type, stride_type,
			    length_type, length_type)
  {
  }

  length_type mult_;
};

#define VSIPL_IMPL_PROVIDE(D, I, O, A, E)		\
template <>                                             \
std::auto_ptr<fft::backend<D, I, O, A, E> >		\
create(Domain<D> const &dom,                            \
       fft::backend<D, I, O, A, E>::scalar_type scale,  \
       bool fast)                                       \
{                                                       \
  if (fast)						\
    return std::auto_ptr<fft::backend<D, I, O, A, E> >	\
      (new impl<D, I, O, A, E, true>(dom, scale));	\
  else                                                  \
    return std::auto_ptr<fft::backend<D, I, O, A, E> >	\
      (new impl<D, I, O, A, E, false>(dom, scale));	\
}

VSIPL_IMPL_PROVIDE(1, float, std::complex<float>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<float>, float, 0, 1)
VSIPL_IMPL_PROVIDE(1, std::complex<float>, std::complex<float>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<float>, std::complex<float>, 0, 1)
VSIPL_IMPL_PROVIDE(1, double, std::complex<double>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<double>, double, 0, 1)
VSIPL_IMPL_PROVIDE(1, std::complex<double>, std::complex<double>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<double>, std::complex<double>, 0, 1)

// TODO:
// VSIPL_IMPL_PROVIDE(2, float, std::complex<float>, 0, -1)
// VSIPL_IMPL_PROVIDE(2, float, std::complex<float>, 1, -1)
// VSIPL_IMPL_PROVIDE(2, std::complex<float>, float, 0, 1)
// VSIPL_IMPL_PROVIDE(2, std::complex<float>, float, 1, 1)
VSIPL_IMPL_PROVIDE(2, std::complex<float>, std::complex<float>, 0, -1)
VSIPL_IMPL_PROVIDE(2, std::complex<float>, std::complex<float>, 1, -1)
VSIPL_IMPL_PROVIDE(2, std::complex<float>, std::complex<float>, 0, 1)
VSIPL_IMPL_PROVIDE(2, std::complex<float>, std::complex<float>, 1, 1)

// Not supported by IPP:
// VSIPL_IMPL_PROVIDE(2, double, std::complex<double>, 0, -1)
// VSIPL_IMPL_PROVIDE(2, double, std::complex<double>, 1, -1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, double, 0, 1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, double, 1, 1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, std::complex<double>, 0, -1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, std::complex<double>, 1, -1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, std::complex<double>, 0, 1)
// VSIPL_IMPL_PROVIDE(2, std::complex<double>, std::complex<double>, 1, 1)

#undef VSIPL_IMPL_PROVIDE

#define VSIPL_IMPL_PROVIDE(I, O, A, E)			\
template <>                                             \
std::auto_ptr<fft::fftm<I, O, A, E> >			\
create(Domain<2> const &dom,                            \
       vsip::impl::Scalar_of<I>::type scale,            \
       bool fast)					\
{                                                       \
  if (fast)                                             \
    return std::auto_ptr<fft::fftm<I, O, A, E> >	\
      (new fftm<I, O, A, E, true>(dom, scale));		\
  else                                                  \
    return std::auto_ptr<fft::fftm<I, O, A, E> >	\
      (new fftm<I, O, A, E, false>(dom, scale));       	\
}

VSIPL_IMPL_PROVIDE(float, std::complex<float>, 0, -1)
VSIPL_IMPL_PROVIDE(float, std::complex<float>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<float>, float, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<float>, float, 1, 1)
VSIPL_IMPL_PROVIDE(std::complex<float>, std::complex<float>, 0, -1)
VSIPL_IMPL_PROVIDE(std::complex<float>, std::complex<float>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<float>, std::complex<float>, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<float>, std::complex<float>, 1, 1)

VSIPL_IMPL_PROVIDE(double, std::complex<double>, 0, -1)
VSIPL_IMPL_PROVIDE(double, std::complex<double>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<double>, double, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<double>, double, 1, 1)
VSIPL_IMPL_PROVIDE(std::complex<double>, std::complex<double>, 0, -1)
VSIPL_IMPL_PROVIDE(std::complex<double>, std::complex<double>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<double>, std::complex<double>, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<double>, std::complex<double>, 1, 1)

#undef VSIPL_IMPL_PROVIDE

} // namespace vsip::impl::ipp
} // namespace vsip::impl
} // namespace vsip
