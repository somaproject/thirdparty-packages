/* Copyright (c) 2006, 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/fftw3/fft_impl.cpp
    @author  Stefan Seefeld
    @date    2006-04-10
    @brief   VSIPL++ Library: FFT wrappers and traits to bridge with 
             FFTW3.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <cassert>

#include <vsip/support.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/fft/backend.hpp>
#include <vsip/core/fft/util.hpp>
#include <vsip/opt/fftw3/fft.hpp>
#include <vsip/core/equal.hpp>
#include <vsip/dense.hpp>
#include <fftw3.h>
#include <cstring>


// 070729: FFTW 3.1.2's split in-place complex-to-complex FFT is
// broken on PowerPC and x86.  The plan captures the gap between the
// real and imaginary components.
//
// 070730: FFTW 3.1.2 split out-of-place complex FFT is broken on
// PowerPC.
//
// 070730: FFTW 3.1.2's split real-complex and complex-real FFT
// also appear broken on x86.
//
// Brave souls may set
//   USE_FFTW_SPLIT 1
//   USE_BROKEN_FFTW_SPLIT 0
// to attempt a work-around: copy, then perform transform out-of-place.

// Control whether FFTW split-complex transforms are performed at all.
#define USE_FFTW_SPLIT 0

// Control whether a subset broken FFTW split-complex transforms are
// worked around.
#define USE_BROKEN_FFTW_SPLIT 0


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace fftw3
{

#if USE_FFTW_SPLIT
typedef vsip::impl::dense_complex_type fftw3_complex_type;
#else
typedef Cmplx_inter_fmt fftw3_complex_type;
#endif

template <dimension_type D>
struct Fft_base<D, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE> >
{
  Fft_base(Domain<D> const& dom, int exp, int flags, bool aligned = false)
    VSIP_THROW((std::bad_alloc))
      : in_buffer_(dom.size()),
	out_buffer_(dom.size()),
        aligned_(aligned)
  {
    if (!aligned) flags |= FFTW_UNALIGNED;
    // For multi-dimensional transforms, these plans assume both
    // input and output data is dense, row-major, interleave-complex
    // format.
    
    for(index_type i=0;i<D;i++) size_[i] = dom[i].size();
    plan_in_place_ =
      Create_plan<fftw3_complex_type>
        ::create<FFTW(plan), FFTW(iodim)>
        (in_buffer_.ptr(), in_buffer_.ptr(), exp, flags, dom);
    
    if (!plan_in_place_) VSIP_IMPL_THROW(std::bad_alloc());

    plan_by_reference_ = Create_plan<fftw3_complex_type>
      ::create<FFTW(plan), FFTW(iodim)>
      (in_buffer_.ptr(), out_buffer_.ptr(), exp, flags, dom);

    if (!plan_by_reference_)
    {
      FFTW(destroy_plan)(plan_in_place_);
      VSIP_IMPL_THROW(std::bad_alloc());
    }
  }
  ~Fft_base() VSIP_NOTHROW
  {
    if (plan_in_place_) FFTW(destroy_plan)(plan_in_place_);
    if (plan_by_reference_) FFTW(destroy_plan)(plan_by_reference_);
  }

  Cmplx_buffer<fftw3_complex_type, SCALAR_TYPE> in_buffer_;
  Cmplx_buffer<fftw3_complex_type, SCALAR_TYPE> out_buffer_;
  FFTW(plan) plan_in_place_;
  FFTW(plan) plan_by_reference_;
  int size_[D];
  bool aligned_;
};

template <vsip::dimension_type D>
struct Fft_base<D, SCALAR_TYPE, std::complex<SCALAR_TYPE> >
{
  Fft_base(Domain<D> const& dom, int A, int flags, bool aligned = false)
    VSIP_THROW((std::bad_alloc))
    : in_buffer_(32, dom.size()),
      out_buffer_(dom.size()),
      aligned_(aligned)
  { 
    if (!aligned) flags |= FFTW_UNALIGNED;
    for (vsip::dimension_type i = 0; i < D; ++i) size_[i] = dom[i].size();  
    // FFTW3 assumes A == D - 1.
    // See also query_layout().
    if (A != D - 1) std::swap(size_[A], size_[D - 1]);
    plan_by_reference_ = Create_plan<fftw3_complex_type>::
      create<FFTW(plan), FFTW(iodim)>
      (in_buffer_.get(), out_buffer_.ptr(), A, flags, dom);
    if (!plan_by_reference_) VSIP_IMPL_THROW(std::bad_alloc());
  }
  ~Fft_base() VSIP_NOTHROW
  {
    if (plan_by_reference_) FFTW(destroy_plan)(plan_by_reference_);
  }

  aligned_array<SCALAR_TYPE> in_buffer_;
  Cmplx_buffer<fftw3_complex_type, SCALAR_TYPE> out_buffer_;
  FFTW(plan) plan_by_reference_;
  int size_[D];
  bool aligned_;
};

template <vsip::dimension_type D>
struct Fft_base<D, std::complex<SCALAR_TYPE>, SCALAR_TYPE>
{
  Fft_base(Domain<D> const& dom, int A, int flags, bool aligned = false)
    VSIP_THROW((std::bad_alloc))
    : in_buffer_(dom.size()),
      out_buffer_(32, dom.size()),
      aligned_(aligned)
  {
    if (!aligned) flags |= FFTW_UNALIGNED;
    for (vsip::dimension_type i = 0; i < D; ++i) size_[i] = dom[i].size();
    // FFTW3 assumes A == D - 1.
    // See also query_layout().
    if (A != D - 1) std::swap(size_[A], size_[D - 1]);
    plan_by_reference_ = Create_plan<fftw3_complex_type>::
      create<FFTW(plan), FFTW(iodim)>
      (in_buffer_.ptr(), out_buffer_.get(), A, flags, dom);

    if (!plan_by_reference_) VSIP_IMPL_THROW(std::bad_alloc());
  }
  ~Fft_base() VSIP_NOTHROW
  {
    if (plan_by_reference_) FFTW(destroy_plan)(plan_by_reference_);
  }

  Cmplx_buffer<fftw3_complex_type, SCALAR_TYPE> in_buffer_;
  aligned_array<SCALAR_TYPE>              out_buffer_;
  FFTW(plan) plan_by_reference_;
  int size_[D];
  bool aligned_;
};

// 1D complex -> complex FFT

template <int A, int E>
class Fft_impl<1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE> >,
    public fft::backend<1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>,
			A, E>

{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<1> const &dom, unsigned number)
    : Fft_base<1, ctype, ctype>(dom, E, convert_NoT(number),
                                // Only require aligned arrays if FFTW3 can actually take 
                                // advantage from it by using SIMD kernels.
                                !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(ctype)) &&
                                !((sizeof(ctype) * dom.length()) % VSIP_IMPL_ALLOC_ALIGNMENT))
  {}
  virtual char const* name() { return "fft-fftw3-1D-complex"; }
  virtual void query_layout(Rt_layout<1> &rtl_inout)
  {
    rtl_inout.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_inout.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    rtl_inout.order = tuple<0, 1, 2>();
    // make default based on library
    rtl_inout.complex = Create_plan<fftw3_complex_type>::format;
  }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual void in_place(ctype *inout, stride_type s, length_type l)
  {
    assert(s == 1 && static_cast<int>(l) == this->size_[0]);
    FFTW(execute_dft)(plan_in_place_,
		      reinterpret_cast<FFTW(complex)*>(inout),
		      reinterpret_cast<FFTW(complex)*>(inout));
  }
  virtual void in_place(ztype inout, stride_type s, length_type l)
  {
    assert(s == 1 && static_cast<int>(l) == this->size_[0]);

#if USE_BROKEN_FFTW_SPLIT
    if (E == -1)
      FFTW(execute_split_dft)(plan_in_place_,
			      inout.first, inout.second,
			      inout.first, inout.second);
    else
      FFTW(execute_split_dft)(plan_in_place_,
			      inout.second, inout.first,
			      inout.second, inout.first);
#else
    typedef Storage<fftw3_complex_type, ctype> storage_type;
    rtype* real = storage_type::get_real_ptr(in_buffer_.ptr());
    rtype* imag = storage_type::get_imag_ptr(in_buffer_.ptr());
    memcpy(real, inout.first, l*sizeof(rtype));
    memcpy(imag, inout.second, l*sizeof(rtype));
    if (E == -1)
      FFTW(execute_split_dft)(plan_by_reference_,
			      real, imag,
			      inout.first, inout.second);
    else
      FFTW(execute_split_dft)(plan_by_reference_,
			      imag, real,
			      inout.second, inout.first);
#endif
  }
  virtual void by_reference(ctype *in, stride_type in_stride,
			    ctype *out, stride_type out_stride,
			    length_type length)
  {
    assert(in_stride == 1 && out_stride == 1 &&
	   static_cast<int>(length) == this->size_[0]);
    FFTW(execute_dft)(plan_by_reference_,
		      reinterpret_cast<FFTW(complex)*>(in),
		      reinterpret_cast<FFTW(complex)*>(out));
  }
  virtual void by_reference(ztype in, stride_type in_stride,
			    ztype out, stride_type out_stride,
			    length_type length)
  {
    assert(in_stride == 1 && out_stride == 1 &&
	   static_cast<int>(length) == this->size_[0]);

    if (E == -1)
      FFTW(execute_split_dft)(plan_by_reference_,
			      in.first,  in.second,
			      out.first, out.second);
    else
      FFTW(execute_split_dft)(plan_by_reference_,
			      in.second,  in.first,
			      out.second, out.first);
  }
};

// 1D real -> complex FFT

template <int A, int E>
class Fft_impl<1, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<1, SCALAR_TYPE, std::complex<SCALAR_TYPE> >,
    public fft::backend<1, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<1> const &dom, unsigned number)
    : Fft_base<1, rtype, ctype>(dom, A, convert_NoT(number),
                                // Only require aligned arrays if FFTW3 can actually take 
                                // advantage from it by using SIMD kernels.
                                !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(rtype)) &&
                                !((sizeof(rtype) * dom.length()) % VSIP_IMPL_ALLOC_ALIGNMENT))
  {}
  virtual char const* name() { return "fft-fftw3-1D-real-forward"; }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual void by_reference(rtype *in, stride_type,
			    ctype *out, stride_type,
			    length_type)
  {
    FFTW(execute_dft_r2c)(plan_by_reference_, 
			  in, reinterpret_cast<FFTW(complex)*>(out));
  }
  virtual void by_reference(rtype *in, stride_type is,
			    ztype out, stride_type os,
			    length_type length)
  {
    assert(is == 1);
    assert(os == 1);
#if USE_BROKEN_FFTW_SPLIT
    FFTW(execute_split_dft_r2c)(plan_by_reference_, 
			  in, out.first, out.second);
#else
    typedef Storage<fftw3_complex_type, ctype> storage_type;
    rtype* out_r = storage_type::get_real_ptr(out_buffer_.ptr());
    rtype* out_i = storage_type::get_imag_ptr(out_buffer_.ptr());
    FFTW(execute_split_dft_r2c)(plan_by_reference_, 
				in, out_r, out_i);
    memcpy(out.first,  out_r, (length/2+1)*sizeof(rtype));
    memcpy(out.second, out_i, (length/2+1)*sizeof(rtype));
#endif
  }
};

// 1D complex -> real FFT

template <int A, int E>
class Fft_impl<1, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
  : Fft_base<1, std::complex<SCALAR_TYPE>, SCALAR_TYPE>,
    public fft::backend<1, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<1> const &dom, unsigned number)
    : Fft_base<1, ctype, rtype>(dom, A, convert_NoT(number),
                                // Only require aligned arrays if FFTW3 can actually take 
                                // advantage from it by using SIMD kernels.
                                !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(rtype)) &&
                                !((sizeof(rtype) * dom.length()) % VSIP_IMPL_ALLOC_ALIGNMENT))
  {}

  virtual char const* name() { return "fft-fftw3-1D-real-inverse"; }
  virtual void query_layout(Rt_layout<1> &rtl_in, Rt_layout<1> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }

  virtual bool requires_copy(Rt_layout<1> &) { return true;}

  virtual void by_reference(ctype *in, stride_type,
			    rtype *out, stride_type,
			    length_type)
  {
    FFTW(execute_dft_c2r)(plan_by_reference_,
			  reinterpret_cast<FFTW(complex)*>(in), out);
  }
  virtual void by_reference(ztype in, stride_type is,
			    rtype *out, stride_type os,
			    length_type length)
  {
    assert(is == 1);
    assert(os == 1);
#if USE_BROKEN_FFTW_SPLIT
    FFTW(execute_split_dft_c2r)(plan_by_reference_,
			  in.first, in.second, out);
#else
    typedef Storage<fftw3_complex_type, ctype> storage_type;
    rtype* in_r = storage_type::get_real_ptr(in_buffer_.ptr());
    rtype* in_i = storage_type::get_imag_ptr(in_buffer_.ptr());
    memcpy(in_r, in.first, (length/2+1)*sizeof(rtype));
    memcpy(in_i, in.second, (length/2+1)*sizeof(rtype));
    FFTW(execute_split_dft_c2r)(plan_by_reference_,
			  in_r, in_i, out);
#endif
  }
};

// 2D complex -> complex FFT

template <int A, int E>
class Fft_impl<2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE> >,
    public fft::backend<2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>,
			A, E>

{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<2> const &dom, unsigned number)
    : Fft_base<2, ctype, ctype>(dom, E, convert_NoT(number))
  {}
  virtual char const* name() { return "fft-fftw3-2D-complex"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    rtl_in.order = row2_type();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual void in_place(ctype *inout,
			stride_type r_stride,
			stride_type c_stride,
			length_type /*rows*/, length_type cols)
  {
    // Check that data is dense row-major.
    assert(r_stride == static_cast<stride_type>(cols));
    assert(c_stride == 1);

    FFTW(execute_dft)(plan_in_place_,
		      reinterpret_cast<FFTW(complex)*>(inout),
		      reinterpret_cast<FFTW(complex)*>(inout));
  }
  /// complex (split) in-place
  virtual void in_place(ztype inout,
			stride_type, stride_type,
			length_type, length_type)
  {
    FFTW(execute_split_dft)(plan_in_place_,
		      inout.first, inout.second,
		      inout.first, inout.second);
  }
  virtual void by_reference(ctype *in,
			    stride_type in_r_stride,
			    stride_type in_c_stride,
			    ctype *out,
			    stride_type out_r_stride,
			    stride_type out_c_stride,
			    length_type /*rows*/, length_type cols)
  {
    // Check that data is dense row-major.
    assert(in_r_stride == static_cast<stride_type>(cols));
    assert(in_c_stride == 1);
    assert(out_r_stride == static_cast<stride_type>(cols));
    assert(out_c_stride == 1);

    FFTW(execute_dft)(plan_by_reference_,
		      reinterpret_cast<FFTW(complex)*>(in), 
		      reinterpret_cast<FFTW(complex)*>(out));
  }
  virtual void by_reference(ztype in,
			    stride_type in_r_stride, stride_type in_c_stride,
			    ztype out,
			    stride_type out_r_stride, stride_type out_c_stride,
			    length_type, length_type cols)
  {
    // Check that data is dense row-major.
    assert(in_r_stride == static_cast<stride_type>(cols));
    assert(in_c_stride == 1);
    assert(out_r_stride == static_cast<stride_type>(cols));
    assert(out_c_stride == 1);

    FFTW(execute_split_dft)(plan_by_reference_,
                            in.first, in.second,
                            out.first, out.second);
  }
};

// 2D real -> complex FFT

template <int A, int E>
class Fft_impl<2, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<2, SCALAR_TYPE, std::complex<SCALAR_TYPE> >,
    public fft::backend<2, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<2> const &dom, unsigned number)
    : Fft_base<2, rtype, ctype>(dom, A, convert_NoT(number))
  {}

  virtual char const* name() { return "fft-fftw3-2D-real-forward"; }

  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // FFTW3 assumes A is the last dimension.
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return true;}

  virtual void by_reference(rtype *in,
			    stride_type, stride_type,
			    ctype *out,
			    stride_type, stride_type,
			    length_type, length_type)
  {
    FFTW(execute_dft_r2c)(plan_by_reference_,
			  in, reinterpret_cast<FFTW(complex)*>(out));
  }
  virtual void by_reference(rtype *in,
			    stride_type, stride_type,
			    ztype out,
			    stride_type, stride_type,
			    length_type, length_type)
  {
    FFTW(execute_split_dft_r2c)(plan_by_reference_,
			  in, out.first, out.second);
  }

};

// 2D complex -> real FFT

template <int A, int E>
class Fft_impl<2, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
  : Fft_base<2, std::complex<SCALAR_TYPE>, SCALAR_TYPE>,
    public fft::backend<2, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<2> const &dom, unsigned number)
    : Fft_base<2, ctype, rtype>(dom, A, convert_NoT(number))
  {}

  virtual char const* name() { return "fft-fftw3-2D-real-inverse"; }

  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // FFTW3 assumes A is the last dimension.
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else rtl_in.order = tuple<0, 1, 2>();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return true;}

  virtual void by_reference(ctype *in,
			    stride_type, stride_type,
			    rtype *out,
			    stride_type, stride_type,
			    length_type, length_type)
  {
    FFTW(execute_dft_c2r)(plan_by_reference_, 
			  reinterpret_cast<FFTW(complex)*>(in), out);
  }
  virtual void by_reference(ztype in,
			    stride_type, stride_type,
			    rtype *out,
			    stride_type, stride_type,
			    length_type, length_type)
  {
    FFTW(execute_split_dft_c2r)(plan_by_reference_, 
			  in.first, in.second, out);
  }

};

// 3D complex -> complex FFT

template <int A, int E>
class Fft_impl<3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE> >,
    public fft::backend<3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>,
			A, E>

{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<3> const &dom, unsigned number)
    : Fft_base<3, ctype, ctype>(dom, E, convert_NoT(number))
  {}
  virtual char const* name() { return "fft-fftw3-3D-complex"; }
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    rtl_in.order = row3_type();
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual void in_place(ctype *inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length)
  {
    assert(static_cast<int>(x_length) == this->size_[0]);
    assert(static_cast<int>(y_length) == this->size_[1]);
    assert(static_cast<int>(z_length) == this->size_[2]);

    // Check that data is dense row-major.
    assert(x_stride == static_cast<stride_type>(y_length*z_length));
    assert(y_stride == static_cast<stride_type>(z_length));
    assert(z_stride == 1);

    FFTW(execute_dft)(plan_in_place_,
		      reinterpret_cast<FFTW(complex)*>(inout),
		      reinterpret_cast<FFTW(complex)*>(inout));
  }
  virtual void in_place(ztype inout,
			stride_type x_stride,
			stride_type y_stride,
			stride_type z_stride,
			length_type x_length,
			length_type y_length,
			length_type z_length)
  {
    assert(static_cast<int>(x_length) == this->size_[0]);
    assert(static_cast<int>(y_length) == this->size_[1]);
    assert(static_cast<int>(z_length) == this->size_[2]);

    // Check that data is dense row-major.
    assert(x_stride == static_cast<stride_type>(y_length*z_length));
    assert(y_stride == static_cast<stride_type>(z_length));
    assert(z_stride == 1);

    FFTW(execute_split_dft)(plan_in_place_,
		      inout.first, inout.second,
		      inout.first, inout.second);
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
    assert(static_cast<int>(x_length) == this->size_[0]);
    assert(static_cast<int>(y_length) == this->size_[1]);
    assert(static_cast<int>(z_length) == this->size_[2]);

    // Check that data is dense row-major.
    assert(in_x_stride == static_cast<stride_type>(y_length*z_length));
    assert(in_y_stride == static_cast<stride_type>(z_length));
    assert(in_z_stride == 1);
    assert(out_x_stride == static_cast<stride_type>(y_length*z_length));
    assert(out_y_stride == static_cast<stride_type>(z_length));
    assert(out_z_stride == 1);

    FFTW(execute_dft)(plan_by_reference_,
		      reinterpret_cast<FFTW(complex)*>(in), 
		      reinterpret_cast<FFTW(complex)*>(out));
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
    assert(static_cast<int>(x_length) == this->size_[0]);
    assert(static_cast<int>(y_length) == this->size_[1]);
    assert(static_cast<int>(z_length) == this->size_[2]);

    // Check that data is dense row-major.
    assert(in_x_stride == static_cast<stride_type>(y_length*z_length));
    assert(in_y_stride == static_cast<stride_type>(z_length));
    assert(in_z_stride == 1);
    assert(out_x_stride == static_cast<stride_type>(y_length*z_length));
    assert(out_y_stride == static_cast<stride_type>(z_length));
    assert(out_z_stride == 1);

    FFTW(execute_split_dft)(plan_by_reference_,
                      in.first, in.second,
                      out.first, out.second);
  }
};

// 3D real -> complex FFT

template <int A, int E>
class Fft_impl<3, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<3, SCALAR_TYPE, std::complex<SCALAR_TYPE> >,
    public fft::backend<3, SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<3> const &dom, unsigned number)
    : Fft_base<3, rtype, ctype>(dom, A, convert_NoT(number))
  {}

  virtual char const* name() { return "fft-fftw3-3D-real-forward"; }

  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // FFTW3 assumes A is the last dimension.
    switch (A)
    {
      case 0: rtl_in.order = tuple<2, 1, 0>(); break;
      case 1: rtl_in.order = tuple<0, 2, 1>(); break;
      default: rtl_in.order = tuple<0, 1, 2>(); break;
    }
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual bool requires_copy(Rt_layout<3> &) { return true;}

  virtual void by_reference(rtype *in,
			    stride_type,
			    stride_type,
			    stride_type,
			    ctype *out,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
    FFTW(execute_dft_r2c)(plan_by_reference_,
			  in, reinterpret_cast<FFTW(complex)*>(out));
  }
  virtual void by_reference(rtype *in,
			    stride_type,
			    stride_type,
			    stride_type,
			    ztype out,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
    FFTW(execute_split_dft_r2c)(plan_by_reference_,
			  in, out.first, out.second);
  }

};

// 3D complex -> real FFT

template <int A, int E>
class Fft_impl<3, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
  : Fft_base<3, std::complex<SCALAR_TYPE>, SCALAR_TYPE>,
    public fft::backend<3, std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fft_impl(Domain<3> const &dom, unsigned number)
    : Fft_base<3, ctype, rtype>(dom, A, convert_NoT(number))
  {}

  virtual char const* name() { return "fft-fftw3-3D-real-inverse"; }
  
  virtual void query_layout(Rt_layout<3> &rtl_in, Rt_layout<3> &rtl_out)
  {
    rtl_in.pack = stride_unit_dense;
    // FFTW3 assumes A is the last dimension.
    switch (A)
    {
      case 0: rtl_in.order = tuple<2, 1, 0>(); break;
      case 1: rtl_in.order = tuple<0, 2, 1>(); break;
      default: rtl_in.order = tuple<0, 1, 2>(); break;
    }
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual bool requires_copy(Rt_layout<3> &) { return true;}

  virtual void by_reference(ctype *in,
			    stride_type,
			    stride_type,
			    stride_type,
			    rtype *out,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
    FFTW(execute_dft_c2r)(plan_by_reference_,
			  reinterpret_cast<FFTW(complex)*>(in), out);
  }
  virtual void by_reference(ztype in,
			    stride_type,
			    stride_type,
			    stride_type,
			    rtype *out,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
    FFTW(execute_split_dft_c2r)(plan_by_reference_,
			  in.first, in.second, out);
  }

};

// real -> complex FFTM

template <int A>
class Fftm_impl<SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, -1>
  : private Fft_base<1, SCALAR_TYPE, std::complex<SCALAR_TYPE> >,
    public fft::fftm<SCALAR_TYPE, std::complex<SCALAR_TYPE>, A, -1>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fftm_impl(Domain<2> const &dom, unsigned number)
    : Fft_base<1, SCALAR_TYPE, std::complex<SCALAR_TYPE> >
      (dom[A], 0, convert_NoT(number),
       // Only require aligned arrays if FFTW3 can actually take 
       // advantage from it by using SIMD kernels.
       !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(rtype)) &&
       !((sizeof(rtype) * dom[A].length()) % VSIP_IMPL_ALLOC_ALIGNMENT)),
      mult_(dom[1-A].size())
  {
  }
  virtual char const* name() { return "fftm-fftw3-real-forward"; }
  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else  rtl_in.order = tuple<0, 1, 2>();
    // make default based on library
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual void by_reference(rtype *in,
			    stride_type i_str_0, stride_type i_str_1,
			    ctype *out,
			    stride_type o_str_0, stride_type o_str_1,
			    length_type rows, length_type cols)
  {
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    for (index_type i = 0; i < n_fft; ++i)
    {
      FFTW(execute_dft_r2c)(plan_by_reference_, 
			    in, reinterpret_cast<FFTW(complex)*>(out));
      in  += in_fft_stride;
      out += out_fft_stride;
    }
  }
  virtual void by_reference(
    rtype*      in,
    stride_type i_str_0,
    stride_type i_str_1,
    ztype       out,
    stride_type o_str_0,
    stride_type o_str_1,
    length_type rows,
    length_type cols)
  {
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    rtype* out_r = out.first;
    rtype* out_i = out.second;

#if !USE_BROKEN_FFTW_SPLIT
    typedef Storage<fftw3_complex_type, ctype> storage_type;
    rtype* tmp_out_r = storage_type::get_real_ptr(out_buffer_.ptr());
    rtype* tmp_out_i = storage_type::get_imag_ptr(out_buffer_.ptr());
#endif

    for (index_type i = 0; i < n_fft; ++i)
    {
#if USE_BROKEN_FFTW_SPLIT
      FFTW(execute_split_dft_r2c)(plan_by_reference_,
				  in, out_r, out_i);
#else
      FFTW(execute_split_dft_r2c)(plan_by_reference_,
				  in, tmp_out_r, tmp_out_i);
      memcpy(out_r, tmp_out_r, (size_[0]/2+1)*sizeof(rtype));
      memcpy(out_i, tmp_out_i, (size_[0]/2+1)*sizeof(rtype));
#endif
      in    += in_fft_stride;
      out_r += out_fft_stride;
      out_i += out_fft_stride;
    }
  }

private:
  length_type mult_;
};

// complex -> real FFTM

template <int A>
class Fftm_impl<std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, 1>
  : private Fft_base<1, std::complex<SCALAR_TYPE>, SCALAR_TYPE>,
    public fft::fftm<std::complex<SCALAR_TYPE>, SCALAR_TYPE, A, 1>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fftm_impl(Domain<2> const &dom, unsigned number)
    : Fft_base<1, std::complex<SCALAR_TYPE>, SCALAR_TYPE>
      (dom[A], 0, convert_NoT(number),
       // Only require aligned arrays if FFTW3 can actually take 
       // advantage from it by using SIMD kernels.
       !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(rtype)) &&
       !((sizeof(rtype) * dom[A].length()) % VSIP_IMPL_ALLOC_ALIGNMENT)),
      mult_(dom[1-A].size())
  {
  }

  virtual char const* name() { return "fftm-fftw3-real-inverse"; }

  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else  rtl_in.order = tuple<0, 1, 2>();
    // make default based on library
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }
  virtual bool requires_copy(Rt_layout<2> &) { return true;}

  virtual void by_reference(ctype *in,
			    stride_type i_str_0, stride_type i_str_1,
			    rtype *out,
			    stride_type o_str_0, stride_type o_str_1,
			    length_type rows, length_type cols)
  {
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    for (index_type i = 0; i < n_fft; ++i)
    {
      FFTW(execute_dft_c2r)(plan_by_reference_, 
			    reinterpret_cast<FFTW(complex)*>(in), out);
      in  += in_fft_stride;
      out += out_fft_stride;
    }
  }
  virtual void by_reference(
    ztype       in,
    stride_type i_str_0,
    stride_type i_str_1,
    rtype*      out,
    stride_type o_str_0,
    stride_type o_str_1,
    length_type rows,
    length_type cols)
  {
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    rtype* in_r = in.first;
    rtype* in_i = in.second;

#if !USE_BROKEN_FFTW_SPLIT
    typedef Storage<fftw3_complex_type, ctype> storage_type;
    rtype* tmp_in_r = storage_type::get_real_ptr(in_buffer_.ptr());
    rtype* tmp_in_i = storage_type::get_imag_ptr(in_buffer_.ptr());
#endif

    for (index_type i = 0; i < n_fft; ++i)
    {
#if USE_BROKEN_FFTW_SPLIT
      FFTW(execute_split_dft_c2r)(plan_by_reference_,
				  in_r, in_i, out);
#else
      memcpy(tmp_in_r, in_r, (size_[0]/2+1)*sizeof(rtype));
      memcpy(tmp_in_i, in_i, (size_[0]/2+1)*sizeof(rtype));
      FFTW(execute_split_dft_c2r)(plan_by_reference_,
				  tmp_in_r, tmp_in_i, out);
#endif
      in_r += in_fft_stride;
      in_i += in_fft_stride;
      out  += out_fft_stride;
    }
  }

private:
  length_type mult_;
};

// complex -> complex FFTM

template <int A, int E>
class Fftm_impl<std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, A, E>
  : private Fft_base<1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE> >,
    public fft::fftm<std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, A, E>
{
  typedef SCALAR_TYPE rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  Fftm_impl(Domain<2> const &dom, int number)
    : Fft_base<1, ctype, ctype>
      (dom[A], E, convert_NoT(number),
       // Only require aligned arrays if FFTW3 can actually take 
       // advantage from it by using SIMD kernels.
       !(VSIP_IMPL_ALLOC_ALIGNMENT % sizeof(ctype)) &&
       !((sizeof(ctype) * dom[A].length()) % VSIP_IMPL_ALLOC_ALIGNMENT)),
      mult_(dom[1-A].size())
  {
  }

  virtual char const* name() { return "fftm-fftw3-complex"; }

  virtual void query_layout(Rt_layout<2> &rtl_inout)
  {
    // By default use unit_stride,
    rtl_inout.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    // an ordering that gives unit strides on the axis perpendicular to A,
    if (A == 0) rtl_inout.order = tuple<1, 0, 2>();
    else rtl_inout.order = tuple<0, 1, 2>();
    // make default based on library
    rtl_inout.complex = Create_plan<fftw3_complex_type>::format;
    rtl_inout.align = VSIP_IMPL_ALLOC_ALIGNMENT;
  }

  virtual void query_layout(Rt_layout<2> &rtl_in, Rt_layout<2> &rtl_out)
  {
    rtl_in.pack = this->aligned_ ? stride_unit_align : stride_unit_dense;
    rtl_in.align = VSIP_IMPL_ALLOC_ALIGNMENT;
    if (A == 0) rtl_in.order = tuple<1, 0, 2>();
    else  rtl_in.order = tuple<0, 1, 2>();
    // make default based on library
    rtl_in.complex = Create_plan<fftw3_complex_type>::format;
    rtl_out = rtl_in;
  }

  virtual void in_place(ctype *inout,
			stride_type str_0, stride_type str_1,
			length_type rows, length_type cols)
  {
    assert((Type_equal<fftw3_complex_type, Cmplx_inter_fmt>::value));

    length_type const n_fft       = (A == 1) ? rows : cols;
    stride_type const fft_stride  = (A == 1) ? str_0 : str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);
    assert(((A == 1) ? str_1 : str_0) == 1);

    for (index_type i = 0; i != n_fft; ++i)
    {
      FFTW(execute_dft)(this->plan_in_place_, 
 			reinterpret_cast<FFTW(complex)*>(inout),
 			reinterpret_cast<FFTW(complex)*>(inout));
      inout += fft_stride;
    }
  }

  virtual void in_place(
    ztype       inout,
    stride_type str_0,
    stride_type str_1,
    length_type rows,
    length_type cols)
  {
    assert((Type_equal<fftw3_complex_type, Cmplx_split_fmt>::value));

    length_type const n_fft       = (A == 1) ? rows : cols;
    stride_type const fft_stride  = (A == 1) ? str_0 : str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);
    assert(((A == 1) ? str_1 : str_0) == 1);

    rtype* real_ptr = inout.first;
    rtype* imag_ptr = inout.second;

    for (index_type i = 0; i != n_fft; ++i)
    {
#if USE_BROKEN_FFTW_SPLIT
      if (E == -1)
	FFTW(execute_split_dft)(this->plan_in_place_,
				real_ptr, imag_ptr,
				real_ptr, imag_ptr);
      else
	FFTW(execute_split_dft)(this->plan_in_place_,
				imag_ptr, real_ptr,
				imag_ptr, real_ptr);
#else
      typedef Storage<fftw3_complex_type, ctype> storage_type;
      rtype* tmp_in_r = storage_type::get_real_ptr(in_buffer_.ptr());
      rtype* tmp_in_i = storage_type::get_imag_ptr(in_buffer_.ptr());
      memcpy(tmp_in_r, real_ptr, size_[0]*sizeof(rtype));
      memcpy(tmp_in_i, imag_ptr, size_[0]*sizeof(rtype));
      if (E == -1)
	FFTW(execute_split_dft)(plan_by_reference_,
				tmp_in_r, tmp_in_i,
				real_ptr, imag_ptr);
      else
	FFTW(execute_split_dft)(plan_by_reference_,
				tmp_in_i, tmp_in_r,
				imag_ptr, real_ptr);
#endif

      real_ptr += fft_stride;
      imag_ptr += fft_stride;
    }
  }

  virtual void by_reference(ctype *in,
			    stride_type i_str_0, stride_type i_str_1,
			    ctype *out,
			    stride_type o_str_0, stride_type o_str_1,
			    length_type rows, length_type cols)
  {
    // If the inputs to the Fftm are distributed, the number of FFTs may
    // be less than mult_.
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    for (index_type i = 0; i != n_fft; ++i)
    {
      FFTW(execute_dft)(plan_by_reference_, 
			reinterpret_cast<FFTW(complex)*>(in), 
			reinterpret_cast<FFTW(complex)*>(out));
      in  += in_fft_stride;
      out += out_fft_stride;
    }
  }

  virtual void by_reference(
    ztype       in,
    stride_type i_str_0,
    stride_type i_str_1,
    ztype       out,
    stride_type o_str_0,
    stride_type o_str_1,
    length_type rows,
    length_type cols)
  {
    // If the inputs to the Fftm are distributed, the number of FFTs may
    // be less than mult_.
    length_type const n_fft          = (A == 1) ? rows : cols;
    length_type const in_fft_stride  = (A == 1) ? i_str_0 : i_str_1;
    length_type const out_fft_stride = (A == 1) ? o_str_0 : o_str_1;

    if (A == 1) assert(rows <= mult_ && static_cast<int>(cols) == size_[0]);
    else        assert(cols <= mult_ && static_cast<int>(rows) == size_[0]);

    rtype* in_real  = in.first;
    rtype* in_imag  = in.second;
    rtype* out_real = out.first;
    rtype* out_imag = out.second;

    for (index_type i = 0; i != n_fft; ++i)
    {
      if (E == -1)
	FFTW(execute_split_dft)(plan_by_reference_, 
				in_real, in_imag, out_real, out_imag);
      else
	FFTW(execute_split_dft)(plan_by_reference_, 
				in_imag, in_real, out_imag, out_real);
      in_real  += in_fft_stride;
      in_imag  += in_fft_stride;
      out_real += out_fft_stride;
      out_imag += out_fft_stride;
    }
  }

private:
  length_type mult_;
};

#define VSIPL_IMPL_PROVIDE(D, I, O, A, E)	       \
template <>                                            \
std::auto_ptr<fft::backend<D, I, O, A, E> >	       \
create(Domain<D> const &dom, unsigned number)	       \
{                                                      \
  return std::auto_ptr<fft::backend<D, I, O, A, E> >   \
    (new Fft_impl<D, I, O, A, E>(dom, number));        \
}

VSIPL_IMPL_PROVIDE(1, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 0, 1)
VSIPL_IMPL_PROVIDE(1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(1, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, 1)
VSIPL_IMPL_PROVIDE(2, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(2, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 0, 1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 1, 1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, 1)
VSIPL_IMPL_PROVIDE(2, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, 1)
VSIPL_IMPL_PROVIDE(3, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(3, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(3, SCALAR_TYPE, std::complex<SCALAR_TYPE>, 2, -1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 0, 1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 1, 1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, SCALAR_TYPE, 2, 1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 2, -1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, 1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, 1)
VSIPL_IMPL_PROVIDE(3, std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 2, 1)

#undef VSIPL_IMPL_PROVIDE

#define VSIPL_IMPL_PROVIDE(I, O, A, E)		       \
template <>                                            \
std::auto_ptr<fft::fftm<I, O, A, E> >		       \
create(Domain<2> const &dom, unsigned number)	       \
{                                                      \
  return std::auto_ptr<fft::fftm<I, O, A, E> >	       \
    (new Fftm_impl<I, O, A, E>(dom, number));          \
}

VSIPL_IMPL_PROVIDE(SCALAR_TYPE, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(SCALAR_TYPE, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, SCALAR_TYPE, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, SCALAR_TYPE, 1, 1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, -1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, -1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 0, 1)
VSIPL_IMPL_PROVIDE(std::complex<SCALAR_TYPE>, std::complex<SCALAR_TYPE>, 1, 1)

#undef VSIPL_IMPL_PROVIDE

} // namespace vsip::impl::fftw3
} // namespace vsip::impl
} // namespace vsip
