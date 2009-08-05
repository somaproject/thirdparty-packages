/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fft/no_fft.hpp
    @author  Stefan Seefeld
    @date    2006-05-01
    @brief   VSIPL++ Library: FFT backend.
*/

#ifndef VSIP_CORE_FFT_NO_FFT_HPP
#define VSIP_CORE_FFT_NO_FFT_HPP

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
struct no_fft_base
{
  no_fft_base() 
  {
//     std::cout << "constructing no_fft_base" << std::endl;
  }
};


template <dimension_type D, typename I, typename O, int A, int E> class no_fft;

// 1D complex -> complex FFT
template <typename T, int A, int E>
class no_fft<1, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<1, std::complex<T>, std::complex<T>, A, E>,
    private no_fft_base

{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void in_place(ctype *, stride_type, length_type)
  {
  }
  virtual void in_place(ztype, stride_type, length_type)
  {
  }
  virtual void by_reference(ctype *, stride_type,
			    ctype *, stride_type,
			    length_type)
  {
  }
  virtual void by_reference(ztype, stride_type,
			    ztype, stride_type,
			    length_type)
  {
  }
};

// 1D real -> complex FFT
template <typename T, int A, int E>
class no_fft<1, T, std::complex<T>, A, E>
  : public fft::backend<1, T, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(rtype *, stride_type,
			    ctype *, stride_type,
			    length_type)
  {
  }
  virtual void by_reference(rtype *, stride_type,
			    ztype, stride_type,
			    length_type)
  {
  }
};

// 1D complex -> real FFT
template <typename T, int A, int E>
class no_fft<1, std::complex<T>, T, A, E>
  : public fft::backend<1, std::complex<T>, T, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(ctype *, stride_type,
			    rtype *, stride_type,
			    length_type)
  {
  }
  virtual void by_reference(ztype, stride_type,
			    rtype *, stride_type,
			    length_type)
  {
  }
};

// 2D complex -> complex FFT
template <typename T, int A, int E>
class no_fft<2, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<2, std::complex<T>, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void in_place(ctype *,
			stride_type, stride_type,
			length_type, length_type)
  {
  }
  virtual void in_place(ztype,
			stride_type, stride_type,
			length_type, length_type)
  {
  }
  virtual void by_reference(ctype *,
			    stride_type, stride_type,
			    ctype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
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
template <typename T, int A, int E>
class no_fft<2, T, std::complex<T>, A, E>
  : public fft::backend<2, T, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(rtype *,
			    stride_type, stride_type,
			    ctype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
  virtual void by_reference(rtype *,
			    stride_type, stride_type,
			    ztype,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }

};

// 2D complex -> real FFT
template <typename T, int A, int E>
class no_fft<2, std::complex<T>, T, A, E>
  : public fft::backend<2, std::complex<T>, T, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(ctype *,
			    stride_type, stride_type,
			    rtype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
  virtual void by_reference(ztype,
			    stride_type, stride_type,
			    rtype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }

};

// 3D complex -> complex FFT
template <typename T, int A, int E>
class no_fft<3, std::complex<T>, std::complex<T>, A, E>
  : public fft::backend<3, std::complex<T>, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void in_place(ctype *,
			stride_type,
			stride_type,
			stride_type,
			length_type,
			length_type,
			length_type)
  {
  }
  virtual void in_place(ztype,
			stride_type,
			stride_type,
			stride_type,
			length_type,
			length_type,
			length_type)
  {
  }
  virtual void by_reference(ctype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    ctype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }
  virtual void by_reference(ztype,
			    stride_type,
			    stride_type,
			    stride_type,
			    ztype,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }
};

// 3D real -> complex FFT
template <typename T, int A, int E>
class no_fft<3, T, std::complex<T>, A, E>
  : public fft::backend<3, T, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(rtype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    ctype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }
  virtual void by_reference(rtype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    ztype,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }

};

// 3D complex -> real FFT
template <typename T, int A, int E>
class no_fft<3, std::complex<T>, T, A, E>
  : public fft::backend<3, std::complex<T>, T, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(ctype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    rtype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }
  virtual void by_reference(ztype,
			    stride_type,
			    stride_type,
			    stride_type,
			    rtype *,
			    stride_type,
			    stride_type,
			    stride_type,
			    length_type,
			    length_type,
			    length_type)
  {
  }

};

template <typename I, typename O, int A, int E> class no_fftm;

// real -> complex FFTM
template <typename T, int A>
class no_fftm<T, std::complex<T>, A, -1>
  : public fft::fftm<T, std::complex<T>, A, -1>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(rtype *,
			    stride_type, stride_type,
			    ctype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
  virtual void by_reference(rtype *,
			    stride_type, stride_type,
			    ztype,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
};

// complex -> real FFTM
template <typename T, int A>
class no_fftm<std::complex<T>, T, A, 1>
  : public fft::fftm<std::complex<T>, T, A, 1>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void by_reference(ctype *,
			    stride_type, stride_type,
			    rtype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
  virtual void by_reference(ztype,
			    stride_type, stride_type,
			    rtype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
};

// complex -> complex FFTM
template <typename T, int A, int E>
class no_fftm<std::complex<T>, std::complex<T>, A, E>
  : public fft::fftm<std::complex<T>, std::complex<T>, A, E>,
    private no_fft_base
{
  typedef T rtype;
  typedef std::complex<rtype> ctype;
  typedef std::pair<rtype*, rtype*> ztype;

public:
  virtual void in_place(ctype *,
			stride_type, stride_type,
			length_type, length_type)
  {
  }

  virtual void in_place(ztype,
			stride_type, stride_type,
			length_type, length_type)
  {
  }

  virtual void by_reference(ctype *,
			    stride_type, stride_type,
			    ctype *,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
  virtual void by_reference(ztype,
			    stride_type, stride_type,
			    ztype,
			    stride_type, stride_type,
			    length_type, length_type)
  {
  }
};

struct No_FFT_tag;

template <dimension_type D,
	  typename I,
	  typename O,
	  int S,
	  vsip::return_mechanism_type R,
	  unsigned N>
struct evaluator<D, I, O, S, R, N, No_FFT_tag>
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<D> const &) { return true;}
  static std::auto_ptr<backend<D, I, O,
 			       axis<I, O, S>::value,
 			       exponent<I, O, S>::value> >
  create(Domain<D> const &/*dom*/, typename Scalar_of<I>::type /*scale*/)
  {
    static int const A = axis<I, O, S>::value;
    static int const E = exponent<I, O, S>::value;
    return std::auto_ptr<backend<D, I, O, A, E> >(new no_fft<D, I, O, A, E>());
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
struct evaluator<I, O, A, E, R, N, fft::No_FFT_tag>
{
  static bool const ct_valid = true;
  static bool rt_valid(Domain<2> const &/*dom*/) { return true;}
  static std::auto_ptr<fft::fftm<I, O, A, E> > 
  create(Domain<2> const &/*dom*/, typename Scalar_of<I>::type /*scale*/)
  {
    return std::auto_ptr<fft::fftm<I, O, A, E> >
      (new fft::no_fftm<I, O, A, E>());
  }
};

} // namespace vsip::impl::fftm
} // namespace vsip::impl
} // namespace vsip

#endif

