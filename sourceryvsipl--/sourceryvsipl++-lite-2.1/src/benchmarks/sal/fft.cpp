/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/fft_sal.cpp
    @author  Jules Bergmann
    @date    2006-08-02
    @brief   VSIPL++ Library: Benchmark for SAL FFT.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/signal.hpp>

#include "benchmarks.hpp"



/***********************************************************************
  Definitions
***********************************************************************/

using namespace vsip;

using vsip::impl::Ext_data;
using vsip::impl::Cmplx_inter_fmt;
using vsip::impl::Cmplx_split_fmt;
using vsip::impl::Stride_unit_dense;
using vsip::impl::SYNC_IN;
using vsip::impl::SYNC_OUT;


// Wrapper class for SAL FFTs.

template <typename T,
	  typename ComplexFmt>
struct sal_fft;



// SAL FFT for interleaved complex (INCOMPLETE).

template <>
struct sal_fft<float, Cmplx_inter_fmt>
{
  typedef COMPLEX ctype;
};



// SAL FFT for split complex.

template <>
struct sal_fft<complex<float>, Cmplx_split_fmt>
{
  typedef COMPLEX_SPLIT type;

  static type to_ptr(std::pair<float*, float*> const& ptr)
  {
    type ret = { ptr.first, ptr.second };
    return ret;
  }

  static void fftop(
    FFT_setup& setup,
    type in,
    type out,
    type tmp,
    int  size,
    int  dir)
  {
    long eflag = 0;  // no caching hints
    fft_zoptx(&setup, &in, 1, &out, 1, &tmp, size, dir, eflag);
  }

  static void fftip(
    FFT_setup& setup,
    type inout,
    type tmp,
    int  size,
    int  dir)
  {
    long eflag = 0;  // no caching hints
    fft_ziptx(&setup, &inout, 1, &tmp, size, dir, eflag);
  }

  static void scale(
    type        data,
    length_type size,
    float       s)
  {
    vsmulx(data.realp, 1, &s, data.realp, 1, size, 0);
    vsmulx(data.imagp, 1, &s, data.imagp, 1, size, 0);
  }

};



template <>
struct sal_fft<complex<float>, Cmplx_inter_fmt>
{
  typedef COMPLEX* type;

  static type to_ptr(std::complex<float>* ptr)
  {
    return (type)ptr;
  }

  static void fftop(
    FFT_setup& setup,
    type in,
    type out,
    type tmp,
    int  size,
    int  dir)
  {
    long eflag = 0;  // no caching hints
    fft_coptx(&setup, in, 2, out, 2, tmp, size, dir, eflag);
  }

  static void fftip(
    FFT_setup& setup,
    type inout,
    type tmp,
    int  size,
    int  dir)
  {
    long eflag = 0;  // no caching hints
    fft_ciptx(&setup, inout, 2, tmp, size, dir, eflag);
  }

  static void scale(
    type        data,
    length_type size,
    float       s)
  {
    float *d = reinterpret_cast<float*>(data);
    vsmulx(d, 1, &s, d, 1, 2 * size, 0);
  }

};



inline unsigned long
ilog2(length_type size)    // assume size = 2^n, != 0, return n.
{
  unsigned int n = 0;
  while (size >>= 1) ++n;
  return n;
}


int
fft_ops(length_type len)
{
  return int(5 * std::log((float)len) / std::log(2.f));
}


template <typename T,
	  typename ComplexFmt>
struct t_fft_op : Benchmark_base
{
  typedef typename impl::Scalar_of<T>::type scalar_type;

  char* what() { return "t_fft_op"; }
  int ops_per_point(length_type len)  { return fft_ops(len); }
  int riob_per_point(length_type) { return -1*(int)sizeof(T); }
  int wiob_per_point(length_type) { return -1*(int)sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    typedef sal_fft<T, ComplexFmt> traits;

    typedef impl::Layout<1, row1_type, Stride_unit_dense, ComplexFmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A  (size, T());
    Vector<T, block_type>   tmp(size, T());
    Vector<T, block_type>   Z  (size);

    unsigned long log2N = ilog2(size);

    FFT_setup     setup;
    unsigned long nbytes  = 0;
    long          options = 0;
    long          dir     = FFT_FORWARD;
    scalar_type   factor  = scalar_type(1) / size;

    fft_setup(log2N, options, &setup, &nbytes);

    A = T(1);
    
    vsip::impl::profile::Timer t1;

    {
      Ext_data<block_type> ext_A(A.block(), SYNC_IN);
      Ext_data<block_type> ext_tmp(tmp.block(), SYNC_IN);
      Ext_data<block_type> ext_Z(Z.block(), SYNC_OUT);

      typename traits::type A_ptr = traits::to_ptr(ext_A.data());
      typename traits::type tmp_ptr = traits::to_ptr(ext_tmp.data());
      typename traits::type Z_ptr = traits::to_ptr(ext_Z.data());
    
      if (!scale_)
      {
	t1.start();
	for (index_type l=0; l<loop; ++l)
	  traits::fftop(setup, A_ptr, Z_ptr, tmp_ptr, log2N, dir);
	t1.stop();
      }
      else
      {
	t1.start();
	for (index_type l=0; l<loop; ++l)
	{
	  traits::fftop(setup, A_ptr, Z_ptr, tmp_ptr, log2N, dir);
	  traits::scale(Z_ptr, size, factor); 
	}
	t1.stop();
      }
    }
    
    if (!equal(Z.get(0), T(scale_ ? 1 : size)))
    {
      std::cout << "t_fft_op: ERROR" << std::endl;
      std::cout << "  got     : " << Z.get(0) << std::endl;
      std::cout << "  expected: " << T(scale_ ? 1 : size) << std::endl;
      abort();
    }
    
    time = t1.delta();
  }

  t_fft_op(bool scale) : scale_(scale) {}

  // Member data
  bool scale_;
};



template <typename T,
	  typename ComplexFmt>
struct t_fft_ip : Benchmark_base
{
  typedef typename impl::Scalar_of<T>::type scalar_type;

  char* what() { return "t_fft_ip"; }
  int ops_per_point(length_type len)  { return fft_ops(len); }
  int riob_per_point(length_type) { return -1*(int)sizeof(T); }
  int wiob_per_point(length_type) { return -1*(int)sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    typedef sal_fft<T, ComplexFmt> traits;

    typedef impl::Layout<1, row1_type, Stride_unit_dense, ComplexFmt> LP;
    typedef impl::Fast_block<1, T, LP, Local_map> block_type;

    Vector<T, block_type>   A  (size, T());
    Vector<T, block_type>   tmp(size, T());

    unsigned long log2N = ilog2(size);

    FFT_setup     setup;
    unsigned long nbytes  = 0;
    long          options = 0;
    long          dir     = FFT_FORWARD;
    scalar_type   factor  = scalar_type(1) / size;

    fft_setup(log2N, options, &setup, &nbytes);

    A = T(0);

    vsip::impl::profile::Timer t1;

    {
      Ext_data<block_type> ext_A(A.block(), SYNC_IN);
      Ext_data<block_type> ext_tmp(tmp.block(), SYNC_IN);

      typename traits::type A_ptr = traits::to_ptr(ext_A.data());
      typename traits::type tmp_ptr = traits::to_ptr(ext_tmp.data());
    
      if (!scale_)
      {
	t1.start();
	for (index_type l=0; l<loop; ++l)
	  traits::fftip(setup, A_ptr, tmp_ptr, log2N, dir);
	t1.stop();
	// Check answer
	A = T(1);
	traits::fftip(setup, A_ptr, tmp_ptr, log2N, dir);
      }
      else
      {
	t1.start();
	for (index_type l=0; l<loop; ++l)
	{
	  traits::fftip(setup, A_ptr, tmp_ptr, log2N, dir);
	  traits::scale(A_ptr, size, factor); 
	}
	t1.stop();
	// Check answer
	A = T(1);
	traits::fftip(setup, A_ptr, tmp_ptr, log2N, dir);
	traits::scale(A_ptr, size, factor); 
      }
    }
    
    if (!equal(A.get(0), T(scale_ ? 1 : size)))
    {
      std::cout << "t_fft_ip: ERROR" << std::endl;
      std::cout << "  got     : " << A.get(0) << std::endl;
      std::cout << "  expected: " << T(scale_ ? 1 : size) << std::endl;
      abort();
    }
    
    time = t1.delta();
  }

  t_fft_ip(bool scale) : scale_(scale) {}

  // Member data
  bool scale_;
};



void
defaults(Loop1P& loop)
{
  loop.start_ = 4;
}



int
test(Loop1P& loop, int what)
{
  switch (what)
  {
  case  1: loop(t_fft_op<complex<float>, Cmplx_split_fmt>(false)); break;
  case  2: loop(t_fft_ip<complex<float>, Cmplx_split_fmt>(false)); break;
  case  5: loop(t_fft_op<complex<float>, Cmplx_split_fmt>(true)); break;
  case  6: loop(t_fft_ip<complex<float>, Cmplx_split_fmt>(true)); break;

  case 11: loop(t_fft_op<complex<float>, Cmplx_inter_fmt>(false)); break;
  case 12: loop(t_fft_ip<complex<float>, Cmplx_inter_fmt>(false)); break;
  case 15: loop(t_fft_op<complex<float>, Cmplx_inter_fmt>(true)); break;
  case 16: loop(t_fft_ip<complex<float>, Cmplx_inter_fmt>(true)); break;

  default: return 0;
  }
  return 1;
}
