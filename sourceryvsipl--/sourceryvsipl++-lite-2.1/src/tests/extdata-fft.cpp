/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    tests/extdata_fft.cpp
    @author  Jules Bergmann
    @date    2005-04-05
    @brief   VSIPL++ Library: Test for extdata for a SP object.

    This file illustrates how data access may be used to implement
    a signal processing object.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>
#include <cassert>
#include <vsip/support.hpp>
#include <vsip/dense.hpp>
#include <vsip/vector.hpp>
#include <vsip/core/fast_block.hpp>

#include <vsip_csl/test.hpp>
#include <vsip_csl/plainblock.hpp>
#include "extdata-output.hpp"

using namespace std;
using namespace vsip;
using namespace vsip_csl;


/***********************************************************************
  Declarations
***********************************************************************/

// Dummy interleaved FFT function, passes input through to output.

template <typename T>
void
fft_unit_stride(
  T const* in,
  T*       out,
  unsigned size)
{
  for (index_type i=0; i<size; ++i)
    out[i] = in[i];
}



// Dummy split FFT function, passes input through to output.

template <typename T>
void
fft_unit_stride_split(
  T const* in_real,
  T const* in_imag,
  T*       out_real,
  T*       out_imag,
  unsigned size)
{
  for (index_type i=0; i<size; ++i)
  {
    out_real[i] = in_real[i];
    out_imag[i] = in_imag[i];
  }
}



struct FFTFor {};
struct FFTInv {};



/// Dummy FFT object that requires access to interleaved complex.

template <typename T,
	  typename Dir>
class Test_FFT_inter
{
  // Constructors.
public:
  Test_FFT_inter(length_type size)
    : size_  (size),
      buffer_(new T[2*size]),
      verbose_(false)
  {}

  ~Test_FFT_inter()
  { delete[] buffer_; }

  template <typename Block1,
	    typename Block2>
  void operator()(
    char const *            str,
    const_Vector<T, Block1> vin,
    Vector      <T, Block2> vout)
  {
    typedef vsip::impl::Layout<1, row1_type,
      vsip::impl::Stride_unit, vsip::impl::Cmplx_inter_fmt>
		LP;
    typedef vsip::impl::Ext_data<Block1, LP> layout1;
    typedef vsip::impl::Ext_data<Block2, LP> layout2;

    // PROFILE: Check sizes, check layout costs.

    // Important to check size.  If vectors are too large, our
    // buffers will overflow.
    test_assert(vin.size()  == size_);
    test_assert(vout.size() == size_);

    if (verbose_)
    {
      cout << "Test_FFT_inter: " << str << endl;
      cout << "Block_layout<Block1>:\n";
      print_layout<typename vsip::impl::Block_layout<Block1>::layout_type>(cout);
      cout << endl;

      cout << "LP:\n";
      print_layout<LP>(cout);
      cout << endl;

      cout << "  access_type<LP>(vin.block()) = "
	   <<    access_type<LP>(vin.block()) << endl;
      cout << "  access_type<LP>(vout.block()) = "
	   <<    access_type<LP>(vout.block()) << endl;
      
      cout << "  mem_required<LP>(vin.block()) = "
	   <<    vsip::impl::mem_required<LP>(vin.block()) << endl;
      cout << "  mem_required<LP>(vout.block()) = "
	   <<    vsip::impl::mem_required<LP>(vout.block()) << endl;
    }

    test_assert(vsip::impl::mem_required<LP>(vin.block())  <= sizeof(T)*size_);
    test_assert(vsip::impl::mem_required<LP>(vout.block()) <= sizeof(T)*size_);

    layout1 rin (vin.block(),  vsip::impl::SYNC_IN,  buffer_ + 0);
    layout2 rout(vout.block(), vsip::impl::SYNC_OUT, buffer_ + size_);

    test_assert(rin.stride(0) == 1);
    test_assert(rin.size(0) == size_);

    test_assert(rout.stride(0) == 1);
    test_assert(rout.size(0) == size_);

    fft_unit_stride(rin.data(), rout.data(), size_);
  }

private:
  length_type size_;
  T*	      buffer_;
  bool        verbose_;
};



/// Dummy FFT object that requires access to split complex.

template <typename T,
	  typename Dir>
class Test_FFT_split
{
  // Constructors.
public:
  Test_FFT_split(length_type size)
    : size_  (size),
      buffer_(new typename T::value_type[4*size]),
      verbose_(false)
  {}

  ~Test_FFT_split()
  { delete[] buffer_; }

  template <typename Block1,
	    typename Block2>
  void operator()(
    char const *            str,
    const_Vector<T, Block1> vin,
    Vector      <T, Block2> vout)
  {
    typedef vsip::impl::Layout<1, row1_type,
      vsip::impl::Stride_unit, vsip::impl::Cmplx_split_fmt>
		LP;
    typedef vsip::impl::Ext_data<Block1, LP> layout1;
    typedef vsip::impl::Ext_data<Block2, LP> layout2;

    // PROFILE: Check sizes, check layout costs.

    // Important to check size.  If vectors are too large, our
    // buffers will overflow.
    test_assert(vin.size()  == size_);
    test_assert(vout.size() == size_);

    if (verbose_)
    {
      cout << "Test_FFT_split: " << str << endl;
      cout << "Block_layout<Block1>:\n";
      print_layout<typename vsip::impl::Block_layout<Block1>::layout_type>(cout);
      cout << endl;
      
      cout << "LP:\n";
      print_layout<LP>(cout);
      cout << endl;
      
      cout << "  access_type<LP>(vin.block()) = "
	   <<    access_type<LP>(vin.block()) << endl;
      cout << "  access_type<LP>(vout.block()) = "
	   <<    access_type<LP>(vout.block()) << endl;
      
      cout << "  mem_required<LP>(vin.block()) = "
	   <<    vsip::impl::mem_required<LP>(vin.block()) << endl;
      cout << "  mem_required<LP>(vout.block()) = "
	   <<    vsip::impl::mem_required<LP>(vout.block()) << endl;
    }

    test_assert(vsip::impl::mem_required<LP>(vin.block())  <= sizeof(T)*size_);
    test_assert(vsip::impl::mem_required<LP>(vout.block()) <= sizeof(T)*size_);

    layout1 rin (vin.block(),  vsip::impl::SYNC_IN,  
		 make_pair(buffer_ + 0, buffer_ + size_));
    layout2 rout(vout.block(), vsip::impl::SYNC_OUT,
		 make_pair(buffer_ + 2*size_, buffer_ + 3*size_));

    test_assert(rin.stride(0) == 1);
    test_assert(rin.size(0) == size_);

    test_assert(rout.stride(0) == 1);
    test_assert(rout.size(0) == size_);

    fft_unit_stride_split(rin.data().first,  rin.data().second,
			  rout.data().first, rout.data().second,
			  size_);
  }

private:
  length_type             size_;
  typename T::value_type* buffer_;
  bool                    verbose_;
};



/***********************************************************************
  Definitions
***********************************************************************/



// Fill vector with sequence of values.

template <typename T,
	  typename Block>
void
fill_view(
  Vector<complex<T>, Block> view,
  int                       k,
  Index<1>	            offset,
  Domain<1>                 /* dom */)
{
  for (index_type i=0; i<view.size(0); ++i)
    view.put(i, complex<T>(T(k*(i + offset[0])+1),
			   T(k*(i + offset[0])+2)));
}



template <typename T,
	  typename Block>
void
fill_view(
  Vector<T, Block> view,
  int              k)
{
  fill_view(view, k, Index<1>(0), Domain<1>(view.size(0)));
}



// Test values in view against sequence.

template <typename T,
	  typename Block>
void
test_view(const_Vector<complex<T>, Block> vec, int k)
{
  for (index_type i=0; i<vec.size(0); ++i)
  {
    if (!equal(vec.get(i), complex<T>(T(k*i+1), T(k*i+2))))
    {
      cout << "ERROR: i        = " << i << endl
	   << "       Got      = " << vec.get(i) << endl
	   << "       expected = " << vec.get(i) << endl;
    }
    test_assert(equal(vec.get(i), complex<T>(T(k*i+1), T(k*i+2))));
  }
}



template <template <typename, typename> class FFT,
	  typename                            Block>
void
test_fft_1d(length_type size, int k)
{
  typedef typename Block::value_type T;
  Vector<T, Block> in (size);
  Vector<T, Block> out(size);

  FFT<T, FFTFor> fft(size);

  fill_view(in, k);
  fft("vector", in, out);
  test_view(out, k);

  fft("subvector", in(Domain<1>(size)), out(Domain<1>(size)));
}

int
main()
{
  test_fft_1d<Test_FFT_inter, vsip::impl::Fast_block<1, complex<float> > >(256, 3);
  test_fft_1d<Test_FFT_split, vsip::impl::Fast_block<1, complex<float> > >(256, 3);
  return 0;
}
