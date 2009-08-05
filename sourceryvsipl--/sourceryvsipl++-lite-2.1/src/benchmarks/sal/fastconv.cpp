/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/fastconv_sal.cpp
    @author  Jules Bergmann, Don McCoy
    @date    2005-10-28
    @brief   VSIPL++ Library: Benchmark for Fast Convolution.

*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/signal.hpp>
#include <vsip/opt/profile.hpp>

#include <vsip_csl/test.hpp>
#include "loop.hpp"
#include "benchmarks.hpp"
#include "fastconv.hpp"

using namespace vsip;
using namespace vsip_csl;

using impl::Stride_unit_dense;
using impl::Cmplx_inter_fmt;
using impl::Cmplx_split_fmt;



/***********************************************************************
  Common definitions
***********************************************************************/

inline unsigned long
ilog2(length_type size)    // assume size = 2^n, != 0, return n.
{
  unsigned int n = 0;
  while (size >>= 1) ++n;
  return n;
}



template <typename ComplexFmt>
struct Impl1ip;		// in-place, phased fast-convolution (fscm)
template <typename ComplexFmt>
struct Impl2ip;		// in-place, interleaved fast-convolution (fcs loop)



/***********************************************************************
  Impl1ip: in-place, phased fast-convolution
***********************************************************************/

template <>
struct t_fastconv_base<complex<float>, Impl1ip<Cmplx_inter_fmt> >
  : fastconv_ops
{
  static length_type const num_args = 1;

  typedef complex<float> T;

  void fastconv(length_type npulse, length_type nrange,
		length_type loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP1;
    typedef impl::Fast_block<1, T, LP1, Local_map> block1_type;

    typedef impl::Layout<2, row2_type, Stride_unit_dense, Cmplx_inter_fmt> LP2;
    typedef impl::Fast_block<2, T, LP2, Local_map> block2_type;
    
    // Create the data cube.
    Matrix<T, block2_type> data(npulse, nrange);
    
    // Create the pulse replica and temporary buffer
    Vector<T, block1_type> replica(nrange);
    Vector<T, block1_type> tmp(nrange);

    // Initialize
    data    = T();
    replica = T();

    // setup direct access to data buffers
    impl::Ext_data<block1_type> ext_replica(replica.block(), impl::SYNC_IN);
    impl::Ext_data<block1_type> ext_tmp(tmp.block(), impl::SYNC_IN);

    // Create weights array
    unsigned long log2N = ilog2(nrange);
    long flag = FFT_FAST_CONVOLUTION;
    FFT_setup setup;
    unsigned long nbytes = 0;
    fft_setup( log2N, flag, &setup, &nbytes );

    // Create filter
    COMPLEX* filter = reinterpret_cast<COMPLEX *>(ext_replica.data());
    COMPLEX* t = reinterpret_cast<COMPLEX *>(ext_tmp.data());
    float scale = 1;
    long f_conv = FFT_INVERSE;  // does not indicate direction, but rather 
                                // convolution as opposed to correlation
    long eflag = 0;  // no caching hints

    fcf_ciptx( &setup, filter, t, &scale, log2N, f_conv, eflag );

    
    // Set up convolution 
    impl::Ext_data<block2_type> ext_data(data.block(), impl::SYNC_IN);
    COMPLEX* msignal = reinterpret_cast<COMPLEX *>(ext_data.data());
    long jr = ext_data.stride(1);
    long jc = ext_data.stride(0);
    unsigned long M = npulse;
    
    vsip::impl::profile::Timer t1;
    
    // Impl1 ip
    t1.start();
    for (index_type l=0; l<loop; ++l)
    {
      // Perform fast convolution
      fcsm_ciptx( &setup, filter, msignal, jr, jc, t, 
		  log2N, M, f_conv, eflag );
    }
    t1.stop();

    // CHECK RESULT
    time = t1.delta();
    fft_free(&setup);
  }
};



/***********************************************************************
  Impl1ip: SPLIT in-place, phased fast-convolution
***********************************************************************/

template <>
struct t_fastconv_base<complex<float>, Impl1ip<Cmplx_split_fmt> >
  : fastconv_ops
{
  static length_type const num_args = 1;

  typedef complex<float> T;

  void fastconv(length_type npulse, length_type nrange,
		length_type loop, float& time)
  {
    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP1;
    typedef impl::Fast_block<1, T, LP1, Local_map> block1_type;

    typedef impl::Layout<2, row2_type, Stride_unit_dense, Cmplx_split_fmt> LP2;
    typedef impl::Fast_block<2, T, LP2, Local_map> block2_type;

    typedef impl::Ext_data<block1_type>::raw_ptr_type ptr_type;
    
    // Create the data cube.
    Matrix<T, block2_type> data(npulse, nrange);
    
    // Create the pulse replica and temporary buffer
    Vector<T, block1_type> replica(nrange);
    Vector<T, block1_type> tmp(nrange);

    // Initialize
    data    = T();
    replica = T();

    // setup direct access to data buffers
    impl::Ext_data<block1_type> ext_replica(replica.block(), impl::SYNC_IN);
    impl::Ext_data<block1_type> ext_tmp(tmp.block(), impl::SYNC_IN);

    // Create weights array
    unsigned long log2N = ilog2(nrange);
    long flag = FFT_FAST_CONVOLUTION;
    FFT_setup setup;
    unsigned long nbytes = 0;
    fft_setup( log2N, flag, &setup, &nbytes );

    // Create filter
    // COMPLEX* filter = reinterpret_cast<COMPLEX *>(ext_replica.data());
    // COMPLEX* t = reinterpret_cast<COMPLEX *>(ext_tmp.data());
    ptr_type p_replica = ext_replica.data();
    ptr_type p_tmp     = ext_tmp.data();

    COMPLEX_SPLIT filter;
    COMPLEX_SPLIT tmpbuf;

    filter.realp = p_replica.first;
    filter.imagp = p_replica.second;
    tmpbuf.realp = p_tmp.first;
    tmpbuf.imagp = p_tmp.second;

    float scale = 1;
    long f_conv = FFT_INVERSE;  // does not indicate direction, but rather 
                                // convolution as opposed to correlation
    long eflag = 0;  // no caching hints

    fcf_ziptx( &setup, &filter, &tmpbuf, &scale, log2N, f_conv, eflag );

    
    // Set up convolution 
    impl::Ext_data<block2_type> ext_data(data.block(), impl::SYNC_IN);
    // COMPLEX* msignal = reinterpret_cast<COMPLEX *>(ext_data.data());

    ptr_type p_data = ext_data.data();
    COMPLEX_SPLIT msignal;
    msignal.realp = p_data.first;
    msignal.imagp = p_data.second;

    long jr = ext_data.stride(1);
    long jc = ext_data.stride(0);
    unsigned long M = npulse;
    
    vsip::impl::profile::Timer t1;
    
    // Impl1 ip
    t1.start();
    for (index_type l=0; l<loop; ++l)
    {
      // Perform fast convolution
      fcsm_ziptx( &setup, &filter, &msignal, jr, jc, &tmpbuf, 
		  log2N, M, f_conv, eflag );
    }
    t1.stop();

    // CHECK RESULT
    time = t1.delta();
    fft_free(&setup);
  }
};





/***********************************************************************
  Impl2ip: out-of-place (tmp), interleaved fast-convolution
***********************************************************************/

template <>
struct t_fastconv_base<complex<float>, Impl2ip<Cmplx_inter_fmt> >
  : fastconv_ops
{
  static length_type const num_args = 1;

  typedef complex<float> T;

  void fastconv(length_type npulse, length_type nrange,
		length_type loop, float& time)
  {
    typedef impl::Layout<2, row2_type, Stride_unit_dense, Cmplx_inter_fmt> LP2;
    typedef impl::Fast_block<2, T, LP2, Local_map> block2_type;
    typedef Matrix<T, block2_type>            view_type;

    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_inter_fmt> LP1;
    typedef impl::Fast_block<1, T, LP1, Local_map> block1_type;
    typedef Vector<T, block1_type>          replica_view_type;

    // Create the data cube.
    view_type data(npulse, nrange);
    
    // Create the pulse replica
    replica_view_type replica(nrange);
    replica_view_type tmp(nrange);

    // Initialize
    data    = T();
    replica = T();

    // setup direct access to data buffers
    impl::Ext_data<block1_type> ext_replica(replica.block(), impl::SYNC_IN);
    impl::Ext_data<block1_type> ext_tmp(tmp.block(), impl::SYNC_IN);

    // Create weights array
    unsigned long log2N = ilog2(nrange);
    long flag = FFT_FAST_CONVOLUTION;
    FFT_setup setup = 0;
    unsigned long nbytes = 0;
    fft_setup( log2N, flag, &setup, &nbytes );

    // Create filter
    COMPLEX* filter = reinterpret_cast<COMPLEX *>(ext_replica.data());
    COMPLEX* t     = reinterpret_cast<COMPLEX *>(ext_tmp.data());
    float scale = 1;
    long f_conv = FFT_INVERSE;  // does not indicate direction, but rather 
                                // convolution as opposed to correlation
    long eflag = 0;  // no caching hints

    fcf_ciptx( &setup, filter, t, &scale, log2N, f_conv, eflag );

    // Set up convolution 
    impl::Ext_data<block2_type> ext_data(data.block(), impl::SYNC_IN);
    COMPLEX* signal    = reinterpret_cast<COMPLEX *>(ext_data.data());
    long signal_stride = 2*ext_data.stride(1);

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
    {
      for (index_type p=0; p<npulse; ++p)
      {
	// Perform fast convolution
	fcs_ciptx( &setup, filter, signal, signal_stride, t, log2N, f_conv, eflag );
      }
    }
    t1.stop();

    // CHECK RESULT
    time = t1.delta();
  }
};



/***********************************************************************
  Impl2ip: SPLIT out-of-place (tmp), interleaved fast-convolution
***********************************************************************/

template <>
struct t_fastconv_base<complex<float>, Impl2ip<Cmplx_split_fmt> >
  : fastconv_ops
{
  static length_type const num_args = 1;

  typedef complex<float> T;

  void fastconv(length_type npulse, length_type nrange,
		length_type loop, float& time)
  {
    typedef impl::Layout<2, row2_type, Stride_unit_dense, Cmplx_split_fmt> LP2;
    typedef impl::Fast_block<2, T, LP2, Local_map> block2_type;
    typedef Matrix<T, block2_type>            view_type;

    typedef impl::Layout<1, row1_type, Stride_unit_dense, Cmplx_split_fmt> LP1;
    typedef impl::Fast_block<1, T, LP1, Local_map> block1_type;
    typedef Vector<T, block1_type>          replica_view_type;

    typedef impl::Ext_data<block1_type>::raw_ptr_type ptr_type;

    // Create the data cube.
    view_type data(npulse, nrange);
    
    // Create the pulse replica
    replica_view_type replica(nrange);
    replica_view_type tmp(nrange);

    // Initialize
    data    = T();
    replica = T();

    // setup direct access to data buffers
    impl::Ext_data<block1_type> ext_replica(replica.block(), impl::SYNC_IN);
    impl::Ext_data<block1_type> ext_tmp(tmp.block(), impl::SYNC_IN);

    // Create weights array
    unsigned long log2N = ilog2(nrange);
    long flag = FFT_FAST_CONVOLUTION;
    FFT_setup setup = 0;
    unsigned long nbytes = 0;
    fft_setup( log2N, flag, &setup, &nbytes );

    // Create filter
    ptr_type p_replica = ext_replica.data();
    ptr_type p_tmp     = ext_tmp.data();

    COMPLEX_SPLIT filter;
    COMPLEX_SPLIT tmpbuf;

    filter.realp = p_replica.first;
    filter.imagp = p_replica.second;
    tmpbuf.realp = p_tmp.first;
    tmpbuf.imagp = p_tmp.second;


    float scale = 1;
    long f_conv = FFT_INVERSE;  // does not indicate direction, but rather 
                                // convolution as opposed to correlation
    long eflag = 0;  // no caching hints

    fcf_ziptx(&setup, &filter, &tmpbuf, &scale, log2N, f_conv, eflag );

    // Set up convolution 
    impl::Ext_data<block2_type> ext_data(data.block(), impl::SYNC_IN);

    ptr_type p_data = ext_data.data();
    COMPLEX_SPLIT signal;
    signal.realp = p_data.first;
    signal.imagp = p_data.second;

    long signal_stride = 1*ext_data.stride(1);

    vsip::impl::profile::Timer t1;
    
    t1.start();
    for (index_type l=0; l<loop; ++l)
    {
      for (index_type p=0; p<npulse; ++p)
      {
	// Perform fast convolution
	fcs_ziptx(&setup, &filter, &signal, signal_stride, &tmpbuf,
		  log2N, f_conv, eflag );
      }
    }
    t1.stop();

    // CHECK RESULT
    time = t1.delta();
  }
};



void
defaults(Loop1P& loop)
{
  loop.cal_        = 4;
  loop.start_      = 4;
  loop.stop_       = 16;
  loop.loop_start_ = 10;
  loop.user_param_ = 64;
}



int
test(Loop1P& loop, int what)
{
  typedef complex<float> C;

  typedef Cmplx_inter_fmt Cif;
  typedef Cmplx_split_fmt Csf;

  length_type param1 = loop.user_param_;
  switch (what)
  {
  case   2: loop(t_fastconv_pf<C, Impl1ip<Cif> >(param1)); break;
  case   6: loop(t_fastconv_pf<C, Impl2ip<Cif> >(param1)); break;

  case  12: loop(t_fastconv_rf<C, Impl1ip<Cif> >(param1)); break;
  case  16: loop(t_fastconv_rf<C, Impl2ip<Cif> >(param1)); break;

  case 102: loop(t_fastconv_pf<C, Impl1ip<Csf> >(param1)); break;
  case 106: loop(t_fastconv_pf<C, Impl2ip<Csf> >(param1)); break;
   
  case 112: loop(t_fastconv_rf<C, Impl1ip<Csf> >(param1)); break;
  case 116: loop(t_fastconv_rf<C, Impl2ip<Csf> >(param1)); break;

  default: return 0;
  }
  return 1;
}
