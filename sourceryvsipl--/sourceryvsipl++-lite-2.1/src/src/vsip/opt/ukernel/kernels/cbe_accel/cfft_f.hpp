/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/cfft_f.hpp
    @author  Jules Bergmann
    @date    2008-06-12
    @brief   VSIPL++ Library: UKernel to compute interleaved-complex float FFT's.
*/

#ifndef VSIP_OPT_UKERNEL_KERNELS_CBE_ACCEL_CFFT_F_HPP
#define VSIP_OPT_UKERNEL_KERNELS_CBE_ACCEL_CFFT_F_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <utility>
#include <complex>
#include <cassert>
#include <spu_intrinsics.h>
#include <cml.h>
#include <cml_core.h>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/kernels/params/fft_param.hpp>

#define MIN_FFT_1D_SIZE	  32
#define MAX_FFT_1D_SIZE	  4096

#define FFT_BUF1_SIZE_BYTES (2*MAX_FFT_1D_SIZE*sizeof(float))
#define FFT_BUF2_SIZE_BYTES (1*MAX_FFT_1D_SIZE*sizeof(float)+128)



/***********************************************************************
  Definitions
***********************************************************************/

struct Fft_kernel : Spu_kernel
{
  typedef std::complex<float>* in0_type;
  typedef std::complex<float>* out0_type;

  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;
  typedef Uk_fft_params param_type;

  Fft_kernel()
    : size(0)
  {}

  void init(param_type& params)
  {
#if DEBUG
    printf("uk_ccfft_f(%d): init:\n", params.size);
#endif

    size  = params.size;
    dir   = params.dir;
    scale = params.scale;

    int rt = cml_fft1d_setup_f(&fft, CML_FFT_CC, size, buf2);
    assert(rt && fft != NULL);
  }

  void compute(
    in0_type  const& in,
    out0_type const& out,
    Pinfo const&     p_in,
    Pinfo const&     p_out)
  {
#if DEBUG
    printf("uk_ccfft_f(%d): compute -- start %05x %05x\n", size, in, out);
#endif
    // Handle inverse FFT explicitly so that shuffle and scale can happen
    // in single step.
    cml_core_ccfft1d_op_mi_f(fft, (float*)in, (float*)out, CML_FFT_FWD);

    if (dir == -1)
    {
      if (scale != 1.f)
	cml_core_rcsvmul1_f(scale, (float*)out, (float*)out, size);
    }
    else
    {
      // Code for the inverse FFT taken from the CBE SDK Libraries
      // Overview and Users Guide, sec. 8.1.
      int const vec_size = 4;
      vector float* start = (vector float*)out;
      vector float* end   = start + 2 * size / vec_size;
      vector float  s0, s1, e0, e1;
      vector unsigned int mask = (vector unsigned int){-1, -1, 0, 0};
      vector float vscale = spu_splats(scale);
      unsigned int i;
      
      // Scale the output vector and swap the order of the outputs.
      // Note: there are two float values for each of 'n' complex values.
      s0 = e1 = *start;
      for (i = 0; i < size / vec_size; ++i) 
      {
	s1 = *(start + 1);
	e0 = *(--end);
	
	*start++ = spu_mul(spu_sel(e0, e1, mask), vscale);
	*end     = spu_mul(spu_sel(s0, s1, mask), vscale);
	s0 = s1;
	e1 = e0;
      }
    }
  }

  // Member data
  size_t      size;
  int         dir;
  float       scale;

  fft1d_f*    fft;

  static char buf1[FFT_BUF1_SIZE_BYTES];
  static char buf2[FFT_BUF2_SIZE_BYTES];
};

#endif // VSIP_OPT_UKERNEL_KERNELS_CBE_ACCEL_CFFT_F_HPP
