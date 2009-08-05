/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/olay1/fft_z.c
    @author  Jules Bergmann
    @date    2008-10-22
    @brief   VSIPL++ Library: Kernel to compute split-complex float FFT's.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <stdio.h>
#include <alf_accel.h>
#include <assert.h>
#include <spu_mfcio.h>
#include <cml.h>
#include <cml_core.h>

#include <vsip/core/acconfig.hpp>
#include <vsip/opt/cbe/fft_params.h>

#define _ALF_MAX_SINGLE_DT_SIZE 16*1024



/***********************************************************************
  Definitions
***********************************************************************/

int kernel_zfft_f(
  void*        context,
  void*        params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int iter_max)
{
  static fft1d_f* obj;
  static char*    buf;
  static size_t   current_size = 0;

  Fft_split_params* fftp = (Fft_split_params *)params;
  unsigned int      size   = fftp->size;
  unsigned int      chunks = fftp->chunks_per_wb;
  unsigned int      i;
  int dir = fftp->direction == fwd_fft ? CML_FFT_FWD : CML_FFT_INV;

  if (iter == iter_max-1 && iter_max * chunks > fftp->chunks_per_spe)
    chunks = fftp->chunks_per_spe % chunks;

  assert(size >= MIN_FFT_1D_SIZE);
  assert(size <= MAX_FFT_1D_SIZE);

  if (iter == 0 && size != current_size)
  {
    if (obj)
    {
      free(buf);
      cml_fft1d_destroy_f_alloc(obj);
    }
    int rt = cml_fft1d_setup_f_alloc(&obj, CML_FFT_CC, size);
    assert(rt && obj != NULL);
    buf = (char*)memalign(16, cml_zzfft1d_buf_size_f(obj));
    assert(buf != NULL);
    current_size = size;
  }

  float*        in_re  = (float*)inout + 0        * size;
  float*        in_im  = (float*)inout + 1*chunks * size;
  float*        out_re = (float*)inout + 2*chunks * size;
  float*        out_im = (float*)inout + 3*chunks * size;

  for (i=0; i<chunks; ++i)
  {
    cml_zzfft1d_op_f(obj,
		     (float*)in_re  + i*size, (float*)in_im  + i*size,
		     (float*)out_re + i*size, (float*)out_im + i*size,
		     dir, buf);
  }

  if (fftp->scale != (double)1.f)
  {
    // Instead of regular split svmul:
    // cml_core_rzsvmul1_f(fftp->scale, out_re,out_im,out_re,out_im,size);
    // Take advantage of real and imag being contiguous:
    cml_core_svmul1_f(fftp->scale, out_re, out_re, 2*size*chunks);
  }

  return 0;
}
