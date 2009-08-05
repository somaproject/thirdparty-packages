/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/alf_fconvm_c.c
    @author  Don McCoy
    @date    2007-03-28
    @brief   VSIPL++ Library: Kernel to compute multiple fast convolutions
               (with distinct convolutions for each row of the matrix).
*/

#include <sys/time.h>
#include <spu_mfcio.h>
#include <alf_accel.h>
#include <assert.h>
#include <cml.h>
#include <cml_core.h>

#include <vsip/opt/cbe/fconv_params.h>

#define _ALF_MAX_SINGLE_DT_SIZE 16*1024



int 
input(
    void*        context,
    void*        params,
    void*        entries,
    unsigned int current_count,
    unsigned int total_count)
{
  unsigned int const FP = 2; // Complex data: 2 floats per point.

  Fastconv_params* fc = (Fastconv_params *)params;
  assert(fc->elements * FP *sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;

#if PPU_IS_32BIT
  // Transfer input.
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, 0);
  ea = fc->ea_input + 
    current_count * FP * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, FP * fc->elements, ALF_DATA_FLOAT, ea);

  // Transfer kernel.
  ea = fc->ea_kernel + 
    current_count * FP * fc->kernel_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, FP * fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);
#else
  // Transfer input.
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, 0);
  ea = fc->ea_input + 
    current_count * FP * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, FP * fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);

  // Transfer kernel.
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, FP*fc->elements*sizeof(float));
  ea = fc->ea_kernel + 
    current_count * FP * fc->kernel_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, FP * fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);
#endif

  return 0;
}



int 
output(
    void*        context,
    void*        params,
    void*        entries,
    unsigned int current_count,
    unsigned int total_count)
{
  unsigned int const FP = 2; // Complex data: 2 floats per point.

  Fastconv_params* fc = (Fastconv_params *)params;
  assert(fc->elements * FP *sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;

  // Transfer output.
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OUT, 0);
  ea = fc->ea_output + 
    current_count * FP * fc->output_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, FP * fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);

  return 0;
}



int 
kernel(
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

  Fastconv_params* fc = (Fastconv_params *)params;
  unsigned int fft_size = fc->elements;
  assert(fft_size <= VSIP_IMPL_MAX_FCONV_SIZE);

  // Reinitialize the FFT object if the fft size changes.
  if (iter == 0 && fft_size != current_size)
  {
    if (obj)
    {
      free(buf);
      cml_fft1d_destroy_f_alloc(obj);
    }
    int rt = cml_fft1d_setup_f_alloc(&obj, CML_FFT_CC, fft_size);
    assert(rt && obj != NULL);
    buf = (char*)memalign(16, cml_zzfft1d_buf_size_f(obj));
    assert(buf != NULL);
    current_size = fft_size;
  }

  float* in    = (float*)input;
  float* coeff = (float*)input + 2 * fft_size;
  float* out   = (float*)output;

  // Perform the forward FFT on the kernel, in place, but
  // only if requested (this step is often done in advance).
  if (fc->transform_kernel)
  {
    cml_ccfft1d_ip_f(obj, coeff, CML_FFT_FWD, buf);
  }

  cml_ccfft1d_ip_f(obj, in, CML_FFT_FWD, buf);
  cml_cvmul1_f(coeff, in, out, fft_size);
  cml_ccfft1d_ip_f(obj, out, CML_FFT_INV, buf);
  cml_core_rcsvmul1_f(1.f / fft_size, out, out, fft_size);

  return 0;
}

ALF_ACCEL_EXPORT_API_LIST_BEGIN
  ALF_ACCEL_EXPORT_API ("input", input);
  ALF_ACCEL_EXPORT_API ("output", output); 
  ALF_ACCEL_EXPORT_API ("kernel", kernel);
ALF_ACCEL_EXPORT_API_LIST_END
