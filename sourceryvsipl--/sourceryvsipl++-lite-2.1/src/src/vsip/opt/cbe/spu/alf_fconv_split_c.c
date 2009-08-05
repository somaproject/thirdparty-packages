/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/alf_fconv_split_c.c
    @author  Don McCoy, Jules Bergmann
    @date    2007-02-28
    @brief   VSIPL++ Library: Kernel to compute fast convolution using
             split-complex.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#define PERFMON 0

#include <sys/time.h>
#include <spu_mfcio.h>
#include <alf_accel.h>
#include <assert.h>
#include <cml.h>
#include <cml_core.h>

#include <vsip/opt/cbe/fconv_params.h>

#define _ALF_MAX_SINGLE_DT_SIZE 16*1024

#if PERFMON
#  include "timer.h"
#  define START_TIMER(x) start_timer(x)
#  define STOP_TIMER(x)  stop_timer(x)
#else
#  define START_TIMER(x)
#  define STOP_TIMER(x)
#endif



/***********************************************************************
  Definitions
***********************************************************************/

// Convolution coefficients.
static float kernel_re[VSIP_IMPL_MAX_FCONV_SPLIT_SIZE] 
       __attribute__ ((aligned (128)));
static float kernel_im[VSIP_IMPL_MAX_FCONV_SPLIT_SIZE] 
       __attribute__ ((aligned (128)));

// Instance-id.  Used to determine when new coefficients must be loaded.
static unsigned int instance_id = 0;



void initialize(
  Fastconv_split_params* fc,
  float*                 p_kernel_re,
  float*                 p_kernel_im,
  fft1d_f*               obj,
  void*                  buf)
{
  unsigned int n    = fc->elements;

  // The kernel matches the input and output size
  mfc_get(p_kernel_re, fc->ea_kernel_re, n*sizeof(float), 31, 0, 0);
  mfc_get(p_kernel_im, fc->ea_kernel_im, n*sizeof(float), 31, 0, 0);
  mfc_write_tag_mask(1<<31);
  mfc_read_tag_status_all();

  if (fc->transform_kernel)
  {
    // Perform the forward FFT on the kernel, in place.  This only need 
    // be done once -- subsequent calls will utilize the same kernel.
    cml_zzfft1d_ip_f(obj,
		     (float*)kernel_re, (float*)kernel_im,
		     CML_FFT_FWD, buf);
  }
}



int input(
  void*        context,
  void*        params,
  void*        entries,
  unsigned int current_count,
  unsigned int total_count)
{
  (void)context;
  (void)total_count;

  Fastconv_split_params* fc = (Fastconv_split_params *)params;
  assert(fc->elements * sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;

  // Transfer input.
#if PPU_IS_32BIT
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, 0);
  ea = fc->ea_input_re + current_count * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);

  ea = fc->ea_input_im + current_count * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);

  ALF_ACCEL_DTL_END(entries);
#else
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, 0);
  ea = fc->ea_input_re + current_count * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);

  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_IN, fc->elements*sizeof(float));
  ea = fc->ea_input_im + current_count * fc->input_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);
#endif

  return 0;
}



int output(
  void*        context,
  void*        params,
  void*        entries,
  unsigned int current_count,
  unsigned int total_count)
{
  (void)context;
  (void)total_count;

  Fastconv_split_params* fc = (Fastconv_split_params *)params;
  assert(fc->elements * sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;

  // Transfer output.
#if PPU_IS_32BIT
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OUT, 0);
  ea = fc->ea_output_re + current_count * fc->output_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);

  ea = fc->ea_output_im + current_count * fc->output_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);

  ALF_ACCEL_DTL_END(entries);
#else
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OUT, 0);
  ea = fc->ea_output_re + current_count * fc->output_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);

  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OUT, fc->elements*sizeof(float));
  ea = fc->ea_output_im + current_count * fc->output_stride * sizeof(float);
  ALF_ACCEL_DTL_ENTRY_ADD(entries, fc->elements, ALF_DATA_FLOAT, ea);
  ALF_ACCEL_DTL_END(entries);
#endif

  return 0;
}



int kernel(
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

  Fastconv_split_params* fc = (Fastconv_split_params *)params;
  unsigned int fft_size = fc->elements;
  assert(fft_size <= VSIP_IMPL_MAX_FCONV_SPLIT_SIZE);

  (void)context;
  (void)iter;
  (void)iter_max;

#if PERFMON
  static acc_timer_t t1;
  static acc_timer_t t2;
  static acc_timer_t t3;
  static acc_timer_t t4;
#endif

  // Initialization establishes the weights (kernel) for the
  // convolution step and the twiddle factors for the FFTs.
  // These are loaded once per task by checking a unique
  // ID passed down from the caller.
  if (iter == 0)
  {
    if (fft_size != current_size)
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
    if (instance_id != fc->instance_id)
    {
      instance_id = fc->instance_id;
      initialize(fc, kernel_re, kernel_im, obj, buf);
#if PERFMON
      t1 = init_timer();
      t2 = init_timer();
      t3 = init_timer();
      t4 = init_timer();
#endif
    }
  }

  float*        in_re  = (float*)input  + 0 * fft_size;
  float*        in_im  = (float*)input  + 1 * fft_size;
  float *       out_re = (float*)output + 0 * fft_size;
  float *       out_im = (float*)output + 1 * fft_size;

  // Switch to frequency space
  START_TIMER(&t1);
  cml_zzfft1d_ip_f(obj,
		   (float*)in_re, (float*)in_im,
		   CML_FFT_FWD, buf);
  STOP_TIMER(&t1);

  // Perform convolution -- now a straight multiplication
  START_TIMER(&t2);
  cml_zvmul1_f(kernel_re, kernel_im, in_re, in_im, out_re, out_im, fft_size);
  STOP_TIMER(&t2);

  // Revert back the time domain
  START_TIMER(&t3);
  cml_zzfft1d_ip_f(obj,
		   (float*)out_re, (float*)out_im,
		   CML_FFT_INV, buf);
  STOP_TIMER(&t3);

  // Scale by 1/n.
  START_TIMER(&t4);
  cml_core_rzsvmul1_f(1.f / fft_size, out_re, out_im, out_re, out_im,
		      fft_size);
  STOP_TIMER(&t4);

#if PERFMON
  if (0 && iter == iter_max-1)
  {
    double total1 = timer_total(&t1);
    double total2 = timer_total(&t2);
    double total3 = timer_total(&t3);
    double total4 = timer_total(&t4);
    double fft_flops = (double)t1.count * 5 * n * log2n;
    double cvm_flops = (double)t1.count * 6 * n;
    double sca_flops = (double)t1.count * 2 * n;
    double fwd_mflops = fft_flops / (total1 * 1e6);
    double cvm_mflops = cvm_flops / (total2 * 1e6);
    double inv_mflops = fft_flops / (total3 * 1e6);
    double sca_mflops = sca_flops / (total4 * 1e6);
    printf("fwd fft: %f s  %f MFLOP/s\n", total1, fwd_mflops);
    printf("cvmul  : %f s  %f MFLOP/s\n", total2, cvm_mflops);
    printf("inv fft: %f s  %f MFLOP/s\n", total3, inv_mflops);
    printf("scale  : %f s  %f MFLOP/s\n", total4, sca_mflops);
  }
#endif

  return 0;
}

ALF_ACCEL_EXPORT_API_LIST_BEGIN
  ALF_ACCEL_EXPORT_API ("input", input);
  ALF_ACCEL_EXPORT_API ("output", output); 
  ALF_ACCEL_EXPORT_API ("kernel", kernel);
ALF_ACCEL_EXPORT_API_LIST_END
