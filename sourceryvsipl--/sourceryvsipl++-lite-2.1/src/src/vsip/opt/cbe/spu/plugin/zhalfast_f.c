/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/plugin/zfft_f.c
    @author  Jules Bergmann
    @date    2009-02-20
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
#include <vsip/opt/cbe/vmmul_params.h>
#include <vsip/opt/cbe/spu/plugin_functions.h>
#include <vsip/opt/cbe/overlay_params.h>

#define _ALF_MAX_SINGLE_DT_SIZE 16*1024


int current_size = 0;
float buf[2*VSIP_IMPL_MAX_VMMUL_SIZE + MAX_FFT_1D_SIZE + 128/sizeof(float)];



/***********************************************************************
  zfft_f Definitions
***********************************************************************/

static int
input_zfft_f(
  Plugin_functions* pf,
  void*             context,
  void*             params,
  void*             entries,
  unsigned int      iter,
  unsigned int      iter_max)
{
  Fft_split_params* fftp   = (Fft_split_params *)params;
  unsigned int      size   = fftp->size;
  unsigned int      chunks = fftp->chunks_per_wb;
  unsigned int      cur_chunks;
  unsigned int      i;
    
  assert(size * sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;

  if (iter == iter_max-1 && iter_max * chunks > fftp->chunks_per_spe)
    cur_chunks = fftp->chunks_per_spe % chunks;
  else
    cur_chunks = chunks;
    
  if (size == fftp->in_blk_stride)
  {
    size *= cur_chunks;
    cur_chunks = 1;
  }

  // Transfer input.
#if PPU_IS_32BIT
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
    
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_input_re +
      (iter * chunks + i) * fftp->in_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
    
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_input_im +
      (iter * chunks + i) * fftp->in_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
  
  (pf->f_dtl_end)(entries);
#else
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_input_re +
      (iter * chunks + i) * fftp->in_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
  (pf->f_dtl_end)(entries);
    
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, size*cur_chunks*sizeof(float));
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_input_im +
      (iter * chunks + i) * fftp->in_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
  (pf->f_dtl_end)(entries);
#endif

  return 0;
}



static int
output_zfft_f(
  Plugin_functions* pf,
  void*             context,
  void*             params,
  void*             entries,
  unsigned int      iter,
  unsigned int      iter_max)
{
  Fft_split_params* fftp   = (Fft_split_params *)params;
  unsigned int      size   = fftp->size;
  unsigned int      chunks = fftp->chunks_per_wb;
  unsigned int      cur_chunks;
  unsigned int      i;
    
  assert(size * sizeof(float) <= _ALF_MAX_SINGLE_DT_SIZE);
  alf_data_addr64_t ea;
    
  if (iter == iter_max-1 && iter_max * chunks > fftp->chunks_per_spe)
    cur_chunks = fftp->chunks_per_spe % chunks;
  else
    cur_chunks = chunks;
  
  if (size == fftp->out_blk_stride)
  {
    size *= cur_chunks;
    cur_chunks = 1;
  }
    
    // Transfer output.
#if PPU_IS_32BIT
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 2*cur_chunks*size*sizeof(float));
    
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_output_re +
      (iter * chunks + i) * fftp->out_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
    
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_output_im +
      (iter * chunks + i) * fftp->out_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
    
  (pf->f_dtl_end)(entries);
#else
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 2*cur_chunks*size*sizeof(float));
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_output_re +
      (iter * chunks + i) * fftp->out_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
  (pf->f_dtl_end)(entries);
    
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 3*cur_chunks*size*sizeof(float));
  for (i=0; i<cur_chunks; ++i)
  {
    ea = fftp->ea_output_im +
      (iter * chunks + i) * fftp->out_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, size, ALF_DATA_FLOAT, ea);
  }
  (pf->f_dtl_end)(entries);
#endif

  return 0;
}



static int
kernel_zfft_f(
  Plugin_functions* pf,
  void*             context,
  void*             params,
  void*             inout,
  unsigned int      iter,
  unsigned int      iter_max)
{
  static fft1d_f* obj;

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
    int rt = cml_fft1d_setup_f(&obj, CML_FFT_CC, size, buf + 16384);
    assert(rt && obj != NULL);
    current_size = size;
  }

  float* in_re  = (float*)inout + 0        * size;
  float* in_im  = (float*)inout + 1*chunks * size;
  float* out_re = (float*)inout + 2*chunks * size;
  float* out_im = (float*)inout + 3*chunks * size;

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



/***********************************************************************
  zvmmul_row_f Definitions
***********************************************************************/

static inline void
add_vector_f(
  Plugin_functions* pf,
  void*             entries, 
  alf_data_addr64_t ea,
  unsigned long     length
  )
{
  unsigned long const max_length = _ALF_MAX_SINGLE_DT_SIZE / sizeof(float);

  while (length > 0)
  {
    unsigned long cur_length = (length > max_length) ? max_length : length;
    (pf->f_dtl_entry_add)(entries, cur_length, ALF_DATA_FLOAT, ea);
    length -= cur_length;
    ea     += cur_length * sizeof(float);
  } 
}



static int
input_zvmmul_row_f(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params, 
  void*             entries, 
  unsigned int      current_count, 
  unsigned int      total_count)
{
  unsigned int const  FP     = 1; // Split-complex data: 1 floats per point.
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;

#if PPU_IS_32BIT
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);

  add_vector_f(pf, entries, 
	       params->ea_input_matrix_re +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  add_vector_f(pf, entries, 
	       params->ea_input_matrix_im +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  if (current_count == 0)
  {
    add_vector_f(pf, entries, params->ea_input_vector_re, length);
    add_vector_f(pf, entries, params->ea_input_vector_im, length);
  }

  (pf->f_dtl_end)(entries);
#else
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
  add_vector_f(pf, entries, 
	       params->ea_input_matrix_re +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  (pf->f_dtl_end)(entries);

  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, length*sizeof(float));
  add_vector_f(pf, entries, 
	       params->ea_input_matrix_im +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  (pf->f_dtl_end)(entries);

  if (current_count == 0)
  {
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 2*length*sizeof(float));
    add_vector_f(pf, entries, params->ea_input_vector_re, length);
    (pf->f_dtl_end)(entries);
    
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 3*length*sizeof(float));
    add_vector_f(pf, entries, params->ea_input_vector_im, length);
    (pf->f_dtl_end)(entries);
  }
#endif

  return 0;
}



static int
output_zvmmul_row_f(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params, 
  void*             entries, 
  unsigned int      current_count, 
  unsigned int      total_count)
{
  unsigned int const  FP     = 1; // Split-complex data: 1 floats per point.
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;

#if PPU_IS_32BIT
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);

  add_vector_f(pf, entries, 
	       params->ea_output_matrix_re +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  add_vector_f(pf, entries, 
	       params->ea_output_matrix_im +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);

  (pf->f_dtl_end)(entries);
#else
  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);
  add_vector_f(pf, entries, 
	       params->ea_output_matrix_re +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  (pf->f_dtl_end)(entries);

  (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, length*sizeof(float));
  add_vector_f(pf, entries, 
	       params->ea_output_matrix_im +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  (pf->f_dtl_end)(entries);
#endif

  return 0;
}






static int
kernel_zvmmul_row_f(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params,
  void*             inout,
  unsigned int      iter,
  unsigned int      iter_count)
{
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       size   = params->length;

  assert(size >= VSIP_IMPL_MIN_VMMUL_SIZE);
  assert(size <= VSIP_IMPL_MAX_VMMUL_SIZE);


  float *a_re = (float *)inout + 2 * size;
  float *a_im = (float *)inout + 3 * size;
  float *b_re = (float *)inout + 0 * size;
  float *b_im = (float *)inout + 1 * size;
  float *r_re = (float *)inout + 0 * size;
  float *r_im = (float *)inout + 1 * size;

  float* save_a_re = buf;
  float* save_a_im = buf + size;

  if (iter == 0)
  {
    int i;

    for (i=0; i<size; ++i)
    {
      save_a_re[i] = a_re[i+params->shift];
      save_a_im[i] = a_im[i+params->shift];
    }
  }

  cml_zvmul1_f(save_a_re, save_a_im, b_re, b_im, r_re, r_im, size);

  return 0;
}



int
input(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params, 
  void*             entries, 
  unsigned int      current_count, 
  unsigned int      total_count)
{
  switch (((Vmmul_split_params*)p_params)->cmd)
  {
  case overlay_zvmmul_row_f:
    return input_zvmmul_row_f(pf, p_context, p_params, entries,
			      current_count, total_count);
  case overlay_zfft_f:
    return input_zfft_f(pf, p_context, p_params, entries,
			current_count, total_count);
  }
  return 1;
}



int
output(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params, 
  void*             entries, 
  unsigned int      current_count, 
  unsigned int      total_count)
{
  switch (((Vmmul_split_params*)p_params)->cmd)
  {
  case overlay_zvmmul_row_f:
    return output_zvmmul_row_f(pf, p_context, p_params, entries,
			       current_count, total_count);
  case overlay_zfft_f:
    return output_zfft_f(pf, p_context, p_params, entries,
			 current_count, total_count);
  }
  return 1;
}



int
kernel(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params,
  void*             inout,
  unsigned int      iter,
  unsigned int      iter_count)
{
  switch (((Vmmul_split_params*)p_params)->cmd)
  {
  case overlay_zvmmul_row_f:
    return kernel_zvmmul_row_f(pf, p_context, p_params, inout,
			       iter, iter_count);
  case overlay_zfft_f:
    return kernel_zfft_f(pf, p_context, p_params, inout,
			 iter, iter_count);
  }
  return 1;
}
