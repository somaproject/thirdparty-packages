/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/olayX/zvmmul_row_f.c
    @author  Jules Bergmann
    @date    2008-11-07
    @brief   VSIPL++ Library: Split-complex vector-matrix-multiply kernel.
*/

#include <spu_intrinsics.h>
#include <assert.h>
#include <alf_accel.h>
#include <cml.h>

#include <vsip/opt/cbe/vmmul_params.h>

#define _ALF_MAX_SINGLE_DT_SIZE 16*1024

static inline void
add_vector_f(
  void*             entries, 
  alf_data_addr64_t ea,
  unsigned long     length
  )
{
  unsigned long const max_length = _ALF_MAX_SINGLE_DT_SIZE / sizeof(float);

  while (length > 0)
  {
    unsigned long cur_length = (length > max_length) ? max_length : length;
    ALF_ACCEL_DTL_ENTRY_ADD(entries, cur_length, ALF_DATA_FLOAT, ea);
    length -= cur_length;
    ea     += cur_length * sizeof(float);
  } 
}


int
input_zvmmul_row_f(
  void*        p_context,
  void*        p_params, 
  void*        entries, 
  unsigned int current_count, 
  unsigned int total_count)
{
  unsigned int const  FP     = 1; // Split-complex data: 1 floats per point.
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;

#if PPU_IS_32BIT
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_IN, 0);

  add_vector_f(entries, 
	       params->ea_input_matrix_re +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  add_vector_f(entries, 
	       params->ea_input_matrix_im +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  if (current_count == 0)
  {
    add_vector_f(entries, params->ea_input_vector_re, length);
    add_vector_f(entries, params->ea_input_vector_im, length);
  }

  ALF_ACCEL_DTL_END(entries);
#else
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_IN, 0);
  add_vector_f(entries, 
	       params->ea_input_matrix_re +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  ALF_ACCEL_DTL_END(entries);

  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_IN, length*sizeof(float));
  add_vector_f(entries, 
	       params->ea_input_matrix_im +
	       current_count * FP * params->input_stride * sizeof(float),
	       length);
  ALF_ACCEL_DTL_END(entries);

  if (current_count == 0)
  {
    ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_IN, 2*length*sizeof(float));
    add_vector_f(entries, params->ea_input_vector_re, length);
    ALF_ACCEL_DTL_END(entries);
    
    ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_IN, 3*length*sizeof(float));
    add_vector_f(entries, params->ea_input_vector_im, length);
    ALF_ACCEL_DTL_END(entries);
  }
#endif

  return 0;
}



int
output_zvmmul_row_f(
  void*        p_context,
  void*        p_params, 
  void*        entries, 
  unsigned int current_count, 
  unsigned int total_count)
{
  unsigned int const  FP     = 1; // Split-complex data: 1 floats per point.
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;

#if PPU_IS_32BIT
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_OUT, 0);

  add_vector_f(entries, 
	       params->ea_output_matrix_re +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  add_vector_f(entries, 
	       params->ea_output_matrix_im +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);

  ALF_ACCEL_DTL_END(entries);
#else
  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_OUT, 0);
  add_vector_f(entries, 
	       params->ea_output_matrix_re +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  ALF_ACCEL_DTL_END(entries);

  ALF_ACCEL_DTL_BEGIN(entries, ALF_BUF_OVL_OUT, length*sizeof(float));
  add_vector_f(entries, 
	       params->ea_output_matrix_im +
	       current_count * FP * params->output_stride * sizeof(float),
	       length);
  ALF_ACCEL_DTL_END(entries);
#endif

  return 0;
}



int
kernel_zvmmul_row_f(
  void*        p_context,
  void*        p_params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int iter_count)
{
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;
  static float* save_a_re;
  static float* save_a_im;

  assert(params->length >= VSIP_IMPL_MIN_VMMUL_SIZE);
  assert(params->length <= VSIP_IMPL_MAX_VMMUL_SIZE);


  float *a_re = (float *)inout + 2*length;
  float *a_im = (float *)inout + 3*length;
  float *b_re = (float *)inout + 0 * length;
  float *b_im = (float *)inout + 1 * length;
  float *r_re = (float *)inout + 0 * length;
  float *r_im = (float *)inout + 1 * length;

  if (iter == 0)
  {
    save_a_re = malloc(2*length*sizeof(float));
    save_a_im = save_a_re + length;

    int i;
    for (i=0; i<length; ++i)
    {
      save_a_re[i] = a_re[i+params->shift];
      save_a_im[i] = a_im[i+params->shift];
    }
    // cml_core_vcopy1_f(a_re, save_a_re, length);
    // cml_core_vcopy1_f(a_im, save_a_im, length);
  }

  cml_zvmul1_f(save_a_re, save_a_im, b_re, b_im, r_re, r_im, length);

  if (iter == iter_count-1)
    free(save_a_re);

  return 0;
}
