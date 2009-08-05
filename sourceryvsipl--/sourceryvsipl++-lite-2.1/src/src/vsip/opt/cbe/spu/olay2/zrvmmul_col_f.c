/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/olayX/zrvmmul_col_f.c
    @author  Jules Bergmann
    @date    2008-10-28
    @brief   VSIPL++ Library: Split-complex vector-matrix-multiply kernel.
*/

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
input_zrvmmul_col_f(
  void*        p_context,
  void*        p_params, 
  void*        entries, 
  unsigned int current_count, 
  unsigned int total_count)
{

  unsigned int const  FP     = 1; // Split-complex data: 1 floats per point.
  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  unsigned long       length = params->length;

  //  Set v_len to minimum of 4 since transfers must be a multiple of 16B.
  unsigned long       v_len  = total_count;
  v_len += params->shift;
  if      (v_len == 0)     v_len = 4;
  else if (v_len % 4 != 0) v_len += 4 - (v_len % 4);

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
    add_vector_f(entries, params->ea_input_vector_re, v_len);
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
    add_vector_f(entries, params->ea_input_vector_re, v_len);
    ALF_ACCEL_DTL_END(entries);
  }

#endif

  return 0;
}



int 
output_zrvmmul_col_f(
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

int kernel_zrvmmul_col_f(
  void*        p_context,
  void*        p_params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int iter_count)
{
  static float* save_a_re;
  static float* save_a_im;

  Vmmul_split_params* params = (Vmmul_split_params*)p_params;
  int length = params->length;
  unsigned long       v_len  = iter_count;
  v_len += params->shift;
  if      (v_len == 0)     v_len = 4;
  else if (v_len % 4 != 0) v_len += 4 - (v_len % 4);

  float *a_re = (float *)inout + 2*length;
  float *b_re = (float *)inout + 0 * length;
  float *b_im = (float *)inout + 1 * length;
  float *r_re = (float *)inout + 0 * length;
  float *r_im = (float *)inout + 1 * length;

  if (iter == 0)
  {
    save_a_re = malloc(1*iter_count*sizeof(float));

    int i;
    for (i=0; i<iter_count; ++i)
      save_a_re[i] = a_re[i+params->shift];
  }

#if DEBUG
  printf("zrvmmul_col_f (%d/%d): %f\n",
	 iter, iter_count,
	 save_a_re[iter]);
#endif


  // 081230: CML function incorrect (issue #243, regressions vmmul_rz_col.cpp)
  // cml_core_rzsvmul1_f(save_a_re[iter], b_re, b_im, r_re, r_im, length);
  int i;
  for (i=0; i<length; ++i)
  {
    r_re[i] = save_a_re[iter] * b_re[i];
    r_im[i] = save_a_re[iter] * b_im[i];
  }

  if (iter == iter_count-1)
    free(save_a_re);

  return 0;
}
