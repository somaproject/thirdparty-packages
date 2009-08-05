/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/plugin/xvmul_f.c
    @author  Jules Bergmann
    @date    2009-02-25
    @brief   VSIPL++ Library: Plugin based vector-mul kernels.
*/

#include <alf_accel.h>
#include <cml.h>

#include <vsip/opt/cbe/vmul_params.h>
#include <vsip/opt/cbe/spu/plugin_functions.h>
#include <vsip/opt/cbe/overlay_params.h>

int input(
  Plugin_functions* pf,
  void*             context,
  void*             params,
  void*             entries,
  unsigned int      iter,
  unsigned int      iter_max)
{
  alf_data_addr64_t ea;

  switch (((Vmul_params*)params)->cmd)
  {
  case overlay_vmul_f:
  case overlay_cvmul_f:
  {
    // words-per-point
    int wpp = (((Common_params*)params)->cmd == overlay_vmul_f) ? 1 : 2;

    Vmul_params* p = (Vmul_params*)params;
    alf_data_addr64_t ea;

#if PPU_IS_32BIT
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);

    // Transfer input A
    ea = p->a_ptr + iter * p->a_blk_stride * wpp * sizeof(float);
    (pf->f_dtl_entry_add)(entries, wpp*p->length, ALF_DATA_FLOAT, ea);
    
    // Transfer input B.
    ea = p->b_ptr + iter * p->b_blk_stride * wpp * sizeof(float);
    (pf->f_dtl_entry_add)(entries, wpp*p->length, ALF_DATA_FLOAT, ea);

    (pf->f_dtl_end)(entries);
#else
    // Transfer input A real
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
    ea = p->a_ptr + iter * p->a_blk_stride * wpp * sizeof(float);
    (pf->f_dtl_entry_add)(entries, wpp*p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);

    // Transfer input B.
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, wpp*p->length*sizeof(float));
    ea = p->b_ptr + iter * p->b_blk_stride * wpp * sizeof(float);
    (pf->f_dtl_entry_add)(entries, wpp*p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
#endif
  }
  break;

  case overlay_zvmul_f:
  {
    Vmul_split_params* p = (Vmul_split_params*)params;
#if PPU_IS_32BIT
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
    
    // Transfer input A real
    ea = p->a_re_ptr + iter * p->a_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    
    ea = p->a_im_ptr + iter * p->a_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    
    // Transfer input B.
    ea = p->b_re_ptr + iter * p->b_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);

    ea = p->b_im_ptr + iter * p->b_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    
    (pf->f_dtl_end)(entries);
#else
    // Transfer input A real
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
    ea = p->a_re_ptr + iter * p->a_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
    
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 1*p->length*sizeof(float));
    ea = p->a_im_ptr + iter * p->a_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
    
    // Transfer input B.
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 2*p->length*sizeof(float));
    ea = p->b_re_ptr + iter * p->b_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
    
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 3*p->length*sizeof(float));
    ea = p->b_im_ptr + iter * p->b_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
#endif
  }
  break;
  }

  return 0;
}



int output(
  Plugin_functions* pf,
  void*             context,
  void*             params,
  void*             entries,
  unsigned int      iter,
  unsigned int      iter_max)
{
  alf_data_addr64_t ea;

  switch (((Vmul_params*)params)->cmd)
  {
  case overlay_vmul_f:
  case overlay_cvmul_f:
    {
      int wpp = (((Common_params*)params)->cmd == overlay_vmul_f) ? 1 : 2;
      
      Vmul_params* p = (Vmul_params*)params;
      alf_data_addr64_t ea;
      
      // Transfer output R.
      (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);
      ea = p->r_ptr + iter *  p->r_blk_stride * wpp * sizeof(float);
      (pf->f_dtl_entry_add)(entries, wpp*p->length, ALF_DATA_FLOAT, ea);
      (pf->f_dtl_end)(entries);
    }
    break;

  case overlay_zvmul_f:
  {
    Vmul_split_params* p = (Vmul_split_params*)params;
    // Transfer output R.
#if PPU_IS_32BIT
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);
    
    ea = p->r_re_ptr + iter *  p->r_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    
    ea = p->r_im_ptr + iter *  p->r_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    
    (pf->f_dtl_end)(entries);
#else
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);
    ea = p->r_re_ptr + iter *  p->r_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
  
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, p->length*sizeof(float));
    ea = p->r_im_ptr + iter *  p->r_blk_stride * sizeof(float);
    (pf->f_dtl_entry_add)(entries, p->length, ALF_DATA_FLOAT, ea);
    (pf->f_dtl_end)(entries);
#endif
  }
    break;
  }

  return 0;
}

int kernel(
  Plugin_functions* pf,
  void*             p_context,
  void*             p_params,
  void*             inout,
  unsigned int      iter,
  unsigned int      n)
{
  Vmul_split_params* p = (Vmul_split_params*)p_params;
  int length = p->length;

  switch (p->cmd)
  {
  case overlay_vmul_f:
    {
      float *a = (float *)inout + 0 * length;
      float *b = (float *)inout + 1 * length;
      float *r = (float *)inout + 0 * length;
      
      cml_vmul1_f(a, b, r, length);
    }
    break;

  case overlay_cvmul_f:
    {
      float *a = (float *)inout + 0 * length;
      float *b = (float *)inout + 2 * length;
      float *r = (float *)inout + 0 * length;
      
      cml_cvmul1_f(a, b, r, length);
    }
    break;

  case overlay_zvmul_f:
    {
      float *a_re = (float *)inout + 0 * length;
      float *a_im = (float *)inout + 1 * length;
      float *b_re = (float *)inout + 2 * length;
      float *b_im = (float *)inout + 3 * length;
      float *r_re = (float *)inout + 0 * length;
      float *r_im = (float *)inout + 1 * length;
      
      cml_zvmul1_f(a_re, a_im, b_re, b_im, r_re, r_im, length);
    }
    break;
  }

  return 0;
}


