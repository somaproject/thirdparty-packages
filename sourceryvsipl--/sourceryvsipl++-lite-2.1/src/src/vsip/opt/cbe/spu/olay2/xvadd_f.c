/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/olay2/xvadd_f.c
    @author  Jules Bergmann
    @date    2008-10-23
    @brief   VSIPL++ Library: Overlay based vector-multiply kernels.
*/

#include <alf_accel.h>
#include <cml.h>

#include <vsip/opt/cbe/vmul_params.h>

static void
zvadd1_f(
  float const* a_re,
  float const* a_im,
  float const* b_re,
  float const* b_im,
  float*       r_re,
  float*       r_im,
  int          length)
{
  while (length >= 4)
  {
    vector float tr0 = *(vector float*)a_re + *(vector float*)b_re;
    vector float ti0 = *(vector float*)a_im + *(vector float*)b_im;
    *(vector float*)r_re = tr0;
    *(vector float*)r_im = ti0;
    length -= 4;
    r_re += 4; r_im += 4;
    a_re += 4; a_im += 4;
    b_re += 4; b_im += 4;
  }

  while (length > 0)
  {
    *r_re = *a_re + *b_re;
    *r_im = *a_im + *b_im;
    --length;
    ++r_re; ++r_im;
    ++a_re; ++a_im;
    ++b_re; ++b_im;
  }
}



static void
vadd1_f(
  float const* a,
  float const* b,
  float*       r,
  int          length)
{
  while (length >= 4)
  {
    vector float tr0 = *(vector float*)a + *(vector float*)b;
    *(vector float*)r = tr0;
    length -= 4;
    r += 4;
    a += 4;
    b += 4;
  }

  while (length > 0)
  {
    *r = *a + *b;
    --length;
    ++r; ++a; ++b;
  }
}



int kernel_vadd_f(
  void*        p_context,
  void*        p_params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int n)
{
  Vmul_params* params = (Vmul_params*)p_params;
  int length = params->length;

  float *a = (float *)inout + 0 * length;
  float *b = (float *)inout + 1 * length;
  float *r = (float *)inout + 0 * length;

  vadd1_f(a, b, r, length);

  return 0;
}

int kernel_cvadd_f(
  void*        p_context,
  void*        p_params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int n)
{
  Vmul_params* params = (Vmul_params*)p_params;
  int length = params->length;

  float *a = (float *)inout + 0 * length;
  float *b = (float *)inout + 2 * length;
  float *r = (float *)inout + 0 * length;

  vadd1_f(a, b, r, 2*length);

  return 0;
}

int kernel_zvadd_f(
  void*        p_context,
  void*        p_params,
  void*        input,
  void*        output,
  void*        inout,
  unsigned int iter,
  unsigned int n)
{
  Vmul_split_params* params = (Vmul_split_params*)p_params;
  int length = params->length;

  float *a_re = (float *)inout + 0 * length;
  float *a_im = (float *)inout + 1 * length;
  float *b_re = (float *)inout + 2 * length;
  float *b_im = (float *)inout + 3 * length;
  float *r_re = (float *)inout + 0 * length;
  float *r_im = (float *)inout + 1 * length;

  zvadd1_f(a_re, a_im, b_re, b_im, r_re, r_im, length);

  return 0;
}
