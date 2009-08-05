/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/olay2/xvmul_f.c
    @author  Jules Bergmann
    @date    2008-10-23
    @brief   VSIPL++ Library: overlay based vector-multiply kernels.
*/

#include <alf_accel.h>
#include <cml.h>

#include <vsip/opt/cbe/vmul_params.h>



int kernel_vmul_f(
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

  cml_vmul1_f(a, b, r, length);

  return 0;
}

int kernel_cvmul_f(
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

  cml_cvmul1_f(a, b, r, length);

  return 0;
}

int kernel_zvmul_f(
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

  cml_zvmul1_f(a_re, a_im, b_re, b_im, r_re, r_im, length);

  return 0;
}
