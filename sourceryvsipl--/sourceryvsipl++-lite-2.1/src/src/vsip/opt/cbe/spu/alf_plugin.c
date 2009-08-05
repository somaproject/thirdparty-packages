/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/spu/alf_plugin.c
    @author  Jules Bergmann
    @date    2009-02-10
    @brief   VSIPL++ Library: Plugin kernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <stdio.h>
#include <alf_accel.h>
#include <assert.h>
#include <spu_mfcio.h>

#include <cml.h>

#include <vsip/core/acconfig.hpp>
#include <vsip/opt/cbe/overlay_params.h>
#include <vsip/opt/cbe/vmul_params.h>
#include <vsip/opt/cbe/spu/plugin_functions.h>

#define DEBUG_LOAD 0
#define DEBUG_CALL 0

// Actual space available is about 160 KB under default circumstances
// However, some kernels require more ALF buffer space (interp_f),
// which will reduce amount available here.
char code_buffer[128+1*1024] __attribute__((section(".mybss")));
alf_data_addr64_t code_ea = 0;
int               code_loading = 0;
Plugin_functions  functions;

kernel_function* p_kernel;
input_function*  p_input;
output_function* p_output;




/***********************************************************************
  Definitions
***********************************************************************/

int input(
  void*        context,
  void*        params,
  void*        entries,
  unsigned int iter,
  unsigned int iter_max)
{
  if (code_ea != ((Common_params*)params)->code_ea)
  {
#if DEBUG_LOAD
    printf("SPU-resident: loading code: %16llx %d (%16llx)\n",
	   ((Common_params*)params)->code_ea,
	   ((Common_params*)params)->code_size, code_ea);
    int t_code_size = ((Common_params*)params)->code_size;
#endif
    code_ea       = ((Common_params*)params)->code_ea;
    int code_size = ((Common_params*)params)->code_size;
    alf_data_addr64_t ea = code_ea;
    char* buf = code_buffer;
    while (code_size > 0)
    {
      int size = code_size > 16384 ? 16384 : code_size;
      mfc_get(buf, ea, size, 31, 0, 0);
      buf += size; ea += size; code_size -= size;
    }

    functions.f_printf        = &printf;
    functions.f_dtl_begin     = &ALF_ACCEL_DTL_BEGIN;
    functions.f_dtl_entry_add = &ALF_ACCEL_DTL_ENTRY_ADD;
    functions.f_dtl_end       = &ALF_ACCEL_DTL_END;

    mfc_write_tag_mask(1<<31);
    mfc_read_tag_status_all();
    p_kernel = *(kernel_function**)(code_buffer+0);
    p_input  = *(input_function**)(code_buffer+4);
    p_output = *(output_function**)(code_buffer+8);
#if DEBUG_LOAD
    printf("SPU-resident: kernel %08x\n", p_kernel);
    printf("SPU-resident: input  %08x\n", p_input);
    printf("SPU-resident: output %08x\n", p_output);
    int i, chksum = 0;
    for (i=0; i<t_code_size/4; ++i)
      chksum = chksum ^ ((int*)code_buffer)[i];
    printf("SPU-resident: chksum %8x\n", chksum);
#endif
  }

#if DEBUG_CALL
  printf("SPU-resident: input %d/%d\n", iter, iter_max);
#endif

  return p_input(&functions, context, params, entries, iter, iter_max);
}



int output(
  void*        context,
  void*        params,
  void*        entries,
  unsigned int iter,
  unsigned int iter_max)
{
#if DEBUG_CALL
  printf("SPU-resident: output %d/%d\n", iter, iter_max);
#endif
  int rv = p_output(&functions, context, params, entries, iter, iter_max);
#if DEBUG_CALL
  printf("SPU-resident: output %d/%d return, rv %d\n", iter, iter_max, rv);
#endif
  return rv;
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
  (void)input; (void)output;

#if DEBUG_CALL
  printf("SPU-resident: kernel %d/%d\n", iter, iter_max);
#endif
  return p_kernel(&functions, context, params, inout, iter, iter_max);
}

ALF_ACCEL_EXPORT_API_LIST_BEGIN
  ALF_ACCEL_EXPORT_API ("input", input);
  ALF_ACCEL_EXPORT_API ("output", output); 
  ALF_ACCEL_EXPORT_API ("kernel", kernel);
ALF_ACCEL_EXPORT_API_LIST_END
