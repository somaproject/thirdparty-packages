/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/overlay_params.h
    @author  Jules Bergmann
    @date    2008-10-23
    @brief   VSIPL++ Library: Parameters for the Overlay kernel.
*/

#ifndef VSIP_OPT_CBE_OVERLAY_PARAMS_H
#define VSIP_OPT_CBE_OVERLAY_PARAMS_H

#include <vsip/opt/cbe/vmul_params.h>
#include <vsip/opt/cbe/fft_params.h>
#include <vsip/opt/cbe/vmmul_params.h>

/***********************************************************************
  Definitions
***********************************************************************/

typedef enum
{
  overlay_cfft_f,
  overlay_zfft_f,

  overlay_vmul_f,
  overlay_cvmul_f,
  overlay_zvmul_f,
  overlay_vadd_f,
  overlay_cvadd_f,
  overlay_zvadd_f,
  overlay_mtrans_f,

  overlay_vmmul_row_f,
  overlay_cvmmul_row_f,
  overlay_zvmmul_row_f,

  overlay_zvmmul_col_f,
  overlay_zrvmmul_row_f,
  overlay_zrvmmul_col_f,
} overlay_op_type;

typedef struct
{
  unsigned long long code_ea;
  int                code_size;
  int                cmd;
} Common_params;

typedef union
{
  Common_params      c;
  Fft_params         cfft;
  Fft_split_params   zfft;
  Vmul_params        vbin;	// non-split element-wise binary vector op
  Vmul_split_params  zvbin;	//     split element-wise binary vector op
  Vmmul_params       vmmul;
  Vmmul_split_params zvmmul;
  char pad[512];
} Overlay_params;



// All overlay kernel uses should use these sizes, so that caching does
// not reload the kernel.

#define VSIP_IMPL_OVERLAY_STACK_SIZE 4096
#define VSIP_IMPL_OVERLAY_BUFFER_SIZE 65536
#define VSIP_IMPL_OVERLAY_DTL_SIZE 128



#endif // VSIP_OPT_CBE_OVERLAY_PARAMS_H
