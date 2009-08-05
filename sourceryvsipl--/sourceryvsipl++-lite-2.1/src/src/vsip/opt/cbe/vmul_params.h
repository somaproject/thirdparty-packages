/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/vmul_params.h
    @author  Jules Bergmann
    @date    2007-02-25
    @brief   VSIPL++ Library: Parameters for vmul_c and vmul_s kernels.
*/

#ifndef VSIP_OPT_CBE_VMUL_PARAMS_H
#define VSIP_OPT_CBE_VMUL_PARAMS_H

/***********************************************************************
  Definitions
***********************************************************************/

#ifdef _cplusplus
namespace vsip
{
namespace impl
{
namespace cbe
{
#endif

// Structures used in DMAs should be sized in multiples of 128-bits

typedef struct
{
  unsigned long long code_ea;
  int                code_size;
  int                cmd;

  unsigned int       length;
  unsigned int       a_blk_stride;
  unsigned int       b_blk_stride;
  unsigned int       r_blk_stride;
  unsigned long long a_ptr; // input
  unsigned long long b_ptr; // input
  unsigned long long r_ptr; // result = A * B
} Vmul_params;

typedef struct
{
  unsigned long long code_ea;
  int                code_size;
  int                cmd;

  unsigned int       length;
  unsigned int       a_blk_stride;
  unsigned int       b_blk_stride;
  unsigned int       r_blk_stride;

  unsigned long long a_im_ptr;
  unsigned long long a_re_ptr;
  unsigned long long b_im_ptr;
  unsigned long long b_re_ptr;

  unsigned long long r_im_ptr;
  unsigned long long r_re_ptr;
  unsigned int       command;
} Vmul_split_params;

#ifdef _cplusplus
} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip
#endif

#endif // VSIP_OPT_CBE_VMUL_PARAMS_H
