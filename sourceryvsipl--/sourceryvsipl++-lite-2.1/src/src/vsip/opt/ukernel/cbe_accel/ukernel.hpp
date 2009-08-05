/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/cbe_accel/ukernel.hpp
    @author  Jules Bergmann
    @date    2008-06-10
    @brief   VSIPL++ Library: User-defined Kernel.
*/

#ifndef VSIP_OPT_UKERNEL_CBE_ACCEL_UKERNEL_HPP
#define VSIP_OPT_UKERNEL_CBE_ACCEL_UKERNEL_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/ukernel_params.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

struct Pinfo
{
  unsigned int dim;            // dimensions in this sub-block
  unsigned int l_total_size;   // total elements for this (local) iteration
  unsigned int l_offset[3];    // offset to beginning of data (if overlap
                               //  is requested, or alignment is required 
                               //  for DMA)
  unsigned int l_size[3];      // elements per dimension for this iteration
  signed int   l_stride[3];    // next-element stride in each dimension
  signed int   g_offset[3];    // local chunk's offset in global view
  signed int   o_leading[3];   // leading overlap
  signed int   o_trailing[3];  // trailing overlap
};



struct Spu_kernel
{
  static unsigned int const in_argc  = 1;
  static unsigned int const out_argc = 1;
  static unsigned int const pre_argc = 0;
  typedef Empty_params param_type;

  static bool const in_place        = false;

  void init_rank(int /*rank*/, int /*nspe*/) {}

  template <typename ParamT>
  void init(ParamT&) {}

  void fini() {}

  template <typename T>
  void pre_compute(T, Pinfo const&) {}
};





#endif // VSIP_OPT_UKERNEL_CBE_ACCEL_UKERNEL_HPP
