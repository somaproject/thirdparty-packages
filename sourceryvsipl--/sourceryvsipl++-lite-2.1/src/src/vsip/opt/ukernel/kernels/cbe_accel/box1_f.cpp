/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/box1_f.cpp
    @author  Jules Bergmann
    @date    2008-08-08
    @brief   VSIPL++ Library: 1D Box filter ukernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/kernels/cbe_accel/box1_f.hpp>

typedef Box1_kernel kernel_type;

#include <vsip/opt/ukernel/cbe_accel/alf_base.hpp>
