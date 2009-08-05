/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/kernels/cbe_accel/zfft_f.hpp
    @author  Jules Bergmann
    @date    2008-08-08
    @brief   VSIPL++ Library: Split-complex fastconv ukernel.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/ukernel/kernels/cbe_accel/zfft_f.hpp>

typedef Fft_kernel kernel_type;

char Fft_kernel::buf1[2*MAX_FFT_1D_SIZE*sizeof(float)]
     __attribute((aligned(128)));
char Fft_kernel::buf2[1*MAX_FFT_1D_SIZE*sizeof(float)+128]
     __attribute((aligned(128)));

#include <vsip/opt/ukernel/cbe_accel/alf_base.hpp>
