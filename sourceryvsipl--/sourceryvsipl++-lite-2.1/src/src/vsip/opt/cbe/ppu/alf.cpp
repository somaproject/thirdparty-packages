/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/alf.cpp
    @author  Stefan Seefeld
    @date    2008-04-18
    @brief   VSIPL++ Library: Wrappers and traits to bridge with IBMs ALF.
*/

#include <alf_cache.h>
#include "alf.hpp"
#include <iostream>

namespace vsip
{
namespace impl
{
namespace cbe
{

Task::Task(
  char const* library,
  char const* image, length_type ssize, length_type psize,
  length_type isize, length_type osize, length_type iosize,
  length_type tsize, alf_handle_t alf, unsigned int spes)
  : image_(image),
    ssize_(ssize),
    psize_(psize),
    isize_(isize),
    osize_(osize),
    iosize_(iosize),
    tsize_(tsize)
{
  (void)spes;
  task_desc desc;

  cached_alf_task_desc_init(desc);
  desc.wb_parm_ctx_buf_size   = psize;
  desc.wb_in_buf_size         = isize;
  desc.wb_out_buf_size        = osize;
  desc.wb_inout_buf_size      = iosize;
  desc.num_dtl_entries        = tsize;
  desc.max_stack_size         = ssize;
  desc.accel_library_ref_l    = library;
  desc.accel_image_ref_l      = image_;
  desc.accel_kernel_ref_l     = "kernel";
  desc.accel_input_dtl_ref_l  = "input";
  desc.accel_output_dtl_ref_l = "output";
  desc.tsk_ctx_data_type      = ALF_DATA_BYTE;

  int status = cached_alf_task_create(alf, &desc, &task_);
  assert(status >= 0);
}

void Task::destroy()
{
  if (task_)
    cached_alf_task_destroy(task_);
}

} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip
