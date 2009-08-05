/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/plugin.hpp
    @author  Jules Bergmann
    @date    2009-02-11
    @brief   VSIPL++ Library: Plugin support.
*/

#ifndef VSIP_OPT_CBE_PPU_PLUGIN_HPP
#define VSIP_OPT_CBE_PPU_PLUGIN_HPP

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cbe
{

void
load_plugin(
  char*&      code_ea,
  int&        code_size,
  char const* library_subdir,
  char const* kernel);

} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CBE_PPU_PLUGIN_HPP
