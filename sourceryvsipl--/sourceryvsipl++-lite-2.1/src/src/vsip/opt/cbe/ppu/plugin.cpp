/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/plugin.cpp
    @author  Jules Bergmann
    @date    2009-02-11
    @brief   VSIPL++ Library: Plugin support.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#define DEBUG 0

#if DEBUG
# include <iostream>
#endif
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <cstdlib>
#include <map>

#include <vsip/opt/cbe/ppu/plugin.hpp>
#include <vsip/core/allocation.hpp>




/***********************************************************************
  Definitions
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace cbe
{

typedef  std::pair<char*, int> plugin_type;
std::map<std::string, plugin_type> plugin_map;

// Load plugin code into PPU memory.

void
load_plugin(
  char*&      code_ea,
  int&        code_size,
  char const* library_subdir,
  char const* kernel)
{
  std::string key = std::string(library_subdir) + std::string(kernel);
  char path[1024], file[1024];
  struct stat sb;

  if (plugin_map.find(key) != plugin_map.end())
  {
#if DEBUG
    std::cout << "load_plugin(" << library_subdir << ", " << kernel << ")\n";
#endif
    plugin_type entry = plugin_map[key];
    code_ea   = entry.first;
    code_size = entry.second;
    return;
  }

  char* env = getenv("ALF_LIBRARY_PATH");
  if (!env) VSIP_IMPL_THROW(std::runtime_error("ALF_LIBRARY_PATH is not defined."));
  strcpy(path, env);

  char* dir = strtok(path, ":");
  while(dir)
  {
    strcpy(file, dir);            strcat(file, "/");
    strcat(file, library_subdir); strcat(file, "/");
    strcat(file, kernel);         strcat(file, ".img");

#if DEBUG
    std::cout << "load_plugin(" << library_subdir << ", " << kernel 
	      << ") try " << file << std::endl;
#endif

    int rt = stat(file, &sb);

    if (rt == 0 && S_ISREG(sb.st_mode))
    {
      int real_size = sb.st_size;
      // pad code size up to 128 B (for better DMA performance)
      code_size = ((real_size) + (-(unsigned)(real_size) % 128));
#if DEBUG
    std::cout << "load_plugin(" << library_subdir << ", " << kernel 
	      << "): size " << code_size << " (" << real_size << ") bytes\n";
#endif
      code_ea = alloc_align<char>(128, code_size);

      std::ifstream stream(file, std::ios::in | std::ios::binary);
      stream.read(code_ea, real_size);
      stream.close();

      plugin_map[key] = plugin_type(code_ea, code_size);
      return;
    }

    dir = strtok(NULL, ":");
  }

  std::string error = "Can't find ukernel plugin: ";
  error += library_subdir;
  error += ' ';
  error += kernel;

  VSIP_IMPL_THROW(std::runtime_error(error));
}

} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip
