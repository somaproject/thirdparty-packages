/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip_csl/png.hpp
    @author  Stefan Seefeld
    @date    2006-03-21
    @brief   VSIPL++ Library: 
*/

#ifndef VSIP_CSL_PNG_HPP
#define VSIP_CSL_PNG_HPP

#include <png.h>
#include <string>

namespace vsip_csl
{
namespace png
{

enum color_type
{ 
  gray = PNG_COLOR_TYPE_GRAY,
  grayalpha = PNG_COLOR_TYPE_GRAY_ALPHA,
  palette = PNG_COLOR_TYPE_PALETTE,
  rgb = PNG_COLOR_TYPE_RGB,
  rgba = PNG_COLOR_TYPE_RGB_ALPHA,
  maskpalette = PNG_COLOR_MASK_PALETTE,
  maskcolor = PNG_COLOR_MASK_COLOR,
  maskalpha = PNG_COLOR_MASK_ALPHA
};

enum interlace_type
{
  none = PNG_INTERLACE_NONE,
  adam7 = PNG_INTERLACE_ADAM7,
  last = PNG_INTERLACE_LAST
};
  
struct info
{
  unsigned long width;
  unsigned long height;
  unsigned long rowbytes;
  unsigned short depth;
  color_type colortype;
  unsigned short compression;
  unsigned short filter;
  interlace_type interlace;
};

class decoder 
{
  static size_t const magic_ = 8;
public:
  decoder(std::streambuf *, info &);
  ~decoder();
  void decode(unsigned char *, unsigned long size);

private:
  std::streambuf *input_;
  png_struct     *png_;
  png_info       *info_;
  png_info       *end_;
};

class encoder 
{
public:
  encoder(std::streambuf *, info const &);
  ~encoder();
  void encode(unsigned char const *, unsigned long size);

private:
  std::streambuf *output_;
  png_struct     *png_;
  png_info       *info_; 
};

} // namespace vsip_csl::png
} // namespace vsip_csl

#endif
