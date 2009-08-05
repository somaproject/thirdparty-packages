/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip_csl/png.cpp
    @author  Stefan Seefeld
    @date    2006-03-21
    @brief   VSIPL++ Library: 
*/
 
#include <vsip_csl/png.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <streambuf>
#include <stdexcept>
#include <cassert>

namespace vsip_csl
{
namespace png
{
namespace
{
void decoder_read(png_structp png_ptr, png_bytep image, png_size_t length) 
{
  std::streambuf *input = static_cast<std::streambuf *>(png_ptr->io_ptr);
  input->sgetn((char *)image, (size_t)length);
}

void decoder_warning(png_structp, png_const_charp msg)
{
  std::cout << "PNG Warning : " << msg << std::endl;
}

void decoder_error(png_structp, png_const_charp msg)
{
  std::cerr << "PNG Error : " << msg << std::endl;
}

void encoder_write(png_structp png_ptr, png_bytep image, png_size_t length) 
{
  std::streambuf *sbuf = static_cast<std::streambuf *>(png_ptr->io_ptr);
  sbuf->sputn((char*)image, (size_t)length);
}

void encoder_flush(png_structp png_ptr) 
{
  std::streambuf *sbuf = static_cast<std::streambuf *>(png_ptr->io_ptr);
  sbuf->pubsync();
}

void encoder_warning(png_structp, png_const_charp msg)
{
  std::cout << "PNG Warning : " << msg << std::endl;
}

void encoder_error(png_structp, png_const_charp msg)
{
  std::cerr << "PNG Error : " << msg << std::endl;
}
} // namespace <anonymous>

decoder::decoder(std::streambuf *sbuf, info &i)
  : input_(sbuf),
    png_(png_create_read_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0)),
    info_(png_create_info_struct(png_)),
    end_(png_create_info_struct(png_))
{
  png_byte header[magic_];
  input_->sgetn((char*)header, magic_);
  if (png_sig_cmp(header, 0, magic_)) throw std::runtime_error("Not a PNG file");

  png_set_sig_bytes(png_, magic_);
  png_set_read_fn(png_, input_, decoder_read);
  png_set_error_fn(png_, input_, decoder_error, decoder_warning);
  png_set_read_status_fn(png_, 0);
  png_read_info(png_, info_);
  png_uint_32 w, h;
  int d, c, in, co, f;
  png_get_IHDR(png_, info_, &w, &h, &d, &c, &in, &co, &f);  
  i.width = w;
  i.height = h;
  i.rowbytes = png_get_rowbytes(png_, info_);
  i.depth = d;
  i.colortype = static_cast<color_type>(c);
  i.compression = co;
  i.filter = f;
  i.interlace = static_cast<interlace_type>(in);
}

decoder::~decoder()
{
  png_destroy_read_struct(&png_, &info_, &end_);
}

void decoder::decode(unsigned char *data, unsigned long size)
{
  png_uint_32 height = png_get_image_height(png_, info_);
  png_uint_32 rowbytes = png_get_rowbytes(png_, info_);
  assert(size >= height * rowbytes);
  std::vector<unsigned char *> rows(height);
  for (png_uint_32 i = 0; i < height; i++)
    rows[i] = data + i * rowbytes;

  png_read_image(png_, &*rows.begin());
  png_read_end(png_, end_);
}

encoder::encoder(std::streambuf *sb, info const &i)
  : output_(sb),
    png_(png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0)),
    info_(png_create_info_struct(png_))
{
  png_set_write_fn(png_, output_, encoder_write, encoder_flush);
  png_set_error_fn(png_, output_, encoder_error, encoder_warning);
  png_set_write_status_fn(png_, 0);

  png_set_IHDR(png_, info_,
	       i.width, i.height, i.depth, i.colortype,
	       i.interlace, i.compression, i.filter);
}

encoder::~encoder()
{
  png_destroy_write_struct(&png_, &info_);
}

void encoder::encode(unsigned char const *data, unsigned long size)
{
  png_uint_32 height = png_get_image_height(png_, info_);
  png_uint_32 rowbytes = png_get_rowbytes(png_, info_);
  assert(size >= height * rowbytes);
  std::vector<unsigned char const *> rows(height);
  for (png_uint_32 i = 0; i < height; i++)
    rows[i] = data + i * rowbytes;
  png_write_info(png_, info_);
  png_write_image(png_, const_cast<unsigned char **>(&*rows.begin()));
  png_write_end(png_, 0);
}

} // namespace vsip_csl::png
} // namespace vsip_csl
