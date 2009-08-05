/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    scripting/signal.cpp
    @author  Stefan Seefeld
    @date    2006-09-20
    @brief   VSIPL++ Library: Python bindings for signal module.

*/
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include <vsip/signal.hpp>
#include "types.hpp"

namespace bpl = boost::python;

namespace
{
typedef vsip::Fft<vsip::Vector, double, dcomplex, 0, vsip::by_reference> Fft;
typedef vsip::Fft<vsip::Vector, dcomplex, dcomplex, 0, vsip::by_reference> CFft;
typedef vsip::Convolution<vsip::const_Vector,
                          vsip::nonsym, vsip::support_full, double>
  Convolution;
typedef vsip::Convolution<vsip::const_Vector,
                          vsip::nonsym, vsip::support_full, dcomplex>
  CConvolution;

CVector fft_by_ref(Fft &fft, Vector input)
{
  CVector output(fft.output_size().length());
  fft.operator()(input, output);
  return output;
}

Vector convolute(Convolution &conv, Vector input)
{
  Vector output(conv.output_size().size());
  conv.operator()(input, output);
  return output;
}

}

BOOST_PYTHON_MODULE(signal)
{
  bpl::class_<Fft, boost::noncopyable>
    fft_fwd("fft_fwd", bpl::init<vsip::length_type, double>());
  fft_fwd.def("__call__", fft_by_ref);

  bpl::class_<Convolution, boost::noncopyable>
    conv("convolution", bpl::init<Vector, long>());
  conv.def("__call__", convolute);
}
