/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    scripting/math.cpp
    @author  Stefan Seefeld
    @date    2006-09-20
    @brief   VSIPL++ Library: Python bindings for math module.

*/
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include <vsip/math.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include "types.hpp"

namespace bpl = boost::python;

namespace
{
template <typename T> T cos(T a) { return vsip::cos(a);}
template <typename T> T mag(T a) { return vsip::mag(a);}
}

#define DEF_UNARY_FUNCTION(name, func)\
  bpl::def(name, func<double>);       \
  bpl::def(name, func<Vector>);       \
  bpl::def(name, func<Matrix>);


BOOST_PYTHON_MODULE(math)
{
  DEF_UNARY_FUNCTION("cos", cos)
  DEF_UNARY_FUNCTION("mag", mag)
}
