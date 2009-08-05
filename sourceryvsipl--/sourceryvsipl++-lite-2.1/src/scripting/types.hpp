/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    scripting/types.hpp
    @author  Stefan Seefeld
    @date    2006-09-20
    @brief   VSIPL++ Library: Common types for all Python bindings.

*/
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>

typedef std::complex<double> dcomplex;

typedef vsip::Vector<double> Vector;
typedef vsip::Vector<dcomplex> CVector;

typedef vsip::Matrix<double> Matrix;
typedef vsip::Matrix<dcomplex> CMatrix;

