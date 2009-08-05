/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    scripting/_base.cpp
    @author  Stefan Seefeld
    @date    2006-09-20
    @brief   VSIPL++ Library: Python bindings for common functionality.

*/
#include <boost/python.hpp>
#include <boost/noncopyable.hpp>
#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include "types.hpp"
#include <stdexcept>

namespace bpl = boost::python;

namespace
{
std::auto_ptr<Vector> construct_vector(bpl::list l)
{
  long len = bpl::extract<long>(l.attr("__len__")());
  std::auto_ptr<Vector> vector(new Vector(len));
  for (long i = 0; i != len; ++i) vector->put(i, bpl::extract<double>(l[i]));
  return vector;
}

void assign_to_vector(Vector &vector, bpl::list l)
{
  long len = bpl::extract<long>(l.attr("__len__")());
  if (vector.size() != len)
    throw std::runtime_error("Attempting to assign vectors of incompatible sizes.");
  for (long i = 0; i != len; ++i) vector.put(i, bpl::extract<double>(l[i]));
}
}

BOOST_PYTHON_MODULE(_base)
{
  // Disambiguate from overload set.
  double (Vector::*v_get)(vsip::index_type) const = &Vector::get;
  vsip::length_type (Matrix::*m_size)(vsip::dimension_type) const 
    = &Matrix::size;
  double (Matrix::*m_get)(vsip::index_type, vsip::index_type) const 
    = &Matrix::get;

  dcomplex (CVector::*cv_get)(vsip::index_type) const = &CVector::get;
  vsip::length_type (CMatrix::*cm_size)(vsip::dimension_type) const
    = &CMatrix::size;
  dcomplex (CMatrix::*cm_get)(vsip::index_type, vsip::index_type) const
    = &CMatrix::get;

  bpl::class_<vsip::vsipl, boost::noncopyable> vsipl("_library");

  bpl::class_<Vector> vector("vector", bpl::init<vsip::length_type>());
  vector.def("__init__", bpl::make_constructor(construct_vector));
  vector.def("assign", assign_to_vector);
  vector.def("length", &Vector::length);
  vector.def("get", v_get);
  vector.def("put", &Vector::put);

  bpl::class_<CVector> cvector("cvector", bpl::init<vsip::length_type>());
  cvector.def("length", &CVector::length);
  cvector.def("get", cv_get);
  cvector.def("put", &CVector::put);

  bpl::class_<Matrix> matrix("matrix",
                             bpl::init<vsip::length_type, vsip::length_type>());
  matrix.def("size", m_size);
  matrix.def("get", m_get);
  matrix.def("put", &Matrix::put);

  bpl::class_<CMatrix> cmatrix("cmatrix",
                               bpl::init<vsip::length_type, vsip::length_type>());
  cmatrix.def("size", cm_size);
  cmatrix.def("get", cm_get);
  cmatrix.def("put", &CMatrix::put);
}
