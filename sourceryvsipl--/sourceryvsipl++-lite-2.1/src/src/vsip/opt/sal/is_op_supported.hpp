/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/is_op_supported.hpp
    @author  Jules Bergmann
    @date    2006-10-26
    @brief   VSIPL++ Library: Mercury SAL ops supported for dispatch.
*/

#ifndef VSIP_OPT_SAL_IS_OP_SUPPORTED_HPP
#define VSIP_OPT_SAL_IS_OP_SUPPORTED_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/view_cast.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace sal
{

/// Traits class to help determine when SAL supports a given
/// binary operation.
///
/// Requirements:
///   OPERATOR is the Operator template class for the operation
///      (from the vsip::impl::op namespace).
///   LTYPE is the type of the left operand.
///   RTYPE is the type of the right operand.
///   RTYPE is the type of the result.
///
/// For LTYPE, RTYPE, and DSTTYPE, vector operands should be represented
/// by a pointer type (for example: float*, complex<float>*, and
/// std::pair<float*, float*>).  scalar operand should be represented
/// by a value type, such as float, complex<float>, std::pair<float, float>.

template <template <typename> class Operator,
	  typename                  SrcType,
	  typename                  DstType>
struct Is_op1_supported
{
  static bool const value = false;
};

template <template <typename, typename> class Operator,
	  typename                            LType,
	  typename                            RType,
	  typename                            DstType>
struct Is_op2_supported
{
  static bool const value = false;
};

template <template <typename, typename, typename> class Operator,
	  typename                            Type1,
	  typename                            Type2,
	  typename                            Type3,
	  typename                            DstType>
struct Is_op3_supported
{
  static bool const value = false;
};


// Tokens for ops not mapping directly to functors.

template <typename> struct copy_token;

template <typename, typename> struct veq_token;
template <typename, typename> struct vne_token;
template <typename, typename> struct vgt_token;
template <typename, typename> struct vge_token;
template <typename, typename> struct vlt_token;
template <typename, typename> struct vle_token;

template <typename, typename, typename> struct cma_token;


#define VSIP_IMPL_OP1SUP(Op, T1, DT)					\
  template <> struct Is_op1_supported<Op, T1, DT >			\
  { static bool const value = true; }

#define VSIP_IMPL_OP2SUP(Op, LT, RT, DT)				\
  template <> struct Is_op2_supported<Op, LT, RT, DT >			\
  { static bool const value = true; }

#define VSIP_IMPL_OP3SUP(Op, T1, T2, T3, DT)				\
  template <> struct Is_op3_supported<Op, T1, T2, T3, DT >		\
  { static bool const value = true; }

typedef std::pair<float*, float*>   split_float;
typedef std::pair<double*, double*> split_double;


/***********************************************************************
  Unary operators and functions provided by SAL
***********************************************************************/

VSIP_IMPL_OP1SUP(magsq_functor, complex<float>*,  float*);
VSIP_IMPL_OP1SUP(magsq_functor, split_float,      float*);

VSIP_IMPL_OP1SUP(op::Minus,     int*,             int*);
VSIP_IMPL_OP1SUP(op::Minus,     float*,           float*);
VSIP_IMPL_OP1SUP(op::Minus,     complex<float>*,  complex<float>*);
VSIP_IMPL_OP1SUP(op::Minus,     split_float,      split_float);

VSIP_IMPL_OP1SUP(mag_functor,   int*,             int*);
VSIP_IMPL_OP1SUP(mag_functor,   float*,           float*);
VSIP_IMPL_OP1SUP(mag_functor,   complex<float>*,  float*);
VSIP_IMPL_OP1SUP(mag_functor,   split_float,      float*);

VSIP_IMPL_OP1SUP(cos_functor,   float*,           float*);

VSIP_IMPL_OP1SUP(sin_functor,   float*,           float*);

VSIP_IMPL_OP1SUP(tan_functor,   float*,           float*);

VSIP_IMPL_OP1SUP(atan_functor,  float*,           float*);

VSIP_IMPL_OP1SUP(log_functor,   float*,           float*);

VSIP_IMPL_OP1SUP(log10_functor, float*,           float*);

VSIP_IMPL_OP1SUP(exp_functor,   float*,           float*);

VSIP_IMPL_OP1SUP(exp10_functor, float*,           float*);

VSIP_IMPL_OP1SUP(sqrt_functor,  float*,           float*);
VSIP_IMPL_OP1SUP(sqrt_functor,  complex<float>*,  complex<float>*);
VSIP_IMPL_OP1SUP(sqrt_functor,  split_float,      split_float);

VSIP_IMPL_OP1SUP(rsqrt_functor, float*,           float*);

VSIP_IMPL_OP1SUP(sq_functor,    float*,           float*);

VSIP_IMPL_OP1SUP(recip_functor, float*,           float*);
VSIP_IMPL_OP1SUP(recip_functor, complex<float>*,  complex<float>*);
VSIP_IMPL_OP1SUP(recip_functor, split_float,      split_float);

VSIP_IMPL_OP1SUP(copy_token,    float*,           float*);
VSIP_IMPL_OP1SUP(copy_token,    complex<float>*,  complex<float>*);
VSIP_IMPL_OP1SUP(copy_token,    split_float,      split_float);

VSIP_IMPL_OP1SUP(copy_token,    split_float,      complex<float>*);
VSIP_IMPL_OP1SUP(copy_token,    complex<float>*,  split_float);

VSIP_IMPL_OP1SUP(copy_token,    unsigned long*,   float*);
VSIP_IMPL_OP1SUP(copy_token,    unsigned short*,  float*);

VSIP_IMPL_OP1SUP(Cast_closure<long          >::Cast, float*, long*);
VSIP_IMPL_OP1SUP(Cast_closure<short         >::Cast, float*, short*);
VSIP_IMPL_OP1SUP(Cast_closure<char          >::Cast, float*, char*);
VSIP_IMPL_OP1SUP(Cast_closure<unsigned long >::Cast, float*, unsigned long*);
VSIP_IMPL_OP1SUP(Cast_closure<unsigned short>::Cast, float*, unsigned short*);
VSIP_IMPL_OP1SUP(Cast_closure<unsigned char >::Cast, float*, unsigned char*);

VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, long*, float*);
VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, short*, float*);
VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, char*, float*);
VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, unsigned long*, float*);
VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, unsigned short*, float*);
VSIP_IMPL_OP1SUP(Cast_closure<float>::Cast, unsigned char*, float*);


#if VSIP_IMPL_HAVE_SAL_DOUBLE
VSIP_IMPL_OP1SUP(magsq_functor, complex<double>*, double*);
VSIP_IMPL_OP1SUP(magsq_functor, split_double,     double*);

VSIP_IMPL_OP1SUP(op::Minus,     double*,          double*);
VSIP_IMPL_OP1SUP(op::Minus,     complex<double>*, complex<double>*);
VSIP_IMPL_OP1SUP(op::Minus,     split_double,     split_double);

VSIP_IMPL_OP1SUP(mag_functor,   double*,          double*);

VSIP_IMPL_OP1SUP(cos_functor,   double*,          double*);

VSIP_IMPL_OP1SUP(sin_functor,   double*,          double*);

VSIP_IMPL_OP1SUP(atan_functor,  double*,          double*);

VSIP_IMPL_OP1SUP(log_functor,   double*,          double*);

VSIP_IMPL_OP1SUP(log10_functor, double*,          double*);

VSIP_IMPL_OP1SUP(exp_functor,   double*,          double*);

VSIP_IMPL_OP1SUP(exp10_functor, double*,          double*);

VSIP_IMPL_OP1SUP(sqrt_functor,  double*,          double*);

VSIP_IMPL_OP1SUP(rsqrt_functor, double*,          double*);

VSIP_IMPL_OP1SUP(sq_functor,    double*,          double*);

// recip: no scalar double
VSIP_IMPL_OP1SUP(recip_functor, complex<double>*, complex<double>*);
VSIP_IMPL_OP1SUP(recip_functor, split_double,     split_double);

VSIP_IMPL_OP1SUP(copy_token,    complex<double>*, complex<double>*);
VSIP_IMPL_OP1SUP(copy_token,    split_double,     split_double);

VSIP_IMPL_OP1SUP(copy_token,    split_double,     complex<double>*);
VSIP_IMPL_OP1SUP(copy_token,    complex<double>*, split_double);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE


/***********************************************************************
  Binary operators and functions provided by SAL
***********************************************************************/

// straight-up vector add
VSIP_IMPL_OP2SUP(op::Add, int*,             int*,            int*);
VSIP_IMPL_OP2SUP(op::Add, float*,           float*,          float*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>*,  complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, split_float,      split_float,     split_float);

VSIP_IMPL_OP2SUP(op::Add, float*,           complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>*,  float*,          complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, float*,           split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Add, split_float,      float*,          split_float);

// scalar-vector vector add
VSIP_IMPL_OP2SUP(op::Add, int,              int*,            int*);
VSIP_IMPL_OP2SUP(op::Add, int*,             int,             int*);
VSIP_IMPL_OP2SUP(op::Add, float,            float*,          float*);
VSIP_IMPL_OP2SUP(op::Add, float*,           float,           float*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>,   complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>*,  complex<float>,  complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>,   split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Add, split_float,      complex<float>,  split_float);

VSIP_IMPL_OP2SUP(op::Add, float,            complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Add, complex<float>*,  float,           complex<float>*);

VSIP_IMPL_OP2SUP(op::Add, float,            split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Add, split_float,      float,           split_float);


// straight-up vector sub
VSIP_IMPL_OP2SUP(op::Sub, int*,             int*,            int*);
VSIP_IMPL_OP2SUP(op::Sub, float*,           float*,          float*);
VSIP_IMPL_OP2SUP(op::Sub, complex<float>*,  complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Sub, split_float,      split_float,     split_float);

VSIP_IMPL_OP2SUP(op::Sub, complex<float>*,  float*,          complex<float>*);
VSIP_IMPL_OP2SUP(op::Sub, split_float,      float*,          split_float);

// scalar-vector vector sub
VSIP_IMPL_OP2SUP(op::Sub, int*,             int,             int*);
// not in sal   (op::Sub, float,            float*,          float*);
VSIP_IMPL_OP2SUP(op::Sub, float*,           float,           float*);
// not in sal   (op::Sub, complex<float>,   complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Sub, complex<float>*,  complex<float>,  complex<float>*);
// not in sal   (op::Sub, complex<float>,   split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Sub, split_float,      complex<float>,  split_float);

VSIP_IMPL_OP2SUP(op::Sub, complex<float>*,  float,           complex<float>*);
VSIP_IMPL_OP2SUP(op::Sub, split_float,      float,           split_float);


// straight-up vector multiply
VSIP_IMPL_OP2SUP(op::Mult, int*,            int*,            int*);
VSIP_IMPL_OP2SUP(op::Mult, float*,          float*,          float*);
VSIP_IMPL_OP2SUP(op::Mult, complex<float>*, complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, split_float,     split_float,     split_float);

// real-complex vector multiply
VSIP_IMPL_OP2SUP(op::Mult, complex<float>*, float*,          complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, float*,          complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, split_float,     float*,          split_float);
VSIP_IMPL_OP2SUP(op::Mult, float*,          split_float,     split_float);

// scalar-vector vector multiply
VSIP_IMPL_OP2SUP(op::Mult, int,             int*,            int*);
VSIP_IMPL_OP2SUP(op::Mult, int*,            int,             int*);
VSIP_IMPL_OP2SUP(op::Mult, float,           float*,          float*);
VSIP_IMPL_OP2SUP(op::Mult, float*,          float,           float*);
VSIP_IMPL_OP2SUP(op::Mult, complex<float>,  complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<float>*, complex<float>,  complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<float>,  split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Mult, split_float,     complex<float>,  split_float);

VSIP_IMPL_OP2SUP(op::Mult, float,           complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<float>*, float,           complex<float>*);

VSIP_IMPL_OP2SUP(op::Mult, float,           split_float,     split_float);
VSIP_IMPL_OP2SUP(op::Mult, split_float,     float,           split_float);



// straight-up vector division
VSIP_IMPL_OP2SUP(op::Div, int*,             int*,            int*);
VSIP_IMPL_OP2SUP(op::Div, float*,           float*,          float*);
VSIP_IMPL_OP2SUP(op::Div, complex<float>*,  complex<float>*, complex<float>*);
VSIP_IMPL_OP2SUP(op::Div, split_float,      split_float,     split_float);

VSIP_IMPL_OP2SUP(op::Div, complex<float>*,  float*,          complex<float>*);
VSIP_IMPL_OP2SUP(op::Div, split_float,      float*,          split_float);

// scalar-vector vector division
// not in sal  (op::Div, int,             int*,            int*);
#if VSIP_IMPL_SAL_HAS_VSDIVIX
VSIP_IMPL_OP2SUP(op::Div, int*,            int,             int*);
#endif
VSIP_IMPL_OP2SUP(op::Div, float,           float*,          float*);
VSIP_IMPL_OP2SUP(op::Div, float*,          float,           float*);
// not in sal   (op::Div, complex<float>,  complex<float>*, complex<float>*);
// not in sal   (op::Div, complex<float>*, complex<float>,  complex<float>*);


// Logical

VSIP_IMPL_OP2SUP(band_functor, int*,             int*,            int*);
VSIP_IMPL_OP2SUP(bor_functor,  int*,             int*,            int*);


// vector comparisons

VSIP_IMPL_OP2SUP(max_functor, float*,             float*,            float*);

VSIP_IMPL_OP2SUP(min_functor, float*,             float*,            float*);


// Vector comparisons to 1/0
VSIP_IMPL_OP2SUP(veq_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(veq_token, int*,    int*,    int*);
VSIP_IMPL_OP2SUP(vne_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(vne_token, int*,    int*,    int*);
VSIP_IMPL_OP2SUP(vgt_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(vgt_token, int*,    int*,    int*);
VSIP_IMPL_OP2SUP(vge_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(vge_token, int*,    int*,    int*);
VSIP_IMPL_OP2SUP(vlt_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(vlt_token, int*,    int*,    int*);
VSIP_IMPL_OP2SUP(vle_token, float*,  float*,  float*);
VSIP_IMPL_OP2SUP(vle_token, int*,    int*,    int*);


#if VSIP_IMPL_HAVE_SAL_DOUBLE
// straight-up vector add
VSIP_IMPL_OP2SUP(op::Add, double*,          double*,         double*);
VSIP_IMPL_OP2SUP(op::Add, complex<double>*, complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Add, split_double,     split_double,    split_double);

// no crvadddx in SAL

// scalar-vector vector add
VSIP_IMPL_OP2SUP(op::Add, double,           double*,         double*);
VSIP_IMPL_OP2SUP(op::Add, double*,          double,          double*);
VSIP_IMPL_OP2SUP(op::Add, complex<double>,  complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Add, complex<double>*, complex<double>, complex<double>*);
VSIP_IMPL_OP2SUP(op::Add, complex<double>,  split_double,    split_double);
VSIP_IMPL_OP2SUP(op::Add, split_double,     complex<double>, split_double);

VSIP_IMPL_OP2SUP(op::Add, double,           complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Add, complex<double>*, double,          complex<double>*);

VSIP_IMPL_OP2SUP(op::Add, double,           split_double,    split_double);
VSIP_IMPL_OP2SUP(op::Add, split_double,     double,          split_double);


// straight-up vector sub
VSIP_IMPL_OP2SUP(op::Sub, double*,          double*,         double*);
VSIP_IMPL_OP2SUP(op::Sub, complex<double>*, complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Sub, split_double,     split_double,    split_double);

// scalar-vector vector sub
// not in sal   (op::Sub, double,           double*,         double*);
VSIP_IMPL_OP2SUP(op::Sub, double*,          double,          double*);
// not in sal   (op::Sub, complex<double>,  complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Sub, complex<double>*, complex<double>, complex<double>*);
// not in sal   (op::Sub, complex<double>,  split_double,    split_double);
VSIP_IMPL_OP2SUP(op::Sub, split_double,     complex<double>, split_double);

VSIP_IMPL_OP2SUP(op::Sub, complex<double>*, double,          complex<double>*);
VSIP_IMPL_OP2SUP(op::Sub, split_double,     double,          split_double);


// straight-up vector multiply
VSIP_IMPL_OP2SUP(op::Mult, double*,         double*,         double*);
VSIP_IMPL_OP2SUP(op::Mult, complex<double>*,complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Mult, split_double,    split_double,    split_double);

// real-complex vector multiply

// scalar-vector vector multiply
VSIP_IMPL_OP2SUP(op::Mult, double,          double*,         double*);
VSIP_IMPL_OP2SUP(op::Mult, double*,         double,          double*);
VSIP_IMPL_OP2SUP(op::Mult, complex<double>, complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<double>*,complex<double>, complex<double>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<double>, split_double,    split_double);
VSIP_IMPL_OP2SUP(op::Mult, split_double,    complex<double>, split_double);

VSIP_IMPL_OP2SUP(op::Mult, double,          complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Mult, complex<double>*,double,          complex<double>*);

VSIP_IMPL_OP2SUP(op::Mult, double,          split_double,    split_double);
VSIP_IMPL_OP2SUP(op::Mult, split_double,    double,          split_double);



// straight-up vector division
VSIP_IMPL_OP2SUP(op::Div, double*,          double*,         double*);
VSIP_IMPL_OP2SUP(op::Div, complex<double>*, complex<double>*,complex<double>*);
VSIP_IMPL_OP2SUP(op::Div, split_double,     split_double,    split_double);

// scalar-vector vector division
// not in sal   (op::Div, double,          double*,         double*);
VSIP_IMPL_OP2SUP(op::Div, double*,         double,          double*);
// not in sal   (op::Div, complex<double>, complex<double>*,complex<double>*);
// not in sal   (op::Div, complex<double>*,complex<double>, complex<double>*);

// vector min/max
VSIP_IMPL_OP2SUP(max_functor, double*,            double*,           double*);
VSIP_IMPL_OP2SUP(min_functor, double*,            double*,           double*);


// Vector comparisons to 1/0
VSIP_IMPL_OP2SUP(veq_token, double*, double*, double*);
VSIP_IMPL_OP2SUP(vne_token, double*, double*, double*);
VSIP_IMPL_OP2SUP(vgt_token, double*, double*, double*);
VSIP_IMPL_OP2SUP(vge_token, double*, double*, double*);
VSIP_IMPL_OP2SUP(vlt_token, double*, double*, double*);
VSIP_IMPL_OP2SUP(vle_token, double*, double*, double*);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE



/***********************************************************************
  Ternary operators and functions provided by SAL.
***********************************************************************/

// Multiply-add

VSIP_IMPL_OP3SUP(ma_functor, float,   float*,  float*, float*);
VSIP_IMPL_OP3SUP(ma_functor, float*,  float,   float*, float*);
VSIP_IMPL_OP3SUP(ma_functor, float*,  float*,  float,  float*);
VSIP_IMPL_OP3SUP(ma_functor, float*,  float*,  float*, float*);

VSIP_IMPL_OP3SUP(ma_functor, complex<float>*, complex<float>, complex<float>*,
		 complex<float>*);
VSIP_IMPL_OP3SUP(ma_functor, complex<float>, complex<float>*, complex<float>*,
		 complex<float>*);
VSIP_IMPL_OP3SUP(ma_functor, split_float, complex<float>, split_float,
		 split_float);
VSIP_IMPL_OP3SUP(ma_functor, complex<float>, split_float, split_float,
		 split_float);

#if VSIP_IMPL_HAVE_SAL_DOUBLE
VSIP_IMPL_OP3SUP(ma_functor, double,   double*,  double*, double*);
VSIP_IMPL_OP3SUP(ma_functor, double*,  double,   double*, double*);
VSIP_IMPL_OP3SUP(ma_functor, double*,  double*,  double,  double*);
VSIP_IMPL_OP3SUP(ma_functor, double*,  double*,  double*, double*);

VSIP_IMPL_OP3SUP(ma_functor,
		 complex<double>*, complex<double>, complex<double>*,
		 complex<double>*);
VSIP_IMPL_OP3SUP(ma_functor,
		 complex<double>, complex<double>*, complex<double>*,
		 complex<double>*);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE


// Multiply-subtract

VSIP_IMPL_OP3SUP(msb_functor, float*,  float,   float*, float*);
// not in sal   (msb_functor, float*,  float*,  float,  float*);
VSIP_IMPL_OP3SUP(msb_functor, float*,  float*,  float*, float*);

#if VSIP_IMPL_HAVE_SAL_DOUBLE
VSIP_IMPL_OP3SUP(msb_functor, double*,  double,   double*, double*);
// not in sal   (msb_functor, double*,  double*,  double,  double*);
VSIP_IMPL_OP3SUP(msb_functor, double*,  double*,  double*, double*);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE

// no complex msb in SAL


// Add-multiply

// not in SAL   (am_functor, float,   float*,  float*, float*);
// not in SAL   (am_functor, float*,  float,   float*, float*);
VSIP_IMPL_OP3SUP(am_functor, float*,  float*,  float,  float*);
VSIP_IMPL_OP3SUP(am_functor, float*,  float*,  float*, float*);

#if VSIP_IMPL_HAVE_SAL_DOUBLE
// not in SAL   (am_functor, double,  double*, double*,double*);
// not in SAL   (am_functor, double*, double,  double*,double*);
// not in SAL   (am_functor, double*, double*, double, double*);
VSIP_IMPL_OP3SUP(am_functor, double*, double*, double*,double*);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE


// Subtract-multiply
VSIP_IMPL_OP3SUP(sbm_functor, float*,  float*,  float,  float*);
VSIP_IMPL_OP3SUP(sbm_functor, float*,  float*,  float*, float*);
#if VSIP_IMPL_HAVE_SAL_DOUBLE
VSIP_IMPL_OP3SUP(sbm_functor, double*, double*, double*,double*);
#endif // VSIP_IMPL_HAVE_SAL_DOUBLE


// Conjugate(multiply)-add

VSIP_IMPL_OP3SUP(cma_token, complex<float>*, complex<float>*, complex<float>*,
		 complex<float>*);
VSIP_IMPL_OP3SUP(cma_token, split_float*, split_float*, split_float*,
		 split_float*);



#undef VSIP_IMPL_OP1SUP
#undef VSIP_IMPL_OP2SUP
#undef VSIP_IMPL_OP3SUP

} // namespace vsip::impl::sal
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_SAL_IS_OP_SUPPORTED_HPP
