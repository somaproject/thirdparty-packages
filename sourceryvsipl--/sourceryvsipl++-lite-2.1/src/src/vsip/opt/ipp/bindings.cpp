/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ipp/bindings.hpp
    @author  Stefan Seefeld
    @date    2005-08-10
    @brief   VSIPL++ Library: Wrappers and traits to bridge with Intel's IPP.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/math.hpp>
#include <vsip/signal.hpp>
#include <vsip/opt/ipp/bindings.hpp>
#include <ipps.h>
#include <ippi.h>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace ipp
{


// Generate a vector from scratch
#define VSIP_IMPL_IPP_VGEN0(FCN, T, IPPFCN, IPPT)			\
void									\
FCN(									\
  T*       Z,								\
  length_type len)							\
{									\
  VSIP_IMPL_COVER_FCN("IPP_VGEN0", IPPFCN)				\
  if (len > 0)								\
  {									\
    IppStatus status = IPPFCN(						\
      reinterpret_cast<IPPT*>(Z),					\
      static_cast<int>(len));						\
    assert(status == ippStsNoErr);					\
  }									\
}



#define VSIP_IMPL_IPP_V(FCN, T, IPPFCN, IPPT)				\
void									\
FCN(									\
  T const* A,								\
  T*       Z,								\
  length_type len)							\
{									\
  VSIP_IMPL_COVER_FCN("IPP_V", IPPFCN)					\
  if (len > 0)								\
  {									\
    IppStatus status = IPPFCN(						\
      reinterpret_cast<IPPT const*>(A),					\
      reinterpret_cast<IPPT*>(Z),					\
      static_cast<int>(len));						\
    assert(status == ippStsNoErr);					\
  }									\
}

// Complex->Real unary function.

#define VSIP_IMPL_IPP_V_CR(FCN, T, IPPFCN, IPPCT, IPPT)			\
void									\
FCN(									\
  complex<T> const* A,							\
  T*                Z,							\
  length_type len)							\
{									\
  VSIP_IMPL_COVER_FCN("IPP_V_CR", IPPFCN)				\
  if (len > 0)								\
  {									\
    IppStatus status = IPPFCN(						\
      reinterpret_cast<IPPCT const*>(A),				\
      reinterpret_cast<IPPT*>(Z),					\
      static_cast<int>(len));						\
    assert(status == ippStsNoErr);					\
  }									\
}

#define VSIP_IMPL_IPP_VV(FCN, T, IPPFCN, IPPT)				\
void									\
FCN(									\
  T const* A,								\
  T const* B,								\
  T*       Z,								\
  length_type len)							\
{									\
  if (len > 0)								\
  {									\
    IppStatus status = IPPFCN(						\
      reinterpret_cast<IPPT const*>(A),					\
      reinterpret_cast<IPPT const*>(B),					\
      reinterpret_cast<IPPT*>(Z),					\
      static_cast<int>(len));						\
    assert(status == ippStsNoErr);					\
  }									\
}

#define VSIP_IMPL_IPP_VV_R(FCN, T, IPPFCN, IPPT)			\
void									\
FCN(									\
  T const* A,								\
  T const* B,								\
  T*       Z,								\
  length_type len)							\
{									\
  if (len > 0)								\
  {									\
    IppStatus status = IPPFCN(						\
      reinterpret_cast<IPPT const*>(B),					\
      reinterpret_cast<IPPT const*>(A),					\
      reinterpret_cast<IPPT*>(Z),					\
      static_cast<int>(len));						\
    assert(status == ippStsNoErr);					\
  }									\
}

// Copy
VSIP_IMPL_IPP_V(vcopy, float,           ippsCopy_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vcopy, double,          ippsCopy_64f,  Ipp64f)
VSIP_IMPL_IPP_V(vcopy, complex<float>,  ippsCopy_32fc, Ipp32fc)
VSIP_IMPL_IPP_V(vcopy, complex<double>, ippsCopy_64fc, Ipp64fc)

// Magnitude (aka abs)
// VSIP_IMPL_IPP_V(vmag, int16,           ippsAbs_16s,  Ipp16s)
// VSIP_IMPL_IPP_V(vmag, int32,           ippsAbs_32s,  Ipp32s)
VSIP_IMPL_IPP_V(vmag, float,           ippsAbs_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vmag, double,          ippsAbs_64f,  Ipp64f)
VSIP_IMPL_IPP_V_CR(vmag, float,  ippsMagnitude_32fc, Ipp32fc, Ipp32f)
VSIP_IMPL_IPP_V_CR(vmag, double, ippsMagnitude_64fc, Ipp64fc, Ipp64f)

// Mag-sq
VSIP_IMPL_IPP_V(vmagsq,    float,  ippsSqr_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vmagsq,    double, ippsSqr_64f,  Ipp64f)
// IPP 4.x does not have magsq for complex floating-point

// Square
VSIP_IMPL_IPP_V(vsq,  float,           ippsSqr_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vsq,  double,          ippsSqr_64f,  Ipp64f)
VSIP_IMPL_IPP_V(vsq,  complex<float>,  ippsSqr_32fc, Ipp32fc)
VSIP_IMPL_IPP_V(vsq,  complex<double>, ippsSqr_64fc, Ipp64fc)

// Square-root
VSIP_IMPL_IPP_V(vsqrt, float,           ippsSqrt_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vsqrt, double,          ippsSqrt_64f,  Ipp64f)
VSIP_IMPL_IPP_V(vsqrt, complex<float>,  ippsSqrt_32fc, Ipp32fc)
VSIP_IMPL_IPP_V(vsqrt, complex<double>, ippsSqrt_64fc, Ipp64fc)

// Exponent
VSIP_IMPL_IPP_V(vexp, float,           ippsExp_32f,  Ipp32f)
VSIP_IMPL_IPP_V(vexp, double,          ippsExp_64f,  Ipp64f)

// Addition
VSIP_IMPL_IPP_VV(vadd, float,           ippsAdd_32f,  Ipp32f)
VSIP_IMPL_IPP_VV(vadd, double,          ippsAdd_64f,  Ipp64f)
VSIP_IMPL_IPP_VV(vadd, complex<float>,  ippsAdd_32fc, Ipp32fc)
VSIP_IMPL_IPP_VV(vadd, complex<double>, ippsAdd_64fc, Ipp64fc)


// Subtraction
// Note: IPP subtract arguments are in reverse order:

VSIP_IMPL_IPP_VV_R(vsub, float,           ippsSub_32f,  Ipp32f)
VSIP_IMPL_IPP_VV_R(vsub, double,          ippsSub_64f,  Ipp64f)
VSIP_IMPL_IPP_VV_R(vsub, complex<float>,  ippsSub_32fc, Ipp32fc)
VSIP_IMPL_IPP_VV_R(vsub, complex<double>, ippsSub_64fc, Ipp64fc)


// Multiplication

VSIP_IMPL_IPP_VV(vmul, float,           ippsMul_32f,  Ipp32f)
VSIP_IMPL_IPP_VV(vmul, double,          ippsMul_64f,  Ipp64f)
VSIP_IMPL_IPP_VV(vmul, complex<float>,  ippsMul_32fc, Ipp32fc)
VSIP_IMPL_IPP_VV(vmul, complex<double>, ippsMul_64fc, Ipp64fc)

// Division

// Note: IPP division arguments are in reverse order:
//
//   ippsDiv(X, Y, Z, ...);
//
// Computes: Z = Y / X.
//
// We swap arguments to IPP so that "A / B" corresponds to vdiv(A, B).

VSIP_IMPL_IPP_VV_R(vdiv, float,           ippsDiv_32f,  Ipp32f)
VSIP_IMPL_IPP_VV_R(vdiv, double,          ippsDiv_64f,  Ipp64f)
VSIP_IMPL_IPP_VV_R(vdiv, complex<float>,  ippsDiv_32fc, Ipp32fc)
VSIP_IMPL_IPP_VV_R(vdiv, complex<double>, ippsDiv_64fc, Ipp64fc)

// Zero
VSIP_IMPL_IPP_VGEN0(vzero, float,           ippsZero_32f,  Ipp32f)
VSIP_IMPL_IPP_VGEN0(vzero, double,          ippsZero_64f,  Ipp64f)
VSIP_IMPL_IPP_VGEN0(vzero, complex<float>,  ippsZero_32fc, Ipp32fc)
VSIP_IMPL_IPP_VGEN0(vzero, complex<double>, ippsZero_64fc, Ipp64fc)



// Scalar-vector functions

void svadd(float A, float const* B, float* Z, length_type len)
{
  if (len > 0)
  {
    IppStatus status = ippsAddC_32f(B, A, Z, static_cast<int>(len));
    assert(status == ippStsNoErr);
  }
}

void svadd(double A, double const* B, double* Z, length_type len)
{
  IppStatus status = ippsAddC_64f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svadd(std::complex<float> A, std::complex<float> const* B,
	   std::complex<float>* Z, length_type len)
{
  Ipp32fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsAddC_32fc(
    reinterpret_cast<Ipp32fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp32fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svadd(std::complex<double> A, std::complex<double> const* B,
	   std::complex<double>* Z, length_type len)
{
  Ipp64fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsAddC_64fc(
    reinterpret_cast<Ipp64fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp64fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}


/// scalar-vector subtraction: scalar - vector

void svsub(float A, float const* B, float* Z, length_type len)
{
  IppStatus status = ippsSubCRev_32f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(double A, double const* B, double* Z, length_type len)
{
  IppStatus status = ippsSubCRev_64f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(std::complex<float> A, std::complex<float> const* B,
	   std::complex<float>* Z, length_type len)
{
  Ipp32fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsSubCRev_32fc(
    reinterpret_cast<Ipp32fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp32fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(std::complex<double> A, std::complex<double> const* B,
	   std::complex<double>* Z, length_type len)
{
  Ipp64fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsSubCRev_64fc(
    reinterpret_cast<Ipp64fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp64fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}



/// scalar-view subtraction: vector - scalar

void svsub(float const* A, float B, float* Z, length_type len)
{
  IppStatus status = ippsSubC_32f(A, B, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(double const* A, double B, double* Z, length_type len)
{
  IppStatus status = ippsSubC_64f(A, B, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(std::complex<float> const* A, std::complex<float> B,
	   std::complex<float>* Z, length_type len)
{
  Ipp32fc ipp_B = { B.real(), B.imag()};
  IppStatus status = ippsSubC_32fc(
    reinterpret_cast<Ipp32fc const*>(A),
    ipp_B,
    reinterpret_cast<Ipp32fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svsub(std::complex<double> const* A, std::complex<double> B,
	   std::complex<double>* Z, length_type len)
{
  Ipp64fc ipp_B = { B.real(), B.imag()};
  IppStatus status = ippsSubC_64fc(
    reinterpret_cast<Ipp64fc const*>(A),
    ipp_B,
    reinterpret_cast<Ipp64fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}



void svmul(float A, float const* B, float* Z, length_type len)
{
  IppStatus status = ippsMulC_32f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svmul(double A, double const* B, double* Z, length_type len)
{
  IppStatus status = ippsMulC_64f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svmul(std::complex<float> A, std::complex<float> const* B,
	   std::complex<float>* Z, length_type len)
{
  Ipp32fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsMulC_32fc(
    reinterpret_cast<Ipp32fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp32fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svmul(std::complex<double> A, std::complex<double> const* B,
	   std::complex<double>* Z, length_type len)
{
  Ipp64fc ipp_A = { A.real(), A.imag()};
  IppStatus status = ippsMulC_64fc(
    reinterpret_cast<Ipp64fc const*>(B),
    ipp_A,
    reinterpret_cast<Ipp64fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}



/// scalar-vector division: scalar / vector

void svdiv(float A, float const* B, float* Z, length_type len)
{
  IppStatus status = ippsDivCRev_32f(B, A, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}



/// scalar-view division: vector / scalar

void svdiv(float const* A, float B, float* Z, length_type len)
{
  IppStatus status = ippsDivC_32f(A, B, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svdiv(double const* A, double B, double* Z, length_type len)
{
  IppStatus status = ippsDivC_64f(A, B, Z, static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svdiv(std::complex<float> const* A, std::complex<float> B,
	   std::complex<float>* Z, length_type len)
{
  Ipp32fc ipp_B = { B.real(), B.imag()};
  IppStatus status = ippsDivC_32fc(
    reinterpret_cast<Ipp32fc const*>(A),
    ipp_B,
    reinterpret_cast<Ipp32fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}

void svdiv(std::complex<double> const* A, std::complex<double> B,
	   std::complex<double>* Z, length_type len)
{
  Ipp64fc ipp_B = { B.real(), B.imag()};
  IppStatus status = ippsDivC_64fc(
    reinterpret_cast<Ipp64fc const*>(A),
    ipp_B,
    reinterpret_cast<Ipp64fc*>(Z),
    static_cast<int>(len));
  assert(status == ippStsNoErr);
}




// Convolution
void conv(float* coeff, length_type coeff_size,
	  float* in,    length_type in_size,
	  float* out)
{
  ippsConv_32f(coeff, coeff_size, in, in_size, out);
}

void conv(double* coeff, length_type coeff_size,
	  double* in,    length_type in_size,
	  double* out)
{
  ippsConv_64f(coeff, coeff_size, in, in_size, out);
}

void
conv_full_2d(
  short*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  short*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  short*      out,
  length_type out_row_stride)
{
  IppiSize coeff_size = { coeff_cols, coeff_rows };
  IppiSize in_size    = { in_cols,    in_rows };
  
  ippiConvFull_16s_C1R(
    coeff, coeff_row_stride*sizeof(short), coeff_size,
    in,    in_row_stride   *sizeof(short),    in_size,
    out,   out_row_stride  *sizeof(short),
    1);
}

void
conv_full_2d(
  float*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  float*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  float*      out,
  length_type out_row_stride)
{
  IppiSize coeff_size = { coeff_cols, coeff_rows };
  IppiSize in_size    = { in_cols,    in_rows };
  
  ippiConvFull_32f_C1R(
    coeff, coeff_row_stride*sizeof(float), coeff_size,
    in,    in_row_stride   *sizeof(float),    in_size,
    out,   out_row_stride  *sizeof(float));
}

void
conv_valid_2d(
  float*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  float*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  float*      out,
  length_type out_row_stride)
{
  IppiSize coeff_size = { coeff_cols, coeff_rows };
  IppiSize in_size    = { in_cols,    in_rows };
  
  ippiConvValid_32f_C1R(
    coeff, coeff_row_stride*sizeof(float), coeff_size,
    in,    in_row_stride   *sizeof(float),    in_size,
    out,   out_row_stride  *sizeof(float));
}

void
conv_valid_2d(
  short*      coeff,
  length_type coeff_rows,
  length_type coeff_cols,
  length_type coeff_row_stride,
  short*      in,
  length_type in_rows,
  length_type in_cols,
  length_type in_row_stride,
  short*      out,
  length_type out_row_stride)
{
  IppiSize coeff_size = { coeff_cols, coeff_rows };
  IppiSize in_size    = { in_cols,    in_rows };
  
  ippiConvValid_16s_C1R(
    coeff, coeff_row_stride*sizeof(short), coeff_size,
    in,    in_row_stride   *sizeof(short),    in_size,
    out,   out_row_stride  *sizeof(short),
    1);
}

} // namespace vsip::impl::ipp
} // namespace vsip::impl
} // namespace vsip
