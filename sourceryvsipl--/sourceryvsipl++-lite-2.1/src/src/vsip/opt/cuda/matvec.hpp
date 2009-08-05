/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/matvec.hpp
    @author  Don McCoy
    @date    2009-02-05
    @brief   VSIPL++ Library: CUDA-based BLAS evaluators 
*/

#ifndef VSIP_OPT_CUDA_MATVEC_HPP
#define VSIP_OPT_CUDA_MATVEC_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/dispatch.hpp>
#include <vsip/opt/cuda/bindings.hpp>
#include <vsip/opt/cuda/device_memory.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{
namespace dispatcher
{

/// CUDA evaluator for vector-vector dot-product (non-conjugated).
template <typename T,
          typename Block0,
          typename Block1>
struct Evaluator<Op_prod_vv_dot, Cuda_tag,
                 T(Block0 const&, Block1 const&)>
{
  static bool const ct_valid =
    impl::cuda::Cuda_traits<T>::valid &&
    Type_equal<T, typename Block0::value_type>::value &&
    Type_equal<T, typename Block1::value_type>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    // check that format is interleaved.
    !Is_split_block<Block0>::value &&
    !Is_split_block<Block1>::value;

  static bool rt_valid(Block0 const& a, Block1 const& b) 
  { 
    Ext_data<Block0> ext_a(const_cast<Block0&>(a));
    Ext_data<Block1> ext_b(const_cast<Block1&>(b));
    
    // check that data is unit stride
    return ((ext_a.stride(0) == 1) && (ext_b.stride(0) == 1));
  }

  static T exec(Block0 const& a, Block1 const& b)
  {
    assert(a.size(1, 0) == b.size(1, 0));

    cuda::Device_memory<Block0 const> dev_a(a);
    cuda::Device_memory<Block1 const> dev_b(b);

    T r = cuda::dot(a.size(1, 0),
                    dev_a.data(), a.impl_stride(1, 0),
                    dev_b.data(), b.impl_stride(1, 0));
    ASSERT_CUBLAS_OK();

    return r;
  }
};


/// CUDA evaluator for vector-vector dot-product (conjugated).
template <typename T,
          typename Block0,
          typename Block1>
struct Evaluator<Op_prod_vv_dot, Cuda_tag,
                 std::complex<T>(Block0 const&, 
                   Unary_expr_block<1, conj_functor, Block1, std::complex<T> > const&)>
{
  static bool const ct_valid = 
    impl::cuda::Cuda_traits<complex<T> >::valid &&
    Type_equal<complex<T>, typename Block0::value_type>::value &&
    Type_equal<complex<T>, typename Block1::value_type>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    // check that format is interleaved.
    !Is_split_block<Block0>::value &&
    !Is_split_block<Block1>::value;

  static bool rt_valid(
    Block0 const& a,
    Unary_expr_block<1, impl::conj_functor, Block1, complex<T> > const& b)
  {
    Ext_data<Block0> ext_a(const_cast<Block0&>(a));
    Ext_data<Block1> ext_b(const_cast<Block1&>(b.op()));
    
    // check that data is unit stride
    return ((ext_a.stride(0) == 1) && (ext_b.stride(0) == 1));
  }

  static complex<T> exec(
    Block0 const& a, 
    Unary_expr_block<1, impl::conj_functor, Block1, complex<T> > const& b)
  {
    assert(a.size(1, 0) == b.size(1, 0));

    cuda::Device_memory<Block0 const> dev_a(a);
    cuda::Device_memory<Block1 const> dev_b(b.op());

    complex<T> r = cuda::dotc(a.size(1, 0),
                              dev_b.data(), b.op().impl_stride(1, 0), 
                              dev_a.data(), a.impl_stride(1, 0) * 2);
    ASSERT_CUBLAS_OK();
    // Note:
    //   BLAS    cdotc(x, y)  => conj(x) * y, while 
    //   VSIPL++ cvjdot(x, y) => x * conj(y)

    return r;
  }
};


} // namespace vsip::impl::dispatcher
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_MATVEC_HPP
