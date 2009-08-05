/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/cml/matvec.hpp
    @author  Don McCoy
    @date    2008-05-07
    @brief   VSIPL++ Library: CML matrix product evaluators.
*/

#ifndef VSIP_OPT_CBE_CML_MATVEC_HPP
#define VSIP_OPT_CBE_CML_MATVEC_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/opt/dispatch.hpp>
#include <vsip/opt/cbe/cml/prod.hpp>
#include <vsip/opt/cbe/cml/traits.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace dispatcher 
{

/// CML evaluator for vector dot products (non-conjugated)
template <typename T,
          typename Block0,
          typename Block1>
struct Evaluator<Op_prod_vv_dot, Cml_tag,
                 T(Block0 const&, Block1 const&)>
{
  typedef typename Block_layout<Block0>::complex_type complex1_type;
  typedef typename Block_layout<Block1>::complex_type complex2_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    // check that all data types are equal
    Type_equal<T, typename Block0::value_type>::value &&
    Type_equal<T, typename Block1::value_type>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    // check complex layout is consistent
    Is_split_block<Block0>::value == Is_split_block<Block1>::value;

  static bool rt_valid(Block0 const&, Block1 const&) { return true; }

  static T exec(Block0 const& a, Block1 const& b)
  {
    assert(a.size(1, 0) == b.size(1, 0));

    Ext_data<Block0> ext_a(const_cast<Block0&>(a));
    Ext_data<Block1> ext_b(const_cast<Block1&>(b));

    T r = T();
    cml::dot( ext_a.data(), ext_a.stride(0),
              ext_b.data(), ext_b.stride(0),
              &r, a.size(1, 0) );

    return r;
  }
};


/// CML evaluator for vector dot products (conjugated)
template <typename T,
          typename Block0,
          typename Block1>
struct Evaluator<Op_prod_vv_dot, Cml_tag,
                 std::complex<T>(Block0 const&, 
                   Unary_expr_block<1, conj_functor, Block1, std::complex<T> > const&)>
{
  typedef typename Block_layout<Block0>::complex_type complex1_type;
  typedef typename Block_layout<Block1>::complex_type complex2_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    // check that types are complex
    Is_complex<typename Block0::value_type>::value &&
    Is_complex<typename Block1::value_type>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    // check complex layout is consistent
    Is_split_block<Block0>::value == Is_split_block<Block1>::value;

  static bool rt_valid(
    Block0 const&, 
    Unary_expr_block<1, impl::conj_functor, Block1, complex<T> > const&)
  { return true; }

  static complex<T> exec(
    Block0 const& a, 
    Unary_expr_block<1, impl::conj_functor, Block1, complex<T> > const& b)
  {
    assert(a.size(1, 0) == b.size(1, 0));

    Ext_data<Block0> ext_a(const_cast<Block0&>(a));
    Ext_data<Block1> ext_b(const_cast<Block1&>(b.op()));

    complex<T> r = complex<T>();
    cml::dotc( ext_a.data(), ext_a.stride(0),
               ext_b.data(), ext_b.stride(0),
               &r, a.size(1, 0) );

    return r;
  }
};


/// CML evaluator for outer products
template <typename T1,
          typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_vv_outer, Cml_tag,
                 void(Block0&, T1, Block1 const&, Block2 const&)>
{
  typedef typename Block_layout<Block0>::order_type order0_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that the output is row-major 
    Type_equal<order0_type, row2_type>::value &&
    // check that all data types are equal
    Type_equal<T1, typename Block0::value_type>::value &&
    Type_equal<T1, typename Block1::value_type>::value &&
    Type_equal<T1, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& r, T1, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    bool unit_stride =
      (ext_r.stride(1) == 1) &&
      (ext_a.stride(0) == 1) && 
      (ext_b.stride(0) == 1);

    return unit_stride;
  }

  static void exec(Block0& r, T1 alpha, Block1 const& a, Block2 const& b)
  {
    assert(a.size(1, 0) == r.size(2, 0));
    assert(b.size(1, 0) == r.size(2, 1));

    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // CML does not support a scaling parameter, so it is built into the
    // wrapper function.
    cml::outer( alpha, 
                ext_a.data(), ext_a.stride(0),
                ext_b.data(), ext_b.stride(0),
                ext_r.data(), ext_r.stride(0),
                a.size(1, 0), b.size(1, 0) );
  }
};


template <typename T1,
          typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_vv_outer, Cml_tag,
                 void(Block0&, std::complex<T1>, Block1 const&, Block2 const&)>
{
  typedef typename Block_layout<Block0>::order_type order0_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that the output is row-major 
    Type_equal<order0_type, row2_type>::value &&
    // check that all data types are equal
    Type_equal<std::complex<T1>, typename Block0::value_type>::value &&
    Type_equal<std::complex<T1>, typename Block1::value_type>::value &&
    Type_equal<std::complex<T1>, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& r, std::complex<T1>, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    bool unit_stride =
      (ext_r.stride(1) == 1) &&
      (ext_a.stride(0) == 1) && 
      (ext_b.stride(0) == 1);

    return unit_stride;
  }

  static void exec(Block0& r, std::complex<T1> alpha, Block1 const& a, Block2 const& b)
  {
    assert(a.size(1, 0) == r.size(2, 0));
    assert(b.size(1, 0) == r.size(2, 1));

    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // CML does not support a scaling parameter, so it is built into the
    // wrapper function.
    cml::outer( alpha, 
                ext_a.data(), ext_a.stride(0),
                ext_b.data(), ext_b.stride(0),
                ext_r.data(), ext_r.stride(0),
                a.size(1, 0), b.size(1, 0) );
  }
};


/// CML evaluator for matrix-vector products
template <typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_mv, Cml_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  typedef typename Block0::value_type T;
  typedef typename Block_layout<Block1>::order_type order1_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that all data types are equal
    Type_equal<T, typename Block1::value_type>::value &&
    Type_equal<T, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& /*r*/, Block1 const& a, Block2 const& /*b*/)
  {
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));

    // For 'a', the dimension with the smallest stride must be one,
    // which depends on whether it is row- or column-major.
    bool is_a_row = Type_equal<order1_type, row2_type>::value;
    stride_type a_stride = is_a_row ? ext_a.stride(1) : ext_a.stride(0);

    return (a_stride == 1);
  }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // Either row- or column-major layouts are supported for the input 
    // matrix by using the identity:
    //   trans(r) = trans(b) * trans(a)
    // or just
    //   r = b * trans(a)  (since r and b are vectors)
    if (Type_equal<order1_type, row2_type>::value)
    {
      cml::mvprod(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(0),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 0),   // M
        a.size(2, 1) ); // N
    }
    else if (Type_equal<order1_type, col2_type>::value)
    {
      cml::vmprod(
        ext_b.data(), ext_b.stride(0),
        ext_a.data(), ext_a.stride(1),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 1),   // N
        a.size(2, 0) ); // M
    }
    else
      assert(0);
  }
};


/// CML evaluator for vector-matrix products
template <typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_vm, Cml_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  typedef typename Block0::value_type T;
  typedef typename Block_layout<Block2>::order_type order2_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that all data types are equal
    Type_equal<T, typename Block1::value_type>::value &&
    Type_equal<T, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& /*r*/, Block1 const& /*a*/, Block2 const& b)
  {
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // For 'b', the dimension with the smallest stride must be one,
    // which depends on whether it is row- or column-major.
    bool is_b_row = Type_equal<order2_type, row2_type>::value;
    stride_type b_stride = is_b_row ? ext_b.stride(1) : ext_b.stride(0);

    return (b_stride == 1);
  }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // Either row- or column-major layouts are supported for the input 
    // matrix by using the identity:
    //   trans(r) = trans(b) * trans(a)
    // or just
    //   r = b * trans(a)  (since r and b are vectors)
    if (Type_equal<order2_type, row2_type>::value)
    {
      cml::vmprod(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(0),
        ext_r.data(), ext_r.stride(0),
        b.size(2, 0),   // M
        b.size(2, 1) ); // N
    }
    else if (Type_equal<order2_type, col2_type>::value)
    {
      cml::mvprod(
        ext_b.data(), ext_b.stride(1),
        ext_a.data(), ext_a.stride(0),
        ext_r.data(), ext_r.stride(0),
        b.size(2, 1),   // N
        b.size(2, 0) ); // M
    }
    else
      assert(0);
  }
};


/// CML evaluator for matrix-matrix products
template <typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_mm, Cml_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  typedef typename Block0::value_type T;
  typedef typename Block_layout<Block0>::order_type order0_type;
  typedef typename Block_layout<Block1>::order_type order1_type;
  typedef typename Block_layout<Block2>::order_type order2_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that all data types are equal
    Type_equal<T, typename Block1::value_type>::value &&
    Type_equal<T, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that the layout is row-major for the first input and the output
    Type_equal<order0_type, row2_type>::value && 
    Type_equal<order1_type, row2_type>::value && 
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // For 'b', the dimension with the smallest stride must be one,
    // which depends on whether it is row- or column-major.
    bool is_b_row = Type_equal<order2_type, row2_type>::value;
    stride_type b_stride = is_b_row ? ext_b.stride(1) : ext_b.stride(0);

    return
      // ensure the data is unit-stide
      ( ext_r.stride(1) == 1 &&
        ext_a.stride(1) == 1 &&
        b_stride        == 1 );
  }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // Either row- or column-major layouts are supported for
    // the second input by mapping them to the normal product
    // or transpose product respectively.
    if (Type_equal<order2_type, row2_type>::value)
    {
      cml::mprod(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(0),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 0),   // M
        a.size(2, 1),   // N
        b.size(2, 1) ); // P
    }
    else if (Type_equal<order2_type, col2_type>::value)
    {
      cml::mprodt(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(1),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 0),   // M
        a.size(2, 1),   // N
        b.size(2, 1) ); // P
    }
    else
      assert(0);
  }
};


/// CML evaluator for matrix-matrix conjugate products
template <typename Block0,
          typename Block1,
          typename Block2>
struct Evaluator<Op_prod_mm_conj, Cml_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  typedef typename Block0::value_type T;
  typedef typename Block_layout<Block0>::order_type order0_type;
  typedef typename Block_layout<Block1>::order_type order1_type;
  typedef typename Block_layout<Block2>::order_type order2_type;

  static bool const ct_valid = 
    // check that CML supports this data type and/or layout
    impl::cml::Cml_supports_block<Block0>::valid &&
    impl::cml::Cml_supports_block<Block1>::valid &&
    impl::cml::Cml_supports_block<Block2>::valid &&
    // check that all data types are equal
    Type_equal<T, typename Block1::value_type>::value &&
    Type_equal<T, typename Block2::value_type>::value &&
    // check that the complex layouts are equal
    Is_split_block<Block0>::value == Is_split_block<Block1>::value &&
    Is_split_block<Block0>::value == Is_split_block<Block2>::value &&
    // check that the layout is row-major for the first input and the output
    Type_equal<order0_type, row2_type>::value && 
    Type_equal<order1_type, row2_type>::value && 
    // check that direct access is supported
    Ext_data_cost<Block0>::value == 0 &&
    Ext_data_cost<Block1>::value == 0 &&
    Ext_data_cost<Block2>::value == 0;

  static bool rt_valid(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // For 'b', the dimension with the smallest stride must be one,
    // which depends on whether it is row- or column-major.
    bool is_b_row = Type_equal<order2_type, row2_type>::value;
    stride_type b_stride = is_b_row ? ext_b.stride(1) : ext_b.stride(0);

    return
      // ensure the data is unit-stide
      ( ext_r.stride(1) == 1 &&
        ext_a.stride(1) == 1 &&
        b_stride        == 1 );
  }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    Ext_data<Block0> ext_r(const_cast<Block0&>(r));
    Ext_data<Block1> ext_a(const_cast<Block1&>(a));
    Ext_data<Block2> ext_b(const_cast<Block2&>(b));

    // Either row- or column-major layouts are supported for
    // the second input by mapping them to the normal product
    // or transpose product respectively.
    if (Type_equal<order2_type, row2_type>::value)
    {
      cml::mprodj(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(0),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 0),   // M
        a.size(2, 1),   // N
        b.size(2, 1) ); // P
    }
    else if (Type_equal<order2_type, col2_type>::value)
    {
      cml::mprodh(
        ext_a.data(), ext_a.stride(0),
        ext_b.data(), ext_b.stride(1),
        ext_r.data(), ext_r.stride(0),
        a.size(2, 0),   // M
        a.size(2, 1),   // N
        b.size(2, 1) ); // P
    }
    else
      assert(0);
  }
};


} // namespace vsip::impl::dispatcher

} // namespace vsip::impl

} // namespace vsip

#endif // VSIP_OPT_CBE_CML_MATVEC_HPP
