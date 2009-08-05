/* Copyright (c) 2007 by CodeSourcery, Inc.  All rights reserved. */

/** @file    vsip_csl/stencil/kernel_builder.hpp
    @author  Stefan Seefeld
    @date    2007-05-07
    @brief   VSIPL++ Library: Harness to generate kernel from expression.

*/

#ifndef VSIP_CSL_STENCIL_KERNEL_BUILDER_HPP
#define VSIP_CSL_STENCIL_KERNEL_BUILDER_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip_csl/stencil/expr.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip_csl
{
namespace stencil
{

//
// Kernel builder: initialize kernel coefficients.
//
template <typename T>
struct Kernel
{
public:
  Kernel(vsip::index_type y, vsip::index_type x,
         vsip::length_type h, vsip::length_type w)
    : y_(y), x_(x), height_(h), width_(w), data_(new T[h*w])
  {
    for (int i = h*w - 1; i >= 0; --i) data_[i] = T(0);
  }
  Kernel(Kernel const &k)
    : y_(k.y_), x_(k.x_), height_(k.height_), width_(k.width_),
      data_(new T[height_*width_])
  {
    for (int i = height_*width_ - 1; i >= 0; --i) data_[i] = T(0);
  }
  ~Kernel() { delete [] data_;}
  Kernel &operator+=(Kernel const &k)
  {
    for (int i = height_*width_ - 1; i >= 0; --i) data_[i] += k.data_[i];
    return *this;
  }
  Kernel &operator-=(Kernel const &k)
  {
    for (int i = height_*width_ - 1; i >= 0; --i) data_[i] -= k.data_[i];
    return *this;
  }
  Kernel &operator*=(T s)
  {
    for (int i = height_*width_ - 1; i >= 0; --i) data_[i] *= s;
    return *this;
  }
  Kernel &operator/=(T s)
  {
    for (int i = height_*width_ - 1; i >= 0; --i) data_[i] /= s;
    return *this;
  }

  vsip::index_type origin(vsip::dimension_type d) const 
  { return d == 0 ? y_ : x_;}
  vsip::length_type size(vsip::dimension_type d) const 
  { return d == 0 ? height_ : width_;}

  T const &operator() (vsip::index_type y, vsip::index_type x) const
  { return data_[x + y * width_];}
  T &operator() (vsip::index_type y, vsip::index_type x)
  { return data_[x + y * width_];}
private:
  vsip::index_type y_, x_;
  vsip::length_type height_, width_;
  T *data_;
};

template <typename E, typename T> 
struct Kernel_builder
{
  static void apply(E e, Kernel<T>& k);
};

template <typename L, typename R, typename T>
void build_kernel(Binary_expr<L, R, vsip::impl::op::Add> e, Kernel<T>& k)
{
  Kernel_builder<L, T>::apply(e.left(), k);
  Kernel_builder<R, T>::apply(e.right(), k);
}

template <typename L, typename R, typename T>
void build_kernel(Binary_expr<L, R, vsip::impl::op::Sub> e, Kernel<T>& k)
{
  Kernel_builder<L, T>::apply(e.left(), k);
  Kernel<T> k2(k);
  Kernel_builder<R, T>::apply(e.right(), k2);
  k -= k2;
}

// Assume S to be a scalar (the only case in which multiplication is defined).
template <typename L, typename R, template <typename, typename> class O,
          typename S, typename T>
void build_kernel(Binary_expr<Binary_expr<L, R, O>, S, vsip::impl::op::Mult> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Binary_expr<L, R, O>, T>::apply(e.left(), k1);
  k1 *= T(e.right());
  k += k1;
}

// Assume S to be a scalar (the only case in which multiplication is defined).
template <typename S, 
          typename L, typename R, template <typename, typename> class O,
          typename T>
void build_kernel(Binary_expr<S, Binary_expr<L, R, O>, vsip::impl::op::Mult> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Binary_expr<L, R, O>, T>::apply(e.right(), k1);
  k1 *= T(e.left());
  k += k1;
}

// Assume S to be a scalar (the only case in which division is defined).
template <typename L, typename R, template <typename, typename> class O,
          typename S, typename T>
void build_kernel(Binary_expr<Binary_expr<L, R, O>, S, vsip::impl::op::Div> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Binary_expr<L, R, O>, T>::apply(e.left(), k1);
  k1 /= T(e.right());
  k += k1;
}

// Assume S to be a scalar (the only case in which multiplication is defined).
template <typename V, typename I, typename J,
          typename S, typename T>
void build_kernel(Binary_expr<Call<V, I, J>, S, vsip::impl::op::Mult> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Call<V, I, J>, T>::apply(e.left(), k1);
  k1 *= T(e.right());
  k += k1;
}

// Assume S to be a scalar (the only case in which multiplication is defined).
template <typename S, 
          typename V, typename I, typename J,
          typename T>
void build_kernel(Binary_expr<S, Call<V, I, J>, vsip::impl::op::Mult> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Call<V, I, J>, T>::apply(e.right(), k1);
  k1 *= T(e.left());
  k += k1;
}

// Assume S to be a scalar (the only case in which division is defined).
template <typename V, typename I, typename J,
          typename S, typename T>
void build_kernel(Binary_expr<Call<V, I, J>, S, vsip::impl::op::Div> e,
                  Kernel<T>& k)
{
  Kernel<T> k1(k);
  Kernel_builder<Call<V, I, J>, T>::apply(e.left(), k1);
  k1 /= T(e.right());
  k += k1;
}

template <vsip::dimension_type D>
vsip::stride_type offset(Iterator<D>) { return 0;}

template <vsip::dimension_type D>
vsip::stride_type offset(Offset<D> o) { return o.i;}

template <typename V, typename I, typename J, typename T>
void build_kernel(Call<V, I, J> e, Kernel<T>& k)
{
  vsip::index_type y = k.origin(0) + offset(e.i());
  vsip::index_type x = k.origin(1) + offset(e.j());
  k(y, x) += T(1);
};

template <typename E, typename T>
void 
Kernel_builder<E, T>::apply(E e, Kernel<T>& k) 
{ build_kernel(e, k);}

} // namespace vsip_csl::stencil
} // namespace vsip_csl

#endif
