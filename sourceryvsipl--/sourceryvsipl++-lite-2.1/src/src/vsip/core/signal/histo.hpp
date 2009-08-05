/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/signal-histo.hpp
    @author  Don McCoy
    @date    2005-11-29
    @brief   VSIPL++ Library: Histogram functions [signal.histo]
*/

#ifndef VSIP_CORE_SIGNAL_HISTO_HPP
#define VSIP_CORE_SIGNAL_HISTO_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{


template <template <typename, typename> class const_View = const_Vector,
          typename T = VSIP_DEFAULT_VALUE_TYPE>
class Histogram
{
public:
  // Constructor and destructor [signal.histo.constructors]
  Histogram(T min_value, T max_value, length_type num_bin)
    VSIP_THROW((std::bad_alloc))
    : min_(min_value),
      max_(max_value),
      num_bin_(num_bin),
      hist_(num_bin, 0)
  {
    assert(min_ < max_);
    assert(num_bin_ >= 3);
    
    delta_ = (max_ - min_) / (num_bin_ - 2);
  }

  /// This constructor, albeit not required by the VSIPL++ spec, is
  /// necessary to implement the C-VSIPL bindings. They allow an external
  /// view to be used to accumulate the histogram values over multiple
  /// calls outside any persistent Histogram<> object / state.
  template <typename Block>
  Histogram(T min_value, T max_value, const_Vector<int, Block> hist)
    VSIP_THROW((std::bad_alloc))
    : min_(min_value),
      max_(max_value),
      num_bin_(hist.size()),
      hist_(num_bin_)
  {
    assert(min_ < max_);
    assert(num_bin_ >= 3);
    delta_ = (max_ - min_) / (num_bin_ - 2);
    hist_ = hist;
  }

  ~Histogram() VSIP_NOTHROW
  {}

  // Histogram operators [signal.histo.operators]
  template <typename Block>
  const_Vector<scalar_i>
  operator()(const_Vector<T, Block> data,
             bool accumulate = false)
    VSIP_NOTHROW
  {
    if (accumulate == false)
      hist_ = 0;
    
    for (index_type i = 0; i < data.size(); ++i)
      hist_(impl_bin(data(i)))++;

    return hist_;
  }

  template <typename Block>
  const_Vector<scalar_i>
  operator()(const_Matrix<T, Block> data,
             bool accumulate = false)
    VSIP_NOTHROW
  {
    if (accumulate == false)
      hist_ = 0;

    for (index_type i = 0; i < data.size(0); ++i)
      for (index_type j = 0; j < data.size(1); ++j)
        hist_(impl_bin(data(i, j)))++;

    return hist_;
  }

  inline index_type
  impl_bin(T value)
  {
    if (value < min_)
      return 0;
    else if (value >= max_)
      return num_bin_ - 1;
    else
      return (index_type)(((value - min_) / delta_) + 1);
  } 

private:
  T min_;
  T max_;
  T delta_;
  length_type const num_bin_;
  Vector<scalar_i> hist_;
};




} // namespace vsip

#endif // VSIP_CORE_SIGNAL_HISTO_HPP
