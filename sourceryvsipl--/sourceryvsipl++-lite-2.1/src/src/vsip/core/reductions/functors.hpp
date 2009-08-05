/* Copyright (c) 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/reductions/functors.hpp
    @author  Jules Bergmann
    @date    2006-05-31
    @brief   VSIPL++ Library: Reduction functors.
	     [math.fns.reductions].

*/

#ifndef VSIP_CORE_REDUCTIONS_FUNTORS_HPP
#define VSIP_CORE_REDUCTIONS_FUNTORS_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/reductions/types.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/***********************************************************************
  Reduction Functors
***********************************************************************/

// Evaluator OpTag for value reductions.

namespace dispatcher
{
template <template <typename> class ReduceT>
struct Op_reduce;
}

template <typename T>
struct All_true
{
  static reduction_type const rtype = reduce_all_true;

  typedef T result_type;
  typedef T accum_type;

  static accum_type initial() { return ~(accum_type()); }

  static accum_type update(accum_type state, T new_value)
    { return band(state, new_value); }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type state) { return (state == T()); }
};


template <>
struct All_true<bool>
{
  static reduction_type const rtype = reduce_all_true_bool;

  typedef bool result_type;
  typedef bool accum_type;

  static bool initial() { return true; }

  static bool update(bool state, bool new_value)
    { return land(state, new_value); }

  static bool value(bool state, length_type)
    { return state; }

  static bool done(bool state) { return state == false; }
};



template <typename T>
struct Any_true
{
  static reduction_type const rtype = reduce_any_true;

  typedef T result_type;
  typedef T accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return bor(state, new_value); }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type state) { return (state == ~(T())); }
};


template <>
struct Any_true<bool>
{
  static reduction_type const rtype = reduce_any_true_bool;

  typedef bool result_type;
  typedef bool accum_type;

  static bool initial() { return false; }

  static bool update(bool state, bool new_value)
    { return lor(state, new_value); }

  static bool value(bool state, length_type)
    { return state; }

  static bool done(bool state) { return state == true; }
};


template <typename T>
struct Mean_value
{
  static reduction_type const rtype = reduce_sum;

  typedef T result_type;
  typedef T accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return state + new_value; }

  static accum_type value(accum_type state, length_type size)
    { return state / static_cast<accum_type>(size); }

  static bool done(accum_type) { return false; }
};



template <typename T>
struct Mean_magsq_value
{
  static reduction_type const rtype = reduce_sum;

  typedef typename Scalar_of<T>::type result_type;
  typedef typename Scalar_of<T>::type accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return state + magsq(new_value); }

  static accum_type value(accum_type state, length_type size)
    { return state / size; }

  static bool done(accum_type) { return false; }
};



template <typename T>
struct Sum_value
{
  static reduction_type const rtype = reduce_sum;

  typedef T result_type;
  typedef T accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return state + new_value; }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type) { return false; }
};



template <typename T>
struct Sum_magsq_value
{
  static reduction_type const rtype = reduce_sum;

  typedef typename Scalar_of<T>::type result_type;
  typedef typename Scalar_of<T>::type accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return state + magsq(new_value); }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type) { return false; }
};



/// Specialization for 'bool': return the number of true values.

template <>
struct Sum_value<bool>
{
  static reduction_type const rtype = reduce_sum;

  typedef length_type result_type;
  typedef length_type accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, bool new_value)
    { return state + (new_value ? 1 : 0); }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type) { return false; }
};



template <typename T>
struct Sum_sq_value
{
  static reduction_type const rtype = reduce_sum;

  typedef T result_type;
  typedef T accum_type;

  static accum_type initial() { return accum_type(); }

  static accum_type update(accum_type state, T new_value)
    { return state + new_value * new_value; }

  static accum_type value(accum_type state, length_type)
    { return state; }

  static bool done(accum_type) { return false; }
};



/***********************************************************************
  Reduction-Index Functors
***********************************************************************/

// Evaluator OpTag for reductions returning index.

namespace dispatcher
{
template <template <typename> class ReduceT>
struct Op_reduce_idx;
}


template <typename T>
class Max_value
{
public:
  typedef T result_type;

  Max_value(T init) : value_(init) {}

  bool next_value(T value)
  {
    if (value > value_)
    {
      value_ = value;
      return true;
    }
    else return false;
  }

  T value() { return value_; }

private:
  T value_;
};



template <typename T>
class Min_value
{
public:
  typedef T result_type;

  Min_value(T init) : value_(init) {}

  bool next_value(T value)
  {
    if (value < value_)
    {
      value_ = value;
      return true;
    }
    else return false;
  }

  T value() { return value_; }

private:
  T value_;
};



template <typename T>
class Max_mag_value
{
public:
  typedef typename Scalar_of<T>::type result_type;

  Max_mag_value(T init) : value_(mag(init)) {}

  bool next_value(T value)
  {
    result_type tmp = mag(value);
    if (tmp > value_)
    {
      value_ = tmp;
      return true;
    }
    else return false;
  }

  result_type value() { return value_; }

private:
  result_type value_;
};



template <typename T>
class Min_mag_value
{
public:
  typedef typename Scalar_of<T>::type result_type;

  Min_mag_value(T init) : value_(mag(init)) {}

  bool next_value(T value)
  {
    result_type tmp = mag(value);
    if (tmp < value_)
    {
      value_ = tmp;
      return true;
    }
    else return false;
  }

  result_type value() { return value_; }

private:
  result_type value_;
};



template <typename T>
class Max_magsq_value
{
public:
  typedef typename Scalar_of<T>::type result_type;

  Max_magsq_value(T init) : value_(magsq(init)) {}

  bool next_value(T value)
  {
    result_type tmp = magsq(value);
    if (tmp > value_)
    {
      value_ = tmp;
      return true;
    }
    else return false;
  }

  result_type value() { return value_; }

private:
  result_type value_;
};



template <typename T>
class Min_magsq_value
{
public:
  typedef typename Scalar_of<T>::type result_type;

  Min_magsq_value(T init) : value_(magsq(init)) {}

  bool next_value(T value)
  {
    result_type tmp = magsq(value);
    if (tmp < value_)
    {
      value_ = tmp;
      return true;
    }
    else return false;
  }

  result_type value() { return value_; }

private:
  result_type value_;
};

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_REDUCTIONS_FUNTORS_HPP
