/* Copyright (c) 2005, 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/reductions/reductions.hpp
    @author  Jules Bergmann
    @date    2005-07-12
    @brief   VSIPL++ Library: Reduction functions.
	     [math.fns.reductions].

*/

#ifndef VSIP_CORE_REDUCTIONS_REDUCTIONS_HPP
#define VSIP_CORE_REDUCTIONS_REDUCTIONS_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/tensor.hpp>
#include <vsip/core/reductions/functors.hpp>
#include <vsip/core/parallel/services.hpp>
#include <vsip/core/impl_tags.hpp>
#include <vsip/core/dispatch.hpp>
#if !VSIP_IMPL_REF_IMPL
#  include <vsip/opt/dispatch.hpp>
#  ifdef VSIP_IMPL_HAVE_SAL
#    include <vsip/opt/sal/eval_reductions.hpp>
#  endif
#endif
#if VSIP_IMPL_HAVE_CVSIP
#  include <vsip/core/cvsip/eval_reductions.hpp>
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace dispatcher
{

/***********************************************************************
  Generic evaluators.
***********************************************************************/

#ifndef VSIP_IMPL_REF_IMPL
template<template <typename> class ReduceT>
struct List<Op_reduce<ReduceT> >
{
  typedef Make_type_list<Cvsip_tag, Mercury_sal_tag, Generic_tag>::type type;
};
#endif


/// Generic evaluator for vector reductions.
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag, 
                 void(T&, Block const&, row1_type, Int_type<1>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, row1_type, Int_type<1>)
  { return true; }

  static void exec(T& r, Block const& a, row1_type, Int_type<1>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<VT>::accum_type state = ReduceT<VT>::initial();

    length_type length = a.size(1, 0);
    PRAGMA_IVDEP
    for (index_type i=0; i<length; ++i)
    {
      state = ReduceT<VT>::update(state, a.get(i));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length);
  }
};


/// Generic evaluator for matrix reductions (tuple<0, 1, 2>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, row2_type, Int_type<2>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, row2_type, Int_type<2>)
  { return true; }

  static void exec(T& r, Block const& a, row2_type, Int_type<2>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_i = a.size(2, 0);
    length_type length_j = a.size(2, 1);

    for (index_type i=0; i<length_i; ++i)
      PRAGMA_IVDEP
      for (index_type j=0; j<length_j; ++j)
      {
	state = ReduceT<VT>::update(state, a.get(i, j));
	if (ReduceT<VT>::done(state)) break;
      }

    r = ReduceT<VT>::value(state, length_i*length_j);
  }
};


/// Generic evaluator for matrix reductions (tuple<2, 1, 0>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, col2_type, Int_type<2>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, col2_type, Int_type<2>)
  { return true; }

  static void exec(T& r, Block const& a, col2_type, Int_type<2>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_i = a.size(2, 0);
    length_type length_j = a.size(2, 1);

    for (index_type j=0; j<length_j; ++j)
    PRAGMA_IVDEP
    for (index_type i=0; i<length_i; ++i)
    {
      state = ReduceT<VT>::update(state, a.get(i, j));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_i*length_j);
  }
};


/// Generic evaluator for tensor reductions (tuple<0, 1, 2>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<0, 1, 2>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<0, 1, 2>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<0, 1, 2>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type i=0; i<length_0; ++i)
    for (index_type j=0; j<length_1; ++j)
    for (index_type k=0; k<length_2; ++k)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};


/// Generic evaluator for tensor reductions (tuple<0, 2, 1>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<0, 2, 1>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<0, 2, 1>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<0, 2, 1>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type i=0; i<length_0; ++i)
    for (index_type k=0; k<length_2; ++k)
    for (index_type j=0; j<length_1; ++j)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};


/// Generic evaluator for tensor reductions (tuple<1, 0, 2>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<1, 0, 2>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<1, 0, 2>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<1, 0, 2>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type j=0; j<length_1; ++j)
    for (index_type i=0; i<length_0; ++i)
    for (index_type k=0; k<length_2; ++k)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};


/// Generic evaluator for tensor reductions (tuple<1, 2, 0>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<1, 2, 0>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<1, 2, 0>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<1, 2, 0>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type j=0; j<length_1; ++j)
    for (index_type k=0; k<length_2; ++k)
    for (index_type i=0; i<length_0; ++i)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};


/// Generic evaluator for tensor reductions (tuple<2, 0, 1>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<2, 0, 1>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<2, 0, 1>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<2, 0, 1>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type k=0; k<length_2; ++k)
    for (index_type i=0; i<length_0; ++i)
    for (index_type j=0; j<length_1; ++j)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};


/// Generic evaluator for tensor reductions (tuple<2, 1, 0>).
template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block>
struct Evaluator<Op_reduce<ReduceT>, Generic_tag,
                 void(T&, Block const&, tuple<2, 1, 0>, Int_type<3>)>
{
  static bool const ct_valid = true;
  static bool rt_valid(T&, Block const&, tuple<2, 1, 0>, Int_type<3>)
  { return true; }

  static void exec(T& r, Block const& a, tuple<2, 1, 0>, Int_type<3>)
  {
    typedef typename Block::value_type VT;
    typename ReduceT<T>::accum_type state = ReduceT<VT>::initial();

    length_type length_0 = a.size(3, 0);
    length_type length_1 = a.size(3, 1);
    length_type length_2 = a.size(3, 2);

    for (index_type k=0; k<length_2; ++k)
    for (index_type j=0; j<length_1; ++j)
    for (index_type i=0; i<length_0; ++i)
    {
      state = ReduceT<VT>::update(state, a.get(i, j, k));
      if (ReduceT<VT>::done(state)) break;
    }

    r = ReduceT<VT>::value(state, length_0*length_1*length_2);
  }
};



/***********************************************************************
  Parallel evaluators.
***********************************************************************/

#if !VSIP_IMPL_REF_IMPL
template <template <typename> class ReduceT,
          typename                  T,
	  typename                  Block,
	  typename                  OrderT,
	  int                       Dim>
struct Par_reduction_eval_base
{
  static bool const ct_valid = 
    !Is_local_map<typename Block::map_type>::value &&
    Reduction_supported<ReduceT<typename Block::value_type>::rtype,
                        typename Block::value_type>::value;

  static bool rt_valid(T&, Block const&, OrderT, Int_type<Dim>)
  { return true; }
};


template <typename                  T,
	  typename                  Block,
	  typename                  OrderT,
	  int                       Dim>
struct Evaluator<Op_reduce<Mean_value>, Parallel_tag,
                 void(T&, Block const&, OrderT, Int_type<Dim>)>
  : Par_reduction_eval_base<Mean_value, T, Block, OrderT, Dim>
{
  static void exec(T& r, Block const& a, OrderT, Int_type<Dim>)
  {
    typedef typename Block::value_type VT;
    T l_r;
    typedef typename Distributed_local_block<Block>::type local_block_type;
    typedef typename Block_layout<local_block_type>::order_type order_type;
    typedef Int_type<Dim>                                       dim_type;
    typedef Mean_value<VT>                                      reduce_type;
    typedef typename Block::map_type                            map_type;

    dispatch<Op_reduce<Sum_value>, void,
             typename Sum_value<VT>::result_type&,
             local_block_type const&,
             order_type,
             dim_type>
      (l_r, get_local_block(a), order_type(), dim_type());

    if (!Type_equal<map_type, Global_map<Block::dim> >::value)
	r = a.map().impl_comm().allreduce(reduce_type::rtype, l_r);
    else
      r = l_r;

    r /= static_cast<typename reduce_type::accum_type>(a.size());
  }
};



template <typename                  T,
	  typename                  Block,
	  typename                  OrderT,
	  int                       Dim>
struct Evaluator<Op_reduce<Mean_magsq_value>, Parallel_tag,
                 void(T&, Block const&, OrderT, Int_type<Dim>)>
  : Par_reduction_eval_base<Mean_magsq_value, T, Block, OrderT, Dim>
{
  static void exec(T& r, Block const& a, OrderT, Int_type<Dim>)
  {
    typedef typename Block::value_type VT;
    T l_r;
    typedef typename Distributed_local_block<Block>::type local_block_type;
    typedef typename Block_layout<local_block_type>::order_type order_type;
    typedef Int_type<Dim>                                       dim_type;
    typedef Mean_magsq_value<VT>                                reduce_type;
    typedef typename Block::map_type                            map_type;

    dispatch<Op_reduce<Sum_magsq_value>, void,
             typename Sum_magsq_value<VT>::result_type&,
             local_block_type const&,
             order_type,
             dim_type>
      (l_r, get_local_block(a), order_type(), dim_type());

    if (!Type_equal<map_type, Global_map<Block::dim> >::value)
      r = a.map().impl_comm().allreduce(reduce_type::rtype, l_r);
    else
      r = l_r;

    r /= static_cast<typename reduce_type::accum_type>(a.size());
  }
};


template <template <typename> class ReduceT,
	  typename                  T,
	  typename                  Block,
	  typename                  OrderT,
	  int                       Dim>
struct Evaluator<Op_reduce<ReduceT>, Parallel_tag, 
                 void(T&, Block const&, OrderT, Int_type<Dim>)>
  : Par_reduction_eval_base<ReduceT, T, Block, OrderT, Dim>
{
  static void exec(T& r, Block const& a, OrderT, Int_type<Dim>)
  {
    typedef typename Block::value_type VT;
    T l_r = T();
    typedef typename Distributed_local_block<Block>::type local_block_type;
    typedef typename Block_layout<local_block_type>::order_type order_type;
    typedef Int_type<Dim>                                       dim_type;
    typedef ReduceT<VT>                                         reduce_type;
    typedef typename Block::map_type                            map_type;

    dispatch<Op_reduce<ReduceT>, void,
             typename ReduceT<VT>::result_type&,
	     local_block_type const&,
             order_type,
             dim_type>
      (l_r, get_local_block(a), order_type(), dim_type());

    if (!Type_equal<map_type, Global_map<Block::dim> >::value)
      r = a.map().impl_comm().allreduce(ReduceT<T>::rtype, l_r);
    else
      r = l_r;
  }
};
#endif

} // namespace vsip::impl::dispatcher



template <template <typename> class ReduceT,
	  typename                  ViewT>
typename ReduceT<typename ViewT::value_type>::result_type
reduce(ViewT v)
{
  typedef typename ViewT::value_type T;
  typedef typename ReduceT<T>::result_type result_type;
  typedef typename ViewT::block_type block_type;
  typedef typename Block_layout<typename ViewT::block_type>::order_type
    order_type;
  typedef Int_type<ViewT::dim> dim_type;

  result_type r;

#if VSIP_IMPL_REF_IMPL
  dispatcher::Evaluator<dispatcher::Op_reduce<ReduceT>, Cvsip_tag, 
    void(result_type&, block_type const&, order_type, dim_type)>::
    exec(r, v.block(), order_type(), dim_type());
#else

  typedef dispatcher::Make_type_list<Parallel_tag, Cvsip_tag, Mercury_sal_tag, 
    Generic_tag>::type list_type;

  dispatcher::Dispatcher<dispatcher::Op_reduce<ReduceT>, 
    void(result_type&, block_type const&, order_type, dim_type), list_type>::
    dispatch(r, v.block(), order_type(), dim_type());
#endif

  return r;
}

} // namespace vsip::impl



/***********************************************************************
  API Reduction Functions
***********************************************************************/

template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::All_true<T>::result_type
alltrue(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::All_true>(v);
}



template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::Any_true<T>::result_type
anytrue(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::Any_true>(v);
}



template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::Mean_value<T>::result_type
meanval(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::Mean_value>(v);
}



// Note: meansqval computes the mean of the magnitude square

template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::Mean_magsq_value<T>::result_type
meansqval(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::Mean_magsq_value>(v);
}



template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::Sum_value<T>::result_type
sumval(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::Sum_value>(v);
}



template <typename                            T,
	  template <typename, typename> class ViewT,
	  typename                            BlockT>
typename impl::Sum_sq_value<T>::result_type
sumsqval(ViewT<T, BlockT> v)
{
  typedef typename impl::Block_layout<BlockT>::order_type order_type;
  return impl::reduce<impl::Sum_sq_value>(v);
}



} // namespace vsip

#endif // VSIP_CORE_REDUCTIONS_REDUCTIONS_HPP
