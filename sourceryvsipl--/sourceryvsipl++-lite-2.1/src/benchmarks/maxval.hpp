/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/maxval.hpp
    @author  Assem Salama
    @date    2007-07-22
    @brief   VSIPL++ Library: Helper file for maxval benchmark

*/
#ifndef BENCHMARKS_MAXVAL_HPP
#define BENCHMARKS_MAXVAL_HPP

// Make a structure to run maxval based on tag.
template <template <typename> class ReduceT,
	  typename T,
          typename Block,
	  vsip::dimension_type dim,
          typename Tag>
struct reduction_op_eval
{

  typedef typename vsip::impl::Block_layout<Block>::order_type order_type;
  typedef vsip::impl::dispatcher::Evaluator<
    vsip::impl::dispatcher::Op_reduce_idx<ReduceT>, Tag,
    void(typename ReduceT<T>::result_type&,
        Block const&,
        vsip::Index<dim>&,
        order_type)> evaluator;

  static void exec(T& r, Block const& a, Index<dim>& idx)
  {
    evaluator::exec(r,a,idx,order_type());
  }
};

// structure to help us create test vectors
template <typename MapT,
          template <typename,typename> class ViewT,
	  typename T,
	  typename Block>
struct create_test_vector_helper {};

// override for when destination is a local view
template<template <typename,typename> class ViewT,
         typename T,
	 typename Block>
struct create_test_vector_helper<vsip::Local_map,ViewT,T,Block>
{
  typedef ViewT<T,Block> dst_view;

  // because this is a local to local, we do a normal assign
  template <typename ViewT1>
  void assign_view(ViewT1 sv)
  { view=sv; };

  // the constructor is very simple too
  create_test_vector_helper(vsip::length_type size) : view(size) {};

  dst_view view;


};

// override for when destination is a distributed view
template<template <typename,typename> class ViewT,
         typename T,
	 typename Block>
struct create_test_vector_helper<vsip::Map<>,ViewT,T,Block>
{
  static vsip::dimension_type const dim = ViewT<T,Block>::dim;
  typedef vsip::Dense<dim,T,typename vsip::impl::Row_major<dim>::type,
    vsip::Map<> > dst_block;
  typedef ViewT<T,dst_block>                 dst_view;

  template <typename ViewT1>
  void assign_view(ViewT1 sv)
  {
    // processor 0 will distribute data to all other procs
    vsip::Vector<vsip::processor_type> pvec_in(1);pvec_in(0)=(0);
    vsip::Map<>                  root_map(pvec_in);
    dst_block                    root_block(sv.size(),root_map);
    dst_view                     root_view(root_block);

    // Ok, now move the vector to the distributed view
    vsip::impl::assign_local(root_view,sv);
    view = root_view;

  };

  create_test_vector_helper(vsip::length_type size) :
    view(size, vsip::Map<>(vsip::num_processors())) {};

  dst_view view;

};

#endif // BENCHMARKS_MAXVAL_HPP
