/* Copyright (c) 2007 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/fft/return_functor.cpp
    @author  Jules Bergmann
    @date    2007-03-09
    @brief   VSIPL++ Library: FFT functor for Return_expr_blocks.
*/

#ifndef VSIP_OPT_FFT_RETURN_FUNCTOR_HPP
#define VSIP_OPT_FFT_RETURN_FUNCTOR_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/block_traits.hpp>
#include <vsip/core/parallel/support_block.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace fft
{

/// Return functor class for Fft.

/// Captures invocation of Fft object on an input block for later
/// evaluation, once the destination block is known.

template <dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
class Fft_return_functor
{
  // Compile-time typedefs.
public:
  typedef BlockT                                              block_type;
  typedef typename block_type::map_type                       map_type;
  typedef typename View_block_storage<block_type>::plain_type block_ref_type;

  typedef Fft_return_functor<Dim, T, 
			     typename Distributed_local_block<BlockT>::type,
			     BackendT, WorkspaceT> local_type;

  // Constructors.
public:
  template <typename GivenViewT>
  Fft_return_functor(
    GivenViewT         in_view,
    Domain<Dim> const& output_size,
    BackendT&          backend,
    WorkspaceT&        workspace)
  : in_block_   (in_view.block()),
    output_size_(output_size),
    backend_    (backend),
    workspace_  (workspace)
  {}

  Fft_return_functor(
    block_ref_type     in_block,
    Domain<Dim> const& output_size,
    BackendT&          backend,
    WorkspaceT&        workspace)
  : in_block_   (in_block),
    output_size_(output_size),
    backend_    (backend),
    workspace_  (workspace)
  {}

  Fft_return_functor(Fft_return_functor const& rhs)
    : in_block_   (rhs.in_block_),
      output_size_(rhs.output_size_),
      backend_    (rhs.backend_),
      workspace_  (rhs.workspace_)
  {}

  // Accessors
public:
  template <typename ResBlockT>
  void apply(ResBlockT& result) const
  {
    workspace_.by_reference_blk(&backend_, in_block_, result);
  }

  length_type size() const
  {
    return output_size_.size();
  }

  length_type size(dimension_type block_dim, dimension_type d) const
  {
    assert(block_dim == Dim);
    return output_size_[d].size();
  }

  local_type local() const
  {
    // The local output size is the same as the global output size
    // along the dimension of the FFT.  Its size along the other
    // dimension matches that of the input local block.
    length_type rows = output_size_[0].size();
    length_type cols = output_size_[1].size();
    if (BackendT::axis == 0)
    {
      cols = block_subblock_domain<2>(in_block_)[1].size();
      rows = (cols == 0) ? 0 : rows;
    }
    else
    {
      rows = block_subblock_domain<2>(in_block_)[0].size();
      cols = (rows == 0) ? 0 : cols;
    }
    Domain<2> l_output_size(rows, cols);

    return local_type(get_local_block(in_block_),
		      l_output_size,
		      backend_,
		      workspace_);
  }

  map_type   const& map()       const { return in_block_.map(); }
  block_type const& block()     const { return in_block_; }
  BackendT   const& backend()   const { return backend_; }
  WorkspaceT const& workspace() const { return workspace_; }

// Member data.
private:
  block_ref_type in_block_;
  Domain<Dim>    output_size_;
  BackendT&      backend_;
  WorkspaceT&    workspace_;
};

} // namespace vsip::impl::fft



// Specliaze parallel block traits for Fft_return_functor.

template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
struct Combine_return_type<CombineT,
	fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> const>
{
  typedef fft::Fft_return_functor<Dim, T,
		typename Combine_return_type<CombineT, BlockT>::tree_type,
		BackendT, WorkspaceT> const
          tree_type;
  typedef tree_type type;
};



template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
struct Combine_return_type<CombineT,
	fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> >
{
  typedef fft::Fft_return_functor<Dim, T,
		typename Combine_return_type<CombineT, BlockT>::tree_type,
		BackendT, WorkspaceT> const
          tree_type;
  typedef tree_type type;
};



template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
typename Combine_return_type<CombineT,
	fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> >::type
apply_combine(
  CombineT const&                                                      combine,
  fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> const& rf)
{
  typedef typename
    Combine_return_type<
		CombineT,
		fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT>
	>::type rf_type;

  return rf_type(apply_combine(combine, rf.in_block()),
		 rf.backend(), rf.workspace());
}



template <typename       VisitorT,
	  dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
void
apply_leaf(
  VisitorT const&                                                      visitor,
  fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> const& rf)
{
  apply_leaf(visitor, rf.in_block());
}



template <dimension_type MapDim,
	  typename       MapT,
	  dimension_type Dim,
	  typename       T,
	  typename       BlockT,
	  typename       BackendT,
	  typename       WorkspaceT>
struct Is_par_same_map<MapDim, MapT,
	fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> const>
{
  typedef fft::Fft_return_functor<Dim, T, BlockT, BackendT, WorkspaceT> const
          rf_type;

  static bool value(MapT const& map, rf_type& rf)
  {
    return Is_par_same_map<MapDim, MapT, BlockT>::value(map, rf.in_block());
  }
};



} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_FFT_RETURN_FUNCTOR_HPP
