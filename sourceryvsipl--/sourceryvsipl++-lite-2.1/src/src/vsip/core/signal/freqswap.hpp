/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/signal/freqswap.hpp
    @author  Don McCoy
    @date    2005-11-29
    @brief   VSIPL++ Library: Frequency swap functions [signal.freqswap]
*/

#ifndef VSIP_CORE_SIGNAL_FREQSWAP_HPP
#define VSIP_CORE_SIGNAL_FREQSWAP_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/domain_utils.hpp>
#include <vsip/core/block_traits.hpp>
#ifndef VSIP_IMPL_REF_IMPL
# include <vsip/opt/expr/return_block.hpp>
#endif


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

template <typename B, dimension_type D = B::dim> class Freqswap_functor;

template <typename B>
class Freqswap_functor<B, 1>
{
public:
  typedef typename B::map_type   map_type;
  typedef typename B::value_type value_type;
  typedef typename View_block_storage<B>::plain_type block_ref_type;

  typedef Freqswap_functor<typename Distributed_local_block<B>::type, 1>
		local_type;

  Freqswap_functor(block_ref_type in) : in_(in) {}

  length_type size(dimension_type block_dim, dimension_type d) const
  {
    assert(block_dim == 1);
    return in_.size(block_dim, d);
  }
  length_type size() const { return in_.size();}

  template <typename B1>
  void apply(B1 &out) const
  {
    // equiv. to r[i] = a[(M/2 + i) mod M],  where i = 0 --> M - 1

    length_type const M = in_.size();

    index_type const ia = M % 2;  // adjustment to get source index

    value_type mid = in_.get(M / 2);

    for (index_type i = 0, ii = M / 2; i < M / 2; ++i, ++ii)
    {
      // Be careful to allow 'out' to alias 'in'
      value_type tmp = in_.get(ii + ia);
      out.put(ii, in_.get(i));
      out.put(i, tmp);
    }

    // if odd, fill in the last row/column(s)
    if (ia)
      out.put(M-1, mid);
  }

  local_type local() const
  {
    return local_type(get_local_block(in_));
  }

  map_type const& map() const { return in_.map();}

private:
  block_ref_type in_;
};

template <typename B>
class Freqswap_functor<B, 2>
{
public:
  typedef typename B::map_type   map_type;
  typedef typename B::value_type value_type;
  typedef typename View_block_storage<B>::plain_type block_ref_type;

  typedef Freqswap_functor<typename Distributed_local_block<B>::type, 2>
		local_type;

  Freqswap_functor(block_ref_type in) : in_(in) {}

  length_type size(dimension_type block_dim, dimension_type d) const
  {
    assert(block_dim == 2);
    return in_.size(block_dim, d);
  }
  length_type size() const { return in_.size();}

  template <typename B1>
  void apply(B1 & out) const
  {
    length_type const M = in_.size(2, 0);
    length_type const N = in_.size(2, 1);

    index_type const ia = M % 2;  // adjustment to get source index
    index_type const ja = N % 2;

    if (is_same_block(in_, out))
    {
      // If in-place, use algorithm that trades O(1) temporary storage
      // for extra copies.

      // First swap rows.
      for (index_type i=0; i < M; ++i)
      {
	value_type mid = in_.get(i, N / 2);

	for (index_type j = 0, jj = N / 2; j < N / 2; ++j, ++jj)
	{
	  // Be careful to allow 'out' to alias 'in'
	  value_type tmp = in_.get(i, jj + ja);
	  out.put(i, jj, in_.get(i, j));
	  out.put(i, j, tmp);
	}

	// if odd, fill in the last row/column(s)
	if (ja) out.put(i, N-1, mid);
      }

      // Second, swap columns.
      for (index_type j=0; j < N; ++j)
      {
	value_type mid = in_.get(M / 2, j);

	for (index_type i = 0, ii = M / 2; i < M / 2; ++i, ++ii)
	{
	  // Be careful to allow 'out' to alias 'in'
	  value_type tmp = in_.get(ii + ia, j);
	  out.put(ii, j, in_.get(i, j));
	  out.put(i,  j, tmp);
	}

	// if odd, fill in the last row/column(s)
	if (ia) out.put(M-1, j, mid);
      }
    }
    else
    {
      // equiv. to out[i,j] = in[(M/2 + i) mod M,(N/2 + i) mod N], 
      //   where i = 0 --> M - 1 and j = 0 --> N - 1
      for (index_type i = 0, ii = M / 2; i < M / 2; ++i, ++ii)
      {
	for (index_type j = 0, jj = N / 2; j < N / 2; ++j, ++jj)
	{
	  value_type tmp = in_.get(ii + ia, jj + ja);
	  out.put(ii, jj, in_.get(i, j));
	  out.put(i, j, tmp);
	  tmp = in_.get(ii + ia, j);
	  out.put(ii, j, in_.get(i, jj + ja));
	  out.put(i, jj, tmp);
	}
      }

      // if odd, fill in the last row/column(s)
      if (ia)
      {
	index_type i = M / 2;
	index_type ii = M - 1;
	for (index_type j = 0, jj = N / 2; j < N / 2; ++j, ++jj)
	{
	  out.put(ii, jj, in_.get(i, j));
	  out.put(ii,  j, in_.get(i,jj + ja));
	}
      }
      if (ja)
      {
	index_type j = N / 2;
	index_type jj = N - 1;
	for (index_type i = 0, ii = M / 2; i < M / 2; ++i, ++ii)
	{
	  out.put(ii, jj, in_.get(i      , j));
	  out.put( i, jj, in_.get(ii + ia, j));
	}
      }
      if (ia && ja) out.put(M - 1, N - 1, in_.get(M / 2, N / 2));
    }
  }

  local_type local() const
  {
    return local_type(get_local_block(in_));
  }

  map_type const& map() const { return in_.map();}

private:
  block_ref_type in_;
};

template <typename B>
struct Freqswap
{
#ifndef VSIP_IMPL_REF_IMPL
  typedef Return_expr_block<B::dim,
                            typename B::value_type,
                            Freqswap_functor<B> >
    block_type;
  typedef typename B::value_type value_type;
  typedef typename View_of_dim<B::dim, value_type, block_type const>::type view_type;
#else
  typedef Dense<B::dim, typename B::value_type> block_type;
  typedef typename B::value_type value_type;
  typedef typename View_of_dim<B::dim, value_type, block_type>::type view_type;
#endif

  static view_type 
  create(typename View_of_dim<B::dim, value_type, B>::const_type v)
  {
    Freqswap_functor<B> func(v.block());
#ifndef VSIP_IMPL_REF_IMPL
    block_type return_block(func);
    return view_type(return_block);
#else
    view_type result(v);
    func.apply(result.block());
    return result;
#endif
  }
};

template <typename CombineT, typename BlockT>
struct Combine_return_type<CombineT, Freqswap_functor<BlockT> const>
{
  typedef Freqswap_functor<BlockT> const tree_type;
  typedef tree_type type;
};

template <typename CombineT, typename BlockT>
struct Combine_return_type<CombineT, Freqswap_functor<BlockT> >
{
  typedef Freqswap_functor<BlockT> const tree_type;
  typedef tree_type type;
};

template <typename CombineT, typename BlockT>
typename Combine_return_type<CombineT, Freqswap_functor<BlockT> >::type
apply_combine(CombineT const& combine, Freqswap_functor<BlockT> const& rf)
{
  typedef 
    typename Combine_return_type<CombineT, Freqswap_functor<BlockT> >::type
    rf_type;

  return rf_type(apply_combine(combine, rf.in_block()),
		 rf.backend(), rf.workspace());
}

template <typename VisitorT, typename BlockT>
void
apply_leaf(VisitorT const& visitor, Freqswap_functor<BlockT> const& rf)
{
  apply_leaf(visitor, rf.in_block());
}

template <dimension_type MapDim, typename MapT, typename BlockT>
struct Is_par_same_map<MapDim, MapT, Freqswap_functor<BlockT> const>
{
  typedef Freqswap_functor<BlockT> const rf_type;

  static bool value(MapT const& map, rf_type& rf)
  {
    return Is_par_same_map<MapDim, MapT, BlockT>::value(map, rf.in_block());
  }
};

} // namespace impl

/// Swaps halves of a vector, or quadrants of a matrix, to remap zero 
/// frequencies from the origin to the middle.

template <template <typename, typename> class const_View,
          typename T,
          typename B>
typename impl::Freqswap<B>::view_type
freqswap(const_View<T, B> in) VSIP_NOTHROW
{
  return impl::Freqswap<B>::create(in);
}

} // namespace vsip

#endif // VSIP_CORE_SIGNAL_FREQSWAP_HPP
