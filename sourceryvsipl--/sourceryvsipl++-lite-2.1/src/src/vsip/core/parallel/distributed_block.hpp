/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/parallel/distributed-block.hpp
    @author  Jules Bergmann
    @date    2005-03-22
    @brief   VSIPL++ Library: Distributed block class.

*/

#ifndef VSIP_IMPL_DISTRIBUTED_BLOCK_HPP
#define VSIP_IMPL_DISTRIBUTED_BLOCK_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/map_fwd.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/domain_utils.hpp>
#include <vsip/core/parallel/get_local_view.hpp>
#include <vsip/core/parallel/proxy_local_block.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

template <typename Block,
	  typename Map>
class Distributed_block
  : public impl::Ref_count<Distributed_block<Block, Map> >
{
  // Compile-time values and types.
public:
  static dimension_type const dim = Block::dim;

  typedef typename Block::value_type           value_type;
  typedef typename Block::reference_type       reference_type;
  typedef typename Block::const_reference_type const_reference_type;

  typedef typename Block_layout<Block>::complex_type impl_complex_type;
  typedef Storage<impl_complex_type, value_type>     impl_storage_type;
  typedef typename impl_storage_type::type           impl_data_type;
  typedef typename impl_storage_type::const_type     impl_const_data_type;

  typedef Map                                  map_type;

  // Non-standard typedefs:
  typedef Block                                local_block_type;
  typedef typename Block_layout<local_block_type>::layout_type
                                               local_LP;
  typedef Proxy_local_block<dim, value_type, local_LP>  proxy_local_block_type;

  // Private compile-time values and types.
private:
  enum private_type {};
  typedef typename impl::Complex_value_type<value_type, private_type>::type uT;
  typedef typename impl::Complex_value_type<value_type, private_type>::ptr_type
          uT_ptr;

  // Constructors and destructor.
public:
  Distributed_block(Domain<dim> const& dom, Map const& map = Map())
    : map_           (map),
      proc_          (local_processor()),
      sb_            (map_.subblock(proc_)),
      subblock_      (NULL)
  {
    map_.impl_apply(dom);
    for (dimension_type d=0; d<dim; ++d)
      size_[d] = dom[d].length();

    Domain<dim> sb_dom = 
      (sb_ != no_subblock) ? map_.template impl_subblock_domain<dim>(sb_)
                           : empty_domain<dim>();

    subblock_ = new Block(sb_dom);
  }

  Distributed_block(Domain<dim> const& dom, value_type value,
		    Map const& map = Map())
    : map_           (map),
      proc_          (local_processor()),
      sb_            (map_.subblock(proc_)),
      subblock_      (NULL)
  {
    map_.impl_apply(dom);
    for (dimension_type d=0; d<dim; ++d)
      size_[d] = dom[d].length();

    Domain<dim> sb_dom = 
      (sb_ != no_subblock) ? map_.template impl_subblock_domain<dim>(sb_)
                           : empty_domain<dim>();

    subblock_ = new Block(sb_dom, value);
  }

  Distributed_block(
    Domain<dim> const& dom, 
    value_type* const  ptr,
    Map const&         map = Map())
    : map_           (map),
      proc_          (local_processor()),
      sb_            (map_.subblock()),
      subblock_      (NULL)
  {
    map_.impl_apply(dom);
    for (dimension_type d=0; d<dim; ++d)
      size_[d] = dom[d].length();

    Domain<dim> sb_dom = 
      (sb_ != no_subblock) ? map_.template impl_subblock_domain<dim>(sb_)
                           : empty_domain<dim>();

    subblock_ = new Block(sb_dom, ptr);
  }

  Distributed_block(
    Domain<dim> const& dom, 
    uT_ptr const       pointer,
    Map const&         map = Map())
    : map_           (map),
      proc_          (local_processor()),
      sb_            (map_.subblock()),
      subblock_      (NULL)
  {
    map_.impl_apply(dom);
    for (dimension_type d=0; d<dim; ++d)
      size_[d] = dom[d].length();

    Domain<dim> sb_dom = 
      (sb_ != no_subblock) ? map_.template impl_subblock_domain<dim>(sb_)
                           : empty_domain<dim>();

    subblock_ = new Block(sb_dom, pointer);
  }

  Distributed_block(
    Domain<dim> const& dom, 
    uT_ptr const       real_pointer,
    uT_ptr const       imag_pointer,
    Map const&         map = Map())
    : map_           (map),
      proc_          (local_processor()),
      sb_            (map_.subblock()),
      subblock_      (NULL)
  {
    map_.impl_apply(dom);
    for (dimension_type d=0; d<dim; ++d)
      size_[d] = dom[d].length();

    Domain<dim> sb_dom = 
      (sb_ != no_subblock) ? map_.template impl_subblock_domain<dim>(sb_)
                           : empty_domain<dim>();

    subblock_ = new Block(sb_dom, real_pointer, imag_pointer);
  }


  ~Distributed_block()
  {
    if (subblock_)
    {
      // PROFILE: issue a warning if subblock is captured.
      subblock_->decrement_count();
    }
  }
    
  // Data accessors.
public:
  // get() on a distributed_block is a broadcast.  The processor
  // owning the index broadcasts the value to the other processors in
  // the data parallel group.
  value_type get(index_type idx) const VSIP_NOTHROW
  {
    // Optimize uni-processor and replicated cases.
    if (    Type_equal<Map, Global_map<1> >::value
        ||  vsip::num_processors() == 1
	|| (   Type_equal<Map, Replicated_map<1> >::value
	    && map_.subblock() != no_subblock))
      return subblock_->get(idx);

    index_type     sb = map_.impl_subblock_from_global_index(Index<1>(idx));
    processor_type pr = *(map_.processor_begin(sb));
    value_type     val = value_type(); // avoid -Wall 'may not be initialized'

    if (pr == proc_)
    {
      assert(map_.subblock() == sb);
      index_type lidx = map_.impl_local_from_global_index(0, idx);
      val = subblock_->get(lidx);
    }

    map_.impl_comm().broadcast(pr, &val, 1);

    return val;
  }

  value_type get(index_type idx0, index_type idx1) const VSIP_NOTHROW
  {
    // Optimize uni-processor and replicated cases.
    if (    Type_equal<Map, Global_map<2> >::value
        ||  vsip::num_processors() == 1
	|| (   Type_equal<Map, Replicated_map<2> >::value
	    && map_.subblock() != no_subblock))
      return subblock_->get(idx0, idx1);

    index_type     sb = map_.impl_subblock_from_global_index(
				Index<2>(idx0, idx1));
    processor_type pr = *(map_.processor_begin(sb));
    value_type     val = value_type(); // avoid -Wall 'may not be initialized'

    if (pr == proc_)
    {
      assert(map_.subblock() == sb);
      index_type l_idx0 = map_.impl_local_from_global_index(0, idx0);
      index_type l_idx1 = map_.impl_local_from_global_index(1, idx1);
      val = subblock_->get(l_idx0, l_idx1);
    }

    map_.impl_comm().broadcast(pr, &val, 1);

    return val;
  }

  value_type get(index_type idx0, index_type idx1, index_type idx2)
    const VSIP_NOTHROW
  {
    // Optimize uni-processor and replicated cases.
    if (    Type_equal<Map, Global_map<3> >::value
        ||  vsip::num_processors() == 1
	|| (   Type_equal<Map, Replicated_map<3> >::value
	    && map_.subblock() != no_subblock))
      return subblock_->get(idx0, idx1, idx2);

    index_type     sb = map_.impl_subblock_from_global_index(
				Index<3>(idx0, idx1, idx2));
    processor_type pr = *(map_.processor_begin(sb));
    value_type     val = value_type(); // avoid -Wall 'may not be initialized'

    if (pr == proc_)
    {
      assert(map_.subblock() == sb);
      index_type l_idx0 = map_.impl_local_from_global_index(0, idx0);
      index_type l_idx1 = map_.impl_local_from_global_index(1, idx1);
      index_type l_idx2 = map_.impl_local_from_global_index(2, idx2);
      val = subblock_->get(l_idx0, l_idx1, l_idx2);
    }

    map_.impl_comm().broadcast(pr, &val, 1);

    return val;
  }


  // put() on a distributed_block is executed only on the processor
  // owning the index.
  void put(index_type idx, value_type val) VSIP_NOTHROW
  {
    index_type     sb = map_.impl_subblock_from_global_index(Index<1>(idx));

    if (map_.subblock() == sb)
    {
      index_type lidx = map_.impl_local_from_global_index(0, idx);
      subblock_->put(lidx, val);
    }
  }

  void put(index_type idx0, index_type idx1, value_type val) VSIP_NOTHROW
  {
    index_type sb = map_.impl_subblock_from_global_index(Index<2>(idx0, idx1));

    if (map_.subblock() == sb)
    {
      index_type l_idx0 = map_.impl_local_from_global_index(0, idx0);
      index_type l_idx1 = map_.impl_local_from_global_index(1, idx1);
      subblock_->put(l_idx0, l_idx1, val);
    }
  }

  void put(index_type idx0, index_type idx1, index_type idx2, value_type val)
    VSIP_NOTHROW
  {
    index_type     sb = map_.impl_subblock_from_global_index(
				Index<3>(idx0, idx1, idx2));

    if (map_.subblock() == sb)
    {
      index_type l_idx0 = map_.impl_local_from_global_index(0, idx0);
      index_type l_idx1 = map_.impl_local_from_global_index(1, idx1);
      index_type l_idx2 = map_.impl_local_from_global_index(2, idx2);
      subblock_->put(l_idx0, l_idx1, l_idx2, val);
    }
  }


  // Support Direct_data interface.
public:
  impl_data_type       impl_data()       VSIP_NOTHROW
  { return subblock_->impl_data(); }

  impl_const_data_type impl_data() const VSIP_NOTHROW
  { return subblock_->impl_data(); }

  stride_type impl_stride(dimension_type D, dimension_type d)
    const VSIP_NOTHROW
  { return subblock_->impl_stride(D, d); }


  // Accessors.
public:
  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type, dimension_type) const VSIP_NOTHROW;
  map_type const& map() const VSIP_NOTHROW 
    { return map_; }

  // length_type num_local_blocks() const { return num_subblocks_; }

  Block& get_local_block() const
  {
    assert(subblock_ != NULL);
    return *subblock_;
  }

  index_type subblock() const { return sb_; }

  void assert_local(index_type sb) const
    { assert(sb == sb_ && subblock_ != NULL); }

  // User storage functions.
public:
  void admit(bool update = true) VSIP_NOTHROW
    { subblock_->admit(update); }

  void release(bool update = true) VSIP_NOTHROW
    { subblock_->release(update); }

  void release(bool update, value_type*& pointer) VSIP_NOTHROW
    { subblock_->release(update, pointer); }
  void release(bool update, uT*& pointer) VSIP_NOTHROW
    { subblock_->release(update, pointer); }
  void release(bool update, uT*& real_pointer, uT*& imag_pointer) VSIP_NOTHROW
    { subblock_->release(update, real_pointer, imag_pointer); }

  void find(value_type*& pointer) VSIP_NOTHROW
    { subblock_->find(pointer); }
  void find(uT*& pointer) VSIP_NOTHROW
    { subblock_->find(pointer); }
  void find(uT*& real_pointer, uT*& imag_pointer) VSIP_NOTHROW
    { subblock_->find(real_pointer, imag_pointer); }

  void rebind(value_type* pointer) VSIP_NOTHROW
    { subblock_->rebind(pointer); }
  void rebind(uT* pointer) VSIP_NOTHROW
    { subblock_->rebind(pointer); }
  void rebind(uT* real_pointer, uT* imag_pointer) VSIP_NOTHROW
    { subblock_->rebind(real_pointer, imag_pointer); }

  enum user_storage_type user_storage() const VSIP_NOTHROW
  { return subblock_->user_storage(); }

  bool admitted() const VSIP_NOTHROW
    { return subblock_->admitted(); }

  // Member data.
public:
  map_type              map_;
  processor_type	proc_;			// This processor in comm.
  index_type   		sb_;
  Block*		subblock_;
  length_type	        size_[dim];
};



/// Specialize block layout trait for Distributed_blocks.
template <typename BlockT,
	  typename MapT>
struct Block_layout<Distributed_block<BlockT, MapT> >
{
  static dimension_type const dim = Block_layout<BlockT>::dim;

  typedef typename Block_layout<BlockT>::access_type  access_type;
  typedef typename Block_layout<BlockT>::order_type   order_type;
  typedef typename Block_layout<BlockT>::pack_type    pack_type;
  typedef typename Block_layout<BlockT>::complex_type complex_type;

  typedef Layout<dim, order_type, pack_type, complex_type> layout_type;
};



/// Specialize Distributed_local_block traits class for Distributed_block.
template <typename Block,
	  typename Map>
struct Distributed_local_block<Distributed_block<Block, Map> >
{
  typedef Block type;
};



/// Specialize Is_simple_distributed_block traits class for Distributed_block.
template <typename Block,
	  typename Map>
struct Is_simple_distributed_block<Distributed_block<Block, Map> >
{
  static bool const value = true;
};



#if VSIP_IMPL_USE_GENERIC_VISITOR_TEMPLATES==0

/// Specialize Combine_return_type for Distributed_block leaves.
template <typename CombineT,
	  typename BlockT,
	  typename MapT>
struct Combine_return_type<CombineT, Distributed_block<BlockT, MapT> >
{
  typedef Distributed_block<BlockT, MapT> block_type;
  typedef typename CombineT::template return_type<block_type>::type
		type;
  typedef typename CombineT::template tree_type<block_type>::type
		tree_type;
};



/// Specialize apply_combine for Distributed_block leaves.
template <typename CombineT,
	  typename BlockT,
	  typename MapT>
typename Combine_return_type<CombineT, Distributed_block<BlockT, MapT> >::type
apply_combine(
  CombineT const&                        combine,
  Distributed_block<BlockT, MapT> const& block)
{
  return combine.apply(block);
}

#endif



/***********************************************************************
  Definitions
***********************************************************************/


/// Return the total size of the block.
template <typename Block,
	  typename Map>
inline length_type
Distributed_block<Block, Map>::size()
  const VSIP_NOTHROW
{
  length_type size = 1;
  for (dimension_type d=0; d<dim; ++d)
    size *= size_[d];
  return size;
}



/// Return the size of the block in a specific dimension.
///
/// Parameters:
///   :block_dim: selects which block-dimensionality
///               (BLOCK_DIM == Block::dim).
///   :d: is the dimension whose length to return (0 <= d < block_dim).
///
/// Returns:
///   The size of dimension d.
template <typename Block,
	  typename Map>
inline
length_type
Distributed_block<Block, Map>::size(
  dimension_type block_dim,
  dimension_type d)
  const VSIP_NOTHROW
{
  assert(block_dim == dim);
  assert(d < dim);
  return size_[d];
}



/// Return the local block for a given subblock.

#if 0
// For now, leave this undefined catches unhandled distributed cases at
// compile-time.
template <typename BlockT>
typename Distributed_local_block<BlockT>::type&
get_local_block(
  BlockT const& /*block*/)
{
  // In general case, we should assume block is not distributed and
  // just return it.
  //
  // For now, through exception to catch unhandled distributed cases.
  VSIP_IMPL_THROW(impl::unimplemented("get_local_block()"));
}
#endif



template <typename BlockT,
	  typename MapT>
BlockT&
get_local_block(
  Distributed_block<BlockT, MapT> const& block)
{
  return block.get_local_block();
}



#if 0
// For now, leave this undefined catches unhandled distributed cases at
// compile-time.
template <typename BlockT>
void
assert_local(
  BlockT const& /*block*/,
  index_type    sb)
{
  // In general case, we should assume block is not distributed and
  // just return it.
  //
  // For now, through exception to catch unhandled distributed cases.
  VSIP_IMPL_THROW(impl::unimplemented("assert_local()"));
}
#endif



template <typename BlockT,
	  typename MapT>
void
assert_local(
  Distributed_block<BlockT, MapT> const& block,
  index_type                             sb)
{
  block.assert_local(sb);
}



} // namespace vsip::impl
} // namespace vsip


#endif // VSIP_IMPL_DISTRIBUTED_BLOCK_HPP
