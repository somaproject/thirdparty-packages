/* Copyright (c) 2006, 2007 by CodeSourcery.  All rights reserved. */

/** @file    vsip/opt/expr/return_block.hpp
    @author  Jules Bergmann
    @date    2006-09-02
    @brief   VSIPL++ Library: Expression return block.

*/

#ifndef VSIP_OPT_EXPR_RETURN_BLOCK_HPP
#define VSIP_OPT_EXPR_RETURN_BLOCK_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <memory>

#include <vsip/core/block_traits.hpp>
#include <vsip/core/view_traits.hpp>
#include <vsip/core/domain_utils.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

/// Expression template block for return block optimization.

/// Requires:
///   DIM to be a dimension.
///   T to be a value type.
///   RETURNFUNCTOR to be ...

template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
class Return_expr_block : public Non_assignable
{
public:
  static dimension_type const dim = Dim;

  typedef T                                value_type;
  typedef value_type&                      reference_type;
  typedef value_type const&                const_reference_type;
  typedef typename ReturnFunctor::map_type map_type;

  typedef Dense<Dim, T> block_type;

  Return_expr_block(ReturnFunctor& rf)
    : rf_   (rf)
    , block_()
  {}

  // Necessary to provide copy constructor, because auto_ptr cannot
  // be copied.
  Return_expr_block(Return_expr_block const& rhs)
    : rf_   (rhs.rf_),
      block_()
  {}

  length_type size() const VSIP_NOTHROW { return rf_.size(); }
  length_type size(dimension_type block_dim, dimension_type d)
    const VSIP_NOTHROW
  { return rf_.size(block_dim, d); }


  void increment_count() const VSIP_NOTHROW {}
  void decrement_count() const VSIP_NOTHROW {}
  map_type const& map() const VSIP_NOTHROW { return rf_.map(); }

  value_type get(index_type i) const
  { 
    assert(block_.get());
    return block_->get(i);
  }

  value_type get(index_type i, index_type j) const
  { 
    assert(block_.get());
    return block_->get(i, j);
  }

  value_type get(index_type i, index_type j, index_type k) const
  { 
    assert(block_.get());
    return block_->get(i, j, k);
  }

  template <typename DstBlock>
  void apply(DstBlock& dst) const
  {
    rf_.apply(dst);
  }

  void loop_fusion_init()
  {
    block_.reset(new block_type(block_domain<Dim>(*this)));
    rf_.apply(*(block_.get()));
  }

  void loop_fusion_fini()
  {
    block_.reset(0);
  }

  ReturnFunctor const& functor() const { return rf_; }

private:
  ReturnFunctor             rf_;
  std::auto_ptr<block_type> block_;
};



/// Specialize traits for Return_expr_block.

template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Is_expr_block<Return_expr_block<Dim, T, ReturnFunctor> >
{ static bool const value = true; };

template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct View_block_storage<Return_expr_block<Dim, T, ReturnFunctor> const>
  : By_value_block_storage<Return_expr_block<Dim, T, ReturnFunctor> const>
{};

template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct View_block_storage<Return_expr_block<Dim, T, ReturnFunctor> >
  : By_value_block_storage<Return_expr_block<Dim, T, ReturnFunctor> >
{};

template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Distributed_local_block<Return_expr_block<Dim, T, ReturnFunctor> const>
{
  typedef Return_expr_block<Dim, T, typename ReturnFunctor::local_type> const
          type;
  typedef Return_expr_block<Dim, T, typename ReturnFunctor::local_type> const
          proxy_type;
};



template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Combine_return_type<CombineT,
			   Return_expr_block<Dim, T, ReturnFunctor> const>
{
  typedef typename Combine_return_type<CombineT, ReturnFunctor>::tree_type
          rf_type;
  typedef Return_expr_block<Dim, T, rf_type> const tree_type;
  typedef tree_type type;
};



template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Combine_return_type<CombineT,
			   Return_expr_block<Dim, T, ReturnFunctor> >
{
  typedef typename Combine_return_type<CombineT, ReturnFunctor>::tree_type
          rf_type;
  typedef Return_expr_block<Dim, T, rf_type> const tree_type;
  typedef tree_type type;
};


  
template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
Return_expr_block<Dim, T, typename ReturnFunctor::local_type>
get_local_block(
  Return_expr_block<Dim, T, ReturnFunctor> const& g_blk)
{
  typedef Return_expr_block<Dim, T, typename ReturnFunctor::local_type>
    block_type;
  typename ReturnFunctor::local_type rf_local(g_blk.functor().local());
  block_type l_blk(rf_local);
  return l_blk;
}



template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
typename Combine_return_type<CombineT,
			     Return_expr_block<Dim, T, ReturnFunctor> const>
		::type
apply_combine(
  CombineT const&                                 combine,
  Return_expr_block<Dim, T, ReturnFunctor> const& block)
{
  typedef typename Combine_return_type<
    CombineT,
    Return_expr_block<Dim, T, ReturnFunctor> const>::type
		block_type;

  return block_type(apply_combine(combine, block.functor()));
}



template <typename       VisitorT,
	  dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
void
apply_leaf(
  VisitorT const&                                 visitor,
  Return_expr_block<Dim, T, ReturnFunctor> const& block)
{
  apply_leaf(visitor, block.functor());
}

template <dimension_type MapDim,
	  typename       MapT,
	  dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Is_par_same_map<MapDim, MapT,
		       Return_expr_block<Dim, T, ReturnFunctor> const>
{
  typedef Return_expr_block<Dim, T, ReturnFunctor> const block_type;

  static bool value(MapT const& map, block_type& block)
  {
    return Is_par_same_map<MapDim, MapT, ReturnFunctor>
      ::value(map, block.functor());
  }
};



template <dimension_type Dim,
	  typename       T,
	  typename       ReturnFunctor>
struct Is_par_reorg_ok<Return_expr_block<Dim, T, ReturnFunctor> const>
{
  static bool const value = false;
};

} // namespace vsip::impl

} // namespace vsip

#endif // VSIP_OPT_EXPR_RETURN_BLOCK_HPP
