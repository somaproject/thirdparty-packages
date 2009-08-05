/* Copyright (c) 2005, 2008 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/setup_assign.hpp
    @author  Jules Bergmann
    @date    2005-08-26
    @brief   VSIPL++ Library: Early binding of an assignment.

*/

#ifndef VSIP_CORE_SETUP_ASSIGN_HPP
#define VSIP_CORE_SETUP_ASSIGN_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/noncopyable.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/parallel/expr.hpp>
#include <vsip/core/parallel/assign_chain.hpp>
#include <vsip/core/dispatch_assign.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/profile.hpp>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

class Setup_assign
   : impl::Non_copyable
{
private:
  class Holder_base
  {
  public:
    virtual ~Holder_base() {}
    virtual void exec() = 0;
    virtual char const * type() = 0;
  };

  class Null_holder : public Holder_base
  {
  public:
    Null_holder() {}
    ~Null_holder() {}
    void exec() {}
    char const * type() { return "Null_holder"; }
  };

  template <dimension_type Dim,
	    typename       DstBlock,
	    typename       SrcBlock>
  class Par_expr_holder : public Holder_base
  {
    typedef typename DstBlock::value_type value1_type;
    typedef typename SrcBlock::value_type value2_type;
  public:
    Par_expr_holder(
      typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst,
      typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src)
      : par_expr_(dst, src)
      {}

    ~Par_expr_holder()
      {}

    void exec() { par_expr_();}
    char const * type() { return "Par_expr_holder"; }


    // Member data
  private:
    vsip::impl::Par_expr<Dim, DstBlock, SrcBlock> par_expr_;
  };



  template <dimension_type Dim,
	    typename       DstBlock,
	    typename       SrcBlock,
	    typename       ParAssignImpl>
  class Par_assign_holder : public Holder_base
  {
    typedef typename DstBlock::value_type value1_type;
    typedef typename SrcBlock::value_type value2_type;
  public:
    Par_assign_holder(
      typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst,
      typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src)
      : par_assign_(dst, src)
      {}

    ~Par_assign_holder()
      {}

    void exec() { par_assign_();}
    char const * type() { return "Par_assign_holder"; }


    // Member data
  private:
    vsip::impl::Par_assign<Dim, value1_type, value2_type,
			   DstBlock, SrcBlock, ParAssignImpl>
		par_assign_;
  };



  template <dimension_type Dim,
	    typename       DstBlock,
	    typename       SrcBlock>
  class Simple_par_expr_holder : public Holder_base
  {
    typedef typename DstBlock::value_type value1_type;
    typedef typename SrcBlock::value_type value2_type;
  public:
    Simple_par_expr_holder(
      typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst,
      typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src)
      : dst_(dst), src_(src)
      {}

    ~Simple_par_expr_holder()
      {}

    void exec() { par_expr_simple(dst_, src_);}
    char const * type() { return "Simple_par_expr_holder"; }

    // Member data
  private:
    typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst_;
    typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src_;
  };



  template <dimension_type Dim,
	    typename       DstBlock,
	    typename       SrcBlock>
  class Ser_expr_holder : public Holder_base
  {
    typedef typename DstBlock::value_type value1_type;
    typedef typename SrcBlock::value_type value2_type;
  public:
    Ser_expr_holder(
      typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst,
      typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src)
      : dst_(dst), src_(src)
      {}

    ~Ser_expr_holder()
      {}

    void exec() { dst_ = src_;}
    char const * type() { return "Ser_expr_holder"; }


    // Member data
  private:
    typename impl::View_of_dim<Dim, value1_type, DstBlock>::type       dst_;
    typename impl::View_of_dim<Dim, value2_type, SrcBlock>::const_type src_;
  };

  template <dimension_type Dim,
	    typename       View1,
	    typename       View2>
  void
  create_holder(View1 dst, View2 src, impl::Tag_serial_expr)
  {
    typedef typename View1::block_type block1_type;
    typedef typename View2::block_type block2_type;
    holder_ = new Ser_expr_holder<Dim, block1_type, block2_type>(dst, src);
  }

  template <dimension_type Dim,
	    typename       View1,
	    typename       View2,
	    typename       ParAssignImpl>
  void
  create_holder(View1 dst, View2 src, impl::Tag_par_assign<ParAssignImpl>)
  {
    typedef typename View1::value_type           value1_type;
    typedef typename View2::value_type           value2_type;
    typedef typename View1::block_type           block1_type;
    typedef typename View2::block_type           block2_type;
    typedef typename View1::block_type::map_type map1_type;
    typedef typename View2::block_type::map_type map2_type;

    if (impl::Is_par_same_map<Dim, map1_type, block2_type>::value(
					dst.block().map(),
					src.block()))
    {
      typedef typename impl::Distributed_local_block<block1_type>::type
	dst_lblock_type;
      typedef typename impl::Distributed_local_block<block2_type>::type
	src_lblock_type;

      typedef typename impl::View_of_dim<Dim, value1_type, dst_lblock_type>::type
	l_dst_view_type;
      typedef typename impl::View_of_dim<Dim, value2_type, src_lblock_type>::type
	l_src_view_type;

      l_dst_view_type l_dst = get_local_view(dst);
      l_src_view_type l_src = get_local_view(src);

      holder_ = new Ser_expr_holder<Dim, dst_lblock_type, src_lblock_type>
	(l_dst, l_src);
    }
    else
    {
      holder_ = new Par_assign_holder<Dim, block1_type, block2_type, ParAssignImpl>(dst, src);
    }
  }

  template <dimension_type Dim,
	    typename       View1,
	    typename       View2>
  void
  create_holder(View1 dst, View2 src, impl::Tag_par_expr)
  {
    typedef typename View1::value_type           value1_type;
    typedef typename View2::value_type           value2_type;
    typedef typename View1::block_type           block1_type;
    typedef typename View2::block_type           block2_type;
    typedef typename View1::block_type::map_type map1_type;
    typedef typename View2::block_type::map_type map2_type;

    if (impl::Is_par_same_map<Dim, map1_type, block2_type>::value(
					dst.block().map(),
					src.block()))
    {
      typedef typename impl::Distributed_local_block<block1_type>::type
	dst_lblock_type;
      typedef typename impl::Distributed_local_block<block2_type>::type
	src_lblock_type;

      typedef typename impl::View_of_dim<Dim, value1_type, dst_lblock_type>::type
	l_dst_view_type;
      typedef typename impl::View_of_dim<Dim, value2_type, src_lblock_type>::type
	l_src_view_type;

      l_dst_view_type l_dst = get_local_view(dst);
      l_src_view_type l_src = get_local_view(src);

      holder_ = new Ser_expr_holder<Dim, dst_lblock_type, src_lblock_type>
	(l_dst, l_src);
    }
    else
    {
      holder_ = new Par_expr_holder<Dim, block1_type, block2_type>(dst, src);
    }
  }


  // Constructors.
public:
   template <template <typename, typename> class View1,
	     template <typename, typename> class View2,
	     typename                            T1,
	     typename                            Block1,
	     typename                            T2,
	     typename                            Block2>
  Setup_assign(
    View1<T1, Block1> dst,
    View2<T2, Block2> src)
  {
    using vsip::impl::ITE_Type;
    using vsip::impl::As_type;
    using vsip::impl::Type_equal;

    dimension_type const dim = View1<T1, Block1>::dim;

    typedef typename Block1::map_type map1_type;
    typedef typename Block2::map_type map2_type;

    typedef typename
      impl::Dispatch_assign_helper<dim, Block1, Block2, true>::type
      dispatch_type;

    create_holder<dim>(dst, src, dispatch_type());
  }

  ~Setup_assign() 
  { delete holder_; }

  void operator()()
  { 
    impl::profile::Scope<impl::profile::par> scope(impl_type());
    holder_->exec();
  }

  char const * impl_type()
  {
    return holder_->type();
  }
  
// Member Data
private:
  Holder_base* holder_;

};

} // namespace vsip

#endif // VSIP_CORE_SETUP_ASSIGN_HPP
