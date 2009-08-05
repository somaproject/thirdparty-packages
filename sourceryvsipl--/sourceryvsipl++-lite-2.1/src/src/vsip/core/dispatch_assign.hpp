/* Copyright (c) 2005, 2006 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/dispatch_assign.hpp
    @author  Jules Bergmann
    @date    2005-05-10
    @brief   VSIPL++ Library: Assignment dispatch.

*/

#ifndef VSIP_CORE_DISPATCH_ASSIGN_HPP
#define VSIP_CORE_DISPATCH_ASSIGN_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/static_assert.hpp>
#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/parallel/map_traits.hpp>
#include <vsip/core/parallel/expr.hpp>
#ifndef VSIP_IMPL_REF_IMPL
# include <vsip/opt/expr/serial_dispatch.hpp>
#endif
#include <vsip/core/parallel/assign.hpp>
#include <vsip/core/dispatch_assign_decl.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// Apply a small meta-program to determine the proper tag for handling
// an assignment:
//
//    if (LHS and RHS are serial)
//      use Tag_serial
//    else if (LHS and RHS are simple distributed)
//      use Tag_par_assign
//    else // (LHS and RHS are distributed expression)
//      use Tag_par_expr

template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2,
	  bool           EarlyBinding>
struct Dispatch_assign_helper
{
  typedef typename Block1::map_type map1_type;
  typedef typename Block2::map_type map2_type;

  // Cannot mix local and distributed data in expressions.
  static bool const is_illegal    =
    !((Is_local_map<map1_type>::value && Is_local_map<map2_type>::value) ||
      (Is_global_map<map1_type>::value && Is_global_map<map2_type>::value));

  static bool const is_local      = Is_local_map<map1_type>::value &&
                                    Is_local_map<map2_type>::value;
  static bool const is_rhs_expr   = Is_expr_block<Block2>::value;
  static bool const is_rhs_simple = Is_simple_distributed_block<Block2>::value;
  static bool const is_rhs_reorg  = Is_par_reorg_ok<Block2>::value;

  static bool const is_lhs_split  = Is_split_block<Block1>::value;
  static bool const is_rhs_split  = Is_split_block<Block2>::value;

  static int const  lhs_cost      = Ext_data_cost<Block1>::value;
  static int const  rhs_cost      = Ext_data_cost<Block2>::value;

  typedef typename
    Choose_par_assign_impl<Dim, Block1, Block2, EarlyBinding>::type
    par_assign_type;

  typedef typename
  ITE_Type<is_illegal,          As_type<Tag_illegal_mix_of_local_and_global_in_assign>,
  ITE_Type<is_local,            As_type<Tag_serial_expr>,
  ITE_Type<is_rhs_simple,       As_type<Tag_par_assign<par_assign_type> >,
  ITE_Type<is_rhs_reorg,        As_type<Tag_par_expr>,
	                        As_type<Tag_par_expr_noreorg> > > > >
		::type type;
};



// Specialization for serial, 1-dimensional assignment.

template <typename       Block1,
	  typename       Block2>
struct Dispatch_assign<1, Block1, Block2, Tag_serial_expr>
{
  typedef typename Block_layout<Block1>::layout_type LP1;
  typedef typename Block_layout<Block2>::layout_type LP2;

  static void exec(Block1& dst, Block2 const& src)
  {
#ifdef VSIP_IMPL_REF_IMPL
    length_type const size = dst.size(1, 0);
    for (index_type i=0; i<size; ++i) dst.put(i, src.get(i));
#else
    Serial_dispatch<1, Block1, Block2, LibraryTagList>::exec(dst, src);
#endif
  }
};



// Specialization for serial, 2-dimensional assignment.

template <typename       Block1,
	  typename       Block2>
struct Dispatch_assign<2, Block1, Block2, Tag_serial_expr>
{
  static void exec(Block1& dst, Block2 const& src)
  {
#ifdef VSIP_IMPL_REF_IMPL
    length_type const rows = dst.size(2, 0);
    length_type const cols = dst.size(2, 1);
    for (index_type i=0; i<rows; ++i)
      for (index_type j=0; j<cols; ++j)
	dst.put(i, j, src.get(i, j));
#else
    Serial_dispatch<2, Block1, Block2, LibraryTagList>::exec(dst, src);
#endif
  }
};



// Specialization for serial, 3-dimensional assignment.

template <typename       Block1,
	  typename       Block2>
struct Dispatch_assign<3, Block1, Block2, Tag_serial_expr>
{
  static void exec(Block1& dst, Block2 const& src)
  {
#ifdef VSIP_IMPL_REF_IMPL
    length_type const size0 = dst.size(3, 0);
    length_type const size1 = dst.size(3, 1);
    length_type const size2 = dst.size(3, 2);

    for (index_type i=0; i<size0; ++i)
      for (index_type j=0; j<size1; ++j)
        for (index_type k=0; k<size2; ++k)
          dst.put(i, j, k, src.get(i, j, k));
#else
    Serial_dispatch<3, Block1, Block2, LibraryTagList>::exec(dst, src);
#endif
  }
};
  


// Specialization for parallel assignment, RHS is simple (A = B)

template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2,
	  typename       ParAssignImpl>
struct Dispatch_assign<Dim, Block1, Block2, Tag_par_assign<ParAssignImpl> >
{
  typedef typename Block1::map_type map1_type;

  typedef typename View_of_dim<Dim, typename Block1::value_type,
			     Block1>::type dst_type;
  typedef typename View_of_dim<Dim, typename Block2::value_type,
			     Block2>::const_type src_type;

  static void exec(Block1& blk1, Block2 const& blk2)
  {
    if (Is_par_same_map<Dim, map1_type, Block2>::value(blk1.map(), blk2))
    {
      // Maps are same, no communication required.
      typedef typename Distributed_local_block<Block1>::type block1_t;
      typedef typename Distributed_local_block<Block2>::type block2_t;
      typedef typename View_block_storage<block1_t>::type::equiv_type stor1_t;
      typedef typename View_block_storage<block2_t>::type::equiv_type stor2_t;

      stor1_t l_blk1 = get_local_block(blk1);
      stor2_t l_blk2 = get_local_block(blk2);

      Dispatch_assign<Dim, block1_t, block2_t>::exec(l_blk1, l_blk2);
    }
    else
    {
      dst_type dst(blk1);
      src_type src(const_cast<Block2&>(blk2));
      
      Par_assign<Dim,
	typename Block1::value_type,
	typename Block2::value_type,
	Block1,
	Block2,
	ParAssignImpl> pa(dst, src);

      pa();
    }
  }
};



// Specialization for distributed expressions where the RHS can be
// reorganized.

template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2>
struct Dispatch_assign<Dim, Block1, Block2, Tag_par_expr>
{
  typedef typename Block1::map_type map1_type;

  typedef typename View_of_dim<Dim, typename Block1::value_type,
			     Block1>::type dst_type;
  typedef typename View_of_dim<Dim, typename Block2::value_type,
			     Block2>::const_type src_type;

  static void exec(Block1& blk1, Block2 const& blk2)
  {
    if (Is_par_same_map<Dim, map1_type, Block2>::value(blk1.map(), blk2))
    {
      // Maps are same, no communication required.
      typedef typename Distributed_local_block<Block1>::type block1_t;
      typedef typename Distributed_local_block<Block2>::type block2_t;
      typedef typename View_block_storage<block1_t>::type::equiv_type stor1_t;
      typedef typename View_block_storage<block2_t>::type::equiv_type stor2_t;

      stor1_t l_blk1 = get_local_block(blk1);
      stor2_t l_blk2 = get_local_block(blk2);

      Dispatch_assign<Dim, block1_t, block2_t>::exec(l_blk1, l_blk2);
    }
    else
    {
      // Maps are different, fall out to general expression.
      dst_type dst(blk1);
      src_type src(const_cast<Block2&>(blk2));
      par_expr(dst, src);
    }
  }
};



// Specialization for distributed expressions that cannot be reorganized.
// 
// Types of expressions that cannot be reorganized:
//  - Return_expr_blocks
//  - Vmmul_expr_blocks

template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2>
struct Dispatch_assign<Dim, Block1, Block2, Tag_par_expr_noreorg>
{
  typedef typename Block1::map_type map1_type;

  typedef typename View_of_dim<Dim, typename Block1::value_type,
			     Block1>::type dst_type;
  typedef typename View_of_dim<Dim, typename Block2::value_type,
			     Block2>::const_type src_type;

  static void exec(Block1& blk1, Block2 const& blk2)
  {
    if (Is_par_same_map<Dim, map1_type, Block2>::value(blk1.map(), blk2))
    {
      // Maps are same, no communication required.
      typedef typename Distributed_local_block<Block1>::type block1_t;
      typedef typename Distributed_local_block<Block2>::type block2_t;
      typedef typename View_block_storage<block1_t>::type::equiv_type stor1_t;
      typedef typename View_block_storage<block2_t>::type::equiv_type stor2_t;

      stor1_t l_blk1 = get_local_block(blk1);
      stor2_t l_blk2 = get_local_block(blk2);

      Dispatch_assign<Dim, block1_t, block2_t>::exec(l_blk1, l_blk2);
    }
    else
    {
      VSIP_IMPL_THROW(impl::unimplemented("Expression cannot be reorganized"));
    }
  }
};

template <dimension_type D, typename DstBlock, typename SrcBlock>
inline void 
assign(DstBlock& dst, SrcBlock& src)
{
  Dispatch_assign<D, DstBlock, SrcBlock>::exec(dst, src);
}

} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_CORE_DISPATCH_ASSIGN_HPP
