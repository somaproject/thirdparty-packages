/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/diag/eval.hpp
    @author  Jules Bergmann
    @date    2006-10-26
    @brief   VSIPL++ Library: Diagnostics for evaluation.
*/

#ifndef VSIP_OPT_DIAG_EVAL_HPP
#define VSIP_OPT_DIAG_EVAL_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/support.hpp>
#include <vsip/core/impl_tags.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>
#include <vsip/opt/expr/serial_dispatch_fwd.hpp>
#include <vsip/core/dispatch_assign_decl.hpp>

#include <iostream>
#include <iomanip>
#include <typeinfo>

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

// Forward decl.

template <typename       TagList,
	  dimension_type Dim,
	  typename       DstBlock,
          typename       SrcBlock>
void
diagnose_eval_list_blk(
  DstBlock&       dst,
  SrcBlock const& src);



namespace diag_detail
{

// Helper class to return the name corresponding to a dispatch tag.

template <typename T> 
struct Dispatch_name
{
  static std::string name() { return "unknown"; }
};

#define VSIP_IMPL_DISPATCH_NAME(TYPE)			\
  template <>						\
  struct Dispatch_name<TYPE> {				\
    static std::string name() { return "" # TYPE; }	\
  };

#define VSIP_IMPL_DISPATCH_NAME_AS(TYPE, ASTYPE)	\
  template <>						\
  struct Dispatch_name<TYPE> {				\
    static std::string name() { return "" # ASTYPE; }	\
  };

VSIP_IMPL_DISPATCH_NAME(Intel_ipp_tag)
VSIP_IMPL_DISPATCH_NAME(Transpose_tag)
VSIP_IMPL_DISPATCH_NAME(Mercury_sal_tag)
VSIP_IMPL_DISPATCH_NAME(Cbe_sdk_tag)
VSIP_IMPL_DISPATCH_NAME(Cuda_tag)
VSIP_IMPL_DISPATCH_NAME(Cml_tag)
VSIP_IMPL_DISPATCH_NAME(Dense_expr_tag)
VSIP_IMPL_DISPATCH_NAME(Copy_tag)
VSIP_IMPL_DISPATCH_NAME(Op_expr_tag)
VSIP_IMPL_DISPATCH_NAME(Simd_builtin_tag)
VSIP_IMPL_DISPATCH_NAME(Simd_loop_fusion_tag)
VSIP_IMPL_DISPATCH_NAME_AS(Simd_unaligned_loop_fusion_tag, Simd_ulf_tag)
VSIP_IMPL_DISPATCH_NAME(Fc_expr_tag)
VSIP_IMPL_DISPATCH_NAME(Rbo_expr_tag)
VSIP_IMPL_DISPATCH_NAME(Mdim_expr_tag)
VSIP_IMPL_DISPATCH_NAME(Loop_fusion_tag)
VSIP_IMPL_DISPATCH_NAME(Cvsip_tag)
VSIP_IMPL_DISPATCH_NAME(Opt_tag)
VSIP_IMPL_DISPATCH_NAME(Generic_tag)

VSIP_IMPL_DISPATCH_NAME(Tag_illegal_mix_of_local_and_global_in_assign)
VSIP_IMPL_DISPATCH_NAME(Tag_serial_expr)
// VSIP_IMPL_DISPATCH_NAME(Tag_par_assign)
VSIP_IMPL_DISPATCH_NAME(Tag_par_expr)
VSIP_IMPL_DISPATCH_NAME(Tag_par_expr_noreorg)



// Represent an unknown type.  Used when a non-view is the RHS of an
// assignment.

struct Unknown_type {};



// Helper class to determine the block type of a view.  Handles non-view
// types (which occur in simple expressions, such as A = 5).

template <dimension_type Dim,
	  typename       T>
struct Block_of
{
  typedef Scalar_block<Dim, T> type;
  static type block(T val) { return type(val); }
};

template <dimension_type Dim,
	  template <typename, typename> class View,
	  typename T,
	  typename BlockT>
struct Block_of<Dim, View<T, BlockT> >
{
  typedef BlockT type;
  static type& block(View<T, BlockT> view) { return view.block(); }
};



// Helper class to conditionally call to Serial_expr_evaluator's
// rt_valid() method, when ct_valid is true.  If ct_valid is false,
// the call may not be valid.

// Primary definition, covers ct_valid == false case.

template <typename SeeT,
	  typename DstBlockT,
	  typename SrcBlockT,
	  bool     CtValid = SeeT::ct_valid>
struct Check_rt_valid
{
  static char const* name()
  {
    return 0;
  }

  static bool rt_valid(DstBlockT&, SrcBlockT const&)
  {
    return false;
  }
};

// Specialization for ct_valid == true case.

template <typename SeeT,
	  typename DstBlockT,
	  typename SrcBlockT>
struct Check_rt_valid<SeeT, DstBlockT, SrcBlockT, true>
{
  static char const* name()
  {
    return SeeT::name();
  }

  static bool rt_valid(DstBlockT& dst, SrcBlockT const& src)
  {
    return SeeT::rt_valid(dst, src);
  }
};



/***********************************************************************
  See_summary --  summarize Serial_expr_eval
***********************************************************************/

// Summarize the evaluation of expression 'dst = src' for dispatch
// tag Tag.  Single line form of diagnose_eval_tag().

template <dimension_type Dim,
	  typename       Tag,
	  typename       DstBlockT,
          typename       SrcBlockT>
struct See_summary
{
  static void exec(
    DstBlockT&       dst,
    SrcBlockT const& src)
  {
    using std::cout;
    using std::setw;
    using std::endl;

    typedef Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Tag>
      see_type;

    bool rt_valid = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
      ::rt_valid(dst, src);

    cout << "  - " << setw(20) << Dispatch_name<Tag>::name()
	 << "  ct: " << setw(5) << (see_type::ct_valid ? "true" : "false")
	 << "  rt: " << setw(5) << (rt_valid ? "true" : "false");

    if (see_type::ct_valid)
    {
      char const* name = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
	::name();
      cout << "  (" << name << ")";
    }
    cout << endl;
  }
};

// Specialization for Transpose_tag

template <typename       DstBlockT,
          typename       SrcBlockT>
struct See_summary<2, Transpose_tag, DstBlockT, SrcBlockT>
{
  typedef Transpose_tag Tag;
  static dimension_type const Dim = 2;

  static void exec(
    DstBlockT&       dst,
    SrcBlockT const& src)
  {
    using std::cout;
    using std::setw;
    using std::endl;

    typedef Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Tag>
      see_type;

    bool rt_valid = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
      ::rt_valid(dst, src);

    cout << "  - " << setw(20) << Dispatch_name<Tag>::name()
	 << "  ct: " << setw(5) << (see_type::ct_valid ? "true" : "false")
	 << "  rt: " << setw(5) << (rt_valid ? "true" : "false");

    if (see_type::ct_valid)
    {
      char const* name = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
	::name();
      cout << "  (" << name << ")";
    }

    cout << " ["
	 << see_type::is_rhs_expr << ", "
	 << see_type::lhs_cost << ", "
	 << see_type::rhs_cost << ", "
	 << see_type::is_rhs_split << ", "
	 << see_type::is_rhs_split << "]";
    cout << endl;
  }
};


// Specialization for Cml_tag

#ifdef VSIP_IMPL_CBE_SDK
template <typename       DstBlockT,
          typename       SrcBlockT>
struct See_summary<2, Cml_tag, DstBlockT, SrcBlockT>
{
  typedef Cml_tag Tag;
  static dimension_type const Dim = 2;

  static void exec(
    DstBlockT&       dst,
    SrcBlockT const& src)
  {
    using std::cout;
    using std::setw;
    using std::endl;

    typedef Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Tag>
      see_type;

    bool rt_valid = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
      ::rt_valid(dst, src);

    cout << "  - " << setw(20) << Dispatch_name<Tag>::name()
	 << "  ct: " << setw(5) << (see_type::ct_valid ? "true" : "false")
	 << "  rt: " << setw(5) << (rt_valid ? "true" : "false");

    if (see_type::ct_valid)
    {
      char const* name = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
	::name();
      cout << "  (" << name << ")";
    }

    cout << " ["
	 << see_type::is_rhs_expr << ", "
	 << see_type::lhs_cost << ", "
	 << see_type::rhs_cost << ", "
	 << see_type::is_rhs_split << ", "
	 << see_type::is_rhs_split << "]";
    cout << endl;
  }
};
#endif



// See_summary specialization for Mdim_expr

template <dimension_type Dim,
	  typename       DstBlockT,
          typename       SrcBlockT,
	  bool           CtValid =
	     Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Mdim_expr_tag>
                                  ::ct_valid>
struct See_summary_mdim_expr
{
  typedef Mdim_expr_tag Tag;
  typedef Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Tag> see_type;

  static void exec(
    DstBlockT&       dst,
    SrcBlockT const& src)
  {
    using std::cout;
    using std::setw;
    using std::endl;

    bool rt_valid = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
      ::rt_valid(dst, src);

    cout << "  - " << setw(20) << Dispatch_name<Tag>::name()
	 << "  ct: " << setw(5) << (see_type::ct_valid ? "true" : "false")
	 << "  rt: " << setw(5) << (rt_valid ? "true" : "false");

    cout << endl;
  }
};



template <dimension_type Dim,
	  typename       DstBlockT,
          typename       SrcBlockT>
struct See_summary_mdim_expr<Dim, DstBlockT, SrcBlockT, true>
{
  typedef Mdim_expr_tag Tag;
  typedef Serial_expr_evaluator<Dim, DstBlockT, SrcBlockT, Tag> see_type;

  static void exec(
    DstBlockT&       dst,
    SrcBlockT const& src)
  {
    using std::cout;
    using std::setw;
    using std::endl;

    bool rt_valid = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
      ::rt_valid(dst, src);

    cout << "  - " << setw(20) << Dispatch_name<Tag>::name()
	 << "  ct: " << setw(5) << (see_type::ct_valid ? "true" : "false")
	 << "  rt: " << setw(5) << (rt_valid ? "true" : "false");

    char const* name = Check_rt_valid<see_type, DstBlockT, SrcBlockT>
	::name();
    cout << "  (" << name << ")";
    cout << endl;
    cout << " ------------\n";
    typedef typename see_type::new_dst_type new_dst_type;
    typedef typename see_type::new_src_type new_src_type;
    new_dst_type n_dst = see_type::diag_helper_dst(dst);
    new_src_type n_src = see_type::diag_helper_src(src);
    diagnose_eval_list_blk<LibraryTagList, 1,
      new_dst_type, new_src_type const>(n_dst, n_src);
    cout << " ------------\n";
  }
};

template <dimension_type Dim,
	  typename       DstBlockT,
          typename       SrcBlockT>
struct See_summary<Dim, Mdim_expr_tag, DstBlockT, SrcBlockT>
  : See_summary_mdim_expr<Dim, DstBlockT, SrcBlockT>
{};



template <dimension_type Dim,
	  typename       Tag,
	  typename       DstBlockT,
          typename       SrcBlockT>
void
see_summary(
  DstBlockT&       dst,
  SrcBlockT const& src)
{
  See_summary<Dim, Tag, DstBlockT, SrcBlockT>::exec(dst, src);
}



// Helper class for diagnose_eval_list() to traverse list of
// dispatch tags.

// Primary definition for list traversal.

template <dimension_type Dim,
	  typename       DstBlock,
	  typename       SrcBlock,
	  typename       TagList,
	  typename       Tag = typename TagList::first,
	  typename       Rest = typename TagList::rest>
struct Diag_eval_list_helper
{
  static void exec(
    DstBlock&       dst,
    SrcBlock const& src)
  {
    see_summary<Dim, Tag, DstBlock, SrcBlock>(dst, src);
    Diag_eval_list_helper<Dim, DstBlock, SrcBlock, Rest>::exec(dst, src);
  }
};

// Specialization for list end.

template <dimension_type Dim,
	  typename       DstBlock,
	  typename       SrcBlock,
	  typename       TagList,
	  typename       Tag>
struct Diag_eval_list_helper<Dim, DstBlock, SrcBlock, TagList, Tag, None_type>
{
  static void exec(
    DstBlock&       dst,
    SrcBlock const& src)
  {
    see_summary<Dim, Tag, DstBlock, SrcBlock>(dst, src);
  }
};

// Specialization for Unknown_type RHS.

template <dimension_type Dim,
	  typename       DstBlock,
	  typename       TagList,
	  typename       Tag,
	  typename       Rest>
struct Diag_eval_list_helper<Dim, DstBlock, Unknown_type, TagList, Tag, Rest>
{
  static void exec(
    DstBlock&           /*dst*/,
    Unknown_type const& /*src*/)
  {
    std::cout << "Delh: unknown type\n";
  }
};



// Helper class for diag_eval_dispatch.

template <dimension_type Dim,
	  typename       DstBlock,
	  typename       SrcBlock,
	  typename       DaTag>
struct Diag_eval_dispatch_helper
{
  static void info(
    DstBlock&       /*dst*/,
    SrcBlock const& /*src*/)
  {
    std::cout << "Diag_eval_dispatch_helper: DaTag not handled" << std::endl
	      << "  DaTag: " << Dispatch_name<DaTag>::name() << std::endl
      ;
  }
};



template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2>
struct Diag_eval_dispatch_helper<Dim, Block1, Block2, Tag_serial_expr>
{
  static void info(
    Block1&       blk1,
    Block2 const& blk2)
  {
    // Equivalent to:
    //   diagnose_eval_list_std(dst, src);
    std::cout << "diagnose_eval_list" << std::endl
	      << "  dst expr: " << typeid(Block1).name() << std::endl
	      << "  src expr: " << typeid(Block2).name() << std::endl;
    Diag_eval_list_helper<Dim, Block1, Block2, vsip::impl::LibraryTagList>
	::exec(blk1, blk2);
  }
};



template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2>
struct Diag_eval_dispatch_helper<Dim, Block1, Block2, Tag_par_expr_noreorg>
{
  typedef typename Block1::map_type map1_type;

  typedef typename View_of_dim<Dim, typename Block1::value_type,
			     Block1>::type dst_type;
  typedef typename View_of_dim<Dim, typename Block2::value_type,
			     Block2>::const_type src_type;

  static void info(
    Block1&       blk1,
    Block2 const& blk2)
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

      std::cout << "LHS and RHS have same map -- local assignment\n";

      // Equivalent to:
      //   diagnose_eval_list_std(dst, src);
      std::cout << "diagnose_eval_list" << std::endl
		<< "  dst expr: " << typeid(stor1_t).name() << std::endl
		<< "  src expr: " << typeid(stor2_t).name() << std::endl;
      Diag_eval_list_helper<Dim, block1_t, block2_t,
	                    vsip::impl::LibraryTagList>
	::exec(l_blk1, l_blk2);
    }
    else
    {
      std::cout << "LHS and RHS have different maps\n";
      std::cout << "error: expr cannot be reorganized\n";
    }
  }
};



template <dimension_type Dim,
	  typename       Block1,
	  typename       Block2>
struct Diag_eval_dispatch_helper<Dim, Block1, Block2, Tag_par_expr>
{
  typedef Dispatch_assign<Dim, Block1, Block2, Tag_par_expr> da_type;

  typedef typename da_type::map1_type map1_type;
  typedef typename da_type::dst_type  dst_type;
  typedef typename da_type::src_type  src_type;

  static void info(
    Block1&       blk1,
    Block2 const& blk2)
  {
    if (Is_par_same_map<Dim, map1_type, Block2>::value(blk1.map(), blk2))
    {
      // Maps are same, no communication required.
      typedef typename Distributed_local_block<Block1>::type block1_t;
      typedef typename Distributed_local_block<Block2>::type block2_t;
      typedef typename View_block_storage<block1_t>::type::equiv_type stor1_t;
      typedef typename View_block_storage<block2_t>::type::equiv_type stor2_t;

      std::cout << "  parallel dim : " << Dim << "  (" << blk1.size(Dim, 0);
      for (dimension_type i=1; i<Dim; ++i)
	std::cout << ", " << blk1.size(Dim, i) ;
      std::cout << ")\n";

      stor1_t l_blk1 = get_local_block(blk1);
      stor2_t l_blk2 = get_local_block(blk2);

      std::cout << "  local dim    : " << Dim << "  (" << l_blk1.size(Dim, 0);
      for (dimension_type i=1; i<Dim; ++i)
	std::cout << ", " << l_blk1.size(Dim, i) ;
      std::cout << ")\n";

      std::cout << "LHS and RHS have same map -- local assignment\n";

      // Equivalent to:
      //   diagnose_eval_list_std(dst, src);
      std::cout << "diagnose_eval_list" << std::endl
		<< "  dst expr: " << typeid(stor1_t).name() << std::endl
		<< "  src expr: " << typeid(stor2_t).name() << std::endl;
      Diag_eval_list_helper<Dim, block1_t, block2_t,
	                    vsip::impl::LibraryTagList>
	::exec(l_blk1, l_blk2);
    }
    else
    {
      std::cout << "LHS and RHS have different maps\n";
      std::cout << "(diagnostics not implemented yet)\n";
    }
  }
};

} // namespace vsip::impl::diag_detail



// Diagnose evaluation of an expression 'dst = src' with dispatch
// tag EvalTag.
//
// Requires
//   EVALTAG to be a dispatch tag for Serial_expr_evaluator.
//   DST to be a LHS view
//   SRC to be a RHS view
//
// Description:
//   Writes diagnostic info on how EvalTag's Serial_expr_evaluator
//   would handle 'dst = src' to cout, including ct_valid and rt_valid.
//
// Example:
//   To determine how the loop fusion evaluator would handle 'A = B + C':
//
//      diagnose_eval_tag<vsip::impl::Loop_fusion_tag>(A, B + C)
//
//   This will produce output like so:
//
//      diagnose_eval_tag:
//        name: Expr_Loop
//        tag: Loop_fusion_tag
//        DstBlockT: ...
//        SrcBlockT: ...
//        ct_valid: true
//        rt_valid: true
//

template <typename       EvalTag,
	  typename       DstViewT,
          typename       SrcViewT>
void
diagnose_eval_tag(
  DstViewT dst,
  SrcViewT src)
{
  using vsip::impl::Serial_expr_evaluator;
  using vsip::impl::Block_layout;
  using vsip::impl::diag_detail::Dispatch_name;

  typedef typename DstViewT::block_type dst_block_type;
  typedef typename SrcViewT::block_type src_block_type;
  
  dimension_type const dim = DstViewT::dim;

  typedef Serial_expr_evaluator<dim, dst_block_type, src_block_type, EvalTag>
    see_type;

  typedef typename Block_layout<dst_block_type>::order_type dst_order_type;

  std::cout << "diagnose_eval_tag:" << std::endl;
  std::cout << "  name: " << see_type::name() << std::endl;
  std::cout << "  tag: " << Dispatch_name<EvalTag>::name() << std::endl;
  std::cout << "  DstBlockT: " << typeid(dst_block_type).name() << std::endl;
  std::cout << "  SrcBlockT: " << typeid(src_block_type).name() << std::endl;
  std::cout << "  ct_valid: " << (see_type::ct_valid ? "true" : "false")
	    << std::endl;
  std::cout << "  rt_valid: "
	    << (see_type::rt_valid(dst.block(), src.block()) ? "true"
		                                             : "false")
	    << std::endl;
}



// Diagnose evaluation of an expression 'dst = src' with a list of
// dispatch tags.
//
// Requires
//   TAGLIST to be a list of dispatch tags for the Serial_expr_evaluator.
//   DST to be a LHS view
//   SRC to be a RHS view
//
// Description:
//   Writes diagnostic info on how each dispatch tag would handle
//   'dst = src' to cout, including ct_valid and rt_valid.
//
// Example:
//   To determine how the standard list of evaluators would handle
//   'A = B + C':
//
//      diagnose_eval_list<vsip::impl::LibraryTagList>(A, B + C)
//
//   This will produce output like so:
//
//      -        Intel_ipp_tag  ct:  true  rt:  true
//      -        Transpose_tag  ct: false  rt: false
//      -      Mercury_sal_tag  ct: false  rt: false
//      -     Simd_builtin_tag  ct:  true  rt:  true
//      -       Dense_expr_tag  ct: false  rt: false
//      -             Copy_tag  ct: false  rt: false
//      -          Op_expr_tag  ct: false  rt: false
//      - Simd_loop_fusion_tag  ct: false  rt: false
//      -      Loop_fusion_tag  ct:  true  rt:  true
//

// Block argument version.

template <typename       TagList,
	  dimension_type Dim,
	  typename       DstBlock,
          typename       SrcBlock>
void
diagnose_eval_list_blk(
  DstBlock&       dst,
  SrcBlock const& src)
{
  using vsip::impl::diag_detail::Diag_eval_list_helper;

  std::cout << "diagnose_eval_list" << std::endl
	    << "  dst expr: " << typeid(DstBlock).name() << std::endl
	    << "  src expr: " << typeid(SrcBlock).name() << std::endl;
  Diag_eval_list_helper<Dim, DstBlock, SrcBlock, TagList>::exec(dst, src);
}



template <typename  TagList,
	  typename  DstViewT,
          typename  SrcViewT>
void
diagnose_eval_list(
  DstViewT dst,
  SrcViewT src)
{
  dimension_type const dim = DstViewT::dim;

  typedef typename diag_detail::Block_of<dim, DstViewT>::type dst_block_type;
  typedef typename diag_detail::Block_of<dim, SrcViewT>::type src_block_type;
  using vsip::impl::diag_detail::Diag_eval_list_helper;

  std::cout << "diagnose_eval_list" << std::endl
	    << "  dst expr: " << typeid(dst_block_type).name() << std::endl
	    << "  src expr: " << typeid(src_block_type).name() << std::endl;
  Diag_eval_list_helper<dim, dst_block_type, src_block_type, TagList>
    ::exec(diag_detail::Block_of<dim, DstViewT>::block(dst),
	   diag_detail::Block_of<dim, SrcViewT>::block(src));
}



// Convenience function for 'diagnose_eval_list<vsip::impl::LibraryTagList>()'

template <typename  DstViewT,
          typename  SrcViewT>
void
diagnose_eval_list_std(
  DstViewT dst,
  SrcViewT src)
{
  diagnose_eval_list<LibraryTagList>(dst, src);
}



// Diagnose Dispatch_assign.

template <typename       DstViewT,
          typename       SrcViewT>
void
diagnose_eval_dispatch(
  DstViewT dst,
  SrcViewT src)
{
  using std::cout;
  using std::endl;
  using std::flush;
  using vsip::impl::diag_detail::Dispatch_name;

  dimension_type const dim = DstViewT::dim;

  typedef typename diag_detail::Block_of<dim, DstViewT>::type dst_block_type;
  typedef typename diag_detail::Block_of<dim, SrcViewT>::type src_block_type;

  typedef Dispatch_assign_helper<dim, dst_block_type, src_block_type, false>
    dah;

  typedef typename dah::type dispatch_type;

  cout << "--------------------------------------------------------\n";
  cout << "diagnose_eval_dispatch:" << std::endl;

  cout << "  dim          : " << dim << "  (" << dst.size(0);
  for (dimension_type i=1; i<dim; ++i)
    cout << ", " << dst.size(i) ;
  cout << ")\n";

  cout << "  DstBlockT    : " << typeid(dst_block_type).name() << endl
       << "  SrcBlockT    : " << typeid(src_block_type).name() << endl
       << "  is_illegal   : " << (dah::is_illegal ? "true" : "false") << endl
       << "  is_rhs_expr  : " << (dah::is_rhs_expr ? "true" : "false") << endl
       << "  is_rhs_simple: " << (dah::is_rhs_simple ? "true" : "false") <<endl
       << "  is_rhs_reorg : " << (dah::is_rhs_reorg ? "true" : "false") << endl
       << "  is_lhs_split : " << (dah::is_lhs_split ? "true" : "false") << endl
       << "  is_rhs_split : " << (dah::is_rhs_split ? "true" : "false") << endl
       << "  lhs_cost     : " << dah::lhs_cost << endl
       << "  rhs_cost     : " << dah::rhs_cost << endl
       << "  TYPE         : " << Dispatch_name<dispatch_type>::name() << endl
    ;
  cout << "--------------------------------------------------------\n";

  diag_detail::Diag_eval_dispatch_helper<dim, dst_block_type, src_block_type,
    dispatch_type>
    ::info(diag_detail::Block_of<dim, DstViewT>::block(dst),
	   diag_detail::Block_of<dim, SrcViewT>::block(src));

  cout << "--------------------------------------------------------\n";
  cout << flush;
}




} // namespace vsip::impl::diag_detail
} // namespace vsip

#endif // VSIP_OPT_DIAG_EVAL_HPP
