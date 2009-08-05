/* Copyright (c) 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/impl_tags.hpp
    @author  Jules Bergmann
    @date    2006-12-06
    @brief   VSIPL++ Library: Implementation Tags.

*/

#ifndef VSIP_CORE_IMPL_TAGS_HPP
#define VSIP_CORE_IMPL_TAGS_HPP

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// Implementation tags.
//
// Each implementation (generic, BLAS, IPP, etc) has a unique
// implementation tag.
//
// Tags are shared by the different evaluators (Serial-expr_evaluator,
// General_dispatch, and Dispatch).
//
// For the serial expression evaluator, these are placed in
// LibraryTagList and determine the order in which expression
// evaluators are tried.  (The order here approximates LibraryTagList)

struct Intel_ipp_tag {};	// Intel IPP Library
struct Transpose_tag {};	// Optimized Matrix Transpose
struct Mercury_sal_tag {};	// Mercury SAL Library
struct Cbe_sdk_tag {};          // IBM CBE SDK.
struct Cml_tag {};              // IBM Cell Math Library
struct Simd_builtin_tag {};	// Builtin SIMD routines (non loop fusion)
struct Dense_expr_tag {};	// Dense multi-dim expr reduction
struct Copy_tag {};		// Optimized Copy
struct Op_expr_tag {};		// Special expr handling (vmmul, etc)
struct Simd_loop_fusion_tag {};	// SIMD Loop Fusion.
struct Simd_unaligned_loop_fusion_tag {};
struct Fc_expr_tag {};		// Fused Fastconv RBO evaluator.
struct Rbo_expr_tag {};		// Return-block expression evaluator.
struct Mdim_expr_tag {};	// Multi-dim expr reduction
struct Loop_fusion_tag {};	// Generic Loop Fusion (base case).

struct Blas_tag {};		// BLAS implementation (ATLAS, MKL, etc)
struct Lapack_tag {};		// LAPACK implementation (ATLAS, MKL, etc)
struct Generic_tag {};		// Generic implementation.
struct Parallel_tag {};		// Parallel implementation.
struct Cvsip_tag {};		// C-VSIPL library.
struct Cuda_tag {};             // NVidia CUDA GPU library

struct Opt_tag {};		// Optimized Tag.


namespace dispatcher
{

// Operation Tags.
//
// Each operation (dot-product, matrix-matrix product, etc) has a 
// unique operation tag.

struct Op_prod_vv_dot;    // vector-vector dot-product
struct Op_prod_vv_outer;  // vector-vector outer-product
struct Op_prod_mm;        // matrix-matrix product
struct Op_prod_mm_conj;   // matrix-matrix conjugate product
struct Op_prod_mv;        // matrix-vector product
struct Op_prod_vm;        // vector-matrix product
struct Op_prod_gemp;      // generalized matrix-matrix product


} // namespace vsip::impl::dispatcher

} // namespace vsip::impl
} // namespace vsip


#endif // VSIP_CORE_IMPL_TAGS_HPP
