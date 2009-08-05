/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/expr/ops_info.hpp
    @author  Jules Bergmann
    @date    2006-08-04
    @brief   VSIPL++ Library: Determine the number of ops per point for
                              an expression template.
*/

#ifndef VSIP_OPT_EXPR_OPS_INFO_HPP
#define VSIP_OPT_EXPR_OPS_INFO_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/metaprogramming.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/expr/scalar_block.hpp>
#include <vsip/core/expr/unary_block.hpp>
#include <vsip/core/expr/binary_block.hpp>
#include <vsip/core/expr/ternary_block.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/coverage.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

// Forward declaration
template <dimension_type VecDim,
	  typename       Block0,
	  typename       Block1>
class Vmmul_expr_block;


/// These generate char tags for given data types, defaulting to int
/// with specializations for common floating point types.  These use
/// BLAS/LAPACK convention.

template <typename T> 
struct Type_name    { static char const value = 'I'; };

#define VSIP_IMPL_TYPE_NAME(T, VALUE)		\
template <>					\
struct Type_name<T > { static char const value = VALUE; };

VSIP_IMPL_TYPE_NAME(float,                'S');
VSIP_IMPL_TYPE_NAME(double,               'D');
VSIP_IMPL_TYPE_NAME(std::complex<float>,  'C');
VSIP_IMPL_TYPE_NAME(std::complex<double>, 'Z');

#undef VSIP_IMPL_TYPE_NAME


template <typename T> 
struct Scalar_type_name    { static char const value = 'i'; };

#define VSIP_IMPL_SCALAR_TYPE_NAME(T, VALUE)	\
template <>					\
struct Scalar_type_name<T > { static char const value = VALUE; };

VSIP_IMPL_SCALAR_TYPE_NAME(float,                's');
VSIP_IMPL_SCALAR_TYPE_NAME(double,               'd');
VSIP_IMPL_SCALAR_TYPE_NAME(std::complex<float>,  'c');
VSIP_IMPL_SCALAR_TYPE_NAME(std::complex<double>, 'z');

#undef VSIP_IMPL_SCALAR_TYPE_NAME



/// Traits classes to determine the ops for a particular operation.

template <template <typename> class UnaryOp,
	  typename                  T1>
struct Unary_op_count
{
  static unsigned const value = 0;
}; 

template <template <typename, 
                    typename> class BinaryOp,
	  typename                  T1,
	  typename                  T2>
struct Binary_op_count
{
  static unsigned const value = 0;
}; 

template <template <typename, 
                    typename, 
                    typename> class TernaryOp,
	  typename                  T1,
	  typename                  T2,
	  typename                  T3>
struct Ternary_op_count
{
  static unsigned const value = 0;
}; 


/// Specializations for Unary types

#define VSIP_IMPL_UNARY_OPS_FUNCTOR(OP, TYPE, VALUE)		\
template <typename T>					\
struct Unary_op_count<OP##_functor, TYPE >		\
{							\
  static unsigned const value = VALUE;			\
}; 

//VSIP_IMPL_UNARY_OPS_FUNCTOR(acos)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(arg)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(asin)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(atan)
VSIP_IMPL_UNARY_OPS_FUNCTOR(bnot,  T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(ceil,  T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(conj,  complex<T>,   1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(cos,   T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(cos,   complex<T>,  12)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(cosh)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(euler)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(exp)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(exp10)
VSIP_IMPL_UNARY_OPS_FUNCTOR(floor, T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(imag,  complex<T>,   0)
VSIP_IMPL_UNARY_OPS_FUNCTOR(lnot,  T,            1)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(log)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(log10)
VSIP_IMPL_UNARY_OPS_FUNCTOR(mag,   T,            0)
VSIP_IMPL_UNARY_OPS_FUNCTOR(mag,   complex<T>,  13)
VSIP_IMPL_UNARY_OPS_FUNCTOR(magsq, T,            3)
VSIP_IMPL_UNARY_OPS_FUNCTOR(neg,   T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(real,  complex<T>,   0)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(recip)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(rsqrt)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sin,   T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sin,   complex<T>,  12)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(sinh)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sq,    T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sq,    complex<T>,   5)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sqrt,  T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(sqrt,  complex<T>,  10)
VSIP_IMPL_UNARY_OPS_FUNCTOR(tan,   T,            1)
VSIP_IMPL_UNARY_OPS_FUNCTOR(tan,   complex<T>,  14)
//VSIP_IMPL_UNARY_OPS_FUNCTOR(tanh)

#undef VSIP_IMPL_UNARY_OPS_FUNCTOR


/// Specializations for Binary types

#define VSIP_IMPL_BINARY_OPS(OP, TYPE1, TYPE2, VALUE)		\
template <typename T1,					\
          typename T2>					\
struct Binary_op_count<OP, TYPE1, TYPE2 >		\
{							\
  static unsigned const value = VALUE;			\
}; 

#define VSIP_IMPL_BINARY_OPS_FUNCTOR(OP, TYPE1, TYPE2, VALUE)	\
        VSIP_IMPL_BINARY_OPS(OP##_functor, TYPE1, TYPE2, VALUE)

VSIP_IMPL_BINARY_OPS(op::Add,  T1,          T2,          1)
VSIP_IMPL_BINARY_OPS(op::Add,  T1,          complex<T2>, 1)
VSIP_IMPL_BINARY_OPS(op::Add,  complex<T1>, T2,          1)
VSIP_IMPL_BINARY_OPS(op::Add,  complex<T1>, complex<T2>, 2)
VSIP_IMPL_BINARY_OPS(op::Sub,  T1,          T2,          1)
VSIP_IMPL_BINARY_OPS(op::Sub,  T1,          complex<T2>, 1)
VSIP_IMPL_BINARY_OPS(op::Sub,  complex<T1>, T2,          1)
VSIP_IMPL_BINARY_OPS(op::Sub,  complex<T1>, complex<T2>, 2)
VSIP_IMPL_BINARY_OPS(op::Mult, T1,          T2,          1)
VSIP_IMPL_BINARY_OPS(op::Mult, T1,          complex<T2>, 2)
VSIP_IMPL_BINARY_OPS(op::Mult, complex<T1>, T2,          2)
VSIP_IMPL_BINARY_OPS(op::Mult, complex<T1>, complex<T2>, 6)
VSIP_IMPL_BINARY_OPS(op::Div,  T1,          T2,          1)
VSIP_IMPL_BINARY_OPS(op::Div,  T1,          complex<T2>, 2)
VSIP_IMPL_BINARY_OPS(op::Div,  complex<T1>, T2,          2)
VSIP_IMPL_BINARY_OPS(op::Div,  complex<T1>, complex<T2>, 6)

VSIP_IMPL_BINARY_OPS_FUNCTOR(add,     T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(add,     T1,          complex<T2>,  1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(add,     complex<T1>, T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(add,     complex<T1>, complex<T2>,  2)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(atan2)
VSIP_IMPL_BINARY_OPS_FUNCTOR(band,    T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(bor,     T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(bxor,    T1,          T2,           1)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(div)
VSIP_IMPL_BINARY_OPS_FUNCTOR(eq,      T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(eq,      complex<T1>, complex<T2>,  2)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(fmod)
VSIP_IMPL_BINARY_OPS_FUNCTOR(ge,      T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(gt,      T1,          T2,           1)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(hypot)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(jmul)
VSIP_IMPL_BINARY_OPS_FUNCTOR(land,    T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(le,      T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(lt,      T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(lor,     T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(lxor,    T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(max,     T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(maxmg,   T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(maxmg,   complex<T1>, complex<T2>, 27)
VSIP_IMPL_BINARY_OPS_FUNCTOR(maxmgsq, T1,          T2,           3)
VSIP_IMPL_BINARY_OPS_FUNCTOR(maxmgsq, complex<T1>, complex<T2>,  7)
VSIP_IMPL_BINARY_OPS_FUNCTOR(min,     T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(minmg,   T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(minmg,   complex<T1>, complex<T2>, 27)
VSIP_IMPL_BINARY_OPS_FUNCTOR(minmgsq, T1,          T2,           3)
VSIP_IMPL_BINARY_OPS_FUNCTOR(minmgsq, complex<T1>, complex<T2>,  7)
VSIP_IMPL_BINARY_OPS_FUNCTOR(ne,      T1,          T2,           1)
VSIP_IMPL_BINARY_OPS_FUNCTOR(ne,      complex<T1>, complex<T2>,  2)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(pow)
//VSIP_IMPL_BINARY_OPS_FUNCTOR(sub)

#undef VSIP_IMPL_BINARY_OPS_FUNCTOR
#undef VSIP_IMPL_BINARY_OPS



/// Specializations for Ternary types

#define VSIP_IMPL_TERNARY_OPS(OP, TYPE1, TYPE2, TYPE3, VALUE)	\
template <typename T1,						\
          typename T2,						\
          typename T3>						\
struct Ternary_op_count<OP, TYPE1, TYPE2, TYPE3 >		\
{								\
  static unsigned const value = VALUE;				\
}; 

#define VSIP_IMPL_TERNARY_OPS_FUNCTOR(OP, TYPE1, TYPE2, TYPE3, VALUE)	\
        VSIP_IMPL_TERNARY_OPS(OP##_functor, TYPE1, TYPE2, TYPE3, VALUE)

// Short synonym for above.
#define VSIP_IMPL_TOF(OP, T1, T2, T3, VALUE) \
    VSIP_IMPL_TERNARY_OPS_FUNCTOR(OP, T1, T2, T3, VALUE)

#define VSIP_IMPL_TERNARY_OPS_RRR(OP, VALUE) \
    VSIP_IMPL_TOF(OP, T1,          T2,          T3,          VALUE)
#define VSIP_IMPL_TERNARY_OPS_RRC(OP, VALUE) \
    VSIP_IMPL_TOF(OP, T1,          T2,          complex<T3>, VALUE)
#define VSIP_IMPL_TERNARY_OPS_RCR(OP, VALUE) \
    VSIP_IMPL_TOF(OP, T1,          complex<T2>, T3,          VALUE)
#define VSIP_IMPL_TERNARY_OPS_RCC(OP, VALUE) \
    VSIP_IMPL_TOF(OP, T1,          complex<T2>, complex<T3>, VALUE)
#define VSIP_IMPL_TERNARY_OPS_CRR(OP, VALUE) \
    VSIP_IMPL_TOF(OP, complex<T1>, T2,          T3,          VALUE)
#define VSIP_IMPL_TERNARY_OPS_CRC(OP, VALUE) \
    VSIP_IMPL_TOF(OP, complex<T1>, T2,          complex<T3>, VALUE)
#define VSIP_IMPL_TERNARY_OPS_CCR(OP, VALUE) \
    VSIP_IMPL_TOF(OP, complex<T1>, complex<T2>, T3,          VALUE)
#define VSIP_IMPL_TERNARY_OPS_CCC(OP, VALUE) \
    VSIP_IMPL_TOF(OP, complex<T1>, complex<T2>, complex<T3>, VALUE)


// The cost for ternary functions is computed by adding the costs for 
// pure real, mixed real-complex and pure complex adds and multiples 
// for the given equation:

//  (t1 + t2) * t3
//                            <  adds  >    <   muls   >
//                            R   M   C     R   M     C
VSIP_IMPL_TERNARY_OPS_RRR(am, 1 + 0 + 0*2 + 1 + 0*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RRC(am, 1 + 0 + 0*2 + 0 + 0*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RCR(am, 0 + 1 + 0*2 + 0 + 1*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RCC(am, 0 + 1 + 0*2 + 0 + 0*2 + 1*6)
VSIP_IMPL_TERNARY_OPS_CRR(am, 0 + 1 + 0*2 + 0 + 1*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_CRC(am, 0 + 1 + 0*2 + 0 + 0*2 + 1*6)
VSIP_IMPL_TERNARY_OPS_CCR(am, 0 + 0 + 1*2 + 0 + 1*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_CCC(am, 0 + 0 + 1*2 + 0 + 0*2 + 1*6)

//  t1 * t2 + (T1(1) - t1) * t3
//                                 <  adds  >    <   muls   >
//                                 R   M   C     R   M     C
VSIP_IMPL_TERNARY_OPS_RRR(expoavg, 2 + 0 + 0*2 + 2 + 0*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RRC(expoavg, 1 + 1 + 0*2 + 1 + 1*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RCR(expoavg, 1 + 1 + 0*2 + 1 + 1*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_RCC(expoavg, 1 + 0 + 1*2 + 0 + 2*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_CRR(expoavg, 0 + 0 + 2*2 + 0 + 2*2 + 0*6)
VSIP_IMPL_TERNARY_OPS_CRC(expoavg, 0 + 0 + 2*2 + 0 + 1*2 + 1*6)
VSIP_IMPL_TERNARY_OPS_CCR(expoavg, 0 + 0 + 2*2 + 0 + 1*2 + 1*6)
VSIP_IMPL_TERNARY_OPS_CCC(expoavg, 0 + 0 + 2*2 + 0 + 0*2 + 2*6)

//VSIP_IMPL_TERNARY_OPS_FUNCTOR(ma)
//VSIP_IMPL_TERNARY_OPS_FUNCTOR(msb)
//VSIP_IMPL_TERNARY_OPS_FUNCTOR(sbm)
//VSIP_IMPL_TERNARY_OPS_FUNCTOR(ite)

#undef VSIP_IMPL_TERNARY_OPS_RRR
#undef VSIP_IMPL_TERNARY_OPS_RRC
#undef VSIP_IMPL_TERNARY_OPS_RCR
#undef VSIP_IMPL_TERNARY_OPS_RCC
#undef VSIP_IMPL_TERNARY_OPS_CRR
#undef VSIP_IMPL_TERNARY_OPS_CRC
#undef VSIP_IMPL_TERNARY_OPS_CCR
#undef VSIP_IMPL_TERNARY_OPS_CCC
#undef VSIP_IMPL_TOF
#undef VSIP_IMPL_TERNARY_OPS_FUNCTOR
#undef VSIP_IMPL_TERNARY_OPS


/// Reduction to count the number operations per point of an expression.

struct Reduce_expr_ops_per_point
{
public:
  template <typename BlockT>
  struct leaf_node
  {
    typedef Int_type<0> type;
  };

  template <dimension_type Dim0,
	    typename       T>
  struct leaf_node<Scalar_block<Dim0, T> >
  {
    typedef Int_type<0> type;
  };

  template <dimension_type            Dim0,
	    template <typename> class Op,
	    typename                  NewBlockT,
	    typename                  NewT>
  struct unary_node
  {
    typedef Int_type<Unary_op_count<Op, NewT>::value +
                     NewBlockT::value> type;
  };

  template <dimension_type                Dim0,
	    template <typename, typename> class Op,
	    typename                      NewLBlock,
	    typename                      NewLType,
	    typename                      NewRBlock,
	    typename                      NewRType>
  struct binary_node
  {
    typedef Int_type<Binary_op_count<Op, NewLType, NewRType>::value +
                     NewLBlock::value +
                     NewRBlock::value> type;
  };

  template <dimension_type                          Dim0,
	    template <typename, typename, typename> class Op,
	    typename                                NewBlock1,
	    typename                                NewType1,
	    typename                                NewBlock2,
	    typename                                NewType2,
	    typename                                NewBlock3,
	    typename                                NewType3>
  struct ternary_node
  {
    typedef Int_type<
      Ternary_op_count<Op, NewType1, NewType2, NewType3>::value +
      NewBlock1::value +
      NewBlock2::value +
      NewBlock3::value> type;
  };

  template <typename BlockT>
  struct transform
  {
    typedef typename leaf_node<BlockT>::type type;
  };

  template <typename BlockT>
  struct transform<BlockT const> : public transform<BlockT> {};

  template <dimension_type            Dim0,
	    template <typename> class Op,
	    typename                  BlockT,
	    typename                  T>
  struct transform<Unary_expr_block<Dim0, Op, BlockT, T> >
  {
    typedef typename unary_node<Dim0, Op,
				typename transform<BlockT>::type,
				T>::type type;
  };

  template <dimension_type                Dim0,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  struct transform<Binary_expr_block<Dim0, Op, LBlock, LType,
				     RBlock, RType> >
  {
    typedef typename binary_node<Dim0, Op,
				typename transform<LBlock>::type, LType,
				typename transform<RBlock>::type, RType>
				::type type;
  };

  template <dimension_type                Dim0,
	    typename                      LBlock,
	    typename                      RBlock>
  struct transform<Vmmul_expr_block<Dim0, LBlock, RBlock> >
  {
    typedef typename binary_node<Dim0, op::Mult,
              typename transform<LBlock>::type, typename LBlock::value_type,
              typename transform<RBlock>::type, typename RBlock::value_type>
              ::type type;
  };

  template <dimension_type                Dim0,
	    template <typename, typename, typename> class Op,
	    typename                      Block1,
	    typename                      Type1,
	    typename                      Block2,
	    typename                      Type2,
	    typename                      Block3,
	    typename                      Type3>
  struct transform<Ternary_expr_block<Dim0, Op, Block1, Type1,
				     Block2, Type2, Block3, Type3> >
  {
    typedef typename ternary_node<Dim0, Op,
				typename transform<Block1>::type, Type1,
				typename transform<Block2>::type, Type2,
				typename transform<Block3>::type, Type3>
				::type type;
  };
};


/// This generates the total number of operations per point in a given 
/// expression.  It also computes the total number of points, as this
/// information is needed to calculate the total number of operations.

template <typename BlockT>
struct Expr_ops_per_point
{
  static length_type size(BlockT const& src)
  {
    length_type size = src.size(BlockT::dim, 0);
    if ( BlockT::dim > 1 )
      size *= src.size(BlockT::dim, 1);
    if ( BlockT::dim > 2 )
      size *= src.size(BlockT::dim, 2);
    return size;
  }

  static unsigned const value =
    Reduce_expr_ops_per_point::template transform<BlockT>::type::value;
};




/// Reduction to generate a tag for the entire expression tree

struct Reduce_expr_op_name
{
public:

  template <typename BlockT>
  struct transform
  {
    static std::string tag() 
    {
      std::string st;
      st = Type_name<typename BlockT::value_type>::value;
      return st;
    }
  };

  template <typename BlockT>
  struct transform<BlockT const> : public transform<BlockT> 
  {};

  template <dimension_type            Dim,
	    typename                  T>
  struct transform<Scalar_block<Dim, T> >
  {
    static std::string tag() 
    {
      std::string st;
      st = Scalar_type_name<T>::value;
      return st;
    }
  };

  template <dimension_type            Dim0,
	    template <typename> class Op,
	    typename                  Block,
	    typename                  Type>
  struct transform<Unary_expr_block<Dim0, Op, 
                                    Block, Type> >
  {
    static std::string tag()
    {
      return Op<Type>::name() + std::string("(") 
        + transform<Block>::tag() + std::string(")");
    } 
  };

  template <dimension_type                Dim0,
	    template <typename, typename> class Op,
	    typename                      LBlock,
	    typename                      LType,
	    typename                      RBlock,
	    typename                      RType>
  struct transform<Binary_expr_block<Dim0, Op, 
                                     LBlock, LType,
                                     RBlock, RType> >
  {
    static std::string tag()
    {
      return Op<LType, RType>::name() + std::string("(")
        + transform<LBlock>::tag() + std::string(",")
        + transform<RBlock>::tag() + std::string(")"); 
    } 
  };

  template <dimension_type                Dim0,
	    typename                      LBlock,
	    typename                      RBlock>
  struct transform<Vmmul_expr_block<Dim0, LBlock, RBlock> >
  {
    static std::string tag()
    {
      return std::string("vmmul") + std::string("(")
        + transform<LBlock>::tag() + std::string(",")
        + transform<RBlock>::tag() + std::string(")"); 
    } 
  };

  template <dimension_type                                Dim0,
	    template <typename, typename, typename> class Op,
	    typename                                      Block1,
	    typename                                      Type1,
	    typename                                      Block2,
	    typename                                      Type2,
	    typename                                      Block3,
	    typename                                      Type3>
  struct transform<Ternary_expr_block<Dim0, Op, 
                                     Block1, Type1,
                                     Block2, Type2,
                                     Block3, Type3> >
  {
    static std::string tag()
    {
      return Op<Type1, Type2, Type3>::name() + std::string("(")
        + transform<Block1>::tag() + std::string(",")
        + transform<Block2>::tag() + std::string(",")
        + transform<Block3>::tag() + std::string(")"); 

    } 
  };

};


/// This generates a tag for an expression in standard prefix notation
/// where the operator is shown, followed by the list of operands in
/// parenthesis.  The operator may be one of  the common binary 
/// operators +-*/ or simply the name of the function.  User-defined
/// expression evaluators will use one of 'unary', 'binary' or 'ternary' 
/// for the function name.  The operand will be one of the letters 'S', 
/// 'D', 'C', and 'Z' for views of those types (using the BLAS convention).
/// Scalar operands use lower-case equivalents of the same letters.
/// For example for matrices of single-precision values where A and B 
/// are real and C is complex:
///
///   A * B + C    -->  +(*(S,S),C)
///   A * 5.f + C  -->  +(*(S,s),C)
///   A * (B + C)  -->  *(S,+(S,C))

template <typename EvalExpr,
          typename BlockT>
struct Expr_op_name
{
  static std::string tag(BlockT const& src)
  {
#if VSIP_IMPL_PROFILE_USE_TYPEID
    return typeid(BlockT).name();
#else
    std::ostringstream  tag;
    tag << EvalExpr::name() << " "
        << BlockT::dim << "D "
        << Reduce_expr_op_name::template transform<BlockT>::tag() << " "
        << src.size(BlockT::dim, 0);
    if ( BlockT::dim > 1 )
      tag << "x" << src.size(BlockT::dim, 1);
    if ( BlockT::dim > 2 )
      tag << "x" << src.size(BlockT::dim, 2);

    return tag.str();
#endif
  }
};


} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_IMPL_EXPR_OPS_INFO_HPP
