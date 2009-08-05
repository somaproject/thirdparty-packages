/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/matvec_prod.hpp
    @author  Jules Bergmann
    @date    2005-09-12
    @brief   VSIPL++ Library: Matrix-Vector product operations

*/

#ifndef VSIP_CORE_MATVEC_PROD_HPP
#define VSIP_CORE_MATVEC_PROD_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/vector.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/matvec.hpp>
#if VSIP_IMPL_CVSIP_FFT
# include <vsip/core/cvsip/matvec.hpp>
#endif
#ifndef VSIP_IMPL_REF_IMPL
# ifdef VSIP_IMPL_CBE_SDK
#  include <vsip/opt/cbe/cml/matvec.hpp>
# endif
# ifdef VSIP_IMPL_HAVE_BLAS
#  include <vsip/opt/lapack/matvec.hpp>
# endif
# ifdef VSIP_IMPL_HAVE_SAL
#  include <vsip/opt/sal/eval_misc.hpp>
# endif
#endif



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace dispatcher
{

#ifndef VSIP_IMPL_REF_IMPL
template<>
struct List<Op_prod_mm>
{
  typedef Make_type_list<Cml_tag, Blas_tag, Mercury_sal_tag, 
    Cvsip_tag, Generic_tag>::type type;
};
#endif

/// Generic evaluator for matrix-matrix products.
template <typename Block0,
	  typename Block1,
	  typename Block2>
struct Evaluator<Op_prod_mm, Generic_tag, 
                 void(Block0&, Block1 const&, Block2 const&)>
{
  static bool const ct_valid = true;
  static bool rt_valid(Block0&, Block1 const&, Block2 const&)
  { return true; }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    typedef typename Block0::value_type RT;

    for (index_type i=0; i<r.size(2, 0); ++i)
      for (index_type j=0; j<r.size(2, 1); ++j)
      {
	RT sum = RT();
	for (index_type k=0; k<a.size(2, 1); ++k)
	{
	  sum += a.get(i, k) * b.get(k, j);
	}
	r.put(i, j, sum);
    }
  }
};


#ifndef VSIP_IMPL_REF_IMPL
template<>
struct List<Op_prod_mm_conj>
{
  typedef Make_type_list<Cml_tag, Blas_tag, Mercury_sal_tag, 
    Cvsip_tag, Generic_tag>::type type;
};
#endif

/// Generic evaluator for matrix-matrix conjugate products.
template <typename Block0,
	  typename Block1,
	  typename Block2>
struct Evaluator<Op_prod_mm_conj, Generic_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  static bool const ct_valid = true;
  static bool rt_valid(Block0&, Block1 const&, Block2 const&)
  { return true; }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    typedef typename Block0::value_type RT;

    for (index_type i=0; i<r.size(2, 0); ++i)
      for (index_type j=0; j<r.size(2, 1); ++j)
      {
	RT sum = RT();
	for (index_type k=0; k<a.size(2, 1); ++k)
	{
	  sum += a.get(i, k) * conj(b.get(k, j));
	}
	r.put(i, j, sum);
    }
  }
};


#ifndef VSIP_IMPL_REF_IMPL
template<>
struct List<Op_prod_mv>
{
  typedef Make_type_list<Cml_tag, Blas_tag, Mercury_sal_tag, 
    Cvsip_tag, Generic_tag>::type type;
};
#endif

/// Generic evaluator for matrix-vector products.
template <typename Block0,
	  typename Block1,
	  typename Block2>
struct Evaluator<Op_prod_mv, Generic_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  static bool const ct_valid = true;
  static bool rt_valid(Block0&, Block1 const&, Block2 const&)
  { return true; }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    typedef typename Block0::value_type RT;

    for (index_type i=0; i<r.size(1, 0); ++i)
    {
      RT sum = RT();
      for (index_type k=0; k<a.size(2, 1); ++k)
      {
        sum += a.get(i, k) * b.get(k);
      }
      r.put(i, sum);
    }
  }
};


#ifndef VSIP_IMPL_REF_IMPL
template<>
struct List<Op_prod_vm>
{
  typedef Make_type_list<Cml_tag, Blas_tag, Mercury_sal_tag, 
    Cvsip_tag, Generic_tag>::type type;
};
#endif

/// Generic evaluator for vector-matrix products.
template <typename Block0,
	  typename Block1,
	  typename Block2>
struct Evaluator<Op_prod_vm, Generic_tag,
                 void(Block0&, Block1 const&, Block2 const&)>
{
  static bool const ct_valid = true;
  static bool rt_valid(Block0&, Block1 const&, Block2 const&)
  { return true; }

  static void exec(Block0& r, Block1 const& a, Block2 const& b)
  {
    typedef typename Block0::value_type RT;

    for (index_type i=0; i<r.size(); ++i)
    {
      RT sum = RT();
      for (index_type k=0; k<b.size(2, 0); ++k)
      {
        sum += a.get(k) * b.get(k, i);
      }
      r.put(i, sum);
    }
  }
};

} // namespace vsip::impl::dispatcher


/// Generic matrix-matrix product.
template <typename T0,
	  typename T1,
	  typename T2,
	  typename Block0,
	  typename Block1,
	  typename Block2>
void
generic_prod(
  const_Matrix<T0, Block0> a,
  const_Matrix<T1, Block1> b,
  Matrix<T2, Block2>       r)
{
  assert(r.size(0) == a.size(0));
  assert(r.size(1) == b.size(1));
  assert(a.size(1) == b.size(0));

#ifdef VSIP_IMPL_REF_IMPL
  dispatcher::Evaluator<dispatcher::Op_prod_mm, Cvsip_tag,
    void(Block2&, Block0 const&, Block1 const&)>::exec
    (r.block(), a.block(), b.block());
#else
  impl::dispatch<dispatcher::Op_prod_mm, void,
    Block2&, Block0 const&, Block1 const&>
    (r.block(), a.block(), b.block());
#endif
}


/// Generic matrix-matrix conjugate product.
template <typename T0,
	  typename T1,
	  typename T2,
	  typename Block0,
	  typename Block1,
	  typename Block2>
void
generic_prodj(
  const_Matrix<T0, Block0> a,
  const_Matrix<T1, Block1> b,
  Matrix<T2, Block2>       r)
{
  assert(r.size(0) == a.size(0));
  assert(r.size(1) == b.size(1));
  assert(a.size(1) == b.size(0));

#ifdef VSIP_IMPL_REF_IMPL
  impl::generic_prod(a, conj(b), r);
#else
  impl::dispatch<dispatcher::Op_prod_mm_conj, void,
    Block2&, Block0 const&, Block1 const&>
    (r.block(), a.block(), b.block());
#endif
}


/// Generic matrix-vector product.
template <typename T0,
	  typename T1,
	  typename T2,
	  typename Block0,
	  typename Block1,
	  typename Block2>
void
generic_prod(
  const_Matrix<T0, Block0> a,
  const_Vector<T1, Block1> b,
  Vector<T2, Block2>       r)
{
  assert(r.size() == a.size(0));
  assert(a.size(1) == b.size());

#ifdef VSIP_IMPL_REF_IMPL
  dispatcher::Evaluator<dispatcher::Op_prod_mv, Cvsip_tag,
    void(Block2&, Block0 const&, Block1 const&)>::exec
    (r.block(), a.block(), b.block());
#else
  impl::dispatch<dispatcher::Op_prod_mv, void,
    Block2&, Block0 const&, Block1 const&>
    (r.block(), a.block(), b.block());
#endif
}


/// Generic vector-matrix product.
template <typename T0,
	  typename T1,
	  typename T2,
	  typename Block0,
	  typename Block1,
	  typename Block2>
void
generic_prod(
  const_Vector<T0, Block0> a,
  const_Matrix<T1, Block1> b,
  Vector<T2, Block2>       r)
{
  assert(r.size() == b.size(1));
  assert(a.size() == b.size(0));

#ifdef VSIP_IMPL_REF_IMPL
  dispatcher::Evaluator<dispatcher::Op_prod_vm, Cvsip_tag,
    void(Block2&, Block0 const&, Block1 const&)>::exec
    (r.block(), a.block(), b.block());
#else
  impl::dispatch<dispatcher::Op_prod_vm, void,
    Block2&, Block0 const&, Block1 const&>
    (r.block(), a.block(), b.block());
#endif
}

} // namespace vsip::impl



/// Matrix-matrix product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Matrix<typename Promotion<T0, T1>::type>
prod(
  const_Matrix<T0, Block0> a,
  const_Matrix<T1, Block1> b)
{
  typedef typename Promotion<T0, T1>::type return_type;

  Matrix<return_type> r(a.size(0), b.size(1));

  impl::generic_prod(a, b, r);

  return r;
}


/// Matrix-vector product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Vector<typename Promotion<T0, T1>::type>
prod(
  const_Matrix<T0, Block0> a,
  const_Vector<T1, Block1> b)
{
  typedef typename Promotion<T0, T1>::type return_type;

  Vector<return_type> r(a.size(0));

  impl::generic_prod(a, b, r);

  return r;
}


/// Vector-Matrix product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Vector<typename Promotion<T0, T1>::type>
prod(
  const_Vector<T0, Block0> a,
  const_Matrix<T1, Block1> b)
{
  typedef typename Promotion<T0, T1>::type return_type;

  Vector<return_type> r(b.size(1));

  impl::generic_prod(a, b, r);

  return r;
}


/// [3x3] Matrix-matrix product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Matrix<typename Promotion<T0, T1>::type>
prod3(
  const_Matrix<T0, Block0> a,
  const_Matrix<T1, Block1> b)
{
  assert( a.size(0) == 3 );
  assert( a.size(1) == 3 );
  typedef typename Promotion<T0, T1>::type return_type;

  Matrix<return_type> r(a.size(0), b.size(1));

  impl::generic_prod(a, b, r);

  return r;
}


/// [3x3] Matrix-vector product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Vector<typename Promotion<T0, T1>::type>
prod3(
  const_Matrix<T0, Block0> a,
  const_Vector<T1, Block1> b)
{
  assert( a.size(0) == 3 );
  assert( a.size(1) == 3 );
  typedef typename Promotion<T0, T1>::type return_type;

  Vector<return_type> r(a.size(0));

  impl::generic_prod(a, b, r);

  return r;
}


/// [4x4] Matrix-matrix product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Matrix<typename Promotion<T0, T1>::type>
prod4(
  const_Matrix<T0, Block0> a,
  const_Matrix<T1, Block1> b)
{
  assert( a.size(0) == 4 );
  assert( a.size(1) == 4 );
  typedef typename Promotion<T0, T1>::type return_type;

  Matrix<return_type> r(a.size(0), b.size(1));

  impl::generic_prod(a, b, r);

  return r;
}


/// [4x4] Matrix-vector product dispatch.
template <typename T0,
	  typename T1,
	  typename Block0,
	  typename Block1>
const_Vector<typename Promotion<T0, T1>::type>
prod4(
  const_Matrix<T0, Block0> a,
  const_Vector<T1, Block1> b)
{
  assert( a.size(0) == 4 );
  assert( a.size(1) == 4 );
  typedef typename Promotion<T0, T1>::type return_type;

  Vector<return_type> r(a.size(0));

  impl::generic_prod(a, b, r);

  return r;
}


/// Matrix-Matrix (with hermitian) product dispatch.
template <typename T0,
          typename T1,
          typename Block0,
          typename Block1>
const_Matrix<typename Promotion<complex<T0>, complex<T1> >::type>
prodh(
  const_Matrix<complex<T0>, Block0> m0,
  const_Matrix<complex<T1>, Block1> m1) 
    VSIP_NOTHROW
{
  typedef typename Promotion<complex<T0>, complex<T1> >::type return_type;

  Matrix<return_type> r(m0.size(0), m1.size(0));

  impl::generic_prodj(m0, trans(m1), r);

  return r;
}


/// Matrix-Matrix (with complex conjugate) product dispatch.
template <typename T0,
          typename T1,
          typename Block0,
          typename Block1>
const_Matrix<typename Promotion<complex<T0>, complex<T1> >::type>
prodj(
  const_Matrix<complex<T0>, Block0> m0,
  const_Matrix<complex<T1>, Block1> m1)
    VSIP_NOTHROW
{
  typedef typename Promotion<complex<T0>, complex<T1> >::type return_type;

  Matrix<return_type> r(m0.size(0), m1.size(1));
  
  impl::generic_prodj(m0, m1, r);

  return r;
}


/// Matrix-Matrix (with transpose) product dispatch.
template <typename T0,
          typename T1,
          typename Block0,
          typename Block1>
const_Matrix<typename Promotion<T0, T1>::type>
prodt(
  const_Matrix<T0, Block0> m0,
  const_Matrix<T1, Block1> m1)
    VSIP_NOTHROW
{
  typedef typename Promotion<T0, T1>::type return_type;

  Matrix<return_type> r(m0.size(0), m1.size(0));

  impl::generic_prod(m0, trans(m1), r);

  return r;
}


} // namespace vsip

#endif // VSIP_IMPL_MATVEC_PROD_HPP
