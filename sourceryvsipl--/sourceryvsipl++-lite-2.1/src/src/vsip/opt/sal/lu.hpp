/* Copyright (c) 2006 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/sal/lu.hpp
    @author  Assem Salama
    @date    2006-04-04
    @brief   VSIPL++ Library: LU linear system solver using SAL.

*/

#ifndef VSIP_IMPL_SAL_SOLVER_LU_HPP
#define VSIP_IMPL_SAL_SOLVER_LU_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/support.hpp>
#include <vsip/matrix.hpp>
#include <vsip/core/math_enum.hpp>
#include <vsip/core/temp_buffer.hpp>
#include <vsip/core/working_view.hpp>
#include <vsip/core/fns_elementwise.hpp>
#include <vsip/core/solver/common.hpp>

#include <sal.h>

// This switch chooses between two sets of LUD-related functions provided 
// by SAL.  Setting to '1' will select the newer mat_lud_sol/dec() variants 
// and setting it to '0' will select the older matlud() and matfbs() pair.

#define VSIP_IMPL_SAL_USE_MAT_LUD  1



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

// A structure that tells us if sal lud impl is available
// for certain types
template <>
struct Is_lud_impl_avail<Mercury_sal_tag, float>
{
  static bool const value = true;
};

template <>
struct Is_lud_impl_avail<Mercury_sal_tag, complex<float> >
{
  static bool const value = true;
};

#if !VSIP_IMPL_SAL_USE_MAT_LUD
template <>
struct Is_lud_impl_avail<Mercury_sal_tag, double>
{
  static bool const value = true;
};

template <>
struct Is_lud_impl_avail<Mercury_sal_tag, complex<double> >
{
  static bool const value = true;
};
#endif

// SAL LUD decomposition functions
#define VSIP_IMPL_SAL_LUD_DEC( T, D_T, SAL_T, SALFCN ) \
inline bool                          \
sal_mat_lud_dec(                     \
  T *c, int ctcols,                  \
  D_T *d, int n)                     \
{                                    \
  return (SALFCN(                     \
   (SAL_T*) c, ctcols,               \
   (D_T*) d, n,                      \
   NULL, NULL, SAL_COND_EST_NONE) == SAL_SUCCESS); \
}

#define VSIP_IMPL_SAL_LUD_DEC_SPLIT( T, D_T, SAL_T, SALFCN ) \
inline bool                          \
sal_mat_lud_dec(                     \
  std::pair<T*,T*> c, int ctcols,    \
  D_T *d, int n)                     \
{                                    \
  return (SALFCN(                    \
   (SAL_T*) &c, ctcols,              \
   (D_T*) d, n,                      \
   NULL, NULL, SAL_COND_EST_NONE) == SAL_SUCCESS); \
}
// Declare LUD decomposition functions
VSIP_IMPL_SAL_LUD_DEC(float,                   int,float,        mat_lud_dec)
VSIP_IMPL_SAL_LUD_DEC(complex<float>,          int,COMPLEX,      cmat_lud_dec)
VSIP_IMPL_SAL_LUD_DEC_SPLIT(float,             int,COMPLEX_SPLIT,zmat_lud_dec)

#undef VSIP_IMPL_SAL_LUD_DEC
#undef VSIP_IMPL_SAL_LUD_DEC_SPLIT


// SAL LUD solver functions
#define VSIP_IMPL_SAL_LUD_SOL( T, D_T, SAL_T, SALFCN ) \
inline bool                          \
sal_mat_lud_sol(                     \
  T *a, int atcols,                  \
  D_T *d, T *b, T *w,                \
  int n,int flag)                    \
{                                    \
  return (SALFCN(                    \
   (SAL_T*) a, atcols,               \
   (D_T*) d,(SAL_T*)b,(SAL_T*)w,     \
   n,flag) == SAL_SUCCESS);          \
}

#define VSIP_IMPL_SAL_LUD_SOL_SPLIT( T, D_T, SAL_T, SALFCN ) \
inline bool                                        \
sal_mat_lud_sol(                                   \
  std::pair<T*,T*> a, int atcols,                  \
  D_T *d, std::pair<T*,T*> b, std::pair<T*,T*> w,  \
  int n,int flag)                                  \
{                                                  \
  return (SALFCN(                                  \
   (SAL_T*) &a, atcols,                            \
   (D_T*) d,(SAL_T*)&b,(SAL_T*)&w,                 \
   n,flag) == SAL_SUCCESS);                        \
}
// Declare LUD solver functions
VSIP_IMPL_SAL_LUD_SOL(float,         int,float,        mat_lud_sol)
VSIP_IMPL_SAL_LUD_SOL(complex<float>,int,COMPLEX,      cmat_lud_sol)
VSIP_IMPL_SAL_LUD_SOL_SPLIT(float,   int,COMPLEX_SPLIT,zmat_lud_sol)

#undef VSIP_IMPL_SAL_LUD_SOL
#undef VSIP_IMPL_SAL_LUD_SOL_SPLIT



// "Legacy" SAL functions - The single-precision versions are listed
// in the Appendix of the SAL Reference manual.  Although the double-
// precision ones are still part of the normal API, we refer to both 
// sets of functions as legacy functions just for ease of naming.

// This function is not provided by SAL but is similar to
// vrecip() which works on floats.
inline void 
vrecipd( double *A, int I, double *C, int K, int N )
{
  while ( N-- )
  {
    *C = 1.0 / *A;
    A += I;
    C += K;
  }
}

// Legacy SAL LUD decomposition functions
// Note that the stride may be passed to the reciprocal function,
// however, the decomposition functions work only with unit strides.

#define VSIP_IMPL_SAL_LUD_DEC( T, SAL_T, SALRECP, SALFCN ) \
inline void                  \
sal_matlud(                  \
  T *r,                      \
  T *c,                      \
  int *d, int n)             \
{                            \
  SALFCN((SAL_T*) c, d, n);  \
  SALRECP((SAL_T*) c, n+1,   \
          (SAL_T*) r, 1, n); \
}

#define VSIP_IMPL_SAL_LUD_DEC_CPLX( T, SAL_T, SALRECP, SALFCN ) \
inline void                    \
sal_matlud(                    \
  T *r,                        \
  T *c,                        \
  int *d, int n)               \
{                              \
  SALFCN((SAL_T*) c, d, n);    \
  SALRECP((SAL_T*) c, 2*(n+1), \
          (SAL_T*) r, 2, n);   \
}

#define VSIP_IMPL_SAL_LUD_DEC_SPLIT( T, SAL_T, SALRECP, SALFCN ) \
inline void                   \
sal_matlud(                   \
  std::pair<T*,T*> r,         \
  std::pair<T*,T*> c,         \
  int *p, int n)              \
{                             \
  SALFCN((SAL_T*) &c, p, n);  \
  SALRECP((SAL_T*) &c, n+1,   \
          (SAL_T*) &r, 1, n); \
}
// Declare LUD decomposition functions
VSIP_IMPL_SAL_LUD_DEC(float,           float,                vrecip,  matlud)
VSIP_IMPL_SAL_LUD_DEC(double,          double,               vrecipd, matludd)
VSIP_IMPL_SAL_LUD_DEC(complex<float>,  COMPLEX,              cvrcip,  cmatlud)
VSIP_IMPL_SAL_LUD_DEC(complex<double>, DOUBLE_COMPLEX,       cvrcipd, cmatludd)
VSIP_IMPL_SAL_LUD_DEC_SPLIT(float,     COMPLEX_SPLIT,        zvrcip,  zmatlud)
VSIP_IMPL_SAL_LUD_DEC_SPLIT(double,    DOUBLE_COMPLEX_SPLIT, zvrcipd, zmatludd)

#undef VSIP_IMPL_SAL_LUD_DEC
#undef VSIP_IMPL_SAL_LUD_DEC_CPLX
#undef VSIP_IMPL_SAL_LUD_DEC_SPLIT


// Legacy LUD solver functions
// Note that the stride may be passed when using complex types, 
// but not with scalar types.  As a result, complex types are 
// restricted to "dense" equivalents (2 for complex interleaved
// and 1 for complex split in SAL terms).

#define VSIP_IMPL_SAL_LUD_SOL( T, SAL_T, SALFCN ) \
inline void                           \
sal_matfbs(                           \
  T *a, T *b, int *p,                 \
  T *c, T *d, int n )                 \
{                                     \
  SALFCN( (SAL_T*) a, (SAL_T*) b, p,  \
          (SAL_T*) c, (SAL_T*) d, n); \
}
#define VSIP_IMPL_SAL_LUD_SOL_CPLX( T, SAL_T, SALFCN ) \
inline void                              \
sal_matfbs(                              \
  T *a, T *b, int *p,                    \
  T *c, T *d, int n )                    \
{                                        \
  SALFCN( (SAL_T*) a, (SAL_T*) b, p,     \
          (SAL_T*) c, (SAL_T*) d, 2, n); \
}
#define VSIP_IMPL_SAL_LUD_SOL_SPLIT( T, SAL_T, SALFCN ) \
inline void                                        \
sal_matfbs(                                        \
  std::pair<T*,T*> a, std::pair<T*,T*> b, int *p,  \
  std::pair<T*,T*> c, std::pair<T*,T*> d, int n )  \
{                                                  \
  SALFCN( (SAL_T*) &a, (SAL_T*) &b, p,             \
          (SAL_T*) &c, (SAL_T*) &d, 1, n);         \
}
// Declare LUD solver functions
VSIP_IMPL_SAL_LUD_SOL(float,                float,                matfbs)
VSIP_IMPL_SAL_LUD_SOL(double,               double,               matfbsd)
VSIP_IMPL_SAL_LUD_SOL_CPLX(complex<float>,  COMPLEX,              cmatfbs)
VSIP_IMPL_SAL_LUD_SOL_CPLX(complex<double>, DOUBLE_COMPLEX,       cmatfbsd)
VSIP_IMPL_SAL_LUD_SOL_SPLIT(float,          COMPLEX_SPLIT,        zmatfbs)
VSIP_IMPL_SAL_LUD_SOL_SPLIT(double,         DOUBLE_COMPLEX_SPLIT, zmatfbsd)

#undef VSIP_IMPL_SAL_LUD_SOL
#undef VSIP_IMPL_SAL_LUD_SOL_CPLX
#undef VSIP_IMPL_SAL_LUD_SOL_SPLIT




/// LU factorization implementation class.  Common functionality
/// for lud by-value and by-reference classes. SAL only supports floats. There
/// are specializations of lud for double farther bellow

template <typename T>
class Lud_impl<T,Mercury_sal_tag>
{
  // The input matrix must be in ROW major form. We want the b matrix
  // to be in COL major form because we need to call the SAL function for
  // each column in the b matrix. SAL supports split and interleaved complex
  // formats. Complex_type tells us which format we will end up using.
  typedef vsip::impl::dense_complex_type complex_type;
  typedef Storage<complex_type, T>       storage_type;
  typedef typename storage_type::type    ptr_type;

  typedef Layout<2, row2_type, Stride_unit_dense, complex_type> data_LP;
  typedef Fast_block<2, T, data_LP> data_block_type;

  typedef Layout<2, col2_type, Stride_unit_dense, complex_type> b_data_LP;
  typedef Fast_block<2, T, b_data_LP> b_data_block_type;

  typedef Dense<1, T> reciprocals_block_type;

  // Constructors, copies, assignments, and destructors.
public:
  Lud_impl(length_type)
    VSIP_THROW((std::bad_alloc));
  Lud_impl(Lud_impl const&)
    VSIP_THROW((std::bad_alloc));

  Lud_impl& operator=(Lud_impl const&) VSIP_NOTHROW;
  ~Lud_impl() VSIP_NOTHROW;

  // Accessors.
public:
  length_type length()const VSIP_NOTHROW { return length_; }

  // Solve systems.
public:
  template <typename Block>
  bool decompose(Matrix<T, Block>) VSIP_NOTHROW;

protected:
  template <mat_op_type tr,
            typename    Block0,
            typename    Block1>
  bool impl_solve(const_Matrix<T, Block0>, Matrix<T, Block1>)
    VSIP_NOTHROW;

  length_type max_decompose_size();

  // Member data.
private:
  typedef std::vector<int, Aligned_allocator<int> > vector_type;

  length_type  length_;                      // Order of A.
  vector_type  ipiv_;                        // Pivot table for Q. This gets
                                             // generated from the decompose and
                                             // gets used in the solve
  Vector<T, reciprocals_block_type> recip_;  // Vector of reciprocals used
                                             // with legacy solvers
  Matrix<T, data_block_type> data_;          // Factorized matrix (A)
};




/***********************************************************************
  Definitions
***********************************************************************/


template <typename T>
Lud_impl<T,Mercury_sal_tag>::Lud_impl(
  length_type length
  )
VSIP_THROW((std::bad_alloc))
  : length_ (length),
    ipiv_   (length_),
    recip_  (length_),
    data_   (length_, length_)
{
  assert(length_ > 0);
}



template <typename T>
Lud_impl<T,Mercury_sal_tag>::Lud_impl(Lud_impl const& lu)
VSIP_THROW((std::bad_alloc))
  : length_ (lu.length_),
    ipiv_   (length_),
    recip_  (length_),
    data_   (length_, length_)
{
  data_ = lu.data_;
  for (index_type i=0; i<length_; ++i)
    ipiv_[i] = lu.ipiv_[i];
  recip_ = lu.recip_;
}



template <typename T>
Lud_impl<T,Mercury_sal_tag>::~Lud_impl()
  VSIP_NOTHROW
{
}


template <typename T>
length_type
Lud_impl<T, Mercury_sal_tag>::max_decompose_size()
{
  return (impl::Is_complex<T>::value ? 512 : 1024);
}



/// Form LU factorization of matrix A
///
/// Requires
///   A to be a square matrix, either
///
/// FLOPS:
///   real   : UPDATE
///   complex: UPDATE

template <typename T>
template <typename Block>
bool
Lud_impl<T,Mercury_sal_tag>::decompose(Matrix<T, Block> m)
  VSIP_NOTHROW
{
  assert(m.size(0) == length_ && m.size(1) == length_);

  assign_local(data_, m);

  Ext_data<data_block_type> a_ext(data_.block());
  Ext_data<reciprocals_block_type> r_ext(recip_.block());
  bool success;

  if(length_ > 1) 
  {
#if VSIP_IMPL_SAL_USE_MAT_LUD
    success = sal_mat_lud_dec(
                     a_ext.data(),a_ext.stride(0),
                     &ipiv_[0], length_);
#else
    if (length_ > max_decompose_size())
      VSIP_IMPL_THROW(unimplemented(
        "Lud_impl<T,Mercury_sal_tag>::decompose - exceeds maximum size"));
    success = true;
    sal_matlud(r_ext.data(), a_ext.data(), &ipiv_[0], length_);
#endif
  }
  else 
  {
    success = true;
  }
  return success;
}


/// Solve Op(A) x = b (where A previously given to decompose)
///
/// Op(A) is
///   A   if tr == mat_ntrans
///   A^T if tr == mat_trans
///   A'  if tr == mat_herm (valid for T complex only)
///
/// Requires
///   B to be a (length, P) matrix
///   X to be a (length, P) matrix
///
/// Effects:
///   X contains solution to Op(A) X = B

template <typename T>
template <mat_op_type tr,
          typename    Block0,
          typename    Block1>
bool
Lud_impl<T,Mercury_sal_tag>::impl_solve(
  const_Matrix<T, Block0> b,
  Matrix<T, Block1>       x)
  VSIP_NOTHROW
{
  assert(b.size(0) == length_);
  assert(b.size(0) == x.size(0) && b.size(1) == x.size(1));

  int trans;
  // We want X matrix to be same layout as B matrix to make it easier to
  // store result.
  Matrix<T, b_data_block_type> b_int(b.size(0),b.size(1));// local copy of b
  Matrix<T, b_data_block_type> x_int(b.size(0),b.size(1));// local copy of x
  Matrix<T, data_block_type> data_int(length_,length_);   // local copy of data

  assign_local(b_int, b);

  if (tr == mat_ntrans)
    trans = SAL_NORMAL_SOLVER;
  else if (tr == mat_trans)
    trans = SAL_TRANSPOSE_SOLVER;
  else if (tr == mat_herm)
  {
    assert(Is_complex<T>::value);
    trans = SAL_TRANSPOSE_SOLVER;
  }

  if(length_ > 1) 
  {
    Ext_data<b_data_block_type> b_ext(b_int.block());
    Ext_data<b_data_block_type> x_ext(x_int.block());
    if(tr == mat_trans) 
    {
      assign_local(data_int,data_);
      data_int=impl::impl_conj(data_int);
    }
    Ext_data<data_block_type>   a_ext((tr == mat_trans)?
                                        data_int.block():data_.block());
    Ext_data<reciprocals_block_type>  r_ext(recip_.block());

    // sal_mat_lud_sol only takes vectors, so, we have to do this for each
    // column in the matrix
    ptr_type b_ptr = b_ext.data();
    ptr_type x_ptr = x_ext.data();
    for(index_type i=0;i<b.size(1);i++) 
    {
#if VSIP_IMPL_SAL_USE_MAT_LUD
      sal_mat_lud_sol(a_ext.data(), a_ext.stride(0),
                      &ipiv_[0],
                      storage_type::offset(b_ptr,i*length_),
                      storage_type::offset(x_ptr,i*length_),
                      length_,trans);
#else
      if (tr == mat_ntrans)
        sal_matfbs(a_ext.data(), r_ext.data(), &ipiv_[0],
                   storage_type::offset(b_ptr, i*length_),
                   storage_type::offset(x_ptr, i*length_),
                   length_);
      else
        VSIP_IMPL_THROW(unimplemented(
          "Lud_impl<mat_op_type!=mat_ntrans>::impl_solve - unimplemented"));
#endif
    }

    assign_local(x, x_int);
  }
  else 
  {
    for(index_type i=0;i<b.size(1);i++)
      if(tr == mat_herm) 
      {
        T result = b.get(0,i)/impl::impl_conj(data_.get(0,0));
        x.put(0,i,result);
      }
      else
        x.put(0,i,b.get(0,i)/data_.get(0,0));
  }


  return true;
}



} // namespace vsip::impl
} // namespace vsip


#endif // VSIP_IMPL_SAL_SOLVER_LU_HPP
