/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    benchmarks/cuda/copy.cpp
    @author  Don McCoy
    @date    2009-03-10
    @brief   VSIPL++ Library: Benchmark for CUDA memory copy.
*/


/***********************************************************************
  Included Files
***********************************************************************/

#include <iostream>

#include <cuda_runtime.h>

#include <vsip/initfin.hpp>
#include <vsip/support.hpp>
#include <vsip/math.hpp>
#include <vsip/core/profile.hpp>
#include <vsip/opt/cuda/kernels.hpp>

#include <vsip_csl/test.hpp>
#include "loop.hpp"

using namespace vsip;
using namespace vsip_csl;


#define DEBUG 0

/***********************************************************************
  Declarations
***********************************************************************/

// ImplTags
struct Impl_host2dev;
struct Impl_dev2host;
struct Impl_host2host;
struct Impl_dev2shared;
struct Impl_dev2dev;
struct Impl_zeroes2dev;

template <typename T,
          typename ImplTag>
struct t_vcopy;



/***********************************************************************
  Vector copy - host to device
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_host2dev> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_host2dev>"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    Vector<T, dst_block_t>   Z(size, T());
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      vsip::impl::Ext_data<dst_block_t> ext_z(Z.block());
      T const* pA = ext_a.data();
      T* pZ = ext_z.data();

      // Allocate device memory
      T* dev;
      cudaMalloc((void**)&dev, size*sizeof(T));

    
      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
        cudaMemcpy(dev, pA, size*sizeof(T), cudaMemcpyHostToDevice);
      t1.stop();
      time = t1.delta();


      // Pull data back, free device memory.
      cudaMemcpy(pZ, dev, size*sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(dev);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
      if (!equal(Z.get(i), A.get(i)))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << A.get(i) << std::endl
		  << "       got     : " << Z.get(i) << std::endl;
      }
      test_assert(equal(Z.get(i), A.get(i)));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};


/***********************************************************************
  Vector copy - device to host
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_dev2host> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_dev2host>"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    Vector<T, dst_block_t>   Z(size, T());
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      vsip::impl::Ext_data<dst_block_t> ext_z(Z.block());
      T* pA = ext_a.data();
      T* pZ = ext_z.data();

      // Allocate device memory
      T* dev;
      cudaMalloc((void**)&dev, size*sizeof(T));

      // Push data to device memory
      cudaMemcpy(dev, pA, size*sizeof(T), cudaMemcpyHostToDevice);

    
      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
        cudaMemcpy(pZ, dev, size*sizeof(T), cudaMemcpyDeviceToHost);
      t1.stop();
      time = t1.delta();


      // Free device memory.
      cudaFree(dev);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
      if (!equal(Z.get(i), A.get(i)))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << A.get(i) << std::endl
		  << "       got     : " << Z.get(i) << std::endl;
      }
      test_assert(equal(Z.get(i), A.get(i)));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};


/***********************************************************************
  Vector copy - host to host
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_host2host> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_host2host>"; }
  int ops_per_point(length_type)  { return 2; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    Vector<T, dst_block_t>   Z(size, T());
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      vsip::impl::Ext_data<dst_block_t> ext_z(Z.block());
      T* pA = ext_a.data();
      T* pZ = ext_z.data();

      // Allocate device memory
      T* dev;
      cudaMalloc((void**)&dev, size*sizeof(T));


      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
      {
        cudaMemcpy(dev, pA, size*sizeof(T), cudaMemcpyHostToDevice);
        cudaMemcpy(pZ, dev, size*sizeof(T), cudaMemcpyDeviceToHost);
      }
      t1.stop();
      time = t1.delta();


      // Free device memory.
      cudaFree(dev);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
      if (!equal(Z.get(i), A.get(i)))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << A.get(i) << std::endl
		  << "       got     : " << Z.get(i) << std::endl;
      }
      test_assert(equal(Z.get(i), A.get(i)));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};


/***********************************************************************
  Vector copy - device to shared (on-chip fast RAM)
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_dev2shared> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_dev2shared>"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    Vector<T, dst_block_t>   Z(size, T());
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      vsip::impl::Ext_data<dst_block_t> ext_z(Z.block());
      T* pA = ext_a.data();
      T* pZ = ext_z.data();

      // Allocate device memory, copy data from host
      T* dev;
      cudaMalloc((void**)&dev, size*sizeof(T));
      cudaMemcpy(dev, pA, size*sizeof(T), cudaMemcpyHostToDevice);

    
      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
        vsip::impl::cuda::copy_device_to_shared(dev, size);
      t1.stop();
      time = t1.delta();


      // Pull data back, free device memory.
      cudaMemcpy(pZ, dev, size*sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(dev);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
#if DEBUG
      if (!equal(Z.get(i), A.get(i)))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << A.get(i) << std::endl
		  << "       got     : " << Z.get(i) << std::endl;
      }
#endif
      test_assert(equal(Z.get(i), A.get(i)));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};


/***********************************************************************
  Vector copy - device to device
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_dev2dev> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_dev2dev>"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    Vector<T, dst_block_t>   Z(size, T());
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      vsip::impl::Ext_data<dst_block_t> ext_z(Z.block());
      T* pA = ext_a.data();
      T* pZ = ext_z.data();

      // Allocate device memory, copy data from host
      T* dev_a;
      cudaMalloc((void**)&dev_a, size*sizeof(T));
      cudaMemcpy(dev_a, pA, size*sizeof(T), cudaMemcpyHostToDevice);
      T* dev_z;
      cudaMalloc((void**)&dev_z, size*sizeof(T));

    
      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
      {
        vsip::impl::cuda::copy_device_to_device(dev_a, dev_z, size);
        cudaThreadSynchronize();
      }
      t1.stop();
      time = t1.delta();


      // Pull data back, free device memory.
      cudaMemcpy(pZ, dev_z, size*sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(dev_z);
      cudaFree(dev_a);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
#if DEBUG
      if (!equal(Z.get(i), A.get(i)))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << A.get(i) << std::endl
		  << "       got     : " << Z.get(i) << std::endl;
      }
#endif
      test_assert(equal(Z.get(i), A.get(i)));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};


/***********************************************************************
  Vector copy - device fill with zero
***********************************************************************/

template <typename T>
struct t_vcopy<T, Impl_zeroes2dev> : Benchmark_base
{
  typedef Dense<1, T, row1_type>  src_block_t;
  typedef Dense<1, T, row1_type>  dst_block_t;
  
  char const* what() { return "CUDA t_vcopy<..., Impl_zero2dev>"; }
  int ops_per_point(length_type)  { return 1; }
  int riob_per_point(length_type) { return 1*sizeof(T); }
  int wiob_per_point(length_type) { return 1*sizeof(T); }
  int mem_per_point(length_type)  { return 1*sizeof(T); }

  void operator()(length_type size, length_type loop, float& time)
  {
    Vector<T, src_block_t>   A(size);
    for (index_type i = 0; i < A.size(); ++i)
      A.put(i, T(i));

    // Scoping is used to control the lifetime of Ext_data<> objects.  These 
    // must be destroyed before accessing data through the view again.
    {     
      vsip::impl::Ext_data<src_block_t> ext_a(A.block());
      T* pA = ext_a.data();

      // Allocate device memory, copy data from host
      T* dev_a;
      cudaMalloc((void**)&dev_a, size*sizeof(T));
      cudaMemcpy(dev_a, pA, size*sizeof(T), cudaMemcpyHostToDevice);

    
      // Benchmark the operation
      vsip::impl::profile::Timer t1;
      t1.start();
      for (index_type l = 0; l < loop; ++l)
      {
        vsip::impl::cuda::copy_zeroes_to_device(dev_a, size);
        cudaThreadSynchronize();
      }
      t1.stop();
      time = t1.delta();


      // Pull data back, free device memory.
      cudaMemcpy(pA, dev_a, size*sizeof(T), cudaMemcpyDeviceToHost);
      cudaFree(dev_a);
    }

    // validate results
    for (index_type i = 0; i < A.size(); ++i)
    {
#if DEBUG
      if (!equal(A.get(i), T()))
      {
	std::cout << "ERROR: at location " << i << std::endl
		  << "       expected: " << T() << std::endl
		  << "       got     : " << A.get(i) << std::endl;
      }
#endif
      test_assert(equal(A.get(i), T()));
    }
    
  }

  t_vcopy()
  {}

  // Member data.
};




void
defaults(Loop1P&)
{
}



int
test(Loop1P& loop, int what)
{
  typedef float F;


  switch (what)
  {
  case  1: loop(t_vcopy<F, Impl_host2dev>()); break;
  case  2: loop(t_vcopy<F, Impl_dev2host>()); break;
  case  3: loop(t_vcopy<F, Impl_host2host>()); break;
  case  4: loop(t_vcopy<F, Impl_dev2shared>()); break;
  case  5: loop(t_vcopy<F, Impl_dev2dev>()); break;
  case  6: loop(t_vcopy<F, Impl_zeroes2dev>()); break;

  case 0:
    std::cout
      << "CUDA copy -- vector copy\n"
      << "   -1 -- host to device copy\n"
      << "   -2 -- device to host copy\n"
      << "   -3 -- host->device->host copy (A = B)\n"
      << "   -4 -- device to shared copy\n"
      << "   -5 -- device to device copy\n"
      << "   -6 -- device fill with zeroes\n"
      ;

  default:
    return 0;
  }
  return 1;
}
