/* Copyright (c) 2009 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/device_memory.hpp
    @author  Don McCoy
    @date    2009-04-05
    @brief   VSIPL++ Library: Device (on GPU) memory class for interfacing
               with CUDA-allocated memory.
*/

#ifndef VSIP_OPT_CUDA_DEVICE_MEMORY_HPP
#define VSIP_OPT_CUDA_DEVICE_MEMORY_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <complex>

#include <vsip/support.hpp>
#include <vsip/opt/cuda/bindings.hpp>
#include <vsip/opt/cuda/device_storage.hpp>
#include <vsip/opt/cuda/gpu_block.hpp>
#include <vsip/opt/cuda/kernels.hpp>
#include <vsip/opt/expr/serial_evaluator.hpp>


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{

/// Helper class used to allocate/deallocate memory and copy data
/// to and from GPU memory space.
template <typename Block>
class Device_memory
{
  typedef typename Block::value_type  T;
  typedef typename Block_layout<Block>::order_type order_type;

public:
  /// Declares GPU buffer used for input or output (to or from GPU respectively).  
  /// Device memory  is allocated.  Copy from host to device occurs if SYNC_IN
  /// or SYNC_INOUT is used.  Copy from device to host memory occurs upon
  /// destruction if SYNC_OUT or SYNC_INOUT is used.
  Device_memory(Block& src_block, sync_action_type sync = SYNC_INOUT)
    : host_(Ext_data<Block>(src_block).data()), block_(src_block), sync_(sync), 
      storage_(block_.size() * sizeof(T))
  {
    if (sync_ & SYNC_IN)
      //      copy_block<Block::dim, order_type, Block>().
      //	host_to_dev(block_, host_, storage_.data());
      copy_host_to_dev(block_, host_, storage_.data());
  }

  /// Host/GPU memory helper class destructor.  Frees all resources after
  /// copying device data back to the host (if SYNC_OUT or SYNC_INOUT is 
  /// specified in the constructor).
  ~Device_memory()
  {
    if (sync_ & SYNC_OUT)
      //      copy_block<Block::dim, order_type, Block>().
      //	dev_to_host(block_, storage_.data(), host_);
      copy_dev_to_host(block_, storage_.data(), host_);
  }

  /// Returns the buffer address (in GPU memory space).
  T* data() { return storage_.data(); }
    
private:
  T* host_;
  Block& block_;
  sync_action_type sync_;
  Device_storage<T> storage_;
};


/// Specialization for const blocks
template <typename Block>
class Device_memory<Block const>
{
  typedef typename Block::value_type  T;
  typedef typename Block_layout<Block>::order_type order_type;

public:
  /// Declares GPU buffer used for input only (to the GPU).  Device memory 
  /// is allocated and data is copied from the host.
  Device_memory(Block const& src_block)
    : host_(Ext_data<Block>(src_block).data()), block_(src_block), sync_(SYNC_IN),
      storage_(block_.size() * sizeof(T))
  {
    //    copy_block<Block::dim, order_type, Block>().
    //      host_to_dev(block_, host_, storage_.data());
    copy_host_to_dev(block_, host_, storage_.data());
  }
  
  /// Returns the buffer address (in GPU memory space).
  T* data() { return storage_.data(); }
    
private:
  T const* host_;
  Block const& block_;
  sync_action_type sync_;
  Device_storage<T> storage_;
};



/// Specialization for Gpu_block 
template <dimension_type Dim>
class Device_memory<Gpu_block<Dim> >
{
  typedef Gpu_block<Dim> block_type;
  typedef typename block_type::value_type T;

public:
  /// Receives a GPU block to be used for input or output.
  Device_memory(block_type& src_block, sync_action_type sync = SYNC_INOUT)
    : block_(src_block), 
      sync_(sync)
  {}

  /// Host/GPU memory helper class destructor.  Copies device data back to
  /// the host (if SYNC_OUT or SYNC_INOUT is specified in the constructor).
  ~Device_memory()
  {
    if (sync_ & SYNC_OUT)
      block_.device_flush();
  }

  /// Returns the buffer address (in GPU memory space).  Data is first copied 
  /// to the GPU if the host copy has been updated since the block was created.
  T* data() { return block_.device_data(); }
    
private:
  block_type& block_;
  sync_action_type sync_;
};


/// Specialization for Gpu_block const
template <dimension_type Dim>
class Device_memory<Gpu_block<Dim> const>
{
  typedef Gpu_block<Dim> const block_type;
  typedef typename block_type::value_type  T;

public:
  /// Receives a GPU block used for input only.  Data is copied if and
  /// only if the host copy has been updated since the block was created.
  Device_memory(block_type& src_block)
    : block_(src_block), 
      sync_(SYNC_IN)
  {}
  
  /// Host/GPU memory helper class destructor.  No data is copied back.
  ~Device_memory() 
  {}

  /// Returns the buffer address (in GPU memory space).
  T const* data() { return block_.device_data(); }
    
private:
  block_type& block_;
  sync_action_type sync_;
};


} // namespace vsip::impl::cuda
} // namespace vsip::impl
} // namespace vsip

#endif // VSIP_OPT_CUDA_DEVICE_MEMORY_HPP
