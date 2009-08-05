/* Copyright (c) 2009 by CodeSourcery.  All rights reserved. 

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cuda/gpu_block.hpp
    @author  Don McCoy
    @date    2009-04-05
    @brief   VSIPL++ Library: GPU block class for use with CUDA
*/

#ifndef VSIP_OPT_GPU_BLOCK_HPP
#define VSIP_OPT_GPU_BLOCK_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/
#include <stdexcept>
#include <string>
#include <utility>

#include <vsip/support.hpp>
#include <vsip/domain.hpp>

#include <vsip/core/refcount.hpp>
#include <vsip/core/layout.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/metaprogramming.hpp>

#include <vsip/opt/cuda/bindings.hpp>
#include <vsip/opt/cuda/device_storage.hpp>



/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{

namespace cuda
{

enum memory_state_type
{
  cpu_and_device_memory_valid,
  device_memory_invalid,
  host_memory_invalid
};

} // namespace cuda



/// Implementation class for Gpu_block.
///
/// Gpu_block_impl implements a multi-dimiensional block that is
/// adapted by Gpu_block specializations for each dimension.
///
/// Requires:
///   Dim to be block dimension,
///   T to be a value type,
///   LP to be the block layout policy, encapsulating dimension-ordering,
///      packing format, and complex format,
///   Map to be a map type.
template <dimension_type Dim   = 1,
	  typename       T     = VSIP_DEFAULT_VALUE_TYPE,
	  typename       LP    = Layout<Dim,
					typename Row_major<Dim>::type,
					Stride_unit_dense,
					Cmplx_inter_fmt>,
	  typename       Map   = Local_map>
class Gpu_block_impl
  : public impl::Ref_count<Gpu_block_impl<Dim, T, LP, Map> >
{
  // Compile-time values and types.
public:
  static dimension_type const dim = Dim;

  typedef T        value_type;
  typedef T&       reference_type;
  typedef T const& const_reference_type;

  typedef typename LP::order_type order_type;
  typedef Map                     map_type;

  // Enable Direct_data access to data.
  template <typename, typename, typename>
  friend class impl::data_access::Low_level_data_access;

  // Implementation types.
public:
  typedef LP                                    layout_type;
  typedef impl::Applied_layout<layout_type>     applied_layout_type;
  typedef Allocated_storage<typename LP::complex_type, T> storage_type;
  typedef cuda::Device_storage<T>               device_storage_type;
  typedef cuda::memory_state_type               memory_state_type;

  // Constructors and destructor.
public:
  Gpu_block_impl(Domain<Dim> const& dom, Map const& = Map());
  Gpu_block_impl(Domain<Dim> const& dom, T value, Map const& = Map());
  ~Gpu_block_impl()
    { storage_.deallocate(layout_.total_size()); }

  // Data accessors.
public:

  /// 1-dimensional accessors
  T get(index_type idx) const
  {
    assert(idx < size());
    return storage_.get(idx);
  }

  void put(index_type idx, T val)
  {
    assert(idx < size());
    storage_.put(idx, val);
    invalidate_device_memory();
  }

protected:
  /// Dim-dimensional accessors
  T get(Index<Dim> idx) const
  {
    for (dimension_type d=0; d<Dim; ++d)
      assert(idx[d] < layout_.size(d));
    return storage_.get(layout_.index(idx));
  }

  void put(Index<Dim> idx, T val)
  {
    for (dimension_type d=0; d<Dim; ++d)
      assert(idx[d] < layout_.size(d));
    storage_.put(layout_.index(idx), val);
    invalidate_device_memory();
  }


  // Accessors.
public:
  length_type size() const;
  length_type size(dimension_type D, dimension_type d) const;
  Map const& map() const { return map_;}

  // Support Direct_data interface.
public:
  typedef typename storage_type::type       data_type;
  typedef typename storage_type::const_type const_data_type;

  data_type       impl_data()       { return storage_.data(); }
  const_data_type impl_data() const { return storage_.data(); }

  stride_type impl_stride(dimension_type D, dimension_type d) const;

  // GPU memory interface
protected:
  T*       impl_device_data()       { return device_storage_.data(); }
  T const* impl_device_data() const { return device_storage_.data(); }

  void invalidate_device_memory()       { memory_state_ = cuda::device_memory_invalid; }
  void invalidate_host_memory()         { memory_state_ = cuda::host_memory_invalid; }
  void sync_device_and_host_memory()    { memory_state_ = cuda::cpu_and_device_memory_valid; }
  bool device_memory_is_invalid() const { return memory_state_ == cuda::device_memory_invalid; }
  bool host_memory_is_invalid() const   { return memory_state_ == cuda::host_memory_invalid; }


  // Hidden copy constructor and assignment.
private:
  Gpu_block_impl(Gpu_block_impl const&);
  Gpu_block_impl& operator=(Gpu_block_impl const&);

  // Member Data
private:
  applied_layout_type layout_;
  storage_type        storage_;
  device_storage_type device_storage_;
  map_type            map_;
  memory_state_type   memory_state_;
};



/// General template declaraiton for Gpu_block block class.
///
/// Requires:
///   Dim to be block dimension,
///   T to be a value type,
///   LP to be the block layout policy, encapsulating dimension-ordering,
///      packing format, and complex format,
///   Map to be a map type.
template <dimension_type Dim   = 1,
	  typename       T     = VSIP_DEFAULT_VALUE_TYPE,
	  typename       LP    = Layout<Dim,
					typename Row_major<Dim>::type,
					Stride_unit_dense,
					Cmplx_inter_fmt> >
class Gpu_block;



/// Gpu_block specialization for 1-dimension block.
template <typename       T,
	  typename       LP>
class Gpu_block<1, T, LP>
  : public impl::Gpu_block_impl<1, T, LP, Local_map>
{
  typedef impl::Gpu_block_impl<1, T, LP, Local_map> base_t;
  typedef Gpu_block<1, T, LP> block_type;
  typedef typename block_type::order_type order_type;

  // Constructors.
public:
  Gpu_block(Domain<1> const& dom, Local_map const& map = Local_map())
    : base_t(dom, map)
  {
    // The base class does not initialize the host memory.
    // GPU memory is marked as invalid.
  }

  Gpu_block(Domain<1> const& dom, T value, Local_map const& map = Local_map())
    : base_t(dom, value, map)
  {
    // The base class initializes the host memory and marks
    // the GPU memory as invalid.
  }

  /// Updates host memory from GPU memory (if necessary).
  void device_flush() 
  { 
    if (base_t::device_memory_is_invalid())
    {
      cuda::copy_host_to_dev(*this, base_t::impl_data(), base_t::impl_device_data());
    }
  }

  /// Returns a pointer to GPU memory, but only after updating it
  /// if necessary.
  T* device_data() 
  { 
    if (base_t::device_memory_is_invalid())
    {
      cuda::copy_host_to_dev(*this, base_t::impl_data(), base_t::impl_device_data());
    }
    // Since a non-const pointer is returned, the GPU data must be
    // presumed altered, even if it is not.  Set the flag that 
    // indicates data must be copied back.
    base_t::invalidate_host_memory();
    return base_t::impl_device_data(); 
  }

  T const* device_data() const 
  { 
    if (base_t::device_memory_is_invalid())
    {
      cuda::copy_host_to_dev(*this, base_t::impl_data(), const_cast<T*>(base_t::impl_device_data()));
      const_cast<block_type*>(this)->base_t::sync_device_and_host_memory();
    }
    return base_t::impl_device_data(); 
  }

  /// 1-dimensional data accessors
  T get(Index<1> idx) const
  {
    // Verify data is on the CPU is valid before returning values
    if (base_t::host_memory_is_invalid())
    {
      cuda::copy_dev_to_host(*this, base_t::impl_device_data(), const_cast<T*>(base_t::impl_data()));
      const_cast<block_type*>(this)->base_t::sync_device_and_host_memory();
    }

    return base_t::get(idx);
  }

  // put is inherited
};



/// Gpu_block specialization for 2-dimension block.
template <typename       T,
	  typename       LP>
class Gpu_block<2, T, LP>
  : public impl::Gpu_block_impl<2, T, LP, Local_map>
{
  typedef impl::Gpu_block_impl<2, T, LP, Local_map> base_t;
  typedef Gpu_block<2, T, LP> block_type;
  typedef typename block_type::order_type order_type;

  // Constructors.
public:
  Gpu_block(Domain<2> const& dom, Local_map const& map = Local_map())
    : base_t(dom, map)
  {}

  Gpu_block(Domain<2> const& dom, T value, Local_map const& map = Local_map())
    : base_t(dom, value, map)
  {}

  /// Updates host memory from GPU memory (if necessary).
  void device_flush() 
  { 
    if (base_t::device_memory_is_invalid())
    {
      //      cuda::copy_block<2, order_type, block_type>().
      //        host_to_dev(*this, base_t::impl_data(), base_t::impl_device_data());
      cuda::copy_host_to_dev(*this, base_t::impl_data(), base_t::impl_device_data());
    }
  }

  /// Returns a pointer to GPU memory, but only after updating it
  /// if necessary.
  T* device_data() 
  { 
    if (base_t::device_memory_is_invalid())
    {
      //      cuda::copy_block<1, order_type, block_type>().
      //        host_to_dev(*this, base_t::impl_data(), const_cast<T*>(base_t::impl_device_data()));
      cuda::copy_host_to_dev(*this, base_t::impl_data(), const_cast<T*>(base_t::impl_device_data()));
      base_t::sync_device_and_host_memory();
    }
    // Since a non-const pointer is returned, the GPU data must be
    // presumed altered, even if it is not.  Set the flag that 
    // indicates data must be copied back.
    base_t::invalidate_host_memory();
    return base_t::impl_device_data(); 
  }

  T const* device_data() const 
  { 
    if (base_t::device_memory_is_invalid())
    {
      //      cuda::copy_block<2, order_type, block_type>().
      //        host_to_dev(*this, base_t::impl_data(), const_cast<T*>(base_t::impl_device_data()));
      cuda::copy_host_to_dev(*this, base_t::impl_data(), const_cast<T*>(base_t::impl_device_data()));
      const_cast<block_type*>(this)->base_t::sync_device_and_host_memory();
    }
    return base_t::impl_device_data(); 
  }

  /// 2-dimensional data accessors.
  T get(Index<2> idx) const
  {
    // Verify data is on the CPU is valid before returning values
    if (base_t::host_memory_is_invalid())
    {
      //      cuda::copy_block<2, order_type, block_type>().
      //        dev_to_host(*this, base_t::impl_device_data(), const_cast<T*>(base_t::impl_data()));
      cuda::copy_dev_to_host(*this, base_t::impl_device_data(), const_cast<T*>(base_t::impl_data()));
      base_t::sync_device_and_host_memory();
    }
    return base_t::get(idx);
  }

  T get(index_type idx0, index_type idx1) const
  { 
    // Verify data is on the CPU is valid before returning values
    if (base_t::host_memory_is_invalid())
    {
      //      cuda::copy_block<2, order_type, block_type>().
      //        dev_to_host(*this, base_t::impl_device_data(), const_cast<T*>(base_t::impl_data()));
      cuda::copy_dev_to_host(*this, base_t::impl_device_data(), const_cast<T*>(base_t::impl_data()));
      const_cast<block_type*>(this)->base_t::sync_device_and_host_memory();
    }
    return base_t::get(Index<2>(idx0, idx1)); 
  }

  void put(index_type idx0, index_type idx1, T val)
  { 
    base_t::put(Index<2>(idx0, idx1), val); 
  }
};


/// Specialize block layout trait for Gpu_blocks.
template <dimension_type Dim,
	  typename       T,
	  typename       LP>
struct Block_layout<Gpu_block<Dim, T, LP> >
{
  static dimension_type const dim = Dim;

  typedef Direct_access_tag         access_type;
  typedef typename LP::order_type   order_type;
  typedef typename LP::pack_type    pack_type;
  typedef typename LP::complex_type complex_type;

  typedef LP layout_type;
};



template <dimension_type Dim,
	  typename       T,
	  typename       LP>
struct Is_modifiable_block<Gpu_block<Dim, T, LP> >
{
  static bool const value = true;
};



/***********************************************************************
  Definitions
***********************************************************************/

/// Construct a Gpu_block_impl.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
Gpu_block_impl<Dim, T, LP, Map>::Gpu_block_impl(
  Domain<Dim> const& dom,
  Map const&         map)
  : layout_ (dom),
    storage_(layout_.total_size()),
    device_storage_(layout_.total_size()),
    map_    (map),
    memory_state_(cuda::device_memory_invalid)
{
  // These checks ensure that only supported template paramters are 
  // passed in from the base class constructor.  As support for these
  // are added, these checks may be removed.
  typedef typename LP::pack_type    pack_type;
  typedef typename LP::complex_type complex_type;
  typedef Local_map                 map_type;
  VSIP_IMPL_STATIC_ASSERT((Type_equal<pack_type, Stride_unit_dense>::value));
  VSIP_IMPL_STATIC_ASSERT((Type_equal<complex_type, Cmplx_inter_fmt>::value));
  VSIP_IMPL_STATIC_ASSERT((Type_equal<map_type, Local_map>::value));
}



/// Construct a Gpu_block_impl with initialized data.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
Gpu_block_impl<Dim, T, LP, Map>::Gpu_block_impl(
  Domain<Dim> const& dom,
  T                  val,
  Map const&         map)
  : layout_ (dom),
    storage_(layout_.total_size(), val),
    device_storage_(layout_.total_size()),
    map_    (map),
    memory_state_(cuda::device_memory_invalid)
{
  // These checks ensure that only supported template paramters are 
  // passed in from the base class constructor.  As support for these
  // are added, these checks may be removed.
  typedef typename LP::pack_type    pack_type;
  typedef typename LP::complex_type complex_type;
  typedef Local_map                 map_type;
  VSIP_IMPL_STATIC_ASSERT((Type_equal<pack_type, Stride_unit_dense>::value));
  VSIP_IMPL_STATIC_ASSERT((Type_equal<complex_type, Cmplx_inter_fmt>::value));
  VSIP_IMPL_STATIC_ASSERT((Type_equal<map_type, Local_map>::value));
}



/// Return the total size of block.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
length_type
Gpu_block_impl<Dim, T, LP, Map>::size() const
{
  length_type retval = layout_.size(0);
  for (dimension_type d=1; d<Dim; ++d)
    retval *= layout_.size(d);
  return retval;
}



/// Return the size of the block in a specific dimension.
///
/// Requires:
///   block_dim selects which block-dimensionality (block_dim <= 2).
///   d is the dimension of interest (0 <= d < block_dim).
/// Returns:
///   The size of dimension d.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
length_type
Gpu_block_impl<Dim, T, LP, Map>::size(
  dimension_type block_dim,
  dimension_type d) const
{
  assert((block_dim == 1 || block_dim == Dim) && (d < block_dim));

  if (block_dim == 1)
    return size();
  else
    return layout_.size(d);
}



/// Return the stride of the block in a specific dimension.
/// Requires:
///   block_dim selects which block-dimensionality (block_dim == 1 or 2).
///   d is the dimension of interest (0 <= d < block_dim).
/// Returns:
///   The stride of dimension d.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
stride_type
Gpu_block_impl<Dim, T, LP, Map>::impl_stride(
  dimension_type block_dim,
  dimension_type d) const
{
  assert(block_dim == 1 || block_dim == Dim);
  assert(d < Dim);

  if (block_dim == 1)
    return 1;
  else
    return layout_.stride(d);
}



} // namespace impl
} // namespace vsip

#endif // VSIP_OPT_GPU_BLOCK_HPP
