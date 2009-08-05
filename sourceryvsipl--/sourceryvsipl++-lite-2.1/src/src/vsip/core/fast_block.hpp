/* Copyright (c) 2005 by CodeSourcery, LLC.  All rights reserved. */

/** @file    vsip/core/fast-block.hpp
    @author  Jules Bergmann
    @date    2005-04-12
    @brief   VSIPL++ Library: Fast block class.

*/

#ifndef VSIP_CORE_FAST_BLOCK_HPP
#define VSIP_CORE_FAST_BLOCK_HPP

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


using vsip::Index;

/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

namespace impl
{



/// Implementation class for Fast_block.
///
/// Fast_block_impl implements a multi-dimiensional block that is
/// adapted by Fast_block specializations for each dimension.
///
/// Template parameters:
///   :Dim: to be block dimension,
///   :T: to be a value type,
///   :LP: to be the block layout policy, encapsulating dimension-ordering,
///        packing format, and complex format,
///   :Map: to be a map type.
template <dimension_type Dim   = 1,
	  typename       T     = VSIP_DEFAULT_VALUE_TYPE,
	  typename       LP    = Layout<Dim,
					typename Row_major<Dim>::type,
					Stride_unit_dense,
					Cmplx_inter_fmt>,
	  typename       Map   = Local_map>
class Fast_block_impl
  : public impl::Ref_count<Fast_block_impl<Dim, T, LP, Map> >
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
  typedef impl::Applied_layout<layout_type>      applied_layout_type;
  typedef Allocated_storage<typename LP::complex_type, T> storage_type;

  // Constructors and destructor.
public:
  Fast_block_impl(Domain<Dim> const& dom, Map const& = Map())
    VSIP_THROW((std::bad_alloc));

  Fast_block_impl(Domain<Dim> const& dom, T value, Map const& = Map())
    VSIP_THROW((std::bad_alloc));

  ~Fast_block_impl() VSIP_NOTHROW
    { storage_.deallocate(layout_.total_size()); }

  // Data accessors.
public:

  // 1-dimensional accessors
  T get(index_type idx) const VSIP_NOTHROW
  {
    assert(idx < size());
    return storage_.get(idx);
  }

  void put(index_type idx, T val) VSIP_NOTHROW
  {
    assert(idx < size());
    storage_.put(idx, val);
  }

  // reference_type       operator()(index_type idx) VSIP_NOTHROW;
  // const_reference_type operator()(index_type idx) const VSIP_NOTHROW;

  // Dim-dimensional accessors
protected:
  T    get(Index<Dim> idx) const VSIP_NOTHROW
  {
    for (dimension_type d=0; d<Dim; ++d)
      assert(idx[d] < layout_.size(d));
    return storage_.get(layout_.index(idx));
  }

  void put(Index<Dim> idx, T val) VSIP_NOTHROW
  {
    for (dimension_type d=0; d<Dim; ++d)
      assert(idx[d] < layout_.size(d));
    storage_.put(layout_.index(idx), val);
  }

  // reference_type       operator()(Point<Dim> idx) VSIP_NOTHROW;
  // const_reference_type operator()(Point<Dim> idx) const VSIP_NOTHROW;

  // Accessors.
public:
  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type D, dimension_type d) const VSIP_NOTHROW;
  Map const& map() const VSIP_NOTHROW { return map_;}

  // Support Direct_data interface.
public:
  typedef typename storage_type::type       data_type;
  typedef typename storage_type::const_type const_data_type;

  data_type       impl_data()       VSIP_NOTHROW { return storage_.data(); }
  const_data_type impl_data() const VSIP_NOTHROW { return storage_.data(); }

  stride_type impl_stride(dimension_type D, dimension_type d)
    const VSIP_NOTHROW;

  // Hidden copy constructor and assignment.
private:
  Fast_block_impl(Fast_block_impl const&);
  Fast_block_impl& operator=(Fast_block_impl const&);

  // Member Data
private:
  applied_layout_type layout_;
  storage_type        storage_;
  map_type            map_;
};



/// General template declaraiton for Fast_block block class.
///
/// Template parameters:
///   :Dim: to be block dimension,
///   :T: to be a value type,
///   :LP: to be the block layout policy, encapsulating dimension-ordering,
///        packing format, and complex format,
///   :Map: to be a map type.
template <dimension_type Dim   = 1,
	  typename       T     = VSIP_DEFAULT_VALUE_TYPE,
	  typename       LP    = Layout<Dim,
					typename Row_major<Dim>::type,
					Stride_unit_dense,
					Cmplx_inter_fmt>,
	  typename       Map   = Local_map>
class Fast_block;



/// Fast_block specialization for 1-dimension block.
template <typename       T,
	  typename       LP,
	  typename       Map>
class Fast_block<1, T, LP, Map>
  : public impl::Fast_block_impl<1, T, LP, Map>
{
  typedef impl::Fast_block_impl<1, T, LP, Map> bast_t;

  // Constructors.
public:
  Fast_block(Domain<1> const& dom, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, map)
  {}

  Fast_block(Domain<1> const& dom, T value, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, value, map)
  {}

  // Use inherited 1-dimensional get/put.
};



/// Fast_block specialization for 2-dimension block.
template <typename       T,
	  typename       LP,
	  typename       Map>
class Fast_block<2, T, LP, Map>
  : public impl::Fast_block_impl<2, T, LP, Map>
{
  typedef impl::Fast_block_impl<2, T, LP, Map> bast_t;

  // Constructors.
public:
  Fast_block(Domain<2> const& dom, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, map)
  {}

  Fast_block(Domain<2> const& dom, T value, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, value, map)
  {}

  // Use inherited 1-dimensional get/put.

  // 2-dimensional data accessors.
public:
  T get(index_type idx0, index_type idx1) const VSIP_NOTHROW
    { return bast_t::get(Index<2>(idx0, idx1)); }
  void put(index_type idx0, index_type idx1, T val) VSIP_NOTHROW
    { bast_t::put(Index<2>(idx0, idx1), val); }
};



/// Fast_block specialization for 3-dimension block.
template <typename T,
	  typename LP,
	  typename Map>
class Fast_block<3, T, LP, Map>
  : public impl::Fast_block_impl<3, T, LP, Map>
{
  typedef impl::Fast_block_impl<3, T, LP, Map> bast_t;

  // Constructors.
public:
  Fast_block(Domain<3> const& dom, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, map)
  {}

  Fast_block(Domain<3> const& dom, T value, Map const& map = Map())
    VSIP_THROW((std::bad_alloc))
  : bast_t(dom, value, map)
  {}

  // Use inherited 1-dimensional get/put.

  // 3-dimensional data accessors.
public:
  T get(index_type idx0, index_type idx1, index_type idx2) const VSIP_NOTHROW
    { return bast_t::get(Index<3>(idx0, idx1, idx2)); }
  void put(index_type idx0, index_type idx1, index_type idx2, T val)
    VSIP_NOTHROW
    { bast_t::put(Index<3>(idx0, idx1, idx2), val); }
};



/// Specialize block layout trait for Fast_blocks.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
struct Block_layout<Fast_block<Dim, T, LP, Map> >
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
	  typename       LP,
	  typename       Map>
struct Block_layout<Fast_block<Dim, T, LP, Map> const>
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
	  typename       LP,
	  typename       Map>
struct Is_modifiable_block<Fast_block<Dim, T, LP, Map> >
{
  static bool const value = true;
};



/***********************************************************************
  Definitions
***********************************************************************/

/// Construct a Fast_block_impl.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
Fast_block_impl<Dim, T, LP, Map>::Fast_block_impl(
  Domain<Dim> const& dom,
  Map const&         map)
  VSIP_THROW((std::bad_alloc))
  : layout_ (dom),
    storage_(layout_.total_size()),
    map_    (map)
{
}



/// Construct a Fast_block_impl with initialized data.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
Fast_block_impl<Dim, T, LP, Map>::Fast_block_impl(
  Domain<Dim> const& dom,
  T                  val,
  Map const&         map)
  VSIP_THROW((std::bad_alloc))
  : layout_ (dom),
    storage_(layout_.total_size(), val),
    map_    (map)
{
}



/// Return the total size of block.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
length_type
Fast_block_impl<Dim, T, LP, Map>::size() const VSIP_NOTHROW
{
  length_type retval = layout_.size(0);
  for (dimension_type d=1; d<Dim; ++d)
    retval *= layout_.size(d);
  return retval;
}



/// Return the size of the block in a specific dimension.
///
/// Parameters:
///   :block_dim: is the dimension whose length to return (0 <= DIM < BLOCK_DIM).
///   :d: selects which block-dimensionality (BLOCK_DIM <= 2).
///
/// Returns:
///   The size of dimension `d`.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
length_type
Fast_block_impl<Dim, T, LP, Map>::size(
  dimension_type block_dim,
  dimension_type d)
  const VSIP_NOTHROW
{
  assert((block_dim == 1 || block_dim == Dim) && (d < block_dim));

  if (block_dim == 1)
    return size();
  else
    return layout_.size(d);
}



/// Parameters:
///   :block_dim: is a valid dimensionality supported by block (DIM == 1 or 2)
///   :d: is a dimension, less than `block_dim`.
///
/// Returns:
///   The stride in dimension D, for dimensionality DIM.
template <dimension_type Dim,
	  typename       T,
	  typename       LP,
	  typename       Map>
inline
stride_type
Fast_block_impl<Dim, T, LP, Map>::impl_stride(
  dimension_type block_dim,
  dimension_type d)
  const VSIP_NOTHROW
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

#endif // VSIP_CORE_FAST_BLOCK_HPP
