/* Copyright (c) 2005, 2006, 2008 by CodeSourcery.  All rights reserved. */

/** @file    vsip/dense.hpp
    @author  Jules Bergmann
    @date    2005-01-20
    @brief   VSIPL++ Library: Dense block class.

*/

#ifndef VSIP_DENSE_HPP
#define VSIP_DENSE_HPP

/***********************************************************************
  Included Files
***********************************************************************/
#include <stdexcept>
#include <string>

#include <vsip/support.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/dense_fwd.hpp>
#include <vsip/core/refcount.hpp>
#include <vsip/core/parallel/local_map.hpp>
#include <vsip/core/layout.hpp>
#include <vsip/core/extdata.hpp>
#include <vsip/core/block_traits.hpp>
#include <vsip/core/parallel/choose_dist_block.hpp>
#include <vsip/domain.hpp>
#include <vsip/core/memory_pool.hpp>

/// Complex storage format for dense blocks.
#if VSIP_IMPL_PREFER_SPLIT_COMPLEX
#  define VSIP_IMPL_DENSE_CMPLX_FMT vsip::impl::Cmplx_split_fmt
#else
#  define VSIP_IMPL_DENSE_CMPLX_FMT vsip::impl::Cmplx_inter_fmt
#endif


/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{

enum user_storage_type 
{
  no_user_format = 0,
  array_format,
  interleaved_format,
  split_format
};



namespace impl
{ 

typedef VSIP_IMPL_DENSE_CMPLX_FMT dense_complex_type;



/// If provided type is complex, extract the component type,
/// otherwise use the provided bogus type.
template <typename T,
	  typename BogusT>
struct Complex_value_type
{
  typedef BogusT type; 
  typedef BogusT ptr_type;
};

template <typename T,
	  typename BogusT>
struct Complex_value_type<complex<T>, BogusT>
{
  typedef T  type; 
  typedef T* ptr_type;
};



/// Class to hold storage that has been given to us by the user.
/// For complex data, this can be in three different formats (array
/// of complex, array of interleaved real, or arrays of split real).
/// The block value type does not determine/restrict which type of
/// data the block can be constructed with.
///
/// General case intended for non-complex types.
template <typename T>
class User_storage
{
  // Constructors.
public:
  User_storage() : format_(no_user_format), data_(0) {}

  User_storage(user_storage_type format, T* data)
    : format_(format), data_(data)
  { assert(format == array_format); }

  /// This constructor is provided so that `User_storage<T>` and
  /// `User_storage<complex<T> >` can be interchanged, however it
  /// should not be called.
  User_storage(user_storage_type format, T* real, T* /*imag*/)
    : format_(format), data_(real)
  { assert(0); }

  // Accessors.
public:
  user_storage_type format() const { return format_; }

  T    get(index_type i) { return this->data_[i]; }
  void put(index_type i, T val) { this->data_[i] = val; }

  // Return the user storage in a format acceptable for initializing
  // a Storage class.
  typename Storage<Cmplx_inter_fmt, T>::type as_storage(Cmplx_inter_fmt) const
  {
    assert(this->format_ == array_format);
    return this->data_;
  }

  // For scalar types, there is no distinction between interleaved and
  // split.
  typename Storage<Cmplx_split_fmt, T>::type as_storage(Cmplx_split_fmt) const
  {
    assert(this->format_ == array_format);
    return this->data_;
  }

  void find(T*& pointer)
  { pointer = this->data_; }

  void rebind(T* pointer)
  { this->data_ = pointer; }

  // Default copy-constructor is OK.
  // Default assignment is OK.

  // Member data.
private:
  user_storage_type format_;
  T* data_;
};



/// User_storage specialization for complex types.
///
/// Can store user-storage in array format, interleaved format, or
/// split format.
template <typename T>
class User_storage<complex<T> >
{
  // Constructors.
public:
  User_storage()
    : format_(no_user_format)
  {
    // Zero everything out so that find() returns NULL values.
    this->u_.data_ = 0;
    this->u_.split_.real_ = 0;
    this->u_.split_.imag_ = 0;
  }

  User_storage(user_storage_type format, complex<T>* data)
    : format_(format)
  {
    assert(format == array_format);
    this->u_.data_ = data;
  }

  User_storage(user_storage_type format, T* real, T* imag)
    : format_(format)
  { 
    assert(format == interleaved_format || format == split_format);
    this->u_.split_.real_ = real;
    this->u_.split_.imag_ = imag;
  }

  // Accessors.
public:
  user_storage_type format() const { return this->format_; }

  complex<T> get(index_type i)
  {
    assert(this->format_ != no_user_format);

    if (this->format_ == array_format)
      return this->u_.data_[i];
    else if (this->format_ == interleaved_format)
      return complex<T>(this->u_.split_.real_[2*i+0],
			this->u_.split_.real_[2*i+1]);
    else // if (format_ == split_format)
      return complex<T>(this->u_.split_.real_[i],
			this->u_.split_.imag_[i]);
  }

  void put(index_type i, complex<T> val)
  {
    assert(this->format_ != no_user_format);

    if (this->format_ == array_format)
      this->u_.data_[i] = val;
    else if (this->format_ == interleaved_format)
    {
      this->u_.split_.real_[2*i+0] = val.real();
      this->u_.split_.real_[2*i+1] = val.imag();
    }
    else // if (format_ == split_format)
    {
      this->u_.split_.real_[i] = val.real();
      this->u_.split_.imag_[i] = val.imag();
    }
  }

  /// Return the user storage in a format acceptable for initializing
  /// a Storage class with array/interleaved format, or return NULL.
  typename Storage<Cmplx_inter_fmt, complex<T> >::type
  as_storage(Cmplx_inter_fmt) const
  {
    assert(this->format_ != no_user_format);

    if (this->format_ == array_format)
      return this->u_.data_;
    else if (this->format_ == interleaved_format)
      return (complex<T>*)this->u_.split_.real_;
    else // if (format_ == split_format)
      return NULL;
  }

  /// Return the user storage in a format acceptable for initializing
  /// a Storage class with split format, or return effective NULL.
  typename Storage<Cmplx_split_fmt, complex<T> >::type
  as_storage(Cmplx_split_fmt) const
  {
    assert(this->format_ != no_user_format);

    if (this->format_ == array_format)
      return std::pair<T*, T*>(0, 0);
    else if (this->format_ == interleaved_format)
      return std::pair<T*, T*>(0, 0);
    else // if (format_ == split_format)
      return std::pair<T*, T*>(this->u_.split_.real_, this->u_.split_.imag_);
  }

  void find(complex<T>*& pointer)
  { pointer = this->u_.data_; }

  void find(T*& pointer)
  { pointer = this->u_.split_.real_; }

  void find(T*& real_pointer, T*& imag_pointer)
  {
    real_pointer = this->u_.split_.real_;
    imag_pointer = this->u_.split_.imag_;
  }

  void rebind(complex<T>* pointer)
  {
    this->u_.data_ = pointer;
    this->format_  = array_format;
  }

  void rebind(T* pointer)
  {
    this->u_.split_.real_ = pointer; 
    this->format_         = interleaved_format;
  }

  void rebind(T* real_pointer, T* imag_pointer)
  {
    this->u_.split_.real_ = real_pointer;
    this->u_.split_.imag_ = imag_pointer;
    this->format_         = split_format;
  }

  // Default copy-constructor is OK.
  // Default assignment is OK.

  // Member data.
private:
  user_storage_type format_;
  union {
    complex<T>*       data_;
    struct {
      T* real_;
      T* imag_;
    } split_;
  } u_;

};


template <typename ComplexFmt,
	  typename T>
class Dense_storage
{
  // Compile-time values and types.
public:
  typedef T*       type;
  typedef T const* const_type;

  // Constructors and destructor.
public:
  Dense_storage(Memory_pool*  pool,
		length_type   size,
		type          buffer = NULL)
    VSIP_THROW((std::bad_alloc))
    : alloc_data_(buffer == NULL),
      data_      (alloc_data_ ? pool->allocate<T>(size) : (T*)buffer)
  {}

  Dense_storage(Memory_pool*  pool,
		length_type   size,
		T             val,
		type          buffer = NULL)
  VSIP_THROW((std::bad_alloc))
    : alloc_data_(buffer == NULL),
      data_      (alloc_data_ ? pool->allocate<T>(size) : (T*)buffer)
  {
    for (index_type i=0; i<size; ++i)
      data_[i] = val;
  }

  ~Dense_storage()
  {
    // user's responsiblity to call deallocate().
    if (alloc_data_)
      assert(data_ == 0);
  }

  // Accessors.
protected:
  void impl_rebind(Memory_pool* pool, length_type size, type buffer);

  void deallocate(Memory_pool* pool, length_type size)
  {
    if (alloc_data_)
    {
      pool->deallocate(data_, size);
      data_ = 0;
    }
  }

  bool is_alloc() const { return alloc_data_; }

  // Accessors.
public:
  T    get(index_type idx) const { return data_[idx]; }
  void put(index_type idx, T val){ data_[idx] = val; }

  T&       impl_ref(index_type idx)       { return data_[idx]; }
  T const& impl_ref(index_type idx) const { return data_[idx]; }

  type       impl_data()       { return data_; }
  const_type impl_data() const { return data_; }

  // Member data.
private:
  bool   alloc_data_;
  T*     data_;
};



template <typename T>
class Dense_storage<Cmplx_split_fmt, vsip::complex<T> >
{
  // Compile-time values and types.
public:
  typedef std::pair<T*, T*>             type;
  typedef std::pair<T const*, T const*> const_type;

  // Constructors and destructor.
public:
  Dense_storage(Memory_pool*  pool,
		length_type   size,
		type          buffer    = type(0, 0))
    VSIP_THROW((std::bad_alloc))
    : alloc_data_(buffer.first == NULL || buffer.second == NULL),
      real_data_ (alloc_data_ ? pool->allocate<T>(size) : buffer.first),
      imag_data_ (alloc_data_ ? pool->allocate<T>(size) : buffer.second)
  {}

  Dense_storage(Memory_pool*     pool,
		length_type      size,
		vsip::complex<T> val,
		type buffer = type(0, 0))
    VSIP_THROW((std::bad_alloc))
    : alloc_data_(buffer.first == NULL || buffer.second == NULL),
      real_data_ (alloc_data_ ? pool->allocate<T>(size) : buffer.first),
      imag_data_ (alloc_data_ ? pool->allocate<T>(size) : buffer.second)
  {
    for (index_type i=0; i<size; ++i)
      real_data_[i] = val.real();
    for (index_type i=0; i<size; ++i)
      imag_data_[i] = val.imag();
  }

  ~Dense_storage()
  {
    // user's responsiblity to call deallocate().
    if (alloc_data_)
      assert(real_data_ == 0 && imag_data_ == 0);
  }

  // Accessors.
protected:
  void impl_rebind(Memory_pool* pool, length_type size, type buffer);

  void deallocate(Memory_pool* pool, length_type size)
  {
    if (alloc_data_)
    {
      pool->deallocate(real_data_, size);
      pool->deallocate(imag_data_, size);
      real_data_ = 0;
      imag_data_ = 0;
    }
  }

  bool is_alloc() const { return alloc_data_; }

  // Accessors.
public:
  vsip::complex<T> get(index_type idx) const
    { return vsip::complex<T>(real_data_[idx], imag_data_[idx]); }

  void put(index_type idx, vsip::complex<T> val)
  {
    real_data_[idx] = val.real();
    imag_data_[idx] = val.imag();
  }

  vsip::complex<T>&       impl_ref(index_type)
    { VSIP_IMPL_THROW(unimplemented(
	"Dense_storage<Cmplx_split_fmt>::impl_ref - unimplemented")); }
  vsip::complex<T> const& impl_ref(index_type) const
    { VSIP_IMPL_THROW(unimplemented(
	"Dense_storage<Cmplx_split_fmt>::impl_ref - unimplemented")); }

  type       impl_data()       { return type(real_data_, imag_data_); }
  const_type impl_data() const { return const_type(real_data_, imag_data_); }

  // Member data.
private:
  bool   alloc_data_;
  T*     real_data_;
  T*     imag_data_;
};



/// Dense implementation.
///
/// Dense blocks will derive from Dense_impl and provide:
///  - Constructors
///  - (Dim > 1)-dimension get/put
///
/// Note:
///   This class has several user-storage functions that are only
///   valid when the value_type is complex.  When the block value_type
///   is complex, these functions accept pointers to the complex
///   component type (i.e. the T type of a complex<T>).  When the
///   block value_type is not complex, these functions accept 
///   pointers to a private type hidden within Dense_impl, preventing
///   their use.
template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT   = Local_map>
class Dense_impl
  : public impl::Dense_storage<ComplexFmt, T>,
    public impl::Ref_count<Dense_impl<Dim, T, OrderT, ComplexFmt, MapT> >
{
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Compile-time values and types.
public:
  static dimension_type const dim = Dim;

  typedef T        value_type;
  typedef T&       reference_type;
  typedef T const& const_reference_type;

  typedef OrderT order_type;
  typedef MapT   map_type;

  // Enable Direct_data access to data.
  template <typename, typename, typename>
  friend class impl::data_access::Low_level_data_access;

  // Implementation types.
public:
  typedef ComplexFmt complex_type;
  typedef impl::Layout<Dim, order_type, impl::Stride_unit_dense,
		       complex_type> layout_type;
  typedef impl::Applied_layout<layout_type>   applied_layout_type;
  typedef impl::Dense_storage<complex_type, T> storage_type;

  // Constructors and destructor.
public:
  Dense_impl(Domain<Dim> const& dom, MapT const& = MapT())
    VSIP_THROW((std::bad_alloc));

  Dense_impl(Domain<Dim> const& dom, T value, MapT const& = MapT())
    VSIP_THROW((std::bad_alloc));

  // User storage constructor.
  Dense_impl(Domain<Dim> const& dom,
	     User_storage<T> const&  data,
	     MapT const& = MapT())
    VSIP_THROW((std::bad_alloc));

  ~Dense_impl() VSIP_NOTHROW
    { storage_type::deallocate(map_.impl_pool(), layout_.total_size()); }

public:
  using storage_type::get;
  using storage_type::put;

protected:
  // Dim-dimensional get/put
  T    get(Index<Dim> const& idx) const VSIP_NOTHROW;
  void put(Index<Dim> const& idx, T val) VSIP_NOTHROW;

  // 2-diminsional get/put
  T    impl_get(index_type idx0, index_type idx1) const VSIP_NOTHROW
    { return this->get(layout_.index(idx0, idx1)); }
  void impl_put(index_type idx0, index_type idx1, T val) VSIP_NOTHROW
    { this->put(layout_.index(idx0, idx1), val); }

  // 3-diminsional get/put
  T    impl_get(index_type idx0, index_type idx1, index_type idx2)
    const VSIP_NOTHROW
    { return this->get(layout_.index(idx0, idx1, idx2)); }
  void impl_put(index_type idx0, index_type idx1, index_type idx2, T val)
    VSIP_NOTHROW
    { this->put(layout_.index(idx0, idx1, idx2), val); }

public:
  using storage_type::impl_ref;

protected:
  // Dim-dimensional lvalue.
  reference_type       impl_ref(Index<Dim> const& idx) VSIP_NOTHROW;
  const_reference_type impl_ref(Index<Dim> const& idx) const VSIP_NOTHROW;

  // Accessors.
public:
  length_type size() const VSIP_NOTHROW;
  length_type size(dimension_type D, dimension_type d) const VSIP_NOTHROW;
  MapT const& map() const VSIP_NOTHROW { return map_;}

  // User storage functions.
public:
  void admit(bool update = true) VSIP_NOTHROW;
  void release(bool update = true) VSIP_NOTHROW;
  void release(bool update, T*& pointer) VSIP_NOTHROW;
  void release(bool update, uT*& pointer) VSIP_NOTHROW;
  void release(bool update, uT*& real_pointer, uT*& imag_pointer) VSIP_NOTHROW;

  void find(T*& pointer) VSIP_NOTHROW;
  void find(uT*& pointer) VSIP_NOTHROW;
  void find(uT*& real_pointer, uT*& imag_pointer) VSIP_NOTHROW;

  void rebind(T* pointer) VSIP_NOTHROW;
  void rebind(uT* pointer) VSIP_NOTHROW;
  void rebind(uT* real_pointer, uT* imag_pointer) VSIP_NOTHROW;

public:
  enum user_storage_type user_storage() const VSIP_NOTHROW;
  bool admitted() const VSIP_NOTHROW;

  // Support Direct_data interface.
public:
  typedef typename storage_type::type       data_type;
  typedef typename storage_type::const_type const_data_type;

  // data_type       impl_data()       VSIP_NOTHROW { return storage_.data(); }
  // const_data_type impl_data() const VSIP_NOTHROW { return storage_.data(); }
  stride_type impl_stride(dimension_type D, dimension_type d)
    const VSIP_NOTHROW;

  stride_type impl_offset() VSIP_NOTHROW
  { return 0; }

  // Hidden copy constructor and assignment.
private:
  Dense_impl(Dense_impl const&);
  Dense_impl& operator=(Dense_impl const&);

  // Member Data
private:
  applied_layout_type layout_;
  User_storage<T>     user_data_;
  map_type            map_;
  bool                admitted_;
};

} // namespace vsip::impl



/// Partial specialization of Dense class template for 1-dimension.

/// Dense block, as defined in standard [view.dense].
///
/// A Dense block is a modifiable, allocatable 1-dimensional block
/// or 1,x-dimensional block, for a fixed x, that explicitly stores
/// one value for each Index in its domain.
///
/// :Requirements:
///   :T: value_type
///   :OrderT: order_type
template <typename T,
	  typename OrderT>
class Dense<1, T, OrderT, Local_map>
  : public impl::Dense_impl<1, T, OrderT, impl::dense_complex_type, Local_map>
{
  typedef impl::Dense_impl<1, T, OrderT, impl::dense_complex_type, Local_map>
   	  base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Compile-time values and types.
  // These are defined in Dense_impl, but redefined for convenience.
public:
  typedef typename base_type::map_type map_type;

  // Constructors.
public:
  Dense(Domain<1> const& dom, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<1> const& dom, T value, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<1> const& dom,
	T*const          pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(array_format, pointer), map)
    {}

  Dense(Domain<1> const& dom,
	uT*const         pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(interleaved_format, pointer, 0),
		  map)
    {}

  Dense(Domain<1> const& dom,
	uT*const         real_pointer,
	uT*const         imag_pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(split_format,
					     real_pointer, imag_pointer), map)
    {}

  // Internal user storage constructor.
  Dense(Domain<1> const&             dom,
	impl::User_storage<T> const& data,
	map_type const&              map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, data, map)
    {}
};



/// Partial specialization of Dense class template for distributed 1-dim.

template <typename T,
	  typename OrderT,
	  typename MapT>
class Dense<1, T, OrderT, MapT>
  : public impl::Choose_dist_block<1, T, OrderT, MapT>::type
{
  typedef typename impl::Choose_dist_block<1, T, OrderT, MapT>::type base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Constructors.
public:
  Dense(Domain<1> const& dom, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<1> const& dom, T value, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<1> const& dom, T* const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<1> const& dom, uT*const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<1> const& dom,
	uT*const real_pointer, uT*const imag_pointer,
	MapT const& map = MapT())
      : base_type(dom, real_pointer, imag_pointer, map)
    {}
};



/// Partial specialization of Dense class template for 1,2-dimension.

template <typename T,
	  typename OrderT>
class Dense<2, T, OrderT, Local_map>
  : public impl::Dense_impl<2, T, OrderT, impl::dense_complex_type, Local_map>
{
  typedef impl::Dense_impl<2, T, OrderT, impl::dense_complex_type, Local_map>
          base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Compile-time values and types, redefined for implementation convenience.
public:
  typedef typename base_type::map_type             map_type;
  typedef typename base_type::reference_type       reference_type;
  typedef typename base_type::const_reference_type const_reference_type;

  // Constructors.
public:
  Dense(Domain<2> const& dom, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<2> const& dom, T value, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<2> const& dom,
	T*const          pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(array_format, pointer), map)
    {}

  Dense(Domain<2> const& dom,
	uT*const         pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(interleaved_format, pointer, 0),
		  map)
    {}
  Dense(Domain<2> const& dom,
	uT*const         real_pointer,
	uT*const         imag_pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(split_format,
					     real_pointer, imag_pointer), map)
    {}

  // Internal user storage constructor.
  Dense(Domain<2> const&             dom,
	impl::User_storage<T> const& data,
	map_type const&              map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, data, map)
    {}

  // 1-dim Data Accessors.
  using base_type::get;
  using base_type::put;

  // 2-dim Data Accessors.
public:
  T get(index_type idx0, index_type idx1) const VSIP_NOTHROW
    { return base_type::impl_get(idx0, idx1); }

  void put(index_type idx0, index_type idx1, T val) VSIP_NOTHROW
    { return base_type::impl_put(idx0, idx1, val); }

  reference_type impl_ref(index_type idx0, index_type idx1)
    VSIP_NOTHROW
    { return base_type::impl_ref(Index<2>(idx0, idx1)); }

  const_reference_type impl_ref(index_type idx0, index_type idx1)
    const VSIP_NOTHROW
    { return base_type::impl_ref(Index<2>(idx0, idx1)); }
};



/// Partial specialization of Dense class template for distributed 1,2-dim.

template <typename T,
	  typename OrderT,
	  typename MapT>
class Dense<2, T, OrderT, MapT>
  : public impl::Choose_dist_block<2, T, OrderT, MapT>::type
{
  typedef typename impl::Choose_dist_block<2, T, OrderT, MapT>::type base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Constructors.
public:
  Dense(Domain<2> const& dom, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<2> const& dom, T value, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<2> const& dom, T* const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<2> const& dom, uT*const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<2> const& dom,
	uT*const real_pointer, uT*const imag_pointer,
	MapT const& map = MapT())
      : base_type(dom, real_pointer, imag_pointer, map)
    {}

};



/// Partial specialization of Dense class template for 1,3-dimension.

template <typename T,
	  typename OrderT>
class Dense<3, T, OrderT, Local_map>
  : public impl::Dense_impl<3, T, OrderT, impl::dense_complex_type, Local_map>
{
  typedef impl::Dense_impl<3, T, OrderT, impl::dense_complex_type, Local_map>
          base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Compile-time values and types, redefined for implementation convenience.
public:
  typedef typename base_type::map_type             map_type;
  typedef typename base_type::reference_type       reference_type;
  typedef typename base_type::const_reference_type const_reference_type;

  // Constructors.
public:
  Dense(Domain<3> const& dom, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<3> const& dom, T value, map_type const& map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<3> const& dom,
	T*const          pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(array_format, pointer), map)
    {}

  Dense(Domain<3> const& dom,
	uT*const         pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(interleaved_format, pointer, 0),
		  map)
    {}
  Dense(Domain<3> const& dom,
	uT*const         real_pointer,
	uT*const         imag_pointer,
	map_type const&  map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, impl::User_storage<T>(split_format,
					     real_pointer, imag_pointer), map)
    {}

  // Internal user storage constructor.
  Dense(Domain<3> const&             dom,
	impl::User_storage<T> const& data,
	map_type const&              map = map_type())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, data, map)
    {}

  // 1-dim Data Accessors.
  using base_type::get;
  using base_type::put;

  // 3-dim Data Accessors.
public:
  T get(index_type idx0, index_type idx1, index_type idx2)
    const VSIP_NOTHROW
    { return base_type::impl_get(idx0, idx1, idx2); }

  void put(index_type idx0, index_type idx1, index_type idx2, T val)
    VSIP_NOTHROW
    { base_type::impl_put(idx0, idx1, idx2, val); }

  reference_type impl_ref(index_type idx0, index_type idx1, index_type idx2)
    VSIP_NOTHROW
    { return base_type::impl_ref(Index<3>(idx0, idx1, idx2)); }

  const_reference_type impl_ref(index_type idx0, index_type idx1,
				  index_type idx2)
    const VSIP_NOTHROW
    { return base_type::impl_ref(Index<3>(idx0, idx1, idx2)); }
};



/// Partial specialization of Dense class template for distributed 1,3-dim.

template <typename T,
	  typename OrderT,
	  typename MapT>
class Dense<3, T, OrderT, MapT>
  : public impl::Choose_dist_block<3, T, OrderT, MapT>::type
{
  typedef typename impl::Choose_dist_block<3, T, OrderT, MapT>::type base_type;
  enum private_type {};
  typedef typename impl::Complex_value_type<T, private_type>::type uT;

  // Constructors.
public:
  Dense(Domain<3> const& dom, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, map)
    {}

  Dense(Domain<3> const& dom, T value, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, value, map)
    {}

  Dense(Domain<3> const& dom, T* const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<3> const& dom, uT*const pointer, MapT const& map = MapT())
    VSIP_THROW((std::bad_alloc))
      : base_type(dom, pointer, map)
    {}

  Dense(Domain<3> const& dom,
	uT*const real_pointer, uT*const imag_pointer,
	MapT const& map = MapT())
      : base_type(dom, real_pointer, imag_pointer, map)
    {}
};



namespace impl
{

/// Specialize block layout trait for Dense_impl blocks.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
struct Block_layout<Dense_impl<Dim, T, OrderT, ComplexFmt, MapT> >
{
  static dimension_type const dim = Dim;

  typedef Direct_access_tag access_type;
  typedef OrderT            order_type;
  typedef Stride_unit_dense pack_type;
  typedef ComplexFmt        complex_type;

  typedef Layout<dim, order_type, pack_type, complex_type> layout_type;
};

/// Specialize block layout trait for Dense blocks.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Block_layout<Dense<Dim, T, OrderT, MapT> >
  : Block_layout<Dense_impl<Dim, T, OrderT, dense_complex_type, MapT> >
{};
#if 0
template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Block_layout<Dense<Dim, T, OrderT, MapT> const>
  : Block_layout<Dense_impl<Dim, T, OrderT, dense_complex_type, MapT> >
{};
#endif


/// Specialize lvalue accessor trait for Dense blocks.
/// Dense provides direct lvalue accessors via impl_ref.

template <typename BlockT,
	  bool     use_proxy = Is_split_block<BlockT>::value>
struct Dense_lvalue_factory_type;

template <typename BlockT>
struct Dense_lvalue_factory_type<BlockT, false>
{
  typedef True_lvalue_factory<BlockT> type;
  template <typename OtherBlock>
  struct Rebind {
    typedef True_lvalue_factory<OtherBlock> type;
  };
};

template <typename BlockT>
struct Dense_lvalue_factory_type<BlockT, true>
{
  typedef Proxy_lvalue_factory<BlockT> type;
  template <typename OtherBlock>
  struct Rebind {
    typedef Proxy_lvalue_factory<OtherBlock> type;
  };
};

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
struct Lvalue_factory_type<Dense<Dim, T, OrderT, Local_map> >
  : public Dense_lvalue_factory_type<Dense<Dim, T, OrderT, Local_map> >
{};



/// Specialize Distributed_local_block traits class for Dense.

/// For a serial map, distributed block and local block are the same.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
struct Distributed_local_block<Dense<Dim, T, OrderT, Local_map> >
{
  typedef Dense<Dim, T, OrderT, Local_map> type;
  typedef Dense<Dim, T, OrderT, Local_map> proxy_type;
};



/// For a distributed map, local block has a serial map.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Distributed_local_block<Dense<Dim, T, OrderT, MapT> >
{
  // We could determine the local block by just changing the map
  // to serial:
  //   typedef Dense<Dim, T, OrderT, Local_map> type;

  // However, to be safe, we'll extract it from the block itself:
  // (local_block is set in the base class Distributed_block.)
  typedef typename Dense<Dim, T, OrderT, MapT>::local_block_type type;
  typedef typename Dense<Dim, T, OrderT, MapT>::proxy_local_block_type
    proxy_type;
};



/// Overload of get_local_block for Dense with serial map.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
Dense<Dim, T, OrderT, Local_map>&
get_local_block(
  Dense<Dim, T, OrderT, Local_map>& block)
{
  return block;
}

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
Dense<Dim, T, OrderT, Local_map> const&
get_local_block(
  Dense<Dim, T, OrderT, Local_map> const& block)
{
  return block;
}



/// Overload of get_local_block for Dense with distributed map.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
inline typename Dense<Dim, T, OrderT, MapT>::local_block_type&
get_local_block(
  Dense<Dim, T, OrderT, MapT> const& block)
{
  return block.get_local_block();
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
inline typename Dense<Dim, T, OrderT, MapT>::proxy_local_block_type
get_local_proxy(
  Dense<Dim, T, OrderT, MapT> const& block,
  index_type                         sb)
{
  return block.impl_proxy_block(sb);
}



/// Assert that subblock is local to block (overload).

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
void
assert_local(
  Dense<Dim, T, OrderT, Local_map> const& /*block*/,
  index_type                              sb)
{
  assert(sb == 0);
}



/// Assert that subblock is local to block (overload).

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
void
assert_local(
  Dense<Dim, T, OrderT, MapT> const& block,
  index_type                         sb)
{
  block.assert_local(sb);
}



/// Specialize Is_simple_distributed_block traits class for Dense.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Is_simple_distributed_block<Dense<Dim, T, OrderT, MapT> >
{
  static bool const value = true;
};



#if VSIP_IMPL_USE_GENERIC_VISITOR_TEMPLATES==0

/// Specialize Combine_return_type for Dense block leaves.

template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Combine_return_type<CombineT, Dense<Dim, T, OrderT, MapT> >
{
  typedef Dense<Dim, T, OrderT, MapT> block_type;
  typedef typename CombineT::template return_type<block_type>::type
		type;
  typedef typename CombineT::template tree_type<block_type>::type
		tree_type;
};



/// Specialize apply_combine for Dense block leaves.

template <typename       CombineT,
	  dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
typename Combine_return_type<CombineT, Dense<Dim, T, OrderT, MapT> >::type
apply_combine(
  CombineT const&                    combine,
  Dense<Dim, T, OrderT, MapT> const& block)
{
  return combine.apply(block);
}
#endif



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT>
struct Is_pas_block<Dense<Dim, T, OrderT, Local_map> >
{
  static bool const value = false;
};

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Is_pas_block<Dense<Dim, T, OrderT, MapT> >
  : Is_pas_block<typename impl::Choose_dist_block<Dim, T, OrderT, MapT>::type>
{};



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       MapT>
struct Is_modifiable_block<Dense<Dim, T, OrderT, MapT> >
{
  static bool const value = true;
};



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
struct Is_modifiable_block<Dense_impl<Dim, T, OrderT, ComplexFmt, MapT> >
{
  static bool const value = true;
};



/***********************************************************************
  Definitions - Dense_storage
***********************************************************************/

/// Rebind the memory referred to by Dense_storage object
///
/// Requires:
///   SIZE to be size object was constructed with.

template <typename ComplexFmt,
	  typename T>
void
Dense_storage<ComplexFmt, T>::impl_rebind(
  Memory_pool* pool,
  length_type  size,
  type         buffer)
{
  if (buffer != NULL)
  {
    if (alloc_data_)
      pool->deallocate<T>(data_, size);
    
    alloc_data_ = false;
    data_       = buffer;
  }
  else // (buffer == NULL
  {
    if (!alloc_data_)
    {
      alloc_data_ = true;
      data_ = pool->allocate<T>(size);
    }
    /* else do nothing - we already own our data */
  }
}



/// Rebind the memory referred to by Dense_storage object
///
/// Requires:
///   SIZE to be size object was constructed with.

template <typename T>
void
Dense_storage<Cmplx_split_fmt, vsip::complex<T> >::impl_rebind(
  Memory_pool* pool,
  length_type  size,
  type         buffer)
{
  if (buffer.first != NULL && buffer.second != NULL)
  {
    if (alloc_data_)
    {
      pool->deallocate(real_data_, size);
      pool->deallocate(imag_data_, size);
    }
    
    alloc_data_ = false;
    real_data_  = buffer.first;
    imag_data_  = buffer.second;
  }
  else // (buffer == NULL)
  {
    if (!alloc_data_)
    {
      alloc_data_ = true;
      real_data_ = pool->allocate<T>(size);
      imag_data_ = pool->allocate<T>(size);
    }
    /* else do nothing - we already own our data */
  }
}



/***********************************************************************
  Definitions - Dense_impl
***********************************************************************/

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::Dense_impl(
  Domain<Dim> const& dom,
  MapT const&        map)
VSIP_THROW((std::bad_alloc))
  : storage_type(map.impl_pool(), applied_layout_type(dom).total_size()),
    layout_     (dom),
    map_        (map),
    admitted_   (true)
{
  map_.impl_apply(dom);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::Dense_impl(
  Domain<Dim> const& dom,
  T                  val,
  MapT const&        map)
VSIP_THROW((std::bad_alloc))
  : storage_type(map.impl_pool(), applied_layout_type(dom).total_size(), val),
    layout_     (dom),
    map_        (map),
    admitted_   (true)
{
  map_.impl_apply(dom);
}



// User-storage constructor.

// The buffer given to storage_ is generated from user_data.
//
// If user_data is compatible with the storage_ type
//   (i.e. user_storage is array_format or interleaved_format and
//   storage_ is Cmplx_inter_fmt, or
//   user_storage is split_format and storage_ is Cmplx_split_fmt)
// then storage_ will use the user provided storage,
// else storage_ will allocate its own memory of the same size.
//
// On admit and release, storage_.is_alloc() is used to determine if
// storage_ and user_data_ are the same (in which case no copy is
// necessary to update), or different (copy necessary).

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::Dense_impl(
  Domain<Dim> const&     dom,
  User_storage<T> const& user_data,
  MapT const&            map)
VSIP_THROW((std::bad_alloc))
  : storage_type(map.impl_pool(), applied_layout_type(dom).total_size(),
		 user_data.as_storage(complex_type())),
    layout_     (dom),
    user_data_  (user_data),
    map_        (map),
    admitted_   (false)
{
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
T
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::get(
  Index<Dim> const& idx)
  const VSIP_NOTHROW
{
  for (dimension_type d=0; d<Dim; ++d)
    assert(idx[d] < layout_.size(d));
  return this->get(layout_.index(idx));
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::put(
  Index<Dim> const& idx,
  T                 val)
  VSIP_NOTHROW
{
  for (dimension_type d=0; d<Dim; ++d)
    assert(idx[d] < layout_.size(d));
  return this->put(layout_.index(idx), val);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
typename Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::reference_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::impl_ref(
  Index<Dim> const& idx) VSIP_NOTHROW
{
  for (dimension_type d=0; d<Dim; ++d)
    assert(idx[d] < layout_.size(d));
  return this->impl_ref(layout_.index(idx));
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
typename Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::const_reference_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::impl_ref(
  Index<Dim> const& idx) const VSIP_NOTHROW
{
  for (dimension_type d=0; d<Dim; ++d)
    assert(idx[d] < layout_.size(d));
  return this->impl_.ref(layout_.index(idx));
}



/// Return the total size of the block.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
length_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::size() const VSIP_NOTHROW
{
  length_type retval = layout_.size(0);
  for (dimension_type d=1; d<Dim; ++d)
    retval *= layout_.size(d);
  return retval;
}



/// Return the size of the block in a specific dimension.

/// Requires:
///   BLOCK_DIM selects which block-dimensionality (BLOCK_DIM == 1).
///   DIM is the dimension whose length to return (0 <= DIM < BLOCK_DIM).
/// Returns:
///   The size of dimension DIM.
template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
length_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::size(
  dimension_type block_dim,
  dimension_type d)
  const VSIP_NOTHROW
{
  assert((block_dim == 1 || block_dim == Dim) && (d < block_dim));
  return (block_dim == 1) ? this->size() : this->layout_.size(d);
}



// Requires:
//   DIM is a valid dimensionality supported by block (DIM == 1 or 2)
//   D is a dimension, less than DIM.
// Returns
//   The stride in dimension D, for dimensionality DIM.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
stride_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::impl_stride(
  dimension_type block_dim, dimension_type d)
  const VSIP_NOTHROW
{
  assert((block_dim == 1 || block_dim == Dim) && (d < block_dim));

  return (block_dim == 1) ? 1 : layout_.stride(d);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::admit(bool update)
  VSIP_NOTHROW
{
  if (!this->admitted_ && this->is_alloc() && update)
  {
    for (index_type i=0; i<this->size(); ++i)
      this->put(i, this->user_data_.get(i));
  }

  this->admitted_ = true;
}



/// Release user-storage from VSIPL++ control to user control.
///
/// Note:
///  - It is not an error to release a block multiple times,
///    but it may signify an application programming error.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::release(bool update)
  VSIP_NOTHROW
{
  if (this->user_data_.format() == no_user_format)
    return;
  if (this->admitted_ && this->is_alloc() && update)
  {
    for (index_type i=0; i<this->size(); ++i)
      this->user_data_.put(i, this->get(i));
  }

  this->admitted_ = false;
}



/// Release user-storage and return pointer (array format).
///
/// Requires:
///   THIS to be a block with either array_format user storage,
///      or no_user_format user storage.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::release(bool update, T*& pointer)
  VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == array_format);
  this->release(update);
  this->user_data_.find(pointer);
}



/// Release user-storeage and return pointer (interleaved format).
///
/// Requires:
///   THIS to be a block with either interleaved_format user storage,
///      or no_user_format user storage.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::release(bool update, uT*& pointer)
  VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == interleaved_format);
  this->release(update);
  this->user_data_.find(pointer);
}



/// Release user-storeage and return pointers (split format).
///
/// Requires:
///   THIS to be a block with either split_format user storage,
///      or no_user_format user storage.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::release(
  bool update,
  uT*& real_pointer,
  uT*& imag_pointer)
VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == split_format);
  this->release(update);
  this->user_data_.find(real_pointer, imag_pointer);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::find(T*& pointer)
  VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == array_format);
  this->user_data_.find(pointer);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::find(uT*& pointer)
  VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == interleaved_format);
  this->user_data_.find(pointer);
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::find(uT*& real_pointer, uT*& imag_pointer)
  VSIP_NOTHROW
{
  assert(this->user_storage() == no_user_format ||
	 this->user_storage() == split_format);
  this->user_data_.find(real_pointer, imag_pointer);
}



/// Rebind user-storage to a new array.
///
/// Requires
///   THIS must be a block with array_format user storage that is
///      currently released.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::rebind(T* pointer)
  VSIP_NOTHROW
{
  assert(!this->admitted() && this->user_storage() == array_format);
  this->user_data_.rebind(pointer);
  this->impl_rebind(
		map_.impl_pool(), layout_.total_size(),
		this->user_data_.as_storage(complex_type()));
}



/// Rebind user-storage to a new interleaved array.
///
/// Requires
///   THIS must be a block with interleaved_format or split_format
///      user storage that is currently released.
///
/// Note:
///   When changing user storage from INTERLEAVED to ARRAY and
///    visa-versa, a rebind() will allocate/deallocate memory.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::rebind(uT* pointer)
  VSIP_NOTHROW
{
  assert(!this->admitted() &&
	 (this->user_storage() == split_format ||
	  this->user_storage() == interleaved_format));
  this->user_data_.rebind(pointer);
  this->impl_rebind(map_.impl_pool(), layout_.total_size(),
		    this->user_data_.as_storage(complex_type()));
}



/// Rebind user-storage to new split arrays.
///
/// Requires
///   THIS must be a block with interleaved_format or split_format
///      user storage that is currently released.
///
/// Note:
///   When changing user storage from INTERLEAVED to ARRAY and
///    visa-versa, a rebind() will allocate/deallocate memory.

template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline void
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::rebind(
  uT* real_pointer,
  uT* imag_pointer)
  VSIP_NOTHROW
{
  assert(!this->admitted() &&
	 (this->user_storage() == split_format ||
	  this->user_storage() == interleaved_format));
  this->user_data_.rebind(real_pointer, imag_pointer);
  this->impl_rebind(map_.impl_pool(),
		    layout_.total_size(),
		    this->user_data_.as_storage(complex_type()));
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
user_storage_type
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::user_storage()
  const VSIP_NOTHROW
{
  return this->user_data_.format();
}



template <dimension_type Dim,
	  typename       T,
	  typename       OrderT,
	  typename       ComplexFmt,
	  typename       MapT>
inline
bool
Dense_impl<Dim, T, OrderT, ComplexFmt, MapT>::admitted()
  const VSIP_NOTHROW
{
  return this->admitted_;
}

} // namespace vsip::impl

} // namespace vsip

#undef VSIP_IMPL_DENSE_CMPLX_FMT

#endif // VSIP_DENSE_HPP
