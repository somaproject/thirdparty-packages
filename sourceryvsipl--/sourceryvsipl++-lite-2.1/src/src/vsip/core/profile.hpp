/* Copyright (c) 2006 by CodeSourcery.  All rights reserved. */

/** @file    vsip/core/profile.hpp
    @author  Stefan Seefeld
    @date    2006-11-24
    @brief   VSIPL++ Library: Profiling routines & classes.
*/

#ifndef VSIP_CORE_PROFILE_HPP
#define VSIP_CORE_PROFILE_HPP

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/noncopyable.hpp>
#include <vsip/core/ops_info.hpp>
#ifndef VSIP_IMPL_REF_IMPL
# include <vsip/opt/profile.hpp>
#endif

namespace vsip
{
namespace impl
{
namespace profile
{
unsigned int const mask =
#ifdef VSIP_IMPL_PROFILER
  VSIP_IMPL_PROFILER
#else
  0
#endif
  ;

#define VSIP_IMPL_PROFILE_MASK_SIGNAL  1
#define VSIP_IMPL_PROFILE_MASK_MATVEC  2
#define VSIP_IMPL_PROFILE_MASK_FNS     4
#define VSIP_IMPL_PROFILE_MASK_USER    8
#define VSIP_IMPL_PROFILE_MASK_PAR     16
#define VSIP_IMPL_PROFILE_MASK_FNS_INT 32

#ifndef VSIP_IMPL_PROFILE
# define VSIP_IMPL_PROFILE(X)
#endif

/// Different operations that may be profiled, each is referred to
/// as a 'feature'.
enum Feature
{
  none,
  signal  = VSIP_IMPL_PROFILE_MASK_SIGNAL, // signal processing (FFT, FIR, etc)
  matvec  = VSIP_IMPL_PROFILE_MASK_MATVEC, // matrix-vector (prod, qr, etc)
  fns     = VSIP_IMPL_PROFILE_MASK_FNS,    // elementwise dispatch (+, -, etc)
  user    = VSIP_IMPL_PROFILE_MASK_USER,   // user defined tag
  par     = VSIP_IMPL_PROFILE_MASK_PAR,    // parallel comms
  fns_int = VSIP_IMPL_PROFILE_MASK_FNS_INT // internal serial_dispatch
};

#if defined(VSIP_IMPL_REF_IMPL)
template <bool>
class Accumulator_base
{
public:
  struct Scope
  {
    Scope(Accumulator_base &) {}
  };
  Accumulator_base(std::string const &, unsigned int = 0) {}
  unsigned int ops() const { return 0;}
  float total() const { return 0.;}
  int   count() const { return 0;}
  float  mops() const { return 0.;}
};

template <bool>
class Scope_base : Non_copyable
{
public:
  Scope_base(std::string const &, int=0) {}
};

enum profiler_mode
{
  pm_trace,
  pm_accum,
  pm_none
};

class Profile
{
public:
  Profile(std::string const &, profiler_mode = pm_accum) {}
};

#endif

template <unsigned int Feature>
class Accumulator : public Accumulator_base<Feature & mask>
{
  typedef Accumulator_base<Feature & mask> base;
public:
  Accumulator(std::string const &n, unsigned int c = 0) : base(n, c) {}
};

template <unsigned int Feature>
class Scope : public Scope_base<Feature & mask>
{
  typedef Scope_base<Feature & mask> base;
public:
  Scope(std::string const &n, int id = 0) : base(n, id) {}
};

} // namespace vsip::impl::profile
} // namespace vsip::impl
} // namespace vsip

#endif
