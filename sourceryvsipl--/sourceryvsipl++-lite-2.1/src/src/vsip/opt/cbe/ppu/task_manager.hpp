/* Copyright (c) 2007, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/task_manager.hpp
    @author  Stefan Seefeld
    @date    2007-01-31
    @brief   VSIPL++ Library: Wrappers and traits to bridge with IBMs CBE SDK.
*/

#ifndef VSIP_OPT_CBE_PPU_TASK_MANAGER_HPP
#define VSIP_OPT_CBE_PPU_TASK_MANAGER_HPP

#if VSIP_IMPL_REF_IMPL
# error "vsip/opt files cannot be used as part of the reference impl."
#endif

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/expr/operations.hpp>
#include <vsip/core/argv_utils.hpp>
#include <vsip/opt/cbe/ppu/alf.hpp>
extern "C"
{
#include <libspe2.h>
}
/***********************************************************************
  Declarations
***********************************************************************/

namespace vsip
{
namespace impl
{

struct Plugin_tag;
struct Fastconv_tag;
struct Fastconvm_tag;
struct Pwarp_tag;


namespace cbe
{

template <typename O, typename S> struct Task_map;

class Task_manager
{
public:
  static Task_manager *instance() { return instance_;}

  static void initialize(int& argc, char**&argv);
  static void finalize() { delete instance_; instance_ = 0;}

  // Return a task for operation O (with signature S).
  // An additional hint parameter may be passed, eventually,
  // indicating how many SPEs to use, etc.
  template <typename O, typename S>
  Task reserve(length_type ssize, // max stack size
	       length_type psize, // parameter buffer size
	       length_type isize, // input buffer size
	       length_type osize, // output buffer size
	       length_type tsize) // number of DMA transfers
  {
    return alf_->create_task("svpp_kernels.so", Task_map<O, S>::image(),
			     ssize, psize, isize, osize, 0, tsize);
  }

  template <typename O, typename S>
  Task reserve_iobuf(
    length_type ssize,  // max stack size
    length_type psize,  // parameter buffer size
    length_type iosize, // input buffer size
    length_type tsize)  // number of DMA transfers
  {
    return alf_->create_task("svpp_kernels.so", Task_map<O, S>::image(),
			     ssize, psize, 0, 0, iosize, tsize);
  }

  length_type num_spes() { return num_spes_; }
  
  ALF* alf_handle() { return alf_; }

private:
  Task_manager(unsigned int num_spes);
  Task_manager(Task_manager const &);
  ~Task_manager();
  Task_manager &operator= (Task_manager const &);

  static Task_manager *instance_;

  ALF*        alf_;
  length_type num_spes_;
};

} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip

// O: operation tag
// S: operation signature
// K: SPE kernel image name
# define DEFINE_TASK(O, S, I)			   \
namespace vsip { namespace impl { namespace cbe {          \
template <>                                                \
struct Task_map<O, S>                                      \
{                                                          \
  static char const *image() { return "alf_" # I "_spu";}  \
};                                                         \
}}}

typedef std::pair<float*,float*> split_float_type;

DEFINE_TASK(Fastconv_tag, void(std::complex<float>, std::complex<float>), fconv_c)
DEFINE_TASK(Fastconv_tag, void(split_float_type, split_float_type), fconv_split_c)
DEFINE_TASK(Fastconvm_tag, void(std::complex<float>, std::complex<float>), fconvm_c)
DEFINE_TASK(Fastconvm_tag, void(split_float_type, split_float_type), fconvm_split_c)
DEFINE_TASK(Pwarp_tag, void(unsigned char, unsigned char), pwarp_ub)
DEFINE_TASK(Plugin_tag, void, plugin)
#endif // VSIP_OPT_CBE_PPU_TASK_MANAGER_HPP
