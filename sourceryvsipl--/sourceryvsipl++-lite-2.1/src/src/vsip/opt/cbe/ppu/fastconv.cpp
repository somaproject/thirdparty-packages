/* Copyright (c) 2007 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/fastconv.cpp
    @author  Don McCoy
    @date    2007-02-23
    @brief   VSIPL++ Library: Wrapper for fast convolution on the SPEs.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <vsip/core/config.hpp>
#include <vsip/core/fns_scalar.hpp>
#include <vsip/core/static_assert.hpp>
#include <vsip/math.hpp>
#include <vsip/opt/cbe/fconv_params.h>
#include <vsip/opt/cbe/ppu/fastconv.hpp>
#include <vsip/opt/cbe/ppu/task_manager.hpp>
#include <vsip/opt/cbe/ppu/util.hpp>
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
namespace cbe
{

// Fast convolution binding for interleaved complex data.


template <dimension_type D,
          typename       T,
	  typename       ComplexFmt>
void
Fastconv_base<D, T, ComplexFmt>::fconv
  (T* in, T* kernel, T* out, length_type rows, length_type length)
{
  Fastconv_params params;

  params.instance_id        = this->instance_id_;
  params.elements           = length;
  params.transform_kernel   = this->transform_kernel_;
  params.ea_kernel          = ea_from_ptr(kernel);
  params.ea_input           = ea_from_ptr(in);
  params.ea_output          = ea_from_ptr(out);
  params.kernel_stride      = length;
  params.input_stride       = length;
  params.output_stride      = length;

  length_type psize = sizeof(params);
  // The stack size takes into account two temporary buffers used
  // to hold the real and imaginary parts of the complex input data.
  length_type stack_size = 4096 + 
    2 * sizeof(T) * cbe::Fastconv_traits<dim, T, complex_type>::max_size;
  typedef typename cbe::Fastconv_traits<dim, T, complex_type>::tag_type tag_type;

  // In the case of a matrix of coefficients (dim == 2), there are two inputs,
  // one row of source data and one row of coefficients.  In the normal case, 
  // where the coeffcients are the same for each row, they are sent in 
  // advance, leaving only one input.
  length_type const num_inputs = (dim == 1 ? 1 : 2);

  Task_manager *mgr = Task_manager::instance();
  Task task = mgr->reserve<tag_type, void(T,T)>
    (stack_size, psize, sizeof(T)*length*num_inputs, sizeof(T)*length, 8);

  length_type spes         = mgr->num_spes();
  length_type rows_per_spe = rows / spes;

  for (index_type i=0; i<spes && i<rows; ++i)
  {
    // If rows don't divide evenly, give the first SPEs one extra.
    length_type my_rows = (i < rows % spes) ? rows_per_spe + 1 : rows_per_spe;

    Workblock block = task.create_workblock(my_rows);
    block.set_parameters(params);
    block.enqueue();
    // Note: for a matrix of coefficients, unique rows are transferred.
    // For the normal case, the address is constant because the same
    // vector is sent once and used repeatedly.
    params.ea_kernel += (dim == 1 ? 0 : sizeof(T) * my_rows * length);
    params.ea_input  += sizeof(T) * my_rows * length;
    params.ea_output += sizeof(T) * my_rows * length;
  }

  task.sync();
}


// Fast convolution binding for split complex data.

template <dimension_type D,
          typename       T,
	  typename       ComplexFmt>
void
Fastconv_base<D, T, ComplexFmt>::fconv(
  std::pair<uT*,uT*> in,
  std::pair<uT*,uT*> kernel,
  std::pair<uT*,uT*> out,
  length_type      rows,
  length_type      length)
{
  Fastconv_split_params params;

  params.instance_id        = this->instance_id_;
  params.elements           = length;
  params.transform_kernel   = this->transform_kernel_;
  params.ea_kernel_re       = ea_from_ptr(kernel.first);
  params.ea_kernel_im       = ea_from_ptr(kernel.second);
  params.ea_input_re        = ea_from_ptr(in.first);
  params.ea_input_im        = ea_from_ptr(in.second);
  params.ea_output_re       = ea_from_ptr(out.first);
  params.ea_output_im       = ea_from_ptr(out.second);
  params.kernel_stride      = length;
  params.input_stride       = length;
  params.output_stride      = length;

  length_type psize = sizeof(params);
  // The split complex version has a smaller stack size requirement,
  // but experimentation indicates that performance is hurt by 
  // reducing this value.
  length_type stack_size = 4096 + 
    sizeof(T)*cbe::Fastconv_traits<dim, T, complex_type>::max_size;
  typedef typename cbe::Fastconv_traits<dim, T, complex_type>::tag_type tag_type;
  length_type const num_inputs = (dim == 1 ? 1 : 2);

  Task_manager *mgr = Task_manager::instance();
  Task task = mgr->reserve<tag_type, void(std::pair<uT*,uT*>,
					  std::pair<uT*,uT*>)>
    (stack_size, psize, sizeof(T)*length*num_inputs, sizeof(T)*length, 8);

  length_type spes         = mgr->num_spes();
  length_type rows_per_spe = rows / spes;

  for (index_type i=0; i<spes && i<rows; ++i)
  {
    // If rows don't divide evenly, give the first SPEs one extra.
    length_type my_rows = (i < rows % spes) ? rows_per_spe + 1 : rows_per_spe;

    Workblock block = task.create_workblock(my_rows);
    block.set_parameters(params);
    block.enqueue();
    params.ea_kernel_re += (dim == 1 ? 0 : sizeof(uT) * my_rows * length);
    params.ea_kernel_im += (dim == 1 ? 0 : sizeof(uT) * my_rows * length);
    params.ea_input_re  += sizeof(uT) * my_rows * length;
    params.ea_input_im  += sizeof(uT) * my_rows * length;
    params.ea_output_re += sizeof(uT) * my_rows * length;
    params.ea_output_im += sizeof(uT) * my_rows * length;
  }

  task.sync();
}



typedef std::complex<float> ctype;
typedef std::pair<float*,float*> ztype;

#define INSTANTIATE_FASTCONV(COEFFS_DIM)                                       \
template void                                                                  \
Fastconv_base<COEFFS_DIM, ctype, Cmplx_inter_fmt>::fconv(                      \
  ctype* in, ctype* kernel, ctype* out, length_type rows, length_type length); \
template<> unsigned int                                                        \
Fastconv_base<COEFFS_DIM, ctype, Cmplx_inter_fmt>::instance_id_counter_ = 0;   \
template void                                                                  \
Fastconv_base<COEFFS_DIM, ctype, Cmplx_split_fmt>::fconv(                      \
  ztype in, ztype kernel, ztype out, length_type rows, length_type length);    \
template<> unsigned int                                                        \
Fastconv_base<COEFFS_DIM, ctype, Cmplx_split_fmt>::instance_id_counter_ = 0; 

INSTANTIATE_FASTCONV(1);
INSTANTIATE_FASTCONV(2);

#undef INSTANTIATE_FASTCONV

          
} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip
