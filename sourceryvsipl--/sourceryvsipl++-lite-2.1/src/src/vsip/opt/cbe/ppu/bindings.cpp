/* Copyright (c) 2006, 2007, 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/cbe/ppu/bindings.cpp
    @author  Stefan Seefeld
    @date    2006-12-29
    @brief   VSIPL++ Library: Wrappers and traits to bridge with IBMs CBE SDK.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#include <algorithm>

#include <vsip/core/config.hpp>
#include <vsip/math.hpp>
#include <vsip/opt/cbe/ppu/bindings.hpp>
#include <vsip/opt/cbe/ppu/task_manager.hpp>
#include <vsip/opt/cbe/ppu/util.hpp>
#include <vsip/opt/cbe/ppu/plugin.hpp>
#include <vsip/opt/cbe/overlay_params.h>
#include <vsip/opt/cbe/vmmul_params.h>
#include <vsip/opt/cbe/vmul_params.h>
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



/***********************************************************************
  General plugin helper class
***********************************************************************/

// Generic elementwise binary vector operation using plugins.

template <typename T>
struct Vbin_plugin
{
  Vbin_plugin() {}

  Vbin_plugin(
    overlay_op_type          op,
    char*                    code_ea,
    int                      code_size,
    std::pair<T*, T*> const& A,
    std::pair<T*, T*> const& B,
    std::pair<T*, T*> const& R,
    length_type              len)
  {
    length_type const chunk_size = 1024;

    assert(is_dma_addr_ok(A.first) && is_dma_addr_ok(A.second));
    assert(is_dma_addr_ok(B.first) && is_dma_addr_ok(B.second));
    assert(is_dma_addr_ok(R.first) && is_dma_addr_ok(R.second));

    Overlay_params params;
    Vmul_split_params* vp = &params.zvbin;

    vp->code_ea   = (uintptr_t)code_ea;
    vp->code_size = code_size;
    vp->cmd       = op;

    vp->length       = chunk_size;
    vp->a_blk_stride = chunk_size;
    vp->b_blk_stride = chunk_size;
    vp->r_blk_stride = chunk_size;

    Task_manager *mgr = Task_manager::instance();

    task = mgr->reserve_iobuf<Plugin_tag, void>
      (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

    assert(sizeof(T)*4*chunk_size <= VSIP_IMPL_OVERLAY_BUFFER_SIZE);
    assert(8 < VSIP_IMPL_OVERLAY_DTL_SIZE);

    length_type chunks = len / chunk_size;
    length_type spes   = mgr->num_spes();
    length_type chunks_per_spe = chunks / spes;
    assert(chunks_per_spe * spes <= chunks);
    
    T const* a_re_ptr = A.first;
    T const* a_im_ptr = A.second;
    T const* b_re_ptr = B.first;
    T const* b_im_ptr = B.second;
    T*       r_re_ptr = R.first;
    T*       r_im_ptr = R.second;

    for (index_type i=0; i<spes && i<chunks; ++i)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type my_chunks = (i < chunks % spes) ? chunks_per_spe + 1
                                                  : chunks_per_spe;
      
      Workblock block = task.create_workblock(my_chunks);
      vp->a_re_ptr = (uintptr_t)a_re_ptr;
      vp->a_im_ptr = (uintptr_t)a_im_ptr;
      vp->b_re_ptr = (uintptr_t)b_re_ptr;
      vp->b_im_ptr = (uintptr_t)b_im_ptr;
      vp->r_re_ptr = (uintptr_t)r_re_ptr;
      vp->r_im_ptr = (uintptr_t)r_im_ptr;
      block.set_parameters(params);
      block.enqueue();

      a_re_ptr += my_chunks*chunk_size;
      a_im_ptr += my_chunks*chunk_size;
      b_re_ptr += my_chunks*chunk_size;
      b_im_ptr += my_chunks*chunk_size;
      r_re_ptr += my_chunks*chunk_size;
      r_im_ptr += my_chunks*chunk_size;
      len -= my_chunks * chunk_size;
    }

    // Cleanup leftover data that doesn't fit into a full chunk.

    // First, handle data that can be DMA'd to the SPEs.  DMA granularity
    // is 16 bytes for sizes 16 bytes or larger.

    length_type const granularity = VSIP_IMPL_CBE_DMA_GRANULARITY / sizeof(T);
    
    if (len >= granularity)
    {
      vp->length = (len / granularity) * granularity;
      assert(is_dma_size_ok(vp->length*sizeof(T)));
      Workblock block = task.create_workblock(1);
      vp->a_re_ptr = (uintptr_t)a_re_ptr;
      vp->a_im_ptr = (uintptr_t)a_im_ptr;
      vp->b_re_ptr = (uintptr_t)b_re_ptr;
      vp->b_im_ptr = (uintptr_t)b_im_ptr;
      vp->r_re_ptr = (uintptr_t)r_re_ptr;
      vp->r_im_ptr = (uintptr_t)r_im_ptr;
      block.set_parameters(params);
      block.enqueue();
      len -= vp->length;
    }

    leftover = len;
  }

  Vbin_plugin(
    overlay_op_type op,
    char*           code_ea,
    int             code_size,
    T const*        A,
    T const*        B,
    T*              R,
    length_type     len)
  {
    length_type const chunk_size = 1024;

    assert(is_dma_addr_ok(A));
    assert(is_dma_addr_ok(B));
    assert(is_dma_addr_ok(R));

    Overlay_params params;
    Vmul_params* vp = &params.vbin;

    vp->code_ea   = (uintptr_t)code_ea;
    vp->code_size = code_size;
    vp->cmd       = op;

    vp->length       = chunk_size;
    vp->a_blk_stride = chunk_size;
    vp->b_blk_stride = chunk_size;
    vp->r_blk_stride = chunk_size;

    Task_manager *mgr = Task_manager::instance();

    task = mgr->reserve_iobuf<Plugin_tag, void>
      (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

    assert(sizeof(T)*4*chunk_size <= VSIP_IMPL_OVERLAY_BUFFER_SIZE);
    assert(8 < VSIP_IMPL_OVERLAY_DTL_SIZE);

    length_type chunks = len / chunk_size;
    length_type spes   = mgr->num_spes();
    length_type chunks_per_spe = chunks / spes;
    assert(chunks_per_spe * spes <= chunks);
    
    T const* a_ptr = A;
    T const* b_ptr = B;
    T*       r_ptr = R;

    for (index_type i=0; i<spes && i<chunks; ++i)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type my_chunks = (i < chunks % spes) ? chunks_per_spe + 1
                                                  : chunks_per_spe;
      
      Workblock block = task.create_workblock(my_chunks);
      vp->a_ptr = (uintptr_t)a_ptr;
      vp->b_ptr = (uintptr_t)b_ptr;
      vp->r_ptr = (uintptr_t)r_ptr;
      block.set_parameters(params);
      block.enqueue();

      a_ptr += my_chunks*chunk_size;
      b_ptr += my_chunks*chunk_size;
      r_ptr += my_chunks*chunk_size;
      len -= my_chunks * chunk_size;
    }

    // Cleanup leftover data that doesn't fit into a full chunk.

    // First, handle data that can be DMA'd to the SPEs.  DMA granularity
    // is 16 bytes for sizes 16 bytes or larger.

    length_type const granularity = VSIP_IMPL_CBE_DMA_GRANULARITY / sizeof(T);
    
    if (len >= granularity)
    {
      vp->length = (len / granularity) * granularity;
      assert(is_dma_size_ok(vp->length*sizeof(T)));
      Workblock block = task.create_workblock(1);
      vp->a_ptr = (uintptr_t)a_ptr;
      vp->b_ptr = (uintptr_t)b_ptr;
      vp->r_ptr = (uintptr_t)r_ptr;
      block.set_parameters(params);
      block.enqueue();
      len -= vp->length;
    }

    leftover = len;
  }

  ~Vbin_plugin() { task.sync(); }

  Task        task;
  length_type leftover;
};



/***********************************************************************
  Vector multiply
***********************************************************************/

template <typename T>
void vmul(T const* A, T const* B, T* R, length_type len)
{
  static char* code = 0;
  static int   size;
  if (code == 0) load_plugin(code, size, "plugin", "xvmul_f");

  overlay_op_type op = Type_equal<T, float>::value ? overlay_vmul_f :
                                                     overlay_cvmul_f;

  Vbin_plugin<T> plugin(op, code, size, A, B, R, len);

  if (plugin.leftover)
    for (index_type i=len-plugin.leftover; i<len; ++i)
      R[i] = A[i] * B[i];
}

template void vmul(float const* A, float const* B, float* R, length_type len);
template void vmul(complex<float> const* A, complex<float> const* B, complex<float>* R, length_type len);

template <typename T>
void vmul(
  std::pair<T*, T*> const& A,
  std::pair<T*, T*> const& B,
  std::pair<T*, T*> const& R,
  length_type              len)
{
  static char* code = 0;
  static int   size;
  if (code == 0) load_plugin(code, size, "plugin", "xvmul_f");
  
  Vbin_plugin<T> plugin(overlay_zvmul_f, code, size, A, B, R, len);

  if (plugin.leftover)
  {
    for (index_type i=len-plugin.leftover; i<len; ++i)
    {
      T t         = A.first[i] * B.first[i]  - A.second[i] * B.second[i];
      R.second[i] = A.first[i] * B.second[i] + A.second[i] * B.first[i];
      R.first[i]  = t;
    }
  }
}

template 
void vmul(
  std::pair<float*, float*> const& A,
  std::pair<float*, float*> const& B,
  std::pair<float*, float*> const& R,
  length_type              len);



/***********************************************************************
  Vector add
***********************************************************************/

template <typename T>
void vadd(T const* A, T const* B, T* R, length_type len)
{
  static char* code = 0;
  static int   size;
  if (code == 0) load_plugin(code, size, "plugin", "xvadd_f");

  overlay_op_type op = Is_complex<T>::value ? overlay_cvadd_f :
                                              overlay_vadd_f;

  Vbin_plugin<T> plugin(op, code, size, A, B, R, len);

  if (plugin.leftover)
    for (index_type i=len-plugin.leftover; i<len; ++i)
      R[i] = A[i] + B[i];
}

template void vadd(float const* A, float const* B, float* R, length_type len);
template void vadd(complex<float> const* A, complex<float> const* B, complex<float>* R, length_type len);



template <typename T>
void vadd(
  std::pair<T*, T*> const& A,
  std::pair<T*, T*> const& B,
  std::pair<T*, T*> const& R,
  length_type              len)
{
  static char* code = 0;
  static int   size;
  if (code == 0) load_plugin(code, size, "plugin", "xvadd_f");

  Vbin_plugin<T> plugin(overlay_zvadd_f, code, size, A, B, R, len);

  if (plugin.leftover)
  {
    for (index_type i=len-plugin.leftover; i<len; ++i)
    {
      R.first[i]  = A.first[i] + B.first[i];
      R.second[i] = A.second[i] + B.second[i];
    }
  }
}

template 
void vadd(
  std::pair<float*, float*> const& A,
  std::pair<float*, float*> const& B,
  std::pair<float*, float*> const& R,
  length_type              len);



/***********************************************************************
  Vector-matrix multiply
***********************************************************************/

// Scalar row-wise vmmul
//
// Computes: R(r, c) = V(c) * M(r, c)
//
// Arguments
//   - m_stride is the distance between input rows (or cols)
//   - r_stride is the distance between output rows (or cols)
//   - lines expresses the number of rows (cols)
//   - length is the size of the input vector 

template <typename T>
void vmmul_row(
  T const* V, T const* M, T* R, 
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length)
{
  assert(length >= VSIP_IMPL_MIN_VMMUL_SIZE);

  static char* code_ea = 0;
  static int   code_size;
  Overlay_params params;

  if (code_ea == 0) load_plugin(code_ea, code_size, "plugin", "chalfast_f");
  Vmmul_params* vp = &params.vmmul;

  vp->code_ea   = (uintptr_t)code_ea;
  vp->code_size = code_size;
  vp->cmd       = Is_complex<T>::value ? overlay_cvmmul_row_f
                                       : overlay_vmmul_row_f;
  vp->input_stride  = m_stride;
  vp->output_stride = r_stride;

  Task_manager* mgr = Task_manager::instance();

  Task task = mgr->reserve_iobuf<Plugin_tag, void>
  (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

  length_type spes            = mgr->num_spes();
  length_type vectors_per_spe = lines / spes;
  assert(vectors_per_spe * spes <= lines);
  length_type const granularity = VSIP_IMPL_CBE_DMA_GRANULARITY / sizeof(T);

  // break large vmmul into smaller vmmuls of MAX size
  index_type pos = 0;
  while (pos<length && (length-pos > VSIP_IMPL_MIN_VMMUL_SIZE))
  {
    length_type my_length = std::min<length_type>(
      length - pos, VSIP_IMPL_MAX_VMMUL_SIZE);

    if (my_length % granularity != 0)
      my_length -= my_length % granularity;

    vp->length = my_length;
    assert(my_length >= VSIP_IMPL_MIN_VMMUL_SIZE);

    vp->ea_input_vector  = ea_from_ptr(V + pos);
    vp->ea_input_matrix  = ea_from_ptr(M + pos);
    vp->ea_output_matrix = ea_from_ptr(R + pos);
    vp->shift = 0;

    for (index_type i=0; i<spes && i<lines; ++i)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type lines_per_spe = (i < lines % spes) ? vectors_per_spe + 1
                                                     : vectors_per_spe;

      Workblock block = task.create_workblock(lines_per_spe);
      block.set_parameters(params);
      block.enqueue();

      vp->ea_input_matrix  += sizeof(T) * lines_per_spe * m_stride;
      vp->ea_output_matrix += sizeof(T) * lines_per_spe * r_stride;
    }

    pos += my_length;
  }

  // Cleanup
  if (length-pos > 0)
  {
    for (index_type r=0; r<lines; ++r)
    {
      for (index_type c=pos; c<length; ++c)
      {
	R[c] = V[c] * M[c];
      }
      M += m_stride;
      R += r_stride;
    }
  }
  task.sync();
}



//template void vmmul_row(
//  float const* V, float const* M, float* R, 
//  stride_type m_stride, stride_type r_stride, length_type length, length_type lines);
template void vmmul_row(
  std::complex<float> const* V, std::complex<float> const* M, std::complex<float>* R, 
  stride_type m_stride, stride_type r_stride, length_type length, length_type lines);



// Split-complex row-wise vmmul
//
// Computes: R(r, c) = V(c) * M(r, c)
//
// Arguments
//   - m_stride is the distance between input rows (or cols)
//   - r_stride is the distance between output rows (or cols)
//   - lines expresses the number of rows (cols)
//   - length is the size of the input vector 

template <typename T>
void vmmul_row(
  std::pair<T*, T*> const& V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const& R,
  stride_type              m_stride,
  stride_type              r_stride, 
  length_type              lines,
  length_type              length)
// Result R = V * M
//   where
//   - m_stride is the distance between input rows (or cols)
//   - r_stride is the distance between output rows (or cols)
//   - lines expresses the number of rows (cols)
//   - length is the size of the input vector 
{
  static char* code_ea = 0;
  static int   code_size;
  Overlay_params params;

  if (code_ea == 0) load_plugin(code_ea, code_size, "plugin", "zhalfast_f");
  Vmmul_split_params* vp = &params.zvmmul;

  vp->code_ea       = (uintptr_t)code_ea;
  vp->code_size     = code_size;
  vp->cmd           = overlay_zvmmul_row_f;
  vp->input_stride  = m_stride;
  vp->output_stride = r_stride;

  Task_manager* mgr = Task_manager::instance();

  Task task = mgr->reserve_iobuf<Plugin_tag, void>
  (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

  length_type spes            = mgr->num_spes();
  length_type vectors_per_spe = lines / spes;
  assert(vectors_per_spe * spes <= lines);
  length_type const granularity = VSIP_IMPL_CBE_DMA_GRANULARITY / sizeof(T);

  // break large vmmul into smaller vmmuls of MAX size
  index_type pos = 0;
  while (pos<length && (length-pos > VSIP_IMPL_MIN_VMMUL_SIZE))
  {
    length_type my_length = std::min<length_type>(
      length - pos, VSIP_IMPL_MAX_VMMUL_SIZE);

    if (my_length % granularity != 0)
      my_length -= my_length % granularity;

    vp->length = my_length;
    assert(my_length >= VSIP_IMPL_MIN_VMMUL_SIZE);

    vp->ea_input_vector_re  = ea_from_ptr(V.first  + pos);
    vp->ea_input_vector_im  = ea_from_ptr(V.second + pos);
    vp->ea_input_matrix_re  = ea_from_ptr(M.first  + pos);
    vp->ea_input_matrix_im  = ea_from_ptr(M.second + pos);
    vp->ea_output_matrix_re = ea_from_ptr(R.first  + pos);
    vp->ea_output_matrix_im = ea_from_ptr(R.second + pos);
    vp->shift = 0;

    for (index_type i=0; i<spes && i<lines; ++i)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type lines_per_spe = (i < lines % spes) ? vectors_per_spe + 1
                                                     : vectors_per_spe;

      Workblock block = task.create_workblock(lines_per_spe);
      block.set_parameters(params);
      block.enqueue();

      vp->ea_input_matrix_re  += sizeof(T) * lines_per_spe * m_stride;
      vp->ea_input_matrix_im  += sizeof(T) * lines_per_spe * m_stride;
      vp->ea_output_matrix_re += sizeof(T) * lines_per_spe * r_stride;
      vp->ea_output_matrix_im += sizeof(T) * lines_per_spe * r_stride;
    }

    pos += my_length;
  }
  if (length-pos > 0)
  {
    T* Vr = V.first;
    T* Vi = V.second;
    T* Mr = M.first;
    T* Mi = M.second;
    T* Rr = R.first;
    T* Ri = R.second;

    for (index_type r=0; r<lines; ++r)
    {
      for (index_type c=pos; c<length; ++c)
      {
	float tmp = Vr[c] * Mr[c] - Vi[c] * Mi[c];
	Ri[c]     = Vr[c] * Mi[c] + Vi[c] * Mr[c];
	Rr[c]     = tmp;
      }
      Mr += m_stride;
      Mi += m_stride;
      Rr += r_stride;
      Ri += r_stride;
    }
  }
  task.sync();
}

template
void vmmul_row(
  std::pair<float*, float*> const& V,
  std::pair<float*, float*> const& M,
  std::pair<float*, float*> const&             R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);



// Split-complex real-scalar column-wise vmmul
//
// Computes: R(r, c) = V(r) * M(r, c)
//
// Arguments
//   - m_stride is the distance between input rows (or cols)
//   - r_stride is the distance between output rows (or cols)
//   - lines expresses the number of rows (cols)
//   - length is the size of the input vector 

template <typename T>
void vmmul_col(
  T const*                 V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const& R,
  stride_type              m_stride,
  stride_type              r_stride, 
  length_type              lines,
  length_type              length)
{
#if USE_FUNCTIONAL_VERSION
  T* Mr = M.first;
  T* Mi = M.second;
  T* Rr = R.first;
  T* Ri = R.second;

  for (index_type r=0; r<lines; ++r)
  {
    for (index_type c=0; c<length; ++c)
    {
      Rr[c] = V[r] * Mr[c];
      Ri[c] = V[r] * Mi[c];
    }
    Mr += m_stride;
    Mi += m_stride;
    Rr += r_stride;
    Ri += r_stride;
  }
#else
  static char* code_ea = 0;
  static int   code_size;
  Overlay_params params;

  if (code_ea == 0) load_plugin(code_ea, code_size, "plugin", "zrvmmul_col_f");
  Vmmul_split_params* vp = &params.zvmmul;

  vp->code_ea   = (uintptr_t)code_ea;
  vp->code_size = code_size;
  vp->cmd       = overlay_zrvmmul_col_f;

  Task_manager* mgr = Task_manager::instance();

  Task task = mgr->reserve_iobuf<Plugin_tag, void>
  (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

  assert(4*sizeof(T)*VSIP_IMPL_MAX_VMMUL_SIZE/2 <= VSIP_IMPL_OVERLAY_BUFFER_SIZE);
  assert(8 < VSIP_IMPL_OVERLAY_DTL_SIZE);

  length_type spes          = mgr->num_spes();
  length_type lines_per_spe = lines / spes;
  if (lines_per_spe == 0 || lines_per_spe % 4 != 0)
    lines_per_spe += 4 - (lines_per_spe % 4);

  for (length_type pos=0; pos<length; pos += 2048)
  {
    length_type this_length = std::min<length_type>(length-pos, 2048);

    vp->length = this_length;
    vp->ea_input_vector_re  = ea_from_ptr(V);
    vp->ea_input_vector_im  = 0;

    if (vp->ea_input_vector_re % 16 != 0)
    {
      int shift = (vp->ea_input_vector_re % 16);
      vp->ea_input_vector_re -= shift;
      vp->shift = shift / sizeof(float);
    }
    else
      vp->shift = 0;

    vp->input_stride        = m_stride;
    vp->output_stride       = r_stride;
    vp->ea_input_matrix_re  = ea_from_ptr(M.first  + pos);
    vp->ea_input_matrix_im  = ea_from_ptr(M.second + pos);
    vp->ea_output_matrix_re = ea_from_ptr(R.first  + pos);
    vp->ea_output_matrix_im = ea_from_ptr(R.second + pos);

    for (length_type rem_lines=lines; rem_lines > 0;)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type my_lines = (lines_per_spe < rem_lines) ? lines_per_spe
                                                         : rem_lines;

      Workblock block = task.create_workblock(my_lines);
      block.set_parameters(params);
      block.enqueue();

      vp->ea_input_vector_re  += sizeof(T) * my_lines;
      vp->ea_input_matrix_re  += sizeof(T) * my_lines * m_stride;
      vp->ea_input_matrix_im  += sizeof(T) * my_lines * m_stride;
      vp->ea_output_matrix_re += sizeof(T) * my_lines * r_stride;
      vp->ea_output_matrix_im += sizeof(T) * my_lines * r_stride;
      rem_lines -= my_lines;
    }
  }
  task.sync();
#endif
}

template
void vmmul_col(
  float const*                     V,
  std::pair<float*, float*> const& M,
  std::pair<float*, float*> const& R,
  stride_type                      m_stride,
  stride_type                      r_stride, 
  length_type                      lines, 
  length_type                      length);



// Split-complex real-scalar column-wise vmmul
//
// Computes R(r, c) = V(r) * M(r, c)
//
// Arguments
//   - m_stride is the distance between input rows (or cols)
//   - r_stride is the distance between output rows (or cols)
//   - lines expresses the number of rows (cols)
//   - length is the size of the input vector 

template <typename T>
void vmmul_col(
  std::pair<T*, T*> const& V,
  std::pair<T*, T*> const& M,
  std::pair<T*, T*> const& R,
  stride_type              m_stride,
  stride_type              r_stride, 
  length_type              lines,
  length_type              length)
{
#if USE_FUNCTIONAL_VERSION
  T* Vr = V.first;
  T* Vi = V.second;
  T* Mr = M.first;
  T* Mi = M.second;
  T* Rr = R.first;
  T* Ri = R.second;

  for (index_type r=0; r<lines; ++r)
  {
    for (index_type c=0; c<length; ++c)
    {
      T tmp = Vr[r] * Mr[c] - Vi[r] * Mi[c];
      Ri[c] = Vr[r] * Mi[c] + Vi[r] * Mr[c];
      Rr[c] = tmp;
    }
    Mr += m_stride;
    Mi += m_stride;
    Rr += r_stride;
    Ri += r_stride;
  }
#else
  static char* code_ea = 0;
  static int   code_size;
  Overlay_params params;

  if (code_ea == 0) load_plugin(code_ea, code_size, "plugin", "zvmmul_col_f");
  Vmmul_split_params* vp = &params.zvmmul;

  vp->code_ea   = (uintptr_t)code_ea;
  vp->code_size = code_size;
  vp->cmd       = overlay_zvmmul_col_f;

  Task_manager* mgr = Task_manager::instance();

  Task task = mgr->reserve_iobuf<Plugin_tag, void>
  (VSIP_IMPL_OVERLAY_STACK_SIZE, sizeof(Overlay_params),
       VSIP_IMPL_OVERLAY_BUFFER_SIZE,
       VSIP_IMPL_OVERLAY_DTL_SIZE);

  assert(4*sizeof(T)*VSIP_IMPL_MAX_VMMUL_SIZE/2 <= VSIP_IMPL_OVERLAY_BUFFER_SIZE);
  assert(8 < VSIP_IMPL_OVERLAY_DTL_SIZE);

  length_type spes          = mgr->num_spes();
  length_type lines_per_spe = lines / spes;
  if (lines_per_spe == 0 || lines_per_spe % 4 != 0)
    lines_per_spe += 4 - (lines_per_spe % 4);

  for (length_type pos=0; pos<length; pos += 2048)
  {
    length_type this_length = std::min<length_type>(length-pos, 2048);

    vp->length = this_length;
    vp->ea_input_vector_re  = ea_from_ptr(V.first);
    vp->ea_input_vector_im  = ea_from_ptr(V.second);

    assert(vp->ea_input_vector_re % 16 == vp->ea_input_vector_im % 16);
    if (vp->ea_input_vector_re % 16 != 0)
    {
      int shift = (vp->ea_input_vector_re % 16);
      vp->ea_input_vector_re -= shift;
      vp->ea_input_vector_im -= shift;
      vp->shift = shift / sizeof(float);
    }
    else
      vp->shift = 0;

    vp->input_stride        = m_stride;
    vp->output_stride       = r_stride;
    vp->ea_input_matrix_re  = ea_from_ptr(M.first  + pos);
    vp->ea_input_matrix_im  = ea_from_ptr(M.second + pos);
    vp->ea_output_matrix_re = ea_from_ptr(R.first  + pos);
    vp->ea_output_matrix_im = ea_from_ptr(R.second + pos);

    for (length_type rem_lines=lines; rem_lines > 0;)
    {
      // If chunks don't divide evenly, give the first SPEs one extra.
      length_type my_lines = (lines_per_spe < rem_lines) ? lines_per_spe
                                                         : rem_lines;

      Workblock block = task.create_workblock(my_lines);
      block.set_parameters(params);
      block.enqueue();

      vp->ea_input_vector_re  += sizeof(T) * my_lines;
      vp->ea_input_vector_im  += sizeof(T) * my_lines;
      vp->ea_input_matrix_re  += sizeof(T) * my_lines * m_stride;
      vp->ea_input_matrix_im  += sizeof(T) * my_lines * m_stride;
      vp->ea_output_matrix_re += sizeof(T) * my_lines * r_stride;
      vp->ea_output_matrix_im += sizeof(T) * my_lines * r_stride;
      rem_lines -= my_lines;
    }
  }
  task.sync();
#endif
}

template
void vmmul_col(
  std::pair<float*, float*> const& V,
  std::pair<float*, float*> const& M,
  std::pair<float*, float*> const& R,
  stride_type m_stride, stride_type r_stride, 
  length_type lines, length_type length);

} // namespace vsip::impl::cbe
} // namespace vsip::impl
} // namespace vsip
