/* Copyright (c) 2008 by CodeSourcery.  All rights reserved.

   This file is available for license from CodeSourcery, Inc. under the terms
   of a commercial license and under the GPL.  It is not part of the VSIPL++
   reference implementation and is not available under the BSD license.
*/
/** @file    vsip/opt/ukernel/cbe_accel/alf_base.hpp
    @author  Jules Bergmann
    @date    2008-06-10
    @brief   VSIPL++ Library: Ukernel base functionality.
*/

/***********************************************************************
  Included Files
***********************************************************************/

#define DEBUG_ALF_BASE 0

#include <utility>
#include <cassert>
#include <complex>

#if DEBUG_ALF_BASE
#  include <cstdio>
#endif
#include <alf_accel.h>
#include <vsip/opt/ukernel/cbe_accel/ukernel.hpp>
#include <vsip/opt/ukernel/ukernel_params.hpp>
#include <vsip/opt/cbe/spu/plugin_functions.h>

// Cell specific
#define DMA_ALIGNMENT 16
#define DMA_SIZE_QUANTUM 16
#define DMA_SIZE_IN_FLOATS (DMA_SIZE_QUANTUM / sizeof(float))
#define DMA_ALIGNMENT_OF(A) ((uintptr_t)(A) & (DMA_ALIGNMENT - 1))
#define IS_DMA_ALIGNED(A) (DMA_ALIGNMENT_OF(A) == 0)
#define IS_DMA_SIZE(S) ((S) % DMA_SIZE_QUANTUM == 0)
#define INCREASE_TO_DMA_SIZE_IN_FLOATS(S) \
  ((S) + (-(unsigned)(S) % DMA_SIZE_IN_FLOATS))


/// Compare two types for "equality".

template <typename T1,
	  typename T2>
struct Type_equal
{
  static bool const value = false;
};

template <typename T>
struct Type_equal<T, T>
{
  static bool const value = true;
  typedef T type;
};



/***********************************************************************
  Definitions
***********************************************************************/

kernel_type ukobj;

typedef Ukernel_params<kernel_type::pre_argc, kernel_type::in_argc,
		       kernel_type::out_argc, kernel_type::param_type>
		param_type;



template <typename T>
inline void
add_entry(
  Plugin_functions* pf,
  void*             entries,
  unsigned int      size,
  alf_data_addr64_t ea)
{
  while (size > 0)
  {
    unsigned int this_size = size*sizeof(T) > 16384 ? 16384/sizeof(T) : size;
    assert(IS_DMA_ALIGNED(ea));
    assert(IS_DMA_SIZE(this_size*sizeof(T)));
    (pf->f_dtl_entry_add)(entries, this_size, ALF_DATA_FLOAT, ea);
    size -= this_size;
    ea   += this_size * sizeof(T);
  }
}



/// Helper functor, converts void buffer pointer to appropriate type.
///
/// The 'off' parameter is a byte offset, while 'size' is in elements.
/// This is necessary because the first is calculated from the amount
/// data previously allocated, which may or may not have the same data
/// type.  Conversely, the second parameter refers to the amount of
/// data for the current segment and it is therefore easier to use
/// pointer arithmetic since the type is known.
///
template <typename T>
struct To_ptr
{
  static T offset(void* data, size_t off, size_t) 
    { return (T)((size_t)data + off); }
};

template <typename T>
struct To_ptr<std::pair<T*, T*> >
{
  static std::pair<T*, T*> offset(void* data, size_t off, size_t size)
    { return std::pair<T*, T*>((T*)((size_t)data + off), 
                               (T*)((size_t)data + off) + size); }
};


/// Converts a size in number of elements (or index value) into an offset 
/// based on the type referenced by the pointer.
///
template <typename T>
struct Byte_offset;

template <typename T>
struct Byte_offset<T*>
{
  static size_t index(size_t size) { return sizeof(T) * size; }
};

template <typename T>
struct Byte_offset<std::pair<T*, T*> >
{
  static size_t index(size_t size) { return 2 * sizeof(T) * size; }
};


void
stream_buffer_size(
  Uk_stream&    stream,
  unsigned int  iter, 
  unsigned int  iter_count,
  unsigned int& num_lines,
  unsigned int& line_size,
  int&          offset,
  char          ptype)
{
  unsigned int chunk_idx;
  unsigned int chunk_idx0;
  unsigned int chunk_idx1;

  if (stream.dim == 3)
  {
    // Currently there is a restriction on the third dimension of a view
    // being limited to a whole distribution.  To handle this, the third
    // dimension is folded into the second (similar to the way a dense
    // 2-D view can be recast as a 1-D view)

    assert(stream.num_chunks2 == 1);

    chunk_idx  = stream.chunk_offset + iter;

    chunk_idx0 = chunk_idx / stream.num_chunks1;
    chunk_idx1 = chunk_idx % stream.num_chunks1;
    offset     = 
        chunk_idx0 * stream.chunk_size0 * stream.stride0 * sizeof(float)
      + chunk_idx1 * stream.chunk_size1 * stream.stride1 * sizeof(float);

    num_lines = stream.chunk_size0;
    line_size = stream.chunk_size1 * stream.chunk_size2;
  }
  else
  {
    chunk_idx  = stream.chunk_offset + iter;

    chunk_idx0 = chunk_idx / stream.num_chunks1;
    chunk_idx1 = chunk_idx % stream.num_chunks1;
    offset     = 
        chunk_idx0 * stream.chunk_size0 * stream.stride0 * sizeof(float)
      + chunk_idx1 * stream.chunk_size1 * stream.stride1 * sizeof(float);

    num_lines = stream.chunk_size0;
    line_size = stream.chunk_size1;
  }

  // Handle last chunk in row/column (if odd sized)
  if (chunk_idx0 == stream.num_chunks0-1 && stream.chunk_size0_extra)
    num_lines = stream.chunk_size0_extra;
  if (chunk_idx1 == stream.num_chunks1-1 && stream.chunk_size1_extra)
    line_size = stream.chunk_size1_extra;

  offset    -= stream.align_shift * sizeof(float);
  line_size += stream.align_shift;

  // Handle overlap.
  if (stream.leading_overlap0 && 
      (chunk_idx0 != 0 || !stream.skip_first_overlap0))
  {
    num_lines += stream.leading_overlap0;
    offset    -= stream.leading_overlap0 * stream.stride0 * sizeof(float);
  }
  if (stream.trailing_overlap0 &&
      (chunk_idx0 != stream.num_chunks0-1 || !stream.skip_last_overlap0))
  {
    num_lines += stream.trailing_overlap0;
  }

  if (stream.leading_overlap1 &&
      (chunk_idx1 != 0 || !stream.skip_first_overlap1))
  {
    unsigned int leading1 =
      INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
    line_size += leading1;
    offset    -= leading1 * sizeof(float);
  }
  if (stream.trailing_overlap1 &&
      (chunk_idx1 != stream.num_chunks1-1 || !stream.skip_last_overlap1))
  {
    unsigned int trailing1 =
      INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.trailing_overlap1);
    line_size += trailing1;
  }

  line_size = INCREASE_TO_DMA_SIZE_IN_FLOATS(line_size);

#if DEBUG_ALF_BASE
  printf("add_stream: type: %c  chunk: %d (%d/%d, %d/%d)  size: %d/%d x %d/%d  stride: %d, %d\n",
    ptype, chunk_idx,
    chunk_idx0, stream.num_chunks0,
    chunk_idx1, stream.num_chunks1, 
    stream.chunk_size0, num_lines, stream.chunk_size1, line_size,
    stream.stride0, stream.stride1);
#endif

}


template <typename PtrT>
void
add_stream(
  Plugin_functions* pf,
  void*             entries,
  Uk_stream&        stream,
  unsigned int      iter, 
  unsigned int      iter_count)
{
  alf_data_addr64_t ea;
  int offset;
  unsigned int num_lines;
  unsigned int line_size;

  char ptype =
    Type_equal<PtrT, float*>::value ? 'S' :
      Type_equal<PtrT, std::complex<float>*>::value ? 'C' :
        Type_equal<PtrT, std::pair<float*, float*> >::value ? 'Z' :
          Type_equal<PtrT, unsigned int*>::value ? 'I' :
            '?';

  stream_buffer_size(stream, iter, iter_count, num_lines, line_size, offset,
		     ptype);

  if (Type_equal<PtrT, float*>::value)
  {
    ea = stream.addr + offset;

    for (int i=0; i<num_lines; ++i)
    {
      alf_data_addr64_t eax = ea + i*stream.stride0 * sizeof(float);
      add_entry<float>(pf, entries, line_size, eax);
    }
  }
  else if (Type_equal<PtrT, std::complex<float>*>::value)
  {
    ea = stream.addr + 2 * offset;

    for (int i=0; i<num_lines; ++i)
      add_entry<float>(pf, entries, 2*line_size,
		       ea + 2*i*stream.stride0*sizeof(float));
  }
  else if (Type_equal<PtrT, std::pair<float*, float*> >::value)
  {
    ea = stream.addr + offset;

    for (int i=0; i<num_lines; ++i)
      add_entry<float>(pf, entries, line_size,
		       ea + i*stream.stride0*sizeof(float));

    ea = stream.addr_split + offset;

    for (int i=0; i<num_lines; ++i)
      add_entry<float>(pf, entries, line_size,
		       ea + i*stream.stride0*sizeof(float));
  }
  else if (Type_equal<PtrT, unsigned int*>::value)
  {
    ea = stream.addr + offset;

    for (int i=0; i<num_lines; ++i)
    {
      alf_data_addr64_t eax = ea + i*stream.stride0 * sizeof(unsigned int);
      add_entry<unsigned int>(pf, entries, line_size, eax);
    }
  }
  else { assert(0); }
}



void
set_chunk_info(
  Uk_stream& stream,
  Pinfo&     pinfo,
  int        iter)
{
#if 0
  {                                                                     
    register volatile vector unsigned int get_r1 asm("1");
    unsigned int stack_pointer   = spu_extract(get_r1, 0);
    unsigned int p_stack_pointer = *(unsigned int*)stack_pointer;
    unsigned int return_addr     = *(unsigned int*)(p_stack_pointer+16);
    printf("STACK(set_chunk_info): cur %05x  pre %05x  ra %05x\n", 
           stack_pointer, p_stack_pointer, return_addr);
  }
#endif

  unsigned int chunk_idx  = stream.chunk_offset + iter;

  pinfo.dim = stream.dim;
  if (stream.dim == 1)
  {
    pinfo.dim         = 1;
    pinfo.g_offset[0] = chunk_idx * stream.chunk_size1;

    if (chunk_idx == stream.num_chunks1-1 && stream.chunk_size1_extra)
      pinfo.l_size[0] = stream.chunk_size1_extra;
    else
      pinfo.l_size[0] = stream.chunk_size1;

    pinfo.l_total_size = pinfo.l_size[0];

    pinfo.l_stride[0] = 1;

    pinfo.l_offset[0] = stream.align_shift;

    if (stream.leading_overlap1 &&
	(chunk_idx != 0 || !stream.skip_first_overlap1))
    {
      pinfo.o_leading[0] = stream.leading_overlap1;
      pinfo.l_offset[0] +=
	INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
    }
    else
    {
      pinfo.o_leading[0] = 0;
    }

    if (stream.trailing_overlap1 &&
	(chunk_idx != stream.num_chunks1-1 || !stream.skip_last_overlap1))
      pinfo.o_trailing[0] = stream.trailing_overlap1;
    else
      pinfo.o_trailing[0] = 0;

  }
  else if (stream.dim == 2)
  {
    unsigned int chunk_idx0 = chunk_idx / stream.num_chunks1;
    unsigned int chunk_idx1 = chunk_idx % stream.num_chunks1;
    pinfo.g_offset[0]   = chunk_idx0 * stream.chunk_size0;
    pinfo.g_offset[1]   = chunk_idx1 * stream.chunk_size1;

    if (chunk_idx0 == stream.num_chunks0-1 && stream.chunk_size0_extra)
      pinfo.l_size[0] = stream.chunk_size0_extra;
    else
      pinfo.l_size[0] = stream.chunk_size0;

    if (chunk_idx1 == stream.num_chunks1-1 && stream.chunk_size1_extra)
      pinfo.l_size[1] = stream.chunk_size1_extra;
    else
      pinfo.l_size[1] = stream.chunk_size1;

    pinfo.l_total_size = pinfo.l_size[0] * pinfo.l_size[1];

    pinfo.l_stride[0]   = pinfo.l_size[1];
    pinfo.l_stride[1]   = 1;
    pinfo.o_leading[0]  = stream.leading_overlap0;
    pinfo.o_leading[1]  = stream.leading_overlap1;
    pinfo.o_trailing[0] = stream.trailing_overlap0;
    pinfo.o_trailing[1] = stream.trailing_overlap1;

    if (stream.leading_overlap0 &&
	(chunk_idx0 != 0 || !stream.skip_first_overlap0))
    {
      pinfo.o_leading[0] = stream.leading_overlap0;
      pinfo.l_offset[0]  = stream.leading_overlap0;
    }
    else
    {
      pinfo.o_leading[0] = 0;
      pinfo.l_offset[0]  = 0;
    }

    if (stream.trailing_overlap0 &&
	(chunk_idx0 != stream.num_chunks0-1 || !stream.skip_last_overlap0))
      pinfo.o_trailing[0] = stream.trailing_overlap0;
    else
      pinfo.o_trailing[0] = 0;

    if (stream.leading_overlap1 &&
	(chunk_idx1 != 0 || !stream.skip_first_overlap1))
    {
      pinfo.o_leading[1] = stream.leading_overlap1;
      pinfo.l_offset[1] =
	stream.align_shift + 
	INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
      pinfo.l_stride[0] += pinfo.l_offset[1];
    }
    else
    {
      pinfo.o_leading[1] = 0;
      pinfo.l_offset[1] = stream.align_shift;
    }

    if (stream.trailing_overlap1 &&
	(chunk_idx1 != stream.num_chunks1-1 || !stream.skip_last_overlap1))
    {
      pinfo.o_trailing[1] = stream.trailing_overlap1;
      pinfo.l_stride[0]  +=
	INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
    }
    else
      pinfo.o_trailing[1] = 0;
  }
  else if (stream.dim == 3)
  {
    unsigned int chunk_idx0 = chunk_idx / stream.num_chunks1;
    unsigned int chunk_idx1 = chunk_idx % stream.num_chunks1;
    pinfo.g_offset[0]   = chunk_idx0 * stream.chunk_size0;
    pinfo.g_offset[1]   = chunk_idx1 * stream.chunk_size1;
    pinfo.g_offset[2]   = 0;

    if (chunk_idx0 == stream.num_chunks0-1 && stream.chunk_size0_extra)
      pinfo.l_size[0] = stream.chunk_size0_extra;
    else
      pinfo.l_size[0] = stream.chunk_size0;

    if (chunk_idx1 == stream.num_chunks1-1 && stream.chunk_size1_extra)
      pinfo.l_size[1] = stream.chunk_size1_extra;
    else
      pinfo.l_size[1] = stream.chunk_size1;

    assert(stream.num_chunks2 == 1);
    pinfo.l_size[2] = stream.chunk_size2;

    pinfo.l_total_size = pinfo.l_size[0] * pinfo.l_size[1] * pinfo.l_size[2];

    pinfo.l_stride[0]   = pinfo.l_size[1] * pinfo.l_size[2];
    pinfo.l_stride[1]   = 1;
    pinfo.l_stride[2]   = 1;
    pinfo.o_leading[0]  = stream.leading_overlap0;
    pinfo.o_leading[1]  = stream.leading_overlap1;
    pinfo.o_leading[2]  = stream.leading_overlap2;
    pinfo.o_trailing[0] = stream.trailing_overlap0;
    pinfo.o_trailing[1] = stream.trailing_overlap1;
    pinfo.o_trailing[2] = stream.trailing_overlap2;

    if (stream.leading_overlap0 &&
	(chunk_idx0 != 0 || !stream.skip_first_overlap0))
    {
      pinfo.o_leading[0] = stream.leading_overlap0;
      pinfo.l_offset[0]  = stream.leading_overlap0;
    }
    else
    {
      pinfo.o_leading[0] = 0;
      pinfo.l_offset[0]  = 0;
    }

    if (stream.trailing_overlap0 &&
        (chunk_idx0 != stream.num_chunks0-1 || !stream.skip_last_overlap0))
      pinfo.o_trailing[0] = stream.trailing_overlap0;
    else
      pinfo.o_trailing[0] = 0;

    if (stream.leading_overlap1 &&
	(chunk_idx1 != 0 || !stream.skip_first_overlap1))
    {
      pinfo.o_leading[1] = stream.leading_overlap1;
      pinfo.l_offset[1] =
	stream.align_shift + 
	INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
      pinfo.l_stride[0] += pinfo.l_offset[1];
    }
    else
    {
      pinfo.o_leading[1] = 0;
      pinfo.l_offset[1] = stream.align_shift;
    }

    if (stream.trailing_overlap1 &&
	(chunk_idx1 != stream.num_chunks1-1 || !stream.skip_last_overlap1))
    {
      pinfo.o_trailing[1] = stream.trailing_overlap1;
      pinfo.l_stride[0]  +=
	INCREASE_TO_DMA_SIZE_IN_FLOATS(stream.leading_overlap1);
    }
    else
      pinfo.o_trailing[1] = 0;

    // overlap for the third dimension is not supported
    assert(stream.leading_overlap2 == 0);
    assert(stream.trailing_overlap2 == 0);
  }
  else
    assert(0);
}


/***********************************************************************
  Z = f(X)
***********************************************************************/

template <typename      KernelT,
	  unsigned int PreArgc = KernelT::pre_argc,
	  unsigned int InArgc  = KernelT::in_argc,
	  unsigned int OutArgc = KernelT::out_argc>
struct Kernel_helper;

template <typename KernelT>
struct Kernel_helper<KernelT, 0, 0, 0>
{
  static void
  input(
    Plugin_functions* /*pf*/,
    param_type*  /*ukp*/,
    void*        /*entries*/,
    unsigned int /*iter*/, 
    unsigned int /*iter_count*/)
  {
  }

  static void
  kernel(
    KernelT&     ukobj,
    param_type*  /*ukp*/,
    void*        /*inout*/,
    unsigned int /*iter*/, 
    unsigned int /*iter_count*/)
  {
    ukobj.compute();
  }
};



template <typename KernelT>
struct Kernel_helper<KernelT, 0, 1, 1>
{
  static void
  input(
    Plugin_functions* pf,
    param_type*       ukp,
    void*             entries,
    unsigned int      iter, 
    unsigned int      iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_out;
    set_chunk_info(ukp->out_stream[0], p_out, iter);
    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);

    // Transfer input A.
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, off1);

    add_stream<in0_type>(pf, entries, ukp->in_stream[0], iter, iter_count);

    (pf->f_dtl_end)(entries);
  }

  static void
  kernel(
    KernelT&     ukobj,
    param_type*  ukp,
    void*        inout,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_in, p_out;

    set_chunk_info(ukp->in_stream[0], p_in,   iter);
    set_chunk_info(ukp->out_stream[0], p_out, iter);

    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);

    ukobj.compute(
      To_ptr<in0_type >::offset(inout, off1, p_in.l_total_size),
      To_ptr<out0_type>::offset(inout,    0, p_out.l_total_size),
      p_in, p_out);
  }
};

template <typename KernelT>
struct Kernel_helper<KernelT, 0, 2, 1>
{
  static void
  input(
    Plugin_functions* pf,
    param_type*  ukp,
    void*        entries,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_out;
    set_chunk_info(ukp->out_stream[0], p_out, iter);
    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);

    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, off1);

    add_stream<in0_type>(pf, entries, ukp->in_stream[0], iter, iter_count);
    add_stream<in1_type>(pf, entries, ukp->in_stream[1], iter, iter_count);

    (pf->f_dtl_end)(entries);
  }

  static void
  kernel(
    KernelT&     ukobj,
    param_type*  ukp,
    void*        inout,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_in0, p_in1, p_out;

    set_chunk_info(ukp->in_stream[0],  p_in0, iter);
    set_chunk_info(ukp->in_stream[1],  p_in1, iter);
    set_chunk_info(ukp->out_stream[0], p_out, iter);

    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);
    size_t off2 = Byte_offset<in0_type>::index(p_in0.l_total_size) + off1;

    ukobj.compute(
      To_ptr<in0_type >::offset(inout, off1, p_in0.l_total_size),
      To_ptr<in1_type >::offset(inout, off2, p_in1.l_total_size),
      To_ptr<out0_type>::offset(inout, 0,    p_out.l_total_size),
      p_in0, p_in1, p_out);
  }
};



template <typename KernelT>
struct Kernel_helper<KernelT, 1, 1, 1>
{
  static void
  input(
    Plugin_functions* pf,
    param_type*  ukp,
    void*        entries,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::out0_type out0_type;

    if (iter < ukp->pre_chunks)
    {
      (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, 0);
      add_stream<in0_type>(pf, entries, ukp->in_stream[0], iter, iter_count);
    }
    else
    {
      Pinfo p_out;
      set_chunk_info(ukp->out_stream[0], p_out, iter);
      size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);

      (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, off1);
      add_stream<in1_type>(pf, entries, ukp->in_stream[1],
			   iter - ukp->pre_chunks, 
			   iter_count - ukp->pre_chunks);
    }

    (pf->f_dtl_end)(entries);
  }

  static void
  kernel(
    KernelT&     ukobj,
    param_type*  ukp,
    void*        inout,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_in0, p_in1, p_out;

    if (iter < ukp->pre_chunks)
    {
      set_chunk_info(ukp->in_stream[0],  p_in0, iter);
      ukobj.pre_compute(
	To_ptr<in0_type >::offset(inout,  0, p_in0.l_total_size),
	p_in0);
    }
    else
    {
      // the iteration count must be adjusted to account for the
      // one used above
      set_chunk_info(ukp->in_stream[1],  p_in1, iter - 1);
      set_chunk_info(ukp->out_stream[0], p_out, iter - 1);
      size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);
      ukobj.compute(
	To_ptr<in1_type >::offset(inout, off1, p_in1.l_total_size),
	To_ptr<out0_type>::offset(inout, 0,    p_out.l_total_size),
	p_in1, p_out);
    }
  }
};
  

template <typename KernelT>
struct Kernel_helper<KernelT, 0, 3, 1>
{
  static void
  input(
    Plugin_functions* pf,
    param_type*  ukp,
    void*        entries,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::in2_type  in2_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_out;
    set_chunk_info(ukp->out_stream[0], p_out, iter);
    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);

    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_IN, off1);

    add_stream<in0_type>(pf, entries, ukp->in_stream[0], iter, iter_count);
    add_stream<in1_type>(pf, entries, ukp->in_stream[1], iter, iter_count);
    add_stream<in2_type>(pf, entries, ukp->in_stream[2], iter, iter_count);

    (pf->f_dtl_end)(entries);
  }

  static void
  kernel(
    KernelT&     ukobj,
    param_type*  ukp,
    void*        inout,
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::in0_type  in0_type;
    typedef typename KernelT::in1_type  in1_type;
    typedef typename KernelT::in2_type  in2_type;
    typedef typename KernelT::out0_type out0_type;

    Pinfo p_in0, p_in1, p_in2, p_out;

    set_chunk_info(ukp->in_stream[0],  p_in0, iter);
    set_chunk_info(ukp->in_stream[1],  p_in1, iter);
    set_chunk_info(ukp->in_stream[2],  p_in2, iter);
    set_chunk_info(ukp->out_stream[0], p_out, iter);

    // Pointers must be extracted from knowledge of the stream sizes as ALF
    // transfers all the input data into one contiguous space.

    size_t off1 = Byte_offset<out0_type>::index(p_out.l_total_size);
    size_t off2 = Byte_offset<in0_type>::index(p_in0.l_total_size) + off1;
    size_t off3 = Byte_offset<in1_type>::index(p_in1.l_total_size) + off2;

    // The To_ptr<> struct calculates the correct offset for a given
    // pointer type (scalar, interleaved complex or split complex).  The 
    // first size passes refers to the previous data segments.  The second 
    // size pertains to the current segment and is only needed to calculate 
    // offsets in the case of split complex.

    ukobj.compute(
      To_ptr<in0_type >::offset(inout, off1, p_in0.l_total_size),
      To_ptr<in1_type >::offset(inout, off2, p_in1.l_total_size),
      To_ptr<in2_type >::offset(inout, off3, p_in2.l_total_size),
      To_ptr<out0_type>::offset(inout, 0,    p_out.l_total_size),
      p_in0, p_in1, p_in2, p_out);
  }
};



template <typename     KernelT,
	  unsigned int OutArgc = KernelT::out_argc>
struct Output_helper
{
  static void
  output(
    Plugin_functions* pf,
    param_type*  ukp,
    void*        entries, 
    unsigned int iter, 
    unsigned int iter_count)
  {}
};

template <typename KernelT>
struct Output_helper<KernelT, 1>
{
  static void
  output(
    Plugin_functions* pf,
    param_type*  ukp,
    void*        entries, 
    unsigned int iter, 
    unsigned int iter_count)
  {
    typedef typename KernelT::out0_type out0_type;

    if (iter < ukp->pre_chunks)
      return;
    else
    {
      iter       -= ukp->pre_chunks;
      iter_count -= ukp->pre_chunks;
    }
    
    // Transfer output Z.
    (pf->f_dtl_begin)(entries, ALF_BUF_OVL_OUT, 0);
    
    add_stream<out0_type>(pf, entries, ukp->out_stream[0], iter, iter_count);
    
    (pf->f_dtl_end)(entries);
  }
};



extern "C"
int
input(
  Plugin_functions* pf,
  void*        p_context,
  void*        params, 
  void*        entries, 
  unsigned int iter, 
  unsigned int iter_count)
{
  param_type*  ukp  = (param_type*)params;

  Kernel_helper<kernel_type>::input(pf, ukp, entries, iter, iter_count);
  return 0;
}



extern "C"
int
output(
  Plugin_functions* pf,
  void*        p_context,
  void*        p_params, 
  void*        entries, 
  unsigned int iter, 
  unsigned int iter_count)
{
  param_type* ukp = (param_type*)p_params;

  Output_helper<kernel_type>::output(pf, ukp, entries, iter, iter_count);

  return 0;
}



extern "C"
int
kernel(
  Plugin_functions* pf,
  void*             p_context,
  void*             params,
  void*             inout,
  unsigned int      iter,
  unsigned int      iter_count)
{
  param_type*  ukp  = (param_type *)params;

  if (iter == 0)
  {
    ukobj.init_rank(ukp->rank, ukp->nspe);
    ukobj.init(ukp->kernel_params);
  }

  Kernel_helper<kernel_type>::kernel(ukobj, ukp, inout, iter, iter_count);

  if (iter == iter_count-1)
    ukobj.fini();

  return 0;
}

#undef DMA_ALIGNMENT
#undef DMA_SIZE_QUANTUM
#undef DMA_SIZE_IN_FLOATS
#undef DMA_ALIGNMENT_OF
#undef IS_DMA_ALIGNED
#undef IS_DMA_SIZE
#undef INCREASE_TO_DMA_SIZE_IN_FLOATS
