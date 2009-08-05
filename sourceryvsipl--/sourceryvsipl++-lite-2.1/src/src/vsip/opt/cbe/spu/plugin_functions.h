#include <alf_accel.h>

typedef int printf_func(const char*, ...);
typedef void dtl_begin_func(void*, ALF_BUF_TYPE_T, unsigned int);
typedef void dtl_entry_add_func(void*, unsigned int, ALF_DATA_TYPE_T, alf_data_addr64_t);
typedef void dtl_end_func(void*);

typedef struct
{
  printf_func*		f_printf;
  dtl_begin_func*	f_dtl_begin;
  dtl_entry_add_func*	f_dtl_entry_add;
  dtl_end_func*		f_dtl_end;
} Plugin_functions;

typedef int kernel_function(
  Plugin_functions* p_functions,
  void*        p_context,
  void*        p_params,
  void*        inout,
  unsigned int iter,
  unsigned int n);

typedef int input_function(
  Plugin_functions* pf,
  void*        context,
  void*        params,
  void*        entries,
  unsigned int iter,
  unsigned int iter_max);

typedef int output_function(
  Plugin_functions* pf,
  void*        context,
  void*        params,
  void*        entries,
  unsigned int iter,
  unsigned int iter_max);

