dnl Copyright (c) 2007 by CodeSourcery, Inc.  All rights reserved.
dnl
dnl File:   cuda.m4
dnl Author: Stefan Seefeld
dnl Date:   2009-02-12
dnl
dnl Contents: CUDA configuration for Sourcery VSIPL++
dnl

AC_DEFUN([SVXX_CHECK_CUDA],
[
  # There are two cuda libraries supported: CUBLAS and CUFFT.  They are
  # configured using --with-cuda and --enable-fft=cuda.  Use of the FFT
  # library requires use of the --with-cuda option.


  # CUBLAS
  if test "$with_cuda" != "no"; then

    CUDA_CPPFLAGS="-I/usr/local/cuda/include"
    save_CPPFLAGS="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $CUDA_CPPFLAGS"
  
    # Find header file
    vsipl_cublas_h_name="not found"
    AC_CHECK_HEADER([cublas.h], [vsipl_cublas_h_name='<cublas.h>'],, [// no prerequisites])
    if test "$vsipl_cublas_h_name" = "not found"; then
      AC_MSG_ERROR([CUDA enabled, but no cublas.h detected])
      CPPFLAGS="$save_CPPFLAGS"
    else

      # Find the library.
      save_LDFLAGS="$LDFLAGS"
      LDFLAGS="$LDFLAGS -L/usr/local/cuda/lib"
      AC_SEARCH_LIBS(cublasGetError, [cublas],
	             [cublas_lib=$ac_lib
		      cuda_found="yes"],
 		     [cuda_found="no"])

      if test "$cuda_found" = "no"; then
        AC_MSG_ERROR([CUDA BLAS library not found])
        CPPFLAGS=$save_CPPFLAGS
        LDFLAGS=$save_LDFLAGS
      else

        # Declare it as found.	
        AC_SUBST(VSIP_IMPL_HAVE_CUDA, 1)
        if test "$neutral_acconfig" = 'y'; then
          CPPFLAGS="$CPPFLAGS -DVSIP_IMPL_HAVE_CUDA=1"
        else
          AC_DEFINE_UNQUOTED(VSIP_IMPL_HAVE_CUDA, 1,
            [Define to set whether or not to use NVIDIA's CUDA libraries.])
        fi
      fi
    fi

    # CUFFT
    if test $enable_cuda_fft != "no"; then

      # Find header file
      vsipl_cufft_h_name="not found"
      AC_CHECK_HEADER([cufft.h], [vsipl_cufft_h_name='<cufft.h>'],, [// no prerequisites])
      if test "$vsipl_cufft_h_name" = "not found"; then
        AC_MSG_ERROR([CUDA FFT enabled, but no cufft.h detected])
      else

        # Find the library file
        AC_SEARCH_LIBS(cufftPlan1d, [cufft],
                       [cufft_lib=$ac_lib
		        cuda_fft_found="yes"],
 		       [cuda_fft_found="no"])

        if test "$cuda_fft_found" = "no"; then
          AC_MSG_ERROR([CUDA FFT library not found])
        else

          # Declare it as found.
	  provide_fft_float=1
          AC_SUBST(VSIP_IMPL_CUDA_FFT, 1)
          if test "$neutral_acconfig" = 'y'; then
            CPPFLAGS="$CPPFLAGS -DVSIP_IMPL_CUDA_FFT=1"
          else
            AC_DEFINE_UNQUOTED(VSIP_IMPL_CUDA_FFT, 1,
              [Define to set whether or not to use NVIDIA's CUDA FFT library.])
          fi
        fi
      fi
    fi
  fi
])
