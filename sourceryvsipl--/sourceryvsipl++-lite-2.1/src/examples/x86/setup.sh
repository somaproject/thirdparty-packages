#! /bin/sh

#########################################################################
# examples/x86/setup.sh -- setup script to configure Sourcery VSIPL++
#			   for use on x86 systems
#
# (29 Oct 08) Jules Bergmann, CodeSourcery, Inc.
#########################################################################

#########################################################################
# Instructions:
#  - Modify variables below to control where and how VSIPL++ is built.
#  - Run setup.sh
#  - Run 'make'
#  - Run 'make install'
#
#    Variables can either be uncommented and edited below, or placed in a
#    separate file that sources this one.
#
# Variables:
#  - src_dir	-- source directory
#  - ipp_dir	-- IPP directory
#  - mkl_dir	-- MKL directory
#  - svpp_prefix-- VSIPL++ installation prefix (VSIPL++ will be installed
#                  here).
#  - comm="ser"	-- set to (ser)ial or (par)allel.
#########################################################################

if test "x$src_dir" = x; then
  src_dir=/opt/sourceryvsipl++-2.0/src
fi

if test "x$comm" = x; then
  comm="ser"			# set to (ser)ial or (par)allel.
fi

if test "x$ipp_dir" = x; then
  ipp_dir=/opt/intel/ipp5/ipp/5.0/em64t
fi

if test "x$mkl_dir" = x; then
  mkl_dir=/opt/intel/mkl721
fi

if test "x$svpp_prefix" = x; then
  svpp_prefix=/opt/sourceryvsipl++-2.0
fi


width="-m64"
opt="-g -O2 -funswitch-loops -fgcse-after-reload -DNDEBUG --param max-inline-insns-single=2000 --param large-function-insns=6000 --param large-function-growth=800 --param inline-unit-growth=300"

export CFLAGS="$width $opt"
export CXXFLAGS="$width $opt"
export LDFLAGS="$width"
export CC=gcc
export CXX=g++
export LD=ld

if test "$comm" == "ser"; then
  args="$args --disable-parallel"
else
  args="$args --enable-parallel"
fi

$src_dir/configure							\
	$args								\
	--with-lapack=mkl						\
	--with-ipp-prefix=$ipp_dir					\
	--with-mkl-prefix=$mkl_dir					\
	--enable-fft=ipp						\
	--with-builtin-simd-routines=generic				\
	--enable-simd-loop-fusion					\
	--with-test-level=1						\
	--prefix=$svpp_prefix						\
	--enable-timer=x86_64_tsc
