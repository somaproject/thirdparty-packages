#! /bin/sh

#########################################################################
# examples/cell/setup.sh -- setup script to configure Sourcery VSIPL++
#			    for use on Cell/B.E. systems
#
# (25 Aug 08) Jules Bergmann, CodeSourcery, Inc.
#########################################################################

#########################################################################
# Instructions:
#  - Modify flags below to control where and how VSIPL++ is built.
#  - Run setup.sh
#  - Run 'make'
#  - Run 'make install'
#
# Variables:
#  - src_dir	-- source directory
#  - sdk_dir	-- SDK install directory (default is for SDK 3.0)
#  - cml_dir	-- CML install directory
#  - svpp_prefix-- VSIPL++ installation prefix (VSIPL++ will be installed
#                  here).
#########################################################################

src_dir=../sourceryvsipl++-20080825
sdk_dir=/opt/ibm/cell-sdk/prototype
cml_dir=/scratch/jules/opt/cml
svpp_prefix=/scratch/jules/build-test/install

bits="-m32"
opt="-mcpu=cell -maltivec -g -O2 -funswitch-loops -fgcse-after-reload -DNDEBUG --param max-inline-insns-single=2000 --param large-function-insns=6000 --param large-function-growth=800 --param inline-unit-growth=300"

export CFLAGS="$bits $opt"
export CXXFLAGS="$bits $opt"
export LDFLAGS="$bits"
export CC=ppu-gcc
export CXX=ppu-g++
export LD=ppu-ld

$src_dir/configure							\
	--with-cbe-sdk							\
	--with-cbe-sdk-prefix=$sdk_dir					\
	--disable-fft-long-double					\
	--disable-parallel						\
	--with-lapack=atlas						\
	--with-atlas-include=/usr/include/atlas				\
	--with-atlas-libdir=/usr/lib/altivec				\
	--enable-fft=cbe_sdk,fftw3					\
	--with-builtin-simd-routines=generic				\
	--with-complex=split						\
	--with-test-level=1						\
	--prefix=$svpp_prefix						\
	--with-cml-prefix=$cml_dir 			 		\
	--enable-timer=power_tb
