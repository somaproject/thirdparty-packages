#! /bin/sh

#########################################################################
# mcoe-setup.sh -- setup script to configure SourceryVSIPL++ for use on
#                  Mercury systems
#
# (27 Feb 05) Jules Bergmann, CodeSourcery, Inc.
#########################################################################

#########################################################################
# Instructions:
#  - Modify flags below to control where and how VSIPL++ is built.
#  - Run mcoe-setup.sh
#  - Run 'make'
#  - Run 'make install'
#
#    Flags can either be uncommented and edited below, or placed in a
#    separate file that sources this one.
#
# Flags:
#
#   dir="."			# Sourcery VSIPL++ source directory.
#   comm="ser"			# set to (ser)ial or (par)allel.
#   fmt="split"			# set to (inter)leaved or (split).
#   opt="y"			# (y) for optimized flags, (n) for debug flags.
#   compiler="def"		# (def)ault, (GHS) - GreenHills, (GNU) - GNU
#   simd_loop_fusion="y"	# (y) for SIMD loop fusion, (n) for not.
#   builtin_simd="y"		# (y) for builtin SIMD routines, (n) for not.
#   pflags="-t ppc7447"		# processor architecture
#   fft="sal,builtin"		# FFT backend(s)
#   testlevel="0"		# Test level
#   prefix="/opt/vsipl++"	# Installation prefix.
#
# Notes:
#  - For compiler, if value is "def", GHS is assumed, but -compiler
#    flag is not given.
#
#########################################################################

# 'dir' is the directory containing SourceryVSIPL++
if test "x$dir" = x; then
  dir="."
fi

if test "x$comm" = x; then
  comm="ser"			# set to (ser)ial or (par)allel.
fi

if test "x$fmt" = x; then
  fmt="split"			# set to (inter)leaved or (split).
fi

if test "x$opt" = x; then
  opt="y"			# (y) for optimized flags, (n) for debug flags.
fi

if test "x$compiler" = x; then
  compiler="def"		# (def), (GHS) for GreenHills, (GNU) for GNU
fi

if test "x$simd_loop_fusion" = x; then
  simd_loop_fusion="y"		# (y) for SIMD loop fusion, (n) for not.
fi

if test "x$builtin_simd" = x; then
  builtin_simd="y"		# (y) for builtin SIMD, (n) for not.
fi

if test "x$sal" = x; then
  sal="y"			# (y) to use SAL, (n) to not.
fi

if test "x$exceptions" = x; then
  exceptions="n"		# (y) for exceptions, (n) for not.
fi

if test "x$pflags" = x; then
  pflags="-t ppc7447"		# processor architecture
fi


# FFT backend.  This controls which backend or backends are used for
# FFTs.
#
# Possible Values:
#   sal,builtin	- Use SAL and builtin Sourcery VSIPL++ FFTW3 (recommended).
#   sal,fftw3	- Use SAL and pre-built FFTW3.
#   sal		- Use SAL only (only power-of-2 sizes are supported).

if test "x$fft" = x; then
  fft="sal,builtin"	# FFT backend.
fi


# Test level.  This controls how hard Sourcery VSIPL++'s test-suite
# tries to test the system.
#
# Values:
#   0 - low-level (avoids long-running and long-compiling tests).
#   1 - default
#   2 - high-level (enables additional long-running tests).

if test "x$testlevel" = x; then
  testlevel="0"		# Test level
fi



# Set 'prefix' the directory where SourceryVSIPL++ Should be installed
if test "x$prefix" = x; then
  prefix="/opt/sourcery-vsipl++"
fi



#########################################################################
cfgflags=""

if test $comm = "pas"; then
  cfg_flags="$cfg_flags --enable-pas"
elif test $comm = "mpi" -o $comm = "par"; then
  cfg_flags="$cfg_flags --enable-mpi=mpipro"
else
  cfg_flags="$cfg_flags --disable-mpi"
fi

# If compiler = "def", assume GHS, but do not set '-compiler GHS' flag.
if test $compiler = "GHS" -o $compiler = "def"; then
  if test $compiler = "GHS"; then
    toolset_flag="-compiler GHS"
  fi
  cxxflags="$pflags $toolset_flag --no_implicit_include"
 
  opt_flags="-Ospeed -Onotailrecursion --max_inlining"
  opt_flags="$opt_flags -DNDEBUG --diag_suppress 175,177,550"
  dbg_flags="-g"

  ex_off_flags="--no_exceptions"
  ex_on_flags="--exceptions"

  fftw3_cflags="-Ospeed $pflags $toolset_flag"
else
  toolset_flag="-compiler GNU"
  cxxflags="$pflags $toolset_flag"

  opt_flags="-Otime -DNDEBUG -w"
  dbg_flags="-g"

  ex_off_flags="-fno-exceptions"
  ex_o_flags=""				# exceptions enabled by default.

  fftw3_cflags="-Otime $pflags $toolset_flag"
fi

if test $opt = "y"; then
  cxxflags="$cxxflags $opt_flags"
else
  cxxflags="$cxxflags $dbg_flags"
fi

if test $builtin_simd = "y"; then
  cfg_flags="$cfg_flags --with-builtin-simd-routines=generic"
fi

if test $simd_loop_fusion = "y"; then
  cfg_flags="$cfg_flags --enable-simd-loop-fusion"
else
  cfg_flags="$cfg_flags --disable-simd-loop-fusion"
fi

if test $sal = "y"; then
  cfg_flags="$cfg_flags --with-sal"
fi

if test $exceptions = "n"; then
  cxxflags="$cxxflags $ex_off_flags"
  cfg_flags="$cfg_flags --disable-exceptions"
else
  cxxflags="$cxxflags $ex_on_flags"
  cfg_flags="$cfg_flags --enable-exceptions"
fi

if test "x$extra_args" != "x"; then
  cfg_flags="$cfg_flags $extra_args"
fi


# select timer
if test "x$timer" = "x"; then
  # timer=realtime
  timer=mcoe_tmr
fi

#########################################################################
# export environment variables

CC=ccmc
CFLAGS="$toolset_flag"
CXX=ccmc++
CXXFLAGS=$cxxflags
AR=armc
AR_FLAGS="$toolset_flag cr"	# armc doesn't support 'u'pdate
LDFLAGS="$pflags $toolset_flag"

export CC
export CXX
export CFLAGS
export CXXFLAGS
export AR
export AR_FLAGS
export LDFLAGS


#########################################################################
# run configure

echo "$dir/configure"
$dir/configure						\
	--prefix=$prefix				\
	--host=powerpc					\
	--enable-fft=$fft				\
	--with-fftw3-cflags="$fftw3_cflags"		\
	--with-fftw3-cfg-opts="--with-our-malloc16"	\
	--with-complex=$fmt				\
	--with-lapack=no				\
	$cfg_flags					\
	--with-test-level=$testlevel			\
	--enable-timer=$timer
