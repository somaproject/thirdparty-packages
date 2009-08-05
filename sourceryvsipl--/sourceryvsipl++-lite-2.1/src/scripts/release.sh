#! /bin/sh

########################################################################
#
# File:   release.sh
# Author: Jules Bergmann
# Date:   2005-01-11
#
# Contents:
#   Script to build Sourcery VSIPL++ release.
#
########################################################################

# SYNOPSIS
#   release.sh -w <command> -d <dir>
#
# OPTIONS
#   -w <command>
#      Specify action of script.  Valid choices are:
#         'src'  - build source package.
#         'bin'  - build binary packages.
#         'test' - test binary packages.
#         'all'  - combination of src, bin, and test (default).
#
#   -d <dir>
#      Specify source directory for 'package' command.
#      Path to command will be <dir>/scripts/package.
#
# DESCRIPTION
#   This script automates the building of source and binary packages
#   for Sourcery VSIPL++.  It uses the 'package' script to perform
#   the following steps:
#     1. Check out sources from SVN,
#     2. Build a source package,
#     3. Build binary packages from source package,
#     4. Test binary packages.

what="all"
dir=$HOME/csl/src/vpp/SVN-HEAD
cfgdir=default
cfgfile=default
pkgs="SerialBuiltin32 SerialBuiltin64 ParallelBuiltin32 ParallelBuiltin64 SerialIntel32 SerialIntel64 ParallelIntel32 ParallelIntel64"
svn_srcdir="svn_srcdir"
src_builddir="vpp-src-build"
test_srcdir="default"
distdir="vpp-dist"
debug="yes"
pkg_opts=""
version="1.4"
host=`hostname`

while getopts "w:c:d:p:C:t:D:T:sS:v:" arg; do
    case $arg in
	w)
	    what=$OPTARG
	    ;;
	c)
	    cfgfile=$OPTARG
	    ;;
	C)
	    svn_srcdir=$OPTARG
	    ;;
	d)
	    dir=$OPTARG
	    ;;
	D)
	    distdir=$OPTARG
	    ;;
	p)
	    pkgs=$OPTARG
	    ;;
	t)
	    test_srcdir=$OPTARG
	    ;;
	s)
	    pkg_opts="$pkg_opts --snapshot"
            version=`date +%Y%m%d`
	    ;;
	v)
	    version=$OPTARG
	    ;;
	\?)
            error "usage: release.sh [-v VERSION]"
	    ;;
    esac
done

srcdir="sourceryvsipl++-$version"
srcpkg="$srcdir.tar.bz2"
prefix="/opt/sourceryvsipl++-$version"
prefix_not_in_tarball="/opt"

pkg_opts="$pkg_opts --prefix=$prefix"
pkg_opts="$pkg_opts --prefix-not-in-tarball=$prefix_not_in_tarball"

package=$dir/scripts/package.py
if test "$cfgdir" = "default"; then
  cfgdir=$dir/scripts
fi

if test "x$cfgfile" = "x"; then
  echo "ERROR: Must specify a config file -c <file>"
  exit -1
fi

if test "$test_srcdir" = "default"; then
  test_srcdir=$srcdir
fi



########################################################################
# Step 0: Set environment variables for PATH and LD_LIBRARY_PATH
#         to known values.
########################################################################

TOOL_DIR=/usr/local/tools/vpp-1.0
GCCTOOL_DIR=/usr/local/tools/gcc-3.4.0
# GC_DIR=/opt/gc6.6/lib
# DOT_DIR=/usr/local/graphviz-2.6
GC_DIR=$HOME/build-cugel/gc6.6/lib
DOT_DIR=$HOME/local/`arch`

ipp_dir=/opt/intel/ipp
mkl_dir=/opt/intel/mkl
pas_dir=$TOOL_DIR/pas

PATH=$TOOL_DIR/sourceryg++/bin
PATH=$PATH:$TOOL_DIR/bin
PATH=$PATH:$GCCTOOL_DIR/bin
PATH=$PATH:/usr/bin
PATH=$PATH:/bin
PATH=$PATH:/usr/local/bin
PATH=$PATH:$DOT_DIR/bin
PATH=$PATH:/opt/renderx/xep
PATH=$PATH:$pas_dir/bin
if test `hostname` = "gannon.codesourcery.com"; then
  PATH=$PATH:/home/jules/local/sun4/bin
fi

LD_LIBRARY_PATH=$TOOL_DIR/sourceryg++/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TOOL_DIR/sourceryg++/lib64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TOOL_DIR/lib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TOOL_DIR/lib64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GC_DIR
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DOT_DIR/lib
# LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DOT_DIR/lib/graphviz
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ipp_dir/em64t/sharedlib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ipp_dir/em64t/sharedlib/linuxem64t
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ipp_dir/ia32_itanium/sharedlib
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ipp_dir/ia32_itanium/sharedlib/linux32
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mkl_dir/lib/em64t
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mkl_dir/lib/32
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$pas_dir/lib
if test `hostname` = "gannon.codesourcery.com"; then
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GCCTOOL_DIR/lib
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GCCTOOL_DIR/lib/sparcv9
fi

if test `hostname` = "gillette" -o `hostname` = "wesleysnipes"; then
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/tools/sdk/lib
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/jules/cell-sdk/sysroot/usr/lib
fi

export PATH
export LD_LIBRARY_PATH

if test "$debug" = "yes"; then
  echo "host            : $host"
  echo "which g++       : " `which g++`
  echo "which dot       : " `which dot`
  echo "which pkg-config: " `which pkg-config`
  echo "configure file  : $cfgfile"
  echo "svn_srcdir      : $svn_srcdir"

  save_IFS=$IFS; IFS=":"
  for d in $LD_LIBRARY_PATH
  do
    echo "LD_LIBRARY_PATH: $d"
  done
  IFS=$save_IFS
fi

# 1. Build/unpack source package.
if test "$what" = "src" -o "$what" = "all"; then
  echo "#####################################################################"
  echo "# build source package                                              #"
  echo "#####################################################################"
  if test -f "$srcpkg"; then
    echo "Source package ($srcpkg) found ... using it"
  else
    echo "Source package ($srcpkg) not found ... creating it"

    if test -d "$svn_srcdir"; then
      echo "ERROR: ($svn_srcdir already exists)"
      exit
    fi

    # 1a. Build source package (includes checkout and patch)
    echo "Build SDist (from $svn_srcdir)"
    $package build_sdist --verbose			\
	--srcdir=$svn_srcdir				\
	--configfile=$cfgfile				\
	--builddir=$src_builddir			\
	--configdir="$dir/scripts"			\
	$pkg_opts					\
	2>&1 > log-src-build

    # 1b. Untar source package.  Use this to build binary packages.
    mv $src_builddir/$srcpkg .
  fi
fi


# 2. Build binary packages.
if test "$what" = "bin" -o "$what" = "all"; then
  echo "#####################################################################"
  echo "# build binary packages                                             #"
  echo "#####################################################################"
  if test -f "$srcpkg"; then
    tar xfj $srcpkg
  else
    echo "No source package found ($srcpkg)."
    exit
  fi

  for pkg in $pkgs; do
    echo "Build: $pkg"
    builddir=vpp-build-$pkg
    $package build_bdist --verbose --srcdir=$srcdir		\
        --no-maintainer-mode					\
	--configfile=$cfgfile					\
	--configdir="$dir/scripts"				\
	--builddir=$builddir					\
	$pkg_opts						\
	--package=$pkg 2>&1 > log-build-$pkg
  done
fi

# 3. Test binary packages.
if test "$what" = "test" -o "$what" = "all"; then
  echo "#####################################################################"
  echo "# test binary packages                                              #"
  echo "#####################################################################"
  for pkg in $pkgs; do
    builddir=vpp-build-$pkg

    pkgfile=`ls vpp-build-$pkg/sourceryvsipl++-*.tar.bz2`

    echo "Test: $pkg ($pkgfile)"

    rm -rf $distdir
    $package test_bdist --verbose				\
	--packagefile=$pkgfile					\
	--distdir=$distdir					\
	--srcdir=$test_srcdir					\
	--configfile=$cfgfile					\
	--configdir="$dir/scripts"				\
	--builddir=$builddir					\
	$pkg_opts						\
	--package=$pkg 2>&1 > log-test-$pkg
  done
fi
