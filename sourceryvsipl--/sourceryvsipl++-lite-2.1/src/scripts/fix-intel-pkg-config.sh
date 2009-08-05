#! /bin/sh

########################################################################
#
# File:   fix-intel-pkg-config.sh
# Author: Jules Bergmann
# Date:   2005-01-19
#
# Contents:
#   Edit pkg-config files to put install prefixes for Intel libraries
#   (IPP and MKL) into variables.
#
########################################################################

# SYNOPSIS
#   fix-intel-pkg-config.sh -p PCFILE [-i IPPDIR] [-m MKLDIR]
#

# IPP Prefix
ipp_prefix="/opt/intel/ipp"

# MKL Prefix
mkl_prefix="/opt/intel/mkl"

drop_ipp_arch="no"

# .pc file
pcfile=""

while getopts "p:i:m:d" arg; do
    case $arg in
	p)
	    pcfile=$OPTARG
	    ;;
	i)
	    ipp_prefix=$OPTARG
	    ;;
	m)
	    mkl_prefix=$OPTARG
	    ;;
	d)
	    drop_ipp_arch="yes";
	    ;;
	\?)
            error "usage: fix-intel-pkg-config.sh -p PCFILE [-i IPPDIR] [-m MKLDIR]"
	    ;;
    esac
done

if test ! -f "$pcfile"; then
  error "error: fix-intel-pkg-config.sh -p PCFILE option required"
fi

if test "$drop_ipp_arch" == "yes"; then
  ipp_prefix=`dirname $ipp_prefix`
fi

echo "ipp_prefix=$ipp_prefix" >  $pcfile.new
echo "mkl_prefix=$mkl_prefix" >> $pcfile.new

cat $pcfile | sed -e "s|$ipp_prefix/|\${ipp_prefix}/|g"	\
            | sed -e "s|$mkl_prefix/|\${mkl_prefix}/|g" >> $pcfile.new

if test -f "$pcfile.new"; then
  mv $pcfile.new $pcfile
fi
