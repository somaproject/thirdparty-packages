######################################################### -*-Makefile-*-
#
# File:   apps/vadd/Makefile
# Author: Mike LeBlanc
# Date:   2009-04-07
#
# Contents: Makefile for vector addition parallel application.
#
########################################################################

########################################################################
# Variables
########################################################################

# This should point to the directory where Sourcery VSIPL++ is installed.
# If installed in /usr/local, or if PKG_CONFIG_PATH is set correctly,
# then this may be left blank, as pkg-config will obtain the path
# from the .pc file.  If you set prefix here, it overrides the above
# environment variable.
prefix =

# This selects the desired library, which will vary depending on the
# platform.  Consult the pkgconfig directory for a complete list of 
# packages (vsipl++-*.pc).  If a symbolic link is used for vsipl++.pc,
# and it is in the path that pkg-config normally searches, then
# this may be left blank (preferred).  Append '-debug' for building a 
# version suitable for debugging or omit to use the optimized version.  
suffix =

ifeq ($(strip $(prefix)),)
pkgcommand := pkg-config vsipl++$(suffix)
fmt-profile-command := $(subst /lib/pkgconfig,,$(PKG_CONFIG_PATH))/bin/fmt-profile.pl
else
pkgcommand := PKG_CONFIG_PATH=$(prefix)/lib/pkgconfig 	\
                     pkg-config vsipl++$(suffix) 	\
                     --define-variable=prefix=$(prefix)
fmt-profile-command := $(prefix)/bin/fmt-profile.pl
endif


CXX      = $(shell ${pkgcommand} --variable=cxx)
CXXFLAGS = $(shell ${pkgcommand} --cflags) \
	   $(shell ${pkgcommand} --variable=cxxflags) \
	   -DVSIP_IMPL_PROFILER=11
LIBS     = $(shell ${pkgcommand} --libs)
 

########################################################################
# Rules
########################################################################

all: vadd

show:
	@echo "pkgcommand: " $(pkgcommand)
	@echo "CXX       : " $(CXX)
	@echo "CXXFLAGS  : " $(CXXFLAGS)
	@echo "LIBS      : " $(LIBS)

clean: 
	rm -f *.o
	rm -f vadd
	rm -f viewtoraw
	rm -f diffview

check: all
	@echo "Running VADD application..."
	./vadd -loop 10

vadd.o: vadd.cpp

vadd: vadd.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

profile: vadd
	@echo "Profiling VADD application ..."
	$(RUN) ./vadd -loop 1000 -N $N --vsipl++-profile-mode=accum --vsipl++-profile-output=vadd.out
	@echo "Formatting profiler output..."
	${fmtprof}  -sec -o vadd-$(NP).txt vadd.out