######################################################### -*-Makefile-*-
#
# File:   apps/ssar/Makefile
# Author: Don McCoy
# Date:   2006-10-28
#
# Contents: Makefile for Scalable SAR (SSAR) benchmark program.
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

# The default precision is single (double may also be used)
precision = single

ifeq ($(precision),double)
ref_image_base := ref_image_dp
ssar_type := SSAR_BASE_TYPE=double
else
ref_image_base := ref_image_sp
ssar_type := SSAR_BASE_TYPE=float
endif

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
	   -DVSIP_IMPL_PROFILER=11 -D$(ssar_type)
LIBS     = $(shell ${pkgcommand} --libs)
 

########################################################################
# Rules
########################################################################

all: ssar viewtopng diffview

show:
	@echo "pkgcommand: " $(pkgcommand)
	@echo "CXX       : " $(CXX)
	@echo "CXXFLAGS  : " $(CXXFLAGS)
	@echo "LIBS      : " $(LIBS)

clean: 
	rm -f *.o
	rm -f ssar
	rm -f viewtoraw
	rm -f diffview

check: all
	@echo "Running SSAR application..."
	./ssar data3
	@echo
	@echo "Comparing output to reference view (should be less than -100)"
	./diffview -r data3/image.view data3/$(ref_image_base).view 756 1144
	@echo
	@echo "Creating viewable image of output"
	./viewtopng -s data3/image.view data3/image.png 1144 756
	@echo "Creating viewable image of reference view"
	./viewtopng -s data3/$(ref_image_base).view data3/$(ref_image_base).png 1144 756

profile1: ssar viewtopng
	@echo "Profiling SSAR application (SCALE = 1)..."
	./ssar data1 -loop 10 --vsipl++-profile-mode=accum --vsipl++-profile-output=profile.out
	@echo "Formatting profiler output..."
	${fmt-profile-command}  -sec -o profile1.txt data1/profile.out
	./make_images.sh data1 438 160 382 266

profile3: ssar viewtopng
	@echo "Profiling SSAR application (SCALE = 3)..."
	./ssar data3 -loop 10 --vsipl++-profile-mode=accum --vsipl++-profile-output=profile.out
	@echo "Formatting profiler output..."
	${fmt-profile-command}  -sec -o profile3.txt data3/profile.out
	./make_images.sh data3 1072 480 1144 756

ssar.o: ssar.cpp kernel1.hpp

ssar: ssar.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

viewtopng: viewtopng.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS) -lpng

diffview: diffview.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
