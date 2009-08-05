########################################################################
#
# File:   GNUmakefile
# Author: Jules Bergmann
# Date:   2007-07-13
#
# Contents: Standalone Makefile for VSIPL++ apps/utils.
#
########################################################################

PKG = vsipl++

CXX      := $(shell pkg-config --variable=cxx $(PKG))
OPT      := $(shell pkg-config --variable=cxxflags $(PKG))
CXXFLAGS := $(shell pkg-config --cflags $(PKG)) $(OPT)
CFLAGS   := $(CXXFLAGS)
LIBS     := $(shell pkg-config --libs   $(PKG))

OBJEXT   := o

cxx_sources := $(wildcard *.cpp)
c_sources   := $(wildcard *.c)

objects     := $(patsubst %.cpp, %.$(OBJEXT), $(cxx_sources))	\
	       $(patsubst %.c, %.$(OBJEXT), $(c_sources))
deps        := $(patsubst %.cpp, %.d, $(cxx_sources))



########################################################################
# Standard Targets
########################################################################

all::

depend:: $(deps)

clean::

ifeq (,$(filter $(MAKECMDGOALS), depend doc clean))
include $(deps)
endif

ifneq (,$(findstring depend, $(MAKECMDGOALS)))
$(deps): %.d:	.FORCE
endif



########################################################################
# Application Targets
########################################################################

all:: $(APPS)

lsmat: lsmat.o
	$(CXX) $(OPT) -o $@ $^ $(LIBS)

clean::
	rm -f $(objects)
	rm -f $(deps)
	rm -f *.raw



########################################################################
# Implicit Rules
########################################################################

# Generate a dependency Makefile fragment for a C++ source file.
# (This recipe is taken from the GNU Make manual.)
%.d: %.cpp
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) \
		      $(call dir_var,$(dir $<),CXXFLAGS) $< \
		      | sed "s|$(*F)\\.$(OBJEXT)[ :]*|$*\\.d $*\\.$(OBJEXT) : |g" > $@'

