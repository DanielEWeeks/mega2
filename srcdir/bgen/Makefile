#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================
# Makefile for BGEN as a library for Mega2

# Put the library in the Mega2 srcdir above this...
# NOTE: -lz is needed when the library is linked...
ifeq ($(OSTYPE),windows)
UP = ..\\
WINDIR=\\
#WINO = -Fo$(@D)/$*.o
WINO = -Fo$*.o
else
UP = ../
WINDIR=/
WINO = -o $@
endif

#FLAGS += -std=c++11 -Igenfile/include 
CFLAGS += -Igenfile/include 

OUT = $(UP)libbgen.a

#SRC = example/bgen_to_vcf.cpp src/MissingValue.cpp src/bgen.cpp src/zlib.cpp
BSRC = src/MissingValue.cpp src/bgen.cpp src/zlib.cpp
BOBJ = $(BSRC:.cpp=.o)
WOBJ = $(subst /,$(WINDIR),$(BOBJ))

# include directories '.' for other vcftools headers, '..' for Mega2 headers... and zlib headers
LOCFLAGS += -I..

# CC, CFLAGS are defined in Mega2 Makedefs

.SUFFIXES: .cpp
.cpp.o:
	$(CXX) $(INCLUDES) $(CFLAGS) $(LOCFLAGS) $(FLAGS) -c $< $(WINO)

all: $(OUT)

$(OUT): $(BOBJ)
	$(AR) $(ARFLAGS) $(AROUT)$(OUT) $(BOBJ)

clean:
	$(RM) $(WOBJ) $(OUT) Makefile.bak

distclean: clean

depend:
	makedepend -Y. -- $(CFLAGS) -- $(INCLUDES) $(BSRC)

# DO NOT DELETE

example/bgen_to_vcf.o: genfile/include/genfile/bgen/bgen.hpp
example/bgen_to_vcf.o: genfile/include/genfile/zlib.hpp
example/bgen_to_vcf.o: genfile/include/genfile/types.hpp
example/bgen_to_vcf.o: genfile/include/genfile/MissingValue.hpp
src/MissingValue.o: genfile/include/genfile/MissingValue.hpp
src/bgen.o: genfile/include/genfile/types.hpp
src/bgen.o: genfile/include/genfile/bgen/bgen.hpp
src/bgen.o: genfile/include/genfile/zlib.hpp
src/bgen.o: genfile/include/genfile/MissingValue.hpp
src/zlib.o: genfile/include/genfile/zlib.hpp
