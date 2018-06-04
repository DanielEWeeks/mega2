#   Mega2: Manipulation Environment for Genetic Analysis
#  
#   Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.
#  
#   Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   and Daniel E. Weeks.
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

L0xFLAGS= -std=c++11
CFLAGS += -D__Mega2__ -Igenfile/include 

ZSTDDIR = zstd-1.1.0/lib
ZSTD = bgenzstd.a
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
	$(CXX) -I$(ZSTDDIR) $(INCLUDES) $(CFLAGS) $(LOCFLAGS) $(FLAGS) -c $< $(WINO)

all: $(OUT)

more: $(OUT) example/bgen_to_vcf

example/bgen_to_vcf: example/bgen_to_vcf.o
	$(CXX) example/bgen_to_vcf.o $(WINO) $(OUT) -lz

example/bgen_to_vcf.o: example/bgen_to_vcf.cpp
example/bgen_to_vcf.o: CFLAGS += $(L0xFLAGS)


$(OUT): $(ZSTDDIR)$(WINDIR)$(ZSTD) $(BOBJ)
	$(AR) $(ARFLAGS) $(AROUT)$(OUT) $(BOBJ) $(subst /,$(WINDIR),$(wildcard $(ZSTDDIR)/*/*.o))
#	$(AR) $(ARFLAGS) $(AROUT)$(OUT) $(BOBJ) $(ZSTDDIR)\decompress\huf_decompress.o $(ZSTDDIR)\decompress\zbuff_decompress.o $(ZSTDDIR)\decompress\zstd_decompress.o $(ZSTDDIR)\common\entropy_common.o  $(ZSTDDIR)\common\fse_decompress.o  $(ZSTDDIR)\common\xxhash.o  $(ZSTDDIR)\common\zstd_common.o $(ZSTDDIR)\dictBuilder\divsufsort.o  $(ZSTDDIR)\dictBuilder\zdict.o $(ZSTDDIR)\legacy\zstd_v01.o  $(ZSTDDIR)\legacy\zstd_v02.o $(ZSTDDIR)\legacy\zstd_v03.o $(ZSTDDIR)\legacy\zstd_v04.o $(ZSTDDIR)\legacy\zstd_v05.o $(ZSTDDIR)\legacy\zstd_v06.o $(ZSTDDIR)\legacy\zstd_v07.o foo

$(ZSTDDIR)$(WINDIR)$(ZSTD):
	$(MAKE) -C $(ZSTDDIR) $(ZSTD)

clean:
	$(MAKE) -C $(ZSTDDIR) $@
	$(RM) $(WOBJ) $(OUT) Makefile.bak
	$(RM) example$(WINDIR)bgen_to_vcf

distclean: clean

depend:
	makedepend -Y. -- $(CFLAGS) -- $(INCLUDES) $(BSRC)

# DO NOT DELETE

src/MissingValue.o: genfile/include/genfile/MissingValue.hpp
src/bgen.o: genfile/include/genfile/types.hpp
src/bgen.o: genfile/include/genfile/bgen/bgen.hpp
src/bgen.o: genfile/include/genfile/zlib.hpp
src/bgen.o: genfile/include/genfile/MissingValue.hpp
src/zlib.o: genfile/include/genfile/zlib.hpp
src/zlib.o: genfile/include/genfile/types.hpp
