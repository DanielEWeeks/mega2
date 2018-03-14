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
# Makefile for VCFTOOLS as a library for Mega2
# Requires gnu make (sometimes called 'gmake')
#

# NOTE: the '-DNDEBUG' is needed to turn off the assert macros in the VCFtools code...
LOCFLAGS += -DNDEBUG

# Put the library in the Mega2 srcdir above this...
# NOTE: -lz is needed when the library is linked...
ifeq ($(OSTYPE),windows)
UP = ..\\
WINO = -Fo$*.o
else
UP = ../
WINO = 
endif

OUT = $(UP)libvcftools.a

# source files.
SRC =  	mega2_vcftools_interface.cpp \
	parameters.cpp header.cpp \
	entry.cpp entry_getters.cpp entry_setters.cpp bcf_entry.cpp bcf_entry_setters.cpp \
	vcf_entry.cpp vcf_entry_setters.cpp \
	variant_file.cpp variant_file_filters.cpp variant_file_index.cpp bcf_file.cpp vcf_file.cpp

OBJ = $(SRC:.cpp=.o)

# include directories '.' for other vcftools headers, '..' for Mega2 headers... and zlib headers
LOCFLAGS += -I..

# CC, CFLAGS are defined in Mega2 Makedefs

.SUFFIXES: .cpp

default: $(OUT)

.cpp.o:
	$(CXX) $(INCLUDES) $(CFLAGS) $(LOCFLAGS) $(FLAGS) -c $< $(WINO)

.c.o:
	$(CC) $(INCLUDES) $(CFLAGS) $(LOCFLAGS) $(FLAGS) -c $< $(WINO)

$(OUT): $(OBJ) #bgzf.o
	$(AR) $(ARFLAGS) $(AROUT)$(OUT) $(OBJ) #bgzf.o

depend:
	makedepend -Y. -- $(CFLAGS) $(LOCFLAGS) -- $(INCLUDES) $(SRC) bgzf.c

clean:
	$(RM) $(OBJ) bgzf.o $(OUT) Makefile.bak

distclean: clean


# DO NOT DELETE

mega2_vcftools_interface.o: ../common.h ../types.hh ../typedefs.h
mega2_vcftools_interface.o: ../errorno.h ../allutil.h ../linkage.h
mega2_vcftools_interface.o: ../linkage_external.h ../makeped.h ../pedtree.h
mega2_vcftools_interface.o: ../batch_input.h ../input_check.h ../hwe.h
mega2_vcftools_interface.o: ../compress_ext.h ../tod.hh ../utils_ext.h
mega2_vcftools_interface.o: ../genetic_utils_ext.h ../error_messages_ext.h
mega2_vcftools_interface.o: ../annotated_ped_file.h ../read_files_ext.h
mega2_vcftools_interface.o: ../marker_lookup_ext.h ../phe_lookup_ext.h
mega2_vcftools_interface.o: ../reorder_loci_ext.h output_log.h parameters.h
mega2_vcftools_interface.o: variant_file.h entry.h bgzf.h vcf_entry.h
mega2_vcftools_interface.o: bcf_entry.h header.h bcf_file.h vcf_file.h
mega2_vcftools_interface.o: mega2_vcftools_interface.h
mega2_vcftools_interface.o: ../annotated_ped_file.h
parameters.o: parameters.h output_log.h
header.o: header.h entry.h bgzf.h output_log.h
entry.o: entry.h bgzf.h output_log.h
entry_getters.o: entry.h bgzf.h output_log.h
entry_setters.o: entry.h bgzf.h output_log.h
bcf_entry.o: bcf_entry.h output_log.h entry.h bgzf.h header.h
bcf_entry_setters.o: bcf_entry.h output_log.h entry.h bgzf.h header.h
vcf_entry.o: vcf_entry.h entry.h bgzf.h output_log.h
vcf_entry_setters.o: vcf_entry.h entry.h bgzf.h output_log.h
variant_file.o: variant_file.h parameters.h output_log.h entry.h bgzf.h
variant_file.o: vcf_entry.h bcf_entry.h header.h
variant_file_filters.o: variant_file.h parameters.h output_log.h entry.h
variant_file_filters.o: bgzf.h vcf_entry.h bcf_entry.h header.h
variant_file_index.o: variant_file.h parameters.h output_log.h entry.h bgzf.h
variant_file_index.o: vcf_entry.h bcf_entry.h header.h
bcf_file.o: bcf_file.h output_log.h parameters.h variant_file.h entry.h
bcf_file.o: bgzf.h vcf_entry.h bcf_entry.h header.h
vcf_file.o: vcf_file.h output_log.h vcf_entry.h entry.h bgzf.h parameters.h
vcf_file.o: variant_file.h bcf_entry.h header.h
bgzf.o: bgzf.h
