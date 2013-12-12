#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2013 Robert Baron, Charles P. Kollar,
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
# Makefile for VCFTOOLS as a library for Mega2
# Requires gnu make (sometimes called 'gmake')
#

OUT = ../libvcftools.a

# source files.
SRC =  	mega2_vcftools_interface.cpp \
	parameters.cpp header.cpp \
	entry.cpp entry_getters.cpp entry_setters.cpp bcf_entry.cpp bcf_entry_setters.cpp vcf_entry.cpp vcf_entry_setters.cpp \
	variant_file.cpp variant_file_filters.cpp variant_file_index.cpp bcf_file.cpp vcf_file.cpp

OBJ = $(SRC:.cpp=.o)

# include directories '.' for other vcftools headers, '..' for Mega2 headers...
INCLUDES = -I. -I..

# compiler flags
CCFLAGS = -g -D_FILE_OFFSET_BITS=64

# compiler variables are set by calling makefile

# library paths
LIBS = -lz

.SUFFIXES: .cpp

default: $(OUT)

.cpp.o:
	$(CXX) $(INCLUDES) $(CCFLAGS) -c $< -o $@

bgzf: bgzf.c bgzf.h

$(OUT): $(OBJ)
	$(CC) -c $(CCFLAGS) -m64 bgzf.c $(LIB) -o bgzf.o
	ar rcs $(OUT) $(OBJ) bgzf.o

depend:
	makedepend -Y. -- $(CFLAGS) -- $(INCLUDES) $(SRC) bgzf.c

clean:
	$(RM) $(OBJ) bgzf.o $(OUT) Makefile.bak

# DO NOT DELETE

mega2_vcftools_interface.o: ../common.h ../typedefs.h ../errorno.h
mega2_vcftools_interface.o: ../allutil.h ../linkage.h ../linkage_external.h
mega2_vcftools_interface.o: ../makeped.h ../pedtree.h ../analysis.h
mega2_vcftools_interface.o: ../analysis_ext.h ../write_fbat_ext.h
mega2_vcftools_interface.o: ../write_plink_ext.h ../plink_core_ext.h
mega2_vcftools_interface.o: ../write_eigenstrat_ext.h ../write_pangaea_ext.h
mega2_vcftools_interface.o: ../write_beagle_ext.h ../write_structure_ext.h
mega2_vcftools_interface.o: ../write_pseq_ext.h ../batch_input.h
mega2_vcftools_interface.o: ../input_check.h ../hwe.h ../compress_ext.h
mega2_vcftools_interface.o: ../utils_ext.h ../error_messages_ext.h
mega2_vcftools_interface.o: output_log.h parameters.h variant_file.h entry.h
mega2_vcftools_interface.o: bgzf.h vcf_entry.h bcf_entry.h header.h
mega2_vcftools_interface.o: bcf_file.h vcf_file.h
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
