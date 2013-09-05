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
# Makefile for compiling Mega2
# Requires gnu make (sometimes called 'gmake')
#

TOP = .
MODULE = mega2
SRCDIR = .
REMOVECVS1 = "../remove_cvs_stuff1.gawk"
REMOVECVS = "../remove_cvs_stuff.gawk"

ifeq ($(OSTYPE), )
    ifeq "$(OS)" "Windows_NT"
	ifeq ($(shell where uname), )
	    OSTYPE = windows
	else
	    OSTYPE := $(shell uname -s)
	endif
    else
        OSTYPE := $(shell uname -s)
    endif
    ifeq "$(OSTYPE)" "Darwin"
	OSTYPE = darwin
    else ifeq "$(OSTYPE)" "Linux"
	OSTYPE = linux
    else ifeq "$(OSTYPE)" ""
        $(error OSTYPE is not specified)
    else ifeq "$(findstring CYGWIN,$(OSTYPE))" "CYGWIN"
	OSTYPE = cygwin
    else ifeq "$(findstring MINGW32,$(OSTYPE))" "MINGW32"
	OSTYPE = mingw
    else
        $(warning OSTYPE is $(OSTYPE))
    endif
endif

EXEC = $(MODULE)_$(OSTYPE)

#CLIBS = -lm
HDRS =  typedefs.h  common.h errorno.h \
	allutil.h   hwe.h  liable.h  linkage.h  pedtree.h \
	R_output.h  makeped.h linkage_ext.h annotated_ped_file.h \
	input_check.h batch_input.h error_sim.h mrecode.h \
	write_mfiles.h

SRCS =  input_check.cpp     select_index.cpp \
	aspex.cpp	    liability_dist.cpp   select_individuals.cpp \
	linkage.cpp	    stamp_time.cpp       output_file_names.cpp \
	cal_const.cpp	    list.cpp	         test_switch.cpp \
	cal_n.cpp	    ln_p_value.cpp	 user_input.cpp \
	cal_prob.cpp	    log_factorial.cpp    utils.cpp \
	check.cpp	    mega2.cpp	         write_apm.cpp \
	check_file.cpp	    output_routines.cpp  write_files.cpp \
	do_switch.cpp	    pedtree.cpp	         write_ghfiles.cpp \
	fcmap.cpp	    print_data.cpp	 write_liability.cpp \
	random_choose.cpp   write_mfiles.cpp     mrecode.cpp \
	gen_user_input.cpp  write_sage_files.cpp splink.cpp \
	genetic_utils.cpp   read_data.cpp	 write_simulate_files.cpp \
	homozygosity.cpp    read_files.cpp	 write_solar_files.cpp \
	hwe.cpp		    reorder_loci.cpp     write_vitesse.cpp \
	hwe_user_input.cpp  scripts.cpp	         cw_routines.cpp \
	makeped1.cpp        omit_ped.cpp         slink.cpp \
	R_output.cpp        write_premakeped.cpp create_summary.cpp \
	write_prest.cpp	    write_pap.cpp        batch_input.cpp \
	makenucs.cpp        append_locus_array.cpp error_sim.cpp \
	write_loki.cpp      error_messages.cpp   annotated_ped_file.cpp \
	                    write_SUP_files.cpp  write_mendel7_files.cpp \
	write_plink_files.cpp cranefoot.cpp      write_mega2annot.cpp \
	grow_string.cpp     marker_lookup.cpppp  write_IQLS.cpp \
	net.cpp             phe_lookup.cpp       analysis.cpp \
	write_fbat.cpp      entry.cpp            loop.cpp \
	write_pangaea.cpp   write_beagle.cpp     plink_core_files.cpp \
	write_eigenstrat.cpp write_structure.cpp batch.cpp \
	compress.cpp

OBJS =  input_check.o     select_index.o \
	aspex.o		  liability_dist.o   select_individuals.o \
	linkage.o	  stamp_time.o       output_file_names.o \
	cal_const.o	  list.o	     test_switch.o \
	cal_n.o		  ln_p_value.o	     user_input.o \
	cal_prob.o	  log_factorial.o    utils.o \
	check.o		  mega2.o	     write_apm.o \
	check_file.o	  output_routines.o  write_files.o \
	do_switch.o	  pedtree.o	     write_ghfiles.o \
	fcmap.o		  print_data.o	     write_liability.o \
	random_choose.o   write_mfiles.o     mrecode.o\
	gen_user_input.o  write_sage_files.o splink.o \
	genetic_utils.o   read_data.o	     write_simulate_files.o \
	homozygosity.o	  read_files.o	     write_solar_files.o \
	hwe.o		  reorder_loci.o     write_vitesse.o \
	hwe_user_input.o  scripts.o          cw_routines.o \
	makeped1.o        omit_ped.o         slink.o \
	R_output.o  	  write_premakeped.o create_summary.o \
	write_prest.o     write_pap.o        batch_input.o \
	makenucs.o        append_locus_array.o error_sim.o \
	write_loki.o      error_messages.o   annotated_ped_file.o \
	                  write_SUP_files.o  write_mendel7_files.o \
	write_plink_files.o cranefoot.o      write_mega2annot.o \
	grow_string.o     marker_lookup.o    write_IQLS.o \
	net.o             phe_lookup.o       analysis.o \
	write_fbat.o      entry.o            loop.o \
	write_pangaea.o   write_beagle.o     plink_core_files.o \
	write_eigenstrat.o write_structure.o batch.o \
	compress.o

# The flags for including the libraries:

IFLAGS = -I$(TOP)

ifneq ($(MINGW),)
	RM=del 
else
	RM=rm -f
endif

$(info including $(TOP)/Makedefs_$(OSTYPE))
-include $(TOP)/Makedefs_$(OSTYPE)
ifndef MAKEDEFS
$(info including $(TOP)/Makedefs_$(OSTYPE) failed.  Trying include $(TOP)/Makedefs_$(OSTYPE).template)
include $(TOP)/Makedefs_$(OSTYPE).template
endif

# CC, CFLAGS are defined in Makedefs

DFLAGS += -DCFREE -DRUNSHELL_SETUP

# Used for getting a direct match with the genetic map files with running Beagle Analysis...
#DFLAGS += -DDONT_CONVERT_TO_HALDANE_FOR_BEAGLE

# The not yet ready for prime time object oriented batch file interface....
#DFLAGS += -DNEW_BATCH

mega2.o utils.o: DFLAGS += -DEXPIRE

ifdef TEST
    DtFLAGS += -DTEST
endif
ifdef HIDE
    DdFLAGS += -DHIDEDATE
    DpFLAGS += -DHIDEPATH
    DfFLAGS += -DHIDEFILE
    DsFLAGS += -DHIDESTATUS
endif

# -DTEST
R_output.o hwe_user_input.o mega2.o mrecode.o utils.o read_files.o write_sage_files.o: DFLAGS += $(DtFLAGS)

# -DDATE
stamp_time.o: DFLAGS += $(DdFLAGS)

# -DHIDESTATUS
aspex.o annotated_ped_file.o batch_input.o genetic_utils.o makeped1.o mega2.o mrecode.o: DFLAGS += $(DsFLAGS)
output_routines.o read_files.o reorder_loci.o user_input.o utils.o write_files.o: DFLAGS += $(DsFLAGS)
write_plink_files.o: DFLAGS += $(DsFLAGS)

# -DHIDEPATH
aspex.o cranefoot.o mega2.o R_output.o write_SUP_files.o write_ghfiles.o write_premakeped.o: DFLAGS += $(DpFLAGS)

# -DHIDEFILE
annotated_ped_file.o makeped1.o mega2.o read_files.o utils.o: DFLAGS += $(DfFLAGS)
write_IQLS.o write_plink_files.o: DFLAGS += $(DfFLAGS)

# -DDATE -DHIDEPATH
write_mfiles.o write_vitesse.o utils.o: DFLAGS += $(DdFLAGS) $(DpFLAGS)

ifeq  "$(HIDE)" "2"
    DaFLAGS += -DSHOWSTATUS
endif
# -DSHOWSTATUS
annotated_ped_file.o mrecode.o: DFLAGS += $(DaFLAGS)

ifeq  "$(TEST)" "2"
    DlFLAGS += -DLOGDIR
    DlFLAGS += -DORDER_HETEROZYGOTE
endif
output_file_names.o: DFLAGS += $(DlFLAGS)
read_files.o compress.o: DFLAGS += $(DlFLAGS)

.c.o:
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CC) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(IFLAGS) $*.c

.cpp.o:	
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CXX) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(IFLAGS) $*.cpp


all: $(MODULE)_$(OSTYPE)

gcc_info.o: $(OBJS)

LINK ?= $(CXX)
LTARGET ?= -o $(EXEC)
$(MODULE)_$(OSTYPE): $(OBJS) $(HDRS) gcc_info.o
	$(RM) $(MODULE)
	$(LINK) $(LFLAGS) $(LTARGET) $(OBJS) $(FLAGS) gcc_info.o $(CLIBS)

ifeq ($(OSTYPE),mingw)
	cp $(MODULE)_$(OSTYPE).exe $(MODULE).exe
else ifeq ($(OSTYPE),windows)
	copy $(MODULE)_$(OSTYPE).exe $(MODULE).exe
else
	ln -s $(MODULE)_$(OSTYPE) $(MODULE)
endif

clean:
	$(RM) $(OBJS) gcc_info.o $(EXEC) $(MODULE) $(EXEC).exe $(MODULE).exe

clean1:
	$(RM) *~

srcdist:	
	for i in $(SRCS) $(HDRS); do $(REMOVECVS1) $$i > ../srcdir/$$i; done
	for i in `ls Makedefs*`; do $(REMOVECVS) $$i > ../srcdir/$$i.template; done
	$(REMOVECVS) Makefile > ../srcdir/Makefile

depend:
	makedepend -Y. -- $(CFLAGS) -- $(SRCS)

# DO NOT DELETE

input_check.o: common.h typedefs.h errorno.h allutil.h linkage.h
input_check.o: linkage_external.h makeped.h pedtree.h analysis.h
input_check.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
input_check.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
input_check.o: write_beagle_ext.h write_structure_ext.h batch_input.h
input_check.o: input_check.h hwe.h compress_ext.h check_ext.h
input_check.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
input_check.o: genetic_utils_ext.h grow_string_ext.h input_check_ext.h
input_check.o: pedtree_ext.h utils_ext.h write_pap_ext.h batch_input_ext.h
select_index.o: common.h typedefs.h errorno.h allutil.h linkage.h
select_index.o: linkage_external.h makeped.h pedtree.h analysis.h
select_index.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
select_index.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
select_index.o: write_beagle_ext.h write_structure_ext.h batch_input.h
select_index.o: input_check.h hwe.h compress_ext.h random_choose_ext.h
aspex.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
aspex.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
aspex.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
aspex.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
aspex.o: batch_input.h input_check.h hwe.h compress_ext.h batch_input_ext.h
aspex.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
aspex.o: genetic_utils_ext.h grow_string_ext.h linkage_ext.h makenucs_ext.h
aspex.o: omit_ped_ext.h output_file_names_ext.h output_routines_ext.h
aspex.o: scripts_ext.h user_input_ext.h utils_ext.h
liability_dist.o: common.h typedefs.h errorno.h allutil.h linkage.h
liability_dist.o: linkage_external.h makeped.h pedtree.h analysis.h
liability_dist.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
liability_dist.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
liability_dist.o: write_beagle_ext.h write_structure_ext.h batch_input.h
liability_dist.o: input_check.h hwe.h compress_ext.h liable.h
liability_dist.o: error_messages_ext.h genetic_utils_ext.h
liability_dist.o: liability_dist_ext.h linkage_ext.h omit_ped_ext.h
liability_dist.o: output_file_names_ext.h utils_ext.h
select_individuals.o: common.h typedefs.h errorno.h allutil.h linkage.h
select_individuals.o: linkage_external.h makeped.h pedtree.h analysis.h
select_individuals.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
select_individuals.o: plink_core_ext.h write_eigenstrat_ext.h
select_individuals.o: write_pangaea_ext.h write_beagle_ext.h
select_individuals.o: write_structure_ext.h batch_input.h input_check.h hwe.h
select_individuals.o: compress_ext.h batch_input_ext.h error_messages_ext.h
select_individuals.o: fcmap_ext.h utils_ext.h
linkage.o: common.h typedefs.h errorno.h allutil.h linkage.h
linkage.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
linkage.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
linkage.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
linkage.o: write_structure_ext.h batch_input.h input_check.h hwe.h
linkage.o: compress_ext.h error_messages_ext.h list_ext.h makeped_ext.h
linkage.o: pedtree_ext.h utils_ext.h user_input_ext.h
linkage.o: append_locus_array_ext.h
output_file_names.o: common.h typedefs.h errorno.h allutil.h linkage.h
output_file_names.o: linkage_external.h makeped.h pedtree.h analysis.h
output_file_names.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
output_file_names.o: plink_core_ext.h write_eigenstrat_ext.h
output_file_names.o: write_pangaea_ext.h write_beagle_ext.h
output_file_names.o: write_structure_ext.h batch_input.h input_check.h hwe.h
output_file_names.o: compress_ext.h batch_input_ext.h error_messages_ext.h
output_file_names.o: fcmap_ext.h genetic_utils_ext.h output_file_names_ext.h
output_file_names.o: utils_ext.h
cal_const.o: common.h typedefs.h errorno.h allutil.h linkage.h
cal_const.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cal_const.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cal_const.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cal_const.o: write_structure_ext.h batch_input.h input_check.h hwe.h
cal_const.o: compress_ext.h hwe_ext.h
list.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
list.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
list.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
list.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
list.o: batch_input.h input_check.h hwe.h compress_ext.h error_messages_ext.h
test_switch.o: common.h typedefs.h errorno.h allutil.h linkage.h
test_switch.o: linkage_external.h makeped.h pedtree.h analysis.h
test_switch.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
test_switch.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
test_switch.o: write_beagle_ext.h write_structure_ext.h batch_input.h
test_switch.o: input_check.h hwe.h compress_ext.h
cal_n.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
cal_n.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
cal_n.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
cal_n.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
cal_n.o: batch_input.h input_check.h hwe.h compress_ext.h
ln_p_value.o: common.h typedefs.h errorno.h allutil.h linkage.h
ln_p_value.o: linkage_external.h makeped.h pedtree.h analysis.h
ln_p_value.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
ln_p_value.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
ln_p_value.o: write_beagle_ext.h write_structure_ext.h batch_input.h
ln_p_value.o: input_check.h hwe.h compress_ext.h hwe_ext.h
user_input.o: common.h typedefs.h errorno.h allutil.h linkage.h
user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
user_input.o: write_beagle_ext.h write_structure_ext.h batch_input.h
user_input.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
user_input.o: cw_routines_ext.h error_messages_ext.h fcmap_ext.h
user_input.o: genetic_utils_ext.h grow_string_ext.h output_file_names_ext.h
user_input.o: plink_ext.h utils_ext.h
cal_prob.o: common.h typedefs.h errorno.h allutil.h linkage.h
cal_prob.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cal_prob.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cal_prob.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cal_prob.o: write_structure_ext.h batch_input.h input_check.h hwe.h
cal_prob.o: compress_ext.h do_switch_ext.h test_switch_ext.h utils_ext.h
utils.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
utils.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
utils.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
utils.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
utils.o: batch_input.h input_check.h hwe.h compress_ext.h batch_input_ext.h
utils.o: error_messages_ext.h fcmap_ext.h net_ext.h user_input_ext.h
utils.o: utils_ext.h version.h
check.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
check.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
check.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
check.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
check.o: batch_input.h input_check.h hwe.h compress_ext.h utils_ext.h
check.o: error_messages_ext.h grow_string_ext.h
mega2.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
mega2.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
mega2.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
mega2.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
mega2.o: batch_input.h input_check.h hwe.h compress_ext.h R_output.h
mega2.o: annotated_ped_file_ext.h aspex_ext.h batch_input_ext.h
mega2.o: cranefoot_ext.h create_summary_ext.h error_messages_ext.h
mega2.o: error_sim_ext.h fcmap_ext.h genetic_utils_ext.h hwe_user_input_ext.h
mega2.o: input_check_ext.h makenucs_ext.h makeped_ext.h marker_lookup_ext.h
mega2.o: mega2annot_ext.h output_file_names_ext.h output_routines_ext.h
mega2.o: pedtree_ext.h read_files_ext.h reorder_loci_ext.h scripts_ext.h
mega2.o: slink_ext.h splink_ext.h user_input_ext.h utils_ext.h
mega2.o: write_IQLS_ext.h write_SUP_files_ext.h write_apm_ext.h
mega2.o: write_files_ext.h write_ghfiles_ext.h write_loki_ext.h
mega2.o: write_mfiles_ext.h write_pap_ext.h write_premakeped_ext.h
mega2.o: write_prest_ext.h write_sage_files_ext.h write_simulate_files_ext.h
mega2.o: write_solar_files_ext.h write_vitesse_ext.h version.h
write_apm.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_apm.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
write_apm.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
write_apm.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
write_apm.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_apm.o: compress_ext.h error_messages_ext.h fcmap_ext.h
write_apm.o: genetic_utils_ext.h input_check_ext.h linkage_ext.h
write_apm.o: omit_ped_ext.h output_file_names_ext.h pedtree_ext.h
write_apm.o: reorder_loci_ext.h user_input_ext.h utils_ext.h
output_routines.o: common.h typedefs.h errorno.h allutil.h linkage.h
output_routines.o: linkage_external.h makeped.h pedtree.h analysis.h
output_routines.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
output_routines.o: plink_core_ext.h write_eigenstrat_ext.h
output_routines.o: write_pangaea_ext.h write_beagle_ext.h
output_routines.o: write_structure_ext.h batch_input.h input_check.h hwe.h
output_routines.o: compress_ext.h error_messages_ext.h grow_string_ext.h
output_routines.o: output_routines_ext.h utils_ext.h
write_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_files.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_files.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_files.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_files.o: fcmap_ext.h genetic_utils_ext.h grow_string_ext.h
write_files.o: linkage_ext.h makeped_ext.h omit_ped_ext.h
write_files.o: output_file_names_ext.h reorder_loci_ext.h user_input_ext.h
write_files.o: utils_ext.h write_ghfiles_ext.h plink_ext.h
do_switch.o: common.h typedefs.h errorno.h allutil.h linkage.h
do_switch.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
do_switch.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
do_switch.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
do_switch.o: write_structure_ext.h batch_input.h input_check.h hwe.h
do_switch.o: compress_ext.h
pedtree.o: common.h typedefs.h errorno.h allutil.h linkage.h
pedtree.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
pedtree.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
pedtree.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
pedtree.o: write_structure_ext.h batch_input.h input_check.h hwe.h
pedtree.o: compress_ext.h error_messages_ext.h linkage_ext.h
write_ghfiles.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_ghfiles.o: linkage_external.h makeped.h pedtree.h analysis.h
write_ghfiles.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_ghfiles.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_ghfiles.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_ghfiles.o: input_check.h hwe.h compress_ext.h R_output.h R_output_ext.h
write_ghfiles.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_ghfiles.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
write_ghfiles.o: output_file_names_ext.h output_routines_ext.h
write_ghfiles.o: user_input_ext.h utils_ext.h write_files_ext.h
fcmap.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
fcmap.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
fcmap.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
fcmap.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
fcmap.o: batch_input.h input_check.h hwe.h compress_ext.h
fcmap.o: error_messages_ext.h utils_ext.h
print_data.o: common.h typedefs.h errorno.h allutil.h linkage.h
print_data.o: linkage_external.h makeped.h pedtree.h analysis.h
print_data.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
print_data.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
print_data.o: write_beagle_ext.h write_structure_ext.h batch_input.h
print_data.o: input_check.h hwe.h compress_ext.h
write_liability.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_liability.o: linkage_external.h makeped.h pedtree.h analysis.h
write_liability.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_liability.o: plink_core_ext.h write_eigenstrat_ext.h
write_liability.o: write_pangaea_ext.h write_beagle_ext.h
write_liability.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_liability.o: compress_ext.h liable.h error_messages_ext.h fcmap_ext.h
write_liability.o: genetic_utils_ext.h output_file_names_ext.h
write_liability.o: user_input_ext.h
random_choose.o: common.h typedefs.h errorno.h allutil.h linkage.h
random_choose.o: linkage_external.h makeped.h pedtree.h analysis.h
random_choose.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
random_choose.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
random_choose.o: write_beagle_ext.h write_structure_ext.h batch_input.h
random_choose.o: input_check.h hwe.h compress_ext.h utils_ext.h
write_mfiles.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_mfiles.o: linkage_external.h makeped.h pedtree.h analysis.h
write_mfiles.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_mfiles.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_mfiles.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_mfiles.o: input_check.h hwe.h compress_ext.h R_output.h write_mfiles.h
write_mfiles.o: batch_input_ext.h R_output_ext.h create_summary_ext.h
write_mfiles.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_mfiles.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_mfiles.o: output_routines_ext.h user_input_ext.h utils_ext.h
write_mfiles.o: write_mendel7_files_ext.h write_mfiles_ext.h
mrecode.o: common.h typedefs.h errorno.h allutil.h linkage.h
mrecode.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
mrecode.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
mrecode.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
mrecode.o: write_structure_ext.h batch_input.h input_check.h hwe.h
mrecode.o: compress_ext.h mrecode.h error_messages_ext.h fcmap_ext.h
mrecode.o: grow_string_ext.h linkage_ext.h makeped_ext.h marker_lookup_ext.h
mrecode.o: omit_ped_ext.h read_files_ext.h select_individuals_ext.h
mrecode.o: utils_ext.h write_files_ext.h
gen_user_input.o: common.h typedefs.h errorno.h allutil.h linkage.h
gen_user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
gen_user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
gen_user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
gen_user_input.o: write_beagle_ext.h write_structure_ext.h batch_input.h
gen_user_input.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
gen_user_input.o: fcmap_ext.h genetic_utils_ext.h output_file_names_ext.h
gen_user_input.o: user_input_ext.h utils_ext.h
write_sage_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_sage_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_sage_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_sage_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_sage_files.o: write_pangaea_ext.h write_beagle_ext.h
write_sage_files.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_sage_files.o: compress_ext.h create_summary_ext.h error_messages_ext.h
write_sage_files.o: fcmap_ext.h genetic_utils_ext.h linkage_ext.h
write_sage_files.o: makenucs_ext.h omit_ped_ext.h output_file_names_ext.h
write_sage_files.o: output_routines_ext.h read_files_ext.h user_input_ext.h
write_sage_files.o: utils_ext.h
splink.o: common.h typedefs.h errorno.h allutil.h linkage.h
splink.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
splink.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
splink.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
splink.o: write_structure_ext.h batch_input.h input_check.h hwe.h
splink.o: compress_ext.h create_summary_ext.h error_messages_ext.h
splink.o: fcmap_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
splink.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
splink.o: utils_ext.h
genetic_utils.o: common.h typedefs.h errorno.h allutil.h linkage.h
genetic_utils.o: linkage_external.h makeped.h pedtree.h analysis.h
genetic_utils.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
genetic_utils.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
genetic_utils.o: write_beagle_ext.h write_structure_ext.h batch_input.h
genetic_utils.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
genetic_utils.o: grow_string_ext.h utils_ext.h user_input_ext.h
read_data.o: common.h typedefs.h errorno.h allutil.h linkage.h
read_data.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
read_data.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
read_data.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
read_data.o: write_structure_ext.h batch_input.h input_check.h hwe.h
read_data.o: compress_ext.h
write_simulate_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_simulate_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_simulate_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_simulate_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_simulate_files.o: write_pangaea_ext.h write_beagle_ext.h
write_simulate_files.o: write_structure_ext.h batch_input.h input_check.h
write_simulate_files.o: hwe.h compress_ext.h error_messages_ext.h fcmap_ext.h
write_simulate_files.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_simulate_files.o: output_file_names_ext.h user_input_ext.h utils_ext.h
write_simulate_files.o: write_files_ext.h write_ghfiles_ext.h
homozygosity.o: common.h typedefs.h errorno.h allutil.h linkage.h
homozygosity.o: linkage_external.h makeped.h pedtree.h analysis.h
homozygosity.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
homozygosity.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
homozygosity.o: write_beagle_ext.h write_structure_ext.h batch_input.h
homozygosity.o: input_check.h hwe.h compress_ext.h genetic_utils_ext.h
read_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
read_files.o: linkage_external.h makeped.h pedtree.h analysis.h
read_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
read_files.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
read_files.o: write_beagle_ext.h write_structure_ext.h batch_input.h
read_files.o: input_check.h hwe.h compress_ext.h mrecode.h
read_files.o: annotated_ped_file.h error_sim.h error_messages_ext.h
read_files.o: fcmap_ext.h grow_string_ext.h linkage_ext.h list_ext.h
read_files.o: makeped_ext.h marker_lookup_ext.h mrecode_ext.h omit_ped_ext.h
read_files.o: plink_ext.h read_files_ext.h reorder_loci_ext.h utils_ext.h
read_files.o: user_input_ext.h batch_input_ext.h write_files_ext.h
write_solar_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_solar_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_solar_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_solar_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_solar_files.o: write_pangaea_ext.h write_beagle_ext.h
write_solar_files.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_solar_files.o: compress_ext.h error_messages_ext.h fcmap_ext.h
write_solar_files.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_solar_files.o: output_file_names_ext.h user_input_ext.h utils_ext.h
hwe.o: hwe.h check_file_ext.h hwe_ext.h read_data_ext.h
reorder_loci.o: common.h typedefs.h errorno.h allutil.h linkage.h
reorder_loci.o: linkage_external.h makeped.h pedtree.h analysis.h
reorder_loci.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
reorder_loci.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
reorder_loci.o: write_beagle_ext.h write_structure_ext.h batch_input.h
reorder_loci.o: input_check.h hwe.h compress_ext.h annotated_ped_file_ext.h
reorder_loci.o: batch_input_ext.h create_summary_ext.h error_messages_ext.h
reorder_loci.o: fcmap_ext.h genetic_utils_ext.h grow_string_ext.h
reorder_loci.o: linkage_ext.h output_routines_ext.h output_file_names_ext.h
reorder_loci.o: user_input_ext.h utils_ext.h write_files_ext.h
write_vitesse.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_vitesse.o: linkage_external.h makeped.h pedtree.h analysis.h
write_vitesse.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_vitesse.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_vitesse.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_vitesse.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
write_vitesse.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_vitesse.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_vitesse.o: user_input_ext.h utils_ext.h write_files_ext.h
hwe_user_input.o: common.h typedefs.h errorno.h allutil.h linkage.h
hwe_user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
hwe_user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
hwe_user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
hwe_user_input.o: write_beagle_ext.h write_structure_ext.h batch_input.h
hwe_user_input.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
hwe_user_input.o: error_messages_ext.h fcmap_ext.h gen_user_input_ext.h
hwe_user_input.o: genetic_utils_ext.h grow_string_ext.h homozygosity_ext.h
hwe_user_input.o: hwe_ext.h linkage_ext.h output_file_names_ext.h
hwe_user_input.o: output_routines_ext.h select_individuals_ext.h
hwe_user_input.o: user_input_ext.h utils_ext.h write_mendel7_files_ext.h
hwe_user_input.o: write_mfiles_ext.h
scripts.o: common.h typedefs.h errorno.h allutil.h linkage.h
scripts.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
scripts.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
scripts.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
scripts.o: write_structure_ext.h batch_input.h input_check.h hwe.h
scripts.o: compress_ext.h error_messages_ext.h fcmap_ext.h
scripts.o: genetic_utils_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
scripts.o: output_file_names_ext.h user_input_ext.h utils_ext.h
scripts.o: write_files_ext.h write_ghfiles_ext.h
makeped1.o: common.h typedefs.h errorno.h allutil.h linkage.h
makeped1.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
makeped1.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
makeped1.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
makeped1.o: write_structure_ext.h batch_input.h input_check.h hwe.h
makeped1.o: compress_ext.h error_messages_ext.h fcmap_ext.h grow_string_ext.h
makeped1.o: linkage_ext.h read_files_ext.h utils_ext.h
omit_ped.o: common.h typedefs.h errorno.h allutil.h linkage.h
omit_ped.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
omit_ped.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
omit_ped.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
omit_ped.o: write_structure_ext.h batch_input.h input_check.h hwe.h
omit_ped.o: compress_ext.h error_messages_ext.h user_input_ext.h utils_ext.h
slink.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
slink.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
slink.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
slink.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
slink.o: batch_input.h input_check.h hwe.h compress_ext.h
slink.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h linkage_ext.h
slink.o: omit_ped_ext.h output_file_names_ext.h pedtree_ext.h
slink.o: user_input_ext.h utils_ext.h write_ghfiles_ext.h
R_output.o: common.h typedefs.h errorno.h allutil.h linkage.h
R_output.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
R_output.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
R_output.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
R_output.o: write_structure_ext.h batch_input.h input_check.h hwe.h
R_output.o: compress_ext.h R_output.h batch_input_ext.h error_messages_ext.h
R_output.o: fcmap_ext.h grow_string_ext.h output_file_names_ext.h utils_ext.h
R_output.o: plink_ext.h user_input_ext.h
write_premakeped.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_premakeped.o: linkage_external.h makeped.h pedtree.h analysis.h
write_premakeped.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_premakeped.o: plink_core_ext.h write_eigenstrat_ext.h
write_premakeped.o: write_pangaea_ext.h write_beagle_ext.h
write_premakeped.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_premakeped.o: compress_ext.h R_output.h R_output_ext.h
write_premakeped.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
write_premakeped.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_premakeped.o: output_file_names_ext.h output_routines_ext.h
write_premakeped.o: user_input_ext.h utils_ext.h write_files_ext.h
write_premakeped.o: write_ghfiles_ext.h write_mfiles_ext.h
create_summary.o: common.h typedefs.h errorno.h allutil.h linkage.h
create_summary.o: linkage_external.h makeped.h pedtree.h analysis.h
create_summary.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
create_summary.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
create_summary.o: write_beagle_ext.h write_structure_ext.h batch_input.h
create_summary.o: input_check.h hwe.h compress_ext.h liable.h
create_summary.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
create_summary.o: grow_string_ext.h input_check_ext.h liability_dist_ext.h
create_summary.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
create_summary.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
create_summary.o: reorder_loci_ext.h user_input_ext.h utils_ext.h
write_prest.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_prest.o: linkage_external.h makeped.h pedtree.h analysis.h
write_prest.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_prest.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_prest.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_prest.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_prest.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_prest.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
write_prest.o: utils_ext.h write_files_ext.h write_ghfiles_ext.h
write_pap.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_pap.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
write_pap.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
write_pap.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
write_pap.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_pap.o: compress_ext.h error_messages_ext.h fcmap_ext.h linkage_ext.h
write_pap.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_pap.o: utils_ext.h
batch_input.o: common.h typedefs.h errorno.h allutil.h linkage.h
batch_input.o: linkage_external.h makeped.h pedtree.h analysis.h
batch_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
batch_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
batch_input.o: write_beagle_ext.h write_structure_ext.h batch_input.h
batch_input.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
batch_input.o: error_messages_ext.h genetic_utils_ext.h reorder_loci_ext.h
batch_input.o: user_input_ext.h plink_ext.h utils_ext.h
makenucs.o: common.h typedefs.h errorno.h allutil.h linkage.h
makenucs.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
makenucs.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
makenucs.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
makenucs.o: write_structure_ext.h batch_input.h input_check.h hwe.h
makenucs.o: compress_ext.h error_messages_ext.h fcmap_ext.h linkage_ext.h
makenucs.o: makenucs_ext.h omit_ped_ext.h output_file_names_ext.h
makenucs.o: output_routines_ext.h pedtree_ext.h user_input_ext.h utils_ext.h
makenucs.o: write_files_ext.h write_ghfiles_ext.h append_locus_array_ext.h
append_locus_array.o: common.h typedefs.h errorno.h allutil.h linkage.h
append_locus_array.o: linkage_external.h makeped.h pedtree.h analysis.h
append_locus_array.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
append_locus_array.o: plink_core_ext.h write_eigenstrat_ext.h
append_locus_array.o: write_pangaea_ext.h write_beagle_ext.h
append_locus_array.o: write_structure_ext.h batch_input.h input_check.h hwe.h
append_locus_array.o: compress_ext.h error_messages_ext.h linkage_ext.h
error_sim.o: common.h typedefs.h errorno.h allutil.h linkage.h
error_sim.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
error_sim.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
error_sim.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
error_sim.o: write_structure_ext.h batch_input.h input_check.h hwe.h
error_sim.o: compress_ext.h error_sim.h batch_input_ext.h
error_sim.o: error_messages_ext.h genetic_utils_ext.h grow_string_ext.h
error_sim.o: list_ext.h output_file_names_ext.h output_routines_ext.h
error_sim.o: reorder_loci_ext.h utils_ext.h
write_loki.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_loki.o: linkage_external.h makeped.h pedtree.h analysis.h
write_loki.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_loki.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_loki.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_loki.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_loki.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_loki.o: output_file_names_ext.h user_input_ext.h utils_ext.h
write_loki.o: write_premakeped_ext.h
error_messages.o: common.h typedefs.h errorno.h allutil.h linkage.h
error_messages.o: linkage_external.h makeped.h pedtree.h analysis.h
error_messages.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
error_messages.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
error_messages.o: write_beagle_ext.h write_structure_ext.h batch_input.h
error_messages.o: input_check.h hwe.h compress_ext.h utils_ext.h
annotated_ped_file.o: common.h typedefs.h errorno.h allutil.h linkage.h
annotated_ped_file.o: linkage_external.h makeped.h pedtree.h analysis.h
annotated_ped_file.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
annotated_ped_file.o: plink_core_ext.h write_eigenstrat_ext.h
annotated_ped_file.o: write_pangaea_ext.h write_beagle_ext.h
annotated_ped_file.o: write_structure_ext.h batch_input.h input_check.h hwe.h
annotated_ped_file.o: compress_ext.h mrecode.h annotated_ped_file.h
annotated_ped_file.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h
annotated_ped_file.o: linkage_ext.h list_ext.h makeped_ext.h
annotated_ped_file.o: marker_lookup_ext.h mrecode_ext.h omit_ped_ext.h
annotated_ped_file.o: phe_lookup_ext.h read_files_ext.h reorder_loci_ext.h
annotated_ped_file.o: select_individuals_ext.h user_input_ext.h utils_ext.h
annotated_ped_file.o: write_files_ext.h plink_ext.h
write_SUP_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_SUP_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_SUP_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_SUP_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_SUP_files.o: write_pangaea_ext.h write_beagle_ext.h
write_SUP_files.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_SUP_files.o: compress_ext.h error_messages_ext.h fcmap_ext.h
write_SUP_files.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_SUP_files.o: user_input_ext.h utils_ext.h write_files_ext.h
write_mendel7_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_mendel7_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_mendel7_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_mendel7_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_mendel7_files.o: write_pangaea_ext.h write_beagle_ext.h
write_mendel7_files.o: write_structure_ext.h batch_input.h input_check.h
write_mendel7_files.o: hwe.h compress_ext.h write_mfiles.h
write_mendel7_files.o: create_summary_ext.h error_messages_ext.h
write_mendel7_files.o: genetic_utils_ext.h utils_ext.h user_input_ext.h
write_mendel7_files.o: write_mendel7_files_ext.h
write_plink_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_plink_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_plink_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_plink_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_plink_files.o: write_pangaea_ext.h write_beagle_ext.h
write_plink_files.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_plink_files.o: compress_ext.h loop.h file_ops.h error_messages_ext.h
write_plink_files.o: utils_ext.h entry.h create_summary_ext.h
write_plink_files.o: output_routines_ext.h sh_util.h linkage_ext.h
write_plink_files.o: fcmap_ext.h omit_ped_ext.h output_file_names_ext.h
write_plink_files.o: user_input_ext.h
cranefoot.o: common.h typedefs.h errorno.h allutil.h linkage.h
cranefoot.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cranefoot.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cranefoot.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cranefoot.o: write_structure_ext.h batch_input.h input_check.h hwe.h
cranefoot.o: compress_ext.h create_summary_ext.h error_messages_ext.h
cranefoot.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h output_file_names_ext.h
cranefoot.o: user_input_ext.h utils_ext.h
write_mega2annot.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_mega2annot.o: linkage_external.h makeped.h pedtree.h analysis.h
write_mega2annot.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_mega2annot.o: plink_core_ext.h write_eigenstrat_ext.h
write_mega2annot.o: write_pangaea_ext.h write_beagle_ext.h
write_mega2annot.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_mega2annot.o: compress_ext.h annotated_ped_file_ext.h
write_mega2annot.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_mega2annot.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_mega2annot.o: output_routines_ext.h user_input_ext.h utils_ext.h
grow_string.o: common.h typedefs.h errorno.h allutil.h linkage.h
grow_string.o: linkage_external.h makeped.h pedtree.h analysis.h
grow_string.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
grow_string.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
grow_string.o: write_beagle_ext.h write_structure_ext.h batch_input.h
grow_string.o: input_check.h hwe.h compress_ext.h utils_ext.h
write_IQLS.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_IQLS.o: linkage_external.h makeped.h pedtree.h analysis.h
write_IQLS.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_IQLS.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_IQLS.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_IQLS.o: input_check.h hwe.h compress_ext.h create_summary_ext.h
write_IQLS.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_IQLS.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
write_IQLS.o: user_input_ext.h utils_ext.h loop.h file_ops.h entry.h
write_IQLS.o: sh_util.h
net.o: common.h error_messages_ext.h
phe_lookup.o: plink_ext.h phe_lookup_ext.h
analysis.o: common.h typedefs.h errorno.h allutil.h linkage.h
analysis.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
analysis.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
analysis.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
analysis.o: write_structure_ext.h batch_input.h input_check.h hwe.h
analysis.o: compress_ext.h utils_ext.h error_messages_ext.h fcmap_ext.h
analysis.o: output_file_names_ext.h aspex_ext.h cranefoot_ext.h
analysis.o: create_summary_ext.h hwe_user_input_ext.h makenucs_ext.h
analysis.o: mega2annot_ext.h scripts_ext.h slink_ext.h splink_ext.h
analysis.o: write_IQLS_ext.h write_SUP_files_ext.h write_apm_ext.h
analysis.o: write_files_ext.h write_ghfiles_ext.h write_loki_ext.h
analysis.o: write_mfiles_ext.h write_pap_ext.h write_premakeped_ext.h
analysis.o: write_prest_ext.h write_sage_files_ext.h
analysis.o: write_simulate_files_ext.h write_solar_files_ext.h
analysis.o: write_vitesse_ext.h
write_fbat.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_fbat.o: linkage_external.h makeped.h pedtree.h analysis.h
write_fbat.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_fbat.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_fbat.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_fbat.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_fbat.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_fbat.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
write_fbat.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_fbat.o: fbat.R.h
entry.o: file_ops.h common.h typedefs.h errorno.h allutil.h linkage.h
entry.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
entry.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
entry.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
entry.o: write_structure_ext.h batch_input.h input_check.h hwe.h
entry.o: compress_ext.h error_messages_ext.h utils_ext.h user_input_ext.h
entry.o: entry.h create_summary_ext.h output_routines_ext.h
loop.o: loop.h file_ops.h common.h typedefs.h errorno.h allutil.h linkage.h
loop.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
loop.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
loop.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
loop.o: write_structure_ext.h batch_input.h input_check.h hwe.h
loop.o: compress_ext.h error_messages_ext.h utils_ext.h entry.h
loop.o: create_summary_ext.h output_routines_ext.h sh_util.h linkage_ext.h
write_pangaea.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_pangaea.o: linkage_external.h makeped.h pedtree.h analysis.h
write_pangaea.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_pangaea.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_pangaea.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_pangaea.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_pangaea.o: error_messages_ext.h utils_ext.h entry.h
write_pangaea.o: create_summary_ext.h output_routines_ext.h sh_util.h
write_pangaea.o: linkage_ext.h fcmap_ext.h genetic_utils_ext.h omit_ped_ext.h
write_pangaea.o: output_file_names_ext.h user_input_ext.h
write_beagle.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_beagle.o: linkage_external.h makeped.h pedtree.h analysis.h
write_beagle.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_beagle.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_beagle.o: write_beagle_ext.h write_structure_ext.h batch_input.h
write_beagle.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_beagle.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_beagle.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
write_beagle.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_beagle.o: genetic_utils_ext.h
plink_core_files.o: common.h typedefs.h errorno.h allutil.h linkage.h
plink_core_files.o: linkage_external.h makeped.h pedtree.h analysis.h
plink_core_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
plink_core_files.o: plink_core_ext.h write_eigenstrat_ext.h
plink_core_files.o: write_pangaea_ext.h write_beagle_ext.h
plink_core_files.o: write_structure_ext.h batch_input.h input_check.h hwe.h
plink_core_files.o: compress_ext.h loop.h file_ops.h error_messages_ext.h
plink_core_files.o: utils_ext.h entry.h create_summary_ext.h
plink_core_files.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
plink_core_files.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_eigenstrat.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_eigenstrat.o: linkage_external.h makeped.h pedtree.h analysis.h
write_eigenstrat.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_eigenstrat.o: plink_core_ext.h write_eigenstrat_ext.h
write_eigenstrat.o: write_pangaea_ext.h write_beagle_ext.h
write_eigenstrat.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_eigenstrat.o: compress_ext.h loop.h file_ops.h error_messages_ext.h
write_eigenstrat.o: utils_ext.h entry.h create_summary_ext.h
write_eigenstrat.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
write_eigenstrat.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_structure.o: common.h typedefs.h errorno.h allutil.h linkage.h
write_structure.o: linkage_external.h makeped.h pedtree.h analysis.h
write_structure.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_structure.o: plink_core_ext.h write_eigenstrat_ext.h
write_structure.o: write_pangaea_ext.h write_beagle_ext.h
write_structure.o: write_structure_ext.h batch_input.h input_check.h hwe.h
write_structure.o: compress_ext.h loop.h file_ops.h error_messages_ext.h
write_structure.o: utils_ext.h entry.h create_summary_ext.h
write_structure.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
write_structure.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
batch.o: common.h typedefs.h errorno.h allutil.h linkage.h linkage_external.h
batch.o: makeped.h pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
batch.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
batch.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
batch.o: batch_input.h input_check.h hwe.h compress_ext.h batch.h
batch.o: reorder_loci_ext.h batch_input_ext.h error_messages_ext.h
batch.o: genetic_utils_ext.h user_input_ext.h plink_ext.h utils_ext.h
batch.o: cw_routines_ext.h
compress.o: common.h typedefs.h errorno.h allutil.h linkage.h
compress.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
compress.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
compress.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
compress.o: write_structure_ext.h batch_input.h input_check.h hwe.h
compress.o: compress_ext.h error_messages_ext.h utils_ext.h
