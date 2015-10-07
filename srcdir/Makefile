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
    else ifeq "$(findstring MSYS_NT,$(OSTYPE))" "MSYS_NT"
	OSTYPE = msys2
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
	grow_string.cpp     marker_lookup.cpp    write_IQLS.cpp \
	net.cpp             phe_lookup.cpp       analysis.cpp \
	write_fbat.cpp      entry.cpp            loop.cpp \
	write_pangaea.cpp   write_beagle.cpp     plink_core_files.cpp \
	write_eigenstrat.cpp write_structure.cpp batch.cpp \
	compress.cpp        write_pseq_files.cpp menu_value_missing.cpp \
	str_utils.cpp       read_impute.cpp	 read_impute_bgen.cpp \
	read_impute_bgen2.cpp

OBJS =  $(SRCS:.cpp=.o)

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

user_input.o read_impute_bgen2.o: IFLAGS += -Ibgen/genfile/include $(stdC11)

write_IQLS.o write_beagle.o write_fbat.o write_loki.o write_mfiles.o: DFLAGS += -DRUNSHELL_SETUP
write_pangaea.o write_plink_files.o write_premakeped.o: DFLAGS += -DRUNSHELL_SETUP
write_structure.o: DFLAGS += -DRUNSHELL_SETUP

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
aspex.o annotated_ped_file.o batch_input.o genetic_utils.o makeped1.o mega2.o: DFLAGS += $(DsFLAGS)
menu_value_missing.o mrecode.o: DFLAGS += $(DsFLAGS)
output_routines.o read_files.o reorder_loci.o user_input.o utils.o write_files.o: DFLAGS += $(DsFLAGS)
write_plink_files.o: DFLAGS += $(DsFLAGS)

# -DHIDEPATH
aspex.o cranefoot.o mega2.o R_output.o write_SUP_files.o write_ghfiles.o write_premakeped.o: DFLAGS += $(DpFLAGS)

# -DHIDEFILE
annotated_ped_file.o makeped1.o mega2.o read_files.o utils.o: DFLAGS += $(DfFLAGS)
write_IQLS.o write_beagle.o write_eigenstrat.o write_fbat.o: DFLAGS += $(DfFLAGS)
write_pangaea.o write_plink_files.o write_pseq_files.o write_structure.o: DFLAGS += $(DfFLAGS)

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

# export CFLAGS so that the called Makefile can find it...
export FLAGS
export OSTYPE
export CFLAGS
export CXX
export CC
export AR
export ARFLAGS
export AROUT

LINK ?= $(CXX)
LTARGET ?= -o $(EXEC)
AR ?= ar
ARFLAGS ?= rcs
AROUT ?=

ZLIBS ?= libz.a
CLIBS += libvcftools.a libbgen.a $(ZLIBS)

stdC11  ?= -std=c++0x
WINO    ?=

.c.o:
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CC) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(IFLAGS) $*.c $(WINO)

.cpp.o:	
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CXX) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(IFLAGS) $*.cpp $(WINO)

all: $(MODULE)_$(OSTYPE)

dbg:
	$(MAKE) $(MAKEFLAGS) TEST=2 HIDE=1 OPT=-g all

tod:
	$(MAKE) $(MAKEFLAGS) TEST=2 HIDE=1 OPT=-g FLAGS=-DTOD all

gcc_info.o: $(OBJS)

SUBDIRS = vcftools bgen zlib-1.2.8
.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) $(MFLAGS) -C $@
vcftools: zlib-1.2.8
libvcftools.a: vcftools
libz.a: zlib-1.2.8

libbgen.a: bgen
bgen: zlib-1.2.8

$(MODULE)_$(OSTYPE): $(SUBDIRS) $(OBJS) $(HDRS) gcc_info.o libvcftools.a libbgen.a libz.a
	$(RM) $(MODULE)
	$(LINK) $(LFLAGS) $(LTARGET) $(OBJS) $(FLAGS) gcc_info.o $(CLIBS)

ifeq ($(OSTYPE),mingw)
	cp $(MODULE)_$(OSTYPE).exe $(MODULE).exe
else ifeq ($(OSTYPE),windows)
	copy $(MODULE)_$(OSTYPE).exe $(MODULE).exe
else ifneq ($(LN),)
	ln -s $(MODULE)_$(OSTYPE) $(MODULE)
endif

clean:
	$(RM) $(OBJS) gcc_info.o $(EXEC) $(MODULE) $(EXEC).exe $(MODULE).exe
#	$(foreach dir,$(SUBDIRS),$(MAKE) -C $(dir) RM="$(RM)" clean &&) echo on
	$(MAKE) $(MFLAGS) -C vcftools   RM="$(RM)" clean
	$(MAKE) $(MFLAGS) -C bgen       RM="$(RM)" clean
	$(MAKE) $(MFLAGS) -C zlib-1.2.8 RM="$(RM)" clean

distclean:
	$(RM) $(OBJS) gcc_info.o $(EXEC) $(MODULE) $(EXEC).exe $(MODULE).exe
#	$(foreach dir,$(SUBDIRS),$(MAKE) -C $(dir) RM="$(RM)" distclean &&) echo on
	$(MAKE) $(MFLAGS) -C vcftools   RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C bgen       RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C zlib-1.2.8 RM="$(RM)" distclean

clean1:
	$(RM) *~

srcdist:
	for i in $(SRCS) $(HDRS); do $(REMOVECVS1) $$i > ../srcdir/$$i; done
	for i in `ls Makedefs*`; do $(REMOVECVS) $$i > ../srcdir/$$i.template; done
	$(REMOVECVS) Makefile > ../srcdir/Makefile

depend:
	makedepend -Y. -- $(CFLAGS) -- $(SRCS)
	$(MAKE) -C vcftools depend
	$(MAKE) -C bgen     depend

# DO NOT DELETE

input_check.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
input_check.o: linkage_external.h makeped.h pedtree.h analysis.h
input_check.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
input_check.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
input_check.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
input_check.o: batch_input.h input_check.h hwe.h compress_ext.h tod.hh
input_check.o: check_ext.h create_summary_ext.h error_messages_ext.h
input_check.o: fcmap_ext.h genetic_utils_ext.h grow_string_ext.h
input_check.o: user_input_ext.h input_ops.hh str_utils.hh input_check_ext.h
input_check.o: pedtree_ext.h utils_ext.h write_pap_ext.h batch_input_ext.h
select_index.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
select_index.o: linkage_external.h makeped.h pedtree.h analysis.h
select_index.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
select_index.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
select_index.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
select_index.o: batch_input.h input_check.h hwe.h compress_ext.h
select_index.o: random_choose_ext.h
aspex.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
aspex.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
aspex.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
aspex.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
aspex.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
aspex.o: hwe.h compress_ext.h batch_input_ext.h create_summary_ext.h
aspex.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
aspex.o: grow_string_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
aspex.o: output_file_names_ext.h output_routines_ext.h scripts_ext.h
aspex.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
aspex.o: write_files_ext.h
liability_dist.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
liability_dist.o: linkage_external.h makeped.h pedtree.h analysis.h
liability_dist.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
liability_dist.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
liability_dist.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
liability_dist.o: batch_input.h input_check.h hwe.h compress_ext.h liable.h
liability_dist.o: error_messages_ext.h genetic_utils_ext.h
liability_dist.o: liability_dist_ext.h linkage_ext.h omit_ped_ext.h
liability_dist.o: output_file_names_ext.h utils_ext.h
select_individuals.o: common.h types.hh typedefs.h errorno.h allutil.h
select_individuals.o: linkage.h linkage_external.h makeped.h pedtree.h
select_individuals.o: analysis.h analysis_ext.h write_fbat_ext.h
select_individuals.o: write_plink_ext.h plink_core_ext.h
select_individuals.o: write_eigenstrat_ext.h write_pangaea_ext.h
select_individuals.o: write_beagle_ext.h write_structure_ext.h
select_individuals.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
select_individuals.o: compress_ext.h mrecode.h batch_input_ext.h
select_individuals.o: error_messages_ext.h fcmap_ext.h utils_ext.h
linkage.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
linkage.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
linkage.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
linkage.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
linkage.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
linkage.o: hwe.h compress_ext.h error_messages_ext.h list_ext.h makeped_ext.h
linkage.o: pedtree_ext.h utils_ext.h user_input_ext.h input_ops.hh
linkage.o: str_utils.hh append_locus_array_ext.h
output_file_names.o: common.h types.hh typedefs.h errorno.h allutil.h
output_file_names.o: linkage.h linkage_external.h makeped.h pedtree.h
output_file_names.o: analysis.h analysis_ext.h write_fbat_ext.h
output_file_names.o: write_plink_ext.h plink_core_ext.h
output_file_names.o: write_eigenstrat_ext.h write_pangaea_ext.h
output_file_names.o: write_beagle_ext.h write_structure_ext.h
output_file_names.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
output_file_names.o: compress_ext.h batch_input_ext.h error_messages_ext.h
output_file_names.o: fcmap_ext.h genetic_utils_ext.h output_file_names_ext.h
output_file_names.o: utils_ext.h
cal_const.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_const.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cal_const.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cal_const.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cal_const.o: write_structure_ext.h write_pseq_ext.h batch_input.h
cal_const.o: input_check.h hwe.h compress_ext.h hwe_ext.h
list.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
list.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
list.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
list.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
list.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
list.o: hwe.h compress_ext.h error_messages_ext.h
test_switch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
test_switch.o: linkage_external.h makeped.h pedtree.h analysis.h
test_switch.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
test_switch.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
test_switch.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
test_switch.o: batch_input.h input_check.h hwe.h compress_ext.h
cal_n.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_n.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cal_n.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cal_n.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cal_n.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
cal_n.o: hwe.h compress_ext.h
ln_p_value.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
ln_p_value.o: linkage_external.h makeped.h pedtree.h analysis.h
ln_p_value.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
ln_p_value.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
ln_p_value.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
ln_p_value.o: batch_input.h input_check.h hwe.h compress_ext.h hwe_ext.h
user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
user_input.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
user_input.o: batch_input.h input_check.h hwe.h compress_ext.h input_ops.hh
user_input.o: str_utils.hh input.hh read_impute.hh annotated_ped_file.h
user_input.o: read_impute_bgen.hh zlib.h zconf.h read_impute_bgen2.hh
user_input.o: batch_input_ext.h cw_routines_ext.h error_messages_ext.h
user_input.o: fcmap_ext.h genetic_utils_ext.h reorder_loci_ext.h
user_input.o: grow_string_ext.h output_file_names_ext.h output_routines_ext.h
user_input.o: plink_ext.h utils_ext.h user_input_ext.h
user_input.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
cal_prob.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_prob.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cal_prob.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cal_prob.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cal_prob.o: write_structure_ext.h write_pseq_ext.h batch_input.h
cal_prob.o: input_check.h hwe.h compress_ext.h do_switch_ext.h
cal_prob.o: test_switch_ext.h utils_ext.h
utils.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
utils.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
utils.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
utils.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
utils.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
utils.o: hwe.h compress_ext.h batch_input_ext.h error_messages_ext.h
utils.o: fcmap_ext.h net_ext.h user_input_ext.h input_ops.hh str_utils.hh
utils.o: utils_ext.h version.h
check.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
check.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
check.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
check.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
check.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
check.o: hwe.h compress_ext.h tod.hh utils_ext.h error_messages_ext.h
check.o: grow_string_ext.h write_files_ext.h
mega2.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
mega2.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
mega2.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
mega2.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
mega2.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
mega2.o: hwe.h compress_ext.h tod.hh R_output.h annotated_ped_file.h
mega2.o: annotated_ped_file_ext.h aspex_ext.h batch_input_ext.h
mega2.o: cranefoot_ext.h create_summary_ext.h error_messages_ext.h
mega2.o: error_sim_ext.h fcmap_ext.h genetic_utils_ext.h hwe_user_input_ext.h
mega2.o: input_check_ext.h makenucs_ext.h makeped_ext.h marker_lookup_ext.h
mega2.o: mega2annot_ext.h menu_value_missing_ext.h output_file_names_ext.h
mega2.o: output_routines_ext.h pedtree_ext.h plink_ext.h read_files_ext.h
mega2.o: reorder_loci_ext.h scripts_ext.h slink_ext.h splink_ext.h
mega2.o: user_input_ext.h input_ops.hh str_utils.hh
mega2.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h utils_ext.h
mega2.o: write_IQLS_ext.h write_SUP_files_ext.h write_apm_ext.h
mega2.o: write_files_ext.h write_ghfiles_ext.h write_loki_ext.h
mega2.o: write_mfiles_ext.h write_pap_ext.h write_premakeped_ext.h
mega2.o: write_prest_ext.h write_sage_files_ext.h write_simulate_files_ext.h
mega2.o: write_solar_files_ext.h write_vitesse_ext.h version.h
write_apm.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_apm.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
write_apm.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
write_apm.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
write_apm.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_apm.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_apm.o: fcmap_ext.h genetic_utils_ext.h input_check_ext.h linkage_ext.h
write_apm.o: omit_ped_ext.h output_file_names_ext.h pedtree_ext.h
write_apm.o: reorder_loci_ext.h user_input_ext.h input_ops.hh str_utils.hh
write_apm.o: utils_ext.h
output_routines.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
output_routines.o: linkage_external.h makeped.h pedtree.h analysis.h
output_routines.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
output_routines.o: plink_core_ext.h write_eigenstrat_ext.h
output_routines.o: write_pangaea_ext.h write_beagle_ext.h
output_routines.o: write_structure_ext.h write_pseq_ext.h batch_input.h
output_routines.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
output_routines.o: grow_string_ext.h output_routines_ext.h utils_ext.h
write_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_files.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_files.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_files.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_files.o: grow_string_ext.h linkage_ext.h makeped_ext.h omit_ped_ext.h
write_files.o: output_file_names_ext.h reorder_loci_ext.h user_input_ext.h
write_files.o: input_ops.hh str_utils.hh utils_ext.h write_ghfiles_ext.h
write_files.o: plink_ext.h
do_switch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
do_switch.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
do_switch.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
do_switch.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
do_switch.o: write_structure_ext.h write_pseq_ext.h batch_input.h
do_switch.o: input_check.h hwe.h compress_ext.h
pedtree.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
pedtree.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
pedtree.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
pedtree.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
pedtree.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
pedtree.o: hwe.h compress_ext.h error_messages_ext.h linkage_ext.h
write_ghfiles.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_ghfiles.o: linkage_external.h makeped.h pedtree.h analysis.h
write_ghfiles.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_ghfiles.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_ghfiles.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_ghfiles.o: batch_input.h input_check.h hwe.h compress_ext.h R_output.h
write_ghfiles.o: R_output_ext.h error_messages_ext.h fcmap_ext.h
write_ghfiles.o: genetic_utils_ext.h linkage_ext.h makenucs_ext.h
write_ghfiles.o: omit_ped_ext.h output_file_names_ext.h output_routines_ext.h
write_ghfiles.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
write_ghfiles.o: write_files_ext.h
fcmap.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
fcmap.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
fcmap.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
fcmap.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
fcmap.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
fcmap.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h
print_data.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
print_data.o: linkage_external.h makeped.h pedtree.h analysis.h
print_data.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
print_data.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
print_data.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
print_data.o: batch_input.h input_check.h hwe.h compress_ext.h
write_liability.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_liability.o: linkage_external.h makeped.h pedtree.h analysis.h
write_liability.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_liability.o: plink_core_ext.h write_eigenstrat_ext.h
write_liability.o: write_pangaea_ext.h write_beagle_ext.h
write_liability.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_liability.o: input_check.h hwe.h compress_ext.h liable.h
write_liability.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_liability.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_liability.o: str_utils.hh
random_choose.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
random_choose.o: linkage_external.h makeped.h pedtree.h analysis.h
random_choose.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
random_choose.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
random_choose.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
random_choose.o: batch_input.h input_check.h hwe.h compress_ext.h utils_ext.h
write_mfiles.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_mfiles.o: linkage_external.h makeped.h pedtree.h analysis.h
write_mfiles.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_mfiles.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_mfiles.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_mfiles.o: batch_input.h input_check.h hwe.h compress_ext.h R_output.h
write_mfiles.o: write_mfiles.h batch_input_ext.h R_output_ext.h
write_mfiles.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
write_mfiles.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_mfiles.o: output_file_names_ext.h output_routines_ext.h
write_mfiles.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
write_mfiles.o: write_mendel7_files_ext.h write_files_ext.h
write_mfiles.o: write_mfiles_ext.h
mrecode.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
mrecode.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
mrecode.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
mrecode.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
mrecode.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
mrecode.o: hwe.h compress_ext.h tod.hh mrecode.h mrecode_ext.h
mrecode.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
mrecode.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h linkage_ext.h
mrecode.o: makeped_ext.h marker_lookup_ext.h omit_ped_ext.h read_files_ext.h
mrecode.o: select_individuals_ext.h utils_ext.h write_files_ext.h
gen_user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
gen_user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
gen_user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
gen_user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
gen_user_input.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
gen_user_input.o: batch_input.h input_check.h hwe.h compress_ext.h
gen_user_input.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
gen_user_input.o: output_file_names_ext.h user_input_ext.h input_ops.hh
gen_user_input.o: str_utils.hh utils_ext.h
write_sage_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_sage_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_sage_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_sage_files.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
write_sage_files.o: write_pangaea_ext.h write_beagle_ext.h
write_sage_files.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_sage_files.o: input_check.h hwe.h compress_ext.h create_summary_ext.h
write_sage_files.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_sage_files.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
write_sage_files.o: output_file_names_ext.h output_routines_ext.h
write_sage_files.o: read_files_ext.h user_input_ext.h input_ops.hh
write_sage_files.o: str_utils.hh utils_ext.h write_files_ext.h
splink.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
splink.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
splink.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
splink.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
splink.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
splink.o: hwe.h compress_ext.h create_summary_ext.h error_messages_ext.h
splink.o: fcmap_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
splink.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
splink.o: input_ops.hh str_utils.hh utils_ext.h write_files_ext.h
genetic_utils.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
genetic_utils.o: linkage_external.h makeped.h pedtree.h analysis.h
genetic_utils.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
genetic_utils.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
genetic_utils.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
genetic_utils.o: batch_input.h input_check.h hwe.h compress_ext.h
genetic_utils.o: error_messages_ext.h grow_string_ext.h utils_ext.h
genetic_utils.o: user_input_ext.h input_ops.hh str_utils.hh
read_data.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
read_data.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
read_data.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
read_data.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
read_data.o: write_structure_ext.h write_pseq_ext.h batch_input.h
read_data.o: input_check.h hwe.h compress_ext.h
write_simulate_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_simulate_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_simulate_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_simulate_files.o: write_plink_ext.h plink_core_ext.h
write_simulate_files.o: write_eigenstrat_ext.h write_pangaea_ext.h
write_simulate_files.o: write_beagle_ext.h write_structure_ext.h
write_simulate_files.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
write_simulate_files.o: compress_ext.h error_messages_ext.h fcmap_ext.h
write_simulate_files.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_simulate_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_simulate_files.o: str_utils.hh utils_ext.h write_files_ext.h
write_simulate_files.o: write_ghfiles_ext.h
homozygosity.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
homozygosity.o: linkage_external.h makeped.h pedtree.h analysis.h
homozygosity.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
homozygosity.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
homozygosity.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
homozygosity.o: batch_input.h input_check.h hwe.h compress_ext.h
homozygosity.o: genetic_utils_ext.h
read_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
read_files.o: linkage_external.h makeped.h pedtree.h analysis.h
read_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
read_files.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
read_files.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
read_files.o: batch_input.h input_check.h hwe.h compress_ext.h mrecode.h
read_files.o: annotated_ped_file.h error_sim.h error_messages_ext.h
read_files.o: fcmap_ext.h grow_string_ext.h linkage_ext.h list_ext.h
read_files.o: makeped_ext.h marker_lookup_ext.h mrecode_ext.h
read_files.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
read_files.o: omit_ped_ext.h plink_ext.h read_files_ext.h reorder_loci_ext.h
read_files.o: utils_ext.h user_input_ext.h input_ops.hh str_utils.hh
read_files.o: batch_input_ext.h write_files_ext.h
write_solar_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_solar_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_solar_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_solar_files.o: write_plink_ext.h plink_core_ext.h
write_solar_files.o: write_eigenstrat_ext.h write_pangaea_ext.h
write_solar_files.o: write_beagle_ext.h write_structure_ext.h
write_solar_files.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
write_solar_files.o: compress_ext.h error_messages_ext.h fcmap_ext.h
write_solar_files.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_solar_files.o: output_file_names_ext.h output_routines_ext.h
write_solar_files.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
write_solar_files.o: write_files_ext.h
hwe.o: hwe.h check_file_ext.h hwe_ext.h read_data_ext.h
reorder_loci.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
reorder_loci.o: linkage_external.h makeped.h pedtree.h analysis.h
reorder_loci.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
reorder_loci.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
reorder_loci.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
reorder_loci.o: batch_input.h input_check.h hwe.h compress_ext.h
reorder_loci.o: annotated_ped_file.h annotated_ped_file_ext.h
reorder_loci.o: batch_input_ext.h create_summary_ext.h error_messages_ext.h
reorder_loci.o: fcmap_ext.h genetic_utils_ext.h grow_string_ext.h
reorder_loci.o: linkage_ext.h output_routines_ext.h output_file_names_ext.h
reorder_loci.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
reorder_loci.o: write_files_ext.h
write_vitesse.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_vitesse.o: linkage_external.h makeped.h pedtree.h analysis.h
write_vitesse.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_vitesse.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_vitesse.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_vitesse.o: batch_input.h input_check.h hwe.h compress_ext.h
write_vitesse.o: batch_input_ext.h error_messages_ext.h fcmap_ext.h
write_vitesse.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_vitesse.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_vitesse.o: str_utils.hh utils_ext.h write_files_ext.h
hwe_user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
hwe_user_input.o: linkage_external.h makeped.h pedtree.h analysis.h
hwe_user_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
hwe_user_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
hwe_user_input.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
hwe_user_input.o: batch_input.h input_check.h hwe.h compress_ext.h mrecode.h
hwe_user_input.o: batch_input_ext.h error_messages_ext.h fcmap_ext.h
hwe_user_input.o: gen_user_input_ext.h genetic_utils_ext.h grow_string_ext.h
hwe_user_input.o: homozygosity_ext.h hwe_ext.h linkage_ext.h
hwe_user_input.o: output_file_names_ext.h output_routines_ext.h
hwe_user_input.o: select_individuals_ext.h user_input_ext.h input_ops.hh
hwe_user_input.o: str_utils.hh utils_ext.h write_mendel7_files_ext.h
hwe_user_input.o: write_mfiles_ext.h write_files_ext.h
scripts.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
scripts.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
scripts.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
scripts.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
scripts.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
scripts.o: hwe.h compress_ext.h error_messages_ext.h fcmap_ext.h
scripts.o: genetic_utils_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
scripts.o: output_file_names_ext.h user_input_ext.h input_ops.hh str_utils.hh
scripts.o: utils_ext.h write_files_ext.h write_ghfiles_ext.h
makeped1.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
makeped1.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
makeped1.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
makeped1.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
makeped1.o: write_structure_ext.h write_pseq_ext.h batch_input.h
makeped1.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
makeped1.o: fcmap_ext.h grow_string_ext.h linkage_ext.h read_files_ext.h
makeped1.o: utils_ext.h
omit_ped.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
omit_ped.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
omit_ped.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
omit_ped.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
omit_ped.o: write_structure_ext.h write_pseq_ext.h batch_input.h
omit_ped.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
omit_ped.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
slink.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
slink.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
slink.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
slink.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
slink.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
slink.o: hwe.h compress_ext.h error_messages_ext.h fcmap_ext.h
slink.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
slink.o: output_file_names_ext.h pedtree_ext.h user_input_ext.h input_ops.hh
slink.o: str_utils.hh utils_ext.h write_ghfiles_ext.h
R_output.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
R_output.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
R_output.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
R_output.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
R_output.o: write_structure_ext.h write_pseq_ext.h batch_input.h
R_output.o: input_check.h hwe.h compress_ext.h R_output.h batch_input_ext.h
R_output.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h
R_output.o: output_file_names_ext.h utils_ext.h plink_ext.h user_input_ext.h
R_output.o: input_ops.hh str_utils.hh
write_premakeped.o: common.h types.hh typedefs.h errorno.h allutil.h
write_premakeped.o: linkage.h linkage_external.h makeped.h pedtree.h
write_premakeped.o: analysis.h analysis_ext.h write_fbat_ext.h
write_premakeped.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
write_premakeped.o: write_pangaea_ext.h write_beagle_ext.h
write_premakeped.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_premakeped.o: input_check.h hwe.h compress_ext.h R_output.h
write_premakeped.o: R_output_ext.h create_summary_ext.h error_messages_ext.h
write_premakeped.o: fcmap_ext.h genetic_utils_ext.h linkage_ext.h
write_premakeped.o: omit_ped_ext.h output_file_names_ext.h
write_premakeped.o: output_routines_ext.h user_input_ext.h input_ops.hh
write_premakeped.o: str_utils.hh utils_ext.h write_files_ext.h
write_premakeped.o: write_ghfiles_ext.h write_mfiles_ext.h
create_summary.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
create_summary.o: linkage_external.h makeped.h pedtree.h analysis.h
create_summary.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
create_summary.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
create_summary.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
create_summary.o: batch_input.h input_check.h hwe.h compress_ext.h liable.h
create_summary.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
create_summary.o: grow_string_ext.h input_check_ext.h liability_dist_ext.h
create_summary.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
create_summary.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
create_summary.o: reorder_loci_ext.h user_input_ext.h input_ops.hh
create_summary.o: str_utils.hh utils_ext.h
write_prest.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_prest.o: linkage_external.h makeped.h pedtree.h analysis.h
write_prest.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_prest.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_prest.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_prest.o: batch_input.h input_check.h hwe.h compress_ext.h
write_prest.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_prest.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
write_prest.o: input_ops.hh str_utils.hh utils_ext.h write_files_ext.h
write_prest.o: write_ghfiles_ext.h
write_pap.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_pap.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
write_pap.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
write_pap.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
write_pap.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_pap.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_pap.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_pap.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
write_pap.o: write_files_ext.h
batch_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
batch_input.o: linkage_external.h makeped.h pedtree.h analysis.h
batch_input.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
batch_input.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
batch_input.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
batch_input.o: batch_input.h input_check.h hwe.h compress_ext.h input_ops.hh
batch_input.o: str_utils.hh batch_input_ext.h error_messages_ext.h
batch_input.o: menu_value_missing_ext.h genetic_utils_ext.h
batch_input.o: reorder_loci_ext.h user_input_ext.h plink_ext.h utils_ext.h
makenucs.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
makenucs.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
makenucs.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
makenucs.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
makenucs.o: write_structure_ext.h write_pseq_ext.h batch_input.h
makenucs.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
makenucs.o: fcmap_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
makenucs.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
makenucs.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
makenucs.o: write_files_ext.h write_ghfiles_ext.h append_locus_array_ext.h
append_locus_array.o: common.h types.hh typedefs.h errorno.h allutil.h
append_locus_array.o: linkage.h linkage_external.h makeped.h pedtree.h
append_locus_array.o: analysis.h analysis_ext.h write_fbat_ext.h
append_locus_array.o: write_plink_ext.h plink_core_ext.h
append_locus_array.o: write_eigenstrat_ext.h write_pangaea_ext.h
append_locus_array.o: write_beagle_ext.h write_structure_ext.h
append_locus_array.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
append_locus_array.o: compress_ext.h error_messages_ext.h linkage_ext.h
error_sim.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
error_sim.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
error_sim.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
error_sim.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
error_sim.o: write_structure_ext.h write_pseq_ext.h batch_input.h
error_sim.o: input_check.h hwe.h compress_ext.h error_sim.h batch_input_ext.h
error_sim.o: error_messages_ext.h genetic_utils_ext.h grow_string_ext.h
error_sim.o: list_ext.h output_file_names_ext.h output_routines_ext.h
error_sim.o: reorder_loci_ext.h utils_ext.h
write_loki.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_loki.o: linkage_external.h makeped.h pedtree.h analysis.h
write_loki.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_loki.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_loki.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_loki.o: batch_input.h input_check.h hwe.h compress_ext.h
write_loki.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_loki.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_loki.o: str_utils.hh utils_ext.h write_premakeped_ext.h
write_loki.o: write_files_ext.h
error_messages.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
error_messages.o: linkage_external.h makeped.h pedtree.h analysis.h
error_messages.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
error_messages.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
error_messages.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
error_messages.o: batch_input.h input_check.h hwe.h compress_ext.h
error_messages.o: utils_ext.h
annotated_ped_file.o: common.h types.hh typedefs.h errorno.h allutil.h
annotated_ped_file.o: linkage.h linkage_external.h makeped.h pedtree.h
annotated_ped_file.o: analysis.h analysis_ext.h write_fbat_ext.h
annotated_ped_file.o: write_plink_ext.h plink_core_ext.h
annotated_ped_file.o: write_eigenstrat_ext.h write_pangaea_ext.h
annotated_ped_file.o: write_beagle_ext.h write_structure_ext.h
annotated_ped_file.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
annotated_ped_file.o: compress_ext.h tod.hh mrecode.h annotated_ped_file.h
annotated_ped_file.o: annotated_ped_file_ext.h error_messages_ext.h
annotated_ped_file.o: fcmap_ext.h grow_string_ext.h linkage_ext.h list_ext.h
annotated_ped_file.o: makeped_ext.h marker_lookup_ext.h mrecode_ext.h
annotated_ped_file.o: vcftools/mega2_vcftools_interface.h
annotated_ped_file.o: annotated_ped_file.h omit_ped_ext.h phe_lookup_ext.h
annotated_ped_file.o: read_files_ext.h reorder_loci_ext.h
annotated_ped_file.o: select_individuals_ext.h user_input_ext.h input_ops.hh
annotated_ped_file.o: str_utils.hh utils_ext.h write_files_ext.h plink_ext.h
annotated_ped_file.o: read_impute.hh
write_SUP_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_SUP_files.o: linkage_external.h makeped.h pedtree.h analysis.h
write_SUP_files.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_SUP_files.o: plink_core_ext.h write_eigenstrat_ext.h
write_SUP_files.o: write_pangaea_ext.h write_beagle_ext.h
write_SUP_files.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_SUP_files.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_SUP_files.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_SUP_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_SUP_files.o: str_utils.hh utils_ext.h write_files_ext.h
write_mendel7_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_mendel7_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_mendel7_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_mendel7_files.o: write_plink_ext.h plink_core_ext.h
write_mendel7_files.o: write_eigenstrat_ext.h write_pangaea_ext.h
write_mendel7_files.o: write_beagle_ext.h write_structure_ext.h
write_mendel7_files.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
write_mendel7_files.o: compress_ext.h write_mfiles.h create_summary_ext.h
write_mendel7_files.o: error_messages_ext.h genetic_utils_ext.h
write_mendel7_files.o: output_routines_ext.h utils_ext.h user_input_ext.h
write_mendel7_files.o: input_ops.hh str_utils.hh write_files_ext.h
write_mendel7_files.o: write_mendel7_files_ext.h
write_plink_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_plink_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_plink_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_plink_files.o: write_plink_ext.h plink_core_ext.h
write_plink_files.o: write_eigenstrat_ext.h write_pangaea_ext.h
write_plink_files.o: write_beagle_ext.h write_structure_ext.h
write_plink_files.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
write_plink_files.o: compress_ext.h tod.hh loop.h file_ops.h zlib.h zconf.h
write_plink_files.o: error_messages_ext.h utils_ext.h entry.h
write_plink_files.o: create_summary_ext.h output_routines_ext.h sh_util.h
write_plink_files.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_plink_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_plink_files.o: str_utils.hh vcftools/mega2_vcftools_interface.h
write_plink_files.o: annotated_ped_file.h
cranefoot.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cranefoot.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
cranefoot.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
cranefoot.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
cranefoot.o: write_structure_ext.h write_pseq_ext.h batch_input.h
cranefoot.o: input_check.h hwe.h compress_ext.h create_summary_ext.h
cranefoot.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
cranefoot.o: output_file_names_ext.h user_input_ext.h input_ops.hh
cranefoot.o: str_utils.hh utils_ext.h write_files_ext.h
write_mega2annot.o: common.h types.hh typedefs.h errorno.h allutil.h
write_mega2annot.o: linkage.h linkage_external.h makeped.h pedtree.h
write_mega2annot.o: analysis.h analysis_ext.h write_fbat_ext.h
write_mega2annot.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
write_mega2annot.o: write_pangaea_ext.h write_beagle_ext.h
write_mega2annot.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_mega2annot.o: input_check.h hwe.h compress_ext.h annotated_ped_file.h
write_mega2annot.o: annotated_ped_file_ext.h error_messages_ext.h fcmap_ext.h
write_mega2annot.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_mega2annot.o: output_file_names_ext.h output_routines_ext.h
write_mega2annot.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
write_mega2annot.o: write_files_ext.h
grow_string.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
grow_string.o: linkage_external.h makeped.h pedtree.h analysis.h
grow_string.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
grow_string.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
grow_string.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
grow_string.o: batch_input.h input_check.h hwe.h compress_ext.h utils_ext.h
marker_lookup.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
marker_lookup.o: linkage_external.h makeped.h pedtree.h analysis.h
marker_lookup.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
marker_lookup.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
marker_lookup.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
marker_lookup.o: batch_input.h input_check.h hwe.h compress_ext.h
write_IQLS.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_IQLS.o: linkage_external.h makeped.h pedtree.h analysis.h
write_IQLS.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_IQLS.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_IQLS.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_IQLS.o: batch_input.h input_check.h hwe.h compress_ext.h
write_IQLS.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
write_IQLS.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_IQLS.o: output_routines_ext.h pedtree_ext.h user_input_ext.h
write_IQLS.o: input_ops.hh str_utils.hh utils_ext.h loop.h file_ops.h zlib.h
write_IQLS.o: zconf.h entry.h sh_util.h
net.o: common.h types.hh error_messages_ext.h
phe_lookup.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
phe_lookup.o: linkage_external.h makeped.h pedtree.h analysis.h
phe_lookup.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
phe_lookup.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
phe_lookup.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
phe_lookup.o: batch_input.h input_check.h hwe.h compress_ext.h utils_ext.h
phe_lookup.o: genetic_utils_ext.h error_messages_ext.h plink_ext.h
phe_lookup.o: phe_lookup_ext.h
analysis.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
analysis.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
analysis.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
analysis.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
analysis.o: write_structure_ext.h write_pseq_ext.h batch_input.h
analysis.o: input_check.h hwe.h compress_ext.h utils_ext.h
analysis.o: menu_value_missing.h error_messages_ext.h fcmap_ext.h
analysis.o: output_file_names_ext.h aspex_ext.h cranefoot_ext.h
analysis.o: create_summary_ext.h hwe_user_input_ext.h makenucs_ext.h
analysis.o: mega2annot_ext.h scripts_ext.h slink_ext.h splink_ext.h
analysis.o: write_IQLS_ext.h write_SUP_files_ext.h write_apm_ext.h
analysis.o: write_files_ext.h write_ghfiles_ext.h write_loki_ext.h
analysis.o: write_mfiles_ext.h write_pap_ext.h write_premakeped_ext.h
analysis.o: write_prest_ext.h write_sage_files_ext.h
analysis.o: write_simulate_files_ext.h write_solar_files_ext.h
analysis.o: write_vitesse_ext.h
write_fbat.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_fbat.o: linkage_external.h makeped.h pedtree.h analysis.h
write_fbat.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_fbat.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_fbat.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_fbat.o: batch_input.h input_check.h hwe.h compress_ext.h loop.h
write_fbat.o: file_ops.h zlib.h zconf.h error_messages_ext.h utils_ext.h
write_fbat.o: entry.h create_summary_ext.h output_routines_ext.h sh_util.h
write_fbat.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_fbat.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_fbat.o: str_utils.hh fbat.R.h
entry.o: file_ops.h zlib.h zconf.h common.h types.hh typedefs.h errorno.h
entry.o: allutil.h linkage.h linkage_external.h makeped.h pedtree.h
entry.o: analysis.h analysis_ext.h write_fbat_ext.h write_plink_ext.h
entry.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
entry.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
entry.o: batch_input.h input_check.h hwe.h compress_ext.h
entry.o: error_messages_ext.h utils_ext.h user_input_ext.h input_ops.hh
entry.o: str_utils.hh entry.h create_summary_ext.h output_routines_ext.h
entry.o: write_files_ext.h
loop.o: loop.h file_ops.h zlib.h zconf.h common.h types.hh typedefs.h
loop.o: errorno.h allutil.h linkage.h linkage_external.h makeped.h pedtree.h
loop.o: analysis.h analysis_ext.h write_fbat_ext.h write_plink_ext.h
loop.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
loop.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
loop.o: batch_input.h input_check.h hwe.h compress_ext.h error_messages_ext.h
loop.o: utils_ext.h entry.h create_summary_ext.h output_routines_ext.h
loop.o: sh_util.h linkage_ext.h tod.hh
write_pangaea.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_pangaea.o: linkage_external.h makeped.h pedtree.h analysis.h
write_pangaea.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_pangaea.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_pangaea.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_pangaea.o: batch_input.h input_check.h hwe.h compress_ext.h loop.h
write_pangaea.o: file_ops.h zlib.h zconf.h error_messages_ext.h utils_ext.h
write_pangaea.o: entry.h create_summary_ext.h output_routines_ext.h sh_util.h
write_pangaea.o: linkage_ext.h fcmap_ext.h genetic_utils_ext.h omit_ped_ext.h
write_pangaea.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_pangaea.o: str_utils.hh
write_beagle.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_beagle.o: linkage_external.h makeped.h pedtree.h analysis.h
write_beagle.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_beagle.o: plink_core_ext.h write_eigenstrat_ext.h write_pangaea_ext.h
write_beagle.o: write_beagle_ext.h write_structure_ext.h write_pseq_ext.h
write_beagle.o: batch_input.h input_check.h hwe.h compress_ext.h loop.h
write_beagle.o: file_ops.h zlib.h zconf.h error_messages_ext.h utils_ext.h
write_beagle.o: entry.h create_summary_ext.h output_routines_ext.h sh_util.h
write_beagle.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_beagle.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_beagle.o: str_utils.hh genetic_utils_ext.h write_files_ext.h
plink_core_files.o: common.h types.hh typedefs.h errorno.h allutil.h
plink_core_files.o: linkage.h linkage_external.h makeped.h pedtree.h
plink_core_files.o: analysis.h analysis_ext.h write_fbat_ext.h
plink_core_files.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
plink_core_files.o: write_pangaea_ext.h write_beagle_ext.h
plink_core_files.o: write_structure_ext.h write_pseq_ext.h batch_input.h
plink_core_files.o: input_check.h hwe.h compress_ext.h tod.hh loop.h
plink_core_files.o: file_ops.h zlib.h zconf.h error_messages_ext.h
plink_core_files.o: utils_ext.h entry.h create_summary_ext.h
plink_core_files.o: output_routines_ext.h sh_util.h linkage_ext.h fcmap_ext.h
plink_core_files.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
plink_core_files.o: input_ops.hh str_utils.hh
write_eigenstrat.o: common.h types.hh typedefs.h errorno.h allutil.h
write_eigenstrat.o: linkage.h linkage_external.h makeped.h pedtree.h
write_eigenstrat.o: analysis.h analysis_ext.h write_fbat_ext.h
write_eigenstrat.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
write_eigenstrat.o: write_pangaea_ext.h write_beagle_ext.h
write_eigenstrat.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_eigenstrat.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_eigenstrat.o: zlib.h zconf.h error_messages_ext.h utils_ext.h entry.h
write_eigenstrat.o: create_summary_ext.h output_routines_ext.h sh_util.h
write_eigenstrat.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_eigenstrat.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_eigenstrat.o: str_utils.hh
write_structure.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_structure.o: linkage_external.h makeped.h pedtree.h analysis.h
write_structure.o: analysis_ext.h write_fbat_ext.h write_plink_ext.h
write_structure.o: plink_core_ext.h write_eigenstrat_ext.h
write_structure.o: write_pangaea_ext.h write_beagle_ext.h
write_structure.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_structure.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_structure.o: zlib.h zconf.h error_messages_ext.h utils_ext.h entry.h
write_structure.o: create_summary_ext.h output_routines_ext.h sh_util.h
write_structure.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_structure.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_structure.o: str_utils.hh
batch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
batch.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
batch.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
batch.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
batch.o: write_structure_ext.h write_pseq_ext.h batch_input.h input_check.h
batch.o: hwe.h compress_ext.h batch.h reorder_loci_ext.h batch_input_ext.h
batch.o: error_messages_ext.h genetic_utils_ext.h user_input_ext.h
batch.o: input_ops.hh str_utils.hh plink_ext.h utils_ext.h cw_routines_ext.h
compress.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
compress.o: linkage_external.h makeped.h pedtree.h analysis.h analysis_ext.h
compress.o: write_fbat_ext.h write_plink_ext.h plink_core_ext.h
compress.o: write_eigenstrat_ext.h write_pangaea_ext.h write_beagle_ext.h
compress.o: write_structure_ext.h write_pseq_ext.h batch_input.h
compress.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
compress.o: utils_ext.h annotated_ped_file.h
write_pseq_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_pseq_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_pseq_files.o: analysis.h analysis_ext.h write_fbat_ext.h
write_pseq_files.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
write_pseq_files.o: write_pangaea_ext.h write_beagle_ext.h
write_pseq_files.o: write_structure_ext.h write_pseq_ext.h batch_input.h
write_pseq_files.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_pseq_files.o: zlib.h zconf.h error_messages_ext.h utils_ext.h entry.h
write_pseq_files.o: create_summary_ext.h output_routines_ext.h sh_util.h
write_pseq_files.o: linkage_ext.h fcmap_ext.h omit_ped_ext.h
write_pseq_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_pseq_files.o: str_utils.hh
menu_value_missing.o: common.h types.hh typedefs.h errorno.h allutil.h
menu_value_missing.o: linkage.h linkage_external.h makeped.h pedtree.h
menu_value_missing.o: analysis.h analysis_ext.h write_fbat_ext.h
menu_value_missing.o: write_plink_ext.h plink_core_ext.h
menu_value_missing.o: write_eigenstrat_ext.h write_pangaea_ext.h
menu_value_missing.o: write_beagle_ext.h write_structure_ext.h
menu_value_missing.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
menu_value_missing.o: compress_ext.h batch_input_ext.h error_messages_ext.h
menu_value_missing.o: fcmap_ext.h plink_ext.h user_input_ext.h input_ops.hh
menu_value_missing.o: str_utils.hh utils_ext.h
str_utils.o: common.h types.hh error_messages_ext.h str_utils.hh
read_impute.o: common.h types.hh tod.hh error_messages_ext.h typedefs.h
read_impute.o: errorno.h allutil.h linkage.h linkage_external.h makeped.h
read_impute.o: pedtree.h analysis.h analysis_ext.h write_fbat_ext.h
read_impute.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
read_impute.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
read_impute.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
read_impute.o: compress_ext.h fcmap_ext.h batch_input_ext.h mrecode.h
read_impute.o: mrecode_ext.h vcftools/mega2_vcftools_interface.h
read_impute.o: annotated_ped_file.h reorder_loci_ext.h plink_ext.h
read_impute.o: annotated_ped_file_ext.h annotated_ped_file.h read_files_ext.h
read_impute.o: str_utils.hh input_ops.hh read_impute.hh
read_impute_bgen.o: zlib.h zconf.h common.h types.hh tod.hh
read_impute_bgen.o: error_messages_ext.h typedefs.h errorno.h allutil.h
read_impute_bgen.o: linkage.h linkage_external.h makeped.h pedtree.h
read_impute_bgen.o: analysis.h analysis_ext.h write_fbat_ext.h
read_impute_bgen.o: write_plink_ext.h plink_core_ext.h write_eigenstrat_ext.h
read_impute_bgen.o: write_pangaea_ext.h write_beagle_ext.h
read_impute_bgen.o: write_structure_ext.h write_pseq_ext.h batch_input.h
read_impute_bgen.o: input_check.h hwe.h compress_ext.h fcmap_ext.h
read_impute_bgen.o: batch_input_ext.h mrecode.h mrecode_ext.h
read_impute_bgen.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
read_impute_bgen.o: reorder_loci_ext.h plink_ext.h read_files_ext.h
read_impute_bgen.o: annotated_ped_file_ext.h annotated_ped_file.h
read_impute_bgen.o: input_ops.hh str_utils.hh read_impute.hh
read_impute_bgen.o: read_impute_bgen.hh
read_impute_bgen2.o: zlib.h zconf.h common.h types.hh tod.hh
read_impute_bgen2.o: error_messages_ext.h typedefs.h errorno.h allutil.h
read_impute_bgen2.o: linkage.h linkage_external.h makeped.h pedtree.h
read_impute_bgen2.o: analysis.h analysis_ext.h write_fbat_ext.h
read_impute_bgen2.o: write_plink_ext.h plink_core_ext.h
read_impute_bgen2.o: write_eigenstrat_ext.h write_pangaea_ext.h
read_impute_bgen2.o: write_beagle_ext.h write_structure_ext.h
read_impute_bgen2.o: write_pseq_ext.h batch_input.h input_check.h hwe.h
read_impute_bgen2.o: compress_ext.h fcmap_ext.h batch_input_ext.h mrecode.h
read_impute_bgen2.o: mrecode_ext.h vcftools/mega2_vcftools_interface.h
read_impute_bgen2.o: annotated_ped_file.h reorder_loci_ext.h plink_ext.h
read_impute_bgen2.o: read_files_ext.h annotated_ped_file_ext.h
read_impute_bgen2.o: annotated_ped_file.h input_ops.hh str_utils.hh
read_impute_bgen2.o: read_impute.hh read_impute_bgen2.hh
