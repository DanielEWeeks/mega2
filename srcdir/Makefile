#   Mega2: Manipulation Environment for Genetic Analysis.
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
            GCCVERSION := $(shell gcc -dumpversion)
	endif
    else
        OSTYPE := $(shell uname -s)
    endif
    $(info type $(OSTYPE) gcc $(GCCVERSION))
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
    else ifeq "$(findstring MINGW64,$(OSTYPE))" "MINGW64"
	OSTYPE = msys2
    else ifeq "$(findstring MSYS_NT,$(OSTYPE))" "MSYS_NT"
        ifeq ($(GCCVERSION),7.3.0)
		OSTYPE = msys2_7
	else
		OSTYPE = msys2
	endif
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
	write_SUP_files.cpp write_mendel7_files.cpp \
	write_plink_files.cpp cranefoot.cpp      write_mega2annot.cpp \
	grow_string.cpp     marker_lookup.cpp    write_IQLS.cpp \
	net.cpp             phe_lookup.cpp       str_utils.cpp \
	entry.cpp           loop.cpp             compress.cpp \
        analysis.cpp        class_old.cpp        \
	write_pangaea.cpp   write_beagle.cpp     plink_core_files.cpp \
	write_eigenstrat.cpp write_structure.cpp batch.cpp \
	write_fbat.cpp      write_pseq_files.cpp menu_value_missing.cpp \
	read_impute.cpp	    read_impute_bgen.cpp read_impute_bgen2.cpp \
	write_shapeit_files.cpp                  write_roadtrips.cpp \
	write_mach.cpp      write_minimacshapeit.cpp        write_vcf.cpp \
	dbmega2.cpp         dbmisc.cpp           dbbatch.cpp \
	dbpedigree.cpp      dballele.cpp         dbloci.cpp \
	dbmap.cpp           dbgenotype.cpp       dbrefallele.cpp \
	read_bcfs.cpp 	    mega2_bcftools_interface.cpp \
        write_MQLS.cpp \

OBJS =  $(SRCS:.cpp=.o)

# The flags for including the libraries:

iFLAGS = -D_mega2_ -I$(TOP)
IFLAGS = -D_mega2_ -I$(shell pwd)

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

dbmega2.o dbmisc.o dbbatch.o dbpedigree.o dballele.o dbloci.o dbmap.o dbgenotype.o: iFLAGS += -Ilib/sqlite3
dbrefallele.o write_vcf.o : iFLAGS += -Ilib/sqlite3

user_input.o read_impute_bgen2.o: iFLAGS += -Ibgen/genfile/include -Ibgen/zstd-1.1.0/lib $(stdC11)

write_IQLS.o write_beagle.o write_fbat.o write_loki.o write_mfiles.o: DFLAGS += -DRUNSHELL_SETUP
write_pangaea.o write_plink_files.o write_premakeped.o: DFLAGS += -DRUNSHELL_SETUP
write_structure.o: DFLAGS += -DRUNSHELL_SETUP

BCFTOOLS_VERSION=1.6
mega2_bcftools_interface.o read_bcfs.o write_vcf.o: iFLAGS += -Ilib/bcftools-$(BCFTOOLS_VERSION) -Ilib/bcftools-$(BCFTOOLS_VERSION)/htslib-$(BCFTOOLS_VERSION) -Ilib/bcftools-$(BCFTOOLS_VERSION)/htslib-$(BCFTOOLS_VERSION)/htslib

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
R_output.o dbmisc.o hwe_user_input.o mega2.o mrecode.o utils.o: DFLAGS += $(DtFLAGS)
read_files.o write_sage_files.o: DFLAGS += $(DtFLAGS)

# -DDATE
dbmega2.o stamp_time.o: DFLAGS += $(DdFLAGS)

# -DHIDESTATUS
aspex.o annotated_ped_file.o batch_input.o error_messages.o genetic_utils.o: DFLAGS += $(DsFLAGS)
makeped1.o mega2.o menu_value_missing.o mrecode.o: DFLAGS += $(DsFLAGS)
output_routines.o read_files.o reorder_loci.o user_input.o utils.o write_files.o: DFLAGS += $(DsFLAGS)
write_plink_files.o: DFLAGS += $(DsFLAGS)

# -DHIDEPATH
aspex.o cranefoot.o mega2.o R_output.o write_SUP_files.o write_ghfiles.o write_premakeped.o: DFLAGS += $(DpFLAGS)
dbbatch.o: DFLAGS += $(DpFLAGS)

# -DHIDEFILE
annotated_ped_file.o makeped1.o mega2.o read_files.o utils.o: DFLAGS += $(DfFLAGS)
write_IQLS.o write_beagle.o write_eigenstrat.o write_fbat.o: DFLAGS += $(DfFLAGS)
write_impute2_files.o write_mach.o write_minimacshapeit.o: DFLAGS += $(DfFLAGS)
write_pangaea.o write_plink_files.o: DFLAGS += $(DfFLAGS)
write_pseq_files.o write_roadtrips.o: DFLAGS += $(DfFLAGS)
write_shapeit_files.o write_structure.o: DFLAGS += $(DfFLAGS)
write_MQLS.o: DFLAGS += $(DfFLAGS)
dbmega2.o: DFLAGS += $(DfFLAGS)

# -DDATE -DHIDEPATH
write_mfiles.o write_vcf.o write_vitesse.o utils.o: DFLAGS += $(DdFLAGS) $(DpFLAGS)

ifeq  "$(HIDE)" "2"
    DaFLAGS += -DSHOWSTATUS
endif
# -DSHOWSTATUS
annotated_ped_file.o mrecode.o: DFLAGS += $(DaFLAGS)

ifeq  "$(TEST)" "2"
    DlFLAGS += -DLOGDIR
endif
output_file_names.o: DFLAGS += $(DlFLAGS)
mega2.o: DFLAGS += -DSORT_HETEROZYGOUS

# export CFLAGS so that the called Makefile can find it...
export FLAGS
export OSTYPE
export CFLAGS
export IFLAGS
export CXX
export CC
export AR
export ARFLAGS
export AROUT
export RANLIB
export SQLxDEFS

LINK ?= $(CXX)
LTARGET ?= -o $(EXEC)
AR ?= ar
ARFLAGS ?= rcs
AROUT ?=
RANLIB = ranlib

ZLIBS ?= libz.a
CLIBS := libvcftools.a libbgen.a $(ZLIBS) libsqlite3.a libbcftools.a $(CLIBS)

##BCFFLAGS =  -lz -lm -lbz2 # –llzma
stdC11  ?= -std=c++0x
WINO    ?=

.c.o:
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CC) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(iFLAGS) $*.c $(WINO)

.cpp.o:	
ifeq ($(MINGW),)
	$(RM) $*.o
endif
	$(CXX) $(CFLAGS) $(DFLAGS) $(FLAGS) -c $(iFLAGS) $*.cpp $(WINO)

all: $(MODULE)_$(OSTYPE)


dbg:
	$(MAKE) $(MAKEFLAGS) TEST=2 HIDE=1 OPT=-g all

tod:
	$(MAKE) $(MAKEFLAGS) TEST=2 HIDE=1 OPT=-g FLAGS=-DTOD all

gcc_info.o: $(OBJS)

SUBDIRS = bgen lib/sqlite3 vcftools zlib-1.2.8 lib/bcftools-$(BCFTOOLS_VERSION)

.PHONY: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) $(MFLAGS) -C $@

libbgen.a: bgen
bgen: zlib-1.2.8

libsqlite3.a: lib/sqlite3

vcftools: zlib-1.2.8
libvcftools.a: vcftools

libz.a: zlib-1.2.8

libbcftools.a: lib/bcftools-$(BCFTOOLS_VERSION)

CFG = bcftool-$(BCFTOOLS_VERSION)scfg

bcftool-$(BCFTOOLS_VERSION)scfg:
	$(MAKE) $(MFLAGS) -C lib/bcftools-$(BCFTOOLS_VERSION) config

$(MODULE)_$(OSTYPE): $(CFG) $(SUBDIRS) $(OBJS) $(HDRS) gcc_info.o libvcftools.a libbgen.a libz.a libsqlite3.a libbcftools.a
	$(RM) $(MODULE)
	$(LINK) $(LFLAGS) $(LTARGET) $(OBJS) $(FLAGS) gcc_info.o $(CLIBS) $(BCFFLAGS)

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
	$(MAKE) $(MFLAGS) -C lib/bcftools-$(BCFTOOLS_VERSION) RM="$(RM)" clean

distclean:
	$(RM) $(OBJS) gcc_info.o $(EXEC) $(MODULE) $(EXEC).exe $(MODULE).exe
#	$(foreach dir,$(SUBDIRS),$(MAKE) -C $(dir) RM="$(RM)" distclean &&) echo on
	$(MAKE) $(MFLAGS) -C vcftools    RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C bgen        RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C zlib-1.2.8  RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C lib/sqlite3 RM="$(RM)" distclean
	$(MAKE) $(MFLAGS) -C lib/bcftools-$(BCFTOOLS_VERSION) RM="$(RM)" distclean

clean1:
	$(RM) *~

srcdist:
	for i in $(SRCS) $(HDRS); do $(REMOVECVS1) $$i > ../srcdir/$$i; done
	for i in `ls Makedefs*`; do $(REMOVECVS) $$i > ../srcdir/$$i.template; done
	$(REMOVECVS) Makefile > ../srcdir/Makefile

depend:
	makedepend -Y. -- $(CFLAGS) -- $(SRCS)
	$(MAKE) -C vcftools    depend
	$(MAKE) -C bgen        depend
	$(MAKE) -C lib/sqlite3 depend

# DO NOT DELETE

input_check.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
input_check.o: linkage_external.h makeped.h pedtree.h batch_input.h
input_check.o: input_check.h hwe.h compress_ext.h class_old.h analysis.h
input_check.o: write_plink_ext.h plink_core_ext.h tod.hh check_ext.h
input_check.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
input_check.o: genetic_utils_ext.h grow_string_ext.h user_input_ext.h
input_check.o: input_ops.hh str_utils.hh annotated_ped_file.h
input_check.o: input_check_ext.h pedtree_ext.h utils_ext.h write_pap_ext.h
input_check.o: batch_input_ext.h
select_index.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
select_index.o: linkage_external.h makeped.h pedtree.h batch_input.h
select_index.o: input_check.h hwe.h compress_ext.h random_choose_ext.h
aspex.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
aspex.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
aspex.o: hwe.h compress_ext.h batch_input_ext.h create_summary_ext.h
aspex.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
aspex.o: grow_string_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
aspex.o: output_file_names_ext.h output_routines_ext.h scripts_ext.h
aspex.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
aspex.o: utils_ext.h write_files_ext.h class_old.h analysis.h
aspex.o: write_plink_ext.h plink_core_ext.h
liability_dist.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
liability_dist.o: linkage_external.h makeped.h pedtree.h batch_input.h
liability_dist.o: input_check.h hwe.h compress_ext.h liable.h
liability_dist.o: error_messages_ext.h genetic_utils_ext.h
liability_dist.o: liability_dist_ext.h linkage_ext.h omit_ped_ext.h
liability_dist.o: output_file_names_ext.h utils_ext.h class_old.h analysis.h
liability_dist.o: write_plink_ext.h plink_core_ext.h
select_individuals.o: common.h types.hh typedefs.h errorno.h allutil.h
select_individuals.o: linkage.h linkage_external.h makeped.h pedtree.h
select_individuals.o: batch_input.h input_check.h hwe.h compress_ext.h
select_individuals.o: mrecode.h batch_input_ext.h error_messages_ext.h
select_individuals.o: fcmap_ext.h utils_ext.h
linkage.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
linkage.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
linkage.o: hwe.h compress_ext.h analysis.h error_messages_ext.h list_ext.h
linkage.o: makeped_ext.h pedtree_ext.h utils_ext.h user_input_ext.h
linkage.o: input_ops.hh str_utils.hh annotated_ped_file.h
linkage.o: append_locus_array_ext.h
output_file_names.o: common.h types.hh typedefs.h errorno.h allutil.h
output_file_names.o: linkage.h linkage_external.h makeped.h pedtree.h
output_file_names.o: batch_input.h input_check.h hwe.h compress_ext.h
output_file_names.o: batch_input_ext.h error_messages_ext.h fcmap_ext.h
output_file_names.o: genetic_utils_ext.h output_file_names_ext.h utils_ext.h
output_file_names.o: class_old.h analysis.h write_plink_ext.h
output_file_names.o: plink_core_ext.h
cal_const.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_const.o: linkage_external.h makeped.h pedtree.h batch_input.h
cal_const.o: input_check.h hwe.h compress_ext.h hwe_ext.h
list.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
list.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
list.o: hwe.h compress_ext.h error_messages_ext.h
test_switch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
test_switch.o: linkage_external.h makeped.h pedtree.h batch_input.h
test_switch.o: input_check.h hwe.h compress_ext.h
cal_n.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_n.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
cal_n.o: hwe.h compress_ext.h
ln_p_value.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
ln_p_value.o: linkage_external.h makeped.h pedtree.h batch_input.h
ln_p_value.o: input_check.h hwe.h compress_ext.h hwe_ext.h
user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
user_input.o: linkage_external.h makeped.h pedtree.h batch_input.h
user_input.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
user_input.o: cw_routines_ext.h error_messages_ext.h fcmap_ext.h
user_input.o: genetic_utils_ext.h reorder_loci_ext.h grow_string_ext.h
user_input.o: output_file_names_ext.h output_routines_ext.h plink_ext.h
user_input.o: utils_ext.h user_input_ext.h input_ops.hh str_utils.hh
user_input.o: annotated_ped_file.h vcftools/mega2_vcftools_interface.h
user_input.o: annotated_ped_file.h write_files_ext.h class_old.h analysis.h
user_input.o: write_plink_ext.h plink_core_ext.h database_dump_ext.h input.hh
user_input.o: read_impute.hh read_impute_bgen.hh read_impute_bgen2.hh
user_input.o: read_bcfs.h
cal_prob.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cal_prob.o: linkage_external.h makeped.h pedtree.h batch_input.h
cal_prob.o: input_check.h hwe.h compress_ext.h do_switch_ext.h
cal_prob.o: test_switch_ext.h utils_ext.h
utils.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
utils.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
utils.o: hwe.h compress_ext.h batch_input_ext.h error_messages_ext.h
utils.o: fcmap_ext.h net_ext.h user_input_ext.h input_ops.hh str_utils.hh
utils.o: annotated_ped_file.h utils_ext.h version.h
check.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
check.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
check.o: hwe.h compress_ext.h tod.hh utils_ext.h error_messages_ext.h
check.o: grow_string_ext.h write_files_ext.h class_old.h analysis.h
check.o: write_plink_ext.h plink_core_ext.h
mega2.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
mega2.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
mega2.o: hwe.h compress_ext.h tod.hh R_output.h annotated_ped_file.h
mega2.o: annotated_ped_file_ext.h aspex_ext.h batch_input_ext.h
mega2.o: cranefoot_ext.h create_summary_ext.h error_messages_ext.h
mega2.o: error_sim_ext.h fcmap_ext.h genetic_utils_ext.h hwe_user_input_ext.h
mega2.o: input_check_ext.h makenucs_ext.h makeped_ext.h marker_lookup_ext.h
mega2.o: mega2annot_ext.h menu_value_missing_ext.h mrecode.h
mega2.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
mega2.o: plink_ext.h read_files_ext.h reorder_loci_ext.h
mega2.o: select_individuals_ext.h scripts_ext.h slink_ext.h splink_ext.h
mega2.o: user_input_ext.h input_ops.hh str_utils.hh
mega2.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h utils_ext.h
mega2.o: write_IQLS_ext.h write_SUP_files_ext.h write_apm_ext.h
mega2.o: write_files_ext.h write_ghfiles_ext.h write_loki_ext.h
mega2.o: write_mfiles_ext.h write_pap_ext.h write_premakeped_ext.h
mega2.o: write_prest_ext.h write_sage_files_ext.h write_simulate_files_ext.h
mega2.o: write_solar_files_ext.h write_vitesse_ext.h version.h class_old.h
mega2.o: analysis.h write_plink_ext.h plink_core_ext.h write_shapeit_ext.h
mega2.o: database_dump_ext.h
write_apm.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_apm.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_apm.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_apm.o: fcmap_ext.h genetic_utils_ext.h input_check_ext.h linkage_ext.h
write_apm.o: omit_ped_ext.h output_file_names_ext.h pedtree_ext.h
write_apm.o: reorder_loci_ext.h user_input_ext.h input_ops.hh str_utils.hh
write_apm.o: annotated_ped_file.h utils_ext.h
output_routines.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
output_routines.o: linkage_external.h makeped.h pedtree.h batch_input.h
output_routines.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
output_routines.o: grow_string_ext.h output_routines_ext.h utils_ext.h
write_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_files.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_files.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_files.o: fcmap_ext.h genetic_utils_ext.h grow_string_ext.h
write_files.o: linkage_ext.h makeped_ext.h omit_ped_ext.h
write_files.o: output_file_names_ext.h reorder_loci_ext.h user_input_ext.h
write_files.o: input_ops.hh str_utils.hh annotated_ped_file.h utils_ext.h
write_files.o: write_ghfiles_ext.h plink_ext.h class_old.h analysis.h
write_files.o: write_plink_ext.h plink_core_ext.h
do_switch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
do_switch.o: linkage_external.h makeped.h pedtree.h batch_input.h
do_switch.o: input_check.h hwe.h compress_ext.h
pedtree.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
pedtree.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
pedtree.o: hwe.h compress_ext.h tod.hh input_check_ext.h error_messages_ext.h
pedtree.o: linkage_ext.h
write_ghfiles.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_ghfiles.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_ghfiles.o: input_check.h hwe.h compress_ext.h R_output.h R_output_ext.h
write_ghfiles.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_ghfiles.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
write_ghfiles.o: output_file_names_ext.h output_routines_ext.h
write_ghfiles.o: user_input_ext.h input_ops.hh str_utils.hh
write_ghfiles.o: annotated_ped_file.h utils_ext.h write_files_ext.h
write_ghfiles.o: class_old.h analysis.h write_plink_ext.h plink_core_ext.h
fcmap.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
fcmap.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
fcmap.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h
print_data.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
print_data.o: linkage_external.h makeped.h pedtree.h batch_input.h
print_data.o: input_check.h hwe.h compress_ext.h
write_liability.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_liability.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_liability.o: input_check.h hwe.h compress_ext.h liable.h
write_liability.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_liability.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_liability.o: str_utils.hh annotated_ped_file.h
random_choose.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
random_choose.o: linkage_external.h makeped.h pedtree.h batch_input.h
random_choose.o: input_check.h hwe.h compress_ext.h utils_ext.h
write_mfiles.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_mfiles.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_mfiles.o: input_check.h hwe.h compress_ext.h R_output.h write_mfiles.h
write_mfiles.o: batch_input_ext.h R_output_ext.h create_summary_ext.h
write_mfiles.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_mfiles.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_mfiles.o: output_routines_ext.h user_input_ext.h input_ops.hh
write_mfiles.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_mfiles.o: write_mendel7_files_ext.h write_files_ext.h
write_mfiles.o: write_mfiles_ext.h class_old.h analysis.h write_plink_ext.h
write_mfiles.o: plink_core_ext.h
mrecode.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
mrecode.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
mrecode.o: hwe.h compress_ext.h tod.hh mrecode.h mrecode_ext.h
mrecode.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
mrecode.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h linkage_ext.h
mrecode.o: makeped_ext.h marker_lookup_ext.h omit_ped_ext.h read_files_ext.h
mrecode.o: select_individuals_ext.h utils_ext.h write_files_ext.h class_old.h
mrecode.o: analysis.h write_plink_ext.h plink_core_ext.h
gen_user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
gen_user_input.o: linkage_external.h makeped.h pedtree.h batch_input.h
gen_user_input.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
gen_user_input.o: fcmap_ext.h genetic_utils_ext.h output_file_names_ext.h
gen_user_input.o: user_input_ext.h input_ops.hh str_utils.hh
gen_user_input.o: annotated_ped_file.h utils_ext.h
write_sage_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_sage_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_sage_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_sage_files.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
write_sage_files.o: genetic_utils_ext.h linkage_ext.h makenucs_ext.h
write_sage_files.o: omit_ped_ext.h output_file_names_ext.h
write_sage_files.o: output_routines_ext.h read_files_ext.h user_input_ext.h
write_sage_files.o: input_ops.hh str_utils.hh annotated_ped_file.h
write_sage_files.o: utils_ext.h write_files_ext.h class_old.h analysis.h
write_sage_files.o: write_plink_ext.h plink_core_ext.h
splink.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
splink.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
splink.o: hwe.h compress_ext.h create_summary_ext.h error_messages_ext.h
splink.o: fcmap_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
splink.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
splink.o: input_ops.hh str_utils.hh annotated_ped_file.h utils_ext.h
splink.o: write_files_ext.h
genetic_utils.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
genetic_utils.o: linkage_external.h makeped.h pedtree.h batch_input.h
genetic_utils.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
genetic_utils.o: grow_string_ext.h utils_ext.h user_input_ext.h input_ops.hh
genetic_utils.o: str_utils.hh annotated_ped_file.h
read_data.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
read_data.o: linkage_external.h makeped.h pedtree.h batch_input.h
read_data.o: input_check.h hwe.h compress_ext.h
write_simulate_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_simulate_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_simulate_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_simulate_files.o: analysis.h error_messages_ext.h fcmap_ext.h
write_simulate_files.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
write_simulate_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_simulate_files.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_simulate_files.o: write_files_ext.h write_ghfiles_ext.h
homozygosity.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
homozygosity.o: linkage_external.h makeped.h pedtree.h batch_input.h
homozygosity.o: input_check.h hwe.h compress_ext.h genetic_utils_ext.h
read_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
read_files.o: linkage_external.h makeped.h pedtree.h batch_input.h
read_files.o: input_check.h hwe.h compress_ext.h mrecode.h
read_files.o: annotated_ped_file.h error_sim.h error_messages_ext.h
read_files.o: fcmap_ext.h grow_string_ext.h linkage_ext.h list_ext.h
read_files.o: makeped_ext.h marker_lookup_ext.h mrecode_ext.h
read_files.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
read_files.o: omit_ped_ext.h plink_ext.h read_files_ext.h reorder_loci_ext.h
read_files.o: utils_ext.h user_input_ext.h input_ops.hh str_utils.hh
read_files.o: batch_input_ext.h write_files_ext.h
write_solar_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_solar_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_solar_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_solar_files.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_solar_files.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_solar_files.o: output_routines_ext.h user_input_ext.h input_ops.hh
write_solar_files.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_solar_files.o: write_files_ext.h class_old.h analysis.h
write_solar_files.o: write_plink_ext.h plink_core_ext.h
hwe.o: hwe.h check_file_ext.h hwe_ext.h read_data_ext.h
reorder_loci.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
reorder_loci.o: linkage_external.h makeped.h pedtree.h batch_input.h
reorder_loci.o: input_check.h hwe.h compress_ext.h annotated_ped_file.h
reorder_loci.o: annotated_ped_file_ext.h batch_input_ext.h
reorder_loci.o: create_summary_ext.h error_messages_ext.h fcmap_ext.h
reorder_loci.o: genetic_utils_ext.h grow_string_ext.h linkage_ext.h
reorder_loci.o: output_routines_ext.h output_file_names_ext.h
reorder_loci.o: user_input_ext.h input_ops.hh str_utils.hh utils_ext.h
reorder_loci.o: write_files_ext.h class_old.h analysis.h write_plink_ext.h
reorder_loci.o: plink_core_ext.h
write_vitesse.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_vitesse.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_vitesse.o: input_check.h hwe.h compress_ext.h batch_input_ext.h
write_vitesse.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_vitesse.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_vitesse.o: user_input_ext.h input_ops.hh str_utils.hh
write_vitesse.o: annotated_ped_file.h utils_ext.h write_files_ext.h
write_vitesse.o: class_old.h analysis.h write_plink_ext.h plink_core_ext.h
hwe_user_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
hwe_user_input.o: linkage_external.h makeped.h pedtree.h batch_input.h
hwe_user_input.o: input_check.h hwe.h compress_ext.h mrecode.h
hwe_user_input.o: batch_input_ext.h error_messages_ext.h fcmap_ext.h
hwe_user_input.o: gen_user_input_ext.h genetic_utils_ext.h grow_string_ext.h
hwe_user_input.o: homozygosity_ext.h hwe_ext.h linkage_ext.h
hwe_user_input.o: output_file_names_ext.h output_routines_ext.h
hwe_user_input.o: select_individuals_ext.h user_input_ext.h input_ops.hh
hwe_user_input.o: str_utils.hh annotated_ped_file.h utils_ext.h
hwe_user_input.o: write_mendel7_files_ext.h write_mfiles_ext.h
hwe_user_input.o: write_files_ext.h class_old.h analysis.h write_plink_ext.h
hwe_user_input.o: plink_core_ext.h
scripts.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
scripts.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
scripts.o: hwe.h compress_ext.h error_messages_ext.h fcmap_ext.h
scripts.o: genetic_utils_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
scripts.o: output_file_names_ext.h user_input_ext.h input_ops.hh str_utils.hh
scripts.o: annotated_ped_file.h utils_ext.h write_files_ext.h
scripts.o: write_ghfiles_ext.h class_old.h analysis.h write_plink_ext.h
scripts.o: plink_core_ext.h
makeped1.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
makeped1.o: linkage_external.h makeped.h pedtree.h batch_input.h
makeped1.o: input_check.h hwe.h compress_ext.h analysis.h makeped_ext.h
makeped1.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h linkage_ext.h
makeped1.o: read_files_ext.h utils_ext.h batch_input_ext.h
omit_ped.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
omit_ped.o: linkage_external.h makeped.h pedtree.h batch_input.h
omit_ped.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
omit_ped.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
omit_ped.o: utils_ext.h
slink.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
slink.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
slink.o: hwe.h compress_ext.h error_messages_ext.h fcmap_ext.h
slink.o: genetic_utils_ext.h linkage_ext.h omit_ped_ext.h
slink.o: output_file_names_ext.h pedtree_ext.h user_input_ext.h input_ops.hh
slink.o: str_utils.hh annotated_ped_file.h utils_ext.h write_ghfiles_ext.h
R_output.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
R_output.o: linkage_external.h makeped.h pedtree.h batch_input.h
R_output.o: input_check.h hwe.h compress_ext.h R_output.h batch_input_ext.h
R_output.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h
R_output.o: output_file_names_ext.h utils_ext.h plink_ext.h user_input_ext.h
R_output.o: input_ops.hh str_utils.hh annotated_ped_file.h class_old.h
R_output.o: analysis.h write_plink_ext.h plink_core_ext.h
write_premakeped.o: common.h types.hh typedefs.h errorno.h allutil.h
write_premakeped.o: linkage.h linkage_external.h makeped.h pedtree.h
write_premakeped.o: batch_input.h input_check.h hwe.h compress_ext.h
write_premakeped.o: R_output.h R_output_ext.h create_summary_ext.h
write_premakeped.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_premakeped.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_premakeped.o: output_routines_ext.h user_input_ext.h input_ops.hh
write_premakeped.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_premakeped.o: write_files_ext.h write_ghfiles_ext.h write_mfiles_ext.h
write_premakeped.o: class_old.h analysis.h write_plink_ext.h plink_core_ext.h
create_summary.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
create_summary.o: linkage_external.h makeped.h pedtree.h batch_input.h
create_summary.o: input_check.h hwe.h compress_ext.h liable.h
create_summary.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
create_summary.o: grow_string_ext.h input_check_ext.h liability_dist_ext.h
create_summary.o: linkage_ext.h makenucs_ext.h omit_ped_ext.h
create_summary.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
create_summary.o: reorder_loci_ext.h user_input_ext.h input_ops.hh
create_summary.o: str_utils.hh annotated_ped_file.h utils_ext.h class_old.h
create_summary.o: analysis.h write_plink_ext.h plink_core_ext.h
write_prest.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_prest.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_prest.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_prest.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_prest.o: output_file_names_ext.h output_routines_ext.h user_input_ext.h
write_prest.o: input_ops.hh str_utils.hh annotated_ped_file.h utils_ext.h
write_prest.o: write_files_ext.h write_ghfiles_ext.h
write_pap.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_pap.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_pap.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_pap.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_pap.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
write_pap.o: utils_ext.h write_files_ext.h class_old.h analysis.h
write_pap.o: write_plink_ext.h plink_core_ext.h
batch_input.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
batch_input.o: linkage_external.h makeped.h pedtree.h batch_input.h
batch_input.o: input_check.h hwe.h compress_ext.h analysis.h input_ops.hh
batch_input.o: str_utils.hh annotated_ped_file.h batch_input_ext.h
batch_input.o: error_messages_ext.h menu_value_missing_ext.h
batch_input.o: genetic_utils_ext.h reorder_loci_ext.h user_input_ext.h
batch_input.o: plink_ext.h utils_ext.h
makenucs.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
makenucs.o: linkage_external.h makeped.h pedtree.h batch_input.h
makenucs.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
makenucs.o: fcmap_ext.h linkage_ext.h makenucs_ext.h omit_ped_ext.h
makenucs.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
makenucs.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
makenucs.o: utils_ext.h write_files_ext.h write_ghfiles_ext.h
makenucs.o: append_locus_array_ext.h class_old.h analysis.h write_plink_ext.h
makenucs.o: plink_core_ext.h
append_locus_array.o: common.h types.hh typedefs.h errorno.h allutil.h
append_locus_array.o: linkage.h linkage_external.h makeped.h pedtree.h
append_locus_array.o: batch_input.h input_check.h hwe.h compress_ext.h
append_locus_array.o: error_messages_ext.h linkage_ext.h
error_sim.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
error_sim.o: linkage_external.h makeped.h pedtree.h batch_input.h
error_sim.o: input_check.h hwe.h compress_ext.h error_sim.h batch_input_ext.h
error_sim.o: error_messages_ext.h genetic_utils_ext.h grow_string_ext.h
error_sim.o: list_ext.h output_file_names_ext.h output_routines_ext.h
error_sim.o: reorder_loci_ext.h utils_ext.h
write_loki.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_loki.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_loki.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_loki.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_loki.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_loki.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_loki.o: write_premakeped_ext.h write_files_ext.h class_old.h analysis.h
write_loki.o: write_plink_ext.h plink_core_ext.h
error_messages.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
error_messages.o: linkage_external.h makeped.h pedtree.h batch_input.h
error_messages.o: input_check.h hwe.h compress_ext.h utils_ext.h
annotated_ped_file.o: common.h types.hh typedefs.h errorno.h allutil.h
annotated_ped_file.o: linkage.h linkage_external.h makeped.h pedtree.h
annotated_ped_file.o: batch_input.h input_check.h hwe.h compress_ext.h tod.hh
annotated_ped_file.o: mrecode.h annotated_ped_file.h annotated_ped_file_ext.h
annotated_ped_file.o: error_messages_ext.h fcmap_ext.h grow_string_ext.h
annotated_ped_file.o: linkage_ext.h list_ext.h makeped_ext.h
annotated_ped_file.o: marker_lookup_ext.h mrecode_ext.h
annotated_ped_file.o: vcftools/mega2_vcftools_interface.h
annotated_ped_file.o: annotated_ped_file.h omit_ped_ext.h phe_lookup_ext.h
annotated_ped_file.o: read_files_ext.h reorder_loci_ext.h
annotated_ped_file.o: select_individuals_ext.h user_input_ext.h input_ops.hh
annotated_ped_file.o: str_utils.hh utils_ext.h write_files_ext.h plink_ext.h
annotated_ped_file.o: read_impute.hh class_old.h analysis.h write_plink_ext.h
annotated_ped_file.o: plink_core_ext.h
write_SUP_files.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_SUP_files.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_SUP_files.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
write_SUP_files.o: fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_SUP_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_SUP_files.o: str_utils.hh annotated_ped_file.h utils_ext.h
write_SUP_files.o: write_files_ext.h class_old.h analysis.h write_plink_ext.h
write_SUP_files.o: plink_core_ext.h
write_mendel7_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_mendel7_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_mendel7_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_mendel7_files.o: write_mfiles.h create_summary_ext.h
write_mendel7_files.o: error_messages_ext.h genetic_utils_ext.h
write_mendel7_files.o: output_routines_ext.h utils_ext.h user_input_ext.h
write_mendel7_files.o: input_ops.hh str_utils.hh annotated_ped_file.h
write_mendel7_files.o: write_files_ext.h write_mendel7_files_ext.h
write_mendel7_files.o: class_old.h analysis.h write_plink_ext.h
write_mendel7_files.o: plink_core_ext.h
write_plink_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_plink_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_plink_files.o: batch_input.h input_check.h hwe.h compress_ext.h tod.hh
write_plink_files.o: loop.h file_ops.h error_messages_ext.h utils_ext.h
write_plink_files.o: entry.h create_summary_ext.h output_routines_ext.h
write_plink_files.o: linkage_ext.h sh_util.h fcmap_ext.h omit_ped_ext.h
write_plink_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_plink_files.o: str_utils.hh annotated_ped_file.h
write_plink_files.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
write_plink_files.o: write_plink_ext.h plink_core_ext.h analysis.h
cranefoot.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
cranefoot.o: linkage_external.h makeped.h pedtree.h batch_input.h
cranefoot.o: input_check.h hwe.h compress_ext.h create_summary_ext.h
cranefoot.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
cranefoot.o: output_file_names_ext.h user_input_ext.h input_ops.hh
cranefoot.o: str_utils.hh annotated_ped_file.h utils_ext.h write_files_ext.h
cranefoot.o: class_old.h analysis.h write_plink_ext.h plink_core_ext.h
write_mega2annot.o: common.h types.hh typedefs.h errorno.h allutil.h
write_mega2annot.o: linkage.h linkage_external.h makeped.h pedtree.h
write_mega2annot.o: batch_input.h input_check.h hwe.h compress_ext.h
write_mega2annot.o: annotated_ped_file.h annotated_ped_file_ext.h
write_mega2annot.o: error_messages_ext.h fcmap_ext.h genetic_utils_ext.h
write_mega2annot.o: linkage_ext.h omit_ped_ext.h output_file_names_ext.h
write_mega2annot.o: output_routines_ext.h user_input_ext.h input_ops.hh
write_mega2annot.o: str_utils.hh utils_ext.h write_files_ext.h class_old.h
write_mega2annot.o: analysis.h write_plink_ext.h plink_core_ext.h
grow_string.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
grow_string.o: linkage_external.h makeped.h pedtree.h batch_input.h
grow_string.o: input_check.h hwe.h compress_ext.h utils_ext.h
marker_lookup.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
marker_lookup.o: linkage_external.h makeped.h pedtree.h batch_input.h
marker_lookup.o: input_check.h hwe.h compress_ext.h
write_IQLS.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_IQLS.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_IQLS.o: input_check.h hwe.h compress_ext.h create_summary_ext.h
write_IQLS.o: error_messages_ext.h fcmap_ext.h linkage_ext.h omit_ped_ext.h
write_IQLS.o: output_file_names_ext.h output_routines_ext.h pedtree_ext.h
write_IQLS.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
write_IQLS.o: utils_ext.h loop.h file_ops.h entry.h sh_util.h class_old.h
write_IQLS.o: analysis.h write_plink_ext.h plink_core_ext.h
net.o: common.h types.hh error_messages_ext.h
phe_lookup.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
phe_lookup.o: linkage_external.h makeped.h pedtree.h batch_input.h
phe_lookup.o: input_check.h hwe.h compress_ext.h utils_ext.h
phe_lookup.o: genetic_utils_ext.h error_messages_ext.h plink_ext.h
phe_lookup.o: phe_lookup_ext.h
str_utils.o: common.h types.hh error_messages_ext.h str_utils.hh
entry.o: file_ops.h common.h types.hh typedefs.h errorno.h allutil.h
entry.o: linkage.h linkage_external.h makeped.h pedtree.h batch_input.h
entry.o: input_check.h hwe.h compress_ext.h error_messages_ext.h utils_ext.h
entry.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
entry.o: entry.h create_summary_ext.h output_routines_ext.h write_files_ext.h
loop.o: loop.h file_ops.h common.h types.hh typedefs.h errorno.h allutil.h
loop.o: linkage.h linkage_external.h makeped.h pedtree.h batch_input.h
loop.o: input_check.h hwe.h compress_ext.h error_messages_ext.h utils_ext.h
loop.o: entry.h create_summary_ext.h output_routines_ext.h linkage_ext.h
loop.o: tod.hh analysis.h output_file_names_ext.h
compress.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
compress.o: linkage_external.h makeped.h pedtree.h batch_input.h
compress.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
compress.o: utils_ext.h annotated_ped_file.h
analysis.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
analysis.o: linkage_external.h makeped.h pedtree.h batch_input.h
analysis.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
analysis.o: utils_ext.h analysis.h class_old.h write_plink_ext.h
analysis.o: plink_core_ext.h write_fbat_ext.h write_eigenstrat_ext.h
analysis.o: write_pangaea_ext.h write_beagle_ext.h write_structure_ext.h
analysis.o: write_pseq_ext.h write_shapeit_ext.h write_roadtrips_ext.h
analysis.o: write_mach_ext.h write_minimacshapeit_ext.h database_dump_ext.h
analysis.o: write_vcf_ext.h write_MQLS_ext.h menu_value_missing.h
class_old.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
class_old.o: linkage_external.h makeped.h pedtree.h batch_input.h
class_old.o: input_check.h hwe.h compress_ext.h analysis.h utils_ext.h
class_old.o: error_messages_ext.h fcmap_ext.h output_file_names_ext.h
class_old.o: aspex_ext.h cranefoot_ext.h create_summary_ext.h
class_old.o: hwe_user_input_ext.h makenucs_ext.h mega2annot_ext.h
class_old.o: scripts_ext.h slink_ext.h splink_ext.h write_IQLS_ext.h
class_old.o: write_SUP_files_ext.h write_apm_ext.h write_files_ext.h
class_old.o: write_ghfiles_ext.h write_loki_ext.h write_mfiles_ext.h
class_old.o: write_pap_ext.h write_premakeped_ext.h write_prest_ext.h
class_old.o: write_sage_files_ext.h write_simulate_files_ext.h
class_old.o: write_solar_files_ext.h write_vitesse_ext.h class_old.h
class_old.o: write_plink_ext.h plink_core_ext.h
write_pangaea.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_pangaea.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_pangaea.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_pangaea.o: error_messages_ext.h utils_ext.h entry.h
write_pangaea.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_pangaea.o: sh_util.h fcmap_ext.h genetic_utils_ext.h omit_ped_ext.h
write_pangaea.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_pangaea.o: str_utils.hh annotated_ped_file.h write_pangaea_ext.h
write_pangaea.o: analysis.h
write_beagle.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_beagle.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_beagle.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_beagle.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_beagle.o: output_routines_ext.h linkage_ext.h sh_util.h fcmap_ext.h
write_beagle.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_beagle.o: input_ops.hh str_utils.hh annotated_ped_file.h
write_beagle.o: genetic_utils_ext.h write_files_ext.h write_beagle_ext.h
write_beagle.o: analysis.h
plink_core_files.o: common.h types.hh typedefs.h errorno.h allutil.h
plink_core_files.o: linkage.h linkage_external.h makeped.h pedtree.h
plink_core_files.o: batch_input.h input_check.h hwe.h compress_ext.h tod.hh
plink_core_files.o: loop.h file_ops.h error_messages_ext.h utils_ext.h
plink_core_files.o: entry.h create_summary_ext.h output_routines_ext.h
plink_core_files.o: linkage_ext.h sh_util.h fcmap_ext.h omit_ped_ext.h
plink_core_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
plink_core_files.o: str_utils.hh annotated_ped_file.h plink_core_ext.h
plink_core_files.o: analysis.h
write_eigenstrat.o: common.h types.hh typedefs.h errorno.h allutil.h
write_eigenstrat.o: linkage.h linkage_external.h makeped.h pedtree.h
write_eigenstrat.o: batch_input.h input_check.h hwe.h compress_ext.h loop.h
write_eigenstrat.o: file_ops.h error_messages_ext.h utils_ext.h entry.h
write_eigenstrat.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_eigenstrat.o: sh_util.h plink_core_ext.h analysis.h fcmap_ext.h
write_eigenstrat.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_eigenstrat.o: input_ops.hh str_utils.hh annotated_ped_file.h
write_eigenstrat.o: write_eigenstrat_ext.h
write_structure.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_structure.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_structure.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_structure.o: error_messages_ext.h utils_ext.h entry.h
write_structure.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_structure.o: sh_util.h fcmap_ext.h omit_ped_ext.h
write_structure.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_structure.o: str_utils.hh annotated_ped_file.h write_structure_ext.h
write_structure.o: analysis.h
batch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
batch.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
batch.o: hwe.h compress_ext.h batch.h analysis.h reorder_loci_ext.h
batch.o: batch_input_ext.h error_messages_ext.h genetic_utils_ext.h
batch.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
batch.o: plink_ext.h utils_ext.h cw_routines_ext.h class_old.h
batch.o: write_plink_ext.h plink_core_ext.h
write_fbat.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_fbat.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_fbat.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_fbat.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_fbat.o: output_routines_ext.h linkage_ext.h sh_util.h fcmap_ext.h
write_fbat.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_fbat.o: input_ops.hh str_utils.hh annotated_ped_file.h write_fbat_ext.h
write_fbat.o: analysis.h fbat.R.h
write_pseq_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_pseq_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_pseq_files.o: batch_input.h input_check.h hwe.h compress_ext.h loop.h
write_pseq_files.o: file_ops.h error_messages_ext.h utils_ext.h entry.h
write_pseq_files.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_pseq_files.o: sh_util.h fcmap_ext.h omit_ped_ext.h
write_pseq_files.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_pseq_files.o: str_utils.hh annotated_ped_file.h plink_core_ext.h
write_pseq_files.o: analysis.h write_plink_ext.h write_pseq_ext.h
menu_value_missing.o: common.h types.hh typedefs.h errorno.h allutil.h
menu_value_missing.o: linkage.h linkage_external.h makeped.h pedtree.h
menu_value_missing.o: batch_input.h input_check.h hwe.h compress_ext.h
menu_value_missing.o: batch_input_ext.h error_messages_ext.h fcmap_ext.h
menu_value_missing.o: plink_ext.h user_input_ext.h input_ops.hh str_utils.hh
menu_value_missing.o: annotated_ped_file.h utils_ext.h class_old.h analysis.h
menu_value_missing.o: write_plink_ext.h plink_core_ext.h write_pseq_ext.h
read_impute.o: common.h types.hh tod.hh error_messages_ext.h typedefs.h
read_impute.o: errorno.h allutil.h linkage.h linkage_external.h makeped.h
read_impute.o: pedtree.h batch_input.h input_check.h hwe.h compress_ext.h
read_impute.o: fcmap_ext.h batch_input_ext.h mrecode.h mrecode_ext.h
read_impute.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
read_impute.o: reorder_loci_ext.h plink_ext.h annotated_ped_file_ext.h
read_impute.o: annotated_ped_file.h read_files_ext.h str_utils.hh
read_impute.o: input_ops.hh read_impute.hh
read_impute_bgen.o: common.h types.hh tod.hh error_messages_ext.h typedefs.h
read_impute_bgen.o: errorno.h allutil.h linkage.h linkage_external.h
read_impute_bgen.o: makeped.h pedtree.h batch_input.h input_check.h hwe.h
read_impute_bgen.o: compress_ext.h fcmap_ext.h batch_input_ext.h mrecode.h
read_impute_bgen.o: mrecode_ext.h vcftools/mega2_vcftools_interface.h
read_impute_bgen.o: annotated_ped_file.h reorder_loci_ext.h plink_ext.h
read_impute_bgen.o: read_files_ext.h annotated_ped_file_ext.h
read_impute_bgen.o: annotated_ped_file.h input_ops.hh str_utils.hh
read_impute_bgen.o: read_impute.hh read_impute_bgen.hh
read_impute_bgen2.o: common.h types.hh tod.hh error_messages_ext.h typedefs.h
read_impute_bgen2.o: errorno.h allutil.h linkage.h linkage_external.h
read_impute_bgen2.o: makeped.h pedtree.h batch_input.h input_check.h hwe.h
read_impute_bgen2.o: compress_ext.h fcmap_ext.h batch_input_ext.h mrecode.h
read_impute_bgen2.o: mrecode_ext.h vcftools/mega2_vcftools_interface.h
read_impute_bgen2.o: annotated_ped_file.h reorder_loci_ext.h plink_ext.h
read_impute_bgen2.o: read_files_ext.h annotated_ped_file_ext.h
read_impute_bgen2.o: annotated_ped_file.h input_ops.hh str_utils.hh
read_impute_bgen2.o: read_impute.hh read_impute_bgen2.hh
write_shapeit_files.o: common.h types.hh typedefs.h errorno.h allutil.h
write_shapeit_files.o: linkage.h linkage_external.h makeped.h pedtree.h
write_shapeit_files.o: batch_input.h input_check.h hwe.h compress_ext.h
write_shapeit_files.o: loop.h file_ops.h error_messages_ext.h utils_ext.h
write_shapeit_files.o: entry.h create_summary_ext.h output_routines_ext.h
write_shapeit_files.o: linkage_ext.h sh_util.h batch_input_ext.h fcmap_ext.h
write_shapeit_files.o: omit_ped_ext.h output_file_names_ext.h
write_shapeit_files.o: user_input_ext.h input_ops.hh str_utils.hh
write_shapeit_files.o: annotated_ped_file.h write_shapeit_ext.h
write_shapeit_files.o: write_plink_ext.h plink_core_ext.h analysis.h
write_shapeit_files.o: write_shapeit_args.h
write_roadtrips.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_roadtrips.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_roadtrips.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_roadtrips.o: error_messages_ext.h utils_ext.h entry.h
write_roadtrips.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_roadtrips.o: sh_util.h fcmap_ext.h genetic_utils_ext.h omit_ped_ext.h
write_roadtrips.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_roadtrips.o: str_utils.hh annotated_ped_file.h write_roadtrips_ext.h
write_roadtrips.o: analysis.h
write_mach.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_mach.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_mach.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_mach.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_mach.o: output_routines_ext.h linkage_ext.h sh_util.h fcmap_ext.h
write_mach.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_mach.o: input_ops.hh str_utils.hh annotated_ped_file.h write_mach_ext.h
write_mach.o: analysis.h write_plink_ext.h plink_core_ext.h
write_minimacshapeit.o: common.h types.hh typedefs.h errorno.h allutil.h
write_minimacshapeit.o: linkage.h linkage_external.h makeped.h pedtree.h
write_minimacshapeit.o: batch_input.h input_check.h hwe.h compress_ext.h
write_minimacshapeit.o: loop.h file_ops.h error_messages_ext.h utils_ext.h
write_minimacshapeit.o: entry.h create_summary_ext.h output_routines_ext.h
write_minimacshapeit.o: linkage_ext.h sh_util.h fcmap_ext.h omit_ped_ext.h
write_minimacshapeit.o: output_file_names_ext.h user_input_ext.h input_ops.hh
write_minimacshapeit.o: str_utils.hh annotated_ped_file.h
write_minimacshapeit.o: write_minimacshapeit_ext.h analysis.h
write_minimacshapeit.o: write_plink_ext.h plink_core_ext.h
write_vcf.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
write_vcf.o: linkage_external.h makeped.h pedtree.h batch_input.h
write_vcf.o: input_check.h hwe.h compress_ext.h loop.h file_ops.h
write_vcf.o: error_messages_ext.h utils_ext.h entry.h create_summary_ext.h
write_vcf.o: output_routines_ext.h linkage_ext.h sh_util.h fcmap_ext.h
write_vcf.o: omit_ped_ext.h output_file_names_ext.h user_input_ext.h
write_vcf.o: input_ops.hh str_utils.hh annotated_ped_file.h write_files_ext.h
write_vcf.o: read_files_ext.h write_vcf_ext.h analysis.h zlib-1.2.8/zlib.h
write_vcf.o: zlib-1.2.8/zconf.h dbrefallele.h dblite.hh
write_vcf.o: mega2_bcftools_interface.h lib/bcftools-1.6/filter.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_defs.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_log.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/vcf.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/kstring.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_endian.h
write_vcf.o: ./lib/bcftools-1.6/htslib-1.6/htslib/tbx.h
write_vcf.o: lib/bcftools-1.6/vcfview.h
write_vcf.o: lib/bcftools-1.6/htslib-1.6/htslib/vcf.h
dbmega2.o: types.hh common.h typedefs.h errorno.h allutil.h linkage.h
dbmega2.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
dbmega2.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h
dbmega2.o: write_files_ext.h tod.hh dblite.hh dbbatch.hh dbmisc.hh
dbmega2.o: dbpedigree.hh dballele.hh dbloci.hh dbmap.hh dbgenotype.hh
dbmega2.o: dbrefallele.h
dbmisc.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbmisc.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
dbmisc.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h input_ops.hh
dbmisc.o: str_utils.hh annotated_ped_file.h dblite.hh dbmisc.hh tod.hh
dbbatch.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbbatch.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
dbbatch.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h input_ops.hh
dbbatch.o: str_utils.hh annotated_ped_file.h dblite.hh dbbatch.hh tod.hh
dbpedigree.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbpedigree.o: linkage_external.h makeped.h pedtree.h batch_input.h
dbpedigree.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
dbpedigree.o: utils_ext.h dblite.hh dbpedigree.hh tod.hh
dballele.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dballele.o: linkage_external.h makeped.h pedtree.h batch_input.h
dballele.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
dballele.o: utils_ext.h dblite.hh dbmisc.hh dballele.hh tod.hh
dbloci.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbloci.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
dbloci.o: hwe.h compress_ext.h read_files_ext.h error_messages_ext.h
dbloci.o: utils_ext.h dblite.hh dbloci.hh dbmisc.hh tod.hh
dbmap.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbmap.o: linkage_external.h makeped.h pedtree.h batch_input.h input_check.h
dbmap.o: hwe.h compress_ext.h error_messages_ext.h utils_ext.h dblite.hh
dbmap.o: dbmisc.hh dbmap.hh tod.hh
dbgenotype.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbgenotype.o: linkage_external.h makeped.h pedtree.h batch_input.h
dbgenotype.o: input_check.h hwe.h compress_ext.h error_messages_ext.h
dbgenotype.o: utils_ext.h dblite.hh dbgenotype.hh tod.hh
dbrefallele.o: common.h types.hh typedefs.h errorno.h allutil.h linkage.h
dbrefallele.o: linkage_external.h makeped.h pedtree.h batch_input.h
dbrefallele.o: input_check.h hwe.h compress_ext.h dbrefallele.h dblite.hh
dbrefallele.o: fcmap_ext.h user_input_ext.h input_ops.hh str_utils.hh
dbrefallele.o: annotated_ped_file.h error_messages_ext.h read_files_ext.h
dbrefallele.o: zlib-1.2.8/zlib.h zlib-1.2.8/zconf.h
read_bcfs.o: common.h types.hh tod.hh error_messages_ext.h
read_bcfs.o: lib/bcftools-1.6/filter.h lib/bcftools-1.6/vcfview.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_defs.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_log.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/vcf.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/kstring.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_endian.h
read_bcfs.o: ./lib/bcftools-1.6/htslib-1.6/htslib/tbx.h typedefs.h errorno.h
read_bcfs.o: allutil.h linkage.h linkage_external.h makeped.h pedtree.h
read_bcfs.o: batch_input.h input_check.h hwe.h compress_ext.h fcmap_ext.h
read_bcfs.o: batch_input_ext.h mrecode.h mrecode_ext.h
read_bcfs.o: vcftools/mega2_vcftools_interface.h annotated_ped_file.h
read_bcfs.o: reorder_loci_ext.h plink_ext.h annotated_ped_file_ext.h
read_bcfs.o: annotated_ped_file.h read_files_ext.h utils_ext.h
read_bcfs.o: phe_lookup_ext.h str_utils.hh read_bcfs.h input_ops.hh
read_bcfs.o: mega2_bcftools_interface.h
mega2_bcftools_interface.o: common.h types.hh typedefs.h errorno.h allutil.h
mega2_bcftools_interface.o: linkage.h linkage_external.h makeped.h pedtree.h
mega2_bcftools_interface.o: batch_input.h input_check.h hwe.h compress_ext.h
mega2_bcftools_interface.o: tod.hh utils_ext.h genetic_utils_ext.h
mega2_bcftools_interface.o: error_messages_ext.h annotated_ped_file.h
mega2_bcftools_interface.o: read_files_ext.h marker_lookup_ext.h
mega2_bcftools_interface.o: phe_lookup_ext.h reorder_loci_ext.h
mega2_bcftools_interface.o: lib/bcftools-1.6/filter.h
mega2_bcftools_interface.o: lib/bcftools-1.6/vcfview.h
mega2_bcftools_interface.o: lib/bcftools-1.6/htslib-1.6/htslib/khash_str2int.h
mega2_bcftools_interface.o: lib/bcftools-1.6/htslib-1.6/htslib/hts.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_defs.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_log.h
mega2_bcftools_interface.o: lib/bcftools-1.6/htslib-1.6/htslib/vcf.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/vcf.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/kstring.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/hts_endian.h
mega2_bcftools_interface.o: ./lib/bcftools-1.6/htslib-1.6/htslib/tbx.h
mega2_bcftools_interface.o: lib/bcftools-1.6/htslib-1.6/htslib/vcfutils.h
mega2_bcftools_interface.o: mega2_bcftools_interface.h
write_MQLS.o: write_MQLS_ext.h analysis.h common.h types.hh typedefs.h
write_MQLS.o: errorno.h allutil.h linkage.h linkage_external.h makeped.h
write_MQLS.o: pedtree.h batch_input.h input_check.h hwe.h compress_ext.h
write_MQLS.o: loop.h file_ops.h error_messages_ext.h utils_ext.h entry.h
write_MQLS.o: create_summary_ext.h output_routines_ext.h linkage_ext.h
write_MQLS.o: sh_util.h fcmap_ext.h omit_ped_ext.h output_file_names_ext.h
write_MQLS.o: user_input_ext.h input_ops.hh str_utils.hh annotated_ped_file.h
write_MQLS.o: write_files_ext.h
