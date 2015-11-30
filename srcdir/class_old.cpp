/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  Daniel E. Weeks, and University of Pittsburgh

  This file is part of the Mega2 program, which is free software; you
  can redistribute it and/or modify it under the terms of the GNU
  General Public License as published by the Free Software Foundation;
  either version 3 of the License, or (at your option) any later
  version.

  Mega2 is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

*/

#include <stdio.h>
#include "common.h"
#include "typedefs.h"
#include "analysis.h"

#include "ctype.h"
#include "utils_ext.h"

#include <string.h>

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "output_file_names_ext.h"

#include "aspex_ext.h"
#include "cranefoot_ext.h"
#include "create_summary_ext.h"
#include "hwe_user_input_ext.h"
#include "makenucs_ext.h"
#include "mega2annot_ext.h"
#include "scripts_ext.h"
#include "slink_ext.h"
#include "splink_ext.h"
#include "write_IQLS_ext.h"
#include "write_SUP_files_ext.h"
#include "write_apm_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
#include "write_loki_ext.h"
#include "write_mfiles_ext.h"
#include "write_pap_ext.h"
#include "write_premakeped_ext.h"
#include "write_prest_ext.h"
#include "write_sage_files_ext.h"
#include "write_simulate_files_ext.h"
#include "write_solar_files_ext.h"
#include "write_vitesse_ext.h"

#include "class_old.h"


void CLASS_ALLEGRO::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
        
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_ALLEGRO::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "al_ped.%s", num);
    sprintf(file_names[1], "al_dat.%s", num);
    sprintf(file_names[2], "al_in.%s", num);
    sprintf(file_names[3], "allegro.%s.sh", num);
}

void CLASS_ALLEGRO::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_ALLELE_FREQ::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_ALLELE_FREQ::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "freq.%s", num);
}

void CLASS_ALLELE_FREQ::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_APM::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    char *pedfl_name = NULL, *locusfl_name = NULL, *omitfl_name = NULL; // assignment to supress warnings
    int disease_locus = 0;
    create_APM_file(pedfl_name, locusfl_name, disease_locus,
                    "2", LPedTreeTop, PedTreeTop, &infl_type, &outfl_type,
                    numchr, mapfl_name, file_names,
                    analysis, omitfl_name, UntypedPedOpt);
}

void CLASS_APM::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "kin_ml.%s", num);
}

void CLASS_APM::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_APM_MULT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    char *pedfl_name = NULL, *locusfl_name = NULL; // assignment to supress warnings
    int disease_locus = 0;

    create_APMULT(pedfl_name, locusfl_name, disease_locus,
                  "2", &LPedTreeTop, PedTreeTop, &infl_type, &outfl_type,
                  numchr, mapfl_name, analysis, file_names,
                  UntypedPedOpt);
}

void CLASS_APM_MULT::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "kin_mult");
    sprintf(file_names[3], "apmmult.%2s.sh", num);
    sprintf(file_names[14], "apmmult_sum.%2s", num);
}

void CLASS_APM_MULT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[14], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_ASPEX::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "Sib-ibd");      break;
    case 2:  strcpy(subprog, "Sib-tdt");      break;
    case 3:  strcpy(subprog, "Sib-phase");    break;
    case 4:  strcpy(subprog, "Sib-map");      break;
    default:                                  break;
    }
}

void CLASS_ASPEX::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (*sub_prog_name == 0) {
        (*analysis)->_suboption = 0;
        return;
    }
    switch(tolower((unsigned char)sub_prog_name[strlen(sub_prog_name) - 1])) {
    case 'd':
        (*analysis)->_suboption = 1; break;
    case 't':
        (*analysis)->_suboption = 2; break;
    case 'e':
        (*analysis)->_suboption = 3; break;
    case 'p':
        (*analysis)->_suboption = 4; break;
    default :
        (*analysis)->_suboption = 0; break;
    }
}

void CLASS_ASPEX::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_aspex_files(&LPedTreeTop,
                      *numchr, file_names, UntypedPedOpt, Top2);
}

void CLASS_ASPEX::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[4], "asp_in.%s", num);
    sprintf(file_names[5], "asp_dat.%s", num);
}

void CLASS_ASPEX::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_CRANEFOOT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_CRANEFOOT_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_CRANEFOOT::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "crnft_ped.%s", num);
    sprintf(file_names[1], "crnft_control.%s", num);
    sprintf(file_names[2], "crnft_shell.%s.sh", num);
}

void CLASS_CRANEFOOT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_CREATE_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_CREATE_SUMMARY::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "seg_sum.%s", num);
    sprintf(file_names[1], "cnt_sum.%s", num);
    sprintf(file_names[2], "sib_sum.%s", num);
}

void CLASS_CREATE_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
}


void CLASS_GENEHUNTER::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GENEHUNTER::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "gh_ped.%s", num);
    sprintf(file_names[1], "gh_dat.%s", num);
    sprintf(file_names[2], "gh_in.%s", num);
    sprintf(file_names[3], "gh.%s.sh", num);
}

void CLASS_GENEHUNTER::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_GENEHUNTERPLUS::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GENEHUNTERPLUS::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "ghp_ped.%s", num);
    sprintf(file_names[1], "ghp_dat.%s", num);
    sprintf(file_names[2], "ghp_in.%s", num);
    sprintf(file_names[3], "ghp.%s.sh", num);
}

void CLASS_GENEHUNTERPLUS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_GENOTYPING_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_GENOTYPING_SUMMARY::gen_file_names(char *file_names[], char *num) {
    strcpy(file_names[0], "genotyping_rate");
}

void CLASS_GENOTYPING_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_GHMLB::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GHMLB::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "mlb_ped.%s", num);
    sprintf(file_names[1], "mlb_dat.%s", num);
    sprintf(file_names[2], "mlb_in.%s", num);
    sprintf(file_names[3], "mlb.%s.sh", num);
}

void CLASS_GHMLB::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_HAPLOTYPE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_HAPLOTYPE::gen_file_names(char *file_names[], char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_haplo.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_HAPLOTYPE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}


void CLASS_HWETEST::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "Gen");      break;
    case 2:  strcpy(subprog, "HWE");      break;
    case 3:  strcpy(subprog, "Chi-sq");   break;
    case 4:  strcpy(subprog, "Exact");    break;
    case 5:  strcpy(subprog, "Mendel");   break;
    default:                              break;
    }
}

void CLASS_HWETEST::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'g':
        (*analysis)->_suboption = 1; break;
    case 'h':
        (*analysis)->_suboption = 2; break;
    case 'c':
        (*analysis)->_suboption = 3; break;
    case 'e':
        (*analysis)->_suboption = 4; break;
    case 'm':
        (*analysis)->_suboption = 5; break;
    default:
        (*analysis)->_suboption = 0; break;
    }
}

void CLASS_HWETEST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    hwe_user_input(LPedTreeTop, numchr, file_names);
}

void CLASS_HWETEST::gen_file_names(char *file_names[], char *num) {
}

void CLASS_HWETEST::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_IBD_EST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_IBD_EST::gen_file_names(char *file_names[], char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_ibd.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}


void CLASS_IBD_EST::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_IQLS::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    int format_checksum = 1; //XX
    if (format_checksum == 0) {
        errorvf("Pedigree, names and map file appear to be in LINKAGE format.\n");
        mssgvf("IQLS requires the input data be in Mega2 format.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    create_IQLS_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_IQLS::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "IQLS_pedigree.%s", num);
    sprintf(file_names[1], "IQLS_marker.%s", num);
    sprintf(file_names[2], "IQLS_parameter.%s", num);
    sprintf(file_names[3], "IQLS.%s.sh", num);
    sprintf(file_names[4], "Idcoefs_pedigree.%s", num);
    sprintf(file_names[5], "Idcoefs_study.%s", num);
}

void CLASS_IQLS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_LIABLE_FREQ::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_LIABLE_FREQ::gen_file_names(char *file_names[], char *num) {
    strcpy(file_names[0], "liabililty_freq");
}

void CLASS_LIABLE_FREQ::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_LINKAGE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_linkage_files(&LPedTreeTop, numchr, file_names, UntypedPedOpt);
}

void CLASS_LINKAGE::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "Lpedin.%s", num);
    sprintf(file_names[1], "Ldatain.%s", num);
}

void CLASS_LINKAGE::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_LOCATION::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_LOCATION::gen_file_names(char *file_names[], char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_loc.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_LOCATION::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}
void CLASS_LOD2::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_LOD2_format_files(&LPedTreeTop, PedTreeTop,
                             analysis, numchr, file_names,
                             UntypedPedOpt);
}

void CLASS_LOD2::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "lod2ped.%s", num);
    sprintf(file_names[1], "lod2data.%s", num);
    sprintf(file_names[3], "lod2.%s.sh", num);
}

void CLASS_LOD2::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_LOKI::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_LOKI_files(&LPedTreeTop, numchr, &infl_type, file_names,
                      UntypedPedOpt, *analysis);
}

void CLASS_LOKI::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "Loki_ped.%s", num);
    sprintf(file_names[1], "Loki_freq.%s", num);
    sprintf(file_names[2], "Loki_map.%s", num);
    sprintf(file_names[3], "Loki_locus.%s", num);
    sprintf(file_names[4], "Loki_link.%s", num);
    sprintf(file_names[5], "Loki_control.%s", num);
    sprintf(file_names[6], "Loki_param.%s", num);
    sprintf(file_names[7], "Loki.%s.sh", num);
}

void CLASS_LOKI::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
}

void CLASS_MEGA2ANNOT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_mega2annot_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_MEGA2ANNOT::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "pedin.%s.mega2", num);
    sprintf(file_names[1], "names.%s.mega2", num);
    sprintf(file_names[2], "map.%s.mega2", num);
    sprintf(file_names[3], "frequency.%s.mega2", num);
    sprintf(file_names[4], "penetrance.%s.mega2", num);
}

void CLASS_MEGA2ANNOT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_MENDEL4::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_mega2annot_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_MENDEL4::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[4], "mendel_map.%s", num);
    sprintf(file_names[0], "mendel_ped.%s", num);
    sprintf(file_names[1], "mendel_locus.%s", num);
    sprintf(file_names[6], "mendel_pen.%s", num);
    sprintf(file_names[3], "mendel_control.%s", num);
    sprintf(file_names[7], "mendel_variable.%s", num);
}

void CLASS_MENDEL4::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[7], numchr);
}

void CLASS_MENDEL7_CSV::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MENDEL7_CSV::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[4], "mendel_map.%s", num);
    sprintf(file_names[0], "mendel_ped.%s", num);
    sprintf(file_names[1], "mendel_locus.%s", num);
    sprintf(file_names[6], "mendel_pen.%s", num);
    sprintf(file_names[3], "mendel_control.%s", num);
    sprintf(file_names[7], "mendel_variable.%s", num);
}

void CLASS_MENDEL7_CSV::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[7], numchr);
}

void CLASS_MENDEL::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MENDEL::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "pedm.%s", num);
    sprintf(file_names[1], "locus.%s", num);
    sprintf(file_names[6], "pen.%s", num);
    sprintf(file_names[2], "batch.%s", num);
    sprintf(file_names[3], "m13bat.%s", num);
}

void CLASS_MENDEL::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_MERLIN::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_merlin_files(&LPedTreeTop, numchr, &infl_type,
                        file_names, UntypedPedOpt, *analysis);
    mssgf("------------------------------------");
    create_mendel_file(&LPedTreeTop, mapfl_name,
                       PedTreeTop, &infl_type, &outfl_type,
                       numchr, analysis, file_names,
                       UntypedPedOpt);
}

void CLASS_MERLIN::gen_file_names(char *file_names[], char *num) {
    /* set simwalk2 file names */
    sprintf(file_names[0], "PEDIGREE.%s", num);
    sprintf(file_names[1], "LOCUS.%s", num);
    sprintf(file_names[2], "BATCH2.%s", num);
    sprintf(file_names[3], "npl.%s.sh", num);
    sprintf(file_names[6], "PEN.%s", num);
    sprintf(file_names[4], "SW2_MAP.%s", num);
    /* Merlin file names */
    sprintf(file_names[9], "sw2merlin_freq.%s", num);
    sprintf(file_names[10], "sw2merlin_ped.%s", num);
    sprintf(file_names[11], "sw2merlin_data.%s", num);
    sprintf(file_names[12], "sw2merlin_map.%s", num);
    sprintf(file_names[13], "sw2merlin_order.%s", num);
    strcpy(file_names[14], "merlin2sw2.pl");
}

void CLASS_MERLIN::replace_chr_number(char *file_names[], int numchr) {
    /* simwalk2 files */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[4], numchr);
    /* Merlin files */
    change_output_chr(file_names[9], numchr);
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[12], numchr);
    change_output_chr(file_names[13], numchr);
}

void CLASS_MERLINONLY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_merlin_files(&LPedTreeTop, numchr, &infl_type,
                        file_names, UntypedPedOpt, *analysis);
}

void CLASS_MERLINONLY::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[9], "merlin_freq.%s", num);
    sprintf(file_names[10], "merlin_ped.%s", num);
    sprintf(file_names[11], "merlin_data.%s", num);
    sprintf(file_names[12], "merlin_map.%s", num);
    sprintf(file_names[13], "merlin.%s.sh", num);
    sprintf(file_names[14], "merlin_model");
}

void CLASS_MERLINONLY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[9], numchr);
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[12], numchr);
    change_output_chr(file_names[13], numchr);
}

void CLASS_MISTYPING::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MISTYPING::gen_file_names(char *file_names[], char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_mis.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_MISTYPING::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_NONPARAMETRIC::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_NONPARAMETRIC::gen_file_names(char *file_names[], char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_npl.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_NONPARAMETRIC::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_NUKE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_nuclear_families(&LPedTreeTop, file_names[0], file_names[1],
                            analysis, &infl_type, UntypedPedOpt, Top2);
}

void CLASS_NUKE::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[1], "nuke_data.%s", num);
    sprintf(file_names[0], "nuke_ped.%s", num);
}

void CLASS_NUKE::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_PAP::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_pap_files(&LPedTreeTop, file_names, UntypedPedOpt, numchr, &infl_type);
}

void CLASS_PAP::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "trip.%s", num);
    sprintf(file_names[1], "header.%s", num);
    sprintf(file_names[2], "phen.%s", num);
    sprintf(file_names[3], "popln.%s", num);
    sprintf(file_names[4], "pap.%s.sh", num);
}

void CLASS_PAP::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_PREMAKEPED::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_premakeped_files(&LPedTreeTop, numchr, &infl_type,
                            file_names, UntypedPedOpt);
}

void CLASS_PREMAKEPED::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "Ppedin.%s", num);
    sprintf(file_names[1], "Pdatain.%s", num);
}

void CLASS_PREMAKEPED::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_PREST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_prest_files(LPedTreeTop, *analysis,
                       file_names, UntypedPedOpt, numchr,
                       &infl_type);
}

void CLASS_PREST::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "prest_ped.%s", num);
    sprintf(file_names[1], "prest_loc.%s", num);
    sprintf(file_names[2], "prest_geno.%s", num);
    sprintf(file_names[3], "prest_chrom.%s", num);
    sprintf(file_names[4], "prest.%s.sh", num);
}

void CLASS_PREST::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_QUANT_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_QUANT_SUMMARY::gen_file_names(char *file_names[], char *num) {
    strcpy(file_names[0], "phenotyping_rate");
}

void CLASS_QUANT_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_SAGE4::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SAGE4_file(&LPedTreeTop, PedTreeTop, &infl_type,
                      numchr, analysis, file_names, UntypedPedOpt);
}

void CLASS_SAGE4::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "sage4_ped.%s", num);
    sprintf(file_names[1], "sage4_dat.%s", num);
    sprintf(file_names[6], "sage4_par.%s", num);
    sprintf(file_names[7], "sage4_map.%s", num);
}

void CLASS_SAGE4::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
}

void CLASS_SAGE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
//  char *mapfl_name = NULL; // assignment to supress warning
    int trait_locus_first;

    create_SAGE_file(&LPedTreeTop, mapfl_name, &trait_locus_first,
                     PedTreeTop, &infl_type, &outfl_type,
                     numchr, analysis, file_names, UntypedPedOpt);
}

void CLASS_SAGE::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "sage_ped.%s", num);
    sprintf(file_names[1], "sage_loc.%s", num);
    sprintf(file_names[2], "sage_par.%s", num);
    sprintf(file_names[3], "sage.%s.sh", num);
    sprintf(file_names[6], "sage_cnt.%s", num);
    sprintf(file_names[7], "sage_sibpal.%s", num);
    sprintf(file_names[10], "sage_cntpar.%s", num);
    sprintf(file_names[11], "sage_cnt.%s.sh", num);
}

void CLASS_SAGE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[0], numchr);
}

void CLASS_SIMULATE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SIMULATE_format_files(&LPedTreeTop, PedTreeTop, analysis,
                                 &infl_type, &outfl_type,
                                 numchr, mapfl_name, file_names,
                                 UntypedPedOpt);
}

void CLASS_SIMULATE::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "simped.%s", num);
    sprintf(file_names[1], "simdata.%s", num);
    sprintf(file_names[2], "problem.%s", num);
    sprintf(file_names[3], "simulate.%s.sh", num);
}

void CLASS_SIMULATE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_SIMWALK2::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int  selection;
    char select[10];

    if (batchANALYSIS) {
        selection = (*analysis)->_suboption;
    } else {
        selection=0;
        while (selection < 1 || selection > 5) {
            draw_line();
            printf("SimWalk2 program options:\n");
            printf("1) Haplotype analysis\n");
            printf("2) Parametric Linkage analysis\n");
            printf("3) Non-Parametric Linkage analysis\n");
            printf("4) IBD analysis\n");
            printf("5) Mistyping analysis\n");
            printf("Enter selection 1-5 > ");
            fcmap(stdin, "%s", select); newline;
            selection = 0;
            sscanf(select, "%d", &selection);
            if (selection < 1 || selection > 5)
                warn_unknown(select);
        }
        draw_line();
    }
    switch (selection) {
    case 1:
        *analysis = HAPLOTYPE;
        break;
    case 2:
        *analysis = LOCATION;
        break;
    case 3:
        *analysis = NONPARAMETRIC;
        break;
    case 4:
        *analysis = IBD_EST;
        break;
    case 5:
        *analysis = MISTYPING;
        break;
    }
}

void CLASS_SIMWALK2::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'h':
        *analysis = HAPLOTYPE;
        break;
    case 'i':
        *analysis = IBD_EST;
        break;
    case 'n':
        *analysis = NONPARAMETRIC;
        break;
    case 'p':
        *analysis = LOCATION;
        break;
    case 'm':
        *analysis = MISTYPING;
        break;
    default:
        *analysis = NULL;
        break;
    }
}

void CLASS_SLINK::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SLINK_format_files(LPedTreeTop, PedTreeTop, analysis,
                              numchr, file_names,
                              UntypedPedOpt);
}

void CLASS_SLINK::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "simped.%s", num);
    sprintf(file_names[1], "simdata.%s", num);
    sprintf(file_names[2], "slinkin.%s", num);
    sprintf(file_names[3], "slink.%s.sh", num);
}

void CLASS_SLINK::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_SOLAR::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SOLAR_files(&LPedTreeTop, PedTreeTop, analysis,
                       &infl_type, &outfl_type,
                       numchr, mapfl_name, file_names, UntypedPedOpt);
}

void CLASS_SOLAR::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "solar_ped.%s", num);
    sprintf(file_names[1], "solar_freq.%s", num);
    sprintf(file_names[6], "solar_phen.%s", num);
    sprintf(file_names[8], "solar_marker.%s", num);
    sprintf(file_names[9], "solar_map.%s", num);
    sprintf(file_names[10],"solar_load.%s.tcl", num);
}

void CLASS_SOLAR::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[8], numchr);
    change_output_chr(file_names[9], numchr);
    /*    change_output_chr(file_names[10], numchr); */
}

void CLASS_SPLINK::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SPLINK(&LPedTreeTop, analysis, numchr,
                  file_names, UntypedPedOpt);
}

void CLASS_SPLINK::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "splink_ped.%s", num);
    sprintf(file_names[3], "splink.%s.sh", num);
}

void CLASS_SPLINK::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_SUP::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SUP_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_SUMMARY::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int  selection;
    char select[10];

    if (batchANALYSIS) {
        selection = (*analysis)->_suboption;
    } else {
        selection = 0;
        while(!(selection == 1 ||
                selection == 2 ||
                selection == 3 ||
                selection == 4 ||
                selection == 5)) {
            draw_line();
            printf("Selection Menu: Summary file options\n");
            printf("1) Create segregation and relative count summary files.\n");
            printf("2) Create allele frequency summary table.\n");
            printf("3) Count alleles and genotypes within groups.\n");
            printf("4) Create genotyping success rate summary.\n");
            printf("5) Create quantitative phenotype summary.\n");
            printf("Enter selection: 1 - 5 > ");
            fcmap(stdin,"%s", select); newline;
            selection=0;
            sscanf(select, "%d", &selection);
            if (selection < 1 || selection > 5)
                warn_unknown(select);
        }
    }

    if (selection == 1) {
        *analysis = CREATE_SUMMARY;
    } else if (selection  == 2) {
        *analysis = TO_ALLELE_FREQ;
    } else if (selection == 3 ) {
        *analysis = TO_LIABLE_FREQ;
    } else if (selection == 4 ) {
        *analysis = GENOTYPING_SUMMARY;
    } else if (selection == 5 ) {
        *analysis = QUANT_SUMMARY;
    }
}

void CLASS_SUMMARY::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'a':
        *analysis = TO_ALLELE_FREQ;
        break;
    case 's':
        *analysis = CREATE_SUMMARY;
        break;
    case 'l':
    case 'c':
        *analysis = TO_LIABLE_FREQ;
        break;
    case 'g':
        *analysis = GENOTYPING_SUMMARY;
        break;
    case 'q':
        *analysis = QUANT_SUMMARY;
        break;
    default:
        *analysis = NULL;
	break;
    }
}

void CLASS_SUP::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "sup_simped.%s", num);
    sprintf(file_names[1], "sup_simdata.%s", num);
    sprintf(file_names[2], "sup_locus.%s", num);
    sprintf(file_names[3], "slinkin.%s", num);
    sprintf(file_names[4], "sup.%s.sh", num);
    sprintf(file_names[5], "sup_outpeds.%s", num);
    sprintf(file_names[6], "sup_mega2_locus.%s", num);
    sprintf(file_names[7], "sup_mega2_map.%s", num);
    sprintf(file_names[8], "sup_mega2_batch.%s", num);
}

void CLASS_SUP::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_TDTMAX::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_TDTMAX(&LPedTreeTop, analysis, numchr, file_names,
                  UntypedPedOpt, Top2);
}


void CLASS_TDTMAX::gen_file_names(char *file_names[], char *num) {
    strcpy(file_names[0], "tdtmax_data");
    sprintf(file_names[3], "tdtmax.%s.sh", num);
}

void CLASS_TDTMAX::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_VITESSE::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "LINKMAP");      break;
    case 2:  strcpy(subprog, "MLINK");        break;
    default:                                  break;
    }
}

void CLASS_VITESSE::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (tolower((unsigned char)sub_prog_name[0]) == 'l') {
        (*analysis)->_suboption = 1;
    } else if (tolower((unsigned char)sub_prog_name[0]) == 'm') {
        (*analysis)->_suboption = 2;
    }
}

void CLASS_VITESSE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_vitesse_files(&LPedTreeTop, numchr, file_names, UntypedPedOpt);
}

void CLASS_VITESSE::gen_file_names(char *file_names[], char *num) {
    sprintf(file_names[0], "vpedin.%s", num);
    sprintf(file_names[1], "vdatain.%s", num);
    sprintf(file_names[3], "vitesse.%s.sh", num);
    sprintf(file_names[4], "vout.%s", num);
    sprintf(file_names[5], "vstrm.%s", num);
}

void CLASS_VITESSE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}
