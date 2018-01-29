/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

===========================================================================
*/

//NrvB: ped_top->Locus ONLY contains markers and is a copy of the Top->LocusTop->Locus data !!

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"
#include "class_old.h"
#include "tod.hh"

#include "input_check.h"

#include "check_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "user_input_ext.h"
#include "input_check_ext.h"
#include "pedtree_ext.h"
#include "utils_ext.h"
#include "write_pap_ext.h"
#include "batch_input_ext.h"
/*
              check_ext.h:  check_locus check_ped_data check_ped_relations check_unique_ids
     create_summary_ext.h:  count_alleles
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  exclaim
        grow_string_ext.h:  grow
        input_check_ext.h:  clear_ped_status mito_transmission_report
            pedtree_ext.h:  renumber_ped
              utils_ext.h:  EXIT draw_line imax log_line summary_time_stamp
          write_pap_ext.h:  set_pap_unique_id
*/


int Display_Errors, Display_Messages;
extern int lastautosome;
char          *zero;
/* prototype definitions */
void            full_check(ped_top *Top, linkage_ped_top *LPedTop,
			   analysis_type analysis);
int             UnMappedLociCheck(linkage_ped_top *Top);
int             input_observed_freq_check(linkage_ped_top *Top, double threshold);
void            check_renumber_ped(ped_tree *Ped);
/*
static  void    write_reset_summary(ped_top *PTop, loc_list *invalidp,
				    person_loc_list *halftypedp,
				    person_loc_list *outofboundsp,
				    int imend, int ht, int aexceed,
				    int uniqueids);
*/
void     create_unique_ids(linkage_ped_top *Top, analysis_type analysis);

/* end of prototypes */

/* Compute the sum of squared differences between input and output alleles
   for all numbered loci, warn if this is above the thershold value
*/
int  input_observed_freq_check(linkage_ped_top *Top, double threshold)
{

    int i, all, num_alleles, found_problem=0;
    double *sum_squared, diff;
    allele_freq_struct ar;
    int first_time=1, flush=0;

    sum_squared = CALLOC((size_t) Top->LocusTop->LocusCnt, double);

    SECTION_LOG_INIT(freq_check);
    for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
        if (Top->LocusTop->Locus[i].Marker->chromosome == MISSING_CHROMO) continue;
        sum_squared[i] = 0.0;
        if (Top->LocusTop->Locus[i].Type == NUMBERED ||
            Top->LocusTop->Locus[i].Type == BINARY) {

            ar.TotalAlleles = 0;
            ar.TotalPeople = 0;
            ar.HalfTyped = 0;
            ar.TotalKnownAlleles = 0;

            num_alleles = Top->LocusTop->Locus[i].AlleleCnt;
            ar.Allele_Bin = CALLOC((size_t) num_alleles, int);
            ar.allele_freq = CALLOC((size_t) num_alleles, double);
            for (all = 0; all < num_alleles; all++) {
                ar.Allele_Bin[all]=0;
            }
            count_alleles(Top, &ar, i, 1, 0);
            sum_squared[i]=0.0;

            for (all = 0; all < num_alleles; all++) {
                ar.allele_freq[all] =
                    ((double) (ar.Allele_Bin[all])) / ((double) (ar.TotalKnownAlleles));
                diff = ar.allele_freq[all] - Top->LocusTop->Locus[i].Allele[all].Frequency;
                sum_squared[i] += diff * diff;
            }

            if (sum_squared[i] > threshold) {
                found_problem=1;

                if (first_time) {
                    /* First error */
                    SECTION_LOG(freq_check);
                    log_line(mssgf);
                    mssgf("For these markers, sum of squared differences between ");
                    sprintf(err_msg,
                            "input and observed allele frequencies exceeded threshold %9.7f.",
                            threshold);
                    mssgf(err_msg);
                    SECTION_LOG_HEADER(freq_check);
                    strcpy(err_msg, "  ");
                    first_time = 0;
                }
                grow(err_msg, "%s: SSD=%9.7f ",
                     Top->LocusTop->Locus[i].LocusName, sum_squared[i]);

                flush++;
                if (flush == 2) {
                    SECTION_LOG(freq_check);
                    mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
                }
            }
            free(ar.Allele_Bin); ar.Allele_Bin= NULL;
            free(ar.allele_freq); ar.allele_freq=NULL;
        }
    }

    if (flush) {
        SECTION_LOG(freq_check);
        mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
    }

    if (SECTION_ERR_COUNT(freq_check) > 0) {
        printf("Check error logs for full list of reset genotypes.\n");
        draw_line();
    }
    SECTION_LOG_FINI(freq_check);

    free(sum_squared);

    return found_problem;
}


/*---------------------------------------------------------------+
  | Given the pedtree Top, check all the pedigree and locus data. |
  | Print a summary if an error was detected.                     |
  |                                                               |
  | infl_type is used to determine if sex errors should be        |
  | ignored (since APM files do not contain sex data).            |
  +---------------------------------------------------------------*/

SECTION_ERR_INIT(check_locus);
SECTION_ERR_INIT(biallele);
SECTION_ERR_INIT(check_ped_relations);
SECTION_ERR_INIT(check_half_type);
SECTION_ERR_INIT(check_oob);
SECTION_ERR_INIT(check_inheritance);
SECTION_ERR_INIT(check_sibship_alleles);

int      set_uniq_check(linkage_ped_top *LPedTop, analysis_type analysis)
{
    int nonuniq;
    int set_uniq = 0;

    if (LPedTop->UniqueIds && analysis != TO_PAP && analysis != IQLS) {
        nonuniq = check_unique_ids(LPedTop); 
/*     if (nonuniq) { */
/*       set_uniq=1; */
/*     } */
    } else {
        /* Mega2 created unique IDs */
        nonuniq = 0;
    }

    if (analysis == TO_PAP || analysis == IQLS) {
        mssgf("PAP/IQLS unique \"ID\" values generated by Mega2.\n");
        set_uniq = 1;
    } else  if (nonuniq) {
        if (analysis == CRANEFOOT) {
            mssgf("New unique \"ID\" values generated by Mega2.\n");
            set_uniq = 1;
        } else {
            set_uniq = 0;
        }
    }

    if (AnalyInputMode == BATCH_FILE_INPUTMODE && BatchValueRead("Default_Set_Uniq")) {
        char ch = (char)set_uniq;
        BatchValueIfSet(ch, "Default_Set_Uniq");
        set_uniq = ch;
    } else if (AnalyInputMode == INTERACTIVE_INPUTMODE && nonuniq &&
               analysis != TO_PAP && analysis != CRANEFOOT && analysis != IQLS) {
        char select_[10];
        int  select;
        int menu_item = 0;
        draw_line();
        printf("0) Done with this menu - please proceed\n"); menu_item++;
        printf(" %d) Generate new unique IDs? [%s]\n", menu_item++, yorn[set_uniq]);

        printf("Select from options 0-%d > ", menu_item-1);
        fcmap(stdin, "%s", select_); newline;
        sscanf(select_, "%d", &select);

        switch(select) {
        case 0:
            break;
        case 1:
            if (select == 1)
                set_uniq = (set_uniq ? 0 : 1);
            BatchValueSet(yorn[set_uniq][0], "Default_Set_Uniq");
            batchf("Default_Set_Uniq");
            break;
        default:
            printf("Unknown option %s\n", select_);
            break;
        }
    }
        
    return set_uniq;
}

void allelecnt_check(linkage_ped_top *Top, analysis_type analysis)
{
    Tod tod_ac("AlleleCnt check");
    linkage_locus_top *LTop = Top->LocusTop;
    int num_nbi = num_reordered;
    int first_time=1;
    int locus = 0;

    SECTION_ERR_INIT(not_bialleleic);

    strcpy(err_msg, "");
//  printf("%d\n", locus_num);
    for (int i = 0; i < num_reordered; i++) {
        locus = reordered_marker_loci[i];
//      printf("%d %d: %d %s\n",
//             i, locus, LTop->Locus[locus].AlleleCnt, LTop->Locus[locus].LocusName);
        if (LTop->Locus[locus].Type == AFFECTION || LTop->Locus[locus].Type == QUANT) {
            num_nbi--;
        } else if (LTop->Locus[locus].AlleleCnt > 2) {
            num_nbi--;

            if (strlen(err_msg) >= 67) {
                SECTION_ERR(not_bialleleic);
                warnf(err_msg);
                strcpy(err_msg, LTop->Locus[locus].LocusName);
            } else if (first_time) {

                SECTION_ERR(not_bialleleic);
                warnf("Excluding these markers because they have more than 2 alleles.");
                SECTION_ERR_HEADER(not_bialleleic);

                strcpy(err_msg, LTop->Locus[locus].LocusName);
                first_time = 0;
            } else {
                grow(err_msg, ", %s", LTop->Locus[locus].LocusName);
            }
        }
    }
    if (strlen(err_msg) > 0) {
        SECTION_ERR(not_bialleleic);
        warnf(err_msg);
    }
    SECTION_ERR_FINI(not_bialleleic);

    tod_ac();

    if (analysis == TO_PLINK && num_nbi == 0) {
        errorvf("%d (of %d) valid marker loci left for PLINK to analyze, terminating Mega2!\n",
                num_nbi, num_reordered);
        EXIT(OUTPUT_FORMAT_ERROR);
    }

    if (analysis == IQLS && num_nbi == 0) {
        errorvf("%d (of %d) valid marker loci left for IQLS to analyze, terminating Mega2!\n",
                num_nbi, num_reordered);
        EXIT(OUTPUT_FORMAT_ERROR);
    }
}

#ifdef DELAY_ZERO
genozero delay_zero;
void dozero(ped_top *Top, linkage_ped_top *LPedTop, analysis_type analysis)
{
   linkage_locus_top *LLTop = LPedTop->LocusTop;

    for (genozerop gp = delay_zero.begin(); gp != delay_zero.end(); gp++) {
        int locus = gp->locus;
        int ped   = gp->ped;
        int per   = gp->per;

        if (per != -1) {
                /*
                set_2alleles(Top->PedTree[ped].Entry[per].Marker, locus,
                             &LLTop->Locus[locus], 0, 0);
                */
                set_2alleles(LPedTop->PedBroken[ped].Entry[per].Marker, locus,
                             &LLTop->Locus[locus], 0, 0);
        } else {
            for (int entry=0; entry < LPedTop->PedBroken[ped].EntryCnt; entry ++) {

                    set_2alleles(LPedTop->PedBroken[ped].Entry[entry].Marker, locus,
                                 &LLTop->Locus[locus], 0, 0);
            }
        }
    }
}
#endif

ped_status      PedStat;
void full_check(ped_top *Top, linkage_ped_top *LPedTop, analysis_type analysis)
{
    int chr;
#ifndef DELAY_ZERO
    int             entry;
#endif
    int             ped, locus, select=-1, menu_item=0;
    int             stat, abortl=0, abortf=0, loc_err=0;
    /* below are the error flags */
    int             aexceed, imend, hmend, freq_mis = 0;
    int             chk_aexceed, chk_imend, chk_hmend;
    int             set_uniq=0, nonuniq;
    linkage_locus_top *LLTop = LPedTop->LocusTop;
    char            select_[10];
    char            toggle_str[50];

    /* These are set if the corresponding menu item is toggled */
    int             exit_select=-1, halftyped_select=-1;
    int             exceedall_select = -1, invalid_select=-1;
    int             uniq_select = -1;
//  int             plink_locus_num;
    int             num_mito_hetero = 0, num_mito_non_maternal = 0;

    zero = REC_UNKNOWN;

    clear_ped_status(&PedStat);

    /* changed the locus checking loop to skip numbered loci
       for QUANT summary. This can be modified to accommodate other
       cases where marker loci need not be changed */

#ifdef need_recode_alleles
    /*   if (num_traits > 0) { */
    /*     printf("Checking locus integrity...\n"); */
    /*     for (locus = 0; locus < num_traits; locus++) { */
    /*       if (global_trait_entries[locus] >= 0) { */
    /* 	stat=check_locus(&(LTop->Locus[global_trait_entries[locus]])); */
    /* 	abortl = imax(abortl, stat); */
    /*       }	 */
    /*     } */
    Tod tod_cl("check_locus");

    plink_locus_num = LLTop->MarkerCnt;
    if (plink_locus_num > 0) {
        for (locus = LLTop->PhenoCnt; locus < LLTop->LocusCnt; locus++) {
            if (LLTop->Locus[locus].Marker->chromosome == MISSING_CHROMO) continue;
            stat=check_locus(LLTop->Locus + locus, analysis, &plink_locus_num);
            abortl = imax(abortl, stat);
        }
        SECTION_ERR_FINI(check_locus);
        SECTION_ERR_FINI(biallele);
    }

    tod_cl();
#endif

    if (! database_dump) {
        if ((analysis == TO_PLINK || analysis == IQLS) /* && plink_locus_num < LTop->LocusCnt */) {
            allelecnt_check(LPedTop, analysis);
        }
        if (abortl > 0 && abortl <= 4) {
            loc_err=1;
        }
    }

    printf("Done checking locus integrity.\n");
    if (Display_Errors == 0) {
        printf("Check error logs for full list of locus data problems.\n");
        draw_line();
        Display_Errors=1;
    }

    if (LPedTop->UniqueIds && analysis != TO_PAP && analysis != IQLS) {
        nonuniq = check_unique_ids(LPedTop); 
    } else {
        /* Mega2 created unique IDs */
        nonuniq = 0;
    }

    if (! database_dump) {
        if (analysis == TO_PAP || analysis == IQLS) {
            set_uniq = 1;
        } else {
            if (nonuniq) {
//              printf("Warning: Duplicate values in the \"ID\" column.\n");
                if (analysis == CRANEFOOT) {
//                  printf("    Unique IDs will be automatically generated as required by %s\n",
//                         ProgName);
                    set_uniq = 1;
                } else {
                    set_uniq = 0;
                }
            }
        }
    }

    printf("Checking pedigree integrity...\n");
    /* separated out the checking of pedigree ids and relationships from
       checking genotypes, so that we can skip the latter, if there are no
       marker loci.
    */
    Tod tod_cpr("check_ped_relations");
    for (ped = 0; ped < Top->PedCnt; ped++)   {
        stat = check_ped_relations(&(Top->PedTree[ped]), &PedStat);
        abortf = imax(abortf, stat);
    }
    SECTION_ERR_FINI(check_ped_relations); // does not use MSSG but should
    tod_cpr();

    Tod tod_mito("mito_transmission_report");
    for (chr=0; chr < MaxChromo; chr++) {
        if (global_chromo_entries[chr] == MITO_CHROMOSOME) {
            mito_transmission_report(Top,
                                     LLTop->PedRecDataType == Postmakeped || 
                                       LLTop->PedRecDataType == Premakeped,
                                     &num_mito_hetero, &num_mito_non_maternal);
            if (num_mito_hetero > 0 || num_mito_non_maternal > 0)
                abortf = imax(abortf, 1);
            break;
        }
    }
    if (Display_Errors == 0) {
        printf("Check ERROR or LOG file for a complete list of errors.\n");
        Display_Errors=1;
    }
    printf("Done checking pedigree integrity.\n");
    tod_mito();

#ifdef need_recoded_alleles
    Tod tod_iofc("input_observed_freq_check");
    if (LLTop->MarkerCnt > 0 /* && analysis != TO_PLINK */) {
        freq_mis = input_observed_freq_check(LPedTop, FreqMismatchThreshold);
        if (freq_mis > 0)
            abortf = imax(abortf, 1);
    } else {
        freq_mis=0;
    }
    tod_iofc();
#endif
    Tod tod_cepi("check epilog code: do reset");
    /* set imend and hmend to 1 if errors are present */
    imend = chk_imend   = 1;
    hmend = chk_hmend   = (Input_Format == in_format_binary_PED) ? 0 : 1; // PLINK binary is NEVER halftyped
    aexceed = chk_aexceed = (Input_Format == in_format_linkage) ? 1 : 0;    // only necessary for Linkage
    if (loc_err == 1 || PedStat.entry_unconnected > 0 || freq_mis || nonuniq ||
        num_mito_hetero > 0 || num_mito_non_maternal > 0) {
        draw_line();  exclaim();
        warnvf("Thus far found these problems/errors in input data (see MEGA2.ERR for details): \n");
        if (loc_err == 1) {
            warnvf(" -> Invalid markers seen: unacceptible allele count or frequencies don't sum to 1.0\n");
        }

        if (PedStat.entry_unconnected > 0) {
            warnvf(" -> Unconnected entries\n");
        }

        if (freq_mis) {
            warnvf(" -> Input and observed allele frequencies do not match.\n");
        }

        if (nonuniq) {
            warnvf(" -> Duplicate values in the \"ID\" column.\n");
        }

        if (num_mito_hetero > 0) {
            warnvf(" ->  Found one or more heterozygous mitochondrial genotypes.\n");
        }

        if (num_mito_non_maternal > 0) {
            warnvf(" ->  Found one or more non-maternal mitochondrial transmissions.\n");
        }

        /*    depending on which menu-items will be displayed,
              set the three _select items which denote which flag to
              toggle or to exit */
        exclaim(); draw_line();
    }

    // Default_Reset_Invalid
    // This option defines how invalid genotypes should be handled without pausing for user-input
    // via the invalid-genotypes menu (which is skipped). If set to yes the genotypes will be reset
    // to unknowns, and if set to no invalid genotypes will not be reset.
    if (InputMode != INTERACTIVE_INPUTMODE && (
                   BatchValueRead("Default_Reset_Halftype") ||
                   BatchValueRead("Default_Reset_Mendelerr") ||
                   BatchValueRead("Default_Reset_Alleleerr") ||
                   BatchValueRead("Default_Set_Uniq")  )) {

        char ch;
        BatchValueGet(ch, "Default_Reset_Halftype");
        hmend = ch;
        BatchValueGet(ch, "Default_Reset_Mendelerr");
        imend = ch;
        BatchValueGet(ch, "Default_Reset_Alleleerr");
        aexceed = ch;
        BatchValueGet(ch, "Default_Set_Uniq");
        set_uniq = ch;

    } else if (InputMode != INTERACTIVE_INPUTMODE &&
        Mega2BatchItems[/* 26 */ Default_Reset_Invalid].items_read) {
        int mask = Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt;
        if (mask == 'n' || mask ==  'N')
            hmend = imend = aexceed = 0;

    } else {
        int imendf = 0, hmendf = 0, aexceedf = 0;
        int mask = 0;
        printf("Specify whether to reset Mendelianly inconsistent loci to missing: \n");
        while (select != 0) {
            strcpy(toggle_str, "");
            menu_item=0;
            draw_line();
            printf("0) Done with this menu - please proceed\n"); menu_item++;
            if (chk_hmend /* PedStat.halftyped */) {
/*
                if (analysis == TO_MENDEL || SIMWALK2(analysis)) {
                    printf(" NOTE: Proceeding will zero out all half-typed individuals.\n");
                } else 
*/
                {
                    printf(" %d) Set half-typed genotypes to unknown [%s].\n", menu_item,
                           yorn[hmendf]);
                    if (!hmendf) {
                        printf("    If \"no\" is indicated, half-typed genotypes\n");
                        printf("    will not be looked for.\n");
                    }
                    sprintf(toggle_str, "; options %d", menu_item);
                    halftyped_select = menu_item++;
                }
            }
            if (chk_imend /* PedStat.genotype_invalid */) {
                printf(" %d) Set all genotypes to unknown within entire pedigrees\n",
                       menu_item);
                printf("    at each Mendelianly-inconsistent locus? [%s]\n",
                       yorn[imendf]);
                if (!imendf) {
                    printf("    If \"no\" is indicated, Mendelianly-inconsistent loci\n");
                    printf("    will not be looked for.\n");
                }
                grow(toggle_str, ", %d", menu_item);
                invalid_select = menu_item++;
            }
            if (chk_aexceed /* PedStat.exceed_allcnt */) {
                printf(" %d) Set out-of-bound genotypes to unknown? [%s]\n",
                       menu_item, yorn[aexceedf]);
                if (!aexceedf) {
                    printf("    If \"no\" is indicated, out-of-bound genotypes\n");
                    printf("    will not be looked for.\n");
                }
                grow(toggle_str, ", %d", menu_item);
                exceedall_select = menu_item++;
            }

            // NOTE: There seems to be no batch file item to cover this case...

            if (nonuniq /* && analysis != TO_PAP && analysis != CRANEFOOT && analysis != IQLS */) {
                printf(" %d) Generate new unique IDs? [%s]\n", menu_item, yorn[set_uniq]);
                grow(toggle_str, ", %d", menu_item);
                uniq_select = menu_item++;
            }

            if (freq_mis) {
                printf("Don't know how to correct for frequency mismatches.\n");
            }
//          if (PedStat.genotype_invalid || PedStat.halftyped || PedStat.exceed_allcnt) {
                strcat(toggle_str, " toggle");
//          }
            printf(" %d) EXIT Mega2.\n", menu_item);
            exit_select=menu_item;

            printf("Select from options 0-%d%s > ", menu_item, toggle_str);
            fcmap(stdin, "%s", select_); newline;
            sscanf(select_, "%d", &select);
            switch(select) {
            case 0:
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
                if (select == halftyped_select) {
                    hmendf = (hmendf ? 0 : 1);
                    mask = 1;
                    break;
                } else if (select == invalid_select) {
                    imendf = (imendf ? 0 : 1);
                    mask = 1;
                    break;
                } else if (select == exceedall_select) {
                    aexceedf = (aexceedf ? 0 : 1);
                    mask = 1;
                    break;
                } else if (select == uniq_select) {
                    set_uniq = (set_uniq ? 0 : 1);
                    mask = 1;
                    break;
                } else if (select == exit_select) {
                    printf("Terminating Mega2.\n");
                    EXIT(EARLY_TERMINATION);
                }
            default:
                printf("Unknown option %s\n", select_);
                break;
            }
        }
        if (mask == 0) {
            BatchValueSet(yorn[mask][0], "Default_Reset_Invalid");
            batchf("Default_Reset_Invalid");
        } else {
            BatchValueSet(yorn[hmendf][0], "Default_Reset_Halftype");
            batchf("Default_Reset_Halftype");

            BatchValueSet(yorn[imendf][0], "Default_Reset_Mendelerr");
            batchf("Default_Reset_Mendelerr");

            BatchValueSet(yorn[aexceedf][0], "Default_Reset_Alleleerr");
            batchf("Default_Reset_Alleleerr");
                                    
            BatchValueSet(yorn[set_uniq][0], "Default_Set_Uniq");
            batchf("Default_Set_Uniq");
        }
        hmend = hmendf;
        imend = imendf;
        aexceed = aexceedf;
    }

    tod_cepi();
    if (abortf == 2 || abortl == 3) {
        exclaim();
        if (abortf == 2)   {
            errorvf("Fatal errors in pedigree file %s.\n", mega2_input_files[0]);
        } else if (abortl == 3)   {
            errorvf("Fatal errors in locus file %s.\n", mega2_input_files[1]);
        }
        exclaim();
        EXIT(DATA_INCONSISTENCY);
    }

    // There are three cases in which genotypes are set to unknown...
    if (Mega2BatchItems[/* 26 */ Default_Reset_Invalid].items_read == 0 &&
        (hmend != 0 || imend != 0 || aexceed != 0)) {
        
        // It is possible to have reached the above menu from batch or interactive mode.
        // Here we note it in the batch file.
        // NOTE: The manual says "Yes" and "No" are the values, but 'copt' is a char.
//      Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt = 'Y';
//      batchf(Default_Reset_Invalid);
    }

    /* else abort = 1*/
    Display_Errors=1;
    bool first = true;

    draw_line();
    Tod tod_hmend_all("reset all half typed");
    Tod tod_hmend(20);
    /* NB
       chk_hmend is true iff there is reason to believe the checking is necessary
       hmend is true iff the user said yes please reset the bad data
       using chk_hmend in the conditional lets the code look for the bad data and report it
          if hmend is also set it is reset and written to .RESET
       using hmend in the conditional does everything iff the user requested a reset
    */
    if (/*chk_*/hmend) {
        FILE *reset_fp = NULL;
        mssgf("Setting any half-typed genotypes to unknowns for the indicated pedigree/person/locus combinations.");
        PedStat.halftyped = 0;
        HalfTypedReset=0;

        for (ped = 0; ped < Top->PedCnt; ped++) {
            tod_hmend.reset();
            for (locus = LLTop->PhenoCnt; locus < LLTop->LocusCnt; locus++) {
                if (LLTop->Locus[locus].Marker->chromosome == MISSING_CHROMO) continue;
                stat = check_half_type(&(Top->PedTree[ped]), &PedStat, LLTop,
                                       ped, locus, LPedTop->UniqueIds, ped,
                                       &reset_fp, &first, hmend);
                abortf = imax(abortf, stat);
                if (stat) {
                    HalfTypedReset=1;
                }
            }
            tod_hmend("reset half typed");
        }
        SECTION_ERR_FINI(check_half_type);
        /* set invalid genos to 0 0 */
        if (reset_fp != NULL)
            fclose(reset_fp);

        msgvf("Status: Reset %d half-typed markers.\n", PedStat.halftyped);
        if (SECTION_ERR_COUNT(check_half_type) > 0)
            msgvf(" Check error logs for full list of half typed genotypes.\n");
        else
            msgvf("\n");

        tod_hmend_all();

    } else {
        HalfTypedReset=0;
    }

    Tod tod_xmend_all("reset ALL out-of-bound genotypes");
    Tod tod_xmend(20);
//  int OOBReset=0;
    if (/*chk_*/aexceed) {
        FILE *reset_fp = NULL;
        mssgf("Setting any genotypes with out-of-bounds alleles to unknown ...");
        PedStat.exceed_allcnt = 0;

        for (ped = 0; ped < Top->PedCnt; ped++) {
            tod_xmend.reset();
            for (locus = LLTop->PhenoCnt; locus < LLTop->LocusCnt; locus++) {
                if (LLTop->Locus[locus].Marker->chromosome == MISSING_CHROMO) continue;
                stat = check_out_of_bounds(&(Top->PedTree[ped]), &PedStat, LLTop,
                                           ped, locus, LPedTop->UniqueIds, ped,
                                           &reset_fp, &first, aexceed);
                abortf = imax(abortf, stat);
                if (stat) {
//                  OOBReset=1;
                }
            }
            tod_xmend("reset out-of-bound genotypes");
        }
        /* set invalid genos to 0 0 */
        SECTION_ERR_FINI(check_oob);
        if (reset_fp != NULL)
            fclose(reset_fp);

        msgvf("Status: Reset %d out-of-bound allele values.\n", PedStat.exceed_allcnt);
        if (SECTION_ERR_COUNT(check_oob) > 0)
            msgvf(" Check error logs for full list of out-of-bound genotypes.\n");
        else 
            msgvf("\n");

        tod_xmend_all();
    } else {
//      OOBReset=0;
    }

    Tod tod_imend_all("reset ALL non mendelian");
    Tod tod_imend(20);
    if (/*chk_*/imend) {
#ifndef DELAY_ZERO
        int recode = LLTop->PedRecDataType == Postmakeped || LLTop->PedRecDataType == Premakeped;
#endif
        FILE *reset_fp = NULL;
        mssgf("Setting any inconsistent genotypes to unknowns for the indicated pedigree/locus combinations.");
        int rm;
        NonMendelianReset=0;
        PedStat.genotype_invalid = 0;
        for (ped = 0; ped < Top->PedCnt; ped++) {
            ped_rec **Sibs = CALLOC((size_t) Top->PedTree[ped].EntryCnt, ped_rec *);
            tod_imend.reset();
            for (locus = LLTop->PhenoCnt; locus < LLTop->LocusCnt; locus++) {
                if (LLTop->Locus[locus].Marker->chromosome > lastautosome) continue;
                if (LLTop->Marker[locus].chromosome == MITO_CHROMOSOME) {
                    continue;           /* Continue without inheritance checks */
                }
                rm = check_invalid_fam(&(Top->PedTree[ped]), &PedStat, LLTop,
                                       ped, locus, LPedTop->UniqueIds,
                                       Sibs);
                if (rm) {
                    if (reset_fp == NULL) {
                        if (first) {
                            reset_fp = fopen(Mega2ResetRun, "w");
                            summary_time_stamp(mega2_input_files, reset_fp, "");
                            first = false;
                        } else
                            reset_fp = fopen(Mega2ResetRun, "a");
                        fprintf(reset_fp, "Mendelianly inconsistent pedigrees %s:\n",
                                imend ? "RESET" : "ALLOWED");
                        fprintf(reset_fp, "Pedigree   Person    Marker\n");
                    }
                    fprintf(reset_fp, "%s      All   %s\n", 
                            Top->PedTree[ped].Name, 
                            LLTop->Locus[locus].LocusName);

                    NonMendelianReset = 1;
                    if (imend) {
#ifdef DELAY_ZERO
                        delay_zero.push_back(geno(locus, ped, -1));
#else
                        for (entry=0; entry < Top->PedTree[ped].EntryCnt; entry ++) {
                            /* set invalid genos to 0 0 */
                            if (recode)
                                set_2alleles(Top->PedTree[ped].Entry[entry].Marker, locus,
                                             &LLTop->Locus[locus], 0, 0);
                            else
                                set_2Ralleles(Top->PedTree[ped].Entry[entry].Marker, locus,
                                              &LLTop->Locus[locus], zero, zero);
                        }
#endif
                    }
                }
            }
            tod_imend("reset non mendelian");
            free(Sibs);
        }
        SECTION_ERR_FINI(check_inheritance);
        SECTION_ERR_FINI(check_sibship_alleles);

        if (reset_fp != NULL) {
            fclose(reset_fp);
        }

        msgvf("Status: Reset %d Mendelianly inconsistent markers values for pedigrees.\n",
              PedStat.genotype_invalid);
        if (SECTION_ERR_COUNT(check_inheritance) > 0 || 
            SECTION_ERR_COUNT(check_sibship_alleles) > 0) {
            msgvf(" Check error logs for full list of pedigrees with reset genotypes.\n");
        }

        tod_imend_all();
    } else {
        NonMendelianReset=0;
    }

    /* Store the resets instead of freeing */
    Tod tod_wrs("write_reset_summary");
/*
    write_reset_summary(Top, invalid_genos, halftyped, exceed_allcnt,
                        PedStat.genotype_invalid, PedStat.halftyped,
                        PedStat.exceed_allcnt, LPedTop->UniqueIds);
*/
    tod_wrs();
    if (Display_Errors == 0) {
        printf("Check error logs for full list of reset genotypes.\n");
        draw_line();
        Display_Errors=1;
    }

    /* Create unique IDs */
    Tod tod_uniq("create_unique_id");
    if (set_uniq)
        create_unique_ids(LPedTop, analysis);
    tod_uniq();

    if (PedStat.genotype_invalid > 0 || PedStat.halftyped > 0 || PedStat.exceed_allcnt > 0) {
        log_line(mssgf);
        if (InputMode != INTERACTIVE_INPUTMODE && Mega2BatchItems[/* 26 */ Default_Reset_Invalid].items_read) {
            mssgf(
                "DEFAULT HANDLING MODE FOR INVALID GENOTYPES SPECIFIED IN BATCH FILE:");
            int copt = Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt;
            if ((copt == 'y' || copt == 'Y')) 
                msgvf(" Reset allele to 0/0\n\n");
            else
                msgvf(" Leave allele alone\n\n");
        }
    }
/*
    if (PedStat.genotype_invalid) {
        if (imend) {
            mssgf("Mendelianly-inconsistent genotypes set to unknown as requested.");
        } else {
            mssgf("User chose NOT to set Mendelianly-inconsistent genotypes to unknowns.");
        }
    }
    if (PedStat.halftyped) {
        if (hmend) {
            mssgf("Half-typed individuals' genotypes set to unknown as requested.");
        } else {
            mssgf("User chose NOT to set half-typed genotypes to unknowns.");
        }
    }
    if (PedStat.exceed_allcnt) {
        if (aexceed) {
            mssgf("Genotypes with out-of-bounds alleles set to unknowns as requested.");
        } else {
            mssgf("User chose NOT to set genotypes with out-of-bounds alleles to unknowns.");
        }
    }
    if (nonuniq  && ! database_dump) {
        if (set_uniq) {
            mssgf("New unique \"ID\" values generated by Mega2.\n");
        } else {
            mssgf("User chose not to generate new \"ID\" values.\n");
        }
    }
*/
    if (PedStat.genotype_invalid > 0 || PedStat.halftyped > 0 ||
        PedStat.exceed_allcnt > 0 || nonuniq) {
        log_line(mssgf);
        warnvf("Mega2 resolved these problems/errors in input data (see MEGA2.ERR for details): \n");
        if (PedStat.halftyped > 0) {
            warnvf(" -> %d Half-typed genotypes\n", PedStat.halftyped);
            /* if (analysis == TO_MENDEL || SIMWALK2(analysis)) */ {
                warnvf("        Half-typed genotypes will be set to unknown \n");
//              warnvf("        as required by %s\n", ProgName);
            }
        }
        if (PedStat.genotype_invalid > 0) {
            warnvf(" -> %d Mendelian-inconsistent genotypes\n", PedStat.genotype_invalid);
            warnvf("        All family members for each inconsistent Family/Marker will be set to unknown \n");
        }

        if (PedStat.exceed_allcnt > 0) {
            warnvf(" -> %d Genotype allele numbers exceed locus allele count\n",
                   PedStat.exceed_allcnt);
        }
        printf("Check %s for details.\n", Mega2Err);
    }

    if (HasMarkers && PedStat.checked_sibship > 0) {
        /* Inheritance checks may have been inadequate */
        draw_line();
        warnf("Inheritance checks may have been inadequate, ");
        warnf("Reason: One or more sibships with untyped parents.");
        warnf("Please verify with a program for checking Mendelian inconsistencies");
        warnf("such as Pedcheck.");
        log_line(mssgf);

    }
    Display_Messages=1;
    Display_Errors=1;
}

void recode_check(linkage_ped_top *LPedTop, analysis_type analysis)
{
    linkage_locus_top *LLTop = LPedTop->LocusTop;
    int stat, abortl = 0;
    Tod tod_cl("check_locus");

    // check_locus() needs recode_locus_top() to have been run.
    // it is not looking at alleles.
    int plink_locus_num = LLTop->MarkerCnt;
    if (plink_locus_num > 0) {
        for (int locus = LLTop->PhenoCnt; locus < LLTop->LocusCnt; locus++) {
            if (LLTop->Locus[locus].Marker->chromosome == MISSING_CHROMO) continue;
            stat=check_locus(LLTop->Locus + locus, analysis, &plink_locus_num);
            abortl = imax(abortl, stat);
        }
        SECTION_ERR_FINI(check_locus);
        SECTION_ERR_FINI(biallele);
    }
    tod_cl();

    int abortf = 0, freq_mis = 0;

    // input_observed_freq_check also needs recode_locus_top() to have been run.
    Tod tod_iofc("input_observed_freq_check");
    if (LLTop->MarkerCnt > 0 /* && analysis != TO_PLINK */) {
        freq_mis = input_observed_freq_check(LPedTop, FreqMismatchThreshold);
        if (freq_mis > 0)
            abortf = imax(abortf, 1);
    } else {
        freq_mis=0;
    }
    tod_iofc();
}

void show_reset_input() {
    if (PedStat.genotype_invalid > -1)
        msgvf("\t%d mendelian-inconsistent markers.\n", PedStat.genotype_invalid);

    if (PedStat.halftyped > -1)
        msgvf("\t%d half-typed markers.\n", PedStat.halftyped);

    if (PedStat.exceed_allcnt > -1)
        msgvf("\t%d out-of-bound allele values.\n", PedStat.exceed_allcnt);

    msgvf("\n");
}


/*---------------------------------------------------------------+
  | Check if the pedigree should be renumbered (if any of the     |
  | members have a parent with a greater ID).                     |
  | Remember to call this function AFTER calls                    |
  | to affected_by_status and reassign_affecteds.                 |
  +---------------------------------------------------------------*/
void            check_renumber_ped(ped_tree *Ped)

{
    int             entry;
    int             ok;

    ped_rec *Entry;

    ok = 1;
    for (entry = 0; (ok) && (entry < Ped->EntryCnt); entry++)
        {
            Entry = &(Ped->Entry[entry]);
            if ((Entry->Mother != NULL) && (Entry->ID < Entry->Mother->ID))
                ok = 0;
            if ((Entry->Father != NULL) && (Entry->ID < Entry->Father->ID))
                ok = 0;
        }

    if (ok) return;

    /*  printf("Renumbering pedigree %s.\n", Ped->Name); */
    if (renumber_ped(Ped) < 0)
        {
            errorvf("Fatal error renumbering the pedigree.\n");
            EXIT(OUTPUT_FORMAT_ERROR);
        }
}

void     create_unique_ids(linkage_ped_top *Top, analysis_type analysis)
{
    int ped, per;
    linkage_ped_rec *entry;

    if (analysis == TO_PAP || analysis == IQLS) {
        set_pap_unique_id(Top);
    } else {
        for (ped=0; ped < Top->PedCnt; ped++) {
            for (per=0; per < Top->Ped[ped].EntryCnt; per++) {
                entry = &(Top->Ped[ped].Entry[per]);
                sprintf(entry->UniqueID, "%d_%d", Top->Ped[ped].Num, entry->ID);
            }
        }
    }
    return;
}
