/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

===========================================================================
*/

//NrvB: ped_top->Locus ONLY contains markers and is a copy of the Top->LocusTop->Locus data !!

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"
#include "tod.hh"

#include "input_check.h"

#include "check_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
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


/* prototype definitions */
void            full_check(ped_top *Top, linkage_ped_top *LPedTop,
			   analysis_type analysis);
int             UnMappedLociCheck(linkage_ped_top *Top);
int             input_observed_freq_check(linkage_ped_top *Top, double threshold);
void            check_renumber_ped(ped_tree *Ped);

static  void    write_reset_summary(ped_top *PTop, loc_list *invalidp,
				    person_loc_list *halftypedp,
				    person_loc_list *outofboundsp,
				    int imend, int ht, int aexceed,
				    int uniqueids);

static void     create_unique_ids(linkage_ped_top *Top, analysis_type analysis);

/* end of prototypes */

/* ==================================================================
   ,  check for  existence for all map positions*/
/* will assume 1st locus is 'disease loci' and skip it in the following check...
   POOR ASSUMPTION: NEED TO FIX THIS */

int             UnMappedLociCheck(linkage_ped_top * Top)
{
    int             i, unmapped = 0;


    for (i = 0; i < num_reordered; i++) {
        if (Top->LocusTop->Marker[reordered_marker_loci[i]].pos_avg < 0.000 &&
            (Top->LocusTop->Locus[reordered_marker_loci[i]].Type == NUMBERED ||
             Top->LocusTop->Locus[reordered_marker_loci[i]].Type == BINARY)) {
            printf("%s ",
                   Top->LocusTop->Locus[reordered_marker_loci[i]].Name);
            unmapped = 1;
        }
    }
    return unmapped;
}

/* Compute the sum of squared differences between input and output alleles
   for all numbered loci, warn if this is above the thershold value
*/
int  input_observed_freq_check(linkage_ped_top *Top, double threshold)

{

    int j, i, all, num_alleles, found_problem=0;
    double *sum_squared, diff;
    allele_freq_struct ar;
    int displayed_errors, first_time=1, flush=0;

    sum_squared = CALLOC((size_t) num_reordered, double);

    Display_Errors=1;
    Display_Messages=1;
    displayed_errors = 0;

    for (j = 0; j < num_reordered; j++) {
        i = reordered_marker_loci[j];
        sum_squared[j] = 0.0;
        if (Top->LocusTop->Locus[i].Type == NUMBERED ||
            Top->LocusTop->Locus[i].Type == BINARY) {

            ar.TotalAlleles = 0;
            ar.TotalPeople = 0;
            ar.HalfTyped = 0;
            ar.TotalKnownAlleles = 0;

            num_alleles=Top->LocusTop->Locus[i].AlleleCnt;
            ar.Allele_Bin = CALLOC((size_t) num_alleles, int);
            ar.allele_freq = CALLOC((size_t) num_alleles, double);
            for (all = 0; all < num_alleles; all++) {
                ar.Allele_Bin[all]=0;
            }
            count_alleles(Top, &ar, i, 1, 0);
            sum_squared[j]=0.0;

            for (all = 0; all < num_alleles; all++) {
                ar.allele_freq[all] =
                    ((double) (ar.Allele_Bin[all])) / ((double) (ar.TotalKnownAlleles));
                diff = ar.allele_freq[all] - Top->LocusTop->Locus[i].Allele[all].Frequency;
                sum_squared[j] += diff * diff;
            }

            if (sum_squared[j] > threshold) {
                found_problem=1;

                if (first_time) {
                    /* First error */
                    log_line(mssgf);
                    warnf("For these markers, sum of squared differences between ");
                    sprintf(err_msg,
                            "input and observed allele frequencies exceeded threshold %9.7f.",
                            threshold);
                    warnf(err_msg);
                    strcpy(err_msg, "  ");
                    first_time = 0;
                }
                grow(err_msg, "%s: SSD=%9.7f ",
                     Top->LocusTop->Locus[i].Name, sum_squared[j]);

                flush++;
                if (flush == 2) {
                    SUPPRESS_MSSG(displayed_errors)
                        if (displayed_errors > MAX_PED_ERRORS) {
                            Display_Messages = 0;
                        }
                    mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
                    displayed_errors++;
                }
            }
            free(ar.Allele_Bin); ar.Allele_Bin= NULL;
            free(ar.allele_freq); ar.allele_freq=NULL;
        }
    }

    if (flush) {
        mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
    }

    if (Display_Messages == 0) {
        printf("Check error logs for full list of reset genotypes.\n");
        draw_line();
        Display_Messages=1;
    }

    free(sum_squared);

    return found_problem;
}

static void write_and_free1(FILE *fp, loc_list *invalidp, ped_top *PedTop,
			   locus_top *LTop)

{
    loc_list *next;
    for (; invalidp != NULL; invalidp = next) {
       next = invalidp->next;
       fprintf(fp, "%s      All   %s\n", 
               PedTop->PedTree[invalidp->ped].Name, 
               LTop->Locus[invalidp->locus].Name);
       free(invalidp);
    }
}

static void write_and_free2(FILE *fp, person_loc_list *invalidp, 
			    ped_top *PedTop,  locus_top *LTop, int uniqueids)

{
    person_loc_list *next;

    for (; invalidp != NULL; invalidp = next) {
        next = invalidp->next;
        fprintf(fp, "%s   %s   %s\n", PedTop->PedTree[invalidp->ped].Name,
                ((uniqueids == 1)? 
                 PedTop->PedTree[invalidp->ped].Entry[invalidp->person].LEntry->UniqueID :
                 PedTop->PedTree[invalidp->ped].Entry[invalidp->person].LEntry->OrigID),
                LTop->Locus[invalidp->locus].Name);
        free(invalidp);
    }
}


static  void  write_reset_summary(ped_top *PedTop,
				  loc_list *invalidp,
				  person_loc_list *halftypedp,
				  person_loc_list *outofboundsp,
				  int imend, int ht, int aexceed,
				  int uniqueids)

{

    FILE *reset_fp;

    if (!imend && !ht && !aexceed) {
        return;
    }

    reset_fp = fopen(Mega2ResetRun, "w");
    summary_time_stamp(mega2_input_files, reset_fp, "");

    if (imend) {
        Display_Messages=1;

        fprintf(reset_fp, "Mendelianly inconsistent pedigrees:\n");
        fprintf(reset_fp, "Pedigree   Person    Marker\n");

        write_and_free1(reset_fp, invalidp, PedTop, PedTop->LocusTop);
    }

    if (ht) {
        Display_Messages=1;
        fprintf(reset_fp, "Half-typed genotypes:\n");
        fprintf(reset_fp, "Pedigree   Person    Marker\n");
        write_and_free2(reset_fp, halftypedp, PedTop, PedTop->LocusTop, uniqueids);
    }

    if (aexceed) {
        Display_Messages=1;
        fprintf(reset_fp, "Genotypes with out-of-bounds alleles.\n");
        fprintf(reset_fp, "Pedigree   Person    Marker\n");
        write_and_free2(reset_fp, outofboundsp, PedTop, PedTop->LocusTop, uniqueids);
    }

    fclose(reset_fp);
    return;
}

/*---------------------------------------------------------------+
  | Given the pedtree Top, check all the pedigree and locus data. |
  | Print a summary if an error was detected.                     |
  |                                                               |
  | infl_type is used to determine if sex errors should be        |
  | ignored (since APM files do not contain sex data).            |
  +---------------------------------------------------------------*/

void      full_check(ped_top *Top, linkage_ped_top *LPedTop,
		     analysis_type analysis)

{
    int chr;
    int             entry, lloc, ped, locus, select=-1, menu_item=0;
    ped_status      PedStat;
    int             stat, abortl=0, abortf=0, loc_err=0;
    /* below are the error flags */
    int             aexceed, imend, hmend, freq_mis;
    int             set_uniq=0, nonuniq;

    register locus_top *LTop = Top->LocusTop;
    loc_list        *invalid_genos, *invalidp;

    person_loc_list *halftyped, *halftypedp, *exceed_allcnt, *outofboundsp;
    char            select_[10];
    char            toggle_str[20];

    /* These are set if the corresponding menu item is toggled */
    int             exit_select=-1, halftyped_select=-1;
    int             exceedall_select = -1, invalid_select=-1;
    int             uniq_select = -1;

    int             displayed_errors=0;
    int             plink_locus_num;

    int             num_mito_hetero = 0, num_mito_non_maternal = 0;

    clear_ped_status(&PedStat);
    HalfTypedReset=-1;
    NonMendelianReset=-1;

    /* changed the locus checking loop to skip numbered loci
       for QUANT summary. This can be modified to accommodate other
       cases where marker loci need not be changed */

    /*   if (num_traits > 0) { */
    /*     printf("Checking locus integrity...\n"); */
    /*     for (locus = 0; locus < num_traits; locus++) { */
    /*       if (global_trait_entries[locus] >= 0) { */
    /* 	stat=check_locus(&(LTop->Locus[global_trait_entries[locus]])); */
    /* 	abortl = imax(abortl, stat); */
    /*       }	 */
    /*     } */
    Tod tod_cl("check_locus");
    Display_Errors=1;
    plink_locus_num = LTop->LocusCnt;
    if (num_reordered > 0) {
        for (locus = 0; locus < LTop->LocusCnt; locus++) {
// LTop is only markers
            stat=check_locus(&(LTop->Locus[locus]),
                             &displayed_errors, analysis,
                             &plink_locus_num);
            abortl = imax(abortl, stat);
        }
    }
    Display_Errors=1;
    displayed_errors = 0;
    tod_cl();

    Tod tod_ac("AlleleCnt check");
    if ((analysis == TO_PLINK || analysis == IQLS) && plink_locus_num < LTop->LocusCnt) {
        int first_time=1;
        warnf("Excluding these markers because they have more than 2 alleles.");
        strcpy(err_msg, "");
        for (locus = 0; locus < LTop->LocusCnt; locus++) {
            if (LTop->Locus[locus].AlleleCnt > 2) {
                SUPPRESS_MSSG(displayed_errors)
                    if (displayed_errors > MAX_PED_ERRORS) {
                        Display_Errors = 0;
                    }

                if (strlen(err_msg) >= 67 || first_time) {
                    warnf(err_msg);
                    displayed_errors++;
                    strcpy(err_msg, LTop->Locus[locus].Name);
                    first_time = 0;
                } else {
                    grow(err_msg, ", %s", LTop->Locus[locus].Name);
                }
            }
        }
        Display_Errors = 1;
	displayed_errors = 0;
        if (strlen(err_msg) > 0) {
            warnf(err_msg);
        }
    }
    tod_cl();

    if (abortl <= 4) {
        loc_err=1;
    }

    if (analysis == TO_PLINK && plink_locus_num <= 0) {
        errorvf("No valid marker loci for PLINK to analyze, terminating Mega2!\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }

    if (analysis == IQLS && plink_locus_num <= 0) {
        errorvf("No valid marker loci for IQLS to analyze, terminating Mega2!\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }


    printf("Done checking locus integrity.\n");
    if (Display_Errors == 0) {
        printf("Check error logs for full list of locus data problems.\n");
        draw_line();
        Display_Errors=1;
    }

    displayed_errors = 0;
    printf("Checking pedigree integrity...\n");
    /* separated out the checking of pedigree ids and relationships from
       checking genotypes, so that we can skip the latter, if there are no
       marker loci.
    */
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
        set_uniq = 1;
    } else {
        if (nonuniq) {
            printf("Warning: Duplicate values in the \"ID\" column.\n");
            if (analysis == CRANEFOOT) {
                printf(
                    "    Unique IDs will be automatically generated as required by %s\n",
                    ProgName);
                set_uniq = 1;
            } else {
                set_uniq = 0;
            }
        }
    }

    Tod tod_cpr("check_ped_relations");
    for (ped = 0; ped < Top->PedCnt; ped++)   {
        stat = check_ped_relations(&(Top->PedTree[ped]), &PedStat);
        abortf = imax(abortf, stat);
    }
    tod_cpr();

    Tod tod_cpd("check_ped_data");
    if (num_reordered > 0 && HasMarkers) {
        /* init the check_sibship flag before checking the pedigrees*/
        PedStat.checked_sibship=0;
        invalid_genos = NULL;
        halftyped = NULL;
        exceed_allcnt = NULL;

        invalidp = invalid_genos;
        halftypedp = halftyped;
        outofboundsp = exceed_allcnt;

        Display_Errors=1;
        displayed_errors = 0;
        printf("Checking for half-typed, Mendelian inconsistencies and invalid alleles ...\n");
        for (ped = 0; ped < Top->PedCnt; ped++)   {
            stat = check_ped_data(&(Top->PedTree[ped]), &PedStat, Top->LocusTop,
                                  ped, &displayed_errors, LPedTop->UniqueIds);
            abortf = imax(abortf, stat);
            if (PedStat.invalid_genotypes != NULL) {
                if (invalidp == NULL) {
                    invalid_genos = PedStat.invalid_genotypes;
                    invalidp = invalid_genos;
                } else {
                    while(invalidp->next != NULL) {
                        invalidp = invalidp->next;
                    }
                    invalidp->next = PedStat.invalid_genotypes;
                }
                PedStat.invalid_genotypes = NULL;
            }
            if (PedStat.half_types != NULL) {
                if (halftypedp == NULL) {
                    halftyped = PedStat.half_types;
                    halftypedp=halftyped;
                } else {
                    while(halftypedp->next != NULL) {
                        halftypedp = halftypedp->next;
                    }
                    halftypedp->next = PedStat.half_types;
                }
                PedStat.half_types = NULL;
            }
            if (PedStat.allele_outof_bounds != NULL) {
                if (outofboundsp == NULL) {
                    exceed_allcnt = PedStat.allele_outof_bounds;
                    outofboundsp=exceed_allcnt;
                } else {
                    while(outofboundsp->next != NULL) {
                        outofboundsp = outofboundsp->next;
                    }
                    outofboundsp->next = PedStat.allele_outof_bounds;
                }
                PedStat.allele_outof_bounds = NULL;
            }
        }
        Display_Errors=1;
        displayed_errors = 0;
    }
    tod_cpd();

    Tod tod_mito("mito_transmission_report");
    for (chr=0; chr < main_chromocnt; chr++) {
        if (global_chromo_entries[chr] == MITO_CHROMOSOME) {
            mito_transmission_report(Top, &num_mito_hetero, &num_mito_non_maternal);
            break;
        }
    }
    if (Display_Errors == 0) {
        printf("Check ERROR or LOG file for a complete list of errors.\n");
        Display_Errors=1;
    }
    printf("Done checking pedigree integrity.\n");
    displayed_errors = 0;
    tod_mito();

    Tod tod_iofc("input_observed_freq_check");
    if (num_reordered > 0 && analysis != TO_PLINK) {
        freq_mis = input_observed_freq_check(LPedTop, FreqMismatchThreshold);
    } else {
        freq_mis=0;
    }
    tod_iofc();

    Tod tod_cepi("check epilog code: do reset");
    /* set imend and hmend to 1 if errors are present */
    imend=((PedStat.genotype_invalid > 0)? 1 : 0);
    hmend=((PedStat.halftyped > 0)? 1 : 0);
    aexceed = ((PedStat.exceed_allcnt > 0)? 1 : 0);

    if ((abortl == 1 || abortf) && (abortl != 2) && (abortf != 2)) {
        draw_line();  exclaim();
        warnvf("Found these problems/errors in input data (see MEGA2.ERR for details): \n");
        if (loc_err == 1) {
            warnvf(" -> Invalid marker allele number or frequencies \n");
        }

        if (PedStat.entry_unconnected > 0) {
            warnvf(" -> Unconnected entries\n");
        }
        if (PedStat.halftyped) {
            warnvf(" -> Half-typed genotypes\n");
            if (analysis == TO_MENDEL || SIMWALK2(analysis)) {
                warnvf("       Half-typed genotypes will be set to unknown \n");
                warnvf("       as required by %s\n", ProgName);
            }
        }
        if (PedStat.genotype_invalid > 0) {
            warnvf(" -> Mendelian-inconsistent genotypes\n");
        }

        if (PedStat.exceed_allcnt > 0) {
            warnvf(" -> Genotype allele numbers exceed locus allele count\n");
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

	// Default_Reset_Invalid
	// This option defines how invalid genotypes should be handled without pausing for user-input
	// via the invalid-genotypes menu (which is skipped). If set to yes the genotypes will be reset
	// to unknowns, and if set to no invalid genotypes will not be reset.
        if (InputMode != INTERACTIVE_INPUTMODE &&
            Mega2BatchItems[/* 26 */ Default_Reset_Invalid].item_read == 1
            ) {
            if (Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt == 'n' ||
                Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt == 'N')
                hmend = imend = aexceed = 0;
        } else {
            printf("How shall we proceed?\n");
            while (select != 0) {
                strcpy(toggle_str, "");
                menu_item=0;
                draw_line();
                printf("0) Done with this menu - please proceed\n"); menu_item++;
                if (PedStat.halftyped) {
                    if (analysis == TO_MENDEL || SIMWALK2(analysis)) {
                        printf(" NOTE: Proceeding will zero out all half-typed individuals.\n");
                    } else {
                        printf(" %d) Set half-typed genotypes to unknown [%s].\n", menu_item,
                               yorn[hmend]);
                        sprintf(toggle_str, ", %d", menu_item);
                        halftyped_select = menu_item++;
                    }
                }
                if (PedStat.genotype_invalid) {
                    printf(" %d) Set all genotypes to unknown within entire pedigrees\n",
                           menu_item);
                    printf("    at each Mendelianly-inconsistent locus? [%s]\n",
                           yorn[imend]);
                    if (toggle_str != NULL) {
                        grow(toggle_str, ", (%d", menu_item);
                    } else {
                        sprintf(toggle_str, ", %d", menu_item);
                    }
                    invalid_select = menu_item++;
                }
                if (PedStat.exceed_allcnt) {
                    printf(" %d) Set out-of-bound genotypes to unknown? [%s]\n",
                           menu_item, yorn[aexceed]);
                    sprintf(toggle_str, ", %d", menu_item);
                    exceedall_select = menu_item++;
                }

		// NOTE: There seems to be no batch file item to cover this case...
                if (nonuniq && analysis != TO_PAP && analysis != CRANEFOOT && analysis != IQLS) {
                    printf(" %d) Generate new unique IDs? [%s]\n", menu_item, yorn[set_uniq]);
                    sprintf(toggle_str, ", %d", menu_item);
                    uniq_select = menu_item++;
                }

                if (freq_mis) {
                    printf("Don't know how to correct for frequency mismatches.\n");
                }
                if (PedStat.genotype_invalid || PedStat.halftyped || PedStat.exceed_allcnt) {
                    strcat(toggle_str, " to toggle)");
                }
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
                        hmend = (hmend ? 0 : 1);
                        break;
                    } else if (select == invalid_select) {
                        imend = (imend ? 0 : 1);
                        break;
                    } else if (select == exceedall_select) {
                        aexceed = (aexceed ? 0 : 1);
                        break;
                    } else if (select == uniq_select) {
                        set_uniq = (set_uniq ? 0 : 1);
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
        }
    }
    tod_cepi();
    if (abortf == 2 || abortl == 3) {
        exclaim();
        if (abortf == 2)   {
            printf("Fatal errors in pedigree file %s.\n", mega2_input_files[0]);
        } else if (abortl == 3)   {
            printf("Fatal errors in locus file %s.\n", mega2_input_files[1]);
        }
        exclaim();
        EXIT(DATA_INCONSISTENCY);
    }

    // There are three cases in which genotypes are set to unknown...
    if (Mega2BatchItems[/* 26 */ Default_Reset_Invalid].item_read != 1 &&
        (hmend != 0 || imend != 0 || aexceed != 0)) {
        
        // It is possible to have reached the above menu from batch or interactive mode.
        // Here we note it in the batch file.
        // NOTE: The manual says "Yes" and "No" are the values, but 'copt' is a char.
        Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt = 'Y';
        batchf(Default_Reset_Invalid);
    }

    /* else abort = 1*/
    displayed_errors=0;
    Tod tod_imend("reset non mendelian");
    if (imend) {
        loc_list *curr, *next;
        int flush=0;

        /* set invalid genos to 0 0 */
        NonMendelianReset=1;

        /* First take out repeated elements */
        curr = invalid_genos;
        next = curr->next;
        while(next != NULL && curr != NULL) {
            if (curr->ped == next->ped && curr->locus == next->locus) {
                curr->next = next->next;
                free(next);
            } else {
                curr = curr->next;
            }
            if (curr->next != NULL) {
                next = curr->next;
            } else {
                break;
            }
        }
        mssgf("Setting Mendelianly-inconsistent genotypes to unknowns for these pedigree/Locus combinations:");
        strcpy(err_msg, "  ");
        invalidp = invalid_genos;
        while(invalidp != NULL) {
            ped = invalidp->ped;
            lloc = Top->LocusTop->Locus[invalidp->locus].linkage_loc_num;
            for (entry=0; entry < Top->PedTree[ped].EntryCnt; entry ++) {
                set_2alleles(Top->PedTree[ped].Entry[entry].LEntry->Marker, lloc, 
                             Top->LocusTop->Locus[invalidp->locus].linkage_loc_rec, 0, 0);
            }
            grow(err_msg, " %s: %s  ",
                 Top->PedTree[ped].Name,
                 LTop->Locus[invalidp->locus].Name);

            flush++;
            if (flush == 3) {
                SUPPRESS_MSSG(displayed_errors)
                    if (displayed_errors > MAX_PED_ERRORS) {
                        Display_Messages = 0;
                    }
                displayed_errors++;
                mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
            }
            invalidp = invalidp->next;
        }

        if (flush) {
            mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
        }
    } else {
        NonMendelianReset=0;
    }
    if (Display_Messages == 0) {
        printf("Check error logs for full list of pedigrees with reset genotypes.\n");
        draw_line();
        Display_Messages=1;
    }

    displayed_errors = 0;
    tod_imend();

    Tod tod_hmend("reset half typed");
    if (hmend) {
        int flush=0;
        mssgf("Setting half-typed genotypes to unknowns for these pedigree/person/locus combinations:");
        HalfTypedReset=1;
        strcpy(err_msg, "  ");
        /* set invalid genos to 0 0 */
        halftypedp = halftyped;
        while (halftypedp != NULL) {
            lloc = Top->LocusTop->Locus[halftypedp->locus].linkage_loc_num;
            entry = halftypedp->person;
            ped = halftypedp->ped;
            set_2alleles(LPedTop->Ped[ped].Entry[entry].Marker, lloc, 
                         Top->LocusTop->Locus[halftypedp->locus].linkage_loc_rec, 0, 0);
            grow(err_msg,
                 " %s, %s: %s  ",
                 Top->PedTree[ped].Name,
                 Top->PedTree[ped].Entry[entry].LEntry->OrigID,
                 LTop->Locus[halftypedp->locus].Name);
            flush++;
            if (flush == 2) {
                SUPPRESS_MSSG(displayed_errors)
                    if (displayed_errors > MAX_PED_ERRORS) {
                        Display_Messages=0;
                    }
                mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
                displayed_errors++;
            }
            halftypedp = halftypedp->next;
        }
        if (flush) {
            mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
        }

        abortf=1;
    } else {
        HalfTypedReset=0;
    }

    if (Display_Messages == 0) {
        printf("Check error logs for full list of reset genotypes.\n");
        draw_line();
        Display_Messages=1;
    }
    tod_hmend();

    Tod tod_oob("reset out of bounds");
    if (aexceed) {
        int flush=0;
        displayed_errors = 0;

        mssgf("Setting genotypes with out-of-bounds alleles to unknown ...");
        strcpy(err_msg, "  ");

        /* set invalid genos to 0 0 */
        outofboundsp = exceed_allcnt;
        while(outofboundsp != NULL) {
            lloc = Top->LocusTop->Locus[outofboundsp->locus].linkage_loc_num;
            entry = outofboundsp->person;
            ped = outofboundsp->ped;
            set_2alleles(LPedTop->Ped[ped].Entry[entry].Marker, lloc, 
                         Top->LocusTop->Locus[outofboundsp->locus].linkage_loc_rec, 0, 0);
            grow(err_msg, "Ped %s, Entry %s:Locus %s  ",
                 Top->PedTree[ped].Name,
                 Top->PedTree[ped].Entry[entry].LEntry->OrigID,
                 LTop->Locus[outofboundsp->locus].Name);
            flush++;
            if (flush == 2) {
                SUPPRESS_MSSG(displayed_errors)
                    if (displayed_errors > MAX_PED_ERRORS) {
                        Display_Messages = 0;
                    }
                mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
                displayed_errors++;

            }
            outofboundsp = outofboundsp->next;
        }

        if (flush) {
            mssgf(err_msg); flush=0; strcpy(err_msg, "  ");
        }
        abortf=1;
    }
    tod_oob();

    Tod tod_uniq("create_unique_id");
    if (set_uniq) {
        /* Create unique IDs */
        create_unique_ids(LPedTop, analysis);
    }
    tod_uniq();

    /* Store the resets instead of freeing */
    Tod tod_wrs("write_reset_summary");
    write_reset_summary(Top, invalid_genos, halftyped, exceed_allcnt,
                        PedStat.genotype_invalid, PedStat.halftyped,
                        PedStat.exceed_allcnt, LPedTop->UniqueIds);
    tod_wrs();
    if (Display_Messages == 0) {
        printf("Check error logs for full list of reset genotypes.\n");
        draw_line();
        Display_Messages=1;
    }
    if (PedStat.genotype_invalid || PedStat.halftyped || PedStat.exceed_allcnt) {
        log_line(mssgf);
        if (InputMode != INTERACTIVE_INPUTMODE && Mega2BatchItems[/* 26 */ Default_Reset_Invalid].item_read) {
            mssgf(
                "DEFAULT HANDLING MODE FOR INVALID GENOTYPES SPECIFIED IN BATCH FILE.\n");
        }
    }

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
    if (nonuniq) {
        if (set_uniq) {
            mssgf("New unique \"ID\" values generated by Mega2.\n");
        } else {
            mssgf("User chose not to generate new \"ID\" values.\n");
        }
    }

    if (PedStat.genotype_invalid || PedStat.halftyped ||
        PedStat.exceed_allcnt || nonuniq) {
        log_line(mssgf);
        printf("Check %s for details.\n", Mega2Err);
    }

    if (HasMarkers && PedStat.checked_sibship) {
        /* Inheritance checks may have been inadequate */
        draw_line();
        warnf("Inheritance checks may have been inadequate, ");
        warnf("Reason: One or more sibships with untyped parents.");
        warnf("Please verify with a program for checking Mendelian inconsistencies");
        warnf("such as Pedcheck.");
        log_line(mssgf);

    }
    Display_Messages=1;
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

    register ped_rec *Entry;

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

static void     create_unique_ids(linkage_ped_top *Top, analysis_type analysis)

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
