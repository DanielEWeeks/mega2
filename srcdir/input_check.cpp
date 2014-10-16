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

#if 0
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
#endif

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
    int             chk_aexceed, chk_imend, chk_hmend;
    int             set_uniq=0, nonuniq;

    register locus_top *LTop = Top->LocusTop;
    char            select_[10];
    char            toggle_str[50];

    /* These are set if the corresponding menu item is toggled */
    int             exit_select=-1, halftyped_select=-1;
    int             exceedall_select = -1, invalid_select=-1;
    int             uniq_select = -1;

    int             displayed_errors=0;
    int             plink_locus_num;

    int             num_mito_hetero = 0, num_mito_non_maternal = 0;

    clear_ped_status(&PedStat);

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

    Tod tod_mito("mito_transmission_report");
    for (chr=0; chr < main_chromocnt; chr++) {
        if (global_chromo_entries[chr] == MITO_CHROMOSOME) {
            mito_transmission_report(Top, &num_mito_hetero, &num_mito_non_maternal);
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
    displayed_errors = 0;
    tod_mito();

    Tod tod_iofc("input_observed_freq_check");
    if (num_reordered > 0 && analysis != TO_PLINK) {
        freq_mis = input_observed_freq_check(LPedTop, FreqMismatchThreshold);
        if (freq_mis > 0)
            abortf = imax(abortf, 1);
    } else {
        freq_mis=0;
    }
    tod_iofc();

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
            warnvf(" -> Invalid marker allele number or frequencies \n");
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
            int imendf = 0, hmendf = 0, aexceedf = 0;
            printf("Specify whether to reset poorly typed individuals and families:\n");
            while (select != 0) {
                strcpy(toggle_str, "");
                menu_item=0;
                draw_line();
                printf("0) Done with this menu - please proceed\n"); menu_item++;
                if (chk_hmend /* PedStat.halftyped */) {
                    if (analysis == TO_MENDEL || SIMWALK2(analysis)) {
                        printf(" NOTE: Proceeding will zero out all half-typed individuals.\n");
                    } else {
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
                if (nonuniq && analysis != TO_PAP && analysis != CRANEFOOT && analysis != IQLS) {
                    printf(" %d) Generate new unique IDs? [%s]\n", menu_item, yorn[set_uniq]);
                    grow(toggle_str, ", %d", menu_item);
                    uniq_select = menu_item++;
                }

                if (freq_mis) {
                    printf("Don't know how to correct for frequency mismatches.\n");
                }
//              if (PedStat.genotype_invalid || PedStat.halftyped || PedStat.exceed_allcnt) {
                    strcat(toggle_str, " toggle");
//              }
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
                        break;
                    } else if (select == invalid_select) {
                        imendf = (imendf ? 0 : 1);
                        break;
                    } else if (select == exceedall_select) {
                        aexceedf = (aexceedf ? 0 : 1);
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
            hmend = hmendf;
            imend = imendf;
            aexceed = aexceedf;
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
    Display_Errors=1;
    displayed_errors=0;
    bool first = true;

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
        mssgf("Setting any half-typed genotypes to unknowns for the indicated pedigree/person/locus combinations:");
        HalfTypedReset=0;

        for (ped = 0; ped < Top->PedCnt; ped++) {
            tod_hmend.reset();
            for (locus = 0; locus < LTop->LocusCnt; locus++) {
                stat = check_half_type(&(Top->PedTree[ped]), &PedStat, LTop,
                                       ped, locus, &displayed_errors, LPedTop->UniqueIds,
                                       &reset_fp, &first, hmend);
                abortf = imax(abortf, stat);
                if (stat) {
                    HalfTypedReset=1;
                }
            }
            tod_hmend("reset half typed");
        }
        /* set invalid genos to 0 0 */
        if (reset_fp != NULL)
            fclose(reset_fp);
        msgvf(" Status: Reset %d half-typed markers.\n", PedStat.halftyped);
        tod_hmend_all();
    } else {
        HalfTypedReset=0;
    }

    if (Display_Errors == 0) {
        printf("Check error logs for full list of half typed genotypes.\n");
        draw_line();
        Display_Errors=1;
    }
    displayed_errors = 0;

    Tod tod_xmend_all("reset ALL out-of-bound genotypes");
    Tod tod_xmend(20);
    int OOBReset=0;
    if (/*chk_*/aexceed) {
        FILE *reset_fp = NULL;
        mssgf("Setting any genotypes with out-of-bounds alleles to unknown ...");

        for (ped = 0; ped < Top->PedCnt; ped++) {
            tod_xmend.reset();
            for (locus = 0; locus < LTop->LocusCnt; locus++) {
                stat = check_out_of_bounds(&(Top->PedTree[ped]), &PedStat, LTop,
                                           ped, locus, &displayed_errors, LPedTop->UniqueIds,
                                           &reset_fp, &first, aexceed);
                abortf = imax(abortf, stat);
                if (stat) {
                    OOBReset=1;
                }
            }
            tod_xmend("reset out-of-bound genotypes");
        }
        /* set invalid genos to 0 0 */
        if (reset_fp != NULL)
            fclose(reset_fp);
        msgvf(" Status: Reset %d out-of-bound allele values.\n", PedStat.exceed_allcnt);
        tod_xmend_all();
    } else {
        OOBReset=0;
    }

    if (Display_Errors == 0) {
        printf("Check error logs for full list of out-of-bound genotypes.\n");
        draw_line();
        Display_Errors=1;
    }
    displayed_errors = 0;

    Tod tod_imend_all("reset ALL non mendelian");
    Tod tod_imend(20);
    if (/*chk_*/imend) {
        FILE *reset_fp = NULL;
        int rm;
        NonMendelianReset=0;
        for (ped = 0; ped < Top->PedCnt; ped++) {
            ped_rec **Sibs = CALLOC((size_t) Top->PedTree[ped].EntryCnt, ped_rec *);
            tod_imend.reset();
            for (locus = 0; locus < LTop->LocusCnt; locus++) {
                if (LTop->Locus[locus].chromosome == MITO_CHROMOSOME) {
                    continue;           /* Continue without inheritance checks */
                }
                rm = check_invalid_fam(&(Top->PedTree[ped]), &PedStat, LTop,
                                       ped, locus, &displayed_errors, LPedTop->UniqueIds,
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
                            LTop->Locus[locus].Name);

                    NonMendelianReset = 1;
                    if (imend) {
                        lloc = LTop->Locus[locus].linkage_loc_num;

                        for (entry=0; entry < Top->PedTree[ped].EntryCnt; entry ++) {
                            /* set invalid genos to 0 0 */
                            set_2alleles(Top->PedTree[ped].Entry[entry].LEntry->Marker, lloc, 
                                         Top->LocusTop->Locus[locus].linkage_loc_rec, 0, 0);
                        }
                    }
                }
            }
            tod_imend("reset non mendelian");
            free(Sibs);
        }
        if (reset_fp != NULL) {
            fclose(reset_fp);
        }
        tod_imend_all();
    } else {
        NonMendelianReset=0;
    }

    if (Display_Errors == 0) {
        printf("Check error logs for full list of pedigrees with reset genotypes.\n");
        draw_line();
        Display_Errors=1;
    }
    displayed_errors = 0;

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

    Tod tod_uniq("create_unique_id");
    if (set_uniq) {
        /* Create unique IDs */
        create_unique_ids(LPedTop, analysis);
    }
    tod_uniq();

    if (PedStat.genotype_invalid || PedStat.halftyped || PedStat.exceed_allcnt) {
        log_line(mssgf);
        if (InputMode != INTERACTIVE_INPUTMODE && Mega2BatchItems[/* 26 */ Default_Reset_Invalid].item_read) {
            mssgf(
                "DEFAULT HANDLING MODE FOR INVALID GENOTYPES SPECIFIED IN BATCH FILE:");
            int copt = Mega2BatchItems[/* 26 */ Default_Reset_Invalid].value.copt;
            if ((copt == 'y' || copt == 'Y')) 
                msgvf(" Reset allele to 0/0\n\n");
            else
                msgvf(" Leave allele alone\n\n");
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
        warnvf("Mega2 also found these problems/errors in input data (see MEGA2.ERR for details): \n");
        if (PedStat.halftyped > 0) {
            warnvf(" -> %d Half-typed genotypes\n", PedStat.halftyped);
            if (analysis == TO_MENDEL || SIMWALK2(analysis)) {
                warnvf("       Half-typed genotypes will be set to unknown \n");
                warnvf("       as required by %s\n", ProgName);
            }
        }
        if (PedStat.genotype_invalid > 0) {
            warnvf(" -> %d Mendelian-inconsistent genotypes\n", PedStat.genotype_invalid);
        }

        if (PedStat.exceed_allcnt > 0) {
            warnvf(" -> %d Genotype allele numbers exceed locus allele count\n",
                   PedStat.exceed_allcnt);
        }
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
    Display_Errors=1;
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
