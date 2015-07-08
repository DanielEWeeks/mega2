/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

/* Changes needed when adding a new output option:

   There are various functions in this file that do analysis-specific
   locus selection checks:

   1) no_trait_allowed():
   This function returns FALSE(0), if the analysis option requires at
   least one trait, TRUE (1) otherwise.  Add your KEYWORD to the first
   block of the case statement, if it needs to have traits.

   2) no_aff_trait_allowed(): Ignore this for the time being, this is not
   used.

   3) allow_trait_combination(): This function return TRUE if traits are
   allowed to be combined in the output as a single set of files, and
   traits do not need to come first in the locus list. Add
   your KEYWORD to this list if your analysis option allows this.

   4) allow_covariates(): If QTLs can be specially designated as
   covariates for analysis by your program, add your KEYWORD to the first
   block of the case.

   #ifdef DEFUNCT
   5) get_trait_positions(): For some analysis options, a single trait can be
   placed between markers, using the reordering menu option 2, and this
   trait can have an unknown position. This function checks to see that
   there is only one such unmapped locus. Add your KEYWORD to the case
   block at the top, if this is appropriate.
   #endif

   6) verify_markers(): If your analysis option has some special
   requirements besides 1-5 on locus selection, then add a special
   checking routine to this function. The variable reset is a flag which
   is to be set by this checking routine if the selection is not
   satisfactory.

*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "common.h"
#include "typedefs.h"

#include "annotated_ped_file.h"
#include "annotated_ped_file_ext.h"
#include "batch_input_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "output_routines_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
 annotated_ped_file_ext.h:  copy_exmap_locmap
        batch_input_ext.h:  batchf
     create_summary_ext.h:  log_quant_selections
     error_messages_ext.h:  errorf mssgf my_calloc my_realloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  check_map_positions loc_type_name log_marker_selections
        grow_string_ext.h:  grow
            linkage_ext.h:  clean_reordered_markers free_all_from_llocusrec
    output_routines_ext.h:  locus_name_width
         user_input_ext.h:  ReOrderMenu define_affection_labels set_missing_quant_input
              utils_ext.h:  EXIT chomp draw_line log_line press_return print_only invalid_value_field
        write_files_ext.h:  write_quant_stats
*/


/* ----------exports-------------- */

linkage_ped_top *ReOrderLoci(linkage_ped_top *Top, int *numchr,
			     analysis_type *analysis);
int             linkage_set_locus_order_new(int num_order, int *order);
int             get_chromosome_list(linkage_locus_top *LTop, int *list1,
				    int *counts, int exit_upon_no_chr,
				    analysis_type analysis);
int             display_selections(linkage_ped_top *Top, int *entries,
				   int num_entries, char *asterisks,
				   char *scanned_str, int traits_only,
				   int map_num);
static int       *check_for_duplicates(int *entrycount, int *number);
int              parse_string(char *entrydummy, int *number,
                              int *num_entries, int num_total_loci, int require_e);

/*---------------internal functions ---------*/
int *init_trp (void);
static int      ReOrderLociByChromosome(linkage_ped_top *Top, int *numchr,
					int get_or_set);
static int      ReOrderLociByPositionNumber(linkage_ped_top *Top,
					    int numchr, int *num_loci,
					    int **selected_loci,
					    analysis_type analysis,
					    int get_or_set, int map_num);
static int     chromo_and_loci_selection(linkage_ped_top *Top,
					 int *num_chromo, int **selected_chromo,
					 int get_ot_set);
static int      select_trait_loci(linkage_ped_top *Top, analysis_type analysis);
#ifdef DEFUNCT
int             ReOrderMappedLoci(linkage_ped_top *Top, int *numchr);
int             linkage_set_locus_order(linkage_ped_top *Top, int *order);
static void     linkage_swap_loci(linkage_locus_rec *locus1,
				  linkage_locus_rec *locus2);
#endif
static int    no_trait_allowed(analysis_type analysis, int mssg);
/* static int    no_aff_trait_allowed(analysis_type analysis, int mssg); */
static int    quant_allowed(linkage_locus_rec Locus, analysis_type analysis);
static int    check_trait_selection(linkage_ped_top *Top,
				    analysis_type analysis, int num_select,
				    int *trait_order, int marker_item);
static int    check_ghmlb_affs(linkage_ped_top *Top, int *traits, int num_affs,
			       int mrk_pos);
static int    all_same_type(linkage_ped_top *Top, int *traits, int num_affs,
			    int markers_item);
static int check_gh_quants(linkage_ped_top *Top, int *traits,
			   int num_trts, int mrk_pos);
static int check_simwalk2_affs(analysis_type anal, int num_select,
			       int *slct,
			       int marker_pos, int sw2check);
static int add_selections(int num_new, int *newer, int *num_old,
			  int *old, int total);
static int allow_trait_combination(analysis_type analysis);
static int allow_covariates(analysis_type analysis);
static int index_selections(int num_select, int *selections,
			    int num_values, int *values,
			    int reverse, int marker_item,
			    int base);
static void display_trait_names(int num_tr, int *trs,
				linkage_locus_top *LTop,
				int offset, int *indexes,
				int base, int marker_item);
static void remove_from_covar(int num_select, int *slct,
			      int *num_covar, int *cvariates);

#ifdef DEFUNCT
static void get_trait_positions(linkage_locus_top *LTop,
				analysis_type analysis);
#endif
static int   ReOrderMappedLoci_new(linkage_ped_top *Top,
				   int loc_type, int **locus_order,
				   int numchr, int *chromolist);
#if 0
static void get_map_num(int *map_num, ext_linkage_locus_top *EXLTop,
			analysis_type analysis);
#endif
static void new_chromo_loci_counts(int num_chromo_selected, int *chromosome_list);

static void reset_ychrom_positions(linkage_ped_top *Top);

typedef struct _tmpnum {
    int number;
    double phys;
    double index;
} tmpnum;

#define INDX_INC(num, indx, i, marker_item)                             \
    for (i=0; i < num; i++) if (indx[i] != marker_item) {indx[i]++;}

/*==============================*/

static int compare_index(const void *ptr1, const void *ptr2)

{
    const tmpnum *t1, *t2;

    t1 = (const tmpnum *) ptr1;
    t2 = (const tmpnum *) ptr2;

    if (t1->index < t2->index) {
        return -1;
    }
    else if (t1->index > t2->index) {
        return 1;
    }
    else {
        return (t1->phys < t2->phys) ? -1 : ((t1->phys == t2->phys) ? 0 : 1);
    }
}

#ifdef DEFUNCT
/* linkage_swap_loci()
 *
 * Slave routine for set_locus_order(). Swaps data in locus records.
 */

/*-----------------------------------------------------------------------
  If two loci on the same chromosome have the same distances, then adjust
  distances by a very small amount
  -------------------------------------------------------------------------*/

static void     linkage_swap_loci(linkage_locus_rec *locus1,
				  linkage_locus_rec *locus2)
{
    register linkage_locus_rec tmplocus;

    tmplocus = *locus1;
    *locus1 = *locus2;
    *locus2 = tmplocus;
}
#endif

/* changed this routine, so that limit on the list size is set by the
   maximum number of chromosome present, sets a global MaxChromo,
   see common.h
   The above comment doesn't make sense as MaxChromo is already set prior
   to this routine, and list and counts are allocated outside as well.
*/
int get_chromosome_list(linkage_locus_top *LTop, int *local_list,
                        int *counts, int exit_upon_no_chr,
                        analysis_type analysis)
{
    int i, j, ui, total_list = 0;
    int new_chr, chr_index;
    SECTION_LOG_INIT(unmapped_errors);
    SECTION_LOG_INIT(invalid_chromosome);
    ui = 0;

    for (j = 0; j < LTop->LocusCnt; j++) {
        if (LTop->Locus[j].Class == TRAIT)
            continue;

        if (LTop->Marker[j].chromosome == UNKNOWN_CHROMO) {

            SECTION_LOG(unmapped_errors);
            warnvf("Locus %s is unmapped\n", LTop->Marker[j].Name);
            unmapped_markers[ui++] = j;
        } else if (LTop->Marker[j].chromosome == MISSING_CHROMO) {
	    // The entry for this locus (marker assumed) was not found in the map file.
            // unmapped_markers[ui++] = j;
            continue;
        } else if (LTop->Marker[j].chromosome < 1) {
            SECTION_LOG(invalid_chromosome);
            warnvf("Invalid chromosome number %d on locus %s.\n",
                    LTop->Marker[j].chromosome, LTop->Marker[j].Name);
            if (exit_upon_no_chr) {
                errorvf("Please correct and restart mega2!\n");
                EXIT(DATA_TYPE_ERROR);
            } else {
                /* 	if (ALLOW_NO_CHR(analysis)) { */
                /* 	  warnf("Chromosome number set to UNKNOWN."); */
                /* 	  ui++; */
                /* 	} */
                /* 	else { */
                warnf("Locus will be omitted from analysis.");
                /* 	} */
                /* 	LTop->Marker[j].chromosome = UNKNOWN_CHROMO; */
            }
        } else {
            /* fill chromosome list and locus counts in input order */
            new_chr = 1;
            if (total_list > 0) {
                for (i = 0; i < total_list; i++) {
                    if (LTop->Marker[j].chromosome == local_list[i]) {
                        chr_index = i;
                        new_chr = 0;
                        break;
                    }
                }
            }
            if (new_chr) {
                local_list[total_list] = LTop->Marker[j].chromosome;
                counts[total_list] = 1;
                total_list++;
            } else {
                counts[chr_index]++;
            }
        }
    }
    SECTION_LOG_FINI(unmapped_errors);
    SECTION_LOG_FINI(invalid_chromosome);
    /* Now add counts for any markers with unknown chromosome on the end */
    if (ui > 0) {
        counts[total_list] = ui;
        local_list[total_list] = UNKNOWN_CHROMO;
        total_list++;
    }

    /* Now convert the counts to cumulative counts */
    for (i = 1; i < total_list; i++) {
        counts[i] += counts[i - 1];
    }

    return total_list;
}

/**
   @brief Determine if an integer is in a string of integers

   Given a string 'str' that contains integers separated by white space,
   determine if the integer 'i' exists in that string.

   @param[in] str    string that conrains integers separated by white space,
   @param[in] i      an integer
   @return           1 if yes, 0 if no
 */
static int is_int_in_string(const char *str, const int i)
{
    char *entry = (char *)str;
    char *endptr;
    int num, num2;

    while (*entry != '\0' && *entry != 'e') {
        num  = strtol(entry, &endptr, 10);
        entry = endptr;
        while ((*entry != '\0' && *entry != 'e') && (isspace((int)*entry))) entry++;
        if (*entry == ',') {
            if (num == i) return 1;
            entry++;
        } else if (*entry == '-' && *(entry+1) != 0) {
            num2  = strtol(entry + 1, &endptr, 10);
            if (i >= num && i <= num2) return 1;
            entry = endptr;
        } else if (num == i) return 1;
    }
    return 0;
}

void get_trait_list(linkage_locus_top *LocusTop, const int check_traits_combine)
{
    int i, j, *affec = NULL, *quant = NULL;
    int   num_affec=0, num_quant=0;
    int *covp, *trp, *affp = NULL, *quantp = NULL;

    if (Mega2Status == INSIDE_ANALYSIS) {
        if (num_traits > 0) {
            free(global_trait_entries);
        }
        if (num_covariates > 0) {
            free(covariates);
        }
    }

    if (Mega2Status != TRAIT_SELECTED) {
        num_traits=0; num_covariates=0;
        /* check the total number of traits */
        for (i=0; i<LocusTop->LocusCnt; i++) {
            if (LocusTop->Locus[i].Class == COVARIATE) {
                num_covariates++;
            }
            else if (LocusTop->Locus[i].Type == AFFECTION) {
                num_traits++;
                num_affec++;
            }
            else if (LocusTop->Locus[i].Type == QUANT) {
                num_quant++;
                num_traits++;
            }
        }

        if (num_traits == 0) {
            return;
        }

        global_trait_entries = CALLOC((size_t) num_traits, int);
        trp = &(global_trait_entries[0]);
        if (num_affec > 0) {
            affec = CALLOC((size_t) num_affec, int);
            affp = &(affec[0]);
        }
        else {
            affp=NULL;
        }
        if (num_quant > 0) {
            quant = CALLOC((size_t) num_quant, int); quantp = &(quant[0]);
        }
        else {
            quantp = NULL;
        }
        if (num_covariates > 0) {
            covariates = CALLOC((size_t) num_covariates, int); covp = &(covariates[0]);
        }
        else {
            covp = NULL;
        }

        for (i=0; i<LocusTop->LocusCnt; i++) {
            if (LocusTop->Locus[i].Class == COVARIATE) {
                *covp = i; covp++;
            }
            else if (LocusTop->Locus[i].Type == AFFECTION) {
                *trp = i; trp++;
                *affp = i; affp++;
            }
            else if (LocusTop->Locus[i].Type == QUANT) {
                *trp = i; trp++;
                *quantp = i; quantp++;
            }
        }
    } else {
        /* global_trait_entries now contains selected traits
           and num_traits is the number of traits selected
        */

        affec = CALLOC((size_t) num_traits, int);
        quant = CALLOC((size_t) num_traits, int);

        /* fill in the quants and affs */

        for (i=0; i < num_traits; i++) {
            j=global_trait_entries[i];
            if (j == -1) {
                continue;
            }
            if (LocusTop->Locus[j].Type == AFFECTION &&
                LocusTop->Locus[j].Class != COVARIATE) {
                affec[num_affec] = j; num_affec++;
            }
            else if (LocusTop->Locus[j].Type == QUANT &&
                     LocusTop->Locus[j].Class != COVARIATE) {
                quant[num_quant] = j; num_quant++;
            }
        }
        mssgf("After selecting traits and covariates");
    }

    if (num_traits > 0) {
        int numt;
        if (Mega2Status != TRAIT_SELECTED) {
            numt = num_traits;
        }
        else {
            numt = ((LoopOverTrait == 0)? num_traits-1 : num_traits);
        }
        sprintf(err_msg, "%d trait %s ", numt, ((numt > 1)? "loci" : "locus"));
        
        // The sentinal associated with 'Traits_Combine' is confusing.
        // If it is not present in the list, then the user will find that no loci are selected
        // and not really be sure why.
        // NOTE: This is a different issue from the 'trait number N is out of the bounds' error
        // message that is generated elsewhere in this code.
        if (check_traits_combine &&
            Mega2BatchItems[/* 15 */ Traits_Combine].items_read &&
            is_int_in_string(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, numt+1) == 0) {
            errorf("The batch file item 'Traits_Combine' must contain an item numbered N+1");
            errorvf("where N is the number of traits (here N+1=%d). For additional information\n", numt+1);
            errorf("please see the Mega2 manual section describing 'Details on batch file items'.");
            errorvf("Currently: Traits_Combine = %s\n", Mega2BatchItems[/* 15 */ Traits_Combine].value.name);
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
    }

    if (num_covariates > 0) {
        grow(err_msg, " and %d covariate(s) ",
             num_covariates);
    }

    mssgf(err_msg);

    if (num_affec > 0) {
        sprintf(err_msg, "      %d Affection status %s", num_affec,
                ((num_affec > 1)? "loci: " : "locus: "));
        mssgf(err_msg);
        display_trait_names(num_affec, affec, LocusTop,
                            16, NULL, 0, LocusTop->PhenoCnt+1);
        newline;
    }

    if (num_quant > 0) {
        sprintf(err_msg, "      %d Quantitative %s", num_quant,
                ((num_quant > 1)? "loci: " : "locus: "));
        mssgf(err_msg);
        display_trait_names(num_quant, quant, LocusTop,
                            16, NULL, 0, LocusTop->PhenoCnt+1);
    }

    /*  log_line(mssgf); */

    if (num_quant > 0) {
	HasQuant=1;
        num_quantitative = num_quant;
        free(quant);
    }
    else {
        HasQuant=0;
        num_quantitative = 0;
    }

    if (num_affec > 0) {
        HasAff=1;      free(affec);
        num_affection = num_affec;
    }
    else {
        HasAff = 0; num_affection = 0;
    }
    strcpy(err_msg, "");
}

static void display_trait_names(int num_tr, int *trs,
				linkage_locus_top *LTop,
				int offset, int *indexes,
				int base,
				int marker_item)

{
    int i, j, k;
    int col=0;

    for (i=0; i < offset; i++) {
        msgvf(" ");
    }
    for (i=0; i < num_tr; i++) {
        if (trs[i] < marker_item) {
            if (indexes == NULL) {
                k = trs[i]-base;
            }
            else {
                k = indexes[trs[i]-base];
            }
            msgvf("%s", LTop->Locus[k].Name);
            col += (int) strlen(LTop->Locus[k].Name);
        }
        else {
            msgvf("[MARKERS]");
            col += 9;
        }

        if ((col + offset) >= 70) {
            msgvf("\n"); col=0;
            for (j=0; j < offset; j++) logvf(" ");
        }
        else {
            if ( i < (num_tr - 1)) {
                msgvf(" ");
                col +=1;
            }
        }
    }
    msgvf("\n");
    return;
}

//
// Parse the chromosome string, returning a positive integer greater then or equal to
// zero for a valid chromosome. Returns -1 on error.
int STR_CHR(const char *dummy)
{
    int chr;
    
    if (!strcasecmp(dummy, "X")         || !strcasecmp(dummy, SEX_CHROMOSOME_STR)) {
        chr=SEX_CHROMOSOME;
    } else if  (!strcasecmp(dummy, "Y") || !strcasecmp(dummy, MALE_CHROMOSOME_STR)) {
        chr=MALE_CHROMOSOME;
    } else if (!strcasecmp(dummy, "XY") || !strcasecmp(dummy, PSEUDO_X_STR)) {
        chr=PSEUDO_X;
    } else if (!strcasecmp(dummy, "MT") || !strcasecmp(dummy, MITO_CHROMOSOME_STR)) {
        chr=MITO_CHROMOSOME;
        // PLINK allows the value of '0' to represent "unplaced"
    } else if (!strcasecmp(dummy, "U") || !strcasecmp(dummy, "0") || !strcasecmp(dummy, UNKNOWN_CHROMO_STR)) {
        chr=UNKNOWN_CHROMO;
    } else {
        // Error out if there is something remaining after parsing what should be an integer,
        // the integer is less than zero or greater than the maximum chromosome that is permitted.
        char *end;
        chr = (int)strtol((const char *)dummy, &end, 10);
        if (strlen(dummy) == 0 || strlen(end) != 0 || errno == ERANGE || chr < 0 || chr > lastautosome+4) chr = -1; // ERROR
        //if (strlen(dummy) == 0 || strlen(end) != 0 || errno == ERANGE || chr < 0) chr = -1; // ERROR
        //chr=atoi(dummy);
    }
    
    return chr;
}

linkage_ped_top *ReOrderLoci(linkage_ped_top *Top, int *numchr,
			     analysis_type *analysis)
{
    int selection, option, map_num = 0; //silly compiler
    int loc_type = 1;
    int i, chr_valid = 0;
    int num_chromo, *selected_chromosomes = NULL;
    int num_loci, *selected_loci = NULL, has_quant;

    ManualReorder = 0;

    // The user chooses the gentic_distance_index in the routine
    // user_input.c:get_genetic_distance_index(ext_linkage_locus_top *EXLTop)
    // which is run in mega2.c:main() before this function is called...
    if (genetic_distance_index != -2) {
        map_num = genetic_distance_index;
    } else {
        errorvf("For the analysis type specified, a Genetic Map was required, but none was chosen.\n");
        EXIT(EARLY_TERMINATION);
    }
    
    if (*analysis == QUANT_SUMMARY && InputMode == BATCH_FILE_INPUTMODE) {
        has_quant = 0; /* 0 = FALSE */
        for (i = 0; i < num_traits; i++) {
            if (Top->LocusTop->Locus[global_trait_entries[i]].Type == QUANT) {
                has_quant = 1; /* 1 = TRUE */
                break;
            }
        }
        if (has_quant) {
            /* has_quant == TRUE if there is at least one quantitative trait locus
             * in the set of traits in the input files.
             */
            /*	printf("This option requires only trait loci to be selected.\n");
                printf("         Use options 1 or 2 to select trait loci\n"); */
            loc_type = -1; /* loc_type == -1: allow user to select only QUANT loci */
        } else {
            errorf("No quantitative loci in input data, terminating Mega2.\n");
            EXIT(EARLY_TERMINATION);
        }
        /* Mega2BatchItems[/ * 15 * / Traits_Combine] is 'Traits_Combine'; Mega2BatchItems[/ * 12 * / Trait_Single] is 'Trait_Single' */
        if (Mega2BatchItems[/* 15 */ Traits_Combine].items_read || Mega2BatchItems[/* 12 */ Trait_Single].items_read) {
            batchREORDER = 1; /* non-zero is true */
            /*  user-specified marker order on current chromosome */
            ReOrderLociByPositionNumber(Top, -1, &num_loci, &(selected_loci), *analysis, 1, 0);
            option = 2;
            ManualReorder = 1;
            *numchr = -1;
        } else {
            batchREORDER = 0; /* zero is false */
        }
    }

    if (*analysis != QUANT_SUMMARY) {
        // SexLinked == 0  did not find m/f-chromosome...
        // SexLinked == 1  found m/f but no autosomes...
        // SexLinked == 2  found m/f and autosomes...
//SL
//        Top->LocusTop->SexLinked = x_linked_check(NumChromo, global_chromo_entries, *analysis);
    }

    /* either not quant summary or not batch Reorder */
    if (batchREORDER) {
        if (NumChromo == 0 || *numchr <= 0) {
            ManualReorder = 1;
            if (!(*analysis == QUANT_SUMMARY && InputMode == BATCH_FILE_INPUTMODE)) {
                if (*analysis == QUANT_SUMMARY)
                    loc_type = -1;
                else
                    loc_type = 0;
                /* Batch Item 11 is 'Loci_Selected'; Batch Item 14 is 'Traits_Loop_over' */
                if (Mega2BatchItems[/* 11 */ Loci_Selected].items_read || Mega2BatchItems[/* 14 */ Traits_Loop_Over].items_read) {
                    /* keep selected trait loci */
                    ReOrderLociByPositionNumber(Top, 0, &num_loci, &selected_loci, *analysis, 1, 0);
                } else {
                    /* keep all trait loci */
                    ReOrderMappedLoci_new(Top, loc_type, &selected_loci, 0, NULL);
                }
            }
            HasMarkers = 0;
        } else {
            if (Top->EXLTop != NULL) {
#ifdef USEOLDMAPCODE
                if (ITEM_READ(41)) {
                    if (Mega2BatchItems[/* 41 */ Output_Map_Num].value.option >= 1
                        && Mega2BatchItems[/* 41 */ Output_Map_Num].value.option
                        <= Top->EXLTop->MapCnt) {
                        copy_exmap_locmap(Top->LocusTop, Top->EXLTop,
                                          Mega2BatchItems[/* 41 */ Output_Map_Num].value.option - 1
					  );
                        map_num = Mega2BatchItems[/* 41 */ Output_Map_Num].value.option - 1;
                    } else {
                        invalid_value_field(41);
                    }
                } else {
#endif /* USEOLDMAPCODE */
                    // The user chooses the gentic_distance_index in the routine
                    // user_input.c:get_genetic_distance_index(ext_linkage_locus_top *EXLTop)
                    // which is run in mega2.c:main() before this function is called...
                    if (genetic_distance_index != -2) {
                        map_num = genetic_distance_index;
                        copy_exmap_locmap(Top->LocusTop, Top->EXLTop, map_num);
                    } else {
                        errorvf("For the analysis type specified, a Genetic Map was required, but none was chosen.\n");
                        EXIT(EARLY_TERMINATION);
                    }
#ifdef USEOLDMAPCODE
                }
#endif /* USEOLDMAPCODE */


#ifndef HIDESTATUS
                // Now we break this down into specifics of the type (possibly auto, sex, female).
                mssgvf("Selected map %s.\n", Top->EXLTop->MapNames[map_num]);
#endif
            } else {
                // if there is no EXLTop structure, just assume one map, the first?
                map_num = 0;
            }

            if (Mega2BatchItems[/* 7 */ Chromosome_Single].items_read) {
                /* single chromosome option */
                /* Find out if chromosome is valid only if
                   there is at least one chromosome*/
                option = 1;
                if (Mega2BatchItems[/* 7 */ Chromosome_Single].value.option == UNKNOWN_CHROMO
                    && ALLOW_NO_CHR(*analysis)) {
                    chr_valid = 1;
                    *numchr = UNKNOWN_CHROMO;
                } else {
                    for (i = 0; i < NumChromo; i++) {
                        if (global_chromo_entries[i] == Mega2BatchItems[/* 7 */ Chromosome_Single].value.option) {
                            chr_valid = 1;
                            *numchr = global_chromo_entries[i];
                            break;
                        }
                    }
                }

                if (!chr_valid) {
                    invalid_value_field(7); /* This terminates if invalid */
                }

                mssgvf("Selected chromosome %d\n", *numchr);

                if (Mega2BatchItems[/* 11 */ Loci_Selected].items_read) {
                    /*  user-specified marker order  on current chromosome */
                    ReOrderLociByPositionNumber(Top, *numchr, &num_loci,
                                                &(selected_loci), *analysis, 1, map_num);
                    option = 2;
                    ManualReorder = 1;
                    mssgf("Locus selections read in from batch file.");
                } else {
                    ReOrderLociByChromosome(Top, numchr, 1);
                    ManualReorder = 0;
                }
            } else if (Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].items_read) {
                char prchr[4];
                /* Multiple chromosomes */
                option = 3;
                ManualReorder = 0;
                if (chromo_and_loci_selection(Top, &num_chromo,
                                              &selected_chromosomes, 1)) {
                    invalid_value_field(8);
                }
                if (main_chromocnt == 1)
                    *numchr = global_chromo_entries[0];
                mssgf("You have selected the following chromosome(s) :");
                strcpy(err_msg, " ");
                for (i = 0; i < main_chromocnt; i++) {
                    CHR_STR(global_chromo_entries[i], prchr);
                    grow(err_msg, " %s", prchr);
                }
                mssgf(err_msg);
                log_line(mssgf);
            } else {
                /* Both Chromosomes_single and Chromosomes_Multiple_Num
                   are missing */
                errorvf("Both %s and %s are missing from batch file.\n",
			C(Mega2BatchItems[/* 7 */ Chromosome_Single].keyword),
			C(Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].keyword)
			);
                EXIT(BATCH_FILE_ITEM_ERROR);
            }
        }
    } else { /* not batch reorder */
        if (*analysis == QUANT_SUMMARY || NumChromo == 0) {
            if (*analysis == QUANT_SUMMARY) {
                has_quant = 0;
                for (i = 0; i < num_traits; i++) {
                    if (Top->LocusTop->Locus[global_trait_entries[i]].Type == QUANT) {
                        has_quant = 1;
                        break;
                    }
                }
                if (has_quant) {
                    printf(
                        "This option requires only quantitative trait loci to be selected.\n");
                    printf("         Use options 1 or 2 to select trait loci\n");
                    loc_type = -1;
                } else {
                    errorf("No quantitative loci in input data, terminating Mega2.\n");
                    EXIT(EARLY_TERMINATION);
                }
            } else {
                /* else if there are no chromosomes */
                warnf("No markers found with an assigned chromosome number.\n         Use options 1 or 2 to select trait loci\n");
                loc_type = 0;
            }
            /* Only keep all trait loci */
            selection = 1;
            option = 1;
            num_loci = ReOrderMappedLoci_new(Top, loc_type, &selected_loci, 0, NULL);
            while (option != 0) {
                ReOrderMenu(0, NULL, &option);
                switch (option) {
                case 0:
                    break;
                case 1:
                    num_loci = ReOrderMappedLoci_new(Top, loc_type,
                                                     &selected_loci, 0, NULL);
                    /* Notes on values of *numchr,
                       -2 = only affection - not used currently,
                       -1 = only QUANT,
                       0 = all traits
                       >0 = markers and traits
                    */
                    selection = 1;
                    break;
                case 2:
                    /* Only selected traits */
                    ReOrderLociByPositionNumber(Top, loc_type, &num_loci,
                                                &selected_loci, *analysis, 0, 0);
                    if (num_loci == 0) {
                        printf("ERROR: No loci selected!\n");
                    }
                    selection = 2;
                    break;
                default:
                        break;
                }
            }

            /* set the batch items correctly */
            /*	Mega2BatchItems[/ * 7 * / Chromosome_Single].value.option = loc_type; */
            /* 	strcpy(Mega2BatchItems[/ * 11 * / Loci_Selected].value.name, ""); */
            strcpy(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, "");
            for (i = 0; i < num_loci; i++) {
                grow(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, " %d", selected_loci[i] + 1);
            }
            grow(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, " %d", num_traits + 1);
            strcat(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, " e");

            ReOrderLociByPositionNumber(Top, loc_type, &num_loci,
					&(selected_loci), *analysis, 1, 0);
            /*      } */
            option = 2;
            ManualReorder = 1;
            HasMarkers = 0;
        } else { /* analysis not quant summary, chr-num non-zero */
            /* initialize chromosome number to the first chromosome */
            selection = 1;
            option = 1;
            if (selected_chromosomes != NULL) {
                free(selected_chromosomes);
                selected_chromosomes = NULL;
            }
            selected_chromosomes = CALLOC((size_t) NumChromo, int);
            num_chromo = 1;
            for (i = 0; i < NumChromo; i++)
                selected_chromosomes[i] = global_chromo_entries[i];
            
            //map_num = 0;
            if (Top->EXLTop != NULL) {
#if 0
                if (*analysis == TO_PLINK || *analysis == IQLS) {
                    /* Set the default map to the first physical map found */
                    for (i = 0; i < Top->EXLTop->MapCnt; i++) {
                        if (Top->EXLTop->map_functions[i] == 'p') {
                            map_num = i;
                            break;
                        }
                    }
                }
#endif /* 0 */
                if (*analysis == TO_PLINK || *analysis == IQLS) {
                    // in these analysis, a physical map is required...
                    if (base_pair_position_index != -2) {
                        // -2 == None
                        // The user chooses the base_pair_position_index in the routine
                        // user_input.c:get_base_pair_position_index(ext_linkage_locus_top *EXLTop)
                        // which is run in mega2.c:main() before this function is called.
                        // for PLINK and IQLS this function will force them to choose a position.
                        map_num = base_pair_position_index;
                    } else {
                        errorvf("For the analysis type specified, a Physical Map was required, but none was chosen.\n");
                        EXIT(EARLY_TERMINATION);
                    }
                }
            }

	  if (map_num >= 0) {
            while (selection != 0) {
                // ==========================================================
                // Locus Reordering Menu
                // 0) Done with this menu - please proceed.
                // *1) Select all loci in map order on chromosome: 1
                // 2) Select by locus number.
                // 3) Select loci on multiple chromosomes[1 ]
                // Enter 0 - 3 >
                ReOrderMenu(num_chromo, selected_chromosomes, &selection);
                switch (selection) {
                case 0: // take the default (1)
                    break;
                case 1: // take the first chromosome...
                    ReOrderLociByChromosome(Top, &(selected_chromosomes[0]), 0);
                    option = 1;
                    break;
                case 2: // Selection by locus number...
                    *numchr = selected_chromosomes[0];
                    ReOrderLociByPositionNumber(Top, *numchr, &num_loci,
                                                &(selected_loci), *analysis, 0, map_num);
                    option = 2;
                    if (num_loci > 0) {
                        /* this is the only case where markers can be out of map order */
                        if (*analysis == TO_ASPEX || *analysis == TO_SAGE4
                            || SIMWALK2(*analysis)) {
                            check_map_positions(Top, *numchr, 1);
                        }
                    }
                    ManualReorder = 1;
                    break;
                case 3: // Select loci on multiple chromosomes...
                    chromo_and_loci_selection(Top, &(num_chromo),
                                              &(selected_chromosomes), 0);
                    option = 3;
                    break;
                default:
                    break;
                }
            } //while (selection != 0) {

            if (Top->EXLTop != NULL) {
                copy_exmap_locmap(Top->LocusTop, Top->EXLTop, map_num);
#ifdef USEOLDMAPCODE
                if (InputMode == INTERACTIVE_INPUTMODE) {
                    Mega2BatchItems[/* 41 */ Output_Map_Num].value.option = map_num + 1;
                    batchf(Output_Map_Num);
                }
#endif /* USEOLDMAPCODE */
            }

          switch (option) {
              case 1:
                  ReOrderLociByChromosome(Top, &(selected_chromosomes[0]), 1);
                  ManualReorder = 0;
                  break;
              case 2:
                  *numchr = selected_chromosomes[0];
                  ReOrderLociByPositionNumber(Top, *numchr, &num_loci,
                                              &(selected_loci), *analysis, 1, map_num);
                  ManualReorder = 1;
                  break;
              case 3:
                  chromo_and_loci_selection(Top, &(num_chromo),
                                            &(selected_chromosomes), 1);
                  ManualReorder = 0;
                  
                  break;
              default:
                  break;
          }
	    }
	}

    } /* end of menu-based reordering */

    /* Take out markers that do not have any usable map positions */
    /*   Now switch XLINKED and YLINKED to NUMBERED, and set HasMarkers */

    /* get_trait_list(Top, 0); */

    if (*analysis == TO_ALLELE_FREQ || *analysis == TO_LIABLE_FREQ || *analysis
        == TO_HWETEST || *analysis == GENOTYPING_SUMMARY || *analysis
        == TO_PREST || *analysis == QUANT_SUMMARY) {

        if (*analysis == QUANT_SUMMARY) {
            if (InputMode == INTERACTIVE_INPUTMODE) {
                log_quant_selections(Top);
            }
            set_missing_quant_input(Top, *analysis);
            if (write_quant_stats(Top, *analysis) != 0)
                set_missing_quant_output(Top, *analysis);
        } else if (*analysis == TO_PREST || *analysis == TO_HWETEST
                   || *analysis == GENOTYPING_SUMMARY) {
            /* get rid of the traits */
            free(global_trait_entries);
            global_trait_entries = CALLOC((size_t) 1, int);
            global_trait_entries[0] = -1;
            num_traits = 1;
            if ((option == 1 || option == 3) && HasMarkers) {
                clean_reordered_markers(Top->LocusTop, *analysis);
            }
        } else if (*analysis == TO_ALLELE_FREQ || *analysis == TO_LIABLE_FREQ) {
            /* special case because this option counts alleles for
               affection status loci */
            int *tr = CALLOC((size_t) num_traits+1, int);
            for (i = 0; i < num_traits; i++) {
                tr[i] = global_trait_entries[i];
            }
            tr[num_traits] = -1;
            free(global_trait_entries);
            global_trait_entries = tr;
            num_traits++;
        }
        if ((option == 1 || option == 3) && HasMarkers) {
            clean_reordered_markers(Top->LocusTop, *analysis);
        }
        Mega2Status = LOCI_REORDERED;
        /*    if (*analysis == TO_PREST) { */
        return Top;
    }

    /* If Y chromosome was selected, set all positions to 0, now that
       we have reordered markers */
    reset_ychrom_positions(Top);

    Mega2Status = LOCI_REORDERED;

    /* If not one of the options that do not require trait selection */
    if ((option == 1 || option == 3) && main_chromocnt > 0) {
        select_trait_loci(Top, *analysis);
    }

    /* At this time we should have the following:
       1) num_reordered = total number of marker loci selected
       2) reordered_loci = list of loci reordered by map and chromosome
       3) num_traits = number of traits selected,
       4) global_trait_entries = list and order of traits selected
       this may contain an item = -1 denoting the position of markers
    */
    /* define_affection labels  use the global_trait_entries */
    define_affection_labels(Top, *analysis);

    /* changed set_missing_quant_input and write_quant_stats to use global_trait_entries as well */
    set_missing_quant_input(Top, *analysis);
    if (write_quant_stats(Top, *analysis) != 0)
        set_missing_quant_output(Top, *analysis);

    if ((option == 1 || option == 3) && HasMarkers) {
        clean_reordered_markers(Top->LocusTop, *analysis);
    }
    /*   if (UnMappedLociCheck(Top))    { */
    /*     draw_line(); */
    /*     printf("WARNING: Some selected loci do not have genetic distance values.\n"); */
    /*     draw_line(); */
    /*   } */

    /*   printf("REorder locus option, %d\n", option); */
    /*   if (option == 2) { */
    /*     get_trait_positions(Top->LocusTop, *analysis); */
    /*   } */
    return Top;
}

void Free_reorder_loci(void) {
/*  might be alloc'ed in linkage.c */
    free(ChrLoci);
}

/*===========================================================================*/

static int cmp(const void *a, const void *b)
{
    int ai = *(int *)a;
    int bi = *(int *)b;
    return (ai < bi) ? -1 : ((ai == bi) ? 0 : 1);
}

static int      ReOrderLociByChromosome(linkage_ped_top * Top,
					int *numchr,
					int get_or_set)
{

    int             i, num_markers=0, selectionmade,  mm;
    linkage_locus_top *LTop;
    int             chromo_cnt = 1, num_cols=3;
    int            *locus_list, *display_chrom;
    char cselect[6];
    char prchr[4];

    LTop = Top->LocusTop;
    chromo_cnt = NumChromo;
    if (get_or_set == 0) {
#if 1
        display_chrom=CALLOC((size_t) NumChromo, int);
        for (i=0; i < NumChromo; i++)
            display_chrom[i] = global_chromo_entries[i];
        qsort(display_chrom, NumChromo, sizeof (int), cmp);
#else
        curr_min_chrom=1e6;
        for (i=0; i < NumChromo; i++) {
            if (global_chromo_entries[i] < curr_min_chrom) {
                curr_min_chrom = global_chromo_entries[i];
            }
        }
        display_chrom[0] = curr_min_chrom;
        for (k=1; k < NumChromo; k++) {
            curr_min_chrom=1e6;
            for (i=0; i < NumChromo; i++) {
                if (global_chromo_entries[i] < curr_min_chrom &&
                    global_chromo_entries[i] > display_chrom[k-1]) {
                    curr_min_chrom = global_chromo_entries[i];
                }
            }
            display_chrom[k] = curr_min_chrom;
        }
#endif
        selectionmade=0;
        printf("  Available chromosomes:\n  ");
        for (mm = 0; mm < chromo_cnt; mm++) {

            if (Ignore_Unmapped(display_chrom[mm])) continue;

            if (display_chrom[mm] < 10) num_cols += 3;
            else num_cols += 4;
            if (num_cols > 69) {
                printf("\n  "); num_cols=3;
            }
            CHR_STR(display_chrom[mm], prchr);
            printf("%s, ", prchr);
        }
        printf("\n  Enter the new chromosome number ");
        if (AllowUnmapped && NumUnmapped > 0)
            printf("(U for unmapped markers) ");
        printf("> ");
        fcmap(stdin, "%s", cselect);
        selectionmade=-1;
/*     if (cselect[0] == 'U') { */
/*       selectionmade = UNKNOWN_CHROMO; */
/*       num_markers = NumUnmapped; */
/*     } */
/*     else { */
        selectionmade = STR_CHR(cselect); // returns -1 on error
        if (selectionmade == UNKNOWN_CHROMO) {
            if (AllowUnmapped) 
                num_markers = NumUnmapped;
            else {
                printf("Unknown Chromosome can not be selected.\n");
                num_markers = -1;
            }
        } else if (selectionmade >= 1 && selectionmade <= MaxChromo) {
            for (mm = LTop->PhenoCnt; mm < LTop->LocusCnt; mm++) {
                if (LTop->Marker[mm].chromosome == selectionmade) {
                    num_markers++;
                }
            }
        } else {
            printf("A chromosome must be a positive integer less than %d or one of the following:\n", lastautosome+4);
            printf("X, Y, XY, MT, U.\n");
        }

        if (num_markers == -1) {
        } else if (num_markers == 0) {
            draw_line();
            printf("No loci found on chromosome %s\n", cselect);
        } else {
            *numchr = selectionmade;
            printf("Chromosome number is now %s.\n", cselect);
        }
        free(display_chrom);
    } else {
        if (InputMode == INTERACTIVE_INPUTMODE) {
            Mega2BatchItems[/* 7 */ Chromosome_Single].value.option = *numchr;
            batchf(Chromosome_Single);
        }
        if (*numchr == UNKNOWN_CHROMO) {
            num_reordered = NumUnmapped;
            reordered_marker_loci = unmapped_markers;
        } else {
            num_markers = ReOrderMappedLoci_new(Top, 1, &locus_list, 1, numchr);
            linkage_set_locus_order_new(num_markers, locus_list);
            new_chromo_loci_counts(1, numchr);
            free(locus_list);
        }

        free(global_chromo_entries);
        global_chromo_entries = CALLOC((size_t) 1, int);
        global_chromo_entries[0] = *numchr;
        main_chromocnt = 1;
    }
    return 1;
}

static void reselect_all_loci(int *num_selected, int *selections)

{
    int i;
    for (i=0; i<*num_selected; i++) {
        selections[i]=0;
    }
    *num_selected=0;
    return;

}


static int verify_markers(int entrycount, int *number,
			  linkage_ped_top *Top,
			  analysis_type analysis)


{
    int naff, *traits, i, sw2check, reset=0;

    /*
     * If your analysis option has some special
     * requirements besides 1-5 on locus selection, then add a special
     * checking routine to this function. The variable reset is a flag which
     * is to be set by this checking routine if the selection is not
     * satisfactory.
     */

    if (entrycount == 0) {
        warnvf("No loci selected.\n");
        return 1;
    }
    traits = CALLOC((size_t) entrycount, int);

    naff=0;
    sw2check=1;
    for (i=0; i < entrycount; i++) {
        if (number[i] < 0 || number[i] >= Top->LocusTop->LocusCnt) {
            warnvf("Locus %d out of bounds!\n", number[i]+1);
            reset=1;
            break;
        }
        /* check for quant loci */
        else if (!quant_allowed(Top->LocusTop->Locus[number[i]], analysis)) {
            reset=1;
            break;
        }
        /* check for non-quantitative loci in QUANT_SUMMARY */
        else if (Top->LocusTop->Locus[number[i]].Type != QUANT &&
                 analysis == QUANT_SUMMARY) {
            errorvf("Locus %s is not a quantitative trait locus. Only quantitative loci allowed for %s.\n",
                    Top->LocusTop->Locus[number[i]].Name,
		    ProgName);
            EXIT(EARLY_TERMINATION);
        }

        /* check if a marker has a usable position */
        if ((Top->LocusTop->Locus[number[i]].Type == XLINKED ||
             Top->LocusTop->Locus[number[i]].Type == YLINKED ||
             Top->LocusTop->Locus[number[i]].Type == NUMBERED) &&
            ((genetic_distance_sex_type_map == SEX_AVERAGED_GDMT &&
	      Top->LocusTop->Marker[number[i]].pos_avg < 0.0) ||
             ((genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
	       genetic_distance_sex_type_map == FEMALE_GDMT) &&
	      Top->LocusTop->Marker[number[i]].pos_female < 0.0))
	    ) {
	    // it's a marker with a sex-averaged map and no position, or a sex-specific map and no female position...
            if (!ALLOW_NO_MAP(analysis)) {
                warnvf("Marker number %d does not have any genetic map positions.\n",
		       number[i]+1);
                reset=1;
                break;
            }
        }

        if (Top->LocusTop->Locus[number[i]].Type == AFFECTION ||
            Top->LocusTop->Locus[number[i]].Type == QUANT) {
            if (analysis == TO_SAGE && i > naff) {
                /* was there a marker locus in the middle? */
                warnvf("SAGE requires all trait loci to be first in the order.\n");
                reset=1;
                break;
            }
            /* error if first aff not found at pos 0
               or more than one aff */
            if ((SIMWALK2(analysis)) && !naff && (i > 0)) {
                warnvf("Simwalk2 options require trait locus to be first in the order.\n");
                naff=1;
                reset=1; sw2check=0;
                break;
            }
            traits[naff] = number[i]+1;   naff++;
        }
    }

    reset +=
        ((analysis == TO_GeneHunter)?
         1 - check_gh_quants(Top, number, entrycount, -1) : 0) +
        ((analysis == TO_GHMLB)?
         1 - check_ghmlb_affs(Top, traits, naff, -1): 0) +
        ((analysis == TO_SAGE)?
         1 - all_same_type(Top, traits, naff, -1): 0) +
        ((SIMWALK2(analysis))?
         1 - check_simwalk2_affs(analysis, naff, traits, -1, sw2check) : 0) +
        ((naff == 0)?
         1 - no_trait_allowed(analysis, 1) : 0);

    if (reset == 0) {
        if ((naff > 1 && allow_trait_combination(analysis)) || naff <= 1) {
            LoopOverTrait = 0;
        } else {
            LoopOverTrait = 1;
        }
    }

    free(traits);
    traits = NULL;

    return(reset);

}

static int *check_for_duplicates(int *entrycount, int *number)
{
    int *number1, i, j;
    int DuplicateCount = 0;
    int *DuplicateNumber;


    DuplicateNumber = CALLOC((size_t) *entrycount, int);
    // CALLOC makes this unnecessary
    //for (i = 0; i < *entrycount; i++)
    //  DuplicateNumber[i]=0;

    // Yes, this is an N^2 algorithm and should be replaced by a hash table
    for (i = 0; i < *entrycount-1; i++)  {
        for (j = i + 1; j < *entrycount; j++)   {
            if (number[i] == number[j]) {
                DuplicateNumber[j] = 1;
                DuplicateCount++;
            }
        }
    }
    if (DuplicateCount > 0)
        warnvf("duplicate locus numbers, removing duplicates.\n");

    number1 = CALLOC((size_t) (*entrycount - DuplicateCount + 1), int);
    j=0; /* only increment if not a duplicate */
    for (i = 0; i < *entrycount; i++) {
        if (DuplicateNumber[i] != 1) {
            number1[j]=number[i]; j++;
        }
    }
    number1[j] = -99;
    free(DuplicateNumber);
    *entrycount = *entrycount - DuplicateCount;
    return(number1);
}

int display_selections(linkage_ped_top *Top, int *entries,
		       int num_entries, char *asterisks, char *scanned_str,
		       int traits_only, int map_num)

{
    int MaxInScreen, LowerLim, UpperLim;
    int k, entry_count, *number1;
    int has_err=0, inp, scanned_text=0;
    int locnamewidth, double_col=0;
    char namestr[10], posstr[10];
    char user_input[100];
    double pos;

    linkage_locus_top *LocusTop = Top->LocusTop;

    if (entries == NULL) {
        /* ignore num_entries, just display all loci */
        number1 = CALLOC((size_t) LocusTop->LocusCnt, int);
        entry_count = LocusTop->LocusCnt;
        for (k=0; k < entry_count; k++) number1[k]=k;
    } else {
        number1 = CALLOC((size_t) num_entries, int);
        entry_count = num_entries;
        for (k=0; k < entry_count; k++) number1[k]=entries[k];
    }

    locus_name_width(LocusTop, &(locnamewidth));

    /* now display the loci from Top */
    MaxInScreen = MaxLinesPerScreen;
    if (ErrorSimOpt == 0) {
        has_err = 0;
    } else {
        for (k = 0; k < entry_count; k++)   {
            if (LocusTop->Marker[number1[k]].error_prob > 0.0) {
                has_err=1;
                break;
            }
        }
    }
    if (locnamewidth <= 15 && has_err == 0 &&
       entry_count > MaxLinesPerScreen) {
        double_col=1;
    } else if (locnamewidth <= 30 && traits_only==1 &&
            entry_count > MaxLinesPerScreen) {
        double_col = 1;
    }

    sprintf(namestr, "%%%ds ", locnamewidth);

    if (double_col) {
        printf("       ");  printf(namestr, "Name");
        if (traits_only == 1) {
            printf("  Type     |");
        } else {
            printf("  Position  Type     |");
        }
        printf("       ");  printf(namestr, "Name");
        if (traits_only == 1) {
            printf("  Type     |");
        } else {
            printf("  Position  Type\n");
        }
    } else {
        printf(namestr, "         Name");
        if (traits_only == 1) {
            printf("Type\n");
        } else {
            if (has_err ==0) {
                printf("Position   Type\n");
            } else {
                printf("Position    Type   Error\n");
            }
        }
    }
    LowerLim=1; UpperLim=entry_count;
    inp = 0;

    if (double_col) {
        for (k = 0; k < entry_count; k+=2)   {
            if (LowerLim == 1 && UpperLim < MaxLinesPerScreen)
                MaxInScreen = UpperLim;
            if (LowerLim > MaxInScreen) {
                inp=press_return();
                if (inp == 1) {
                    break;
                } else if (inp == -1 && scanned_str != NULL) {
                    fflush(stdout);
                    (void)fgets(user_input, 99, stdin);
                    strcat(scanned_str, user_input);
                    scanned_text=1;
                }
                LowerLim=1;
            }

            if (asterisks != NULL) {
                printf("%c", asterisks[k]);
            } else {
                printf(" ");
            }
            printf("%4d) ", k+1);
            printf(namestr, LocusTop->Locus[number1[k]].Name);
            if (traits_only == 1) {
                printf(" %-9s ", loc_type_name(LocusTop->Locus[number1[k]].Type));
            } else {
                if (LocusTop->Locus[number1[k]].Class == TRAIT) {
                    pos = UNKNOWN_POSITION;
                } else if (Top->EXLTop == NULL) {
		  // What about SEX_SPECIFIC???
                    pos = LocusTop->Marker[number1[k]].pos_avg;
                } else {
                    if (LocusTop->Marker[number1[k]].chromosome == SEX_CHROMOSOME) {
                        pos = ((Top->EXLTop->EXLocus[number1[k]].pos_female[map_num] >= 0.0)?
                               Top->EXLTop->EXLocus[number1[k]].pos_female[map_num] :
                               Top->EXLTop->EXLocus[number1[k]].positions[map_num]);
                    } else if (LocusTop->Marker[number1[k]].chromosome == MALE_CHROMOSOME) {
                        pos = UNKNOWN_POSITION;
                    } else {
                        pos = Top->EXLTop->EXLocus[number1[k]].positions[map_num];
                    }
                }

                if (pos >= 0.0) {
                    sprintf(posstr, "%8.5f", pos);
                } else {
                    sprintf(posstr, "   -    ");
                }
                printf(" %9s %-9s ",
                       posstr,
                       loc_type_name(LocusTop->Locus[number1[k]].Type));
            }
            if ((k+1) < entry_count) {
                printf("|");
                if (asterisks != NULL) {
                    printf("%c", asterisks[k+1]);
                } else {
                    printf(" ");
                }
                printf("%4d) ", k+2);
                printf(namestr, LocusTop->Locus[number1[k+1]].Name);
                if (traits_only == 1) {
                    printf(" %-9s ", loc_type_name(LocusTop->Locus[number1[k+1]].Type));
                } else {
                    if (LocusTop->Locus[number1[k+1]].Class == TRAIT) {
                        pos = UNKNOWN_POSITION;
                    } else if	 (Top->EXLTop == NULL) {
                        pos = LocusTop->Marker[number1[k+1]].pos_avg;
                    } else {
                        if (LocusTop->Marker[number1[k+1]].chromosome == SEX_CHROMOSOME) {
                            pos = ((Top->EXLTop->EXLocus[number1[k+1]].pos_female[map_num] >= 0.0)?
                                   Top->EXLTop->EXLocus[number1[k+1]].pos_female[map_num] :
                                   Top->EXLTop->EXLocus[number1[k+1]].positions[map_num]);
                        } else if (LocusTop->Marker[number1[k+1]].chromosome == MALE_CHROMOSOME) {
                            pos = UNKNOWN_POSITION;
                        } else {
                            pos = Top->EXLTop->EXLocus[number1[k+1]].positions[map_num];
                        }
                    }

                    if (pos >= 0.0) {
                        sprintf(posstr, "%8.5f", pos);
                    } else {
                        sprintf(posstr, "   -    ");
                    }
                    printf(" %9s %-9s \n",
                           posstr,
                           loc_type_name(LocusTop->Locus[number1[k+1]].Type));
                }
            }
            LowerLim++; UpperLim--;
        }
    } else {
        for (k = 0; k < entry_count; k++)   {
            if (LowerLim == 1 && UpperLim < MaxLinesPerScreen)
                MaxInScreen = UpperLim;
            if (LowerLim > MaxInScreen) {
                inp=press_return();
                if (inp == 1) {
                    break;
                } else if (inp == -1 && scanned_str != NULL) {
                    fflush(stdout);
                    (void)fgets(user_input, 99, stdin);
                    strcat(scanned_str, user_input);
                    scanned_text=1;
                }
                LowerLim=1;
            }

            if (asterisks != NULL) {
                printf("%c ", asterisks[k]);
            }
            printf("%-4d ", k+1);
            printf(namestr, LocusTop->Locus[number1[k]].Name);
            if (traits_only == 1) {
                printf("%-9s ", loc_type_name(LocusTop->Locus[number1[k]].Type));
            } else {
                if (LocusTop->Locus[number1[k]].Class == TRAIT) {
                    pos = UNKNOWN_POSITION;
                } else if (Top->EXLTop == NULL) {
                    pos = LocusTop->Marker[number1[k]].pos_avg;
                } else {
                    /* not trait and has a extended locus top */
                    if (LocusTop->Marker[number1[k]].chromosome == SEX_CHROMOSOME) {
                        pos = ((Top->EXLTop->EXLocus[number1[k]].pos_female[map_num] >= 0.0)?
                               Top->EXLTop->EXLocus[number1[k]].pos_female[map_num] :
                               Top->EXLTop->EXLocus[number1[k]].positions[map_num]);
                    } else if (LocusTop->Marker[number1[k]].chromosome == MALE_CHROMOSOME) {
                        pos = UNKNOWN_POSITION;
                    } else {
                        pos = Top->EXLTop->EXLocus[number1[k]].positions[map_num];
                    }
                }
                if (pos >= 0.0) {
                    sprintf(posstr, "%8.5f", pos);
                } else {
                    sprintf(posstr, "   -    ");
                }
                printf("%s  %-9s ",
                       posstr,
                       loc_type_name(LocusTop->Locus[number1[k]].Type));
                if (LocusTop->Locus[number1[k]].Type == NUMBERED &&
                   ErrorSimOpt == 1 &&
                   LocusTop->Marker[number1[k]].error_prob >= 0.0) {
                    printf("%4.3f", LocusTop->Marker[number1[k]].error_prob);
                }
            }
            printf("\n");
            LowerLim++; UpperLim--;
        }
    }
    printf("\n");
    free(number1);
    return scanned_text;
}

/* here the user interface will permit the user to enter a string of
   position numbers (position of loci in the LOCUS DATA FILE file).
   The order user enters the numbers is the new loci number...
   eg.  '1 52 7 8 22' means loci #1 is first, #52 is second, and so on.
   Assume that ReOrderLociByChromosomeNumber has been called at least once,
   so that LocusOrder and NumMarkers have been set.
*/

#define BAD_FORM -1
#define CONTINUE_INPUT 1
#define FINISH_INPUT 2

int parse_string(char *entrydummy, int *number,
		 int *num_entries, int num_total_loci, int require_e)

{
    int rflag=0, nflag=1, num=0;
    int i, entrycount = *num_entries;
    char *entry;
    int lhs=0;

    entry = &(entrydummy[0]);

    while(*entry != '\0') {
        if ((!isdigit((int) *entry)) && (*entry != ' ') &&
           (*entry != '\n') && (*entry != '-') &&
           (*entry != 'e') && (*entry != 'E') &&
           (*entry != 'm') && (*entry != 'M')) {
            printf("ERROR: The character '%c' is bad (not numeric, '-' or whitespace)!\n",*entry);
            return BAD_FORM;
        }

        if (isdigit((int)*entry)) {
            nflag=0; /* reading a number */
            while(isdigit((int)*entry)) {
                num=10*num+ ((int)*entry - 48);
                entry++;
            }
        }

        if (*entry == 'e' || *entry == 'E' || *entry == '\n' ||
            *entry == 'm' || *entry == 'M' ||
            *entry == ' ' || *entry == '-' || *entry == '\0') {
            if (*entry == '-') {
                if (nflag==0) {
                    /* did we read the beginning of a range? */
                    lhs=num;
                    if (number != NULL) {
                        number[entrycount] = num - 1;
                    }
                    entrycount++;
                }
                if (entrycount==0) {
                    /* no numbers read  */
                    printf("ERROR: '-' encountered without a preceding number.\n");
                    return BAD_FORM;
                } else if (rflag == 1) {
                    printf("ERROR: Encountered invalid range %s       should be only two numbers I-J!\n",
                           entrydummy);
                    return BAD_FORM;
                } else {
                    rflag=1;
                }
                num=0;
                nflag=1;
                entry++;
            } else {
                /* if entry is not '-' */
                if (nflag==0)  {
                    /* finished reading a number */
/* This change in response to Bug 240 breaks covariate selection:
   if (num > num_total_loci) {
   printf("ERROR: Invalid selection %s      Number exceeds maximum value allowed (%d)!\n",
   entrydummy, num_total_loci);
   return BAD_FORM;
   }
*/
                    if (rflag==1) {
                        /* right side of a range */
                        /*	    for (i=number[entrycount-1] + 1; i < num; i++) { */
                        for (i=lhs+1; i <= num; i++) {
                            if (number != NULL) {
                                number[entrycount] = i-1;
                            }
                            entrycount++;
                        }
                        rflag=0;
                    } else {
                        /* just a single locus */
                        if (number != NULL) {
                            number[entrycount]  = num - 1;
                        }
                        entrycount++;
                    }
                }
                nflag=1;
                num=0;
                if (*entry == ' ') {
                    while(*entry == ' ') {
                        entry++;
                    }
                } else if (*entry == 'M' || *entry == 'M') {
                    number[entrycount] = -1;
                    entry++;
                } else if (*entry == 'e' || *entry == 'E') {
                    *num_entries = entrycount;
                    return FINISH_INPUT;
                } else if (*entry == '\n' || *entry == '\0') {
                    *num_entries = entrycount;
                    if (require_e) {
                        return CONTINUE_INPUT;
                    } else {
                        return FINISH_INPUT;
                    }
                }
            }
        }
    }
    /* This will never happen, as the string has to end with a newline,
       unless the string is very long */
    if (require_e) {
        return CONTINUE_INPUT;
    } else {
        *num_entries = entrycount;
        return FINISH_INPUT;
    }
}


/* #ifdef __MWERKS__ */
/* #pragma global_optimizer on */
/* #pragma optimization_level 2 */
/* #endif */

/* this works on a list of indexes that are locus numbers of loci
   that are present on the selected chromosome. Traits are also included
   in this list, unless 'm' is specified.
*/

//
// numchr == -1 select only the quantitative traits from the global_trait_entries.

static void invalid_value_in_list(int batch_item, char *list_item)

{
    errorvf("Invalid value-list item %s for keyword %s.\n",
            list_item, C(Mega2BatchItems[batch_item].keyword));
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static int       ReOrderLociByPositionNumber(linkage_ped_top *Top,
					     int numchr, int *num_loci,
					     int **selected_loci,
					     analysis_type analysis,
					     int get_or_set,
					     int map_num)
{
    char            *MarkSelection, *entry;
    char            entrycopy[FILENAME_LENGTH], entrydummy[FILENAME_LENGTH];
    int             *number1=NULL, *number, num_chr_loci=0, *chromo_loci;
    int             i, j, entrycount, reset=0;
    int             HasBeenReOrdered = 0, num;
    int             press_retval;
    int             *old_number, old_entrycount;
    int             *display_loci = NULL;
    tmpnum          *tmpNumber;
    int (*comp_index)(const void *, const void *) = compare_index;
    linkage_locus_top *LTop = Top->LocusTop;

    if (numchr < 1) {
        if (numchr == -1) {
            for (i = 0; i < num_traits; i++) {
                if (LTop->Locus[global_trait_entries[i]].Type == QUANT) {
                    num_chr_loci++;
                }
            }
        } else if (numchr == 0) {
            num_chr_loci = num_traits;
        }

        chromo_loci = CALLOC((size_t) num_chr_loci + num_traits, int);
        number = CALLOC((size_t) num_chr_loci + num_traits, int);
        if (numchr == -1) {
            for (i = 0, j = 0; i < num_traits; i++) {
                if (LTop->Locus[global_trait_entries[i]].Type == QUANT) {
                    chromo_loci[j++] = global_trait_entries[i];
                }
            }
        } else {
            for (i = 0; i < num_traits; i++) {
                chromo_loci[i] = global_trait_entries[i];
            }
        }
    } else {
        if (numchr == UNKNOWN_CHROMO) {
            num_chr_loci = NumUnmapped;
            tmpNumber = CALLOC((size_t) num_chr_loci, tmpnum);
            for (i = 0; i < num_chr_loci; i++) {
                tmpNumber[i].index = UNKNOWN_POSITION;
                tmpNumber[i].number = unmapped_markers[i];
            }
        } else {
            // Find the entry for this chromosome...
            for (i = 0; i < NumChromo; i++) {
                if (numchr == global_chromo_entries[i]) {
                    break;
                }
            }
	    // The number of loci at a particular chromosome....
	    // chromo_loci_count is a cumulative count...
            num_chr_loci = chromo_loci_count[i] - ((i > 0) ? chromo_loci_count[i - 1] : 0);
            tmpNumber = CALLOC((size_t) num_chr_loci, tmpnum);
            // Reorder the loci by position if on a valid chromosome.
	    // Remember to use the appropriate position based on the input map selection...
            for (i = LTop->PhenoCnt, j = 0; i < LTop->LocusCnt; i++) {
                if (LTop->Marker[i].chromosome == numchr) {
		    if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                        tmpNumber[j].index = LTop->Marker[i].pos_avg;
		    } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
			       genetic_distance_sex_type_map == FEMALE_GDMT) {
                        tmpNumber[j].index = LTop->Marker[i].pos_female;
		    } else {
		      errorvf("You must have specified a genetic distance sex map type.\n");
		      EXIT(DATA_TYPE_ERROR);
		    }
                    if (base_pair_position_index >= 0)
                        tmpNumber[j].phys = Top->EXLTop->EXLocus[i].positions[base_pair_position_index];
                    else
                        tmpNumber[j].phys = i;
                    tmpNumber[j++].number = i;
                }
            }
            qsort((void *) tmpNumber, (size_t) num_chr_loci, sizeof(tmpnum), comp_index);
        }

        number = CALLOC((size_t) num_chr_loci+num_traits, int);

	// an array containing the traits loci, followed by the 
        chromo_loci = CALLOC((size_t) num_chr_loci+num_traits, int);
        for (i = 0; i < num_traits; i++) {
            chromo_loci[i] = global_trait_entries[i];
        }
        for (i = 0, j = num_traits; i < num_chr_loci; i++) {
            chromo_loci[j++] = tmpNumber[i].number;
        }
        num_chr_loci += num_traits;

        free(tmpNumber);
    }

    if (batchREORDER) {
        if (numchr == -1) {
            if (Mega2BatchItems[/* 15 */ Traits_Combine].items_read) {
                strcpy(entrydummy, Mega2BatchItems[/* 15 */ Traits_Combine].value.name);
            } else if (Mega2BatchItems[/* 12 */ Trait_Single].items_read) {
                sprintf(entrydummy, "%d", Mega2BatchItems[/* 12 */ Trait_Single].value.option);
            }
        } else {
            if (Mega2BatchItems[/* 11 */ Loci_Selected].items_read) {
                strcpy(entrydummy, Mega2BatchItems[/* 11 */ Loci_Selected].value.name);
            } else if (Mega2BatchItems[/* 14 */ Traits_Loop_Over].items_read) {
                strcpy(entrydummy, Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name);
           }
        }
        entry = &(entrydummy[0]);
        if (*entry == ' ') {
            while (isspace((int) *entry)) {
                entry++;
            }
        }
        if (*entry == 'M' || *entry == 'm') {
            if (numchr < 1) {
                errorf("No marker loci to select.");
                invalid_value_in_list(11, entry);
            }
            entrycount = 0;
            for (j = 0; j < num_chr_loci; j++) {
                if (LTop->Locus[chromo_loci[j]].Type == NUMBERED) {
                    number[entrycount] = chromo_loci[j];
                    entrycount++;
                }
            }
        } else {
            entrycount = 0;
            num = parse_string(entry, number, &entrycount, num_chr_loci, 0);
            if (num == BAD_FORM) {
                invalid_value_in_list(11, &(entry[0]));
            }
            if (numchr >= 0) {
                for (i = 0; i < entrycount; i++) {
                    if (number[i] < 0 || number[i] >= num_chr_loci) {
                        printf("Number %d out of bounds\n", number[i]);
                        invalid_value_in_list(11, &(entry[0]));
                    }
                }
            } else {
                /* Just check against entire list of numbers */
                for (i = 0; i < entrycount; i++) {
                    if (number[i] < 0 || number[i] >= Top->LocusTop->LocusCnt) {
                        printf("Number %d out of bounds\n", number[i]);
                        invalid_value_in_list(11, &(entry[0]));
                    }
                }
            }

        }

        /* Usual locus selection checks */
        reset = verify_markers(entrycount, number, Top, analysis);
        if (reset > 0) {
            invalid_value_field(11); /* Terminate mega2 */
        }
        /* now check for duplicates... */
        number1 = check_for_duplicates(&entrycount, number);

        *selected_loci = CALLOC((size_t) entrycount, int);
        *num_loci = entrycount;

        for (i = 0; i < entrycount; i++) {
            if (numchr >= 0) {
                (*selected_loci)[i] = chromo_loci[number1[i]]; /* loci are reorered */
            } else {
                (*selected_loci)[i] = number1[i];
            }
        }
    } else if (get_or_set == 0) {
        MarkSelection = CALLOC((size_t) num_chr_loci, char);
        old_number = CALLOC((size_t) num_chr_loci, int);
        old_entrycount = 0;

        if (numchr == -1) {
            printf("       Select quantitative trait loci        \n");
        } else if (numchr == 0) {
            printf("       Select affection status loci\n");
        } else {
            printf("       Select marker loci on chromosome %d         \n",
                   numchr);
        }
        printf(
            "Select a set of loci by entering their numbers in the desired order.\n");
        printf("Separate each number by a space.\n");
        printf(
            "  Enter 'v' to view the original list of loci and their numbers.\n");
        if (numchr > 0) {
            printf(
                "  Enter 'm' to select all marker loci on current chromosome.\n");
        }
        printf("  Enter 'o' to list the ordered loci you have selected.\n");
        printf("  Enter 'e' to terminate the selection process.\n");

        strcpy(entrydummy, "v");
        strcpy(entrycopy, "");
        press_retval = 0;

        while (1) {
            entrycount = 0;
            if (press_retval == 0) {
                printf("Enter 'v', 'm', 'o', 'e' or a set of locus numbers > ");
                fflush(stdout);
                (void)fgets(entrydummy, FILENAME_LENGTH - 1, stdin);
                newline;
            }
            entry = &(entrydummy[0]);
            while (isspace((int) *entry)) {
                entry++;
            }

            switch (*entry) {
            case 'm':
            case 'M':
                if (numchr == -1) {
                    printf("Marker loci not handled by analysis type.\n");
                    break;
                } else if (numchr == 0) {
                    printf("No markers to select.\n");
                    break;
                }
            press_retval = 0;
            strcpy(entrycopy, "M");
            /* Select all codominant markers */
            printf("Selecting all marker loci...\n");
            old_entrycount = 0;
            for (j = 0; j < num_chr_loci; j++) {
                if (LTop->Locus[chromo_loci[j]].Type == NUMBERED) {
/*old						old_number[old_entrycount] = chromo_loci[j];
  new*/						old_number[old_entrycount] = j;
                    MarkSelection[j] = '*';
                    old_entrycount++;
                }
            }
            break;

            case 'o':
            case 'O':
                if (old_entrycount > 0) {
                    if (display_loci != NULL) {
                        free(display_loci);
                    }
                    display_loci = CALLOC((size_t) old_entrycount, int);
                    for (i = 0; i < old_entrycount; i++) {
                        display_loci[i] = chromo_loci[old_number[i]];
                    }
                    /* here  are current selections (in NEW order!) */
                    printf("Selected loci (in new order)\n");
                    press_retval = display_selections(Top, display_loci,
                                                      old_entrycount, NULL, NULL, 0, map_num);
                } else {
                    printf("No loci have been selected yet.\n");
                }
            break;

            case 'v':
            case 'V':
                for (i = 0; i < num_chr_loci; i++)
                    MarkSelection[i] = ' ';
            if (old_entrycount > 0) {
                for (i = 0; i < old_entrycount; i++)
                    MarkSelection[old_number[i]] = '*';
            }

            strcpy(entrydummy, "");
            printf("Original list of loci\n");
            press_retval = display_selections(Top, chromo_loci,
                                              num_chr_loci, MarkSelection, entrydummy, 0, map_num);
            break;

            default:
                press_retval = 0;
                /* get the next number */
                strcat(entrycopy, entry);
                chomp(entrycopy);
                num = parse_string(entry, number, &entrycount, num_chr_loci, 1);

                if (num == BAD_FORM) {
                    printf("            Re-select  all loci.\n");
                    reselect_all_loci(&entrycount, number);
                    reselect_all_loci(&old_entrycount, old_number);
                    strcpy(entrydummy, "v");
                    strcpy(entrycopy, "");
                    entry = &(entrydummy[0]);
                } else if (num == CONTINUE_INPUT || num == FINISH_INPUT) {
                    if (entrycount > 0) {
                        /* now check for duplicates... */
                        number1 = check_for_duplicates(&entrycount, number);
                    }
                    if (num == CONTINUE_INPUT) {
                        for (i = 0; i < entrycount; i++) {
                            old_number[old_entrycount + i] = number1[i];
                        }
                        old_entrycount += entrycount;
                    } else {

                        /* allocate the selected_loci and copy these numbers over */
                        /* first copy over the old loci, then the new ones */

                        if (*selected_loci != NULL) {
                            free(*selected_loci);
                        }
                        *selected_loci
                            = CALLOC((size_t) old_entrycount + entrycount, int);
                        *num_loci = old_entrycount + entrycount;

                        for (i = 0; i < old_entrycount; i++) {
                            (*selected_loci)[i] = chromo_loci[old_number[i]];
                        }

                        reset = 0;
                        for (j = 0; j < entrycount; j++) {
                            if (number1[j] < 0 || number1[j] >= num_chr_loci) {
                                printf("Number %d out of bounds\n", number1[j]
                                       + 1);
                                reset = 1;
                                break;
                            }
                            (*selected_loci)[i + j] = chromo_loci[number1[j]];
                        }

                        if (!reset) {
                            reset = verify_markers(entrycount + old_entrycount,
                                                   *selected_loci, Top, analysis);
                        }
                        if (reset > 0) {
                            strcpy(entrydummy, "v");
                            entry = &(entrydummy[0]);
                            strcpy(entrycopy, "");
                            printf("            Re-select  all loci.\n");
                            reselect_all_loci(&entrycount, number);
                            reselect_all_loci(&old_entrycount, old_number);
                            break;
                        }

                        /* loci are reordered */
                        /* display  selections */
                        print_only("\nSelected loci (in final order)");
                        number1
                            = check_for_duplicates(num_loci, *selected_loci);
                        for (i = 0; i < *num_loci; i++) {
                            (*selected_loci)[i] = number1[i];
                        }
                        log_marker_selections(Top, *selected_loci, *num_loci,
                                              NULL, (void (*)(const char *)) print_only, ErrorSimOpt);
                        /* sleep(2); */
                        /*  Write the selections into a batch file */
                        if (InputMode == INTERACTIVE_INPUTMODE) {
                            if (numchr == -1) {
                                strcpy(Mega2BatchItems[/* 15 */ Traits_Combine].value.name,
                                       entrycopy);
                            } else {
                                /* 	      printf("entrycopy %s\n", entrycopy); */
                                Mega2BatchItems[/* 7 */ Chromosome_Single].value.option = numchr;
                                /* 	      batchf(Chromosome_Single); */
                                strcpy(Mega2BatchItems[/* 11 */ Loci_Selected].value.name,
                                       entrycopy);
                                /*	      batchf(Loci_Selected); */
                            }
                        }

                        /* we are done */
                        free(MarkSelection);
                        free(number);
                        free(number1);
                        free(old_number);
                        if (display_loci != NULL) {
                            free(display_loci);
                        }
                        return HasBeenReOrdered = 1;
                    }
                }
                break;
            }
        }
    }

    if (batchREORDER || get_or_set) {
        free(global_chromo_entries);
        global_chromo_entries = CALLOC((size_t) 1, int);
        main_chromocnt = 1;
        global_chromo_entries[0] = numchr;
        linkage_set_locus_order_new(*num_loci, *selected_loci);
        mssgf("\nSelected loci (in order)");
        log_marker_selections(Top, *selected_loci, *num_loci, NULL, (void (*)(const char *)) mssgf,
                              ErrorSimOpt);
        if (InputMode == INTERACTIVE_INPUTMODE) {
            if (numchr <= 0) {
                batchf(Traits_Combine);
            } else {
                batchf(Chromosome_Single);
                batchf(Loci_Selected);
            }
        }
        NumChrLoci = *num_loci;
        ChrLoci = *selected_loci;
        /* reset the global_trait_entries data structure */
        for (i = 0, j = 0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Class == TRAIT) {
                j++;
            }
        }
        num_traits = j;
	global_trait_entries
	  = REALLOC(global_trait_entries, (size_t) num_traits + ((LoopOverTrait == 0) ? 2 : 1), int);

        HasMarkers = 0; // do any of the loci have markers?
        for (i = 0, j = 0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Class == TRAIT) {
                // gather the indices of the traits into global_trait_entries...
                global_trait_entries[j++] = ChrLoci[i];
            } else if (LTop->Locus[ChrLoci[i]].Type == XLINKED
                       || LTop->Locus[ChrLoci[i]].Type == YLINKED
                       || LTop->Locus[ChrLoci[i]].Type == NUMBERED) {
	      // 
                LTop->Locus[ChrLoci[i]].Type = NUMBERED;
                HasMarkers = 1;
            }
        }
        if (LoopOverTrait == 0) {
            num_traits++;
            global_trait_entries[j] = -1;
            global_trait_entries[j + 1] = -99;
        } else {
            global_trait_entries[j] = -99;
        }
        /* sleep(2); */
        Mega2Status = TRAIT_SELECTED;
        get_trait_list(LTop, 0);
    }

    return 1;
}
/* After selecting chromosomes, do the cumulative counts again, so that it
   is easy to subset loci on each chromosome during output */

static void new_chromo_loci_counts(int num_chromo_selected, int *chromosome_list)
{
    int i,j, locus_count = 0;
    int *new_counts = CALLOC((size_t) num_chromo_selected, int);

    for (j=0; j < num_chromo_selected; j++) {
        for (i=0; i < NumChromo; i++) {
            if (chromosome_list[j] == global_chromo_entries[i]) {
                if (i==0) {
                    locus_count = chromo_loci_count[i];
                } else {
                    locus_count =  chromo_loci_count[i] - chromo_loci_count[i-1];
                }
                break;
            }
        }
        if (j == 0) {
            new_counts[j] = locus_count;
        } else {
            new_counts[j] = new_counts[j-1] + locus_count;
        }
    }

    free(chromo_loci_count);
    chromo_loci_count = new_counts;
}

/*-----------------------------------------------------------*/
/* local routine to verify that the chromosomes specified
   by user actually exist in the data read in */

static int verify_chromosomes(int *chromosome_list, int entry_index)

{
    int i,j;
    int found;
    char intstr[10];
    for (i = 0; i < entry_index; i++) {
        if (chromosome_list[i] == UNKNOWN_CHROMO) {
            if (NumUnmapped > 0 && AllowUnmapped) {
                continue;
            } else {
                if (InputMode == INTERACTIVE_INPUTMODE) {
                    printf("ERROR: No unmapped loci or output option does not allow unmapped loci.\n");
                    return 1;
                } else {
                    sprintf(intstr, "%d", chromosome_list[i]);
                    invalid_value_in_list(9, intstr);
                }
            }
        } else {
            found = 0;
            for (j = 0; j < NumChromo; j++) {
                if (chromosome_list[i] == global_chromo_entries[j]) {
                    found = 1; break;
                }
            }
            if (!found) {
                if (InputMode == INTERACTIVE_INPUTMODE) {
                    printf("ERROR: Selected  chromosome %d not present in data.\n",
                           chromosome_list[i]);
                    printf("       Please reselect all chromosomes.\n");
                    return(1);
                } else {
                    char prchr[4];
                    CHR_STR(chromosome_list[i], prchr);
//                  sprintf(intstr, "%d", chromosome_list[i]);
                    strcpy(intstr, prchr);
                    invalid_value_in_list(9, intstr);
                }
            }
        }
    }
    return 0;
}


/*=================================================================
  here the user may specify chromsomes to be included for analysis */

static int    chromo_and_loci_selection(linkage_ped_top *LPedTreeTop,
                                        int *num_chromo, int **chromosome_list,
                                        int get_or_set)
{

    int            do_again = 0, i, k;
    int            entry_index;
    char           reply[4], *chromo_string_p, chromo_string[201];
    char           entry_[MAX_CHROMO_NUM][3];
    char           choice[10];
    int            chr_digit, num_choice=0;
    int            *locus_order, loc_type = 1;
    int            *display_chrom;
    int            item, only_autosomes = 0;
    int            auto_it, all_it, sp_it, un_all_it;
    int            num_auto=NumChromo;
    char           prchr[4];

    /*  locus_order = CALLOC((size_t) LPedTreeTop->LocusTop->LocusCnt+1, int); */

    if (!get_or_set) {
        display_chrom=CALLOC((size_t) NumChromo, int);
#if 1
        for (i=0; i < NumChromo; i++) {
            if (LPedTreeTop->LocusTop->SexLinked == 2 &&
                (global_chromo_entries[i] == SEX_CHROMOSOME ||
                 global_chromo_entries[i] == MALE_CHROMOSOME ||
                 global_chromo_entries[i] == PSEUDO_X ||
                 global_chromo_entries[i] == MITO_CHROMOSOME ||
                 global_chromo_entries[i] == UNKNOWN_CHROMO )) {
                num_auto--;
                only_autosomes = 1;
            }
            display_chrom[i] = global_chromo_entries[i];
        }
        qsort(display_chrom, NumChromo, sizeof (int), cmp);
#else
        curr_min_chrom=1e6;
        for (i=0; i < NumChromo; i++) {
            if (LPedTreeTop->LocusTop->SexLinked == 2 &&
                (global_chromo_entries[i] == SEX_CHROMOSOME ||
                 global_chromo_entries[i] == MALE_CHROMOSOME ||
                 global_chromo_entries[i] == PSEUDO_X)) {
                num_auto--;
                only_autosomes = 1;
            }
            if (global_chromo_entries[i] < curr_min_chrom) {
                curr_min_chrom = global_chromo_entries[i];
            }
        }
        display_chrom[0] = curr_min_chrom;
        for (k=1; k < NumChromo; k++) {
            curr_min_chrom=1e6;
            for (i=0; i < NumChromo; i++) {
                if (global_chromo_entries[i] < curr_min_chrom &&
                    global_chromo_entries[i] > display_chrom[k-1]) {
                    curr_min_chrom = global_chromo_entries[i];
                }
            }
            if (curr_min_chrom != UNKNOWN_CHROMO || AllowUnmapped)
                display_chrom[k] = curr_min_chrom;
        }
#endif
        auto_it = all_it = sp_it = un_all_it = -1;
        item=1;
        if (only_autosomes) {
            printf("%d) Select all autosomes. \n", item); auto_it = item;
            item++;
        }
        printf("%d) Select all chromosomes. \n", item); all_it = item; item++;
        printf("%d) Select markers on specific chromosomes. \n", item); sp_it = item;
        if (NumUnmapped > 0 && AllowUnmapped) {
            item++;
            printf("%d) Select all marker loci (mapped and unmapped). \n", item); un_all_it=item;
        }
        printf("Select from 1-%d > ", item);
        fcmap(stdin, "%s", choice); newline;
        num_choice = -1;
        sscanf(choice, "%d", &num_choice);

        /* Selecting all markers in order on all chromosomes */

        if (num_choice == all_it) {
            printf("Selecting all chromosomes: ");
            *num_chromo = 0;
            /* chromosome_list needs to be allocated because of
               the subsequent code, even though it is simply a copy
               of global_chromo_entries */
            *chromosome_list = CALLOC((size_t) NumChromo, int);
            for (i=0; i< NumChromo; i++) {
                if (display_chrom[i] != UNKNOWN_CHROMO) {
                    (*num_chromo)++;
                    (*chromosome_list)[i] = display_chrom[i];
                }
            }
        } else if (num_choice == auto_it) {
            printf("Selecting all autosomes: ");
            *chromosome_list = CALLOC((size_t) num_auto, int);
            *num_chromo = num_auto;
            k=0;
            for (i=0; i< NumChromo; i++) {
                if (display_chrom[i] != SEX_CHROMOSOME &&
                    display_chrom[i] != MALE_CHROMOSOME &&
                    display_chrom[i] != PSEUDO_X &&
                    display_chrom[i] != MITO_CHROMOSOME &&
                    display_chrom[i] != UNKNOWN_CHROMO) {
                    (*chromosome_list)[k] = display_chrom[i];
                    k++;
                }
            }
        } else if (num_choice == un_all_it) {
            printf("Selecting all chromosomes and unmapped loci: ");
            *num_chromo = NumChromo;
            *chromosome_list = CALLOC((size_t) NumChromo, int);
            for (i=0; i< NumChromo; i++) {
                (*chromosome_list)[i] = display_chrom[i];
            }
        } else if (num_choice == sp_it) {
            do_again=1;
            while (do_again) {
                /* selection phase */
                printf("Select from the following chromosomes:\n");
                for (i = 0; i < NumChromo; i++) {
                    if (Ignore_Unmapped(display_chrom[i])) continue;
                    CHR_STR(display_chrom[i], prchr);
                    printf("%s ", prchr);
                    if (((i+1) % 20 == 0) && (i > 0))
                        printf("\n");
                }
                if (((i+1) % 20) != 0) newline;

                printf("Enter string of chromosome numbers, \n");
                printf("( press 'RETURN' or 'Enter' to terminate selection) > ");
                strcpy(chromo_string, "");
                fflush(stdout);
                /* fgets will store the newline */
                (void)fgets(chromo_string, 200, stdin); printf("\n");
                entry_index = 0;
                chromo_string_p = &chromo_string[0];
                while (*chromo_string_p != '\n' && *chromo_string_p != 'e') {
                    int chr_str = tolower(*chromo_string_p);
                    if (!isdigit((int)chr_str) && chr_str != 'e' &&
                        chr_str != ' ' && chr_str != '\n' &&
                        chr_str != 'x' && chr_str != 'y' &&
                        chr_str != 'm' && chr_str != 't' &&
                        chr_str != 'u') {
                        printf("Invalid characters in input.\n");
                        entry_index=0; do_again=1;
                        break;
                    }
                    while(*chromo_string_p == ' ') chromo_string_p++;
                    if (*chromo_string_p == '\0' || *chromo_string_p == 'e' ||
                       *chromo_string_p == '\n') break;
//                  if (isdigit((int)*chromo_string_p))
                    {
                        if (entry_index >= NumChromo) {
                            /* break out of selection phase */
                            printf("Exceeded number of chromosomes present in data, %d\n",
                                   NumChromo);
                            printf(
                                "Continue, ignoring the rest (enter 'n' to re-select)?(y/n) > ");
                            fcmap(stdin, "%s", reply); newline;
                            if (tolower((unsigned char)*reply) == 'n') entry_index = 0;
                            break;
                        }
                        /* if we are still under NumChromo, read in next chromosome */
                        chr_digit=0; strcpy(entry_[entry_index], " ");
                        while(!isspace((int)*chromo_string_p) && *chromo_string_p != '\n') {
                            entry_[entry_index][chr_digit]= *chromo_string_p;
                            chr_digit++;
                            entry_[entry_index][chr_digit]='\0';
                            chromo_string_p++;
                        }
                        entry_index++;
                    }
                }
                /* verification  phase */
                if (entry_index <= 0) {
                    printf("... No chromosomes specified, select again.\n ");
                    do_again=1;
                } else {
                    *chromosome_list = CALLOC((size_t) entry_index, int);
                    for (i=0; i < entry_index; i++)
                        (*chromosome_list)[i] = STR_CHR(entry_[i]);
                    do_again = verify_chromosomes(*chromosome_list, entry_index);
                }
            }
            *num_chromo = entry_index;
        } else {
            warn_unknown(choice);
        }
        for (i=0; i< *num_chromo; i++) {
            CHR_STR((*chromosome_list)[i], prchr);
            printf("%s ", prchr);
        }
        printf("\n");

        free(display_chrom);
    } else {
        if (batchREORDER) {
            /* Read the number of chromosomes */
            entry_index =  *num_chromo = Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].value.option;
            *chromosome_list=CALLOC((size_t) entry_index, int);
            for (i=0; i < entry_index; i++)
                (*chromosome_list)[i] = Mega2BatchItems[/* 9 */ Chromosomes_Multiple].value.mult_opts[i];
            do_again = verify_chromosomes(*chromosome_list, entry_index);
            if (do_again) {
                /* in input mode we loop over the chromo selection until valid chromosomes
                   are specified */
                return -1;
            }
        } else {
            entry_index = *num_chromo;
        }

        if (InputMode == INTERACTIVE_INPUTMODE) {
            /* write the selections into the batch file */
            Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].value.option = entry_index;
            batchf(Chromosomes_Multiple_Num);

            Mega2BatchItems[/* 9 */ Chromosomes_Multiple].value.mult_opts = CALLOC((size_t) entry_index+1, int);
            for(i=0; i < entry_index; i++) {
                Mega2BatchItems[/* 9 */ Chromosomes_Multiple].value.mult_opts[i] = (*chromosome_list)[i];
            }
            Mega2BatchItems[/* 9 */ Chromosomes_Multiple].value.mult_opts[entry_index] = -99;
            batchf(Chromosomes_Multiple);
            mssgf("You have selected the following chromosome(s) :");
            strcpy(err_msg, " ");
            for(i=0; i < entry_index; i++) {
                CHR_STR((*chromosome_list)[i], prchr);
                grow(err_msg, " %s", prchr);
            }
            mssgf(err_msg); log_line(mssgf);
        }

        /* now set the global_chromo_entries and main_chromocnt */

        k=ReOrderMappedLoci_new(LPedTreeTop, loc_type, &locus_order,
                                entry_index, *chromosome_list);
        linkage_set_locus_order_new(k, locus_order);
        new_chromo_loci_counts(entry_index, *chromosome_list);
        main_chromocnt  = entry_index;
        free(global_chromo_entries);
        global_chromo_entries = *chromosome_list;
        /* Now copy markers on selected chromosomes
           into the global locus list  */
        free(locus_order);
    }
    return 0;
}

static void remove_marker_item(int *num_trs, int *trait_select,
                               int marker_item)

{
    int i=0;
    int mrk_pos=-1;

    while(i < *num_trs) {
        if (trait_select[i] == marker_item) {
            mrk_pos = i;
            break;
        }
        i++;
    }
    if (mrk_pos > -1) {
        while (i < (*num_trs-1)) {
            trait_select[i]=trait_select[i+1];
            i++;
        }
        (*num_trs)--;
    }
}

static void create_trait_string(char *str, int number_traits, int *traits)

{

    int i;
    size_t curr_size = FILENAME_LENGTH;

    strcpy(str, "");

    for (i=0; i < number_traits; i++) {
        if (strlen(str) >= curr_size) {
            str = REALLOC(str, curr_size+FILENAME_LENGTH, char);
        }
        grow(str, " %d", traits[i]);
    }
}

/*---------------------------------------------------------------*/

/* New implementation of trait selection:
   a) Input mode:
   User sees a list of trait loci numbered 1 to N, where N = number of
   trait loci. In combination mode, user also sees a marker item labelled
   N+1. Thus, each item (decremented by 1) is an index into the
   global_trait_entries array.

   b) Batch mode: same as input mode, numbers through N represent trait indexes
   into the global_trait_entries, and any number not between 1 and N
   represents the marker item.

   Selection of trait locus numbers is maintained as indexes and translated
   into actual locus numbers for display purposes, until the selection process
   is over, then global_trait_entries is reallocated to contain only selected
   traits as per these indexes.

*/
static int select_trait_loci(linkage_ped_top *LTop, analysis_type analysis)

{
    int i, j, k, opt=-1, done=0, item; /* flags and counters */
    int *number1, *traits, *trait_select, *traits_left; /* trait lists */
    int num_allowed=0, num_select, num_trait_select, num_traits_left, naff=0;
    int default_trait=0, invalid;
    int             retval;
    char str[FILENAME_LENGTH], strcopy[FILENAME_LENGTH], reply[10];
    char LocusTypeName[5][17] = {
        "UNKNOWN         ",
        "Quantitative    ",
        "Affection Status",
        "Binary  ",
        "Numbered"
    };

    char loop_or_combine[2][40] = {
        "Combine, ordered as in item 2",
        "Trait-first, one trait at a time"
    };

    char one_trait_case[2][30] = {
        "Order as in item 2",
        "Trait-first"
    };

    int base=0, cbase=0, marker_item;
    int overwrite;
    int *covar_items;

    number1 = traits = trait_select = traits_left = covar_items = NULL; /* trait lists */
    /* These following steps are followed with and without a batch file */
    /* No traits, and our analysis requires a trait */

    if (num_traits == 0) {
        if (no_trait_allowed(analysis, 0) == 0) {
            errorvf("%s requires at least one trait locus.\n", ProgName);
            EXIT(OUTPUT_FORMAT_ERROR);
        } else {
            return 1;
        }
    }

    /* Count the number of AFFECTION traits and
       set the default trait to the first affection trait */
    for (i = num_traits-1; i >= 0; i--) {
        if (LTop->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
            default_trait = i;
            naff++;
        }
    }

    if (naff==0 && no_trait_allowed(analysis, 0) == 0) {
        errorvf("%s requires at least one affection trait locus.\n", ProgName);
        EXIT(OUTPUT_FORMAT_ERROR);
    }

    /* If there is exactly one affection trait locus, select it automatically
       if the analysis requires a trait, ignore batch file entry as well */

    if (naff == 1 && no_trait_allowed(analysis, 0) == 0) {
        marker_item = num_traits + 1;
        printf("Found one affection trait locus %s\n",
               LTop->LocusTop->Locus[global_trait_entries[default_trait]].Name);
        printf("%s requires at least one trait, selecting this locus for analysis.\n",
               ProgName);
        /* sleep(2); */
        LoopOverTrait=0;

/*     for (i=0; i < num_reordered; i++) { */
/*       main_order[j] = reordered_marker_loci[i]; j++; */
/*     } */
/*     linkage_set_locus_order_new(j, main_order); */

        /* reset the global_trait_entries and num_traits */
        global_trait_entries = REALLOC(global_trait_entries, (size_t) 2, int);
        global_trait_entries[0] = global_trait_entries[default_trait];
        global_trait_entries[1] = -1;
        num_traits = 2;
        Mega2Status = TRAIT_SELECTED;
        if (InputMode == INTERACTIVE_INPUTMODE) {
            sprintf(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, "%d %d e",
                    default_trait+1, marker_item);
            batchf(Traits_Combine);
        }
        return 1;
    }

    traits = CALLOC((size_t) num_traits+1, int);
    trait_select = CALLOC((size_t) num_traits+1, int); /* selected traits */
    traits_left =  CALLOC((size_t) num_traits+1, int); /* traits left for covariates */
    covariates = CALLOC((size_t) num_traits+1, int); /* selected covariates */
    covar_items = CALLOC((size_t) num_traits, int);

    /* Now we have a choice  of traits */
    if (!batchTRAIT) {
        /*  marker_item = LTop->LocusTop->LocusCnt+1; */
        /*     base=1; */
        marker_item = global_trait_entries[num_traits-1]+2;
        done=0;
        /* No covariates */
        num_covariates = 0;
        /* default single_trait is first trait */
        if (allow_trait_combination(analysis)) {
            LoopOverTrait=0;
            num_trait_select = num_select = 2;
            traits[0] = default_trait;       traits[1] = marker_item;
            sprintf(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, "%d %d e",
                    default_trait+1, marker_item);
            strcpy(strcopy, Mega2BatchItems[/* 15 */ Traits_Combine].value.name);

        } else {
            LoopOverTrait = 1;
            num_trait_select = num_select = 1;    traits[0] = default_trait;
            sprintf(Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name, "%d e", default_trait+1);
            strcpy(strcopy, Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name);
        }
        while (opt != 0) {
            /* Fill in the allowable covariate list according to
               traits selected */
            num_traits_left=0;
            for (i=0; i < num_traits; i++) {
                invalid = 0;
                for (j=0; j < num_trait_select; j++) {
                    if (traits[j] == marker_item) continue;
                    if (global_trait_entries[i] == global_trait_entries[traits[j]-base]
                       || LTop->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
                        invalid = 1;
                    }
                }
                if (!invalid) {
                    traits_left[num_traits_left] = global_trait_entries[i];
                    covar_items[num_traits_left]=i;
                    num_traits_left++;
                }
            }

            item=1;
            if (analysis == TO_SAGE && num_traits > 2 && LoopOverTrait == 0) {
                printf("                      SAGE:\n");
                printf("If multiple affection status loci are selected,\n");
                printf("      FCOR files will be set up in trait sub-directories.\n");
                printf("      SIBPAL files will be set up in the main directory.\n");
                printf("If multiple QTLs are selected,\n");
                printf("      FCOR and SIBPAL files will be set up in the main directory.\n");
            }
	    draw_line();
            printf("Trait and covariate selection menu:\n");
            if (allow_trait_combination(analysis) && LoopOverTrait == 1) {
                printf("    (To select markers only, first toggle item 1.)\n");
            }
            printf("0) Done with this menu - please proceed.\n");

            if (allow_trait_combination(analysis) && num_traits >= 1) {
                if (num_traits > 1) {
                    printf(" %d) Create trait-specific files or combine: [%s]\n",
                           item, loop_or_combine[LoopOverTrait]);
                    if (LoopOverTrait == 0) {
                        printf("    Note: To include markers, select the [MARKERS] entry in item 2.\n");
                    } else {
                        printf("    Note: markers already selected will be automatically included.\n");
                    }
                } else {
                    printf(" %d) Order trait vs. markers: [%s]\n", item,
                           one_trait_case[LoopOverTrait]);
                }

                item++;
            } else {
                if (!allow_trait_combination(analysis)) {
                    printf(" Traits will be analyzed one at a time.\n");
                }
            }
            printf(" %d) Trait%s loci selected: [", item,
                   ((LoopOverTrait == 0)? " and marker" : ""));
            display_trait_names(num_trait_select, traits,
                                LTop->LocusTop, 0, global_trait_entries,
                                base, num_traits+1);
            printf("]\n");

            if (allow_covariates(analysis) && num_traits > 1) {
                item++;
                printf(" %d) Covariates selected: [", item);
                display_trait_names(num_covariates, covariates,
                                    LTop->LocusTop, 0, global_trait_entries,
                                    cbase, num_traits+1);
                printf("]\n");
            }

            printf("Enter option 0 - %d > ", item);
            fcmap(stdin, "%s", reply); newline;
            sscanf(reply, "%d", &opt);
            /*      printf("%d, %d\n", traits[0], marker_item); */
            switch(opt) {
            case 0:
                break;
            case 1:
                if (LoopOverTrait == 0 || allow_trait_combination(analysis)) {
                    LoopOverTrait = TOGGLE(LoopOverTrait);
                    if (LoopOverTrait == 1) {
                        remove_marker_item(&num_trait_select, traits, marker_item);
                        create_trait_string(Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name,
                                            num_trait_select, traits);
                    } else {
                        create_trait_string(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, num_trait_select,
                                            traits);
                    }
                    if (base == 0) {
                        INDX_INC(num_trait_select, traits, i, marker_item);
                        base=1;
                    }

                    retval = check_trait_selection(LTop, analysis, num_trait_select,
                                                   traits, marker_item);
                    /* revert back to default selection */

                    if (retval == 0) {
                        num_select=0;
                        if (LoopOverTrait == 0) {
                            num_trait_select = 2;
                            traits[0] = default_trait;
                            traits[1] = marker_item;
                        } else {
                            num_trait_select = 1;
                            traits[0] = default_trait;
                        }
                    } else if (retval == -1 && LoopOverTrait == 1) {
                        printf("Reason: in Trait Loop over mode, \n");
                        printf("        at least one trait has to be selected.\n");
                        printf("To select 0 traits, first turn on Combine-traits mode\n.");
                        num_trait_select = num_select = 1;
                        traits[0] = default_trait;
                        /* revert back to default selection */
                    } else {
                        num_trait_select=num_select=retval;
                    }
                    break;
                } /* if not (ALLOW_OPT3) fall through to default */
                opt++;
            case 2:
            case 3:
                num_select=0; /* default choice is single locus */
                if (opt == 2 || allow_covariates(analysis)) {
                    if (opt == 2) {
                        base = 1;
                        printf("Select traits from the following list:\n");
                        for (i=0; i < num_traits; i++) {
                            int trait = global_trait_entries[i];
                            printf("%d)     %-15s      %s\n", i+1,
                                   LTop->LocusTop->Locus[trait].Name,
                                   LocusTypeName[LTop->LocusTop->Locus[trait].Type]);
                        }

                        if (LoopOverTrait == 0) {
                            printf("%d)     [MARKERS]\n", num_traits+1);
                            num_allowed = num_traits + 1;
                        } else {
                            num_allowed = num_traits;
                        }
                    } else if (allow_covariates(analysis)) {
                        cbase = 1;

                        printf("Select covariate(s) from the following list:\n");
                        for (i=0; i < num_traits_left; i++) {
                            printf("%d)     %-15s      %s\n", covar_items[i]+1,
                                   LTop->LocusTop->Locus[traits_left[i]].Name,
                                   LocusTypeName[LTop->LocusTop->Locus[traits_left[i]].Type]);
                        }
                        num_allowed = num_traits_left;
                    }


                    num_select = 0; done = 0;
                    overwrite=1;

                    while(!done) {
                        if (opt == 2) {
                            if (no_trait_allowed(analysis,0) &&
                               LoopOverTrait == 0) {
                                printf("Enter trait locus numbers, 'e' to terminate\n");
                                printf( "     or only marker item for 0 traits.\n");
                            } else {
                                printf("Enter trait locus numbers, 'e' to terminate > ");
                            }
                        } else {
                            printf("Enter trait locus numbers, 'return' to terminate\n");
                            printf("      or only 'e' for 0 traits > ");
                        }

                        fflush(stdout);
                        (void)fgets(str, FILENAME_LENGTH-1, stdin);    newline;
                        parse_string(str, trait_select, &num_select, num_allowed, 1);
                        INDX_INC(num_select, trait_select, i, marker_item);

                        if (num_select == 0) break;

                        strcpy(strcopy, str); chomp(strcopy);
                        done=1;

                        /* verify after reading in the string and store loci
                           entered up to now */
                        if (opt == 2) {
/* 	      if (LoopOverTrait == 1) { */
/* 		retval = check_trait_selection(LTop, analysis, num_select, */
/* 				      trait_select, -1); */
/* 	      } */
/* 	      else { */
                            retval = check_trait_selection(LTop, analysis, num_select,
                                                           trait_select, marker_item);
                            /*	      } */

                            if (retval == 0) {
                                num_select=0;
                                if (LoopOverTrait == 0) {
                                    num_trait_select = 2;
                                    traits[0] = default_trait+base;
                                    traits[1] = marker_item;
                                } else {
                                    num_trait_select = 1;
                                    traits[0] = default_trait+base;
                                }
                                /* revert back to default selection */
                            } else if (retval == -1 && LoopOverTrait == 1) {
                                printf("Reason: in Trait Loop over mode, \n");
                                printf("        at least one trait has to be selected.\n");
                                printf("To select 0 traits, first turn on Combine-traits mode\n.");
                                num_select=0;
                                num_trait_select = 1;
                                traits[0] = default_trait+base;
                                /* revert back to default selection */
                            } else {
                                if (overwrite) {
                                    num_trait_select=0;
                                    overwrite=0;
                                }
                                num_select=retval;
                                add_selections(num_select, trait_select,
                                               &num_trait_select, traits,
                                               ((LoopOverTrait == 0)? marker_item : num_traits));
                                remove_from_covar(num_select, trait_select,
                                                  &num_covariates, covariates);
                                if (LoopOverTrait == 0) {
                                    strcpy(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, strcopy);
                                } else {
                                    strcpy(Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name, strcopy);
                                }
                                num_select=0;
                            }
                        } else {
                            /* checking covariates is to simply see if the indexes
                               are within bounds of allowable number of traits */

                            for(j=0; j < num_select; j++) {
                                retval=-1;
                                for (i=0; i < num_traits_left; i++) {
                                    if ((global_trait_entries[trait_select[j] - 1]) == traits_left[i]) {
                                        retval = i;
                                        break;
                                    }
                                }
                                if (retval == -1) {
                                    printf("Invalid selection, previous covariate list unchanged.\n");
                                    num_select=0;
                                    break;
                                }
                            }

                            if (overwrite) {
                                num_covariates=0;
                                overwrite=0;
                            }
                            if (num_select > 0) {
                                add_selections(num_select, trait_select, &num_covariates,
                                               covariates, num_traits);
                                strcpy(Mega2BatchItems[/* 32 */ Covariates_Selected].value.name, strcopy);

                            }
                        }
                    }
                    break;
                }
                /* else opt == 2  AND ALLOW_COVAR is false, simply drop down */
            default:
                warn_unknown(reply);
                break;
            }
            draw_line();
        }

        /* Log covariates */
        if (num_covariates > 0) {
            if (InputMode == INTERACTIVE_INPUTMODE) {
                batchf(Covariates_Selected);
            }
            index_selections(num_covariates, covariates,
                             num_traits, global_trait_entries, 1,
                             marker_item, cbase);
        }

        /* Have to log selections in the batch file */
        if (LoopOverTrait == 1 && num_trait_select == 1) {
            if (InputMode == INTERACTIVE_INPUTMODE) {
                Mega2BatchItems[/* 12 */ Trait_Single].value.option = traits[0] +
                    ((base == 1)? 0 : 1);
                batchf(Trait_Single);
            }
            index_selections(num_trait_select, traits,
                             num_traits, global_trait_entries, 1, -1, base);
        } else {
            if (LoopOverTrait == 1) {
                if (InputMode == INTERACTIVE_INPUTMODE) {
                    batchf(Traits_Loop_Over);
                }
                index_selections(num_trait_select, traits,
                                 num_traits, global_trait_entries, 1, -1, base);
            } else {
                if (InputMode == INTERACTIVE_INPUTMODE) {
                    batchf(Traits_Combine);
                }
                marker_item=index_selections(num_trait_select, traits,
                                             num_traits, global_trait_entries, 1,
                                             marker_item,base);
            }
        }
    }
    /* Batch file input */
    else {
        /*  marker_item = LTop->LocusTop->LocusCnt+1; */
        marker_item = global_trait_entries[num_traits-1]+2;
        /* We assume that the values for traits and covariates are indexes into
           the trait list, not the entire locus list */
        if (ITEM_READ(12)) {
            /*  Single trait locus, test to see if this selection is correct */
            traits = CALLOC((size_t) 1, int);
            num_trait_select = 1;
            traits[0] = Mega2BatchItems[/* 12 */ Trait_Single].value.option;
            LoopOverTrait = 1;
            if (check_trait_selection(LTop,
                                     analysis, 1,
                                     traits,
                                     marker_item) == 0) {
                invalid_value_field(12);
            }
            index_selections(num_trait_select, traits,
                             num_traits, global_trait_entries, 1, -1,1);

/*       Reverse map trait numbers onto indexes for  */
/* 	 global_trait_entries  - check all selections   */
        } else if (ITEM_READ(14)) {
            /* LoopOverTraits traits option */
            LoopOverTrait=1;
            num_trait_select=0;
            /*      num_trait_select=Mega2BatchItems[/ * 13 * / Traits_Num].value.option; */
            traits = CALLOC((size_t) num_traits+1, int);
            parse_string(Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name, traits,
                         &num_trait_select, num_traits, 0);
            number1 = check_for_duplicates(&num_trait_select, traits);
            free(traits); traits = CALLOC((size_t) num_trait_select, int);
            for (i=0; i < num_trait_select; i++) {
                traits[i] = number1[i];
            }
            free(number1);
            INDX_INC(num_trait_select, traits, i, -1);
            remove_marker_item(&num_trait_select, traits, num_traits+1);
            remove_marker_item(&num_trait_select, traits, marker_item);
            if (check_trait_selection(LTop, analysis, num_trait_select,
                                      traits,
                                     -1) == 0){
                invalid_value_field(14);
            }

            index_selections(num_trait_select, traits,
                             num_traits, global_trait_entries, 1, -1,1);

            /* Reverse map trait numbers onto indexes for
               global_trait_entries  - check all selections   */
        } else if (ITEM_READ(15)) {
            /* Combine traits option
               Trait list should also include a marker position =
               num_traits+1
            */
            LoopOverTrait = 0;
            /*      num_trait_select=Mega2BatchItems[/ * 13 * / Traits_Num].value.option; */
            traits = CALLOC((size_t) num_traits+1, int);
            num_trait_select = 0;
            parse_string(Mega2BatchItems[/* 15 */ Traits_Combine].value.name, traits,
                         &num_trait_select, num_traits + 1, 0);
            /*       number1 = check_for_duplicates(&num_trait_select, traits); */
            /*       free(traits); traits = number1; */
            number1 = check_for_duplicates(&num_trait_select, traits);
            free(traits);
            traits = CALLOC((size_t) num_trait_select, int);
            for (i = 0; i < num_trait_select; i++) {
                traits[i] = number1[i];
                if (traits[i] < 0 || traits[i] >= LTop->LocusTop->LocusCnt) {
                    errorvf("The trait number %d is out of the bounds [1, %d].\n",
			    traits[i]+1,
			    LTop->LocusTop->LocusCnt);
                    EXIT(EARLY_TERMINATION);
                }

            }
            free(number1);
            INDX_INC(num_trait_select, traits, i, marker_item);

            if ((num_trait_select = check_trait_selection(LTop, analysis, num_trait_select,
                                                          traits, marker_item)) == 0) {
                invalid_value_field(15);
            }

            /* Reverse map trait numbers onto indexes for
               global_trait_entries  - check all selections   */
            marker_item = index_selections(num_trait_select, traits,
                                           num_traits, global_trait_entries, 1, marker_item, 1);
        }

        if (ITEM_READ(32)) {
            /* Store covariates */
            num_covariates=0;
            invalid=0;
            parse_string(Mega2BatchItems[/* 32 */ Covariates_Selected].value.name,
                         covariates, &num_covariates, num_traits, 0);
            INDX_INC(num_covariates, covariates, i, marker_item);
            index_selections(num_covariates, covariates,
                             num_traits,
                             global_trait_entries, 1,
                             marker_item,1);

            for (k=0; k < num_covariates; k++) {
                for (j=0; j < num_trait_select; j++) {
                    if (traits[j] == covariates[k] ||
                        LTop->LocusTop->Locus[covariates[k]].Type == AFFECTION) {
                        invalid=1; break;
                    }
                }
                if (invalid) {
                    invalid_value_field(32);
                }
            }
        }
    }

    if (num_trait_select == 1 && traits[0] != marker_item) {
        LoopOverTrait = 1;
    }

    /* Log trait and covariate selections */
    if (LoopOverTrait == 0) {
        if (global_chromo_entries[0] > 0 && HasMarkers) {
            mssgf("Output will combine markers and the following selected traits:");
        } else {
            mssgf("Output will combine the following selected traits:");
        }
    } else {
        mssgf("Output will be created with one trait locus ");
        mssgf("at a time with the following traits:");
        /*    marker_item=-1; */
    }

    display_trait_names(num_trait_select, traits, LTop->LocusTop,
                        16, NULL, 0, marker_item);
    newline;

    if ((analysis == TO_SAGE && num_trait_select > 1)
       && (LoopOverTrait == 0)) {
        mssgf("                      SAGE: Combine-traits");
        mssgf("If multiple affection status loci are selected,");
        mssgf("      FCOR files will be set up in trait sub-directories.");
        mssgf("      SIBPAL files will be set up in the main directory.");
        mssgf("If multiple QTLs are selected,");
        mssgf("      FCOR and SIBPAL files will be set up in the main directory.");
    }

    if (num_covariates > 0) {
        mssgf("Output will include the following covariates:");
        display_trait_names(num_covariates, covariates, LTop->LocusTop,
                            16, NULL, 0, marker_item);
        newline;
        /* Change the class of covariates to COVARIATES */

        for (i=0; i < num_covariates; i++) {
            LTop->LocusTop->Locus[covariates[i]].Class = COVARIATE;
        }
    }
    /* NM - removed code to re-create a combined set of loci,
       instead keep these as list of traits, list of covariates and
       a combined locus list.
       GLOBAL global_trait_entries : traits
       GLOBAL covariates : covariates
       GLOBAL reordered_marker_loci : combined
       If LoopOverTraits is set, then reordered_marker_loci will
       not contain traits.

    */

    /* first switch the traits to global_trait_entries */

    num_traits = num_trait_select;
    if (HasMarkers == 0) {
        linkage_set_locus_order_new(num_trait_select, traits);
        NumChrLoci=num_reordered;
        ChrLoci=reordered_marker_loci;
        ManualReorder=1;
    }

    global_trait_entries = REALLOC(global_trait_entries, (size_t) num_traits+1, int);
    for (i = 0; i < num_trait_select; i++) {
        global_trait_entries[i] = (traits[i] != marker_item) ? traits[i] : -1;
    }
    global_trait_entries[num_traits] = -99;

    Mega2Status = TRAIT_SELECTED;
    get_trait_list(LTop->LocusTop, 0);

    /*  cfree(main_order);  */
    free(traits); free(trait_select);
    return 1;
}

/* see if the option allows 0 traits, removed genehunter from the list */
/* This is meant for at least one affection trait */

static int no_trait_allowed(analysis_type analysis, int mssg)

{
    /* This function returns FALSE(0) if the analysis option requires at
     * least one trait, TRUE (1) otherwise.  Add your KEYWORD to the first
     * block of the case statement, if it needs to have traits.
     */

    if (analysis->allow_no_trait()) {
        if (mssg) {
            sprintf(err_msg, "%s requires at least one trait.", ProgName);
            errorf(err_msg);
        }
        return 0;
    } else
        return 1;
}

/*
  static int no_aff_trait_allowed(analysis_type analysis, int mssg)
  {
    if (analysis->allow_no_aff_trait()) {
      if (mssg) {
        sprintf(err_msg, "%s requires at least one trait.", ProgName);
        errorf(err_msg);
      }
    return 0;
    } else
  return 1;
  }
*/

static int allow_trait_combination(analysis_type analysis)

{
/*
 * This function return TRUE if traits are
 * allowed to be combined in the output as a single set of files, and
 * traits do not need to come first in the locus list. Add
 * your KEYWORD to this list if your analysis option allows this.
 */

    if (analysis->allow_trait_combination())
        return 1;
    else
        return 0;
}

static int allow_covariates(analysis_type analysis)

{
/*
 * If QTLs can be specially designated as
 * covariates for analysis by your program, add your KEYWORD to the first
 * block of the case.
 */
    if (analysis->allow_covariates())
        return 1;
    else
        return 0;
}


/*
   Check that trait loci were properly selected, some analysis options do
   not allow quant loci. These are: ASPEX, GH-plus, apm, apm-mult, splink,
   simulate, segregation summary, liability groups, tdtmax, and simwalk2.
*/
static int  quant_allowed(linkage_locus_rec Locus, analysis_type analysis)
{
    /* now check for QTLs against analysis type */
    if (QTL_DISALLOW  && Locus.Type == QUANT) {
        printf("ERROR: %s does not allow QTLs.\n", ProgName);
        return 0;
    }
    return 1;
}

#ifdef DEFUNCT
static void get_trait_positions(linkage_locus_top *LTop,
				analysis_type analysis)
{
    int i, i1, j, *traits = CALLOC((size_t) LTop->LocusCnt, int);
    int non_monotonic;
    char str[FILENAME_LENGTH];
    char *c;
    double lastval;

    /*
     * For some analysis options, a single trait can be
     * placed between markers, using the reordering menu option 2, and this
     * trait can have an unknown position. This function checks to see that
     * there is only one such unmapped locus. Add your KEYWORD to the case
     * block at the top, if this is appropriate.
     */

    switch(analysis) {
    case TO_SLINK:
    case TO_SIMULATE:
        /* Check thetas */
        j=0;
        for (i1=0; i1 < num_reordered; i1++) {
            i = reordered_marker_loci[i1];
            if ((LTop->Locus[i].Type == AFFECTION ||
                 LTop->Locus[i].Type == QUANT) &&
                 (i > 0) &&
		 (i < LTop->LocusCnt-1)) {

	      if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                non_monotonic=((LTop->Marker[i].pos_avg < 0)? 1: 0);
                if (!non_monotonic) {
                    if ((LTop->Marker[i].pos_avg >= LTop->Marker[i-1].pos_avg) &&
                        (LTop->Marker[i+1].pos_avg >= LTop->Marker[i].pos_avg)) {
                        non_monotonic=1;
                    }
                }
	      } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
			 genetic_distance_sex_type_map == FEMALE_GDMT) {
                non_monotonic=((LTop->Marker[i].pos_female < 0)? 1: 0);
                if (!non_monotonic) {
                    if ((LTop->Marker[i].pos_female >= LTop->Marker[i-1].pos_female) &&
                        (LTop->Marker[i+1].pos_female >= LTop->Marker[i].pos_female)) {
                        non_monotonic=1;
                    }
                }
	      }

                if (non_monotonic) {
                    traits[j]=i; j++;
                }
            }
        }
        if (j > 0) {
            printf("Locus order contains multiple unmapped trait loci.\n");
            printf("%s requires all loci but one to have genetic positions.\n",
                   ProgName);
            while(1) {
                if (j == 1) {
                    printf("Enter position for %s (2nd unmapped locus) in cM > ",
                           LTop->Locus[traits[0]].Name);
                } else {
                    printf("Enter upto %d trait locus positions in cM,\n", j);
                    printf("starting with 2nd unmapped locus %s.\n",
                           LTop->Locus[traits[0]].Name);
                    printf("if fewer than %d positions are provided, \n", j);
                    printf("the last position will be assigned to the remainder of the unmapped loci > ");
                }
                fflush(stdout);
                (void)fgets(str, FILENAME_LENGTH-1, stdin);
                newline;
                if (strcmp(str, "")) {
                    c = strtok(str, " "); lastval=atof(c);
                    for(i=0; i < j; i++) {
                        LTop->Marker[traits[i]].pos_avg = lastval;
                        c = strtok(NULL, " ");
                        if (c != NULL) {
                            lastval=atof(c);
                        }
                    }
                    break;
                }
            }
            /* set chromocome numbers */
            for(i=0; i < j; i++) {
                LTop->Marker[traits[i]].chromosome = LTop->Marker[traits[i]-1].chromosome;
            }
        }
        break;
    default:
        break;
    }
}
#endif

static int     check_trait_selection(linkage_ped_top *Top,
				     analysis_type analysis,
				     int num_select,
				     int *trait_order,
				     int marker_item)
{
    int j, tr;
    int marker_position=-1;

    if (num_select == 0) {
        return -1;
    }


    /* are numbers within bounds and marker item present ? */
    if (marker_item > 0) {
        for (j = num_select-1 ; j >= 0; j--) {
            if (trait_order[j] == num_traits+1) {
                marker_item = trait_order[j] = marker_item;
                marker_position=j;
                break;
            }
        }
    }

    /* at least one trait,
       this means at least one trait */
    if (LoopOverTrait == 0 && num_select < 2 &&
        trait_order[0] == marker_item && !no_trait_allowed(analysis, 1)) {
        return 0;
    }

    for (j=0; j < num_select; j++) {
        if (trait_order[j] != marker_item &&
            (trait_order[j] < 1 || trait_order[j] > num_traits)) {
            printf("Invalid selection %d, previous list unchanged.\n",
                   trait_order[j]);
            return 0;
        }
    }

    if (LoopOverTrait == 1 && marker_position >= 0) {
        /* take out the marker item */
        tr=0;
        for (j=marker_position; j < (num_select-1); j++) {
            trait_order[j] = trait_order[j+1];
        }
        num_select--;
    }
/*   else if (LoopOverTrait == 0 && marker_position < 0) {  */
/*     /\* add marker item at the end*\/     */
/*     warnf("User did not include [MARKERS], appended this item to selected traits."); */
/*     num_select++;  */
/*     trait_order[num_select-1] = marker_item;  */
/*   }  */
    else if (LoopOverTrait == 0 && marker_position < 0) {
        HasMarkers=0;
    }

    if (analysis == TO_GHMLB && LoopOverTrait == 0) {
        if (!check_ghmlb_affs(Top, trait_order, num_select, marker_item)) {
            return 0;
        }
    }


    /* quant loci allowed */
    for (j=0; j < num_select; j++) {
        if (LoopOverTrait == 0 && global_chromo_entries[0] > 0 &&
           trait_order[j] == marker_item) ;
        else {
            tr = global_trait_entries[trait_order[j]-1];
            if (quant_allowed(Top->LocusTop->Locus[tr], analysis) == 0) {
                printf("Please re-select all loci.\n");
                return 0;
            }
        }
    }

    if (analysis == TO_SAGE && LoopOverTrait == 0 && num_select > 0) {
        if (!all_same_type(Top, trait_order, num_select, marker_item)) {
            printf("Please re-select all loci.\n");
            return 0;
        }
    }

    if (analysis == TO_GeneHunter && LoopOverTrait == 0) {
        if (!check_gh_quants(Top, trait_order,
                            num_select, marker_item)) {
            printf("Please re-select all loci.\n");
            return 0;
        }
    }

    if ((SIMWALK2(analysis)) && LoopOverTrait == 0) {
        if (!check_simwalk2_affs(analysis, num_select, trait_order,
                                marker_item,1)) {
            printf("Please re-select all loci.\n");
            return 0;
        }
    }

    return num_select;

}

static int check_ghmlb_affs(linkage_ped_top *Top, int *traits,
			    int num_affs, int mrk_pos)

{
    /* check to see if there is at most 1 quant + at most 1 affection
       for mlb and that the order is right
    */
    switch(num_affs) {
    case 0:
        printf("ERROR: At least one trait required.\n");
        return 0;
    case 1:
        return 1;
    case 2:
        if (mrk_pos == -1) {
            /* no MARKERS item in the trait order - entry through reordring opt 2*/
            if (Top->LocusTop->Locus[traits[0]-1].Type == QUANT &&
                Top->LocusTop->Locus[traits[1]-1].Type == AFFECTION)
                return 1;
            else {
                printf("ERROR: Selection of multiple trait loci requires\n");
                printf("       locus order to be quantitative locus followed by affection status.\n");
                return 0;
            }
        } else {
            /* trait order may contain a MARKERS item */
            if (traits[1] == mrk_pos) {
                /* don't care what traits[0]'s type is,
                   may have to add a dummy affction status if
                   it is of type quant. */
                return 1;
            } else {
                if (Top->LocusTop->Locus[traits[0]-1].Type == QUANT &&
                    Top->LocusTop->Locus[traits[1]-1].Type == AFFECTION){
                    /* ideal case, no dummy required later on */
                    return 1;
                } else {
                    /* anything else is wrong order */
                    printf("ERROR: Locus order invalid\n");
                    printf("       locus order to be quantitative locus followed by affection status.\n");
                    return 0;
                }
            }
        }
    case 3:
        /* the third item has to be the markers item */
        if (Top->LocusTop->Locus[traits[0]-1].Type == QUANT &&
            Top->LocusTop->Locus[traits[1]-1].Type == AFFECTION &&
            traits[2] == mrk_pos && mrk_pos != -1)
            return 1;
        else {
            printf("ERROR: Selection of multiple trait loci requires\n");
            printf("       locus order to be quantitative locus followed by affection status.\n");
            return 0;
        }
    default:
        printf("Too many trait loci selected, maximum 2 allowed.\n");
        return 0;
    }
/*    return 1; can not get here; all cases above return */
}

static int all_same_type(linkage_ped_top *Top, int *traits, int num_affs,
			 int markers_item)

{
    int i;
    linkage_locus_type loc_type;

    /* check that the last item is markers */
    if (markers_item != -1) {
        for (i=0; i < num_affs; i++) {
            if (traits[i] == markers_item && i != (num_affs-1)) {
                printf("ERROR: SAGE requires traits to be first in locus order.\n");
                return 0;
            }
        }
    }

    loc_type = Top->LocusTop->Locus[traits[0]].Type;

    for (i=1; i < num_affs-1; i++) {
        if (Top->LocusTop->Locus[traits[i]-1].Type != loc_type) {
            printf("ERROR: SAGE requires all trait loci to be of the same type.\n");
            return 0;
        }
    }

    return 1;
}

static int check_gh_quants(linkage_ped_top *Top, int *traits,
			   int number_traits, int mrk_pos)

{
    int l, tr, affec=0;
    int found_quant=0;

    /* first check that all quants are at the end of the list
     */

    for(l=0; l < number_traits && traits[l]; l++) {
        tr=((mrk_pos == -1)? traits[l] : traits[l]-1);
        if (traits[l] == mrk_pos) ;
        else if (Top->LocusTop->Locus[tr].Type == AFFECTION) {
            affec++;
        }
    }
    if (affec > 1) {
        printf("ERROR: GeneHunter requires exactly one affection\n");
        printf("       status locus.\n");
        return 0;
    }

    if (affec==1) {
        if (traits[0]  == mrk_pos) {
            printf("ERROR: GeneHunter requires the affection status\n");
            printf("       locus to be first in the order.\n");
            return 0;
        }
        tr = ((mrk_pos == -1)? traits[0] : traits[0]-1);
        if (Top->LocusTop->Locus[tr].Type != AFFECTION) {
            printf("ERROR: GeneHunter requires the affection status\n");
            printf("       locus to be first in the order.\n");
            return 0;
        }
    }

    for (l=0; l < number_traits; l++) {
        tr = ((mrk_pos == -1)? traits[l] : traits[l]-1);

        if (traits[l] == mrk_pos && found_quant) {
            printf("ERROR: GeneHunter requires QTLs to be last in the order.\n");
            return 0;
        } else if (Top->LocusTop->Locus[tr].Type == QUANT &&
                found_quant) ;
        else if (Top->LocusTop->Locus[tr].Type == QUANT &&
                 !found_quant)
            found_quant=1;
        else if (Top->LocusTop->Locus[tr].Type != QUANT &&
                 found_quant) {
            printf("ERROR: GeneHunter requires QTLs to be last in the order.\n");
            return 0;
        }
    }
    return 1;
}

static int check_simwalk2_affs(analysis_type analysis,
			       int num_select, int *slct,
			       int marker_pos, int sw2check)

{
    /* All simwalk2 options can handle at most one trait, so there should
       be at most two items in select, and the second one should
       be marker_pos */


    if (!sw2check) {
        /* found a previous error */
        return 1;
    }

    if (marker_pos == -1 &&  num_select > 1) {
        /* No marker item inside select */
        printf("ERROR: SimWalk2 allows at most one trait locus,\n");
        printf("       Selection contains multiple traits.\n");
        return 0;
    }
    if (num_select > 2) {
        printf("ERROR: SimWalk2 allows at most one trait locus,\n");
        printf("       Selection contains multiple traits.\n");
        return 0;
    }
    if (num_select == 2) {
        if (slct[1] != marker_pos) {
            printf("ERROR: SimWalk2 requires trait locus to be first in the order.\n");
            printf("       Option 3 may only be used to \"select markers only\".\n");
            return 0;
        }
    }
    return 1;
}

#ifdef DEFUNCT
int ReOrderMappedLoci(linkage_ped_top *Top, int *numchr)
{

    char            reply[4];
    int             *tmpIndex, tmpswap, maxIndex, *tmpOrder;
    int             *tmpNumber, i, j, count;
    linkage_locus_top *LTop;

    strcpy(reply, "No");
    LTop = Top->LocusTop;
    tmpNumber=CALLOC((size_t) LTop->LocusCnt, int);

    if (*numchr <= 0) {
        count=0;
        switch(*numchr) {
        case -2:
            for (i=0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == AFFECTION) {
                    tmpNumber[count] = i+1; count++;
                }
            }
            break;

        case -1:
            for (i=0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == QUANT) {
                    tmpNumber[count] = i+1; count++;
                }
            }
            break;

        case 0:
            for (i=0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == AFFECTION ||
                    LTop->Locus[i].Type == QUANT) {
                    tmpNumber[count] = i+1; count++;
                }
            }
            break;
        default:
            break;
        }
        tmpNumber[count]= -99;
        if (InputMode == INTERACTIVE_INPUTMODE) {
            Mega2BatchItems[/* 7 */ Chromosome_Single].value.option = 0;
            batchf(Chromosome_Single);
            /*      Mega2BatchItems[/ * 10 * / Loci_Selected_Num].value.option = count;
                    batchf(Loci_Selected_Num);
            */
                        strcpy(Mega2BatchItems[/* 11 */ Loci_Selected].value.name, "");
            for(i=0; i < (count-1); i++) {
                grow(Mega2BatchItems[/* 11 */ Loci_Selected].value.name, " %d",
                     tmpNumber[i]);
            }
            strcat(Mega2BatchItems[/* 11 */ Loci_Selected].value.name, " e");
            batchf(Loci_Selected);
        }

        linkage_set_locus_order(Top, tmpNumber);
        free(tmpNumber);
        return 1;
    }

    tmpIndex=CALLOC((size_t) LTop->LocusCnt, int);
    tmpOrder=CALLOC((size_t) LTop->LocusCnt, int);

    for (i=0; i<LTop->LocusCnt; i++) {
        tmpNumber[i]=-1; tmpIndex[i]= -99;
        tmpOrder[i]=-99;
    }
    count=0;
    for (i = 0; i < LTop->LocusCnt; i++) {
        if ( (LTop->Locus[i].Class == MARKER && LTop->Marker[i].chromosome == *numchr) ||
            LTop->Locus[i].Type == AFFECTION ||
            LTop->Locus[i].Type == QUANT) {
            /* only indexes that are mapped */
            tmpIndex[count] = i;
            tmpNumber[count] = LTop->Locus[i].number;
            count++;
        }
    }  /* index as in the locus order */

    /* We don't have to reorder loci by map, this has already been done */
    if (Mega2Status == INSIDE_ANALYSIS) {
        for (i=0; i < LTop->LocusCnt; i++) {
            if (tmpIndex[i] != -99) {
                tmpOrder[i] = tmpIndex[i]+1;
            }
        }
    } else {
        /* use locus.number to order now, leave it up completely to the user to
           specify the correct order inside the map file */
        for (i = 0; i < count; i++) {
            maxIndex = i;
            for (j = i ; j < count; j++)  {
                if (tmpNumber[j] < tmpNumber[maxIndex]) {
                    maxIndex = j;
                    /* do a swap between maxIndex and count-i for both numbers and Index*/
                }
            }
            tmpOrder[i]=tmpIndex[maxIndex]+1;
            tmpswap = tmpNumber[i];
            tmpNumber[i] = tmpNumber[maxIndex];
            tmpNumber[maxIndex] = tmpswap;
            tmpswap = tmpIndex[i];
            tmpIndex[i] = tmpIndex[maxIndex];
            tmpIndex[maxIndex] = tmpswap;
        }
    }

    linkage_set_locus_order(Top, tmpOrder);
    /* adjust_distances(Top); */

    free(tmpOrder);
    free(tmpNumber); free(tmpIndex);
    if (Mega2Status == INSIDE_ANALYSIS) {
       get_trait_list(Top->LocusTop, 0);
    }
    return 1;
}
#endif

/* ReOrderMappedLoci_new looks at a list of selected chromosomes,
   main_chromocnt and chromo_loci_count,
   and fills in loci into locus_order in map order
   This should be called before new_chromo_loci_counts
*/

static int ReOrderMappedLoci_new(linkage_ped_top *Top, int locus_type,
				 int **locus_order, int numchr, int *chromosomes)
{

    int             i, j, k, count, total_count;
    tmpnum          *tmpNumber;

    int (*comp_index)(const void *, const void *) = compare_index;

    linkage_locus_top *LTop = Top->LocusTop;

    /* copy one chromosome worth of indexes and marker numbers
       into temporary arrays, then sort the indexes by their
       corresponding numbers */

    if (locus_type <= 0) {
        count = 0;
        *locus_order = CALLOC((size_t) num_traits, int);
        switch (locus_type) {
        case -2: /* locus_type == -2: select only AFFECTION status trait loci */
            for (i = 0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == AFFECTION) {
                    (*locus_order)[count] = i;
                    count++;
                }
            }
            break;

        case -1: /* locus_type == -1: select only QUANT trait loci */
            for (i = 0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == QUANT) {
                    (*locus_order)[count] = i;
                    count++;
                }
            }
            break;
        case 0: /* locus_type == 0: select only AFFECTION or QUANT trait loci */
            for (i = 0; i < LTop->LocusCnt; i++) {
                if (LTop->Locus[i].Type == AFFECTION || LTop->Locus[i].Type
                    == QUANT) {
                    (*locus_order)[count] = i;
                    count++;
                }
            }
            break;
        default:
            break;
        }
        return count;
    } else {
        total_count = 0;
        for (j = 0; j < numchr; j++) {
            /* find out the total number of marker loci */
            for (i = 0; i < NumChromo; i++) {
                if (chromosomes[j] == global_chromo_entries[i]) {
                    if (i == 0) {
                        total_count += chromo_loci_count[i];
                    } else {
                        total_count += chromo_loci_count[i]
                            - chromo_loci_count[i - 1];
                    }
                    break;
                }
            }
        }

        *locus_order = CALLOC((size_t) total_count, int);

        k = 0;
        for (j = 0; j < numchr; j++) {
            for (i = 0; i < NumChromo; i++) {
                if (chromosomes[j] == global_chromo_entries[i]) {
                    break;
                }
            }
            if (i == 0) {
                count = chromo_loci_count[i];
            } else {
                count = chromo_loci_count[i] - chromo_loci_count[i - 1];
            }
            tmpNumber = CALLOC((size_t) count, tmpnum);
            count = 0;

            /* find loci on chromosome */
            for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
                if (LTop->Marker[i].chromosome != UNKNOWN_CHROMO
                    && LTop->Marker[i].chromosome == chromosomes[j]) {

		    if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
		      tmpNumber[count].index = LTop->Marker[i].pos_avg;
		    } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
			       genetic_distance_sex_type_map == FEMALE_GDMT) {
                        tmpNumber[count].index = LTop->Marker[i].pos_female;
		    } else {
		      errorvf("You must have specified a genetic distance sex map type.\n");
		      EXIT(DATA_TYPE_ERROR);
		    }
                    if (base_pair_position_index >= 0)
                        tmpNumber[count].phys = Top->EXLTop->EXLocus[i].positions[base_pair_position_index];
                    else
                        tmpNumber[count].phys = i;
                    tmpNumber[count].number = i;
                    count++;
                }
            }
            qsort((void *) tmpNumber, (size_t) count, sizeof(tmpnum), comp_index);
            for (i = 0; i < count; i++) {
                (*locus_order)[k + i] = tmpNumber[i].number;
            }
            free(tmpNumber);
            k += count;
        }

        return total_count;
    }
}

/* Simply copies over the ordered marker loci into the global list */
/* does not change the linkage-ped-top */

int linkage_set_locus_order_new(int num_loci, int *order)

{
    int i;

    num_reordered = num_loci;

    if (reordered_marker_loci != NULL) {
        free(reordered_marker_loci);
    }
    reordered_marker_loci = CALLOC((size_t) num_loci, int);
    for (i=0; i < num_loci; i++) {
        reordered_marker_loci[i] = order[i];

    }
    return 1;
}

#ifdef DEFUNCT
/* linkage_set_locus_order()  (Original by M. Schroeder; modified by DEW)
 * (also note that after this function is executed,
 *  the original data structure linkage_ped_top is
 *  corrupted (i.e. only the newly selected subset of loci is valid).
 * Given the pedigree top and an array of integers of ID's in
 * which the loci should be ordered, order them.
 *
 * Ie.: If order[] = { 1, 3, 2 }, order the loci so that
 * the first remains first, the third is moved into the
 * second place, and the second one is moved to the end.
 *
 * Only existing loci may be in order[], and they may each
 * only be in order[] once. Not all need to be present,
 * however; those that are present will be first in the
 * Locus[] and Allele_*[] arrays, with those that were not
 * present moved to the end and put into some random order.
 *
 * order[] must have as many entries as there are loci,
 * but if not all loci are present then the extra entries
 * in order[] must be set to -99
 *
 * Example: You have 3 loci. You wish to save the first and
 * the third, in that order, but you don't care about the second.
 * In that case, order[] should be { 1, 3, -99 }.
 *
 * Note that not only must the order of the Locus[] be changed,
 * but also the Allele_1[] and Allele_2[] of all members of
 * all pedigrees. Yuck!!!
 *
 * Return non-zero on error.
 * Modified by nandita: take out the pedigree record modifications until
 * later.
 */
/* the terminating order number is -99 */

int             linkage_set_locus_order(linkage_ped_top *Top, int *order)
{
    int             i, j;
    register int    tmpi;
    register linkage_locus_top *LTop = Top->LocusTop;
    pheno_pedrec_data tmpdatapp;
    register linkage_ped_rec *Entry;
    register person_node_type *PEntry;
    int             ped, entry_;
    int            *place = CALLOC((size_t) LTop->LocusCnt, int);

//xx rt
    /*
     * THIS IS VERY WRONG AND VERY BAD -- rvb
     */

    /* set up place to initial order which is always... */
    for (i = 0; i < LTop->LocusCnt; i++)
        place[i] = i + 1;

    /* There's no really good way to do this. */
    for (i = 0; (i < LTop->LocusCnt) && (order[i] != -99); i++) {
        if (place[i] != order[i])  {
//xx rt
            abort();
            if (i == LTop->LocusCnt - 1)   {
                /* Last locus, and they don't agree? */
                free(place);
                errorvf("Internal error setting locus order!\n");
                EXIT(SYSTEM_ERROR);
            } else   {
                for (j = i + 1; (j < LTop->LocusCnt) && (place[j] != order[i]); j++);
                linkage_swap_loci(&(LTop->Locus[i]), &(LTop->Locus[j]));
                /* Now for all pedigree members... */
                for (ped = 0; ped < Top->PedCnt; ped++)
                    if (Top->pedfile_type == POSTMAKEPED_PFT) {
                        for (entry_ = 0; entry_ < Top->Ped[ped].EntryCnt; entry_++) {
                            Entry = &(Top->Ped[ped].Entry[entry_]);

                            tmpdatapp = Entry->Pheno[i];
                            Entry->Pheno[i] = Entry->Pheno[j];
                            Entry->Pheno[j] = tmpdatapp;
                        }
                    } else {
                        /* pre makeped format */
                        for (entry_ = 0; entry_ < Top->PTop[ped].num_persons; entry_++) {
                            PEntry=&(Top->PTop[ped].persons[entry_]);
                            tmpdatapp= PEntry->pheno[i];
                            PEntry->pheno[i] = PEntry->pheno[j];
                            PEntry->pheno[j] = tmpdatapp;
                        }
                    }

                /* don't forget to change place */
                tmpi = place[i];
                place[i] = place[j];
                place[j] = tmpi;
            }
        }
    }
    free(place);

    /* Set locus count to the new count */
    j = 0;
    for (i = 0; i < LTop->LocusCnt && order[i] != -99; i++)
        j += 1;

    /* we should free up locus records
       and entry data which are not mapped any more */
    for(i=j; i < LTop->LocusCnt; i++) {
        free_all_from_llocusrec(&(LTop->Locus[i]));
    }

    Top->LocusTop->LocusCnt = j;
    /* The Data size is modified as necessary, otherwise we leave
       unfreed memory */

    for(ped=0; ped < Top->PedCnt; ped++) {
        if (Top->pedfile_type == POSTMAKEPED_PFT) {
            for(entry_=0; entry_ < Top->Ped[ped].EntryCnt; entry_++) {
                Top->Ped[ped].Entry[entry_].Pheno =
                    REALLOC(Top->Ped[ped].Entry[entry_].Pheno,
                            (size_t) Top->LocusTop->PhenoCnt, pheno_pedrec_data);
            }
        } else {
            for(entry_=0; entry_ < Top->PTop[ped].num_persons; entry_++) {
                Top->PTop[ped].persons[entry_].pheno =
                    REALLOC(Top->PTop[ped].persons[entry_].pheno,
                            (size_t) Top->LocusTop->PhenoCnt, pheno_pedrec_data);
            }
        }
    }

    return 0;
}
#endif

static int index_selections(int num_select, int *selections,
			    int num_values, int *values,
			    int reverse, int marker_item,
			    int base)

{

    /*  if reverse is 0, then selections are a subset of
        values,
        -> replace selections with corresponding indexes into
        values,

        if reverse is 1, then selections are indexes
        into values,
        -> replace selections with the values
        assume that selections are 1-based, but arrays are 0-based.
    */

    int i, j;
    int found;

    if (reverse) {
        /* first set the Marker item */
        if (marker_item > -1) {
            for (j = 0; j < num_select; j++) {
                if (selections[j] >= marker_item) {
                    selections[j] = values[num_values-1]+base+1;
                    HasMarkers=1;
                    marker_item = values[num_values-1]+base+1;
                    /*	  printf("%d %d\n", j+1, marker_item); */
                    break;
                }
            }
        }

        for (j = 0; j < num_select; j++) {
            if (marker_item > -1 && selections[j] == marker_item) {
                continue;
            } else {
                /*	printf("%d, %d\n", j, selections[j]); */
                if ((selections[j] > num_values) || (selections[j] < base)) {
                    printf("Cannot reverse index, index out of bounds.\n");
                    return 0;
                } else {
                    selections[j] = values[selections[j] - base];
                }
            }
        }
    } else {
        for (j = 0; j < num_select; j++) {
            found = -1;
            for (i = 0; i < num_values; i++) {
                if (values[i] == (selections[j] - 1)) {
                    selections[j] = i;
                    found = 1;
                    break;
                }
            }
            if (found == -1) {
                if (marker_item > -1) {
                    HasMarkers=0;
                    if (selections[j] == marker_item) {
                        HasMarkers=1;
                    }
                } else {
                    printf("Could not find selection in value list.\n");
                    return 0;
                }
            }
        }
    }
    return marker_item;
}

static int add_selections(int num_new, int *newer, int *num_old,
			  int *old, int total)

{

    /* assume that *old terminates with a -99, and start
       filling after *num_old elements, after checking for
       duplicates */

    int exists, i, j;

    for (j=0; j < num_new; j++) {
        exists = 0;
        for (i=0; i < *num_old; i++) {
            if (old[i] == newer[j]) { exists=1; break; }
        }
        if (!exists) {
            if (*num_old >= total) {
                printf("Cannot add to old list, not enough space.\n");
                return 0;
            }
            old[*num_old] = newer[j];
            (*num_old)++;
        }
    }
    return 1;
}

/* remove selections from covariate list (if they exist) */
static void remove_from_covar(int num_select, int *slct,
                              int *num_covar, int *local_covariates)

{

    int *tmp;

    int i, j, exists, num_new=0;

    if (*num_covar == 0) {
        return;
    }
    tmp = CALLOC((size_t) *num_covar, int);
    for (j=0; j < *num_covar; j++) {
        exists = 0;
        for (i = 0; i < num_select; i++) {
            if (local_covariates[j] == slct[i]) {
                exists=1; break;
            }
        }
        if (!exists) {
            tmp[num_new]=local_covariates[j];
            num_new++;
        }
    }

    for (j=0; j < num_new; j++) {
        local_covariates[j] = tmp[j];
    }
    if (*num_covar != num_new) {
        warnvf("previous selected covariates selected as traits.\n");
        warnvf("Removing these from covariate list.\n");
    }

    *num_covar = num_new;
    free(tmp);
}


//
// Determine if there are m/f or autosomial chromosomes in our list
//
int x_linked_check(int chromocnt, int *chromosomes, analysis_type analysis)
{
/*   char creply[10], stars[4]; */
/*   int i, choice; */
    int i;
    int autosomal_chr=0, sex_chr=0;
    int sex_linked = 0;

    // Look for one sex chromosome...
    for (i=0; i < chromocnt; i++) {
        if (chromosomes[i] == SEX_CHROMOSOME ||
            chromosomes[i] == MALE_CHROMOSOME) {
            sex_chr = 1;
            break;
        }
    }

    // Since we will only have autosomes...
    if (sex_chr == 0) return 0;

    // look for one autosome...
    for (i=0; i < chromocnt; i++) {
        if (chromosomes[i] != SEX_CHROMOSOME &&
            chromosomes[i] != MALE_CHROMOSOME) {
            autosomal_chr = 1;
            break;
        }
    }

    if (autosomal_chr == 1 && sex_chr == 1) {
        sex_linked = 2;
    } else if (autosomal_chr == 0 && sex_chr == 1) {
        sex_linked = 1;
    } else if (sex_chr == 0) {
        sex_linked = 0;
    }

    /*   if (Mega2BatchItems[/ * 30 * / Xlinked_Analysis_Mode].items_read == 1) { */
    /*     sex_linked = Mega2BatchItems[/ * 30 * / Xlinked_Analysis_Mode].value.option;  */
/*   } */
/*   else { */
/*     strcpy(stars, "  *"); */
/*     choice=3; sex_linked = 2; */
/*     while(choice) { */
/*       printf("Data includes sex-chromosome (chr %d), ", SEX_CHROMOSOME); */
/*       printf("Select analysis mode from the menu below:\n"); */
/*       printf("0) Done with this menu - please proceed.\n"); */
/*       printf("%c1) Treat ALL markers as autosomal.\n", stars[0]); */
/*       printf("%c2) Treat ALL markers as sex-linked.\n", stars[1]); */
/*       printf("%c3) Treat ONLY markers on X(%d) or Y(%d) as sex-linked.\n",  */
/* 	     stars[2], SEX_CHROMOSOME, MALE_CHROMOSOME); */

/*       printf("Enter selection 0 - 3 > "); */
/*       fcmap(stdin, "%s", creply); newline; */
/*       sscanf(creply, "%d", &choice); */
/*       switch(choice) { */
/*       case 0:  */
/* 	break; */
/*       case 1:  */
/* 	strcpy(stars, "*  ");	 */
/* 	sex_linked= choice - 1; */
/* 	break; */
/*       case 2: */
/* 	strcpy(stars, " * "); */
/* 	sex_linked= choice - 1; */
/* 	break; */
/*       case 3: */
/* 	strcpy(stars, "  *");	       */
/* 	sex_linked= choice - 1; */
/* 	break; */
/*       default: */
/* 	printf("Unknown option %s.\n", creply); */
/* 	break; */
/*       } */
/*     } */
/*   }     */

    switch(sex_linked) {
    case 0: // did not find m/f-chromosome...
        warnf("Treating ALL markers as autosomal.");
        break;

    case 1: // found m/f but no autosomes...
        if (SIMWALK2(analysis) || analysis->forbid_sex_linked_loci()){
            errorvf("%s does not allow sex-linked loci.\n", ProgName);
            EXIT(EARLY_TERMINATION);
        } else {
            warnf("Treating ALL markers as sex-linked.");
        }
        break;
    case 2: // found m/f and autosomes...
        if (SIMWALK2(analysis) || analysis->forbid_sex_linked_loci()){
            if (autosomal_chr) {
                warnvf("Loci on sex-chromosome(s) will be skipped since\n%s allows only autosomal loci.\n", 
                        ProgName);
            } else {
                errorvf("%s allows only autosomal loci.\nNo autosomal loci to analyze.\n",
                        ProgName);
                EXIT(EARLY_TERMINATION);
            }
        } else {
            warnvf("Only markers on chromosome %d will be analyzed as sex-linked.\n",
                   SEX_CHROMOSOME);
        }
        break;
    }

/*   if (InputMode == INTERACTIVE_INPUTMODE) { */
/*     Mega2BatchItems[/ * 30 * / Xlinked_Analysis_Mode].value.option = sex_linked; */
/*     batchf(Xlinked_Analysis_Mode); */
/*   } */

    return sex_linked;
}

/* New function to return pointer to the first trait */
int *init_trp (void)

{
    int *i = &(global_trait_entries[0]);
    if (*i == -1) {
        i++;
    }
    return i;
}

#if 0
/* New locus selection menu */
static void get_map_num(int *map_num, ext_linkage_locus_top *EXLTop,
			analysis_type analysis)

{
    int m;
    char mapnums[10];

    printf("List of maps:\n");
    for (m=0; m < EXLTop->MapCnt; m++) {
        if (m == *map_num) {
            printf("*%d) %s\n", m+1, EXLTop->MapNames[m]);
        } else {
            printf(" %d) %s\n", m+1, EXLTop->MapNames[m]);
        }
    }
    printf("Enter new map number 1 - %d > ", EXLTop->MapCnt);
    fcmap(stdin, "%s", mapnums); newline;
    sscanf(mapnums, "%d", &m);
    if (m < 1 || m > EXLTop->MapCnt) {
        warn_unknown(mapnums);
    } else {
/*     if (EXLTop->map_functions[m-1] == 'h' || EXLTop->map_functions[m-1] == 'k') { */
/*       if (analysis != TO_PLINK) { */
/* 	*map_num = m-1; */
/*       } */
/*       else { */
/* 	printf("Only physical maps allowed for option %s\n", ProgName); */
/*       } */
/*     } */
/*     else if (EXLTop->map_functions[m-1] == 'p') { */
/*       if (analysis == TO_PLINK) { */
        *map_num = m-1;
/*       } */
/*       else { */
/* 	printf("Only genetic map allowed for option %s\n", ProgName); */
/*       } */
    }
    return;

}
#endif

/* Set y-chromosome positions to 0, so that output map does not have
   genetic map positions or recombination fractions */
static void reset_ychrom_positions(linkage_ped_top *Top)
{
    int has_y=-1, i, loc, upper, lower;

    for (i=0; i < main_chromocnt; i++) {
        if (global_chromo_entries[i] == MALE_CHROMOSOME) {
            has_y = i; break;
        }
    }

    if (has_y == -1  || Top->LocusTop->map_distance_type == 'p') {
        return;
    }

    printf("Now setting genetic positions of Y-linked loci to zero.\n");
    lower = ((has_y > 0)? chromo_loci_count[has_y-1] : 0);
    upper = chromo_loci_count[has_y];

    for (loc = lower; loc < upper; loc++) {
        Top->LocusTop->Marker[reordered_marker_loci[loc]].pos_avg =
            Top->LocusTop->Marker[reordered_marker_loci[loc]].pos_male =
            Top->LocusTop->Marker[reordered_marker_loci[loc]].pos_female = 0.0;
    }
}
