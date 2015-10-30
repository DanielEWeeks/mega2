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

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "makeped_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "reorder_loci_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_ghfiles_ext.h"
#include "plink_ext.h"

#include "class_old.h"

#include <string>
#include <map>

/*
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  chrom_num_to_name kosambi_to_haldane
        grow_string_ext.h:  grow
            linkage_ext.h:  count_lgenotypes get_loci_on_chromosome get_unmapped_loci switch_penetrances
            makeped_ext.h:  count_pgenotypes
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
       reorder_loci_ext.h:  get_trait_list
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line global_ou_files imax log_line
      write_ghfiles_ext.h:  write_gh_locus_file
*/



void write_quantitative_data(FILE *filep, int locusnm,
			     linkage_locus_rec *locus,
			     linkage_ped_rec *entry);
void write_affection_data(FILE *filep, int locusnm,
			  linkage_locus_rec *locus,
			  linkage_ped_rec *entry);
void write_binary_data(FILE *filep, int locusnm,
		       linkage_locus_rec *locus,
		       linkage_ped_rec *entry);

void write_simulate_numbered_data(FILE *filep, int locusnm,
                                  linkage_locus_rec *locus,
                                  linkage_ped_rec *entry);
void write_affection_tdtmax(FILE *filep, int locusnm,
                            linkage_locus_rec *locus,
                            linkage_ped_rec *entry);

int save_linkage_peds(char *pedfl, linkage_ped_top *Top,
		      analysis_type analysis, int num_loci,
		      int *loci);

void  create_linkage_files(linkage_ped_top **Top, int *numchr, char *file_names[],
			   int untyped_ped_opt);

int write_quant_stats(linkage_ped_top *Top, analysis_type analysis);


/*============= prototypes */

//
// In the past vitesse quantitative data was treated separately. However, the code was
// the same, so this routine is used to write vitesse quantitative data now.
//
// This routine is also used by makenucs.cpp:write_nuclear_peds()
//
// The analysis types that use this need to be setup to use the 'Value_Missing_Quant_On_Output'
// batch file item. To do this they should have an method defined in the analysis class
// output_quant_can_define_missing_value() which returns 'true'. If the missing output quant
// value should be numeric the analysis class method output_quant_must_be_numeric() should
// be defined returning 'true'.
//
// http://linkage.rockefeller.edu/soft/linkage/
// LINKAGE User's Guide (version 5.2) Section 2.5 'Quantitative Variables' States that...
// "Unknown phenotypes are entered as "0.0". (The code for unknown quantitative values
// is a program constant that can be changed.)"
// However, the Mega2 code (below) that writes the quantitative data uses MissingQuant and
// the outputing of "0.0" has been commented out?! The reason for this is that "0.0" is not
// a good value for a missing quant (as it may tend to overlap with real data values) and so
// people generally recompile LINKAGE to use something more reasonable.
void write_quantitative_data(FILE *filep, int locusnm,
			     linkage_locus_rec *locus,
			     linkage_ped_rec *entry)
{
    fprintf(filep, "  ");
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON) {
        if (ITEM_READ(Value_Missing_Quant_On_Output)) {
            fprintf(filep, "%s ",
		  Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        } else {
            // This should have been caught before this point, but just in case...
            errorvf("The batch file item 'Value_Missing_Quant_On_Output' must be specified because\n");
            errorvf("one or more QTL has been found to be missing in the input data.\n");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
    /* fprintf(filep, "     0.0   "); */
    } else {
        fprintf(filep, "%10.5f", entry->Pheno[locusnm].Quant);
    }
}

void write_affection_tdtmax(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    if (entry->Pheno[locusnm].Affection.Status != UNDEF)
        {
            if (entry->Pheno[locusnm].Affection.Status == 2)
                fprintf(filep,"y ");
            else fprintf(filep,"x ");
        }
    else fprintf(filep, "x ");
    if (locus->Pheno->Props.Affection.ClassCnt > 1) {
        errorvf("tdtmax doesn't handle loci with more than one liability class\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }
}

void write_affection_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    fprintf(filep, "  ");
    if (entry->Pheno[locusnm].Affection.Status != UNDEF)
        fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Status);
    else fprintf(filep, "0");
    if (locus->Pheno->Props.Affection.ClassCnt > 1) {
        fprintf(filep, " ");
        if (entry->Pheno[locusnm].Affection.Class != UNDEF)
            fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Class);
        else
            fprintf(filep, "0");
    }
}


void write_binary_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)

{
    int allele;
    fprintf(filep, "  ");
    for (allele = 0; allele < locus->AlleleCnt; allele++) {
        int a1, a2;
        get_2alleles(entry->Marker, locusnm, &a1, &a2);

        if ((allele == a1 - 1) || (allele == a2 - 1))
	  fprintf(filep, "1 ");
        else
	  fprintf(filep, "0 ");
    }
}


void write_simulate_numbered_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)

{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);

    if (a1 > 0)
        fprintf(filep, " 1 ");
    else
        fprintf(filep, " 0 ");
}


//
// Holds a cache of the string value of numerical alleles...
static std::map<int, char *> numbered_format_string_map;
/**
   @brief Format the allele appropriately.

   This routine provides a way of returning the string representation
   of the allele. If certain criteria hold then the character/string name
   is returned, otherwise a string representation of the number is returned.

   @return a string constant representation of the allele.
 */
const char *format_allele(const linkage_locus_rec *locus, const int allele)
{
    // Big 'if' predicate....
    // First see if we can/should return the character allele representation...

    // Explaining the 'force_numeric_alleles' flag...
    // For microsatalites that give you alleles as number that are not contiguous,
    // in this potentially spaces numeric case we also recode. A command line argument
    // has been added to Mega2 which says to always pick the recode (e.g., sequential numberic values).
    if ((force_numeric_alleles == 0) &&
        AnalysisOpt->allele_data_use_name_if_available() &&
        allele > 0 &&
        allele <= locus->AlleleCnt &&
        // When the input is linkage there will be no name array
        // entry for this allele.
        locus->Allele[allele-1].AlleleName != NULL &&
        // When the allele is not present for any of the markers
        // it is the case that Frequency == 0 and so we will use
        // numeric representation instead.
        locus->Allele[allele-1].Frequency != 0 &&
        // seem mrecode.cpp
        strncmp("dummy", locus->Allele[allele-1].AlleleName, 5) != 0) {

        return locus->Allele[allele-1].AlleleName;
    }

    // Here, the character allele is not used, or not available.
    //
    // At this point we look for the string representation of the numeric allele
    // in the map. If it is found, then we return it. If it is not found then
    // we add it, and then return it.
    std::map<int, char *>::iterator it = numbered_format_string_map.find(allele);
    if (it == numbered_format_string_map.end()) {
        char allele_buf[1024], *str;
        snprintf(allele_buf, 1024, "%d", allele);
        numbered_format_string_map[allele] = str = strdup(allele_buf);
        return str;
    }

    return it->second;
}

/**
 @brief Print the alleles associated with a person (entry) at a locus
*/
void write_numbered_data(FILE *filep, const int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);
    fprintf(filep, "  %2s %2s", format_allele(locus, a1), format_allele(locus, a2));
}

/* assign consecutive entry ids to entries,
   and renumber the links correctly */

typedef struct _kv {
    int key;
    int val;
} kv;

static
int LinkPerIdSort(const void *av, const void *bv)
{
    const kv *a = (const kv *)av;
    const kv *b = (const kv *)bv;
    return a->key < b->key ? -1 : (a->key == b->key ? 0 : 1) ;
}

static int LinkPerIdLookup(int k, kv *indexes, int len)
{
    kv kvk, *kvp;

    kvk.key = k;

    if (k == 0) return 0;
    kvp = (kv *)bsearch(&kvk, indexes, len, sizeof (kv), LinkPerIdSort);
    if (kvp != NULL) {
/*      printf("lookup %d,   ans = %p, k/v %d/%d\n", k, kvp, kvp->key, kvp->val);*/
        return kvp->val;
    } else {
/*      printf("lookup %d,   failed\n", k);*/
        return 0;
    }
}

static void renumber_linkage_ped(linkage_ped_top *Top)

{
    int ped, entry;
    int len;
    kv *indexes, *kvp;

    for (ped=0; ped < Top->PedCnt; ped++) {
        /* use the max entry id to assign as many elements,
           so that the entry id can be directly hashed into
           the entry number
        */
        len = Top->Ped[ped].EntryCnt;

        indexes = CALLOC((size_t) len, kv);
        /* reverse indexes */
        for (entry=0, kvp = indexes; entry < Top->Ped[ped].EntryCnt; entry++, kvp++) {
            kvp->key = Top->Ped[ped].Entry[entry].ID;
            kvp->val = entry+1;
        }
        qsort(indexes, (size_t)len, sizeof (kv), LinkPerIdSort);
            

        for(entry=0; entry < Top->Ped[ped].EntryCnt; entry++) {
            Top->Ped[ped].Entry[entry].ID = LinkPerIdLookup(Top->Ped[ped].Entry[entry].ID, indexes, len);
            Top->Ped[ped].Entry[entry].Father =
                LinkPerIdLookup(Top->Ped[ped].Entry[entry].Father, indexes, len);
            Top->Ped[ped].Entry[entry].Mother =
                LinkPerIdLookup(Top->Ped[ped].Entry[entry].Mother, indexes, len);
            Top->Ped[ped].Entry[entry].Next_PA_Sib =
                LinkPerIdLookup(Top->Ped[ped].Entry[entry].Next_PA_Sib, indexes, len);
            Top->Ped[ped].Entry[entry].Next_MA_Sib =
                LinkPerIdLookup(Top->Ped[ped].Entry[entry].Next_MA_Sib, indexes, len);
            Top->Ped[ped].Entry[entry].First_Offspring =
                LinkPerIdLookup(Top->Ped[ped].Entry[entry].First_Offspring, indexes, len);
        }

        free(indexes);
    }
    return;
}

/*
  Write linkage format pedigree file.
  For some of the analysis options, a single trait locus
  has to be written first.
*/
int save_linkage_peds(char *outfl_name, linkage_ped_top *Top,
		      analysis_type analysis, int num_loci, int *loci)
{
    int tr, nloop, num_affec=num_traits;
    int ped, entry, locus, non_numeric_peds = 0;
    int *trp, num_loc, locus_id;
    register linkage_ped_rec *Entry;
    linkage_locus_rec *Loc;
    linkage_locus_top *LTop;
    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];

    /* see if all pedigree ids are numeric */
    for (ped = 0; ped < Top->PedCnt; ped++) {
        if (atoi(Top->Ped[ped].Name) < 1) {
            non_numeric_peds = 1;
            break;
        }
    }

    NLOOP;

    LTop = Top->LocusTop;
    if (num_affec > 0) trp = &(global_trait_entries[0]);
    else {
        trp = NULL;
        if (LoopOverTrait == 1) {
            LoopOverTrait = 0;
        }
    }
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(pedfl, "%s/%s", output_paths[tr], outfl_name);
        if ((filep = fopen(pedfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_loci == 0) num_loc=NumChrLoci;
        else {
            num_loc= num_loci;
        }
        if (analysis == TO_SIMULATE) {
            fprintf(filep, "%d ", NumTypedPeds);
            for (ped = 0; ped < Top->PedCnt; ped++) {
                if (UntypedPeds != NULL) {
                    if (UntypedPeds[ped]) {
                        continue;
                    }
                }
                fprintf(filep,"%d ",Top->Ped[ped].EntryCnt);
            }
            fprintf(filep,"\n");
            renumber_linkage_ped(Top);
        }
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }

            /* find the greatest allele number for each locus */
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* write the pedigree and entry numbers */
                if (non_numeric_peds) {
                    fprintf(filep, "%6d %3d ", Top->Ped[ped].Num, Entry->ID);
                } else {
                    fprintf(filep, "%6s %3d ", Top->Ped[ped].Name, Entry->ID);
                }
                /* write the parents */
                fprintf(filep, " %3d %3d ", Entry->Father, Entry->Mother);
                /* write the first sib */
                fprintf(filep, " %3d", Entry->First_Offspring);
                /* write the next sibs */
                if (Entry->Father == 0) fprintf(filep, "   0");
                else fprintf(filep, " %3d", Entry->Next_PA_Sib);
                if (Entry->Mother == 0) fprintf(filep, "   0");
                else fprintf(filep, " %3d", Entry->Next_MA_Sib);
                /* write the sex */
                if (Entry->Sex == MALE_ID) fprintf(filep, " %d", MALE_ID);
                else if (Entry->Sex == FEMALE_ID) fprintf(filep, " %d", FEMALE_ID);
                else fprintf(filep, " 0");
                /* write the proband */
                fprintf(filep, " %d ",Entry->OrigProband);
                /* write the affection locus for some analysis types */
                if (LoopOverTrait == 1 && analysis != TO_VITESSE) {
                    switch(LTop->Locus[*trp].Type) {
                    case AFFECTION:
                        write_affection_data(filep, *trp, &(LTop->Locus[*trp]),
                                             Entry);
                        break;
                    case QUANT:
                        write_quantitative_data(filep, *trp,
                                                &(LTop->Locus[*trp]),
                                                Entry);
                        break;
                    default:
                        break;
                    }
                }
                /* Now write the genotype data */
                for (locus = 0; locus < num_loc; locus++) {
                    if (num_loci == 0) { locus_id=ChrLoci[locus]; }
                    else {
                        locus_id=loci[locus];
                        if (locus_id == -9) locus_id =  *trp;
                    }
                    Loc = &(LTOP->Locus[locus_id]);
                    switch(Loc->Type) {
                    case QUANT:
                        if (analysis == TO_VITESSE)
                            write_quantitative_data(filep, locus_id, Loc, Entry);
                        else if (LoopOverTrait == 0)
                            write_quantitative_data(filep, locus_id, Loc, Entry);
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0 || analysis == TO_VITESSE)
                            write_affection_data(filep, locus_id, Loc, Entry);
                        break;
                    case NUMBERED:
                        if (Top->pedfile_type==PREMAKEPED_PFT || analysis != TO_SIMULATE) {
                            write_numbered_data(filep, locus_id, Loc, Entry);
                        } else {
                            write_simulate_numbered_data(filep, locus_id, Loc, Entry);
                        }
                        break;
                    case BINARY:
                        write_binary_data(filep, locus_id, Loc, Entry);
                        break;
                    default:
                        fprintf(stderr, "unknown locus type (locus %d) while writing output\n",
                                locus+1);
                        break;
                    }
                }
                /* finally, write the silly ped# and per# */
                fprintf(filep, "  Ped: %s   Per: %s", Top->Ped[ped].Name, Entry->OrigID);
                if (Top->UniqueIds && analysis  == TO_LINKAGE) {
                    fprintf(filep, "  ID: %s ",  Entry->UniqueID);
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }
    return non_numeric_peds;
}


/*
 * Write the locus data file.
 */


static void linkage_file_names(char *file_names[], int glob_file)

{
    /* if there are multiple chromosomes, just put all of them inside one file */

    int select=-1;
    char stem[5], cselect[10], fl_stat[12];

    strcpy(stem, "");

    if (main_chromocnt > 1) {
        if (glob_file) {
            change_output_chr(file_names[0], 0);
            change_output_chr(file_names[1], 0);
        } else {
            change_output_chr(file_names[0], -9);
            change_output_chr(file_names[1], -9);
            strcpy(stem ,"stem");
        }
    }

    if (DEFAULT_OUTFILES) return;

    while(select != 0) {
        draw_line();
        print_outfile_mssg();
        printf("Linkage option file name selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");

        printf(" 1) Pedigree file name %s:    %-15s", stem, file_names[0]);
        if (glob_file==1 || main_chromocnt == 1)
            printf("\t%s", file_status(file_names[0], fl_stat));
        newline;

        printf(" 2) Locus file name %s:       %-15s", stem, file_names[1]);
        if (glob_file==1 || main_chromocnt == 1)
            printf("\t%s", file_status(file_names[1], fl_stat));
        newline;

        printf("Select options 0-2 > ");

        fcmap(stdin, "%s", cselect); newline;
        sscanf(cselect, "%d", &select);
        test_modified(select);

        if (select==1) {
            printf("Enter new pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]); newline;
        }
        if (select == 2) {
            printf("Enter new locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]); newline;
        }

        if (select < 0 || select > 2) {
            warn_unknown(cselect);
        }
    }
}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void  create_linkage_files(linkage_ped_top **LPedTop, int *numchr, char *file_names[],
			   int untyped_ped_opt)

{

    int i, glob_files;
    int no_orig_ped=0;
    int sex_linked, unmapped_markers_p=0;

    analysis_type analysis=TO_LINKAGE;

    if (main_chromocnt > 1)
        glob_files=global_ou_files(1);
    else
        glob_files=0;
//SL
    if (Top->LocusTop->SexLinked == 2 && glob_files == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }

    // Linkage option file name selection menu:
    // ==========================================================
    // 0) Done with this menu - please proceed
    // 1) Pedigree file name stem:    Lpedin         
    // 2) Locus file name stem:       Ldatain        
    // Select options 0-2 > 
    linkage_file_names(file_names, glob_files);

    /* for creating subsets of data */
    /*   Copy=new_lpedtop(); */
    /*   copy_linkage_ped_top(Top, Copy, 1); */

    for (i=0; i < main_chromocnt; i++) {
        if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
            unmapped_markers_p = 1;
            break;
        }
    }
    create_mssg(TO_LINKAGE);
    // Process each user selected chromosome in tern...
    for (i=0; i < main_chromocnt; i++) {
        if (!glob_files) {
            *numchr=global_chromo_entries[i];
            change_output_chr(file_names[0], *numchr);
            change_output_chr(file_names[1], *numchr);
            if (*numchr == UNKNOWN_CHROMO) {
                get_unmapped_loci(2);
            } else {
                get_loci_on_chromosome(*numchr);
            }
        } else {
            *numchr = 0;
            get_loci_on_chromosome(0);
            if (unmapped_markers_p) {
                /* append the unmapped loci */
                get_unmapped_loci(1);
            }
        }
        sex_linked =
            (((*numchr == SEX_CHROMOSOME && Top->LocusTop->SexLinked == 2) ||
              (Top->LocusTop->SexLinked == 1))? 1 : 0);

        omit_peds(untyped_ped_opt, Top);
        no_orig_ped = save_linkage_peds(file_names[0],
                                        Top, analysis, 0, NULL);
        write_gh_locus_file(file_names[1], Top->LocusTop, analysis, sex_linked);
        sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "        Locus file:           %s", file_names[1]);
        mssgf(err_msg);
        draw_line();
        if (glob_files) {
            break;
        }
    }

    if (no_orig_ped) {
        warnf("Could not use original pedigree ids in output pedigree file.");
        warnf("One or more original pedigree ids are non-numeric.");
        draw_line();
    }

/*   free_all_including_lpedtop(&Copy); */
}

#undef Top

/**
   This routine should be called after MissingQuant has been set by a call to
   user_input.cpp:set_missing_quant_input().

   This is BEFORE makeped, hence, Top->pedfile_type and pedfile_type are both
   set to the input file type, which works correctly.

   Returns true (!= 0) if a missing quantitative trait is found.
*/
int write_quant_stats(linkage_ped_top *Top,
                      analysis_type analysis)

{
    int i,j,k,l,q;
    int num_qloc=0, *qloc, *qlocp;
    double meanq, minq, maxq, var, stdev, sumq, skew, curt;
    int missingq;
    double *ped_mean, *ped_max, *ped_min, *ped_stdev;
    double *loc_mean, *loc_min, *loc_max, *loc_stdev;
    double *loc_skew, *loc_curt, dev;

    linkage_ped_rec *Entry;
    person_node_type *person;
    int ZeroQuant=0, missing_between_min_and_max, found_missing_quant=0;
    int *num_phenotyped_peds, *num_phenos, this_ped_phenotyped, *loc_missing;

    // Does this analysis type allow quantitative traits?
    if (!(analysis)->qtl_allow()) return 0;

    /* count the number of quant loci */
    for (i=0; i < num_traits; i++) {
        if (global_trait_entries[i] == -1) continue;
        if (Top->LocusTop->Locus[global_trait_entries[i]].Type == QUANT) num_qloc++;
    }

    if ((num_qloc+num_covariates) > 0) {

        num_qloc += num_covariates;
        ped_mean = CALLOC((size_t) Top->PedCnt, double);
        ped_min = CALLOC((size_t) Top->PedCnt, double);
        ped_max = CALLOC((size_t) Top->PedCnt, double);
        ped_stdev = CALLOC((size_t) Top->PedCnt, double);

        loc_mean = CALLOC((size_t) num_qloc, double);
        loc_min = CALLOC((size_t) num_qloc, double);
        loc_max = CALLOC((size_t) num_qloc, double);
        loc_stdev = CALLOC((size_t) num_qloc, double);
        loc_skew = CALLOC((size_t) num_qloc, double);
        loc_curt = CALLOC((size_t) num_qloc, double);

        num_phenotyped_peds = CALLOC((size_t) num_qloc, int);
        num_phenos = CALLOC((size_t) num_qloc, int);
        loc_missing = CALLOC((size_t) num_qloc, int);

        qloc=CALLOC((size_t)(num_traits+num_covariates), int);
        qlocp=&(qloc[0]);
        for (l=0; l < num_traits; l++) {
            *qlocp++ = global_trait_entries[l];
        }
        for (l=0; l < num_covariates; l++) {
            *qlocp++ = covariates[l];
        }

        SECTION_LOG_INIT(quant_stat);
        SECTION_LOG(quant_stat);
        mssgf("------------------------------------------------------------");
        mssgf("Per-pedigree quantitative phenotype summary:");
        mssgf("Pedigree     Mean      Std Dev     Minimum    Maximum  #Phenotypes");
        mssgf("------------------------------------------------------------");
//      SECTION_LOG_HEADER(quant_stat)
        q=0;
        for (l=0; l < num_traits+num_covariates; l++) {
            if (qloc[l] == -1) continue;

            if (Top->LocusTop->Locus[qloc[l]].Type == QUANT) {
                i = qloc[l];
                SECTION_LOG_HEADER(quant_stat);
                SECTION_LOG(quant_stat);
                mssgf(Top->LocusTop->Locus[i].LocusName);
                mssgf("------------------------------------------------------------");
                /* initialize across-ped variables */
                num_phenotyped_peds[q]=0;
                sumq=0.0; num_phenos[q]=0;
                /* minq = LARGE; maxq = -LARGE; */
                /* set to Missing for each trait locus */
                minq = MissingQuant; maxq = MissingQuant;
                missingq = 0;
                var = stdev = skew = curt = 0.0;
                for (j=0; j < Top->PedCnt; j++) {
                    /* Initialize per-per variables */
                    ped_mean[j]=0.0;
                    /* 	  ped_min[j]=LARGE; */
                    /* 	  ped_max[j]=-LARGE; */
                    ped_min[j]=ped_max[j]=MissingQuant;
                    this_ped_phenotyped=0;

                    if (Top->pedfile_type == POSTMAKEPED_PFT) {
                        /* post-makeped */
                        for (k=0; k < Top->Ped[j].EntryCnt; k++) {
                            Entry = &(Top->Ped[j].Entry[k]);

                            if (fabs(Entry->Pheno[i].Quant - MissingQuant) >  EPSILON) {
                                this_ped_phenotyped++;
                                if (MissingQuant != 0.0 && analysis == TO_LINKAGE && Entry->Pheno[i].Quant == 0.0) {
                                    ZeroQuant=1;
                                }
                                sumq += Entry->Pheno[i].Quant;
                                ped_mean[j] += Entry->Pheno[i].Quant;

                                /* Have the ped min and max been set yet?
                                   If not set to the first available phenotype in pedigree ? */
                                if (fabs(ped_max[j] - MissingQuant) <= EPSILON) {
                                    ped_max[j] = Entry->Pheno[i].Quant;
                                    ped_min[j] = Entry->Pheno[i].Quant;
                                } else {
                                    /* Now do the usual comparison */
                                    if (Entry->Pheno[i].Quant < ped_min[j]) ped_min[j] = Entry->Pheno[i].Quant;
                                    if (Entry->Pheno[i].Quant > ped_max[j]) ped_max[j] = Entry->Pheno[i].Quant;
                                }

                                /* Have the loc min and max been set yet?
                                   If not set to the first available phenotype. */
                                if (fabs(minq - MissingQuant) <= EPSILON) {
                                    minq = Entry->Pheno[i].Quant;
                                    maxq = Entry->Pheno[i].Quant;
                                } else {
                                    if (Entry->Pheno[i].Quant < minq) minq = Entry->Pheno[i].Quant;
                                    if (Entry->Pheno[i].Quant > maxq) maxq = Entry->Pheno[i].Quant;
                                }
                                /* ditto here */
                            } else {
                                // indicate that a MissingQuant was found in the input QTL data...
                                missingq++;
                            }
                        }
                    } else {
                        /* pre-makeped */
                        for (k=0; k < Top->PTop[j].num_persons; k++) {
                            person = &(Top->PTop[j].persons[k]);

                            if (fabs(person->pheno[i].Quant - MissingQuant) > EPSILON) {
                                this_ped_phenotyped++;
                                sumq += person->pheno[i].Quant;
                                ped_mean[j] += person->pheno[i].Quant;
                                ped_stdev[j] += (person->pheno[i].Quant)*(person->pheno[i].Quant);

                                /* Have the ped min and max been set yet?
                                   If not set to the first available phenotype in pedigree ? */
                                if (fabs(ped_max[j] - MissingQuant) <= EPSILON) {
                                    /* set both the values */
                                    ped_max[j] = person->pheno[i].Quant;
                                    ped_min[j] = person->pheno[i].Quant;
                                } else {
                                    if (person->pheno[i].Quant < ped_min[j]) ped_min[j] = person->pheno[i].Quant;
                                    if (person->pheno[i].Quant > ped_max[j]) ped_max[j] = person->pheno[i].Quant;
                                }

                                /* Have the loc min and max been set yet?
                                   If not set to the first available phenotype. */
                                if (fabs(minq - MissingQuant) <= EPSILON) {
                                    minq = person->pheno[i].Quant;
                                    maxq = person->pheno[i].Quant;
                                } else {
                                    /* Now do the usual comparison */
                                    if (person->pheno[i].Quant < minq) minq = person->pheno[i].Quant;
                                    if (person->pheno[i].Quant > maxq) maxq = person->pheno[i].Quant;
                                }
                            } else {
                                // indicate that a MissingQuant was found in the input QTL data...
                                missingq++;
                            }
                        }
                    }

                    if (this_ped_phenotyped == 0) {
                        ped_mean[j] = ped_max[j] = ped_min[j] = 0.0;
                    } else {
                        ped_mean[j] =
                            ped_mean[j]/(double)this_ped_phenotyped;
                        num_phenotyped_peds[q] +=
                            ((this_ped_phenotyped >0)? 1 : 0);
                        num_phenos[q] += this_ped_phenotyped;
                    }
                }

                if (num_phenos[q])
                    meanq = sumq/(double)(num_phenos[q]);
                else
                    meanq = 0.0;

                /* Second pass to set the skew and curtosis after computing mean.
                   Sums of devs raised to 2,3,4 power while looping over individuals.
                */
                SECTION_LOG_HEADER(quant_stat);
                SECTION_LOG(quant_stat);
                mssgf("------------------------------------------------------------");
                for (j=0; j < Top->PedCnt; j++) {
                    this_ped_phenotyped=0;
                    ped_stdev[j]=0.0;
                    if (Top->pedfile_type == POSTMAKEPED_PFT) {
                        /* post-makeped */
                        for (k=0; k < Top->Ped[j].EntryCnt; k++) {
                            Entry = &(Top->Ped[j].Entry[k]);
                            if (fabs(Entry->Pheno[i].Quant - MissingQuant) >  EPSILON) {
                                this_ped_phenotyped++;
                                dev = Entry->Pheno[i].Quant - meanq;
                                var += dev*dev;
                                skew += pow(dev, 3.0);
                                curt += pow(dev, 4.0);
                                ped_stdev[j] +=
                                    fabs(Entry->Pheno[i].Quant - ped_mean[j])*
                                    fabs(Entry->Pheno[i].Quant - ped_mean[j]);
                            }
                        }
                    } else {
                        for (k=0; k < Top->PTop[j].num_persons; k++) {
                            person = &(Top->PTop[j].persons[k]);
                            if (fabs(person->pheno[i].Quant - MissingQuant) > EPSILON) {
                                this_ped_phenotyped++;
                                dev = person->pheno[i].Quant - meanq;
                                var += dev*dev;
                                skew += pow(dev, 3.0);
                                curt += pow(dev, 4.0);
                                ped_stdev[j] +=
                                    fabs(person->pheno[i].Quant - ped_mean[j])*
                                    fabs(person->pheno[i].Quant - ped_mean[j]);
                            }
                        }
                    }

                    if (this_ped_phenotyped >= 2) {
                        ped_stdev[j] = sqrt(ped_stdev[j]/(double)(this_ped_phenotyped - 1));
                    } else {
                        ped_stdev[j] = 0.0;
                    }
                    /* write the quant stats for pedigree */
                    SECTION_LOG(quant_stat);
                    if (Top->pedfile_type == PREMAKEPED_PFT) {
                        msgvf("%-10d ", Top->PTop[j].ped);
                    } else {
                        msgvf("%-10s ", Top->Ped[j].Name);
                    }
                    msgvf("%10.5f %10.5f %10.5f %10.5f %8d\n",
                            ped_mean[j], ped_stdev[j], ped_min[j],
                            ped_max[j], this_ped_phenotyped);

                }

                /* Now divide each sum by the appropriate denom */
                if (num_phenos[q] > 1) {
                    var /= (double)(num_phenos[q] - 1);
                    stdev = sqrt(var);
                    skew /= (pow(stdev, 3.0)*((double)(num_phenos[q] - 1)));
                    curt /= (pow(stdev, 4.0)*((double)(num_phenos[q] - 1)));
                    curt -= 3.0;
                } else {
                    var=stdev = skew = curt = 0.0;
                }
                loc_mean[q] = meanq;
                loc_stdev[q] = stdev;
                loc_skew[q] = skew;
                loc_curt[q] = curt;
                loc_min[q] = minq;
                loc_max[q] = maxq;
                loc_missing[q] = missingq;
                // flag returned by this routine indicating that a missing QTL value was found...
                if (missingq > 0) found_missing_quant = 1;
                q++;
            }
        }
        SECTION_LOG_FINI(quant_stat);

        /* Now display combined statistics on screen as well as log */
        SECTION_LOG_INIT(quant_stat_sum);
        SECTION_LOG(quant_stat_sum);
        mssgf("               Quantitative trait phenotype statistics");
        mssgf("-------------------------------------------------------------------------------");
        mssgf("QTL            Missing  Minimum   Maximum   Total      Pedigrees    Total");
        mssgf("                                            pedigrees  phenotyped   phenotypes");
        mssgf("-------------------------------------------------------------------------------");
        SECTION_LOG_HEADER(quant_stat_sum);

        #define FLAG_USED_TO_INDICATE_MISSING_BETWEEN_MIN_AND_MAX                 '*'
        missing_between_min_and_max = 0;
        for (q=0,l=0; l < num_traits+num_covariates; l++) {
            if (qloc[l] == -1) continue;
            
            if (Top->LocusTop->Locus[qloc[l]].Type == QUANT) {
                char missing_between_min_and_max_flag;
                i = qloc[l];
                
                missing_between_min_and_max_flag = ' ';
                if (loc_missing[q] > 0 && loc_min[q] <= MissingQuant && MissingQuant <= loc_max[q]) {
                    missing_between_min_and_max_flag = FLAG_USED_TO_INDICATE_MISSING_BETWEEN_MIN_AND_MAX;
                    missing_between_min_and_max = 1;
                }
                SECTION_LOG(quant_stat_sum);
                mssgvf("%-15s %6d%c %7.3f  %8.3f   %9d   %9d    %9d\n",
                       Top->LocusTop->Locus[i].LocusName,
                       loc_missing[q], missing_between_min_and_max_flag,
                       loc_min[q], loc_max[q],
                       Top->PedCnt, num_phenotyped_peds[q], num_phenos[q]);
                q++;
            }
        }
        if (missing_between_min_and_max != 0) {
            mssgvf("Missing value counts marked with a '%c' have a value that falls between\n",
                   FLAG_USED_TO_INDICATE_MISSING_BETWEEN_MIN_AND_MAX);
            mssgvf("the Minimum and Maximum.\n");
        }
	if (fabs(MissingQuant - QMISSING) <= EPSILON)
            mssgvf("NOTE: The Missing QTL value on input has been assigned as 'NA'.\n");
	else
            mssgvf("NOTE: The Missing QTL value on input has been assigned as '%4.2f'.\n",  MissingQuant);

        SECTION_LOG_HEADER(quant_stat_sum)
        SECTION_LOG(quant_stat_sum)
        mssgf("-------------------------------------------------------------------------");
        mssgf("QTL                 Mean   Std Dev  Skewness  Kurtosis");
        mssgf("-------------------------------------------------------------------------");
        SECTION_LOG_HEADER(quant_stat_sum);
        q=0;
        for (l=0; l < num_traits+num_covariates; l++) {
            if (qloc[l] == -1) {
                continue;
            }
            if (Top->LocusTop->Locus[qloc[l]].Type == QUANT) {
                i = qloc[l];
                sprintf(err_msg,
                        "%-15s %8.3f  %8.3f  %8.3f  %8.3f",
                        Top->LocusTop->Locus[i].LocusName,
                        loc_mean[q], loc_stdev[q], loc_skew[q], loc_curt[q]);
                SECTION_LOG(quant_stat_sum);
                mssgf(err_msg);
                q++;
            }
        }

        log_line(mssgf);
        SECTION_LOG_FINI(quant_stat_sum);
        if (ZeroQuant) {
            warnf("Found KNOWN quantitative phenotypes with value 0.0");
            warnf("Linkage-format UNKNOWN quantitative phenotype is assigned 0.0");
            log_line(mssgf);
        }

        free(ped_mean);
        free(ped_min);
        free(ped_max);
        free(ped_stdev);
        free(loc_mean);
        free(loc_min);
        free(loc_max);
        free(loc_stdev);
        free(num_phenotyped_peds);
        free(num_phenos);
        free(loc_missing);
        free(qloc);
    }
    return found_missing_quant;
}

void write_locus_stats(linkage_locus_top *LTop, file_format locus_file_type)

{

    int i, k, marker_count;
    char chrom_name[4];
    int curr_min_chrom, curr_min_chrom_index;
    int *display_chrom = NumChromo ? CALLOC((size_t) NumChromo, int) : NULL;

    /* display_chrom will hold indexes of the chromosome
       such that they are displayed in order */
    if (NumChromo) {
        curr_min_chrom=MAX_CHROMO_NUM;
        curr_min_chrom_index=0; //compiler: if loop fails to assign
        for (i=0; i < NumChromo; i++) {
            if (global_chromo_entries[i] < curr_min_chrom) {
                curr_min_chrom = global_chromo_entries[i];
                curr_min_chrom_index = i;
            }
        }
        display_chrom[0] = curr_min_chrom_index;
    }
    for (k=1; k < NumChromo; k++) {
        curr_min_chrom=MAX_CHROMO_NUM + UNKNOWN_CHROMO;
        for (i=0; i < NumChromo; i++) {
            if (global_chromo_entries[i] < curr_min_chrom &&
                global_chromo_entries[i] > global_chromo_entries[display_chrom[k-1]]) {
                curr_min_chrom = global_chromo_entries[i];
                curr_min_chrom_index = i;
            }
        }
        display_chrom[k] = curr_min_chrom_index;
    }

    marker_count=0;

#ifndef HIDESTATUS
    log_line(mssgf);
    switch (locus_file_type) {
    case NAMES:
        mssgf("Locus file is in NAMES format.");
        break;

    case LINKAGE:
        mssgf("Locus file is in LINKAGE format.");
        break;

    default:
        break;
    }
#endif

    mssgvf("Total number of loci =  %d\n", LTop->LocusCnt);
    get_trait_list(LTop, 1);
    marker_count = LTop->LocusCnt - num_affection - num_quantitative;
    mssgvf("      %d Marker %s\n", marker_count,
           ((marker_count > 1)? "loci " : "locus "));

    if (NumUnmapped > 1) {
        mssgvf("      %d Unmapped loci\n", NumUnmapped);
    } else if (NumUnmapped == 1) {
        mssgf("      1 Unmapped locus");
    }

    mssgf("Number of loci found per chromosome (chromosome:number)");

    strcpy(err_msg, "   ");
    for (k=0; k < NumChromo - 1; k++) {
        i = display_chrom[k];
        if (strlen(err_msg) >= 60) {
            mssgf(err_msg);
            strcpy(err_msg, "   ");
        }
        chrom_num_to_name(global_chromo_entries[i], &(chrom_name[0]));
        grow(err_msg, "%s:%d, ",
             chrom_name,
             ((i==0)? chromo_loci_count[i] :
              chromo_loci_count[i] - chromo_loci_count[i-1]));
    }

    if (NumChromo) {
        i = display_chrom[NumChromo - 1];
        chrom_num_to_name(global_chromo_entries[i],
                          &(chrom_name[0]));
        if (i==0) {
            grow(err_msg, "%s:%d",
                 chrom_name,  chromo_loci_count[i]);
        } else {
            grow(err_msg, "%s:%d",
                 chrom_name,  chromo_loci_count[i]  - chromo_loci_count[i - 1]);
        }
    } else {
        grow(err_msg, "no markers specified");
    }
    mssgf(err_msg);
    strcpy(err_msg, "");
    free(display_chrom);
    log_line(mssgf);
}

void write_ped_stats(linkage_ped_top *Top)

{
    size_t individual_count, male_count, female_count;
    size_t typed, untyped, male_typed, female_typed, half_typed;
    size_t peds_typed;
/* void (*log_func)() = mssgf; */

    if (Mega2Status > LOCI_REORDERED && num_reordered == 0) {
        return;
    }

    if (Top->pedfile_type == PREMAKEPED_PFT) {
        count_pgenotypes(Top, &individual_count, &male_count, &female_count,
                         &untyped, &typed, &peds_typed, &male_typed, &female_typed,
                         &half_typed);
#ifndef HIDESTATUS
        if (PLINK.plink == PED_format || PLINK.plink == binary_PED_format) {
            mssgf("Input pedigree file is in PLINK-fam format. ");
        } else {
            mssgf("Input pedigree file is in pre-makeped format. ");
        }
#endif /* HIDESTATUS */
    } else {
        count_lgenotypes(Top, &individual_count, &male_count, &female_count,
                         &untyped, &typed, &peds_typed, &male_typed, &female_typed,
                         &half_typed);
#ifndef HIDESTATUS
        if (Top->pedfile_type == POSTMAKEPED_PFT) {
            mssgf("Input pedigree file is in post-makeped format.");
        }
#endif /* HIDESTATUS */
    }

    /* Leave the first part as it is */
    mssgf("                                                Marker Genotypes");
    mssgf("                                                Fully    Half");
    mssgf("     Pedigrees   People   Males   Females       Typed    Typed     Total");

    sprintf(err_msg,
#ifdef MS_PRINTF
            "TOTAL   %6d %8Iu %7Iu   %7Iu    %8Iu %8Iu %9Iu",
#else
            "TOTAL   %6d %8zu %7zu   %7zu    %8zu %8zu %9zu",
#endif
            Top->PedCnt, individual_count, male_count, female_count,
            typed, half_typed, typed+untyped+half_typed);
    mssgf(err_msg);
    sprintf(err_msg,
#ifdef MS_PRINTF
            "Typed   %6Iu %8Iu %7Iu   %7Iu",
#else
            "Typed   %6zu %8zu %7zu   %7zu",
#endif
            peds_typed, male_typed+female_typed, male_typed, female_typed);
    mssgf(err_msg);
    sprintf(err_msg,
#ifdef MS_PRINTF
            "Untyped %6Iu %8Iu %7Iu   %7Iu",
#else
            "Untyped %6zu %8zu %7zu   %7zu",
#endif
            Top->PedCnt - peds_typed,
            individual_count - male_typed - female_typed,
            male_count - male_typed,
            female_count - female_typed);
    mssgf(err_msg);
    log_line(mssgf);
}

/* This is the same as write_gh_locus_file, except it takes
   a vector of selected loci.
   Used by VITESSE and SUP.
   For VITESSE, trait position is fixed for all trait loci, and
   indicated by a -9 inside the locus list <loci>.
   For SUP, trait needs to be inserted into the correct position as
   user indicated. Positions for each trait are stored in SLINK
   params.  Recombination fractions are replaced by intermarker
   distances in Haldane cM.

   Tried to avoid creating two almost identical routines for the
   two analysis options.
*/
static void write_linkage_locfile_inorder_sex_averaged(linkage_locus_top *LTop,
						       char *loutfl_name,
						       int locus_num,
						       int *loci, int sex_linked,
						       analysis_type analysis)
{
    int tr, nloop, num_affec=num_traits, locus, allele, tmpi, tmpi2;
    int *trp, trnum, num_loci = 0,  *locus_order = NULL, trait_locus = 0; //compiler: too hard
    char fl[2*FILENAME_LENGTH];
    register linkage_locus_rec *Locus;
    FILE *filep;

    double *positions = NULL; //compiler: too hard

    if (analysis == TO_VITESSE) {
        locus_order = &(loci[0]);
        num_loci = locus_num;
        positions=NULL;
    } else if (analysis == TO_SUP) {
        if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
            locus_order = CALLOC((size_t) locus_num+1, int);
            positions = CALLOC((size_t) locus_num+1 , double);
        } else {
            locus_order = CALLOC((size_t) locus_num, int);
            positions = CALLOC((size_t) locus_num , double);
            LTop->Run.slink.trait_positions[0]= -99;
        }
    }

    NLOOP;

    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        trp=&(global_trait_entries[0]);
    } else {
        trp=NULL;
    }
    trnum=0;

    for (tr=0; tr<= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(fl, "%s/%s", output_paths[tr], loutfl_name);
        if ((filep=fopen(fl, "w")) == NULL) {
            errorvf("Could not open %s for writing.\n", fl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (analysis == TO_SUP) {
            /* create the locus_order from loci, inserting -9 only if this
               trait is to be simulated on specified chromosome */
            tmpi=0; tmpi2=0;
            num_loci = locus_num+1;
	    // What about SEX_SPECIFIC???
            while((tmpi < locus_num) &&
                  (LTop->Marker[loci[tmpi]].pos_avg <
                   LTop->Run.slink.trait_positions[trnum])) {
                locus_order[tmpi2] = loci[tmpi];
                positions[tmpi2] = LTop->Marker[loci[tmpi]].pos_avg;
                tmpi++; tmpi2++;
            }
            if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
                /* insert the trait locus indicator, and the position */
                locus_order[tmpi2]=-9;
                positions[tmpi2] = LTop->Run.slink.trait_positions[trnum];
                tmpi2++;
            }
            /* copy over the rest of the marker loci  and their positions */
            while(tmpi < locus_num) {
                locus_order[tmpi2] = loci[tmpi];
                positions[tmpi2] = LTop->Marker[loci[tmpi]].pos_avg;
                tmpi++; tmpi2++;
            }
        }

        fprintf(filep, "%d %d %d %d\n%d %.3f %.3f %d\n",
                num_loci, LTop->RiskLocus, sex_linked,
                LTop->Program, LTop->MutLocus, LTop->MutMale,
                LTop->MutFemale, LTop->Haplotype);

        for (tmpi=0; tmpi < num_loci; tmpi++) {
            /* write the order in the locus file */
            fprintf(filep, "%d ", tmpi+1);
        }
        fputc('\n', filep);

        if (analysis == TO_VITESSE || analysis == TO_SUP) // hopefully to make compiler happy
        for (locus = 0; locus < num_loci; locus++) {
            if (locus_order[locus] == -9) {
                Locus=&(LTop->Locus[*trp]);
                trait_locus=locus;
            } else  {
                Locus=&(LTop->Locus[locus_order[locus]]);
            }
            fprintf(filep, "%d %d", Locus->Type - 1, Locus->AlleleCnt);
            if (Locus->LocusName != NULL)
                fprintf(filep, " #%s\n", Locus->LocusName);
            for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
            }
            fputc('\n', filep);
            switch(Locus->Type) {
            case QUANT:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Quant.ClassCnt);
                for (tmpi = Locus->Pheno->Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                    for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                        fprintf(filep," %.6f", Locus->Pheno->Props.Quant.Mean[tmpi-1][tmpi2]);
                fprintf(filep, " << GENOTYPE MEANS\n");
                /* ASSUMES only one trait per qtl */
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Variance[0][0]);
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Multiplier);
                break;
            case AFFECTION:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < Locus->Pheno->Props.Affection.ClassCnt; tmpi++) {
                    if (sex_linked) {
                        for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                        fputc('\n', filep);
                        for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                        fputc('\n', filep);
                    } else {
                        for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                        fputc('\n', filep);
                    }
                }
                break;
            case BINARY:
                fprintf(filep, "%d\n", Locus->Marker->Props.Binary.FactorCnt);
                for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++) {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(filep, " %d", (int) Locus->Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', filep);
                }
                break;
            case NUMBERED:
                /* nothing to do  */
                break;
            default:
	        errorvf("Unknown locus type %d\n", (int) LTop->Locus[locus].Type);
                EXIT(DATA_TYPE_ERROR);
                break;
            }
        }

        fprintf(filep, "%d %d\n", LTop->SexDiff, LTop->Interference);
        /* now for the recombination information */
        if (analysis == TO_VITESSE) {
            if (LTop->Program == MLINK) {
                fprintf(filep, "%5.4f\n", LTop->Run.mlink.start_theta);
                fprintf(filep, "1 %5.4f %5.4f\n", LTop->Run.mlink.increment,
                        LTop->Run.mlink.stop_theta);
            } else {
                for (tmpi=0; tmpi < num_loci-1; tmpi++)
                    fprintf(filep, "%9.8f ", LTop->Run.linkmap.recomb_frac[tmpi]);
                fprintf(filep, "\n");
                fprintf(filep, "%d %9.8f %d \n", trait_locus+1,
                        LTop->Run.linkmap.stop_theta,
                        LTop->Run.linkmap.num_evals);
            }
        } else if (analysis == TO_SUP) {
            /* From the SUP documentation:
             * "The intermarker distances must be entered in units of centiMorgans,
             * according to the Haldane map (SUP will not simulate interference)."
             */
            for (tmpi = 0; tmpi < num_loci - 1; tmpi++) {
                if (LTop->map_distance_type == 'k') {
                    fprintf(filep, "%7.4f ",
			    kosambi_to_haldane(positions[tmpi + 1] - positions[tmpi]));
                } else {
                    fprintf(filep, "%7.4f ", positions[tmpi + 1] - positions[tmpi]);
                }
            }
        }
        fprintf(filep, "\n");
        fprintf(filep,"1 0.10000 0.45000\n");
        fclose(filep);
        trnum++;
        trp++;
        if (nloop == 1) break;
    }

    if (analysis == TO_SUP) {
        free(positions); free(locus_order);
    }
}

static void write_linkage_locfile_inorder_sex_specific(linkage_locus_top *LTop,
						       char *loutfl_name,
						       int locus_num,
						       int *loci, int sex_linked,
						       analysis_type analysis)
{
    int tr, nloop, num_affec=num_traits, locus, allele, tmpi, tmpi2;
    int *trp, trnum, num_loci = 0,  *locus_order = NULL, trait_locus = 0; //compiler: too hard
    char fl[2*FILENAME_LENGTH];
    register linkage_locus_rec *Locus;
    FILE *filep;

    double *positions_male = NULL, *positions_female = NULL; //compiler: too hard

    if (analysis == TO_VITESSE) {
        locus_order = &(loci[0]);
        num_loci = locus_num;
        positions_male = positions_female = NULL;
    } else if (analysis == TO_SUP) {
        if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
            locus_order = CALLOC((size_t) locus_num+1, int);
            positions_male = CALLOC((size_t) locus_num+1 , double);
            positions_female = CALLOC((size_t) locus_num+1 , double);
        } else {
            locus_order = CALLOC((size_t) locus_num, int);
            positions_male = CALLOC((size_t) locus_num , double);
            positions_female = CALLOC((size_t) locus_num , double);
            LTop->Run.slink.trait_positions[0]= -99;
        }
    }

    NLOOP;

    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        trp=&(global_trait_entries[0]);
    } else {
        trp=NULL;
    }
    trnum=0;

    for (tr=0; tr<= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(fl, "%s/%s", output_paths[tr], loutfl_name);
        if ((filep=fopen(fl, "w")) == NULL) {
            errorvf("Could not open %s for writing.\n", fl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (analysis == TO_SUP) {
            /* create the locus_order from loci, inserting -9 only if this
               trait is to be simulated on specified chromosome */
            tmpi=0; tmpi2=0;
            num_loci = locus_num+1;

            while((tmpi < locus_num) &&
                  (LTop->Marker[loci[tmpi]].pos_avg <
                   LTop->Run.slink.trait_positions[trnum])) {
                locus_order[tmpi2] = loci[tmpi];
                positions_male[tmpi2] = LTop->Marker[loci[tmpi]].pos_male;
                positions_female[tmpi2] = LTop->Marker[loci[tmpi]].pos_female;
                tmpi++; tmpi2++;
            }
            if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
                /* insert the trait locus indicator, and the position */
                locus_order[tmpi2]=-9;
                positions_male[tmpi2] = LTop->Run.slink.trait_positions[trnum];
                positions_female[tmpi2] = LTop->Run.slink.trait_positions[trnum];
                tmpi2++;
            }
            /* copy over the rest of the marker loci  and their positions */
            while(tmpi < locus_num) {
                locus_order[tmpi2] = loci[tmpi];
                positions_male[tmpi2] = LTop->Marker[loci[tmpi]].pos_male;
                positions_female[tmpi2] = LTop->Marker[loci[tmpi]].pos_female;
                tmpi++; tmpi2++;
            }
        }

        fprintf(filep, "%d %d %d %d\n%d %.3f %.3f %d\n",
                num_loci, LTop->RiskLocus, sex_linked,
                LTop->Program, LTop->MutLocus, LTop->MutMale,
                LTop->MutFemale, LTop->Haplotype);

        for (tmpi=0; tmpi < num_loci; tmpi++) {
            /* write the order in the locus file */
            fprintf(filep, "%d ", tmpi+1);
        }
        fputc('\n', filep);

        if (analysis == TO_VITESSE || analysis == TO_SUP) // hopefully to make compiler happy
        for (locus = 0; locus < num_loci; locus++) {
            if (locus_order[locus] == -9) {
                Locus=&(LTop->Locus[*trp]);
                trait_locus=locus;
            } else  {
                Locus=&(LTop->Locus[locus_order[locus]]);
            }
            fprintf(filep, "%d %d", Locus->Type - 1, Locus->AlleleCnt);
            if (Locus->LocusName != NULL)
                fprintf(filep, " #%s\n", Locus->LocusName);
            for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
            }
            fputc('\n', filep);
            switch(Locus->Type) {
            case QUANT:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Quant.ClassCnt);
                for (tmpi = Locus->Pheno->Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                    for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                        fprintf(filep," %.6f", Locus->Pheno->Props.Quant.Mean[tmpi-1][tmpi2]);
                fprintf(filep, " << GENOTYPE MEANS\n");
                /* ASSUMES only one trait per qtl */
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Variance[0][0]);
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Multiplier);
                break;
            case AFFECTION:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < Locus->Pheno->Props.Affection.ClassCnt; tmpi++) {
                    if (sex_linked) {
                        for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                        fputc('\n', filep);
                        for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                        fputc('\n', filep);
                    } else {
                        for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                        fputc('\n', filep);
                    }
                }
                break;
            case BINARY:
                fprintf(filep, "%d\n", Locus->Marker->Props.Binary.FactorCnt);
                for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++) {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(filep, " %d", (int) Locus->Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', filep);
                }
                break;
            case NUMBERED:
                /* nothing to do  */
                break;
            default:
                errorvf("Unknown locus type %d\n", (int) LTop->Locus[locus].Type);
                EXIT(DATA_TYPE_ERROR);
                break;
            }
        }

        fprintf(filep, "%d %d\n", LTop->SexDiff, LTop->Interference);
        /* now for the recombination information */
        if (analysis == TO_VITESSE) {
            if (LTop->Program == MLINK) {
                fprintf(filep, "%5.4f\n",
			LTop->Run.mlink.start_theta);
                fprintf(filep, "1 %5.4f %5.4f\n",
			LTop->Run.mlink.increment,
                        LTop->Run.mlink.stop_theta);
            } else {
                for (tmpi=0; tmpi < num_loci-1; tmpi++)
                    fprintf(filep, "%9.8f ", LTop->Run.linkmap.recomb_frac[tmpi]);
                fprintf(filep, "\n");
                fprintf(filep, "%d %9.8f %d \n",
			trait_locus+1,
                        LTop->Run.linkmap.stop_theta,
                        LTop->Run.linkmap.num_evals);
            }
        } else if (analysis == TO_SUP) {
            /* From the SUP documentation:
             * "The intermarker distances must be entered in units of centiMorgans,
             * according to the Haldane map (SUP will not simulate interference)."
             */
            for (tmpi = 0; tmpi < num_loci - 1; tmpi++) {
                if (LTop->map_distance_type == 'k') {
                    fprintf(filep, "%7.4f ",
			    kosambi_to_haldane(positions_male[tmpi + 1] - positions_male[tmpi]));
                } else {
                    fprintf(filep, "%7.4f ",
			    positions_male[tmpi + 1] - positions_male[tmpi]);
                }
            }
	    fprintf(filep, "\n");
            for (tmpi = 0; tmpi < num_loci - 1; tmpi++) {
                if (LTop->map_distance_type == 'k') {
                    fprintf(filep, "%7.4f ",
			    kosambi_to_haldane(positions_female[tmpi + 1] - positions_female[tmpi]));
                } else {
                    fprintf(filep, "%7.4f ",
			    positions_female[tmpi + 1] - positions_female[tmpi]);
                }
            }
	    fprintf(filep, "\n");
        }
        fprintf(filep,"1 0.10000 0.45000\n");
        fclose(filep);
        trnum++;
        trp++;
        if (nloop == 1) break;
    }

    if (analysis == TO_SUP) {
        free(positions_male); free(positions_female);
	free(locus_order);
    }
}

void write_linkage_locfile_inorder(linkage_locus_top *LTop,
				   char *loutfl_name,
				   int locus_num,
				   int *loci, int sex_linked,
				   analysis_type analysis)
{
    if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
        write_linkage_locfile_inorder_sex_averaged(LTop, loutfl_name, locus_num, loci, sex_linked, analysis);
    } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
        write_linkage_locfile_inorder_sex_specific(LTop, loutfl_name, locus_num, loci, sex_linked, analysis);
    }
}
