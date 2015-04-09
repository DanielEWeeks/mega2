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
#include <ctype.h>

#include "common.h"
#include "typedefs.h"
#include "tod.hh"

#include "mrecode.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "makeped_ext.h"
#include "marker_lookup_ext.h"
#include "omit_ped_ext.h"
#include "read_files_ext.h"
#include "select_individuals_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "vcftools/mega2_vcftools_interface.h"

/*
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
        grow_string_ext.h:  grow
            linkage_ext.h:  clear_laffclass new_llocustop
            makeped_ext.h:  check_pre_makeped read_pre_makeped
      marker_lookup_ext.h:  allele_prop_reset search_allele
           omit_ped_ext.h:  omit_peds
         read_files_ext.h:  premakeped_omit_file read_linkage_ped_file
 select_individuals_ext.h:  founders_or_everyone get_count_option random_ped_member
              utils_ext.h:  EXIT draw_line log_line summary_time_stamp
        write_files_ext.h:  write_ped_stats
*/


linkage_locus_top *read_common_marker_data(int all_loci, int num_markers, char **names, char *types, int annotated, double penetrances_read, double freq, double pen[5]);

void recode_liability_class(pheno_type *pheno_list, linkage_locus_top *LTop, int locus);

/* function for inserting a raw allele into a marker allele list,
   three cases: new allele_list (first non-zero allele)
   new allele
   existing allele
*/

void free_marker_item(allele_list_type *marker_item)
{
    if (marker_item == NULL) {
        return;
    } else if (marker_item->next == NULL) {
        free(marker_item);
        marker_item = NULL;
        return;
    } else if (marker_item->next != NULL) {
        free_marker_item(marker_item->next);
        free(marker_item);
    }
    return;
}

void free_marker_item_class(class_list_type *liability_class)
{
    if (liability_class == NULL) {
        return;
    } else if (liability_class->next == NULL) {
        free(liability_class);
        liability_class = NULL;
        return;
    } else if (liability_class->next != NULL) {
        free_marker_item_class(liability_class->next);
        free(liability_class);
    }
    return;
}

/* traverse the list until we find an allele value
   that is greater than or equal to the value we
   are inserting. If we find a value greater than the
   insert-value, then we have to insert this allele between
   the last two members.

*/

/* NM: 9/17/04
   count_allele is a flag, skip=0, count=1
*/

/* function for inserting a raw allele into a marker allele list,
   three cases: new allele_list (first non-zero allele)
   new allele
   existing allele
*/

static void count_this_allele(allele_list_type *marker_item,
                              const char *all_val, int opt)

{
    allele_list_type *marker_itemp = marker_item;

    if (allelecmp(all_val, REC_UNKNOWN) == 0) {
        return;
    }

    marker_itemp = (allele_list_type *)allele2allele_prop_prop(all_val);
/*
    while(allelecmp(marker_itemp->allele_freq.name,
                 all_val)) {
        marker_itemp = marker_itemp->next;
    }
*/
    if (opt != 4) {
        if (opt == 1) {
            (marker_itemp->allele_freq.unique_count)++;
            (marker_itemp->allele_freq.random_count)++;
            (marker_itemp->allele_freq.founder_count)++;
        }

        if (opt == 2) {
            (marker_itemp->allele_freq.random_count)++;
            (marker_itemp->allele_freq.unique_count)++;
        }

        if (opt == 3) {
            (marker_itemp->allele_freq.unique_count)++;
        }
    } else {
/*
 * handled elsewhere
 *        (marker_itemp->allele_freq.everyone_count)++;
 */
    }
    return;

}
/*
  case 1:
            (marker_itemp->allele_freq.founder_count)++;
  case 2:
            (marker_itemp->allele_freq.random_count)++;
  case 3:
            (marker_itemp->allele_freq.unique_count)++;
            break
  case 4:
           (marker_itemp->allele_freq.everyone_count)++;
 */
static allele_list_type *insert_into_allele_list(allele_list_type **marker_item,
				   const char *all_val)

{
    allele_list_type *marker_itemp = *marker_item;
    allele_list_type *marker_itemp1 = *marker_item;
    allele_list_type *new_entry;

    /* all_val comes from {P,}RAllele ... already canonical */
    if (allelecmp(all_val, REC_UNKNOWN) == 0) {
        return (allele_list_type *) 0;
    }

    if (marker_itemp == NULL) {
        marker_itemp = CALLOC((size_t) 1, allele_list_type);
        marker_itemp->next = NULL;
        marker_itemp->allele_freq.count =
            marker_itemp->allele_freq.founder_count =
            marker_itemp->allele_freq.random_count =
            marker_itemp->allele_freq.unique_count =
            marker_itemp->allele_freq.everyone_count = 0;
        marker_itemp->allele_freq.name = all_val;
        marker_itemp->allele_freq.index = 1;
        *marker_item = marker_itemp;
        new_entry = marker_itemp;
    } else {
        while((marker_itemp != NULL) &&
              (strcmp(marker_itemp->allele_freq.name, all_val) <= 0)) {
            marker_itemp1 = marker_itemp;
            marker_itemp = marker_itemp->next;
        }

        if (allelecmp(marker_itemp1->allele_freq.name, all_val)==0) {
            new_entry = marker_itemp1;
            /* do nothing */
        } else {
            new_entry = CALLOC((size_t) 1, allele_list_type);
            new_entry->allele_freq.name = all_val;
            new_entry->allele_freq.count =
                new_entry->allele_freq.founder_count =
                new_entry->allele_freq.random_count =
                new_entry->allele_freq.unique_count =
                new_entry->allele_freq.everyone_count = 0;

            if (marker_itemp == NULL) {
                /* Need to add a new at the end */
                marker_itemp1->next = new_entry;
                new_entry->next=NULL;
                new_entry->allele_freq.index =
                    marker_itemp1->allele_freq.index+1;
                marker_itemp=new_entry->next;
            } else if (marker_itemp == marker_itemp1) {
                /* insert at the beginning of the list */
                new_entry->next=marker_itemp1;
                new_entry->allele_freq.index = 1;
                *marker_item=new_entry;
            } else {
                /* Need to insert between itemp and item1 */
                new_entry->next = marker_itemp1->next;
                new_entry->allele_freq.index =
                    marker_itemp1->allele_freq.index + 1;
                marker_itemp1->next = new_entry;
            }

            while(marker_itemp != NULL) {
                (marker_itemp->allele_freq.index)++;
                marker_itemp = marker_itemp->next;
            }
        }
    }
    return new_entry;
}


/*-------------------end of insert_into_allele_list------------------*/


static void insert_into_class_list(class_list_type **class_list,
				   int class_num)

{
    class_list_type *class_listp = *class_list;
    class_list_type *new_entry = CALLOC((size_t) 1, class_list_type);

    new_entry->next = NULL;
    new_entry->l_class.num = class_num;
    new_entry->l_class.autosomal_pen = CALLOC((size_t) 3, double);

    new_entry->l_class.autosomal_pen[0] = DOM_G11;
    new_entry->l_class.autosomal_pen[1] = DOM_G12;
    new_entry->l_class.autosomal_pen[2] = DOM_G22;
    if (class_listp == NULL) {
        class_listp = new_entry;
        *class_list = class_listp;
    } else {
        while (class_listp->next != NULL) {
            class_listp = class_listp->next;
        }
        class_listp->next = new_entry;
    }
    return;
}

#ifndef num_alleles
#define num_alleles data.num_alleles
#endif

#ifndef num_classes
#define num_classes data.num_classes
#endif

void count_classes(pheno_type *pheno_list, linkage_locus_top *LTop)

{
    int m, count;
    class_list_type *classp, *classlp;

    for (m=0; m < LTop->LocusCnt; m++) {
        if (LTop->Locus[m].Type == AFFECTION) {
            classp  = pheno_list[m].class_list;
            if (classp == NULL) continue;
            classlp = NULL;
            count=0;
            while(classp != NULL) {
                count++;
                classlp = classp;
                classp  = classp->next;
            }
            count = classlp == NULL ? 0 : classlp->l_class.num;
            pheno_list[m].num_classes = count;
            LTop->Pheno[m].Props.Affection.ClassCnt = count;
        }
    }
}

void load_classes(pheno_type *pheno, linkage_locus_top *LTop)
{
    int m, count = 0, error = 0;
    class_list_type *classp, *classlp;

    for (m=0; m < LTop->PhenoCnt; m++) {
        if (LTop->Locus[m].Type == AFFECTION) {
            classp  = pheno[m].class_list;
            if (classp == NULL) {
                if (LTop->Pheno[m].Props.Affection.ClassCnt > 1) {
                    errorvf("%s has multiple liability classes, but no defined penetrance.\n Please use a penetrance file to specify them.\n",
                            LTop->Locus[m].Name);
                    error++;
                    continue;
                } else if (LTop->Pheno[m].Props.Affection.ClassCnt == 1) {
                    msgvf("Trait '%s' will be assigned the default penetrance: (%.4f %.4f %.4f)\n",
                          LTop->Pheno[m].Name, DOM_G11, DOM_G12, DOM_G22);
                    count = 1;
                }
            } else {
                classlp = NULL;
                count=0;
                while(classp != NULL) {
                    count++;
                    classlp = classp;
                    classp  = classp->next;
                }
                count = classlp == NULL ? 0 : classlp->l_class.num;
            }
            pheno[m].num_classes = count;
            LTop->Pheno[m].Props.Affection.ClassCnt = count;
            recode_liability_class(pheno, LTop, m);
        }
    }

    if (error) {
        EXIT(INPUT_DATA_ERROR);
    }

    return;
}

/* Function to traverse an allele_list and just count the
   length of the list */

void count_raw_alleles(marker_type *markers, linkage_locus_top *LTop)

{
    int m, count;
    allele_list_type *allelep;

    for (m = LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
/*     if (markers[m].estimate_frequencies == 0) { */
/*       continue; */
/*     } */
        allelep = markers[m].first_allele;
        count=0;
        while(allelep != NULL) {
            count++;
            allelep = allelep->next;
        }
        markers[m].num_alleles = count;
    }
    return;
}

/*-------------------end of count_alleles------------------*/
/*---- function sort_allele_list to sort allele names
  if these do not need to be recoded */

struct kv {
    int key;
    allele_list_type *val;
} kvstore[100];

static int kvsort(const void *va, const void *vb)
{
    struct kv *a = (struct kv *)va;
    struct kv *b = (struct kv *)vb;
    return a->key - b->key;
}

static void order_numeric_alleles(marker_type *marker)
{
    int i, nalleles = marker->num_alleles;
    struct kv *order;
    allele_list_type *allelep;
    int allele_index;

    if (nalleles <= 0) {
        return;
    }
    order = (nalleles < 100) ? kvstore : CALLOC((size_t) nalleles, struct kv);
    allelep = marker->first_allele;

    for (i=0; i < nalleles; i++) {
        allele_index = atoi(allelep->allele_freq.name);
        if (allele_index <= 0) {
            marker->recode_alleles = 1;
            if (order != kvstore)
                free(order);
            return;
        }
        order[i].key = allele_index;
        order[i].val = allelep;
        allelep = allelep->next;
    }

    qsort(order, nalleles, sizeof(struct kv), kvsort);

    allelep = order[0].val;
    for (i=0; i < nalleles-1; i++) {
        order[i].val->next = order[i+1].val;
        order[i].val->allele_freq.index = i+1;
        order[i].val = 0;
    }
    order[i].val->next = NULL;
    order[i].val->allele_freq.index = i+1;
    if (order[i].key != i+1)
        marker->recode_alleles = 1;
    order[i].val = 0;
    marker->first_allele = allelep;

/* show list for debugging
    allelep = marker->first_allele;
    for (i=0; i < marker->num_alleles; i++) {
        printf("sort #%d %d: %d %s\n",
               marker->num_alleles, i, allelep->allele_freq.index, allelep->allele_freq.name);
        allelep = allelep->next;
    }
    printf("\n");
*/

    if (order != kvstore)
        free(order);
}

/* function to check if the alleles have to be recoded */

void need_recoding(marker_type *marker)
{
    order_numeric_alleles(marker);
}



/* Function to add up allele counts for all 4 options.
   Then assign the correct count depending on count_option, and
   also compute the frequencies by dividing each count by the total
   for that option
*/

void convert_to_freq(marker_type *marker_list,
		     linkage_locus_top *LTop,
		     int count_option1,
		     analysis_type analysis)

{

    int count_option, m;
    int total_count;
    allele_list_type *allelep;

    for (m=0; m < LTop->LocusCnt; m++) {
        if (LTop->Locus[m].Type == NUMBERED ||
                LTop->Locus[m].Type == XLINKED ||
                LTop->Locus[m].Type == YLINKED) {

            if (analysis == TO_PLINK && marker_list[m].num_alleles <= 2) {
                continue;
            }

            if (analysis == TO_PLINK) {
                marker_list[m].estimate_frequencies = 1;
            }

            if (marker_list[m].estimate_frequencies) {
                count_option=count_option1;
            } else {
                count_option=0;
            }


            /*      LTop->Locus[m].Type = NUMBERED; */
            switch(count_option) {
            case 1:
                marker_list[m].num_people = marker_list[m].num_founders;
                break;
            case 2:
                marker_list[m].num_people = marker_list[m].num_random;
                break;
            case 3:
                marker_list[m].num_people = marker_list[m].num_unique;
                break;
            case 4:
                marker_list[m].num_people = marker_list[m].num_everyone;
                break;
            default:
                break;
            }

            total_count=0;
            allelep = marker_list[m].first_allele;
            while(allelep != NULL) {
                switch(count_option) {
                case 1:
                    allelep->allele_freq.count = allelep->allele_freq.founder_count;
                    break;
                case 2:
                    allelep->allele_freq.count = allelep->allele_freq.random_count;
                    break;
                case 3:
                    allelep->allele_freq.count = allelep->allele_freq.unique_count;
                    break;
                case 0:
                    /* This is for the annotated pedigree file option */
                case 4:
                    allelep->allele_freq.count = allelep->allele_freq.everyone_count;
                    break;
                default:
                    break;
                }
                total_count += allelep->allele_freq.count;
                allelep = allelep->next;
            }
            marker_list[m].num_total_alleles = total_count;
            allelep = marker_list[m].first_allele;
            while(allelep != NULL) {
                if (!marker_list[m].recode_alleles) {
                    allelep->allele_freq.index = atoi(allelep->allele_freq.name);
                }

                /* compute frequency only if non-annotated */
                if (count_option > 0) {
                    if (total_count > 0) {
                        allelep->allele_freq.freq =
                            (double)(allelep->allele_freq.count)/(double)total_count;
                    } else {
                        if (allelep->next == NULL) {
                            if (LTop->Locus[m].number != -1) {
                                warnf("No genotyped individuals found among those selected.");
                                sprintf(err_msg,
                                        "All allele frequencies will be set to 0 for marker %s.",
                                        LTop->Locus[m].Name);
                                warnf(err_msg);
                            }
                        }
                        allelep->allele_freq.freq = 0.0;
                    }
                }
                allelep = allelep->next;
            }
        }
    }
}

/*----------------------end of convert_to_freq----------------------*/

/*--------------assign_dummy_alleles----------------- */
/* Function to assign one or two dummy alleles to markers
   that have less than 2 alleles */


void assign_dummy_alleles(marker_type *marker_list,
			  linkage_locus_top *LTop)
{
    int m;
    allele_list_type *all;


    for(m = LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
        if (marker_list[m].estimate_frequencies == 0 ) continue;

        if (LTop->Locus[m].Type == NUMBERED ||
           LTop->Locus[m].Type == XLINKED ||
           LTop->Locus[m].Type == YLINKED) {

            if (marker_list[m].num_alleles == 1) {
                marker_list[m].first_allele->allele_freq.freq=1.0;
                all=CALLOC((size_t) 1, allele_list_type);
                all->allele_freq.freq  = 0.0;
                all->allele_freq.name  = canonical_allele("dummy");
                all->allele_freq.index = 2;
                all->next = NULL;
                marker_list[m].first_allele->next = all;
                marker_list[m].num_alleles=2;
            } else if (marker_list[m].num_alleles == 0) {
                all = CALLOC((size_t) 1, allele_list_type);
                all->allele_freq.freq  = 0.5;
                all->allele_freq.name  = canonical_allele("dummy1");
                all->allele_freq.index = 1;
                marker_list[m].first_allele = all;

                all = CALLOC((size_t) 1, allele_list_type);
                all->allele_freq.freq  = 0.5;
                all->allele_freq.name  = canonical_allele("dummy2");
                all->allele_freq.index = 2;
                all->next = NULL;
                marker_list[m].first_allele->next = all;

                marker_list[m].num_alleles=2;

            }
        }
    }
    return;
}

void write_recode_summary(marker_type *markers,
			  linkage_locus_top *LTop, int inc_ht)

{

    int m, found_freq_est = 0;
    int found_recode_alleles = 0;
    char count_info[FILENAME_LENGTH];
    char recode_info[FILENAME_LENGTH];
    double founder_freq, random_freq, unique_freq, everyone_freq;
    int tot_founder, tot_random, tot_unique, tot_everyone;

    FILE *fp;
    allele_list_type *allelep;

    for (m=LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
        if (markers[m].estimate_frequencies == 1) {
            found_freq_est++;
        }
        if (markers[m].recode_alleles) {
            found_recode_alleles++;
        }
    }
    if (found_freq_est == 0 && found_recode_alleles == 0) return;

    if ((fp = fopen(Mega2RecodeRun, "w")) == NULL) {
      errorvf("Could not open %s for writing.\n", Mega2RecodeRun);
      EXIT(FILE_WRITE_ERROR);
    }

    if (found_freq_est) {
        sprintf(recode_info,
                "Allele frequency computed from: \n");

        switch(SelectIndividuals) {
        case 2:
            strcat(recode_info,
                   "     Genotyped founders or a randomly chosen individual");

            break;
        case 1:
            strcat(recode_info, "     Genotyped founders only");


            break;
        case 3:
            strcat(recode_info,
                   "     Genotyped founders or a randomly chosen individual and \n");
            strcat(recode_info,
                   "     individuals with unique alleles");
            break;
        case 4:
            strcat(recode_info,
                   "     All genotyped individuals");

            break;
        default:
            break;
        }

        strcat(recode_info,
               ((inc_ht)?
                "\n     including half-typed individuals." :
                "\n     excluding half-typed individuals."));
        strcat(recode_info,
               "\n---------------------\nRecoded loci:");

    } else {
        strcpy(recode_info, "No frequency estimates done on this data.");
    }

    summary_time_stamp(mega2_input_files, fp, recode_info);

    for (m=0; m < LTop->LocusCnt; m++) {

        if (LTop->Locus[m].Type == NUMBERED || LTop->Locus[m].Type == XLINKED ||
            LTop->Locus[m].Type == YLINKED) {

            if (markers[m].estimate_frequencies) {
                /* Full summary  is to be written */

                /* compute the other frequencies */
                tot_founder = tot_random = tot_unique = tot_everyone=0;

                allelep = markers[m].first_allele;
                while(allelep != NULL) {
                    tot_founder += allelep->allele_freq.founder_count;
                    tot_random += allelep->allele_freq.random_count;
                    tot_unique += allelep->allele_freq.unique_count;
                    tot_everyone += allelep->allele_freq.everyone_count;
                    allelep = allelep->next;
                }

                allelep = markers[m].first_allele;
                if (LTop->Locus[m].Type == YLINKED) {
                    fprintf(fp, "Heterozygous males have been excluded.\n");
                }
                sprintf(count_info,
                        "%d founders, %d randomly selected individuals",
                        markers[m].num_founders,
                        markers[m].num_random - markers[m].num_founders);

                if (markers[m].num_random < markers[m].num_unique) {
                    grow(count_info,
                         ", %d individuals with unique alleles",
                         markers[m].num_unique - markers[m].num_random);
                }
                fprintf(fp, "Frequencies and counts:\n");
                fprintf(fp,
                        "   1=Founders       2=1+Random       3=2+Unique     4=Everyone\n");
                fprintf(fp,
                        "Marker %d: %s has %d alleles; \n",
                        m+1, LTop->Locus[m].Name, markers[m].num_alleles);
                fprintf(fp, "Option %d frequencies computed using %d %s:",
                        SelectIndividuals,
                        markers[m].num_people,
                        ((SelectIndividuals == 1)? "founders": "individuals"));
                if (markers[m].num_half_typed > 0) {
                    fprintf(fp, " including %d half-typed individuals",
                            markers[m].num_half_typed);
                }
                fprintf(fp, "\n");
                if (SelectIndividuals == 2 || SelectIndividuals == 3) {
                    fprintf(fp, "%s\n", count_info);
                }

                fprintf(fp, "Allele  Code      ");
                fprintf(fp,
                        "Freq1     Count1  Freq2     Count2  Freq3     Count3  Freq4     Count4\n");

                allelep = markers[m].first_allele;
                while(allelep != NULL) {
                    founder_freq =
                        ((tot_founder > 0)?
                         (double)allelep->allele_freq.founder_count/(double)tot_founder : 0);
                    random_freq =
                        ((tot_random > 0)?
                         (double)allelep->allele_freq.random_count/(double)tot_random : 0);
                    unique_freq =
                        ((tot_unique > 0)?
                         (double)allelep->allele_freq.unique_count/(double)tot_unique : 0);
                    everyone_freq =
                        ((tot_everyone > 0)?
                         (double)allelep->allele_freq.everyone_count/(double)tot_everyone: 0);

                    fprintf(fp, "%6d %-10s %7.5f  %7d  %7.5f  %7d  %7.5f  %7d  %7.5f  %7d\n",
                            allelep->allele_freq.index,
                            allelep->allele_freq.name,
                            founder_freq, allelep->allele_freq.founder_count,
                            random_freq,  allelep->allele_freq.random_count,
                            unique_freq,  allelep->allele_freq.unique_count,
                            everyone_freq,  allelep->allele_freq.everyone_count);
                    allelep = allelep->next;
                }
                fprintf(fp, "TOTAL             %16d  %16d  %16d  %16d\n",
                        tot_founder, tot_random, tot_unique, tot_everyone);
                if (inc_ht) {
                    fprintf(fp, "Half-typed        %16d  %16d  %16d  %16d\n",
                            markers[m].ht_founders, markers[m].ht_random,
                            markers[m].ht_unique, markers[m].ht_everyone);
                }
                fprintf(fp, "\n");

            } else {
                if (markers[m].recode_alleles) {
                    /* Only write the allele code mapping */

                    fprintf(fp, "Marker %d: %s has %d alleles; \n",
                            m+1, LTop->Locus[m].Name, markers[m].num_alleles);
                    fprintf(fp, "Allele  Code\n");
                    allelep = markers[m].first_allele;
                    while(allelep != NULL) {
                        fprintf(fp, "%6d  %-10s\n", allelep->allele_freq.index,
                                allelep->allele_freq.name);
                        allelep = allelep->next;
                    }
                }
                /* else skip this marker from the summary */
            }
        }
    }

    fclose(fp);

}

void recode_liability_class(pheno_type *pheno_list, linkage_locus_top *LTop,  int locus)
{
    int all, nc;
    class_list_type *recoded_class;

    LTop->Pheno[locus].Props.Affection.PenCnt = 3;

    if (LTop->Pheno[locus].Props.Affection.Class) {
        free(LTop->Pheno[locus].Props.Affection.Class[0].AutoPen);
        free(LTop->Pheno[locus].Props.Affection.Class[0].FemalePen);
        free(LTop->Pheno[locus].Props.Affection.Class[0].MalePen);
        free(LTop->Pheno[locus].Props.Affection.Class);
    }

    LTop->Pheno[locus].Props.Affection.Class =
        CALLOC((size_t) pheno_list[locus].num_classes, linkage_affection_class);

    nc = pheno_list[locus].num_classes;

    for (all = 0; all < nc; all++) {
        clear_laffclass(&(LTop->Pheno[locus].Props.Affection.Class[all]));
        LTop->Pheno[locus].Props.Affection.Class[all].AutoDef     = 1;
        LTop->Pheno[locus].Props.Affection.Class[all].AutoPen    = CALLOC((size_t) 3, double);
        LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[0] = DOM_G11;
        LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[1] = DOM_G12;
        LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[2] = DOM_G22;

        LTop->Pheno[locus].Props.Affection.Class[all].FemaleDef    = 1;
        LTop->Pheno[locus].Props.Affection.Class[all].FemalePen    = CALLOC((size_t) 3, double);
        LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[0] = DOM_G11;
        LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[1] = DOM_G12;
        LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[2] = DOM_G22;

        LTop->Pheno[locus].Props.Affection.Class[all].MaleDef    = 1;
        LTop->Pheno[locus].Props.Affection.Class[all].MalePen    = CALLOC((size_t) 2, double);
        LTop->Pheno[locus].Props.Affection.Class[all].MalePen[0] = DOM_G11;
        LTop->Pheno[locus].Props.Affection.Class[all].MalePen[1] = DOM_G22;
    }

    recoded_class = &(pheno_list[locus].class_list[0]);

    while (recoded_class != NULL) {

        all = recoded_class->l_class.num - 1;
        if (recoded_class->l_class.autosomal_pen != NULL) {
            LTop->Pheno[locus].Props.Affection.Class[all].AutoDef    = 0;
            LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[0] =
                recoded_class->l_class.autosomal_pen[0];
            LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[1] =
                recoded_class->l_class.autosomal_pen[1];
            LTop->Pheno[locus].Props.Affection.Class[all].AutoPen[2] =
                recoded_class->l_class.autosomal_pen[2];
        }

        if (recoded_class->l_class.female_pen != NULL) {
            LTop->Pheno[locus].Props.Affection.Class[all].FemaleDef    = 0;
            LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[0] =
                recoded_class->l_class.female_pen[0];
            LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[1] =
                recoded_class->l_class.female_pen[1];
            LTop->Pheno[locus].Props.Affection.Class[all].FemalePen[2] =
                recoded_class->l_class.female_pen[2];
        }

        if (recoded_class->l_class.male_pen != NULL) {
            LTop->Pheno[locus].Props.Affection.Class[all].MaleDef    = 0;
            LTop->Pheno[locus].Props.Affection.Class[all].MalePen[0] =
                recoded_class->l_class.male_pen[0];
            LTop->Pheno[locus].Props.Affection.Class[all].MalePen[1] =
                recoded_class->l_class.male_pen[1];
        }

        recoded_class = recoded_class->next;

    }
}

/* Function to fill up the linkage locus top (mostly with default
   values, except for the NUMBERED alleles, where we have frequencies
   and allele_counts */

void recode_locus_top(marker_type *marker_list, pheno_type *pheno_list, linkage_locus_top *LTop)
{
    int m, all;
    allele_list_type *recoded_allele;

    /* First assign number of alleles and frequencies from the
       marker list to the Locus Top */


    for (m=0; m < LTop->LocusCnt; m++) {
        switch(LTop->Locus[m].Type) {
        case NUMBERED:
        case XLINKED:
        case YLINKED:
            LTop->Locus[m].Class = MARKER;
            LTop->Locus[m].AlleleCnt = marker_list[m].num_alleles;
            LTop->Locus[m].Allele
                = CALLOC((size_t) marker_list[m].num_alleles, linkage_allele_rec);

            recoded_allele = marker_list[m].first_allele;

            for (all = 0; all < marker_list[m].num_alleles; all++) {
                LTop->Locus[m].Allele[all].Frequency
                    = recoded_allele->allele_freq.freq;
                LTop->Locus[m].Allele[all].index = recoded_allele->allele_freq.index;
                LTop->Locus[m].Allele[all].name = recoded_allele->allele_freq.name;
                recoded_allele = recoded_allele->next;
            }
            LTop->Marker[m].Props.Numbered.Recoded = marker_list[m].recode_alleles;
            LTop->Marker[m].Props.Numbered.EstimateFrequencies = marker_list[m].estimate_frequencies;
            LTop->Marker[m].Props.Numbered.NumAlleles = marker_list[m].num_total_alleles;
            LTop->Marker[m].Props.Numbered.SelectOpt = SelectIndividuals;
            break;

        case AFFECTION:
            LTop->Locus[m].Class = TRAIT;
            /* This has not yet been set up */
            LTop->Locus[m].AlleleCnt = 2;
            if (LTop->Locus[m].Allele) free(LTop->Locus[m].Allele);
            LTop->Locus[m].Allele = 	CALLOC((size_t) 2, linkage_allele_rec);

            LTop->Locus[m].Allele[0].Frequency = pheno_list[m].first_allele->allele_freq.freq;
            LTop->Locus[m].Allele[1].Frequency = pheno_list[m].first_allele->next->allele_freq.freq;

//          LTop->Pheno[m].Props.Affection.ClassCnt = marker_list[m].num_classes;

            recode_liability_class(pheno_list, LTop, m);

            break;

        case QUANT:
            LTop->Locus[m].Class = TRAIT;
            LTop->Locus[m].AlleleCnt = 2;
            LTop->Locus[m].Allele = CALLOC((size_t) 2, linkage_allele_rec);
            LTop->Locus[m].Allele[0].Frequency =
                pheno_list[m].first_allele->allele_freq.freq;
            LTop->Locus[m].Allele[1].Frequency =
                pheno_list[m].first_allele->next->allele_freq.freq;

            LTop->Pheno[m].Props.Quant.ClassCnt=1;
            LTop->Pheno[m].Props.Quant.Mean = CALLOC((size_t) 1, double *);
            LTop->Pheno[m].Props.Quant.Mean[0] = CALLOC((size_t) 3, double);
            LTop->Pheno[m].Props.Quant.Mean[0][0] = 0.0;
            LTop->Pheno[m].Props.Quant.Mean[0][1] = 0.0;
            LTop->Pheno[m].Props.Quant.Variance = CALLOC((size_t) 1, double *);
            LTop->Pheno[m].Props.Quant.Variance[0] = CALLOC((size_t) 1, double);
            LTop->Pheno[m].Props.Quant.Variance[0][0] = 1.0;
            LTop->Pheno[m].Props.Quant.Multiplier = 1.0;
            break;

        default:
            break;
        }
    }

    LTop->SexDiff = LTop->Interference = 0;
    if (LTop->MarkerCnt > 1) {
        LTop->MaleRecomb = CALLOC((size_t) LTop->LocusCnt - 1, double);
        for (all = 0; all < LTop->LocusCnt - 1; all++)
            LTop->MaleRecomb[all] = 0.0;
    }
    LTop->PedRecDataType = ((LTop->PedRecDataType == Raw_premake)?
                            Premakeped : Postmakeped);
    LTop->FlankRecomb = UNDEF;
    LTop->FMRatio = UNDEF;
    LTop->FemaleRecomb = NULL;
    LTop->Run.linkmap.recomb_frac=NULL; /* leave it to vitesse
                                           to allocate */
    LTop->Run.mlink.start_theta = 1;
    LTop->Run.mlink.increment = 0.05;
    LTop->Run.mlink.stop_theta = 0.0;
    LTop->Run.mlink.num_evals=0; /* for now assume this will be put in later
				    inside Vitesse */
    return;


}

/*-----------------------end of recode_locus_top------------------*/

/* Function to recode genotypes in all pedigrees,
   Iterate over each marker, and each allele, then recode
   each person whose genotype at that marker is equal to the
   allele-name. */

void recode_ped_top(marker_type *marker_list, linkage_ped_top *Top, plink_info_type *plink_info)

{
    int m, ped, per, entrycnt;
    allele_list_type *allele;
    void *marker;
    void *marker_copy;
    int rall;
    const char *all;
    int geno_recoded[2], allele_unrecoded=0;
    int displayed_errors=0;
    const char *all1, *all2;
    int a1, a2;

#ifdef DEBUG
    int num_skip;

    if (Top->LocusTop->LocusCnt <= 10) {
        num_skip=1;
    } else if (Top->LocusTop->LocusCnt <= 100) {
        num_skip = 10;
    } else {
        num_skip = 100;
    }
#endif

    Display_Errors=1;

#ifndef HIDESTATUS
    mssgf("Recoding pedigree genotypes ... ");
#endif
    for(ped=0; ped < Top->PedCnt; ped++) {
        entrycnt = ((pedfile_type == POSTMAKEPED_PFT) ?
                    Top->Ped[ped].EntryCnt :
                    Top->PTop[ped].num_persons);
        for(per=0; per < entrycnt; per++) {

            marker = (pedfile_type == POSTMAKEPED_PFT)?
                      Top->Ped[ped].Entry[per].Marker :
                      Top->PTop[ped].persons[per].marker;
            if (marker != NOTYPED_ALLELES) {
                marker_copy = marker_alloc((size_t) Top->LocusTop->MarkerCnt, Top->LocusTop->PhenoCnt);
                for(m = Top->LocusTop->PhenoCnt; m < Top->LocusTop->LocusCnt; m++) {
                    a1 = a2 = 0;
                    if (UntypedPeds[ped] == 1) {
                        geno_recoded[0]=geno_recoded[1]=1;
                    } else {
                        get_2Ralleles(marker, m, &all1, &all2);
                        geno_recoded[0]=geno_recoded[1]=0;

                        if (allelecmp(all1, REC_UNKNOWN) == 0) {
                            a1 = 0;
                            geno_recoded[0]=1;
                        }

                        if (allelecmp(all2, REC_UNKNOWN) == 0) {
                            a2 = 0;
                            geno_recoded[1]=1;
                        }

                        if (allelecmp(all1, REC_UNKNOWN) ||
                            allelecmp(all2, REC_UNKNOWN)) {
                            allele = marker_list[m].first_allele;
                            while(allele != NULL) {
                                rall = allele->allele_freq.index;
                                 all = allele->allele_freq.name;

                                if (allelecmp(all1, all) == 0) {
                                    a1 = rall;
                                    geno_recoded[0]=1;
                                }

                                if (allelecmp(all2, all) == 0) {
                                    a2 = rall;
                                    geno_recoded[1]=1;
                                }
                                allele = allele->next;
                            }
                        }
                    }
                    set_2alleles(marker_copy, m,
                                 &Top->LocusTop->Locus[m], a1, a2);
                    /* free the raw alleles */
                    if (geno_recoded[0] == 0 || geno_recoded[1] == 0) {
                        if (pedfile_type == POSTMAKEPED_PFT) {
                            sprintf(err_msg,
                                    "Ped %s, person %s, marker %s: ",
                                    Top->Ped[ped].Name,
                                    Top->Ped[ped].Entry[per].UniqueID,
                                    Top->LocusTop->Locus[m].Name);
                        } else {
                            sprintf(err_msg, "Ped %s, person %s, marker %s: ",
                                    Top->PTop[ped].Name,
                                    Top->PTop[ped].persons[per].uniqueid,
                                    Top->LocusTop->Locus[m].Name);
                        }
                        if (geno_recoded[0] == 0 && geno_recoded[1] == 0) {
                            strcat(err_msg, " both alleles unrecognized");
                        } else if (geno_recoded[0] == 0) {
                            strcat(err_msg, " allele 1 unrecognized");
                        } else {
                            strcat(err_msg, " allele 2 unrecognized");
                        }

                        grow(err_msg, " in genotype %s/%s", all1, all2);

                        SUPPRESS_MSSG(displayed_errors)
                            if (displayed_errors > MAX_PED_ERRORS) {
                                Display_Errors = 0;
                            }
                        displayed_errors++;
                        errorf(err_msg);
                        allele_unrecoded=1;
                    }
                }
// put it back
                if (Top->pedfile_type == POSTMAKEPED_PFT) {
                    Top->Ped[ped].Entry[per].Marker = marker_copy;
                } else {
                    Top->PTop[ped].persons[per].marker = marker_copy;
                }

                /* finally someone did it right ! */
                marker_free(marker, Top->LocusTop->PhenoCnt);
            }
        }
    }
    
    if (Display_Errors == 0) {
        printf("Check error logs for full list of unrecoded alleles.\n");
        draw_line();
        Display_Errors=1;
    }


    if (allele_unrecoded) {
        log_line(mssgf);
        errorf("Alleles unrecoded due to one or more of these reasons: ");
        errorf("  1) Unknown allele definition is incorrect.");
        errorf("  2) Allele labels missing from frequency file for these markers.");
        errorf("  3) Alleles not encountered in genotypes selected for frequency estimation.");

        mssgf("Please correct errors and restart Mega2.");
        EXIT(DATA_INCONSISTENCY);
    }

    order_heterozygous_allele(Top);
}

/*-----------------------end of recode_ped_top------------------*/


/*  void copy_pedrec_rawdata(raw_allele_data *d1, raw_allele_data *d2,  */
/*  			 linkage_locus_type Type) */

/*  { */

/*    switch (Type) { */
/*    case QUANT: */
/*      break; */
/*    case AFFECTION: */
/*      break; */
/*    case NUMBERED: */
/*      d2->Allele_1 = d1->Allele_1; */
/*      d2->Allele_2 = d1->Allele_2; */
/*      d2->Recoded_1 = d1->Recoded_1; */
/*      d2->Recoded_2 = d1->Recoded_2; */

/*      break; */
/*    default: */
/*      fprintf(stderr, "ERROR: unknown data type!\n"); */
/*      EXIT(-1); */
/*      break; */
/*    } */
/*    return; */
/*  } */

/*---------------------------------------------------------------*/

/* Function to read in the names file into a preliminary
   linkage_locus_top  structure, without allele counts, frequencies
   etc, only with names and types.
*/
linkage_locus_top *read_marker_data(FILE *fp,
				    int type_col,
				    int name_col)
{
    int i, num_markers, annotated=0;
    char nextline[FILENAME_LENGTH];
    char ltype[3], marker[100];
    char *types, **names;
    int penetrances_read = 0;
    double freq, pen[5];
    int num_read, lch, line_num;
    /* error flags */
    int invalid_marker_type =0, missing_columns = 0;
    int displayed_messages = 0;
    int loci_name_space = 0;

    num_markers=linecount(fp);

    types = CALLOC((size_t) num_markers, char);
    names = CALLOC((size_t) num_markers, char *);
    num_markers = 0;

    Display_Errors = 1;

    if (type_col < 0) {
        type_col = 0;
        name_col = 1;
        annotated=0;
        line_num=0;
    } else {
        /* annotated format */
        annotated=1;
        skip_line(fp, lch);
        line_num=1;
    }
    i=-1;
    while(!feof(fp)) {
        *nextline = 0;
        (void)fgets(nextline, FILENAME_LENGTH - 1, fp);
        if (*nextline == 0) {
            break;
        }
        if (nextline[0] == '#' || nextline[0] == '\n') {
            continue;
        }
        if (type_col == 0) {
            num_read=sscanf(nextline, "%s %s", ltype, marker);
        } else {
            num_read=sscanf(nextline, "%s %s", marker, ltype);
        }
        line_num++;
        if (num_read > 0 && num_read < 2) {
            sprintf(err_msg,
                    "Line %d: Found fewer than 2 columns (at least 2 required).",
                    line_num);
            SUPPRESS_MSSG(displayed_messages);
            if (displayed_messages > MAX_PED_ERRORS) {
                Display_Errors = 0;
            }
            errorf(err_msg);
            displayed_messages++;
            missing_columns++;
        } else {
            if ((!strcmp(ltype, "M")) || (!strcmp(ltype, "A")) ||
               (!strcmp(ltype, "L")) || (!strcmp(ltype, "C")) ||
               (!strcmp(ltype, "X")) || (!strcmp(ltype, "Y")) ||
               (!strcmp(ltype, "T"))) {
                i++;
                types[i]=ltype[0];
                names[i] = CALLOC(strlen(marker)+1, char);
                loci_name_space += strlen(marker)+1;
                strcpy(names[i], marker);
                if ((!strcmp(ltype, "M")) || (!strcmp(ltype, "X")) || (!strcmp(ltype, "Y")) )
                    num_markers++;
                if (!annotated && (!strcmp(ltype, "A"))) {
                    penetrances_read=0;
                    penetrances_read =
                        sscanf(nextline, "%s %s %lf %lf %lf %lf %lf %lf",
                               ltype, marker, &freq,
                               &(pen[0]), &(pen[1]), &(pen[2]),
                               &(pen[3]), &(pen[4]));
                }
            } else {
                sprintf(err_msg,
                        "Marker type %s for marker %s not recognized.",
                        ltype, marker);
                SUPPRESS_MSSG(displayed_messages);
                if (displayed_messages > MAX_PED_ERRORS) {
                    Display_Errors = 0;
                }
                errorf(err_msg);
                displayed_messages++;
                invalid_marker_type++;
            }
        }
    }
    Display_Errors = 1;

    if (missing_columns) {
        printf("Found lines with less than 2 columns, see %s for details.\n",
               Mega2Err);
    }

    if (invalid_marker_type) {
        printf("Invalid marker types, see %s for details.\n", Mega2Err);

    }

    if (missing_columns || invalid_marker_type) {
        errorf("Found fatal errors in locus file.");
        errorf("Please fix errors and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: locus strings: %d MB (%d x %d)\n",
          loci_name_space / 1024 / 1024,
          i+1,  loci_name_space / (i+1) );
#endif

    return read_common_marker_data(i+1, num_markers, names, types, annotated, penetrances_read, freq, pen);
}

linkage_locus_top *read_common_marker_data(int all_loci, int num_markers, char **names, char *types, int annotated, double penetrances_read, double freq, double pen[5])
{
    int i, m, p;
    int cnt_x = 0, cnt_y = 0, cnt_a = 0;
    linkage_locus_rec *Locus;
    linkage_locus_top *LTop = new_llocustop();

    LTop->LocusCnt = all_loci;
    LTop->NumPedigreeCols = 0;
    LTop->RiskLocus=0;
    LTop->SexLinked=0;
    LTop->Program=MLINK;
    LTop->MutLocus=0;
    LTop->MutMale=LTop->MutFemale=0;
    LTop->Haplotype=0;

    LTop->Locus   = CALLOC((size_t) all_loci,    linkage_locus_rec);

    LTop->PhenoCnt = all_loci - num_markers;
    LTop->Pheno  = (LTop->PhenoCnt > 0) ? CALLOC((size_t) LTop->PhenoCnt, pheno_rec) : 0;

    LTop->MarkerCnt = num_markers;
    LTop->Marker = (LTop->MarkerCnt > 0) ? (CALLOC((size_t) LTop->MarkerCnt, marker_rec) - LTop->PhenoCnt) : 0;
#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: LTop->Locus: %d MB (%d x %d)\n",
          all_loci * sizeof(linkage_locus_rec) / 1024 / 1024,
          all_loci,  sizeof(linkage_locus_rec));

    msgvf("ALLOC SPACE: LTop->Pheno: %d MB (%d x %d)\n",
          (all_loci-num_markers) * sizeof(pheno_rec) / 1024 / 1024,
          (all_loci-num_markers),  sizeof(pheno_rec));

    msgvf("ALLOC SPACE: LTop->Marker: %d MB (%d x %d)\n",
          num_markers * sizeof(marker_rec) / 1024 / 1024,
          num_markers,  sizeof(marker_rec));
#endif

    LTop->NumPedigreeCols = 0;
    for (i=0, m=0, p=0; i < all_loci; i++) {
        int idx;
        marker_rec *Marker;
        pheno_rec *Pheno;
        switch(types[i]) {
        case 'M':
        case 'X':
        case 'Y':
            idx = LTop->PhenoCnt + m++;
            Locus = &(LTop->Locus[idx]);
            Locus->number = -1;
            Locus->Name = CALLOC(strlen(names[i])+1, char);
            strcpy(Locus->Name, names[i]);
            Locus->Class = MARKER;
            HasMarkers=1;
            Locus->col_num = LTop->NumPedigreeCols;
            Locus->AlleleCnt=0;

            Marker = &(LTop->Marker[idx]);
            Locus->Marker = Marker;
            Locus->Pheno  = 0;
            Marker->Name = Locus->Name;
            Marker->col_num = LTop->NumPedigreeCols;
            LTop->NumPedigreeCols += 2;
            /*    clear_llocusrec(Locus, TYPE_UNSET); */
            switch(types[i]) {
            case 'M':
                Locus->Type = NUMBERED;
                cnt_a++;
                break;

            case 'X':
                Locus->Type = XLINKED;
                cnt_x++;
                break;

            case 'Y':
                Locus->Type = YLINKED;
                cnt_y++;
                break;
            }
            break;
        case 'A':
        case 'L':
        case 'T':
        case 'C':
            idx = p++;
            Locus = &(LTop->Locus[idx]);
            Locus->number = -1;
            Locus->Name = CALLOC(strlen(names[i])+1, char);
            strcpy(Locus->Name, names[i]);
            Locus->Class = TRAIT;

            Pheno = &(LTop->Pheno[idx]);
            Locus->Marker = 0;
            Locus->Pheno  = Pheno;
            Pheno->Name = Locus->Name;
            Locus->col_num = LTop->NumPedigreeCols;
            Pheno->col_num = LTop->NumPedigreeCols;

            switch(types[i]) {
            case 'A':
            case 'L':
            Locus->Type = AFFECTION;
            if (types[i] == 'L') {
                LTop->NumPedigreeCols += 2;
                Pheno->Props.Affection.ClassCnt = 2;
                Pheno->Props.Affection.Class = NULL;
                Locus->AlleleCnt = 0;
            } else {
                LTop->NumPedigreeCols ++;
                Pheno->Props.Affection.ClassCnt = 1;
                Locus->AlleleCnt = 2;
                Locus->Allele = CALLOC((size_t) 2, linkage_allele_rec);
                Pheno->Props.Affection.Class = CALLOC((size_t) 1, linkage_affection_class);
                Pheno->Props.Affection.PenCnt = 3;
                clear_laffclass(&(Pheno->Props.Affection.Class[0]));
                Pheno->Props.Affection.Class[0].AutoPen   = CALLOC((size_t) 3, double);
                Pheno->Props.Affection.Class[0].FemalePen = CALLOC((size_t) 3, double);
                Pheno->Props.Affection.Class[0].MalePen   = CALLOC((size_t) 2, double);

                if (!annotated) {
                    if (penetrances_read < 6) {
                        sprintf(err_msg, "For affection status locus %s :", Locus->Name);
                        mssgf(err_msg);
                        mssgf("Names file did not contain allele frequency and penetrances.");
                        mssgf("Setting penetrances to default values (A/A=0.05 A/a=0.90 a/a=0.90).");
                        mssgf("and setting the allele frequencies to (0.50, 0.50)");
                        Locus->Allele[0].Frequency = Locus->Allele[1].Frequency = 0.5;
                        Pheno->Props.Affection.Class[0].AutoPen[0] = DOM_G11;
                        Pheno->Props.Affection.Class[0].AutoPen[1] = DOM_G12;
                        Pheno->Props.Affection.Class[0].AutoPen[2] = DOM_G22;

                        Pheno->Props.Affection.Class[0].FemalePen[0] = DOM_G11;
                        Pheno->Props.Affection.Class[0].FemalePen[1] = DOM_G12;
                        Pheno->Props.Affection.Class[0].FemalePen[2] = DOM_G22;

                        Pheno->Props.Affection.Class[0].MalePen[0] = DOM_G11;
                        Pheno->Props.Affection.Class[0].MalePen[1] = DOM_G22;
                        draw_line();
                    } else {
                        Locus->Allele[0].Frequency = 1 - freq;
                        Locus->Allele[1].Frequency = freq;
                        if (penetrances_read == 8) {
                            Pheno->Props.Affection.Class[0].FemalePen[0] = pen[0];
                            Pheno->Props.Affection.Class[0].FemalePen[1] = pen[1];
                            Pheno->Props.Affection.Class[0].FemalePen[2] = pen[2];
                            Pheno->Props.Affection.Class[0].MalePen[0] = pen[3];
                            Pheno->Props.Affection.Class[0].MalePen[1] = pen[4];
                        } else {
                            Pheno->Props.Affection.Class[0].AutoPen[0] = pen[0];
                            Pheno->Props.Affection.Class[0].AutoPen[1] = pen[1];
                            Pheno->Props.Affection.Class[0].AutoPen[2] = pen[2];

                            sprintf(err_msg, "For affection status locus %s :", Locus->Name);
                            mssgf(err_msg);
                            mssgf("Names file did not contain 2 male penetrances.");
                            mssgf("Assuming data is for Autosomes Only.");
                            draw_line();
                        }
                    }
                }
            }
            break;

            case 'T':
            case 'C':
                Locus->Type = QUANT;
                LTop->NumPedigreeCols ++;
                Locus->AlleleCnt=0;
                break;

            default:
                break;
            }
            break;
        default:
            break;
        }

        /* free(names[i]); */
        free(names[i]); // This only works if everyone copies here, and currently they seem to do so
    }
    LTop->SexLinked = (cnt_x + cnt_y) ? (cnt_a ? 2 : 1) : 0;

    // cpk: There is a memory leak here of ltype[] because it is generated by a strcpy,
    // and we only put the first char if it into types[]. 
    free(types);
    free(names);
/*   get_trait_list(LTop); */
    return LTop;
}

/*----------------------end of read_marker_data------------------------*/
//r default_trait_penetrance() below is called by some one that does not know about Locus split
void default_trait_penetrances(linkage_ped_top *Top, int locus,
                               pheno_type *pheno_listi)

{
    int ped, per, entrycount;
    int lclass, i;

    if (Top->LocusTop->Locus[locus].Type == AFFECTION) {
        /* Just store the number of classes,
           ClassCnt = 0 => multiple liability class possible */
        for (ped=0; ped < Top->PedCnt; ped++) {
            if (pedfile_type == POSTMAKEPED_PFT) {
                entrycount = Top->Ped[ped].EntryCnt;
            } else {
                entrycount = Top->PTop[ped].num_persons;
            }

            for(per=0; per < entrycount; per++) {
                if (Top->LocusTop->Pheno[locus].Props.Affection.ClassCnt == 0) {
                    if (pedfile_type == POSTMAKEPED_PFT) {
                        lclass = Top->Ped[ped].Entry[per].Pheno[locus].Affection.Class;
                    } else {
                        lclass = Top->PTop[ped].persons[per].pheno[locus].Affection.Class;
                    }
                    /* store the liability class */
                    pheno_listi->num_classes =
                        ((pheno_listi->num_classes < lclass ) ?
                         lclass : pheno_listi->num_classes) ;

                } else {
                    pheno_listi->num_classes = 1;
                }
            }
        }

        for (i=0; i < pheno_listi->num_classes; i++) {
            insert_into_class_list(&(pheno_listi->class_list), i+1);
        }
    }
    return;
}

/* NM: 7-15-08: set flag whether this marker needs to be recoded,
   8-21-08: For recode, we count everyone, then depending on the chromsome,
   count males once, females twice for X, 0 for Y.

*/

linkage_ped_top *count_allele_freq(linkage_ped_top *Top,
				    int locus,
				    int count_option,
				    marker_type *marker_listi,
                                    allelecnt **member_ids,
				    int count_ht)
{
    int ped, per, ped2, rare, rareped, pedcnt = Top->PedCnt;
    const char  *all1, *all2;
    int sex, entrycount, found, ht;
    record_type rec;
    linkage_locus_type LocType = Top->LocusTop->Locus[locus].Type;
    allele_list_type *all1fp, *all2fp;
    linkage_ped_rec  *pr = NULL, *prp = NULL;
    person_node_type *pn = NULL, *pnp = NULL;
    int homoz = 0, homozper = 0, rand = 0;
    int founder = 0;
    int typed1, typed2;
    int ht_ok  = count_ht ? 1 : 0;
    int xlinked = 0;

    allele_prop_reset();

    if (LocType == AFFECTION || LocType == QUANT) {
        marker_listi->first_allele = CALLOC((size_t) 1, allele_list_type);
        marker_listi->first_allele->next = CALLOC((size_t) 1, allele_list_type);
        marker_listi->first_allele->next->next = NULL;
    } else if ((LocType == NUMBERED) || (LocType == XLINKED) || (LocType == YLINKED)) {

        for (ped2=0; ped2 < 2 * pedcnt; ped2++) {
            if (ped2 < pedcnt) {
                rare = 0;
                ped = ped2;
            } else {
                rare = 1;
                ped = ped2-pedcnt;
            }
/*
            rare = (ped2 < pedcnt) ? 0 : 1;
            ped  = (ped2 < pedcnt) ? ped2 : ped2 - pedcnt;
*/
            if (UntypedPeds[ped] == 1) continue;

            if (pedfile_type == POSTMAKEPED_PFT) {
                entrycount = Top->Ped[ped].EntryCnt;
                pr = Top->Ped[ped].Entry;
                rec = Raw_postmake;
            } else {
                entrycount = Top->PTop[ped].num_persons;
                pn = Top->PTop[ped].persons;
                rec = Raw_premake;
            }
            homozper = 0;
            rand = 0;
            found = 0;
            rareped = 0;
            for(per=0; per <= entrycount && !rand; per++) {
                if (per == entrycount) {
                    rand = 1;
                    if (homozper == 0) {
                        continue;
                    } else if (found == 0 && UntypedPeds[ped] != 1) { //Untyped redundant
                        int idx = randomindex(homozper);
                        per = member_ids[ped][idx].num;
                    } else {
                        continue;
                    }
                }
                member_ids[ped][per].num=0;

            /* get alleles */
                if (pedfile_type == POSTMAKEPED_PFT) {
                    prp = &pr[per];
                    get_2Ralleles(prp->Marker, locus, &all1, &all2);
                    sex  = prp->Sex;
                    founder = IS_LFOUNDER(*prp);
                } else {
                    pnp = &pn[per];
                    get_2Ralleles(pnp->marker, locus, &all1, &all2);
                    sex  = pnp->gender;
                    founder = PFOUNDER(*pnp);
                }

                if (LocType == YLINKED && sex == 2) {
                    /* skip the female entirely */
                    continue;
                }

                if (!allelecmp(all1, REC_UNKNOWN)) {
                    const char *tmp = all1;
                    all1 = all2;
                    all2 = tmp;
                }
            /* get freq struct */
                all1fp = (allele_list_type *) allele2allele_prop_prop(all1);
                if (all1fp == (void *) 0) {
                    all1fp = insert_into_allele_list(&(marker_listi->first_allele), all1);
                    allele2allele_prop_prop(all1) = all1fp;
                }

                all2fp = (allele_list_type *) allele2allele_prop_prop(all2);
                if (all2fp == (void *) 0) {
                    all2fp = insert_into_allele_list(&(marker_listi->first_allele), all2);
                    allele2allele_prop_prop(all2) = all2fp;
                }

             /* typing */
                homoz = all1 == all2;
                if (rec == Raw_postmake) {
                    const char *a1, *a2;
                    get_2Ralleles(prp->Marker, locus, &a1, &a2);
                    typed1=allelecmp(a1, REC_UNKNOWN);
                    typed2=allelecmp(a2, REC_UNKNOWN);
                } else if (rec == Postmakeped) {
                    int a1, a2;
                    get_2alleles(prp->Marker, locus, &a1, &a2);
                    typed1=(a1)? 1 : 0;
                    typed2=(a2)? 1 : 0;
                } else if (rec == Raw_premake) {
                    const char *a1, *a2;
                    get_2Ralleles(pnp->marker, locus, &a1, &a2);
                    typed1=allelecmp(a1, REC_UNKNOWN);
                    typed2=allelecmp(a2, REC_UNKNOWN);
                } else {
                    int a1, a2;
                    get_2alleles(pnp->marker, locus, &a1, &a2);
                    typed1=(a1)? 1 : 0;
                    typed2=(a2)? 1 : 0;
                }

            /* find if this person is genotyped at the locus */
                if (typed1) {
                    ht = (typed2) ? 2 : ht_ok;
//                    has_typed++;
                } else {
                    ht = (typed2) ? ht_ok : 0;
//                    if (typed2) has_typed++;
                }
		// if it's not a male
                if (!(xlinked && sex == 1)) {
                    int homoz_test = ((LocType == YLINKED) && typed1 && typed2) ? homoz : 1;
                    if (homoz_test)
                        member_ids[ped][homozper].num = per;
                    homozper++;
                }

                if (ht >= 1) {

//says all1 & all2 have values that are different Heterozygote && YLINKED
                    if (allelecmp(all1, REC_UNKNOWN) && allelecmp(all2, REC_UNKNOWN) &&
                        allelecmp(all1, all2) && LocType == YLINKED) {
                        continue;
                    }
                    if (ht == 1) {
                        marker_listi->num_half_typed++;
                    }

                    /* Count the first allele for everyone if X-linked or autosomal,
                       only males, if Y-linked */
//                        count_this_allele(marker_listi->first_allele,
//                                          all1, 1);
                    if (rare) {
                        if (all1fp->allele_freq.random_count == 0) {
/*                            all1fp->allele_freq.rare_count++; */
                            rareped = 1;
                        }
                        if (all2fp->allele_freq.random_count == 0) {
/*                            all2fp->allele_freq.rare_count++; */
                            rareped = 1;
                        }
                        if (rareped) break; /* out if per loop */
                    } else {
                        if (founder) {
                            found=1;
                            marker_listi->num_founders++;
                            marker_listi->num_random++;
                            marker_listi->num_unique++;
                            if (ht == 1) {
                                marker_listi->ht_founders++;
                                marker_listi->ht_random++;
                                marker_listi->ht_unique++;
                            }
                            if (LocType == NUMBERED || LocType == XLINKED ||
                                (sex == 1 && LocType == YLINKED)) {
                                all1fp->allele_freq.founder_count++;
                                all1fp->allele_freq.random_count++;
                                all1fp->allele_freq.unique_count++;
                            }
                            if ((LocType == NUMBERED) || (sex == 2 && LocType == XLINKED)) {
                                all2fp->allele_freq.founder_count++;
                                all2fp->allele_freq.random_count++;
                                all2fp->allele_freq.unique_count++;
                            }
                        } else if (rand) {
                            marker_listi->num_random ++;
                            marker_listi->num_unique++;
                            if (ht == 1) {
                                marker_listi->ht_random++;
                                marker_listi->ht_unique++;
                            }
                            if (LocType == NUMBERED || LocType == XLINKED ||
                                (sex == 1 && LocType == YLINKED)) {
                                all1fp->allele_freq.random_count++;
                                all1fp->allele_freq.unique_count++;
                            }
                            if ((LocType == NUMBERED) || (sex == 2 && LocType == XLINKED)) {
                                all2fp->allele_freq.random_count++;
                                all2fp->allele_freq.unique_count++;
                            }
                        }

                        if (ht == 1)
                            marker_listi->ht_everyone++;
                    }
                }
            }
        }
    }
    /*
      xfer rare_count -> random_count
     */
    return Top;
}

//w
int Display_y_female, y_female = 0;
linkage_ped_top *create_allele_list(linkage_ped_top *Top,
				    int locus,
				    marker_type *marker_listi,
                                    allelecnt **member_ids,
				    int count_ht)
{
    int ped, per;
    const char  *all1, *all2;
    int sex, entrycount;
//  record_type rec;
    linkage_locus_type LocType = Top->LocusTop->Locus[locus].Type;
    allele_list_type *al1p, *al2p;

    Display_Errors = 1;
    allele_prop_reset();
    if (LocType == AFFECTION || LocType == QUANT) {
//      marker_listi->first_allele = CALLOC((size_t) 1, allele_list_type);
//      marker_listi->first_allele->next = CALLOC((size_t) 1, allele_list_type);
//      marker_listi->first_allele->next->next = NULL;
//      For Traits, just put the answers in during recode phase
    } else if ((LocType == NUMBERED) || (LocType == XLINKED) || (LocType == YLINKED)) {

        /* first, for this marker, create the list of all new alleles,
           and set the count to 0 */

        marker_listi->num_people=0;
        marker_listi->num_half_typed=0;

        marker_listi->num_founders = marker_listi->num_random = 0;
        marker_listi->num_unique   = marker_listi->num_everyone = 0;

        marker_listi->ht_founders = marker_listi->ht_random = 0;
        marker_listi->ht_unique   = marker_listi->ht_everyone = 0;

        if (pedfile_type == POSTMAKEPED_PFT) {
//          rec = Raw_postmake;
            allelecnt **mp = member_ids;
            for (ped=0; ped < Top->PedCnt; ped++, mp++) {
                entrycount = Top->Ped[ped].EntryCnt;
                linkage_ped_rec *te = &Top->Ped[ped].Entry[0];
                allelecnt *mpp = mp[0];
                for(per=0; per < entrycount; per++, te++, mpp++) {
                    get_2Ralleles(te->Marker, locus, &all1, &all2);
                    mpp->num  = 0;
                    mpp->all1 = all1;
                    mpp->all2 = all2;
                    if (!allelecmp(all1, REC_UNKNOWN)) {
                        if (!allelecmp(all2, REC_UNKNOWN)) {
                            continue; // 0/0
                        }
                        const char *tmp = all1;
                        all1 = all2;
                        all2 = tmp;
                    }
                    al1p = (allele_list_type *)allele2allele_prop_prop(all1);
                    if (al1p == (void *) 0) {
                        al1p = insert_into_allele_list(&(marker_listi->first_allele),
                                                      all1);
                        allele2allele_prop_prop(all1) = al1p;
                    }
                    al2p = (allele_list_type *)allele2allele_prop_prop(all2);
                    if (al2p == (void *) 0) {
                        al2p = insert_into_allele_list(&(marker_listi->first_allele),
                                                      all2);
                        allele2allele_prop_prop(all2) = al2p;
                    }
                    if (Top->LocusTop->Locus[locus].number == -1) continue;
                    if (!allelecmp(all2, REC_UNKNOWN)) {
                        if (!count_ht)
                            continue; // if no half type allowed
                        marker_listi->ht_everyone++;
                    }
                    sex = te->Sex;
                    if (LocType == NUMBERED) {
                        marker_listi->num_everyone++;  // one or the other must be != REC_UNKNOWN
                        if (al1p) al1p->allele_freq.everyone_count++;
                        if (al2p) al2p->allele_freq.everyone_count++;
                    } else if (sex == 1) {
                        if (al1p && (al1p == al2p)) {
                            marker_listi->num_everyone++;
                            al1p->allele_freq.everyone_count++;
                        }
                    } else if (sex == 2) {
                        if (LocType == XLINKED) {
                            marker_listi->num_everyone++;
                            if (al1p) al1p->allele_freq.everyone_count++;
                            if (al2p) al2p->allele_freq.everyone_count++;
                        } else if (LocType == YLINKED) {
                            if (allelecmp(all1, REC_UNKNOWN) || allelecmp(all2, REC_UNKNOWN)) {
                                SUPPRESS_MSSG_NESTED(y_female);
                                warnvf("marker %s (on Y chromosome) observed for female (%s) with value%s %s/%s\n",
                                       Top->LocusTop->Locus[locus].Name,
                                       Top->Ped[ped].Entry[per].UniqueID,
                                       all1 != all2 ? "(s)" : "",
                                       all1,
                                       all1 != all2 ? all2 : "");
                            }
                        }
                    }

                }
            }
        } else {
//          rec = Raw_premake;
            allelecnt **mp = member_ids;
            for (ped=0; ped < Top->PedCnt; ped++, mp++) {
                entrycount = Top->PTop[ped].num_persons;
                person_node_type *tp = &Top->PTop[ped].persons[0];
                allelecnt *mpp = mp[0];
                for(per=0; per < entrycount; per++, tp++, mpp++) {
                    get_2Ralleles(tp->marker, locus, &all1, &all2);
                    mpp->num  = 0;
                    mpp->all1 = all1;
                    mpp->all2 = all2;
                    if (!allelecmp(all1, REC_UNKNOWN)) {
                        if (!allelecmp(all2, REC_UNKNOWN)) {
                            continue; // 0/0
                        }
                        const char *tmp = all1;
                        all1 = all2;
                        all2 = tmp;
                    }
                    al1p = (allele_list_type *)allele2allele_prop_prop(all1);
                    if (al1p == (void *) 0) {
                        al1p = insert_into_allele_list(&(marker_listi->first_allele),
                                                      all1);
                        allele2allele_prop_prop(all1) = al1p;
                    }
                    al2p = (allele_list_type *)allele2allele_prop_prop(all2);
                    if (al2p == (void *) 0) {
                        al2p = insert_into_allele_list(&(marker_listi->first_allele),
                                                      all2);
                        allele2allele_prop_prop(all2) = al2p;
                    }
                    if (Top->LocusTop->Locus[locus].number == -1) continue;
                    if (!allelecmp(all2, REC_UNKNOWN)) {
                        if (!count_ht)
                            continue; // if no half type allowed
                        marker_listi->ht_everyone++;
                    }
                    sex = tp->gender;
                    if (LocType == NUMBERED) {
                        marker_listi->num_everyone++;  // one or the other must be != REC_UNKNOWN
                        if (al1p) al1p->allele_freq.everyone_count++;
                        if (al2p) al2p->allele_freq.everyone_count++;
                    } else if (sex == 1) {
                        if (al1p && (al1p == al2p)) {
                            marker_listi->num_everyone++;
                            al1p->allele_freq.everyone_count++;
                        }
                    } else if (sex == 2) {
                        if (LocType == XLINKED) {
                            marker_listi->num_everyone++;
                            if (al1p) al1p->allele_freq.everyone_count++;
                            if (al2p) al2p->allele_freq.everyone_count++;
                        } else if (LocType == YLINKED) {
                            if (allelecmp(all1, REC_UNKNOWN) || allelecmp(all2, REC_UNKNOWN)) {
                                SUPPRESS_MSSG_NESTED(y_female);
                                warnvf("marker %s (on Y chromosome) observed for female (%s) with value%s %s/%s\n",
                                       Top->LocusTop->Locus[locus].Name,
                                       Top->PTop[ped].persons[per].uniqueid,
                                       all1 != all2 ? "(s)" : "",
                                       all1,
                                       all1 != all2 ? all2 : "");
                            }
                        }
                    }
                }
            }
        }
        Display_Errors = 1;
    }
    return Top;
}

linkage_ped_top *count_allele_list(linkage_ped_top *Top,
                                   int locus,
                                   int count_option,
                                   marker_type *marker_listi,
                                   allelecnt **member_ids,
                                   int count_ht)
{
    int ped, per;
    const char  *all1, *all2;
    int sex, entrycount, found, counted1, counted2, ht;
    int genotyped;
    record_type rec;
    allele_list_type *allelep;
    linkage_locus_type LocType = Top->LocusTop->Locus[locus].Type;
//  allele_list_type *alp;

    if ((LocType == NUMBERED) || (LocType == XLINKED) || (LocType == YLINKED)) {
        marker_listi->num_people=0;
        marker_listi->num_half_typed=0;
        marker_listi->num_founders =  marker_listi->num_random =
            marker_listi->num_unique = marker_listi->num_everyone =
            marker_listi->ht_founders = marker_listi->ht_random =
            marker_listi->ht_unique = marker_listi->ht_everyone = 0;

        /* Now create the count for each allele */

        /* First founders */
//      Tod tod_counting("founders, everyone, count this");
        rec = (pedfile_type == POSTMAKEPED_PFT) ? Raw_postmake : Raw_premake;

        for (ped=0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds[ped] == 1) {
                continue;
            }
            found=0;
            entrycount = ((pedfile_type == POSTMAKEPED_PFT)? Top->Ped[ped].EntryCnt :
                          Top->PTop[ped].num_persons);
            founders_or_everyone(Top, locus, rec, ped,
                                 member_ids[ped],
                                 1, count_ht, 0);

            for(per=0; per < entrycount; per++) {
                if (member_ids[ped][per].num >= 1) {
                    sex = (pedfile_type == POSTMAKEPED_PFT) ?
                               Top->Ped[ped].Entry[per].Sex :
                               Top->PTop[ped].persons[per].gender;
                    all1 = member_ids[ped][per].all1;
                    all2 = member_ids[ped][per].all2;

                    /* printf("1: %d %d %s(%p) %s(%p) %d\n",
                              ped, per, all1, all1, all2, all2, sex); */

                    /* Assumes that female genotypes have been set to unknown,
                       skip heterozygous males */
                    if (allelecmp(all1, REC_UNKNOWN) && allelecmp(all2, REC_UNKNOWN) &&
                        allelecmp(all1, all2) && LocType == YLINKED) {
                        continue;
                    }

                    if (member_ids[ped][per].num == 1) {
                        /* switch the two alleles if the 1st one is unknown,
                           so that we always count the known allele */
                        /* RETURN VALUES The strcmp() and strncmp() return an integer greater than, equal to, or
                           less than 0, according as the string s1 is greater than, equal to, or less than the string s2.

                           Logical expressions connected by && and || have value 1 if true and 0 if false.

                           (!allelecmp(all1, REC_UNKNOWN) is equivalent to (all1 == REC_UNKNOWN)

                        */
                        if (!allelecmp(all1, REC_UNKNOWN)) {
                            const char *tmp = all1;
                            all1 = all2;
                            all2 = tmp;
/*
                            strcpy(all1, all2);
                            strcpy(all2, REC_UNKNOWN);
*/
                        }
                    }

                    /* Count the first allele for everyone if X-linked or autosomal,
                       only males, if Y-linked */
                    if (LocType == NUMBERED || LocType == XLINKED ||
                        (sex == 1 && LocType == YLINKED)) {
                        count_this_allele(marker_listi->first_allele,
                                          all1, 1);
                        found=1;
                        marker_listi->num_founders ++;
                        marker_listi->num_random ++;
                        marker_listi->num_unique ++;

                        if (member_ids[ped][per].num == 1) {
                            (marker_listi->ht_founders)++;
                            (marker_listi->ht_random)++;
                            (marker_listi->ht_unique)++;
                        }
                    }
                    /* Count the second allele for everyone if autosomal,
                       only females if X-linked */
                    if ((LocType == NUMBERED) || (sex == 2 && LocType == XLINKED)) {
                        count_this_allele(marker_listi->first_allele,
                                          all2, 1);
                    }
                }
            }

            /* Next the random */
            if (found == 0 && UntypedPeds[ped] != 1) {
                /* make the selection so that the person returned is always counted */
                if (LocType == XLINKED) {
                    /* For x-linked loci, make sure that a female is selected */
                    per = random_ped_member(Top, locus, ped, rec,
                                            count_ht, &genotyped, 1);
                } else {
                    /* If y-linked, females are not typed, hence, not selected */
                    per = random_ped_member(Top, locus, ped, rec,
                                            count_ht, &genotyped, 0);
                }

                if (per > -1 && genotyped >= 1) {
                    member_ids[ped][per].num=genotyped;
                    /* This can be incremented here because we always count the
                       person index returned */
                    (marker_listi->num_random)++;
                    (marker_listi->num_unique)++;

                    sex = (pedfile_type == POSTMAKEPED_PFT) ?
                               Top->Ped[ped].Entry[per].Sex :
                               Top->PTop[ped].persons[per].gender;
                    all1 = member_ids[ped][per].all1;
                    all2 = member_ids[ped][per].all2;

                    /* printf("2: %d %d %s(%p) %s(%p) %d\n",
                              ped, per, all1, all1, all2, all2, sex); */

                    if (member_ids[ped][per].num == 1) {
                        /* switch the two alleles if the 1st one is unknown,
                           so that we always count the known allele */
                        if (!allelecmp(all1, REC_UNKNOWN)) {
                            const char *tmp = all1;
                            all1 = all2;
                            all2 = tmp;
/*
                            strcpy(all1, all2);
                            strcpy(all2, REC_UNKNOWN);
*/
                        }
                    }

                    /* Count first allele for NUMBERED and XLINKED for everyone,
                       YLINKED for males.
                       Here, count the first allele from males at YLINKED loci if
                       they are coded as homozygous.  But how are YLINKED loci
                       entered into the program? */
                    /*
                     *     (allelecmp(all1, REC_UNKNOWN)) is equivalent to (all1 != REC_UNKNOWN)
                     *
                     */
                    if (LocType == XLINKED || LocType == NUMBERED ||
                        (LocType == YLINKED && sex == 1  &&
                         !((allelecmp(all1, REC_UNKNOWN)) && (allelecmp(all2, REC_UNKNOWN)) && allelecmp(all1, all2)))) {
                        count_this_allele(marker_listi->first_allele,
                                          all1, 2);
                    }

                    /* Count the second alleles for numbered and Xlinked-female */
                    if (LocType == NUMBERED || (LocType == XLINKED && sex == 2)) {
                        count_this_allele(marker_listi->first_allele,
                                          all2, 2);
                    }

                    if (genotyped == 1) {
                        (marker_listi->ht_random)++;
                        (marker_listi->ht_unique)++;
                    }
                }
            }
        }
//      tod_counting();

        /* unique alleles */
//      Tod tod_alleles("while thru alleles and ped & halftype");
        allelep = marker_listi->first_allele;
        while (allelep != NULL) {
            /* If there are still alleles with 0 counts
               insert non-founders with new alleles */
            if (allelep->allele_freq.random_count > 0) {
                allelep = allelep->next;
                continue;
            }
            for (ped=0; ped < Top->PedCnt; ped++) {
                if (UntypedPeds[ped] == 1) {
                    continue;
                }
                /* Else insert at most one person with that allele
                   from each pedigree */
                entrycount =
                    ((pedfile_type == POSTMAKEPED_PFT) ?
                     Top->Ped[ped].EntryCnt : Top->PTop[ped].num_persons);

                for(per=0; per < entrycount; per++) {
                    sex = (pedfile_type == POSTMAKEPED_PFT) ?
                               Top->Ped[ped].Entry[per].Sex :
                               Top->PTop[ped].persons[per].gender;
                    all1 = member_ids[ped][per].all1;
                    all2 = member_ids[ped][per].all2;

                    /* printf("3: %d %d %s(%p) %s(%p) %d\n",
                              ped, per, all1, all1, all2, all2, sex); */
                    ht=counted1=counted2=0;

                    /* ignore females' genotypes */
                    if (LocType == YLINKED && sex == 2) {
                        continue;
                    }

                    if (!allelecmp(all1, REC_UNKNOWN) && !allelecmp(all2, REC_UNKNOWN)) {
                        continue;
                    }
                    if (!allelecmp(all1, REC_UNKNOWN) || !allelecmp(all2, REC_UNKNOWN)) {
                        ht=1;
                    }
                    if (ht && !count_ht) {
                        continue;
                    }

                    /* Check that males are homozygous for Xlinked and Ylinked loci, if
                       both alleles are set to known */
                    if (sex == 1 && ht == 0 &&
                       (LocType == XLINKED || LocType == YLINKED) &&
                        allelecmp(all1, all2)) {
                        continue;
                    }

                    /* For XLINKED or YLINKED loci, and a male, count only one known allele,
                       Females are automatically skipped for Ylinked loci */
                    if (sex == 1 &&
                        (LocType == XLINKED || LocType == YLINKED)) {
                        if (!allelecmp(all1, allelep->allele_freq.name)) {
                            (allelep->allele_freq.unique_count) ++;
                            counted1=counted2=1;
                        } else if (!allelecmp(all2, allelep->allele_freq.name)) {
                            (allelep->allele_freq.unique_count) ++;
                            counted1 = counted2 = 1;
                        }
                        continue;
                    }

                    /* For NUMBERED and everyone and XLINKED and females */

                    if (LocType == NUMBERED || (LocType == XLINKED && sex == 2)) {
                        /* if allele1 is unique, increment allele1, then also insert allele2 */
                        if (!allelecmp(all1, allelep->allele_freq.name)) {
                            (allelep->allele_freq.unique_count) ++;
                            counted1=1;
                            if (allelecmp(all2, REC_UNKNOWN)) {
                                count_this_allele(marker_listi->first_allele,
                                                  all2, 3);

                                member_ids[ped][per].num=2;
                                counted2=1;
                            } else {
                                member_ids[ped][per].num=1;
                            }
                        } else if (!allelecmp(all2, allelep->allele_freq.name)) {
                            /* count allele2 as unique and also insert allele1 */
                            (allelep->allele_freq.unique_count) ++;
                            counted2=1;
                            if (allelecmp(all1, REC_UNKNOWN)) {
                                count_this_allele(marker_listi->first_allele,
                                                  all1, 3);

                                member_ids[ped][per].num=2;
                                counted1 = 1;
                            } else {
                                member_ids[ped][per].num=1;
                            }
                        }
                    }

                    if (counted2 || counted1) {
                        marker_listi->num_unique ++;
                        if (member_ids[ped][per].num == 1) {
                            (marker_listi->ht_unique) ++;
                        }
                        /* This breaks out of the persons loop and goes to next pedigree */
                        break;
                    }
                }
            }
            allelep = allelep->next;
        }
//      tod_alleles();
//      Tod tod_lt4("count_option < 4");
        if (count_option < 4) {
            /* count up the half_typed people for this marker now */
            for(ped=0; ped < Top->PedCnt; ped++) {
                if (UntypedPeds[ped] == 1) continue;

                entrycount =
                    ((pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt : Top->PTop[ped].num_persons);
                for(per=0; per < entrycount; per++) {
                    marker_listi->num_half_typed += ((member_ids[ped][per].num == 1)? 1: 0);
                }
            }
        }
//      tod_lt4();
        /* else leave it for later */

//      Tod tod_more("founder or everyone");
        /* count everyone after setting all member ids to 0 */
        for (ped=0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds[ped] == 1) continue;
            /* Else insert at most one person with that allele
               from each pedigree */
            entrycount =
                ((pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt : Top->PTop[ped].num_persons);
            for(per=0; per < entrycount; per++) {
                member_ids[ped][per].num=0;
            }
            founders_or_everyone(Top, locus, rec, ped,
                                 member_ids[ped],
                                 4, count_ht, 0);

            for(per=0; per < entrycount; per++) {
                if (member_ids[ped][per].num >= 1) {
                    sex = (pedfile_type == POSTMAKEPED_PFT) ?
                               Top->Ped[ped].Entry[per].Sex :
                               Top->PTop[ped].persons[per].gender;
                    all1 = member_ids[ped][per].all1;
                    all2 = member_ids[ped][per].all2;

                    /* printf("4: %d %d %s(%p) %s(%p) %d\n",
                              ped, per, all1, all1, all2, all2, sex); */

                    if (member_ids[ped][per].num == 1) {
                        /* This is a half-typed person */
                        if (count_ht == 0) {
                            continue;
                        } else {
                            /* else switch the two alleles if the 1st one is unknown,
                               so that we always count the known allele */
                            if (!allelecmp(all1, REC_UNKNOWN)) {
                                const char *tmp = all1;
                                all1 = all2;
                                all2 = tmp;
/*
                                strcpy(all1, all2);
                                strcpy(all2, REC_UNKNOWN);
*/
                            }
                        }
                    }

                    if (sex == 1 && member_ids[ped][per].num == 2 && allelecmp(all1, all2) &&
                        (LocType == XLINKED || LocType == YLINKED)) {
                        continue;
                    }

                    if (LocType == NUMBERED || LocType == XLINKED ||
                        (LocType == YLINKED && sex == 1)) {
                        count_this_allele(marker_listi->first_allele,
                                          all1, 4);
                        marker_listi->num_everyone ++;
                    }

                    if ((LocType == NUMBERED) || (LocType == XLINKED && sex == 2)) {
                        count_this_allele(marker_listi->first_allele,
                                          all2, 4);
                    }
                    if (member_ids[ped][per].num == 1) {
                        (marker_listi->ht_everyone) ++;
                    }
                }
            }
        }
//      tod_more();
//      Tod tod_eq4("eq 4");
        if (count_option == 4) {
            /* count up the half_typed people for this marker */
            for(ped=0; ped < Top->PedCnt; ped++) {
                if (UntypedPeds[ped] == 1) continue;
                entrycount =
                    ((pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt : Top->PTop[ped].num_persons);
                for(per=0; per < entrycount; per++) {
                    marker_listi->num_half_typed += ((member_ids[ped][per].num == 1)? 1: 0);
                }
            }
        }
//      tod_eq4();
    }

//m    count_allele_freq(Top, locus, count_option, marker_listi, member_ids, count_ht);

    return Top;
}

extern void set_missing_quant_input(linkage_ped_top *Top, analysis_type analysis);

linkage_ped_top  *create_full_marker_data(
    FILE *pedfilep,
    char *omitfl_name,
    int *untyped_ped_opt,
    linkage_locus_top *LTop,
    int *col2locus,
    analysis_type analysis)
{
    FILE *fp=pedfilep;
    int i, ped_count, numlines;
    int count_option, count_halftyped;
    marker_type *marker_list;
    pheno_type *pheno_list;
    linkage_ped_top *Top;

    Mega2Status = INSIDE_RECODE;
    /*  read_marker_data(lfilep, &LTop); */

    if (LTop->PhenoCnt > 0) {
        pheno_list=CALLOC((size_t) LTop->PhenoCnt, pheno_type);
        for (i=0; i < LTop->PhenoCnt; i++) {
            pheno_type *pheno_listi = &(pheno_list[i]);
            pheno_listi->first_allele = NULL;
            pheno_listi->estimate_frequencies = 1;
            pheno_listi->first_allele = CALLOC((size_t) 1, allele_list_type);
            pheno_listi->first_allele->allele_freq.freq = 0.5;
            pheno_listi->first_allele->next = CALLOC((size_t) 1, allele_list_type);
            pheno_listi->first_allele->next->allele_freq.freq = 0.5;
            pheno_listi->first_allele->next->next = NULL;
            if (LTop->Locus[i].Type == QUANT) {
                pheno_listi->num_classes = 1;
                pheno_listi->num_alleles = 2;
            } else if (LTop->Locus[i].Type == AFFECTION) {
                pheno_listi->num_classes = 0;
                pheno_listi->num_alleles = 2;
            }
        }
    } else 
        pheno_list = 0;

    if (LTop->MarkerCnt > 0) {
        marker_list = (CALLOC((size_t) LTop->MarkerCnt, marker_type)) - LTop->PhenoCnt;
        for(i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
            marker_list[i].first_allele = NULL;
            if (LTop->Locus[i].Type == NUMBERED || LTop->Locus[i].Type == XLINKED ||
               LTop->Locus[i].Type == YLINKED) {

                marker_list[i].num_alleles = 0;
                marker_list[i].num_classes = 0;

            }
            marker_list[i].estimate_frequencies = 1;
            marker_list[i].recode_alleles = 0;
        }
    } else
        marker_list = 0;

    set_missing_quant_input((linkage_ped_top *) NULL, analysis);

    ped_count = check_pre_makeped(fp, &numlines);
    if (ped_count) {
        pedfile_type = PREMAKEPED_PFT;
        LTop->PedRecDataType = Raw_premake;
        rewind(fp);
        Top = read_pre_makeped(fp, ped_count, numlines, NULL, LTop, col2locus);
    } else {
        pedfile_type = POSTMAKEPED_PFT;
        LTop->PedRecDataType = Raw_postmake;
        rewind(fp);
        Top = read_linkage_ped_file(fp, LTop, col2locus);
    }
    log_line(mssgf);
    mssgf("Input pedigree data contains:");

    order_heterozygous_allele(Top);
    write_ped_stats(Top, pedfile_type);
    if (omitfl_name != NULL) {
        premakeped_omit_file(Top, omitfl_name,1);
    }

    if (HasMarkers) {
        count_option =
            get_count_option(1, &count_halftyped,
                             "Select individuals to compute allele frequencies\n       for recoded marker loci:");
    } else {
        count_option=4;
    }

/* omit_peds returns currently VOID so the code below ... */
/*   if (omit_peds(*untyped_ped_opt, Top) == 0) { */
/*     mssgf("terminating Mega2!"); */
/*     EXIT(EARLY_TERMINATION); */
/*   } */
    if (UntypedPeds == NULL) {
        UntypedPeds = CALLOC((size_t) Top->PedCnt, int);
    }
    omit_peds(*untyped_ped_opt, Top);
    if (analysis != TO_PLINK) {
        mssgf("Counting alleles and computing frequencies ...");
    } else {
        mssgf("Counting alleles for each marker ...");
    }
    {
        int ped, entrycount;
        allelecnt **member_ids;
        member_ids=CALLOC((size_t) Top->PedCnt, allelecnt *);
        for (ped=0; ped < Top->PedCnt; ped++) {
            entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt:
                                               Top->PTop[ped].num_persons;
            member_ids[ped]=CALLOC((size_t) entrycount, allelecnt);
        }

        Tod tod_cal1("create allele list");
        Tod tod_fr1(20);
        Tod tod_fr1_x4(20); // 0.000,688
        for (i=0; i < LTop->LocusCnt; i++) {
            tod_fr1.reset();
            tod_fr1_x4.reset();
#ifdef DEBUG
            if (((i+1) % 100) == 0) {
                printf("Marker %s ..", LTop->Locus[i].Name);
            }
#endif
        /* in this function, we also decide whether to recode or not */
            if (LTop->Locus[i].Class == MARKER) {
                create_allele_list(Top, i, &(marker_list[i]),
                                   member_ids, count_halftyped);
            
                tod_fr1("create allele list 1 freq");
                if ( (count_option != 4 || analysis == TO_HWETEST || analysis == TO_SIMULATE) &&
                      LTop->Locus[i].number != -1 ) {
                    count_allele_list(Top, i, count_option, &(marker_list[i]),
                                      member_ids, count_halftyped);
                    tod_fr1_x4("create allele list 4 freq");
                }
            } else
                default_trait_penetrances(Top, i, &(pheno_list[i]));
        }
        for(ped=0; ped < Top->PedCnt; ped++) {
            free(member_ids[ped]);
        }
        free(member_ids);
        tod_cal1();
    }
    SUPPRESS_MSSG_NESTED_FINI(y_female);

    count_raw_alleles(marker_list, Top->LocusTop);

    if (analysis != TO_PLINK) {
        count_classes(pheno_list, Top->LocusTop);
    }

    /* see if alleles are numeric and consecutive, else
       needs recoding, also sort alleles */

    for (i=0; i < LTop->LocusCnt; i++) {
        if (LTop->Locus[i].Type == NUMBERED ||
            LTop->Locus[i].Type == XLINKED ||
            LTop->Locus[i].Type == YLINKED) {

            need_recoding(&(marker_list[i]));

/*       if (LTop->Locus[i].Type == XLINKED || */
/*    	     LTop->Locus[i].Type == YLINKED) { */
/* 	        LTop->Locus[i].Type = NUMBERED; */
/*       } */
        }
    }
    convert_to_freq(marker_list, Top->LocusTop, count_option, analysis);
    assign_dummy_alleles(marker_list, Top->LocusTop);

    write_recode_summary(marker_list, Top->LocusTop, count_halftyped);
    recode_ped_top(marker_list, Top, (plink_info_type *)NULL);
    recode_locus_top(marker_list, pheno_list, Top->LocusTop);
    Mega2Status = DONE_RECODE;
    log_line(mssgf);
    mssgf("Pedigree data summary after recoding:");
    write_ped_stats(Top, pedfile_type);

    for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
        free_marker_item(marker_list[i].first_allele);
    }

    if (marker_list) {
        free(marker_list + LTop->PhenoCnt);
        marker_list=NULL;
    }

    return Top;
}

/* Not in names format if we cannot read
   'M', 'X', 'L', 'A', 'C', or 'T'
*/

//
// Read the first non-whitespace character in the file
// If it is a digit assume that it is in linkage (LINKAGE) format.
// Otherwise, assume that the file is in Mega2 (NAMES) format.
file_format check_locus_file_format(FILE *fp)
{
    char first_char;

    (void)fscanf(fp, "%c", &first_char);
    while(isspace((int)first_char)) {
        (void)fscanf(fp, "%c", &first_char);
    }
    rewind(fp);

    if (isdigit((int)first_char)) {
#ifndef HIDESTATUS
        mssgf("Locus file is in Linkage format");
#endif
        return LINKAGE;
    } else {
#ifndef HIDESTATUS
        mssgf("Locus file is \"Names\" format");
#endif
        return NAMES;
    }
}


#ifdef DEBUG

static void aff_lines(FILE *fp, int nclass, char *name)

{
    int i;

    fprintf(fp, "1   2  # %s \n", name);
    fprintf(fp, "0.500000 0.500000   << GENE FREQUENCIES\n");
    fprintf(fp, "%d  << NO. OF LIABILITY CLASSES\n", nclass);
    for (i=0; i < nclass; i++) {
        fprintf(fp, "0.0000 0.0000 1.0000 << PENETRANCES \n");
    }
    return;
}

static void quant_lines(FILE *fp, char *name)

{

    fprintf(fp, "0   2  # %s\n", name);
    fprintf(fp, "0.500000 0.500000   # GENE FREQUENCIES\n");
    fprintf(fp, "1 # NO. OF TRAITS\n");
    fprintf(fp, "0.000  0.000  0.000 # GENOTYPE MEANS\n");
    fprintf(fp, "1.000  # VARIANCE - COVARIANCE MATRIX\n");
    fprintf(fp, "1.000  # MULTIPLIER FOR VARIANCE IN HETEROZYGOTES\n");
    return;

}

static void marker_lines(FILE *fp, marker_type *marker)
{
    int i;
    allele_list_type *allelep = marker->first_allele;

    fprintf(fp, "3   %d # %s\n", marker->num_alleles,
            marker->marker_name);
    while(allelep != NULL) {
        fprintf(fp, " %6.5f", allelep->allele_freq.freq);
        allelep=allelep->next;
    }
    fprintf(fp, "\n");
}

int main(int argc, char *argv[]) {

    char names_file[FILENAME_LENGTH];
    char ped_file[FILENAME_LENGTH];
    char output_file[FILENAME_LENGTH];
    char *types, **names;
    int num_markers;

    int lines;

    strcpy(names_file, argv[1]);
    strcpy(ped_file, argv[2]);
    strcpy(output_file, argv[3]);

    lines=linecount(names_file);
    types = CALLOC((size_t) lines, char);
    names = CALLOC((size_t) lines, char *);

    read_marker_data(names_file, types, names, &num_markers);
    read_pedigree_data(ped_file, types, names,
                       num_markers, output_file);

    /* free(names);*/

    /* print the header, not sex-linked for now */

    fprintf(fp, "%d  0  0 5\n", num_markers);
    fprintf(fp, " 0  0.0  0.0  0\n");
    for (m=1; m <= num_markers; m++) {
        fprintf(fp, " %d", m);
    }
    fprintf(fp, "\n");

    /* Now write the markers in order */
    for (m=0; m < num_markers; m++) {
        switch(types[m]) {
        case 'L':
        case 'A':
            aff_lines(fp, marker_list[m].num_classes, marker_list[m].marker_name);
        break;

        case 'M':
        case 'X':
            marker_lines(fp, &(marker_list[m]));
        break;

        case 'C':
        case 'T':
            quant_lines(fp, marker_list[m].marker_name);
        break;

        default:
            break;
        }
    }

    /* Now for the footer - dummy values for now */

    fprintf(fp, "  0  0\n");
    for (m=0; m < num_markers; m++) {
        fprintf(fp, " 0.0000");
    }
    fprintf(fp, "\n");
    fprintf(fp, "  1  0.05  0.5\n");

    fclose(fp);
    free(types);

}

#endif

// //////////////////////////////////////////////////////////////// //

/* plink parameters */
extern char *PLINKtrait;
/* end plink parameters */

linkage_locus_top *read_marker_only_data(FILE *fp, int cols, char **phe_names, int *phe_types)
{
    int i, num_markers, annotated=0;
    char nextline[FILENAME_LENGTH];
    char ltype[10], marker[100];
    char *types, **names;
    int penetrances_read = 0;
    int num_read, line_num;
    /* error flags */
    int invalid_marker_type =0, missing_columns = 0;
    int displayed_messages = 0;
    int loci_name_space = 0;

    num_markers = linecount(fp);
    num_markers += cols;

    types = CALLOC((size_t) num_markers, char);
    names = CALLOC((size_t) num_markers, char *);
    num_markers = 0;

    /* annotated format */
    annotated=1;
    line_num=0;

    for (i = 0; i < cols; i++) {
        if (phe_types[i] == 0)
            strcpy(ltype, "A");
        else
            strcpy(ltype, "T");
        types[i]=ltype[0];

        names[i] = CALLOC(strlen(phe_names[i])+1, char);
        strcpy(names[i], phe_names[i]);
    }
    i--;

    Display_Errors = 1;
    while(!feof(fp)) {
        *nextline = 0;
        (void)fgets(nextline, FILENAME_LENGTH - 1, fp);
        if (*nextline == 0) {
            break;
        }
        if (nextline[0] == '#' || nextline[0] == '\n') {
            continue;
        }
        num_read=sscanf(nextline, "%s %s", ltype, marker);
        line_num++;
        if (num_read > 0 && num_read < 2) {
            sprintf(err_msg,
                    "Line %d: Found fewer than 2 columns (at least 2 required).",
                    line_num);
            SUPPRESS_MSSG(displayed_messages);
            if (displayed_messages > MAX_PED_ERRORS) {
                Display_Errors = 0;
            }
            errorf(err_msg);
            displayed_messages++;
            missing_columns++;
        } else {
            if (!strcasecmp(ltype, SEX_CHROMOSOME_STR) || !strcasecmp(ltype, "X"))
                strcpy(ltype, "X");
            else if (!strcasecmp(ltype, MALE_CHROMOSOME_STR) || !strcasecmp(ltype, "Y"))
                strcpy(ltype, "Y");
            else
                strcpy(ltype, "M");
            if ((!strcmp(ltype, "M")) || (!strcmp(ltype, "A")) ||
               (!strcmp(ltype, "L")) || (!strcmp(ltype, "C")) ||
               (!strcmp(ltype, "X")) || (!strcmp(ltype, "Y")) ||
               (!strcmp(ltype, "T"))) {
                i++;
                types[i]=ltype[0];
                names[i] = CALLOC(strlen(marker)+1, char);
                loci_name_space += strlen(marker)+1;
                if ((!strcmp(ltype, "M")) || (!strcmp(ltype, "X")) || (!strcmp(ltype, "Y")) )
                    num_markers++;
                strcpy(names[i], marker);
            } else {
                sprintf(err_msg,
                        "Marker type %s for marker %s not recognized.",
                        ltype, marker);
                SUPPRESS_MSSG(displayed_messages);
                if (displayed_messages > MAX_PED_ERRORS) {
                    Display_Errors = 0;
                }
                errorf(err_msg);
                displayed_messages++;
                invalid_marker_type++;
            }
        }
    }
    Display_Errors = 1;

    if (missing_columns) {
        printf("Found lines with less than 2 columns, see %s for details.\n",
               Mega2Err);
    }

    if (invalid_marker_type) {
        printf("Invalid marker types, see %s for details.\n", Mega2Err);

    }

    if (missing_columns || invalid_marker_type) {
        errorf("Found fatal errors in locus file.");
        errorf("Please fix errors and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: locus strings: %d MB (%d x %d)\n",
          loci_name_space / 1024 / 1024,
          i+1,  loci_name_space / (i+1) );
#endif

    return read_common_marker_data(i+1, num_markers, names, types, annotated, penetrances_read, 0, NULL);
}


//
// Process the genotype 'map' and phenotype 'phe_*' marker data.
void read_m2_map_as_names_file(m2_map map,
                               linkage_locus_top **LTop,
                               const int cols,
                               char **phe_names,
                               int *phe_types)
{
    int i, num_markers = 0, annotated = 1; /* annotated format */
    char ltype[10];
    int penetrances_read = 0;
    unsigned int map_i = 0;
    int loci_name_space = 0;
    char *types, **names;
    size_t names_vector_size = map.size() + cols;

    types = CALLOC(names_vector_size, char);
    names = CALLOC(names_vector_size, char *);

    // The phenotypes are loaded in the beginning of the arrays....
    for (i = 0; i < cols; i++) {
        if (phe_types[i] == 0)
            strcpy(ltype, "A"); // Affection Status
        else
            strcpy(ltype, "T"); // Trait

        types[i]=ltype[0];
        names[i] = CALLOC(strlen(phe_names[i])+1, char);

        strcpy(names[i], phe_names[i]);
    }
    i--;

    Display_Errors = 1;
    // The genotypes are loaded at the end....
    // This is a superset of markers from a Mega2 map file. Because this fakes reading from a
    // names file, the marker information is split into two vectors 'names' and 'types'. We must
    // fill both in here as well as give a tally of the markers...
    for (map_i = 0; map_i < map.size(); map_i++) {
      m2_map_entry map_entry = map.get_entry(map_i);
      string marker_name = map_entry.get_marker_name();
      size_t marker_name_size = marker_name.size() + 1;

      i++;

      // The names file lists marker names by type (listed in parenthesis below).
      // Markers: autosomal numbered (M), X-linked numbered (X), Y-linked numbered (Y).
      // For a vcf file, to determine the type we must look at the CHROM (integer representation).
      //
      // Need to set 'types[]' so that LTop->Locus[].Type is initialized correctly in
      // 'read_common_marker_data()' below. SEX_CHOMOSOME, and MALE_CHROMOSOME are set correctly
      // based on the organism (e.g., human, dog...). These are also used to convert from strings
      // to numbers in 'VCFtools_get_map()' from which 'map' is derived.
      strcpy(ltype, map_entry.get_chr_type_string().c_str());
      types[i]=ltype[0];

      names[i] = CALLOC(marker_name_size, char);
      strcpy(names[i], marker_name.c_str());

      loci_name_space += marker_name_size;
      num_markers++;
    }

#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: locus strings: %d MB (%d x %d)\n",
          loci_name_space / 1024 / 1024,
          i+1,  loci_name_space / (i+1) );
#endif

    *LTop = read_common_marker_data(i+1, num_markers, names, types, annotated, penetrances_read, 0, NULL);
    // NOTE: *types, and **names are freed at the end of read_common_marker_data()...
}
