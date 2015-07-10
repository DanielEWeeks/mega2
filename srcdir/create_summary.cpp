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
#include "liable.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "input_check_ext.h"
#include "liability_dist_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "pedtree_ext.h"
#include "reorder_loci_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc my_malloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  draw_histogram kinship
        grow_string_ext.h:  grow
        input_check_ext.h:  check_renumber_ped
     liability_dist_ext.h:  liability_summary
            linkage_ext.h:  get_loci_on_chromosome get_unmapped_loci is_typed_lentry new_lpedtop
           makenucs_ext.h:  makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
            pedtree_ext.h:  convert_to_pedtree is_typed reassign_affecteds
       reorder_loci_ext.h:  display_selections
         user_input_ext.h:  pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line imax summary_time_stamp
*/


/* Separated out these routines from mega2.c */

/* prototype definitions */

void    create_summary_file(linkage_ped_top * LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr);

void marker_typing_summary(FILE *sum_fp, linkage_ped_top *LPedTreeTop,
			   int numchr);

void    seg_affected_by_status(linkage_ped_top *Top, int locus,
			       char *SegFileName,  char **input_files);
void    affected_by_status(linkage_ped_top *Top, int locus);
int aff_status_entry(int status, int eclass, linkage_locus_rec *Locus);

static void    create_allele_freq_tables(FILE * filep, linkage_ped_top * Top1,
                                         int inc_ht);
static void    create_sex_allele_freq_tables(FILE *filep,
					     linkage_ped_top *Top1,
                                             int inc_ht);
static void create_Y_allele_freq_tables(FILE *filep, linkage_ped_top *Top1,
                                        int inc_ht);
static void    write_sex_allele_freq_table(FILE *filep,
					   linkage_ped_top *Top1,
					   allele_freq_struct *ar,
					   int index);
static void    write_allele_freq_table(FILE *filep, linkage_ped_top *Top1,
				       allele_freq_struct *ar, int index);
static void  write_Y_allele_freq_table(FILE *filep, linkage_ped_top *Top1,
				       allele_freq_struct *ar,
				       int index);
static void    aff_rel_count(ped_top *Top, int *numchr,
			     char *cnt_sum_name,
			     linkage_ped_top *LPedTreeTop,
			     char **input_files, int trait_pos);
static void  aff_sib_count(linkage_ped_top *LPedTreeTop,
			   char *sib_sum_name,
			   int trp);
static void quant_phenotype_summary(linkage_ped_top *Top,
				    char *file_names[]);

/* end of prototypes */

void   create_summary_file(linkage_ped_top * LPedTreeTop,
			   analysis_type   *analysis,
			   char *file_names[],
			   int untyped_ped_opt,
			   int *numchr)
{

    int             *trp, tr, i, nloop, num_affec=num_traits;
    char            seg_sum_name[2*FILENAME_LENGTH],
        cnt_sum_name[2*FILENAME_LENGTH],
        sib_sum_name[2*FILENAME_LENGTH],
        all_sum_name[2*FILENAME_LENGTH];

    FILE           *filep;
    int             p, choice;
    int             inc_ht;
    char            fl_stat[12], **input_files, cchoice[4];
    ped_top         *PedTreeTop;
    linkage_ped_top *NewTop;

/*   input_files[0]=mega2_input_files[0]; */
/*   input_files[1]=mega2_input_files[1]; */
/*   input_files[2]=mega2_input_files[2]; */
/*   input_files[3]=mega2_input_files[3]; */

    input_files = mega2_input_files;

    if (*analysis == TO_LIABLE_FREQ) {
        if (num_traits < 1) {
            errorvf("there are NO trait loci present.\n");
            EXIT(OUTPUT_FORMAT_ERROR);
        }
        liability_summary(input_files, LPedTreeTop, numchr, untyped_ped_opt);
        return;
    }

    /* ALLELE frequencies */

    if (*analysis == TO_ALLELE_FREQ)  {
        if (main_chromocnt > 1) change_output_chr(file_names[0], 0);
        inc_ht = !HalfTypedReset;

        if (InputMode == INTERACTIVE_INPUTMODE || (! DEFAULT_OUTFILES)) {
            do	{
                draw_line();
                print_outfile_mssg();
                printf("Allele frequency summary file name menu:\n");
                printf("0) Done with this menu - please proceed.\n");
                printf(" 1) Allele frequency summary file name:    %s\t[%s]\n",
                       file_names[0],
                       ((access(file_names[0], F_OK) == 0) ? "overwrite" : "new"));
                printf("Enter 0 or 1 > ");

                fcmap(stdin, "%s", cchoice); newline;
                sscanf(cchoice, "%d", &choice);
                test_modified(choice);

                switch(choice) {
                case 0:
                    break;
                case 1 :
                    printf("Enter new allele frequency summary file name: > ");
                    fcmap(stdin, "%s", file_names[0]); newline;
                    break;
                default:
                    warn_unknown(cchoice);
                    break;
                }
            } while (choice != 0);
        }
        sprintf(all_sum_name, "%s/%s", output_paths[0], file_names[0]);
        fclose(fopen(all_sum_name, "w"));

        for (i = 0; i < main_chromocnt; i++)  {
            if (main_chromocnt > 1) {
                *numchr = global_chromo_entries[i];
            }
            omit_peds(untyped_ped_opt, LPedTreeTop);
            if (*numchr == UNKNOWN_CHROMO) {
                get_unmapped_loci(0);
            } else {
                get_loci_on_chromosome(*numchr);
            }
            filep=fopen(all_sum_name, "a");
            if (filep == NULL) {
                errorvf("Could not open %s\n", all_sum_name);
		EXIT(FILE_WRITE_ERROR);
            }

            if (i == 0) {
                summary_time_stamp(input_files, filep, "");
            }
            if (*numchr == SEX_CHROMOSOME) {
                create_sex_allele_freq_tables(filep, LPedTreeTop, inc_ht);
            } else if (*numchr == MALE_CHROMOSOME) {
                create_Y_allele_freq_tables(filep, LPedTreeTop, inc_ht);
            } else {
                create_allele_freq_tables(filep, LPedTreeTop, inc_ht);
            }
            fclose(filep);
        }
        create_mssg(*analysis);
        sprintf(err_msg, "    Summary file: %s", file_names[0]);
        mssgf(err_msg);
        draw_line();
        return;
    }

    if (*analysis == QUANT_SUMMARY) {
/*     NumChrLoci=0; */
/*     ChrLoci = NULL; */
        /* summary of quantitative phenotypes */
        /*    omit_peds(untyped_ped_opt, LPedTreeTop); */
        quant_phenotype_summary(LPedTreeTop, file_names);
    }

    /* Genotyping success rate summary */

    if (*analysis == GENOTYPING_SUMMARY) {
        if (main_chromocnt > 1)  {
            change_output_chr(file_names[0], -9);
        }
        if (!DEFAULT_OUTFILES) {
            do	{
                draw_line();
                print_outfile_mssg();
                printf("Genotyping rate summary file name selection menu:\n");
                printf("0) Done with this menu - please proceed.\n");
                if (main_chromocnt > 1) {
                    printf(" 1) Genotyping rate summary file name stem:    %s\n",
                           file_names[0]);
                } else {
                    printf(" 1) Genotyping rate summary file name:    %s\t[%s]\n",
                           file_names[0],
                           ((access(file_names[0], F_OK) == 0) ? "overwrite" : "new"));
                }

                printf("Enter 0 or 1 > ");
                fcmap(stdin, "%s", cchoice); newline;
                sscanf(cchoice, "%d", &choice);
                test_modified(choice);

                switch(choice) {
                case 0:
                    break;
                case 1:
                    printf("Enter new genotyping rate summary file name: > ");
                    fcmap(stdin, "%s", file_names[0]); newline;
                    break;

                default:
                    warn_unknown(cchoice);
                    break;
                }

            } while (choice != 0);
        }

        create_mssg(*analysis);
        draw_line();
        if (main_chromocnt > 1) {

            /* write the combined file */
            change_output_chr(file_names[0], 0);
            sprintf(all_sum_name, "%s/%s", output_paths[0], file_names[0]);
            if ((filep=fopen(all_sum_name, "w")) == NULL) {
                errorvf("Could not open %s\n", all_sum_name);
		EXIT(FILE_WRITE_ERROR);
            }
            summary_time_stamp(input_files, filep,
                               "IDs are from :Ped, :Per fields if present, or pedigree, person ids.");
            get_loci_on_chromosome(0);
            for (i=0; i < main_chromocnt; i++) {
                if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
                    get_unmapped_loci(1);
                    break;
                }
            }
            marker_typing_summary(filep, LPedTreeTop, 0);
            /*	   fclose(filep); */
            sprintf(err_msg, "Combined summary file:         %s", file_names[0]);
            mssgf(err_msg);
        }

        for (i = 0; i < main_chromocnt; i++)  {
            if (main_chromocnt > 1) {
                change_output_chr(file_names[0], -9);
            }
            omit_peds(untyped_ped_opt, LPedTreeTop);

            /* Now write chromosome specific files */
            *numchr = global_chromo_entries[i];
            if (*numchr == UNKNOWN_CHROMO) {
                get_unmapped_loci(0);
            } else {
                get_loci_on_chromosome(*numchr);
            }
            change_output_chr(file_names[0], *numchr);
            sprintf(all_sum_name, "%s/%s", output_paths[0], file_names[0]);
            if ((filep=fopen(all_sum_name, "w")) == NULL) {
                errorvf("Could not open %s\n", all_sum_name);
		EXIT(FILE_WRITE_ERROR);
            }
            summary_time_stamp(input_files, filep,
                               "IDs are from :Ped, :Per fields if present, or pedigree, person ids.");
            marker_typing_summary(filep, LPedTreeTop, *numchr);
            /*      fclose(filep); */
            sprintf(err_msg, "Chromosome %d summary file:     %s", *numchr, file_names[0]);
            mssgf(err_msg);
        }
        return;

    }

    /* SEGREGATION and allele counts */
    if (*analysis == CREATE_SUMMARY) {
        if (num_traits < 1) {
            errorf("No trait locus present.");
            EXIT(OUTPUT_FORMAT_ERROR);
        }

        /* select file names */
        choice=1;
        if (main_chromocnt > 1) {
            get_loci_on_chromosome(0);
            for (i=0; i < main_chromocnt; i++) {
                if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
                    get_unmapped_loci(1);
                    break;
                }
            }
            change_output_chr(file_names[0], 0);
            change_output_chr(file_names[1], 0);
            change_output_chr(file_names[2], 0);
        } else {
            if (global_chromo_entries[0] == UNKNOWN_CHROMO) {
                get_unmapped_loci(0);
            } else {
                get_loci_on_chromosome(global_chromo_entries[0]);
            }
        }
        /*    LPedTreeTop->OrigIds=1; */
        if (! DEFAULT_OUTFILES) {
            while(choice != 0) {
                print_outfile_mssg();
                printf("Summary file names selection menu:\n");
                draw_line();
                printf("0) Done with this menu - please proceed.\n");
                printf(" 1) Segregation summary file name:    %-15s\t%s\n",
                       file_names[0], file_status(file_names[0], fl_stat));
                printf(" 2) Relative count file name:         %-11s\t%s\n",
                       file_names[1], file_status(file_names[1], fl_stat));
                printf(" 3) Sibship  summary file name:       %-11s\t%s\n",
                       file_names[2], file_status(file_names[2], fl_stat));
                pedigree_id_item(4, CREATE_SUMMARY, OrigIds[1], 29, 2,
                                 LPedTreeTop->OrigIds);

                printf("Select options 0 - 4 > ");
                fcmap(stdin, "%s", cchoice); newline;
                sscanf(cchoice, "%d", &choice);
                test_modified(choice);

                switch (choice)	{
                case 0:
                    break;
                case 1:
                    printf("Enter new Segregation summary file name: > ");
                    fcmap(stdin, "%s", file_names[0]); newline;
                    break;
                case 2:
                    printf("Enter new Relative Count file name: > ");
                    fcmap(stdin, "%s", file_names[1]); newline;
                    break;
                case 3:
                    printf("Enter new Sibship summary file name: > ");
                    fcmap(stdin, "%s", file_names[2]); newline;
                    break;
                case 4:
                    OrigIds[1] =
                        pedigree_id_item(0, CREATE_SUMMARY, OrigIds[1], 0, 1,
                                         LPedTreeTop->OrigIds);
                    break;
                default:
                    printf("Invalid option %s\n", cchoice);
                    break;
                }
            }
            draw_line();
        }

        omit_peds(untyped_ped_opt, LPedTreeTop);
        PedTreeTop = convert_to_pedtree(LPedTreeTop, 0);
        if (PedTreeTop == NULL)    {
            errorvf("PedTreeTop is NULL.\n");
            EXIT(OUTPUT_FORMAT_ERROR);
        }

        NewTop = new_lpedtop();
        /*    sort_peds(PedTreeTop); */
        for(p=0; p < PedTreeTop->PedCnt; p++) {
            check_renumber_ped(&(PedTreeTop->PedTree[p]));
        }

        makenucs1(LPedTreeTop, PedTreeTop, NewTop);
        NLOOP;
        if (LoopOverTrait == 0) {
            trp = global_trait_entries;
            for (tr=0; tr < num_traits; tr++) {
                if (global_trait_entries[tr] < 0) continue;
                if (LPedTreeTop->LocusTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
                    trp = &(global_trait_entries[tr]);
                    break;
                }
            }
        } else {
            trp = &(global_trait_entries[0]);
        }
        for (tr=0; tr <= nloop; tr++) {
            if (nloop > 1 && tr==0) continue;
            sprintf(seg_sum_name, "%s/%s", output_paths[tr], file_names[0]);
            sprintf(cnt_sum_name, "%s/%s", output_paths[tr], file_names[1]);
            sprintf(sib_sum_name, "%s/%s", output_paths[tr], file_names[2]);
            seg_affected_by_status(LPedTreeTop, *trp,
                                   seg_sum_name, input_files);
            /* set the affected count for the current disease locus */
            affected_by_status(LPedTreeTop, *trp);
            reassign_affecteds(PedTreeTop);
            aff_rel_count(PedTreeTop, &(global_chromo_entries[0]), cnt_sum_name,
                          LPedTreeTop, input_files, *trp);
            aff_sib_count(NewTop, sib_sum_name, *trp);
            if (nloop == 1) break;
            trp++;
        }

        create_mssg(*analysis);
        sprintf(err_msg, "    Segregation summary file:      %s",
                file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "    Relative count summary file:   %s",
                file_names[1]);
        mssgf(err_msg);
        sprintf(err_msg, "    Sibship summary file:          %s",
                file_names[2]);
        mssgf(err_msg);
        draw_line();
        return;
    }
}


static void create_sex_allele_freq_tables(FILE *filep, linkage_ped_top *Top1,
                                          int inc_ht)
{

    int j, ks, is;
    int             i, k, k1, allele1, allele2;
    int             num_alleles, genotyped;
    allele_freq_struct    *ar;
    linkage_ped_rec *tpe;

    /* taken out the test for sex-chromosome,
       this is checked before invocation of this routine - Nandita */

    ar = MALLOC(allele_freq_struct);
    ar->PhenoCount=0;

    /* count phenotyped of at affection status */
    for (i=0; i<num_traits; i++) {
        SKIP_TRI(i);
        if (Top1->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION)
            for (j = 0; j < Top1->PedCnt; j++) {
                if (UntypedPeds != NULL) {
                    if (UntypedPeds[j]) {
                        continue;
                    }
                }
                for (ks = 0; ks < Top1->Ped[j].EntryCnt; ks++)
                    if (Top1->Ped[j].Entry[ks].Pheno[i].Affection.Status != 0)
                        ar->PhenoCount++;
            }
    }

    for (k1 = 0; k1 < NumChrLoci; k1++) {
        k= ChrLoci[k1];
        if (Top1->LocusTop->Locus[k].Type != NUMBERED &&
            Top1->LocusTop->Locus[k].Type != BINARY) continue;

        /* initialize allele count bins */
        ar->GenoCount = 0;
        ar->TotalAlleles = 0;
        ar->TotalKnownAlleles = 0;
        ar->TotalMaleAlleles = 0;
        ar->TotalFemaleAlleles = 0;
        ar->TotalPeople = 0;
        ar->HalfTyped=0;
        ar->FemaleCnt = 0;
        ar->MaleCnt = 0;
        ar->FemaleHetzCnt = 0;
        ar->MalesGenotyped = 0;
        ar->FemalesGenotyped = 0;
        ar->FemaleHetzObs=0;
        num_alleles=Top1->LocusTop->Locus[k].AlleleCnt;

        ar->male_allele_freq = CALLOC((size_t) num_alleles, double);
        ar->female_allele_freq = CALLOC((size_t) num_alleles, double);
        ar->allele_freq = CALLOC((size_t) num_alleles, double);
        ar->input_freq = CALLOC((size_t) num_alleles, double);
        for (i = 0; i < num_alleles; i++)
            ar->input_freq[i] =
                Top1->LocusTop->Locus[k].Allele[i].Frequency;

        ar->Male_Allele_Bin = CALLOC((size_t) num_alleles, int);
        ar->Female_Allele_Bin = CALLOC((size_t) num_alleles, int);
        ar->Allele_Bin = CALLOC((size_t) num_alleles, int);

        ar->genotype_matrix = (int **) CALLOC_PTR((size_t) num_alleles, int *);
        for (i = 0; i < num_alleles; i++)
            ar->genotype_matrix[i] = CALLOC((size_t) num_alleles, int);

        /* now tally allele_bins for each locus */
        for (is = 0; is < Top1->PedCnt; is++)  {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[is]) {
                    continue;
                }
            }
            /* Possible error: cast size_t to int */
            ar->TotalPeople += (int) Top1->Ped[is].EntryCnt;
            for (j = 0; j < Top1->Ped[is].EntryCnt; j++)   {
                tpe = &(Top1->Ped[is].Entry[j]);

                if (tpe->Sex == 2) ar->FemaleCnt++;
                else ar->MaleCnt++;

                get_2alleles(tpe->Marker, k, &allele1, &allele2);
                /* Count the alleles */
                if (inc_ht) {
                    genotyped = ((allele1 != 0 || allele2 != 0)? 1: 0);
                } else {
                    genotyped = ((allele1 != 0 && allele2 != 0)? 1: 0);
                }
                if (genotyped) {
                    /* females */
                    if (tpe->Sex == 2) {
                        if (allele1) {
                            ar->Allele_Bin[allele1 - 1]++;
                            ar->Female_Allele_Bin[allele1 - 1]++;
                        }
                        if (allele2) {
                            ar->Allele_Bin[allele2 - 1]++;
                            ar->Female_Allele_Bin[allele2 - 1]++;
                        }
                        if (genotyped && (!allele1 || !allele2)) {
                            ar->HalfTyped++;
                        }
                    } else {
                        /* males */
                        if (allele1 && !allele2) {
                            ar->Allele_Bin[allele1 - 1]++;
                            ar->Male_Allele_Bin[allele1 - 1]++;
                            ar->HalfTyped++;
                        } else if (!allele1 && allele2) {
                            ar->Allele_Bin[allele2 - 1]++;
                            ar->Male_Allele_Bin[allele2 - 1]++;
                            ar->HalfTyped++;
                        } else if (allele1 && allele2) {
                            /* assume that allele1 is the same as allele2,
                               this should have passed the Mendelian checks */
                            ar->Allele_Bin[allele1 - 1]++;
                            ar->Male_Allele_Bin[allele1 - 1]++;
                        }
                    }
                }
                /* count fully typed genotypes */
                if (allele1 != 0 && allele2 != 0) {
                    ar->GenoCount++;
                    (tpe->Sex == 1)? ar->MalesGenotyped++ : ar->FemalesGenotyped++;
                    ar->genotype_matrix[allele1 - 1][allele2 - 1]++;
                    if (allele1 != allele2 && tpe->Sex == 2) {
                        ar->FemaleHetzCnt++;
                    }
                }
            }
        }

        for (i=0; i < num_alleles; i++)     {
            ar->TotalMaleAlleles += ar->Male_Allele_Bin[i];
            ar->TotalFemaleAlleles += ar->Female_Allele_Bin[i];
            ar->TotalKnownAlleles += ar->Allele_Bin[i];
        }

        if (ar->TotalKnownAlleles == 0) {
            sprintf(err_msg,
                    "Marker %s is fully untyped, skipping this marker\n",
                    Top1->LocusTop->Locus[k].Name);
            mssgf(err_msg);
            fprintf(filep, "----------------------------------------\n");
            fprintf(filep, "%s\n", err_msg);
            fprintf(filep, "----------------------------------------\n");
            continue;
        }

        for (i=0; i < num_alleles; i++) {
            ar->allele_freq[i] =
                ((double) ar->Allele_Bin[i]) / ((double) ar->TotalKnownAlleles);
            if (ar->TotalFemaleAlleles > 0) {
                ar->female_allele_freq[i] =
                    ((double) ar->Female_Allele_Bin[i]) /
                    ((double) ar->TotalFemaleAlleles);
            } else {
                ar->female_allele_freq[i] = 0;
            }
            if (ar->TotalMaleAlleles > 0) {
                ar->male_allele_freq[i] =
                    ((double) ar->Male_Allele_Bin[i]) /
                    ((double) ar->TotalMaleAlleles);
            } else {
                ar->male_allele_freq[i] = 0;
            }
            ar->FemaleHetzObs +=
                (ar->female_allele_freq[i] * ar->female_allele_freq[i]);
        }
        ar->FemaleHetzObs = 1 - ar->FemaleHetzObs;

        write_sex_allele_freq_table(filep, Top1, ar, k);
        free(ar->Allele_Bin);  ar->Allele_Bin = NULL;
        free(ar->Male_Allele_Bin); ar->Male_Allele_Bin = NULL;
        free(ar->Female_Allele_Bin); ar->Female_Allele_Bin = NULL;
        free(ar->allele_freq); ar->allele_freq=NULL;
        free(ar->input_freq); ar->input_freq=NULL;
        free(ar->male_allele_freq); ar->male_allele_freq=NULL;
        free(ar->female_allele_freq); ar->male_allele_freq=NULL;
        for (i=0; i< num_alleles; i++)
            free(ar->genotype_matrix[i]);
        free(ar->genotype_matrix); ar->genotype_matrix = NULL;
    }
    return;
}

static void create_Y_allele_freq_tables(FILE *filep, linkage_ped_top *Top1,
                                        int inc_ht)
{
    int j, ks, is;
    int             i, k, k1, allele1, allele2;
    int             num_alleles, genotyped;
    allele_freq_struct    *ar;
    linkage_ped_rec *tpe;

    /* taken out the test for sex-chromosome,
       this is checked before invocation of this routine - Nandita */

    ar = CALLOC((size_t) 1, allele_freq_struct);
    ar->PhenoCount=0;

    /* count phenotyped of at affection status */
    for (i=0; i<num_traits; i++) {
        SKIP_TRI(i);
        if (Top1->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION)
            for (j = 0; j < Top1->PedCnt; j++) {
                if (UntypedPeds != NULL) {
                    if (UntypedPeds[j]) {
                        continue;
                    }
                }
                for (ks = 0; ks < Top1->Ped[j].EntryCnt; ks++)
                    if (Top1->Ped[j].Entry[ks].Pheno[i].Affection.Status != 0 &&
                        Top1->Ped[j].Entry[ks].Sex == 1)
                        ar->PhenoCount++;
            }
    }

    for (k1 = 0; k1 < NumChrLoci; k1++) {
        k= ChrLoci[k1];
        if (Top1->LocusTop->Locus[k].Type != NUMBERED &&
            Top1->LocusTop->Locus[k].Type != BINARY) continue;

        /* initialize allele count bins */
        ar->GenoCount = 0;
        ar->TotalMaleAlleles = 0;
        ar->HalfTyped=0;
        ar->MaleCnt = 0;
        ar->MalesGenotyped = 0;
        ar->Hetz = 0.0;

        num_alleles=Top1->LocusTop->Locus[k].AlleleCnt;

        ar->male_allele_freq = CALLOC((size_t) num_alleles, double);
        ar->input_freq = CALLOC((size_t) num_alleles, double);
        for (i = 0; i < num_alleles; i++)
            ar->input_freq[i] =
                Top1->LocusTop->Locus[k].Allele[i].Frequency;

        ar->Male_Allele_Bin = CALLOC((size_t) num_alleles, int);
        ar->genotype_matrix = CALLOC((size_t) num_alleles, int *);

        for (i = 0; i < num_alleles; i++)
            ar->genotype_matrix[i] = CALLOC((size_t) num_alleles, int);

        /* now tally allele_bins for each locus */
        for (is = 0; is < Top1->PedCnt; is++)  {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[is]) {
                    continue;
                }
            }

            for (j = 0; j < Top1->Ped[is].EntryCnt; j++)   {
                if (Top1->Ped[is].Entry[j].Sex == 2) {
                    continue;
                }

                tpe = &(Top1->Ped[is].Entry[j]);
                ar->MaleCnt++;

                get_2alleles(tpe->Marker, k, &allele1, &allele2);
                /* Count the alleles */
                if (inc_ht) {
                    genotyped = ((allele1 != 0 || allele2 != 0)? 1: 0);
                } else {
                    genotyped = ((allele1 != 0 && allele2 != 0)? 1: 0);
                }
                if (genotyped) {
                    /* females */
                    if (allele1) {
                        ar->Male_Allele_Bin[allele1 - 1]++;
                    }
                    if (allele2 && allele2 != allele1) {
                        ar->Male_Allele_Bin[allele2 - 1]++;
                    }
                    if (genotyped && (!allele1 || !allele2)) {
                        ar->HalfTyped++;
                    }
                }
                /* count fully typed genotypes */
                if (allele1 != 0 && allele2 != 0) {
                    ar->MalesGenotyped++;
                    ar->genotype_matrix[allele1 - 1][allele2 - 1]++;
                    if (allele1 != allele2 && tpe->Sex == 1) {
                        ar->Hetz++;
                    }
                }
            }
        }

        for (i=0; i < num_alleles; i++)     {
            ar->TotalMaleAlleles += ar->Male_Allele_Bin[i];
        }

        if (ar->TotalMaleAlleles == 0) {
            sprintf(err_msg,
                    "Marker %s is fully untyped, skipping this marker\n",
                    Top1->LocusTop->Locus[k].Name);
            mssgf(err_msg);
            fprintf(filep, "----------------------------------------\n");
            fprintf(filep, "%s\n", err_msg);
            fprintf(filep, "----------------------------------------\n");
            continue;
        }

        for (i=0; i < num_alleles; i++) {
            ar->male_allele_freq[i] =
                ((double) ar->Male_Allele_Bin[i]) / ((double) ar->TotalMaleAlleles);
        }

        write_Y_allele_freq_table(filep, Top1, ar, k);
        free(ar->Male_Allele_Bin); ar->Male_Allele_Bin = NULL;
        free(ar->input_freq); ar->input_freq=NULL;
        free(ar->male_allele_freq); ar->male_allele_freq=NULL;

        for (i=0; i< num_alleles; i++) {
            free(ar->genotype_matrix[i]);
        }
        free(ar->genotype_matrix); ar->genotype_matrix = NULL;
    }
    return;
}

void count_alleles(linkage_ped_top *Top1, allele_freq_struct *ar,
		   int loc, int inc_ht, int comp_geno)

{
    int i;
    int j, genotyped, all;
    int allele1, allele2, entry_count;
    void *Marker;

    /* fill in the allele and genotype counts assuming that numbered
       alleles range from 1 to allelecnt */
    for (i = 0; i < Top1->PedCnt; i++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[i]) {
                continue;
            }
        }
        /* Possible error: casting size_t to int */
        entry_count =
            (Top1->pedfile_type == PREMAKEPED_PFT)?
            (int) Top1->PTop[i].num_persons : (int) Top1->Ped[i].EntryCnt;

        ar->TotalPeople += entry_count;

        for (j = 0; j < entry_count; j++)   {
            Marker = (Top1->pedfile_type == POSTMAKEPED_PFT) ?
                Top1->Ped[i].Entry[j].Marker :
                Top1->PTop[i].persons[j].marker;

            get_2alleles(Marker, loc, &allele1, &allele2);

            if (inc_ht) {
                genotyped = (allele1 != 0 || allele2 != 0)? 1: 0;
            } else {
                genotyped = (allele1 != 0 && allele2 != 0)? 1: 0;
            }
            /* count alleles */
            if (genotyped && allele1 &&
                (allele1 <= Top1->LocusTop->Locus[loc].AlleleCnt)) {
                (ar->Allele_Bin[allele1-1])++;
            }
            if (genotyped && allele2 &&
                (allele2 <= Top1->LocusTop->Locus[loc].AlleleCnt)) {
                (ar->Allele_Bin[allele2-1])++;
            }

            if (genotyped && (allele1 == 0 || allele2 == 0)) {
                (ar->HalfTyped)++;
            }

            if (comp_geno) {
                if (allele1 && allele2) {
                    (ar->GenoCount)++;
                    if (allele1 != allele2)
                        (ar->HeteroGenoCount)++;
                    (ar->genotype_matrix[allele1-1][allele2-1])++;
                    if (allele1 != allele2)
                        (ar->genotype_matrix[allele2-1][allele1-1])++;
                }
            }
        }
    }

    for (all = 0; all < Top1->LocusTop->Locus[loc].AlleleCnt; all++) {
        ar->TotalKnownAlleles +=  ar->Allele_Bin[all];
    }

    return;
}



static void  create_allele_freq_tables(FILE *filep, linkage_ped_top *Top1,
                                       int inc_ht)
{

    int  j, ks;
    int             i, k, k1;
    int             num_alleles;
    allele_freq_struct   ar;


    ar.PhenoCount=0;

    for (i = 0; i < num_traits; i++) {
        SKIP_TRI(i);
        if (Top1->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
            for (j = 0; j < Top1->PedCnt; j++) {
                if (UntypedPeds != NULL) {
                    if (UntypedPeds[j]) {
                        continue;
                    }
                }

                for (ks = 0; ks < Top1->Ped[j].EntryCnt; ks++)
                    if (Top1->Ped[j].Entry[ks].Pheno[global_trait_entries[i]].Affection.Status != 0)
                        ar.PhenoCount++;
            }
        }
    }

    for (k1 = 0; k1 < NumChrLoci; k1++) {
        k = ChrLoci[k1];
        if (Top1->LocusTop->Locus[k].Type != NUMBERED &&
            Top1->LocusTop->Locus[k].Type != BINARY) continue;

        /* initialize the ar datastructure */
        ar.HeteroGenoCount = 0;
        ar.GenoCount = 0;
        ar.TotalAlleles = 0;
        ar.TotalPeople = 0;
        ar.HalfTyped = 0;
        ar.TotalKnownAlleles = 0;
        ar.Hetz = 0.0;
        num_alleles=Top1->LocusTop->Locus[k].AlleleCnt;
        ar.Allele_Bin = CALLOC((size_t) num_alleles, int);
        ar.allele_freq = CALLOC((size_t) num_alleles, double);
        ar.input_freq = CALLOC((size_t) num_alleles, double);
        ar.genotype_matrix = CALLOC((size_t) num_alleles, int *);
        for (i = 0; i < num_alleles; i++) {
            ar.genotype_matrix[i] = CALLOC((size_t) num_alleles, int);
            ar.input_freq[i] = Top1->LocusTop->Locus[k].Allele[i].Frequency;
/*       printf("%d %f\n", i+1, ar.input_freq[i]); */
        }
        /* fill in the allele and genotype counts assuming that numbered
           alleles range from 1 to allelecnt */
        count_alleles(Top1, &ar, k, inc_ht, 1);

        if (ar.TotalKnownAlleles == 0) {
            sprintf(err_msg,
                    "Marker %s is fully untyped, skipping this marker.",
                    Top1->LocusTop->Locus[k].Name);
            mssgf(err_msg);
            fprintf(filep, "----------------------------------------\n");
            fprintf(filep, "%s\n", err_msg);
            fprintf(filep, "----------------------------------------\n");
        } else {
            /* Calculate frequencies, and heterozygosity (1- homozygosity) */
            for (i = 0; i < num_alleles; i++) {
                ar.allele_freq[i] =
                    ((double)ar.Allele_Bin[i]) / ((double)ar.TotalKnownAlleles);
                ar.Hetz = ar.Hetz + (ar.allele_freq[i] * ar.allele_freq[i]);
            }

            ar.Hetz = 1 - ar.Hetz;

            /* write the file */
            write_allele_freq_table(filep, Top1, &ar, k);
        }
        /* free ar */
        free(ar.Allele_Bin); ar.Allele_Bin= NULL;
        free(ar.allele_freq); ar.allele_freq=NULL;
        free(ar.input_freq); ar.input_freq=NULL;
        for (i=0; i< num_alleles; i++)
            free(ar.genotype_matrix[i]);
        free(ar.genotype_matrix); ar.genotype_matrix=NULL;
    }

    return;
}

static void  write_sex_allele_freq_table(FILE *filep, linkage_ped_top *Top1,
					 allele_freq_struct *ar,
					 int index)
{

    int             i, j, k, p,  width, allele_width, num_alleles;
    int             max_width;
    double           expected_geno;
    char            blanks[10], format_str[6], float_format[7];

    strcpy(blanks, "         "); /* 9 blanks */

    num_alleles=Top1->LocusTop->Locus[index].AlleleCnt;

    width= -1;
    for (i = 0; i < Top1->LocusTop->Locus[index].AlleleCnt; i++) {
        for (j = 0; j < Top1->LocusTop->Locus[index].AlleleCnt; j++) {
            if (ar->genotype_matrix[j][i] > width)
                width = ar->genotype_matrix[j][i];
        }
    }
    width = (int) floor(log10((double) width)) + 1;
    allele_width = (int) floor(log10((double) num_alleles)) + 1;

    if (width > allele_width)    max_width = width;
    else   max_width = allele_width;

    fprintf(filep, "\n---------------------------------------------------------\n");

    if (Top1->LocusTop->Locus[index].Type == AFFECTION)    {
        fprintf(filep, "Disease locus '%s'\n", Top1->LocusTop->Locus[index].Name);
        fprintf(filep, "There are %d pedigrees containing %d individuals,\nof whom %d are phenotyped.\n",
                NumTypedPeds, ar->TotalPeople, ar->PhenoCount);
        return;
    }

    fprintf(filep, "\nChromosome %d : locus %s with %d alleles\n",
            Top1->LocusTop->Marker[index].chromosome,
            Top1->LocusTop->Locus[index].Name,
            Top1->LocusTop->Locus[index].AlleleCnt);

    fprintf(filep, "There are %d pedigrees containing %d individuals, of whom \n",
            NumTypedPeds, ar->TotalPeople);
    if (ar->HalfTyped) {
        fprintf(filep, " %d are half-typed and ", ar->HalfTyped);
    }
    fprintf(filep, " %d are fully typed.\n",
            ar->GenoCount);
    fprintf(filep, "Breakdown of genotypes:\n");
    fprintf(filep, "Male         %6d  out of     %6d\n",
            ar->MalesGenotyped, ar->MaleCnt);
    fprintf(filep, "Female       %6d  out of     %6d\n",
            ar->FemalesGenotyped, ar->FemaleCnt);
    fprintf(filep, "Total        %6d  out of     %6d\n\n",
            ar->GenoCount, ar->TotalPeople);

    fprintf(filep, "Heterozygous Females   Frequency  Count\n");
    fprintf(filep, "             Observed     %5.4f %7d\n",
            (double) (ar->FemaleHetzCnt) / (ar->FemalesGenotyped),
            (int) ar->FemaleHetzCnt);
    fprintf(filep, "             Expected     %5.4f %7.0f\n\n",
            ar->FemaleHetzObs,
            (double) (ar->FemalesGenotyped * ar->FemaleHetzObs));

    fprintf(filep, "Allele Frequency Table:\n");
    fprintf(filep, "                      Total              Females              Males \n");
    fprintf(filep, "Allele#   Input    Observed  Count   Observed  Count    Observed  Count\n");
    for (i = 0; i < num_alleles; i++) {
        fprintf(filep, "%6d    %6.4f", i+1,
                ar->input_freq[i]);
        fprintf(filep, "    %6.4f   %5d", ar->allele_freq[i], ar->Allele_Bin[i]);
        fprintf(filep, "   %6.4f    %5d",
                ar->female_allele_freq[i], ar->Female_Allele_Bin[i]);
        fprintf(filep, "     %6.4f   %5d",
                ar->male_allele_freq[i], ar->Male_Allele_Bin[i]);
        fprintf(filep, "\n");
    }

    fprintf(filep, "Total Alleles observed = %d\n\n", ar->TotalKnownAlleles);

    if (ar->MalesGenotyped > 0) {
        draw_histogram(filep, num_alleles, ar, "Male Allele Frequency Table", "male");
    } else {
        fprintf(filep, "-----------------------------------------\n");
        fprintf(filep, "    No genotyped males for this marker\n\n");
    }
    if (ar->FemalesGenotyped > 0) {
        draw_histogram(filep, num_alleles, ar,
                       "Female Allele Frequency Table",  "female");
    } else {
        fprintf(filep, "-----------------------------------------\n");
        fprintf(filep, "    No genotyped females for this marker\n\n");
    }


    /* now print out lower tridiagonal matrix of genotype distribution: */
    fprintf(filep, "Female Genotype Distribution Table (fully-typed females):\n");
    fprintf(filep, "         Observed \n");
    fprintf(filep, "         Expected \n\n");
    fprintf(filep, "allele #\n");

    /* now let's determine whether a line will be printed or not */

    for (j = 0; j < num_alleles; j++) {
        if (ar->Female_Allele_Bin[j] > 0)  {
            /* print this allele */
            sprintf(format_str, " %%%dd)", max_width);
            fprintf(filep, "\n");
            fprintf(filep, format_str, j+1);

            for (p = 0; p <= j; p++) {
                if (ar->Female_Allele_Bin[p] > 0) {
                    if (ar->genotype_matrix[j][p] == 0) {
                        for (k=0; k < max_width; k++) fprintf(filep, " ");
                        fprintf(filep, "-");
                    } else {
                        sprintf(format_str, " %%%dd", max_width);
                        fprintf(filep, format_str, ar->genotype_matrix[j][p]);
                    }
                }
            }
            fprintf(filep, "\n");
            for (k=0; k<max_width+2; k++) fprintf(filep, " ");

            for (i = 0; i <= j; i++)   {
                if (ar->Female_Allele_Bin[i] > 0) {
                    expected_geno = (double) ar->FemalesGenotyped *
                        ar->female_allele_freq[j] * ar->female_allele_freq[i];
                    if (i != j) expected_geno = 2.0 * expected_geno;
                    sprintf(float_format, " %%%d.0f", max_width);
                    fprintf(filep, float_format, expected_geno);
                }
            }
            fprintf(filep, "\n");
        }
    }
    fprintf(filep, "\n");
    /* for the allele numbers at the bottom */
    for (k=0; k < max_width+2; k++) fprintf(filep, " ");

    for (k = 0; k < num_alleles; k++) {
        if (ar->Female_Allele_Bin[k] > 0) {
            sprintf(format_str, " %%%dd", max_width);
            fprintf(filep, format_str, k + 1);
        }
    }

    fprintf(filep, " <- allele #\n");
    return;

}

static void  write_Y_allele_freq_table(FILE *filep, linkage_ped_top *Top1,
				       allele_freq_struct *ar,
				       int index)
{

    int             i, j,  width, num_alleles;
    char            blanks[10];

    strcpy(blanks, "         "); /* 9 blanks */

    num_alleles=Top1->LocusTop->Locus[index].AlleleCnt;

    width= -1;
    for (i = 0; i < Top1->LocusTop->Locus[index].AlleleCnt; i++) {
        for (j = 0; j < Top1->LocusTop->Locus[index].AlleleCnt; j++) {
            if (ar->genotype_matrix[j][i] > width)
                width = ar->genotype_matrix[j][i];
        }
    }
    width = (int) floor(log10((double) width)) + 1;

    fprintf(filep, "\n---------------------------------------------------------\n");

    if (Top1->LocusTop->Locus[index].Type == AFFECTION)    {
        fprintf(filep, "Disease locus '%s'\n", Top1->LocusTop->Locus[index].Name);
        fprintf(filep,
                "There are %d pedigrees containing %d males,\nof whom %d are phenotyped.\n",
                NumTypedPeds, ar->MaleCnt, ar->PhenoCount);
        return;
    }

    fprintf(filep, "\nChromosome %d : locus %s with %d alleles\n",
            Top1->LocusTop->Marker[index].chromosome,
            Top1->LocusTop->Locus[index].Name,
            Top1->LocusTop->Locus[index].AlleleCnt);

    fprintf(filep, "There are %d pedigrees containing %d males, of whom \n",
            NumTypedPeds, ar->MaleCnt);
    if (ar->HalfTyped) {
        fprintf(filep, " %d are half-typed and ", ar->HalfTyped);
    }
    fprintf(filep, " %d are fully typed.\n",
            ar->MalesGenotyped);

    if (ar->Hetz != 0.0) {
        fprintf(filep, "Heterozygous Males   Frequency  Count\n");
        fprintf(filep, "             Observed     %5.4f %7d\n",
                (double) (ar->Hetz) / (ar->MalesGenotyped),
                (int) ar->Hetz);
        fprintf(filep, "             Expected     0\n\n");
    }
    fprintf(filep, "Allele Frequency Table:\n");
    fprintf(filep, "Allele#   Input   Observed  Count \n");
    for (i = 0; i < num_alleles; i++) {
        fprintf(filep, "%6d    %6.4f", i+1,
                ar->input_freq[i]);
        fprintf(filep, "     %6.4f   %5d",
                ar->male_allele_freq[i], ar->Male_Allele_Bin[i]);
        fprintf(filep, "\n");
    }

    fprintf(filep, "Total Alleles observed = %d\n\n", ar->TotalMaleAlleles);

    if (ar->MalesGenotyped > 0) {
        draw_histogram(filep, num_alleles, ar,
                       "Male Allele Frequency Table", "male");
    } else {
        fprintf(filep, "-----------------------------------------\n");
        fprintf(filep, "    No genotyped males for this marker\n\n");
    }

    return;

}

/* =======================================================================*/

static void            write_allele_freq_table(FILE *filep,
                                               linkage_ped_top *Top1,
                                               allele_freq_struct *ar, int index)
{

    int             i, j, k, width, allele_width;
    int             num_alleles, max_width;
    double           expected_geno;
    char            blanks[10], format_str[10], float_format[10];

    strcpy(blanks, "         "); /* 9 blanks */

    fprintf(filep, "\n-----------------------------------------------------------\n");

    if (Top1->LocusTop->Locus[index].Type == AFFECTION)    {
        fprintf(filep, "Disease locus '%s'\n", Top1->LocusTop->Locus[index].Name);
        fprintf(filep, "There are %d pedigrees containing %d individuals,\n",
                NumTypedPeds, ar->TotalPeople);
        fprintf(filep, "of whom %d are phenotyped.\n", ar->PhenoCount);
        return;
    }

    num_alleles=Top1->LocusTop->Locus[index].AlleleCnt;
    width=0;
    for (i = 0; i < num_alleles; i++) {
        for (j = 0; j <= i; j++)
            if (ar->genotype_matrix[i][j] > width)
                width = ar->genotype_matrix[i][j];
    }
    if (width > 0) {
        width = (int) floor(log10((double) width)) + 1;
    }

    allele_width =  (int) floor(log10((double) num_alleles)) + 1;

    if (width > allele_width)    max_width = width;
    else  max_width = allele_width;

    fprintf(filep, "\nChromosome %d : locus %s with %d alleles\n",
            Top1->LocusTop->Marker[index].chromosome,
            Top1->LocusTop->Locus[index].Name,
            Top1->LocusTop->Locus[index].AlleleCnt);

    fprintf(filep, "There are %d pedigrees containing %d individuals, of whom\n",
            NumTypedPeds, ar->TotalPeople);

    if (ar->HalfTyped) {
        fprintf(filep, " %d are half-typed and \n", ar->HalfTyped);
    }
    fprintf(filep, " %d are fully typed.\n", ar->GenoCount);
    fprintf(filep, "\n");
    fprintf(filep, "Heterozygosity Frequency   Count\n");
    fprintf(filep, "     Observed  %5.4f    %7d\n",
            ((double) (ar->HeteroGenoCount)) /
            ((double) (ar->GenoCount)), (int) (ar->HeteroGenoCount));
    fprintf(filep, "     Expected  %5.4f    %7.0f\n",
            ar->Hetz, ((double) ar->GenoCount * ar->Hetz));

    fprintf(filep, "Total Alleles observed   = %5d", ar->TotalKnownAlleles);
    fprintf(filep, "\n\n");

    draw_histogram(filep, num_alleles, ar,
                   "Allele Frequency Histogram", "all");

    /* now print out lower tridiagonal matrix of genotype distribution: */
    fprintf(filep, "Genotype Distribution Table (fully-typed individuals):\n");
    fprintf(filep, "         Observed \n");
    fprintf(filep, "         Expected \n\n");
    fprintf(filep, "allele #\n");

    for (i=0; i < num_alleles; i++) {
        if (ar->Allele_Bin[i] > 0) {
            /* print <allele>) */
            sprintf(format_str, " %%%dd)", max_width);
            fprintf(filep, "\n");
            fprintf(filep, format_str, i+1);

            /* if total allele count is positive, there is at least one non-zero genotype */
            for (j=0; j<=i; j++) {
                /* print the observed genotype frequencies */
                if (ar->Allele_Bin[j] >0) {
                    if (ar->genotype_matrix[i][j] == 0)   {
                        /* print a '-' properly aligned */
                        for(k=0; k < max_width; k++) fprintf(filep, " ");
                        fprintf(filep, "-");
                    } else {
                        /* print the allele count */
                        sprintf(format_str, " %%%dd", max_width);
                        fprintf(filep, format_str, ar->genotype_matrix[i][j]);
                    }
                }
            }
            fprintf(filep, "\n");
            for(k=0; k < max_width+2; k++) fprintf(filep, " ");

            for (j=0; j<=i; j++)   {
                /* print the expected frequency */
                if (ar->Allele_Bin[j] > 0) {
                    if (j!=i)
                        expected_geno = 2.0 * ((double) ar->GenoCount)*(ar->allele_freq[j])*
                            (ar->allele_freq[i]);
                    else
                        expected_geno = ((double) ar->GenoCount)*(ar->allele_freq[i])*
                            (ar->allele_freq[i]);
                    sprintf(float_format, " %%%d.0f", max_width);
                    fprintf(filep, float_format, expected_geno);
                }
            }
            fprintf(filep, "\n");
        }
    }
    fprintf(filep, "\n");
    /* for the allele numbers at the bottom */
    for(k=0; k < max_width+2; k++) fprintf(filep, " ");

    for (i = 0; i < num_alleles; i++)  {
        if (ar->Allele_Bin[i] > 0) {
            sprintf(format_str, " %%%dd", max_width);
            fprintf(filep, format_str, i + 1);
        }
    }
    fprintf(filep, " <- allele #\n");

    return;
}

//
// status and class are from linkage.h:linkage_pedrec_data is a union (Affection(int Status, int Class)
// return 0 for unknown; otherwise 1 or 2 for the allele.
int aff_status_entry(int status, int eclass, linkage_locus_rec *Locus)
{
    int i, aff;

    // NOTE:
    // The ANALYSIS method 'allow_affection_liability_class()' returning 'false' for your
    // analysis class causes 'Labels' to have a value of NULL, 'NumLabel' of 0, and
    // 'liability_multiplier' of 0 (see user_input.cpp:define_affection_labels()).
    // This implies that 'aff' will not be set to '2' since eclass is never '2'.
    NumLabels = Locus->Pheno->Props.Affection.NumLabels;
    Labels = Locus->Pheno->Props.Affection.Labels;

    if (status == 0)
        aff=0;
    else {
        aff= eclass + liability_multiplier*status;
        for (i=0; i < NumLabels; i++) {
            if (aff == Labels[i]) {
                aff=2;
                break;
            }
        }
        if (aff != 2) { aff = 1; }
    }

    return aff;
}

/* try a new version of the routine based only on the string */

/* -----------------------------------------------------------------------------*/
/*  here is original (CHAPM version) affected_by_status , modified as noted,
    will be used at this time, exclusively by the create_apmult() */

void            affected_by_status(linkage_ped_top *Top, int locus)
{

    int             ped, entry;
    int             AffEntryCnt, i;
    int             aff;

    AffEntryCnt = 0;
    for (i=0; i< num_traits; i++) {
        if (global_trait_entries[i] == locus) {
            break;
        }
    }
    if (i >= num_traits) {
        errorf("Trait locus invalid.");
        EXIT(DATA_TYPE_ERROR);
    }

    for (ped = 0; ped < Top->PedCnt; ped++) {
        for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
            aff = aff_status_entry(Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Status,
                                   Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Class,
                                   &(Top->LocusTop->Locus[locus]));
            Top->Ped[ped].Entry[entry].Orig_status = aff;
            /* Allocate space to store an integer in TmpData */
            Top->Ped[ped].Entry[entry].TmpData = MALLOC(int);
/*       if (aff == 0) { */
/* 	continue; */
/*       } */
            if (aff != 2) {
                aff=1;
                /* TmpData is a pointer to void, so to store an integer in it,
                 * we need to first cast it to point to an integer and then
                 * dereference that.  Hence, the * (int *) syntax here:
                 */
                *(int*)Top->Ped[ped].Entry[entry].TmpData = 0;
            } else {
                *(int*)Top->Ped[ped].Entry[entry].TmpData = 1;
                AffEntryCnt++;
            }
            /*      Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Status=aff; */
        }
    }
#ifdef DEBUG
    sprintf(err_msg,
            "Found %d members with affection status 2 at locus %s.\n",
            AffEntryCnt, Top->LocusTop->Locus[locus].Name);
    mssgf(err_msg);
#endif
    return;
}

/* Modified by Nandita 11/18/02 to use the defined constant
   LIABILITYMULTIPLIER.
   This wasn't really necessary, because all we are doing here
   is matching up numbers, we could have used the old code.
   However, this works just as well.
*/

/*---------------------------------------------------------------+
  | Prompt (if necessary) for the status number to use to         |
  | label affecteds. Give the user as much information as         |
  | possible. Then label them.                                    |
  |                                                               |
  | This is used instead of affected_by_class() when there is     |
  | only one class.                                               |
  +---------------------------------------------------------------*/
void            seg_affected_by_status(linkage_ped_top *Top,
				       int locus,
				       char *SegFileName,
				       char **input_files)
{

#define ENTRY  Top->Ped[ped].Entry[entry]
#define STATUS i3

    int          ped, entry, entry1, mo, fa;
    int             UnAffCnt, AffEntryCnt;
    int             i1, i2, i3;
    int             segmat[4][4][4];
    int             sum12, sum;
    FILE           *segfp;
    char           message[FILENAME_LENGTH];

    for (i1=0; i1< num_traits; i1++) {
        if (global_trait_entries[i1] == locus) {
            break;
        }
    }
    if (i1 >= num_traits) {
        errorf("Trait locus is invalid.");
        EXIT(DATA_TYPE_ERROR);
    }

    if ((segfp = fopen(SegFileName, "w")) == NULL) {
        errorvf("Unable to open file '%s' for writing\n", SegFileName);
        EXIT(FILE_WRITE_ERROR);
    }

    sprintf(message, "Affection status locus %s",
            Top->LocusTop->Locus[locus].Name);

    summary_time_stamp(input_files, segfp, message);
    for (i1 = 0; i1 < 4; i1++)
        for (i2 = 0; i2 < 4; i2++)
            for (i3 = 0; i3 < 4; i3++)
                segmat[i1][i2][i3] = 0;

    AffEntryCnt = 0;
    UnAffCnt = 0;

    /* NM- Mar 2009, this routine used a two-pass method to first set all
       the affection status, and then create the seg-matrix in the second pass.
       This involves destroying the original status information inside the
       pedigree records. So, this has been changed to a one-pass, where we compute
       the affection status and create the matrix in one go, so that affection status
       values do not have to be stored.
    */
    /* Comment out first loop */

/*   for (ped = 0; ped < Top->PedCnt; ped++) { */
/*     if (UntypedPeds != NULL) { */
/*       if (UntypedPeds[ped]) { */
/* 	continue; */
/*       } */
/*     } */

/*     for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  { */
/*       /\* Decide who is affected *\/ */

/*       if (Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Status == 0) */
/* 	aff=0; */
/*       else { */
/* 	aff=Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Class + */
/* 	  LIABILITYMULTIPLIER * Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Status; */
/* 	for (i1=0; i1<num_labels; i1++) { */
/* 	  if (aff == labels[i1]) {aff=2; break;} */
/* 	} */
/* 	if (aff != 2) aff=1; */
/*       } */
/*       Top->Ped[ped].Entry[entry].Pheno[locus].Affection.Status=aff; */
/*     } */
/*   } */

    /* set status and increment matrix in one go */
    for (ped = 0; ped < Top->PedCnt; ped++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[ped]) {
                continue;
            }
        }

        for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
            /* Segregation analysis table */
            /* always compute the entry's status, and copy to Orig_status */
            i3 = aff_status_entry(ENTRY.Pheno[locus].Affection.Status,
                                  ENTRY.Pheno[locus].Affection.Class,
                                  &(Top->LocusTop->Locus[locus]));

            /*      if (i3  != 2) { i3 = 1;} */

            /* find parents and do the rest of the table */
            if (Top->Ped[ped].Entry[entry].Mother > 0) {
                mo=fa= -1;
                for(entry1 = 0; entry1 < Top->Ped[ped].EntryCnt; entry1++) {
                    if (mo != -1 && fa != -1) break;

                    if (Top->Ped[ped].Entry[entry1].ID == ENTRY.Mother) mo=entry1;
                    if (Top->Ped[ped].Entry[entry1].ID == ENTRY.Father) fa=entry1;
                }


                i1 = aff_status_entry(Top->Ped[ped].Entry[mo].Pheno[locus].Affection.Status,
                                      Top->Ped[ped].Entry[mo].Pheno[locus].Affection.Class,
                                      &(Top->LocusTop->Locus[locus]));
                /*	if (i1 != 2) { i1 = 1;} */

                i2 = aff_status_entry(Top->Ped[ped].Entry[fa].Pheno[locus].Affection.Status,
                                      Top->Ped[ped].Entry[fa].Pheno[locus].Affection.Class,
                                      &(Top->LocusTop->Locus[locus]));

                /*	if (i2 != 2) { i2 = 1;} */

                /* redefined macro to do the opposite */
                /*	i3 = STATUS; */
                segmat[i1][i2][i3]++;
            }

            /* Store original status for transfer to Pedtree */
            Top->Ped[ped].Entry[entry].Orig_status = STATUS;

            if (STATUS == 1) {
                UnAffCnt++;
                Top->Ped[ped].Entry[entry].TmpData = (void *) 0;
            }
            if (STATUS == 2)      {
                AffEntryCnt++;
                Top->Ped[ped].Entry[entry].TmpData = (void *) 1;
            }/* affected are counted here */
        }
    }
    sprintf(err_msg, "Summary for affection status locus %s:",
            Top->LocusTop->Locus[locus].Name);
    mssgf(err_msg);
    sprintf(err_msg, "%d affected members with status 2.\n", AffEntryCnt);
    mssgf(err_msg);
    sprintf(err_msg, "%d unaffected members with status 1.\n", UnAffCnt);
    mssgf(err_msg);
    fprintf(segfp, "%d affected members with status 2.\n", AffEntryCnt);
    fprintf(segfp, "%d unaffected members with status 1.\n\n", UnAffCnt);
    fprintf(segfp, " Segregation Mating Table Counts:\nMo Fa              0           1           2 Seg Ratio Aff/(Unaff + Aff)\n");
    fprintf(segfp, "------           ---         ---         ---  -------- --------\n");
    for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++)  {
            sum = 0;
            sum12 = 0;
            fprintf(segfp, "%d x %d -> ", i1, i2);
            for (i3 = 0; i3 < 3; i3++)   {
                fprintf(segfp, " %10d ", segmat[i1][i2][i3]);
                sum += segmat[i1][i2][i3];
                if (i3 > 0)
                    sum12 += segmat[i1][i2][i3];
            }
            if (sum > 0)
                fprintf(segfp, " %8.3f", 100.0 * (float) segmat[i1][i2][2] / ((float) sum));
            else
                fprintf(segfp, "        ");
            if (sum12 > 0)
                fprintf(segfp, " %8.3f", 100.0 * (float) segmat[i1][i2][2] / ((float) sum12));
            else
                fprintf(segfp, "        ");
            fprintf(segfp, "\n");
        }
    }
    fclose(segfp);
}


static void            aff_rel_count(ped_top *Top, int *numchr,
				     char *cnt_sum_name,
				     linkage_ped_top *LPedTreeTop,
				     char **input_files, int trait_pos)


{   /* user may now specify cnt_sum_file name */
    int is, js;
    int             ped, entry, child, i, j, l, HasHalfSibs = 0;
    int             tind = 0, taff = 0;
    int ttyped = 0, ttyped_all = 0, taff_typed = 0,taff_typed_all = 0;
    int fwid, tpeds_2typed = 0;
    int             tpairs = 0, typed = 0, allele1, allele2;
    int             tpairs_typed = 0; /* Number of genotyped affected sib pairs */
    int             thalfsibpairs = 0;	/* Number of affected half-sib pairs */
    int             thalfsibs_typed = 0;	/* Number of genotyped affected half-sib pairs */
    int             tparchildpairs = 0;	/* Number of affected parent-child
					 * pairs */
    int             tparchild_typed = 0;	/* Number of genotyped affected parent-child
                                                 * pairs */
    int             tavuncularpairs = 0;	/* Number of affected avuncular pairs */
    int             tavunc_typed = 0;	/* Number of genotyped affected avuncular pairs */
    int             trel[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    boolean         trel_print[9];
    double           k;
    double          int_part, frac_part;
    ped_tree       *PedTree;
    ped_rec        *PedEntry;
    FILE           *cntfp;
    char           format[10];

    field_widths(LPedTreeTop, NULL, &fwid, NULL, NULL, NULL);
    fwid = imax(8, fwid);
    create_formats(fwid, 0, &(format[0]), NULL);


    if ((cntfp = fopen(cnt_sum_name, "w")) == NULL) {
        errorvf("Unable to open file '%s' for writing\n", cntfp);
        EXIT(FILE_WRITE_ERROR);
    }

    summary_time_stamp(input_files, cntfp, "");
    for (ped = 0; ped < Top->PedCnt; ped++)    {
        /* initialize before we skip pedigrees as per omit list */
        PedTree = &(Top->PedTree[ped]);
        for (i = 0; i < 33; i++)
            PedTree->rel[i] = 0;
        PedTree->Naff = 0;	/* Number of affecteds in pedigree */
        PedTree->Ntyped = 0;	/* Number of people typed at at least one locus */
        PedTree->Ntyped_all = 0;	/* Number of people typed at ALL loci */
        PedTree->Naff_typed = 0;	/* Number of affecteds typed at at least one
                                         * locus */
        PedTree->Naff_typed_all = 0;	/* Number of affecteds typed at ALL loci */
        PedTree->Nsibships = 0;	/* Number of sibships with at least one
				 * affected */
        PedTree->Nsibpairs = 0;	/* Number of affected sib pairs */
        PedTree->Nsibs_typed = 0;	/* Number of genotyped affected sib pairs */
        PedTree->Nhalfsibpairs = 0;	/* Number of affected half-sib pairs */
        PedTree->Nhalfsibs_typed = 0;	/* Number of genotyped affected half-sib pairs */
        PedTree->Nparchildpairs = 0;	/* Number of affected parent-child pairs */
        PedTree->Nparchild_typed = 0;	/* Number of genotyped affected parent-child pairs */
        PedTree->Navuncularpairs = 0;	/* Number of affected avuncular pairs */
        PedTree->Navunc_typed = 0;	/* Number of genotyped affected avuncular pairs */
        PedTree->min_off_aff = 0;	/* Min # of affecteds in a sibship */
        PedTree->max_off_aff = 0;	/* Max # of affecteds in a sibship */
        /* end of initialization */

        if (UntypedPeds != NULL) {
            if (UntypedPeds[ped]) {
                continue;
            }
        }
        /* do the rest for included pedigrees */
        PedTree->min_off_aff = 99;	/* Min # of affecteds in a sibship */
        PedTree->max_off_aff = -99;	/* Max # of affecteds in a sibship */
        for (entry = 0; entry < PedTree->EntryCnt; entry++)	{
            PedEntry = &(PedTree->Entry[entry]);
            PedEntry->Mother_aff = 0;
            PedEntry->Father_aff = 0;
            PedEntry->Off_aff = 0;
            PedEntry->Dau_nrm = 0;
            PedEntry->Son_nrm = 0;
            PedEntry->Dau_aff = 0;
            PedEntry->Son_aff = 0;
            PedEntry->N_son = 0;
            PedEntry->N_dau = 0;
        }
        if (PedTree->AffectedCnt > 1) {
            for (i = 0; i < PedTree->AffectedCnt; i++) {
                for (j = i + 1; j < PedTree->AffectedCnt; j++)  {
                    /* For each pair of affecteds, figure out their kinship.  If it is k/32,
                       then increment the k-th entry in the rel array. */
                    /* If it isn't k/32 where k is an integer, then need to print it out */
                    /* Need to distinguish parent-child pairs from sib-sib pairs - I think it would
                       suffice to simply subtract the number of sibpairs from the number with
                       kinship 1/4. */

                    k = 32 * kinship(PedTree, PedTree->Affected[i], PedTree->Affected[j]);
                    frac_part = modf(k, &int_part);

                    if ((int_part == 2 || int_part == 4 || int_part == 6 || int_part == 8) &&
                        frac_part == 0.0 ) {
                        /*  	   printf("%d, %d : %d\n", PedTree->Affected[i]->ID, */
                        /*  		  PedTree->Affected[j]->ID, (int)k); */

                        PedTree->rel[(int) int_part] += 1;
                        if ((int) int_part == 4)  {
                            /*
                             * If the mother of one is the grandmother of the other AND the father
                             * of one is the grandfather of the other, then they must be an
                             * avuncular pair Assumption: If one parent is present, both must be
                             */
                            if ((PedTree->Affected[i]->Mother != NULL) &&
                                (PedTree->Affected[j]->Mother != NULL)) {
                                if (PedTree->Affected[i]->Mother->Mother != NULL) {
                                    if ((PedTree->Affected[i]->Mother->Mother->ID ==
                                         PedTree->Affected[j]->Mother->ID) &&
                                        (PedTree->Affected[i]->Mother->Father->ID ==
                                         PedTree->Affected[j]->Father->ID)) {
                                        PedTree->Navuncularpairs += 1;
                                        if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                            is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                            /* Count affected avuncular pairs where both members are genotyped */
                                            PedTree->Navunc_typed += 1;
                                        }
                                    }
                                }
                            }
                            if ((PedTree->Affected[i]->Father != NULL) &&
                                (PedTree->Affected[j]->Mother != NULL)) {
                                if (PedTree->Affected[i]->Father->Mother != NULL) {
                                    if ((PedTree->Affected[i]->Father->Mother->ID ==
                                         PedTree->Affected[j]->Mother->ID) &&
                                        (PedTree->Affected[i]->Father->Father->ID ==
                                         PedTree->Affected[j]->Father->ID)) {
                                        PedTree->Navuncularpairs += 1;
                                        if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                            is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                            /* Count affected avuncular pairs where both members are genotyped */
                                            PedTree->Navunc_typed += 1;
                                        }
                                    }
                                }
                            }
                            if ((PedTree->Affected[j]->Mother != NULL) &&
                                (PedTree->Affected[i]->Mother != NULL)) {
                                if (PedTree->Affected[j]->Mother->Mother != NULL) {
                                    if ((PedTree->Affected[j]->Mother->Mother->ID ==
                                         PedTree->Affected[i]->Mother->ID) &&
                                        (PedTree->Affected[j]->Mother->Father->ID ==
                                         PedTree->Affected[i]->Father->ID)) {
                                        PedTree->Navuncularpairs += 1;
                                        if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                            is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                            /* Count affected avuncular pairs where both members are genotyped */
                                            PedTree->Navunc_typed += 1;
                                        }
                                    }
                                }
                            }
                            if ((PedTree->Affected[j]->Father != NULL) &&
                                (PedTree->Affected[i]->Mother != NULL)) {
                                if (PedTree->Affected[j]->Father->Mother != NULL) {
                                    if ((PedTree->Affected[j]->Father->Mother->ID ==
                                         PedTree->Affected[i]->Mother->ID) &&
                                        (PedTree->Affected[j]->Father->Father->ID ==
                                         PedTree->Affected[i]->Father->ID)) {
                                        PedTree->Navuncularpairs += 1;
                                        if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                            is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                            /* Count affected avuncular pairs where both members are genotyped */
                                            PedTree->Navunc_typed += 1;
                                        }
                                    }
                                }
                            }

                            if ((PedTree->Affected[j]->Father != NULL) &&
                                (PedTree->Affected[i]->Mother != NULL) &&
                                (PedTree->Affected[j]->Father != NULL) &&
                                (PedTree->Affected[j]->Mother != NULL)) {
                                if (((PedTree->Affected[i]->Mother->ID != PedTree->Affected[j]->Mother->ID)
                                     && (PedTree->Affected[i]->Father->ID == PedTree->Affected[j]->Father->ID))
                                    || ((PedTree->Affected[i]->Mother->ID == PedTree->Affected[j]->Mother->ID)
                                        && (PedTree->Affected[i]->Father->ID !=
                                            PedTree->Affected[j]->Father->ID)))  {
                                    PedTree->Nhalfsibpairs += 1;
                                    if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                        is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                        /* Count affected halfsib pairs where both members are genotyped */
                                        PedTree->Nhalfsibs_typed += 1;
                                    }
                                }
                            }
                        }
                        if ((int) int_part == 8) {
                            if (PedTree->Affected[i]->Mother != NULL) {
                                if ((PedTree->Affected[i]->Mother->ID == PedTree->Affected[j]->ID) ||
                                    (PedTree->Affected[i]->Father->ID == PedTree->Affected[j]->ID)) {
                                    PedTree->Nparchildpairs += 1;
                                    if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                        is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                        /* Count affected parent-child pairs where both members are genotyped */
                                        PedTree->Nparchild_typed += 1;
                                    }
                                }
                            }
                            if (PedTree->Affected[j]->Mother != NULL) {
                                if ((PedTree->Affected[j]->Mother->ID == PedTree->Affected[i]->ID) ||
                                    (PedTree->Affected[j]->Father->ID == PedTree->Affected[i]->ID)) {
                                    PedTree->Nparchildpairs += 1;
                                    if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                        is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                        /* Count affected parent-child pairs where both members are genotyped */
                                        PedTree->Nparchild_typed += 1;
                                    }
                                }
                            }
                            if ((PedTree->Affected[i]->Mother != NULL) &&
                                (PedTree->Affected[j]->Mother != NULL) &&
                                (PedTree->Affected[i]->Father != NULL) &&
                                (PedTree->Affected[j]->Father != NULL)) {
                                /*
                                 * Full sibs only if both parents are the same and kinship = 8/32 =
                                 * 1/4
                                 */
                                if ((PedTree->Affected[i]->Mother->ID == PedTree->Affected[j]->Mother->ID) &&
                                    (PedTree->Affected[i]->Father->ID == PedTree->Affected[j]->Father->ID))  {
                                    /* Only count full sibs as sibpairs */
                                    PedTree->Nsibpairs += 1;
                                    if (is_typed(PedTree->Affected[i], Top->LocusTop, TYPED_AT_ANY_LOCUS) &&
                                        is_typed(PedTree->Affected[j], Top->LocusTop, TYPED_AT_ANY_LOCUS)) {
                                        /* Count affected sib pairs where both members are genotyped */
                                        PedTree->Nsibs_typed += 1;
                                    }
                                }
                            }
                        }
                    }
/* 	  if (frac_part != 0.0) { */
/* 	    printf("DEBUG: %8s %d %d ", */
/* 		   PedTree->Name,PedTree->Affected[i]->ID,PedTree->Affected[j]->ID); */
/* 	    printf(" %7.5f %7.5f %7.5f %s\n",  */
/* 		   k, frac_part, int_part,(frac_part == 0.0 ? "" : "Alert")); */
/* 	  } */
                }   /* Loop on affecteds */
            }	/* Loop on affecteds */
        }
        for (entry = 0; entry < PedTree->EntryCnt; entry++)  {
            PedEntry = &(PedTree->Entry[entry]);
            if (is_typed(PedEntry, Top->LocusTop, TYPED_AT_ANY_LOCUS))
                PedTree->Ntyped += 1;

            if (is_typed(PedEntry, Top->LocusTop, TYPED_AT_ALL_LOCI))
                PedTree->Ntyped_all += 1;
            if (PedEntry->Status == 2) {
                PedTree->Naff += 1;
                if (is_typed(PedEntry, Top->LocusTop, TYPED_AT_ANY_LOCUS))
                    PedTree->Naff_typed += 1;
                if (is_typed(PedEntry, Top->LocusTop, TYPED_AT_ALL_LOCI))
                    PedTree->Naff_typed_all += 1;
            }
            if (PedEntry->Mother != NULL)
                if (PedEntry->Mother->Status == 2)
                    PedEntry->Mother_aff += 1;
            if (PedEntry->Father != NULL)
                if (PedEntry->Father->Status == 2)
                    PedEntry->Father_aff += 1;
            for (child = 0; child < PedEntry->OffspringCnt; child++)
                if (PedEntry->Offspring[child] != NULL)  {
                    if (PedEntry->Offspring[child]->Sex == MALE_ID)
                        PedEntry->N_son += 1;
                    else
                        PedEntry->N_dau += 1;
                    if (PedEntry->Offspring[child]->Status == 2) {
                        PedEntry->Off_aff += 1;
                        if (PedEntry->Offspring[child]->Sex == MALE_ID)
                            PedEntry->Son_aff += 1;
                        else
                            PedEntry->Dau_aff += 1;
                    }
                    if (PedEntry->Offspring[child]->Orig_status == 1)  {
                        if (PedEntry->Offspring[child]->Sex == MALE_ID)
                            PedEntry->Son_nrm += 1;
                        else
                            PedEntry->Dau_nrm += 1;
                    }
                }	/* loop on child */
        } /* loop on entry */
        for (entry = 0; entry < PedTree->EntryCnt; entry++)  {
            PedEntry = &(PedTree->Entry[entry]);
            if ((PedEntry->Off_aff > 0) && (PedEntry->Off_aff > PedTree->max_off_aff))
                PedTree->max_off_aff = PedEntry->Off_aff;
            if ((PedEntry->Off_aff > 0) && (PedEntry->Off_aff < PedTree->min_off_aff))
                PedTree->min_off_aff = PedEntry->Off_aff;
            /*
             * Can count the sibships by counting all mothers who have at least one
             * affected child
             */
            /* NOTE: This doesn't work with half-sibs */

            if ((PedEntry->Off_aff > 0) && (PedEntry->Sex != MALE_ID))  {
                PedTree->Nsibships += 1;
            }
        } /* loop on entry */
    }  /* loop on ped */
    /* Count number of typed pedigrees */
    for (is = 0; is < LPedTreeTop->PedCnt; is++) {
        LPedTreeTop->Ped[is].IsTyped = 0;
        for (js = 0; js < LPedTreeTop->Ped[is].EntryCnt; js++)  {
            for (l = 0; l < NumChrLoci; l++)	   {
                if (LPedTreeTop->LocusTop->Locus[ChrLoci[l]].Type == NUMBERED) {
                    get_2alleles(LPedTreeTop->Ped[is].Entry[js].Marker, ChrLoci[l], &allele1, &allele2);
                    if (allele1 != 0 && allele2 != 0) {
                        LPedTreeTop->Ped[is].IsTyped = 1;
                        break;
                    }
                }
            }
            if (LPedTreeTop->Ped[is].IsTyped ==1) break;
        }
    }

    typed = 0;
    for (is = 0; is < LPedTreeTop->PedCnt; is++)
        if (LPedTreeTop->Ped[is].IsTyped >= 1) typed++;

    for (ped = 0; ped < Top->PedCnt; ped++)  {
        PedTree = &(Top->PedTree[ped]);
        for (i = 0; i <= 8; i++)  {
            /* The kinship 8/32 = 1/4 category includes full sibs, parent-child pairs. */
            trel[i] += PedTree->rel[i];
        }
    }
    for (i = 0; i <= 8; i++)
        if (trel[i] > 0)
            trel_print[i] = true;
        else
            trel_print[i] = false;
    fprintf(cntfp, "   Summary Counts for affection status locus %s:\n",
            LPedTreeTop->LocusTop->Locus[trait_pos].Name);
    fprintf(cntfp, "There are %d pedigrees, of which %d are typed at \n",
            Top->PedCnt, typed);
    if (main_chromocnt == 1)
        fprintf(cntfp, "one or more of %d markers on chromosome %d. \n",
                ((LoopOverTrait == 0) ?
                 NumChrLoci - (num_traits-1) : NumChrLoci - num_traits),
                *numchr);
    else {
        fprintf(cntfp, "one or more of %d markers on chromosomes",
                ((LoopOverTrait == 0) ?
                 NumChrLoci - (num_traits-1) : NumChrLoci - num_traits));

        for (i=0; i<main_chromocnt-1; i++)
            fprintf(cntfp, " %d,", global_chromo_entries[i]);
        fprintf(cntfp, " and %d.\n", global_chromo_entries[main_chromocnt-1]);
    }
    fprintf(cntfp, "             All members  |  Affecteds   |      Sibships  Aff  | Affected pairs with\n");
    fprintf(cntfp, "                 Typed    |    Typed     |         Naff   Sib  | kinship k/32\n");
    fprintf(cntfp, "Pedigree      #   >1  All |   #   >1  All|   #  Min  Max  Pairs|");
    for (i = 0; i <= 8; i++)
        if (trel_print[i])
            fprintf(cntfp, "%4d ", i);
    fprintf(cntfp, "\n");


    for (ped = 0; ped < Top->PedCnt; ped++)   {
        PedTree = &(Top->PedTree[ped]);
        if (PedTree->min_off_aff == 99) PedTree->min_off_aff = 0;
        if (PedTree->max_off_aff == -99) PedTree->max_off_aff = 0;

        prID_ped(cntfp, ped, format, &LPedTreeTop->Ped[ped], "");
        fprintf(cntfp, ":  %4d %4d %4d|%4d %4d %4d|%4d %4d %4d  %5d|",
                PedTree->EntryCnt, PedTree->Ntyped,
                PedTree->Ntyped_all, PedTree->Naff, PedTree->Naff_typed,
                PedTree->Naff_typed_all, PedTree->Nsibships,
                PedTree->min_off_aff, PedTree->max_off_aff, PedTree->Nsibpairs);
        for (i = 0; i <= 8; i++) {
            if (trel_print[i])
                fprintf(cntfp, "%4d ", PedTree->rel[i]);
        }
        fprintf(cntfp, "\n");
        tpairs += PedTree->Nsibpairs;
        tpairs_typed += PedTree->Nsibs_typed;
        thalfsibpairs += PedTree->Nhalfsibpairs;
        thalfsibs_typed += PedTree->Nhalfsibs_typed;
        tparchildpairs += PedTree->Nparchildpairs;
        tparchild_typed += PedTree->Nparchild_typed;
        tavuncularpairs += PedTree->Navuncularpairs;
        tavunc_typed += PedTree->Navunc_typed;
        tind += PedTree->EntryCnt;
        taff += PedTree->Naff;
        ttyped += PedTree->Ntyped;
        ttyped_all += PedTree->Ntyped_all;
        taff_typed += PedTree->Naff_typed;
        if (PedTree->Naff_typed > 1) {
            tpeds_2typed += 1;
        }
        taff_typed_all += PedTree->Naff_typed_all;
    }  /* loop on ped */
    fprintf(cntfp, "   TOTAL    %4d %4d %4d|%4d %4d %4d|                 %4d|",
            tind, ttyped, ttyped_all, taff, taff_typed, taff_typed_all, tpairs);
    for (i = 0; i <= 8; i++)   {
        if (trel_print[i])
            fprintf(cntfp, "%4d ", trel[i]);
    }
    fprintf(cntfp, "\n\n");
    fprintf(cntfp,"Number of families with >=2 genotyped affecteds:  %4d\n",tpeds_2typed);
    fprintf(cntfp,"Number of affecteds:                              %4d\n",taff);
    fprintf(cntfp,"Number of genotyped affecteds:                    %4d\n",taff_typed);
    fprintf(cntfp,"\n");
    fprintf(cntfp, "Type of affected pair      Count    Kinship  Genotyped\n");
    fprintf(cntfp, " sib pairs                  %4d      8/32     %4d\n", tpairs,tpairs_typed);
    fprintf(cntfp, " parent-child pairs         %4d      8/32     %4d\n", tparchildpairs,tparchild_typed);
    fprintf(cntfp, " halfsib pairs              %4d      4/32     %4d\n", thalfsibpairs,thalfsibs_typed);
    fprintf(cntfp, " avuncular pairs            %4d      4/32     %4d\n", tavuncularpairs,tavunc_typed);
    fprintf(cntfp, "\n");
    fprintf(cntfp, "NOTE: Only non-inbred affected relative pairs are counted in this table above\n");
    fprintf(cntfp, "      Kinship does not indicate the exact type of relative pair\n");
    fprintf(cntfp,
            "The kinship is 8/32 = 1/4 for parent-child and full sib pairs\n");
    fprintf(cntfp,
            "The kinship is 4/32 = 1/8 for half-sib, avuncular, and double first cousin pairs\n");
    fprintf(cntfp, "The kinship is 2/32 = 1/16 for first cousin pairs\n");
    fprintf(cntfp,
            "  The following affected relative pairs had kinships that did not fit in the\n");
    fprintf(cntfp, "  table above (person ids may have been renumbered):\n");
    for (ped = 0; ped < Top->PedCnt; ped++)   {
        PedTree = &(Top->PedTree[ped]);
        if (PedTree->AffectedCnt > 1)
            for (i = 0; i < PedTree->AffectedCnt; i++)
                for (j = i + 1; j < PedTree->AffectedCnt; j++)	   {
                    k = 32 * kinship(PedTree, PedTree->Affected[i], PedTree->Affected[j]);
                    frac_part = modf(k, &int_part);
                    if (!(int_part == 2 || int_part == 4 || int_part == 6 || int_part == 8) ||
                        frac_part != 0.0 ) {
                        prID_ped(cntfp, ped, format, &LPedTreeTop->Ped[ped], "");
                        fprintf(cntfp, ": pair %3d %3d has ",
                                PedTree->Affected[i]->ID, PedTree->Affected[j]->ID);
                        fprintf(cntfp, "kinship = %7.5f/32\n", k);
                    }
                }	/* Loop on affecteds */
    }

    /*
     * here the set of pedigrees are checked to see if one or more pedigrees
     * have half-sibs...
     */
    for (ped = 0; ped < Top->PedCnt; ped++)    {
        PedTree = &(Top->PedTree[ped]);
        for (entry = 0; entry < PedTree->EntryCnt; entry++)       {
            PedEntry = &(PedTree->Entry[entry]);
            if ((PedEntry->Off_aff > 0) && (PedEntry->Sex != MALE_ID))
                if (PedEntry->Offspring[0] != NULL)
                    if (PedEntry->Offspring[0]->Father->Off_aff != PedEntry->Off_aff)
                        HasHalfSibs = 1;
        }
    }

    /* if an affected half sib exists in a pedigree, it will be listed here... */
    if (HasHalfSibs)   {
        fprintf(cntfp, "\n* * * * * * * * * * * * * * * * * * * * * * * *");
        fprintf(cntfp, "\nNOTE: The minimum and maximum number of sibs are counted through the mother\n");
        fprintf(cntfp, "      without regard to whether or not they are half-sibs.  The following \n");
        fprintf(cntfp, "      pedigrees have affected half-sibs in them:\n");
        for (ped = 0; ped < Top->PedCnt; ped++) {
            PedTree = &(Top->PedTree[ped]);
            for (entry = 0; entry < PedTree->EntryCnt; entry++)  {
                PedEntry = &(PedTree->Entry[entry]);
                if ((PedEntry->Off_aff > 0) && (PedEntry->Sex != MALE_ID)){
                    if (PedEntry->Offspring[0] != NULL) {
                        if (PedEntry->Offspring[0]->Father->Off_aff !=
                            PedEntry->Off_aff) {
                            prID_ped(cntfp, ped, format, &LPedTreeTop->Ped[ped], "");
                            fprintf(cntfp,
                                    ": Father & mother have different numbers of affected offspring: %d %d\n",
                                    PedEntry->Offspring[0]->Father->Off_aff, PedEntry->Off_aff);
                        }
                    }
                }
            }
        }
        fprintf(cntfp, "\n* * * * * * * * * * * * * * * * * * * * * * * *\n");
    }
    fprintf(cntfp, "\nPedigree     Detailed Summary - ALL Pedigrees\n");
    fprintf(cntfp,
            "            Num Ntyped  all  Naff typed typed_all Sibs Halfsibs ParChld Avunc\n");
    for (ped = 0; ped < Top->PedCnt; ped++)   {
        PedTree = &(Top->PedTree[ped]);
        prID_ped(cntfp, ped, format, &LPedTreeTop->Ped[ped], "");
        fprintf(cntfp, ": %5d %5d %5d %5d %5d %5d   %6d %6d %6d %6d\n",
                PedTree->EntryCnt, PedTree->Ntyped,
                PedTree->Ntyped_all, PedTree->Naff, PedTree->Naff_typed, PedTree->Naff_typed_all,
                PedTree->Nsibpairs, PedTree->Nhalfsibpairs, PedTree->Nparchildpairs,
                PedTree->Navuncularpairs);

        fprintf(cntfp, "     |            Affected Relatives      |  Normal   | Parent\n");
        fprintf(cntfp, " #Off|  Off     Dau       Son     Pa   Ma |  Dau  Son |  ID\n");

        for (entry = 0; entry < PedTree->EntryCnt; entry++)     {
            PedEntry = &(PedTree->Entry[entry]);
            if (PedEntry->Off_aff > 0)
                fprintf(cntfp,
                        "%4d | %4d %4d/%-4d %4d/%-4d %4d %4d | %4d %4d |(%s %s %s)\n",
                        PedEntry->OffspringCnt, PedEntry->Off_aff, PedEntry->Dau_aff,
                        PedEntry->N_dau, PedEntry->Son_aff, PedEntry->N_son,
                        (((PedEntry->Sex == MALE_ID) && (PedEntry->Status == 2)) ? 1 : 0),
                        (((PedEntry->Sex != MALE_ID) && (PedEntry->Status == 2)) ? 1 : 0),
                        PedEntry->Dau_nrm, PedEntry->Son_nrm,
                        PedEntry->LEntry->OrigID,
                        (PedEntry->Sex == MALE_ID ? "Father" : "Mother"),
                        (PedEntry->Status == 2 ? "Affected" : ""));
        }
    }  /* loop on ped */
    fclose(cntfp);
}   /* end of aff_rel_count */


#ifdef CMP_PAIR
#undef CMP_PAIR
#endif

#define CMP_PAIR(a1, a2, b1, b2) a1 && a2 && b1 && b2 &&        \
    ((a1 != b1) || (a1 != b2) ||                                \
     (a2 != b1) || (a2 != b2))

static void aff_sib_count(linkage_ped_top *NewTop, char *sib_sum_name,
			  int tr)

{

    FILE *sibfp;

    int ped, per;
    int any_typed, all_typed, aff, norm, unk;
    int parents_any_typed, parents_all_typed;
    int affected, all_flag, any_flag, loc;

    char stat_str[3][11]= {" Unknown ", " Normal ", " Affected "};
    struct {
        int sum_aff, sum_normal, sum_unknown, sum_sibs;
        int sum_typed_at_any_loc, sum_typed_at_all_loci;
        int trio_2p_typed_any_loc, trio_2p_typed_all_loci;
        int trio_1p_typed_any_loc, trio_1p_typed_all_loci;
        int trio_0p_typed_any_loc, trio_0p_typed_all_loci;
        int *stdt_pairs;
    } counts;

    int non_singleton=0;

    counts.sum_aff= counts.sum_normal= counts.sum_unknown= counts.sum_sibs= 0;
    counts.sum_typed_at_any_loc= counts.sum_typed_at_all_loci= 0;
    counts.trio_2p_typed_any_loc= counts.trio_2p_typed_all_loci= 0;
    counts.trio_1p_typed_any_loc= counts.trio_1p_typed_all_loci= 0;
    counts.trio_0p_typed_any_loc= counts.trio_0p_typed_all_loci= 0;

    if ((sibfp = fopen(sib_sum_name, "w")) == NULL) {
        errorvf("Unable to open summary file %s\n", sib_sum_name);
        EXIT(FILE_WRITE_ERROR);
    }

    summary_time_stamp(mega2_input_files, sibfp, "");
    /* First check if there are any nuclear pedigrees with more than one person
     */
    for(ped=0; ped < NewTop->PedCnt; ped++) {
        if (NewTop->Ped[ped].EntryCnt > 1) {
            non_singleton=1;
            break;
        }
    }
    if (!non_singleton) {
        fprintf(sibfp, "None of the pedigrees have any offspring,\n");
        fprintf(sibfp, "unable to create sibling-pair summary.\n");
        fclose(sibfp); return;
    }

    /* First the sibship section for all nuclear families */

    fprintf(sibfp, "Affection status counts for sibships for %s\n",
            NewTop->LocusTop->Locus[tr].Name);
    fprintf(sibfp,
            "--------------------------------------------------------------------------\n");
    fprintf(sibfp,
            "Ped     Nuclear  Sibship  Typed at    Typed at \n");
    fprintf(sibfp,
            "Number  Fam #    Size     any locus   all loci    Affected  Normal Unknown\n");
    fprintf(sibfp,
            "--------------------------------------------------------------------------\n");
    for(ped=0; ped < NewTop->PedCnt; ped++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[NewTop->Ped[ped].origped] == 1) {
                continue;
            }
        }
        if (NewTop->Ped[ped].EntryCnt < 3) {
            continue;
        }
        any_typed=0; all_typed=0;
        aff=0; norm=0; unk=0;
        parents_any_typed =
            is_typed_lentry(NewTop->Ped[ped].Entry[0],
                            NewTop->LocusTop, TYPED_AT_ANY_LOCUS, NULL) +
            is_typed_lentry(NewTop->Ped[ped].Entry[1],
                            NewTop->LocusTop, TYPED_AT_ANY_LOCUS, NULL);

        parents_all_typed =
            is_typed_lentry(NewTop->Ped[ped].Entry[0],
                            NewTop->LocusTop, TYPED_AT_ALL_LOCI, NULL) +
            is_typed_lentry(NewTop->Ped[ped].Entry[1],
                            NewTop->LocusTop, TYPED_AT_ALL_LOCI, NULL);

        for(per=2; per < NewTop->Ped[ped].EntryCnt; per++) {
            any_flag = is_typed_lentry(NewTop->Ped[ped].Entry[per],
                                       NewTop->LocusTop, TYPED_AT_ANY_LOCUS, NULL);
            all_flag = is_typed_lentry(NewTop->Ped[ped].Entry[per],
                                       NewTop->LocusTop, TYPED_AT_ALL_LOCI, NULL);

            all_typed += all_flag; any_typed += any_flag;

            affected = aff_status_entry(NewTop->Ped[ped].Entry[per].Pheno[tr].Affection.Status,
                                        NewTop->Ped[ped].Entry[per].Pheno[tr].Affection.Class,
                                        &(NewTop->LocusTop->Locus[tr]));

            switch(affected) {
            case 0:
                unk++; break;
            case 1:
                norm++; break;
            case 2:
                aff++;
                /* if affected, then also update the trio count */
                switch(parents_all_typed) {
                case 0:
                    counts.trio_0p_typed_all_loci += all_flag;	  break;
                case 1:
                    counts.trio_1p_typed_all_loci += all_flag;	  break;
                case 2:
                    counts.trio_2p_typed_all_loci += all_flag;      break;
                default:
                    break;
                }

                switch(parents_any_typed) {
                case 0:
                    counts.trio_0p_typed_any_loc += any_flag;      break;
                case 1:
                    counts.trio_1p_typed_any_loc += any_flag;      break;
                case 2:
                    counts.trio_2p_typed_any_loc += any_flag;      break;
                default:
                    break;
                }
                break;
            default:
                /* should not happen */
                break;
            }
        }
        counts.sum_typed_at_any_loc += any_typed;
        counts.sum_typed_at_all_loci += all_typed;
        counts.sum_aff += aff;
        counts.sum_normal += norm;
        counts.sum_unknown += unk;
        /* Possible error: casting size_t to int */
        counts.sum_sibs += (int) NewTop->Ped[ped].EntryCnt - 2,
            fprintf(sibfp,
                    "%6d  %7d  %7d  %9d   %8d    %8d %7d %7d\n",
                    NewTop->Ped[ped].OriginalID,
                    NewTop->Ped[ped].Num,
                    NewTop->Ped[ped].EntryCnt - 2,
                    any_typed, all_typed,
                    aff, norm, unk);

    }

    fprintf(sibfp,
            "--------------------------------------------------------------------------\n");
    fprintf(sibfp,
            "TOTAL   %7d  %7d  %9d   %8d    %8d %7d %7d\n",
            NewTop->PedCnt, counts.sum_sibs, counts.sum_typed_at_any_loc,
            counts.sum_typed_at_all_loci,
            counts.sum_aff, counts.sum_normal, counts.sum_unknown);

    /* Second, the parent status section for all nuclear families */
    fprintf(sibfp,
            "--------------------------------------------------------------------------\n");

    fprintf(sibfp, "\nAffection status of parents for %s\n",
            NewTop->LocusTop->Locus[tr].Name);
    fprintf(sibfp,
            "---------------------------------------------------\n");

    fprintf(sibfp,
            "Ped     Nuclear  #Parents typed at    Parent status:\n");
    fprintf(sibfp,
            "Number  Fam #    Any locus  All loci   Mother X Father\n");
    fprintf(sibfp,
            "---------------------------------------------------\n");

    for(ped=0; ped < NewTop->PedCnt; ped++) {
        if (NewTop->Ped[ped].EntryCnt < 3) {
            continue;
        }
        fprintf(sibfp,
                "%6d  %7d  %9d  %8d   %sX%s \n",
                NewTop->Ped[ped].OriginalID,
                NewTop->Ped[ped].Num,
                (is_typed_lentry(NewTop->Ped[ped].Entry[0], NewTop->LocusTop,
                                 TYPED_AT_ANY_LOCUS, NULL) +
                 is_typed_lentry(NewTop->Ped[ped].Entry[1], NewTop->LocusTop,
                                 TYPED_AT_ANY_LOCUS, NULL)),
                (is_typed_lentry(NewTop->Ped[ped].Entry[0], NewTop->LocusTop,
                                 TYPED_AT_ALL_LOCI, NULL) +
                 is_typed_lentry(NewTop->Ped[ped].Entry[1], NewTop->LocusTop,
                                 TYPED_AT_ALL_LOCI, NULL)),
                stat_str[NewTop->Ped[ped].Entry[0].Pheno[tr].Affection.Status],
                stat_str[NewTop->Ped[ped].Entry[1].Pheno[tr].Affection.Status]);
    }
    fprintf(sibfp,
            "---------------------------------------------------\n");
    fprintf(sibfp, "Trios containing affected offspring:\n");
    fprintf(sibfp,
            "                     2 typed parents  1 parent typed  0 parent typed\n");
    fprintf(sibfp,
            "                     ---------------  --------------  --------------\n");
    fprintf(sibfp,
            "At any locus:            %7d         %7d          %7d\n",
            counts.trio_2p_typed_any_loc, counts.trio_1p_typed_any_loc,
            counts.trio_0p_typed_any_loc);
    fprintf(sibfp,
            "At all loci:             %7d         %7d          %7d\n",
            counts.trio_2p_typed_all_loci, counts.trio_1p_typed_all_loci,
            counts.trio_0p_typed_all_loci);

    /* Sib-pairs suitable for S-TDT analysis,
       For each marker, write number of discordant sib-pairs
       who are also have different genotypes */
    counts.stdt_pairs = CALLOC((size_t) NumChrLoci, int);
    for(loc=0; loc < NumChrLoci; loc++) {
        if (NewTop->LocusTop->Locus[ChrLoci[loc]].Type == NUMBERED) {
            counts.stdt_pairs[loc]=0;
        }
    }

    for(ped=0; ped < NewTop->PedCnt; ped++) {
        int per2;
        int a1, a2, b1, b2;
        if (NewTop->Ped[ped].EntryCnt < 4) {
            continue;
        }
        for(per=2; per < (NewTop->Ped[ped].EntryCnt-1); per++) {
            for(per2=per+1; per2<NewTop->Ped[ped].EntryCnt; per2++) {
                if ((NewTop->Ped[ped].Entry[per].Pheno[tr].Affection.Status
                    == 0) ||
                   (NewTop->Ped[ped].Entry[per2].Pheno[tr].Affection.Status
                    == 0)) {
                    continue;
                }
                if (NewTop->Ped[ped].Entry[per].Pheno[tr].Affection.Status
                    != NewTop->Ped[ped].Entry[per2].Pheno[tr].Affection.Status) {
                    for(loc=0; loc < NumChrLoci; loc++) {
                        if (NewTop->LocusTop->Locus[ChrLoci[loc]].Type != NUMBERED) {
                            continue;
                        }
                        get_2alleles(NewTop->Ped[ped].Entry[per].Marker, ChrLoci[loc], &a1, &a2);
                        get_2alleles(NewTop->Ped[ped].Entry[per2].Marker, ChrLoci[loc], &b1, &b2);

                        if (CMP_PAIR(a1, a2, b1, b2)) {
                            counts.stdt_pairs[loc]++;
                        }
                    }
                }
            }
        }
    }
    fprintf(sibfp,
            "--------------------------------------------------------------------------\n");
    /* Now print the summary */
    fprintf(sibfp, "Discordant sib-pairs with different genotypes:\n");
    fprintf(sibfp, "   Marker         Number\n");
    for(loc=0; loc < NumChrLoci; loc++) {
        if (NewTop->LocusTop->Locus[ChrLoci[loc]].Type != NUMBERED) {
            continue;
        }
        fprintf(sibfp,"    %-15s %d\n",
                NewTop->LocusTop->Locus[ChrLoci[loc]].Name,
                counts.stdt_pairs[loc]);
    }
    fprintf(sibfp,
            "---------------------------------------------------\n");

    fclose(sibfp);
    free(counts.stdt_pairs);
    return;
}


void marker_typing_summary(FILE *sum_fp, linkage_ped_top *LPedTreeTop,
			   int numchr)

{

    /* 2 sections
       1. Number of People per marker (Total, typed-all %success typed-affecteds)
       maximum and minimum over the markers
       2. Number of markers per person (Total, typed, %success rate)
       maximum and minimum over all individuals
    */

    int ped, per, nper=0;
    int nmrk=0, loc, loc1;
    int max_per_marker, min_per_marker;
    int max_per_person, min_per_person;
    int *typed_per_marker, *typed_per_person;
    int wrote_ped;
    int untyped_mssg_disp;
    int a1, a2;

    nmrk=0;
    for(loc1=0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        if ((LPedTreeTop->LocusTop->Locus[loc].Type != AFFECTION) &&
            (LPedTreeTop->LocusTop->Locus[loc].Type != QUANT)) {
            nmrk++;
        }
    }
    for(ped=0; ped<LPedTreeTop->PedCnt; ped++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[ped]) {
                continue;
            }
        }

        nper += LPedTreeTop->Ped[ped].EntryCnt;
        for(per=0; per<LPedTreeTop->Ped[ped].EntryCnt; per++) {
            LPedTreeTop->Ped[ped].Entry[per].Ngeno = 0;
            for(loc1=0; loc1 < NumChrLoci; loc1++) {
                loc = ChrLoci[loc1];
                get_2alleles(LPedTreeTop->Ped[ped].Entry[per].Marker, loc, &a1, &a2);

                if ((LPedTreeTop->LocusTop->Locus[loc].Type != AFFECTION) &&
                    (LPedTreeTop->LocusTop->Locus[loc].Type != QUANT) &&
                    (a1 != 0) &&
                    (a2 != 0)) {
                    (LPedTreeTop->Ped[ped].Entry[per].Ngeno)++;
                }
            }
        }
    }

    typed_per_person = CALLOC((size_t) nper, int);
    typed_per_marker = CALLOC((size_t) NumChrLoci, int);

    for(loc=0; loc < NumChrLoci; loc++) {
        typed_per_marker[loc]=0;
    }

    min_per_person=9999999;
    max_per_person=0;

    for(ped=0; ped < LPedTreeTop->PedCnt; ped++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[ped]) {
                continue;
            }
        }

        for(per=0; per < LPedTreeTop->Ped[ped].EntryCnt; per++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) continue;
            }
            typed_per_person[per]=LPedTreeTop->Ped[ped].Entry[per].Ngeno;
            if (typed_per_person[per] > 0 &&  typed_per_person[per] < min_per_person) {
                min_per_person = typed_per_person[per];
            }
            if (typed_per_person[per] > max_per_person) {
                max_per_person = typed_per_person[per];
            }
            for(loc1=0; loc1 < NumChrLoci; loc1++) {
                loc = ChrLoci[loc1];
                if ((LPedTreeTop->LocusTop->Locus[loc].Type == NUMBERED ||
                    LPedTreeTop->LocusTop->Locus[loc].Type == BINARY)) {
                    /* allow half-types */
                    get_2alleles(LPedTreeTop->Ped[ped].Entry[per].Marker, loc, &a1, &a2);
                    if (a1 != 0 || a2 != 0) {
                        (typed_per_marker[loc1])++;
                    }
                }
            }
        }
    }

    /* print per person summary */
    if (numchr > 0) {
        fprintf(sum_fp, " Genotyping success rate summary for chromosome %d\n",
                numchr);
    } else {
        fprintf(sum_fp, " Genotyping success rate summary for all chromosomes\n");
    }
    fprintf(sum_fp,
            "---------------------------------------------------\n");
    fprintf(sum_fp, "Summary of per-person genotyping rate \n");
    fprintf(sum_fp, "------------------------------------- \n");
    fprintf(sum_fp, "Number of markers: %d\n", nmrk);
    fprintf(sum_fp, "Maximum number of markers typed %d\n", max_per_person);
    fprintf(sum_fp, "Minimum number of markers typed %d\n",
            ((min_per_person == 9999999)? 0 : min_per_person));
    fprintf(sum_fp,
            "---------------------------------------------------\n");
    fprintf(sum_fp,
            "Ped      Person   #Markers typed  Success Rate\n");
    fprintf(sum_fp,
            "---------------------------------------------------\n");
    for(ped=0; ped < LPedTreeTop->PedCnt; ped++) {
        for(per=0; per < LPedTreeTop->Ped[ped].EntryCnt; per++) {
            if (LPedTreeTop->Ped[ped].Entry[per].Ngeno > 0) {
                fprintf(sum_fp,
                        "%-8s  %-8s %12d   %10.2f\n",
                        LPedTreeTop->Ped[ped].Name,
                        LPedTreeTop->Ped[ped].Entry[per].OrigID,
                        LPedTreeTop->Ped[ped].Entry[per].Ngeno,
                        (double)LPedTreeTop->Ped[ped].Entry[per].Ngeno/(double)nmrk);
            }
        }

    }

    fprintf(sum_fp, "---------------------------------------------------\n\n");

    min_per_marker=9999999;
    max_per_marker=0;
    for(loc1=0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        if ((LPedTreeTop->LocusTop->Locus[loc].Type == NUMBERED) ||
           (LPedTreeTop->LocusTop->Locus[loc].Type == BINARY)) {

            if (typed_per_marker[loc1] < min_per_marker) {
                min_per_marker = typed_per_marker[loc1];
            }

            if (typed_per_marker[loc1] > max_per_marker) {
                max_per_marker = typed_per_marker[loc1];
            }
        }
    }

    fprintf(sum_fp, "Summary of per-marker genotyping rate \n");
    fprintf(sum_fp, "------------------------------------- \n");
    fprintf(sum_fp, "Total number of individuals: %d\n", nper);
    fprintf(sum_fp, "Maximum number of individuals typed %d\n",
            max_per_marker);
    fprintf(sum_fp, "Minimum number of individuals typed %d\n",
            min_per_marker);

    fprintf(sum_fp,
            "---------------------------------------------------\n");

    fprintf(sum_fp,
            "Marker      Alleles   #Individuals    Percentage  \n");
    fprintf(sum_fp,
            "                       Typed          Typed       \n");
    fprintf(sum_fp,
            "---------------------------------------------------\n");

    for(loc1=0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        if ((LPedTreeTop->LocusTop->Locus[loc].Type == NUMBERED) ||
           (LPedTreeTop->LocusTop->Locus[loc].Type == BINARY)) {
            fprintf(sum_fp,
                    "%-10s  %7d   %12d          %5.2f\n",
                    LPedTreeTop->LocusTop->Locus[loc].Name,
                    LPedTreeTop->LocusTop->Locus[loc].AlleleCnt,
                    typed_per_marker[loc1],
                    100.0 * ((double)typed_per_marker[loc1]/max_per_marker));
        }
    }

    fprintf(sum_fp,
            "---------------------------------------------------\n");

    free(typed_per_person);
    free(typed_per_marker);
    untyped_mssg_disp=0;
    for(ped=0; ped < LPedTreeTop->PedCnt; ped++) {
        wrote_ped=0;
        for(per=0; per < LPedTreeTop->Ped[ped].EntryCnt; per++) {
            if (LPedTreeTop->Ped[ped].Entry[per].Ngeno == 0) {
                if (untyped_mssg_disp == 0) {
                    fprintf(sum_fp, "The following persons were fully untyped:\n");
                    untyped_mssg_disp=1;
                }
                if (!wrote_ped) {
                    fprintf(sum_fp, "Ped %s: ", LPedTreeTop->Ped[ped].Name);
                    wrote_ped=1;
                }
                fprintf(sum_fp, "%s ", LPedTreeTop->Ped[ped].Entry[per].OrigID);
            }
        }
        if (wrote_ped) { fprintf(sum_fp, "\n"); }
    }

    fclose(sum_fp);
    return;
}

/*
                           NM - quant summary doesn't have a marker item so don't have
                           to check for marker item in global_trait_entries.
*/
static void quant_phenotype_summary(linkage_ped_top *Top, char *file_names[])

{

    int *num_phenos = CALLOC((size_t) Top->PedCnt, int);
    int *num_founders  = CALLOC((size_t) Top->PedCnt, int);
    int *founders_pheno = CALLOC((size_t) Top->PedCnt, int);
    int i1, i, num_quant=0;
    int ped, ind;
    int qindex = 0; // silly compiler
    char fname[FILENAME_LENGTH];
    FILE *fp;
    int choice = -1, combine_loci = 0;

    linkage_ped_rec Entry;

    for (i=0; i < NumChrLoci; i++) {
        if (Top->LocusTop->Locus[ChrLoci[i]].Type == QUANT) {
            num_quant++;
            if (num_quant==1) {
                qindex=ChrLoci[i];
            }
        }
    }

    if (num_quant == 0) {
        warnf("No quantitative loci selected.");
        free(num_phenos);
        return;
    }

    if (num_quant > 1) {
        change_output_chr(file_names[0], 0);
        combine_loci = 1;
    } else {
        change_output_chr(file_names[0], -9);
        grow(file_names[0], ".%s", Top->LocusTop->Locus[qindex].Name);
    }

    if (! DEFAULT_OUTFILES) {
        do	{
            draw_line();
            print_outfile_mssg();
            printf("Phenotyping summary file name selection menu:\n");
            printf("0) Done with this menu - please proceed.\n");
            if (!combine_loci && num_quant > 1) {
                printf(" 1) Phenotyping rate summary file name stem:    %s\n",
                       file_names[0]);
            } else {
                printf(" 1) Phenotyping rate summary file name:    %s\t[%s]\n",
                       file_names[0],
                       ((access(file_names[0], F_OK) == 0) ? "overwrite" : "new"));
            }

            if (num_quant > 1) {
                printf(" 2) Combine loci into one output file:    %s\n",
                       yorn[combine_loci]);
                printf("Enter 0-2, 2 to toggle > ");
            } else {
                printf("Enter 0 or 1 > ");
            }

            fcmap(stdin, "%d", &choice); newline;
            test_modified(choice);

            if (choice == 1) {
                printf("Enter new phenotyping rate summary file name %s: > ",
                       ((combine_loci==0)? "stem": ""));
                fcmap(stdin, "%s", file_names[0]); newline;
            }

            if (choice == 2) {
                combine_loci = TOGGLE(combine_loci);
                change_output_chr(file_names[0], -9);
            }

        } while (choice != 0);
    }
    if (combine_loci) {
	/* Put the header information only once, at the beginning */
        sprintf(fname, "%s/%s", output_paths[0], file_names[0]);
        fp = fopen(fname, "w");
        summary_time_stamp(mega2_input_files, fp, "");
        fclose(fp);
    }

    create_mssg(QUANT_SUMMARY);

    for (i1=0; i1< NumChrLoci; i1++) {
        i = ChrLoci[i1];
        if (Top->LocusTop->Locus[i].Type != QUANT) continue;
        if (!combine_loci && num_quant > 1) {
            sprintf(fname, "%s/%s.%s", output_paths[0],
                    file_names[0], Top->LocusTop->Locus[i].Name);
        } else {
            sprintf(fname, "%s/%s", output_paths[0], file_names[0]);
        }

        if (combine_loci) {
            fp = fopen(fname, "a");
        } else {
            fp = fopen(fname, "w");
            summary_time_stamp(mega2_input_files, fp, "");
        }
        fprintf(fp, "Phenotype Summary for %s\n",
                Top->LocusTop->Locus[i].Name);
        fprintf(fp, "---------------------------------------\n");
        fprintf(fp,
                "\t\t\tAll\t\t\t\tFounders\n");
        fprintf(fp,
                "Pedigree\tMembers\tPhenotyped\t %%Phenotyped\tTotal\tPhenotyped\t%%\n");

        for(ped=0; ped<Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            num_phenos[ped]=0;
            num_founders[ped]=0;
            founders_pheno[ped]=0;
            for(ind = 0; ind < Top->Ped[ped].EntryCnt; ind++) {
                Entry = Top->Ped[ped].Entry[ind];
                if (IS_LFOUNDER(Entry)) {
                    num_founders[ped]++;
                }
                if (fabs(Entry.Pheno[i].Quant - MissingQuant) > EPSILON) {
                    num_phenos[ped]++;
                    if (IS_LFOUNDER(Entry)) {
                        founders_pheno[ped]++;
                    }
                }
            }
            fprintf(fp, "%s\t\t%d\t\t%d\t\t%4.1f\t%d\t%d\t\t%4.1f\n",
                    Top->Ped[ped].Name, Top->Ped[ped].EntryCnt,
                    num_phenos[ped],
                    100*(double)num_phenos[ped]/(double)Top->Ped[ped].EntryCnt,
                    num_founders[ped], founders_pheno[ped],
                    100*(double)founders_pheno[ped]/(double)num_founders[ped]);
        }
        fprintf(fp, "---------------------------------------\n");
        fclose(fp);
        if (!combine_loci && num_quant > 1) {
            sprintf(err_msg, "   Summary file:     %s.%s", file_names[0],
                    Top->LocusTop->Locus[i].Name);
            mssgf(err_msg);
        }
    }

    if (combine_loci || num_quant == 1) {
        sprintf(err_msg, "     Summary file:     %s", file_names[0]);
        mssgf(err_msg);
    }

    free(num_phenos); free(num_founders); free(founders_pheno);
    return;

}

void log_quant_selections(linkage_ped_top *Top)

{

    int num_select=0;
    int *select;

    num_select = NumChrLoci;
    select = &(ChrLoci[0]);
    display_selections(Top,
                       &(select[0]),
                       num_select, NULL, NULL, 0, 0);
    return;
}
