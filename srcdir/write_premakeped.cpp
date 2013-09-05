/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2013 Robert Baron, Charles P. Kollar,
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
#include "R_output.h"

#include "R_output_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
#include "write_mfiles_ext.h"
/*
           R_output_ext.h:  append_R_commands merlin_R_setup
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  switch_map
            linkage_ext.h:  connect_loops get_loci_on_chromosome get_unmapped_loci switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg shorten_path
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line global_ou_files log_line script_time_stamp strtail
        write_files_ext.h:  write_affection_data write_binary_data write_numbered_data write_quantitative_data
      write_ghfiles_ext.h:  write_gh_locus_file
       write_mfiles_ext.h:  get_aff_status
*/

#define MERLIN_MAX_LOCUS_NAME_LEN                    8


static void init_merlin_opts(merlin_opt_type *merlin_opt);

#define MERLIN_R_SETUP  ((num_traits > 0) && (analysis == TO_MERLINONLY) && \
                         (LoopOverTrait == 1 ||                         \
			  (LoopOverTrait == 0 && num_traits == 2)))

#define NUM_MERLIN_OPTS 61
static void merlin_affection_data(FILE *filep,
				  int locusnm,
				  int num_labels,
				  int *labels,
				  linkage_locus_rec *locus,
				  linkage_ped_rec *entry)
{
    int affected=0;

    if (locus->Pheno->Props.Affection.ClassCnt == 1) {
        affected = entry->Pheno[locusnm].Affection.Status;
    } else {
        affected = aff_status_entry(entry->Pheno[locusnm].Affection.Status,
                                    entry->Pheno[locusnm].Affection.Class,
                                    locus);
    }
    fprintf(filep, "  %d", affected);
    return;
}

void save_premakeped_peds(char *outfl_name, linkage_ped_top *Top,
			  analysis_type analysis, int model_file)
{
    int tr, nloop, num_affec=num_traits;
    int ped, entry, locus, loc1, j;
    int *trp, *aff_status=NULL, *paff_status, *aff_loci;
    register linkage_ped_rec *Entry;
    linkage_locus_rec *Loc;
    linkage_locus_top *LTop;
    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];
    char pformat[6], fformat[6];
    int p, pwid, fwid;

    if (analysis == TO_MERLINONLY) {
        if (num_affection > 0 && LoopOverTrait == 0) {
            aff_status = CALLOC((size_t) num_affection, int);
            aff_loci = CALLOC((size_t) num_affection, int);
            j=0;
            for (tr=0; tr < num_traits; tr++) {
                SKIP_TRI(tr)
                    if (Top->LocusTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
                        aff_loci[j++]=global_trait_entries[tr];
                    }
            }
        } else {
            aff_loci = CALLOC((size_t) 1, int);
            aff_status = CALLOC((size_t) 1, int);
        }
    }

    NLOOP;

    field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);
    create_formats(fwid, pwid, fformat, pformat);

    LTop = Top->LocusTop;
    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(pedfl, "%s/%s", output_paths[tr], outfl_name);
        if ((filep = fopen(pedfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }

        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }

            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* write the pedigree and entry numbers */
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(filep, fformat, Top->Ped[ped].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(filep, fformat, ped+1);
                } else {
                    fprintf(filep, fformat, Top->Ped[ped].Num);
                }
                if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                    fprintf(filep, pformat, Entry->OrigID);
                    if (Entry->Father != 0) {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].OrigID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].OrigID);
                    }
                } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                    fprintf(filep, pformat, Entry->UniqueID);
                    if (Entry->Father != 0) {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].UniqueID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].UniqueID);
                    }
                } else {
                    fprintf(filep, pformat, Entry->ID);
                    if (Entry->Father != 0) {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].ID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].ID);
                    }
                }
                if (Entry->Father == 0){
                    /* create the 0 strings */
                    for(p=0; p < (pwid-1); p++)    fprintf(filep, " ");
                    fprintf(filep, "0 ");
                    for(p=0; p< (pwid-1); p++)    fprintf(filep, " ");
                    fprintf(filep, "0 ");
                }

                /* write the sex */
                if (Entry->Sex == MALE_ID) fprintf(filep, " %d", MALE_ID);
                else if (Entry->Sex == FEMALE_ID) fprintf(filep, " %d", FEMALE_ID);
                else fprintf(filep, " 0");
                /* write the affection locus for some analysis types */
                if (LoopOverTrait == 1 && trp  != NULL) {
                    switch(LTop->Locus[*trp].Type) {
                    case AFFECTION:
                        if (analysis == TO_MERLIN || analysis == TO_MERLINONLY) {
                            merlin_affection_data(filep, *trp,
                                                  NumLabels, Labels,
                                                  &(LTop->Locus[*trp]),
                                                  Entry);
                        } else {
                            write_affection_data(filep, *trp, &(LTop->Locus[*trp]),
                                                 Entry);
                        }

                        if (analysis == TO_MERLINONLY) {
                            aff_status[0] =
                                get_aff_status(&(Top->LocusTop->Locus[*trp]), analysis, Entry, *trp);
                        }
                        break;

                    case QUANT:
                        if (analysis == TO_MERLIN || analysis == TO_MERLINONLY) {
                            fprintf(filep, "  ");
                            // The documentation for Merlin says to use ‘x’ rather than using the command line
                            // argument ‘-x missing_value_code’ to specify a missing quantitative value.
                            // This being the case the batch file parameter 'Value_Missing_Quant_On_Output' is
			    // not used here.
                            if (fabs(Entry->Pheno[*trp].Quant - MissingQuant) <= EPSILON)
                                fprintf(filep, "    x     ");
                            else
                                fprintf(filep, "%10.5f", Entry->Pheno[*trp].Quant);
                        } else {
                            // For the other analysis that use this routine, the missing quant on output
                            // parameter will be used.
                            // This would include analysis: TO_LOKI, TO_PREMAKEPED.
                            write_quantitative_data(filep, *trp, &(LTop->Locus[*trp]), Entry);
                        }
                        break;

                    default:
                        break;
                    }

                }
                /* Now write the genotype data */
                /* added NM: 8/29/05, special handling for Merlin
                   affection status loci with multiple liability
                   classes. */

                if (num_affection > 0) paff_status=&(aff_status[0]);

                for (loc1 = 0; loc1 < NumChrLoci; loc1++) {
                    locus = ChrLoci[loc1];
                    Loc = &(LTop->Locus[locus]);
                    switch(Loc->Type) {
                    case AFFECTION:
                        if (LoopOverTrait == 0){
                            if (analysis == TO_MERLIN || analysis == TO_MERLINONLY) {
                                merlin_affection_data(filep, locus,
                                                      NumLabels, Labels,
                                                      Loc,  Entry);
                            } else {
                                write_affection_data(filep, locus, Loc, Entry);
                            }

                            if (analysis == TO_MERLINONLY) {
                                *paff_status = get_aff_status(Loc, analysis, Entry, locus);
                                paff_status++;
                            }
                        }
                        break;

                    case QUANT:
                        if (LoopOverTrait==0) {
                            if (analysis == TO_MERLIN || analysis == TO_MERLINONLY) {
                                fprintf(filep, "  ");
                                // The documentation for Merlin says to use ‘x’ rather than using the command line
                                // argument ‘-x missing_value_code’ to specify a missing quantitative value.
                                if (fabs(Entry->Pheno[locus].Quant - MissingQuant) <= EPSILON)
                                    fprintf(filep, "    x     ");
                                else
                                    fprintf(filep, "%10.5f", Entry->Pheno[locus].Quant);
                            } else {
                                // This would include analysis: TO_LOKI, TO_PREMAKEPED.
                                write_quantitative_data(filep, locus, Loc, Entry);
                            }
                        }
                        break;
                    case NUMBERED:
                        write_numbered_data(filep, locus, Loc, Entry);
                        break;
                    case BINARY:
                        write_binary_data(filep, locus, Loc, Entry);
                        break;
                    default:
                        fprintf(stderr, "unknown locus type (locus %d) while writing output\n",
                                locus+1);
                        break;
                    }
                }
                if (analysis == TO_PREMAKEPED && Top->UniqueIds == 1) {
                    fprintf(filep, " ID: %s", Entry->UniqueID);
                } else if (analysis == TO_MERLINONLY && model_file) {
                    if (LoopOverTrait == 1 &&
                        Top->LocusTop->Pheno[*trp].Props.Affection.ClassCnt > 1) {
                        fprintf(filep,"%8d", aff_status[0]);
                    } else {
                        for (j=0; j < num_affection; j++) {
                            if (Top->LocusTop->Pheno[aff_loci[j]].Props.Affection.ClassCnt > 1) {
                                fprintf(filep,"%8d ", aff_status[j]);
                            }
                        }
                    }
                }

                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }

    if (analysis == TO_MERLINONLY) {
        free(aff_status); free(aff_loci);
    }
}

/* function to create model files with trait penetrances
   for one or more liability classes */

static void merlin_model_file(char *model_file, linkage_locus_top *LTop)

{

    int a, l;
    FILE *filep;
    char mfl_name[2*FILENAME_LENGTH];
    linkage_locus_rec *Loc;

    if (LoopOverTrait == 0 || num_traits == 1) {
        sprintf(mfl_name, "%s/%s", output_paths[0], model_file);
        if ((filep = fopen(mfl_name, "w")) == NULL) {
            errorvf("Unable to open model file '%s'\n", mfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        /*     fprintf(filep, "AFFECTION        ALLELE_FREQ   PENETRANCES  LABEL\n"); */
    }

    for (l=0; l < num_traits; l++) {
        SKIP_TRI(l)
            if (LTop->Locus[global_trait_entries[l]].Type == AFFECTION) {
                if (LoopOverTrait == 1 && num_traits >= 2) {
                    sprintf(mfl_name, "%s/%s", output_paths[l+1], model_file);
                    if ((filep = fopen(mfl_name, "w")) == NULL) {
                        errorvf("Unable to open model file '%s'\n", mfl_name);
                        EXIT(FILE_WRITE_ERROR);
                    }
                    fprintf(filep, "AFFECTION        ALLELE_FREQ   PENETRANCES  LABEL\n");
                }
                Loc = &(LTop->Locus[global_trait_entries[l]]);
                fprintf(filep, "%15s ", Loc->Name);
                fprintf(filep, "%f ", Loc->Allele[1].Frequency);
                if (Loc->Pheno->Props.Affection.ClassCnt == 1) {
                    fprintf(filep, "%f,%f,%f  ",
                            Loc->Pheno->Props.Affection.Class[0].AutoPen[0],
                            Loc->Pheno->Props.Affection.Class[0].AutoPen[1],
                            Loc->Pheno->Props.Affection.Class[0].AutoPen[2]);
                    fprintf(filep, "%s\n", Loc->Name);
                } else {
                    fprintf(filep, "*       %s\n", Loc->Name);
                    for (a=0; a < Loc->Pheno->Props.Affection.ClassCnt; a++) {
                        fprintf(filep, "%15s_covar = %d    ",
                                Loc->Name, (liability_multiplier + a + 1));
                        /* Use 1-Pen */
                        fprintf(filep, "%f,%f,%f\n",
                                (1.0 - Loc->Pheno->Props.Affection.Class[a].AutoPen[0]),
                                (1.0 - Loc->Pheno->Props.Affection.Class[a].AutoPen[1]),
                                (1.0 - Loc->Pheno->Props.Affection.Class[a].AutoPen[2]));
                        fprintf(filep, "%15s_covar = %d    ",
                                Loc->Name, (liability_multiplier*2 + a + 1));
                        /* Use Pen */
                        fprintf(filep, "%f,%f,%f\n",
                                Loc->Pheno->Props.Affection.Class[a].AutoPen[0],
                                Loc->Pheno->Props.Affection.Class[a].AutoPen[1],
                                Loc->Pheno->Props.Affection.Class[a].AutoPen[2]);
                        /* print the otherwise clause */
                    }
                    fprintf(filep, "     OTHERWISE     0.0,0.0,0.0\n");

                }

                if (LoopOverTrait == 1 && num_traits >= 2) {
                    fclose(filep);
                }
            }
    }
    if (LoopOverTrait == 0 || num_traits == 1) {
        fclose(filep);
    }
    return;

}

static void merlin_file_names(char  *file_names[],
			      int orig_id_opt[2], analysis_type analysis,
			      int has_orig, int has_uniq, int *r_setup,
			      int num_affec, int names_file)

{

    int choice = -1, item;
    int map_it, freq_it, shell_order_it, ped_it, ind_it, r_it;
    char cchoice[10], stem[15];
    char fl_stat[12];
    int input_R_setup = *r_setup;

    if (DEFAULT_OUTFILES)  {
        mssgf("Output file names set to defaults.");
        return;
    }

    map_it= freq_it= shell_order_it= ped_it= ind_it=  r_it= -1;
/*   orig_id_opt[0] = 1; */
/*   orig_id_opt[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2; */

    strcpy(stem, (num_traits > 1 && LoopOverTrait) ? "in each subdir" : "");

    // If more than one chromosome is specified, then always use one file for each.
    if (main_chromocnt > 1) {
        strcpy(stem, "stem");
        // Here -9 is a flag to return only the <file_name> removing <extension> and <rest> if they exist...
        change_output_chr(file_names[9], -9);
        change_output_chr(file_names[10], -9);
        change_output_chr(file_names[11], -9);
        change_output_chr(file_names[12], -9);
        change_output_chr(file_names[13], -9);
    }

    while (choice != 0) {
        item=3;
        draw_line();
        print_outfile_mssg();
        printf("  Merlin file name menu\n");

/*     if ((num_traits > 2 && LoopOverTrait == 0)) { */
/*       mssgf("Npl-plots not implemented for the combined traits option"); */
/*       mssgf("Try running Mega2 with Loop-over-traits option instead."); */
/*     } */

        draw_line();
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt == 1) {
            printf(" 1) Locus file name:                     %-15s\t%s\n",
                   file_names[11], file_status(file_names[11], fl_stat));
            printf(" 2) Pedigree filename:                   %-15s\t%s\n",
                   file_names[10], file_status(file_names[10], fl_stat));
            if (*r_setup || names_file) {
                printf(" %d) Map file name:                       %-15s\t%s\n",
                       item, file_names[12], file_status(file_names[12], fl_stat));
                map_it=item;
                item++;
            }
            if (names_file) {
                freq_it = item;
                printf(" %d) Frequency file name:                 %-15s\t%s\n",
                       item, file_names[9], file_status(file_names[9], fl_stat));
                item++;
            }
            shell_order_it=item;
            if (analysis == TO_MERLIN) {
                printf(" %d) Order file name:                     %-15s\t%s\n",
                       item, file_names[13], file_status(file_names[13], fl_stat));
            } else {
                printf(" %d) C-shell script name:                 %-15s\t%s\n",
                       item, file_names[13], file_status(file_names[13], fl_stat));
            }
            item++;
        } else {
            printf(" 1) Locus file name %s:             %-15s\n", stem,
                   file_names[11]);
            printf(" 2) Pedigree file name %s:          %-15s\n", stem,
                   file_names[10]);
            if (*r_setup || names_file) {
                map_it = item;
                printf(" %d) Map file name %s:               %-15s\n", item,
                       stem, file_names[12]); item++;
            }
            if (names_file) {
                freq_it = item;
                printf(" %d) Frequency file name %s:         %-15s\t\n", item,
                       stem, file_names[9]); item++;
            }
            shell_order_it = item;
            if (analysis == TO_MERLIN) {
                printf(" %d) Order file name %s:             %-15s\n",
                       item, stem, file_names[13]); item++;
            } else {
                printf(" %d) C-shell script name %s:              %-15s\n",
                       item,  stem, file_names[13]);
                item++;
            }
        }

        ind_it = item;
        individual_id_item(item, analysis, orig_id_opt[0], 31, 2, 0, 0);
        item++;
        ped_it = item;
        pedigree_id_item(item, analysis, orig_id_opt[1], 26, 2, 0);
        if (analysis == TO_MERLINONLY) {
            if (input_R_setup) {
                item++;
                printf(" %d) Generate R-plots for Merlin:         [%s]\n",
                       item, yorn[*r_setup]); r_it = item;
                printf("Select 1-%d to enter new values, %d to toggle > ",
                       item-1, item);
            } else {
                *r_setup=0;
                printf("Select 1-%d to enter new values > ", item);
            }
        } else {
            printf("Select options 0-%d > ", item);
        }

        choice = -1;
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        case 1:
            printf("New locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[11]);
            printf("\n");
            break;
        case 2:
            printf("New pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[10]);
            printf("\n");
            break;
        default:
            if (choice == map_it) {
                printf("New map file name %s > ", stem);
                fcmap(stdin, "%s", file_names[12]);
                printf("\n");
            } else if (choice == freq_it) {
                printf("New frequency file name %s > ", stem);
                fcmap(stdin, "%s", file_names[9]);
                printf("\n");
            } else if (choice == shell_order_it) {
                if (analysis == TO_MERLIN) {
                    printf("New order file name %s > ", stem);
                } else {
                    printf("New C-shell script name %s > ", stem);
                }
                fcmap(stdin, "%s", file_names[13]);
                newline;
            } else if (choice == ind_it) {
                orig_id_opt[0] =
                    individual_id_item(0, analysis, orig_id_opt[0],
                                       35, 1, has_orig, has_uniq);
            } else if (choice == ped_it) {
                orig_id_opt[1] = pedigree_id_item(0, analysis, orig_id_opt[1],
                                                  35, 1, has_orig);
            } else if (choice == r_it) {
                if (num_traits > 0) {
                    *r_setup=TOGGLE(*r_setup);
                    break;
                }
            } else {
                warn_unknown(cchoice);
            }
            break;
        }
    }

    individual_id_item(0, analysis, orig_id_opt[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, analysis, orig_id_opt[1], 0, 3, has_orig);
    if (analysis == TO_MERLINONLY && main_chromocnt > 1) {
        change_output_chr(file_names[13], global_chromo_entries[0]);
        strcat(file_names[13], ".sh");
    }

}

static void pmk_file_names(int glob_files, char  *file_names[],
			   int *orig_id_opt, int has_orig, int has_uniq)

{
    int choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12];
    analysis_type analysis = TO_PREMAKEPED;

    ((num_traits > 1 && LoopOverTrait)?
     strcpy(stem, "in each subdir") : strcpy(stem, ""));

    if (main_chromocnt > 1 && glob_files != 1) strcpy(stem, "stem");

    if (glob_files) {
        // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, 0);
    }

    if (!glob_files && main_chromocnt > 1) {
        // returns only the <file_name> removing <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, -9);
    }

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        return;
    }

    while (choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("  Pre-makeped file name menu\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt == 1 || glob_files) {
            printf(" 1) Locus file name:           %-15s\t%s\n", file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 2) Pedigree filename:         %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
        } else {
            printf(" 1) Locus file name %s:       %-15s\n", stem, file_names[1]);
            printf(" 2) Pedigree file name %s:    %-15s\n", stem, file_names[0]);
        }

        individual_id_item(3, TO_PREMAKEPED, orig_id_opt[0], 43, 2, 0, 0);
        pedigree_id_item(4, TO_PREMAKEPED, orig_id_opt[1], 43, 2, 0);

        printf("Select options 0-4 > ");

        choice = -1;
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        case 1:
            printf("New locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]);
            printf("\n");
            break;
        case 2:
            printf("New pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);
            printf("\n");
            break;
        case 3:
            orig_id_opt[0] =
                individual_id_item(0, TO_PREMAKEPED, orig_id_opt[0], 35, 1,
                                   has_orig, has_uniq);
            break;
        case 4:
            orig_id_opt[1] = pedigree_id_item(0, TO_PREMAKEPED,
                                              orig_id_opt[1], 35, 1, has_orig);
            break;
        default:
            warn_unknown(cchoice);
            break;
        }
    }
    individual_id_item(0, TO_PREMAKEPED, orig_id_opt[0], 0, 3, has_orig,
                       has_uniq);
    pedigree_id_item(0, TO_PREMAKEPED, orig_id_opt[1], 0, 3, has_orig);
}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void  create_premakeped_files(linkage_ped_top **LPedTop, int *numchr,
			      file_format  *infl_type,
			      char *file_names[], int  untyped_ped_opt)

{
    int i, glob_files;
    int xlinked;

    linkage_ped_tree *LPed;
    analysis_type analysis=TO_PREMAKEPED;

    if (main_chromocnt > 1)
        glob_files=global_ou_files(1);
    else
        glob_files=0;
//SL
    if (Top->LocusTop->SexLinked == 2 && glob_files == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].item_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }

    pmk_file_names(glob_files, file_names, OrigIds, Top->OrigIds,
                   Top->UniqueIds);

    if (*infl_type == LINKAGE) {
        for (i = 0; i < Top->PedCnt; i++)  {
            LPed = &(Top->Ped[i]);
            if (LPed->Loops != NULL)   {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)    {
                    errorf("Fatal error! Aborting.");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    /* for creating subsets of data */
    /*   Copy=new_lpedtop(); */
    /*   copy_linkage_ped_top(Top, Copy, 1); */

    for (i=0; i < main_chromocnt; i++) {
        if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
            UnmappedSelected = 1;
            break;
        }
    }

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
            get_loci_on_chromosome(0);
            if (UnmappedSelected==1) {
                /* append the unmapped loci */
                get_unmapped_loci(1);
            }
        }
        xlinked = (((Top->LocusTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
                    (Top->LocusTop->SexLinked == 1))? 1 : 0);

        omit_peds(untyped_ped_opt, Top);

        if (i==0)   create_mssg(TO_PREMAKEPED);

        save_premakeped_peds(file_names[0], Top, TO_PREMAKEPED, 0);
        write_gh_locus_file(file_names[1], Top->LocusTop, analysis, xlinked);

        sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "        Locus file:           %s", file_names[1]);
        mssgf(err_msg);
        draw_line();
        if (glob_files) {
            break;
        }
    }

/*   free_all_including_lpedtop(&Copy); */
}

#undef Top


static void merlin_output_ped_order(linkage_ped_top *Top, char *order_fl)

{
    int ped, loc;
    int nloop, tr, num_affec=num_traits;

    FILE *fp;
    char ofl[2*FILENAME_LENGTH];

    NLOOP;

    for(tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(ofl, "%s/%s", output_paths[tr], order_fl);
        if ((fp = fopen(ofl, "w")) == NULL) {
            errorvf("could not open order file %s.\n", ofl);
            EXIT(FILE_WRITE_ERROR);
        }

        for(ped=0; ped < Top->PedCnt; ped++) {
            if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                fprintf(fp, "%s\t", Top->Ped[ped].Name);
            } else if (OrigIds[1] == 3) {
                fprintf(fp, "%d\t", ped + 1);
            } else {
                fprintf(fp, "%d\t", Top->Ped[ped].Num);
            }
        }
        fprintf(fp, "\n");

        for(loc=0; loc < NumChrLoci; loc++) {
            if (Top->LocusTop->Locus[ChrLoci[loc]].Type == NUMBERED)
                fprintf(fp, "%s\t",
                        strtail(Top->LocusTop->Locus[ChrLoci[loc]].Name, MERLIN_MAX_LOCUS_NAME_LEN));
        }
        fprintf(fp, "\n");
        fclose(fp);
        if (nloop == 1) break;
    }
}

/********************************
   Merlin option, for now separate files per chromosome, since
   Merlin simply concatenates results  for each chromosome
   into the same file otherwise
********************************/

static void write_merlin_map(linkage_locus_top *LTop,
			     char *mfl_name, int numchr,
			     analysis_type analysis)
{
    int nloop, tr, locus, locus1, num_affec=num_traits;
    int display_warning_XwAverage=0, display_warning_SexwAverage=0, display_warning_FemaleEverywhere=0;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;
    // Where to store the maps converted to Haldane...
    double *position=NULL, *female_pos=NULL, *male_pos=NULL;

    NLOOP;

    // Convert from Kosambi to Haldane if necessary for all three maps <average, famele, male>
    if (LTop->map_distance_type == 'k') {
        // Assume that the average map is always there...
        if ((position = CALLOC((size_t) NumChrLoci, double)) == NULL) {
            errorf("Unable to allocate memory.\n");
            EXIT(MEMORY_ALLOC_ERROR);
        }
        switch_map(LTop, position, SEX_AVERAGED_MAP);

        if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
	    genetic_distance_sex_type_map == FEMALE_GDMT) {

            if ((female_pos = CALLOC((size_t) NumChrLoci, double)) == NULL) {
                errorf("Unable to allocate memory.\n");
                EXIT(MEMORY_ALLOC_ERROR);
            }
            switch_map(LTop, female_pos, FEMALE_SEX_MAP);

	    if (genetic_distance_sex_type_map != FEMALE_GDMT) {
	      if ((male_pos = CALLOC((size_t) NumChrLoci, double)) == NULL) {
                errorf("Unable to allocate memory.\n");
                EXIT(MEMORY_ALLOC_ERROR);
	      }
	      switch_map(LTop, male_pos, MALE_SEX_MAP);
	    }
        }
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(mfl, "%s/%s", output_paths[tr], mfl_name);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("could not open map file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        // The Merlin map rules:
        // - If on the X, put the female map into all three map columns of the Merlin map file
        // if a sex-specific or female map was specified. However, if the user has specified a
        // sex-averaged map then warn them that it is being used for the X (in a three column map,
        // just in case they want to fudge it like it currently lets them do).
        // - If on an autosome, when sex-specific maps are chosen, all three columns of the Merlin
        // map file "POSITION FEMALE_POSITION MALE_POSITION" must be populated with their correct
        // values from the sex-averaged, female, and map maps, respectively.  So to set up for
        // sex-specific analyses on the autosome, the user must have input a map file that
        // contains valid values for all three of these maps.

        // From PedigreeGlobals::LoadMarkerMap(): three columns is !sexSpecificMap; five columns is sexSpecificMap
        fprintf(fp,  "CHR  MARKER         POSITION");
        if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
            genetic_distance_sex_type_map == FEMALE_GDMT) {
            fprintf(fp, "   FEMALE_POSITION  MALE_POSITION");
        }
        fprintf(fp, "\n");

        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            locus = ChrLoci[locus1];
            if (LTop->Locus[locus].Type == NUMBERED ||
                LTop->Locus[locus].Type == BINARY) {
                
                fprintf(fp, "%2d  %-15s",
                        LTop->Marker[locus].chromosome,
                        ((analysis == TO_MERLIN) ?
                         strtail(LTop->Locus[locus].Name, MERLIN_MAX_LOCUS_NAME_LEN) :
                         LTop->Locus[locus].Name));
                
                if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                    // Even for the X chromosome we use the sex-averaged position to give the user an opportunity
                    // to "fudge it" by copying the female map position into the sex-average map...
                    double apos = (LTop->map_distance_type == 'k') ? position[locus1] : LTop->Marker[locus].pos_avg;
                    fprintf(fp, " %10f\n", apos); // three column map
                    if (LTop->Marker[locus].chromosome == SEX_CHROMOSOME && display_warning_XwAverage == 0) {
                        warnf("Value_Genetic_Distance_SexTypeMap specified a sex-averaged map");
                        warnf("which is being used on the X-Chromosome.");
                        display_warning_XwAverage++;
                    }
                    
                } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
                    // On the X chromosome only use female map positions...
                    if (LTop->Marker[locus].chromosome == SEX_CHROMOSOME) {
                        double fpos = (LTop->map_distance_type == 'k') ? female_pos[locus1] : LTop->Marker[locus].pos_female;
                        fprintf(fp, " %10f    %10f     %10f\n", fpos, fpos, fpos); // five column map
                        if (display_warning_FemaleEverywhere == 0) {
                            warnf("You specified sex-specific analysis on the X chromosome.");
                            warnf("For Merlin to work the female map must be copied to all positions.");
                            display_warning_FemaleEverywhere++;
                        }
                    } else {
                        double apos=0, fpos=0.0, mpos=0.0;
                        if (LTop->map_distance_type == 'k') {
                            apos = position[locus1];
                            fpos = female_pos[locus1];
                            mpos = male_pos[locus1];
                        } else {
                            apos = LTop->Marker[locus].pos_avg;
                            fpos = LTop->Marker[locus].pos_female;
                            mpos = LTop->Marker[locus].pos_male;
                        }
                        fprintf(fp, " %10f    %10f     %10f\n", apos, fpos, mpos); // five column map
                        if (display_warning_SexwAverage == 0) {
                            warnf("You specified sex-specific analysis, but Merlin also requires");
                            warnf("a sex-averaged map, so it will also be used.");
                            warnf("For Merlin to work properly you must have included a sex-averaged");
                            warnf("map in your input map file.");
                            display_warning_SexwAverage++;
                        }
                    }
                    
                } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
                    if (LTop->Marker[locus].chromosome == SEX_CHROMOSOME) {
                        double fpos = (LTop->map_distance_type == 'k') ? female_pos[locus1] : LTop->Marker[locus].pos_female;
                        fprintf(fp, " %10f    %10f     %10f\n", fpos, fpos, fpos);
                    } else {
                        errorvf("A female-only genetic map can only be used on the X chromosome.\n");
                        EXIT(DATA_TYPE_ERROR);
                    }
                } else {
                    errorvf("Either a sex-averaged, sex-specific, or female only genetic map must be provided to write the map file.\n");
                    EXIT(DATA_TYPE_ERROR);
                }
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }

    if (position != NULL) free(position);
    if (female_pos != NULL) free(female_pos);
    if (male_pos != NULL) free(male_pos);
}

static void write_merlin_freq(linkage_locus_top *LTop,
			      char *freqfl, int numchr,
			      analysis_type analysis)

{
    int nloop, tr, locus, locus1, num_affec=num_traits ;
    int all;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], freqfl);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("could not open frequency file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            locus = ChrLoci[locus1];
            if (LTop->Locus[locus].Type == NUMBERED ||
                LTop->Locus[locus].Type == BINARY) {
                fprintf(fp, "M  %s\n",
                        ((analysis == TO_MERLIN)?
                         strtail(LTop->Locus[locus].Name, MERLIN_MAX_LOCUS_NAME_LEN) :
                         LTop->Locus[locus].Name));
                fprintf(fp, "F ");
                for(all=0; all < LTop->Locus[locus].AlleleCnt; all++) {
                    fprintf(fp, "%6f ",
                            LTop->Locus[locus].Allele[all].Frequency);
                }
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
}

static void write_qtdt_locus_file(char *loc_file, linkage_locus_top *LTop,
				  analysis_type analysis, int model_file)

{

    int nloop, tr, *trp, locus, num_affec=num_traits ;
    char locfl[2*FILENAME_LENGTH];
    int locus1;
    FILE *fp;

    NLOOP;

    trp = (num_traits > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(locfl, "%s/%s", output_paths[tr], loc_file);
        if ((fp = fopen(locfl, "w")) == NULL) {
            errorvf("could not open locus file %s.\n", locfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (LoopOverTrait == 1) {
            /* trait locus first */
            if (LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(fp, "A  ");
            } else {
                fprintf(fp, "T  ") ;
            }
            fprintf(fp, "%s\n",
                    ((analysis == TO_MERLIN) ?
                     strtail(LTop->Locus[*trp].Name, MERLIN_MAX_LOCUS_NAME_LEN) :
                     LTop->Locus[*trp].Name));
        }

        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            char *locname;

            locus = ChrLoci[locus1];
            locname =  ((analysis == TO_MERLIN) ?
                        strtail(LTop->Locus[locus].Name, MERLIN_MAX_LOCUS_NAME_LEN) :
                        LTop->Locus[locus].Name);

            if (LTop->Locus[locus].Class == COVARIATE) {
                fprintf(fp, "C  %s\n", locname);
            } else {
                switch(LTop->Locus[locus].Type) {
                case NUMBERED:
                case BINARY:
                    fprintf(fp, "M  %s\n", locname);
                    break;
                case AFFECTION:
                    if (LoopOverTrait == 0) {
                        fprintf(fp, "A  %s\n", locname);
                    }
                    break;
                case QUANT:
                    if (LoopOverTrait == 0) {
                        fprintf(fp, "T  %s\n", locname);
                    }
                    break;
                case TYPE_UNSET:
                default:
                    break;
                }
            }
        }
        /* Write the covariates if necessary */
        if (model_file) {
            if (LoopOverTrait == 1 && LTop->Locus[*trp].Type == AFFECTION
                && LTop->Pheno[*trp].Props.Affection.ClassCnt > 1) {
                fprintf(fp, "C  %s_covar\n", LTop->Locus[*trp].Name);
            } else {
                for (locus=0; locus < num_traits; locus++) {
                    SKIP_TRI(locus)
                        if (LTop->Locus[global_trait_entries[locus]].Type == AFFECTION
                            && LTop->Pheno[global_trait_entries[locus]].Props.Affection.ClassCnt > 1) {
                            fprintf(fp, "C  %s_covar\n", LTop->Locus[global_trait_entries[locus]].Name);
                        }
                }
            }
        }
        fclose(fp);
        if (nloop == 1) break;
        trp++;

    }
}

/* A single overall shell script is created for all traits and all chromosomes
   Then, for each trait and each chromosome, we also create a separate
   shell-script.
*/
static void merlin_shell_script(int numchr, int global,
				char *output_files[],
				char *opts,
				linkage_locus_top *LTop,
				int qtdt,
				merlin_opt_type merlin_opt)
{

    /* pedigree_file = output_files[10]
       locus_file = output_files[11]
       map_file = output_files[12]
       shell_script = output_files[13];
    */

    int nloop, tr, *trp, num_affec=num_traits;
    FILE *filep;
    char mfl[2*FILENAME_LENGTH];
    char *syscmd;
    char chr_str[3];
    char merlinoutfile[FILENAME_LENGTH];
    char merlinpdffile[FILENAME_LENGTH];
    char merlintablefile[FILENAME_LENGTH];
    char global_shell[FILENAME_LENGTH];
    char merlin_prog[7];
    int xlinked;
    int tabulate=0;
#ifdef RUNSHELL_SETUP
    char            merlin_prog_env[8];
#endif /* RUNSHELL_SETUP */

    xlinked=(((LTop->SexLinked == 2 && numchr == SEX_CHROMOSOME) ||
              (LTop->SexLinked == 1))? 1 : 0);

#ifdef RUNSHELL_SETUP
    strcpy(merlin_prog_env, (xlinked ? "_MINX" : "_MERLIN"));
#endif /* RUNSHELL_SETUP */
    strcpy(merlin_prog, (xlinked? "minx" : "merlin"));

    if (strcmp(output_paths[0], ".")) {
        sprintf(global_shell, "%s/%s", output_paths[0], output_files[13]);
    } else {
        strcpy(global_shell, output_files[13]);
    }
    change_output_chr(global_shell, 0);

    NLOOP;

    CHR_STR(numchr, chr_str);
    MERLINOUT(chr_str, merlinoutfile);
    MERLINPDF(chr_str, merlinpdffile);
    MERLINTABLE(chr_str, merlintablefile);

    trp = (num_traits > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && global == 1) {
            /*  create a top level shell script */
            if ((filep = fopen(global_shell, "w")) == NULL) {
                errorvf("could not open shell file %s.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "#!/bin/csh -f\n");
            /* print some identification */
            fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    global_shell
#endif
                );
            script_time_stamp(filep);
            fclose(filep);
            syscmd = CALLOC((strlen(global_shell) + 10), char);
            sprintf(syscmd, "chmod +x %s", global_shell);
            System(syscmd);
            free(syscmd);
        }

        if (nloop > 1 && tr == 0) {
            continue;
        }

        sprintf(mfl, "%s/%s", output_paths[tr], file_names[13]);
        if ((filep = fopen(mfl, "w")) == NULL) {
            errorvf("Unable to open %s for writing.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "#!/bin/csh -f\n");
        /* print some identification */
        fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                NOPATH
#else
                file_names[13]
#endif
            );
        script_time_stamp(filep);
        fprintf(filep, "# Chromosome number:    %d\n", numchr);
#ifdef RUNSHELL_SETUP
        // This handles the environment variable setup to allow the checking
        // functions in 'batch_run' to work correctly...
        fprintf_env_checkset_csh(filep, merlin_prog_env, merlin_prog);
        fprintf_env_checkset_csh(filep, "_MERLIN2SW2", "merlin2sw2.pl");
        fprintf_env_checkset_csh(filep, "_SIMWALK2", "simwalk2");
#endif /* RUNSHELL_SETUP */
        fprintf(filep, " if (-e %s) then\n", merlinoutfile);
        fprintf(filep, "    /bin/rm %s\n", merlinoutfile);
        fprintf(filep, "endif\n");
        fprintf(filep, " if (-e %s) then\n", merlintablefile);
        fprintf(filep, "    /bin/rm %s\n", merlintablefile);
        fprintf(filep, "endif\n");
        fprintf(filep, " if (-e %s) then\n", merlinpdffile);
        fprintf(filep, "    /bin/rm %s\n", merlinpdffile);
        fprintf(filep, "endif\n");

        fprintf(filep,
                "echo Running merlin on %s and  %s \n",
                output_files[10], output_files[11]);
#ifdef RUNSHELL_SETUP
        fprintf(filep, "$%s -p %s -d %s ", merlin_prog_env, output_files[10], output_files[11]);
#else /* RUNSHELL_SETUP */
        fprintf(filep, "%s -p %s -d %s ", merlin_prog, output_files[10], output_files[11]);
#endif /* RUNSHELL_SETUP */

        if (qtdt) {
            fprintf(filep, "-m %s -f %s ", output_files[12], output_files[9]);
        }

        if (opts == NULL) {
            if (trp==NULL) {
                fprintf(filep, "--ibd >> %s\n", merlinoutfile);
            } else if (LTop->Locus[*trp].Type == QUANT) {
                fprintf(filep, " --qtl --tabulate >> %s\n", merlinoutfile);
                tabulate=1;
                trp++;
            } else if (LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "--npl --pairs --tabulate >> %s\n", merlinoutfile);
                tabulate=1;
                trp++;
            }
        } else {
            /* Add tabulate if one of the following options are selected */
            tabulate = merlin_opt.tabulate;

            if ((merlin_opt.npl == 1 || merlin_opt.pairs == 1 || merlin_opt.qtl == 1 ||
                 merlin_opt.parametric == 1 || merlin_opt.vc == 1)) {
                /* tabulate should always be 1 */
                if (merlin_opt.tabulate != 1) {
                    fprintf(filep, "%s --tabulate >> %s\n", opts, merlinoutfile);
                } else {
                    fprintf(filep, "%s >> %s\n", opts, merlinoutfile);
                }
                tabulate = 1;
            } else {
                fprintf(filep, "%s >> %s\n", opts, merlinoutfile);
            }
        }
	fprintf_status_check_csh(filep, "merlin", 0);

        fprintf(filep, "cat %s\n", merlinoutfile);
        fprintf(filep, "echo ===================== >> %s\n",
                merlinoutfile);
        fprintf(filep, "cat %s >> %s\n", output_files[12],
                merlinoutfile);
        fprintf(filep, "echo Created output file %s\n",
                merlinoutfile);

        if (tabulate == 1) {
            if (merlin_opt.pairs == 1 || merlin_opt.npl == 1 ||
                merlin_opt.qtl == 1) {
                fprintf(filep, "if (-e merlin-nonparametric.tbl) then\n");
                fprintf(filep, "   mv merlin-nonparametric.tbl %s\n", merlintablefile);
                fprintf(filep, "   echo Renamed merlin-nonparametric.tbl to %s\n",
                        merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");
            } else if (merlin_opt.parametric == 1) {
                fprintf(filep, "if (-e merlin-parametric.tbl) then\n");
                fprintf(filep, "   mv merlin-parametric.tbl %s\n", merlintablefile);
                fprintf(filep, "   echo Renamed merlin-parametric.tbl to %s\n", merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");
            } else if (merlin_opt.vc == 1) {
                fprintf(filep, "if (-e merlin-vc-chr%2s.tbl) then\n", chr_str);
                fprintf(filep, "   mv merlin-vc-chr%2s.tbl %s\n", chr_str, merlintablefile);
                fprintf(filep, "   echo Renamed merlin-vc-chr%2s.tbl to %s\n", chr_str, merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");

            } else {
                sprintf(err_msg, "Tabulate option not available for selected options %s.\n",
                        opts);
                errorf(err_msg);
            }
        }
        if (merlin_opt.pdf == 1) {
            fprintf(filep, "if (-e merlin.pdf) then \n");
            fprintf(filep, "   mv merlin.pdf %s\n", merlinpdffile);
            fprintf(filep, "   echo Renamed merlin.pdf to %s\n", merlinpdffile);
            fprintf(filep, "else \n");
            fprintf(filep, "   echo Failed to create merlin pdf file!\n");
            fprintf(filep, "endif\n");

        }
        fclose(filep);
        syscmd = CALLOC((strlen(mfl) + 10), char);
        sprintf(syscmd, "chmod +x %s", mfl);
        System(syscmd);
        free(syscmd);

        if (global) {
            if ((filep = fopen(global_shell, "a")) == NULL) {
	        errorvf("Unable to open %s for writing.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "echo Running Merlin shell script %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    mfl
#endif
                );
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd %s\n", trait_paths[tr]);
            }
            fprintf(filep, "./%s\n", file_names[13]);
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd ..\n");
            }
            fclose(filep);
        }

        if (nloop == 1) break;

    }

    return;
}


static int merlin_option_check(char *option,
			       merlin_opt_type *merlin_opt,
			       char *file_names[],
			       int *create_merlin_model)

{
    char opt_copy[FILENAME_LENGTH];
    char opt_output[FILENAME_LENGTH];
    char *nextopt, prevopt[100];
    float num;
    int valid, i;
    char y[10];


    char knownopts[NUM_MERLIN_OPTS][20] =
        {
            "--error", "--information", "--likelihood", "--model",
            "--ibd", "--kinship", "--matrices", "--extended", "--select",
            "--npl", "--pairs", "--qtl", "--deviates", "--exp",
            "--vc", "--useCovariates", "--ascertainment", "--unlinked",
            "--infer", "--assoc", "--fastAssoc", "--filter", "--custom",
            "--best", "--sample", "--all", "--founders", "--horizontal",
            "--zero", "--one", "--two", "--three", "--singlepoint",
            "--steps", "--maxStep", "--minStep", "--grid", "--start", "--stop",
            "--clusters", "--distance", "--rsq", "--cfreq",
            "--bits", "--megabytes", "--minutes",
            "--trim", "--noCoupleBits", "--swap", "--smallSwap",
            "--quiet", "--markerNames", "--frequencies", "--perFamily", "--pdf",
            "--tabulate", "--prefix",
            "--simulate", "--reruns", "--save", "--trait"
        };

    init_merlin_opts(merlin_opt);


    strcpy(opt_copy, option);
    nextopt=strtok(opt_copy, " ");
    strcpy(opt_output, "");

    while (nextopt != NULL) {
        if (!strncmp(nextopt, "--",(size_t) 2)) {
            valid=0;
            for (i=0; i < NUM_MERLIN_OPTS; i++) {
                if (!strcmp(nextopt, knownopts[i])) {
                    valid=1; break;
                }
            }
            if (!valid) { sprintf(err_msg, "Unknown option %s.\n", nextopt); warnf(err_msg);}
        }
        nextopt = strtok(NULL, " ");
    }

    /* Now re-initialize the string */
    strcpy(opt_copy, option);
    nextopt=strtok(opt_copy, " ");

    while(nextopt != NULL) {
        strcat(opt_output, nextopt);
        strcat(opt_output, " ");
        /* Check that these options are followed by another arguments */
        if (!strcmp(nextopt, "--unlinked") || !strcmp(nextopt, "--steps") ||
            !strcmp(nextopt, "--maxStep") || !strcmp(nextopt, "--minStep") ||
            !strcmp(nextopt, "--grid") || !strcmp(nextopt, "--start") ||
            !strcmp(nextopt, "--stop") || !strcmp(nextopt, "--bits") ||
            !strcmp(nextopt, "--rsq") || !strcmp(nextopt, "--distance") ||
            !strcmp(nextopt, "--reruns") ||
            !strcmp(nextopt, "--megabytes") || !strcmp(nextopt, "--minutes")) {

            strcpy(prevopt, nextopt);
            nextopt = strtok(NULL, " ");
            if (nextopt != NULL) {
                if (sscanf(nextopt, "%f", &num) != 1) {
                    printf("ERROR: Option %s needs a numeric argument.\n", prevopt);
                    return 0;
                }
            } else {
                printf("ERROR: Option %s needs a numeric argument.\n", prevopt);
                return 0;
            }
            strcat(opt_output, nextopt); strcat(opt_output, " ");
        } else if (!strcmp(nextopt, "--model") || !strcmp(nextopt, "--custom")) {
            strcpy(prevopt, nextopt);
            /* Before reading in the next option we need to figure out whether
               this is a constant or another option */
            nextopt = strtok(NULL, " ");
            if (nextopt == NULL || ! strncmp(nextopt, "--", (size_t) 2)) {
                printf("ERROR: Option %s needs a file name argument.\n", prevopt);
                return 0;
            }

            if (!strcmp(prevopt, "--model")) {
                shorten_path(nextopt, file_names[14]);
                if (access(file_names[14], F_OK) != 0) {
                    if (*create_merlin_model) {
                        printf("File %s will be created from traits automatically.\n",
                               file_names[14]);
                    } else {
                        printf("WARNING: File %s not found, Merlin may not run correctly.\n",
                               file_names[14]);

                        printf("Create model file(s) %s in output folder? [y/n](default n)",
                               file_names[14]);
                        fflush(stdout);
                        (void)fgets(y, 9, stdin); newline; fflush(stdin);
                        if (y[0] == 'Y' || y[0] == 'y') {
                            *create_merlin_model=1;
                        } else {
                            printf("WARNING: No model file, Merlin may not run correctly.\n");
                        }
                    }
                }
                if (strcmp(file_names[14], "None")) {
                    strcat(opt_output, file_names[14]); strcat(opt_output, " ");
                }
                merlin_opt->parametric = 1;
            } else {
                strcat(opt_output, nextopt); strcat(opt_output, " ");
                if (access(nextopt, F_OK) != 0) {
                    printf("WARNING: file %s provided for %s does not exist or is unreadable.\n",
                           prevopt, nextopt);
                }
            }
        } else if (!strcmp(nextopt, "--prefix")) {
            nextopt = strtok(NULL, " ");
            if (nextopt == NULL || !strncmp(nextopt, "--", (size_t) 2)) {
                printf("ERROR: Option %s needs a string name argument.\n", prevopt);
                return 0;
            }
            strcat(opt_output, nextopt); strcat(opt_output, " ");
            warnf("Mega2 generated shell script recognizes only the default prefix \"merlin\"");
            warnf("Merlin shell script will fail with alternate prefixes.");
            merlin_opt->prefix=1;
            /* otherwise process this string as an option */
        } else {
            if (!strcmp(nextopt, "--npl")) { merlin_opt->npl=1; }
            if (!strcmp(nextopt, "--pairs")) { merlin_opt->pairs=1; }
            if (!strcmp(nextopt, "--qtl")) { merlin_opt->qtl=1; }
            if (!strcmp(nextopt, "--vc")) { merlin_opt->vc=1; }
            if (!strcmp(nextopt, "--markerNames")) { merlin_opt->markernames=1; }
            if (!strcmp(nextopt, "--pdf")) { merlin_opt->pdf=1; }
            if (!strcmp(nextopt, "--tabulate")) { merlin_opt->tabulate=1; }
        }
        nextopt = strtok(NULL, " ");
    }

    strcpy(option, opt_output);

    return 1;

}

static void init_merlin_opts(merlin_opt_type *merlin_opt)
{
    merlin_opt->npl=
        merlin_opt->qtl =
        merlin_opt->pairs =
        merlin_opt->vc =
        merlin_opt->parametric =
        merlin_opt->pdf =
        merlin_opt->markernames =
        merlin_opt->prefix =
        merlin_opt->tabulate = 0;
}

static void merlin_options(char *option,
			   merlin_opt_type *merlin_opt,
			   char *file_names[],
			   int *create_liability_model)

{
    int done_ = 0;
    char cdone_[10];


    if (HasAff) {
        strcpy(option, "--npl --pairs --tabulate");
    } else if (HasQuant) {
        strcpy(option, "--qtl --tabulate");
    } else {
        strcpy(option, "--ibd");
    }

    *create_liability_model = 0;

    if (DEFAULT_OPTIONS) {
        /* This is called only to set the correct flags in merlin_opt structure */
        merlin_option_check(option, merlin_opt, file_names, create_liability_model);
    } else {
        done_ = -1;
        while(done_) {
            /* New merlin run parameters menu */
            draw_line();
            printf("Merlin run options menu:\n");
            printf("0) Done with this menu, please proceed.\n");
            printf(" 1) Select Merlin option [%s].\n", option);
            if (strcmp(file_names[14], "None")) {
                printf(" 2) %s Merlin model file %s from affection trait parameters [%s].\n",
                       ((access(file_names[14], F_OK) == 0)?
                        "Recreate" : "Create"),
                       file_names[14],
                       yorn[*create_liability_model]);
            } else {
                printf(" 2) Create Merlin model file from affection trait parameters [%s].\n",
                       yorn[*create_liability_model]);
            }
            if (*create_liability_model && strcmp(file_names[14], "None")) {
                printf("    To use this model file, use the Merlin analysis option:\n");
                printf("       --model %s\n", file_names[14]);
            }
            printf("Select 0-2 (2 to toggle) > ");
            fcmap(stdin, "%s", cdone_);
            newline;
/*       done_ = atoi(cdone_); */
            done_ = (int)strtol(cdone_, (char **)NULL, 10);
/*       printf("%d\n", done_); */
/*       sleep(1); */
            switch(done_) {
            case 0:
                if (!merlin_option_check(option, merlin_opt, file_names, create_liability_model)) {
                    printf("Error: bad option string, discarding options.\n");
                    strcpy(option, "");
                }
                break;
            case 1:
                draw_line();
                printf("Merlin option selection:\n");
                printf("Enter one or more of the options listed below, ");
                printf("separated by spaces, \n");
                printf("   e.g. --steps 2 --bits 32 \n");
                printf("Some options must be followed by a number as indicated\n");

                printf("For no options press <RETURN> or <Enter>\n\n");

                printf("         General : --error, --information, --likelihood, --model\n");
                printf("      IBD States : --ibd, --kinship, --matrices, --extended, --select\n");
                printf("     NPL Linkage : --npl, --pairs, --qtl, --deviates, --exp\n");
                printf("     VC Linkage  : --vc, --useCovariates, --ascertainment\n");
                printf("     Association : --infer, --assoc, --fastAssoc, --filter, --custom\n");
                printf("     Haplotyping : --best, --sample, --all, --founders, --horizontal\n");
                printf("   Recombination : --zero, --one, --two, --three, --singlepoint\n");
                printf("       Positions : --steps, --maxStep, --minStep, --grid, --start, --stop\n");
                printf(" Marker Clusters : --clusters, --distance, --rsq, --cfreq\n");
                printf("          Limits : --bits, --megabytes, --minutes\n");
                printf("     Performance : --trim, --noCoupleBits, --swap, --cache\n");
                printf("          Output : --quiet, --markerNames, --frequencies, --perFamily, --pdf,\n");
                printf("                   --prefix\n");
                printf("      Simulation : --simulate, --reruns, --save\n");

                printf("These options can be followed by one additional argument each:\n");
                printf(" model, custom, steps, maxStep, minStep, grid, start, stop, \n");
                printf(" clusters, megabytes, minutes, cache, prefix\n");
                printf("\nEnter options as text (at most 200 characters) > ");
                strcpy(option, "");
                fflush(stdout);
                (void)fgets(option, FILENAME_LENGTH-1, stdin); newline;
                option[strlen(option)-1]='\0';
                if (!merlin_option_check(option, merlin_opt, file_names, create_liability_model)) {
                    printf("Error: bad option string, discarding options.\n");
                    strcpy(option, "");
                }
                break;
            case 2:
                printf("Set value to \"yes\" to automatically generate model file\n");
                printf("using affection status trait locus parameters from input files.\n");
                *create_liability_model = TOGGLE(*create_liability_model);
                break;
            default:
                warn_unknown(cdone_);
                break;
            }
        }
    }
}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void create_merlin_files(linkage_ped_top **LPedTop, int *numchr,
			 file_format  *infl_type,
			 char *file_names[], int  untyped_ped_opt,
			 analysis_type analysis)
{
    int i, num_affec=0;
    int R_setup, global;
    char opts[FILENAME_LENGTH];
    int qtdt, xlinked, model_file=0;
    merlin_opt_type merlin_opt;
    /*  linkage_ped_top *Copy; */
    linkage_ped_tree *LPed;


/*    qtdt = ((Top->LocusTop->PedRecDataType == Premakeped || */
/*  	   Top->LocusTop->PedRecDataType == Postmakeped)? 0 : 1); */

    qtdt=1;
    strcpy(opts, "");
    init_merlin_opts(&merlin_opt);

    /* Input merlin options */
    if (analysis == TO_MERLINONLY) {
        merlin_options(&(opts[0]), &(merlin_opt), file_names, &model_file);

        if (!DEFAULT_OPTIONS) {
            if (strcmp(opts, "")) {
                sprintf(err_msg, "User options to Merlin: %s.", opts);
                mssgf(err_msg);
            } else {
                if (HasAff) {
                    sprintf(opts, "--npl --pairs");
                    mssgf("Using default Merlin options --npl and --pairs");
                    merlin_opt.npl=merlin_opt.pairs=1;
                    merlin_opt.markernames=0;
                } else if (HasQuant) {
                    sprintf(opts, "--qtl");
                    merlin_opt.qtl=1;
                    mssgf("Using default Merlin option --qtl.");
                } else {
                    sprintf(opts, "--ibd");
                    mssgf("Using default Merlin option --ibd.");
                }
            }
            log_line(mssgf);
        }
    }

    for (i=0; i < num_traits; i++) {
        SKIP_TRI(i)
            if ((*LPedTop)->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION ||
               (*LPedTop)->LocusTop->Locus[global_trait_entries[i]].Type == QUANT) {
                num_affec++;
            }
    }

    if ((analysis == TO_MERLINONLY) &&
        (merlin_opt.npl==1 ||  merlin_opt.pairs==1 ||
         merlin_opt.qtl==1 ||  merlin_opt.vc==1 ||
         merlin_opt.parametric == 1)) {
        R_setup = 1;
    } else {
        R_setup = 0;
    }

    merlin_file_names(file_names, OrigIds, analysis, Top->OrigIds,
                      Top->UniqueIds, &R_setup, num_affec, qtdt);

    if (*infl_type == LINKAGE) {
        for (i = 0; i < Top->PedCnt; i++)  {
            LPed = &(Top->Ped[i]);
            if (LPed->Loops != NULL)   {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)    {
                    errorf("FATAL: Aborting.");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    global = ((main_chromocnt > 1 || (num_traits > 0 && LoopOverTrait == 1))? 1: 0);

    if (R_setup) {
        R_setup = merlin_R_setup(*numchr, *LPedTop, merlin_opt);
        if (!merlin_opt.markernames) {
            strcat(opts, " --markerNames");
        }
    }
    /* for creating subsets of data */
/*   Copy=new_lpedtop(); */
/*   copy_linkage_ped_top(Top, Copy, 1); */

    for (i=0; i < main_chromocnt; i++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[i];
            analysis->replace_chr_number(file_names, *numchr);
/*       free_all_including_lpedtop(&Top); */
/*       Top=new_lpedtop(); */
/*       copy_linkage_ped_top(Copy, Top,1); */
/*       ReOrderMappedLoci(Top, numchr); */
            draw_line();
        }
        get_loci_on_chromosome(*numchr);
        xlinked = (((Top->LocusTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
                    (Top->LocusTop->SexLinked == 1))? 1 : 0);

        if (xlinked) {
            if (analysis == TO_MERLIN && xlinked) {
                warnf("Files not created for X-Chromosome.");
                continue;
            }
        }

        omit_peds(untyped_ped_opt, Top);
        save_premakeped_peds(file_names[10], Top, analysis, model_file);
        if (qtdt) {
            write_qtdt_locus_file(file_names[11], Top->LocusTop, analysis, model_file);
        } else {
            write_gh_locus_file(file_names[11], Top->LocusTop, analysis, xlinked);
        }

        if (R_setup || qtdt) {
            write_merlin_map(Top->LocusTop, file_names[12], *numchr, analysis);
        }
        if (qtdt) {
            write_merlin_freq(Top->LocusTop, file_names[9], *numchr, analysis);
        }
        if (analysis == TO_MERLIN) {
            merlin_output_ped_order(Top, file_names[13]);
        } else {
            merlin_shell_script(*numchr, global, file_names, opts,
                                Top->LocusTop, qtdt, merlin_opt);

            if (global) global++;
        }

        if (R_setup && !qtdt) {
            log_line(mssgf);
            mssgf("  Merlin map file will be used for generating plots only.");
            mssgf("  Merlin will use thetas from the locus file.");
        }

        if (i==0)       create_mssg(analysis);

        sprintf(err_msg, "        Pedigree file:         %s",
                file_names[10]);
        mssgf(err_msg);
        sprintf(err_msg, "        Locus file:            %s",
                file_names[11]);
        mssgf(err_msg);
        if (R_setup || qtdt) {
            sprintf(err_msg, "        Map file:              %s",
                    file_names[12]);
            mssgf(err_msg);
        }
        if (qtdt) {
            sprintf(err_msg, "        Frequency file:        %s",
                    file_names[9]);
            mssgf(err_msg);
        }
        if (analysis == TO_MERLIN) {
            sprintf(err_msg, "        Order file:            %s",
                    file_names[13]);
        } else {
            sprintf(err_msg, "        C-shell script:        %s",
                    file_names[13]);
        }
        mssgf(err_msg);
    }


    if (analysis != TO_MERLIN) {
        if (global || R_setup || model_file) {
            mssgf("Additional files: ");
        }

        if (model_file) {
            merlin_model_file(file_names[14], Top->LocusTop);
            if (LoopOverTrait == 1 && num_traits >= 2) {
                sprintf(err_msg, "        Trait-specific model file: %s",
                        file_names[14]);
            } else {
                sprintf(err_msg, "        Combined model file:    %s",
                        file_names[14]);
            }
            mssgf(err_msg);
        }
        if (global) {
            change_output_chr(file_names[13], 0);
            sprintf(err_msg, "        Combined shell script: %s", file_names[13]);
            mssgf(err_msg);
        }
        if (R_setup) {
            append_R_commands(file_names[13], Rshell_file);
            sprintf(err_msg, "        R-script file:         %s", Rscript_file);
            mssgf(err_msg);
            sprintf(err_msg, "        Shell file to run R:   %s", Rshell_file);
            mssgf(err_msg);
            if (base_pair_position_index >= 0) {
	         // Needed to produce UCSC Custom Tracks which is only done if a physical map file exists.
                 sprintf(err_msg, "        Physical map for R:    %s", Rmap_file);
                 mssgf(err_msg);
            }
        }

        if (R_setup) {
            log_line(mssgf);
            sprintf(err_msg,
                    "Run %s in %s for all analysis and creation of R-plots",
                    file_names[13],
#ifdef HIDEPATH
                    NOPATH
#else
                    ((strcmp(output_paths[0], ".")?
                      output_paths[0] : "current directory"))
#endif
                );
            mssgf(err_msg);
        }
    }

    if (merlin_opt.prefix == 1) {
        printf("***************\n");
        warnf("Merlin shell script will fail if prefix selected is not the default \"merlin\".");
        printf("***************\n");
    }

    /*  draw_line(); */

    /*  free_all_including_lpedtop(&Copy); */
}

#undef Top
