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

/* Functions for writing out Mega2 annotated formatted files,
   for now only handles integer pedigree and person IDs. Works
   like l2a.py otherwise.
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "annotated_ped_file_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
 annotated_ped_file_ext.h:  write_annotated_aff write_annotated_numbered write_annotated_quant
     error_messages_ext.h:  mssgf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  chrom_num_to_name
            linkage_ext.h:  get_loci_on_chromosome get_unmapped_loci set_link_IDs switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line global_ou_files
*/


void write_annotated_names_file(linkage_locus_top *LTop,
				char *mfl_name);

/* 1) function declarations:
   These are optional but avoid compiler warnings, if functions are
   called before they are declared, and are good for clarity. */

/* 1a) Single entry function that is called from mega2.c */
void create_mega2annot_files(linkage_ped_top **LPedTop, char *file_names[],
                             int untyped_ped_opt);

/* 1b) Local function definitions
   static f1()
   static f2() etc.

*/

/* static void marker_name_width(linkage_locus_top *LTop, int *MrkWidth)  */

/* { */

/*   int m; */

/*   *MrkWidth=0; */
/*   for(m=0; m < LTop->LocusCnt; m++) { */
/*     *MrkWidth =  */
/*       ((strlen(LTop->Locus[m].Name) > *MrkWidth)? */
/*        strlen(LTop->Locus[m].Name) :  *MrkWidth); */
/*   } */
/*   return; */
/* } */


/* This is a function that implements the file-names menu for MEGA2ANNOT,
 */

static void rename_mega2annot_locus(linkage_locus_top *LTop)

{
    int i;
    char *dot;

    for (i=0; i < NumChrLoci; i++) {
        dot = strchr(LTop->Locus[ChrLoci[i]].Name, '.');
        if (dot != NULL) {
            printf("%s\n", dot);
            *dot = '_';
        }
    }
}

static void mega2annot_file_names(char *file_names[], int *combine_chromo,
				  int has_orig, int has_uniq)
{
    int item, choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12];
    analysis_type analysis = MEGA2ANNOT;
    int pen_it, pedid_it, perid_it;

    if (DEFAULT_OUTFILES)  {
        mssgf("Output file names set to defaults.");
        // don't display the menu for the user...
        choice = 0;
        //return;
    }

    // If there is more than one chromosome specified defer to the batch file item
    // or if it is not read, ask the user...
    if (main_chromocnt > 1) {
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].item_read) {
            *combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
        } else {
            // Ask the user if they want chromosome specific files.
            //  0 == One set of output files per chromosome (default).
            //  1 == Combine chromosomes into one output file.
            // In many of the other analysis modes this question is "bundled" in with the menu below.
            if (choice != 0) *combine_chromo=global_ou_files(1);
        }
    }

    pen_it = pedid_it = perid_it = -1;

    strcpy(stem, (num_traits > 1 && LoopOverTrait) ? "in each subdir" : "");

    if (main_chromocnt > 1) {
        if (!(*combine_chromo)) {
            strcpy(stem, "stem");
            // returns only the <file_name> removing <extension> and <rest> if they exist...
            analysis->replace_chr_number(file_names, -9);
        } else {
            // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
            analysis->replace_chr_number(file_names, 0);
        }
    }

    while (choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("  Mega2 annotated file name menu\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt == 1 || *combine_chromo) {
            printf(" 1) Pedigree filename:                   %-15s\t%s\n",
                   file_names[0], file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:                     %-15s\t%s\n",
                   file_names[1], file_status(file_names[1], fl_stat));
            printf(" 3) Map file name:                       %-15s\t%s\n",
                   file_names[2], file_status(file_names[2], fl_stat));
            printf(" 4) Frequency file name:                 %-15s\t%s\n",
                   file_names[3], file_status(file_names[3], fl_stat));
            if (HasAff) {
                printf(" 5) Penetrance file name:                %-15s\t%s\n",
                       file_names[4], file_status(file_names[4], fl_stat));
                pen_it = item = 5;
            } else {
                item=4;
            }
        } else {
            printf(" 1) Pedigree file name %s:          %-15s\n", stem,
                   file_names[0]);
            printf(" 2) Locus file name %s:             %-15s\n", stem,
                   file_names[1]);
            printf(" 3) Map file name %s:               %-15s\n", stem,
                   file_names[2]);
            printf(" 4) Frequency file name %s:         %-15s\t\n", stem,
                   file_names[3]);

            if (HasAff) {
                printf(" 5) Penetrance file name %s:        %-15s\n", stem,
                       file_names[4]);
                pen_it = item = 5;
            } else {
                item=4;
            }
        }

        item++;
        perid_it = item;
        individual_id_item(item, analysis, OrigIds[0], 31, 2, has_orig, has_uniq);

        item++;
        pedid_it = item;
        pedigree_id_item(item, analysis, OrigIds[1], 26, 2, has_orig);

        printf("Select options 0-%d > ", item);

        choice = -1;
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        case 2:
            printf("New locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]);
            printf("\n");
            break;
        case 1:
            printf("New pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);
            printf("\n");
            break;
        case 3:
            printf("New map file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]);
            printf("\n");
            break;
        case 4:
            printf("New frequency file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]);
            printf("\n");
            break;
        default:
            if (choice == perid_it) {
                OrigIds[0] =
                    individual_id_item(0, TO_PREMAKEPED, OrigIds[0], 35, 1,
                                       has_orig, has_uniq);
            } else if (choice == pedid_it) {
                OrigIds[1] = pedigree_id_item(0, TO_PREMAKEPED,
                                              OrigIds[1], 35, 1, has_orig);
            } else if (choice == pen_it) {
                printf("New penetrance file name %s > ", stem);
                fcmap(stdin, "%s", file_names[4]);
                printf("\n");
            } else {
                warn_unknown(cchoice);
            }
            break;
        }
    } // while (choice != 0) {

    individual_id_item(0, TO_PREMAKEPED, OrigIds[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, TO_PREMAKEPED, OrigIds[1], 0, 3, has_orig);
}

void write_annotated_names_file(linkage_locus_top *LTop,
				char *mfl_name)
{

    int nloop, tr, num_affec=num_traits ;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;
    int m1, m;
    int *trp;

    NLOOP;

    trp = (num_traits > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr = 0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], mfl_name);
        if ((fp = fopen(mfl, "w")) == NULL) {
	    errorvf("Could not open names file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        /* print the header */
        fprintf(fp, "Type Name\n");

        if (LoopOverTrait == 1) {
            /* trait locus first */
            if (LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(fp, "A  ");
            } else {
                fprintf(fp, "T  ");
            }
            fprintf(fp, "%s\n", LTop->Locus[*trp].Name);
            trp++;
        }

        for (m1 = 0; m1 < NumChrLoci; m1++) {
            m = ChrLoci[m1];
            switch (LTop->Locus[m].Type) {
            case NUMBERED:
                if (LTop->Locus[m].chromosome == SEX_CHROMOSOME) {
                    fprintf(fp, " X  %s\n", LTop->Locus[m].Name);
                } else if (LTop->Locus[m].chromosome == MALE_CHROMOSOME) {
                    fprintf(fp, " Y  %s\n", LTop->Locus[m].Name);
                } else {
                    fprintf(fp, " M  %s\n", LTop->Locus[m].Name);
                }
                break;
            case AFFECTION:
                if (LoopOverTrait == 0) {
                    if (LTop->Locus[m].Data.Affection.ClassCnt > 1) {
                        fprintf(fp, " L  %s\n", LTop->Locus[m].Name);
                    } else {
                        fprintf(fp, " A  %s\n", LTop->Locus[m].Name);
                    }
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0 || LTop->Locus[m].Class == COVARIATE) {
                    if (LTop->Locus[m].Class == COVARIATE) {
                        fprintf(fp, " C  %s\n", LTop->Locus[m].Name);
                    } else {
                        fprintf(fp, " T  %s\n", LTop->Locus[m].Name);
                    }
                }
                break;
            default:
                break;
            }
        }
        fclose(fp);
        if (nloop == 1)
            break;
    }
    return;
}

static void annotated_ped_file(char *outfl_name, linkage_ped_top *Top,
			       int pedfile_type)

{

    int tr, nloop, num_affec=num_traits;
    int ped, entry, locus;
    int *trp;
    linkage_ped_rec *Entry;
    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];
    char pformat[6], fformat[6];
    int p, pwid, fwid;
    linkage_locus_top *LTop;

    NLOOP;

    field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);
    create_formats(fwid, pwid, fformat, pformat);

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
        /* print the header */
        fprintf(filep, "Pedigree ID Father Mother ");
        if (pedfile_type == POSTMAKEPED_PFT) {
            fprintf(filep, "FirstOff NextMatSib NextPatSib ");
        }
        fprintf(filep, "Sex ");
        if (pedfile_type == POSTMAKEPED_PFT) {
            fprintf(filep, "Proband ");
        }
        if (LoopOverTrait == 1) {
            switch(LTop->Locus[*trp].Type) {
            case AFFECTION:
                ((LTop->Locus[*trp].LAFFDATA.ClassCnt > 1)?
                 fprintf(filep, "%s.L.1 %s.L.2 ", LTop->Locus[*trp].Name,  LTop->Locus[*trp].Name) :
                 fprintf(filep, "%s.A ", LTop->Locus[*trp].Name));
                break;
            case QUANT:
                fprintf(filep, "%s.T ", LTop->Locus[*trp].Name);
                break;
            default:
                break;
            }
        }
        /* Now write the phenotype column headers */
        for (locus = 0; locus < NumChrLoci; locus++) {
            linkage_locus_rec *Loc = &(LTop->Locus[ChrLoci[locus]]);
            char ext;

            switch(Loc->Type) {
            case AFFECTION :
                if (LoopOverTrait == 0) {
                    (Loc->LAFFDATA.ClassCnt > 1)?
                        fprintf(filep, "%s.L.1 %s.L.2 ",
                                Loc->Name,  Loc->Name) :  fprintf(filep, "%s.A ", Loc->Name);
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0) {
                    fprintf(filep, "%s.T ", Loc->Name);
                }
                break;
            case NUMBERED:
                if (Loc->chromosome == SEX_CHROMOSOME) {
                    ext = 'X';
                } else if (Loc->chromosome == MALE_CHROMOSOME) {
                    ext = 'Y';
                } else {
                    ext = 'M';
                }
                fprintf(filep, "%s.%c.1 %s.%c.2 ", Loc->Name,  ext, Loc->Name, ext);
                break;

            default:
                break;
            }
        }

        fprintf(filep, "PedID  PerID\n");

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
                } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                    fprintf(filep, pformat, Entry->UniqueID);
                } else {
                    fprintf(filep, pformat, Entry->ID);
                }

                if (Entry->Father != 0) {
                    if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].OrigID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].OrigID);
                    } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].UniqueID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].UniqueID);
                    } else {
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Father-1].ID);
                        fprintf(filep, pformat,
                                Top->Ped[ped].Entry[Entry->Mother-1].ID);
                    }
                } else {
                    /* create the 0 strings */
                    for(p=0; p < (pwid-1); p++)    fprintf(filep, " ");
                    fprintf(filep, "0 ");
                    for(p=0; p< (pwid-1); p++)    fprintf(filep, " ");
                    fprintf(filep, "0 ");
                }

                if (pedfile_type == POSTMAKEPED_PFT) {
                    if (Entry->First_Offspring != 0) {
                        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->First_Offspring-1].OrigID);
                        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->First_Offspring-1].UniqueID);
                        } else {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->First_Offspring-1].ID);
                        }
                    } else {
                        /* create the 0 strings */
                        for(p=0; p < (pwid-1); p++)    fprintf(filep, " ");
                        fprintf(filep, "0 ");
                    }


                    if (Entry->Next_PA_Sib != 0) {
                        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_PA_Sib-1].OrigID);
                        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_PA_Sib-1].UniqueID);
                        } else {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_PA_Sib-1].ID);
                        }
                    } else {
                        for (p=0; p< (pwid-1); p++)    fprintf(filep, " ");
                        fprintf(filep, "0 ");
                    }

                    if (Entry->Next_MA_Sib != 0) {
                        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_MA_Sib-1].OrigID);
                        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_MA_Sib-1].UniqueID);
                        } else {
                            fprintf(filep, pformat,
                                    Top->Ped[ped].Entry[Entry->Next_MA_Sib-1].ID);
                        }
                    } else {
                        for (p=0; p< (pwid-1); p++)    fprintf(filep, " ");
                        fprintf(filep, "0 ");
                    }
                }

                /* write the sex */
                if (Entry->Sex == MALE_ID) fprintf(filep, " %d", MALE_ID);
                else if (Entry->Sex == FEMALE_ID) fprintf(filep, " %d", FEMALE_ID);
                else fprintf(filep, " 0");

                if (pedfile_type == POSTMAKEPED_PFT) {
                    /* write the proband */
                    fprintf(filep, " %d ", Entry->OrigProband);
                }

                /* write the affection locus for some analysis types */
                if (LoopOverTrait == 1) {
                    switch(LTop->Locus[*trp].Type) {
                    case AFFECTION:
                        write_annotated_aff(filep, *trp, &(LTop->Locus[*trp]),
                                            Entry);
                        break;
                    case QUANT:
                        write_annotated_quant(filep, *trp, Entry);
                        break;
                    default:
                        break;
                    }
                }
                /* Now write the genotype data */
                for (locus = 0; locus < NumChrLoci; locus++) {
                    linkage_locus_rec *Loc = &(LTop->Locus[ChrLoci[locus]]);
                    switch(Loc->Type) {
                    case QUANT:
                        if (LoopOverTrait == 0)  {
                            write_annotated_quant(filep, ChrLoci[locus], Entry);
                        }
                        break;

                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            write_annotated_aff(filep, ChrLoci[locus], Loc, Entry);
                        }
                        break;

                    case NUMBERED:
                        write_annotated_numbered(filep, ChrLoci[locus], Entry);
                        break;

                    default:
                        break;
                    }
                }

                /* finally, write the silly ped# and per# */
                fprintf(filep, "  %s  %s", Top->Ped[ped].Name, Entry->OrigID);
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }

    return;
}

static void annotated_map_file(linkage_locus_top *LTop, char *mfl_name)

{

    int nloop, tr, locus1, num_affec=num_traits ;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;
    char chrom_name[4];

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], mfl_name);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open map file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (LTop->map_distance_type == 'k' ||
            LTop->map_distance_type == 'h') {
	  if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            fprintf(fp,  "Name  Chromosome  Map.%c.a\n",
		    LTop->map_distance_type);
	  } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
		     genetic_distance_sex_type_map == FEMALE_GDMT) {
	    // We need a heading for the male map even though there may not be one...
	    fprintf(fp, "Name  Chromosome  Map.%c.f  Map.%c.m\n",
		    LTop->map_distance_type,
		    LTop->map_distance_type);
	  }
	} else {
            fprintf(fp,  "Name  Chromosome Map.p.a\n");
        }

        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            int locus = ChrLoci[locus1];
            if (LTop->Locus[locus].Type == NUMBERED ||
                LTop->Locus[locus].Type == BINARY) {

                if (LTop->Locus[locus].chromosome < SEX_CHROMOSOME) {
                    fprintf(fp, "%-15s %2d",
                            LTop->Locus[locus].Name,
                            LTop->Locus[locus].chromosome);
                } else {
                    fprintf(fp, "%-15s %s",
                            LTop->Locus[locus].Name,
                            chrom_num_to_name(LTop->Locus[locus].chromosome, &(chrom_name[0])));
                }

		// NOTE: the old version of mega2 would output three columns if there were
		// sex specific maps available. The new mega2 does not allow you to pick
		// all three (a, m, f) even if they exist.
		if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
		    fprintf(fp, " %10f\n",
			    LTop->Locus[locus].position);
		} else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
		    fprintf(fp, "    %10f     %10f\n", 
			    LTop->Locus[locus].pos_female,
			    LTop->Locus[locus].pos_male);
		} else if (genetic_distance_sex_type_map == FEMALE_GDMT &&
			   LTop->Locus[locus].chromosome != SEX_CHROMOSOME) {
		    errorvf("When a female only map is specified only operations on the X chromosome are permitted.\n");
		    EXIT(INPUT_DATA_ERROR);
		}
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
}



static void annotated_frequency_file(linkage_locus_top *LTop, char *freqfl)
{
    int nloop, tr, locus, locus1, num_affec=num_traits ;
    int all;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    for (tr = 0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0)
            continue;
        sprintf(mfl, "%s/%s", output_paths[tr], freqfl);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open frequency file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(fp, "Name Allele Frequency\n");
        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            locus = ChrLoci[locus1];
            for (all = 0; all < LTop->Locus[locus].AlleleCnt; all++) {
                fprintf(fp, "%s %d %8f\n", LTop->Locus[locus].Name, all + 1,
                        LTop->Locus[locus].Allele[all].Frequency);
            }
        }
        fclose(fp);
        if (nloop == 1)
            break;
    }

    return;

}

static void write_loc_pen(FILE *filep, linkage_locus_rec *Locus, linkage_affection_class *clp, int xlinked)
{
    int tmpi2, i;

    for (tmpi2 = 1; tmpi2 <= Locus->LAFFDATA.ClassCnt; tmpi2++) {
        if (xlinked) {
            if (xlinked == 2) {
                if (clp[tmpi2 - 1].AutoDef == 0) {
                    fprintf(filep, "%s %d autosomal", Locus->Name, tmpi2);
                    for (i = 0; i < Locus->LAFFDATA.PenCnt; i++)
                        fprintf(filep, " %f",  clp[tmpi2 - 1].AutoPen[i]);
                    fprintf(filep, "\n");
                }
            }

            if (clp[tmpi2 - 1].FemaleDef == 0) {
                fprintf(filep, "%s %d female", Locus->Name, tmpi2);
                for (i = 0; i < Locus->LAFFDATA.PenCnt; i++)
                    fprintf(filep, " %f",  clp[tmpi2 - 1].FemalePen[i]);
                fprintf(filep, "\n");
            }

            if (clp[tmpi2 - 1].MaleDef == 0) {
                fprintf(filep, "%s %d male", Locus->Name, tmpi2);
                for (i = 0; i < Locus->AlleleCnt; i++)
                    fprintf(filep, " %f",  clp[tmpi2 - 1].MalePen[i]);
                fprintf(filep, "\n");
            }
        } else {
            fprintf(filep, "%s %d", Locus->Name, tmpi2);
            for (i = 0; i < Locus->LAFFDATA.PenCnt; i++)
                fprintf(filep, " %f",  clp[tmpi2 - 1].AutoPen[i]);
            fprintf(filep, "\n");
        }
    }
}

static void annotated_penetrance_file(char *fl_name, linkage_locus_top * LTop, int xlinked)
{
    int       tr;
//    int	tmpi2;
    linkage_locus_rec *Locus;
    FILE      *filep;
    char      poutfl_name[2*FILENAME_LENGTH];
//    int xlinked = 0;
    linkage_affection_class *clp = NULL;


    if (LoopOverTrait == 1 && num_traits > 1) {
        for (tr=0; tr < num_traits; tr++) {
            Locus = &(LTop->Locus[global_trait_entries[tr]]);
            if (Locus->Type == AFFECTION) {
                sprintf(poutfl_name, "%s/%s", output_paths[tr+1], fl_name);
                if ((filep = fopen(poutfl_name, "w")) == NULL) {
                    errorvf("Could not open penetrance file %s.\n", poutfl_name);
                    EXIT(FILE_WRITE_ERROR);
                }

                clp = Locus->LAFFDATA.Class;
/*
                xlinked = 0;
                for (tmpi2 = 1; tmpi2 <= Locus->LAFFDATA.ClassCnt; tmpi2++)
                    xlinked = xlinked || clp[tmpi2-1].MaleDef != 1 || clp[tmpi2-1].FemaleDef != 1;
*/
                fprintf(filep, "Name   Class  %sPen.11    Pen.12    Pen.22\n", xlinked ? "Sex   " : "");

                write_loc_pen(filep, Locus, clp, xlinked);

                fclose(filep);}
            
        }
    } else {
        sprintf(poutfl_name, "%s/%s", output_paths[0], fl_name);
        if ((filep = fopen(poutfl_name, "w")) == NULL) {
            errorvf("Could not open penetrance file %s.\n", poutfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

/*
        xlinked = 0;
        for (tr=0; tr < num_traits; tr++) {
            if (global_trait_entries[tr] < 0) continue;

            Locus = &(LTop->Locus[global_trait_entries[tr]]);
            if (Locus->Type == AFFECTION) {
                clp = Locus->LAFFDATA.Class;
                for (tmpi2 = 1; tmpi2 <= Locus->LAFFDATA.ClassCnt; tmpi2++)
                    xlinked = xlinked || clp[tmpi2-1].MaleDef != 1 || clp[tmpi2-1].FemaleDef != 1;
            }
        }
*/
        fprintf(filep, "Name   Class  %sPen.11    Pen.12    Pen.22\n", xlinked ? "Sex   " : "");
        for (tr=0; tr < num_traits; tr++) {
            if (global_trait_entries[tr] < 0) continue;

            Locus = &(LTop->Locus[global_trait_entries[tr]]);
            if (Locus->Type == AFFECTION) {
                clp = Locus->LAFFDATA.Class;
                write_loc_pen(filep, Locus, clp, xlinked);
            }
        }
        fclose(filep);
    }
}

/* entry function */
#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void create_mega2annot_files(linkage_ped_top **LPedTop, char *file_names[],
			     int untyped_ped_opt)

{
    int i, fwid, pwid, mwid;
    int numchr;
    int xlinked;
    int combine_chromo=0;    // One set of output files per chromosome (default)...
    analysis_type analysis = MEGA2ANNOT;
    linkage_locus_top *LTop = Top->LocusTop;

    /* get the output file names, whether combined output for all chromosomes etc. */
    mega2annot_file_names(file_names, &combine_chromo, Top->OrigIds, Top->UniqueIds);
//SL
/*
    if (Top->LocusTop->SexLinked == 2 && combine_chromo == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/ * 50 * / Loop_Over_Chromosomes].item_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }
*/
    create_mssg(analysis);

    /* Now set the links properly, if input file was a post-makeped format */
    set_link_IDs(Top);

    /* This function sets the formatting for pedigree and individual
       ids, as well as locus Ids */
    field_widths(Top, LTop, &fwid, &pwid, NULL, &mwid);

    /* Add in code for setting other run-realted parameters that are
       specific to this option, and if these are going to remain
       constant across all chromosomes */

    /* Now write the main loop for creating chromosome-specific files:
     */
    for (i = 0; i < main_chromocnt; i++) {
        if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
            UnmappedSelected = 1;
            break;
        }
    }

    for (i = 0; i < main_chromocnt; i++) {
        if (!combine_chromo) {
            numchr = global_chromo_entries[i];
            analysis->replace_chr_number(file_names, numchr);
            if (numchr == UNKNOWN_CHROMO) {
                get_unmapped_loci(2);
            } else {
                get_loci_on_chromosome(numchr);
            }
        } else {
            // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
            analysis->replace_chr_number(file_names, 0);
            get_loci_on_chromosome(0);
            if (UnmappedSelected == 1) {
                /* append the unmapped loci */
                get_unmapped_loci(1);
            }
        }

        xlinked = ( (Top->LocusTop->SexLinked == 2 && numchr == SEX_CHROMOSOME) ) ?
                  1 : Top->LocusTop->SexLinked;

        xlinked = (Top->LocusTop->SexLinked == 2 ) ? 
                       ( (numchr == SEX_CHROMOSOME) ? 1 : combine_chromo ? 2 : 0) :
                       Top->LocusTop->SexLinked;

        omit_peds(untyped_ped_opt, Top);

        rename_mega2annot_locus(Top->LocusTop);
        annotated_ped_file(file_names[0], Top, pedfile_type);
        sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
        mssgf(err_msg);

        write_annotated_names_file(LTop, file_names[1]);
        sprintf(err_msg, "        Names file:           %s", file_names[1]);
        mssgf(err_msg);

        annotated_map_file(LTop, file_names[2]);
        sprintf(err_msg, "        Map file:             %s", file_names[2]);
        mssgf(err_msg);

        annotated_frequency_file(LTop, file_names[3]);
        sprintf(err_msg, "        Frequency file:       %s", file_names[3]);
        mssgf(err_msg);

        if (HasAff) {
            annotated_penetrance_file(file_names[4], LTop, xlinked);
            sprintf(err_msg, "        Penetrance file:      %s", file_names[4]);
            mssgf(err_msg);
        }

        draw_line();
        if (combine_chromo) break;
    }
}
