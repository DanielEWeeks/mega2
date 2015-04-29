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
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_unmapped_loci
           omit_ped_ext.h:  is_typed_at_markers omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
*/


/* Draw pedigrees with user-selected traits and marker genotypes.

   1) Creates a pedigree file, control file and shell script for each
   chromosome.

   2) If a single binary trait is selected as the primary affection
   status, then shade the individuals appropriately, this is automatic
   for loop-over-mode.

   3) Geenotyping status can be used for shading.

   4) Give the user a choice of displaying percentage genotyped as a
   quantitative trait as well.
*/



#define CRANEFOOTOUT(chr, of) sprintf(of, "crnft_chr%d", chr);

/* 1) function declarations:
   These are optional but avoid compiler warnings, if functions are
   called before they are declared, and are good for clarity. */

/* 1a) Single entry function that is called from mega2.c */

void create_CRANEFOOT_files(linkage_ped_top **LPedTop, char *file_names[],
                            int untyped_ped_opt);

/* 1b) Local function definitions  */
static void cranefoot_primary_aff(linkage_locus_top *LTop, int *primary_aff);
static void CRANEFOOT_file_names(char *file_names[], int *combine_chromo,
				 int *disp_genotyped,
				 int *disp_genotyped_perc,
				 linkage_locus_top *LTop, int has_orig);

static void save_CRANEFOOT_pedigrees(char *pedfl_name, linkage_ped_top *Top,
				     int primary_aff,
				     int disp_geno_percent);

static void write_CRANEFOOT_global_shell(char *file_names[], int combine_chromo);
static void write_CRANEFOOT_shell_script(int numchr, char *file_names[],
					 int first_time);
static void write_CRANEFOOT_control(char *file_names[], int numchr,
				    linkage_locus_top *LTop, int primary_aff,
				    int disp_geno_percent);
static void genotyped_percent(linkage_ped_top *Top);

/* primary pedigree file, contains pedigrees, phenotypes, genotypes*/

static void save_CRANEFOOT_pedigrees(char *pedfl_name, linkage_ped_top *Top,
				     int primary_aff,
				     int disp_geno_percent)
{

    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];
    int ped, per;
    int m, locus, loc;
    int tr, nloop, num_affec = num_traits;
    int aff, *trp;
    linkage_ped_rec *tpe;
    int a1, a2;

    NLOOP;

    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        trp = &(global_trait_entries[0]);
    } else {
        trp=NULL;
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        sprintf(pedfl, "%s/%s", output_paths[tr], pedfl_name);
        if ((filep = fopen(pedfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }
        /* First print the column header */
        fprintf(filep, "Pedigree\tName\tFather\tMother\tGender");

        if (LoopOverTrait == 1) {
            fprintf(filep, "\t%s", Top->LocusTop->Locus[*trp].Name);

        }
        /* Now print all the markers */
        for (locus = 0; locus < NumChrLoci; locus++)  {
            switch (Top->LocusTop->Locus[ChrLoci[locus]].Type)   {
            case AFFECTION:
                if (LoopOverTrait == 0) {
                    fprintf(filep, "\t%s", Top->LocusTop->Locus[ChrLoci[locus]].Name);
                }
                break;

            case QUANT:
                if (LoopOverTrait == 0) {
                    fprintf(filep, "\t%s", Top->LocusTop->Locus[ChrLoci[locus]].Name);
                }
                break;

            case BINARY:
            case NUMBERED:
                fprintf(filep, "\t%s", Top->LocusTop->Locus[ChrLoci[locus]].Name);
                break;
            default:
                break;
            }
        }

        if (primary_aff == -2) {
            fprintf(filep, "\tGenotyped");
        }

        if (disp_geno_percent) {
            fprintf(filep, "\tGenotyping-rate");
        }

        fprintf(filep, "\n");

        /* Now the individuals */
        for (ped=0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            for (per = 0; per < Top->Ped[ped].EntryCnt; per++) {
                tpe = &(Top->Ped[ped].Entry[per]);
                /* set the pedigree identifier to name */
                fprintf(filep, "%s", Top->Ped[ped].Name);

                /* set the individual identifier to unique ids */
                fprintf(filep, "\t%s", tpe->UniqueID);
                if (tpe->Father != 0) {
                    fprintf(filep, "\t%s",
                            Top->Ped[ped].Entry[tpe->Father-1].UniqueID);
                    fprintf(filep, "\t%s",
                            Top->Ped[ped].Entry[tpe->Mother-1].UniqueID);
                } else {
                    /* create the 0 strings */
                    fprintf(filep, "\t");
                    fprintf(filep, "\t");
                }

                if (tpe->Sex == 2) {
                    fprintf(filep, "\tF");
                } else if (tpe->Sex == 1) {
                    fprintf(filep, "\tM");
                } else {
                    fprintf(filep, "\t");
                }

                if (LoopOverTrait == 1) {
                    /* trait locus */
                    switch(Top->LocusTop->Locus[*trp].Type) {
                    case AFFECTION:
                        /* print the shape value */
                        aff = aff_status_entry(tpe->Pheno[*trp].Affection.Status,
                                               tpe->Pheno[*trp].Affection.Class,
                                               &(Top->LocusTop->Locus[*trp]));
                        if (aff == 0) {
                            /* unknown = forward hatched */
                            fprintf(filep, "\t11");
                        } else if (aff == 1) {
                            /* unaffected = unshaded */
                            fprintf(filep, "\t1");
                        } else {
                            /* affected = filled */
                            fprintf(filep, "\t99");
                        }
                        break;

                    case QUANT:
                        /* print the quantitative phenotype for display */
                        if (fabs(tpe->Pheno[*trp].Quant - MissingQuant) <= EPSILON) {
                            fprintf(filep, "\tunknown");
                        } else {
                            fprintf(filep, "\t%10f", tpe->Pheno[*trp].Quant);
                        }
                        break;
                    default:
                        break;
                    }
                }
                /* Now print the markers */
                for (m=0; m < NumChrLoci; m++) {
                    loc = ChrLoci[m];
                    switch(Top->LocusTop->Locus[loc].Type) {
                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            if (primary_aff == loc) {
                                /* print shape */
                                aff = aff_status_entry(tpe->Pheno[loc].Affection.Status,
                                                       tpe->Pheno[loc].Affection.Class,
                                                       &(Top->LocusTop->Locus[loc]));
                                if (aff == 0) {
                                    fprintf(filep, "\t11");
                                } else if (aff == 1) {
                                    fprintf(filep, "\t1");
                                } else {
                                    fprintf(filep, "\t99");
                                }
                            } else {
                                fprintf(filep, "\t%1d",
                                        aff_status_entry(tpe->Pheno[loc].Affection.Status,
                                                         tpe->Pheno[loc].Affection.Class,
                                                         &(Top->LocusTop->Locus[loc])));
                            }
                        }
                        break;
                    case QUANT:
                        if (LoopOverTrait == 0) {
                            if (fabs(tpe->Pheno[loc].Quant - MissingQuant) < EPSILON) {
                                fprintf(filep, "\tunknown");
                            } else {
                                fprintf(filep, "\t%10f", tpe->Pheno[loc].Quant);
                            }
                        }
                        break;
                    case NUMBERED:
                    case BINARY:
                        get_2alleles(tpe->Marker, loc, &a1, &a2);
			fprintf(filep, "\t%s/%s",
				format_allele(&(Top->LocusTop->Locus[loc]), a1),
				format_allele(&(Top->LocusTop->Locus[loc]), a2));
                        break;

                    default:
                        break;
                    }
                }

                /* display genotyped if set */
                if (primary_aff == -2) {
                    /* print shape for genotyped == affected */
                    /* IsTyped field for each individual is set by call to
                       omit_peds */
                    aff = tpe->IsTyped;
                    if (aff == 0) {
                        /* unaffected = unshaded */
                        fprintf(filep, "\t1");
                    } else {
                        /* affected = filled */
                        fprintf(filep, "\t99");
                    }
                }
                if (disp_geno_percent) {
                    fprintf(filep, "\t%.5f", tpe->PercentTyped);
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }
    return;

}

/* This has to come before the entry level function, it allows the
   pointer to the LPedTop to be renamed for convinience */

#ifdef Top
#undef Top
#endif

#define Top (*LPedTop)

/* 2) Template for the entry function */

void create_CRANEFOOT_files(linkage_ped_top **LPedTop, char *file_names[],
                            int untyped_ped_opt)
{
    int primary_aff, disp_geno_percent;
    int first_time=1, combine_chromo;
    int i, numchr=0;
    analysis_type analysis = CRANEFOOT;

    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    combine_chromo = main_chromocnt > 1;
    // This is the batch file item that controls whether you wish to comnine the
    // chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
    // This is a derective from the user which will override the default...
    if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
      combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;

    /* get the output file names, also use this function to set output
       related parameters, like whether an overall shell script should
       be created for all chromosome etc. */
    CRANEFOOT_file_names(file_names, &combine_chromo, &primary_aff, &disp_geno_percent,
                         Top->LocusTop, Top->OrigIds);

    /* This displays a line on the screen stating that the following
       files will be created for CRANEFOOT */

    create_mssg(analysis);

    /* This function sets the formatting for pedigree and individual
       ids, as well as locus Ids */

    /* Now write the main loop for creating chromosome-specific files:
     */

    for (i=0; i < main_chromocnt; i++) {
        if (HasMarkers) {
            numchr=global_chromo_entries[i];
            analysis->replace_chr_number(file_names, numchr);
        }

        if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
            get_unmapped_loci(0);
        } else {
            get_loci_on_chromosome(global_chromo_entries[i]);
        }

        /* This sets the IsTyped field for each individual */
        omit_peds(untyped_ped_opt, Top);
        if (disp_geno_percent) {
            genotyped_percent(Top);
        }
        save_CRANEFOOT_pedigrees(file_names[0], Top, primary_aff,
                                 disp_geno_percent);
        sprintf(err_msg, "    CRANEFOOT pedigree file:           %s", file_names[0]);
        mssgf(err_msg);

        write_CRANEFOOT_control(file_names, numchr, Top->LocusTop, primary_aff,
                                disp_geno_percent);

        sprintf(err_msg, "    CRANEFOOT control file:            %s", file_names[1]);
        mssgf(err_msg);

        write_CRANEFOOT_shell_script(numchr, file_names, first_time);

        if (!combine_chromo) {
            sprintf(err_msg, "    C-shell script:                    %s", file_names[2]);
            mssgf(err_msg);
        }

        if (combine_chromo) {
            /* for the next chromosome */
            first_time = 0;
        }
    }

    if (combine_chromo) {
        sprintf(err_msg, "    C-shell script:                    %s", file_names[2]);
        mssgf(err_msg);
    }

    if (LoopOverTrait == 1 && num_traits > 1) {
        write_CRANEFOOT_global_shell(file_names, combine_chromo);
        change_output_chr(file_names[2], 0);
        sprintf(err_msg, "    Combined shell script:           %s", file_names[2]);
        mssgf(err_msg);

        sprintf(err_msg,
                "Run %s in %s for all traits and chromosomes",
                file_names[2],
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




/* If there are multiple traits, then create a top level file for all
   traits in the output folder, global decides whether the
   trait-specific shell scripts contain calls to all chromosomes,
   else these are chromosome-specific.
*/
static void write_CRANEFOOT_global_shell(char *file_names[], int combine_chromo)
{
    char global_shell[2*FILENAME_LENGTH];
    int i, tr;
    FILE *filep;
    char *syscmd;

    sprintf(global_shell, "%s/%s", output_paths[0], file_names[2]);
    change_output_chr(global_shell, 0);

    if (combine_chromo) {
        change_output_chr(file_names[2], 0);
    }

    if ((filep = fopen(global_shell, "w")) == NULL) {
        errorvf("Unable to open %s for writing.\n", global_shell);
        EXIT(FILE_WRITE_ERROR);
    }

    for (tr=1; tr <= num_traits; tr++) {
        fprintf(filep, "cd %s\n", trait_paths[tr]);
        if (combine_chromo) {
            fprintf(filep, "./%s\n", file_names[2]);
        } else {
            for (i=0; i < main_chromocnt; i++ ) {
                change_output_chr(file_names[2], global_chromo_entries[i]);
                fprintf(filep, "./%s\n", file_names[2]);
            }
        }
        fprintf(filep, "cd ..\n");
    }

    fclose(filep);
    syscmd = CALLOC((strlen(global_shell) + 10), char);
    sprintf(syscmd, "chmod +x %s", global_shell);
    System(syscmd);
    free(syscmd);
}

/* global decides whether the shell script should call CRANEFOOT
   for all chromosomes, first_time = whether this is the first chromosome.
   Since one call to CRANEFOOT creates a single PS file, we cannot
   combine chromosomes. If there are no markers selected, numchr = 0,
   and output files will have the suffix "all".
*/

static void write_CRANEFOOT_shell_script(int numchr, char *file_names[],
					 int first_time)
{
    int  nloop, tr,/* *trp, */ num_affec=num_traits;
    FILE *filep;
    char shfl[2*FILENAME_LENGTH];
    char *syscmd;
    char cranefootoutfile[100];

    NLOOP;

    CRANEFOOTOUT(numchr, cranefootoutfile);

/*
    if (num_traits == 0) {
        trp=NULL;
    } else {
        trp= &(global_trait_entries[0]);
    }
*/
    for (tr=0; tr <= nloop; tr++) {

        if (nloop > 1 && tr == 0) {
            continue;
        }

        /* trait specific shell script (if applicable) */
        sprintf(shfl, "%s/%s", output_paths[tr], file_names[2]);

        if (first_time) {
            /* Create a new shell script file */
            if ((filep = fopen(shfl, "w")) == NULL) {
                errorvf("Unable to open %s for writing.\n", shfl);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "#!/bin/csh -f\n");
            /* print some identification */
            fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    shfl
#endif
                );
            script_time_stamp(filep);
        } else {
            /* append to the existing shell script file */
            if ((filep = fopen(shfl, "a")) == NULL) {
                errorvf("Unable to open %s for writing.\n", shfl);
                EXIT(FILE_WRITE_ERROR);
            }
        }

        if (numchr > 0) {
            fprintf(filep, "echo Running CRANEFOOT on chromosome %d\n", numchr);
        }
        fprintf(filep, " if (-e %s.ps) then\n", cranefootoutfile);
        fprintf(filep, "    /bin/rm %s.ps\n", cranefootoutfile);
        fprintf(filep, " endif\n");
        fprintf(filep, " if (-e %s.topology.txt) then\n", cranefootoutfile);
        fprintf(filep, "    /bin/rm %s.topology.txt\n", cranefootoutfile);
        fprintf(filep, " endif\n");

        fprintf(filep, "\n");

        fprintf(filep, " cranefoot %s\n ", file_names[1]);

        fclose(filep);

        if (first_time) {
            syscmd = CALLOC((strlen(shfl) + 10), char);
            sprintf(syscmd, "chmod +x %s", shfl);
            System(syscmd);
            free(syscmd);
        }
        if (nloop == 1) break;

    }

    return;

}

static void write_CRANEFOOT_control(char *file_names[], int numchr,
				    linkage_locus_top *LTop, int primary_aff,
				    int disp_geno_percent)

{

    int nloop, tr, trp, num_affec=num_traits;
    char cfl[2*FILENAME_LENGTH];
    int loc;

    FILE *fp;

    char cranefootoutfile[100];

    CRANEFOOTOUT(numchr, cranefootoutfile);

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        if (LoopOverTrait == 1 && tr > 0) {
            trp = tr - 1;
        } else {
            trp = tr;
        }

        sprintf(cfl, "%s/%s", output_paths[tr], file_names[1]);
        if ((fp = fopen(cfl, "w")) == NULL) {
            errorvf("Could not open control file %s for writing.\n", file_names[1]);
            EXIT(FILE_WRITE_ERROR);
        }

        /* First write down the pedigree variables */
        fprintf(fp, "PedigreeFile %s\n", file_names[0]);
        fprintf(fp, "PedigreeName %s\n", cranefootoutfile);
        fprintf(fp, "SubgraphVariable Pedigree\n");
        fprintf(fp, "NameVariable Name\n");
        fprintf(fp, "FatherVariable Father\n");
        fprintf(fp, "MotherVariable Mother\n");
        fprintf(fp, "GenderVariable Gender\n");
        fprintf(fp, "TextVariable Name\n");
        /* Now write the text variables which are phenotypes and markers */
        if (LoopOverTrait == 1) {
            switch(LTop->Locus[trp].Type) {
            case AFFECTION:
                fprintf(fp, "PatternVariable %s\n", LTop->Locus[trp].Name);
                break;
            case QUANT:
                fprintf(fp, "TextVariable %s \n", LTop->Locus[trp].Name);
                break;
            default:
                break;
            }
        }

        /* Now print the markers */
        for (loc = 0; loc < NumChrLoci; loc++) {
            switch(LTop->Locus[ChrLoci[loc]].Type) {
            case AFFECTION:
                if (LoopOverTrait == 0) {
                    if (primary_aff == ChrLoci[loc]) {
                        fprintf(fp, "PatternVariable %s\n", LTop->Locus[ChrLoci[loc]].Name);
                        fprintf(fp, "PatternInfo Unknown 11\n");
                        fprintf(fp, "PatternInfo Affected 99\n");
                    } else {
                        fprintf(fp, "TextVariable %s\n", LTop->Locus[ChrLoci[loc]].Name);
                    }
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0) {
                    fprintf(fp, "TextVariable %s\n", LTop->Locus[ChrLoci[loc]].Name);
                }
                break;
            case BINARY:
            case NUMBERED:
                fprintf(fp, "TextVariable %s\n", LTop->Locus[ChrLoci[loc]].Name);
                break;

            default:
                break;
            }
        }
        /* Now the primary aff and percent */
        if (primary_aff == -2) {
            fprintf(fp, "PatternVariable Genotyped\n");
            fprintf(fp, "PatternInfo Genotyped 99\n");
        }

        if (disp_geno_percent) {
            fprintf(fp, "TextVariable Genotyping-rate\n");
        }

        fclose(fp);
        if (nloop == 1) break;

    }
}

/**
   This is a function that implements the file-names menu for CRANEFOOT,
   The user is not allowed to change the choice of
   individual IDs, however, output pedigree id choice can be changed
 */
static void CRANEFOOT_file_names(char *file_names[], int *combine_chromo,
				 int *primary_aff,
				 int *disp_geno_percent,
				 linkage_locus_top *LTop,
				 int has_orig)
{
    int i, choice=-1;
    char cchoice[10], stem[5], fl_stat[12];
    analysis_type anal = CRANEFOOT;
    int iped, iprime, igdisp;
    int iglob, ipedf, icont, ishell;

    strcpy(stem, "   ");

    if (*combine_chromo) {
        // returns only the <file_name> removing <extension> and <rest> if they exist...
        anal->replace_chr_number(file_names, -9);
        change_output_chr(file_names[2], 0);
    } else {
        if (!HasMarkers) {
            // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
            anal->replace_chr_number(file_names, 0);
            change_output_chr(file_names[2], 0);
        }
    }

    *primary_aff = -2;

    for (i=0; i < num_traits; i++) {
        SKIP_TRI(i)
            if (LTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
                /* set it directly to the trait locus number, not index of trait list */
                *primary_aff = global_trait_entries[i];
                break;
            }
    }
    /* Display genotyped percentage by default */
    *disp_geno_percent =  (HasMarkers ? 1 : 0);

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        /* set the primary affection locus to the first affection status locus */
        return;
    }

    /* output file name menu */
    iped = iprime = igdisp =
        iglob = ipedf = icont = ishell = -1;

    while (choice != 0) {
        print_outfile_mssg();
        printf("Output file names menu:\n");
        printf("0) Done with this menu - please proceed\n");
        i=1;
        if (main_chromocnt > 1) {
            printf(" %d) Run all chromosomes using one shell script? %s\n",
                   i, yorn[*combine_chromo]);
            iglob=i; i++;

            printf(" %d) Pedigree file name stem                     %-15s\n",
                   i, file_names[0]);
            ipedf = i; i++;

            printf(" %d) Control file name stem                      %-15s\n",
                   i, file_names[1]);
            icont = i; i++;

            printf(" %d) Shell file name %s                          %-15s  %s\n",
                   i, stem, file_names[2],
                   ((*combine_chromo == 1)?
                    file_status(file_names[2], fl_stat) : ""));
            ishell = i;
        } else {
            printf(" %d) Pedigree file name                          %-15s  %s\n",
                   i, file_names[0], file_status(file_names[0], fl_stat));
            ipedf = i; i++;

            printf(" %d) Control file name                           %-15s  %s\n",
                   i, file_names[1], file_status(file_names[1], fl_stat));
            icont = i; i++;

            printf(" %d) Shell file name                             %-15s  %s\n",
                   i, file_names[2], file_status(file_names[2], fl_stat));
            ishell = i;

        }

        if (LoopOverTrait == 0) {
            i++;
            printf(" %d) Affection locus to use for shading          ",
                   i);
            if (*primary_aff == -1) {
                printf("none\n");
            } else if (*primary_aff == -2) {
                printf("genotyped-status\n");
            } else {
                printf("%s\n", LTop->Locus[*primary_aff].Name);
            }

            iprime = i;
        }

        if (HasMarkers) {
            i++;
            printf(" %d) Display genotyping-rate \n", i);
            printf("            for each individual                 [%s]\n",
                   yorn[*disp_geno_percent]);
            igdisp = i;
        }

        if (has_orig) {
            i++;
            pedigree_id_item(i, anal, OrigIds[1], 55, 2, 0);
            iped=i;
        }

        printf("Enter options 0-%d > ", i);

        fcmap(stdin, "%s", cchoice); printf("\n");
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch(choice) {
        case 0:
            break;
        default:
            if (choice == iglob) {
                *combine_chromo = TOGGLE(*combine_chromo);
                if (*combine_chromo) {
                    change_output_chr(file_names[2], 0);
                    strcpy(stem, "    ");
                    strcat(file_names[2], ".sh");
                } else {
                    change_output_chr(file_names[2], -9);
                    strcpy(stem, "stem");
                }
            } else if (choice == ipedf) {
                printf("New pedigree file name %s > ", stem);
                fcmap(stdin, "%s", file_names[0]);
                printf("\n");
            } else if (choice == icont) {
                printf("New control file name %s > ", stem);
                fcmap(stdin, "%s", file_names[1]);
                printf("\n");
            } else if (choice == ishell) {
                printf("Name of C-Shell file name %s > ", stem);
                fcmap(stdin, "%s", file_names[2]);
                printf("\n");
            } else if (choice == iped) {
                OrigIds[1] = pedigree_id_item(0, anal, OrigIds[1], 35, 1,
                                              has_orig);
            } else if (choice == iprime) {
                cranefoot_primary_aff(LTop, primary_aff);
            } else if (choice == igdisp) {
                *disp_geno_percent = TOGGLE(*disp_geno_percent);
            } else {
                warn_unknown(cchoice);
            }
            break;
        }
    }

    pedigree_id_item(0, anal, OrigIds[1], 0, 3, has_orig);

    if (main_chromocnt > 1 && !combine_chromo) {
        anal->replace_chr_number(file_names, global_chromo_entries[0]);
    }
}

static void cranefoot_primary_aff(linkage_locus_top *LTop, int *primary_aff)

{
    int a, choice, loc;
    char cchoice[10];
    int *affs = CALLOC((size_t)num_affection+2, int);

    a=0;
    for (loc=0; loc < num_traits; loc++) {
        SKIP_TRI(loc)
            if (LTop->Locus[global_trait_entries[loc]].Type == AFFECTION) {
                affs[a] = global_trait_entries[loc];
                a++;
            }
    }

    if (HasMarkers) {
        affs[a] = -2;
    }
    affs[a+1] = -1;

    choice = -1;
    while(choice) {
        printf("Select affection status trait to use for shading:\n");
        printf("0) Done with this menu - please proceed.\n");
        for (loc=0; loc < a; loc++) {
            if (affs[loc] == *primary_aff) {
                printf("*%d) %s\n", loc+1, LTop->Locus[affs[loc]].Name);
            } else {
                printf(" %d) %s\n", loc+1, LTop->Locus[affs[loc]].Name);
            }
        }
        if (HasMarkers) {
            affs[a] = -2;
            if (*primary_aff == -2) {
                printf("*");
            } else {
                printf(" ");
            }
            printf("%d) GENOTYPED\n", a+1);
        }
        affs[a+1] = -1;
        if (*primary_aff == -1) {
            printf("*");
        } else {
            printf(" ");
        }
        printf("%d) NONE (do not shade)\n", a+2);

        printf("Select from 1-%d > ", a+2);
        fcmap(stdin, "%s", cchoice);
        sscanf(cchoice, "%d", &choice);
        if (choice == 0) {
            ;
        } else if (choice >= 1 && choice <= a+2) {
            *primary_aff = affs[choice-1];
        } else {
            warn_unknown(cchoice);
        }
    }
    free(affs);
    return;
}

static void genotyped_percent(linkage_ped_top *LPedTreeTop)

{
    entry_type Entry;
    int i, j;
    double percent_typed;

    for (i = 0; i < LPedTreeTop->PedCnt; i++) {
        for (j=0; j < LPedTreeTop->Ped[i].EntryCnt; j++) {
            Entry.LEntry = &(LPedTreeTop->Ped[i].Entry[j]);
            is_typed_at_markers(Entry, LPedTreeTop->pedfile_type,
                                LPedTreeTop->LocusTop, &percent_typed);
            LPedTreeTop->Ped[i].Entry[j].PercentTyped = percent_typed;
        }
    }
    return;

}
