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

/*   Output files for SUP
     sup_simped.<chr>
     sup_simdata.<chr>
     sup_locus.<chr>
     sup_shell.<chr>.sh

*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc my_realloc
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line imax script_time_stamp
        write_files_ext.h:  write_linkage_locfile_inorder
*/


#define DEFAULT_sup_sim_pheno 0
#define DEFAULT_sup_seed 21341
#define  DEFAULT_sup_replicates 1
#define  DEFAULT_sup_proportion  0.00

/* prototypes */
void  create_SUP_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt);

static void SUP_file_names(char *file_names[]);

/* Just copied over from slink.c */
static void sup_slinkin_file(char *slinkin_name, int seed, int repl,
			     double proportion);

static void write_sup_shell(char *file_names[], int xlinked_specific);

static void sup_sim_opts(double *trait_positions, int *seed, int *repl,
                         double *proportion, int *sim_pheno,
                         linkage_locus_top *LTop);

static void get_sim_trait_pos(double *trait_positions,
			      linkage_locus_top *LTop);

static void write_SLINK_locus_file(linkage_locus_top *LTop,
				   char *lfl_name, int num_founders,
				   int sex_linked);

static void save_SUP_slink_peds(char *pedfl_name, linkage_ped_top *Top,
				int *num_founder_alleles, int sim_pheno,
				int xlinked);

static void SUP_file_names(char *file_names[])

{
    char reply[10], fl_stat[12];
    int user_reply = -1;
    char stem[5];
    char *outfl_name,  *slink_loutfl_name,  *sup_loutfl_name;
    char *slinkin_name,  *sup_shell_name,  *sup_output_pedfile;
    char *mega2_locus, *mega2_map, *mega2_batch;


    if (DEFAULT_OUTFILES) {
        return;
    }
    outfl_name = file_names[0];  slink_loutfl_name = file_names[1];
    sup_loutfl_name = file_names[2]; slinkin_name = file_names[3];
    sup_shell_name = file_names[4];   sup_output_pedfile = file_names[5];
    mega2_locus = file_names[6];  mega2_map = file_names[7];
    mega2_batch = file_names[8];

    while (user_reply != 0) {
        draw_line();
        print_outfile_mssg();
        printf("SUP output file name menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        printf(" 1) Pedigree file name:                  %-15s\t%s\n",
               outfl_name,
               file_status(outfl_name, fl_stat));
        printf(" 2) SLINK locus file name:               %-15s\t%s\n",
               slink_loutfl_name,
               file_status(slink_loutfl_name, fl_stat));
        printf(" 3) SLINK parameter file name:           %-15s\t%s\n",
               slinkin_name,
               file_status(slinkin_name, fl_stat));
        if (main_chromocnt <= 1) {
            strcpy(stem, "");
            printf(" 4) SUP locus file name:                 %-15s\t%s\n",
                   sup_loutfl_name,
                   file_status(sup_loutfl_name, fl_stat));
            printf(" 5) SUP output pedigree file name:       %-15s\t%s\n",
                   sup_output_pedfile,
                   file_status(sup_output_pedfile, fl_stat));
        } else {
            strcpy(stem, "stem");
            change_output_chr(sup_loutfl_name, -9);
            change_output_chr(sup_output_pedfile, -9);

            printf(" 4) SUP locus file  name stem:           %-15s\n",
                   sup_loutfl_name);
            printf(" 5) SUP output pedigree file name stem:  %-15s\n",
                   sup_output_pedfile);
        }

        printf(" 6) C-shell file  name:                  %-15s\t%s\n",
               sup_shell_name,
               file_status(sup_shell_name, fl_stat));

        if (HasLoops) {
            printf(" 7) Mega2 locus file  name:              %-15s\t%s\n",
                   mega2_locus, file_status(mega2_locus, fl_stat));

            printf(" 8) Mega2 map file  name:                %-15s\t%s\n",
                   mega2_map, file_status(mega2_map, fl_stat));

            printf(" 9) Mega2 batch file  name:              %-15s\t%s\n",
                   mega2_batch, file_status(mega2_batch, fl_stat));

            printf("Enter options 1-9 to change file names > ");
        } else {
            printf("Enter options 1-6 to change file names > ");
        }
        fcmap(stdin, "%s", reply); newline;
        user_reply=-1;
        sscanf(reply, "%d", &user_reply);
        test_modified(user_reply);

        switch (user_reply)  {
        case 0:
            break;
        case 1:
            printf("Enter new pedigree file name: > ");
            fcmap(stdin, "%s", outfl_name); newline;
            break;
        case 2:
            printf("Enter new SLINK locus file name: > ");
            fcmap(stdin, "%s", slink_loutfl_name); newline;
            break;
        case 3:
            printf("Enter new SLINK parameter file name: > ");
            fcmap(stdin, "%s", slinkin_name); newline;
            break;
        case 4:
            printf("Enter new SUP locus file name %s: > ", stem);
            fcmap(stdin, "%s", sup_loutfl_name); newline;
            break;
        case 5:
            printf("Enter new SUP output pedigree file name %s: > ", stem);
            fcmap(stdin, "%s", sup_output_pedfile); newline;
            break;
        case 6:
            printf("Enter new SUP C-shell file name: > ");
            fcmap(stdin, "%s", sup_shell_name); newline;
            break;
        case 7:
            if (HasLoops) {
                printf("Enter new Mega2 locus file name: > ");
                fcmap(stdin, "%s", mega2_locus); newline;
                break;
            }
        case 8:
            if (HasLoops) {
                printf("Enter new Mega2 map file name: > ");
                fcmap(stdin, "%s", mega2_map); newline;
                break;
            }
        case 9:
            if (HasLoops) {
                printf("Enter new Mega2 batch file name: > ");
                fcmap(stdin, "%s", sup_shell_name); newline;
                break;
            }
        default:
            warn_unknown(reply);
            break;
        }
    }
    return;

}

static void sup_sim_opts(double *trait_positions, int *seed, int *repl,
                         double *proportion, int *sim_pheno,
                         linkage_locus_top *LTop)

{

    char choice[5];
    int choice_;
    int tr, t;

    *seed = DEFAULT_sup_seed;
    *repl = DEFAULT_sup_replicates;
    *proportion = DEFAULT_sup_proportion;
    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        *sim_pheno= DEFAULT_sup_sim_pheno;
    } else {
        *sim_pheno=0;
    }

    if (DEFAULT_OPTIONS) {
        mssgf("SUP simulation parameters set to defaults");
        mssgf("as requested in the batch file.");
        if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
            mssgf("Traits will be set to position 0 on first chromosome.");
        }
        choice_ = 0;
    } else {
        choice_ = -1;
    }
    while (choice_ != 0)  {
        draw_line();
        printf("Simulation paramters menu:\n");
        printf("0) Done with this menu - please proceed.\n");
        printf(" 1) Random seed: %d\n", *seed);
        printf(" 2) Number of replicates: %d\n", *repl);
        if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
            printf(" 3) Proportion of UNLINKED families: %f\n", *proportion);
            printf(" 4) Simulated or original phenotypes: [%s]\n",
                   ((*sim_pheno)? "simulated" : "original"));
            printf(" 5) Change simulated trait location\n");
            printf("Enter selection 0 - 5 > ");
        } else {
            printf("Enter selection 0 - 2 > ");
        }

        fcmap(stdin, "%s", choice); newline;

        sscanf(choice, "%d", &choice_);
        switch (choice_)  {
        case 0:
            break;
        case 1:
            printf("Enter value of random seed ( 25000 -> 30000 ) > ");
            fcmap(stdin, "%d", seed); newline;
            if (*seed < 25000 || *seed > 30000) {
                printf("Seed value out of bounds.\n");
                *seed = DEFAULT_sup_seed;
            }
            break;
        case 2:
            printf("Enter desired number of replicates (minimum 1) > ");
            fcmap(stdin, "%d", repl); newline;
            if (*repl < 1) {
                printf("Invalid number of replicates.\n");
                *repl = DEFAULT_sup_replicates;
            }
            break;
        case 3:
            if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
                printf("Enter proportion (0  - 1) of UNLINKED families\n");
                printf("(e.g., enter '0.1' to have 90%% of your families be linked) > ");
                fcmap(stdin, "%g", proportion); newline;
                if (*proportion < 0.0 || *proportion > 1.0) {
                    printf("Invalid proportion.\n");
                    *proportion = DEFAULT_sup_proportion;
                }
                break;
            }
        case 4:
            if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
                *sim_pheno = TOGGLE(*sim_pheno); break;
            }
        case 5:
            if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
                get_sim_trait_pos(trait_positions, LTop);
                break;
            }
        default:
            warn_unknown(choice);
            break;
        }
        draw_line();
    }
    t=0;
    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        mssgf("Simulated trait positions for SUP.");
        mssgf("     Trait name    Position");
        for (tr=0; tr < num_traits; tr++) {
            if (global_trait_entries[tr] < 0) continue;
            sprintf(err_msg, "%3d) %15s  %7.4f",
                    t+1, LTop->Locus[global_trait_entries[tr]].Name,
                    trait_positions[tr]);
            mssgf(err_msg);
            t++;
        }
    } else {
        mssgf("No traits will be simulated.");
    }
    return;
}

static void get_sim_trait_pos(double *trait_positions,
			      linkage_locus_top *LTop)

{
    int m, m1, mrk_pos=-1;
    int done=0;
    char input_str[10];
    int item;

    draw_line();
    printf("Trait position menu:\n");
    while(!done) {
        printf("0) Done with this menu - please proceed.\n");
        printf("   Locus          Position\n");
        m1=0;
        for (m=0; m < num_traits; m++) {
            if (global_trait_entries[m] < 0) {
                mrk_pos=m;
                continue;
            }

            printf(" %2d)  %-13s   %f\n", m1+1,
                   LTop->Locus[global_trait_entries[m]].Name,
                   trait_positions[m]);
            m1++;
        }

        printf("Enter trait number or 0 to continue > ");
        fcmap(stdin, "%s", input_str); newline;
        sscanf(input_str, "%d", &item);
        if (mrk_pos >= 0 && (item-1) >= mrk_pos) item++;

        switch(item) {
        case 0:
            done=1;
            break;
        default:
            if (item >= 1 && item <= num_traits) {
                printf("Enter position [%7.4f] >", trait_positions[item-1]);
                (void)scanf("%lf", &trait_positions[item-1]); newline;
            } else {
                warn_unknown(input_str);
            }
        }
    }
    return;
}

static void create_locus_order(linkage_locus_top *LTop, int numchr,
			       int **locus_order, int *num_loci)

{
    int i, *loc;

    *num_loci = 0;
    *locus_order = CALLOC((size_t) NumChrLoci, int);
    loc = *locus_order;

    /* may need to add in covariates in the future */
    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Class != TRAIT) {
            *loc = ChrLoci[i]; loc++; (*num_loci)++;
        }
    }
    *locus_order = REALLOC(*locus_order, (size_t) *num_loci, int);
    return;
}

/* This is for SLINK to create the descent patterns */

static void write_SLINK_locus_file(linkage_locus_top *LTop,
				   char *loutfl_name, int num_founders,
				   int sex_linked)
{
    int tr, nloop, num_affec = num_traits;
    int *trp;
    int tmpi, tmpi2;
    FILE *filep;
    char lfl_name[2*FILENAME_LENGTH];
    double hap1_freq, hap2_freq;

    NLOOP;

    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits >= 1)) {
        trp = &(global_trait_entries[0]);
        for (tr=0; tr <= nloop; tr++) {
            if (nloop > 1 && tr==0) continue;
            sprintf(lfl_name, "%s/%s", output_paths[tr], loutfl_name);
            if ((filep = fopen(lfl_name, "w")) == NULL) {
                errorvf("Unable to open locus file '%s'\n", lfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            /* write the headers of the locus file */
            fprintf(filep, "3 0 %d %d\n%d %.3f %.3f 1\n",
                    sex_linked,
                    LTop->Program,
                    LTop->MutLocus, LTop->MutMale,
                    LTop->MutFemale);

            fprintf(filep, "1 2 3\n");
            fprintf(filep, "%d %d", (int) LTop->Locus[*trp].Type - 1,
                    LTop->Locus[*trp].AlleleCnt);
            fprintf(filep, " # %s\n", LTop->Locus[*trp].Name);
            hap1_freq = LTop->Locus[*trp].Allele[0].Frequency/(double)num_founders;
            hap2_freq = LTop->Locus[*trp].Allele[1].Frequency/(double)num_founders;

            /* don't print frequency */
            if (LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "%d\n", LTop->Pheno[*trp].Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < LTop->Pheno[*trp].Props.Affection.ClassCnt; tmpi++) {
                    if (sex_linked) {
                        for (tmpi2 = 0; tmpi2 < LTop->Pheno[*trp].Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", LTop->Pheno[*trp].Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                        fputc('\n', filep);
                        for (tmpi2 = 0; tmpi2 < LTop->Locus[*trp].AlleleCnt; tmpi2++)
                            fprintf(filep, " %.4f", LTop->Pheno[*trp].Props.Affection.Class[tmpi].MalePen[tmpi2]);
                        fputc('\n', filep);
                    } else {
                        for (tmpi2 = 0; tmpi2 < LTop->Pheno[*trp].Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", LTop->Pheno[*trp].Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                        fputc('\n', filep);
                    }
                }
            } else {
                fprintf(filep, "%d\n", LTop->Pheno[*trp].Props.Quant.ClassCnt);
                for (tmpi = LTop->Pheno[*trp].Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                    for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                        fprintf(filep," %7.6f", LTop->Pheno[*trp].Props.Quant.Mean[tmpi-1][tmpi2]);
                fprintf(filep, " << GENOTYPE MEANS\n");
                /* ASSUMES only one trait per qtl */
                fprintf(filep," %7.6f\n", LTop->Pheno[*trp].Props.Quant.Variance[0][0]);
                fprintf(filep," %7.6f\n", LTop->Pheno[*trp].Props.Quant.Multiplier);
            }

            fprintf(filep, "3 2 # Proxy\n");
            fprintf(filep, "3 %d # Descent\n", num_founders);
            /* Now print haplotype frequencies,
               perfect LD between trait and proxy,
               thus only 1-1-X and 2-2-X haplotype frequencies are non-zero.
            */
            /* First num_founders equifrequent haplotypes = major-allele-freq of trait/N */
            for (tmpi=0; tmpi < num_founders; tmpi++) {
                fprintf(filep, "%7.4f ", hap1_freq);
            }
            fprintf(filep, "\n");
            for (tmpi=0; tmpi < num_founders; tmpi++) {
                fprintf(filep, "0.0 ");
            }
            fprintf(filep, "\n");
            for (tmpi=0; tmpi < num_founders; tmpi++) {
                fprintf(filep, "0.0 ");
            }
            fprintf(filep, "\n");
            for (tmpi=0; tmpi < num_founders; tmpi++) {
                fprintf(filep, "%7.4f ", hap2_freq);
            }
            fprintf(filep, "\n");

            fprintf(filep, "%d %d\n", LTop->SexDiff,
                    LTop->Interference);

            fprintf(filep, "0.0 0.0\n");
            fprintf(filep,"1 0.10000 0.45000\n");
            fclose(filep);
            if (nloop == 1) break;
            trp++;

        }
    } else {
        sprintf(lfl_name, "%s/%s", output_paths[0], loutfl_name);
        if ((filep = fopen(lfl_name, "w")) == NULL) {
            errorvf("Unable to open locus file '%s'\n", lfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        /* write the headers of the locus file */
        fprintf(filep, "1 0 0 %d\n%d %.3f %.3f 1\n",
                LTop->Program,
                LTop->MutLocus, LTop->MutMale,
                LTop->MutFemale);

        fprintf(filep, "1\n");
        fprintf(filep, "3 %d # Descent\n", num_founders);
        for (tmpi=0; tmpi < num_founders; tmpi++) {
            fprintf(filep, "%7.4f ", 1.0/(double)num_founders);
        }
        fprintf(filep, "\n");
        fprintf(filep, "%d %d\n", LTop->SexDiff,
                LTop->Interference);

        fprintf(filep, "\n");
        fprintf(filep,"1 0.10000 0.45000\n");
        fclose(filep);
    }
    return;
}

static void save_SUP_slink_peds(char *pedfl_name, linkage_ped_top *Top,
				int *num_founder_alleles, int sim_pheno,
				int xlinked)
{
    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];
    int ped, per;
    int tr, nloop, num_affec = num_traits;
    int *trp;
    linkage_ped_rec *tpe;
    int num_ped_alleles = 0;


    NLOOP;

    if ((num_traits >= 2) || (LoopOverTrait == 1 && num_traits > 0)) {
        trp = &(global_trait_entries[0]);
    } else {
        trp=NULL;
    }
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        *num_founder_alleles = 0;
        sprintf(pedfl, "%s/%s", output_paths[tr], pedfl_name);
        if ((filep = fopen(pedfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }
        for (ped=0; ped < Top->PedCnt; ped++) {
            num_ped_alleles=0;
            for (per = 0; per < Top->Ped[ped].EntryCnt; per++) {
                tpe = &(Top->Ped[ped].Entry[per]);
                fprintf(filep, "%-4d %-4d %-4d %-4d %-4d %-4d %-4d %1d %1d ",
                        Top->Ped[ped].Num,
                        tpe->ID,
                        tpe->Father,
                        tpe->Mother,
                        tpe->First_Offspring,
                        tpe->Next_PA_Sib,
                        tpe->Next_MA_Sib,
                        tpe->Sex,
                        tpe->OrigProband);
                if (trp != NULL) {
                    /* trait locus */
                    switch(Top->LocusTop->Locus[*trp].Type) {
                    case AFFECTION:
                        if (sim_pheno) {
                            fprintf(filep, " 0");
                        } else {
                            fprintf(filep, " %1d", tpe->Pheno[*trp].Affection.Status);
                        }
                        if (Top->LocusTop->Pheno[*trp].Props.Affection.ClassCnt > 1) {
                            fprintf(filep, " %2d", tpe->Pheno[*trp].Affection.Class);
                        }
                        break;

                    case QUANT:
                        if (fabs(tpe->Pheno[*trp].Quant - MissingQuant) < EPSILON || sim_pheno) {
                            fprintf(filep, "    0.0   ");
                        } else {
                            fprintf(filep, " %10.5f", tpe->Pheno[*trp].Quant);
                        }
                        break;
                    default:
                        break;
                    }
                    /* Now the proxy marker */
                    fprintf(filep, "  0  0");
                }
                /* Then the descent marker */

                if (IS_LFOUNDER((*tpe)) && tpe->OrigProband < 2) {
                    if (tpe->Sex == 1 && xlinked) {
                        fprintf(filep, " %3d %3d", num_ped_alleles+1, num_ped_alleles+1);
                        num_ped_alleles ++;
                    } else {
                        fprintf(filep, " %3d %3d", num_ped_alleles+1, num_ped_alleles+2);
                        num_ped_alleles += 2;
                    }
                } else {
                    fprintf(filep, "   0   0");
                }
                /* here is the 'code' column ' */
                if (sim_pheno) {
                    fprintf(filep, " 1\n");
                } else {
                    fprintf(filep, " 2\n");
                }
            }
            *num_founder_alleles = imax(*num_founder_alleles, num_ped_alleles);
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;

    }
    return;
}

static void sup_slinkin_file(char *slinkin_name, int seed, int repl,
			     double proportion)

{
    FILE *fp;
    int   tr, nloop, num_affec = num_traits;
    char  slinkin[2*FILENAME_LENGTH];

    NLOOP;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        sprintf(slinkin, "%s/%s", output_paths[tr], slinkin_name);
        if ((fp = fopen(slinkin, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", slinkin);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(fp, "%d\n%d\n%d\n%f\n", seed, repl, 1, proportion);
        fclose(fp);

        if (nloop == 1) break;
    }
}

/* Used the example run.csh provided with SUP as an example.
   Creates one shell file per trait (if multiple traits), as
   well as a global shell file to run the trait-specific shell files.
   Otherwise only a single shell script.
*/
static void   write_mega2_locus(char *loutfl_name, linkage_locus_top *LTop,
				int num_founders, int sex_linked)

{

    int tr, nloop, num_affec = num_traits;
    int *trp;
    FILE *filep;
    char lfl_name[2*FILENAME_LENGTH];
    /*  int sex_linked=0; */
    int tmpi, tmpi2;

    NLOOP;

    trp = &(global_trait_entries[0]);
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(lfl_name, "%s/%s", output_paths[tr], loutfl_name);
        if ((filep = fopen(lfl_name, "w")) == NULL) {
            errorvf("Unable to open locus file '%s'\n", lfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        /* write the headers of the locus file */
        fprintf(filep, "3 0 %d %d\n%d %.3f %.3f 0\n",
                sex_linked,
                LTop->Program,
                LTop->MutLocus, LTop->MutMale,
                LTop->MutFemale);

        fprintf(filep, "1 2 3\n");
        fprintf(filep, "%d %d", (int) LTop->Locus[*trp].Type - 1,
                LTop->Locus[*trp].AlleleCnt);
        fprintf(filep, " # %s\n", LTop->Locus[*trp].Name);
        for (tmpi2 = 0; tmpi2 < LTop->Locus[*trp].AlleleCnt; tmpi2++) {
            fprintf(filep, " %.6f", LTop->Locus[*trp].Allele[tmpi2].Frequency);
        }
        fprintf(filep, "\n");

        if (LTop->Locus[*trp].Type == AFFECTION) {
            fprintf(filep, "%d\n", LTop->Pheno[*trp].Props.Affection.ClassCnt);
            for (tmpi = 0; tmpi < LTop->Pheno[*trp].Props.Affection.ClassCnt; tmpi++) {
                linkage_locus_rec *Locus = &LTop->Locus[*trp];
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
        } else {
            fprintf(filep, "%d\n", LTop->Pheno[*trp].Props.Quant.ClassCnt);
            for (tmpi = LTop->Pheno[*trp].Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                    fprintf(filep," %7.6f", LTop->Pheno[*trp].Props.Quant.Mean[tmpi-1][tmpi2]);
            fprintf(filep, " << GENOTYPE MEANS\n");
            /* ASSUMES only one trait per qtl */
            fprintf(filep," %7.6f\n", LTop->Pheno[*trp].Props.Quant.Variance[0][0]);
            fprintf(filep," %7.6f\n", LTop->Pheno[*trp].Props.Quant.Multiplier);
        }
        fprintf(filep, "3 2 # Proxy\n");
        fprintf(filep, "0.5  0.5\n");
        fprintf(filep, "3 %d # Descent\n", num_founders);
        for (tmpi=0; tmpi < num_founders; tmpi++) {
            fprintf(filep, "%7.4f", 1.0/(double)num_founders);
        }
        fprintf(filep, "\n");
        fprintf(filep, "%d %d\n", LTop->SexDiff,
                LTop->Interference);

        fprintf(filep, "0.0 0.0\n");
        fprintf(filep,"1 0.10000 0.45000\n");
        fclose(filep);
        if (nloop == 1) break;
        trp++;

    }
    return;

}

static void write_mega2_map(char *mapfl_name, char map_function)

{

    FILE *filep;
    int tr, nloop, num_affec = num_traits;
    char fl_name[2*FILENAME_LENGTH];

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], mapfl_name);

        if ((filep = fopen(fl_name, "w")) == NULL) {
            errorvf("Unable to open map file '%s'\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "Chromosome %s  Marker\n",
                ((map_function == 'k')? "Kosambi" : "Haldane"));
        fprintf(filep, "1      0.0      Proxy\n");
        fprintf(filep, "1      1.0      Descent\n");
        fclose(filep);
        if (nloop == 1) break;

    }
    return;
}

static void write_mega2_batch(char *file_names[], linkage_locus_top *LTop)

{

    FILE *filep;
    int tr, nloop, num_affec = num_traits;
    char fl_name[2*FILENAME_LENGTH];
    int *trp;

    NLOOP;

    trp = &(global_trait_entries[0]);
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], file_names[8]);

        if ((filep = fopen(fl_name, "w")) == NULL) {
            errorvf("Unable to open batch file '%s'\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        script_time_stamp(filep);
        fprintf(filep, "%s=pedfile.dat\n", Mega2BatchItems[/* 0 */ Input_Pedigree_File].keyword);
        fprintf(filep, "%s=%s\n", Mega2BatchItems[/* 1 */ Input_Locus_File].keyword, file_names[6]);
        fprintf(filep, "%s=%s\n", Mega2BatchItems[/* 2 */ Input_Map_File].keyword, file_names[7]);
        fprintf(filep, "%s=2\n", Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].keyword);
        fprintf(filep, "%s=23\n", Mega2BatchItems[/* 5 */ Analysis_Option].keyword);
        fprintf(filep, "%s=1\n", Mega2BatchItems[/* 7 */ Chromosome_Single].keyword);
        fprintf(filep, "%s=1\n", Mega2BatchItems[/* 12 */ Trait_Single].keyword);
        if (LTop->Locus[*trp].Type == QUANT) {
            fprintf(filep, "%s=%7.4f\n", Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].keyword, MissingQuant);
        }
        fclose(filep);
        if (nloop == 1) {
            break;
        }
        trp++;
    }
}

static void   write_sup_shell(char *file_names[], int xlinked_specific)

{

    int i, nloop, tr, num_affec=num_traits;
    FILE *filep;
    char mfl[2*FILENAME_LENGTH];
    char global_shell[FILENAME_LENGTH];
    char *syscmd;

    NLOOP;

    if (nloop > 1) {
        if (strcmp(output_paths[0], ".")) {
            sprintf(global_shell, "%s/%s", output_paths[0], file_names[4]);
        } else {
            strcpy(global_shell, file_names[4]);
        }
    }

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) {
            /*  create  the top level shell script */
            if ((filep = fopen(global_shell, "w")) == NULL) {
                errorvf("Unable to open %s for writing.\n", global_shell);
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

        /* Now the trait_specific shell files */
        if (nloop > 1 && tr == 0) {
            continue;
        }

        sprintf(mfl, "%s/%s", output_paths[tr], file_names[4]);
        if ((filep = fopen(mfl, "w")) == NULL) {
            errorvf("Unable to open %s for writing.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "#!/bin/csh -f\n");
        /* print some identification */
        fprintf(filep, "# C-shell file name: %s\n", file_names[4]);
        script_time_stamp(filep);
        if (!xlinked_specific || main_chromocnt > 1) {
            /* at least one autosome */
            fprintf(filep, "cp %s simped.dat\n", file_names[0]);
            fprintf(filep, "cp %s simdata.dat\n", file_names[1]);
            fprintf(filep, "cp %s slinkin.dat\n", file_names[3]);
            fprintf(filep, "slink\n");
        }

        for (i=0; i < main_chromocnt; i++) {
            if (xlinked_specific && global_chromo_entries[i] == SEX_CHROMOSOME) {
                continue;
            }
            change_output_chr(file_names[2], global_chromo_entries[i]);
            change_output_chr(file_names[5], global_chromo_entries[i]);

            if (HasLoops) {
                fprintf(filep, "%s/mega2 %s << EOF\n", mega2_path, file_names[8]);
                fprintf(filep, "2 \nprepedfile.dat\n0\nEOF\n");
                fprintf(filep, "sup -pre prepedfile.dat -d %s -o %s -nd\n",
                        file_names[2], file_names[5]);
            } else {
                fprintf(filep, "sup -p pedfile.dat -d %s -o %s -nd\n",
                        file_names[2], file_names[5]);
            }
        }

        if (xlinked_specific) {
            change_output_chr(file_names[0], SEX_CHROMOSOME);
            change_output_chr(file_names[1], SEX_CHROMOSOME);
            change_output_chr(file_names[3], SEX_CHROMOSOME);
            fprintf(filep, "cp %s simped.dat\n", file_names[0]);
            fprintf(filep, "cp %s simdata.dat\n", file_names[1]);
            fprintf(filep, "cp %s slinkin.dat\n", file_names[3]);
            fprintf(filep, "slink\n");

            change_output_chr(file_names[2], SEX_CHROMOSOME);
            change_output_chr(file_names[5], SEX_CHROMOSOME);
            if (HasLoops) {
                change_output_chr(file_names[8], SEX_CHROMOSOME);
                fprintf(filep, "%s/mega2 %s << EOF\n", mega2_path, file_names[8]);
                fprintf(filep, "2 \nprepedfile.dat\n0\nEOF\n");
                fprintf(filep, "sup -pre prepedfile.dat -d %s -o %s -nd\n",
                        file_names[2], file_names[5]);
            } else {
                fprintf(filep, "sup -p pedfile.dat -d %s -o %s -nd\n",
                        file_names[2], file_names[5]);
            }
        }

        fclose(filep);
        syscmd = CALLOC((strlen(mfl) + 15), char);
        sprintf(syscmd, "chmod +x %s", mfl);
        System(syscmd);
        free(syscmd);

        if (nloop > 1) {
            if ((filep = fopen(global_shell, "a")) == NULL) {
                errorvf("Unable to open %s for writing.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "echo Running SUP shell script %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    mfl
#endif
                );
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd %s\n", trait_paths[tr]);
            }
            fprintf(filep, "./%s\n", file_names[4]);
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd ..\n");
            }
            fclose(filep);

        }

        if (nloop == 1) break;
    }
}


#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

/* output files:
   file[0] = slink pedigree file with inheritance information
   file[1] = slink_locus file with dummy loci
   file[2] = sup locus file with the actual loci
   file[3] = slink control file
   file[4] = sup shell script
   file[5] = sup output file
   file[6:8] = Mega2-format files for breaking loops

   If there are multiple chromosomes, only file[2] and file[5] have to
   have chromosome-specific extensions, all others are created only once.

   If data has both x-linked and autosomal loci, then we need to create
   two sets of the global files (0,1,3,4,6,7,8), one for x-linked,
   and one for the other chromosomes.

*/
void  create_SUP_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt)

{

    int i, numchr, num_founders, num_loci, *locus_order;
    int sex_linked, seed, repl, sim_pheno, xlinked_specific=0;
    linkage_locus_top *LTop = Top->LocusTop;
    double proportion, *trait_positions;
    analysis_type analysis=TO_SUP;


    if (main_chromocnt > 1) {
        change_output_chr(file_names[0], 0);
        change_output_chr(file_names[1], 0);
        change_output_chr(file_names[3], 0);
        change_output_chr(file_names[4], 0);
        change_output_chr(file_names[6], 0);
        change_output_chr(file_names[7], 0);
        change_output_chr(file_names[8], 0);
    }

    SUP_file_names(file_names);

    trait_positions = CALLOC((size_t) num_traits, double);

    sup_sim_opts(trait_positions, &seed, &repl,
                 &proportion, &sim_pheno, Top->LocusTop);

    LTop->Run.slink.trait_positions = trait_positions;

    create_mssg(analysis);

    /* first create the x-linked-specific files */
    for (i=0; i < main_chromocnt; i++) {
        numchr=global_chromo_entries[i];
        sex_linked =
            (((numchr == SEX_CHROMOSOME && Top->LocusTop->SexLinked == 2) ||
              (Top->LocusTop->SexLinked == 1))? 1 : 0);

        if (sex_linked) {
            get_loci_on_chromosome(numchr);
            omit_peds(untyped_ped_opt, Top);

            change_output_chr(file_names[0], numchr);
            change_output_chr(file_names[1], numchr);
            change_output_chr(file_names[2], numchr);

            save_SUP_slink_peds(file_names[0],
                                Top, &num_founders, sim_pheno, 1);
            sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
            mssgf(err_msg);

            write_SLINK_locus_file(LTop, file_names[1], num_founders, sex_linked);
            sprintf(err_msg, "        SLINK locus file:     %s", file_names[1]);
            mssgf(err_msg);

            create_locus_order(LTop, numchr, &locus_order, &num_loci);
            write_linkage_locfile_inorder(LTop, file_names[2],
                                          num_loci, locus_order, sex_linked,
                                          analysis);

            sprintf(err_msg, "        SUP locus file:       %s", file_names[2]);
            mssgf(err_msg);

            if (HasLoops) {
                change_output_chr(file_names[6], numchr);
                change_output_chr(file_names[7], numchr);
                change_output_chr(file_names[8], numchr);

                write_mega2_locus(file_names[6], LTop, num_founders,1);
                sprintf(err_msg, "        Mega2 locus file:     %s", file_names[6]);
                mssgf(err_msg);

                write_mega2_map(file_names[7], LTop->map_distance_type);
                sprintf(err_msg, "        Mega2 map file:       %s", file_names[7]);
                mssgf(err_msg);

                write_mega2_batch(file_names, LTop);
                sprintf(err_msg, "        Mega2 batch file:     %s", file_names[8]);
                mssgf(err_msg);

                change_output_chr(file_names[6], 0);
                change_output_chr(file_names[7], 0);
                change_output_chr(file_names[8], 0);
            }

            change_output_chr(file_names[0], 0);
            change_output_chr(file_names[1], 0);
            xlinked_specific=1;
            break;
        }
    }

    if (!xlinked_specific || main_chromocnt > 1) {
        /* omit based on all autosomes */
        if (!xlinked_specific) {
            get_loci_on_chromosome(0);
        } else {
            get_loci_on_chromosome(-1);
        }
        omit_peds(untyped_ped_opt, Top);
        save_SUP_slink_peds(file_names[0],
                            Top, &num_founders, sim_pheno, 0);
        sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
        mssgf(err_msg);

        write_SLINK_locus_file(LTop, file_names[1], num_founders, 0);
        sprintf(err_msg, "        SLINK locus file:     %s", file_names[1]);
        mssgf(err_msg);
    }

    for (i=0; i < main_chromocnt; i++) {
        numchr=global_chromo_entries[i];
        change_output_chr(file_names[2], numchr);
        sex_linked =
            (((numchr == SEX_CHROMOSOME && Top->LocusTop->SexLinked == 2) ||
              (Top->LocusTop->SexLinked == 1))? 1 : 0);
        if (!sex_linked) {
            get_loci_on_chromosome(numchr);
            create_locus_order(LTop, numchr, &locus_order, &num_loci);
            write_linkage_locfile_inorder(LTop, file_names[2],
                                          num_loci, locus_order, sex_linked,
                                          analysis);

            sprintf(err_msg, "        SUP locus file:       %s", file_names[2]);
            mssgf(err_msg);

        }
    }

    sup_slinkin_file(file_names[3], seed, repl, proportion);

    sprintf(err_msg, "        SLINK control file:   %s", file_names[3]);
    mssgf(err_msg);

    if (HasLoops && (!xlinked_specific || main_chromocnt > 1)) {
        write_mega2_locus(file_names[6], LTop, num_founders, 0);
        sprintf(err_msg, "        Mega2 locus file:     %s", file_names[6]);
        mssgf(err_msg);

        write_mega2_map(file_names[7], LTop->map_distance_type);
        sprintf(err_msg, "        Mega2 map file:       %s", file_names[7]);
        mssgf(err_msg);

        write_mega2_batch(file_names, LTop);
        sprintf(err_msg, "        Mega2 batch file:     %s", file_names[8]);
        mssgf(err_msg);
    }

    write_sup_shell(file_names, xlinked_specific);
    sprintf(err_msg, "        SUP shell file:       %s", file_names[4]);
    mssgf(err_msg);

    free(trait_positions);
    return;

}
