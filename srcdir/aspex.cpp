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

#include "batch_input_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "scripts_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
        batch_input_ext.h:  batchf
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  check_map_positions delete_file haldane_to_kosambi kosambi_to_haldane sub_prog_name write_time
        grow_string_ext.h:  grow
            linkage_ext.h:  get_loci_on_chromosome new_lpedtop
           makenucs_ext.h:  create_nuked_fids makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
            scripts_ext.h:  create_sib_phase_graph_scripts
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
*/


#ifndef ASPEX_H
#define ASPEX_H
#define TRUE_FALSE(v) ((v)? "true" : "false")
#endif

typedef struct tcl_opts_ {

    int mapopt , boolean[14]; /*set by user */

    int FidWidth, PidWidth, AlleleWidth, LocWidth; /* automatically set */
    int dloc, nloc, sex_linked, *allele_cnt;  /* automatically set */

    double max_step, error_freq, exclusion_level; /*set by user */
    double *dist;
    double **allele_freq;  /* set by program except 1st & last element of dist */

    int *locus_ids;
    char **locus_names, **dlocus_name; /* automatically set */

}tcl_opts_type;


/*----------------------------------------------------------------------------
  Function to free the opt data structure
  ---------------------------------------------------------------------------- */

static void free_opt(tcl_opts_type opt)
{
    int i;

    for(i=0; i<opt.nloc; i++) {
        free(opt.allele_freq[i]);
        free(opt.locus_names[i]);
    }
    for(i=0; i < num_traits; i++) {
        if (opt.dlocus_name[i] != NULL) {
            free(opt.dlocus_name[i]);
        }
    }
    free(opt.dlocus_name);
    free(opt.locus_ids);
    free(opt.dist);
    free(opt.allele_freq);
    free(opt.locus_names);
    free(opt.allele_cnt);
    return;
}


/*----------------------------------------------------------------------------
  Function to get 1st set of options (nulear? and global_cshell?)
  ---------------------------------------------------------------------------- */

static void AspexMainOptions(int *nuked, int *global_cshell,
			     int single_chr)

{
    int option = -1, exit_loop=0;
    char opt[10];
    *nuked=1;
    *global_cshell=(single_chr? 0 : 1);

    if (DEFAULT_OPTIONS) {
        mssgf("Aspex run parameters set to default values");
        mssgf("as requested in the batch file.");
        return;
    }

    while (!exit_loop) {
        draw_line();
        printf("Main analysis options: \n");
        printf("0) Done with this menu - please proceed.\n");
        printf(" 1) Convert all pedigrees to nuclear families [%s]\n",
               ((*nuked)? "yes" : "no"));
        if (single_chr==1) {
            printf("Enter 1 to toggle setting, 0 to proceed > ");
            *global_cshell=0;
        } else {
            printf(" 2) Single C-shell script for all chromosomes [%s]\n",
                   ((*global_cshell)? "yes" : "no"));
            printf("Enter 1 or 2 to toggle settings > ");
        }
        fcmap(stdin, "%s", opt); newline;
        option = -1;
        sscanf(opt, "%d", &option);

        switch (option) {
        case 1:
            *nuked = ((*nuked)? 0 : 1);
            break;
        case 2:
            if (single_chr==1) {
                printf("Unknown option %d\n", option);
            } else {
                *global_cshell = ((*global_cshell)? 0 : 1);
            }
            break;
        case 0:
            exit_loop=1;
            break;

        default:
            printf("Unknown option %s\n", opt);
            break;
        }
    }
    return;

}

/*----------------------------------------------------------------------------
  Function to get risk ratio options from user
  ---------------------------------------------------------------------------- */

static void AspexRiskRatio(int AspexPrg,
                           int *SelectRiskOpt, double *risk, double *InitRisk,
                           double *RiskIncrement, double *FinalRisk)

{

    int choice=-1;
    char choice_[10], Marker[2], risk_str[50];

    Marker[0]='*'; Marker[1]=' ';
    /* Default values */
    *SelectRiskOpt=1;
    *risk = 2.0;
    *InitRisk =   1.5;
    *RiskIncrement =   0.5;
    *FinalRisk =   2.5;

    if (AspexPrg==2 || AspexPrg ==4 || DEFAULT_OPTIONS) return;

    while (choice != 0) {
        printf("Aspex risk ratio selection:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        printf("%c 1) Pre-defined single risk ratio   [%4.2f]\n", Marker[0], *risk);
        printf("%c 2) Pre-defined grid of risk ratios [start=%4.2f increment=%4.2f end=%4.2f]\n",
               Marker[1], *InitRisk, *RiskIncrement, *FinalRisk);

        printf("Enter 1 or 2  to change risk ratios, 0 to proceed > ");
        fcmap(stdin, "%s", choice_); newline;
        choice = -1;
        sscanf(choice_, "%d", &choice);

        switch (choice) {

        case 1:
            *SelectRiskOpt=1; Marker[0]='*'; Marker[1]=' ';
            printf("Enter new risk ratio > ");
            fcmap(stdin, "%g", risk); printf("\n");
            break;

        case 2:
            *SelectRiskOpt=2; Marker[1]='*'; Marker[0]=' ';
            printf("Enter new start, increment and final risk ratios \n(enter \"none\" to use default values) > ");
            fflush(stdout);
            (void)fgets(risk_str, 49, stdin); printf("\n");
            if (strcmp(risk_str, "none")) {
                sscanf(risk_str, "%lg %lg %lg", InitRisk, RiskIncrement, FinalRisk);
            }
            break;

        case 0:
            break;

        default:
            printf("Unknown option %s \n", choice_);
            break;
        }
    }
    return;

}

/*----------------------------------------------------------------------------
  Function to set/modify output file names
  ---------------------------------------------------------------------------- */

static void AspexOutfileNames(int global_cshell, int AspexOption,
			      char *file_names[],
			      int numchr, int nuked,
			      int *orig_id_opt, int has_orig, int has_uniq)

{

    int select_reply=-1, exit_loop=0;
    char select_[10], stem[5], chr_str[3];
    char fl_stat[12];
    char shell_file_names[4][10]= {
        "sib_ibd",
        "sib_tdt",
        "sib_phase",
        "sib_map"
    };

    if (nuked) {
        orig_id_opt[0]=1;
        orig_id_opt[1]=4;
    }

    /* get file_name stems */
    switch (AspexOption) {
    case 1:
    case 2:
    case 3:
    case 4:
        strcpy(file_names[3], shell_file_names[AspexOption-1]); break;
    default:
        errorf("Unknown aspex option.");
        EXIT(INPUT_DATA_ERROR);
        break;
    }

    if (global_cshell)
        strcat(file_names[3], ".all.sh");
    if (main_chromocnt == 1) {
        CHR_STR(numchr, chr_str);
        grow(file_names[3], ".%s.sh", chr_str);
    }

    if (DEFAULT_OUTFILES) {
        if (global_cshell || main_chromocnt == 1) {
            /* initialize the combined shell script */
            FILE *fp = fopen(file_names[3], "w");
            if (fp == NULL) {
                errorvf("Could not open %s for writing\n", file_names[3]);
                EXIT(FILE_WRITE_ERROR);
            }
            fclose(fp);
        }
        if (main_chromocnt > 1 && !global_cshell) {
            strcat(file_names[3], ".01.sh");
        }
        return;
    }
    while (!exit_loop) {
        draw_line();
        /*    if (main_chromocnt>1)
              need multiple asp_in and asp_dat files, therefore set only the name stem */
        print_outfile_mssg();
        printf("ASPEX file name menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt<=1) {
            strcpy(stem, "");
            printf(" 1) Locus file name:           %-15s\t%s\n", file_names[4],
                   file_status(file_names[4], fl_stat));
            printf(" 2) Pedigree datafile name:    %-15s\t%s\n", file_names[5],
                   file_status(file_names[5], fl_stat));
            printf(" 3) C-shell script file name:  %-15s\t%s\n", file_names[3],
                   file_status(file_names[3], fl_stat));
        } else {
            strcpy(stem, "stem");
            printf(" 1) Locus file name %s:           %-15s\n", stem,
                   file_names[4]);
            printf(" 2) Pedigree datafile name %s:    %-15s\n",  stem,
                   file_names[5]);
            if (global_cshell) {
                printf(" 3) C-shell script file name:       %-15s\t%s\n",
                       file_names[3], file_status(file_names[3], fl_stat));
            } else {
                printf(" 3) C-shell script file name %s:  %-15s\n", stem,
                       file_names[3]);
            }
        }

        individual_id_item(4, TO_ASPEX, orig_id_opt[0], 43, 2, 0, 0);
        pedigree_id_item(5, TO_ASPEX, orig_id_opt[1], 43, 2, 0);
        printf("Select options 0-5 to enter new values > ");

        fcmap(stdin, "%s", select_); newline;
        select_reply = -1;
        sscanf(select_, "%d", &select_reply);
        test_modified(select_reply);

        switch (select_reply) {
        case 0:
            exit_loop = 1; break;

        case 1:
            printf("\nEnter NEW TCL-format locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[4]); printf("\n");
            break;

        case 2:
            printf("\nEnter NEW ASPEX-format pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[5]); printf("\n");
            break;

        case 3:
            if (global_cshell || main_chromocnt == 1) {
                printf("\nEnter NEW ASPEX C-Shell script file name > ");
            } else {
                printf("\nEnter NEW ASPEX C-Shell script file name stem > ");
            }
            fcmap(stdin, "%s", file_names[3]); printf("\n");
            break;

        case 4:
            orig_id_opt[0] = individual_id_item(0, TO_ASPEX, orig_id_opt[0], 35, 1,
                                                has_orig, has_uniq);
            break;

        case 5:
            if (nuked) {
                orig_id_opt[1] = pedigree_id_item(0, TO_NUKE,
                                                  orig_id_opt[1], 35, 1,
                                                  has_orig);
            } else {
                orig_id_opt[1] = pedigree_id_item(0, TO_ASPEX,
                                                  orig_id_opt[1], 35, 1,
                                                  has_orig);
            }
            break;

        default:
            printf("Unknown option %s \n", select_);
            break;
        }
    }

    if (global_cshell || main_chromocnt == 1) {
        /* initialize the combined shell script */
        FILE *fp = fopen(file_names[3], "w");
        fclose(fp);
    }
    if (main_chromocnt > 1 && !global_cshell) {
        strcat(file_names[3], ".01.sh");
    }
    return;

}

/*----------------------------------------------------------------------------
  Function to initialize shell script
  --------------------------------------------------------------------------- */

static void AspexCshellFile(char *cshellname, char *aspin_name,
			    int numchr, int first_time,
			    char *prog, int prog_num,
			    char *aspdat_name,
			    int risk_opt, int num_ratios,
			    double *ratios)

{
    int i, tr, nloop, num_affec=num_traits;
    FILE *shell_fp;
    char fl_name[2*FILENAME_LENGTH], mode[2], syscmd[2*FILENAME_LENGTH], chr_str[3];
    int xlinked =  ((numchr == SEX_CHROMOSOME)? 1 : 0);

    CHR_STR(numchr, chr_str);
    if (first_time)  strcpy(mode, "w");
    else strcpy(mode, "a");

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], cshellname);

        if (first_time){
            delete_file(fl_name);
        }

        if ((shell_fp = fopen(fl_name, mode)) == NULL)    {
            errorvf("Unable to open aspex shell file %s.\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        if (first_time) {
            fprintf(shell_fp, "#!/bin/csh -f\n");
            script_time_stamp(shell_fp);
        }
        /* write a single line to execute the appropriate program,
           don't know what the sib_ibd option does */

        if (risk_opt==1)
            fprintf(shell_fp, "%s -f %s %s\n", prog, aspin_name, aspdat_name);

        /* xmgr graph scripts for multiple risk ratios */
        if (risk_opt == 2)  {
            /* set the right chromosome number for the first run */
            change_output_chr(aspin_name, numchr);
            if (prog_num == 3) {
                create_sib_phase_graph_scripts(aspin_name, aspdat_name, 0.00,
                                               shell_fp, chr_str, 0, xlinked);
            } else {
                fprintf(shell_fp, "%s -f %s %s \n",
                        prog, aspin_name, aspdat_name);
            }
            /* take out the chromosome number for the subsequent runs, and insert
               the risk ratios (this is also done inside sib_pohase_graph_scripts) */
            change_output_chr(aspin_name, -9);
            for (i=0; i< num_ratios; i++) {
                if (prog_num  == 3) {
                    create_sib_phase_graph_scripts(aspin_name, aspdat_name, ratios[i],
                                                   shell_fp, chr_str, i+1, xlinked);
                } else {
                    fprintf(shell_fp, "%s -f %s.%4.2f.%s %s \n",
                            prog, aspin_name, ratios[i], chr_str, aspdat_name);
                }
            }
        }
        fclose(shell_fp);
        sprintf(syscmd, "chmod +x %s", fl_name);
        System(syscmd);

        if (nloop == 1) break;
    }
    return;
}

/*----------------------------------------------------------------------------
  Function to get Aspex program name from user
  --------------------------------------------------------------------------- */

static void AspexProgramOption(int *SelectAspexPrg)

{
    int exit_loop=0;
    char select[10], sub_prog[20];

    if (batchANALYSIS) {
        *SelectAspexPrg=TO_ASPEX->_suboption;
        return;
    }

    *SelectAspexPrg = -1;
    while (!exit_loop) {
        draw_line();
        printf("ASPEX program selection menu:\n");
        draw_line();
        printf("1) sib_ibd\n");
        printf("2) sib_tdt\n");
        printf("3) sib_phase\n");
        printf("4) sib_map\n");
        printf("Enter option [1-4] > ");
        fcmap(stdin, "%s",  select); newline;
        *SelectAspexPrg = -1;
        sscanf(select, "%d", SelectAspexPrg);
        if (*SelectAspexPrg >= 1  && *SelectAspexPrg <= 4)
            exit_loop=1;
        else
            printf("Unknown option %s\n", select);

    }
    
    TO_ASPEX->_suboption = *SelectAspexPrg;
    TO_ASPEX->sub_prog_name(*SelectAspexPrg, &sub_prog[0]);
    strcpy(Mega2BatchItems[/* 6 */ Analysis_Sub_Option].value.name, sub_prog);
    batchf(Analysis_Sub_Option);
}

/*----------------------------------------------------------------------------
  Function to get TCL options from user
  ---------------------------------------------------------------------------- */

static void AspexTCLOptions(int AspexPrg, char *AspexPrgName,
                            tcl_opts_type *opt, int risk_opt)

{


    int  select=-1, exit_loop=0;
    int item=0;
    char select_[10];

    /* set TCL parameters to default values : these values have been changed
       according to the ASPEX documentation */

    /* boolean values */
    opt->boolean[0]= 1; /* discard_partial */
    opt->boolean[1]= 1; /* linear model */
    if (risk_opt==1) opt->boolean[2]= 1; /*  most_likely */
    opt->boolean[3]= 0; /* fixed_step_ */
    opt->boolean[4]= 1; /* no_Dv_ */
    opt->boolean[5]= 0; /* fixed_freq_ */
    opt->boolean[6]= 1; /* truncate_sharing_ */
    opt->boolean[7]= 0; /* count_once_   */
    opt->boolean[8]= 0; /* first_pair_ */
    opt->boolean[9]= 0; /*count_unaffected_ */
    opt->boolean[10]= 0; /* count_discordant_ */
    opt->boolean[11]= 0; /* sex_split */
    opt->boolean[12]= 0; /* show_pairs */
    opt->boolean[13]= 0; /* show_freq */

    opt->mapopt=0;

    /* double values */

    opt->dist[0] = 10.0;
    opt->dist[1] = 10.0;
    opt->max_step=  0.01;
    opt->error_freq=  0.00;
    opt->exclusion_level=  0.0;

    if (DEFAULT_OPTIONS) {
        return;
    }

    while (!exit_loop) {
        item=0;
        printf("ASPEX TCL default parameters for %s:\n", AspexPrgName);
        printf("(Genetic distance increases from the left to the right on the chromosome.)\n");
        draw_line();

        printf("\n%d) Done with this menu - please proceed.\n", item); item++;
        printf(" %d) discard_partial:        [%s]\n", item, TRUE_FALSE(opt->boolean[0])); item++;
        printf(" %d) linear_model:           [%s]\n", item, TRUE_FALSE(opt->boolean[1])); item++;
        if (risk_opt==1) {
            printf(" %d) most_likely:            [%s] \n", item, TRUE_FALSE(opt->boolean[2])); item++;
        }
        printf(" %d) fixed_step:             [%s] \n", item, TRUE_FALSE(opt->boolean[3])); item++;
        printf(" %d) no_Dv:                  [%s]\n",  item, TRUE_FALSE(opt->boolean[4])); item++;
        printf(" %d) fixed_freq:             [%s]\n",  item, TRUE_FALSE(opt->boolean[5])); item++;
        printf(" %d) truncate_sharing:       [%s]\n",  item, TRUE_FALSE(opt->boolean[6])); item++;
        printf(" %d) count_once:             [%s]\n",  item, TRUE_FALSE(opt->boolean[7])); item++;
        printf(" %d) first_pair:             [%s]\n",  item, TRUE_FALSE(opt->boolean[8])); item++;
        printf(" %d) count_unaffected:      [%s]\n",  item, TRUE_FALSE(opt->boolean[9])); item++;
        printf(" %d) count_discordant:      [%s]\n",  item, TRUE_FALSE(opt->boolean[10])) ;item++;
        printf(" %d) sex_split:             [%s]\n",  item, TRUE_FALSE(opt->boolean[11])) ;item++;
        printf(" %d) show_pairs:            [%s]\n",  item, TRUE_FALSE(opt->boolean[12])) ;item++;
        if (AspexPrg == 3) {
            printf(" %d) show_freqs:            [%s]\n",  item, TRUE_FALSE(opt->boolean[13])) ;item++;
        }
        printf(" %d) mapping:               [%s]\n", item, ((opt->mapopt==0)? "haldane" : "kosambi"));
        item++;
        printf(" %d) max_step:              [%7.4f cM]\n", item, 100 * opt->max_step); item++;
        printf(" %d) error_freq             [%5.4f]\n", item, opt->error_freq); item++;
        printf(" %d) exclusion_level        [%5.4f]\n", item, opt->exclusion_level); item++;
        printf(" %d) Distance to left end\n", item);  item++;
        printf("     of chromosome          [%7.4f cM]\n", opt->dist[0]);
        printf(" %d) Distance to right end\n", item);
        printf("     of chromosome          [%7.4f cM]\n", opt->dist[1]);
        draw_line();
        printf("(Selecting \"true/false\" options will toggle their values)\n");
        printf("Select from options 0-%d > ", item);

        fcmap(stdin, "%s", select_); newline;
        select = -1;
        sscanf(select_, "%d", &select);

        if (risk_opt==2 && select >= 3)  select++;
        /* no most_likely option, so user's input 3 becomes 4, 4 becomes 5 and so on */
        if (AspexPrg != 3 && select >=14) select++;
        /* no show_freqs option, so user's input 14 becomes 15, and so on */
        switch (select) {
        case 0:
            exit_loop=1; break;
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
            opt->boolean[select-1]=TOGGLE(opt->boolean[select-1]);
            break;
        case 15:
            opt->mapopt=TOGGLE(opt->mapopt);
            break;
        case 16:
            printf("Enter max_step (cM): > ");
            fcmap(stdin, "%g", &(opt->max_step));
            printf("\n");
            break;
        case 17:
            printf("Enter error_freq: (between 0.0 and 1.0) > ");
            fcmap(stdin, "%g", &(opt->error_freq));
            printf("\n");
            if (opt->error_freq < 0.0 || opt->error_freq > 1.0) {
                printf("Warning: error_freq %5.4f out of bounds, using default 0.0\n", opt->error_freq);
                opt->error_freq= 0.0;
            }
            break;
        case 18:
            printf("Enter exclusion_level: (should be less than or equal to 0) > ");
            fcmap(stdin, "%g", &(opt->exclusion_level));
            printf("\n");
            if (opt->exclusion_level > 0.0) {
                printf("Warning: Positive exclusion level %5.4f, using default 0.0\n", opt->exclusion_level);
                opt->exclusion_level=  0.0;
            }
            break;
        case 19:
            printf("Enter distance (cM) to the LEFT end of chromosome > ");
            fcmap(stdin, "%g", &(opt->dist[0])); newline;
            break;
        case 20:
            printf("Enter distance (cM) to the RIGHT end of chromosome > ");
            fcmap(stdin, "%g", &(opt->dist[1])); newline;
            break;
        default:
            printf("Unknown option %s\n", select_);
            break;
        }
    }
    return;

}

/*----------------------------------------------------------------------------
  Function to fill in locus_related stuff into the opt data structure for
  each cromosome.

  NM - aspex does not allow multiple traits to be combined, so don't have
  to check for marker item in global_trait_entries.
  ---------------------------------------------------------------------------- */

static int assign_distances(tcl_opts_type *opt, int numchr,
			    linkage_locus_top *LTop1,
			    int SelectionNumber, char AspexProgram[4][10])

{
    int *missing_index, *ok_index, i, j, nloc, miss=0;
    double left_dist, right_dist, *distances;
    char user_reply[3];


    miss=0; opt->dloc=0;
    nloc=0;
    missing_index=CALLOC((size_t)NumChrLoci, int);
    ok_index=CALLOC((size_t)NumChrLoci,int);

    /* find out if all marker loci have positive genetic distances */
    opt->dlocus_name = CALLOC((size_t)num_traits, char *);
    for (i=0; i<num_traits; i++) {
        opt->dlocus_name[i] = NULL;
    }

    for (i=0,j=0; i<num_traits; i++) {
        if (global_trait_entries[i] < 0) continue;
        opt->dlocus_name[j] = CALLOC((size_t)MAX_NAMELEN, char);
        sprintf(opt->dlocus_name[j++], "%s", LTop1->Locus[global_trait_entries[i]].Name);
    }
    for (i = 0; i < NumChrLoci; i++)  {
        if (LTop1->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop1->Locus[ChrLoci[i]].Type == BINARY) {
            if (LTop1->Marker[ChrLoci[i]].pos_avg < 0.0) {
                missing_index[miss] = ChrLoci[i];
                miss++;
            } else {
                ok_index[nloc] = ChrLoci[i];
                nloc++;
            }
        }
    }

    if (miss>0) {
        if (SelectionNumber == 1 || SelectionNumber == 3) {
            errorf(
                "The following loci don't have genetic distances:\n");

            sprintf(err_msg,
                    "The ASPEX program '%s' needs ALL inter-loci distances",
                    AspexProgram[SelectionNumber - 1]);
            errorf(err_msg);

        } else {
            warnf(
                "The following loci don't have genetic distances:");

        }
        for (i = 0; i < miss; i++) {
            sprintf(err_msg, "%s ", LTop1->Locus[missing_index[i]].Name);
            errorf(err_msg);
        }
        printf("Do you wish to abort mega2 ? (y/n)[y] > ");
        fcmap(stdin, "%s", user_reply); newline;

        if (user_reply[0] != 'n'  && user_reply[0] != 'N') {
            errorvf("Mega2 run aborted.\n");
            EXIT(EARLY_TERMINATION);
        }
        if (SelectionNumber == 1 || SelectionNumber == 3) {
            sprintf(err_msg,
                    "Excluding chromosome %d from analysis.", numchr);
            mssgf(err_msg);
            return 0;
        }
    }

    opt->nloc=nloc;
    left_dist=opt->dist[0]; right_dist=opt->dist[1];
    free(opt->dist);
    opt->dist= CALLOC((size_t)(opt->nloc)+1, double);
    opt->dist[0]=left_dist;
    opt->dist[opt->nloc]=right_dist;

    opt->allele_cnt= CALLOC((size_t)opt->nloc, int);
    opt->allele_freq=(double **) CALLOC_PTR((size_t)opt->nloc, double *);
    opt->locus_names=(char **) CALLOC_PTR((size_t)opt->nloc, char *);
    distances= CALLOC((size_t)opt->nloc, double);

    for (nloc = 0; nloc < opt->nloc; nloc++)  {
        i=ok_index[nloc];
        opt->locus_names[nloc]=
            CALLOC(strlen(LTop1->Locus[i].Name)+1, char);
        strcpy(opt->locus_names[nloc], LTop1->Locus[i].Name);
        opt->allele_cnt[nloc] = (int) LTop1->Locus[i].AlleleCnt;
        opt->allele_freq[nloc]=
            CALLOC((size_t) LTop1->Locus[i].AlleleCnt, double);
        distances[nloc] = LTop1->Marker[i].pos_avg;
        for (j = 0; j < opt->allele_cnt[nloc]; j++)
            opt->allele_freq[nloc][j] = LTop1->Locus[i].Allele[j].Frequency;
    }

    opt->locus_ids= CALLOC((size_t)opt->nloc, int);

    for (i = 0; i < opt->nloc; i++)
        opt->locus_ids[i]=ok_index[i];

    for (i=1; i < opt->nloc; i++) {
        opt->dist[i] = distances[i] - distances[i-1];
        if (opt->mapopt == 0 && LTop1->map_distance_type == 'k') {
            opt->dist[i] = kosambi_to_haldane(opt->dist[i]);
        }
        if (opt->mapopt == 1 && LTop1->map_distance_type == 'h') {
            opt->dist[i] = haldane_to_kosambi(opt->dist[i]);
        }
    }

    free(ok_index);
    free(missing_index); free(distances);

    return 1;
}



/*----------------------------------------------------------------------------
  Function to write TCL options into asp_in.##
  ---------------------------------------------------------------------------- */

static void write_contents(char *aspin_name, char *aspex_program,
			   double risk,  tcl_opts_type opt)

{
    int tr, nloop, *trp, num_affec = num_traits;
    int i, j;
    char fl_name[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;
    trp = &(global_trait_entries[0]);

    for (tr=0; tr <=nloop; tr++) {
        if (tr==0 && nloop > 1) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], aspin_name);
        if ((fp = fopen(fl_name, "w")) == NULL) {
            errorvf("Could not open %s for writing.\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(fp, "#This ASPEX TCL file is: '%s'\n", aspin_name);

        fprintf(fp, "#Locus file :    %s\n",
#ifdef HIDEPATH
                NOPATH
#else
                mega2_input_files[0]
#endif
            );
        fprintf(fp, "#Pedigree file : %s\n",
#ifdef HIDEPATH
                NOPATH
#else
                mega2_input_files[1]
#endif
            );

        fprintf(fp, "#Map file :      %s\n",
#ifdef HIDEPATH
                NOPATH
#else
                mega2_input_files[2]
#endif
            );
        fprintf(fp, "#Date     :      "); write_time(fp);
        fprintf(fp, "#This is a TCL parameter file for ASPEX '%s' program.\n", aspex_program);
        fprintf(fp, "#the set linkage_format must be set true to use these ASPEX output files.\n");
        if (num_affec >= 1) {
            if (tr > 0) {
                fprintf(fp, "#disease locus %s\n", opt.dlocus_name[tr-1]);
            } else {
                fprintf(fp, "#disease locus %s\n", opt.dlocus_name[0]);
            }
        } else {
            fprintf(fp, "#disease locus none\n");
        }
        fprintf(fp, "set nloc %d\n", opt.nloc);
        fprintf(fp, "set loc {\n");
        for (i = 0; i < opt.nloc; i++)
            fprintf(fp, " \"%s\"\n", opt.locus_names[i]);
        fprintf(fp, " }\n");
        fprintf(fp, "set blank 0\n");
        fprintf(fp, "set fid_width %d\n", opt.FidWidth);
        fprintf(fp, "set pid_width %d\n", opt.PidWidth);
        fprintf(fp, "set loc_width %d\n", opt.LocWidth);
        fprintf(fp, "set allele_width %d\n", opt.AlleleWidth);
        fprintf(fp, "set discard_partial %s\n", TRUE_FALSE(opt.boolean[0]));
        fprintf(fp, "set sex_linked %s \n", TRUE_FALSE(opt.sex_linked));
        fprintf(fp, "set dist {\n");
        if (opt.dist[0] < 0.0) { fprintf(fp, " 0.0999\n"); }
        else { fprintf(fp, " %5.4f\n", opt.dist[0]/100.0);  }

        for (i = 1; i < opt.nloc; i++) {
            if (opt.dist[i] < 0.0) {fprintf(fp, " 0.0999\n"); }
            else { fprintf(fp, " %5.4f\n", opt.dist[i]/100.0);}
        }
        fprintf(fp, " %5.4f\n", opt.dist[opt.nloc]/100.0);
        fprintf(fp, " }\n");
        fprintf(fp, "set mapping %s\n", ((opt.mapopt==0)? "haldane": "kosambi"));
        fprintf(fp, "set risk %4.2f\n", risk);
        if (opt.boolean[1]) {
            /* additive model */
            fprintf(fp, "#Identity by descent probabilities : Additive model\n");
            if (!opt.sex_linked)
                fprintf(fp, "set z \"[expr 0.25/$risk] 0.50 [expr 0.50-0.25/$risk]\"\n");
            else
                fprintf(fp, "set z \"[expr 0.50/$risk] [expr (1.0 - (0.50/$risk))]\"\n");
        } else {
            fprintf(fp, "#Identity by descent probabilities : Multiplicative model\n");
            if (opt.sex_linked) {
                printf("Don't know how to compute z for sex-linked in multiplicative model\n");
                printf("Using linear model\n");
                fprintf(fp, "set z \"[expr 0.50/$risk] [expr (1.0 - (0.50/$risk))]\"\n");
            } else {
                fprintf(fp,
                        "set z \"[expr 0.25/$risk], [expr 1/sqrt($risk)-0.5/$risk], [expr 1+0.25/$risk-1/sqrt($risk)]\"\n");
            }
        }

        fprintf(fp, "set max_step %5.4f\n", opt.max_step);
        fprintf(fp, "set fix_step %s\n", TRUE_FALSE(opt.boolean[3]));
        fprintf(fp, "set most_likely %s \n", TRUE_FALSE(opt.boolean[2]));
        fprintf(fp, "set linear_model %s \n", TRUE_FALSE(opt.boolean[1]));
        fprintf(fp, "set no_Dv %s \n", TRUE_FALSE(opt.boolean[4]));
        fprintf(fp, "set fixed_freq %s \n", TRUE_FALSE(opt.boolean[5]));
        fprintf(fp, "set truncate_sharing %s\n", TRUE_FALSE(opt.boolean[6]));
        fprintf(fp, "set exclusion_level %5.4f\n", opt.exclusion_level);
        fprintf(fp, "set count_once %s\n", TRUE_FALSE(opt.boolean[7]));
        fprintf(fp, "set first_pair %s\n", TRUE_FALSE(opt.boolean[8]));
        fprintf(fp, "set count_unaffected %s\n", TRUE_FALSE(opt.boolean[9]));
        fprintf(fp, "set count_discordant %s\n", TRUE_FALSE(opt.boolean[10]));
        fprintf(fp, "set sex_split %s\n", TRUE_FALSE(opt.boolean[11]));
        fprintf(fp, "set show_pairs %s\n", TRUE_FALSE(opt.boolean[12]));
        if (strcmp(aspex_program, "sib_phase")==0)
            fprintf(fp, "set show_freqs %s\n", TRUE_FALSE(opt.boolean[13]));
        fprintf(fp, "set error_freq %f\n", opt.error_freq);

        if (opt.boolean[5])  {
            for (i = 0; i < opt.nloc; i++) {
                fprintf(fp, "freq \"%s\" {\n", opt.locus_names[i]);
                for (j = 0; j < opt.allele_cnt[i]; j++)
                    fprintf(fp, " \"%d\" %f \n", j + 1, opt.allele_freq[i][j]);
                fprintf(fp, " }\n");
            }
        }
        fclose(fp);
        if (nloop == 1) break;
        trp++;

    }

    sprintf(err_msg, "   TCL parameters file:       %s", aspin_name);
    mssgf(err_msg);
    return;
}


static void AspexPedFile(char *aspdat_name, tcl_opts_type opt,
			 linkage_ped_top *TTop)

{

    int i, j;
    int  k, lcount, tr, *trp, nloop, num_affec=num_traits;
    char fl_name[2*FILENAME_LENGTH], fformat[8], pformat[8];
    FILE *fp;
    linkage_ped_rec *tpe;

    NLOOP;

    trp = &(global_trait_entries[0]);

    for (tr=0; tr <=nloop; tr++) {
        if (tr==0 && nloop > 1) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], aspdat_name);

        if ((fp = fopen(fl_name, "w")) == NULL) {
            errorvf("Could not open %s for writing.\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        for (k = 0; k < opt.nloc; k++)
            fprintf(fp, "%s ", opt.locus_names[k]);
        fprintf(fp, "\n");

        create_formats(opt.FidWidth, opt.PidWidth, fformat, pformat);

        for (j = 0; j < TTop->PedCnt; j++)  {
            for (i = 0; i < TTop->Ped[j].EntryCnt; i++) {
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(fp, fformat, TTop->Ped[j].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(fp, fformat, (int) j+1);
                } else {
                    fprintf(fp, fformat, TTop->Ped[j].Num);
                }

                tpe = &(TTop->Ped[j].Entry[i]);

                if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                    fprintf(fp, pformat, tpe->UniqueID);
                    if (tpe->Father>0) {
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Father-1].UniqueID);
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Mother-1].UniqueID);
                    } else {
                        fprintf(fp, pformat, "0");
                        fprintf(fp, pformat, "0");
                    }
                }

                else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    fprintf(fp, pformat, tpe->OrigID);
                    if (tpe->Father>0) {
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Father-1].OrigID);
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Mother-1].OrigID);
                    } else {
                        fprintf(fp, pformat, "0");
                        fprintf(fp, pformat, "0");
                    }
                } else {
                    fprintf(fp, pformat, tpe->ID);
                    if (tpe->Father > 0) {
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Father-1].ID);
                        fprintf(fp, pformat, TTop->Ped[j].Entry[tpe->Mother-1].ID);
                    } else {
                        fprintf(fp, pformat, 0);
                        fprintf(fp, pformat, 0);
                    }
                }
                fprintf(fp, "%-4d", tpe->Sex);

                if (num_affec > 0) {
                    fprintf(fp, " %-4d",
                            aff_status_entry(tpe->Pheno[*trp].Affection.Status, tpe->Pheno[*trp].Affection.Class,
                                             &(TTop->LocusTop->Locus[*trp])));
                }
                /* 	printf("%d nloc\n", chromo_loci_count[numchr]); */
                for (lcount = 0; lcount < opt.nloc; lcount++)  {
                    int allele1, allele2;
                    linkage_locus_rec *locus = &TTop->LocusTop->Locus[*trp];
                    get_2alleles(tpe->Marker, opt.locus_ids[lcount], &allele1, &allele2);
                    fprintf(fp, "  %2s", format_allele(locus, allele1));
                    if (TTop->LocusTop->Marker[opt.locus_ids[lcount]].chromosome == SEX_CHROMOSOME &&
                        tpe->Sex == 1)
                        /* = 1 = male , 2= female */
                        fprintf(fp, " Y");
                    else {
                        fprintf(fp, "  %2s", format_allele(locus, allele2));
                    }
                    /*      max of 10 digits/characters per allele identifier... */
                }
                fprintf(fp, "\n");
            }
        }
        trp++;
        fclose(fp);
        if (nloop == 1) break;
    }

    sprintf(err_msg, "   Pedigree file:             %s", aspdat_name);
    mssgf(err_msg);
    return;
}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)
/*----------------------------------------------------------------------------
  Function to create the Aspex files: asp_in.##, asp_dat.##
  and <aspex_prog_name>.sh
  ---------------------------------------------------------------------------- */

void  create_aspex_files(linkage_ped_top **LPedTop,
			 int numchr, char *file_names[],
			 int untyped_ped_opt, linkage_ped_top **NukeTop)
{

    char            AspexProgram[4][10] = {"sib_ibd",
                                           "sib_tdt",
                                           "sib_phase",
                                           "sib_map"};
    char		  chr_str[3], aspin_name1[FILENAME_LENGTH];

    int             SelectRiskOpt, SelectionNumber;
    int             num_ratios;

    double           this_risk, left, right, risk, *risk_ratios = NULL, InitRisk,
        RiskIncrement, FinalRisk;  //Note: risk_rations is only used if riskopt == 2

    int             nuked, global_cshell;
    int             i, count, first_time;

    linkage_ped_top *NewTop = NULL;

    linkage_ped_top *TTop = NULL; //silly compiler
    linkage_locus_top *LTop1;
    tcl_opts_type opt;

    /* This function (create_aspex_files) is called after
       the loci have been reordered */

    /* Step 1: Main options
       convert to nuclear families?
       have a global cshell file? */
    AspexMainOptions(&nuked, &global_cshell, main_chromocnt);

    /* Create a copy of Top, NewTop which will be nuked */

    /* Step 2: select ASPEX program */
    draw_line();
    AspexProgramOption(&SelectionNumber);
#ifndef HIDESTATUS
    mssgvf("Selected ASPEX program %s.\n", AspexProgram[SelectionNumber-1]);
#endif
    /* Step 3:  check for presence and position of affection status locus
       These are critical or not depending on the aspex program */

    if ((num_traits < 1) || (LoopOverTrait == 0 && num_traits < 2)) {
        /* cannot recover for 1st 3 aspex programs*/
        if (SelectionNumber != 4) {
            errorvf("There are NO trait loci present.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    }

    /* Step 5: read in risk ratio */
    AspexRiskRatio(SelectionNumber, &SelectRiskOpt, &risk,
                   &InitRisk, &RiskIncrement, &FinalRisk);

    /* Step 5.1: create a vector of risk ratios */
    if (SelectRiskOpt==2) {
        this_risk=InitRisk;
        num_ratios=(int) ceil((FinalRisk-InitRisk)/RiskIncrement) + 1;

        risk_ratios= CALLOC((size_t)num_ratios, double);
        num_ratios = 0;
        while (this_risk <= FinalRisk) {
            risk_ratios[num_ratios] = this_risk;
            num_ratios++;
            this_risk = InitRisk + num_ratios * RiskIncrement;
        }
    }

    /* Step 6: read in the Tcl parameters
       main_chromocnt is a global variable */

    /* allocate two elements for distances to the left and right ends of the chromosome
       the inter-loci distances will be filled in later for each chromosome */
    opt.dist= CALLOC((size_t)2, double);
    AspexTCLOptions(SelectionNumber, AspexProgram[SelectionNumber-1], &opt, SelectRiskOpt);

    /* Step 7: get output file names, the risk ratio option is necessary for
       naming files*/
    if (main_chromocnt > 1) {
        change_output_chr(file_names[4], -9);
        change_output_chr(file_names[5], -9);
    }

    AspexOutfileNames(global_cshell, SelectionNumber, file_names, numchr,
                      nuked, &(OrigIds[0]), Top->OrigIds, Top->UniqueIds);


    /* Step 8: Now ready to write the scripts and output files */
    if (nuked) {
        NewTop = new_lpedtop();
        makenucs1(Top, NULL, NewTop);
        create_nuked_fids(NewTop);
    }

    first_time=1;
    for (count = 0; count < main_chromocnt; count++) {
        draw_line();
        if (main_chromocnt > 1) {
            numchr = global_chromo_entries[count];
            change_output_chr(file_names[4], numchr);
            change_output_chr(file_names[5], numchr);
            if (!global_cshell) {
                change_output_chr(file_names[3], numchr);
            }
        }

        get_loci_on_chromosome(numchr);
        /* omit untyped peds */
        omit_peds(untyped_ped_opt, Top);

        if (nuked) {
            TTop = NewTop;
        } else if (!nuked) {
            TTop=Top;
        }

        LTop1=TTop->LocusTop;
        /* Step 9.0: decide on field widths and create formats */
        field_widths(TTop, LTop1,
                     &(opt.FidWidth), &(opt.PidWidth), &(opt.AlleleWidth),
                     &(opt.LocWidth));

        CHR_STR(numchr, chr_str);
        if (first_time) {
            create_mssg(TO_ASPEX);
        }

        /* Step 9.1 : write the TCL parameter file(s) */

        /* Step 9.1.a: assign inter-marker distances exect the first and the last*/

        opt.sex_linked = ((numchr == SEX_CHROMOSOME) ? 1: 0);
        if (assign_distances(&opt, numchr, LTop1, SelectionNumber, AspexProgram) == 0)
            continue;

        /* 1. single risk ratio */
        if (SelectRiskOpt==2) {
            /* the first file should have most_likely set to true */
            opt.boolean[2] = 1;
        }

        write_contents(file_names[4], AspexProgram[SelectionNumber-1],
                       risk, opt);

        /* 2. grid of risk ratios */
        if (SelectRiskOpt==2) {
            opt.boolean[2]= 0;
            change_output_chr(file_names[4], -9);
            for(i=0; i<num_ratios; i++) {
                sprintf(aspin_name1, "%s.%4.2f.%2s", file_names[4], risk_ratios[i],
                        chr_str);
                write_contents(aspin_name1, AspexProgram[SelectionNumber-1],
                               risk_ratios[i],  opt);
            }
            /* change_output_chr(file_names[4], numchr); */
        }

        /* Step 9.2:
         * Create the aspex data file (using the linkage_format parameter set
         * on ), presence of affection data column is assumed here
         */
        AspexPedFile(file_names[5], opt, TTop);
        /* Step 9.3 :
           create the shell scripts
        */

        if (!global_cshell) first_time=1;
        AspexCshellFile(file_names[3], file_names[4], numchr, first_time,
                        AspexProgram[SelectionNumber-1], SelectionNumber,
                        file_names[5], SelectRiskOpt, num_ratios,
                        risk_ratios);

        first_time=0;

        if (!global_cshell) {
            /* close file and set it to be executable */
            sprintf(err_msg, "   C-shell script:            %s", file_names[3]);
            mssgf(err_msg);
        }
        left= opt.dist[0];
        right= opt.dist[opt.nloc];
        free_opt(opt);
        opt.dist= CALLOC((size_t)2, double);
        opt.dist[0]=left; opt.dist[1]=right;

    } /* End of for loop */

    /* Step 10: close file_handles and clear data structures
       if a global shell script has been created, change its mode */

    if (global_cshell) {
        sprintf(err_msg, "   C-shell script:            %s", file_names[3]);
        mssgf(err_msg);
    }

    draw_line();
    if (SetMarkerPosToSpecial) {
        check_map_positions(TTop, numchr, 2);
        draw_line();
    }

    free(opt.dist);
    *NukeTop = NewTop;
}

#undef Top
