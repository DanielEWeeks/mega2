/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

#include "common.h"
#include "typedefs.h"
#include "analysis.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file
            linkage_ext.h:  get_loci_on_chromosome switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line log_line randomnum script_time_stamp
        write_files_ext.h:  save_linkage_peds
      write_ghfiles_ext.h:  write_gh_locus_file
*/


#define DEFAULT_SIMULATE_num_cycles 200
#define SIMUL_PEDOUT(chrstr, of) sprintf(of, "pedfile.%s", chrstr);
#define SIMUL_SEEDOUT(chrstr, of) sprintf(of, "problem.%s", chrstr);

static void  simulate_file_names(char *file_names[], analysis_type *analysis)

{
    int choice=-1;
    char fl_stat[12], stem[5];

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        log_line(mssgf);
        return;
    }
    while(choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("SIMULATE-format output file name selection:\n");
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt > 1) {
            // returns only the <file_name> removing <extension> and <rest> if they exist...
            (*analysis)->replace_chr_number(file_names, -9);
            strcpy(stem, "stem");
            printf(" 1) Pedigree file name stem:       %s\n", file_names[0]);
            printf(" 2) Locus file name stem:          %s\n", file_names[1]);
            printf(" 3) Parameter file name stem:      %s\n", file_names[2]);
            printf(" 4) Shell script file name stem:   %s\n", file_names[3]);
        } else {
            strcpy(stem, "");
            printf(" 1) Pedigree file name:            %-15s\t%s\n",file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:               %-15s\t%s\n",file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 3) Parameter file name:           %-15s\t%s\n",file_names[2],
                   file_status(file_names[2], fl_stat));
            printf(" 4) Shell script file name:        %-15s\t%s\n", file_names[3],
                   file_status(file_names[3], fl_stat));
        }

        printf("Enter option 0-4 > ");
        fcmap(stdin, "%d", &choice); newline;
        test_modified(choice);

        switch (choice) {
        case 1:
            printf("Enter new  pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]); printf("\n"); newline;
            break;
        case 2:
            printf("Enter new  locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]); printf("\n"); newline;
            break;
        case 3:
            printf("Enter new  parameter file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]); printf("\n");  newline;
            break;

        case 4:
            printf("Enter new  shell script file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]); printf("\n");  newline;
            break;

        case 0:
            break;
        default:
            printf("Unknown option %d\n", choice);
            break;
        }
    }
    if (main_chromocnt > 1) {
        /* put back the .sh extension */
        change_output_chr(file_names[3], global_chromo_entries[0]);
        strcat(file_names[3], ".sh");
    }
    return;
}

static int simulate_num_cycles(void)

{
    int choice=-1;
    int num_cycles= DEFAULT_SIMULATE_num_cycles;

    if (DEFAULT_OPTIONS) {
        mssgf("Number of replicates set to default");
        mssgf("as requested in the batch file.");
        return num_cycles;
    }
    draw_line();
    printf("Number of simulation replicates selection: \n");
    while(choice != 0) {
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Number of replicates:   %d\n", num_cycles);
        printf("Enter option 0 or 1 > ");
        fcmap(stdin, "%d", &choice);
        printf("\n");
        if (choice == 1) {
            printf("Enter number of replicates > ");
            fcmap(stdin, "%d", &num_cycles); printf("\n");
        }
    }
    return num_cycles;
}

static void write_simulate_params(char *slinkin, int num_cycles)

{
    int tr, nloop, num_affec= num_traits;
    char sp[2*FILENAME_LENGTH];
    FILE  *slinkin_p;

    NLOOP;

    for(tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;
        sprintf(sp, "%s/%s", output_paths[tr], slinkin);
        if ((slinkin_p = fopen(sp, "w")) == NULL) {
            errorvf("Unable to open file %s for writing\n", sp);
            EXIT(FILE_WRITE_ERROR);
        }
        randomnum();
        fprintf(slinkin_p,"%d %d %d ", seed1, seed2, seed3);
        fprintf(slinkin_p,"%d\n", num_cycles);
        fclose(slinkin_p);
        if (nloop == 1) break;
    }
    return;
}


static void write_simulate_script(char *file_names[], int numchr,
				  int first_time, int write_global)

{

    char shfl[2*FILENAME_LENGTH], *syscmd;
    char global_shell[FILENAME_LENGTH];
    FILE *shellfp;
    int nloop, tr, num_affec=num_traits;
    char simul_pedout[20];
    char simul_seedout[20];
    char chr_str[3];

    NLOOP;

    CHR_STR(numchr, chr_str);
    SIMUL_PEDOUT(chr_str, simul_pedout );
    SIMUL_SEEDOUT(chr_str, simul_seedout);

    for (tr=0; tr <=nloop; tr++) {
        if (tr==0 && nloop > 1) continue;

        sprintf(shfl, "%s/%s", output_paths[tr], file_names[3]);
        syscmd=CALLOC(strlen(shfl)+15, char);
        delete_file(shfl);
        if ((shellfp = fopen(shfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", shfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(shellfp, "#!/bin/csh -f\n");
        script_time_stamp(shellfp);
        fprintf(shellfp, "unalias rm\n");
        fprintf(shellfp, "rm -f simdata.dat\n");
        fprintf(shellfp, "rm -f simped.dat\n");
        fprintf(shellfp, "rm -f problem.dat\n");
        fprintf(shellfp, "rm -f pedfile.dat\n");
        fprintf(shellfp, "cp %s simdata.dat\n", file_names[1]);
        fprintf(shellfp, "cp %s simped.dat\n", file_names[0]);
        fprintf(shellfp, "cp %s problem.dat\n", file_names[2]);
        fprintf(shellfp, "simulate\n");
        fprintf(shellfp, "cp pedfile.dat %s\n", simul_pedout);
        fprintf(shellfp, "cp problem.dat %s\n", simul_seedout);
        fclose(shellfp);
        sprintf(syscmd, "chmod +x %s", shfl);
        System(syscmd);
        if (write_global) {
            strcpy(global_shell, file_names[3]);
            change_output_chr(global_shell, 0);
            sprintf(shfl, "%s/%s", output_paths[tr], global_shell);

            if (first_time) {
                delete_file(shfl);
                if ((shellfp = fopen(shfl, "w")) == NULL) {
                    errorvf("Unable to open file '%s' for writing\n", shfl);
                    EXIT(FILE_WRITE_ERROR);
                }
                script_time_stamp(shellfp);
                sprintf(syscmd, "chmod +x %s", shfl);
                System(syscmd);
            } else {
                if ((shellfp = fopen(shfl, "a")) == NULL) {
                    errorvf("Unable to open file '%s' for writing\n", shfl);
                    EXIT(FILE_WRITE_ERROR);
                }
            }
            fprintf(shellfp, "./%s\n", file_names[3]);
            fclose(shellfp);
        }

        free(syscmd);
        if (nloop == 1) break;
    }
}


/*        change LPedTreeTop to *Top */
#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*Top)

void create_SIMULATE_format_files(linkage_ped_top **Top,
				  ped_top *PedTreeTop,
				  analysis_type *analysis,
				  file_format *infl_type,
				  file_format *outfl_type, int *numchr,
				  char *mapfl_name, char *file_names[],
				  int untyped_ped_opt)
{

    int             i, num_cycles;
    int             xlinked, first_time;

    num_cycles=simulate_num_cycles();
    simulate_file_names(file_names, analysis);

    /* make a copy of the pedigree tree */
/*  LPedTreeTop_copy = new_lpedtop(); */
/*  copy_linkage_ped_top(LPedTreeTop, LPedTreeTop_copy, 1); */
    first_time=1;
    for (i=0; i < main_chromocnt; i++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[i];
            (*analysis)->replace_chr_number(file_names, *numchr);
        }
        get_loci_on_chromosome(*numchr);
        xlinked =
            (((*numchr == SEX_CHROMOSOME && LPedTreeTop->LocusTop->SexLinked == 2) ||
              (LPedTreeTop->LocusTop->SexLinked == 1))? 1 : 0);

        /* omit untyped  peds */
        omit_peds(untyped_ped_opt, LPedTreeTop);
        save_linkage_peds(file_names[0],LPedTreeTop,*analysis, 0, NULL);
        write_simulate_params(file_names[2], num_cycles);
        write_gh_locus_file(file_names[1], LPedTreeTop->LocusTop, *analysis,
                            xlinked);
        if (main_chromocnt > 1) {
            write_simulate_script(file_names, *numchr, first_time, 1);
        } else {
            write_simulate_script(file_names, *numchr, first_time, 0);
        }
        if (i==0) create_mssg(*analysis);
        sprintf(err_msg, "    Pedigree file:     %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "    Locus file:        %s", file_names[1]);
        mssgf(err_msg);
        sprintf(err_msg, "    Parameter file:    %s", file_names[2]);
        mssgf(err_msg);
        sprintf(err_msg, "    Shell script file: %s", file_names[3]);
        mssgf(err_msg);
        first_time=0;
    }

    if (main_chromocnt > 1) {
        change_output_chr(file_names[3], 0);
        sprintf(err_msg, "    Combined shell script file for all chromosomes: %s",
                file_names[3]);
        mssgf(err_msg);
    }
    /*  free_all_including_lpedtop(&LPedTreeTop_copy); */
    return;

}   /* End of create_SIMULATE_format_files  */
