/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

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

/***********************************************************/
/*                gen_user_input.c                         */
/*    produces input files                                 */
/*    Lazzeroni & Lange's GEN program                      */
/*    Author : Nandita Mukhopadhyay, June 1999             */
/***********************************************************/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  genoindx
  output_file_names_ext.h:  CHR_STR print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line fcopy summary_time_stamp
*/



void init_gen_outfile(char **input_files, char *outfile_name, int count_option)

{
    char count_option_str[255];
    FILE *outf;

    sprintf(count_option_str, "Counted individuals consist of: \n");
    switch (count_option) {
    case 2:
        strcat(count_option_str,
               "Genotyped founders or a randomly chosen genotyped person from each pedigree");
        break;
    case 1:
        strcat(count_option_str, "Genotyped founders only");     break;
    case 4:
        strcat(count_option_str, "All genotyped individuals");    break;
    }
/*  outf=fopen(outfile_name, "a");*/
    outf=fopen(outfile_name, "w");
    summary_time_stamp(input_files, outf, count_option_str);
    fclose(outf);
}


void gen_run_options(int *ninditer, int *nmciter, int *ntests)

{
    int parameter_option=9;

    *ninditer=0; *nmciter=5000; *ntests=1;

    if (DEFAULT_OPTIONS) {
        mssgf("GEN run parameters set to default values");
        mssgf("as requested in the batch file.");
        return;
    }

    while (parameter_option !=0) {
        draw_line();
        printf("Simulation parameters\n");
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Number of independent samples    [%4d]\n", *ninditer);
        printf(" 2) Number of Markov Chain steps     [%4d]\n", *nmciter);
        printf(" 3) Number of test runs              [%4d]\n", *ntests);
        printf("Select option 0-3 > ");
        fcmap(stdin, "%d", &parameter_option); newline;
        switch (parameter_option) {
        case 0:
            break;
        case 1:
            printf("Enter number of independent samples > ");
            fcmap(stdin, "%d", ninditer); newline;
            if (*ninditer > 1)
                *nmciter=0;
            break;
        case 2:
            printf("Enter number of Markov Chain steps > ");
            fcmap(stdin, "%d", nmciter); newline;
            if (*ninditer > 1)
                *ninditer=0;
            break;

        case 3:
            printf("Enter number of test runs: > ");
            fcmap(stdin, "%d", ntests); newline;
            break;

        default:
            printf("Unknown option: enter a number between 0 to 3\n");
            break;
        }
    }
    return;
}


char *gen_input_file(int *geno_count, int num_alleles, char *locus_name,
                     int ninditer, int nmciter, int ntests)

{
    /* creates a file named after locus_name: gen_dat.<locus_name>
       with the following format:
       Line 1:  NINDITER NMCITER NTESTS
       Line 2:  T or F (set it to false for now, for genotype-count data
       Line 3:  NROWS  NLOCI MAXALLS
       Line 4+: all1  count
       all2
       :
    */

    int nrows;
    int *all1, *all2;
    int num_genos, ii, jj, row;
    char *infile_name;
    FILE *fin;

    nrows=0;
    /* count non_zero genotypes */
    num_genos=NUMGENOS(num_alleles);
    for (ii=0; ii<num_genos; ii++) {
        if (geno_count[ii])   nrows ++;
    }
    all1= CALLOC((size_t) nrows, int);
    all2= CALLOC((size_t) nrows, int);

    /* fill in the alleles */
    row=0;
    for (ii=1; ii<num_alleles+1; ii++)
        for (jj=1; jj<ii+1; jj++)
            if (geno_count[genoindx(ii,jj)]){
                all1[row]=ii; all2[row]=jj; row++;
            }

    /* start writing the data */
    infile_name = CALLOC((strlen(locus_name)+9), char);
    strncpy(infile_name, "gen_dat.", (size_t) 9);
    strcat(infile_name, locus_name);

    fin=fopen(infile_name, "w");
    fprintf(fin, "%d, %d, %d\n", ninditer, nmciter, ntests);
    fprintf(fin, "F\n");
    fprintf(fin, "%d 1 %d\n", 2*nrows, num_alleles);
    for (ii=0; ii<nrows; ii++) {
        fprintf(fin, "%d, %d\n", all1[ii], geno_count[genoindx(all1[ii], all2[ii])]);
        fprintf(fin, "%d\n", all2[ii]);
    }
    fclose(fin);
    return infile_name;
}


void execute_gen(char *gen_executable_name, char *locus_name, const char *gen_outfile,
		 double homo_o, double homo_e, int num_genos)

{
    /* Runs the Gen program pointed to by <gen_executable_name>
       (1) renames the gen_dat.<locus> file "Gen.in"
       (2) Invokes Gen on Gen.in
       (3) copies the contents of GEN.OUT to gen_results.<chr_num>
       (4) parses gen_results.<> to form a summary table
    */

    FILE *fp;
    char *infile_name, *exec_str;

    infile_name= CALLOC(strlen(locus_name)+9, char);
    strncpy(infile_name, "gen_dat.", (size_t) 8);
    strcat(infile_name, locus_name);
    exec_str= CALLOC((size_t) 50, char);
    sprintf(exec_str, "cp %s Gen.in", infile_name);
    System(exec_str);
    free(exec_str); exec_str=NULL;
    System(gen_executable_name);
    fp=fopen(gen_outfile, "a");
    fprintf(fp, "==============\n Marker %s\n=============\n",
            locus_name);
    fprintf(fp, "Observed homozygosity: %5.4f\n", homo_o);
    fprintf(fp, "Expected homozygosity: %5.4f\n", homo_e);
    fprintf(fp, "Genotypes counted    : %d\n", num_genos);
    fclose(fp);
    if (access("GEN.OUT", F_OK) != 0) {
        /* Gen program did not run */
        printf("Gen did not run successfully, aborting further runs.\n");
        EXIT(EARLY_TERMINATION);
    }
    fcopy("GEN.OUT", gen_outfile, "a");
    return;
}

void gen_outfile_name(int numchr, char *outfile_name1, char *outfile_name2)

{
    int exit_loop, user_reply;
    char *OutFileExist1, *OutFileExist2, num[3];
    char flname[2*FILENAME_LENGTH];

    OutFileExist1 = CALLOC((size_t) FILENAME_LENGTH, char);
    OutFileExist2 = CALLOC((size_t) FILENAME_LENGTH, char);

    if (main_chromocnt > 1) {
        sprintf(outfile_name1, "gen_results.all");
        sprintf(outfile_name2, "gen_table.all");
    } else {
        CHR_STR(numchr, num);
        sprintf(outfile_name1, "gen_results.%s", num);
        sprintf(outfile_name2, "gen_table.%s", num);
    }

    if (!DEFAULT_OUTFILES) {
        exit_loop=0;
        do {
            sprintf(flname, "%s/%s", output_paths[0], outfile_name1);
            ((access(flname, F_OK) == 0)?
             strcpy(OutFileExist1, "append") :
             strcpy(OutFileExist1, "new"));

            sprintf(flname, "%s/%s", output_paths[0], outfile_name2);
            ((access(flname, F_OK) == 0)?
             strcpy(OutFileExist2, "append") :
             strcpy(OutFileExist2, "new"));

            draw_line();
            print_outfile_mssg();
            printf("GEN output file menu\n");

            printf("0) Done with this menu - please proceed\n");
            printf(" 1) GEN long-form output file name: %s\t[%s]\n",
                   outfile_name1, OutFileExist1);
            printf(" 2) GEN table-form output file name: %s\t[%s]\n",
                   outfile_name2, OutFileExist2);
            printf("Select from options 0 - 2 > ");
            fcmap(stdin, "%d", &user_reply); newline;
            test_modified(user_reply);

            if (user_reply == 0)
                exit_loop = 1;

            if (user_reply == 1) {
                printf("Enter NEW long-form output file name > ");
                fcmap(stdin, "%s", outfile_name1); newline;
            }
            if (user_reply == 2) {
                printf("Enter NEW table-form output file name > ");
                fcmap(stdin, "%s", outfile_name2); newline;
            }
            if (user_reply != 0 && user_reply != 1 && user_reply !=2)
                printf("Unknown option, enter 0, 1 or 2.\n");

        } while (!exit_loop);

        /* output type - table or long format */
        draw_line();
    }
    return;
}
