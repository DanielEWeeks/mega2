/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"
#include "mrecode.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "gen_user_input_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "homozygosity_ext.h"
#include "hwe_ext.h"
#include "linkage_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "select_individuals_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_mendel7_files_ext.h"
#include "write_mfiles_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
        batch_input_ext.h:  batchf
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
     gen_user_input_ext.h:  execute_gen gen_input_file gen_outfile_name gen_run_options init_gen_outfile
      genetic_utils_ext.h:  delete_file genoindx init_file
        grow_string_ext.h:  grow
       homozygosity_ext.h:  expected_homozygosity observed_homozygosity
                hwe_ext.h:  hwe_test
            linkage_ext.h:  get_loci_on_chromosome get_unmapped_loci
  output_file_names_ext.h:  CHR_STR create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  field_widths
 select_individuals_ext.h:  get_count_option select_individuals
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line fcopy log_line script_time_stamp summary_time_stamp
write_mendel7_files_ext.h:  csv_write_mendel_locus_file csv_write_mendel_map
       write_mfiles_ext.h:  create_mendel_formats write_mendel_locus_file write_mendel_map
*/


char *get_hwe_option(char prog_names[5][8]);

static void hwe_copyright(FILE *fp)

{

    fprintf(fp, "These results are produced by the HWE program by \n");
    fprintf(fp, "Sun-Wei Guo and Elizabeth Thompson.  If you publish \n");
    fprintf(fp, "any of these results,please cite: \n\n");
    fprintf(fp, "Guo, S.-W., Thompson, E. T. (1992) \n");
    fprintf(fp, "Performing the Exact Test of Hardy-Weinberg Proportion \n");
    fprintf(fp, "for Multiple Alleles. Biometrics 48:361-372\n");
}


static void make_output_table(const char *outfile_name1, char *outfile_name2, int option)

{
    char *exec_str;
    size_t len, mega_path_pos;


    /* length of make_gen_table.pl(=17) + 2 blanks + 1 null character = 20 */
    len=strlen(outfile_name1) + strlen(outfile_name2) + 22;
    len+=FILENAME_LENGTH+6;
    exec_str= CALLOC(len, char);

    /* now copy the name of the perl executable path into exec_str with
       a space at the end */
    if (option==1)
        strcpy(exec_str, perl_pgm("make_gen_table.pl"));
    else
        strcpy(exec_str, perl_pgm("make_hwe_table.pl"));

    mega_path_pos = strlen(exec_str);
    exec_str[mega_path_pos] = ' ';
    mega_path_pos++;

    /* copy the input file name into the exec_str advancing after the space */
    sprintf(&exec_str[mega_path_pos], "%s ", outfile_name1);
    mega_path_pos+=strlen(outfile_name1) + 1;

    /* then the output file_name */
    sprintf(&exec_str[mega_path_pos], "%s\n", outfile_name2);
    System(exec_str);
    free(exec_str);
    return;
}

char *get_hwe_option(char prog_names[5][8])

{
    /* used to be HWE or GEN,
     */
    int loop=1, hwe_type;
    char asterisks[6], cchoice[10];

    /* Guo or the other */

    strcpy(asterisks,  "*    ");
    if (batchANALYSIS) {
        hwe_type = TO_HWETEST->_suboption;
    } else {
        draw_line();
        hwe_type=1;
        while (1){
            printf("   HARDY WEINBERG PROGRAM SELECTION MENU\n");
            printf("Selections 3 and 4 require the R statistical package\n");
            printf("as well as the genetics library for R to be installed.\n");
            printf("\n0) Done with the menu - please proceed.\n");
            printf("%c1) The GEN program by Lazzeroni and Lange\n",
                   asterisks[0]);
            printf("%c2) The HWE program by Guo and Thompson\n",
                   asterisks[1]);
            printf("%c3) The chi-square test from R genetics library\n", asterisks[2]);
            printf("%c4) The exact test for bi-allelic markers from R genetics library\n",
                   asterisks[3]);
            printf("%c5) MENDEL program for testing disequilibrium\n",
                   asterisks[4]);
            printf("Enter selection: 0 - 5 > ");
            fcmap(stdin, "%s", cchoice); newline;
            sscanf(cchoice, "%d", &loop);

            switch(loop) {
            case 0:
                break;
            case 1:
                hwe_type=loop;
                sprintf(asterisks, "*    ");
                break;
            case 2:
                hwe_type=loop;
                sprintf(asterisks, " *   ");
                break;
            case 3:
                hwe_type=loop;
                sprintf(asterisks, "  *  ");
                break;
            case 4:
                hwe_type=loop;
                sprintf(asterisks, "   * ");
                break;
            case 5:
                hwe_type=loop;
                sprintf(asterisks, "    *");
                break;
            default:
                warn_unknown(cchoice);
                hwe_type=1;
                break;
            }
            if (!loop) break;
        }

        sprintf(err_msg, "Hardy-Weinberg program option: %s",
                prog_names[hwe_type - 1]);
        mssgf(err_msg);

        TO_HWETEST->_suboption = hwe_type;
        strcpy(Mega2BatchItems[/* 6 */ Analysis_Sub_Option].value.name, prog_names[hwe_type - 1]);
        batchf(Analysis_Sub_Option);
    }
    return(prog_names[hwe_type-1]);
}

static void init_hwe_outfile(char **input_files, char *outfile_name, char *OutFileExist, int count_option)

{
    char write_mode[3], count_option_str[255];
    FILE *fp;

    sprintf(count_option_str, "Counted individuals consist of: \n");
    if (strcmp(OutFileExist, "append") == 0) strcpy(write_mode, "a");
    else strcpy(write_mode, "w");

    fp = fopen(outfile_name, write_mode);
    if (strcmp(write_mode, "w") == 0) hwe_copyright(fp);
    switch (count_option) {
    case 2:
        strcat(count_option_str,
               "Genotyped founders or a randomly chosen genotyped person from each pedigree");
        break;
    case 1:
        strcat(count_option_str, "Genotyped founders only");     break;
    case 4:
        strcat(count_option_str, "All genotyped individuals");     break;
    }
    summary_time_stamp(input_files, fp,  count_option_str);
    fclose(fp);
}


static void  hwe_outfile_name(int numchr, char *outfile_name1, char *outfile_name2,
                              char **input_files, int count_option)

{
    /* outfile names for the HWE option */

    char num[3], *OutFileExist1, *OutFileExist2;
    int user_reply, exit_loop=0;
    char flname[2*FILENAME_LENGTH];

    OutFileExist1 = CALLOC((size_t) FILENAME_LENGTH, char);
    OutFileExist2 = CALLOC((size_t) FILENAME_LENGTH, char);

    if (main_chromocnt > 1) {
        sprintf(outfile_name1, "hwe_results.all");
        sprintf(outfile_name2, "hwe_table.all");
    } else {
        CHR_STR(numchr, num);
        sprintf(outfile_name1, "hwe_results.%s", num);
        sprintf(outfile_name2, "hwe_table.%s", num);
    }

    if (!DEFAULT_OUTFILES) {

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

            printf("HWE output file menu\n");

            printf("0) Done with this menu - please proceed\n");
            printf(" 1) HWE long-form output file name: %s\t[%s]\n",
                   outfile_name1, OutFileExist1);
            printf(" 2) HWE table-form output file name: %s\t[%s]\n",
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
    }
    /* output type - table or long format */
    sprintf(flname, "%s/%s", output_paths[0], outfile_name1);
    init_hwe_outfile(input_files, flname, OutFileExist1, count_option);
    sprintf(flname, "%s/%s", output_paths[0], outfile_name2);
    init_hwe_outfile(input_files, flname, OutFileExist2, count_option);
    return;
}


static void hwe_run_options(int *scycles, int *samples, int *dcycles)

{
    /* run parameters for the HWE option */
    int user_value, parameter_option=99;
    char cchoice[10];

    *dcycles=4000; *samples=20; *scycles=5000;
    if (DEFAULT_OPTIONS) return;
    /* dcycles: number of dememorization cycles,
       samples: batch
       scycles: number of batches */
    while (parameter_option !=0){
        draw_line();
        printf("\nSimulation parameters for HWE\n");
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Number of batches               [%4d]\n", *scycles);
        printf(" 2) Batch size                      [%4d]\n", *samples);
        printf(" 3) Number of dememorization cycles [%4d]\n", *dcycles);
        printf("Enter item number > ");
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &parameter_option);
        user_value = 0;
        switch (parameter_option) {
        case 0:
            break;
        case 1:
            printf("Enter number of batches > ");
            fcmap(stdin, "%s", cchoice); newline;
            sscanf(cchoice, "%d", &user_value);
            if (user_value>0)  *scycles = user_value;
            break;
        case 2:
            printf("Enter batch size > ");
            fcmap(stdin, "%d", &user_value);
            printf("\n");
            if (user_value>0)  *samples = user_value;
            break;
        case 3:
            printf("Enter number of dememorization cycles > ");
            fcmap(stdin, "%d", &user_value);
            printf("\n");
            if (user_value>0)  *dcycles = user_value;
            break;
        default:
            warn_unknown(cchoice);
            break;
        }
    }
    return;

}

/* Since we are counting genotypes here, we exclude half-typed individuals */

static int count_genotypes(linkage_ped_top *LPedTreeTop, int locus_id,
                           int count_option, int *geno_count, int xlinked)

{
    int ii, jj, allele1, allele2;
    int num_genotypes=0;
    allelecnt *member_ids;

    for (ii=0; ii < LPedTreeTop->PedCnt; ii++) {
        member_ids = CALLOC((size_t) LPedTreeTop->Ped[ii].EntryCnt, allelecnt);
        select_individuals(LPedTreeTop, locus_id, count_option, ii,
                           member_ids, 0, xlinked);
        for(jj=0; jj < LPedTreeTop->Ped[ii].EntryCnt; jj++) {
            if (member_ids[jj].num == 2) {
                get_2alleles(LPedTreeTop->Ped[ii].Entry[jj].Marker, locus_id, &allele1, &allele2);
                num_genotypes++;
                geno_count[genoindx(allele1, allele2)] ++;
            }
        }
        free(member_ids);
    }
    return num_genotypes;
}

char *hwe_input_file(int *geno_count, int num_alleles, char *locus_name,
                     int scycles, int samples, int dcycles)

{
    /* create an input file for the HWE program for each marker */
    int ii, jj;
    char *infile_name;
    FILE *fin;

    infile_name = CALLOC((strlen(locus_name)+9), char);
    sprintf(infile_name, "hwe_dat.%s", locus_name);
    fin = fopen(infile_name, "w");
    fprintf(fin, "%d\n", num_alleles);
    for (ii=1; ii < num_alleles + 1; ii++){
        for (jj=1; jj < ii+1; jj++)
            fprintf(fin, "%d\t", geno_count[genoindx(ii,jj)]);
        fprintf(fin, "\n");
    }
    fprintf(fin, "%d\t", dcycles);
    fprintf(fin, "%d\t%d", scycles, samples);

    fprintf(fin, "\n");
    fclose(fin);
    return infile_name;
}

static void hwe_R_outfile_names(char *geno_file,
				char *outfile_name,
				int numchr,
				char *hwe_option)

{

    /* outfile names for the HWE option */

    char num[3], *exists, copt[10];
    int user_reply, exit_loop=0;
    char flname[2*FILENAME_LENGTH];

    exists = (char *) CALLOC((size_t) 10, char);
    sprintf(geno_file, "hwe_genos");

    if (main_chromocnt > 1) {
        sprintf(outfile_name, "hwe_%s_table.all", hwe_option);
    } else {
        CHR_STR(numchr, num);
        sprintf(outfile_name, "hwe_%s_table.%s", hwe_option, num);
    }

    if (!DEFAULT_OUTFILES) {
        while (!exit_loop) {
            sprintf(flname, "%s/%s", output_paths[0], outfile_name);
            ((access(flname, F_OK) == 0)?
             strcpy(exists, "append") : strcpy(exists, "new"));
            draw_line();
            print_outfile_mssg();
            printf("HWE output file name menu:\n");
            printf("0) Done with this menu - please proceed\n");
            printf(" 1) Genotype file name stem: %s\n", geno_file);
            printf(" 2) Table-form output file name: %s\t[%s]\n",
                   outfile_name, exists);
            printf("Enter selection: 0, 1 or 2 > ");
            fcmap(stdin, "%s", copt); newline;
            sscanf(copt, "%d", &user_reply);
            test_modified(user_reply);
            switch (user_reply) {
            case 0:
                exit_loop = 1;
                break;
            case 1:
                printf("Enter new genotype file name stem > ");
                fcmap(stdin, "%s", geno_file); newline;
                break;
            case 2:
                printf("Enter new table-form output file name > ");
                fcmap(stdin, "%s", outfile_name); newline;
                break;
            default:
                warn_unknown(copt);
                break;
            }
        }
    }
    return;

}

static FILE *init_hwe_R_script(char *R_script)

{

    FILE *Rfl;
    Rfl = fopen(R_script, "w");
    fprintf(Rfl, "# R-script file %s\n", R_script);
    script_time_stamp(Rfl);
    fprintf(Rfl, "library(\"genetics\");\n");
    fprintf(Rfl, "pvalues<-NULL;\n");
    fprintf(Rfl, "g11<-NULL;\n");
    fprintf(Rfl, "g12<-NULL;\n");
    fprintf(Rfl, "g22<-NULL;\n");
    fprintf(Rfl, "chi1<-NULL;\n");
    fprintf(Rfl, "chi2<-NULL;\n");
    fprintf(Rfl, "stat<-NULL;\n");
    return(Rfl);
}

static void Rfl_footer(FILE *Rfl, char *hwe_opt,
		       char *outfile, int all_one)

{
    if (!strcmp(hwe_opt, "EXACT")) {
        if (all_one == 1) {
            fprintf(Rfl ,"output<-cbind(g11, g12, g22, pvalues)\n");
        } else {
            fprintf(Rfl ,"output<-cbind(g22, g12, g11, pvalues)\n");
        }
    } else {
        fprintf(Rfl, "output<-cbind(chi1, chi2, stat, pvalues)\n");
    }
    fprintf(Rfl,
            "ou.table <- data.frame(output, row.names=markers);\n");
    fprintf(Rfl, "names(ou.table) <- \n");
    if (!strcmp(hwe_opt, "EXACT")) {
        fprintf(Rfl, "  c(\"1/1\", \"1/2\", \"2/2\", \"P-value\");\n");
    } else {
        fprintf(Rfl,
                "  c(\"max(exp-obs)\", \"max2(exp-obs)\", \"Chi-sq\", \"P-value\");\n");
    }
    fprintf(Rfl,
            "write.table(format(ou.table), file=\"%s\", quote=F, append=T);\n",
            outfile);
    fclose(Rfl);

}

static void save_mendel_hwe_individuals(linkage_ped_top *Top,
					int count_option, int *loci_indexes,
					int num_loci, char *geno_file,
					int csv_format, int xlinked)


{

    int m, i;
    char fformat[6], pformat[6];
    int mendel_ind_number=0;
    int fwid, pwid, ped, ind;
    allelecnt **select, *ind_selected;
    FILE *hwe_in = fopen(geno_file, "w");

    /* Now print the selected individuals */
    field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);

    if (fwid >= 7) fwid = 8;
    else fwid++;

    if (pwid >= 7) pwid=8;
    else pwid++;

    if (!csv_format) {
        fprintf(hwe_in, "(I2,1x,A%d)\n", fwid);
        fprintf(hwe_in, "(3A%d,1X,2A1,(T%d,3(A8),:))\n", pwid, 3*pwid+4);
        create_mendel_formats(fwid, pwid, fformat, pformat);
    }

    select = CALLOC((size_t) num_loci, allelecnt *);
    /* This is a matrix of selection for each locus */
    for (ped = 0; ped < Top->PedCnt; ped++) {
        ind_selected=CALLOC((size_t) Top->Ped[ped].EntryCnt, allelecnt);
        /* array of flags to denote if individual has been selected
           at any marker */
        for (ind = 0; ind < Top->Ped[ped].EntryCnt; ind++)  {
            ind_selected[ind].num=0;
        }
        /* Now select individuals for each locus */
        for (m=0; m < num_loci; m++) {
            select[m]=CALLOC((size_t) Top->Ped[ped].EntryCnt, allelecnt);
            select_individuals(Top, loci_indexes[m],
                               count_option, ped, select[m], 0, xlinked);
        }
        /* set the ind_selected flags */
        for (ind = 0; ind < Top->Ped[ped].EntryCnt; ind++)  {
            for (m = 0; m < num_loci; m++) {
                if (select[m][ind].num) {
                    ind_selected[ind].num++;
                }
            }
            if (!ind_selected[ind].num) {
                continue;
            }
            /* output only if the individual has been selected */
            mendel_ind_number++;
            /* write mendel ped number */
            if (!csv_format) {
                fprintf(hwe_in, " 1");
                fprintf(hwe_in, fformat, mendel_ind_number);
                fprintf(hwe_in, "\n");

                /* write individual */
                fprintf(hwe_in, pformat, Top->Ped[ped].Entry[ind].ID);
                /* write blanks for parents */
                for (i=0; i < 2*pwid; i++) {
                    fprintf(hwe_in, " ");
                }
                /* write sex */
                if (Top->Ped[ped].Entry[ind].Sex == MALE_ID)    fprintf(hwe_in, " M ");
                if (Top->Ped[ped].Entry[ind].Sex == FEMALE_ID)  fprintf(hwe_in, " F ");
            } else {
                /* ped number */
                fprintf(hwe_in, "%d,", mendel_ind_number);
                /* indiv id */
                fprintf(hwe_in, "%d,", Top->Ped[ped].Entry[ind].ID);
                /* parents and sex */
                if (Top->Ped[ped].Entry[ind].Sex == MALE_ID)
                    fprintf(hwe_in, ",,M,");
                if (Top->Ped[ped].Entry[ind].Sex == FEMALE_ID)
                    fprintf(hwe_in, ",,F,");

            }

            /* write genotypes */
            for (m = 0; m < num_loci; m++) {
                if (!csv_format) {
                    if (m > 0 && (m % 3) == 0) {
                        fprintf(hwe_in, "\n");
                        for(i=0; i < 3*pwid+3; i++) {
                            fprintf(hwe_in, " ");
                        }
                    }
                }

                if (select[m][ind].num) {
                    int a1, a2;
                    linkage_locus_rec *locus = &Top->LocusTop->Locus[m];
                    get_2alleles(Top->Ped[ped].Entry[ind].Marker, loci_indexes[m], &a1, &a2);
                    if (csv_format) {
                        fprintf(hwe_in, ",%s/%s", format_allele(locus, a1), format_allele(locus, a2));
                    } else {
                        fprintf(hwe_in, "%3s/%3s", format_allele(locus, a1), format_allele(locus, a2));
                    }
                } else {
                    fprintf(hwe_in, ",");
                }
            }
            fprintf(hwe_in, "\n");
        }
        free(ind_selected);
        /* Now free the select data structure */
        for (m=0; m < num_loci; m++) {
            free(select[m]);
            select[m]=NULL;
        }
    }
    free(select);
    fclose(hwe_in);
    return;
}

static void hwe_R_setup(char *hwe_option, linkage_ped_top *Top,
			int count_option, int *loci_indexes,
			int num_loci,
			char *geno_file, char *R_script,
			char *outfile, int inc_ht, int xlinked)

{
    int m, ind, ped;
    allelecnt **member_ids;
    FILE *hwe_in, *Rfl;
    char mrkfl[2*FILENAME_LENGTH];
    int num_halftyped, num_fulltyped=0;
    int set_first_allele=1, first_allele=1;

    /*  Step 1: Initialize the R-file */
    Rfl=init_hwe_R_script(R_script);
    fprintf(Rfl, "markers <- c(");
    /* Step 2: Create a file containing list of alleles */
    /* for only valid markers   */

    member_ids = CALLOC((size_t) Top->PedCnt, allelecnt *);
    for (ped=0; ped < Top->PedCnt; ped++) {
        member_ids[ped]= CALLOC((size_t) Top->Ped[ped].EntryCnt, allelecnt);
    }

    for (m = 0; m < num_loci; m++) {
        linkage_locus_rec *locus = &Top->LocusTop->Locus[loci_indexes[m]];
        /* initialize R file for every marker */
        sprintf(mrkfl, "%s.%s", geno_file, locus->LocusName);
        hwe_in = fopen(mrkfl, "w");
        fprintf(hwe_in, "Pedigree Individual allele1\tallele2\n");

        /* initalize selected matrix for all pedigrees and individuals */
        for (ped=0; ped < Top->PedCnt; ped++) {
            for (ind =0; ind < Top->Ped[ped].EntryCnt; ind++) {
                member_ids[ped][ind].num=0;
            }
        }
        /* select individuals for each pedigree, then write out only selected
           individuals */
        num_halftyped=0;
        num_fulltyped=0;
        for (ped=0; ped < Top->PedCnt; ped++) {
            select_individuals(Top, loci_indexes[m], count_option, ped, member_ids[ped],
                               inc_ht, xlinked);
            for (ind =0; ind < Top->Ped[ped].EntryCnt; ind++) {
                int a1, a2;
                if (member_ids[ped][ind].num >= 1) {
                    get_2alleles(Top->Ped[ped].Entry[ind].Marker, loci_indexes[m], &a1, &a2);
                    if (set_first_allele) {
                        first_allele = a1;
                        set_first_allele = 0;
                    }
                    fprintf(hwe_in, "%6d     %5d      ",
                            Top->Ped[ped].Num,
                            Top->Ped[ped].Entry[ind].ID);
                    fprintf(hwe_in, "%3s     %3s \n", format_allele(locus, a1), format_allele(locus, a2));
                    if (member_ids[ped][ind].num == 1) {
                        num_halftyped++;
                    } else {
                        num_fulltyped++;
                    }
                }
            }
        }
        fclose(hwe_in);
        /* Log the number of half-types */
        sprintf(err_msg, "Marker %s: Found %d fully-typed", locus->LocusName, num_fulltyped);
        if (num_halftyped > 0) {
            grow(err_msg, " and %d half-typed individuals.",
                 num_halftyped);
        } else {
            strcat(err_msg, " individuals.");
        }
        mssgf(err_msg);

        /* add the marker name to the list in the R-script */
        if (m < (num_loci-1)) {
            fprintf(Rfl, "\"%s\", ", locus->LocusName);
        } else {
            fprintf(Rfl, "\"%s\");\n", locus->LocusName);
        }
    }

    /* Step 2: Create an R command to process these markers */
    fprintf(Rfl, "unlink(\"%s\")\n", outfile);
    fprintf(Rfl, "for (i in 1:length(markers)) {\n");
    fprintf(Rfl, "  file<-paste(\"%s\", markers[i], sep=\".\");\n",
            geno_file);
    fprintf(Rfl, "  alleles <- read.table(file, header=T);\n");
    fprintf(Rfl, "  genos <- \n");
    fprintf(Rfl,
            "      genotype(alleles$allele1, alleles$allele2, reorder=\"no\");\n");

    if (!strcmp(hwe_option, "EXACT")) {
        /* store pvalue, and observed genotypes */
        fprintf(Rfl, "  hwe <- HWE.exact(genos);\n");
        fprintf(Rfl, "  pvalues[i] <- round(hwe$p.value,3);\n");
        fprintf(Rfl, "  g11[i]     <- round(hwe$statistic[1],0);\n");
        fprintf(Rfl, "  g12[i]     <- round(hwe$statistic[2],0);\n");
        fprintf(Rfl, "  g22[i]     <- round(hwe$statistic[3],0);\n");
    } else {
        /* store pvalue, expected and observed count and chisq value */
#ifdef TEST
        fprintf(Rfl, "  hwe <- HWE.chisq(genos,simulate.p.value=FALSE);\n");
#else
        fprintf(Rfl, "  hwe <- HWE.chisq(genos);\n");
#endif
        fprintf(Rfl, "  pvalues[i] <- round(hwe$p.value,3);\n");
        fprintf(Rfl, "# get the two largest differences values\n");
        fprintf(Rfl, "  diff <- abs(hwe$expected - hwe$observed);\n");
        fprintf(Rfl, "  chi1[i] <- max(diff);\n");
        fprintf(Rfl, "  diff[diff == max(diff)] <- 0;\n");
        fprintf(Rfl, "  chi2[i] <- max(diff);\n");
        fprintf(Rfl, "  stat[i]    <- round(hwe$statistic,3);\n");
    }

    fprintf(Rfl, "}\n");
    Rfl_footer(Rfl, hwe_option, outfile, first_allele);

    return;
}

static void write_mendel_hwe_control(char *file_names[], char *chr_str,
				     int csv_format)

{

    char flname[2*FILENAME_LENGTH];
    FILE *fp;

    sprintf(flname, "%s/%s", output_paths[0], file_names[3]);
    fp=fopen(flname, "w");
    fprintf(fp, "OUTPUT_FILE = Mendel.%s.out\n", chr_str);
    fprintf(fp, "DEFINITION_FILE = %s\n", file_names[1]);
    fprintf(fp, "PEDIGREE_FILE = %s\n", file_names[0]);
    fprintf(fp, "MAP_FILE = %s\n", file_names[2]);
    fprintf(fp, "SUMMARY_FILE=Summary.%s.out\n", chr_str);
    fprintf(fp, "ANALYSIS_OPTION = Genetic_equilibrium\n");
    fprintf(fp, "NUMBER_OF_MARKERS_INCLUDED=1\n");
    if (csv_format) {
        fprintf(fp, "DEFAULT_LIST_READ=TRUE\n");
    }
    fclose(fp);

}


static void mendel_hwe_file_names(char *file_names[], int numchr,
				  int *csv_format)

{

    char chr[3], flname[2*FILENAME_LENGTH], fl_stat[12];
    char cchoice[10];
    int opt, exit_loop=0;

    *csv_format=1;

    if (main_chromocnt > 1) {
        sprintf(file_names[0], "hwe_mendel_ped.all");
        sprintf(file_names[1], "hwe_mendel_loc.all");
        sprintf(file_names[2], "hwe_mendel_map.all");
        sprintf(file_names[3], "hwe_mendel_control.all");
    } else {
        CHR_STR(numchr, chr);
        sprintf(file_names[0], "hwe_mendel_ped.%s", chr);
        sprintf(file_names[1], "hwe_mendel_loc.%s", chr);
        sprintf(file_names[2], "hwe_mendel_map.%s", chr);
        sprintf(file_names[3], "hwe_mendel_control.%s", chr);
    }

    if (DEFAULT_OUTFILES) {
        return;
    }
    while (!exit_loop) {
        draw_line();
        print_outfile_mssg();
        printf("HWE-MENDEL output file name menu:\n");
        printf("0) Done with this menu - still proceed\n");
        printf(" 1) Create CSV-format files:  %s\n",
               yorn[*csv_format]);
        sprintf(flname, "%s/%s", output_paths[0], file_names[0]);
        printf(" 2) Pedigree file name:       %s\t%s\n",
               file_names[0], file_status(flname, fl_stat));
        sprintf(flname, "%s/%s", output_paths[0], file_names[1]);
        printf(" 3) Locus file name:          %s\t%s\n",
               file_names[1], file_status(flname, fl_stat));
        sprintf(flname, "%s/%s", output_paths[0], file_names[2]);
        printf(" 4) Map file name:            %s\t%s\n",
               file_names[2], file_status(flname, fl_stat));
        sprintf(flname, "%s/%s", output_paths[0], file_names[3]);
        printf(" 5) Control file name:        %s\t%s\n",
               file_names[3], file_status(flname, fl_stat));

        printf("Select options 0 - 5, 1 to toggle > ");
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &opt);
        test_modified(opt);

        switch(opt) {
        case 0:
            exit_loop=1;
            break;

        case 1:
            *csv_format = TOGGLE(*csv_format);
            break;

        case 3:
            printf("New locus file name > ");
            fcmap(stdin, "%s", file_names[1]);
            newline;
            break;

        case 2:
            printf("New pedigree file name > ");
            fcmap(stdin, "%s", file_names[0]);
            newline;
            break;

        case 4:
            printf("New map file name > ");
            fcmap(stdin, "%s", file_names[2]);
            newline;
            break;

        case 5:
            printf("Name of control file name > ");
            fcmap(stdin, "%s", file_names[3]);
            newline;
            break;

        default:
            warn_unknown(cchoice);
        }
    }
    return;
}

void hwe_user_input(linkage_ped_top *LPedTreeTop, int *numchr,
                    char *file_names[])

{
    /* selects the analysis method, reorders and selects loci, then runs HWE or GEN
       on the selected loci. After analysis, parses the output appropriately */

    char *infile_name, *outfile_name1, *outfile_name2;
    int num_loci, reject_type, count_option, incl_halftyped;
    int dcycles, samples, scycles;
    int xlinked = 0,  ii, ii1, kk, kkk, dummy;
    int *geno_count, num_genotypes;
    char *hwe_option;
    int *loci_indxs;
    double *homo_e = NULL, *homo_o = NULL;
    linkage_locus_top *ltop;
    FILE *fp1;
    char **input_files, *exec_cp, *gen_executable_name;
    char R_script[30], chr[4];
    char prog_names[5][8] = { "GEN", "HWE", "CHISQ", "EXACT", "MENDEL"};

#define _mrk_name(i)   (ltop->Locus[i].LocusName)
#define _num_all(i)   (ltop->Locus[i].AlleleCnt)

    input_files=mega2_input_files;
    sprintf(R_script, "HWE_test.R");
    /* removed initialization of seed inside mega2.c */

    /* omit untyped pedigrees */
/*   omit_peds(untyped_ped_opt, LPedTreeTop); */
    /* get the analysis option */
    hwe_option = CALLOC((size_t) 16, char);
    strcpy(hwe_option, get_hwe_option(prog_names));

    /* select_loci for analysis */
    ltop=LPedTreeTop->LocusTop;
    num_loci=0;
    /* see if we need to add unmapped loci */

    get_loci_on_chromosome(0);
    for (ii=0; ii < main_chromocnt; ii++) {
        if (global_chromo_entries[ii] == UNKNOWN_CHROMO) {
            get_unmapped_loci(1);
            break;
        }
    }

    loci_indxs = (int *) CALLOC((size_t) NumChrLoci, int);

    for (ii1 = 0; ii1 < NumChrLoci; ii1++) {
        ii=ChrLoci[ii1];
        if (ltop->Marker[ii].chromosome == MALE_CHROMOSOME) {
            warnvf("No HWE test for Y-linked loci, skipping %s.\n", _mrk_name(ii));
            continue;
        }
        xlinked =
            (((ltop->Marker[ii].chromosome == SEX_CHROMOSOME &&
               LPedTreeTop->LocusTop->SexLinked == 2) ||
              (LPedTreeTop->LocusTop->SexLinked == 1))? 1 : 0);

        reject_type=0;
        if (ltop->Locus[ii].Type != NUMBERED) reject_type =1;
        else if (_num_all(ii) < 3) reject_type =2;
        else if (_num_all(ii) > MAX_ALLELE) reject_type = 3;
        else if (!strcmp(hwe_option, "EXACT") && (_num_all(ii) > 2)) reject_type = 4;

        switch (reject_type) {
        case 0:
            loci_indxs[num_loci] = ii;
            num_loci ++;
            /* fills up loci_indxs from 1 to n+1 */
            break;
        case 1:
            warnvf("Excluding locus %s from analysis (not numbered type).\n",
                    _mrk_name(ii));
            break;
        case 2:
            if (strcmp(hwe_option, "EXACT") &&
                strcmp(hwe_option, "CHISQ") &&
                strcmp(hwe_option, "MENDEL")) {
                warnvf("Excluding locus %s from analysis (less than 3 alleles).\n",
                        _mrk_name(ii));
            } else {
                loci_indxs[num_loci] = ii;
                num_loci++;
            }
            break;
        case 3:
            warnvf("Excluding locus %s from analysis (more than %d alleles).\n",
                    _mrk_name(ii), MAX_ALLELE);
            break;
        case 4:
            warnvf("Excluding locus %s from analysis (more than 2 alleles).\n",
                    _mrk_name(ii));
            break;
        default:
            break;
        } // switch (reject_type)
    } // for (ii1 = 0; ii1 < NumChrLoci; ii1++) {
    if (num_loci == 0) {
      errorvf("All locus have been excluded from analysis.\n");
      EXIT(OUTPUT_FORMAT_ERROR);
    }
    /* get count option */
    draw_line();

    count_option=
        get_count_option(0, &dummy,
                         "Select individuals to test Hardy-Weinberg equilibrium on:");
    incl_halftyped=((!HalfTypedReset &&
                     (!strcmp(hwe_option, "EXACT") || !strcmp(hwe_option, "CHISQ")) ) ?
                    1: 0);
    if (count_option == 3) {
        /* To match recode menu */
        count_option=4;
    }

    if (!(strcmp(hwe_option, "CHISQ")) ||
        !(strcmp(hwe_option, "EXACT"))) {
        /* for R options only */
        outfile_name1 = CALLOC((size_t) 2*FILENAME_LENGTH, char);
        outfile_name2 = CALLOC((size_t) 2*FILENAME_LENGTH, char);
        hwe_R_outfile_names(file_names[0], file_names[1], *numchr,
                            hwe_option);

        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[0]);
        sprintf(outfile_name2, "%s/%s", output_paths[0], file_names[1]);

        hwe_R_setup(hwe_option, LPedTreeTop, count_option,
                    loci_indxs, num_loci, outfile_name1, R_script,
                    outfile_name2, incl_halftyped, xlinked);
        printf("Now running R to compute the HWE statistics ...");
        fflush(stdout);
        exec_cp = CALLOC((size_t) 25, char);
        sprintf(exec_cp, "R CMD BATCH %s", R_script);
        System(exec_cp);
        if (access(outfile_name2, F_OK) != 0) {
            warnvf("Unable to create %s, consult %sout for details\n",
		   file_names[1], R_script);
        }

        newline; log_line(mssgf);
        create_mssg(TO_HWETEST);
        for (ii=0; ii < num_loci; ii++) {
            sprintf(err_msg, "    Genotypes file:             %s.%s",
                    file_names[0], LPedTreeTop->LocusTop->Locus[loci_indxs[ii]].LocusName);
            mssgf(err_msg);
        }
        if (access(outfile_name1, F_OK) == 0) {
            sprintf(err_msg, "    Long-form output file:      %s", file_names[2]);
        }
        if (access(outfile_name2, F_OK) == 0) {
            sprintf(err_msg, "    Table-form output file:     %s", file_names[1]);
            mssgf(err_msg);
        }
        /*    delete_file(R_script); */
        free(exec_cp);
        free(outfile_name1); free(outfile_name2);
        OrigIds[0]=5;
        OrigIds[1]= (basefile_type == 0)? 1 : 2;
    } else if (!strcmp(hwe_option, "MENDEL")) {
        int csv_format;
        int sex_linked =
            (((*numchr == SEX_CHROMOSOME && LPedTreeTop->LocusTop->SexLinked == 2) ||
              (LPedTreeTop->LocusTop->SexLinked == 1))? 1 : 0);

        mendel_hwe_file_names(file_names, *numchr, &csv_format);
        /*    rename_loci(LPedTreeTop->LocusTop); */
        CHR_STR(*numchr, chr);

        create_mssg(TO_HWETEST);
        outfile_name1 = CALLOC((size_t) 2*FILENAME_LENGTH, char);
        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[0]);
        init_file(outfile_name1);
        save_mendel_hwe_individuals(LPedTreeTop,
                                    count_option, loci_indxs,
                                    num_loci, outfile_name1, csv_format,
                                    sex_linked);

        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[1]);
        init_file(outfile_name1);
        if (csv_format) {
            csv_write_mendel_locus_file(file_names[1],
                                        LPedTreeTop->LocusTop,
                                        sex_linked);

        } else {
            write_mendel_locus_file(file_names[1],
                                    LPedTreeTop->LocusTop, TO_HWETEST,
                                    sex_linked);
        }

        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[2]);
        init_file(outfile_name1);
        if (csv_format) {
            csv_write_mendel_map(file_names[2], LPedTreeTop->LocusTop, TO_HWETEST,
                                 sex_linked);
        } else {
            write_mendel_map(file_names[2],
                             LPedTreeTop->LocusTop, TO_HWETEST);
        }

        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[3]);
        init_file(outfile_name1);
        write_mendel_hwe_control(file_names, chr, csv_format);
        sprintf(err_msg, "        Pedigree file:        %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "        Control file:         %s", file_names[3]);
        mssgf(err_msg);
        free(outfile_name1);
    } else if (strcmp(hwe_option, "HWE")==0) {

        /* get the output file names and write copyright notice and time stamp */
        outfile_name1 = CALLOC((size_t) 2*FILENAME_LENGTH, char);
        outfile_name2 = CALLOC((size_t) 2*FILENAME_LENGTH, char);
        hwe_outfile_name(*numchr, file_names[0], file_names[1],
                         input_files, count_option);
        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[0]);
        sprintf(outfile_name2, "%s/%s", output_paths[0], file_names[1]);

        /* get run parameters */
        hwe_run_options(&scycles, &samples, &dcycles);

        /* execute hwe */
        /* initialize the temp file */
        fp1 = fopen("hwe_tmp", "w");    fclose(fp1);

        printf("Counting genotypes...\n");
        printf("The code for HWE is included in mega2 by permission of \n");
        printf("       Sun-Wei Guo and Elizabeth Thompson.\n");
        printf("Any publications using results generated by this HWE option \n");
        printf("should cite this paper:\n");
        printf("Guo, S.-W., Thompson, E. T. (1992) \n");
        printf("Performing the Exact Test of Hardy-Weinberg Proportion \n");
        printf("for Multiple Alleles. Biometrics 48:361-372\n");
        log_line(mssgf);

        /* initialize the gen_count and homozygosity data structures */
        homo_e= CALLOC((size_t) num_loci, double);
        homo_o= CALLOC((size_t) num_loci, double);

        for (kk=0; kk < num_loci; kk++){
            kkk = loci_indxs[kk];
            geno_count=(int *) CALLOC((size_t) NUMGENOS(_num_all(kkk)), int);
            num_genotypes=count_genotypes(LPedTreeTop, kkk, count_option,
                                          geno_count, xlinked);
            homo_e[kk]=expected_homozygosity(LPedTreeTop, kkk);
            homo_o[kk]=observed_homozygosity(geno_count, _num_all(kkk));
            infile_name=hwe_input_file(geno_count, _num_all(kkk), _mrk_name(kkk),
                                       scycles, samples, dcycles);

            hwe_test(_mrk_name(kkk), infile_name, "hwe_tmp",
                     homo_o[kk], homo_e[kk], num_genotypes);
            sprintf(err_msg, "    Marker %10s", _mrk_name(kkk));
            mssgf(err_msg);
            delete_file(infile_name);
            free(infile_name); infile_name = NULL;
            free(geno_count); geno_count=NULL;
        }

        /* append hwe_tmp to outfile1, the raw output file,
           have to create a second temp file to use cat */
        exec_cp = CALLOC((size_t) 13 + strlen(outfile_name1), char);
        sprintf(exec_cp, "cp %s hwe_tmp1", outfile_name1);
        System(exec_cp);
        free(exec_cp); exec_cp=NULL;
        exec_cp = CALLOC((size_t) 24 + strlen(outfile_name1), char);
        sprintf(exec_cp, "cat hwe_tmp1 hwe_tmp > %s", outfile_name1);
        System(exec_cp);
        free(exec_cp); exec_cp = NULL;
        draw_line();
        create_mssg(TO_HWETEST);
        make_output_table("hwe_tmp", outfile_name2, 2);
        sprintf(err_msg, "    Long-form output file:      %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "    Table-form output file:     %s", file_names[1]);
        mssgf(err_msg);
        draw_line();
        System("/bin/rm -f hwe_tmp hwe_tmp1");
        free(outfile_name1); free(outfile_name2);
        outfile_name1=NULL; outfile_name2=NULL;
    } else if (!strcmp(hwe_option, "GEN")) {
        /* GEN option */
        /* get the output file names and write copyright notice and time stamp */
        outfile_name1 = CALLOC((size_t) FILENAME_LENGTH, char);
        outfile_name2 = CALLOC((size_t) FILENAME_LENGTH, char);

        gen_outfile_name(*numchr, file_names[0], file_names[1]);

        sprintf(outfile_name1, "%s/%s", output_paths[0], file_names[0]);
        sprintf(outfile_name2, "%s/%s", output_paths[0], file_names[1]);

        init_gen_outfile(input_files, outfile_name1, count_option);
        init_gen_outfile(input_files, outfile_name2, count_option);

        /* get run parameters */
        gen_run_options(&scycles, &samples, &dcycles);

        gen_executable_name= CALLOC((size_t) 4, char);
        strcpy(gen_executable_name, "Gen");
        /* Execute gen */
        /* initialize the temp file */
        fp1 = fopen("gen_tmp", "w");
        fprintf(fp1, "\n");
        fclose(fp1);

        /* initialize the gen_count and homozygosity data structures */
        homo_e= CALLOC((size_t) num_loci, double);
        homo_o= CALLOC((size_t) num_loci, double);

        printf("Counting genotypes...\n");
        for (kk=0; kk < num_loci; kk++){
            kkk = loci_indxs[kk];
            printf("Marker %s ...", _mrk_name(kkk));
            geno_count= CALLOC((size_t) NUMGENOS(_num_all(kkk)), int);
            num_genotypes = count_genotypes(LPedTreeTop, kkk, count_option,
                                            geno_count, xlinked);
            homo_e[kk]=expected_homozygosity(LPedTreeTop, kkk);
            homo_o[kk]=observed_homozygosity(geno_count, _num_all(kkk));

            printf("%5.4f, %5.4f \n", homo_o[kk], homo_e[kk]);
            infile_name=gen_input_file(geno_count, _num_all(kkk), _mrk_name(kkk),
                                       scycles, samples, dcycles);
            execute_gen(gen_executable_name, _mrk_name(kkk), "gen_tmp",
                        homo_o[kk], homo_e[kk], num_genotypes);
            /* append gen_tmp to outfile1, the raw output file */
            delete_file(infile_name);
            free(infile_name); infile_name = NULL;
            free(geno_count); geno_count=NULL;
        }
        draw_line();
        fcopy("gen_tmp", outfile_name1, "a");
        create_mssg(TO_HWETEST);
        make_output_table("gen_tmp", outfile_name2, 1);
        sprintf(err_msg, "    Long-form output file:      %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "    Table-form output file:     %s", file_names[1]);
        mssgf(err_msg);
        draw_line();
        System("/bin/rm -f gen_tmp gen_tmp1");
        free(outfile_name1); free(outfile_name2);
        outfile_name1=NULL; outfile_name2=NULL;
    }

    /* finished processing, free datastructures */
    free(loci_indxs);
    if (!strcmp(hwe_option, "HWE") || !strcmp(hwe_option, "GEN")) {
        free(homo_o); free(homo_e);
    }
#undef _num_all
#undef _mrk_name
    return;
}
