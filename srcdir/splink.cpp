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

/* Logs incorporated from output_routines.c */

/* Revision 1.20.2.11.2.2  2001/03/26 17:49:22  nandita
   1) Since SPLINK does not handle liability classes, added a call to
   define_affection_labels.
   2) Added a Space at the start of every genotype field.
   3) Writing the affection data explicitly instead of using write_affection_data().
   TESTING:  on "bevan" data.

   Revision 1.20.2.11.2.1  2001/03/26 17:03:44  nandita
   In keeping with similar linkage format pedigree output, added a SPACE
   *before* the numbered marker genotype.

   Revision 1.20.2.11  2000/11/22 17:23:32  dweeks
   DEW: Changed some printf("ERROR 's to fprintf("ERROR & tightened some err msgs.

   Revision 1.20.2.7  2000/08/16 13:40:45  nandita
   fixed splink help option

   Revision 1.20.2.6  2000/08/16 13:16:39  nandita
   corrected bug in splink options menu

   Revision 1.20.2.5  2000/08/14 16:12:19  dweeks
   DW: Removed instructions re SPLINK 'help' command, as the help command
   does not work.  Also corrected some of the punctuation re the SPLINK
   options.

*/


#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc my_malloc
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome new_lpedtop
           makenucs_ext.h:  create_nuked_fids makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  ped_name_width person_id_width
         user_input_ext.h:  enter_number test_modified
              utils_ext.h:  EXIT draw_line press_return script_time_stamp
*/


/*-------------- static ------------*/
static void CreateSplinkCShell(char *options, int single_loci);
static void            splink_help(char *opt);
static char            *splink_option_check(char *option);

/*------------- to be exported --------------- */

void create_SPLINK(linkage_ped_top ** Top,
		   analysis_type * analysis,  int *numchr, char *file_names[],
		   int untyped_ped_opt);

/*------------prototypes------- */

/*----------------display help on splink options---------*/

static void splink_help(char *opt)

{

    draw_line();
    printf("Help: SPLINK option descriptions: ");
    switch (*opt)
        {
        case 'a':
            printf(" flag '-asp'\n");
            printf("Include families w/at least an affected sib pair:\n");
            printf("By default other families are included since they contribute information\n");
            printf("about haplotype frequencies.\n");
            break;
        case 'b':
            printf(" flag '-bs#', where '#' is number of bootstrap samples\n");
            printf("Do bootstrap significance test using # bootstrap samples.\n");
            printf("By default bootstrap testing is not performed.\n");
            break;
        case 'c':
            printf(" flag '-c'\n");
            printf("Include only 'complete' families (i.e. families with both parents typed)\n");
            break;
        case 'd':
            printf(" flag '-d#'\n");
            printf("Write locus data to file #. If # is absent, the data are appended\n");
            printf("to the file datafile.dat.");
            break;
        case 'f':
            printf(" flag '-f#'\n");
            printf("Specify MAXIMUM number of nuclear families, (default = 500) \n");
            break;
        case 'h':
            printf(" flag '-h'\n");
            printf("Help (lists command line options).\n");
            break;
        case 'i':
            printf(" flag '-i#'\n");
            printf("Specify # as an identifier string to be used to label output.  \n");
            printf("If not specified, the string ??? is used.");
            break;
        case 'l':
            printf(" flag '-l#'\n");
            printf("Specify number of of marker loci.\n");
            printf("If missing,  it is assumed that the number is given on the input file.");
            break;
        case 'm':
            printf(" flag '-mf'\n");
            printf("Allow multiple nuclear families from one pedigree \n");
            printf("( although the relationship between these families will be\n");
            printf("  ignored). If not set, only the first nuclear family\n");
            printf("  encountered  in each pedigree is used. Default is for this flag to\n");
            printf("  be set, but it may be unset by -nimf or -mf- .");
            break;
        case 'n':
            printf("  flag '-n#'\n");
            printf("Specify maximum of persons on data file, (default = 5000)\n");
            break;
        case 'o':
            printf(" flag '-o#'\n");
            printf("Specify that the IBD assignment scores for each family will be\n");
            printf("written to file #. By default this option is turned off\n");
            break;
        case 'p':
            printf(" flag '-pt'\n");
            printf("Specifies that estimation of z0, z1, z2 will be restricted to values\n");
            printf("compatible with linkage. For an autosomal marker this is the\n");
            printf("'possible triangle' and for an X marker, z1 is constrained to be\n");
            printf("no smaller than z0. Default is for the flag to be set, but it can\n");
            printf("be unset by -nopt or -pt- .\n");
            break;
        case 'r':
            printf(" flag '-rs#'\n");
            printf("Specify an integer seed #, for the pseudo-random number generator.\n");
            printf("This is used in bootstrap testing (see option 2).\n");
            break;
        case 's':
            printf(" flag '-s#'\n");
            printf("Selects score test. Permitted values of # are,\n");
            printf("for autosomal markers:\n");
            printf("       2	2 df test\n");
            printf("       pt	1-2 df 'possible triangle' test\n");
            printf("       1d	1 df test with the constraint z1 = 0.5\n");
            printf("       1r	1 df test with the constraint z1 = 2 * z0\n");
            printf("for X markers (only maternal haplotype scored)\n");
            printf("       m        1 df test of null z0 = z1 = 0.5\n");
            printf("       m1       1 df 1-sided test z1 > z0 ( z2 = 0 )\n");
            printf("More than one -s option may be selected, but the first\n");
            printf("three are only legal for autosomal markers and that the last\n");
            printf("two for X markers. By default, if the -pt flag is set,\n");
            printf("the 'pt' test is selected for an autosomal marker and\n");
            printf("the 'm1' test is selected for an X marker. All other tests\n");
            printf("are unselected. If the -pt flag is unset, the default score tests\n");
            printf("are the '2' and 'm' options.");
            break;
        case 'S':
            printf(" flag '-S+'\n");
            printf("Write brief ONE LINE summary to file #. If absent, the output is appended\n");
            break;
        case 'x':
            printf(" flag '-x#'\n");
            printf("Specify maximum allowable ambiguity for parental haplotyping.\n");
            printf("If there are more than # possible parental haplotypes assignments,\n");
            printf("the family is excluded. Note that, for speed, possible parental haplotypes\n");
            printf("are stored in dynamically allocated memory and the -x option will help\n");
            printf("you if you run out of memory.");
            break;
        case 'X':
            printf(" flag '-X'\n");
            printf("The marker loci are on the X-chromosome.\n");
            printf("All marker loci are on the X-chromosome.\n");
            break;
        case 'w':
            printf(" flag 'wt'\n");
            printf("Weight multiple sibpair pair within the same nuclear family by 2/A,\n");
            printf("where A is the number of affected sibs. Default is for such comparisons\n");
            printf("to be weighted.");
            break;
        case 'v':
            printf(" flag 'v#'\n");
            printf("Specifies variance estimator to be used for score tests.\n");
            printf("Values of # are n, t, or e for native, theoretical and\n");
            printf("empirical estimates respectively. More than one -v option can be\n");
            printf("specified, (Default -vn)");
            break;
        default:
            break;
        }
    press_return();
}

/*----------------------splink_options checking-------------*/
static char *splink_option_check(char *option)

{
    int n;
    char *value;

    value= CALLOC((size_t) FILENAME_LENGTH, char);
    value[0]='\0';
    switch (option[0]) {
    case 'a':
    case 'c':
    case 'h':
    case 'p':
    case 'X':
    case 'm':
    case 'w':
        sprintf(value, " ");
    break;
    case 'b':
    case 'r':
    case 'd':
    case 'f':
    case 'i':
    case 'l':
    case 'n':
    case 's':
    case 'S':
    case 'x':
    case 'v':
        if (option[0] == 'r' || option[0] == 'b' || option[0] == 'S')
            option++;

    if (strlen(option) < 1) {
        printf("No numerical value provided with option %s, please specify > ",
               option);
        enter_number(&n); newline;
        option++;
        sprintf(option, "%d", n);
        sprintf(value, "%d", n);
    }
    else
        sprintf(value, "%s", option);
    break;

    case 'o':
        if (strlen(option) < 2) {
            do {
                printf("\nEnter file NAME to which the IBD assignment for each\n");
                printf("family will be written: > ");
                fcmap(stdin, "%s", value);
            }  while (value[0] == '\0');
            option++;
            sprintf(option, "%s", value);
        }
        else
            sprintf(value, "%s", option); // ??
        break;

    default:
        warn_unknown(option);
        break;
    }
    return &(value[0]);

}


/*--------------------------------------------------------------------*/
static void            CreateSplinkCShell(char *options, int single_loci)
{
    char            choice[101], *value;
    int             done_;
    char            tmp[10], *opt_ptr;
    int             done_parsing, j;
    char            *expression, option_string[100] = {""};

    done_ = 0;

    if (!DEFAULT_OPTIONS) {
        while(!done_) {
            draw_line();
            printf("SPLINK option selection:\n");
            printf("Enter one or more of the option(s) listed below, separated by spaces\n");
            printf("     ( e.g. '-asp -c -mf -S+232 -s23' )\n OR\n");
            printf("Enter 'help' followed by option for a description (e.g. 'help -asp')\n");
            printf("For no options press <RETURN> or <Enter>\n");

            printf("  -asp            -n#\n");
            printf("  -bs#            -o#\n");
            printf("  -c              -pt\n");
            printf("  -d#             -rs#\n");
            printf("  -f#             -s#\n");

            if (!single_loci)  {
                printf("  -h              -S+#\n");
                printf("  -i#             -x#\n");
                printf("  -X              -v#\n");
                printf("  -mf             -wt\n");
            } else  {
                printf("  -h              -x#\n");
                printf("  -i#             -v\n");
                printf("  -X              -wt#\n");
                printf("  -mf\n");
            }

            printf("Enter options as a string of text (at most 100 characters):\n");
            printf("Terminate string by pressing <RETURN> or <Entry> key > ");
            strcpy(choice, "");
            fflush(stdout);
            IgnoreValue(fgets(choice, 100, stdin)); newline;

            if (*choice == 0) break;

            if (strncasecmp("help", choice, (size_t) 4) == 0) {
                /* print the help right here */
                expression = &(choice[4]);
                sscanf(expression, "%s", tmp);
                if (!(tmp[0] == '-' || isalpha((int)tmp[0])))   {
                    printf("No option specified for help!\n");
                } else {
                    if (tmp[0] == '-') splink_help(&tmp[1]);
                    else     	 splink_help(&tmp[0]);
                }
            } else {
                expression = &(choice[0]);
                done_parsing=0; opt_ptr=&option_string[0];
                while(!done_parsing) {
                    while (*expression == ' ') expression++;
                    if (*expression == '-') {
                        expression++;
                        j=0;
                        while (*expression != ' ' && *expression != '\n') {
                            tmp[j]=*expression++; j++;
                        }
                        tmp[j]= '\0';
                        value=splink_option_check(tmp);
                        if (value != NULL) {
                            sprintf(opt_ptr, "-%s ", tmp);
                            opt_ptr+= strlen(tmp)+2;
                        }
                    }
                    if (*expression == '\n')   {
                        done_ = 1;
                        done_parsing = 1;
                        *opt_ptr = '\0';
                    }
                }
            }
        }
    }
    printf("\nSPLINK options selected --> %s\n", option_string);
    strcpy(options, option_string);
    return;
}

#ifdef Top
#undef Top
#endif
#define Top (*LPTop)

/*--------------------------------------------------------------------*/
void            create_SPLINK(linkage_ped_top **LPTop,
			      analysis_type * analysis, int *numchr,
			      char *file_names[],
			      int untyped_ped_opt)
{

#define DEFAULT_SPLINK_ped_opt 2
    int             a1, a2;
    char            *options;
    int             tr, nloop, *trp, num_affec= num_traits, i, m, kk,k;
    linkage_locus_top *LTop;
    linkage_ped_top *Top2; /* *Copy; */
    linkage_ped_rec *tpe;
    FILE           *fp, *cfp;
    char            fl_stat[12], cchoice[10], syscmd[100];
    int             select, choice;
    char            pedfl_name[2*FILENAME_LENGTH+MAX_NAMELEN];
    char            cshfl_name[2*FILENAME_LENGTH];
    char            pfl[FILENAME_LENGTH], fformat[6], pformat[6];
    int             niter, locus_cnt, wid;
    char           stars[3] = " *";
    int            loop=1;

    LTop = Top->LocusTop;
    options = MALLOC_STR(200, char);
    select = DEFAULT_SPLINK_ped_opt;

    if (! DEFAULT_OPTIONS) {
        while (loop) {
            draw_line();
            printf("SPLINK file conversion options:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf("%c 1) Single pedigree file for all markers.\n", stars[0]);
            printf("%c 2) One pedigree file for each marker.\n", stars[1]);
            printf("Enter selection: 1 or 2 to select, 0 to proceed > ");
            fcmap(stdin, "%s", cchoice); newline;
            sscanf(cchoice, "%d", &loop);
            if (loop == 0) {
                ;
            } else if (loop == 1) {
                select=loop;
                strcpy(stars, "* ");
            } else if (loop == 2) {
                select=loop;
                strcpy(stars, " *");
            } else {
                warn_unknown(cchoice);
            }
        }
    }
    if (select==1) {
        if (main_chromocnt > 1)
            change_output_chr(file_names[0], 0);
        else
            change_output_chr(file_names[0], *numchr);
    } else
        change_output_chr(file_names[0], -9);

    if (main_chromocnt > 1) change_output_chr(file_names[3], 0);
    else change_output_chr(file_names[3], *numchr);

    if (!DEFAULT_OUTFILES) {
        choice=-1;

        while (choice != 0)  {
            draw_line();
            print_outfile_mssg();
            printf("SPLINK-format file name selection menu:\n");
            draw_line();
            printf("0) Done with menu - please proceed.\n");
            if (select == 1) {
                printf(" 1) Pedigree file name:         %-15s\t%s\n", file_names[0],
                       file_status(file_names[0], fl_stat));
            } else {
                printf(" 1) Pedigree file name stem:    %-15s\n", file_names[0]);
            }
            printf(" 2) C-shell file name :         %-15s\t%s\n", file_names[3],
                   file_status(file_names[3], fl_stat));
            printf("Enter selection:  0, 1 or 2 > ");
            fcmap(stdin, "%s", cchoice);  newline;
            sscanf(cchoice, "%d", &choice);
            test_modified(choice);

            if (choice == 1)     {
                printf("\nEnter NEW pedigree file name: > "); newline;
                fcmap(stdin, "%s", file_names[0]);
            }
            if (choice == 2)     {
                printf("\nEnter NEW C-shell file name: > "); newline;
                fcmap(stdin, "%s", file_names[3]);
            }
            if (!(choice==0 || choice==1 || choice==2)) {
                warn_unknown(cchoice);
            }
        }
    }

    /* omit pedigrees not typed */
    omit_peds(untyped_ped_opt, Top);
    Top2 = new_lpedtop();
    makenucs1(Top, NULL, Top2);
    create_nuked_fids(Top2);
    CreateSplinkCShell(options, 0);

    NLOOP;
    trp = &(global_trait_entries[0]);

    ped_name_width(Top2, &wid);
    sprintf(fformat, "%%%ds ", wid);
    person_id_width(Top2, &wid);
    sprintf(pformat, "%%%dd ", wid);

    get_loci_on_chromosome(0);
    locus_cnt = ((LoopOverTrait == 1)?
                 NumChrLoci - num_affec : NumChrLoci - num_affec + 1);

    if (select == 1) {
        for (tr=0; tr <= nloop; tr++) {
            if (tr == 0 && nloop > 1) continue;
            /* 1. pedigree file */
            sprintf(pedfl_name, "%s/%s", output_paths[tr], file_names[0]);
            if ((fp = fopen(pedfl_name, "w")) == NULL)  {
                errorvf("Unable to open file %s\n", pedfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            /* set the affection status correctly */
            /*      affected_by_status(Top2, *trp, "2", analysis); */
            for (m = 0; m < Top2->PedCnt; m++)  {
                for (i = 0; i < Top2->Ped[m].EntryCnt; i++)   {
                    tpe = &(Top2->Ped[m].Entry[i]);
                    fprintf(fp, fformat, Top2->Ped[m].Name);
                    fprintf(fp, pformat, tpe->ID);
                    fprintf(fp, pformat, tpe->Father);
                    fprintf(fp, pformat, tpe->Mother);
                    fprintf(fp, "%1d ", tpe->Sex);
                    /* write the affection status locus */
                    /* write_affection_data(fp, *trp, &(Top2->LocusTop->Locus[*trp]), tpe); */
                    fprintf(fp, " %1d", aff_status_entry(tpe->Pheno[*trp].Affection.Status,
                                                         tpe->Pheno[*trp].Affection.Class,
                                                         &(Top2->LocusTop->Locus[*trp])));
                    /* write the genotypes */
                    for (k = 0; k < NumChrLoci; k++) {
                        linkage_locus_rec *locus = &Top2->LocusTop->Locus[ChrLoci[k]];
                        kk=ChrLoci[k];
                        switch (locus->Type) {
                        case AFFECTION:
                            break;
                        case NUMBERED:
                            get_2alleles(tpe->Marker, kk, &a1, &a2);
                            // NOTE: The second allele is left alligned...
                            fprintf(fp, " %2s/%-2s",
                                    format_allele(locus, a1), format_allele(locus, a2));
                            break;
                        case QUANT:
                        case BINARY:
                            errorf("Unsupported locus type in SPLINK files.");
                            EXIT(DATA_TYPE_ERROR);
                            break;
                        default:
                            errorvf("unknown locus type (locus %d) while writing output\n",
                                    kk + 1);
                            EXIT(DATA_TYPE_ERROR);
                        }
                    }
                    fprintf(fp, "\n");
                }
            }
            fclose(fp);
            /* 2. the shell file */
            sprintf(cshfl_name, "%s/%s", output_paths[tr], file_names[3]);
            if ((cfp = fopen(cshfl_name, "w")) == NULL) {
                errorvf("Unable to open file %s.\n", cshfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(cfp, "#!/bin/csh -f\n");
            script_time_stamp(cfp);
            fprintf(cfp, "touch splink.lst\n");
            fprintf(cfp, "date >> splink.lst\n");

            fprintf(cfp, "splink < %s %s -l%d", file_names[0], options, locus_cnt);
            sprintf(syscmd, "chmod +x %s", cshfl_name);
            fclose(cfp);
            System(syscmd);

            if (nloop == 1) break;
            trp++;
        }

        /*    draw_line(); */
        create_mssg(*analysis);
        sprintf(err_msg, "    Pedigree file:       %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "    C-shell script file: %s", file_names[3]);
        mssgf(err_msg);
    } else {
        /* if separate pedigree files */
        niter=1;
        create_mssg(*analysis);
        for (k = 0; k < NumChrLoci; k++)   {
            linkage_locus_rec *locus = &LTop->Locus[ChrLoci[k]];
            kk=ChrLoci[k];
            if ((LTop->Locus[kk].Type == NUMBERED) ||
                (LTop->Locus[kk].Type == BINARY)) {
                sprintf(pfl, "%s.%s", file_names[0], locus->LocusName);
                trp = &(global_trait_entries[0]);
                for (tr = 0; tr <= nloop; tr++) {
                    if (nloop > 1 && tr==0) continue;
                    sprintf(pedfl_name, "%s/%s", output_paths[tr], pfl);
                    if ((fp = fopen(pedfl_name, "w")) == NULL)   {
                        errorvf("Unable to open file '%s'\n", pedfl_name);
                        EXIT(FILE_WRITE_ERROR);
                    }
                    for (m = 0; m < Top2->PedCnt; m++) {
                        for (i = 0; i < Top2->Ped[m].EntryCnt; i++)  {
                            tpe = &(Top2->Ped[m].Entry[i]);
                            fprintf(fp, fformat, Top2->Ped[m].Name);
                            fprintf(fp, pformat, tpe->ID);
                            fprintf(fp, pformat, tpe->Father);
                            fprintf(fp, pformat, tpe->Mother);
                            fprintf(fp, "%1d ", tpe->Sex);
                            /*   write_affection_data(fp, 0,
                                 &(Top2->LocusTop->Locus[*trp]), tpe); */
                            fprintf(fp, " %1d", tpe->Pheno[*trp].Affection.Status);
                            get_2alleles(tpe->Marker, kk, &a1, &a2);
                            // NOTE: The second allele is left alligned...
                            fprintf(fp, " %2s/%-2s\n",
                                    format_allele(locus, a1), format_allele(locus, a2));
                        }
                    }
                    fclose(fp);
                    /* shell  file */
                    sprintf(cshfl_name, "%s/%s", output_paths[tr], file_names[3]);
                    if (niter == 1)	{
                        if ((cfp = fopen(cshfl_name, "w")) == NULL) {
                            errorvf("Unable to open file %s.\n", cshfl_name);
                            EXIT(FILE_WRITE_ERROR);
                        }
                        fprintf(cfp, "#!/bin/csh -f\n");
                        script_time_stamp(cfp);
                        fprintf(cfp, "touch splink.lst\n");
                        fprintf(cfp, "date >> splink.lst\n");
                    } else {
                        if ((cfp = fopen(cshfl_name, "a")) == NULL) {
                            errorvf("Unable to open file %s.\n", cshfl_name);
                            EXIT(FILE_WRITE_ERROR);
                        }
                    }
                    fprintf(cfp, "echo Running SPLINK on the marker %s\n", locus->LocusName);
                    fprintf(cfp, "splink -l1 -i%s -S+ %s < %s > /dev/null \n",
                            locus->LocusName, options, pfl);
                    if (niter == locus_cnt) {
                        fprintf(cfp, "echo SPLINK analyses are completed.\n");
                        fprintf(cfp,
                                "echo Please see splink.lst for a summary of the results.\n");
                    }
                    fclose(cfp);
                    sprintf(syscmd, "chmod +x %s", cshfl_name);
                    System(syscmd);
                    if (nloop == 1) break;
                    trp++;
                }
                sprintf(err_msg, "    Pedigree file:       %s", pfl);
                mssgf(err_msg);
                niter++;

            }
        } /* every marker */
        sprintf(err_msg, "    C-shell script file: %s",
                file_names[3]);
        mssgf(err_msg);
    }
}

#undef Top
