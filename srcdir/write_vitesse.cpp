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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "typedefs.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
        batch_input_ext.h:  batchf
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file haldane_theta kosambi_theta log_marker_selections sub_prog_name
            linkage_ext.h:  get_loci_on_chromosome switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT chmod_X_file draw_line global_ou_files script_time_stamp
        write_files_ext.h:  save_linkage_peds
*/


void create_vitesse_files(linkage_ped_top **LTop, int *numchr,
			  char *file_names[], int untyped_ped_opt);


char vitesse_prog_names[2][8] = {"LINKMAP", "MLINK  "};

static void vitesse_program_option(int *prog_num)

{
    int exit_loop=0;
    int i;
    char sub_prog[20], stars[3] = "* ";

    if (batchANALYSIS) {
        *prog_num = TO_VITESSE->_suboption;
        return;
    }
    *prog_num = 1;

    while (!exit_loop) {
        draw_line();
        printf("Vitesse program selection menu:\n");
        printf("0) Done with this menu - please proceed.\n");
        printf("%c1) %s\n", stars[0], vitesse_prog_names[0]);
        printf("%c2) %s\n", stars[1], vitesse_prog_names[1]);

        printf("Enter selection 0, 1 or 2 > ");

        fcmap(stdin, "%d", &i); newline;
        if (i == 1 || i == 2) {
            *prog_num = i;
            if ( i == 1) {
                strcpy(stars, "* ");
            } else {
                strcpy(stars, " *");
            }
        } else if (i == 0) {
            exit_loop = 1;
        } else {
            printf("Unknown option %d. \n", *prog_num);
        }
    }
    TO_VITESSE->_suboption = *prog_num;
    TO_VITESSE->sub_prog_name(*prog_num, sub_prog);
    strcpy(Mega2BatchItems[/* 6 */ Analysis_Sub_Option].value.name, sub_prog);
    batchf(Analysis_Sub_Option);
}

static void vitesse_file_names(int single_file, char *file_names[])

{
    char stem[5], file_type[15], *fname, fl_stat[12];
    int select=-1;

    change_output_chr(file_names[0], -9);
    change_output_chr(file_names[1], -9);

    if (single_file==1) {
        change_output_chr(file_names[3], 0);
        change_output_chr(file_names[4], 0);
        change_output_chr(file_names[5], 0);
    }

    if (single_file == 0 && main_chromocnt > 1) {
        change_output_chr(file_names[3], -9);
        change_output_chr(file_names[4], -9);
        change_output_chr(file_names[5], -9);
    }

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        select =0;
    }

    while(select != 0) {
        strcpy(stem, "");
        draw_line();
        print_outfile_mssg();
        printf("Vitesse output file name selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Pedigree file name stem:        %-15s\n", file_names[0]);
        printf(" 2) Data file name stem:            %-15s\n", file_names[1]);
        if (single_file == 1 || main_chromocnt == 1) {
            printf(" 3) C-shell script file name:       %-15s\t%s\n",
                   file_names[3], file_status(file_names[3], fl_stat));
            printf(" 4) Vitesse output file name:       %-15s\t%s\n",
                   file_names[4], file_status(file_names[4], fl_stat));
            printf(" 5) Stream file name:               %-15s\t%s\n",
                   file_names[5], file_status(file_names[5], fl_stat));
        } else {
            strcpy(stem, "stem");
            printf(" 3) C-shell script file name stem:  %-15s\n", file_names[3]);
            printf(" 4) Vitesse output file name stem:  %-15s\n", file_names[4]);
            printf(" 5) Stream file name stem:          %-15s\n", file_names[5]);
        }

        printf("Select option 0-5 > ");
        fcmap(stdin, "%d", &select); newline;
        test_modified(select);
        fname = NULL; // silly compiler does not understand switch stmnts
        switch(select) {
        case 0:
            break;
        case 1:
            strcpy(file_type, "pedigree"); fname=file_names[0];
            break;
        case 2:
            strcpy(file_type, "data"); fname=file_names[1];
            break;
        case 3:
            strcpy(file_type, "C-shell script"); fname=file_names[3];
            break;
        case 4:
            strcpy(file_type, "vitesse output"); fname=file_names[4];
            break;
        case 5:
            strcpy(file_type, "stream"); fname=file_names[5];
            break;
        default:
            printf("Unknown option %d\n", select);
            break;
        }
        if (select <= 5 && select >=1) {
            printf("Enter new %s file name %s > ", file_type, stem);
            fcmap(stdin, "%s", fname);
        }
    }
    if (single_file == 0 && main_chromocnt > 1) {
        change_output_chr(file_names[3], 1);
        strcat(file_names[3], ".sh");
    }
}

static void vitesse_linkmap_marker_interval(int *mint_opt, int *intv_size,
					    int **markers,
					    linkage_locus_top *LTop)

{
    int i, select=-1;
    char selected[2], cselect[10];
    int *chr_markers = NULL, num_chr_markers=0, found_invalid=0; // too hard for compiler to decide chr_markers is assigned

    if (main_chromocnt == 1) {
        chr_markers = CALLOC((size_t) num_reordered, int);
        for (i=0; i<num_reordered; i++) {
            if (LTop->Locus[reordered_marker_loci[i]].Type == NUMBERED ||
                LTop->Locus[reordered_marker_loci[i]].Type == BINARY) {
                chr_markers[num_chr_markers]=reordered_marker_loci[i];
                num_chr_markers++;
            }
        }
    }
    *intv_size=2;
    *mint_opt=1;
    *markers=NULL;

    if (DEFAULT_OPTIONS) {
        mssgf("LINKMAP marker interval set to default value");
        mssgf("as requested in the batch file.");
        select=0;
    }

    selected[0]='*'; selected[1]=' ';
    draw_line();
    while (select != 0) {
        printf("Linkmap marker interval selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        printf("%c1) Analyze all intervals [size %d].\n", selected[0],
               *intv_size);
        if (main_chromocnt == 1) {
            printf("%c2) Specify an interval [", selected[1]);
            if (*markers != NULL) {
                for(i=0; i < *intv_size; i++) {
                    printf(" %d", (*markers)[i]);
                }
            }
            printf("]\n");
            printf("Enter option 0-2 > ");
        } else {
            printf("Enter option 0 or 1 > ");
        }
        fcmap(stdin, "%s", cselect); newline;
        sscanf(cselect, "%d", &select);
        switch(select) {
        case 0:
            break;
        case 1:
            printf("Enter new interval size > ");
            fcmap(stdin, "%d", intv_size); newline;
            if (main_chromocnt == 1) {
                if (num_chr_markers < *intv_size) {
                    printf("There are only %d markers, please select again.\n",
                           num_chr_markers);
                    draw_line();
                    *intv_size=2;
                }
            }
            selected[0]='*'; selected[1]=' ';
            *mint_opt=1;
            break;

        case 2:
            if (main_chromocnt == 1) {
                printf("Markers available for selection:\n");
                for(i=0; i < num_chr_markers; i++) {
                    printf(" %d) %s\n",
                           i+1, LTop->Marker[chr_markers[i]].MarkerName);
                }
                printf("Enter size of new interval > ");
                fcmap(stdin, "%d", intv_size); newline;
                if (*markers != NULL) free(*markers);
                *markers =  CALLOC((size_t) *intv_size, int);
                printf("Enter %d loci > ", *intv_size);
                for (i=0; i< *intv_size; i++) {
                    fcmap(stdin, "%d", &((*markers)[i]));
                }
                /* check selection */
                found_invalid=0;
                for (i=0; i< *intv_size; i++) {
                    if (!((*markers)[i] >=1 && (*markers)[i] <= num_chr_markers)) {
                        printf("Invalid selection %d, please select again.\n",
                               (*markers)[i]);
                        found_invalid=1;
                        break;
                    }
                }
                newline;
                if (found_invalid) {
                    free(*markers);
                    *markers=NULL; *intv_size=2;
                    selected[0]='*'; selected[1]=' ';
                    *mint_opt=1;
                } else {
                    selected[1]='*'; selected[0]=' ';
                    *mint_opt=2;
                }
                break;
            }
        default:
            warn_unknown(cselect);
            select=-1;
            break;
        }
    }
    if (*mint_opt == 2) {
        for (i=0; i< *intv_size; i++) {
            (*markers)[i] = chr_markers[(*markers)[i]-1];
        }
    } else {
        draw_line();
        mssgf("For analysis Vitesse LINKMAP:");
        sprintf(err_msg, "Option : Analyze all intervals of size %d.", *intv_size);
        mssgf(err_msg);
    }
    return;
}

static void verify_intv_size(int *intv_size, int num_loci, int numchr)

{
    int choice = -1;
    char cselect[10];


    while (choice) {
        if (*intv_size < 1 || *intv_size > num_loci) {
            printf("There are only %d marker loci on chr %d.\n",
                   num_loci, numchr);
            *intv_size=((*intv_size < 1)? 2: num_loci);
            printf("Enter new interval size for chr %d:\n", numchr);
            draw_line();
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) Linkmap interval size : %d \n", *intv_size);
            printf("Enter 0 or 1 > ");
            fflush(stdout);
            IgnoreValue(fgets(cselect, 9, stdin)); newline;
            choice=atoi(cselect);
            if (choice  != 0 &&  choice  != 1) {
                warn_unknown(cselect);
            } else if (choice == 1) {
                printf("Enter new interval size > ");
                fcmap(stdin, "%d", intv_size);
            }
        } else {
            choice=0;
        }
    }
    return;
}



static void vitesse_linkmap_params(double *stop_theta, int *num_evals)

{
    int select=-1;

    *stop_theta=0.46;
    *num_evals=5;

    if (DEFAULT_OPTIONS) {
        mssgf("Linkmap parameters set to default values");
        mssgf("as requested in the batch file.");
        select =0;
    }

    while (select != 0) {
        draw_line();
        printf("LINKMAP parameter modification menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        /*    printf(" 1) Stop value for theta:      %f\n", *stop_theta); */
        printf(" 1) Number of evaluations:     %d\n", *num_evals);
        printf("Enter 0 or 1 > ");
        fcmap(stdin, "%d", &select); newline;
        switch(select) {
        case 1:
            /*      printf("Enter stop value for theta > ");
                    fcmap(stdin, "%f", stop_theta);
                    if (*stop_theta > 0.5) {
                    printf("WARNING: Stop value of theta cannot be greater than 0.5\n");
                    *stop_theta=0.5;
                    }
                    newline;
                    break;
                    case 2:
            */
            printf("Enter number of evaluations > ");
            fcmap(stdin, "%d", num_evals);
            newline;
            break;
        case 0:
            break;
        default:
            printf("Unknown option %d\n", select);
            break;
        }
    }
    return;
}

static void vitesse_mlink_params(double *start_theta, double *increment,
                                 double *stop_theta, int *num_evals)
{
    /* sex difference, theta to vary, Increment in theta, Stop value,
       num_addition_evals */

    int select=-1;

    *start_theta=0.0;
    *increment = 0.05;
    *stop_theta = 0.46;
    *num_evals=0;

    if (DEFAULT_OPTIONS) {
        mssgf("MLINK parameters set to default values");
        mssgf("as requested in the batch file.");
        select =0;
    }

    while (select != 0) {
        draw_line();
        printf("MLINK parameter menu:\n");
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Start value of theta     [%5.4f]\n", *start_theta);
        printf(" 2) Increment in theta       [%5.4f]\n", *increment);
        printf(" 3) Stop value of theta      [%5.4f]\n", *stop_theta);

        printf("Enter option 0 - 3 > ");
        fcmap(stdin, "%d", &select);
        switch(select) {
        case 1:
            printf("Enter start value of theta > ");
            fcmap(stdin, "%g", start_theta); newline;
            break;
        case 2:
            printf("Enter increment in theta > ");
            fcmap(stdin, "%g", increment); newline;
            break;
        case 3:
            printf("Enter stop value of theta > ");
            fcmap(stdin, "%g", stop_theta); newline;
            break;
        case 0:
            break;
        default:
            printf("Unknown option %d\n", select);
            break;
        }
    }

    return;
}

static void assign_mlink_vals(double start_theta, double increment, double stop_theta,
                              int num_evals, linkage_locus_top *loc_top)

{

    loc_top->Program = MLINK;
    loc_top->SexDiff=0;
    loc_top->Interference=0;
    loc_top->Run.mlink.increment=increment;
    loc_top->Run.mlink.stop_theta=stop_theta;
    loc_top->Run.mlink.start_theta=start_theta;
    loc_top->Run.mlink.num_evals=num_evals;

    return;
}


static void assign_linkmap_vals(double stop_theta, int num_evals,
                                int num_thetas, linkage_locus_top *loc_top)

{

    loc_top->Program = LINKMAP;
    loc_top->SexDiff=0;
    loc_top->Interference=0;
    loc_top->Run.linkmap.trait_marker=1;
    loc_top->Run.linkmap.num_evals=num_evals;
    loc_top->Run.linkmap.stop_theta=stop_theta;
    loc_top->Run.linkmap.recomb_frac= CALLOC((size_t) num_thetas, double);

    return;

}

static void write_vitesse_lspstm(FILE *fstm, char *file_names[],
                                 linkage_locus_top *LTop, int locus_num,
                                 int *loci, int program, int *loci1)

{

    int i;
/*  extern   time_t  time();
    extern   struct  tm  *localtime(); */
#ifdef HIDEDATE
    extern   time_t  NOTIMEval;
#endif
    time_t   time_val;
    struct   tm     *time_stt;
    char     time_str[40];
    char     time_elements[5][20];

/*    char     *month[12] = */
/*    { */
/*      "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec" */
/*    }; */

    fprintf(fstm, "@\n");
    fprintf(fstm, "LSP\n");
    fprintf(fstm, "%s\n", vitesse_prog_names[program-1]);
    fprintf(fstm, "%s\n",
#ifdef HIDEPATH
            NOPATH
#else
            file_names[0]
#endif
        );
    fprintf(fstm, "%s\n",
#ifdef HIDEPATH
            NOPATH
#else
            file_names[1]
#endif
        );
    fprintf(fstm, "%s\n", file_names[2]);
    fprintf(fstm, "%s\n", file_names[3]);
    fprintf(fstm, "lsp.log\n");
    fprintf(fstm, "lsp.stm\n");
    /* instead of write_time()
       stuff copied from mega2_0/write_vitesse.c on watson */
#ifdef HIDEDATE
    time_val = NOTIMEval;
#else
    time(&time_val);
#endif
    time_stt = localtime (&time_val);
    sprintf(time_str, "%s", asctime(time_stt));
    sscanf(time_str, "%s %s %s %s %s",  time_elements[0], time_elements[1],
           time_elements[2], time_elements[3],  time_elements[4]);

    /* now format the time_str */
    fprintf(fstm, "%s-%s-%c%c %s\n", time_elements[2],
            time_elements[1], time_elements[4][2],
            time_elements[4][3], time_elements[3]);

    /* fprintf(fstm, "%2d-%s-%d  %02d:%02d:%02d\n",
       time_stt->tm_mday,month[time_stt->tm_mon],time_stt->tm_year,
       time_stt->tm_hour,time_stt->tm_min,time_stt->tm_sec);
    */
    /*------------------------*/
    fprintf(fstm, "%d\n", locus_num);
    for (i=0; i<locus_num; i++) {
        if (loci[i] == -9)  fprintf(fstm, " %d", LTop->Run.linkmap.trait_marker+1);
        else {
            fprintf(fstm, " %d", loci1[i]);
        }
    }
    fprintf(fstm, "\n");
    fprintf(fstm, "0\n");
    if (LTop->Program==LINKMAP) {
        fprintf(fstm, "%d\n", LTop->Run.linkmap.trait_marker+1);
        fprintf(fstm, "%d\n", LTop->Run.linkmap.num_evals);
        fprintf(fstm, "%9.8f\n", LTop->Run.linkmap.stop_theta);
    } else {
        fprintf(fstm, "1\n");
        fprintf(fstm, "%9.8f\n", LTop->Run.mlink.increment);
        fprintf(fstm, "%9.8f\n", LTop->Run.mlink.stop_theta);

    }
    fprintf(fstm, "0\n");

    if (LTop->Program == LINKMAP) {
        for (i=0; i < locus_num-1; i++)
            fprintf(fstm, " %9.8f", LTop->Run.linkmap.recomb_frac[i]);
        fprintf(fstm, "\n~\n");
    } else {
        fprintf(fstm, " %9.8f",  LTop->Run.mlink.start_theta);
        fprintf(fstm, "\n0\n~\n");
    }


    return;

}

static void write_vitesse_lsplog(FILE *flog, char *file_names[],
                                 linkage_locus_top *LTop, int locus_num,
                                 int *loci, int program, int run_num, int *loci1)

{

    int i;
#ifdef __INTEL__
    long now;
#else
    time_t now;
#endif

    /*  char *ctime(); */

    time (&now);
    fprintf(flog, "****************************************************************\n\n");
    fprintf(flog, "                      %s RUN %d\n\n",
            vitesse_prog_names[program-1], run_num);
    fprintf(flog, "     Pedigree file              : %s\n",
#ifdef HIDEPATH
            NOPATH
#else
            file_names[0]
#endif
        );
    fprintf(flog, "     Parameter file             : %s\n",
#ifdef HIDEPATH
            NOPATH
#else
            file_names[1]
#endif
        );
    fprintf(flog, "     Output Pedigree file       : %s\n", file_names[2]);
    fprintf(flog, "     Output Parameter file      : %s\n", file_names[3]);
    fprintf(flog, "     Log file                   : lsp.log\n");
    fprintf(flog, "     Stream file                : lsp.stm\n");
    fprintf(flog, "     Date Run                   : %s",
#ifdef HIDEDATE
            NOcTIME
#else
            ctime(&now)
#endif
        );
    fprintf(flog, "     Sex Difference             : 0\n");
    if (LTop->Program == LINKMAP) {
        fprintf(flog, "     Test Locus                 : 1\n");
        fprintf(flog, "     Number of iterations       : %d\n", LTop->Run.linkmap.num_evals);
    } else {
        fprintf(flog, "     Recomb. Fraction to vary   : 1\n");
        fprintf(flog, "     Increment Value            : %9.8f\n", LTop->Run.mlink.increment);
        fprintf(flog, "     Stop Value                 : %9.8f\n", LTop->Run.mlink.stop_theta);
    }
    fprintf(flog, "     Locus Order                :");
    /* shift the marker loci so that they are numbered consecutively with
       the trait locus - this is WRONG, changed  back to original locus numbers*/

    for (i=0; i< locus_num;i++) {
        if (loci[i]==-9) {
            fprintf(flog, " %d", LTop->Run.linkmap.trait_marker +1);
        } else {
            fprintf(flog, " %d", loci1[i]);
            /*LTop->Locus[loci[i]].LocusName); */
        }
    }
    fprintf(flog, "\n");
    fprintf(flog, "     Locus Names                :");
    for (i=0; i< locus_num;i++) {
        if (loci[i]==-9) fprintf(flog, " %s",
                                 LTop->Pheno[LTop->Run.linkmap.trait_marker].TraitName);
        else fprintf(flog, " %s", LTop->Locus[loci[i]].LocusName);
    }
    fprintf(flog, "\n");
    if (LTop->Program == LINKMAP) {
        fprintf(flog, "     Male recomb. fractions     :");
        for (i=0; i<locus_num-1; i++) {
            fprintf(flog, " %9.8f", LTop->Run.linkmap.recomb_frac[i]);
        }
        fprintf(flog, "\n");
    } else {
        fprintf(flog, "     Initial recomb. fraction   : %9.8f",
                LTop->Run.mlink.start_theta);
    }

    fprintf(flog, "\n");
    fprintf(flog, "****************************************************************\n\n");

    return;

}

static void write_vitesse_locus_file(char *loutfl_name,
				     linkage_locus_top *LTop,
				     int locus_num, int *loci,
				     int sex_linked)

{
    /* this is the same as write_gh_locus_file, except it takes a vector of selected loci */

    int tr, nloop, num_affec=num_traits, locus, allele, tmpi, tmpi2;
    int *trp, num_loci, trait_locus = 0; //compiler: not obvious
    linkage_locus_rec *Locus;
    char fl[2*FILENAME_LENGTH];
    FILE *filep;

    NLOOP;

    for (tmpi=0; tmpi < locus_num; tmpi++) {
        if (loci[tmpi] == -9) {
            trait_locus = tmpi; break;
        }
    }

    trp=&(global_trait_entries[0]);

    for (tr=0; tr<= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(fl, "%s/%s", output_paths[tr], loutfl_name);
        if ((filep=fopen(fl, "w")) == NULL) {
            errorvf("Could not open %s for writing.\n", fl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(filep, "%d %d %d %d\n%d %.3f %.3f %d\n",
                locus_num, LTop->RiskLocus, sex_linked,
                LTop->Program, LTop->MutLocus, LTop->MutMale,
                LTop->MutFemale, LTop->Haplotype);

        for (tmpi=0; tmpi < locus_num; tmpi++)
            fprintf(filep, "%d ", tmpi+1);

        fputc('\n', filep);
        num_loci = locus_num;

        for (locus = 0; locus <num_loci; locus++) {
            if (locus == trait_locus)	Locus=&(LTop->Locus[*trp]);
            else Locus=&(LTop->Locus[loci[locus]]);
            fprintf(filep, "%d %d", Locus->Type - 1, Locus->AlleleCnt);
            if (Locus->LocusName != NULL)
                fprintf(filep, " #%s", Locus->LocusName);
            fputc('\n', filep);
            for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
            }
            fputc('\n', filep);
            switch(Locus->Type) {
            case QUANT:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Quant.ClassCnt);
                for (tmpi = Locus->Pheno->Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                    for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                        fprintf(filep," %.6f", Locus->Pheno->Props.Quant.Mean[tmpi-1][tmpi2]);
                fprintf(filep, " << GENOTYPE MEANS\n");
                /* ASSUMES only one trait per qtl */
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Variance[0][0]);
                fprintf(filep," %.6f\n", Locus->Pheno->Props.Quant.Multiplier);
                break;
            case AFFECTION:
                fprintf(filep, "%d\n", Locus->Pheno->Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < Locus->Pheno->Props.Affection.ClassCnt; tmpi++) {
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
                break;
            case BINARY:
                fprintf(filep, "%d\n", Locus->Marker->Props.Binary.FactorCnt);
                for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++) {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(filep, " %d", (int) Locus->Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', filep);
                }
                break;
            case NUMBERED:
                /* nothing to do  */
                break;
            default:
                errorvf("Unknown locus type %d\n", (int) LTop->Locus[locus].Type);
                EXIT(DATA_TYPE_ERROR);
                break;
            }
        }
        /* now for the recombination information */
        fprintf(filep, "%d %d\n", LTop->SexDiff, LTop->Interference);
        if (LTop->Program == MLINK) {
            fprintf(filep, "%5.4f\n", LTop->Run.mlink.start_theta);
            fprintf(filep, "1 %5.4f %5.4f\n", LTop->Run.mlink.increment,
                    LTop->Run.mlink.stop_theta);
        } else {
            for (tmpi=0; tmpi < num_loci-1; tmpi++)
                fprintf(filep, "%9.8f ", LTop->Run.linkmap.recomb_frac[tmpi]);
            fprintf(filep, "\n");
            fprintf(filep, "%d %9.8f %d \n", trait_locus+1,
                    LTop->Run.linkmap.stop_theta,
                    LTop->Run.linkmap.num_evals);
        }
        trp++;
        fclose(filep);
        if (nloop == 1) break;
    }
    return;

}

static void write_vitesse_header(FILE *fshell, char *streamfl_name)

{

    fprintf(fshell, "#!/bin/csh -f\n");
    fprintf(fshell, "#\n");
    script_time_stamp(fshell);
    fprintf(fshell, "\n");
    fprintf(fshell, "unalias rm\n");
    fprintf(fshell, "unalias cp\n");
    return;
}

static void write_vitesse_order_run(char *shell, int *perm,
				    linkage_ped_top *Top,
				    int num_markers, int program,
				    char *outfl_name, char *loutfl_name,
				    char *finalfl_name, char *streamfl_name,
				    int run_num, int first_time, int *perm1)

{
    /* num_markers is the length of perm */

    int i, tr, *trp, nloop, num_affec=num_traits;
    FILE *fshell, *fstm, *flog;
    char md[2], *TmpFile, *filenames[4], flname[3*FILENAME_LENGTH];

    NLOOP;
    trp = &(global_trait_entries[0]);
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
/*      for (i=0; i<num_markers; i++) { */
/*        if (perm[i] == -9) { */
/*  	perm[i] = *trp; break; */
/*        } */
/*      } */
        Top->LocusTop->Run.linkmap.trait_marker = *trp;
        sprintf(flname, "%s/%s", output_paths[tr], shell);
        if (first_time) {
            strcpy(md, "w"); delete_file(flname);
        }
        else  strcpy(md, "a");

        if ((fshell=fopen(flname, md)) == NULL) {
            errorvf("Could not open shell file %s.\n", flname);
            EXIT(FILE_WRITE_ERROR);
        }
        if (first_time)  write_vitesse_header(fshell, streamfl_name);

        fprintf(fshell, "echo\necho\n");
        /* echo to stdout */
        fprintf(fshell,
                "echo \"*********************************************************\"\n");
        fprintf(fshell, "echo \"Run   %d - analysis type %s :\" \n",
                run_num, vitesse_prog_names[program-1]);
        fprintf(fshell, "echo Loci: ");
        for (i=0; i<num_markers; i++) {
            if (perm[i] == -9) {
                fprintf(fshell, "%s ", Top->LocusTop->Pheno[*trp].TraitName);
            } else {
                fprintf(fshell, "%s ", Top->LocusTop->Locus[perm[i]].LocusName);
            }
        }
        fprintf(fshell, "\n");
        fprintf(fshell,
                "echo \"*********************************************************\"\n");
        fprintf(fshell, "echo\n");

        /* reproduce the contents of log and stm files */
        TmpFile= CALLOC((size_t) 3*FILENAME_LENGTH, char);
        sprintf(TmpFile, "%s/STM.%s", output_paths[tr], &(outfl_name[7]));
        fstm=fopen(TmpFile, "w");

        filenames[0]=mega2_input_files[0];
        filenames[1]=mega2_input_files[1];
        filenames[2]=outfl_name;
        filenames[3]=loutfl_name;

        write_vitesse_lspstm(fstm, filenames, Top->LocusTop, num_markers, perm,
                             program, perm1);
        fclose(fstm);
        /* create these files */
        if (run_num==1)
            fprintf(fshell, "\t cat STM.%s >! %s\n", &(outfl_name[7]),
                    streamfl_name);
        else
            fprintf(fshell, "\t cat STM.%s >> %s\n", &(outfl_name[7]),
                    streamfl_name);
        sprintf(TmpFile, "%s/LOG.%s", output_paths[tr], &(outfl_name[7]));
        flog=fopen(TmpFile, "w");
        write_vitesse_lsplog(flog, filenames, Top->LocusTop, num_markers, perm,
                             program, run_num, perm1);
        fclose(flog);
        if (run_num == 1)
            fprintf(fshell, "\t cat LOG.%s >! %s\n", &(outfl_name[7]),
                    finalfl_name);
        else
            fprintf(fshell, "\t cat LOG.%s >> %s\n", &(outfl_name[7]),
                    finalfl_name);
        free(TmpFile);
        /*----------------*/
        fprintf(fshell, "\tif (-e lsp.log) then\n\t\trm -f lsp.log\n\tendif\n");
        fprintf(fshell, "\tif (-e vstream.dat) then\n\t\trm -f vstream.dat\n\tendif\n");
        fprintf(fshell, "\tcp %s pedfile.dat\n", outfl_name);
        fprintf(fshell, "\tcp %s datafile.dat\n", loutfl_name);

        fprintf(fshell, "\tvitesse");
        if (ITEM_READ(Value_Missing_Quant_On_Output))
           fprintf(fshell," -U %s", Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        fprintf(fshell, "\n");

        fprintf(fshell, "\t\tcat voutfile.dat >> %s\n", finalfl_name);
        fprintf(fshell, "\t\tif !(-e vstream.dat) then\n");
        fprintf(fshell,
                "\t\t\techo VITESSE failed to run on %s and %s\n",
                outfl_name, loutfl_name);
        fprintf(fshell, "\t\t\techo Please check your files for errors\n");
        fprintf(fshell,
                "\t\t\techo Note VITESSE may also fail to handle pedigrees with loops.\n");
        fprintf(fshell, "\t\texit (1)\n");
        fprintf(fshell, "\t\tendif\n");
        fprintf(fshell, "\tcat vstream.dat >> %s\n",
		streamfl_name);
        if (ITEM_READ(Value_Missing_Quant_On_Output)) {
            fprintf(fshell, "\techo WARNING: The missing value was set to '%s' via the -U command line option.\n",
		    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
            fprintf(fshell, "\techo WARNING: Using the batch file item 'Value_Missing_Quant_On_Output'.\n");
        }
        fclose(fshell);
        chmod_X_file(flname);
        trp++;
        if (nloop == 1) break;
    }
    sprintf(err_msg, "         Stream file:         STM.%s",
            &(outfl_name[7]));
    mssgf(err_msg);
    sprintf(err_msg, "         Log file:            LOG.%s",
            &(outfl_name[7]));
    mssgf(err_msg);
    return;
}

#ifdef LTop
#undef LTop
#endif
#define LTop (*Top)

void create_vitesse_files(linkage_ped_top **Top, int *numchr,
			  char *file_names[], int untyped_ped_opt)
{
    int i, j, k, jj,  single_file, run_num;
    char chr_str[3];
    int mint_opt = 0, intv_size, intv_size_width, *markers = NULL, tmp_intv_size; //compiler: mint_opt,markers are always
    int first_time=1, num_evals;                                                  //assigned in vitesse_linkmap_marker_interval()
    int prog_num, *chr_markers;
    double start_theta, stop_theta, increment, th;
    char *outfl_name, *loutfl_name;
    char *file_ext;
    int *perms, *perms1, num_loci, xlinked;
    analysis_type analysis;

    analysis = TO_VITESSE;

    draw_line();
    printf("Analysis option: Vitesse\n");
    draw_line();

    if (main_chromocnt > 1)
        single_file=global_ou_files(1);
    else
        single_file=0;
//SL
    if ((*Top)->LocusTop->SexLinked == 2 && single_file == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }

    vitesse_program_option(&prog_num);
    vitesse_file_names(single_file, file_names);

    if (prog_num == 1) {
        /* linkmap */
        vitesse_linkmap_marker_interval(&mint_opt, &intv_size, &markers,
                                        LTop->LocusTop);
        vitesse_linkmap_params(&stop_theta, &num_evals);
        if (mint_opt == 2) {
            draw_line();
            mssgf("Selected the following markers for Vitesse LINKMAP:");
            log_marker_selections(LTop, markers, intv_size, NULL, (void (*)(const char *)) mssgf, 0);
        }
    } else {
        vitesse_mlink_params(&start_theta, &increment,  &stop_theta, &num_evals);
    }

    first_time = 1;

    for (i=0; i< main_chromocnt; i++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[i];
            if (!single_file) {
                first_time=1;
                change_output_chr(file_names[3], *numchr);
                change_output_chr(file_names[4], *numchr);
                change_output_chr(file_names[5], *numchr);
            }
        }

        get_loci_on_chromosome(*numchr);
        xlinked = (((LTop->LocusTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
                    (LTop->LocusTop->SexLinked == 1))? 1 : 0);

        if (prog_num == 1)
            assign_linkmap_vals(stop_theta, num_evals, intv_size, LTop->LocusTop);
        else
            assign_mlink_vals(start_theta, increment, stop_theta, num_evals, LTop->LocusTop);

        /* omit untyped pedigrees */
        omit_peds(untyped_ped_opt, LTop);
        chr_markers= CALLOC((size_t) NumChrLoci, int);
        k=0;
        CHR_STR(*numchr, chr_str);
        for (j = 0; j<NumChrLoci; j++) {
            if (LTop->LocusTop->Locus[ChrLoci[j]].Type == NUMBERED ||
                LTop->LocusTop->Locus[ChrLoci[j]].Type == BINARY) {
                chr_markers[k]=ChrLoci[j];
                k++;
            }
        }
        num_loci=k;
        
	    if (InputMode == BATCH_FILE_INPUTMODE && num_loci == 0) {
            errorvf("No markers are selected. A Vitesse Linkmap cannot be run without having any markers selected.\n");
            EXIT(DATA_INCONSISTENCY);
	    }
        
        if (i == 0) create_mssg(TO_VITESSE);
        else
            draw_line();

        /* for each program type */
        if (prog_num == 1) {
            /* LINKMAP */
            if (mint_opt == 2) {
                run_num=1;    	/* user-specified interval */
                intv_size_width = (int) floor(log10((double)intv_size)) + 1;
                perms= CALLOC((size_t) intv_size+1, int);
                perms1= CALLOC((size_t) intv_size+1, int);
                for (j=0; j<=intv_size; j++) {
                    perms[j]=-9; /* this will be filled in with the disease locus */
                    perms1[j]=-9; /* this will be filled in with the disease locus */
                    if (j > 0)
                        for (k=0; k<j; k++) {
                            perms[k]=markers[k];    perms1[k]=k+1;
                        }

                    if (j < intv_size)
                        for (k=j; k<intv_size; k++) {
                            perms[k+1]=markers[k]; perms1[k+1]=k+1;
                        }

                    /* assign recombination fractions and stop_thetas*/
                    /* a. trait is at position 0 */
                    if (j==0) {
                        LTop->LocusTop->Run.linkmap.recomb_frac[0]=0.5;
                        for (k=1; k<intv_size; k++) {
                            th=LTop->LocusTop->Marker[markers[k]].pos_avg -
                                LTop->LocusTop->Marker[markers[k-1]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[k]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.stop_theta=0.0;
                    }
                    /* b. trait is at position intv_size */
                    if (j==intv_size) {
                        LTop->LocusTop->Run.linkmap.recomb_frac[intv_size-1]=0.0;
                        for (k=0; k<intv_size-1; k++) {
                            th=LTop->LocusTop->Marker[markers[k+1]].pos_avg -
                                LTop->LocusTop->Marker[markers[k]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[k]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.stop_theta=0.5;
                    }
                    /* c. trait is betwwen two  marker loci at position j */
                    if (j>0 && j<intv_size) {
                        for (k=1; k<j-1; k++) {
                            th=LTop->LocusTop->Marker[markers[k+1]].pos_avg -
                                LTop->LocusTop->Marker[markers[k]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[k]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.recomb_frac[j]=0.0;
                        for (k=j; k<intv_size; k++) {
                            th=LTop->LocusTop->Marker[markers[k]].pos_avg -
                                LTop->LocusTop->Marker[markers[k-1]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[k]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }

                        th=LTop->LocusTop->Marker[markers[j]].pos_avg -
                            LTop->LocusTop->Marker[markers[j-1]].pos_avg;
                        LTop->LocusTop->Run.linkmap.stop_theta=
                            ((LTop->LocusTop->map_distance_type == 'h')?
                             haldane_theta(th) : kosambi_theta(th));

                    }

                    file_ext = CALLOC(strlen(LTop->LocusTop->Marker[markers[0]].MarkerName) +
                                      strlen(LTop->LocusTop->Marker[markers[intv_size-1]].MarkerName) +
                                      3 + intv_size_width, char);

                    sprintf(file_ext, "%s-%s_%d",
                            LTop->LocusTop->Marker[markers[0]].MarkerName,
                            LTop->LocusTop->Marker[markers[intv_size-1]].MarkerName, j+1);
                    outfl_name = CALLOC(strlen(file_ext)+strlen(file_names[0]) +
                                        strlen(chr_str)+5, char);
                    loutfl_name = CALLOC(strlen(file_ext)+strlen(file_names[1]) +
                                         strlen(chr_str)+5, char);
                    sprintf(outfl_name,  "%s.%s.%s",  file_names[0], file_ext, chr_str);
                    sprintf(loutfl_name, "%s.%s.%s",  file_names[1], file_ext, chr_str);

                    save_linkage_peds(outfl_name, LTop, analysis, intv_size+1, perms);
                    write_vitesse_locus_file(loutfl_name, LTop->LocusTop,
                                             intv_size+1, perms, xlinked);
/* 	  write_linkage_locfile_inorder(loutfl_name, LTop->LocusTop,  */
/* 					intv_size+1, perms, xlinked, */
/* 					TO_VITESSE); */

                    write_vitesse_order_run(file_names[3], perms, LTop, intv_size+1,
                                            prog_num, outfl_name, loutfl_name,
                                            file_names[4], file_names[5],
                                            run_num, first_time, perms1);

                    first_time=0;
                    sprintf(err_msg, "         Pedigree file:       %s", outfl_name);
                    mssgf(err_msg);
                    sprintf(err_msg, "         Locus file:          %s", loutfl_name);
                    mssgf(err_msg);
                    free(file_ext); free(outfl_name); free(loutfl_name);
                    run_num++;
                }
                free(perms); perms = NULL;
                free(perms1); perms1 = NULL;
                continue;
            } /* user-specified interval */

            /* all intervals, sliding window for each marker less interval_size */
            tmp_intv_size=intv_size;
            verify_intv_size(&intv_size, num_loci, *numchr);
            run_num=1;
            intv_size_width = (int) floor(log10((double)intv_size)) + 1;
            markers= CALLOC((size_t) intv_size, int);
            perms =  CALLOC((size_t) intv_size+1, int);
            perms1 =  CALLOC((size_t) intv_size+1, int);
            for (j=0; j< (num_loci-intv_size+1); j++) {
                /* first create the window */
                for (k = 0; k < intv_size; k++)  {
                    markers[k]=chr_markers[j+k];
                }
                /* place disease locus at different intervals */
                for (k=0; k <= intv_size; k++) {
                    perms[k]=-9; 	  perms1[k]=-9;
                    if (k > 0) {
                        for (jj=0; jj<k; jj++) {
                            perms[jj]=markers[jj];   perms1[jj]=j+jj+2;
                        }
                    }

                    if (k < intv_size) {
                        for (jj=k; jj < intv_size; jj++) {
                            perms[jj+1]=markers[jj]; 	      perms1[jj+1]=j+jj+2;
                        }
                    }

                    /* assign recombination fractions and stop_thetas*/
                    /* a. trait is at position 0 */
                    if (k==0) {
                        LTop->LocusTop->Run.linkmap.recomb_frac[0]=0.5;
                        for (jj=1; jj<intv_size; jj++) {
                            th=LTop->LocusTop->Marker[markers[jj]].pos_avg -
                                LTop->LocusTop->Marker[markers[jj-1]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[jj]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.stop_theta=0.0;
                    }
                    /* b. trait is at position intv_size */
                    if (k==intv_size) {
                        LTop->LocusTop->Run.linkmap.recomb_frac[intv_size-1]=0.0;
                        for (jj=0; jj<intv_size-1; jj++) {
                            th=LTop->LocusTop->Marker[markers[jj+1]].pos_avg -
                                LTop->LocusTop->Marker[markers[jj]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[jj]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.stop_theta=0.5;
                    }

                    /* c. trait is betwwen two  marker loci at position j */
                    if (k>0 && k<intv_size) {
                        for (jj=0; jj<k-1; jj++) {
                            th=LTop->LocusTop->Marker[markers[jj+1]].pos_avg -
                                LTop->LocusTop->Marker[markers[jj]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[jj]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }
                        LTop->LocusTop->Run.linkmap.recomb_frac[k-1]=0.0;
                        for (jj=k; jj<intv_size; jj++) {
                            th=LTop->LocusTop->Marker[markers[jj]].pos_avg -
                                LTop->LocusTop->Marker[markers[jj-1]].pos_avg;
                            LTop->LocusTop->Run.linkmap.recomb_frac[jj]=
                                ((LTop->LocusTop->map_distance_type == 'h')?
                                 haldane_theta(th) : kosambi_theta(th));
                        }

                        th=LTop->LocusTop->Marker[markers[k]].pos_avg -
                            LTop->LocusTop->Marker[markers[k-1]].pos_avg;
                        LTop->LocusTop->Run.linkmap.stop_theta=
                            ((LTop->LocusTop->map_distance_type == 'h')?
                             haldane_theta(th) : kosambi_theta(th));

                    }
                    file_ext = CALLOC(strlen(LTop->LocusTop->Marker[markers[0]].MarkerName) +
                                      strlen(LTop->LocusTop->Marker[markers[intv_size-1]].MarkerName) +
                                      3 + intv_size_width, char);
                    sprintf(file_ext, "%s-%s_%d",
                            LTop->LocusTop->Marker[markers[0]].MarkerName,
                            LTop->LocusTop->Marker[markers[intv_size-1]].MarkerName, k+1);
                    outfl_name = CALLOC(strlen(file_ext)+strlen(file_names[0]) +
                                        strlen(chr_str)+5, char);
                    loutfl_name = CALLOC(strlen(file_ext)+strlen(file_names[1]) +
                                         strlen(chr_str)+5, char);

                    sprintf(outfl_name,  "%s.%s.%s",  file_names[0], file_ext, chr_str);
                    sprintf(loutfl_name, "%s.%s.%s",  file_names[1], file_ext, chr_str);

                    save_linkage_peds(outfl_name, LTop, analysis, intv_size+1, perms);
                    write_vitesse_locus_file(loutfl_name, LTop->LocusTop,
                                             intv_size+1, perms, xlinked);

                    write_vitesse_order_run(file_names[3], perms, LTop, intv_size+1,
                                            prog_num, outfl_name, loutfl_name,
                                            file_names[4], file_names[5],
                                            run_num, first_time, perms1);
                    first_time=0;
                    run_num++;
                    sprintf(err_msg, "         Pedigree file:       %s", outfl_name);
                    mssgf(err_msg);
                    sprintf(err_msg, "         Locus file:          %s", loutfl_name);
                    mssgf(err_msg);
                    free(file_ext); free(outfl_name); free(loutfl_name);
                }
            }
            free(markers);
            free(perms); perms = NULL;
            free(perms1); perms1 = NULL;
            /* END LINKMAP- option 1 */
            intv_size=tmp_intv_size;
        } else {
            /* MLINK */
            run_num=1;
            perms= CALLOC((size_t) 2, int);
            perms1= CALLOC((size_t) 2, int);
            perms[0]=-9;       perms1[0]=-9;
            for (j=0; j < num_loci; j++) {
                perms[1]=chr_markers[j];
                perms1[1]=j+2;
                outfl_name = CALLOC(strlen(file_names[0]) +
                                    strlen(LTop->LocusTop->Locus[perms[1]].LocusName)+5, char);
                loutfl_name = CALLOC(strlen(file_names[1]) +
                                     strlen(LTop->LocusTop->Locus[perms[1]].LocusName)+5, char);

                sprintf(outfl_name, "%s.%s", file_names[0], LTop->LocusTop->Locus[perms[1]].LocusName);
                sprintf(loutfl_name, "%s.%s", file_names[1], LTop->LocusTop->Locus[perms[1]].LocusName);

                save_linkage_peds(outfl_name, LTop, analysis, 2, perms);
                write_vitesse_locus_file(loutfl_name, LTop->LocusTop, 2, perms, xlinked);

                sprintf(err_msg, "         Pedigree file:       %s", outfl_name);
                mssgf(err_msg);
                sprintf(err_msg, "         Locus file:          %s", loutfl_name);
                mssgf(err_msg);
                write_vitesse_order_run(file_names[3], perms, LTop, 2, prog_num,
                                        outfl_name, loutfl_name, file_names[4],
                                        file_names[5], run_num, first_time, perms1);
                first_time=0;
                run_num++;
                free(outfl_name); free(loutfl_name);
            }
            free(perms); perms = NULL;
            free(perms1); perms1 = NULL;
        } /* MLINK */

/*      if (prog_num==1) */
/*        free(LTop->LocusTop->Run.linkmap.recomb_frac); */
        if (main_chromocnt > 1) {
            if (single_file != 1) {
                sprintf(err_msg, "         C-shell file:        %s", file_names[3]);
                mssgf(err_msg);
            }
        }
        free(chr_markers); chr_markers=NULL;
    }
    if (single_file || main_chromocnt == 1) {
        sprintf(err_msg, "         C-shell file:        %s", file_names[3]);
        mssgf(err_msg);
    }
}

#undef LTop
