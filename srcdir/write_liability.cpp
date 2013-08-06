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

#include "common.h"
#include "typedefs.h"

#include "liable.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  create_genos safe_divide
  output_file_names_ext.h:  CHR_STR shorten_path
         user_input_ext.h:  test_modified
*/


char **get_liablfl_name(int *numchr, int tab_text);
void get_output_options(int *ou_percent, int *empty_cols,
			int *empty_rows, int *tab_text);
void write_liable_dist(liable_allele_dist *allele_dist, char **liablefl_name,
		       linkage_ped_top *LPedTreeTop, int num_AFFECTION,
		       int num_NUMBERED, int ou_percent,
		       int empty_rows, int empty_cols, int tab_text);


/* function write liability distribution
   File can be 1) Tab-text for export into Excel
   2) Aligned and prettified for visual inspection

*/


char **get_liablfl_name(int *numchr, int tab_text)

{

    /* if there are more than one chormosome, all the results will be compiled into
       one single groups file named "groups.all" by default, otherwise this file will
       have the chromosome number suffix */

    char **liablefl_name, num[4];
    int exit_loop, user_reply;
    int i;
    char flname[2*FILENAME_LENGTH];

    liablefl_name = CALLOC((size_t) 3, char *);
    liablefl_name[0] = CALLOC((size_t) 2*FILENAME_LENGTH, char);
    liablefl_name[1] = CALLOC((size_t) 2*FILENAME_LENGTH, char);
    liablefl_name[2] = CALLOC((size_t) 2*FILENAME_LENGTH, char);

    CHR_STR(*numchr, num);

    if (main_chromocnt > 1)
        strcpy(liablefl_name[0], "groups.all");
    else
        sprintf(liablefl_name[0], "groups.%s", num);

    strcpy(liablefl_name[1], "alleles");
    strcpy(liablefl_name[2], "genotypes");

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults");
    } else {
        exit_loop=0;
        while (!exit_loop) {
            printf("\nCounts file name change options:\n");
            printf("0) Done with this menu - please proceed\n");
            printf(" 1) Group counts summary file name :      %s\t%s\n",
                   liablefl_name[0],
                   ((access(liablefl_name[0], F_OK) == 0)? "[Overwrite]" : "[New]"));
            if (tab_text == 1) {
                printf(" 2) Allele counts file name stem :        %s\n", liablefl_name[1]);
                printf(" 3) Genotype counts file name stem :      %s\n", liablefl_name[2]);
                printf("Enter selection: 0-3 > ");
            } else {
                printf("\nEnter selection: 0 or 1 > ");
            }
            fcmap(stdin, "%d", &user_reply); newline;
            test_modified(user_reply);

            switch (user_reply) {
            case 0:
                exit_loop = 1;
                break;
            case 1:
                printf("\nEnter NEW group counts summary file name: > ");
                fcmap(stdin, "%s", liablefl_name[0]); newline;
                break;
            case 2:
                if (tab_text ==1 ) {
                    printf("\nEnter NEW allele counts file name stem: > ");
                    fcmap(stdin, "%s", liablefl_name[1]); newline;
                    break;
                }
            case 3:
                if (tab_text == 1) {
                    printf("\nEnter NEW genotype counts file name stem: > ");
                    fcmap(stdin, "%s", liablefl_name[2]); newline;
                    break;
                }
            default:
                printf("\nUnknown option %d.\n", user_reply);
                break;
            }
        }
    }
    for (i=0; i < 3; i++) {
        sprintf(flname, "%s/%s", output_paths[0], liablefl_name[i]);
        strcpy(liablefl_name[i], flname);
    }
    return liablefl_name;

}

void get_output_options(int *ou_percent, int *empty_cols, int *empty_rows, int *tab_text)

{
    int user_reply, exit_loop =0;

    if (DEFAULT_OPTIONS) {
        mssgf("Liability summary output options set to defaults.");
        mssgf("as requested in the batch file");
        return;
    }
    while (!exit_loop) {
        printf("\nSummary file output options :\n");
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Display percentages in addition to counts [%s]\n",
               ((*ou_percent) ? "yes" : "no"));
        printf(" 2) Include rows with all zero counts [%s]\n",
               ((*empty_rows) ? "yes" : "no"));
        printf(" 3) Include columns with all zero counts [%s]\n",
               ((*empty_cols) ? "yes" : "no"));
        printf(" 4) Output tab-text file [%s]\n",
               ((*tab_text) ? "yes" : "no"));
        printf("\nEnter selection (0 to proceed, item number to toggle setting) > ");
        fcmap(stdin, "%d", &user_reply); newline;
        switch (user_reply) {
        case 0:
            exit_loop=1; break;
        case 1:
            *ou_percent =((*ou_percent)? 0 : 1); break;
        case 2:
            *empty_rows =((*empty_rows)? 0: 1); break;
        case 3:
            *empty_cols =((*empty_cols)? 0: 1); break;
        case 4:
            *tab_text = ((*tab_text)? 0: 1); break;
        default:
            printf("Unknown option %d, enter a number between 0 and 4\n", user_reply);
            break;
        }
    }
    return;
}

static void create_header_p(int empty_cols, int *header_p,
			    liable_allele_dist allele_dist, int nn, int lclasses)
{
    int mm, ll;

    for (mm=0; mm<lclasses; mm++){
        ll=3*mm;
        if (!empty_cols) {
            if (allele_dist.MARKERTOTAL(nn,mm,0)>0) { header_p[ll] = 1; }
            else { header_p[ll]=0;}
            if (allele_dist.MARKERTOTAL(nn,mm,1)>0) { header_p[ll+1] = 1;}
            else { header_p[ll+1]=0; }
            if (allele_dist.MARKERTOTAL(nn,mm,2)>0) { header_p[ll+2] = 1;}
            else { header_p[ll+2]=0;}
        } else {
            header_p[ll]=1;
            header_p[ll+1]=1;
            header_p[ll+2]=1;
        }
    }
}

static char *create_class_headers(int lclasses, int *header_p,
				  int ou_percent, int tab_text)

{
    char *header_str, *header_pos, lb_st[17];
    int i, l, num_cols=0;
    size_t pos;
    int col_width, len_header=0;

    header_str=NULL;
    if (tab_text) {
        if (lclasses<=1) {
            /* Possible error: casting size_t to int */
            len_header=header_p[0]*(int) strlen("\tUnknown") +
                header_p[1]*(int) strlen("\tNormal") +
                header_p[2]*(int) strlen("\tAffected");
            if (ou_percent)
                len_header += (int) strlen("\t\t\t");
            header_str=CALLOC((size_t) len_header+2, char);
            pos=0;
            header_pos=&header_str[pos];
            if (header_p[0]) {
                strncpy(header_pos, "\tUnknown", strlen("\tUnknown"));
                pos += strlen("\tUnknown");
                header_pos = &header_str[pos];
                if (ou_percent) strncpy(header_pos, "\t", (size_t) 1);
                pos += ou_percent * strlen("\t");
            }
            header_pos=&header_str[pos];

            if (header_p[1]) {
                strncpy(header_pos, "\tNormal", strlen("\tNormal"));
                pos += strlen("\tNormal");
                header_pos=&header_str[pos];
                if (ou_percent) strncpy(header_pos, "\t", (size_t) 1);
                pos += ou_percent;
            }
            header_pos=&header_str[pos];

            if (header_p[2]){
                strncpy(header_pos, "\tAffected", strlen("\tAffected"));
                pos += strlen("\tAffected");
                header_pos=&header_str[pos];
                if (ou_percent) strncpy(header_pos, "\t", (size_t)  1);
                pos += ou_percent;
                header_pos = &(header_str[pos]);
                /* don't put a tab after the last column */
            }
            strcpy(header_pos, "");
            return header_str;
        } /* lclasses <=1 */

        i=0;
        while (i< 3*lclasses) {
            num_cols += header_p[i]+header_p[i+1]+header_p[i+2];
            i += 3;
        }


        col_width=15;
        len_header = num_cols*col_width;
        pos=0;

        header_str=CALLOC((size_t) len_header+10, char);
        if (num_cols>48) printf("Warning: file is at least 50 columns wide\n");

        i=0; l=1; /* liability class */
        while (i < 3*lclasses) {
            if (header_p[i]) {
                header_pos = &header_str[pos];
                if (ou_percent)
                    sprintf(lb_st, "\tUnk(%1d)\t", l);
                else
                    sprintf(lb_st, "\tUnk(%1d)", l);
                strncpy(header_pos, lb_st, (size_t)  col_width);
                pos += strlen(lb_st);
            }
            if (header_p[i+1]) {
                header_pos = &header_str[pos];
                if (ou_percent)
                    sprintf(lb_st, "\tNorm(%1d)\t", l);
                else
                    sprintf(lb_st, "\tNorm(%1d)", l);
                strncpy(header_pos, lb_st, (size_t)  col_width);
                pos +=  strlen(lb_st);
            }
            if (header_p[i+2]) {
                header_pos = &header_str[pos];
                if (ou_percent)
                    sprintf(lb_st, "\tAff(%1d)\t", l);
                else
                    sprintf(lb_st, "\tAff(%1d)", l);
                strncpy(header_pos, lb_st, (size_t)  col_width);
                pos +=  strlen(lb_st);
            }
            i += 3; l++;
        } /* lclasses > 1 */

        header_str[pos] = '\0';
        return header_str;
    } /* if (tab_text) */
    else {
        /* can be justified up to 5 digits of counts */
        col_width=16;
        if (lclasses<=1){
            len_header= col_width*(header_p[0] + header_p[1] + header_p[2]);
            header_str=CALLOC((size_t) len_header+1, char);
            pos=0; header_pos = &(header_str[pos]);
            if (header_p[0]) {
                strncpy(header_pos, "    Unknown     ", (size_t)  col_width);
                pos += col_width;
            }
            if (header_p[1]) {
                header_pos = &(header_str[pos]);
                strncpy(header_pos, "     Normal     ", (size_t) col_width);
                pos += col_width;
            }
            if (header_p[2]) {
                header_pos = &(header_str[pos]);
                strncpy(header_pos, "    Affected    ", (size_t) col_width);
                pos += col_width;
            }
        } /* lclasses <=1 */
        else {
            i=0;
            while (i< 3*lclasses) {
                len_header += header_p[i]*col_width +
                    header_p[i+1]*col_width + header_p[i+2]*col_width;
                i += 3;
            }
            pos=0;
            header_str=CALLOC((size_t) len_header+10, char);
            i = 0; l = 1;
            while (i < 3*lclasses) {
                if (header_p[i]) {
                    header_pos = &(header_str[pos]);
                    sprintf(lb_st,"  (Unknown,%1d)   ", l);
                    strncpy(header_pos, lb_st, (size_t) col_width);
                    pos += strlen(lb_st);
                }
                if (header_p[i+1]) {
                    header_pos = &(header_str[pos]);
                    sprintf(lb_st,"  (Normal,%1d)    ", l);
                    strncpy(header_pos, lb_st, (size_t) col_width);
                    pos += strlen(lb_st);
                }
                if (header_p[i+2]) {
                    header_pos = &(header_str[pos]);
                    sprintf(lb_st, "  (Affected,%1d)  ", l);
                    strncpy(header_pos, lb_st, (size_t) col_width);
                    pos += strlen(lb_st);
                }
                l++; i += 3;
            }
        } /* lclasses > 1 */
        header_str[pos] = '\0';

        return &(header_str[0]);
    } /* if (!tab_text) */
}


void write_liable_dist(liable_allele_dist *allele_dist, char **liablefl_name,
		       linkage_ped_top *LPedTreeTop, int num_AFFECTION,
		       int num_NUMBERED, int ou_percent,
		       int empty_rows, int empty_cols, int tab_text)

{
    char *disease_name, *header, *ou_file_name;
    char *file_ext;
    int *loci_NUMBERED, lclasses, locus_AFFECTION, num_alleles;
    int num_genos, **genotypes, *header_p;
    FILE *fp;
    linkage_locus_top *ltop;
    int kk, nn, mm, ll, pos;
    char shortname[FILENAME_LENGTH];

    /* for indicating whether to write a column */

    if (!tab_text) {
        fp = fopen(liablefl_name[0], "a");
    }

    for (kk=0; kk<num_AFFECTION; kk++){
        locus_AFFECTION=allele_dist[kk].locus_AFFECTION;
        ltop=LPedTreeTop->LocusTop;
        lclasses=ltop->Locus[locus_AFFECTION].Data.Affection.ClassCnt;
        header_p=NULL;
        header_p=CALLOC((size_t) lclasses*3, int);

        disease_name = ltop->Locus[locus_AFFECTION].Name;
        loci_NUMBERED=allele_dist[kk].loci_NUMBERED;
        for (nn=0; nn<num_NUMBERED; nn++){
            if (tab_text) {
                file_ext= CALLOC(strlen(ltop->Locus[loci_NUMBERED[nn]].Name)+
                                 strlen(disease_name)+3, char);
                sprintf(file_ext, "_%s.%s", disease_name, ltop->Locus[loci_NUMBERED[nn]].Name);
                ou_file_name =
                    CALLOC((strlen(liablefl_name[1])+strlen(file_ext)+1), char);
                strncpy(ou_file_name, liablefl_name[1], strlen(liablefl_name[1]));
                strcat(ou_file_name, file_ext);
                fp = fopen(ou_file_name, "w");
            }
/*       marker_str = ltop->Locus[loci_NUMBERED[nn]].Name; */

            num_alleles=ltop->Locus[loci_NUMBERED[nn]].AlleleCnt;
            create_header_p(empty_cols, header_p, allele_dist[kk], nn, lclasses);
            header=create_class_headers(lclasses, header_p, ou_percent, tab_text);

            /* for each marker, write the allele_data first, then the genotype_data */
            if (tab_text) {
                fprintf(fp, "%s%sTotal", ltop->Locus[loci_NUMBERED[nn]].Name, header);
                if (ou_percent) fprintf(fp, "\tTotal Freq");
            } else {
                fprintf(fp, "Trait %s\n", disease_name);
                fprintf(fp, "      Allele counts for marker  %s\n",
                        ltop->Locus[loci_NUMBERED[nn]].Name);
                fprintf(fp, "Allele   %s Total    ", header);
                if (ou_percent) fprintf(fp, "Total Freq");
            }
            fprintf(fp, "\n");
            for (ll=0; ll<num_alleles; ll++){
                if (!empty_rows || (allele_dist[kk].ALLELETOTAL(nn,ll)>0)) {
                    if (tab_text) fprintf(fp, "%5d\t", ll+1);
                    else  fprintf(fp, "%5d    ", ll+1);
                    for (mm=0;mm<lclasses;mm++){
                        pos=3*mm;
                        if (tab_text) {

                            if (header_p[pos]){
                                fprintf(fp, "%5d\t", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 0));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].STATUSPERCENT(nn,ll, mm, 0));
                            }
                            if (header_p[pos+1]){
                                fprintf(fp, "%5d\t", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 1));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].STATUSPERCENT(nn,ll, mm, 1));
                            }
                            if (header_p[pos+2]){
                                fprintf(fp, "%5d\t", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 2));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].STATUSPERCENT(nn,ll, mm, 2));
                            }
                        } /* tab_text */
                        else {
                            if (header_p[pos]){
                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f ",
                                            allele_dist[kk].LIABLESTATUS(nn,ll, mm, 0),
                                            allele_dist[kk].STATUSPERCENT(nn,ll, mm, 0));
                                else
                                    fprintf(fp, "     %5d     ", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 0));
                            }
                            if (header_p[pos+1]){
                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f ",
                                            allele_dist[kk].LIABLESTATUS(nn,ll, mm, 1),
                                            allele_dist[kk].STATUSPERCENT(nn,ll, mm, 1));
                                else
                                    fprintf(fp, "     %5d     ", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 1));
                            }
                            if (header_p[pos+2]){
                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f ",
                                            allele_dist[kk].LIABLESTATUS(nn,ll, mm, 2),
                                            allele_dist[kk].STATUSPERCENT(nn,ll, mm, 2));
                                else
                                    fprintf(fp, "     %5d     ", allele_dist[kk].LIABLESTATUS(nn,ll, mm, 2));
                            }
                        }  /* not tab_text */
                    }
                    if (tab_text) {
                        fprintf(fp, "%5d\t", allele_dist[kk].ALLELETOTAL(nn,ll));
                        if (ou_percent)
                            fprintf(fp, "%5.4f\t",
                                    safe_divide(allele_dist[kk].ALLELETOTAL(nn,ll),
                                                allele_dist[kk].TOTAL(nn)));
                    } else {
                        if (ou_percent)
                            fprintf(fp, "%5d     %5.4f",
                                    allele_dist[kk].ALLELETOTAL(nn,ll),
                                    safe_divide(allele_dist[kk].ALLELETOTAL(nn,ll),
                                                allele_dist[kk].TOTAL(nn)));
                        else
                            fprintf(fp, "%5d", allele_dist[kk].ALLELETOTAL(nn,ll));
                    }
                    fprintf(fp,"\n");
                } /* Row is not all zeros or empty_rows should be written */
            } /* each allele */
            /* write the column totals */
            if (tab_text)       fprintf(fp, "TOTAL\t");
            else fprintf(fp, "TOTAL    ");
            for (mm=0; mm<lclasses; mm++){
                ll=3*mm;
                if (tab_text) {
                    if (header_p[ll]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERTOTAL(nn, mm,0));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                    if (header_p[ll+1]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERTOTAL(nn,mm,1));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                    if (header_p[ll+2]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERTOTAL(nn,mm,2));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                } else {
                    if (header_p[ll]) {
                        if (ou_percent)
                            fprintf(fp, "%5d          ", allele_dist[kk].MARKERTOTAL(nn, mm,0));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERTOTAL(nn, mm,0));

                    }
                    if (header_p[ll+1]) {
                        if (ou_percent)
                            fprintf(fp, "%5d          ", allele_dist[kk].MARKERTOTAL(nn,mm,1));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERTOTAL(nn,mm,1));

                    }
                    if (header_p[ll+2]) {
                        if (ou_percent)
                            fprintf(fp, "%5d          ", allele_dist[kk].MARKERTOTAL(nn,mm,2));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERTOTAL(nn,mm,2));
                    }
                }
            }
            fprintf(fp, "%5d\n", allele_dist[kk].TOTAL(nn));

            if (tab_text) {
                fclose(fp); 	/* close the marker file */
                shorten_path(ou_file_name, shortname);
                sprintf(err_msg, "      Allele count summary file:    %s", shortname);
                mssgf(err_msg);
                free(ou_file_name); /* free output file name */
                ou_file_name = NULL;

                /* reallocate output file name */
                ou_file_name =
                    CALLOC((strlen(liablefl_name[2])+strlen(file_ext)+1), char);
                strncpy(ou_file_name, liablefl_name[2], strlen(liablefl_name[2]));
                strcat(ou_file_name, file_ext);
                /* open the genotype file for writing */
                fp = fopen(ou_file_name, "w");
            }
            /* write the genotypes */

            num_genos=NUMGENOS(num_alleles);  genotypes=create_genos(num_alleles);
            if (tab_text) {
                fprintf(fp, "%s%sTotal", ltop->Locus[loci_NUMBERED[nn]].Name, header);
                if (ou_percent) fprintf(fp, "\tTotal Freq");
            } else {
                fprintf(fp, "\n      Genotype counts for marker  %s\n",
                        ltop->Locus[loci_NUMBERED[nn]].Name);
                fprintf(fp, "Genotype %s Total   ", header);
                if (ou_percent) fprintf(fp, "Total Freq");
            }
            fprintf(fp, "\n");
            for (ll=0; ll<num_genos; ll++){
                if (empty_rows==1 || allele_dist[kk].GENOTOTAL(nn, ll)>0){
                    if (tab_text)
                        fprintf(fp, "%2d/%2d\t", genotypes[ll][0], genotypes[ll][1]);
                    else
                        fprintf(fp, "%2d/%2d    ", genotypes[ll][0], genotypes[ll][1]);
                    for (mm=0;mm<lclasses;mm++){
                        pos=3*mm;
                        if (tab_text) {
                            if (header_p[pos]){
                                fprintf(fp, "%5d\t", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 0));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 0));
                            }
                            if (header_p[pos+1]){
                                fprintf(fp, "%5d\t", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 1));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 1));
                            }
                            if (header_p[pos+2]){
                                fprintf(fp, "%5d\t", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 2));
                                if (ou_percent)
                                    fprintf(fp, "%5.4f\t", allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 2));
                            }
                        } /* tab_text */
                        else {
                            if (header_p[pos]){

                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f",
                                            allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 0),
                                            allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 0));
                                else
                                    fprintf(fp, "     %5d     ", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 0));

                            }
                            if (header_p[pos+1]){
                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f ",
                                            allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 1),
                                            allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 1));
                                else
                                    fprintf(fp, "     %5d    ", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 1));
                            }
                            if (header_p[pos+2]){
                                if (ou_percent)
                                    fprintf(fp, "%5d   %5.4f ",
                                            allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 2),
                                            allele_dist[kk].GSTATUSPERCENT(nn,ll, mm, 2));
                                else
                                    fprintf(fp, "     %5d     ", allele_dist[kk].GLIABLESTATUS(nn,ll, mm, 2));

                            }
                        } /* not tab_text */
                    }
                    if (tab_text) {
                        fprintf(fp, "%5d\t", allele_dist[kk].GENOTOTAL(nn, ll));
                        if (ou_percent)
                            fprintf(fp, "%5.4f",
                                    safe_divide(allele_dist[kk].GENOTOTAL(nn,ll),
                                                allele_dist[kk].GTOTAL(nn)));
                    } else {
                        if (ou_percent)
                            fprintf(fp, "%5d     %5.4f", allele_dist[kk].GENOTOTAL(nn, ll),
                                    safe_divide(allele_dist[kk].GENOTOTAL(nn,ll),
                                                allele_dist[kk].GTOTAL(nn)));
                        else
                            fprintf(fp, "%5d", allele_dist[kk].GENOTOTAL(nn, ll));
                    }
                    fprintf(fp, "\n");
                } /* non-zero row, or output all rows */
            } /* each genotype */
            /* write the column totals */
            if (tab_text)	fprintf(fp, "TOTAL\t");
            else 	fprintf(fp, "TOTAL    ");
            for (mm=0; mm<lclasses; mm++){
                ll=3*mm;
                if (tab_text) {
                    if (header_p[ll]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERGTOTAL(nn, mm, 0));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                    if (header_p[ll+1]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERGTOTAL(nn,mm,1));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                    if (header_p[ll+2]) {
                        fprintf(fp, "%5d\t", allele_dist[kk].MARKERGTOTAL(nn,mm,2));
                        if (ou_percent) fprintf(fp, "\t");
                    }
                } else {
                    if (header_p[ll]) {
                        if (ou_percent)
                            fprintf(fp, "%5d         ", allele_dist[kk].MARKERGTOTAL(nn, mm, 0));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERGTOTAL(nn,mm,0));

                    }
                    if (header_p[ll+1]) {
                        if (ou_percent)
                            fprintf(fp, "%5d          ", allele_dist[kk].MARKERGTOTAL(nn,mm,1));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERGTOTAL(nn,mm,1));
                    }
                    if (header_p[ll+2]) {
                        if (ou_percent)
                            fprintf(fp, "%5d          ", allele_dist[kk].MARKERGTOTAL(nn,mm,2));
                        else
                            fprintf(fp, "     %5d     ", allele_dist[kk].MARKERGTOTAL(nn,mm,2));
                    }
                }
            }
            fprintf(fp, "%5d\n\n\n", allele_dist[kk].GTOTAL(nn));
            free(genotypes);

            if (tab_text) {
                fclose(fp);
                free(file_ext);
                shorten_path(ou_file_name, shortname);
                sprintf(err_msg, "      Genotype count summary file:  %s",
                        shortname);
                mssgf(err_msg);
                free(ou_file_name);
                file_ext = NULL; ou_file_name = NULL;
            }
            if (header != NULL) {
                free(header);
            }
        } /* each marker */
        free(header_p);
    } /* each trait */

    if (!tab_text)  fclose(fp);
}
