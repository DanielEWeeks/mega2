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
#include "genetic_utils_ext.h"
#include "liability_dist_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc my_malloc
      genetic_utils_ext.h:  create_genos genoindx safe_divide
     liability_dist_ext.h:  get_liablfl_name get_output_options write_liable_dist
            linkage_ext.h:  get_loci_on_chromosome get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  create_mssg shorten_path
              utils_ext.h:  draw_line log_line summary_time_stamp
*/




static liable_allele_dist *alloc_allele_dist(int num_affected, int *loci_affected,
                                             int num_numbered, int *loci_numbered,
                                             linkage_locus_top *Llocus_top)

{
    liable_allele_dist *LiableAlleleDist;
    int locus_num, num_classes, num_alleles, num_genos;
    int mm, nn, jj, kk;
    int aff_locus_num;
    int **genos, **matrix;

    LiableAlleleDist = CALLOC((size_t) num_affected, liable_allele_dist);

    /* fill in the affected loci */
    for (mm=0; mm <num_affected; mm++){
        aff_locus_num=loci_affected[mm];
        LiableAlleleDist[mm].locus_AFFECTION = aff_locus_num;
        num_classes=Llocus_top->Locus[aff_locus_num].Data.Affection.ClassCnt;
        LiableAlleleDist[mm].num_liabls=num_classes;
        LiableAlleleDist[mm].MarkerLiability= MALLOC_STR(num_numbered, marker_liabl_type);
        LiableAlleleDist[mm].loci_NUMBERED = MALLOC_STR(num_numbered, int);
        for (nn=0; nn<num_numbered; nn++){
            locus_num=loci_numbered[nn];
            LiableAlleleDist[mm].loci_NUMBERED[nn]=locus_num;
            LiableAlleleDist[mm].MarkerLiability[nn].MarkerLocus=locus_num;
            LiableAlleleDist[mm].MarkerLiability[nn].num_alleles=
                Llocus_top->Locus[locus_num].AlleleCnt;
            num_alleles=Llocus_top->Locus[locus_num].AlleleCnt;
            /* allocate matrix for allele columns sums */
            LiableAlleleDist[mm].MarkerLiability[nn].MarkerTotal=CALLOC((size_t) num_classes, int *);
            for(jj=0; jj<num_classes; jj++) {
                LiableAlleleDist[mm].MarkerLiability[nn].MarkerTotal[jj]=CALLOC((size_t) 3, int);
                LiableAlleleDist[mm].MarkerLiability[nn].MarkerTotal[jj][0] =
                    LiableAlleleDist[mm].MarkerLiability[nn].MarkerTotal[jj][1] =
                    LiableAlleleDist[mm].MarkerLiability[nn].MarkerTotal[jj][2] = 0;
            }
/*       /\* allocate and initialize allele-related fields *\/ */
            LiableAlleleDist[mm].MarkerLiability[nn].AlleleTotal=
                MALLOC_STR(num_alleles, int);
            LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability =
                MALLOC_STR(num_alleles, allele_liabl_type);
            /* allocate and initialize genotype-related fields */
            /* create the genotypes--allele_pairs */
            num_genos=NUMGENOS(num_alleles); genos=create_genos(num_alleles);
            LiableAlleleDist[mm].MarkerLiability[nn].num_genos=num_genos;

            /* allocate matrix for genotype column sums */
            matrix= MALLOC_STR(num_classes, int *);
            for(jj=0; jj<num_classes; jj++) {
                matrix[jj] = MALLOC_STR(3, int);
                matrix[jj][0]=0; matrix[jj][1]=0; matrix[jj][2]=0;
            }
            LiableAlleleDist[mm].MarkerLiability[nn].MarkerGTotal=matrix;
            LiableAlleleDist[mm].MarkerLiability[nn].GenoTotal=
                MALLOC_STR(num_genos, int);
            LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability =
                MALLOC_STR(num_genos, geno_liabl_type);
            /* allocate each allele_record structure */
            for (kk=0; kk<num_alleles; kk++){
                LiableAlleleDist[mm].MarkerLiability[nn].AlleleTotal[kk]=0;
                LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Allele=kk+1;
                LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Liabls=
                    MALLOC_STR(num_classes, liabl_type);
                for (jj=0; jj<num_classes; jj++){
                    LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Liabls[jj].Class=
                        jj+1;
                    LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Liabls[jj].Status[0]=
                        0;
                    LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Liabls[jj].Status[1]=
                        0;
                    LiableAlleleDist[mm].MarkerLiability[nn].AlleleLiability[kk].Liabls[jj].Status[2]=
                        0;
                }
            }
            /* allocate each genotype_record structure */
            for (kk=0; kk<num_genos; kk++){
                LiableAlleleDist[mm].MarkerLiability[nn].GenoTotal[kk]=0;
                LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Allele1=genos[kk][0];
                LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Allele2=genos[kk][1];
                LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Liabls=
                    MALLOC_STR(num_classes, liabl_type);
                for (jj=0; jj<num_classes; jj++){
                    LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Liabls[jj].Class=
                        jj+1;
                    LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Liabls[jj].Status[0]=
                        0;
                    LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Liabls[jj].Status[1]=
                        0;
                    LiableAlleleDist[mm].MarkerLiability[nn].GenoLiability[kk].Liabls[jj].Status[2]=
                        0;
                }
            }
        } /* num_numbered */
    } /* num_affected */
    return LiableAlleleDist;
}

static void free_allele_dist(liable_allele_dist *allele_dist, int num_affection, int num_NUMBERED)

{
    int ii, jj, kk;
    int num_ALLELES;

    for (ii=0; ii<num_affection; ii++) {
        for (jj=0; jj<num_NUMBERED; jj++){
            num_ALLELES=allele_dist[ii].MarkerLiability[jj].num_alleles;
            for (kk=0; kk< num_ALLELES; kk++)
                free(allele_dist[ii].MarkerLiability[jj].AlleleLiability[kk].Liabls);
            free(allele_dist[ii].MarkerLiability[jj].AlleleLiability);
            free(allele_dist[ii].MarkerLiability[jj].GenoLiability);
            free(allele_dist[ii].MarkerLiability[jj].AlleleTotal);
            free(allele_dist[ii].MarkerLiability[jj].GenoTotal);
            free(allele_dist[ii].MarkerLiability[jj].MarkerTotal);
            free(allele_dist[ii].MarkerLiability[jj].MarkerGTotal);
        }
        free(allele_dist[ii].MarkerLiability);
        free(allele_dist[ii].loci_NUMBERED);
    }
    free(allele_dist);
}



static liable_allele_dist *create_allele_dist(linkage_ped_top *LPedTreeTop,
                                              int num_NUMBERED, int *loci_NUMBERED,
                                              int num_AFFECTION, int *loci_AFFECTION)

{
    liable_allele_dist *allele_dist;

    int ii, jj, kk, nn, mm, ll, num_alleles, num_genos;
    int affected_status, affected_class, allele1, allele2;
    int num_classes, t_locus, m_locus;
    linkage_locus_top *ltop;
    linkage_ped_tree ped_tree;


    ltop=LPedTreeTop->LocusTop;
    /* create the allele_distribution data_structure:
       a list of length= number of AFFECTED loci
       each element of type liable_allele_dist */
    allele_dist=
        alloc_allele_dist(num_AFFECTION, loci_AFFECTION, num_NUMBERED, loci_NUMBERED, ltop);
    for (kk=0; kk<num_AFFECTION; kk++){
        /* for each trait locus of affection type */
        t_locus=loci_AFFECTION[kk];
        num_classes=ltop->Locus[t_locus].Data.Affection.ClassCnt;
        for (mm=0; mm<num_NUMBERED; mm++){
            /* for each marker */
            m_locus=loci_NUMBERED[mm];
            allele_dist[kk].MARKERLIABLE(mm).Total=0;
            allele_dist[kk].MARKERLIABLE(mm).GTotal=0;
            num_alleles=ltop->Locus[m_locus].AlleleCnt;
            num_genos=NUMGENOS(num_alleles);
            for(nn=0; nn<num_classes; nn++){
                /* for each lability class */
                for (ii=0; ii<LPedTreeTop->PedCnt; ii++){
                    ped_tree=LPedTreeTop->Ped[ii];
                    for (jj=0; jj<ped_tree.EntryCnt; jj++){
                        /* for each entry in each pedigree */
                        affected_status=ped_tree.Entry[jj].ESTATUS(t_locus);
                        affected_class=ped_tree.Entry[jj].ECLASS(t_locus);
                        if (affected_class==(nn+1)) {
                            /* current_class is nn */
                            /* assign to the allele dist */
                            allele1= ped_tree.Entry[jj].EALLELE1(m_locus);
                            if (allele1>0){
                                ll = allele1 - 1;
                                allele_dist[kk].LIABLESTATUS(mm,ll,nn,affected_status) += 1;
                                allele_dist[kk].MARKERTOTAL(mm, nn, affected_status) += 1;
                                allele_dist[kk].ALLELETOTAL(mm, ll) += 1;
                            }
                            allele2=ped_tree.Entry[jj].EALLELE2(m_locus);
                            if (allele2>0) {
                                ll = allele2 - 1;
                                allele_dist[kk].LIABLESTATUS(mm,ll,nn,affected_status) += 1;
                                allele_dist[kk].MARKERTOTAL(mm, nn, affected_status) += 1;
                                allele_dist[kk].ALLELETOTAL(mm, ll) += 1;
                            }

                            /* assign to the genotype dist */
                            if (allele1>0 && allele2 >0) {
                                ll=genoindx(allele1, allele2);
                                allele_dist[kk].GLIABLESTATUS(mm,ll,nn,affected_status) += 1;
                                allele_dist[kk].MARKERGTOTAL(mm, nn, affected_status) += 1;
                                allele_dist[kk].GENOTOTAL(mm, ll) += 1;
                            }
                        } /* if affected_class */
                    }	/* each entry */
                }     /* each pedigree */
                /* loop for each allele to calculate allele totals*/
                for (ll=0; ll<num_alleles; ll++){
                    allele_dist[kk].STATUSPERCENT(mm,ll,nn,0)=
                        safe_divide(allele_dist[kk].LIABLESTATUS(mm,ll,nn,0),
                                    allele_dist[kk].MARKERTOTAL(mm, nn, 0));
                    allele_dist[kk].STATUSPERCENT(mm,ll,nn,1)=
                        safe_divide(allele_dist[kk].LIABLESTATUS(mm,ll,nn,1),
                                    allele_dist[kk].MARKERTOTAL(mm, nn, 1));
                    allele_dist[kk].STATUSPERCENT(mm,ll,nn,2)=
                        safe_divide(allele_dist[kk].LIABLESTATUS(mm,ll,nn,2),
                                    allele_dist[kk].MARKERTOTAL(mm, nn, 2));
#ifdef DEBUG
                    /*  printf("percent0, %3.2f, percent1 %3.2f, percent2, %3.2f\n",
                        allele_dist[kk].STATUSPERCENT(mm,ll,nn,0),
                        allele_dist[kk].STATUSPERCENT(mm,ll,nn,1),
                        allele_dist[kk].STATUSPERCENT(mm,ll,nn,2)); */
#endif
                } /* for each allele */
                /* Loop for each genotype */
                for (ll=0; ll<num_genos; ll++){
                    allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,0)=
                        safe_divide(allele_dist[kk].GLIABLESTATUS(mm,ll,nn,0),
                                    allele_dist[kk].MARKERGTOTAL(mm, nn, 0));
                    allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,1)=
                        safe_divide(allele_dist[kk].GLIABLESTATUS(mm,ll,nn,1),
                                    allele_dist[kk].MARKERGTOTAL(mm, nn, 1));
                    allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,2)=
                        safe_divide(allele_dist[kk].GLIABLESTATUS(mm,ll,nn,2),
                                    allele_dist[kk].MARKERGTOTAL(mm, nn, 2));
#ifdef DEBUG
                    /*	  printf("gpercent0, %3.2f, gpercent1 %3.2f, gpercent2, %3.2f\n",
                          allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,0),
                          allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,1),
                          allele_dist[kk].GSTATUSPERCENT(mm,ll,nn,2));
                    */
#endif
                }       /* for each genotype */
            }         /* each class */
            for (ll=0; ll<num_alleles; ll++)
                allele_dist[kk].TOTAL(mm) += allele_dist[kk].ALLELETOTAL(mm, ll);
            for (ll=0; ll<num_genos; ll++)
                allele_dist[kk].GTOTAL(mm) += allele_dist[kk].GENOTOTAL(mm, ll);
        }         /* each marker */
    }             /* each ffection locus */
    return allele_dist;
}
/*
  NM - liability summary does not allow multiple traits to be combined, so don't have
  to check for marker item in global_trait_entries.

*/
void liability_summary(char **input_files, linkage_ped_top *LPedTreeTop,
		       int *numchr, int untyped_ped_opt)

{
    int *loci_NUMBERED, *loci_AFFECTION;
    int total_NUMBERED=0, num_NUMBERED=0, num_AFFECTION=0;
    int i_NUM=0, i_AFF=0, kk;
    int ii;
    linkage_locus_top *ltop;
    liable_allele_dist *allele_liabl_dist;
    int ou_percent, empty_rows, empty_cols, tab_text;
    char **liablefl_name, flname[FILENAME_LENGTH];
    FILE *fp;
    /* printf("did nothing\n"); */

    /* get the file name */
    ou_percent = 1;
    empty_cols = 0;   empty_rows=0;
    tab_text = 0;
    printf("\nOption: Count alleles and genotypes within groups: \n");
    draw_line();
    get_output_options(&ou_percent, &empty_cols, &empty_rows, &tab_text);

    printf("\nOption: Count alleles and genotypes within groups: \n");
    draw_line();
    liablefl_name = get_liablfl_name(numchr, tab_text);

    fclose(fopen(liablefl_name[0], "w"));
    ltop=LPedTreeTop->LocusTop;

    num_AFFECTION=0;
    i_AFF=0;
    for (kk=0; kk<num_traits; kk++) {
        if (global_trait_entries[kk] < 0) continue;
        if (ltop->Locus[global_trait_entries[kk]].Type == AFFECTION)
            num_AFFECTION += 1;
    }


    log_line(mssgf);
    if (num_AFFECTION > 1) {
        sprintf(err_msg, "Found %d affection traits:", num_AFFECTION);
        mssgf(err_msg);
    } else {
        mssgf("Found 1 affection trait:");
    }

    loci_AFFECTION = CALLOC((size_t) num_AFFECTION, int);

    for (kk=0; kk < num_traits; kk++) {
        if (global_trait_entries[kk] < 0) continue;
        if (ltop->Locus[global_trait_entries[kk]].Type == AFFECTION) {
            loci_AFFECTION[i_AFF++]=global_trait_entries[kk];
            sprintf(err_msg, "  %s", ltop->Locus[global_trait_entries[kk]].Name);
            mssgf(err_msg);
        }
    }

    if (main_chromocnt > 1) {
        get_loci_on_chromosome(0);
        for (ii=0; ii < main_chromocnt; ii++) {
            if (global_chromo_entries[ii] == UNKNOWN_CHROMO) {
                get_unmapped_loci(1);
                break;
            }
        }
    } else {
        if (*numchr == UNKNOWN_CHROMO) {
            get_unmapped_loci(0);
        } else {
            get_loci_on_chromosome(*numchr);
        }
    }

    num_NUMBERED=0;
    i_NUM=0;
    omit_peds(untyped_ped_opt, LPedTreeTop);
    for (kk=0; kk<NumChrLoci; kk++) {
        if (ltop->Locus[ChrLoci[kk]].Type == NUMBERED ||
            ltop->Locus[ChrLoci[kk]].Type == BINARY) {
            num_NUMBERED += 1;
        }
    }
    loci_NUMBERED = CALLOC((size_t) num_NUMBERED, int);

    for (kk=0; kk<NumChrLoci; kk++) {
        if (ltop->Locus[ChrLoci[kk]].Type == NUMBERED ||
            ltop->Locus[ChrLoci[kk]].Type == BINARY) {
            loci_NUMBERED[i_NUM++]= ChrLoci[kk];
        }
    }


    allele_liabl_dist= create_allele_dist(LPedTreeTop, num_NUMBERED,
                                          loci_NUMBERED,
                                          num_AFFECTION, loci_AFFECTION);

    fp = fopen(liablefl_name[0], "w") ;
    summary_time_stamp(input_files, fp, "");
    fclose(fp);
    create_mssg(TO_LIABLE_FREQ);
    write_liable_dist(allele_liabl_dist, liablefl_name,
                      LPedTreeTop, num_AFFECTION, num_NUMBERED, ou_percent,
                      empty_rows, empty_cols, tab_text);
    /* create the text summary in addition to the tab text file */
    if (tab_text == 1)
        write_liable_dist(allele_liabl_dist, liablefl_name,
                          LPedTreeTop, num_AFFECTION, num_NUMBERED, ou_percent,
                          empty_rows, empty_cols, 0);

    free(loci_NUMBERED);
    loci_NUMBERED= NULL;
    free_allele_dist(allele_liabl_dist, num_AFFECTION, num_NUMBERED);
    total_NUMBERED += num_NUMBERED;
    /*    if (tab_text == 1) {log_line(mssgf); } */


    shorten_path(liablefl_name[0], flname);
    sprintf(err_msg, "      Overall summary file:         %s", flname);
    mssgf(err_msg);

    free(loci_AFFECTION);
    /*  free_all_including_lpedtop(&LPedTreeTop); */
    return;
}
