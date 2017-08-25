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

#include "common.h"
#include "typedefs.h"

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "read_files_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  check_map_positions haldane_theta kosambi_theta
            linkage_ext.h:  connect_loops copy_linkage_ped_top free_all_including_lpedtop get_loci_on_chromosome new_lpedtop
           makenucs_ext.h:  makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         read_files_ext.h:  free_lpedtop_pedinfo
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
*/



/*------------------ static functions ----------------------*/

static void write_SAGE_cntcsh(char *cntcshfl_name, char *cntfl_name,
			      char *pedfl_name,
			      char *cntpar_name, char *parfl_name,
			      linkage_ped_top *Top,
			      int num_fams, int orig_loopover_opt);
static void write_SAGE_csh(char *cshfl_name, char *outfl_name, char *batfl_name, char *sibfl_name,
                           char* loutfl_name,linkage_ped_top *Top,
                           int num_fams);
static void write_SAGE_sibpal(char *poutfl_name, int numchr, linkage_ped_top *LLTop,
                              analysis_type analysis, int pwid, int fwid);
static void write_SAGE_cntpar(char *batfl_name, int numchr, linkage_ped_top *LLTop,
			      analysis_type analysis,
			      int pwid, int fwid);
static void write_SAGE_par(char *batfl_name, int numchr,
			   linkage_ped_top *LLTop, analysis_type analysis,
			   int pwid, int fwid);
static void sagewrite_quantitative_data(FILE *filep, int locusnm,
					linkage_locus_rec *locus,
					linkage_ped_rec *entry);
static void sage4write_quantitative_data(FILE *filep, int locusnm,
					 linkage_locus_rec *locus,
                                         linkage_ped_rec *entry);
static void sagewrite_affection_data(FILE *filep, int locusnm,
                                     linkage_locus_rec *locus,
                                     linkage_ped_rec *entry);
static void sagewrite_binary_data(FILE *filep, int locusnm,
				  linkage_locus_rec  *locus,
				  linkage_ped_rec *entry);
static int save_SAGE_peds(char *outfl_name, linkage_ped_top *Top,
			  analysis_type analysis, int pwid, int fwid);
static void  write_SAGE_locus_file(char *loutfl_name,
				   linkage_locus_top * LTop);
static void sage_file_names(char *file_names[], linkage_ped_top *Top);

static void write_SAGE_cnt_file(char *poutfl_name,
				linkage_ped_top *Top,
				int pwid, int fwid);

static void write_SAGE4_par(char *poutfl_name, linkage_ped_top *Top);
static void write_SAGE4_map(char *mapfl_name, int numchr, linkage_locus_top *Top);


/*------------------exported function(s) -----------------------*/
void create_SAGE_file(linkage_ped_top **LPedTreeTop, char *mapfl_name,
                      int *trait, ped_top *PedTreeTop,
                      file_format *infl_type, file_format *outfl_type,
                      int *numchr, analysis_type *analysis,
                      char *file_names[], int untyped_ped_opt);

void create_SAGE4_file(linkage_ped_top **LPedTreeTop,
                       ped_top *PedTreeTop, file_format *infl_type,
                       int *numchr, analysis_type *analysis,
		       char *file_names[], int untyped_ped_opt);

/*--------------------------------------------------------------*/

//
// The current version of SAGE supports the 'SIBPAL' program, however
// the arguments and data structure of the input files seems to be
// completely different from this. See...
// S.A.G.E.v6.3 UserDoc Chapter 18.3 [SIBPAL] Program Input
static void   write_SAGE_sibpal(char *sibfl_name, int numchr,
				linkage_ped_top *LLTop,
				analysis_type analysis,
				int pwid, int fwid)
{
    FILE           *filep;
    char           trait_name[11], sfl[2*FILENAME_LENGTH];
    int            i, *trp, tr, nloop, num_affec=num_traits;
    int            ntr, nmrk;
/*  char marker_name[MAXLOCI][FILENAME_LENGTH]; */

    nmrk = ((LoopOverTrait == 0)? NumChrLoci - num_affec + 1:
            NumChrLoci - num_affec);
    ntr =  ((LoopOverTrait == 0)? num_affec - 1: num_affec);
    /*  printf("%s", sibfl_name); */
    NLOOP;

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop> 1 && tr==0) continue;
        SKIP_TRI(tr-1);
        strcpy(trait_name, "");
        sprintf(sfl, "%s/%s", output_paths[tr], sibfl_name);
        if ((filep = fopen(sfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", sfl);
            EXIT(FILE_WRITE_ERROR);
        }
        if (LoopOverTrait == 1 || num_traits == 2) {
            /* single trait */
            if (strlen(LLTop->LocusTop->Pheno[*trp].TraitName) > 10) {
                strncpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName, (size_t) 10);
                trait_name[10]='\0';
            } else {
                strcpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName);
            }

            if (LLTop->LocusTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "SIBPAL: %s\n", trait_name);
                fprintf(filep, "    5    5    5    0\n");
                fprintf(filep, "    2%5d    1    1    0    1    0\n", nmrk);

                fprintf(filep, "%10s    0   0     2    1\n", trait_name);
                fprintf(filep, "(A3,1X,I%d,1X,A%d,1X,%dx,F2.0,T%d,%d(A6))\n",
                        fwid, pwid,(2*pwid+4),
                        (4+fwid+1+3*(pwid+1)+2+2+1),
                        nmrk);
            } else  if (LLTop->LocusTop->Locus[*trp].Type == QUANT) {
                fprintf(filep, "SIBPAL: %s\n", trait_name);
                fprintf(filep, "    5    5    5    0\n");
                fprintf(filep, "    1%5d    1    1    0    1    0\n",
                        nmrk);
                fprintf(filep, "%10s    0%4d.    2    1\n", trait_name, (int)MissingQuant);
                fprintf(filep, "(A3,1X,I%d,1X,A%d,1X,%dx,F7.3,1X,T%d,%d(A6))\n",
                        fwid, pwid, (2*pwid+4),
                        (4+fwid+1+3*(pwid+1)+2+8+1),
                        nmrk);
            }
            trp++;
        } else {
            /* multiple traits, assume that reordering made sure
               that all traits are of the same type and at the
               beginning of the order */
            if (num_affec == 0) {
                fprintf(filep, "SIBPAL: MARKERS \n");
                fprintf(filep, "    5    5    5    0\n");
                fprintf(filep, "    4%5d    0    0    1    0\n", nmrk);
                /*	fprintf(filep, "no-trait      0   0.    2    1\n"); */
                fprintf(filep, "(A3,1X,I%d,1X,A%d,1X,%dx,T%d,%d(A6))\n",
                        fwid, pwid, (2*pwid+4),
                        (4+fwid+1+3*(pwid+1)+2+1),
                        nmrk);
            } else {
                fprintf(filep, "SIBPAL: MULTI-TRAIT\n");
                fprintf(filep, "    5    5    5    0\n");
                if (HasAff) {
                    fprintf(filep, "    2%5d%5d    0    1    0\n",
                            nmrk, ntr);
                    if (strlen(LLTop->LocusTop->Pheno[*trp].TraitName) > 10) {
                        strncpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName,
                                (size_t) 10);
                        trait_name[10]='\0';
                    } else {
                        strcpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName);
                    }

                    for (i=0; i < num_affec; i++) {
                        if (global_trait_entries[i] < 0) continue;
                        if (strlen(LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName) > 10) {
                            strncpy(trait_name,
                                    LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName, (size_t) 10);
                            trait_name[10] = '\0';
                        } else {
                            strcpy(trait_name,
                                   LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName);
                        }
                        fprintf(filep, "%10s   0    0     2    1\n", trait_name);
                    }
                    fprintf(filep, "(A3,1X,I%d,1X,A%d,1X,%dx,%d(F2.0),T%d,%d(A6))\n",
                            fwid, pwid, (2*pwid+2+2),
                            ntr,
                            (4+fwid+1+3*(pwid+1)+2+2*ntr+1),
                            nmrk);
                } else {
                    fprintf(filep, "    1%5d%5d    0    1    0\n",
                            nmrk, ntr);
                    for (i=0; i < num_affec; i++) {
                        if (global_trait_entries[i] < 0) continue;
                        if (strlen(LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName) > 10) {
                            strncpy(trait_name, LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName, (size_t) 10);
                            trait_name[10] = '\0';
                        } else {
                            strcpy(trait_name,
                                   LLTop->LocusTop->Pheno[global_trait_entries[i]].TraitName);
                        }
                        fprintf(filep, "%10s   0 %4d.    2    1\n", trait_name, (int)MissingQuant);
                    }
                    fprintf(filep,
                            "(A3,1X,I%d,1X,A%d,1X,%dx,%d(F7.3,1X),T%d,%d(A6))\n",
                            fwid, pwid, (4+2*pwid), num_affec,
                            (4+fwid+1+3*(pwid+1)+2+8*ntr+1), nmrk);
                }
            }
        }
        fclose(filep);
        if (nloop==1) break;

    }
} /*  write_SAGE_sibpal */

static void write_SAGE_cntcsh(char *cshfl_name, char *cntfl_name,
			      char *pedfl_name,
			      char *cntpar_name, char *parfl_name,
			      linkage_ped_top *Top, int num_fams,
			      int multiple_affs)
{
    int *trp, tr, nloop, num_affec=num_traits;
    FILE *csh_unit;
    char     syscmd[2*FILENAME_LENGTH], cfl[2*FILENAME_LENGTH];
    int i, num_inds=0;

    for (i=0; i<Top->PedCnt; i++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[i]) continue;
        }
        num_inds += Top->Ped[i].EntryCnt;
    }

    NLOOP;

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        if (nloop == 1 && LoopOverTrait == 0) {
            sprintf(cfl, "%s/%s", output_paths[0], cshfl_name);
        } else {
            SKIP_TRI(tr-1)
                sprintf(cfl, "%s/%s", output_paths[tr], cshfl_name);
        }
        if ((csh_unit=fopen(cfl,"w")) == NULL)   {
            draw_line();
            errorvf("Unable to open file '%s'\n", cfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(csh_unit,"#!/bin/csh -f\n");
        script_time_stamp(csh_unit);
        fprintf(csh_unit,"unalias rm\n");
        fprintf(csh_unit,"unalias cp\n");
        fprintf(csh_unit,
                "foreach f (\"fort.1\" \"fort.2\" \"fort.11\" \"fort.12\" \"fort.21\" \"fort.22\" \"fort.13\" \"fort.23\")\n");
        fprintf(csh_unit, "    if (-e $f)  then \n");
        fprintf(csh_unit, "       rm $f\n");
        fprintf(csh_unit, "    endif\n");
        fprintf(csh_unit, "end\n");

        /* to run FSP: */
        /* 1. copy  sage_par.## to fort.1 */
        if (multiple_affs == 1) {
            /* copy the parameter file from the top level directory */
            fprintf(csh_unit,"cp ../%s fort.1\n", parfl_name);
        } else {
            /* Copy parameter file  from current directory */
            fprintf(csh_unit,"cp %s fort.1\n", parfl_name);
        }

        /* 2.
           for affs : copy sage_cnt.## fort.11
           for qtls : copy sage_ped.## fort.11
           cntfl_name will be set appropriately outside this routine
        */
        if (Top->LocusTop->Locus[*trp].Type == AFFECTION) {
            fprintf(csh_unit,"cp %s fort.11\n", cntfl_name);
        } else {
            fprintf(csh_unit,"cp %s fort.11\n", pedfl_name);
        }

        /* 3. Run  FSP */
        fprintf(csh_unit,"echo Running fsp on %s and %s\n",
                parfl_name, cntfl_name);
        fprintf(csh_unit,"fsp -f %d -i %d\n", num_fams, num_inds);

        /* To run FCOR */
        /* 1. Copy sage_cntpar.## to fort.1 */
        fprintf(csh_unit,"cp %s fort.1\n", cntpar_name);
        /* 2. Create the fort.11 file:
           if aff then copy sage_cnt.## fort.11
           otherwise copy sage_ped.## into fort.11 */
        if (Top->LocusTop->Locus[*trp].Type == AFFECTION) {
            fprintf(csh_unit,"cp %s fort.11\n", cntfl_name);
        } else {
            fprintf(csh_unit,"cp %s fort.11\n", pedfl_name);
        }
        /* 3. Create the fort.12 file: */
        fprintf(csh_unit,"cp fort.22 fort.12\n");
        /* 4. delete all fort.2* files */
        fprintf(csh_unit,"rm fort.2*\n");
        /* 5. Run fcor */
        fprintf(csh_unit,"echo Running fcor using %s and fsp-generated files\n",
                cntfl_name);

        fprintf(csh_unit,"fcor -i %d\n", num_inds);
        fprintf(csh_unit,
                "echo The file fort.21 contains the detailed output from FCOR\n");
        fprintf(csh_unit,
                "echo The file fort.22 contains the summary counts from FCOR\n");
        fprintf(csh_unit,
                "echo WARNING  Missing count value assumed to be %d\n",
                (int)MissingQuant);
        fclose(csh_unit);

        sprintf(syscmd, "chmod +x %s", cfl);
        System(syscmd);
        if (nloop==1) break;
        trp++;
    }
}

static void     write_SAGE_csh(char *cshfl_name, char *outfl_name,
			       char *batfl_name, char *sibfl_name,
			       char *loutfl_name,
			       linkage_ped_top *Top,
			       int num_fams)

{
    FILE *csh_unit;
    char cfl[2*FILENAME_LENGTH], syscmd[2*FILENAME_LENGTH];
    int num_inds=0;
    int i, tr, nloop, num_affec=num_traits;


    /* first check the number of individuals */

    for (i=0; i<Top->PedCnt; i++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[i]) continue;
        }
        num_inds += Top->Ped[i].EntryCnt;
    }

    NLOOP;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        if (nloop == 1 && LoopOverTrait == 0) {
            sprintf(cfl, "%s/%s", output_paths[0], cshfl_name);
        } else {
            SKIP_TRI(tr);
            sprintf(cfl, "%s/%s", output_paths[tr], cshfl_name);
        }
        if ((csh_unit=fopen(cfl,"w")) == NULL)   {
            errorvf("Unable to open file '%s'\n", cfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(csh_unit,"#!/bin/csh -f\n");
        script_time_stamp(csh_unit);
        fprintf(csh_unit,"unalias rm\n");
        fprintf(csh_unit,"unalias cp\n");
        fprintf(csh_unit,
                "foreach f (\"fort.1\" \"fort.2\" \"fort.11\" \"fort.12\" \"fort.21\" \"fort.22\" \"fort.23\")\n");
        fprintf(csh_unit, "    if (-e $f)  then \n");
        fprintf(csh_unit, "       rm $f\n");
        fprintf(csh_unit, "    endif\n");
        fprintf(csh_unit, "end\n");
        fprintf(csh_unit,"cp %s fort.11\n",outfl_name);
        fprintf(csh_unit,"cp %s fort.1\n",batfl_name);
        fprintf(csh_unit,"echo Running fsp on %s and %s\n",batfl_name,outfl_name);
        fprintf(csh_unit,"fsp -f %d -i %d\n",
                num_fams, num_inds);
        if (num_traits > 1 || (LoopOverTrait == 1 && num_traits >= 1)) {
            fprintf(csh_unit,"mv -f fort.22 fort.13\n");
            fprintf(csh_unit,"mv -f fort.11 fort.12\n");
            fprintf(csh_unit,"if (-e fort.11) then\n");
            fprintf(csh_unit,"    rm fort.11\n");
            fprintf(csh_unit,"endif\n");
            fprintf(csh_unit,"cp %s fort.11\n",loutfl_name);
            fprintf(csh_unit,"if (-e fort.1) then\n");
            fprintf(csh_unit,"    rm fort.1\n");
            fprintf(csh_unit,"endif\n");
            fprintf(csh_unit,"cp %s fort.1\n",sibfl_name);
            fprintf(csh_unit,"rm fort.2*\n");
            fprintf(csh_unit,
                    "echo Running sibpal using %s and fsp-generated files\n",
                    sibfl_name);
            fprintf(csh_unit,
                    "sibpal -l %d -f %d -i %d\n",
                    ((LoopOverTrait == 0)? NumChrLoci :
                     NumChrLoci - num_affec + 1),  num_fams, num_inds);
            fprintf(csh_unit,
                    "echo The file fort.21 contains the detailed output from SIBPAL\n");
            fprintf(csh_unit,
                    "echo The file fort.22 contains the summary output from SIBPAL\n");
            fprintf(csh_unit,
                    "echo The file fort.23 contains the optional output from SIBPAL\n");
            fprintf(csh_unit,
                    "echo WARNING  Missing trait value assumed to be %d\n",
                    (int)MissingQuant);
        }
        fclose(csh_unit);
        sprintf(syscmd, "chmod +x %s", cfl);
        System(syscmd);
        if (nloop==1) break;
    }
}

static void  write_SAGE_cntpar(char  *cbatfl_name, int numchr,
			       linkage_ped_top *LLTop, analysis_type analysis,
			       int pwid, int fwid)
{
    int            i, *trp, tr, nloop, num_affec=num_traits;
    char           cfl[2*FILENAME_LENGTH], trait_name[11];
    FILE           *filep;

    NLOOP;

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop> 1 && tr==0) continue;
        SKIP_TRI(tr-1);
        sprintf(cfl, "%s/%s", output_paths[tr], cbatfl_name);

        if ((filep = fopen(cfl, "w")) == NULL)  {
            errorvf("Opening SAGE parameter file %s.\n", cfl);
            EXIT(FILE_WRITE_ERROR);
        }
        if (LoopOverTrait == 1 || num_traits == 1) {
            if (LLTop->LocusTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "Genotyped affected relative pairs for %s\n",
                        LLTop->LocusTop->Pheno[*trp].TraitName);
                fprintf(filep, "    1    1    0    0    0    0    0\n");
                fprintf(filep, "MF\n");
                fprintf(filep, "(4X,I%d,1X,A%d,1X,T%d,A1,1X,F8.0)\n",
                        fwid, pwid, (4+fwid+1+3*(pwid+1)+1));
                /* length("FSP ") + family_id + 3*per_id + 1 +tab to" */
            } else {
                fprintf(filep, "Correlations for QTL %s\n",
                        LLTop->LocusTop->Pheno[*trp].TraitName);
                fprintf(filep, "    1    1\n");
                fprintf(filep, "MF\n");
                fprintf(filep, "(4X,I%d,1X,A%d,T%d,A1,1X,F7.3)\n",
                        fwid, pwid, (4+fwid+1+3*pwid+3+1));
                /* length("FSP ") + family_id + 1 + 3*(per_id + 1)" */
            }

            /* only the current trait locus */
            if (strlen(LLTop->LocusTop->Pheno[*trp].TraitName) > 10) {
                strncpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName, (size_t) 10);
                trait_name[10] = '\0';
            } else {
                strcpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName);
            }
            fprintf(filep, "%s\n", trait_name);
            if (LLTop->LocusTop->Locus[*trp].Type == QUANT) {
	      // Apparently the older versions of SAGE only support integer values for the
	      // value of a missing quantitative trait.
                fprintf(filep, "%10d\n", (int)MissingQuant);
            } else {
                fprintf(filep, "       0.0\n");
            }
        }

        else if (LoopOverTrait == 0 && num_traits > 0) {
            /* Only possible if all are QTLs */
            fprintf(filep, "QTL analysis for all traits\n");
            fprintf(filep, "%5d    1\n", num_traits);
            fprintf(filep, "MF\n");
            fprintf(filep, "(4X,I%d,1X,A%d,T%d,A1,1X,%d(F7.3,1X))\n",
                    fwid, pwid, (4+fwid+1+3*(pwid+1)+1), num_traits);
            /* length("FSP ") + family_id + 1 + 3*(per_id + 1) */

            /* Now the names and default values */
            for (i=0; i<num_traits; i++) {
                SKIP_TRII(*trp);
                if (strlen(LLTop->LocusTop->Pheno[*trp].TraitName) > 10) {
                    strncpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName,
                            (size_t) 10);
                    trait_name[10] = '\0';
                } else {
                    strcpy(trait_name, LLTop->LocusTop->Pheno[*trp].TraitName);
                }
                fprintf(filep, "%-10s", trait_name);
                trp++;
            }
            fprintf(filep, "\n");
            /* For all trait loci :  missing values */
            trp = &(global_trait_entries[0]);
            for (i=0; i<num_traits; i++) {
                SKIP_TRII(*trp);
                if (LLTop->LocusTop->Locus[*trp].Type == QUANT) {
                    fprintf(filep, "%10d", (int)MissingQuant);
                } else {
                    fprintf(filep, "       0.0");
                }
                trp++;
            }
            fprintf(filep, "\n");
        }

        fclose(filep);
        if (nloop==1) break;
        trp++; /* advance only for multiple iterations of the loop */
    }
} /*  write_SAGE_cntpar */

static void  write_SAGE_par(char *batfl_name, int numchr,
			    linkage_ped_top *LLTop,
			    analysis_type analysis, int pwid, int fwid)
{
    int            tr, nloop, num_affec=num_traits;
    char           pfl[2*FILENAME_LENGTH];
    FILE           *filep;

    NLOOP;

    for(tr=0; tr <=nloop; tr++) {
        if (nloop> 1 && tr==0) continue;
        sprintf(pfl, "%s/%s", output_paths[tr], batfl_name);
        if ((filep=fopen(pfl,"w")) == NULL)   {
            errorvf("Unable to open file '%s'\n", pfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(filep, "FSP parameter file\n");
        fprintf(filep, "    1    0    1    0    1MF\n");
        fprintf(filep, "(A3,1X,I%d,1X,3(A%d,1X),A1,1X)\n",
                fwid, pwid);
        fclose(filep);
        if (nloop==1) break;
    }
    return;
} /*  write_SAGE_par */


static void   sagewrite_quantitative_data(FILE *filep, int locusnm,
					  linkage_locus_rec *locus,
					  linkage_ped_rec *entry)
{
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON) {
        // This is the old sage for which I could not find source or documentation...
        // Should this be truncated to an int?
        if (ITEM_READ(Value_Missing_Quant_On_Output)) {
            fprintf(filep, " %s ",
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        } else {
	  // This should have been caught before this point, but just in case...
	  errorvf("The batch file item 'Value_Missing_Quant_On_Output' must be specified because\n");
	  errorvf("one or more QTL has been found to be missing in the input data.\n");
	  EXIT(BATCH_FILE_ITEM_ERROR);
	}
    } else {
        fprintf(filep, "%7.3f ", entry->Pheno[locusnm].Quant);
    }
}

//
// In 'opensage-read-only/core_code/trunk/src/c++/include/rped/rped.ipp' the
// numeric missing code is set to 'NaN' (not a number), and the string missing
// code to the empty string (e.g., ""). In the old Mega2 code "NA" was used for
// the default value of a missing quant in sage4+ mode. All of the Sage documentation
// that I have been able to find states that a missing quantitative value
// must be a number or a space.
// Since it has been impossible for me to find old SAGE documentation and code, only
// the more recent versions of SAGE have been adapted to use 'Value_Missing_Quant_On_Output'.
static void   sage4write_quantitative_data(FILE *filep, int locusnm,
                                           linkage_locus_rec *locus,
                                           linkage_ped_rec *entry)
{
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON) {
        //fprintf(filep, "   NA   ");
        if (ITEM_READ(Value_Missing_Quant_On_Output)) {
            fprintf(filep, " %s ",
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        } else {
	  // This should have been caught before this point, but just in case...
	  errorvf("The batch file item 'Value_Missing_Quant_On_Output' must be specified because\n");
	  errorvf("one or more QTL has been found to be missing in the input data.\n");
	  EXIT(BATCH_FILE_ITEM_ERROR);
	}
    } else {
        fprintf(filep, "%7.3f ", entry->Pheno[locusnm].Quant);
    }
}

static void  sagewrite_affection_data(FILE *filep, int locusnm,
				      linkage_locus_rec *locus,
				      linkage_ped_rec *entry)
{
    int aff = aff_status_entry(entry->Pheno[locusnm].Affection.Status, entry->Pheno[locusnm].Affection.Class,
                               locus);

    if (aff > 0)
        fprintf(filep, "%1d ", aff);
    else
        fprintf(filep, "  ");
    return;

}


static void  sagewrite_binary_data(FILE *filep,  int locusnm,
				   linkage_locus_rec *locus,
				   linkage_ped_rec *entry)
{
    int a1, a2;
    int allele;
/* Need to check for completely untyped person */
    for (allele = 0; allele < locus->AlleleCnt; allele++)
        {
            get_2alleles(entry->Marker, locusnm, &a1, &a2);
            if ((allele == a1 - 1)
                || (allele == a2 - 1))
                fprintf(filep, "1");
            else
                fprintf(filep, "0");
            for (allele = locus->AlleleCnt; allele < 7; allele++)
                fprintf(filep, " ");
        }
}

static void sagewrite_numbered_data(FILE *filep, const int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);
    // For a single character allele, all of these format strings will use 6 columns...
    if (a1 == 0)
        fputs(" 0/ 0 ", filep);
    else
        fprintf(filep, "%2s/%2s ", format_allele(locus, a1), format_allele(locus, a2));
}


/*
 * Save all the pedigrees to a file in SAGE format with
 * the first 'trait' field containing, for only the affecteds,
 * an indicator variable which is "1" if they are typed at
 * any loci.
 *  Assumes:  only one affection status locus
 */
static void  write_SAGE_cnt_file(char *poutfl_name,
				 linkage_ped_top *Top,
				 int pwid, int fwid)
{
#define LTOP Top->LocusTop

    int tr, *trp, nloop, ped, entry, locus1, num_affec = num_traits;
    FILE     *filep;
    char     pfl[2*FILENAME_LENGTH], fformat[10], pformat[10];
    linkage_ped_rec *Entry;
    int b;

    create_formats(fwid, pwid, fformat, pformat);
    /*  sprintf(fformat, "FSP %%%dd", fwid);
        sprintf(pformat, "%%%dd", pwid);
    */
    NLOOP;

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr<=nloop; tr++) {
        if (nloop>1 && tr==0) continue;

        /* SKIP_TRI(tr-1); */
        sprintf(pfl, "%s/%s", output_paths[tr], poutfl_name);
        if ((filep = fopen(pfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pfl);
            EXIT(FILE_WRITE_ERROR);
        }
        for (ped = 0; ped < Top->PedCnt; ped++)    {
            /* find the greatest allele number for each locus */
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* Write a dummy study id and the pedigree id as a  numeric code */
                fprintf(filep, "FSP ");
                fprintf(filep, fformat, Top->Ped[ped].Num);
                /* write the entry number  */
                fprintf(filep, pformat, entry + 1);
                /* write the parents */
                if (Entry->Father != 0) {
                    fprintf(filep, pformat, Entry->Father);
                    fprintf(filep, pformat, Entry->Mother);
                } else {
                    for(b=0; b < 2*(pwid+1); b++) fprintf(filep, " ");
                }
                /* write the sex */
                if (Entry->Sex == MALE_ID)
                    fprintf(filep, "M ");
                else if (Entry->Sex == FEMALE_ID)
                    fprintf(filep, "F ");
                else
                    fprintf(filep, "? ");
                /* now write the number of genotypes for any affected person */
                if (trp != NULL) {
                    locus1 = *trp;
                    SKIP_TRII(*trp);
                    if (Top->LocusTop->Locus[locus1].Type == AFFECTION) {
                        if (Entry->Pheno[locus1].Affection.Status == 2) {
                            fprintf(filep,"%8d ",Entry->Ngeno);
                        }
                        else
                            fprintf(filep,"         ");
                        sagewrite_affection_data(filep, locus1,
                                                 &(Top->LocusTop->Locus[locus1]), Entry);
                        fprintf(filep, "\n");
                    } else {
                        sagewrite_quantitative_data(filep, locus1,
                                                    &(Top->LocusTop->Locus[locus1]),
                                                    Entry);

                        fprintf(filep, "\n");
                    }
                }
            }
        }
        fclose(filep);
        if (nloop==1) break;
        trp++;
    }
} /* End of write_SAGE_cnt_file */


/*
 * Save all the pedigrees to a file in SAGE format. Return the
 * number of members saved (total).
 */
static int   save_SAGE_peds(char *outfl_name, linkage_ped_top *Top,
			    analysis_type analysis, int pwid, int fwid)
{
#define LTOP Top->LocusTop
    int       *trp, tr, nloop, ped, entry, locus1, num_affec = num_traits;
    char      pedfl[2*FILENAME_LENGTH], fformat[6], pformat[6];
    FILE      *filep;
    linkage_ped_rec *Entry;
    int  b, l1;

    NLOOP;

    create_formats(fwid, pwid, fformat, pformat);

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(pedfl, "%s/%s", output_paths[tr], outfl_name);

        if ((filep = fopen(pedfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }

        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            /* find the greatest allele number for each locus */
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);
                if (analysis == TO_SAGE) {
                    /* Write a dummy study id and the pedigree id as a 5 numeric code */
                    fprintf(filep, "FSP ");
                    fprintf(filep, fformat, Top->Ped[ped].Num);
                } else {
                    /* Write the pedigree id as a "fwid" numeric code */
                    prID_ped(filep, ped, fformat, &Top->Ped[ped]);
                }
                /* write the entry number  */
                if (analysis == TO_SAGE) {
                    fprintf(filep, pformat, entry + 1);
                    /* write the parents */
                    if (Entry->Father != 0) {
                        fprintf(filep, pformat, Entry->Father);
                        fprintf(filep, pformat, Entry->Mother);
                    } else {
                        for(b=0; b<2*(pwid + 1); b++)
                            fprintf(filep, " ");
                    }
                } else {
                    prID_fam(filep, pformat, Entry, Top->Ped[ped].Entry);
                    if (Entry->Father == 0) {
                        for(b=0; b<2*(pwid + 1); b++)
                            fprintf(filep, " ");
                    }
                }
                /* write the sex */
                if (Entry->Sex == MALE_ID)
                    fprintf(filep, "M ");
                else if (Entry->Sex == FEMALE_ID)
                    fprintf(filep, "F ");
                else
                    fprintf(filep, "? ");
                if (LoopOverTrait==1) {
                    /* first write the affection locus */
                    locus1 = *trp;
                    switch (Top->LocusTop->Locus[locus1].Type)   {
                    case QUANT:
                        if (analysis == TO_SAGE) {
                            sagewrite_quantitative_data(filep, locus1,
                                                        &(Top->LocusTop->Locus[locus1]),
                                                        Entry);
                        } else {
                            sage4write_quantitative_data(filep, locus1,
                                                         &(Top->LocusTop->Locus[locus1]),
                                                         Entry);
                        }
                        break;
                    case AFFECTION:
                        sagewrite_affection_data(filep, locus1,
                                                 &(Top->LocusTop->Locus[locus1]), Entry);
                        break;
                    default:
                        printf("INVALID locus type.\n");
                        break;
                    }
                }
                /* now write genotype data */
                for (l1 = 0; l1 < NumChrLoci; l1++)   {
                    locus1 = ChrLoci[l1];
                    switch (Top->LocusTop->Locus[locus1].Type)   {
                    case QUANT:
                        if (LoopOverTrait==0) {
                            if (analysis == TO_SAGE) {
                                sagewrite_quantitative_data(filep, locus1,
                                                            &(Top->LocusTop->Locus[locus1]),
                                                            Entry);
                            } else {
                                sage4write_quantitative_data(filep, locus1,
                                                             &(Top->LocusTop->Locus[locus1]),
                                                             Entry);
                            }
                        }
                        break;
                    case AFFECTION:
                        if (LoopOverTrait==0) {
                            sagewrite_affection_data(filep, locus1,
                                                     &(Top->LocusTop->Locus[locus1]),
                                                     Entry);
                        }
                        break;
                    case NUMBERED:
                        sagewrite_numbered_data(filep, locus1,
                                                &(Top->LocusTop->Locus[locus1]), Entry);
                        break;
                    case BINARY:
                        sagewrite_binary_data(filep, locus1,
                                              &(Top->LocusTop->Locus[locus1]), Entry);
                        break;
                    default:
                        printf("Unknown locus type (locus %d).\n",
                               locus1 + 1);
                    } /* end of switch on locus type */
                } /* End of loop locus1 to Top->LocusTop->LocusCnt */
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop==1) break;
        trp++;
    }
    return 1;
} /* End of save_SAGE_peds */


/*
 * Write the locus data file.
 */
static void  write_SAGE_locus_file(char *loutfl_name,
				   linkage_locus_top * LTop)
{
    int             num, locus1, allele;
    int             tr, nloop, num_affec=num_traits;
    char            locfl[2*FILENAME_LENGTH];
    FILE            *filep;
    linkage_locus_rec *Locus;

    num=0;
    for (locus1 = 0; locus1 < NumChrLoci; locus1++)  {
        if (LTop->Locus[ChrLoci[locus1]].Type == NUMBERED) {
            num++;
        }
    }

    if (num == 0) {
        warnf("No mapped loci: locus file not created.");
        return;
    }

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        /*    SKIP_TRI(tr-1); */
        sprintf(locfl, "%s/%s", output_paths[tr], loutfl_name);

        if ((filep=fopen(locfl, "w"))==NULL) {
            errorvf("Unable to open file '%s' for writing\n", locfl);
            EXIT(FILE_WRITE_ERROR);
        }

        for (locus1 = 0; locus1 < NumChrLoci; locus1++)  {
            Locus = &(LTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED) {
                if (Locus->LocusName != NULL)
                    fprintf(filep, "%-8s", Locus->LocusName);
                else
                    fprintf(filep, "%8d", locus1 + 1);
                fprintf(filep, "\n");
                for (allele = 0; allele < Locus->AlleleCnt; allele++)   {
#ifdef TEST                    
//                  fprintf(filep, "%-4d = %9.6f\n",
//                          allele + 1, Locus->Allele[allele].Frequency);
                    fprintf(filep, "%-4s = %9.6f\n",
                            format_allele(Locus, allele + 1),
                            Locus->Allele[allele].Frequency);
#else
//                  fprintf(filep, "%-4d = %8.5f\n",
//                          allele + 1, Locus->Allele[allele].Frequency);
                    fprintf(filep, "%-4s = %8.5f\n",
                            format_allele(Locus, allele + 1),
                            Locus->Allele[allele].Frequency);
#endif
                }
                fprintf(filep,";\n;\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;

    }
    sprintf(err_msg, "        Locus file:            %s", loutfl_name);
    mssgf(err_msg);
}

static void sage_file_names(char *file_names[], linkage_ped_top *Top)

{

    int i, select=-1;
    analysis_type analysis = TO_SAGE;
    char stem[5], fl_stat[12], shellfl_type[7], cselect[5];
    char filename[FILENAME_LENGTH], *name;

    ((num_traits > 1 || (LoopOverTrait == 1 && num_traits > 0))?
     strcpy(shellfl_type, "SIBPAL") :
     strcpy(shellfl_type, "FSP   "));

    if (main_chromocnt > 1) {
        // returns only the <file_name> removing <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, -9);
        strcpy(stem, "stem");
    }
    else
        strcpy(stem, "");
    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        select=0;
    }

    while(select != 0) {
        name = NULL;
        i=4;
        draw_line();
        print_outfile_mssg();
        printf("SAGE output file name selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt > 1) {
            printf(" 1) Pedigree file name stem:             %s\n", file_names[0]);
            printf(" 2) Locus file name stem:                %s\n", file_names[1]);
            printf(" 3) FSP parameter file name stem:        %s\n", file_names[2]);
            printf(" 4) %s C-shell file name stem:       %s\n",
                   shellfl_type, file_names[3]);
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                i++;
                printf(" %d) SIBPAL parameter file name stem:     %s\n", i,
                       file_names[7]);
                i++;
                if (Top->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                    printf(" %d) Genotype count file name stem:       %s\n", i,
                           file_names[6]);
                    i++;
                }
                printf(" %d) FCOR parameter file name stem:       %s\n", i,
                       file_names[10]);       i++;
                printf(" %d) FCOR C-shell file name stem:         %s\n", i,
                       file_names[11]);
            }
        } else {
            printf(" 1) Pedigree file name:             %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:                %-15s\t%s\n", file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 3) FSP parameter file name:        %-15s\t%s\n", file_names[2],
                   file_status(file_names[2], fl_stat));
            printf(" 4) %s C-shell file name:       %-15s\t%s\n",
                   shellfl_type, file_names[3], file_status(file_names[3], fl_stat));
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                i++;
                printf(" %d) SIBPAL parameter file name:     %-15s\t%s\n", i,
                       file_names[7], file_status(file_names[7], fl_stat));

                i++;
                if (Top->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                    printf(" %d) Genotype count file name:       %-15s\t%s\n", i,
                           file_names[6], file_status(file_names[6], fl_stat)); i++;
                }
                printf(" %d) FCOR parameter file name:       %-15s\t%s\n", i,
                       file_names[10],
                       file_status(file_names[10], fl_stat));        i++;
                printf(" %d) FCOR C-shell file name:         %-15s\t%s\n", i,
                       file_names[11],
                       file_status(file_names[11], fl_stat));
            }
        }

/*      individual_id_item(i, TO_SAGE, OrigIds[0], 43, 2,0, 0); */
/*      i++; */
/*      pedigree_id_item(i, TO_SAGE, OrigIds[1], 43, 2, 0); */

        printf("Enter option 0 - %d > ", i);
        fflush(stdout);
        IgnoreValue(fgets(cselect, 4, stdin)); newline;
        select=-1;
        sscanf(cselect, "%d", &select);
        test_modified(select);
        if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
            if (Top->LocusTop->Locus[global_trait_entries[0]].Type != AFFECTION
               && select > 5) {
                select++;
            }
        }
        switch(select) {
        case 0:
            break;
        case 1:
            strcpy(filename, "pedigree");  name=file_names[0];    break;
        case 2:
            strcpy(filename, "locus"); name=file_names[1];     break;
        case 3:
            strcpy(filename, "FSP parameter"); name=file_names[2];      break;
        case 4:
            sprintf(filename, "%s C-shell", shellfl_type); name=file_names[3];   break;
        case 5:
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                strcpy(filename, "SIBPAL parameter"); name=file_names[7];     break;
            }
        case 6:
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                if (Top->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                    strcpy(filename, "genotype count"); name=file_names[6];
                    break;
                }
            }
        case 7:
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                strcpy(filename, "FCOR parameter"); name=file_names[10];     break;
            }
        case 8:
            if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
                strcpy(filename, "FCOR C-shell"); name=file_names[11];     break;
            }
/*      case 9: */
/*        OrigIds[0] = individual_id_item(0, TO_SAGE, OrigIds[0], 35, 1, */
/*  				      Top->OrigIds, Top->UniqueIds); */
/*        break; */

/*      case 10: */
/*        OrigIds[1] = pedigree_id_item(0, TO_SAGE, OrigIds[1], 35, 1, */
/*  				    Top->OrigIds); */
/*        break; */

        default:
            printf("Unknown option %s\n", cselect);  select=-1;    break;
        }
        if (select > 0) {
//n
            if (name == NULL) {
                printf("Too few traits for selection number\n.  Select again.");
                newline;
            } else {
                printf("Enter NEW %s file name %s > ", filename, stem);
                fcmap(stdin, "%s", name); newline;
            }
        }
        strcpy(cselect,  "");
    }
/*    individual_id_item(0, TO_SAGE, OrigIds[0], 0, 3, */
/*  		     Top->OrigIds, Top->UniqueIds); */
/*    pedigree_id_item(0, TO_SAGE, OrigIds[1], 0, 3, Top->OrigIds); */
}

/* main routine to create sage files */

#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*Top)

void create_SAGE_file(linkage_ped_top **Top, char *mapfl_name,
		      int *trait,
		      ped_top *PedTreeTop, file_format *infl_type,
		      file_format *outfl_type,
		      int *numchr, analysis_type *analysis,
		      char *file_names[], int untyped_ped_opt)
{
    int  i, j, k, pwid, fwid;
    linkage_ped_tree *LPed;

    char loutfl_name[FILENAME_LENGTH], outfl_name[FILENAME_LENGTH];
    char poutfl_name[FILENAME_LENGTH], shfl3[FILENAME_LENGTH];
    char shfl11[FILENAME_LENGTH];
    /*  linkage_ped_top  *LPedTreeTop_copy; */
    linkage_ped_top *NewTop;
    /* int dummy; */
    int old_loop_over_trait, *global_trait_entries_old, num_traits_old;

    /* shall we recount the number of genotypes here? */
    /* Why?
       count_lgenotypes(LPedTreeTop,  &dummy, &dummy, &dummy, &dummy, &dummy,
       &dummy, &dummy, &dummy, &dummy); */
    sage_file_names(file_names, LPedTreeTop);

    field_widths(LPedTreeTop, NULL, &fwid, &pwid, NULL, NULL);
    /* set_missing_quant_input(LPedTreeTop);*/

    /*  LPedTreeTop_copy = new_lpedtop();
        copy_linkage_ped_top(LPedTreeTop, LPedTreeTop_copy, 1);
    */
    NewTop = new_lpedtop();
    /* make nuclear families */
    makenucs1(LPedTreeTop, NULL, NewTop);

    for (j=0; j< main_chromocnt; j++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[j];
            (*analysis)->replace_chr_number(file_names, *numchr);
/*       free_all_including_lpedtop(&LPedTreeTop); */
/*       LPedTreeTop=new_lpedtop(); */
/*       copy_linkage_ped_top(LPedTreeTop_copy, LPedTreeTop, 1); */
/*       ReOrderMappedLoci(LPedTreeTop, numchr); */
            sprintf(shfl3, "%s.sh", file_names[3]);
            sprintf(shfl11, "%s.sh", file_names[11]);

        } else {
            strcpy(shfl3, file_names[3]);
            strcpy(shfl11, file_names[11]);
        }
        get_loci_on_chromosome(*numchr);
        strcpy(loutfl_name, file_names[1]);
        strcpy(outfl_name,file_names[0]);
        strcpy(poutfl_name,file_names[6]);
        if (*infl_type == LINKAGE) {
            for (i = 0; i < LPedTreeTop->PedCnt; i++)  {
                LPed = &(LPedTreeTop->Ped[i]);
                if (LPed->Loops != NULL)   {
                    errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                    if (connect_loops(LPed, LPedTreeTop) < 0)  {
                        errorf("FATAL: Aborting.\n");
                        EXIT(LOOP_CONNECTION_ERROR);
                    }
                }
            }
        }
        /*
          SIPPAL files:
          file_names[0] = sage_ped
          file_names[1] = sage_loc
          file_names[2] = sage_par - fsp
          file_names[3] = sage.sh
          file_names[7] = sage_sibpal
          FCOR files used for two purposes
          (a) Correlations
          (b) Affected relative pair counts
          file_names[2] = sage_par - fsp
          file_names[6] = sage_cnt - rel-pair counts
          file_names[10] =  sage_cntpar - rel-pair counts
          file_names[11] = sage_cnt.sh - rel-pair counts and fcor
        */
        change_output_chr(file_names[3], *numchr);
        change_output_chr(file_names[11], *numchr);

        /* omit untyped peds */
        omit_peds(untyped_ped_opt, LPedTreeTop);
        /* write the SIBPAL related files */
        save_SAGE_peds(outfl_name, LPedTreeTop, TO_SAGE, pwid, fwid);
        write_SAGE_par(file_names[2], *numchr, LPedTreeTop, *analysis, pwid, fwid);
        write_SAGE_csh(shfl3, file_names[0], file_names[2],
                       file_names[7], file_names[1],
                       LPedTreeTop, NewTop->PedCnt);
        write_SAGE_locus_file(loutfl_name, LPedTreeTop->LocusTop);

        old_loop_over_trait = LoopOverTrait;
        if (num_traits > 1 || (num_traits == 1 && global_trait_entries[0] >= 0)) {

            write_SAGE_sibpal(file_names[7], *numchr, LPedTreeTop,
                              *analysis, pwid, fwid);
            /* write the FCOR related files */
            if (LoopOverTrait == 0 && num_traits > 2 &&
                LPedTreeTop->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                LoopOverTrait = 1;
                /* create a new global_trait_entries list */
                global_trait_entries_old = CALLOC((size_t) num_traits, int);
                for (i=0; i < num_traits; i++) {
                    global_trait_entries_old[i] = global_trait_entries[i];
                }
                num_traits_old = num_traits;
                /* modify the global_trait_entries list */
                i=0;
                /* This gets rid of the -1 item from the global_trait_entries list */
                while(global_trait_entries[i] >= 0) i++;
                for (k=i; k < num_traits - 1; k++) {
                    global_trait_entries[k] = global_trait_entries[k+1];
                }
                num_traits--;
                write_SAGE_cntcsh(shfl11, file_names[6], file_names[0],
                                  file_names[10], file_names[2],
                                  LPedTreeTop, NewTop->PedCnt, 1);
                write_SAGE_cnt_file(poutfl_name, LPedTreeTop, pwid, fwid);
                num_traits = num_traits_old;
                free(global_trait_entries);
                global_trait_entries = &(global_trait_entries_old[0]);
            } else {
                write_SAGE_cntcsh(shfl11, file_names[6], file_names[0],
                                  file_names[10], file_names[2],
                                  LPedTreeTop, NewTop->PedCnt, 0);
                if (LPedTreeTop->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                    write_SAGE_cnt_file(poutfl_name, LPedTreeTop, pwid, fwid);
                }
            }

            write_SAGE_cntpar(file_names[10], *numchr, LPedTreeTop, *analysis,
                              pwid, fwid);
        }
        LoopOverTrait = old_loop_over_trait;

        if (j==0) create_mssg(*analysis);
        draw_line();
        sprintf(err_msg, "        Pedigree file:         %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "        FSP parameter file:    %s", file_names[2]);
        mssgf(err_msg);
        sprintf(err_msg, "        Locus file:            %s", loutfl_name);
        mssgf(err_msg);

        if (num_traits > 1 || (num_traits == 1 && global_trait_entries[0] >= 0)) {
            sprintf(err_msg, "        SIBPAL C-shell file:   %s", shfl3);
            mssgf(err_msg);
            sprintf(err_msg, "        SIBPAL parameter file: %s", file_names[7]);
            mssgf(err_msg);
            if (LoopOverTrait == 0 && num_traits > 2) {
                mssgf("Sage files for each trait,  ");
            }
            if (LPedTreeTop->LocusTop->Locus[global_trait_entries[0]].Type == AFFECTION) {
                sprintf(err_msg, "        Genotype count file:   %s", file_names[6]);
                mssgf(err_msg);
            }
            sprintf(err_msg, "        FCOR parameter file:   %s", file_names[10]);
            mssgf(err_msg);
            sprintf(err_msg, "        FCOR C-shell file:     %s", shfl11);
            mssgf(err_msg);
        } else {
            sprintf(err_msg, "        FSP C-shell file:      %s", shfl3);
            mssgf(err_msg);
            warnf("FCOR and SIBPAL files not created (no traits were selected).");
        }
    }

    free_lpedtop_pedinfo(NewTop);

/*   free_all_including_lpedtop(&LPedTreeTop_copy); */
} /* End of create_SAGE_file */

#undef LPedTreeTop

/*************************************
	 SAGE 4.0 routines
**************************************/

static void sage4_file_names(char *file_names[], int has_orig, int has_uniq)

{

    int select=-1;
    analysis_type analysis = TO_SAGE4;
    char stem[5], fl_stat[12];
    char filename[FILENAME_LENGTH], *name = (char *)""; //silly compiler!

    if (main_chromocnt > 1) {
        // returns only the <file_name> removing <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, -9);
        strcpy(stem, "stem");
    }
    else
        strcpy(stem, "");

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        return;
    }

    while(select != 0) {
        draw_line();
        printf("SAGE 4.0 output file names selection menu:\n");
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt > 1) {
            printf(" 1) Pedigree file name stem:         %-15s\n", file_names[0]);
            printf(" 2) Locus file name stem:            %-15s\n", file_names[1]);
            printf(" 3) Parameter file name stem:        %-15s\n", file_names[6]);
            printf(" 4) Map file name stem:              %-15s\n", file_names[7]);
        } else {
            printf(" 1) Pedigree file name:             %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:                %-15s\t%s\n", file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 3) Parameter file name:            %-15s\t%s\n", file_names[6],
                   file_status(file_names[6], fl_stat));
            printf(" 4) Map file name:                  %-15s\t%s\n", file_names[7],
                   file_status(file_names[7], fl_stat));
        }

        individual_id_item(5, TO_SAGE4, OrigIds[0], 36, 2, has_orig, has_uniq);
        pedigree_id_item(6, TO_SAGE4, OrigIds[1], 36, 2, has_orig);

        printf("Enter option 0 - 6 > ");
        /*    printf("Enter option 0 - 4 > "); */
        fcmap(stdin, "%d", &select);
        test_modified(select);
        switch(select) {
        case 0:
            break;
        case 1:
            strcpy(filename, "pedigree");  name=file_names[0];    break;
        case 2:
            strcpy(filename, "locus"); name=file_names[1];     break;
        case 3:
            strcpy(filename, "patameter"); name=file_names[6];      break;
        case 4:
            strcpy(filename, "map"); name=file_names[7];     break;
        case 5:
            OrigIds[0] = individual_id_item(0, TO_SAGE4, OrigIds[0], 36, 1,
                                            has_orig, has_uniq);
            break;

        case 6:
            OrigIds[1] = pedigree_id_item(0, TO_SAGE4, OrigIds[1], 35, 1,
                                          has_orig);
            break;

        default:
            printf("Unknown option %d\n", select);      break;
        }
        if (select > 0 && select < 5) {
            printf("Enter NEW %s file name %s > ", filename, stem);
            fcmap(stdin, "%s", name); newline;
        }
    }
    /* log the id option */
    individual_id_item(0, TO_SAGE4, OrigIds[0], 0, 3,
                       has_orig, has_uniq);
    pedigree_id_item(0, TO_SAGE4, OrigIds[1], 0, 3, has_orig);
}

static void write_SAGE4_par(char *poutfl_name,
			    linkage_ped_top *Top)
{
    FILE *fp;
    int m, *trp;
    int tr, nloop, num_affec=num_traits;
    char pfl[2*FILENAME_LENGTH];
    int pwid, fwid;

    NLOOP;
    trp = &(global_trait_entries[0]);
    field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(pfl, "%s/%s",
                output_paths[tr], poutfl_name);

        if ((fp=fopen(pfl, "w"))==NULL) {
            errorvf("Unable to open file '%s' for writing\n", pfl);
            EXIT(FILE_WRITE_ERROR);
        }

        script_time_stamp(fp);
        fprintf(fp, "#pedigree file format\n");
        fprintf(fp, "pedigree,column{\n");
        fprintf(fp, "  # Field formats\n");
        /* find out the maximum length of the pedigree and id fields
         */
        /* format str looks like
           format="4A5,1X,A1,A3,1X,A1,A2,A3,1X"
        */
        fprintf(fp, "  format=\"A%d,1X,3(A%d,1X),A1,1X", fwid, pwid); /* pedigree data fields */
        if (LoopOverTrait==1) {
            /* write the trait first */
            m = *trp;
            switch (Top->LocusTop->Locus[m].Type)   {
            case QUANT:
                fprintf(fp, ",A7,1X"); /* quant phenotype */
                break;
            case AFFECTION:
                fprintf(fp, ",A1,1X"); /* affection 0 1 2 */
                break;
            default:
                printf("INVALID locus type.\n");
                break;
            }
        }
        /* now write markers */
        for (m = 0; m < NumChrLoci; m++)   {
            switch (Top->LocusTop->Locus[ChrLoci[m]].Type)   {
            case QUANT:
                if (LoopOverTrait==0) {
                    fprintf(fp, ",A7,1X"); /* quant phenotype */
                }
                break;
            case AFFECTION:
                if (LoopOverTrait==0) {
                    fprintf(fp, ",A1,1X"); /* affection 0 1 2 */
                }
                break;
            case BINARY:
            case NUMBERED:
                fprintf(fp, ",A6"); /* xxx/xxx */
                break;
            default:
                break;
            }
        }
        fprintf(fp, "\"\n"); /* end of format string */
        /* now for the field descriptions */
        fprintf(fp, "  # Field encoding parameters\n\n");
        fprintf(fp, "  individual_missing_value=\" \"\n");
        fprintf(fp, "  sex_code,male=\"M\",");
        fprintf(fp, "  female=\"F\",unknown=\"?\"\n\n");
        fprintf(fp, "  # Order in which fields are to be read.\n");

        fprintf(fp, "  # Family Structure\n");

        fprintf(fp, "  pedigree_id\n");
        fprintf(fp, "  individual_id\n");
        fprintf(fp, "  parent_id\n");
        fprintf(fp, "  parent_id\n");
        fprintf(fp, "  sex_field\n");
        fprintf(fp, "\n\n");
        fprintf(fp, "  # Phenotypes and markers\n");
        if (LoopOverTrait==1) {
            /* write the trait first */
            m = *trp;
            switch (Top->LocusTop->Locus[m].Type)   {
            case QUANT:
                // Missing values can be set in the parameter file....
                // Previous Mega2 code set the missing value to NA, but it is clear from
                // the SAGE code and documentation that a 'numeric code' is required for
                // quantitative traits. So we use the batch file value if specified.
                // Since the batch file item must be specified (or is set) if there are missing
                // values then we know that if it is not then there are no missing values...
                //fprintf(fp,"  trait=%s,missing=\"NA\"\n", Top->LocusTop->Pheno[*trp].TraitName);
                fprintf(fp,"  trait=%s", Top->LocusTop->Pheno[*trp].TraitName);
                if (ITEM_READ(Value_Missing_Quant_On_Output))
                    fprintf(fp,",missing=\"%s\"",
                            Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
                fprintf(fp,"\n");
                break;
            case AFFECTION:
                fprintf(fp, "  phenotype=%s,binary,",
                        Top->LocusTop->Pheno[*trp].TraitName);
                fprintf(fp, "affected=\"2\",");
                fprintf(fp, "unaffected=\"1\",");
                fprintf(fp, "missing=\" \"\n"); /* affection 0 1 2 */

                break;
            default:
                printf("INVALID locus type.\n");
                break;
            }
        }
        /* now write markers */
        for (m = 0; m < NumChrLoci; m++)   {
            int m1 = ChrLoci[m];
            switch (Top->LocusTop->Locus[m1].Type)   {
            case QUANT:
                if (LoopOverTrait==0) {
                    /* quant phenotype */
                    //fprintf(fp,"  trait=%s,missing=\"NA\"\n", Top->LocusTop->Pheno[*trp].TraitName);
                    fprintf(fp,"  trait=%s", Top->LocusTop->Pheno[m1].TraitName);
                    if (ITEM_READ(Value_Missing_Quant_On_Output)) {
                        fprintf(fp,",missing=\"%s\"",
                                Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
                    }
                    fprintf(fp,"\n");
                }
                break;
            case AFFECTION:
                if (LoopOverTrait==0) {
                    fprintf(fp, "  phenotype=%s,binary,",
                            Top->LocusTop->Locus[m1].LocusName);
                    fprintf(fp, "affected=\"2\",");
                    fprintf(fp, "unaffected=\"1\",");
                    fprintf(fp, "missing=\" \"\n"); /* affection 0 1 2 */
                }
                break;
            case BINARY:
            case NUMBERED:
                fprintf(fp, "  marker=%s",
                        Top->LocusTop->Locus[m1].LocusName);
                fprintf(fp, ",missing=\"0\"\n");
                break;
            default:
                printf("Unknown locus type (locus %d).\n",
                       m1 + 1);
            }  /* switch */
        } /* loop over loci */

        fprintf(fp, "}\n\n# analysis parameters\n\n");

        fclose(fp);
        if (nloop==1) break;
        trp++;
    }
}

static void   write_SAGE4_map(char *mapfl_name,
			      int numchr,
			      linkage_locus_top *LTop)
{
    FILE *fp;
    int i, m;
    double theta;
    int tr, nloop, num_affec=num_traits;
    char mfl[2*FILENAME_LENGTH];
    int *markers=CALLOC((size_t) NumChrLoci, int);

    NLOOP;

    i=0;
    for(m=0; m<NumChrLoci; m++) {
        if (LTop->Locus[ChrLoci[m]].Type == NUMBERED ||
           LTop->Locus[ChrLoci[m]].Type == BINARY) {
            markers[i]=ChrLoci[m]; i++;
        }
    }

    if (i == 0) {
        warnf("No mapped loci: map file not created");
        return;
    }
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(mfl, "%s/%s",
                output_paths[tr], mapfl_name);

        if ((fp=fopen(mfl, "w"))==NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n",  mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(fp, "genome, map=%s\n{\n",
                ((LTop->map_distance_type=='h')?
                 "Haldane": "Kosambi"));
        fprintf(fp,"  region=chr%d\n  {\n",  numchr);
        fprintf(fp, "    marker = %s\n",
                LTop->Marker[markers[0]].MarkerName);
        for(m=1; m<i; m++) {
            if ((LTop->Marker[markers[m]].pos_avg -
                LTop->Marker[markers[m-1]].pos_avg) < 0.0) {
                fprintf(fp, "    theta = 0.09999\n");
            } else {
                theta=((LTop->map_distance_type == 'h')?
                       haldane_theta(LTop->Marker[markers[m]].pos_avg -
                                     LTop->Marker[markers[m-1]].pos_avg):
                       kosambi_theta(LTop->Marker[markers[m]].pos_avg -
                                     LTop->Marker[markers[m-1]].pos_avg));
                fprintf(fp, "    theta = %0.7g\n", theta);
            }
            fprintf(fp, "    marker = %s\n",
                    LTop->Marker[markers[m]].MarkerName);
        }
        fprintf(fp, "  }\n}\n");
        fclose(fp);
        if (nloop==1) break;
    }
    sprintf(err_msg, "        Map file:              %s", mapfl_name);
    mssgf(err_msg);
}


#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*LPedTop)

void create_SAGE4_file(linkage_ped_top **LPedTop,
                       ped_top *PedTreeTop, file_format *infl_type,
                       int *numchr, analysis_type *analysis,
                       char *file_names[], int untyped_ped_opt)

{
    int  i, j, pwid, fwid;
    linkage_ped_tree *LPed;

    char loutfl_name[FILENAME_LENGTH], outfl_name[FILENAME_LENGTH];
    char mapfl_name[FILENAME_LENGTH], poutfl_name[FILENAME_LENGTH];

    linkage_ped_top  *LPedTreeTop_copy;

    sage4_file_names(file_names, LPedTreeTop->OrigIds, LPedTreeTop->UniqueIds);
    field_widths(LPedTreeTop, NULL, &fwid, &pwid, NULL, NULL);
    /* set_missing_quant_input(LPedTreeTop);*/

    LPedTreeTop_copy = new_lpedtop();
    copy_linkage_ped_top(LPedTreeTop, LPedTreeTop_copy, 1);

    for (j=0; j< main_chromocnt; j++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[j];
            (*analysis)->replace_chr_number(file_names, *numchr);
/*       free_all_including_lpedtop(&LPedTreeTop); */
/*       LPedTreeTop=new_lpedtop(); */
/*       copy_linkage_ped_top(LPedTreeTop_copy, LPedTreeTop, 1); */
/*       ReOrderMappedLoci(LPedTreeTop, numchr); */
        }
        get_loci_on_chromosome(*numchr);
        strcpy(outfl_name,file_names[0]);
        strcpy(loutfl_name, file_names[1]);
        strcpy(poutfl_name,file_names[6]);
        strcpy(mapfl_name,file_names[7]);

        if (*infl_type == LINKAGE) {
            for (i = 0; i < LPedTreeTop->PedCnt; i++)  {
                LPed = &(LPedTreeTop->Ped[i]);
                if (LPed->Loops != NULL)   {
                    errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                    if (connect_loops(LPed, LPedTreeTop) < 0)  {
                        errorf("FATAL: Aborting.");
                        EXIT(LOOP_CONNECTION_ERROR);
                    }
                }
            }
        }

        /* omit untyped peds */
        omit_peds(untyped_ped_opt, LPedTreeTop);
        save_SAGE_peds(outfl_name, LPedTreeTop, TO_SAGE4, pwid, fwid);
        write_SAGE4_par(poutfl_name, LPedTreeTop);
        if (j==0) create_mssg(*analysis);
        draw_line();
        sprintf(err_msg, "        Pedigree file:         %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "        Parameter file:        %s", file_names[6]);
        mssgf(err_msg);
        write_SAGE_locus_file(loutfl_name, LPedTreeTop->LocusTop);
/*      sprintf(err_msg, "        Locus file:            %s", loutfl_name);
        mssgf(err_msg); */
        write_SAGE4_map(mapfl_name, *numchr, LPedTreeTop->LocusTop);
    }
    if (SetMarkerPosToSpecial) {
        check_map_positions(LPedTreeTop, *numchr, 2);
        draw_line();
    }
    /* free_all_including_lpedtop(LPedTreeTop); */
    delink_linkage_ped_top(LPedTreeTop_copy);
    free_all_including_lpedtop(&LPedTreeTop_copy);
} /* End of create_SAGE4_file */
