/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  switch_map
            linkage_ext.h:  connect_loops get_loci_on_chromosome switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
*/


/*==================*/

static void  SOLARwrite_quantitative_data(FILE *filep, int locusnm,
					  linkage_ped_rec *entry)
{
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON) {
        fprintf(filep, "          ");
    } else {
        fprintf(filep, "%10.5f", entry->Pheno[locusnm].Quant);
    }
}


static void            SOLARwrite_affection_data(FILE *filep, int locusnm,
						 linkage_locus_rec *locus,
						 linkage_ped_rec *entry)

{
    if (entry->Pheno[locusnm].Affection.Status > 0)
        fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Status);
    if (locus->Pheno->Props.Affection.ClassCnt > 1)
        {
            if (entry->Pheno[locusnm].Affection.Status > 0)
                fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Class);
        }
}


static void            SOLARwrite_numbered_data(FILE *filep, int locusnm,
						linkage_ped_rec *entry)

{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);

    if (a1 == 0)
        fprintf(filep, " 0/ 0");
    else
        fprintf(filep, "%2d/%2d", a1, a2);
}

static void            SOLARwrite_numbered_xdata(FILE *filep, int locusnm,
                                                 linkage_ped_rec *entry)

{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);

    if (a1 == 0)
        fprintf(filep, "  / 0");
    else {
        if (a1 == a2) {
            fprintf(filep, "  /%2d", a2);
        } else {
            fprintf(filep, "%2d/%2d",  a1, a2);
        }
    }
}


/*
 * Save all the pedigrees to a file in SOLAR format. Return the
 * number of members saved (total).
 */
static int save_SOLAR_peds(char *fl_name, linkage_ped_top *Top)

{
#define LTOP Top->LocusTop
    int            tr, nloop, num_affec=num_traits;
    int            ped, entry;
    int hhid=-1;
    FILE *filep;
    char pedfl[2*FILENAME_LENGTH];
    register linkage_ped_rec *Entry;

    /* first figure out if there is a trait called HHID
     */

    for(tr=0; tr < num_traits; tr++) {
        SKIP_TRI(tr);
        if (!(strcmp(Top->LocusTop->Locus[global_trait_entries[tr]].Name, "HHID"))) {
            hhid=global_trait_entries[tr];
            break;
        }
    }

    NLOOP;

    for (tr=0; tr<=nloop; tr++) {
        if (nloop>1 && tr==0) continue;
        sprintf(pedfl, "%s/%s", output_paths[tr], fl_name);
        if ((filep = fopen(pedfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pedfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (hhid >=0) {
            fprintf(filep,"FAMID,ID,FA,MO,SEX,HHID\n");
        } else {
            fprintf(filep,"FAMID,ID,FA,MO,SEX\n");
        }
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* Write the pedigree id */
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(filep, "%s,", Top->Ped[ped].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(filep, "%d,", ped+1);
                } else {
                    fprintf(filep, "%d,", Top->Ped[ped].Num);
                }
                /* write the entry number  */
                /* write the parents */
                if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                    fprintf(filep, "%s,", Entry->UniqueID);
                    if (Entry->Father) {
                        fprintf(filep, "%s,%s,",
                                Top->Ped[ped].Entry[Entry->Father-1].UniqueID,
                                Top->Ped[ped].Entry[Entry->Mother-1].UniqueID);
                    } else {
                        fprintf(filep, "0,0,");
                    }
                } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    fprintf(filep, "%s,", Entry->OrigID);
                    if (Entry->Father) {
                        fprintf(filep, "%s,%s,",
                                Top->Ped[ped].Entry[Entry->Father-1].OrigID,
                                Top->Ped[ped].Entry[Entry->Mother-1].OrigID);
                    } else {
                        fprintf(filep, "0,0,");
                    }
                } else {
                    fprintf(filep, "%d,", Entry->ID);
                    fprintf(filep, "%d,%d,", Entry->Father, Entry->Mother);
                }

                /* write the sex */
                if (Entry->Sex == MALE_ID)
                    fprintf(filep, "M");
                else if (Entry->Sex == FEMALE_ID)
                    fprintf(filep, "F");
                else
                    fprintf(filep, "?");
                if (hhid >=0) {
                    fprintf(filep, ",");
                    SOLARwrite_quantitative_data(filep, hhid, Entry);
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop==1) break;
    }
    return 1;
} /* End of save_SOLAR_peds */


/*
 * Write the locus data file.
 */
static void  write_SOLAR_locus_file(char *fl_name, linkage_locus_top *LTop)
{
    int        tr, nloop, locus1, allele, num_affec=num_traits;
    FILE *filep;
    char lfl[2*FILENAME_LENGTH];
    register linkage_locus_rec *Locus;

    NLOOP;

    for (tr=0; tr<=nloop; tr++) {
        if (nloop>1 && tr==0) continue;
        sprintf(lfl, "%s/%s", output_paths[tr], fl_name);
        if ((filep = fopen(lfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", lfl);
            EXIT(FILE_WRITE_ERROR);
        }
        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            Locus = &(LTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED || Locus->Type == BINARY) {
                if (Locus->Name != NULL)
                    fprintf(filep, "%-8s ", Locus->Name);
                else
                    fprintf(filep, "%8d ", locus1 + 1);
                for (allele = 0; allele < Locus->AlleleCnt; allele++)   {
                    fprintf(filep, "%d %9f  ", allele + 1,
                            Locus->Allele[allele].Frequency);
                }
                fprintf(filep,"\n");
            }
        }
        fclose(filep);
        if (nloop==1) break;
    }
    return;
}

/* Looping over traits is implemented in here differently than
   all other routines. Since only traits are written, we
   either write all traits, or write one trait at a time
*/
static int write_SOLAR_pheno(char *fl_name, linkage_ped_top *Top)
{
    int             ped, entry, locus1;
    int      cov, *trp, tr, num_affec=num_traits;
    FILE     *filep;
    char     pfl[2*FILENAME_LENGTH];
    register linkage_ped_rec *Entry;
    register linkage_locus_rec *Locus;

    if (num_affec == 0) {
        printf("No trait loci, solar phenotype file not created.\n");
        return 0;
    }

    if ((num_affec == 1) &&
	((global_trait_entries[0] < 0) ||
	 (!strcmp(Top->LocusTop->Locus[global_trait_entries[0]].Name, "HHID")))) {
        printf("No trait loci, solar phenotype file not created.\n");
        return 0;
    }

    if (LoopOverTrait==0 || num_affec==1) {
        /* write all traits */
        sprintf(pfl, "%s/%s", output_paths[0], fl_name);
        if ((filep = fopen(pfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep,"FAMID,ID");
        for(tr=0; tr < num_affec; tr++) {
            SKIP_TRI(tr);
            locus1=global_trait_entries[tr];
            Locus = &(Top->LocusTop->Locus[locus1]);
            if (strcmp(Locus->Name, "HHID")) {
                fprintf(filep, ",%s", Locus->Name);
            }
        }
        for(cov=0; cov < num_covariates; cov++) {
            fprintf(filep, ",%s",
                    Top->LocusTop->Locus[covariates[cov]].Name);
        }

        fprintf(filep, "\n");
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* Write the pedigree id */
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(filep, "%s,", Top->Ped[ped].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(filep, "%d,", ped+1);
                } else {
                    fprintf(filep, "%d,", Top->Ped[ped].Num);
                }
                /* write the entry number  */
                if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                    fprintf(filep, "%s", Entry->UniqueID);
                } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    fprintf(filep, "%s", Entry->OrigID);
                } else {
                    fprintf(filep, "%d", Entry->ID);
                }

                /* now write genotype data  */
                for(tr=0; tr<num_affec; tr++) {
                    SKIP_TRI(tr);
                    locus1=global_trait_entries[tr];
                    switch (Top->LocusTop->Locus[locus1].Type)  {
                    case QUANT:
                        if (strcmp(Top->LocusTop->Locus[locus1].Name, "HHID")) {
                            fprintf(filep, ",");
                            SOLARwrite_quantitative_data(filep, locus1, Entry);
                        }
                        break;
                    case AFFECTION:
                        fprintf(filep, ",");
                        SOLARwrite_affection_data(filep, locus1,
                                                  &(Top->LocusTop->Locus[locus1]),
                                                  Entry);
                        break;
                    default:
                        break;
                    }
                }
                for(tr=0; tr < num_covariates; tr++) {
                    locus1=covariates[tr];
                    switch (Top->LocusTop->Locus[locus1].Type)  {
                    case QUANT:
                        fprintf(filep, ",");
                        SOLARwrite_quantitative_data(filep, locus1, Entry);
                        break;
                    case AFFECTION:
                        fprintf(filep, ",");
                        SOLARwrite_affection_data(filep, locus1,
                                                  &(Top->LocusTop->Locus[locus1]),
                                                  Entry);
                        break;
                    default:
                        break;
                    }
                }
                fprintf(filep,"\n");
            }
        }
        fclose(filep);
        return 1;
    }

    /* loopovertrait = 1  and num_traits > 1*/
    trp = &(global_trait_entries[0]);
    for (tr=1; tr<=num_affec; tr++) {
        sprintf(pfl, "%s/%s", output_paths[tr], fl_name);
        if ((filep = fopen(pfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", pfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep,"FAMID,ID");
        Locus = &(Top->LocusTop->Locus[*trp]);
        fprintf(filep, ",%s", Locus->Name);
        for(cov=0; cov < num_covariates; cov++) {
            fprintf(filep, ",%s",
                    Top->LocusTop->Locus[covariates[cov]].Name);
        }
        fprintf(filep, "\n");
        for (ped = 0; ped < Top->PedCnt; ped++) {
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);

                /* Write the pedigree id */
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(filep, "%s,", Top->Ped[ped].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(filep, "%d,", ped+1);
                } else {
                    fprintf(filep, "%d,", Top->Ped[ped].Num);
                }
                /* write the entry number  */
                if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                    fprintf(filep, "%s", Entry->UniqueID);
                } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    fprintf(filep, "%s", Entry->OrigID);
                } else {
                    fprintf(filep, "%d", Entry->ID);
                }

                /* now write genotype data  */
                switch (Locus->Type)  {
                case QUANT:
                    if (strcmp(Locus->Name, "HHID")) {
                        fprintf(filep, ",");
                        SOLARwrite_quantitative_data(filep, *trp, Entry);
                    }
                    break;
                case AFFECTION:
                    fprintf(filep, ",");
                    SOLARwrite_affection_data(filep, *trp, Locus, Entry);
                    break;

                default:
                    break;
                }
                for(cov=0; cov < num_covariates; cov++) {
                    Locus = &(Top->LocusTop->Locus[covariates[cov]]);
                    switch (Locus->Type)  {
                    case QUANT:
                        if (strcmp(Locus->Name, "HHID")) {
                            fprintf(filep, ",");
                            SOLARwrite_quantitative_data(filep,
                                                         covariates[cov],
                                                         Entry);
                        }
                        break;
                    case AFFECTION:
                        fprintf(filep, ",");
                        SOLARwrite_affection_data(filep,
                                                  covariates[cov],
                                                  Locus, Entry);
                        break;

                    default:
                        break;
                    }
                }
                fprintf(filep,"\n");
            }
        }
        fclose(filep);
        trp++;
    }
    return 1;
}

static int write_SOLAR_geno(char *flname, linkage_ped_top *Top, int sex_linked)
{
    int             ped, entry, locus1;
    int tr, nloop, num_affec=num_traits;
    char    gfl[2*FILENAME_LENGTH];
    FILE    *filep;
    register linkage_ped_rec *Entry;
    register linkage_locus_rec *Locus;


    NLOOP;

    for (tr=0; tr<=nloop; tr++) {
        if (nloop>1 && tr==0) continue;
        sprintf(gfl, "%s/%s", output_paths[tr], flname);
        if ((filep = fopen(gfl, "w")) == NULL) {
            errorvf("Unable to open file %s.\n", gfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep,"FAMID,ID");
        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            Locus = &(Top->LocusTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED)  {
                if (Locus->Name != NULL)
                    fprintf(filep, ",%s", Locus->Name);
                else
                    fprintf(filep, ",%d", locus1 + 1);
                /*     if (locus1 < Top->LocusTop->LocusCnt-1)
                       fprintf(filep,","); */
            }
        }
        fprintf(filep,"\n");
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* Write the pedigree id */
                if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                    fprintf(filep, "%s,", Top->Ped[ped].Name);
                } else if (OrigIds[1] == 3) {
                    fprintf(filep, "%d,", ped+1);
                } else {
                    fprintf(filep, "%d,", Top->Ped[ped].Num);
                }

                /* write the entry number  */
                if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                    fprintf(filep, "%s", Entry->UniqueID);
                } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    fprintf(filep, "%s", Entry->OrigID);
                } else {
                    fprintf(filep, "%d", Entry->ID);
                }

                /* now write genotype data  */
                for (locus1 = 0; locus1 < NumChrLoci; locus1++)  {
                    int l1=ChrLoci[locus1];
                    switch (Top->LocusTop->Locus[l1].Type)   {
                    case QUANT:
                        break;
                    case AFFECTION:
                        break;
                    case BINARY:
                    case NUMBERED:
                        fprintf(filep,",");
                        if (Top->LocusTop->SexLinked && IS_MALE(*Entry))
                            SOLARwrite_numbered_xdata(filep, l1, Entry);
                        else
                            SOLARwrite_numbered_data(filep, l1,  Entry);
                        break;
                    default:
                        fprintf(stderr, "unknown locus type (locus %d).\n", l1 + 1);
                    }
                }    /* End of loop locus1 to Top->LocusTop->LocusCnt */
                fprintf(filep,"\n");

            }
        }
        fclose(filep);
        if (nloop==1) break;
    }
    return 1;
}

static int write_SOLAR_map(char *fl_name, linkage_ped_top *Top, int chr)

{
    int  locus1, tr, nloop, num_affec=num_traits;
    double *position= NULL;
    FILE *filep;
    char mapfl[2*FILENAME_LENGTH];

    linkage_locus_rec *Locus;

    if (Top->LocusTop->map_distance_type == 'h') {
        position = CALLOC((size_t) NumChrLoci, double);
        switch_map(Top->LocusTop, position, SEX_AVERAGED_MAP);
    }

    NLOOP;

    for (tr=0; tr<=nloop; tr++) {
        if (nloop>1 && tr==0) continue;
        sprintf(mapfl, "%s/%s", output_paths[tr], fl_name);
        if ((filep = fopen(mapfl, "w")) == NULL) {
            errorvf("Unable to open file '%s' for writing\n", mapfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(filep,"%d\n", chr); /* Current chromosome number */
        for (locus1 = 0; locus1 < NumChrLoci; locus1++)  {
            Locus = &(Top->LocusTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED || Locus->Type == BINARY)  {
                fprintf(filep, "%10s  %10f\n", Locus->Name,
                        (Top->LocusTop->map_distance_type == 'h') ?
                        position[locus1] : Locus->Marker->pos_avg);
            }
        }
        fclose(filep);
        if (nloop == 1) break;
    }

    if (position != NULL) free(position);
    return 1;
}


static void solar_tcl_script(char *file_names[], int has_markers, int first_time,
			     int numchr, int sex_linked)
{

    /* Create a tck script with a general function that takes chromosome
       number as an argument, then create a bunch of functions that
       simply pass the chromosome number to this function */

    char tcl_fl[2*FILENAME_LENGTH];
    char file_stem[FILENAME_LENGTH];
    int  tr, nloop, num_affec=num_traits;
    FILE *filep;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(tcl_fl, "%s/%s", output_paths[tr], file_names[10]);
        if (first_time) {
            filep = fopen(tcl_fl, "w");
            script_time_stamp(filep);
            fprintf(filep, "proc solar_startup {chromo} {\n");
            fprintf(filep, " if {$chromo < 10 } {\n");
            fprintf(filep, "    if {$chromo == 0} {\n");
            fprintf(filep, "       set chr \"all\"\n");
            fprintf(filep, "    } else {\n");
            fprintf(filep, "     set chr \"0$chromo\"\n");
            fprintf(filep, "    }\n");
            fprintf(filep, " }  else {\n");
            fprintf(filep, "     set chr \"$chromo\"\n");
            fprintf(filep, " }\n");
            fprintf(filep, " puts \"Loading pedigrees\"\n");
            strcpy(file_stem, file_names[0]); change_output_chr(file_stem, -9);
            fprintf(filep, " load ped %s.$chr\n", file_stem);
            strcpy(file_stem, file_names[6]); change_output_chr(file_stem, -9);
            fprintf(filep, " load phenotypes %s.$chr\n", file_stem);
            if (has_markers) {
                strcpy(file_stem, file_names[9]); change_output_chr(file_stem, -9);
                fprintf(filep, " load map %s.$chr\n", file_stem);
                strcpy(file_stem, file_names[1]); change_output_chr(file_stem, -9);
                fprintf(filep, " load freq %s.$chr\n", file_stem);
                strcpy(file_stem, file_names[8]); change_output_chr(file_stem, -9);
                fprintf(filep, " load marker %s.$chr\n", file_stem);
                if (sex_linked) {
                    fprintf(filep, " ibdoption xlinked y\n");
                }
            }
            fprintf(filep, "}\n");
            fprintf(filep, "\n");
        } else {
            filep = fopen(tcl_fl, "a");
        }
        fprintf(filep, "proc load%d {} {\n", numchr);
        fprintf(filep, "   solar_startup %d\n", numchr);
        fprintf(filep, "}\n");
        fclose(filep);
        if (nloop == 1) break;
    }
}

static void get_solar_ou_filenames(analysis_type *anal, char *file_names[],
				   int has_markers, int has_orig,
				   int has_uniq)
{
    int i, choice=-1;
    char stem[5], fl_stat[12];

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        choice = 0;
    }
    if (main_chromocnt > 1) {
        strcpy(stem, "stem");
        (*anal)->replace_chr_number(file_names, -9);
    } else  {
        strcpy(stem, "");
    }

    /* output file name menu */
    while (choice != 0) {
        print_outfile_mssg();
        printf("Output file name change option:\n");
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt > 1) {
            printf(" 1) Pedigree file name stem          %s\n", file_names[0]);
            printf(" 2) Phenotype file name stem         %s\n", file_names[6]);
            if (has_markers) {
                printf(" 3) Locus file name stem             %s\n", file_names[1]);
                printf(" 4) Genotype file name stem          %s\n", file_names[8]);
                printf(" 5) Map file name stem               %s\n", file_names[9]);
                i=6;
            } else {
                i=3;
            }
        } else {
            printf(" 1) Pedigree file name           %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Phenotype file name          %-15s\t%s\n", file_names[6],
                   file_status(file_names[6], fl_stat));
            if (has_markers) {
                printf(" 3) Locus file name              %-15s\t%s\n", file_names[1],
                       file_status(file_names[1], fl_stat));
                printf(" 4) Genotype file name           %-15s\t%s\n", file_names[8],
                       file_status(file_names[8], fl_stat));
                printf(" 5) Map file name                %-15s\t%s\n", file_names[9],
                       file_status(file_names[9], fl_stat));
                i=6;
            } else {
                i=3;
            }
        }
        printf(" %d) TCL script name %s             %-15s\t%s\n", i,
               ((main_chromocnt >1)? "    ": ""),
               file_names[10], file_status(file_names[10], fl_stat));
        i++;
        individual_id_item(i, TO_SOLAR, OrigIds[0], 43, 2,0, 0);
        i++;
        pedigree_id_item(i, TO_SOLAR, OrigIds[1], 43, 2, 0);

        printf("Enter options 0-%d > ", i);
        fcmap(stdin, "%d", &choice); printf("\n");
        test_modified(choice);

        switch(choice) {
        case 0:
            break;
        case 1:
            printf("Enter new pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);    newline;
            break;
        case 2:
            printf("Enter new phenotype file name %s > ", stem);
            fcmap(stdin, "%s", file_names[6]);    newline;
            break;
        case 3:
            if (has_markers) {
                printf("Enter new locus file name %s > ", stem);
                fcmap(stdin, "%s", file_names[1]);    newline;
                break;
            }
            /* else fall through */
        case 4:
            if (has_markers) {
                printf("Enter new genotype file name %s > ", stem);
                fcmap(stdin, "%s", file_names[8]);    newline;
                break;
            }
            /* else fall through */
        case 5:
            if (has_markers) {
                printf("Enter new map file name %s > ", stem);
                fcmap(stdin, "%s", file_names[9]);    newline;
                break;
            }
            /* else fall through */
        case 6:
            printf("Enter new tcl script name > ");
            fcmap(stdin, "%s", file_names[10]); newline;
            break;
        case 7:
            OrigIds[0] = individual_id_item(0, TO_SOLAR, OrigIds[0], 35, 1,
                                            has_orig, has_uniq);
            break;
        case 8:
            OrigIds[1] = pedigree_id_item(0, TO_SOLAR, OrigIds[1], 35, 1,
                                          has_orig);
            break;
        default:
            printf("Unknown option %d\n", choice);
        }
    }

    /* log the id option */
    individual_id_item(0, TO_SOLAR, OrigIds[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, TO_SOLAR, OrigIds[1], 0, 3, has_orig);
    return;

}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void create_SOLAR_files(linkage_ped_top **LPedTop, ped_top *PedTreeTop,
			analysis_type *analysis,
			file_format *infl_type,
			file_format *outfl_type, int *numchr,
			char *mapfl_name, char *file_names[],
			int untyped_ped_opt)
{
    int xlinked, i;

    linkage_ped_tree *LPed;

    int has_markers = HasMarkers;

    if (main_chromocnt > 1) {
        change_output_chr(file_names[10], 0);
    } else {
        change_output_chr(file_names[10], *numchr);
    }
    get_solar_ou_filenames(analysis, file_names, has_markers,
                           Top->OrigIds, Top->UniqueIds);
    /* make a copy of the pedigree tree */

    if (*infl_type == LINKAGE) {
        for (i=0; i < Top->PedCnt; i++) {
            LPed=&(Top->Ped[i]);
            if (LPed->Loops != NULL)  {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)   {
                    errorf("FATAL: Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    /* set_missing_quant_input(Top);*/
    for (i=0; i < main_chromocnt; i++) {
        if (main_chromocnt > 1) {
            *numchr = global_chromo_entries[i];
            (*analysis)->replace_chr_number(file_names, *numchr);
        }
        xlinked =
            (((Top->LocusTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
              (Top->LocusTop->SexLinked == 1))? 1 : 0);
        get_loci_on_chromosome(*numchr);

        /* omit untyped peds */
        omit_peds(untyped_ped_opt, Top);
        write_SOLAR_pheno(file_names[6], Top);
        save_SOLAR_peds(file_names[0], Top);
        if (has_markers) {
            write_SOLAR_locus_file(file_names[1], Top->LocusTop);
            write_SOLAR_geno(file_names[8], Top, xlinked);
            write_SOLAR_map(file_names[9], Top, *numchr);
        }

        if (i==0) {
            solar_tcl_script(file_names, has_markers, 1, *numchr, xlinked);
        } else {
            solar_tcl_script(file_names, has_markers, 0, *numchr, xlinked);
        }

        if (i==0) create_mssg(*analysis);
        sprintf(err_msg, "        Pedigree file:         %s", file_names[0]);
        mssgf(err_msg);
        if ((num_traits > 1) || (LoopOverTrait == 1 && num_traits > 0)) {
            sprintf(err_msg, "        Phenotype file:        %s", file_names[6]);
            mssgf(err_msg);
        }

        if (has_markers) {
            sprintf(err_msg,
                    "        Allele frequency file: %s", file_names[1]);
            mssgf(err_msg);
            sprintf(err_msg,
                    "        Genotype file:         %s", file_names[8]);
            mssgf(err_msg);
            sprintf(err_msg,
                    "        Map file:              %s", file_names[9]);
            mssgf(err_msg);
        }

        if (i == (main_chromocnt - 1)) {
            sprintf(err_msg, "        Tcl script:            %s", file_names[10]);
            mssgf(err_msg);
        }

        if (xlinked) {
            printf("             ibdoption set to xlinked in %s.\n", file_names[10]);
        }
    }
} /* end of create_SOLAR_files */

#undef Top
