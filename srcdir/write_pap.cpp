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

/* write_pap.c
   create pedigree and genotype files for the pap program
   trip.dat, phen.dat and header.dat
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc
              fcmap_ext.h:  fcmap
            linkage_ext.h:  connect_loops get_loci_on_chromosome
           omit_ped_ext.h:  is_phenotyped is_typed_at_markers omit_peds
  output_file_names_ext.h:  create_mssg file_status print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
*/


typedef struct triplet_type_ {
    linkage_ped_rec *Entry, *Father, *Mother;
    int PedNum;
} triplet_type;
int compare_triplets(const void *trp1, const void *trp2);

/* pap always required numeric IDs */
void set_pap_unique_id(linkage_ped_top *Top)

{

    int ped, per;
/*   int pwid, fwid; */
    int  recordnum=1;
/*  int uniqueid; */

    OrigIds[0]=3;

    if (Top->IndivCnt > 100000000) {
        errorf("Too many individuals (> 10,000,000) to assign unique ids");
        errorf("Unique ids must fit within 8 digits.\n");
        EXIT(EXCEEDED_FIELD_MAX);
    } else {
        recordnum=1;
        /* Set the unique ids */
        for(ped=0; ped < Top->PedCnt; ped++) {
            for(per=0; per < Top->Ped[ped].EntryCnt; per++, recordnum++) {
                sprintf(Top->Ped[ped].Entry[per].UniqueID,
                        "%d", recordnum);
            }
        }
    }
}


/*   /\* field_widths(Top, NULL, &fwid, &pwid, NULL, NULL); *\/ */
/*   fwid--; pwid--; */

/*   if ((fwid+pwid) > 7) { */
/*     /\* Find out if there are too many records *\/ */
/*     for(ped=0; ped < Top->PedCnt; ped++) { */
/*       for(per=0; per < Top->Ped[ped].EntryCnt; per++, recordnum++) { */
/* 	; */
/*       } */
/*     } */
/*     if (recordnum > 100000000) { */
/*       errorf("Too many individuals (> 10,000,000) to assign unique ids"); */
/*       errorf("Unique ids must fit within 8 digits.\n"); */
/*       EXIT(-1); */
/*     } */
/*     else { */
/*       recordnum=1; */
/*       /\* Set the unique ids *\/ */
/*       for(ped=0; ped < Top->PedCnt; ped++) { */
/* 	for(per=0; per < Top->Ped[ped].EntryCnt; per++, recordnum++) { */
/* 	  sprintf(Top->Ped[ped].Entry[per].UniqueID, */
/* 		  "%8d", recordnum); */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   else { */
/*     for(ped=0; ped < Top->PedCnt; ped++) { */
/*       for(per=0; per < Top->Ped[ped].EntryCnt; per++) { */
/* 	uniqueid = */
/* 	  (int)pow(10.0, (double)(pwid+1))*Top->Ped[ped].Num +  */
/* 	  Top->Ped[ped].Entry[per].ID; */
/* 	sprintf(Top->Ped[ped].Entry[per].UniqueID, */
/* 		"%8d", uniqueid); */
/*       } */
/*     } */
/*   }        */
/*   return; */
/* } */

/* compare triplets:
   triplets have to be completely sorted, first on pedigree (which
   is true before this routine is called, then on father, then on
   mother, then on offspring.
*/

int compare_triplets(const void *trp1, const void *trp2)

{
    const triplet_type *trip1,  *trip2;

    trip1 = (const triplet_type *) trp1;
    trip2 = (const triplet_type *) trp2;

    if (atoi(trip1->Father->UniqueID) < atoi(trip2->Father->UniqueID)) {
        return -1;
    } else if (atoi(trip1->Father->UniqueID) == atoi(trip2->Father->UniqueID)) {
        if (atoi(trip1->Mother->UniqueID) < atoi(trip2->Mother->UniqueID)) {
            return -1;
        } else if (atoi(trip1->Mother->UniqueID) == atoi(trip2->Mother->UniqueID)) {
            if (atoi(trip1->Entry->UniqueID) < atoi(trip2->Entry->UniqueID)) {
                return -1;
            } else if (atoi(trip1->Entry->UniqueID) == atoi(trip2->Entry->UniqueID)) {
                return 0;
            } else {
                return 1;
            }
        } else {
            return 1;
        }
    }
    else
        return 1;
}

static void save_pap_peds(char *trip_file_name,
			  linkage_ped_top *Top)

{
    int tr, nloop, num_affec=num_traits;
    char tr_flname[2*FILENAME_LENGTH];
    int ped, per, rec, recordnum;
    linkage_ped_rec *Entry, *Father, *Mother;
    FILE *filep;
    triplet_type *triplets;
    int (*comp_trip)(const void *, const void *) = compare_triplets;

    recordnum=0;
    for (ped=0; ped < Top->PedCnt; ped++) {
        if (UntypedPeds != NULL && UntypedPeds[ped]) continue;

        recordnum += Top->Ped[ped].EntryCnt;
    }

    triplets = CALLOC((size_t) recordnum, triplet_type);
    recordnum=0;
    for(ped=0; ped < Top->PedCnt; ped++) {
        if (UntypedPeds != NULL && UntypedPeds[ped]) continue;

        for(per=0; per < Top->Ped[ped].EntryCnt; per++) {
            Entry = &(Top->Ped[ped].Entry[per]);
            if (Entry->Father != 0) {
                Father = &(Top->Ped[ped].Entry[Entry->Father-1]);
                Mother = &(Top->Ped[ped].Entry[Entry->Mother-1]);
                triplets[recordnum].Entry = Entry;
                triplets[recordnum].Father = Father;
                triplets[recordnum].Mother = Mother;
                triplets[recordnum].PedNum = Top->Ped[ped].Num;
                recordnum++;
            }
        }
    }

    /* Now do a quick sort */
    qsort(((void *) triplets), (size_t) recordnum, sizeof(triplet_type), comp_trip);

    NLOOP;

    for(tr=0; tr <= num_traits; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(tr_flname, "%s/%s", output_paths[tr], trip_file_name);
        if ((filep = fopen(tr_flname, "w")) == NULL) {
            errorvf("Unable to open file %s for writing.\n", tr_flname);
            EXIT(FILE_WRITE_ERROR);
        }

        // http://hasstedt.genetics.utah.edu/pap/papmanual70.pdf
        //
        // B.1. Pedigree Structure: trip.dat
        // (see section II.1)
        // Line repeated for each non-founder pedigree member.
        //   LPEDGR Pedigree number (columns: 1-8; type Integer)
        //   JFATH Father's ID number (columns: 9-16; type Integer)
        //   JMOTH Mother's ID number (columns: 17-24; type Integer)
        //   JCHLD Pedigree member's ID number (columns: 25-32; type Integer)

        for(rec=0; rec < recordnum; rec++) {
            fprintf(filep, "%8d%8s%8s%8s\n",
                    triplets[rec].PedNum,
                    triplets[rec].Father->UniqueID,
                    triplets[rec].Mother->UniqueID,
                    triplets[rec].Entry->UniqueID);
        }
        fclose(filep);
        if (nloop == 1) break;
    }
    free(triplets);
}


static void save_pap_phen(char *phen_file_name,
			  linkage_ped_top *Top,
			  int include_sex)
{
    int ped, per, loc, loc1, ncol;
    int tr, *trp, nloop, num_affec=num_traits;
    char phen_fl[2*FILENAME_LENGTH];
    entry_type comb_entry;
    int a1, a2;

    FILE *filep;

    NLOOP;

    if (num_affec > 0) {
        trp = &(global_trait_entries[0]);
    }
    else trp =NULL;

    for(tr=0; tr <= nloop; tr++) {
        if (tr==0 && nloop > 1) continue;

        sprintf(phen_fl, "%s/%s", output_paths[tr], phen_file_name);
        if ((filep = fopen(phen_fl, "w")) == NULL) {
            errorvf("Unable to open file %s for writing.\n", phen_fl);
            EXIT(FILE_WRITE_ERROR);
        }

        for(ped=0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL && UntypedPeds[ped]) continue;

            for(per=0; per < Top->Ped[ped].EntryCnt; per++) {
                comb_entry.LEntry = &(Top->Ped[ped].Entry[per]);
                if (is_typed_at_markers(comb_entry, 0, Top->LocusTop, NULL) ||
                    is_phenotyped(comb_entry, 0, Top->LocusTop)) {
                    
                    // B.4. Phenotype Data: phen.dat and ascer.dat
                    // (see section II.4
                    // Line repeated for each measured individual.
                    // ID ID number corresponding to trip.dat (columns: 1-8; type Integer)
                    // PHEN Individual's phenotype for column I, I = 1, NCOL (columns: 9-16,17-24,etc; type: Real)

                    fprintf(filep, "%8s", Top->Ped[ped].Entry[per].UniqueID);
                    ncol=1;

                    if (include_sex) {
                        fprintf(filep, "%8d", Top->Ped[ped].Entry[per].Sex);
                        ncol++;
                    }

                    if (LoopOverTrait == 1) {
                        /* print the tr-th trait locus */
                        switch(Top->LocusTop->Locus[*trp].Type) {
                        case AFFECTION:
                            fprintf(filep, "%8d",
                                    Top->Ped[ped].Entry[per].Pheno[*trp].Affection.Status);
                            break;

                        case QUANT:
                            {   double q = Top->Ped[ped].Entry[per].Pheno[*trp].Quant;
                                if (fabs(q - MissingQuant) <= EPSILON)
                                    fprintf(filep, "   -9999"); // default missing value code
                                else 
                                    fprintf(filep, "%8.2f", q);
                            }
                            break;
                        default:
                            break;
                        }
                        ncol++;
                    }
                    for(loc1=0; loc1 < NumChrLoci; loc1++) {
                        loc = ChrLoci[loc1];
                        if (ncol==10) {
                            fprintf(filep, "\n"); ncol=1;
                        }
                        switch(Top->LocusTop->Locus[loc].Type) {
                        case NUMBERED:
                        case BINARY:
                            get_2alleles(Top->Ped[ped].Entry[per].Marker, loc, &a1, &a2);
                            fprintf(filep, "%8s%8s",
                                    format_allele(&(Top->LocusTop->Locus[loc]), a1),
                                    format_allele(&(Top->LocusTop->Locus[loc]), a2));
                            ncol += 2;
                            break;

                        case QUANT:
                            if (LoopOverTrait == 0) {
                                double q = Top->Ped[ped].Entry[per].Pheno[loc].Quant;
                                if (fabs(q - MissingQuant) <= EPSILON)
                                    fprintf(filep, "   -9999"); // default missing value code
                                else 
                                    fprintf(filep, "%8.2f", q);
                                ncol++;
                            }
                            break;
                        case AFFECTION:
                            if (LoopOverTrait == 0) {
                                fprintf(filep, "%8d",
                                        Top->Ped[ped].Entry[per].Pheno[loc].Affection.Status);
                                ncol++;
                            }
                            break;
                        default:
                            break;
                        }
                    }
                    fprintf(filep, "\n");
                }
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }
}

//
// http://hasstedt.genetics.utah.edu/pap/papmanual70.pdf
// B.3. Trait/Marker Descriptions: header.dat
//
static int save_pap_hdr(char *hdr_file_name, linkage_locus_top *LocusTop,
                        int include_sex, int use_pop_file)
{
    int loc, num_vars, num_cols, ncol=0;
    int *trp, tr, nloop, num_affec=num_traits;
    int tmpi, num_loc = 0; //silly compiler
    char hdr_fl[2*FILENAME_LENGTH];
    char lname[9];
    FILE *filep;
    linkage_locus_rec *Locus;

    num_cols = ((LoopOverTrait == 0)?
                2*(NumChrLoci - (num_affec-1))+include_sex+num_affec-1 :
                2*(NumChrLoci - num_affec)+1+include_sex);
    num_vars = ((LoopOverTrait == 0)?
                NumChrLoci + include_sex :
                NumChrLoci - num_affec + 1 + include_sex);

    if (num_cols > 38) {
        sprintf(err_msg,
                "Too many variables (%d) needing %d columns in %s file.\n",
                num_vars, num_cols, hdr_file_name);
        errorf(err_msg);
        errorf("PAP allows only 39 columns including individual Id.\n");
        mssgf("Skipping this chromosome.");
        return -1;
    }

    NLOOP;

    if (num_affec > 0) {
        trp = &(global_trait_entries[0]);
    }

    for (tr=0; tr<=nloop; tr++) {
        if (tr==0 && nloop > 1) continue;
        sprintf(hdr_fl, "%s/%s", output_paths[tr], hdr_file_name);
        if ((filep = fopen(hdr_fl, "w"))==NULL) {
            errorvf("Unable to open file %s for writing.\n", hdr_fl);
            EXIT(FILE_WRITE_ERROR);
        }
        if (LoopOverTrait == 0) {
            num_loc=0;
        }

        //
        // Required line(s) describes the names of variables included in phen.dat (section II.4).
        //
        // NDATA      Number of variables in phen.dat (columns: 1-4; type: Integer)
        // NCOL       Number of columns in phen.dat (columns: 5-8; type: Integer)
        // NAME(I)    Name of variable I, I = 1, NDATA (columns: 9-16,17-24, etc; type: Char)
        //
        // NOTE: 9 items to a line, additional entries are on lines.
        fprintf(filep, "%4d%4d", num_vars, num_cols);
        ncol=0;
        if (include_sex) {
            fprintf(filep, "     sex"); ncol=1;
        }
        if (LoopOverTrait == 1) {
            Locus = &(LocusTop->Locus[*trp]);
            if (strlen(Locus->Name) <= 8)
                fprintf(filep, "%8s", Locus->Name);
            else {
                for (tmpi=0 ; tmpi < 8; tmpi++)
                    lname[tmpi]=Locus->Name[strlen(Locus->Name) - 8 + tmpi];
                lname[8]='\0';
                fprintf(filep, "%8s", lname);
                /*	    strcpy(Locus->Name, lname); */
            }
            ncol++;
        }

        for (loc=0; loc < NumChrLoci; loc++) {
            Locus = &(LocusTop->Locus[ChrLoci[loc]]);
            if (ncol == 9) {
                fprintf(filep, "\n"); ncol=1;
            }
            switch(Locus->Type) {
            case AFFECTION:
            case QUANT:
                if (LoopOverTrait == 0)  {
                    if (strlen(Locus->Name) <= 8)
                        fprintf(filep, "%8s", Locus->Name);
                    else {
                        for (tmpi=0 ; tmpi < 8; tmpi++)
                            lname[tmpi]=Locus->Name[strlen(Locus->Name) - 8 + tmpi];
                        lname[8]='\0';
                        fprintf(filep, "%8s", lname);
                        /*	    strcpy(Locus->Name, lname); */
                    }
                    ncol++;
                }
                break;

            case NUMBERED:
            case BINARY:
                if (strlen(Locus->Name) <= 8)
                    fprintf(filep, "%8s", Locus->Name);
                else {
                    for (tmpi=0 ; tmpi < 8; tmpi++)
                        lname[tmpi]=Locus->Name[strlen(Locus->Name) - 8 + tmpi];
                    lname[8]='\0';
                    fprintf(filep, "%8s", lname);
                    /*	    strcpy(Locus->Name, lname); */
                }
                ncol++;
                break;
            default:
                break;
            }
        }
        fprintf(filep, "\n");

        //
        // Required line(s) describes the types of variables included in phen.dat (section II.4).
        //
        // IVARTP(I) Type of variable I, I = 1, NDATA (columns: 9-12,17-20,etc; type: Integer)
        // = 1, discrete
        // = 2, quantitative
        // = 3, autosomal marker
        // = 4, X-linked marker
        // = 5, autosomal codominant marker specified as alleles
        // = 6, X-linked codominant marker specified as alleles
        //
        //  IVARDF(I) Definition of variable I, I = 1, NDATA (columns: 13-16,21-24,etc; type: Integer)
        // = 0, if IVARTP(I) = 1 and prevalence/incidence rates will not be used or if IVARTP(I) = 2
        // = location in popln.dat of prevalence/incidence rates if IVARTP(I) = 1
        // = (-) number of alleles, if IVARTP(I) = 3, 4, 5, or 6 and if no frequency/phenotype information in popln.dat
        // = location in popln.dat of frequency/phenotype information, if IVARTP(I) = 3, 4, 5, or 6
        //
        // NOTE: 9 items to a line, additional entries are on lines.

        ncol=0;
        // skip to the 9th column...
        fprintf(filep, "        ");
        if (include_sex) {
            fprintf(filep, "   1   0"); ncol += 8;
        }

        if (LoopOverTrait == 1) {
            switch(LocusTop->Locus[*trp].Type) {
            case AFFECTION:
                fprintf(filep, "   1"); // discrete
                fprintf(filep, "   0");
                break;
            case QUANT:
                fprintf(filep, "   2"); // quantitative
                fprintf(filep, "   0");
                break;
            default:
                break;
            }
            ncol += 8;
        }
        for (loc=0; loc < NumChrLoci; loc++) {
            switch(LocusTop->Locus[ChrLoci[loc]].Type) {
            case AFFECTION:
                if (LoopOverTrait == 0)	  {
                    fprintf(filep, "   1"); fprintf(filep, "   0");
                    ncol += 8;
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0)    {
                    fprintf(filep, "   2"); fprintf(filep, "   0");
                    ncol += 8;
                }
                break;
            case NUMBERED:
            case BINARY:
                if ((LocusTop->SexLinked == 2 &&
                    LocusTop->Marker[ChrLoci[loc]].chromosome == SEX_CHROMOSOME) ||
                   (LocusTop->SexLinked == 1)) {
                    fprintf(filep, "   6"); // X-linked
                } else {
                    fprintf(filep, "   5"); // autosomal
                }
                ncol += 4;
                if (use_pop_file) {
                    if (LoopOverTrait == 0) {
                        fprintf(filep, "%4d", num_loc+1); num_loc++;
                    } else {
                        fprintf(filep, "%4d", loc);
                    }
                } else {
                    fprintf(filep, "%4d", -LocusTop->Locus[ChrLoci[loc]].AlleleCnt);
                }
                ncol += 4;
                break;
            default:
                break;
            }
            if (ncol==72) {
                fprintf(filep, "\n"); ncol=0;
            }
        }
        if (ncol > 0) fprintf(filep, "\n");
        
        //
        // Optional line(s) specifying missing value code. MISVL(I) = -9999 if not included.
        //
        // -1 Code to indicate line of missing values (columns: 1-8; type: Integer)
        // MISVL(I) Missing value code for variable I, I = 1, NDATA (columns: 9-16,17-24,etc; type: Real)
        //
        // NOTE: 9 items to a line, additional entries are on lines.
        
        ncol=0;
        fprintf(filep, "      -1"); 	  ncol += 8;
        if (include_sex) {
            fprintf(filep, "       0"); 	  ncol += 8;
        }
        if (LoopOverTrait == 1) {
            if (LocusTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "       0");
            } else {
                // THIS IS A FIXED WITH FIELD!!!
               if (ITEM_READ(Value_Missing_Quant_On_Output))
                  fprintf(filep, "   %s",
                          Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
               else
                  fprintf(filep, "   -9999");
            }
            ncol += 8;

        }
        for (loc=0; loc < NumChrLoci; loc++) {
            switch(LocusTop->Locus[ChrLoci[loc]].Type) {
            case AFFECTION:
                if (LoopOverTrait == 0)    {
                    fprintf(filep, "       0"); 	  ncol += 8;
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0)    {
                    // THIS IS A FIXED WITH FIELD!!!
                    if (ITEM_READ(Value_Missing_Quant_On_Output))
                        fprintf(filep, "   %s",
                                Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
                    else
                        fprintf(filep, "   -9999");
                    ncol += 8;
                }
                break;
            case NUMBERED:
            case BINARY:
                fprintf(filep, "       0"); 	  ncol += 8;
                break;
            default:
                break;
            }
            if (ncol==72) {
                fprintf(filep, "\n"); ncol=0;
            }
        }
        if (ncol > 0) fprintf(filep, "\n");

        //
        // Optional line(s). PHSIM(I) = -99999 if not included.
        // -2 Code to indicate line of simulation indicators (columns: 1-8; type: Integer)
        // PHSIM(I) Simulation indicator for variable I, I = 1, NDATA (columns: 9-16,17-24,etc; type: Real)
        //
        // Optional line(s). MEAN(I) = 0 if not included.
        // -3 Code to indicate line of values to subtract (columns: 1-8; type: Integer)
        // MEAN(I) Value to subtract from variable I, I = 1, NDATA (columns: 9-16,17-24,etc; type: Real)
        //
        // Optional line(s). SD(I) = 1 if not included.
        // -4 Code to indicate line of values to divide by (columns: 1-8; type: Integer)
        // SD(I) Value by which to divide variable I, I = 1, NDATA 9-16,17-24,etc Real
        //
        // Optional line(s). PWR(I) = 1 if not included.
        // -5 Code to indicate line of power 1-8 Integer
        // PWR(I) for the transformation function y = r/P[(x/r+1)P-1] [Maclean et al 1976],
        // where x represents the phenotype, p represents the power and r = 6.

        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }

    return 0;
}


static void pap_popln_file(char *pop_file_name, linkage_locus_top *LTop)

{
    /* contains only markers for now */
    int tr, nloop, num_affec=num_traits;
    int loc, loc1, all;
    char pop_fl[2*FILENAME_LENGTH];
    FILE *filep;

    NLOOP;

    for(tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(pop_fl, "%s/%s", output_paths[tr], pop_file_name);
        if ((filep = fopen(pop_fl, "w"))==NULL) {
            errorvf("Unable to open file %s for writing.\n", pop_fl);
            EXIT(FILE_WRITE_ERROR);
        }
        for(loc1=0; loc1 < NumChrLoci; loc1++) {
            loc = ChrLoci[loc1];
            if (LTop->Locus[loc].Type == NUMBERED ||
               LTop->Locus[loc].Type == BINARY) {
                fprintf(filep, "   1 F %c\n",
                        (((LTop->SexLinked == 2 &&
                           LTop->Marker[loc].chromosome == SEX_CHROMOSOME) ||
                          (LTop->SexLinked == 1))? 'T' : 'F'));
                fprintf(filep, "%3d", LTop->Locus[loc].AlleleCnt);
                for(all=0; all < LTop->Locus[loc].AlleleCnt; all++) {
                    fprintf(filep, "%15f", LTop->Locus[loc].Allele[all].Frequency);
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
    }
}

static void pap_shell_file(char *file_names[])
{
    FILE *filep;
    char shfl[2*FILENAME_LENGTH];
    int tr, nloop, num_affec=num_traits;
    char *syscmd;
    NLOOP;

    for(tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(shfl, "%s/%s", output_paths[tr], file_names[4]);
        if ((filep = fopen(shfl, "w"))==NULL) {
            errorvf("Unable to open file %s for writing.\n", shfl);
            EXIT(FILE_WRITE_ERROR);
        }
        syscmd = CALLOC((strlen(shfl) + strlen("chmod +x ") + 1), char);

        fprintf(filep, "#!/bin/csh -f\n");
        fprintf(filep, "# C-shell file name %s\n", file_names[4]);
        script_time_stamp(filep);
        fprintf(filep, "#\n");
        fprintf(filep, "unalias cp\n");
        fprintf(filep, "cp %s popln.dat\n", file_names[3]);
        fprintf(filep, "cp %s phen.dat\n", file_names[2]);
        fprintf(filep, "cp %s header.dat\n", file_names[1]);
        fprintf(filep, "cp %s trip.dat\n", file_names[0]);
        fclose(filep);
        sprintf(syscmd, "chmod +x %s", shfl);
        System(syscmd);
        free(syscmd);
        if (nloop == 1) break;
    }
}

static void pap_file_names(char *file_names[], int *include_sex,
			   int *use_pop_file)
{
    int i=6;
    char cselect[5], fl_stat[12], file_type[11];
    analysis_type analysis = TO_PAP;

    *include_sex=1;
    *use_pop_file = 1;
    if (main_chromocnt > 1) {
        // Here -9 is a flag to return only the <file_name> removing <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, -9);
    }

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults");
        i=0;
    }

    while(i) {
        draw_line();
        print_outfile_mssg();
        printf("PAP output file name selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt <= 1) {
            printf(" 1) Pedigree file name:         %-15s\t%s\n",
                   file_names[0], file_status(file_names[0], fl_stat));
            printf(" 2) Header file name:           %-15s\t%s\n",
                   file_names[1], file_status(file_names[1], fl_stat));
            printf(" 3) Phenotypes file name:       %-15s\t%s\n",
                   file_names[2], file_status(file_names[2], fl_stat));
            printf(" 4) Population info file name:  %-15s\t%s\n",
                   file_names[3], file_status(file_names[3], fl_stat));
            printf(" 5) C-Shell script file name:   %-15s\t%s\n",
                   file_names[4], file_status(file_names[4], fl_stat));
        } else {
            printf(" 1) Pedigree file name stem:         %-15s\n",
                   file_names[0]);
            printf(" 2) Header file name stem:           %-15s\n",
                   file_names[1]);
            printf(" 3) Phenotype file name stem:        %-15s\n",
                   file_names[2]);
            printf(" 4) Population info file name stem:  %-15s\n",
                   file_names[3]);
            printf(" 5) C-shell script file name stem:   %-15s\n",
                   file_names[4]);
        }
        printf(" 6) Include Sex as a phenotype:      [%s]\n",
               yorn[*include_sex]);
        printf(" 7) Read allele frequencies from popln file: [%s]\n",
               yorn[*use_pop_file]);
        printf("Select options 1-5 to enter new file names, 6 and 7 to toggle > ");
        fflush(stdout);
        (void)fgets(cselect, 4, stdin); newline;
        sscanf(cselect, "%d", &i);
        test_modified(i);

        if (i < 6 && i > 0) {
            switch(i) {
            case 1:
                sprintf(file_type, "Pedigree");
                break;
            case 2:
                sprintf(file_type, "Header");
                break;
            case 3:
                sprintf(file_type, "Phenotype");
                break;
            case 4:
                sprintf(file_type, "Population");
                break;
            case 5:
                sprintf(file_type, "C-shell");
                break;
            }
            printf("Enter new %s file name > ", file_type);
            fcmap(stdin, "%s", file_names[i-1]); newline;
        } else if (i == 6) {
            *include_sex = TOGGLE(*include_sex);
        } else if (i == 7) {
            *use_pop_file = TOGGLE(*use_pop_file);
        } else if (i != 0) {
            warn_unknown(cselect);
        }
    }

    // If more than one chromosome is specified, then just used the first one.
    // There is no provision for using multiple chromosomes in this analysis option.
    if (main_chromocnt > 1) {
        analysis->replace_chr_number(file_names, global_chromo_entries[0]);
        strcat(file_names[4], ".sh");
    }
}

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

//
// http://hasstedt.genetics.utah.edu/pap/
//
// PAP is a package of programs for computing likelihoods or simulating phenotypes on pedigree
// members using genetic models. The package may be used for segregation analysis,
// variance components analysis, linkage analysis, measured genotype analysis, transmission
// disequilibrium testing, or genetic model fitting. It is distributed as the FORTRAN source.
// A large suite of examples, keyed to the manual, is included.

void create_pap_files(linkage_ped_top **LPedTop, char *file_names[],
		      int untyped_ped_opt, int *numchr,
		      file_format *infl_type)
{
    int ped, i, include_sex, use_pop_file;
    analysis_type analysis = TO_PAP;
    linkage_ped_tree *LPed;

    if (*infl_type == LINKAGE) {
        for (ped = 0; ped < Top->PedCnt; ped++)  {
            LPed = &(Top->Ped[ped]);
            if (LPed->Loops != NULL)   {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)  {
                    errorvf("FATAL: Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    pap_file_names(file_names, &include_sex, &use_pop_file);
/*   set_unique_id(Top); */

/*   if (main_chromocnt > 1) { */
/*     Copy = new_lpedtop(); */
/*     copy_linkage_ped_top(Top, Copy, 1); */
/*   } */

    for (i=0; i < main_chromocnt; i++) {
        if (main_chromocnt > 1) {
            *numchr = global_chromo_entries[i];
            analysis->replace_chr_number(file_names, *numchr);
/*       free_all_including_lpedtop(&Top); */
/*       Top=new_lpedtop(); */
/*       copy_linkage_ped_top(Copy, Top, 1); */
/*       ReOrderMappedLoci(Top, numchr); */
        }
        get_loci_on_chromosome(*numchr);
        if (i==0) {
            create_mssg(analysis);
        }
        omit_peds(untyped_ped_opt, Top);
        if (save_pap_hdr(file_names[1], Top->LocusTop, include_sex, use_pop_file)==0) {
            save_pap_peds(file_names[0], Top);
            save_pap_phen(file_names[2], Top, include_sex);
            pap_popln_file(file_names[3], Top->LocusTop);
            pap_shell_file(file_names);
            sprintf(err_msg, "\tPedigree file:   %s", file_names[0]);
            mssgf(err_msg);
            sprintf(err_msg, "\tHeader file:     %s", file_names[1]);
            mssgf(err_msg);
            sprintf(err_msg, "\tPhenotype file:  %s", file_names[2]);
            mssgf(err_msg);
            sprintf(err_msg, "\tPopulation file: %s", file_names[3]);
            mssgf(err_msg);
            sprintf(err_msg, "\tC-shell file:    %s", file_names[4]);
            mssgf(err_msg);
            draw_line();
        }
    }

/*   if (main_chromocnt > 1) { */
/*     free_all_including_lpedtop(&Copy); */
/*   } */
}

#undef Top
