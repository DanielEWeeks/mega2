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
#include <math.h>

#include "common.h"
#include "typedefs.h"
#include "write_mfiles.h"

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "genetic_utils_ext.h"
#include "utils_ext.h"
#include "user_input_ext.h"
#define WRITE_MENDEL7_FILES_EXT_H
#include "write_mendel7_files_ext.h"

/*
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  mssgf my_calloc warnf
      genetic_utils_ext.h:  chrom_num_to_name delete_file haldane_theta kosambi_theta
              utils_ext.h:  EXIT
*/


void csv_save_mendel_peds(char *file_name, linkage_ped_top *Top);

void csv_write_mendel_locus_file(char *file_name,
				 linkage_locus_top *LTop,
				 int sex_linked);
void csv_mendel7_pen_file(char *fl_name,
			  linkage_ped_top *Top,
                          int sex_linked);

void csv_write_mendel_map(char *mapfile, linkage_locus_top *LTop,
			  analysis_type analysis, int xlinked);

/* modified frm save_mendel_peds,
   don't need the number of affection  loci, or set formats,
   or shorten names.
   Pedigree name is written on every line.
   Every missing value is a null string (nothing between commas)
   Write the qtl variables at the end.
*/

void csv_save_mendel_peds(char *outfl_name, linkage_ped_top *Top)

{
    int      nloop, tr, ped, entry, locus1, *trp;
    int      /*num_labels, *labels, */ num_affec = num_traits;
    register linkage_ped_rec *Entry;
    int      affected;
    int      a1, a2;

    char     ofl[2*FILENAME_LENGTH];
    FILE     *filep;

    NLOOP;

    if (num_affec > 0) {
        trp = &(global_trait_entries[0]);
    } else {
        trp=NULL;
    }

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;
/*
        if (trp != NULL) {
            num_labels=NumLabels;
            labels = Labels;
        } else {
            labels=NULL; num_labels=0;
        }
*/
        sprintf(ofl, "%s/%s", output_paths[tr], outfl_name);
        if ((filep = fopen(ofl, "a")) == NULL) {
            errorvf("Could not open file %s for writing.\n", ofl);
            EXIT(FILE_WRITE_ERROR);
        }

        // for each pedigree....
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            // for each individual/person in that pedigree...
	    // see: http://www.genetics.ucla.edu/software/Mendel_current_doc.pdf
	    // page#35 (The Pedigree File).
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
                    // The individual's name...
                    fprintf(filep, "%s,", Entry->UniqueID);
                    // Name of parent #1, parent #2
                    if (Entry->Father) {
                        fprintf(filep, "%s,%s,",
                                Top->Ped[ped].Entry[Entry->Father-1].UniqueID,
                                Top->Ped[ped].Entry[Entry->Mother-1].UniqueID);
                    } else {
                        fprintf(filep, ",,");
                    }
                } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                    // The individual's name...
                    fprintf(filep, "%s,", Entry->OrigID);
                    // Name of parent #1, parent #2
                    if (Entry->Father) {
                        fprintf(filep, "%s,%s,",
                                Top->Ped[ped].Entry[Entry->Father-1].OrigID,
                                Top->Ped[ped].Entry[Entry->Mother-1].OrigID);
                    } else {
                        fprintf(filep, ",,");
                    }
                } else {
                    fprintf(filep, "%d,", Entry->ID);
                    if (Entry->Father > 0) {
                        fprintf(filep, "%d,%d,", Entry->Father, Entry->Mother);
                    } else {
                        fprintf(filep, ",,");
                    }
                }

                // Individual's sex...
                if (Entry->Sex == MALE_ID)
                    fprintf(filep, "M,");
                else if (Entry->Sex == FEMALE_ID)
                    fprintf(filep, "F,");
                else
                    fprintf(filep, "?,");

		// Where is Identical twin status????

                /* if LoopOverTrait mode, then write the single trait locus */
                if (LoopOverTrait == 1 && num_affec > 0) {
                    switch(Top->LocusTop->Locus[*trp].Type) {
                    case AFFECTION:
                        fprintf(filep,",");
                        affected=aff_status_entry(Entry->Pheno[*trp].Affection.Status, Entry->Pheno[*trp].Affection.Class,
                                                  &(Top->LocusTop->Locus[*trp]));
/* 	    aff_status= Entry->Pheno[*trp].Affection.Class+ liability_multiplier*Entry->Pheno[*trp].Affection.Status; */
/* 	    for (i=0; i<num_labels; i++) { */
/* 	      if (aff_status==labels[i]) { */
/* 		affected=2; break; */
/* 	      } */
/* 	    } */
/* 	    if (affected != 2) { */
/* 	      if (Entry->Pheno[*trp].Affection.Status == 0) { */
/* 		affected = 0; */
/* 	      } */
/* 	      else { */
/* 		affected = 1; */
/* 	      } */
/* 	    } */
                        if (affected > 0) {
                            fprintf(filep, "%d", affected);
                        }
                        break;
                    default:
                        break;
                    }
                }

                /* now write genotype data  */
                for (locus1 = 0; locus1 < NumChrLoci; locus1++)  {
                    int l1=ChrLoci[locus1];
                    switch (Top->LocusTop->Locus[l1].Type)   {
                    case QUANT:
                        if (LoopOverTrait == 0) {
                            fprintf(filep,",");
                            if (fabs(Entry->Pheno[l1].Quant - MissingQuant) > EPSILON) {
                                fprintf(filep, "%8f ", Entry->Pheno[l1].Quant);
                            }
                        }
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            fprintf(filep,",");
                            affected=aff_status_entry(Entry->Pheno[l1].Affection.Status, Entry->Pheno[l1].Affection.Class,
                                                      &(Top->LocusTop->Locus[l1]));
/* 	      aff_status= Entry->Pheno[l1].Affection.Class+  */
/* 		liability_multiplier*Entry->Pheno[l1].Affection.Status; */
/* 	      for (i=0; i<num_labels; i++) { */
/* 		if (aff_status==labels[i]) { */
/* 		  affected=2; break; */
/* 		} */
/* 	      } */
/* 	      if (affected != 2) { */
/* 		if (Entry->Pheno[l1].Affection.Status == 0) { */
/* 		  affected = 0; */
/* 		} */
/* 		else { */
/* 		  affected = 1; */
/* 		} */
/* 	      }	     */
                            if (affected > 0) {
                                fprintf(filep, "%d", affected);
                            }
                        }
                        break;
                    case BINARY:
                    case NUMBERED:
                        fprintf(filep,",");
                        get_2alleles(Entry->Marker, l1, &a1, &a2);
                        if (a1 > 0 && a2 > 0) {
                            fprintf(filep, "%d/%d", a1, a2);
                        } else if (a1 > 0 && a2 == 0) {
                            fprintf(filep, "%d/", a1);
                        } else if (a1 == 0 && a2 > 0) {
                            fprintf(filep, "%d/", a2);
                        }

                        break;
                    default:
                        fprintf(stderr, "unknown locus type (locus %d).\n", l1 + 1);
                    }
                }
                if (LoopOverTrait == 1 && num_affec > 0) {
                    switch(Top->LocusTop->Locus[*trp].Type) {
                    case QUANT:
                        fprintf(filep,",");
                        if (fabs(Entry->Pheno[*trp].Quant - MissingQuant) > EPSILON) {
                            fprintf(filep, "%8f ", Entry->Pheno[*trp].Quant);
                        }
                        break;
                    default:
                        break;
                    }
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop==1) break;
        trp++;
    }

    sprintf(err_msg, "        Pedigree file:        %s", outfl_name);
    mssgf(err_msg);
    return;

} /* End of csv_save_mendel_peds */

/* Write the qtl variables at the end */

void  csv_write_mendel_locus_file(char *file_name,
				  linkage_locus_top *LTop,
				  int sex_linked)
{

    int      nloop, tr, locus1, allele, tmpi, tmpi2, *trp;
    register linkage_locus_rec *Locus;
    char     chromo[4], loutfl_name[2*FILENAME_LENGTH];
    int      numgen, ipen, num_affec = num_traits;
    FILE     *fp;
    double pen;

    NLOOP;

    if (num_affec > 0) trp = &(global_trait_entries[0]);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(loutfl_name, "%s/%s", output_paths[tr], file_name);
        if ((fp = fopen(loutfl_name, "a")) == NULL) {
            errorvf("Could not open locus file %s.\n", loutfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        /* first write the <tr>th affection locus */
        if (LoopOverTrait == 1 && num_affec > 0) {
            /* this pointer is advanced only if there are multiple
               traits being looped over. */
            Locus = &(LTop->Locus[*trp]);
            if (Locus->Type == AFFECTION) {
                if (Locus->Pheno->Props.Affection.ClassCnt > 1) {
                    MENDEL7_ML_AFFECTION
                        } else {
                    MENDEL7_AFFECTION
                        }
            }
        }

        /* Now write the numbered and binary loci */
        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            Locus = &(LTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED || Locus->Type == BINARY) {
                if (Locus->Name != NULL) {
		  fprintf(fp, "%s,", strtail(Locus->Name, MENDEL7_MAX_LOCUS_NAME_LEN));
                }
                else
                    fprintf(fp, "%d,", locus1 + 1);
                fprintf(fp, "%s, %d\n",	chrom_num_to_name(Locus->Marker->chromosome, &(chromo[0])),
                        ((Locus->Type == NUMBERED)?
                         Locus->AlleleCnt : Locus->Marker->Props.Binary.FactorCnt));
            }
            switch (Locus->Type)  {
            case BINARY:
                for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++)  {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(fp, " %d", (int) Locus->Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', fp);
                }
                break;
            case NUMBERED:
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(fp, "%d,%8f\n", allele + 1,
                            Locus->Allele[allele].Frequency);
                }
                break;
            case AFFECTION:
                if (LoopOverTrait == 0) {
                    if (Locus->Pheno->Props.Affection.ClassCnt > 1) {
                        MENDEL7_ML_AFFECTION
                            } else {
                        MENDEL7_AFFECTION
                            }
                }
                break;
            case QUANT:
                if (LoopOverTrait == 0) {
		  fprintf(fp, "%s,variable,0\n", strtail(Locus->Name, MENDEL7_MAX_VARIABLE_NAME_LEN));
                }
                break;
            default:
                errorvf("Unknown locus type %d\n", (int) Locus->Type);
                EXIT(DATA_TYPE_ERROR);
                break;
            }
        }
        if (LoopOverTrait == 1 && num_affec > 0) {
            /* this pointer is advanced only if there are multiple
               traits being looped over. */
            Locus = &(LTop->Locus[*trp]);
            if (Locus->Type == QUANT) {
	      fprintf(fp, "%s,variable,0\n", strtail(LTop->Locus[*trp].Name, MENDEL7_MAX_VARIABLE_NAME_LEN));
            }
        }
        trp++;

        fclose(fp);
        if (nloop == 1) break;
    }

    sprintf(err_msg, "        Definition file:      %s", file_name);
    mssgf(err_msg);
    return;

}
/*------------------ End of write_mendel_definition_file------------ */

/* csv mendel7 penetrance file,
   In loop-over mode, separate files are created with one trait in each
   In combine traits mode, penetrance records for each traits are written
   sequentially, i.e. all individuals for trait-1, then all individuals
   for trait-2 etc. */

/* NM - Dec 01 2008 , call this function only if there are multiple
   liability classes */

void csv_mendel7_pen_file(char *fl_name, linkage_ped_top *Top, int sex_linked)

{
    int       tr;
    register  linkage_locus_rec *Locus;
    int       trp, ipen, ped, entry;
    FILE      *filep;
    char poutfl_name[2*FILENAME_LENGTH];
    char genostr[3][4] = {"1/1", "1/2", "2/2"};
    double pen;
    linkage_ped_rec *Entry;
    linkage_affection_class *lacp;
    int      sex;
    int      affected;

    for (tr=0; tr <= num_traits; tr++) {
        if (LoopOverTrait == 0 || num_traits == 1) {
            if (tr == 0) {
                sprintf(poutfl_name, "%s/%s", output_paths[0], fl_name);
                if ((filep = fopen(poutfl_name, "w")) == NULL) {
                    errorvf("Could not open file %s for writing.\n", poutfl_name);
                    EXIT(FILE_WRITE_ERROR);
                }
            }
            trp = global_trait_entries[tr];
            if (trp == -1) continue;
            if (tr == num_traits) break;
            if (Top->LocusTop->Locus[trp].Type == QUANT) continue;
/*       if (Top->LocusTop->Locus[trp].Type == AFFECTION && */
/* 	  Top->LocusTop->Pheno[trp].Props.Affection.ClassCnt == 1) continue; */
            /* file is already open for writing */

        } else {
            if (tr == 0) continue;
            trp = global_trait_entries[tr-1];
            if (trp == -1) continue;
            if (Top->LocusTop->Locus[trp].Type == QUANT) continue;
            if (Top->LocusTop->Locus[trp].Type == AFFECTION &&
                Top->LocusTop->Pheno[trp].Props.Affection.ClassCnt == 1) continue;

            sprintf(poutfl_name, "%s/%s", output_paths[tr], fl_name);
            if ((filep = fopen(poutfl_name, "w")) == NULL) {
                errorvf("Could not open file %s for writing.\n", poutfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
        }
        Locus = &(Top->LocusTop->Locus[trp]);
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }

            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);
                lacp  = &(Locus->Pheno->Props.Affection.Class[Entry->Pheno[trp].Affection.Class-1]);
                sex   = Entry->Sex;
                /* don't write the unknown penetrances */
                if (Entry->Pheno[trp].Affection.Status == 0) {
                    continue;
                }

                affected= aff_status_entry(Entry->Pheno[trp].Affection.Status, Entry->Pheno[trp].Affection.Class,  Locus);
                /* This works because we have already skipped 0 status individuals */

                if (affected != 2) {
                    affected = 1;
                }
                for (ipen=0; ipen < 3; ipen++) {
                    /* write up to 3 entries */
                    if (sex_linked == 0)
                        pen = ((affected == 1) ? 1.0 - lacp->AutoPen[ipen] : lacp->AutoPen[ipen]);
                    else if (sex_linked && sex == 1) {
                        pen = ((affected == 1) ? 1.0 - lacp->MalePen[ipen] : lacp->MalePen[ipen]);
                        if (ipen) ipen++;
                    } else if (sex_linked && sex == 2)
                        pen = ((affected == 1) ? 1.0 - lacp->FemalePen[ipen] : lacp->FemalePen[ipen]);
                    else
                        pen = 0;

                    if (pen > 0.0) {
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
                            fprintf(filep, "%s,", Entry->UniqueID);
                        } else if (OrigIds[0] == 1 || OrigIds[0] == 2) {
                            fprintf(filep, "%s,", Entry->OrigID);
                        } else {
                            fprintf(filep, "%d,", Entry->ID);
                        }

                        fprintf(filep, "%s,%s,%8f,%d\n",
                                strtail(Locus->Name, MENDEL7_MAX_LOCUS_NAME_LEN), genostr[ipen], pen, affected);
                    }
                }
            }
        }
        /* in loop-over trait mode, close the current file */
        if (LoopOverTrait == 1)  {
            fclose(filep);
        }
    }

    if (LoopOverTrait == 0) {
        fclose(filep);
    }

    sprintf(err_msg, "      Penetrance file:        %s",
            fl_name);
    mssgf(err_msg);
    return;

}

#if 0
static void csv_write_mendel_map_using_sex_specific(char *mapfile, linkage_locus_top *LTop,
						    analysis_type analysis, int xlinked)
{
    int i, j, tr, nloop=0, *trp, num_affec=num_traits;
    int num_markers, *markers, num_numbered=0;
    char mapfl_name[2*FILENAME_LENGTH];
    double difff = INVALID_POS_DIFF, diffm = INVALID_POS_DIFF;

    FILE *filep;

    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY) {
            num_numbered++;
        }
    }

    if (num_numbered == 0) {
        num_markers = NumChrLoci;
    } else {
        num_markers = num_numbered;
    }

    markers = CALLOC((size_t) num_markers, int);

    j=0;
    if (num_numbered > 0) {
        for (i=0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[i]].Type == BINARY) {
                markers[j] = ChrLoci[i]; j++;
            }
        }
    } else if (LoopOverTrait == 0) {
        for (i=0; i < NumChrLoci; i++) {
            markers[i] = ChrLoci[i];
        }
    }

    NLOOP;

    if (analysis != TO_HWETEST) {
        if (num_affec > 0) trp = &(global_trait_entries[0]);
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(mapfl_name, "%s/%s", output_paths[tr], mapfile);

        if (num_markers <= 0) {
            delete_file(mapfl_name); continue;
        }

        if ((filep = fopen(mapfl_name, "a")) == NULL) {
            errorvf("Could not open map file %s.\n", mapfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_numbered == 0) {
            if (LoopOverTrait == 1 && LTop->Locus[*trp].Type == AFFECTION) {
	      fprintf(filep, "%s\n", strtail(LTop->Locus[*trp].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
                if (xlinked) {
                    fprintf(filep, ",\n");
                } else {
                    fprintf(filep, "\n");
                }
            }
        }

        for(i = 1; i < num_markers; i++) {
            switch(LTop->Locus[markers[i]].Type) {
            case NUMBERED:
            case BINARY:
                switch(LTop->Locus[markers[i-1]].Type) {
                case NUMBERED:
                case BINARY:
                    if (LTop->Marker[markers[i]].chromosome == LTop->Marker[markers[i-1]].chromosome) {
		               diffm = LTop->Marker[markers[i]].pos_male - LTop->Marker[markers[i-1]].pos_male;
		               difff = LTop->Marker[markers[i]].pos_female - LTop->Marker[markers[i-1]].pos_female;
                    } else {
                        diffm = difff = INVALID_POS_DIFF;
                    }
                    break;

                case AFFECTION :
                    if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                        difff = diffm = INVALID_POS_DIFF;
                    }
                    break;
                default:
                    break;
                }
                break;
            case AFFECTION :
                if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                    difff = diffm = INVALID_POS_DIFF;
                }
                break;
            default:
                break;
            }

            if (diffm <= INVALID_POS_DIFF) {
	      fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
                if (xlinked) {
                    fprintf(filep, ",\n");
                } else {
                    fprintf(filep, "\n");
                }
            } else {
                if (diffm < 0.0) {
                    sprintf(err_msg,
                            "%s (%.4g) and %s (%.4g) are not ordered by increasing male map distance!",
                            LTop->Marker[markers[i-1]].Name, LTop->Locus[markers[i-1]].pos_male,
                            LTop->Marker[markers[i]].Name, LTop->Locus[markers[i]].pos_male);
                    warnf(err_msg);
                    warnf("Setting the distance between these to 0.001 cM.");
                    diffm=0.001;
                }
                if (difff < 0.0 && xlinked) {
		    sprintf(err_msg,
			    "%s (%.4g) and %s (%.4g) are not ordered by increasing female map distance!",
			    LTop->Marker[markers[i-1]].Name, LTop->Locus[markers[i-1]].pos_female,
			    LTop->Marker[markers[i]].Name, LTop->Locus[markers[i]].pos_female);
		    warnf(err_msg);
		    warnf("Setting the distance between these to 0.001 cM.");
                    difff=0.001;
                }

                if (LTop->map_distance_type == 'k') {
                    diffm=kosambi_theta(diffm);
                    if (xlinked) {
                        difff=kosambi_theta(difff);
                    }
                } else {
                    diffm=haldane_theta(diffm);
                    if (xlinked) {
                        difff=haldane_theta(difff);
                    }
                }
                fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
                if (xlinked) {
		    fprintf(filep, "%f,%f ! %s\n", difff, diffm,
			    (LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
                } else {
                    fprintf(filep, "%f ! %s\n", diffm,
			    (LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
                }
            }
        }

        fprintf(filep, "%s\n", strtail(LTop->Locus[markers[num_markers-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
        fclose(filep);
        if (nloop == 1) break;
    }

    sprintf(err_msg, "        Map file:             %s", mapfile);
    mssgf(err_msg);
    free(markers);
}
#endif
static void csv_write_mendel_map_using_sex_specific(char *mapfile, linkage_locus_top *LTop,
						    analysis_type analysis, int xlinked)
{
    int i, j, tr, nloop=0, *trp, num_affec=num_traits;
    int num_markers, *markers, num_numbered=0;
    char mapfl_name[2*FILENAME_LENGTH];
    double difff = INVALID_POS_DIFF, diffm = INVALID_POS_DIFF;

    FILE *filep;

    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY) {
            num_numbered++;
        }
    }

    if (num_numbered == 0) {
        num_markers = NumChrLoci;
    } else {
        num_markers = num_numbered;
    }

    markers = CALLOC((size_t) num_markers, int);

    j=0;
    if (num_numbered > 0) {
        for (i=0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[i]].Type == BINARY) {
                markers[j] = ChrLoci[i]; j++;
            }
        }
    } else if (LoopOverTrait == 0) {
        for (i=0; i < NumChrLoci; i++) {
            markers[i] = ChrLoci[i];
        }
    }

    NLOOP;

    if (analysis != TO_HWETEST) {
        if (num_affec > 0) trp = &(global_trait_entries[0]);
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(mapfl_name, "%s/%s", output_paths[tr], mapfile);

        if (num_markers <= 0) {
            delete_file(mapfl_name); continue;
        }

        if ((filep = fopen(mapfl_name, "a")) == NULL) {
            errorvf("Could not open map file %s.\n", mapfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_numbered == 0) {
            if (LoopOverTrait == 1 && LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "%s\n", strtail(LTop->Locus[*trp].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
		fprintf(filep, ",\n");
            }
        }

        for(i = 1; i < num_markers; i++) {
            switch(LTop->Locus[markers[i]].Type) {
            case NUMBERED:
            case BINARY:
            case XLINKED:
            case YLINKED:
                switch(LTop->Locus[markers[i-1]].Type) {
                case NUMBERED:
                case BINARY:
                case XLINKED:
                case YLINKED:
                    if (LTop->Marker[markers[i]].chromosome == LTop->Marker[markers[i-1]].chromosome) {
		               diffm = LTop->Marker[markers[i]].pos_male - LTop->Marker[markers[i-1]].pos_male;
		               difff = LTop->Marker[markers[i]].pos_female - LTop->Marker[markers[i-1]].pos_female;
                    } else {
                        diffm = difff = INVALID_POS_DIFF;
                    }
                    break;

                case AFFECTION :
                    if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                        difff = diffm = INVALID_POS_DIFF;
                    }
                    break;
                default:
                    break;
                }
                break;
            case AFFECTION :
                if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                    difff = diffm = INVALID_POS_DIFF;
                }
                break;
            default:
                break;
            }

            if (diffm <= INVALID_POS_DIFF) {
                fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
		fprintf(filep, ",\n");
            } else {
                if (diffm < 0.0) {
                    warnvf("%s (%.4g) and %s (%.4g) are not ordered by increasing male map distance!\n",
                            LTop->Marker[markers[i-1]].Name, LTop->Marker[markers[i-1]].pos_male,
                            LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_male);
                    warnf("Setting the distance between these to 0.001 cM.");
                    diffm=0.001;
                }
                if (difff < 0.0) {
		    warnvf("%s (%.4g) and %s (%.4g) are not ordered by increasing female map distance!\n",
			    LTop->Marker[markers[i-1]].Name, LTop->Marker[markers[i-1]].pos_female,
			    LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_female);
		    warnf("Setting the distance between these to 0.001 cM.");
                    difff=0.001;
                }

                if (LTop->map_distance_type == 'k') {
                    diffm=kosambi_theta(diffm);
		    difff=kosambi_theta(difff);
                } else {
                    diffm=haldane_theta(diffm);
		    difff=haldane_theta(difff);
                }
                fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
		fprintf(filep, "%f,%f ! %s\n", difff, diffm,
			(LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
            }
        }

        fprintf(filep, "%s\n", strtail(LTop->Locus[markers[num_markers-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
        fclose(filep);
        if (nloop == 1) break;
    }

    sprintf(err_msg, "        Map file:             %s", mapfile);
    mssgf(err_msg);
    free(markers);
}

static void csv_write_mendel_map_using_sex_averaged(char *mapfile, linkage_locus_top *LTop,
						    analysis_type analysis)
{
  int i, j, tr, nloop=0, *trp, num_affec=num_traits, display_warn_x=0;
    int num_markers, *markers, num_numbered=0;
    char mapfl_name[2*FILENAME_LENGTH];
    double diff = INVALID_POS_DIFF;

    FILE *filep;

    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY) {
            num_numbered++;
        }
    }

    num_markers = num_numbered;
    if (num_numbered == 0) {
        num_markers = NumChrLoci;
    }

    markers = CALLOC((size_t) num_markers, int);

    if (num_numbered > 0) {
        for (i=0,j=0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[i]].Type == BINARY) {
                markers[j++] = ChrLoci[i];
            }
        }
    } else if (LoopOverTrait == 0) {
        for (i=0; i < NumChrLoci; i++) {
            markers[i] = ChrLoci[i];
        }
    }

    NLOOP;

    if (analysis != TO_HWETEST) {
        if (num_affec > 0) trp = &(global_trait_entries[0]);
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;

        sprintf(mapfl_name, "%s/%s", output_paths[tr], mapfile);

        if (num_markers <= 0) {
            delete_file(mapfl_name); continue;
        }

        if ((filep = fopen(mapfl_name, "a")) == NULL) {
            errorvf("Could not open map file %s.\n", mapfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_numbered == 0) {
            if (LoopOverTrait == 1 && LTop->Locus[*trp].Type == AFFECTION) {
	      fprintf(filep, "%s\n", strtail(LTop->Locus[*trp].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
		fprintf(filep, "\n");
            }
        }

        for(i = 1; i < num_markers; i++) {
            switch(LTop->Locus[markers[i]].Type) {
                case NUMBERED:
                case BINARY:
                case XLINKED:
                case YLINKED:
                    switch(LTop->Locus[markers[i-1]].Type) {
                        case YLINKED:
                            errorvf("Since sex-averaged genetic map was specified the map cannot be used for markers on the Y chromosomes.\n");
                            EXIT(INPUT_DATA_ERROR);
                            break;
                        case XLINKED:
			    if (display_warn_x == 0) {
                                 warnvf("Since sex-averaged genetic map was specified it will be used for markers on the X chromosomes.\n");
                                 display_warn_x++;
                            }
                            // DONT break...
                        case NUMBERED:
                        case BINARY:
                            if (LTop->Marker[markers[i]].chromosome == LTop->Marker[markers[i-1]].chromosome) {
                                diff = LTop->Marker[markers[i]].pos_avg - LTop->Marker[markers[i-1]].pos_avg;
                            } else {
                                diff = INVALID_POS_DIFF;
                            }
                            break;
                        case AFFECTION :
                            if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                                diff = INVALID_POS_DIFF;
                            }
                            break;
                        default:
                            break;
                    }
                    break;
                case AFFECTION :
                    if (LoopOverTrait == 0 && analysis != TO_HWETEST && num_numbered == 0) {
                        diff = INVALID_POS_DIFF;
                    }
                    break;
                default:
                    break;
            }

            if (diff <= INVALID_POS_DIFF) {
                fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
		fprintf(filep, "\n");
            } else {
                if (diff < 0.0) {
                  warnvf("%s (%.4g) and %s (%.4g) are not ordered by increasing map distance!\n",
                          LTop->Marker[markers[i-1]].Name, LTop->Marker[markers[i-1]].pos_avg,
                          LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_avg);
                  warnf("Setting the distance between these to 0.001 cM.");
                  diff = 0.001;
                }

                if (LTop->map_distance_type == 'k') {
                  diff=kosambi_theta(diff);
                } else {
                  diff=haldane_theta(diff);
                }
                fprintf(filep, "%s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
                fprintf(filep, "%f ! %s\n", diff,
                        (LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
            }
        }

        fprintf(filep, "%s\n", strtail(LTop->Locus[markers[num_markers-1]].Name, MENDEL7_MAX_LOCUS_NAME_LEN));
        fclose(filep);
        if (nloop == 1) break;
    }

    sprintf(err_msg, "        Map file:             %s", mapfile);
    mssgf(err_msg);
    free(markers);
}


void csv_write_mendel_map(char *mapfile, linkage_locus_top *LTop,
                          analysis_type analysis, int xlinked)
{
    if (genetic_distance_index > -1) {
        if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            csv_write_mendel_map_using_sex_averaged(mapfile, LTop, analysis);
        } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
            csv_write_mendel_map_using_sex_specific(mapfile, LTop, analysis, xlinked);
        } else {
            errorvf("Either a sex averaged or sex specific map must be provided to write the map file.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    } else {
        // write an error message if no genetic map was specified (currently one must be specified).
    }
}


/*-----------------end of cvs map file---------------------*/
