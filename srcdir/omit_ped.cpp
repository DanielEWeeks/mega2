/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

/* function omit_ped.c
   Based on the untyped_ped_opt, omit one or more pedigrees
   Should be invoked after read_omit_file
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  mssgf
         user_input_ext.h:  untyped_ped_messg
              utils_ext.h:  EXIT log_line
*/


int is_typed_at_markers(entry_type Entry,
                        linkage_ped_top *Top,
                        linkage_locus_top *LTop,
			double *percent_typed);

int is_phenotyped(entry_type Entry, linkage_ped_top *Top, linkage_locus_top *LTop);

/* static void delete_pedigree(linkage_ped_top *Top, int untyped); */

static void output_untyped_pedlist(int *ungeno, linkage_ped_top *Top, int opt)

{
    int i, untyped=0;
    for (i=0; i < Top->PedCnt; i++) {
        if (ungeno[i] == 1) {
            untyped++;
        }

    }

    if (untyped == Top->PedCnt) {
        errorvf("No typed pedigrees to save, terminating Mega2!\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }

    if (untyped == 0) {
        mssgf("Including all pedigrees (all are typed).");
        return;
    }
    NumTypedPeds = Top->PedCnt - untyped;

    log_line(mssgf);
    mssgf("HALF-TYPED individuals NOT INCLUDED in this step.");
    mssgf("Excluding these pedigrees from analysis because");
    mssgf(err_msg);
    switch(opt) {
    case 1:
        mssgf("they are completely un-genotyped.");
        break;
    case 3:
        mssgf("they are not fully genotyped at one or more markers.");
        break;
    default:
        sprintf(err_msg,
                "they have less than %d individuals who are fully genotyped at any marker.\n",
                opt - 4 + 1);
        mssgf(err_msg);
        break;
    }

    strcpy(err_msg, "");
    if (Top->pedfile_type == POSTMAKEPED_PFT) {
        for (i = 0; i < Top->PedCnt; i++) {
            if (ungeno[i] == 1) {
                if (strlen(err_msg) <= 80) {
                    strcat(err_msg, Top->Ped[i].Name); strcat(err_msg, " ");
                } else {
                    mssgf(err_msg);
                    strcpy(err_msg, Top->Ped[i].Name); strcat(err_msg, " ");
                }
            }
        }
    } else {
        for (i = 0; i < Top->PedCnt; i++) {
            if (ungeno[i] == 1) {
                if (strlen(err_msg) <= 80) {
                    strcat(err_msg, Top->PTop[i].Name); strcat(err_msg, " ");
                } else {
                    mssgf(err_msg);
                    strcpy(err_msg, Top->PTop[i].Name); strcat(err_msg, " ");
                }
            }
        }
    }
    if (strlen(err_msg) > 0) {
        mssgf(err_msg);
    }

    log_line(mssgf);
    return;
}

//
// Omit pedigrees under certain circumstances.
//
// See user_input.c:untyped_ped_messg() for a description of the 'untyped_ped_opt' values.
void omit_peds(int untyped_ped_opt, linkage_ped_top *LPedTreeTop)
{
    int untyped, num_untyped, *untyped_peds;
    int i, j, IsTyped, entry_count;
    char messg[100];
    entry_type Entry;
    double percent_typed;

    if (Mega2Status >= INSIDE_ANALYSIS && !HasMarkers) {
        return;
    }

    untyped_peds = &(UntypedPeds[0]);
    for (i = 0; i < LPedTreeTop->PedCnt; i++) {
        untyped_peds[i] = 0;
    }

    if (Mega2Status != INSIDE_ANALYSIS) {
        sprintf(err_msg, "Pedigree exclusion option : %s.",
                untyped_ped_messg(untyped_ped_opt, messg));
        mssgf(err_msg);
    }

    if (untyped_ped_opt == 2) {
        // Include all pedigrees whether typed or not...
        NumTypedPeds = LPedTreeTop->PedCnt;
        return;
    }

    untyped = 0;

    for (i = 0; i < LPedTreeTop->PedCnt; i++) {
        /* find out the number of typed individuals and set the value for
           IsTyped
        */

        IsTyped=0;
        if (LPedTreeTop->pedfile_type == PREMAKEPED_PFT) {
            entry_count = LPedTreeTop->PTop[i].num_persons;
        } else {
            entry_count = LPedTreeTop->Ped[i].EntryCnt;
        }

        Entry.LEntry = NULL; // silly compiler
        for (j=0; j < entry_count; j++) {
            if (LPedTreeTop->pedfile_type == PREMAKEPED_PFT) {
                Entry.PEntry = &(LPedTreeTop->PTop[i].persons[j]);
            } else if (LPedTreeTop->pedfile_type == POSTMAKEPED_PFT) {
                Entry.LEntry = &(LPedTreeTop->Ped[i].Entry[j]);
            }

            if (is_typed_at_markers(Entry, LPedTreeTop,
                                    LPedTreeTop->LocusTop, &percent_typed) == 1) {
                IsTyped ++;
            }

        }

        if (LPedTreeTop->pedfile_type == POSTMAKEPED_PFT) {
            Entry.LEntry->IsTyped = ((percent_typed > 0)? 1: 0);
            Entry.LEntry->PercentTyped = percent_typed;
        }

        switch(untyped_ped_opt) {
        case 1:
            // Omit any completely untyped pedigrees...
            if (IsTyped == 0) {
                untyped_peds[i]=1;
                untyped++;
            }
            break;
        case 3:
            // Exclude pedigrees not fully typed at one or more markers...
            if (IsTyped < entry_count) {
                untyped_peds[i] = 1;
                untyped++;
            }
            break;

        case 4:
            // Exclude any pedigree with 0 marker-typed people...
        default:
            // Exclude any pedigree with %d or less marker-typed people...
            if (untyped_ped_opt >= 4) {
                num_untyped= untyped_ped_opt - 4;
                if (IsTyped <= num_untyped) {
                    untyped_peds[i] = 1;
                    untyped++;
                }
            } else {
                // untyped_ped_opt =< 0
                printf("Unknown untyping option %d\n", untyped_ped_opt);
                return;
            }
            break;
        }
    }

    output_untyped_pedlist(untyped_peds, LPedTreeTop, untyped_ped_opt);
}



/* had to rewrite this routine because,
   the existing is_typed routine does not work on
   the linkage data structure.
   An individual is typed if genotypes at any marker
   locus is known
*/

int is_typed_at_markers(entry_type Entry, linkage_ped_top *Top,
			linkage_locus_top *LTop, double *percent_typed)

{
    int num_marker_loci=0, num_markers, j;
    int i, is_typed=0;
    void *marker;
    const char *ar1, *ar2;
    int a1, a2;

    if (Top->pedfile_type == POSTMAKEPED_PFT) {
        marker = Entry.LEntry->Marker;
    } else {
        marker = Entry.PEntry->marker;
    }

    if (ChrLoci == NULL) {
        num_markers = LTop->LocusCnt;
    } else {
        num_markers = NumChrLoci;
    }
    for (j=0; j< num_markers; j++) {
        if (ChrLoci != NULL) {
            i = ChrLoci[j];
        } else {
            i = j;
        }
        if (LTop->Locus[i].Type == NUMBERED ||
            LTop->Locus[i].Type == BINARY) {
            if (LTop->PedRecDataType == Premakeped ||
               LTop->PedRecDataType == Postmakeped) {
                num_marker_loci++;
                get_2alleles(marker, i, &a1, &a2);
                is_typed += (a1 > 0 && a2 > 0) ?  1 : 0;
            } else {
                get_2Ralleles(marker, i, &ar1, &ar2);
                is_typed +=
                    ((strcmp(ar1, REC_UNKNOWN) &&
                      strcmp(ar2, REC_UNKNOWN))?
                     1 : 0);
            }
        }
    }

    if (percent_typed != NULL) {
        *percent_typed = (double)is_typed/(double)num_marker_loci;
    }

    if (is_typed > 0) {
        return 1;
    } else {
        return 0;
    }

}

int is_phenotyped(entry_type Entry, linkage_ped_top *Top, linkage_locus_top *LTop)

{


    int i1, i, is_typed=0;

    for (i1=0; i1 < num_traits; i1++) {
        i=global_trait_entries[i1];
        if (i  == -1) continue;

        if (LTop->Locus[i].Type == QUANT) {
            if (Top->pedfile_type == POSTMAKEPED_PFT) {
                /* linkage format  */
                is_typed += ((fabs(Entry.LEntry->Pheno[i].Quant - MissingQuant) > EPSILON)?
                             1 : 0);
            } else {
                /* pre-makeped format */
                is_typed += ((fabs(Entry.PEntry->pheno[i].Quant - MissingQuant) > EPSILON)?
                             1 : 0);
            }
        } else  if (LTop->Locus[i].Type == AFFECTION) {
            if (Top->pedfile_type == POSTMAKEPED_PFT) {
                is_typed += ((Entry.LEntry->Pheno[i].Affection.Status)? 1 : 0);
            } else {
                is_typed += ((Entry.PEntry->pheno[i].Affection.Status)? 1 : 0);
            }
        }
    }

    if (is_typed > 0) {
        return 1;
    } else {
        return 0;
    }

}
