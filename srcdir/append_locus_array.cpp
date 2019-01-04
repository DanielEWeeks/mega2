/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

/* moved this routine out of pedtree.c */

/* support routines for append_locus_array to clean up the code */
/* This is a temporary change */

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "linkage_ext.h"
/*
     error_messages_ext.h:  my_calloc
            linkage_ext.h:  copy_linkage_locus_rec copy_lpedtree1
*/


static void copy_pp_static_data(const marriage_graph_type *From,
				marriage_graph_type *To)
{
    To->ped=From->ped;
    To->max_id=From->max_id;
    To->num_persons=From->num_persons;
    To->num_marriages=0;
    To->total_data_cols=0;
    To->num_selected_loci=0;
    /* This is before breaking loops, so marriages is NULL */
    To->marriages=NULL;
    strcpy(To->Name, From->Name);
}

/* entry data that is not specific to loci */

void copy_le_static_data(const linkage_ped_rec *From, linkage_ped_rec *To)
{
    int i;
    strcpy(To->OrigID, From->OrigID);
    strcpy(To->PerPre, From->PerPre);
    strcpy(To->UniqueID, From->UniqueID);
    /* FamName is not copied over because it is not needed anymore */
    strcpy(To->FamName, From->FamName);
    To->ID = From->ID;
    To->Father = From->Father;
    To->Mother = From->Mother;
    To->First_Offspring = From->First_Offspring;
    To->Next_PA_Sib = From->Next_PA_Sib;
    To->Next_MA_Sib = From->Next_MA_Sib;
    To->Sex = From->Sex;
    To->OrigProband = From->OrigProband;
    To->Orig_status = From->Orig_status;
    To->genocnt     = From->genocnt;
    To->Ngeno = From->Ngeno;
    /* Possible error: Why was this cast to int and then to void?
       To->TmpData = (void *)((int) From->TmpData); */
    To->TmpData = (void *)(From->TmpData);
    To->IsTyped = From->IsTyped;
    i=0;
    if (From->loopbreakers != NULL) {
        To->loopbreakers = CALLOC((size_t) MAXLOOPMEMBERS+1, int);
        while(From->loopbreakers[i] != -99) {
            To->loopbreakers[i] = From->loopbreakers[i];
            i++;
        }
        To->loopbreakers[i] = -99;
    } else {
        To->loopbreakers = NULL;
    }
    To->ext_ped_num = From->ext_ped_num;
    To->ext_per_num = From->ext_per_num;
}

static void copy_pe_static_data(const person_node_type *From, person_node_type *To)
{
    int k;

    To->father = From->father;
    To->mother = From->mother;
    To->gender = From->gender;
    To->indiv  = From->indiv;
    strcpy(To->uniqueid, From->uniqueid);
    To->node_id = From->node_id;
    To->proband = From->proband;
    To->num_to_marriages=0;
    for (k=0; k<MAXMARRIAGES; k++) To->to_marriage_node_id[k] = -1;
    To->from_marriage_node_id=-1;
    To->genocnt     = From->genocnt;
}

/* affection locus data */

static int copy_trait_loc_data(const pheno_pedrec_data *From,
                               pheno_pedrec_data *To,
                               const linkage_locus_type Type)
{
    int copied;

    switch (Type) {
    case NUMBERED:
    case BINARY:
        copied=0;
        break;
    case AFFECTION:
        To->Affection.Status = From->Affection.Status;
        To->Affection.Class = From->Affection.Class;
        copied=1;
        break;
    case QUANT:
        To->Quant = From->Quant;
        copied=1;
        break;
    default:
        copied=0;
        break;
    }
    return copied;
}

/* marker locus data */

static int copy_marker_loc_data(void *From, int from,
                                void *To,   int to,
                                const linkage_locus_type Type)
{
    int copied;

    switch (Type) {
    case NUMBERED:
    case BINARY:
        copy_2alleles(To, From, to, from);
        copied=1;
        break;
    case AFFECTION:
    case QUANT:
        copied=0;
        break;
    default:
        copied=0;
        break;
    }
    return copied;
}

/* function to copy over the static locustop stuff */

static void copy_static_locus_info(const linkage_locus_top *From, linkage_locus_top *To)
{
    int i;

    To->NumPedigreeCols = From->NumPedigreeCols;
    To->PedRecDataType  = From->PedRecDataType;
    To->RiskLocus       = From->RiskLocus;
    To->RiskAllele      = From->RiskAllele;
    To->SexLinked       = From->SexLinked;
    To->Program         = From->Program;
    To->MutLocus        = From->MutLocus;
    To->Haplotype       = From->Haplotype;
    To->MutMale         = From->MutMale;
    To->MutFemale       = From->MutFemale;
    To->SexDiff         = From->SexDiff;
    To->Interference    = From->Interference;
    To->FlankRecomb     = From->FlankRecomb;
    To->map_distance_type = From->map_distance_type;
    for (i = 0; i < From->LocusCnt - 1; i++)
        To->MaleRecomb[i] = From->MaleRecomb[i];
}

/* ==========================================================================*/
/* append_locus_array:  Puts the locus array of *From onto the beginning of
 *  *To.
 *  Pre-conditions:  *To has its locus count set to the final count, and
 enough memory allocated for that purpose.
 start contains the locus number at which to begin
 the insertion.
 *  Post-conditions:  The loci in *From have been placed into *To with the
 *  corresponding Alleles in the Peds having been copied as well.
 *                    *start is incremented by the number of loci
 *                    inserted.
 */
/*  CAUTION:
    Assumes that on the first call to this function, start is zero.
    In this case, this function copies over non-marker  specific
    pedigree and locus information.
    Assumes on subsequent calls, start is not zero.
    NM -  this has been changed, we are specifically passing in the position
    of insertion now.
*/
void append_locus_array(const linkage_ped_top *From,
                        linkage_ped_top *To,
                        int *start)
{
    int   i, j, k;
    int current, saved;
    int num_entries;
    pheno_pedrec_data *pheno1, *pheno2;
    void *marker1, *marker2;

    /* Using the variable start to keep track of where to start appending */
    current=0;
    if (*start == 0) {
        /* first call to append */
        /* Step 1. Copy over pedigree information */
        To->PedCnt = From->PedCnt;
        To->pedfile_type=From->pedfile_type;
        To->OrigIds = From->OrigIds;
        To->UniqueIds = From->UniqueIds;
        /* Code modified from copy routines of linkage.c */

        /* Only copy over the pedigree information once, on the first
           call to append */
        for (i = 0; i < From->PedCnt; i++)  {
            if (To->pedfile_type == POSTMAKEPED_PFT) {
                copy_lpedtree1(&(From->Ped[i]), &(To->Ped[i]));
                /* entry info that doesn't change from one chromosome to another */
                for (j = 0; j < From->Ped[i].EntryCnt; j++)  {
                    copy_le_static_data(&(From->Ped[i].Entry[j]), &(To->Ped[i].Entry[j]));
                    /* affection and quant loci */
                    current=0;
                    for (k = 0; k < From->LocusTop->LocusCnt; k++)  {
                        pheno1=&(From->Ped[i].Entry[j].Pheno[k]);
                        pheno2=&(To->Ped[i].Entry[j].Pheno[current]);
                        current +=
                            copy_trait_loc_data(pheno1, pheno2,
                                                From->LocusTop->Locus[k].Type);
                    }
                }
            } else {
                /* pre-makeped format */
                /* 1. Pedigree info */
                copy_pp_static_data(&(From->PTop[i]),&(To->PTop[i]));
                /* entry info that doesn't change from one chromosome to another */
                for (j = 0; j < From->PTop[i].num_persons; j++)  {
                    copy_pe_static_data(&(From->PTop[i].persons[j]),
                                        &(To->PTop[i].persons[j]));

                    /* affection and quant loci */
                    current=0;
                    for (k = 0; k < From->LocusTop->LocusCnt; k++)  {
                        pheno1=&(From->PTop[i].persons[j].pheno[k]);
                        pheno2=&(To->PTop[i].persons[j].pheno[current]);
                        current +=
                            copy_trait_loc_data(pheno1, pheno2,
                                                From->LocusTop->Locus[k].Type);
                    }
                }
            }
        }

        current=0;
        /* Step 2. locus top related stuff that is performed only once*/
        copy_static_locus_info(From->LocusTop, To->LocusTop);
        for (i = 0; i < From->LocusTop->LocusCnt; i++) {
            if (From->LocusTop->Locus[i].Type == AFFECTION ||
               From->LocusTop->Locus[i].Type == QUANT) {
                To->LocusTop->Locus[i].Allele =
                    CALLOC((size_t) From->LocusTop->Locus[i].AlleleCnt, linkage_allele_rec);
                copy_linkage_locus_rec(&(From->LocusTop->Locus[i]),
                                       &(To->LocusTop->Locus[current]));
                current++;
            }
        }
    }

    saved =  current;
    /* the rest is executed everytime this routine is called
       including copying genotypes, and marker locus information
    */

    for (i = 0; i < From->PedCnt; i++) {
        if (To->pedfile_type == POSTMAKEPED_PFT) num_entries=From->Ped[i].EntryCnt;
        else num_entries=From->PTop[i].num_persons;

        for (j = 0; j <num_entries ; j++)	{

            if (*start > 0) current = *start;
            else current = saved;
            for (k = 0; k < From->LocusTop->LocusCnt; k++)  {
                if (To->pedfile_type == POSTMAKEPED_PFT) {
                    marker1=From->Ped[i].Entry[j].Marker;
                    marker2=To->Ped[i].Entry[j].Marker;
                } else {
                    marker1=From->PTop[i].persons[j].marker;
                    marker2=To->PTop[i].persons[j].marker;
                }
                current +=
                    copy_marker_loc_data(marker1, k, marker2, current,
                                         From->LocusTop->Locus[k].Type);
            }
        }
    }

    if (*start > 0) current = *start;
    else current = saved;

    for (i = 0; i < From->LocusTop->LocusCnt; i++) {
        if (From->LocusTop->Locus[i].Type == NUMBERED ||
            From->LocusTop->Locus[i].Type == XLINKED) {
            To->LocusTop->Locus[current].Allele =
                CALLOC((size_t) From->LocusTop->Locus[i].AlleleCnt, linkage_allele_rec);
            copy_linkage_locus_rec(&(From->LocusTop->Locus[i]),
                                   &(To->LocusTop->Locus[current]));
            current++;
        }
    }
    *start = current;
}
