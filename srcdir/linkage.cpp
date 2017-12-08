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
#include "analysis.h"

#include "linkage.h"

#include "error_messages_ext.h"
#include "list_ext.h"
#include "makeped_ext.h"
#include "pedtree_ext.h"
#include "utils_ext.h"
#include "user_input_ext.h"
#include "append_locus_array_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc my_malloc my_realloc warnf
               list_ext.h:  append_to_list_tail free_all_from_list list_free_all list_iterate new_list pop_first_list_entry
            makeped_ext.h:  copy_pedrec_data copy_preped_peds free_marriage_graph
            pedtree_ext.h:  copy_le_static_data
              utils_ext.h:  EXIT
*/

// http://www.jurgott.org/linkage/LinkagePC.html
// http://www.jurgott.org/linkage/liped.html

void copy_linkage_locus_top1(linkage_locus_top *from, linkage_locus_top *to,
			     int locuscnt, int *locus_inds);
void copy_linkage_ped_top(linkage_ped_top * Original, linkage_ped_top * Copy, int malloc_it);
/* void free_linkage_ped_top(linkage_ped_top * Copy); */
void malloc_linkage_ped_top(linkage_ped_top * Original, linkage_ped_top * Copy,
			    int alloc_alleles);
void delink_linkage_locus_top(linkage_locus_top *Copy);
void malloc_linkage_locus_top(linkage_ped_top *Original,
			      linkage_ped_top *Copy,
			      int alloc_alleles);
void copy_lpedtree1(linkage_ped_tree *From, linkage_ped_tree *To);
void copy_lpedrec(linkage_ped_rec *From, linkage_ped_rec *To,
		  linkage_ped_top *top);
void copy_lpedrec_raw_alleles(linkage_ped_rec *From, linkage_ped_rec *To,
			      linkage_ped_top *top);
void copy_lpedtree(linkage_ped_tree *From, linkage_ped_tree *To);
void clear_lallelerec(linkage_allele_rec *Allele);
void clear_laffclass(linkage_affection_class *Class);
void clear_llocusdata(linkage_locus_data *Data, linkage_locus_type Type);
void clear_llocusrec(linkage_locus_rec *Locus, linkage_locus_type Type);
void clear_llocustop(linkage_locus_top *LTop);
void clear_lpedrec(linkage_ped_rec *Entry);
void clear_lpedtree( linkage_ped_tree *Ped);
void clear_lpedtop(linkage_ped_top *PTop);
void clear_llooprec(linkage_loop_rec *Loop);
linkage_locus_top *new_llocustop(void);
linkage_ped_rec *new_lpedrec(void);
linkage_ped_tree *new_lpedtree(void);
linkage_ped_top *new_lpedtop(void);
linkage_loop_rec *new_llooprec(void);
void free_all_from_llocusrec(linkage_locus_rec *Locus);
void free_all_from_llocustop(linkage_locus_top *LTop);
void free_all_from_lpedrec(linkage_ped_rec *Entry, int offset);
void free_pedrec_raw_alleles(void *Entry, linkage_locus_top *LTop);
void free_all_from_lpedtree(linkage_ped_tree *Ped, int offset);
void free_all_from_lpedtop(linkage_ped_top *PTop);
void free_all_including_lpedtop(linkage_ped_top **PTop);
void free_all_including_llocustop(linkage_locus_top *LTop);
void copy_linkage_locus_top(linkage_locus_top *from, linkage_locus_top *to);
linkage_loop_rec *get_llooprec(listhandle **List, int Num);
void count_lgenotypes(linkage_ped_top *Top, size_t *num_inds, size_t *males,
		      size_t *females, size_t *num_untyped, size_t *num_typed,
		      size_t *peds_typed, size_t *males_typed, size_t *females_typed,
		      size_t *half_typed);
int is_typed_lentry(linkage_ped_rec Entry, linkage_locus_top *LTop,
		    int mode, double *percent_typed);

void clean_reordered_markers(linkage_locus_top *LTop, analysis_type analysis);

#ifdef DEBUG_MEGA2
static linkage_allele_rec *new_lallelerec();
static linkage_locus_rec *new_llocusrec(linkage_locus_type Type);
#endif /* DEBUG_MEGA2 */

static void    *lcollapse_iterate(list_data EntryData, void *It);
static void    *loop_list_iterate(list_data EntryData, void *Arg);
static void    connect_loopbreakers(linkage_ped_rec *S_Entry,
				    linkage_ped_rec *D_Entry,
				    linkage_ped_tree *Ped,
				    linkage_ped_top *Top1);
static    void  copy_llocusdata(linkage_locus_rec *from, linkage_locus_rec *to,
				linkage_locus_type Type, int num_alleles);

/* This function has been modified (NM- Dec 6,2001)
   It simply  used to assign pointers: now, we actually copy
   over the data. This assumes that the *To structure is
   previously allocated  */

void copy_lpedrec(linkage_ped_rec *From, linkage_ped_rec *To,
		  linkage_ped_top *top)

{

    copy_le_static_data(From, To);
/*
    copy_pedrec_data(From->Data, To->Data, top->LocusTop);
*/
/*  Note From->Data must be set to NULL IFF you don't want both structures to SHARE data */
    To->Pheno  = From->Pheno;
    To->Marker = From->Marker;
}

/* This is a special function to be called when we are
   reading in pedigree data for the purpose of recoding, and
   should be only called along with copy_lpedrec. For now,
   try only assigning the pointers, so have to be careful, that
   data.RAlleles are not freed before recoding is completed.

*/

void   copy_lpedtree(linkage_ped_tree *From, linkage_ped_tree *To)

{
    if ((From != NULL) && (To != NULL))
        {
            *To = *From;
// do this instead if you have to:
//   To->Num = From->Num;
//   To->EntryCnt = From->EntryCnt;
//   To->Entry = From->Entry;
//   To->Proband = From->Proband;
//   To->Loops = From->Loops;
        }
}

void            clear_lallelerec(linkage_allele_rec *Allele)

{
    if (Allele == NULL)
        return;
    Allele->Frequency = UNDEF;
}


void            clear_laffclass(linkage_affection_class *Class)

{
    if (Class == NULL)
        return;
    Class->AutoPen = NULL;
}


void            clear_llocusdata(linkage_locus_rec *Locus, linkage_locus_type Type)

{
    if (Locus == NULL)
        return;
    switch (Type)
        {
        case QUANT:
            Locus->Pheno->Props.Quant.ClassCnt = 0;
            Locus->Pheno->Props.Quant.Mean = NULL;
            Locus->Pheno->Props.Quant.Variance = NULL;
            break;
        case AFFECTION:
            Locus->Pheno->Props.Affection.ClassCnt = 0;
            Locus->Pheno->Props.Affection.PenCnt = 0;
            Locus->Pheno->Props.Affection.Class = NULL;
            break;
        case BINARY:
            Locus->Marker->chromosome=0;
            Locus->Marker->Props.Binary.FactorCnt = 0;
            Locus->Marker->Props.Binary.Factor = NULL;
            break;
        case NUMBERED:
            Locus->Marker->chromosome=0;
            Locus->Marker->Props.Numbered.Recoded=0;
            Locus->Marker->Props.Numbered.NumAlleles=0;
            Locus->Marker->Props.Numbered.SelectOpt=0;
            break;
        default:
            break;
        }
}

/* new function  to copy over locus data - NM 11/29/2001*/

static    void  copy_llocusdata(linkage_locus_rec *fromrec, linkage_locus_rec *torec,
				linkage_locus_type Type, int num_alleles)

{
    int  i, j;
    int js;
    pheno_data *p_from = &(fromrec->Pheno->Props), *p_to = &(torec->Pheno->Props);
    marker_data *m_from = &(fromrec->Marker->Props), *m_to = &(torec->Marker->Props);

   switch (Type) {
    case QUANT:
        p_to->Quant.ClassCnt = p_from->Quant.ClassCnt;
        p_to->Quant.Mean     = p_from->Quant.Mean;
        p_to->Quant.Variance = p_from->Quant.Variance;
        break;
    case AFFECTION:
        /* refer to code for reading in affection status locus
           in read_linkage_locus_file(), read_files.c */
        p_to->Affection.ClassCnt = p_from->Affection.ClassCnt;
        p_to->Affection.PenCnt   = p_from->Affection.PenCnt;
        p_to->Affection.Class    = CALLOC((size_t) p_from->Affection.ClassCnt,
                                     linkage_affection_class);
        for(i=0; i < p_from->Affection.ClassCnt; i++) {
            p_to->Affection.Class[i].AutoPen = CALLOC((size_t) p_from->Affection.PenCnt, double);
            for(j=0; j < p_from->Affection.PenCnt; j++) {
                p_to->Affection.Class[i].AutoPen[j] =
                    p_from->Affection.Class[i].AutoPen[j];
            }
            p_to->Affection.Class[i].FemalePen = CALLOC((size_t) p_from->Affection.PenCnt, double);
//xx        for (js= 0; js < num_alleles; js++)
            for (js= 0; js < p_from->Affection.PenCnt; js++) {
                p_to->Affection.Class[i].FemalePen[js] = p_from->Affection.Class[i].FemalePen[js];
            }
            p_to->Affection.Class[i].MalePen = CALLOC((size_t) num_alleles, double);
            for (js= 0; js < num_alleles; js++) {
                p_to->Affection.Class[i].MalePen[js] = p_from->Affection.Class[i].MalePen[js];
            }
        }
        break;
    case BINARY:
        /* refer to code for reading in binary locus
           in read_linkage_locus_file(), read_files.c */
        m_to->Binary.FactorCnt = m_from->Binary.FactorCnt;
        m_to->Binary.Factor = CALLOC((size_t) m_from->Binary.FactorCnt, char *);
        for(i=0; i < m_from->Binary.FactorCnt; i++) {
            m_to->Binary.Factor[i] = CALLOC((size_t) num_alleles, char);
            for(js=0; js < num_alleles; js++) {
                m_to->Binary.Factor[i][js] = m_from->Binary.Factor[i][js];
            }
        }
        break;
    case NUMBERED:
        m_to->Numbered.Recoded=    m_from->Numbered.Recoded;
        m_to->Numbered.NumAlleles= m_from->Numbered.NumAlleles;
        m_to->Numbered.SelectOpt=  m_from->Numbered.SelectOpt;
        break;
    default:
        break;
    }
}


void clear_llocusrec(linkage_locus_rec *Locus, linkage_locus_type Type)
{
    if (Locus == NULL) return;
    Locus->LocusName = NULL;
    Locus->AlleleCnt = 0;
    Locus->Allele = NULL;
    Locus->Type = TYPE_UNSET;
    Locus->Class = CLASS_UNSET;
    clear_llocusdata(Locus, Type);
}


void            clear_llocustop(linkage_locus_top *LTop)
{
    if (LTop == NULL) return;
    LTop->LocusCnt = 0;
    LTop->PhenoCnt = 0;
    LTop->MarkerCnt = 0;
    LTop->NumPedigreeCols=0;
    LTop->Locus = NULL;
    LTop->Pheno = NULL;
    LTop->Marker = NULL;
    LTop->RiskLocus = UNDEF;
    LTop->RiskAllele = UNDEF;
    LTop->SexLinked = UNDEF;
    LTop->Program = UNDEF;
    LTop->MutLocus = UNDEF;
    LTop->Haplotype = UNDEF;
    LTop->MutMale = UNDEF;
    LTop->MutFemale = UNDEF;
    LTop->map_distance_type = 0;
    LTop->SexDiff = UNDEF;
    LTop->Interference = UNDEF;
    LTop->MaleRecomb = NULL;
    LTop->FemaleRecomb = NULL;
    LTop->FlankRecomb = UNDEF;
    LTop->FMRatio = UNDEF;
    LTop->Run.linkmap.recomb_frac=NULL;
}

void            clear_lpedrec(linkage_ped_rec *Entry)
{
    if (Entry == NULL) return;
    Entry->ID = UNDEF;
    Entry->Mother = UNDEF;
    Entry->Father = UNDEF;
    Entry->First_Offspring = UNDEF;
    Entry->Next_PA_Sib = UNDEF;
    Entry->Next_MA_Sib = UNDEF;
    Entry->Sex = UNDEF;
    Entry->Pheno = NULL;
    Entry->Marker = NULL;
    Entry->TmpData = NULL;
    Entry->Orig_status = 0;
    Entry->Ngeno = 0;
    Entry->IsTyped  = 0;
    strcpy(Entry->OrigID, "");
    strcpy(Entry->PerPre, "");
    strcpy(Entry->FamName, "");
    strcpy(Entry->UniqueID, "");
    Entry->loopbreakers=NULL;
    Entry->ext_ped_num=UNDEF;
    Entry->ext_per_num=UNDEF;
    Entry->OrigProband = 0;
}


void            clear_lpedtree( linkage_ped_tree *Ped)
{
    if (Ped == NULL)
        return;
    Ped->Num = UNDEF;
    Ped->EntryCnt = 0;
    Ped->Entry = NULL;
    Ped->OriginalID = 0;
    Ped->IsTyped  = 0;
    Ped->Proband = UNDEF;
    Ped->Loops = NULL;
}

void            clear_lpedtop(linkage_ped_top *PTop)

{
    if (PTop == NULL)
        return;
    PTop->OrigIds=0;
    PTop->UniqueIds=0;
    PTop->IndivCnt=0;
    PTop->pedfile_type = 0;
    PTop->PedCnt = 0;
    PTop->Ped = NULL;
    PTop->PTop = NULL;
    PTop->LocusTop = NULL;
    PTop->EXLTop = NULL;
}


void            clear_llooprec(linkage_loop_rec *Loop)

{
    int i;
    if (Loop == NULL)
        return;
    Loop->Num = UNDEF;
    for (i=0; i<MAXMARRIAGES; i++) {
        Loop->Entry[i] = UNDEF;
    }
}


#ifdef DEBUG_MEGA2
linkage_allele_rec *new_lallelerec()
{
    linkage_allele_rec *newrec = MALLOC(linkage_allele_rec);
    clear_lallelerec(newrec);
    return newrec;
}


linkage_locus_rec *new_llocusrec(linkage_locus_type Type)

{
    linkage_locus_rec *newrec = MALLOC(linkage_locus_rec);
    clear_llocusrec(newrec, Type);
    return newrec;
}
#endif

linkage_locus_top *new_llocustop(void)
{
    linkage_locus_top *newtop = MALLOC(linkage_locus_top);
    clear_llocustop(newtop);
    return newtop;
}


linkage_ped_rec *new_lpedrec(void)
{
    linkage_ped_rec *newrec = MALLOC(linkage_ped_rec);
    clear_lpedrec(newrec);
    return newrec;
}


linkage_ped_tree *new_lpedtree(void)
{
    linkage_ped_tree *newtree = MALLOC(linkage_ped_tree);
    clear_lpedtree(newtree);
    return newtree;
}


linkage_ped_top *new_lpedtop(void)
{
    linkage_ped_top *newtop = MALLOC(linkage_ped_top);
    clear_lpedtop(newtop);
    return newtop;
}


linkage_loop_rec *new_llooprec(void)
{
    linkage_loop_rec *newloop = MALLOC(linkage_loop_rec);
    clear_llooprec(newloop);
    return newloop;
}


void   free_all_from_llocusrec(linkage_locus_rec *Locus)

{
    int    factor;
    int classIdx;

    if (Locus == NULL)
        return;
    if (Locus->LocusName != NULL) {
        free(Locus->LocusName);
        Locus->LocusName = NULL;
    }
    if (Locus->Allele != NULL)     {
        free(Locus->Allele);
        Locus->Allele = NULL;
    }
    switch (Locus->Type)     {
    case AFFECTION:
        if (Locus->Pheno->Props.Affection.Class != NULL)   {
            for (classIdx = 0; classIdx < Locus->Pheno->Props.Affection.ClassCnt; classIdx++) {
                if (Locus->Pheno->Props.Affection.Class[classIdx].AutoPen != NULL)
                    free(Locus->Pheno->Props.Affection.Class[classIdx].AutoPen);
                if (Locus->Pheno->Props.Affection.Class[classIdx].FemalePen != NULL)
                    free(Locus->Pheno->Props.Affection.Class[classIdx].FemalePen);
                if (Locus->Pheno->Props.Affection.Class[classIdx].MalePen != NULL)
                    free(Locus->Pheno->Props.Affection.Class[classIdx].MalePen);
            }
            free(Locus->Pheno->Props.Affection.Class);
            Locus->Pheno->Props.Affection.Class = NULL;
        }

        break;

    case BINARY:
        if (Locus->Marker->Props.Binary.Factor != NULL)   {
            for (factor = 0; factor < Locus->Marker->Props.Binary.FactorCnt; factor++)    {
                if (Locus->Marker->Props.Binary.Factor[factor] != NULL)
                    free(Locus->Marker->Props.Binary.Factor[factor]);
            }
            free(Locus->Marker->Props.Binary.Factor);
        }
        break;
    case QUANT:
    case NUMBERED:
        /* Nothing to do, since quant has 3 atomic variables inside the
           locus->Data, and Numbered  has nothing */
    default:
        break;
    }
}

void   free_all_from_llocustop(linkage_locus_top *LTop)

{
    int             locus;
    if (LTop == NULL)
        return;
    if (LTop->Locus != NULL) {
        for (locus = 0; locus < LTop->LocusCnt; locus++)
            free_all_from_llocusrec(&(LTop->Locus[locus]));
        free(LTop->Locus);
        LTop->Locus = NULL;
    }
    if (LTop->MaleRecomb != NULL)
        free(LTop->MaleRecomb);
    if (LTop->FemaleRecomb != NULL)
        free(LTop->FemaleRecomb);
    if (LTop->Program == LINKMAP &&
       LTop->Run.linkmap.recomb_frac != NULL) {
        free(LTop->Run.linkmap.recomb_frac);
    }
    return;

}

void free_all_from_lpedrec(linkage_ped_rec *Entry, int offset)

{
    if (Entry == NULL)
        return;
/*
    if (Entry->Data != NULL) {
            free(Entry->Data);
            Entry->Data = NULL;
    }
*/
    if (Entry->Pheno != NULL) {
            free(Entry->Pheno);
            Entry->Pheno = NULL;
    }
    if (Entry->Marker != NULL) {
        marker_free(Entry->Marker, offset);
        Entry->Marker = NULL;
    }
    if (Entry->loopbreakers != NULL)
        free(Entry->loopbreakers);
}

void free_all_from_lpedtree(linkage_ped_tree *Ped, int offset)

{
    int entry;
    if (Ped == NULL)
        return;

    if (Ped->Loops != NULL) {
        free_all_from_list(Ped->Loops, free);
    }

    if (Ped->Entry != NULL) {
            for (entry = 0; entry < Ped->EntryCnt; entry++)
                free_all_from_lpedrec(&(Ped->Entry[entry]), offset);
            free(Ped->Entry);
            Ped->Entry = NULL;
    }
    if (Ped->Loops != NULL) {
            list_free_all(Ped->Loops, NULL);
            Ped->Loops = NULL;
    }
}

void            free_lpedtop_pedinfo(linkage_ped_top *PTop)

{

    int ped;

    if (PTop == NULL) {
        return;
    }

    if (PTop->Ped != NULL) {
        for (ped = 0; ped < (PTop->PedCnt-1); ped++) {
            /*      printf("free %d %s\n", ped+1, PTop->Ped[ped].Name); */
            if (PTop->Ped[ped].Loops != NULL) {
                free_all_from_list(PTop->Ped[ped].Loops, free);
            }
            if (PTop->Ped[ped].Entry != NULL) {
                free(PTop->Ped[ped].Entry);
            }
        }
    }
    return;
}

void free_all_from_lpedtop(linkage_ped_top *PTop)

{
    int ped;
    int offset = PTop->LocusTop->PhenoCnt;

    if (PTop == NULL)
        return;
    if (PTop->pedfile_type == 0) {
        if (PTop->Ped != NULL) {
            for (ped = 0; ped < PTop->PedCnt; ped++)
                free_all_from_lpedtree(&(PTop->Ped[ped]), offset);
            free(PTop->Ped);
            PTop->Ped = NULL;
        }
    } else {
        if (PTop->PTop != NULL) {
            for (ped = 0; ped < PTop->PedCnt; ped++)
                free_marriage_graph(&(PTop->PTop[ped]));
            free(PTop->PTop);
            PTop->PTop=NULL;
        }
    }
    return;
}


void            free_all_including_lpedtop(linkage_ped_top **PTop)

{
    if (*PTop == NULL)
        return;
    free_all_from_lpedtop(*PTop);
    free_all_from_llocustop((*PTop)->LocusTop);
    (*PTop)->LocusTop = NULL;
    free(*PTop);
    *PTop=NULL;
}

/* new function copy_linkage_locus_rec to physically copy
   a locus record to another one which has been
   pre-allocated */
void  copy_linkage_locus_rec(linkage_locus_rec *from,
			     linkage_locus_rec *to)

{

    int j;

    if (from->LocusName != NULL)
        strcpy(to->LocusName, from->LocusName);
    to->AlleleCnt = from->AlleleCnt;
    to->Type = from->Type;
//  to->Class = to->Class; broken since v4.5.4
    to->Class = from->Class;
    /* We need new code here to copy over the .Data structure */
    copy_llocusdata(from, to, from->Type, from->AlleleCnt);
    to->number = from->number;
    switch(to->Type) {
    case BINARY:
    case NUMBERED:
        for (j = 0; j < from->AlleleCnt; j++)  {
            to->Allele[j].Frequency = from->Allele[j].Frequency;
            to->Allele[j].index = from->Allele[j].index;
        }
        to->Marker->pos_avg = from->Marker->pos_avg;
        to->Marker->pos_male = from->Marker->pos_male;
        to->Marker->pos_female = from->Marker->pos_female;
        to->Marker->chromosome = from->Marker->chromosome;
        to->Marker->error_prob = from->Marker->error_prob;
        break;
    case AFFECTION:
    case QUANT:
        for (j = 0; j < from->AlleleCnt; j++)  {
            to->Allele[j].Frequency = from->Allele[j].Frequency;
            to->Allele[j].index = from->Allele[j].index;
        }
        break;
    default:
        break;
    }
}

/* NM- June 27, separated out copying over static info,
   that doesn't depend on the number of loci */

static void copy_ltop_static_info(linkage_locus_top *from, linkage_locus_top *to)

{
    to->NumPedigreeCols = from->NumPedigreeCols;
    to->PedRecDataType = from->PedRecDataType;
    to->RiskLocus = from->RiskLocus;
    to->RiskAllele = from->RiskAllele;
    to->SexLinked = from->SexLinked;
    to->Program = from->Program;
    to->MutLocus = from->MutLocus;
    to->Haplotype = from->Haplotype;
    to->MutMale = from->MutMale;
    to->MutFemale = from->MutFemale;
    to->SexDiff = from->SexDiff;
    to->Interference = from->Interference;
    to->FlankRecomb = from->FlankRecomb;
    to->map_distance_type = from->map_distance_type;
    if (to->Program == MLINK) {
        to->Run.mlink.start_theta = from->Run.mlink.start_theta;
        to->Run.mlink.increment = from->Run.mlink.increment;
        to->Run.mlink.stop_theta = from->Run.mlink.stop_theta;
        to->Run.mlink.num_evals = from->Run.mlink.num_evals;
    }
    if (to->Program == LINKMAP) {
        to->Run.linkmap.trait_marker = from->Run.linkmap.trait_marker;
        to->Run.linkmap.stop_theta = from->Run.linkmap.stop_theta;
        to->Run.linkmap.num_evals =  from->Run.linkmap.num_evals;
    }

    if (from->SexDiff == CONSTANT_SEX_DIFF) {
        to->FMRatio = from->FMRatio;
    }
}

void copy_linkage_locus_top(linkage_locus_top * from, linkage_locus_top * to)
{
    int             i;

    to->LocusCnt  = from->LocusCnt;
    to->PhenoCnt  = from->PhenoCnt;
    to->MarkerCnt = from->MarkerCnt;

    copy_ltop_static_info(from, to);

    for (i = 0; i < from->LocusCnt - 1; i++)
        to->MaleRecomb[i] = from->MaleRecomb[i];
    if (to->Program == LINKMAP && from->Run.linkmap.recomb_frac != NULL) {

        to->Run.linkmap.recomb_frac =  CALLOC((size_t) to->LocusCnt-1, double);
        for (i = 0; i < to->LocusCnt - 1; i++)
            to->Run.linkmap.recomb_frac[i] = from->Run.linkmap.recomb_frac[i];
    } else {
        to->Run.linkmap.recomb_frac=NULL;
    }

    if (from->SexDiff == VARIABLE_SEX_DIFF && from->FemaleRecomb != NULL) {
        for (i = 0; i < from->LocusCnt - 1; i++)
            to->FemaleRecomb[i] = from->FemaleRecomb[i];
    }

    /* Since structures can be assigned under ANSI C, then this
       should work: However, that only works for things that have
       values in the arrays, but does not work as desired for pointers */

    /* instead of assigning pointers (to->Locus[i] = from->Locus[i])
       copy over values, so that when one is freed, the other doesn't
       get corrupted as well.*/

    for (i = 0; i < from->LocusCnt; i++) {
        copy_linkage_locus_rec(&(from->Locus[i]), &(to->Locus[i]));
    }
}

/* same as above, except, a list of loci is also specified */
void copy_linkage_locus_top1(linkage_locus_top *from, linkage_locus_top *to,
			     int locuscnt, int *locus_inds)



{

    int ii, i;

    to->LocusCnt = locuscnt;
    copy_ltop_static_info(from, to);

    if (locuscnt > 0) {
        for (ii = 0; ii < locuscnt - 1; ii++) {
            to->MaleRecomb[ii] = from->MaleRecomb[locus_inds[ii]];
        }
    }
    if (to->Program == LINKMAP && from->Run.linkmap.recomb_frac != NULL && locuscnt > 0) {
        to->Run.linkmap.recomb_frac = CALLOC((size_t) locuscnt-1, double);
        for (ii = 0; ii < to->LocusCnt - 1; ii++) {
            to->Run.linkmap.recomb_frac[ii] = from->Run.linkmap.recomb_frac[locus_inds[ii]];
        }
    } else {
        to->Run.linkmap.recomb_frac=NULL;
    }

    if (from->SexDiff == VARIABLE_SEX_DIFF && locuscnt > 0) {
        for (i = 0; i < locuscnt - 1; i++)
            to->FemaleRecomb[i] = from->FemaleRecomb[locus_inds[i]];
    }

    for (i = 0; i < locuscnt; i++) {
        copy_linkage_locus_rec(&(from->Locus[locus_inds[i]]), &(to->Locus[i]));
    }

    return;
}

/* There was a problem with this code, I assumed that
   get_llooprec(List, Num) was grabbing the <Num>th record,
   if fact, that code had been changed by me some time ago, to
   grab the record with Data->Num = <Num>
*/

void  copy_lpedtree1(linkage_ped_tree *From, linkage_ped_tree *To)

{

    linkage_loop_rec *loop, *loop2;
    list_entry     *list_ptr;
    int j;

    To->Num = From->Num;
    strcpy(To->Name, From->Name);
    strcpy(To->PedPre, From->PedPre);
    To->EntryCnt = From->EntryCnt;
    To->Proband = From->Proband;
    To->OriginalID = From->OriginalID;
    To->IsTyped = From->IsTyped;
    /* loops code has been moved from copy_linkage_ped_top to here */
    if (From->Loops != NULL)  {
        /* get_llooprec(List, Num) gets the loop record with data->Num = Num, or,
           if it does not exist, it creates it */
        /* Loop through all the loop entries in Original */
        list_ptr = From->Loops->Head;
        while (list_ptr != NULL)    {
            loop = (linkage_loop_rec *) list_ptr->Data;
            list_ptr = list_ptr->Next;
            /* Create the new loop record in the Copy */
            loop2 = get_llooprec(&(To->Loops), loop->Num);
            loop2->Num = loop->Num;
            loop2->numEntry = loop->numEntry;
            for (j=0; j < MAXMARRIAGES; j++)
                loop2->Entry[j] = loop->Entry[j];
        }/* Loop on 'loop' */
    } /* if LPed->Loops != Null */
    return;
}

/* All allocation parts have been removed from this */
/*
 * It appears that malloc_linkage_ped_top is ONLY called by
 * copy_linkage_ped_top.  The latter seems to hardly be used
 * anymore, but when used it is after premakeped has been converted
 * to POSTMAKEPED_PFT (i.e. 0) so the pedfile_type check is always TRUE
 */
void  copy_linkage_ped_top(linkage_ped_top *Original, linkage_ped_top *Copy,
			   int malloc_it)
{

    int             i, j;

    if (malloc_it == 1) {
        if (Copy == NULL)  {
            /* This pointer has to be allocated outside this routine
               to make it permanent */
            errorvf("In copy_linkage_ped_top: Null pointer\n");
            EXIT(SYSTEM_ERROR);
        } else {
            malloc_linkage_ped_top(Original, Copy, 1);
        }
    }

    Copy->OrigIds = Original->OrigIds;
    Copy->UniqueIds = Original->UniqueIds;
    Copy->PedCnt = Original->PedCnt;
    Copy->pedfile_type = Original->pedfile_type;

    /* added check for pedfile_type - Nandita (Jan 27, 2000) */

    if (Copy->pedfile_type == 0) {
        for (i = 0; i < Original->PedCnt; i++) {
            /*      clear_lpedtree(&(Copy->Ped[i])); */
            copy_lpedtree1(&(Original->Ped[i]), &(Copy->Ped[i]));
            for (j = 0; j < Original->Ped[i].EntryCnt; j++) {
                copy_lpedrec(&(Original->Ped[i].Entry[j]),
                             &(Copy->Ped[i].Entry[j]), Original);
                /* Pheno/Marker ptrs is SHARED.  since we don't dp:
                   Original...->Pheno  = NULL;
                   Original...->Marker = NULL;
                   SO
                   You MUST "delink_linkage_ped_top(Copy)" before you free the Copy!
                   This function does Copy...->Pheno = NULL; Copy...->Marker = NULL;
                 */
            }
        }
    }

    /* added branch for pe-makeped file type */
    else {
       /*
        for (i = 0; i < Original->PedCnt; i++) {
            / *      clear_prepedtree(&(Copy->PTop[i]));      * /
            for (j = 0; j < Original->PTop[i].num_persons; j++) {
                Copy->PTop[i].persons[j].data =
                    CALLOC((size_t) Original->LocusTop->LocusCnt, linkage_pedrec_data);
            }
        }
       */
        /* added a routine copy_prepedrec() to makeped.c */
        copy_preped_peds(Original, Copy);
    } /* else */

    /*   Copy->LocusTop = new_llocustop(); */
    copy_linkage_locus_top(Original->LocusTop, Copy->LocusTop);

    return;

}

/* All allocation parts have been removed from this */
void  delink_linkage_ped_top(linkage_ped_top *Copy)
{
    int             i, j;

    /* added check for pedfile_type - Nandita (Jan 27, 2000) */

    if (Copy->pedfile_type == 0) {
        for (i = 0; i < Copy->PedCnt; i++) {
            for (j = 0; j < Copy->Ped[i].EntryCnt; j++) {
                Copy->Ped[i].Entry[j].Pheno = NULL;
                Copy->Ped[i].Entry[j].Marker = NULL;
            }
        }
    } else {
        for (i = 0; i < Copy->PedCnt; i++) {
            for (j = 0; j < Copy->PTop[i].num_persons; j++) {
                Copy->PTop[i].persons[j].pheno = NULL;
                Copy->PTop[i].persons[j].marker = NULL;
            }
        }
    }
}

void malloc_linkage_locus_top(linkage_ped_top *Original,
                              linkage_ped_top *Copy,
                              int alloc_alleles)

{
    size_t i;
    size_t             locuscnt;

    if (Copy == NULL) {
        errorvf("NULL pointer Copy!\n");
        EXIT(SYSTEM_ERROR);
    }

    locuscnt = (size_t)Original->LocusTop->LocusCnt;
    if (locuscnt > 1) {
        Copy->LocusTop->MaleRecomb = CALLOC(locuscnt - 1, double);
        for (i = 0; i < locuscnt-1; i++)
            Copy->LocusTop->MaleRecomb[i] = 0.0;
        Copy->LocusTop->FemaleRecomb = CALLOC(locuscnt - 1, double);
        for (i = 0; i < locuscnt-1; i++)
            Copy->LocusTop->FemaleRecomb[i] = 0.0;
    }

    Copy->LocusTop->Locus = CALLOC(locuscnt, linkage_locus_rec);
    Copy->LocusTop->Pheno  = CALLOC((size_t) Original->LocusTop->PhenoCnt, pheno_rec);
    // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
    Copy->LocusTop->Marker = (CALLOC((size_t) Original->LocusTop->MarkerCnt, marker_rec)) - Original->LocusTop->PhenoCnt;

    for (i = 0; i < locuscnt; i++) {
        Copy->LocusTop->Locus[i].LocusName = CALLOC((size_t) FILENAME_LENGTH, char);
        linkage_locus_rec *Locus = &Copy->LocusTop->Locus[i];
        switch(Original->LocusTop->Locus[i].Type) {
        case BINARY:
        case NUMBERED:
            Locus->Marker = &Copy->LocusTop->Marker[i];
            Locus->Pheno  = 0;
            if (alloc_alleles) {
                Locus->AlleleCnt = Original->LocusTop->Locus[i].AlleleCnt;
                Locus->Allele =
                    CALLOC((size_t) Original->LocusTop->Locus[i].AlleleCnt, linkage_allele_rec);
            }
            break;
        case AFFECTION:
        case QUANT:
            Locus->Marker = 0;
            Locus->Pheno  = &Copy->LocusTop->Pheno[i];
            if (alloc_alleles) {
                Locus->AlleleCnt = Original->LocusTop->Locus[i].AlleleCnt;
                Locus->Allele =
                    CALLOC((size_t) Original->LocusTop->Locus[i].AlleleCnt, linkage_allele_rec);
            }
            break;
        default:
            break;
        }
    }

    return;

}

/* This is a copy of the previous function except that
   it specifies which loci to copy
*/

void malloc_linkage_locus_top1(linkage_ped_top *Original,
			       linkage_ped_top *Copy,
			       int alloc_alleles,
			       size_t locuscnt,
			       int *locus_inds)

{
    size_t i;

    if (Copy == NULL) {
        errorvf("NULL pointer Copy!\n");
        EXIT(SYSTEM_ERROR);
    }

    if (locuscnt > 1) {
        Copy->LocusTop->MaleRecomb = CALLOC(locuscnt - 1, double);
        for (i = 0; i < locuscnt-1; i++)
            Copy->LocusTop->MaleRecomb[i] = 0.0;
        Copy->LocusTop->FemaleRecomb = CALLOC(locuscnt - 1, double);
        for (i = 0; i < locuscnt-1; i++)
            Copy->LocusTop->FemaleRecomb[i] = 0.0;
    }

    Copy->LocusTop->Locus = CALLOC(locuscnt, linkage_locus_rec);

    for (i = 0; i < locuscnt; i++) {
        Copy->LocusTop->Locus[i].LocusName = CALLOC((size_t) FILENAME_LENGTH, char);
        if (alloc_alleles) {
            Copy->LocusTop->Locus[i].Allele =
                CALLOC((size_t) Original->LocusTop->Locus[locus_inds[i]].AlleleCnt,
                       linkage_allele_rec);
        }
    }

    return;

}

/* this is a new (not complete fcn() stub) which will malloc MOST of a
   linkage_ped_top to make it look like an existing one
   only space allocation and initialization.
   The malloc_alleles flag is necessary because inside
   append_locus_array, we don't know the ordering of markers
   a priori (not the same as the original) */
void  malloc_linkage_ped_top(linkage_ped_top *Original,
			     linkage_ped_top *Copy,
			     int alloc_alleles)
{

    int i;

    if (Copy == NULL)  {
        errorvf("In malloc_linkage_ped_top: Null pointer?\n");
        EXIT(SYSTEM_ERROR);
    }

    Copy->LocusTop=new_llocustop();
    malloc_linkage_locus_top(Original, Copy, alloc_alleles);

    if (Original->pedfile_type == 0)
        Copy->Ped = CALLOC((size_t) Original->PedCnt, linkage_ped_tree);
    else
        Copy->PTop = CALLOC((size_t) Original->PedCnt, marriage_graph_type);

    for (i = 0; i < Original->PedCnt; i++) {
        if (Original->pedfile_type == 0) {
            Copy->Ped[i].Entry = CALLOC((size_t) Original->Ped[i].EntryCnt, linkage_ped_rec);
/*
    No reason to allocate space since copy_lpedrec() just slams the "from ptr" into the "to".
    If a "move" was desired, the caller then sets the "from ptr" to NULL.  No field by field
    copy is ever done.
            for (j = 0; j < Original->Ped[i].EntryCnt; j++) {
                Copy->Ped[i].Entry[j].Pheno  = CALLOC((size_t) Original->LocusTop->PhenoCnt, pheno_pedrec_data);
                Copy->Ped[i].Entry[j].Marker = marker_alloc((size_t) Original->LocusTop->MarkerCnt, Original->LocusTop->PhenoCnt);
            }
*/
        } else {
            Copy->PTop[i].persons =
                CALLOC((size_t) Original->PTop[i].num_persons, person_node_type);
/*
            for (j = 0; j < Original->PTop[i].num_persons; j++) {
                Copy->PTop[i].persons[j].pheno  = CALLOC((size_t) Original->LocusTop->PhenoCnt, pheno_pedrec_data);
 I don't think this is necessary
                Copy->PTop[i].persons[j].marker = marker_alloc((size_t) Original->LocusTop->MarkerCnt, Original->LocusTop->PhenoCnt);
            }
*/
        }
    }

    return;
}

void  malloc_lpedtree_entries(linkage_ped_tree *Copy, size_t entrycount,
			      size_t locuscnt)

{
    /* This will allocate enough space for a pedtree
       given the number of entries and information
       about locus data */
    Copy->Entry = CALLOC(entrycount, linkage_ped_rec);
}


/* helper function for get_llooprec(), used to scan the loop list */
static void    *loop_list_iterate(list_data EntryData, void *Arg)
{
    if (((linkage_loop_rec *) EntryData)->Num == * (int *)Arg) {
        return (void *) EntryData;
    }
    return (void *)NULL;
}


/* get_llooprec()
 *
 * Get the loop record number Num. If it doesn't exist, create it.
 * Note that the list may need to be created, too.
 */
linkage_loop_rec *get_llooprec(listhandle **List, int Num)

{

    linkage_loop_rec *loop = NULL;

    if (*List == NULL)    {
        /* take a short-cut */
        *List = new_list();
        loop = new_llooprec();
        loop->Num = Num;
        loop->numEntry=0;
        append_to_list_tail(*List, (list_data) loop);

    } else {
	/* list_iterate's second argument is a void pointer */
        loop = (linkage_loop_rec *) list_iterate(*List,(void *) &Num,
                                                 loop_list_iterate);
        if (loop == NULL) {
            loop = new_llooprec();
            loop->Num = Num;
            loop->numEntry=0;
            append_to_list_tail(*List, (list_data) loop);
        }
    }
    return loop;
}

static void    *lcollapse_iterate(list_data EntryData, void *It)

{

    id_pair_type *this_id = (id_pair_type *) EntryData;

    if (*(int *)It == this_id->old)
        *(int *)It = this_id->new_linkage;

    return NULL;	/* continue */
}


/* collapse_lpedtree()
 *
 * This routine frees all entries with IDs set to UNDEF. It then reallocates
 * a new Entry[] array, copies all the entries over, and fixes their
 * ID's, parents, sibs, and children. It then frees the old Entry[] array
 * and assigns the new one to the pedigree.
 *
 * This is called, for example, by connect_loops(), since that function
 * has to delete duplicate entries.
 */
static void    collapse_lpedtree(linkage_ped_tree *Ped,
				 linkage_ped_top *Top)
{
    int entry, entry1;
    int newentry;
    int NewEntryCnt;
    linkage_ped_rec *Entry, *NEntry;
    linkage_ped_rec *NewEntry;

    id_pair_type *id_pair;
    listhandle     *Modified = new_list();

    int offset = Top->LocusTop->PhenoCnt;

    /* count the number of entries to be saved */
    /* and construct the list of ID's of deleted entries */
    NewEntryCnt = 0;
    for (entry = 0; entry < Ped->EntryCnt; entry++) {
        if (Ped->Entry[entry].ID != UNDEF)  {
            NewEntryCnt++;
        }
    }
    /* create the new array */
    NewEntry = CALLOC((size_t) NewEntryCnt, linkage_ped_rec);
    newentry = 0;
    for (entry = 0; entry < Ped->EntryCnt; entry++) {
        Entry = &(Ped->Entry[entry]);
        if (Entry->ID == UNDEF)
            continue;
        NEntry = &(NewEntry[newentry]);

        copy_lpedrec(Entry, NEntry, Top);
        Entry->Pheno  = NULL;
        Entry->Marker = NULL;
        free_all_from_lpedrec(Entry, offset);
        /* Change the field numbers of entry to match the modified list */
        /* father, mother, first-off, pa_sib, ma_sib */
        list_iterate(Modified, (void *) &(NEntry->Father), lcollapse_iterate);
        list_iterate(Modified, (void *) &(NEntry->Mother), lcollapse_iterate);
        list_iterate(Modified, (void *) &(NEntry->First_Offspring),
                     lcollapse_iterate);
        list_iterate(Modified, (void *) &(NEntry->Next_PA_Sib),
                     lcollapse_iterate);
        list_iterate(Modified, (void *) &(NEntry->Next_MA_Sib),
                     lcollapse_iterate);
        /* change the fields of all other entries to match the new ID */
        if (NEntry->ID != (newentry+1)) {
            for(entry1=0; entry1 < Ped->EntryCnt; entry1++) {
                /* father */
                if (NewEntry[entry1].Father == NEntry->ID) {
                    NewEntry[entry1].Father = newentry+1;
                }
                /* mother */
                if (NewEntry[entry1].Mother == NEntry->ID) {
                    NewEntry[entry1].Mother = newentry+1;
                }
                /* First_Offspring */
                if (NewEntry[entry1].First_Offspring == NEntry->ID) {
                    NewEntry[entry1].First_Offspring = newentry+1;
                }
                /* Next_PA_Sib */
                if (NewEntry[entry1].Next_PA_Sib == NEntry->ID) {
                    NewEntry[entry1].Next_PA_Sib = newentry+1;
                }
                /* Next_MA_Sib */
                if (NewEntry[entry1].Next_MA_Sib == NEntry->ID) {
                    NewEntry[entry1].Next_MA_Sib = newentry+1;
                }
            }
            id_pair = CALLOC((size_t) 1, id_pair_type);
            id_pair->new_linkage = newentry+1;
            id_pair->old = NEntry->ID;
            append_to_list_tail(Modified, (list_data) id_pair);
        }
        NEntry->ID = newentry+1;
        newentry++;
    }

    /* also change the proband field of the pedigree */
    list_iterate(Modified, (void *) &(Ped->Proband), lcollapse_iterate);
    list_free_all(Modified, free);
    free(Ped->Entry);
    Ped->Entry = NewEntry;
    Ped->EntryCnt = NewEntryCnt;
    return;
}


static void  log_loop_person(linkage_ped_rec *D_entry,
			     linkage_ped_rec *S_entry,
			     linkage_ped_tree *Ped)
{

    sprintf(err_msg, "Ped %s: Loop-breaker ids %s and %s",
            Ped->Name, S_entry->UniqueID, D_entry->UniqueID);
    mssgf(err_msg);
    sprintf(err_msg, "        Saving %s, deleting %s",
            S_entry->UniqueID, D_entry->UniqueID);
    mssgf(err_msg);
    return;
}

/* connect_loops()
 *
 * This routine goes through all the loops and reconnects them.
 * It collapses the pedigree afterward.
 *
 * Return -1 if there was a fatal error, 0 otherwise.
 * ****Modified by Nandita - (2/17/00), to connect multiple loops
 *
 */

int  connect_loops(linkage_ped_tree *Ped, linkage_ped_top *Top1)

{
    linkage_loop_rec *loop1, *Loop1=NULL, *Loop2=NULL;
    linkage_ped_rec *D_Entry = NULL, *S_Entry;
    int lb, i, save, *delete_ppl, found_loop=0;
    int j;

    int offset = Top1->LocusTop->PhenoCnt;

    if (Ped->Loops == NULL) {
        /* now collapse the Entry array */
        collapse_lpedtree(Ped, Top1);
        return 0;
    }


    delete_ppl=CALLOC((size_t) Ped->EntryCnt, int);
    for (j=0; j < Ped->EntryCnt; j++) delete_ppl[j]=0;

    while (loop1 = (linkage_loop_rec *) pop_first_list_entry(Ped->Loops),
           loop1 != NULL) {
        if (loop1->Num == 0 && loop1->numEntry == 1) {
            Loop1 = loop1;
            continue;
        }
        if (loop1->Num == 1 && loop1->numEntry == 1) {
            /* also a singleton */
            Loop2 = loop1;
            continue;
        }

        if (found_loop == 0) {
            sprintf(err_msg, "Connecting loops in pedigree %d...", Ped->Num);
            mssgf(err_msg);
            found_loop=1;
        }

        /* else */
        save = -1;
        for (i=0; i<loop1->numEntry; i++) {
            /* identify which loop breaker to keep */
            if (loop1->Entry[i] <= 0)  {
                errorf("Incomplete loop record! Are all members properly labelled?");
                sprintf(err_msg, "Entry id: %s", Ped->Entry[loop1->Entry[i]-1].UniqueID);
                errorf(err_msg);
                return -1;
            } else  {
                /* one to keep has non_zero parents*/
                /* (by picking it this way we won't have to correct any parents) */
                if ((Ped->Entry[loop1->Entry[i] - 1].Mother > 0)
                    && (Ped->Entry[loop1->Entry[i] - 1].Father > 0))   {
                    if (save > -1) {
                        errorf("More than one loop-breaker with non-zero parents.");
                        EXIT(DATA_INCONSISTENCY);
                    }
                    save=i;
                }
            }
            if (save == -1) {
                /* all clones are founders */
                save=0; /* just save the first entry */
            }
        }

        S_Entry = &(Ped->Entry[loop1->Entry[save] - 1]);
        S_Entry->loopbreakers = CALLOC((size_t) MAXLOOPMEMBERS+1, int);
        lb=0;
        for (i=0; i<loop1->numEntry; i++) {
            /* compare the saved entry with genotypes of other equivalent entries */
            if (save == i) {
                continue;
            }
            delete_ppl[loop1->Entry[i] - 1]=1;
            D_Entry = &(Ped->Entry[loop1->Entry[i] - 1]);
            S_Entry->loopbreakers[lb]=D_Entry->ID; lb++;
            connect_loopbreakers(S_Entry, D_Entry, Ped, Top1);
        }
        S_Entry->loopbreakers[lb]=-99;    /* don't forget to free the loop record */
        free(loop1);
        loop1 = NULL;
    }
    /* Now process the Loop1 and Loop2 */
    if (Loop1 != NULL && Loop2 != NULL) {
        if (found_loop == 0) {
            sprintf(err_msg, "Connecting loops in pedigree %d...", Ped->Num);
            mssgf(err_msg);
        }
        delete_ppl[Loop2->Entry[0] - 1] = 1;
        D_Entry = &(Ped->Entry[Loop2->Entry[0] - 1]);
        S_Entry = &(Ped->Entry[Loop1->Entry[0] - 1]);
        connect_loopbreakers(S_Entry, D_Entry, Ped, Top1);
        /*    S_Entry->OrigProband=0; */
        free(Loop1); Loop1=NULL;
        free(Loop2); Loop2=NULL;
    }

    if (Loop1 != NULL) free(Loop1);

    /* free the memory associated with all deleted entries */
    for (j=0; j<Ped->EntryCnt; j++) {
        if (delete_ppl[j] == 1) { // implies D_Entry was defined via above code
//     free_all_from_lpedrec() must ignore Pheno&Marker
            D_Entry->Pheno = 0;
            D_Entry->Marker = 0;
            free_all_from_lpedrec(D_Entry, offset);
            /* mark it as deleted */
            Ped->Entry[j].ID = UNDEF;
        }
    }

    /* delete the loops list (it has no more entries) */
    free_list(Ped->Loops);
    Ped->Loops = NULL;

    /* now collapse the Entry array */
    collapse_lpedtree(Ped, Top1);
    free(delete_ppl);
    return 0;
}

void count_lgenotypes(linkage_ped_top *Top, size_t *num_inds,
		      size_t *males, size_t *females,
		      size_t *num_untyped, size_t *num_typed,
		      size_t *peds_typed, size_t *males_typed,
		      size_t *females_typed, size_t *num_half_typed)
{
    int i, j;
    size_t k, l;
    size_t male_count=0, female_count=0, individual_count=0;
    size_t typed=0, half_typed=0, untyped=0;
    size_t this_person_typed, this_male_typed, this_female_typed, this_ped_typed;
    size_t numloc;
    int first_time=1;

    *peds_typed = *males_typed = *females_typed = 0;
    *num_half_typed = 0;
    strcpy(err_msg, "");
    SECTION_LOG_INIT(untyped_pedigree);
    for (i = 0; i < Top->PedCnt; i++) {
        this_ped_typed = 0;
        individual_count = individual_count + Top->Ped[i].EntryCnt;
        for (j = 0; j < Top->Ped[i].EntryCnt; j++)  {
            linkage_ped_rec *pp = &(Top->Ped[i].Entry[j]);
            int gender = pp->Sex;
            if (gender == MALE_ID) {
                male_count++;
            } else if (gender == FEMALE_ID) {
                female_count++;
            }
            Top->Ped[i].Entry[j].Ngeno = 0;
            Top->Ped[i].Entry[j].IsTyped = 0;
            this_male_typed=0; this_female_typed=0;

            if (Mega2Status < LOCI_REORDERED || database_dump) {
                numloc = Top->LocusTop->LocusCnt;
            } else {
                numloc = num_reordered;
            }
            for (l = 0; l < numloc; l++)   {
                if (Mega2Status < LOCI_REORDERED || database_dump) {
                    k = l;
                } else {
                    k = reordered_marker_loci[l];
                }
                if (Top->LocusTop->Locus[k].Type != AFFECTION &&
                    Top->LocusTop->Locus[k].Type != QUANT)    {
                    this_person_typed=0;
                    /* why
                    if (Top->LocusTop->Locus[k].Type == YLINKED && gender == FEMALE_ID) {
                        continue;
                    }
                    */
                    if (Mega2Status <= INSIDE_RECODE && Top->LocusTop->PedRecDataType == Raw_postmake) {
                        linkage_ped_rec *pp = Top->PedBroken[i].Entry + j;
                        this_person_typed = num_typed_2Ralleles(pp->Marker, k);

                    } else {
                        linkage_ped_rec *pp = Top->PedBroken[i].Entry + j;
                        this_person_typed = num_typed_2alleles(pp->Marker, k);
                    }
                    Top->Ped[i].Entry[j].IsTyped += this_person_typed;
                    /* If we are in recode, we are counting half-types,
                       otherwise this rouitne is called after half-typed
                       individuals have been reset or not reset */
                    if (this_person_typed) {
                        Top->Ped[i].Entry[j].Ngeno += 1;
                        if (gender == MALE_ID) {
                            this_male_typed = 1;
                        } else if (gender == FEMALE_ID) {
                            this_female_typed = 1;
                        }
                        if (this_person_typed == 2) {
                            typed++;
                        } else {
                            half_typed++;
                        }
                        this_ped_typed = 1;
                    }
                    else
                        untyped++;
                }
            }
            if (Top->Ped[i].Entry[j].IsTyped == (int)(2 * numloc) )
                Top->Ped[i].Entry[j].IsTyped = 2;
            else if (Top->Ped[i].Entry[j].IsTyped)
                Top->Ped[i].Entry[j].IsTyped = 1;

            if (gender == MALE_ID) {
                *males_typed += this_male_typed;
            } else if (gender == FEMALE_ID) {
                *females_typed += this_female_typed;
            }
        }

        if (! this_ped_typed) {
            /* output the pedigree number into log file */
            if (strlen(err_msg) > 70) {
                if (first_time) {
                    SECTION_LOG(untyped_pedigree);
                    mssgf("------------------------------------------------------------");
                    mssgf("Completely untyped pedigrees:");
                    SECTION_LOG_HEADER(untyped_pedigree);
                    first_time=0;
                }
                SECTION_LOG(untyped_pedigree);
                mssgf(err_msg);
                strcpy(err_msg, "");
            }
            strcat(err_msg, Top->Ped[i].Name); strcat(err_msg, " ");
        }
        *peds_typed += this_ped_typed;
    }

    if (strlen(err_msg) > 0) {
        if (first_time) {
            SECTION_LOG(untyped_pedigree);
            mssgf("------------------------------------------------------------");
            mssgf("Completely untyped pedigrees:");
            SECTION_LOG_HEADER(untyped_pedigree);
            first_time=0;
        }
        SECTION_LOG(untyped_pedigree);
        mssgf(err_msg);
        strcpy(err_msg, "");
    }
    SECTION_LOG_FINI(untyped_pedigree);

    *num_inds=individual_count;
    *males = male_count;
    *females = female_count;
    *num_untyped = untyped;
    *num_typed = typed;
    *num_half_typed = half_typed;
}

static void connect_loopbreakers(linkage_ped_rec *S_Entry,
				 linkage_ped_rec *D_Entry,
				 linkage_ped_tree *Ped,
				 linkage_ped_top *Top1)
{
    int j, locus1;
    int sa1, sa2, da1, da2;

    /* Check that genotypes match between loop person and duplicate */
    for (locus1 = 0; locus1 < Top1->LocusTop->LocusCnt; locus1++)
        switch (Top1->LocusTop->Locus[locus1].Type)   {
        case QUANT:
            if (fabs(S_Entry->Pheno[locus1].Quant - D_Entry->Pheno[locus1].Quant) >= EPSILON) {
                sprintf(err_msg, "Ped %d, Trait locus %s:",
                        Ped->Num, Top1->LocusTop->Locus[locus1].LocusName);
                errorf(err_msg);
                sprintf(err_msg, "Loop person %d has phenotype %10.7f",
                        S_Entry->ID, S_Entry->Pheno[locus1].Quant);
                errorf(err_msg);
                sprintf(err_msg, "Duplicate person %d has phenotype %10.7f",
                        D_Entry->ID, D_Entry->Pheno[locus1].Quant);
                errorf(err_msg);
                EXIT(DATA_INCONSISTENCY);
            }
            break;
        case AFFECTION:
            if (S_Entry->Pheno[locus1].Affection.Status != D_Entry->Pheno[locus1].Affection.Status) {
                sprintf(err_msg, "Ped %d, Trait locus %s:",
                        Ped->Num, Top1->LocusTop->Locus[locus1].LocusName);
                errorf(err_msg);
                sprintf(err_msg, "Loop person %d has status %d",
                        S_Entry->ID, S_Entry->Pheno[locus1].Affection.Status);
                errorf(err_msg);
                sprintf(err_msg, "Duplicate person %d has status %d",
                        D_Entry->ID, D_Entry->Pheno[locus1].Affection.Status);
                errorf(err_msg);
                EXIT(DATA_INCONSISTENCY);
            }
            if (Top1->LocusTop->Pheno[locus1].Props.Affection.ClassCnt > 1) {
                if (S_Entry->Pheno[locus1].Affection.Class != D_Entry->Pheno[locus1].Affection.Class) {
                    sprintf(err_msg, "Ped %d, Trait locus %s:",
                            Ped->Num, Top1->LocusTop->Locus[locus1].LocusName);
                    errorf(err_msg);
                    sprintf(err_msg, "Loop person %d has class %d",
                            S_Entry->ID, S_Entry->Pheno[locus1].Affection.Class);
                    errorf(err_msg);
                    sprintf(err_msg, "Duplicate person %d has class %d",
                            D_Entry->ID, D_Entry->Pheno[locus1].Affection.Class);
                    errorf(err_msg);
                    EXIT(DATA_INCONSISTENCY);
                }
            }

            break;
        case BINARY:
        case NUMBERED:
        case XLINKED:
        case YLINKED:
            get_2alleles(S_Entry->Marker, locus1, &sa1, &sa2);
            get_2alleles(D_Entry->Marker, locus1, &da1, &da2);
            if (sa1 != da1 || sa2 != da2) {
                if (Top1->LocusTop->Marker[locus1].chromosome > 0) {
		    Display_Errors=1;
                    sprintf(err_msg, "Ped %d, Chr %d, Locus %d:",
                            Ped->Num,
                            Top1->LocusTop->Marker[locus1].chromosome, locus1);
                    errorf(err_msg);
                    // Since this is linkage we can only have numeric alleles...
                    sprintf(err_msg, "Loop person %d has genotype %d/%d",
                            S_Entry->ID, sa1, sa2);
                    errorf(err_msg);
                    sprintf(err_msg, "Duplicate person %d has genotype %d/%d",
                            D_Entry->ID, da1, da2);
                    errorf(err_msg);
                    EXIT(DATA_INCONSISTENCY);
                }
            }
            break;
        default:
            warnvf("unknown locus type (locus %d) while writing output\n", locus1 + 1);
        }

    /* reset the parent pointers  in the pedigree */
    for (j=0; j < Ped->EntryCnt; j++) {
        if (Ped->Entry[j].Father == D_Entry->ID) {
            Ped->Entry[j].Father = S_Entry->ID;
            if (S_Entry->First_Offspring == 0) {
                S_Entry->First_Offspring = Ped->Entry[j].ID;
            }
        }
        if (Ped->Entry[j].Mother == D_Entry->ID) {
            Ped->Entry[j].Mother = S_Entry->ID;
            if (S_Entry->First_Offspring == 0) {
                S_Entry->First_Offspring = Ped->Entry[j].ID;
            }
        }
    }
    /* Why are we setting the saved entry values to Deleted entry value? */
    /* Nandita - we were deleting the member with parents,
       now we are deleting the member without parents.
       S_Entry->Next_PA_Sib = D_Entry->Next_PA_Sib;
       S_Entry->Next_MA_Sib = D_Entry->Next_MA_Sib;
    */

    if (Ped->Proband == D_Entry->ID) {
        Ped->Proband = S_Entry->ID;
    }
    log_loop_person(D_Entry, S_Entry, Ped);
    return;
}

//
// Used by get_loci_on_chromosome() to indicate the first and last index of the markers
// in the ChrLoci array. The markers are copied together, but will not necessarly be at the
// beginning or the end of the array.
int ChrLociFirstMarkerIndex, ChrLociLastMarkerIndex;

/**
   @brief Find the loci on a chromosome

   Create locus list for a list of chromosomes, setting ChrLoci[], and NumChrLoci.
   Check the global_trait_entries, and include traits.

   @param[in]   the chromocome number (-1 is says all).
   @return      uses globals ChrLoci[], and NumChrLoci.
*/
void get_loci_on_chromosome(const int numchr)
{
    int chromindex = 0, upper, lower; //silly compiler
    int i, j, k;

    if (ManualReorder == 1) return;

    // Determine the number of loci to consider.
    // Traits are included if 'num_traits' has been set.

    NumChrLoci=0;
    NumChrSite=0;
    if (num_traits > 0) {
//c      NumChrLoci = num_traits - ((LoopOverTrait == 0) ? 1 : 0);// This breaks when there is no '-1' in global_trait_entries[].
// It looks like the problem occurs in the batch_run tests 'change_chrom' where:
//        sed -e s/Chromosome_Single=.*/Chromosome_Single=$3/g $file  > $target
//        echo "Trait_Single=1" >> $target
// In this case: global_trait_entries == <0, -99>; LoopOverTrait == 0; num_traits == 1;
//r      NumChrLoci = num_traits - (global_trait_entries[num_traits - 1] == -1 ? 1 : 0);
        // Here we are adjusting NumChrLoci to hold the count of the actual number of traits that
        // we will be working with. In addition to trait indexes, the array global_trait_entries may
        // contain a flag where the value is == -1 within the range 0 < num_traits. The infamous -99
        // flag is a centenal and SHOULD occur at the index num_traits, but not within.
        NumChrLoci = num_traits;
        for (i = 0; i < num_traits; i++)
            if (global_trait_entries[i] < 0) NumChrLoci--;
    }

    // if length of reordered_marker_loci greater than zero...
    if (num_reordered > 0) {
        if (numchr == 0) {
            NumChrLoci += num_reordered;
            NumChrSite += num_reordered;
        } else if (numchr == -1) {
            /* everything but the x-chromosome */
            for (i=0; i < main_chromocnt; i++) {
                if (global_chromo_entries[i] == SEX_CHROMOSOME) {
                    int locuscnt = chromo_loci_count[i] - ((i > 0) ? chromo_loci_count[i-1] : 0);
                    NumChrLoci += (num_reordered - locuscnt);
                    NumChrSite += (num_reordered - locuscnt);
                    chromindex = i;
                    break;
                }
            }
        } else {
            // for everything that the user selected...
            for (i=0; i < main_chromocnt; i++) {
                if (numchr == global_chromo_entries[i]) {
                    // chromo_loci_count is cumulative, so only the 0th entry contains the correct number...
                    int locuscnt = chromo_loci_count[i] - ((i > 0) ? chromo_loci_count[i-1] : 0);
                    NumChrLoci += locuscnt;
                    NumChrSite += locuscnt;
                    chromindex = i;
                    break;
                }
            }
            /* not handling invalid chromo number */
        }
    }

    if (num_covariates > 0 && covariates != NULL) {
        NumChrLoci += num_covariates;
    }

    // Deallocate (if necessary) and reallocate...
    if (ChrLoci != NULL) free(ChrLoci);
    ChrLoci = CALLOC((size_t)NumChrLoci, int);
    
    // Currently only used by Beagle to sort the markers...
    ChrLociFirstMarkerIndex = ChrLociLastMarkerIndex = -1;

    // Populate the loci array...
    // NOTE: num_traits should be the length of global_trait_entries...
    int chr_copied = false;
    for (j=0, i=0; i < num_traits; i++) {
        if (global_trait_entries[i] >= 0) {
            ChrLoci[j++] = global_trait_entries[i];
        } else {
            chr_copied = true;
            if (ChrLociFirstMarkerIndex == -1) ChrLociFirstMarkerIndex = j;
            if (num_reordered > 0) {
                if (numchr == 0) {
                    for (k=0; k < num_reordered; k++) {
                        ChrLoci[j++] = reordered_marker_loci[k];
                    }
                    ChrLociLastMarkerIndex = j-1;
                } else if (numchr == -1) {
                    /* skip the loci on X-chromosome */
                    if (chromindex == 0) {
                        lower=0; upper = chromo_loci_count[chromindex];
                    } else {
                        lower=chromo_loci_count[chromindex-1];
                        upper = chromo_loci_count[chromindex];
                    }
                    for (k=0; k < num_reordered; k++) {
                        if (k < lower || k >= upper) {
                            ChrLoci[j++] = reordered_marker_loci[k];
                        }
                    }
                    ChrLociLastMarkerIndex = j-1;
                } else {
                    if (chromindex == 0) {
                        lower=0; upper = chromo_loci_count[chromindex];
                    } else {
                        lower=chromo_loci_count[chromindex-1];
                        upper = chromo_loci_count[chromindex];
                    }
                    for (k = lower; k < upper; k++) {
                        ChrLoci[j++] = reordered_marker_loci[k];
                    }
                    ChrLociLastMarkerIndex = j-1;
                }
            }
        }
    }

    if (num_reordered > 0) {
//      if (LoopOverTrait != 0 || num_traits == 0)
        if (! chr_copied)
        {
            /* insert the markers now */
            if (ChrLociFirstMarkerIndex == -1) ChrLociFirstMarkerIndex = j;
            if (numchr == 0) {
                for (k=0; k < num_reordered; k++) {
                    ChrLoci[j++] = reordered_marker_loci[k];
                }
                ChrLociLastMarkerIndex = j-1;
            } else {
                if (chromindex == 0) {
                    lower=0; upper = chromo_loci_count[chromindex];
                } else {
                    lower=chromo_loci_count[chromindex-1];
                    upper = chromo_loci_count[chromindex];
                }
                for (k = lower; k < upper; k++) {
                    ChrLoci[j++] = reordered_marker_loci[k];
                }
                ChrLociLastMarkerIndex = j-1;
            }
        }
    } else {
        warnf("No markers to save.");
    }

    if (num_covariates > 0 && covariates != NULL) {
        for (k=0; k < num_covariates; k++) {
            ChrLoci[j++] = covariates[k];
        }
    }

#ifdef DEBUG_MEGA2
    printf("Copied over %d loci\n", NumChrLoci);
#endif
}

void get_unmapped_loci(int append)

{
    int i,j,k;

    if (append == 1) {
        ChrLoci = REALLOC(ChrLoci, (size_t)NumChrLoci+NumUnmapped, int);
        for (i=0; i < NumUnmapped; i++) {
            ChrLoci[NumChrLoci + i] = unmapped_markers[i];
        }
        NumChrLoci += NumUnmapped;
    } else {
        if (ManualReorder == 1) {
            return;
        }

        if (ChrLoci != NULL) {
            free(ChrLoci);
        }
        if (append == 2) {
            /* first store the traits in Combine mode */
            NumChrLoci = (LoopOverTrait == 1 || num_traits == 0)? num_traits : num_traits - 1;
            NumChrLoci += NumUnmapped;
        } else {
            NumChrLoci = NumUnmapped;
        }
        ChrLoci = CALLOC((size_t)NumChrLoci, int);

        if (append == 2) {
            int chr_copied = false;
            j=0;
            for (i=0; i < num_traits; i++) {
                if (global_trait_entries[i] >= 0) {
                    ChrLoci[j++] = global_trait_entries[i];
                } else {
                    chr_copied = true;
                    for (k=0; k < NumUnmapped; k++) {
                        ChrLoci[j++] = unmapped_markers[k];
                    }
                }
            }
            /* Now insert the markers if LoopOverTrait = 1 */
//          if (LoopOverTrait == 1 || num_traits == 0)
            if (! chr_copied)
            {
                for (k=0; k < NumUnmapped; k++) {
                    ChrLoci[j++] = unmapped_markers[k];
                }
            }
        } else {
            for (i=0; i < NumUnmapped; i++) {
                ChrLoci[i] = unmapped_markers[i];
            }
        }
    }

    return;
}

void markers_on_chromosome(int numchr)
{
    if (numchr == -1) {
        int i;
        get_loci_on_chromosome(0);
        for (i=0; i < main_chromocnt; i++) {
            if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
                get_unmapped_loci(1);
                break;
            }
        }
    } else if (numchr == UNKNOWN_CHROMO)
        get_unmapped_loci(0);
    else
        get_loci_on_chromosome(numchr);
}

void free_chr_loci(void)
{
    NumChrLoci=0;
    free(ChrLoci);
}

int is_typed_lentry(linkage_ped_rec Entry, linkage_locus_top *LTop,
		    int mode, double *percent_typed)

{

    int loc, m;
    int *markers, num_markers = 0;
    int typed=0;


    for(loc=0; loc < NumChrLoci; loc++) {
        if (LTop->Locus[ChrLoci[loc]].Type == NUMBERED ||
           LTop->Locus[ChrLoci[loc]].Type == BINARY) {
            num_markers++;
        }
    }
    /* Possible error: Why should this return 1 = 'true' if there are zero markers? */
    if (num_markers <= 0) {
        return 1;
    }

    markers = CALLOC((size_t)num_markers, int);

    m=0;
    for(loc=0; loc < NumChrLoci; loc++) {
        if (LTop->Locus[ChrLoci[loc]].Type == NUMBERED ||
           LTop->Locus[ChrLoci[loc]].Type == BINARY) {
            markers[m]=ChrLoci[loc]; m++;
        }
    }

    for(m=0; m < num_markers; m++) {
        int a1, a2;
        get_2alleles(Entry.Marker, markers[m], &a1, &a2);
        if (a1 != 0 || a2 != 0)
            typed++;
    }

    free(markers);

    if (percent_typed != NULL) {
        *percent_typed = (double) (typed/num_markers);
    }

    /* Logical expressions are defined to have value 1 if true, and 0 if false */

    if (mode == TYPED_AT_ANY_LOCUS) {
        return((typed > 0)? 1 : 0);
    } else {
        return((typed < num_markers)? 0: 1);
    }

}

void clean_reordered_markers(linkage_locus_top *LTop, analysis_type analysis)
{
    int i, j, upper, lower;
    int num_valid, *tmp=CALLOC((size_t)num_reordered, int);
    int *chromo_loci_final_count = CALLOC((size_t)main_chromocnt, int);
    
    HasMarkers=0;
    
    num_valid=0;
    for (i=0; i < main_chromocnt; i++) {
        if (i==0) {
            lower=0;
            chromo_loci_final_count[i] = chromo_loci_count[i];
        } else {
            lower = chromo_loci_count[i-1];
            chromo_loci_final_count[i] = chromo_loci_count[i] - chromo_loci_count[i-1];
        }
        upper = chromo_loci_count[i];
        for (j=lower; j < upper; j++) {
            /* switch types */
            if (LTop->Locus[reordered_marker_loci[j]].Type == XLINKED ||
                LTop->Locus[reordered_marker_loci[j]].Type == YLINKED ||
                LTop->Locus[reordered_marker_loci[j]].Type == NUMBERED) {
                LTop->Locus[reordered_marker_loci[j]].Type = NUMBERED;
                HasMarkers=1;
            }
            /*     Check positions */
            if (LTop->Locus[reordered_marker_loci[j]].Type == NUMBERED) {
                int chr = LTop->Marker[reordered_marker_loci[j]].chromosome;
#ifdef USEOLDMAPCODE
                if (LTop->Marker[reordered_marker_loci[j]].position >= 0.0 ||
                    ALLOW_NO_MAP(analysis)) {
                    tmp[num_valid] = reordered_marker_loci[j];
                    num_valid++;
                } else if (LTop->Marker[reordered_marker_loci[j]].pos_male >= 0.0 &&
                           LTop->Marker[reordered_marker_loci[j]].pos_female >= 0.0) {
                    tmp[num_valid] = reordered_marker_loci[j];
                    num_valid++;
                    // compute the missing sex-averaged position by taking the averaged of
                    // the sex-specific information?!
                    // Likely this means that the code is making an assumption about having
                    // sex-averated information later in the code. That is troubling...
                    LTop->Marker[reordered_marker_loci[j]].pos_avg  =
                    (LTop->Marker[reordered_marker_loci[j]].pos_male +
                     LTop->Marker[reordered_marker_loci[j]].pos_female)/2.0;
                } else {
                    warnvf("Missing map position for marker %s, it will be excluded from analysis.\n",
                            LTop->Marker[reordered_marker_loci[j]].MarkerName);
                    chromo_loci_final_count[i]--;
                }
#else /* USEOLDMAPCODE */
                // Since the user has the option now of explicitly choosing a sex-averaged,
                // sex-specific, or female map even though all three maps (a, m, f) are present,
                // we must take this into consideration. We must also worry about the sequelae of this...
//              asm("int $3");
                if (ALLOW_NO_MAP(analysis) || analysis->allow_no_genetic_map() ||
                    // with a sex-averaged map, we only know about the average position...
                    (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT &&
                     LTop->Marker[reordered_marker_loci[j]].pos_avg >= 0.0) ||
                    // with a sex-specific map, we must have both male an female for this genetic marker...
                    (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT &&
                     // only test or a male position if this is not the X-chromosome...
                     (chr != SEX_CHROMOSOME ? LTop->Marker[reordered_marker_loci[j]].pos_male >= 0.0 : 1) &&
                     LTop->Marker[reordered_marker_loci[j]].pos_female >= 0.0) ||
                    // with a female only map, we must have the femele for this marker, and it 's only valid on the X chromosome...
                    (genetic_distance_sex_type_map == FEMALE_GDMT &&
                     LTop->Marker[reordered_marker_loci[j]].pos_female >= 0.0 &&
                     LTop->Marker[reordered_marker_loci[j]].chromosome == SEX_CHROMOSOME)
                    ) {
                    tmp[num_valid] = reordered_marker_loci[j];
                    num_valid++;
                } else {
                    warnvf("Missing or unusable map position for marker %s, it will be excluded from analysis.\n",
                            LTop->Marker[reordered_marker_loci[j]].MarkerName);
                    chromo_loci_final_count[i]--;
                }
#endif /* USEOLDMAPCODE */
            }
        }
    }

    // ??????? Put a breakpoint here. Make a marker go away from the .MAP file and that marker should
    // not show up in the 'tmp' vector and therefore never copied into the 'reordered_marker-loci' vector.
    
    if (num_valid > 0) {
        reordered_marker_loci = REALLOC(reordered_marker_loci, (size_t)num_valid, int);
        HasMarkers = 1;
        for (i=0; i < num_valid; i++) {
            reordered_marker_loci[i] = tmp[i];
        }
    } else {
        free(reordered_marker_loci);
        reordered_marker_loci = NULL;
    }
    
    /* set the chromo_loci_counts */
    for (i=0; i < main_chromocnt; i++) {
        if (i==0) {
            chromo_loci_count[i] = chromo_loci_final_count[i];
        } else {
            chromo_loci_count[i] = chromo_loci_count[i-1] + chromo_loci_final_count[i];
        }
    }
    
    num_reordered = num_valid;
    /*   printf("%d, num_reordered\n.", num_reordered); */
    free(tmp); free(chromo_loci_final_count);
}

void set_link_IDs(linkage_ped_top *Top)
{
    int ped, entry;
    int entry1, l;
    linkage_ped_rec *Entry;
    linkage_ped_tree Ped;
    int **ids;

    for (ped=0; ped < Top->PedCnt; ped++) {
        ids = CALLOC((size_t) Top->Ped[ped].EntryCnt, int*);
        for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
            ids[entry] = CALLOC((size_t) 5, int);
            for (l=0; l < 5; l++) {
                ids[entry][l] = 0;
            }
        }
        Ped = Top->Ped[ped];
        for (entry = 0; entry < Ped.EntryCnt; entry++) {
            Entry = &(Ped.Entry[entry]);
            for (entry1=0; entry1 < (int) Ped.EntryCnt; entry1++) {
                /* Father */
                if (Entry->Father > 0 && Entry->Father == Ped.Entry[entry1].ID) {
                    ids[entry][0] = entry1+1;
                }
                /* Mother */
                if (Entry->Mother > 0 && Entry->Mother == Ped.Entry[entry1].ID) {
                    ids[entry][1] = entry1+1;
                }
                /* First_Offspring */
                if (Entry->First_Offspring > 0 && Entry->First_Offspring == Ped.Entry[entry1].ID) {
                    ids[entry][2] = entry1+1;
                }
                /* Next_PA_Sib */
                if (Entry->Next_PA_Sib > 0 && Entry->Next_PA_Sib == Ped.Entry[entry1].ID ) {
                    ids[entry][3] = entry1+1;
                }
                /* Next_MA_Sib */
                if (Entry->Next_MA_Sib > 0 && Entry->Next_MA_Sib == Ped.Entry[entry1].ID) {
                    ids[entry][4] = entry1+1;
                }
            }
        }

        for (entry = 0; entry < Ped.EntryCnt; entry++) {
            Entry = &(Ped.Entry[entry]);

            if (ids[entry][0] > 0) {
                Entry->Father = ids[entry][0];
            }

            if (ids[entry][1] > 0) {
                Entry->Mother = ids[entry][1];
            }

            if (ids[entry][2] > 0) {
                Entry->First_Offspring = ids[entry][2];
            }

            if (ids[entry][3] > 0) {
                Entry->Next_PA_Sib = ids[entry][3];
            }

            if (ids[entry][4] > 0) {
                Entry->Next_MA_Sib = ids[entry][4];
            }

            free(ids[entry]);
        }
        free(ids);
    }
    return;
}
