 /*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#include "pedtree.h"

#include "error_messages_ext.h"
#include "linkage_ext.h"
/*
     error_messages_ext.h:  my_calloc my_malloc my_realloc
            linkage_ext.h:  is_typed_lentry
*/


/*        Exported functions           */

int set_locus_orderDEFUNCT(ped_top *Top, int *order);
void clear_ped_tree( ped_tree *Tree);
ped_top *new_ped_top(void);
int is_typed(ped_rec *Entry, locus_top *LTop1, int mode);
ped_top        *convert_to_pedtree(linkage_ped_top *Top, int assign_affs);
void free_all_including_ped_top(ped_top *Top,  void (*PTmpFreeFun)(void *TmpData),
				void (*ETmpFreeFun)(void *TmpData));
void free_all_from_ped_tree(ped_tree *Ped,
                            void (*PTmpFreeFun)(void *TmpData),
                            void (*ETmpFreeFun)(void *TmpData));
void free_all_from_locus_top(locus_top *LTop);
void add_to_offspring(ped_rec *Parent, ped_rec *Child);
int renumber_ped(ped_tree *Ped);
void remove_untyped_affecteds(ped_tree *Ped, locus_top *LTop,
			      int mode, int *order, int n);
void reassign_affecteds(ped_top *PTop);

/********************************************/

/*        internal functions          */
static void swap_lociDEFUNCT(locus_rec *locus1, locus_rec *locus2);
static void clear_locus_top(locus_top *LTop);
static void clear_ped_top(ped_top *Top);
static void free_all_from_ped_rec(ped_rec *Entry,
				  void (*TmpFreeFun)(void *TmpData));
static void free_all_from_locus_rec(locus_rec *Locus);
static void free_all_from_ped_top(ped_top *Top, void (*PTmpFreeFun)(void *TmpData),
                                  void (*ETmpFreeFun)(void *TmpData));
static locus_top *new_locus_top(void);
static void clear_ped_rec(ped_rec *Entry);
static int sort_affecteds(ped_tree *Ped);

/* static void assign_ped_id(ped_rec *Entry, int *current_id); */
/**********************************************/

static void clear_locus_top(locus_top *LTop)
{
    LTop->LocusCnt = 0;
    LTop->Locus = NULL;
    LTop->SexLinked = 0;
}


static locus_top *new_locus_top(void)
{
    locus_top *newtop = MALLOC(locus_top);
    clear_locus_top(newtop);
    return newtop;
}


static void clear_ped_rec(ped_rec *Entry)
{
    Entry->ID = UNDEF;
    /*  strcpy(Entry->OrigID, ""); */
    Entry->Status = UNDEF;
    Entry->Father = NULL;
    Entry->Mother = NULL;
    Entry->OffspringCnt = 0;
    Entry->Offspring = NULL;
    Entry->PA_Sib = NULL;
    Entry->MA_Sib = NULL;
    Entry->Sex = UNDEF;
    Entry->Allele_1 = NULL;
    Entry->Allele_2 = NULL;
    Entry->TmpData = NULL;
    Entry->LEntry = NULL;
}

void clear_ped_tree( ped_tree *Tree)
{
    Tree->Name=NULL;
    Tree->Entry = NULL;
    Tree->EntryCnt = 0;
    Tree->AffectedCnt = 0;
    Tree->Affected = NULL;
    Tree->Proband = NULL;
    Tree->TmpData = NULL;
}


static void clear_ped_top(ped_top *Top)
{
    Top->PedCnt = 0;
    Top->PedTree = NULL;
    Top->LocusTop = NULL;
}


ped_top *new_ped_top(void) {
    ped_top *newtop = MALLOC(ped_top);

    clear_ped_top(newtop);
    return newtop;
}


static void free_all_from_ped_rec(ped_rec *Entry, void (*TmpFreeFun)(void *TmpData))

{
    if (Entry == NULL) return;
    if (Entry->Offspring != NULL)
        free(Entry->Offspring);
    if (Entry->Allele_1 != NULL)
        free(Entry->Allele_1);
    if (Entry->Allele_2 != NULL)
        free(Entry->Allele_2);
    if ((Entry->TmpData != NULL) && (TmpFreeFun != NULL))
        TmpFreeFun(Entry->TmpData);
}


void free_all_from_ped_tree(ped_tree *Ped, void (*PTmpFreeFun)(void *TmpData),
			    void (*ETmpFreeFun)(void *TmpData))

{
    int entry;
    if (Ped == NULL) return;
    if (Ped->Name != NULL) free(Ped->Name);
    if (Ped->Entry != NULL) {
        for (entry = 0; entry < Ped->EntryCnt; entry++)
            free_all_from_ped_rec(&(Ped->Entry[entry]), ETmpFreeFun);
        free(Ped->Entry);  Ped->Entry = NULL;
    }
    if (Ped->Affected != NULL)
        free(Ped->Affected);
    if ((Ped->TmpData != NULL) && (PTmpFreeFun != NULL))
        PTmpFreeFun(Ped->TmpData);
}


static void free_all_from_locus_rec(locus_rec *Locus)
{
    if (Locus == NULL) return;
    if (Locus->LocusName != NULL) free(Locus->LocusName);
}


void free_all_from_locus_top(locus_top *LTop)
{
    int locus;

    if (LTop == NULL) return;
    if (LTop->Locus != NULL) {
        for (locus = 0; locus < LTop->LocusCnt; locus++)
            free_all_from_locus_rec(&(LTop->Locus[locus]));
        free(LTop->Locus);  LTop->Locus = NULL;
    }
}


static void free_all_from_ped_top(ped_top *Top, void (*PTmpFreeFun)(void *TmpData),
                                  void (*ETmpFreeFun)(void *TmpData))


{
    int ped;
    if (Top == NULL) return;
    if (Top->PedTree != NULL) {
        for (ped = 0; ped < Top->PedCnt; ped++)
            free_all_from_ped_tree(&(Top->PedTree[ped]), PTmpFreeFun, ETmpFreeFun);
        free(Top->PedTree);  Top->PedTree = NULL;
    }
    if (Top->LocusTop != NULL) {
        free_all_from_locus_top(Top->LocusTop);
        free(Top->LocusTop);  Top->LocusTop = NULL;
    }
}


void free_all_including_ped_top(ped_top *Top,  void (*PTmpFreeFun)(void *TmpData),
				void (*ETmpFreeFun)(void *TmpData))

{
    if (Top == NULL) return;
    free_all_from_ped_top(Top, PTmpFreeFun, ETmpFreeFun);
    free(Top);
}



void add_to_offspring(ped_rec *Parent, ped_rec *Child)

{
    if ((Parent == NULL) || (Child == NULL)) return;
    if (Parent->OffspringCnt == 0) {
        Parent->OffspringCnt = 1;
        Parent->Offspring = MALLOC(ped_rec *);
        Parent->Offspring[0] = Child;
    } else {
        Parent->OffspringCnt++;
        Parent->Offspring = REALLOC(Parent->Offspring, (size_t) Parent->OffspringCnt, ped_rec *);
        Parent->Offspring[Parent->OffspringCnt - 1] = Child;
    }
}


/* swap_loci()
 *
 * Slave routine for set_locus_order(). Swaps data in locus records.
 */
static void swap_lociDEFUNCT(locus_rec *locus1, locus_rec *locus2)

{
    register locus_rec tmplocus;

    tmplocus = *locus1;
    *locus1 = *locus2;
    *locus2 = tmplocus;
}


/* set_locus_order()
 *
 * Given the pedigree top and an array of integers of ID's in
 * which the loci should be ordered, order them.
 *
 * Ie.: If order[] = { 1, 3, 2 }, order the loci so that
 * the first remains first, the third is moved into the
 * second place, and the second one is moved to the end.
 *
 * Only existing loci may be in order[], and they may each
 * only be in order[] once. Not all need to be present,
 * however; those that are present will be first in the
 * Locus[] and Allele_*[] arrays, with those that were not
 * present moved to the end and put into some random order.
 *
 * order[] must have as many entries as there are loci,
 * but if not all loci are present then the extra entries
 * in order[] must be set to 0.
 *
 * Example: You have 3 loci. You wish to save the first and
 * the third, in that order, but you don't care about the second.
 * In that case, order[] should be { 1, 3, 0 }.
 *
 * Note that not only must the order of the Locus[] be changed,
 * but also the Allele_1[] and Allele_2[] of all members of
 * all pedigrees. Yuck!!!
 *
 * Return non-zero on error.
 */

int set_locus_orderDEFUNCT(ped_top *Top, int *order)
{
    int i, j;
    register int tmpi;
    register locus_top *LTop = Top->LocusTop;
    register int tmpallelep;
    register ped_rec *Entry;
    int ped, entry;
    int *place = CALLOC((size_t) LTop->LocusCnt, int);

    /* set up place to initial order which is always... */
    for (i = 0; i < LTop->LocusCnt; i++)
        place[i] = i+1;

    /* There's no really good way to do this. */
    for (i = 0; (i < LTop->LocusCnt) && (order[i] != 0); i++) {
        if (place[i] != order[i]) {
            if (i == LTop->LocusCnt - 1) {
                /* Last locus, and they don't agree? */
                fprintf(stderr, "ERROR: Internal error setting locus order!\n");
                free(place);
                return -1;
            } else {
                for (j = i+1; (j < LTop->LocusCnt) && (place[j] != order[i]); j++);
                swap_lociDEFUNCT(&(LTop->Locus[i]), &(LTop->Locus[j]));
                /* Now for all pedigree members... */
                for (ped = 0; ped < Top->PedCnt; ped++)
                    for (entry = 0; entry < Top->PedTree[ped].EntryCnt; entry++) {
                        Entry = &(Top->PedTree[ped].Entry[entry]);
                        // I've never seen Allele_1/Allele_2 initialized to a value
                        // but then, of course, we're never called.
                        tmpallelep = Entry->Allele_1[i];
                        Entry->Allele_1[i] = Entry->Allele_1[j];
                        Entry->Allele_1[j] = tmpallelep;
                        tmpallelep = Entry->Allele_2[i];
                        Entry->Allele_2[i] = Entry->Allele_2[j];
                        Entry->Allele_2[j] = tmpallelep;
                    }
                /* don't forget to change place */
                tmpi = place[i];
                place[i] = place[j];
                place[j] = tmpi;
            }
        }
    }
    free(place);
    return 0;
}


/* sort_affecteds()
 *
 * Sort the affecteds so that their IDs are in increasing order.
 *
 * This is a dumb n^2 routine. Maybe later we will want to
 * write something better.
 *
 * Return -1 on error, 0 on no error.
 */

static int sort_affecteds(ped_tree *Ped)
{
    int aff, aff2, lowest_ID;
    ped_rec **lowest_aff = 0; // will be set unless lowest_ID remains 0; but 0 does not allow ref of lowest_aff
    register ped_rec *tmp_aff;

    for (aff = 0; aff < Ped->AffectedCnt; aff++) {
        lowest_ID = 0;
        for (aff2 = aff; aff2 < Ped->AffectedCnt; aff2++) {
            if ((lowest_ID == 0)
                || (lowest_ID > Ped->Affected[aff2]->ID)) {
                lowest_ID = Ped->Affected[aff2]->ID;
                lowest_aff = &(Ped->Affected[aff2]);
            }
        }
        if (lowest_ID == 0) {
            fprintf(stderr, "ERROR: Internal error sorting affecteds\n");
            return -1;
        } else {
            tmp_aff = *lowest_aff;
            *lowest_aff = Ped->Affected[aff];
            Ped->Affected[aff] = tmp_aff;
        }
    }
    return 0;
}


/* renumber_ped()
 *
 * Rearrange the pedigree so that all members have
 * ID's greater than their parents'.
 *
 * We want to allow for the possibility that the entry
 * ID's are bad.
 *
 * Return -1 on error, 0 on no error.
 */
int renumber_ped(ped_tree *Ped)
{
    int entry, child, flag;
    register ped_rec *Entry, from_rec, to_rec;
    int entry_ID, to_index;
    int *place = CALLOC((size_t) Ped->EntryCnt, int);
    int *mother = CALLOC((size_t) Ped->EntryCnt, int);
    int *father = CALLOC((size_t) Ped->EntryCnt, int);
    int *affected = (Ped->AffectedCnt == 0 ? NULL : CALLOC((size_t) Ped->AffectedCnt, int));

    /* place[] maintains a map of old IDs to new IDs */
    /* initialize place[] */
    for (entry = 0; entry < Ped->EntryCnt; entry++)
        {
            place[entry] = 0;
            Entry = &(Ped->Entry[entry]);
            mother[entry] = (Entry->Mother == NULL ? 0 : Entry->Mother->ID);
            father[entry] = (Entry->Father == NULL ? 0 : Entry->Father->ID);
        }
    for (entry = 0; entry < Ped->AffectedCnt; entry++)
        {
            affected[entry] = Ped->Affected[entry]->ID;
        }

    /* Fix parent, child, and affected pointers so that they are ID's */
/*  for (entry = 0; entry < Ped->EntryCnt; entry++) {
    Entry = &(Ped->Entry[entry]);
    if (Entry->Father != NULL)
    Entry->Father = (ped_rec *) Entry->Father->ID;
    else
    Entry->Father = (ped_rec *) 0;
    if (Entry->Mother != NULL)
    Entry->Mother = (ped_rec *) Entry->Mother->ID;
    else
    Entry->Mother = (ped_rec *) 0;
    for (child = 0; child < Entry->OffspringCnt; child++)
    Entry->Offspring[child] = (ped_rec *) Entry->Offspring[child]->ID;
    }
    for (entry = 0; entry < Ped->AffectedCnt; entry++)
    Ped->Affected[entry] = (ped_rec *) Ped->Affected[entry]->ID; */

    /* Now set new IDs until we get a permissible, complete pedigree */
    /* This is the algorithm Dr. Weeks suggested. */
    entry_ID = 1;
    do {
        flag = 0;
        for (entry = 0; entry < Ped->EntryCnt; entry++) {
            Entry = &(Ped->Entry[entry]);
            if (place[entry] != 0)
                continue;
            if (((mother[entry] == 0)
                 || (place[mother[entry] - 1] != 0))
                && ((father[entry] == 0)
                    || (place[father[entry] - 1] != 0))) {
                place[entry] = entry_ID;
                entry_ID++;
                flag = 1;
            }
        }
    } while (flag);

    /* Swap the actual pedigree members around according to their new places */
    flag = 1;
    do {
        /* Find the first one out of place */
        to_index = 0;
        for (entry = 0; (entry < Ped->EntryCnt) && (to_index == 0); entry++) {
            if (place[entry] > 0) {
                if (entry + 1 != place[entry]) {
                    to_index = place[entry];
                    from_rec = Ped->Entry[entry];
                }
                place[entry] = -place[entry]; /* show that we've been here */
            }
        }
        /* entry contains the old ID of the first */
        /* to_index contains its new ID           */
        /* from_rec is a copy of it               */
        if (to_index > 0) {
            /* Worm through the array */
            entry--; /* make the ID an index */
            for (;;) {
                to_index--; /* make it a real index */
                to_rec = Ped->Entry[to_index]; /* copy before overwriting */
                Ped->Entry[to_index] = from_rec; /* overwrite */
                Ped->Entry[to_index].ID = to_index + 1; /* fix the ID */
                if (to_index == entry)
                    break; /* stop when we overwrite the first */
                from_rec = to_rec; /* copy this one to 'from' space */
                place[to_index] = -place[to_index]; /* show that we've been here */
                to_index = -place[to_index]; /* find where to put this one */
            }
        } else
            flag = 0;
    } while (flag);

    /* Make all place[] entries positive again */
    for (entry = 0; entry < Ped->EntryCnt; entry++)
        if (place[entry] < 0) place[entry] = -place[entry];

    /* Now we just need to fix the parent, child, and affected pointers */
    for (entry = 0; entry < Ped->EntryCnt; entry++) {
        Entry = &(Ped->Entry[entry]);
        if ( Entry->Father != NULL)
            Entry->Father = &(Ped->Entry[place[(Entry->Father->ID) - 1] - 1]);
        else
            Entry->Father = NULL;
        if (Entry->Mother != NULL)
            Entry->Mother = &(Ped->Entry[place[(Entry->Mother->ID) - 1] - 1]);
        else
            Entry->Mother = NULL;
        for (child = 0; child < Entry->OffspringCnt; child++) {
            Entry->Offspring[child]
                = &(Ped->Entry[place[(Entry->Offspring[child]->ID) - 1] - 1]);
        }
    }


    for (entry = 0; entry < Ped->AffectedCnt; entry++)
        Ped->Affected[entry] = &(Ped->Entry[place[affected[entry] - 1] - 1]);


    free(place);
    free(mother);
    free(father);
    free(affected);

    /* we probably have to sort the affecteds */
    return sort_affecteds(Ped);

}



/* is_typed()
 *
 * Return 1 if the entry is typed, 0 if not.
 *
 * Whether it's typed depends on whether mode
 * is TYPED_AT_ALL_LOCI or TYPED_AT_ANY_LOCUS.
 */

int is_typed(ped_rec *Entry, locus_top *LTop1, int mode)
{

    int typed;

    typed = is_typed_lentry(*(Entry->LEntry), LTop1->LTop, mode, NULL);

/*   for (locus1 = 0; locus1 < LTop1->LocusCnt; locus1++) { */

/*     if ((Entry->Allele_1[locus1] != 0) && (Entry->Allele_2[locus1] != 0)) */
/*       typed++; */
/*   } */

/*   if (mode & TYPED_AT_ANY_LOCUS) { */
/*     if (typed > 0) return 1; */
/*     else return 0; */
/*   } */
/*   else { */
/*     if (typed == LTop1->LocusCnt) return 1; */
/*     else */
/*       return 0; */
/*   } */
    return typed;
}

static int is_typed_order(ped_rec *Entry,
                          int mode, int *order, int num_loci)
{
    int locus1;
    int typed, i;

    typed = 0;
    for(i=0; i<num_loci; i++) {
        locus1=order[i];
        if ((Entry->Allele_1[locus1] >0) && (Entry->Allele_2[locus1] > 0)) {
            typed++;
        }
    }

    if (mode==TYPED_AT_ANY_LOCUS) {
        if (typed > 0)  return 1;
        else return 0;
    } else {
        if (typed == num_loci) return 1;
        else return 0;
    }
}


/* remove_untyped_affecteds()
 *
 * Remove from the Affected[] array pointers to members that
 * are untyped at all loci.
 *
 * This is easily done by scanning the array, keeping track
 * of which affecteds are typed and how many there are, then
 * creating a new array, copying the typed ones over, freeing
 * the old array, assigning the new array, and updating AffectedCnt.
 *
 * Mode can be TYPED_AT_ALL_LOCI or TYPED_AT_ANY_LOCUS. If it
 * is the former, affecteds must be typed at all loci to be saved.
 * If it is the latter, affecteds must only be typed at any locus.
 * (The former is the default.)
 */

void remove_untyped_affecteds(ped_tree *Ped, locus_top *LTop, int mode,
			      int *order, int num_loci)
{
    int aff, affcnt;
    int newaff;
    int *typed;
    ped_rec **NewAffected = 0;  // stupid msvc: if affcnt > 0 NewAFfected is assigned ... and used

    /* Record which ones are typed and their number. */
    affcnt = 0;
/* typed=CALLOC((size_t) Ped->AffectedCnt, int); */
    typed = (Ped->AffectedCnt == 0 ? NULL : CALLOC((size_t) Ped->AffectedCnt, int));

    for (aff = 0; aff < Ped->AffectedCnt; aff++) {
        if (order != NULL) {
            typed[aff]=is_typed_order(Ped->Affected[aff], mode, order, num_loci);
        } else {
            typed[aff] = is_typed(Ped->Affected[aff], LTop, mode);
        }
        if (typed[aff])  {
            affcnt++;
        }
    }

    if (affcnt > 0) {
        /* Make the new array. */
        NewAffected = CALLOC((size_t) affcnt, ped_rec *);

        /* Copy the affecteds into the new array. */
        newaff = 0;
        for (aff = 0; aff < Ped->AffectedCnt; aff++) {
            if (typed[aff]) {
                NewAffected[newaff] = Ped->Affected[aff];
                newaff++;
            }
        }

    }

    /* Free the old array. */
    free(Ped->Affected);
    Ped->Affected = NULL;

    if (affcnt > 0) {
        /* Assign the new array and change AffectedCnt. */
        Ped->Affected = NewAffected;
        Ped->AffectedCnt = affcnt;
    } else {
        Ped->AffectedCnt = 0;
    }

    free(typed);

    return;
}


/*---------------------------------------------------------------+
  | Convert linkage structures to pedtree structures.             |
  | This routine presumes that the affected members will be       |
  | marked by their TmpData set to (int) 1.                       |
  |                                                               |
  | It also reconnects loops.                                     |
  |                                                               |
  | Return the new pedtree top.                                   |
  +---------------------------------------------------------------*/

ped_top        *convert_to_pedtree(linkage_ped_top *Top,
				   int assign_affs)

{
    int             ped, entry1, entry, locus1;
    int             llocus1, mito_locus, l1;
    ped_top        *PTop;
    register ped_tree *PPed;
    register linkage_ped_tree *LPed;
    register ped_rec *PRec;
    register linkage_ped_rec *LRec;
    locus_top      *PLTop;
    linkage_locus_top *LLTop;
    /*  int plocus, LocusMapSize, *LocusMap=NULL; */

#ifdef DEBUG
    printf("\nConverting to pedigree tree ...\n");
#endif

    LLTop = Top->LocusTop;
    PTop = new_ped_top();

    /* Create the locus map. */
    /*   LocusMapSize = num_reordered; */
    /*   LocusMap = CALLOC((size_t) LocusMapSize, int); */

    /* Convert the locus data first. */
    /* Store the mitochondrial marker numbers in a separate list */
    PedTreeMito = NULL;
    for (l1=0; l1 < main_chromocnt; l1++) {
        if (global_chromo_entries[l1] == MITO_CHROMOSOME) {
            if (l1 == 0) {
                PedTreeMito = CALLOC((size_t) chromo_loci_count[l1]+1, int);
            } else {
                PedTreeMito = CALLOC((size_t) chromo_loci_count[l1] - chromo_loci_count[l1-1] + 2, int);
            }
            break;
        }
    }

    locus1=0;
    for (l1 = 0; l1 < num_reordered; l1++) {
        if (LLTop->Locus[reordered_marker_loci[l1]].Class != TRAIT) {
            locus1++;
        }
    }

    PLTop = (PTop->LocusTop = new_locus_top());

    /* Count the number of loci with specified allelic genotypes. */

    /* Now convert it. */
    PLTop->LocusCnt = locus1;
    PLTop->SexLinked = LLTop->SexLinked;

    if (PLTop->LocusCnt > 0) {
        PLTop->Locus = CALLOC((size_t) PLTop->LocusCnt, locus_rec);
        /*    plocus = 0; */
        locus1=0; mito_locus = 0;
        for (l1 = 0; l1 < num_reordered; l1++) {
            if (LLTop->Locus[reordered_marker_loci[l1]].Class == TRAIT) {
                continue;
            }
            llocus1 = reordered_marker_loci[l1];
            /* Set the name. */
            PLTop->Locus[locus1].LocusName = strdup(LLTop->Locus[llocus1].LocusName);
            PLTop->Locus[locus1].chromosome = LLTop->Marker[llocus1].chromosome;
            /* Allele records are identical. */
            PLTop->Locus[locus1].AlleleCnt = LLTop->Locus[llocus1].AlleleCnt;
            PLTop->Locus[locus1].Allele = (allele_rec *) LLTop->Locus[locus1].Allele;
            PLTop->Locus[locus1].linkage_loc_num = llocus1;
            PLTop->Locus[locus1].linkage_loc_rec = &(LLTop->Locus[llocus1]);
            if (LLTop->Marker[llocus1].chromosome == MITO_CHROMOSOME) {
                PedTreeMito[mito_locus] = locus1;
                mito_locus++;
            }
            locus1++;
        }
        if (PedTreeMito != NULL) {
            PedTreeMito[mito_locus] = -99;
        }
    } else {
        PLTop->Locus=NULL;
    }
    PLTop->LTop = LLTop;

    PTop->PedCnt = Top->PedCnt;
    PTop->PedTree = CALLOC((size_t) PTop->PedCnt, ped_tree);
    /*  save_linkage_peds("TMP1", Top, TO_LINKAGE, 0, NULL); */
    for (ped = 0; ped < PTop->PedCnt; ped++) {
        LPed = &(Top->Ped[ped]);
        PPed = &(PTop->PedTree[ped]);
        clear_ped_tree(PPed);
        PPed->Name = strdup(LPed->Name);
        PPed->PedPre = strdup(LPed->PedPre);
        PPed->EntryCnt = LPed->EntryCnt;

        PPed->Entry = CALLOC((size_t) PPed->EntryCnt, ped_rec);

        /* clear the entries */
        for (entry = 0; entry < PPed->EntryCnt; entry++) {
            clear_ped_rec(&(PPed->Entry[entry]));
            PRec = &(PPed->Entry[entry]);
            LRec = &(LPed->Entry[entry]);

            PRec->Orig_status = LRec->Orig_status;
            PRec->LEntry = LRec;
            PRec->Marker = LRec->Marker;
            PRec->ID = entry+1;
            PRec->Sex = LRec->Sex;
            PRec->connected = 0;
            PRec->sibship_checked = 0;

            /* Initialize the allele arrays. */
            if (FirstTime) {
                /* Set the IsTyped field of both PedTree and LPedTree */
                LRec->IsTyped=0; PRec->IsTyped = 0;
                locus1=0;
                for (l1=0; l1 < num_reordered; l1++) {
                    int a1, a2;
                    if (LLTop->Locus[reordered_marker_loci[l1]].Class == TRAIT) {
                        continue;
                    }
                    get_2alleles(LRec->Marker, reordered_marker_loci[l1], &a1, &a2);
                    if (a1 > 0 && a2 > 0 ) {
                        PRec->IsTyped++;
                    }
                    locus1++;
                }
                if (PRec->IsTyped == locus1) {
                    LRec->IsTyped = 2;
                } else if (PRec->IsTyped > 0) {
                    LRec->IsTyped = 1;
                } else {
                    LRec->IsTyped =0;
                }
            }
        }
        /* Settle parents, siblings and proband. */
        for (entry = 0; entry < PPed->EntryCnt; entry++)  {
            LRec = &(LPed->Entry[entry]);
            PRec = &(PPed->Entry[entry]);
/*       printf("%d %d %d\n", LRec->ID, LRec->Next_PA_Sib, LRec->Next_MA_Sib); */

            if (LRec->Mother) {
                for (entry1=0; entry1 < LPed->EntryCnt; entry1++) {
                    if (LPed->Entry[entry1].ID == LRec->Mother) {
                        PRec->Mother = &(PPed->Entry[entry1]);
                        add_to_offspring(PRec->Mother, PRec);
                        break;
                    }
                }
            }

            if (LRec->Father) {
                for (entry1=0; entry1 < LPed->EntryCnt; entry1++) {
                    if (LPed->Entry[entry1].ID == LRec->Father) {
                        PRec->Father = &(PPed->Entry[entry1]);
                        add_to_offspring(PRec->Father, PRec);
                        break;
                    }
                }
            }

            if (LRec->Next_PA_Sib){
                for (entry1=0; entry1 < LPed->EntryCnt; entry1++) {
                    if (LPed->Entry[entry1].ID == LRec->Next_PA_Sib) {
                        PRec->PA_Sib = &(PPed->Entry[entry1]);
                        break;
                    }
                }
            }

            if (LRec->Next_MA_Sib){
                for (entry1=0; entry1 < LPed->EntryCnt; entry1++) {
                    if (LPed->Entry[entry1].ID == LRec->Next_MA_Sib) {
                        PRec->MA_Sib = &(PPed->Entry[entry1]);
                        break;
                    }
                }
            }
        }

        /* Now go through and resolve offspring. */
        /* Now find the proband. */
        if (LPed->Proband > 0) {
            for (entry1=0; entry1 < PPed->EntryCnt; entry1++) {
                if (LPed->Proband == LPed->Entry[entry1].ID) {
                    PPed->Proband = &(PPed->Entry[entry1]);
                    break;
                }
            }
        }
        else
            PPed->Proband = NULL;

    }

    /* if we have already set the affected status */

    if (assign_affs) {
        reassign_affecteds(PTop);
    }
    return PTop;
}

/*---------------------------------------------------*/

/* reassign_affecteds: after a call to affected_by_status with
   a specific disease locus we have to re-compute the affected status
   of entries in the pedtree, and set the affecteds correctly
*/

void reassign_affecteds(ped_top *PTop)

{
    int ped, entry;
    int affnum;
    linkage_ped_rec *lentry;

    for (ped = 0; ped < PTop->PedCnt; ped++) {
        PTop->PedTree[ped].AffectedCnt = 0;
        for (entry = 0; entry < PTop->PedTree[ped].EntryCnt; entry++) {
            lentry = PTop->PedTree[ped].Entry[entry].LEntry;
            PTop->PedTree[ped].Entry[entry].Orig_status = lentry->Orig_status;
            /*	LTop->Ped[ped].Entry[entry].Orig_status; */

            if (*(int *) (lentry->TmpData) == 1) {
                /* This person is affected. */
                PTop->PedTree[ped].AffectedCnt++;
                PTop->PedTree[ped].Entry[entry].Status = AFFECTED_STATUS;
            } else
                PTop->PedTree[ped].Entry[entry].Status = UNAFFECTED_STATUS;
        }
        if (PTop->PedTree[ped].Affected != NULL) {
            free(PTop->PedTree[ped].Affected);
        }
        /*   PTop->PedTree[ped].Affected =
             CALLOC((size_t) PTop->PedTree[ped].AffectedCnt,
             ped_rec *); */
        PTop->PedTree[ped].Affected
            = (PTop->PedTree[ped].AffectedCnt == 0 ? NULL
               : CALLOC((size_t) PTop->PedTree[ped].AffectedCnt, ped_rec *));

        /* Also resolve affecteds into array.    */
        affnum = 0;
        for (entry = 0; entry < PTop->PedTree[ped].EntryCnt; entry++) {
            if (PTop->PedTree[ped].Entry[entry].Status == AFFECTED_STATUS) {
                PTop->PedTree[ped].Affected[affnum++]
                    = &(PTop->PedTree[ped].Entry[entry]);
            }
        }
    }

    return;
}
