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

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"
#include "tod.hh"

#include "utils_ext.h"
#include "errorno.h"

#include "error_messages_ext.h"
#include "grow_string_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
     error_messages_ext.h:  errorf my_calloc warnf
        grow_string_ext.h:  grow
*/


/* Made major changes to the check_ped routine, and
   check_inheritance routine NM - 3/1/00 */

/* exported functions */
int check_locus(linkage_locus_rec *Locus,
		analysis_type analysis,
		int *plink_locus_num);
void clear_ped_status(ped_status *Stat);

void mito_transmission_report(ped_top *PTop,
			      int *num_hetero, int *num_non_maternal);

/*-----------------------*/

/* static functions */
/*  static int check_ordinal(ped_tree *PedTree, ped_status *PedStatus); */

static int check_inheritance(ped_rec *PedEntry,
			     linkage_locus_top *LTop, const int locus,
                             int A1, int A2, const char *name,
                             const int uniqueids);

static int check_ped_connected(ped_tree *Ped);
static void connected_traverse(ped_rec *Entry, int curr_ped);
static int check_sibship_alleles(ped_tree *PedTree, ped_rec **Sibs, int locus,
				 linkage_locus_top *LTop1, int *sibship_checked,
				 int uniqueids);
static int autosomal_or_female_xlinked(int all1, int all2,
				       int *alleles, int *allelecnt,
				       int maxalleles);
static int male_xlinked(int all1, int *allele, int *allelecnt);
/*----------------*/
int cmp_ids(const void *aa, const void *bb);

void clear_ped_status(ped_status *Stat)
{
    Stat->father_invalid = 0;
    Stat->mother_invalid = 0;
    Stat->parents_not_defined=0;
    Stat->incomplete_entry = 0;
    Stat->offspring_out_of_range = 0;

    Stat->genotype_invalid = -1;
    Stat->halftyped = -1;
    Stat->exceed_allcnt = -1;
    Stat->entry_unconnected = -1;
    Stat->checked_sibship = -1;

    Stat->bad_ID = 0;
    Stat->bad_sex = 0;
    Stat->unknown = 0;
    Stat->invalid_genotypes=NULL;
    Stat->half_types=NULL;
    Stat->allele_outof_bounds = NULL;
}

/* check_ordinal()
 * check that the entry IDs are ordinal
 * if they are not ordinal, warn the user and stop
 */

/* static int check_ordinal(ped_tree *PedTree, ped_status *PedStatus) */


/* check_inheritance()
 *
 * Check the alleles A1 and A2 against the parents'
 * PA* are paternal, MA* are maternal.
 *
 * Return non-zero if the combination A1 / A2 isn't
 * possible given the parents' genotypes.???   non-zero means it is a food inheritance, lwa

 modfied 4/8/97, lwa, to include checks on Xlinked loci */


/* here locus1 is the linkage-locus-top index and ploc is the
   pedtree->locustop index
*/

/* Ignore females' genotypes for y-linked loci */

SECTION_ERR_EXTERN(check_inheritance);
static int check_inheritance(ped_rec *PedEntry,
			     linkage_locus_top *LTop, const int locus,
                             int A1, int A2, const char *name,
			     const int uniqueids)
{
    int PA[2], MA[2];
    ped_rec *Father=PedEntry->Father;
    ped_rec *Mother=PedEntry->Mother;
    /* Changed the definition of sex-linked to specifically x-linked,
       this only considers loci on X chromo to be x-linked. */
    int sex_linked = ((LTop->Marker[locus].chromosome == SEX_CHROMOSOME
                       && LTop->SexLinked >= 1) ? 1 : 0);

    int y_linked = ((LTop->Marker[locus].chromosome == MALE_CHROMOSOME
		     && LTop->SexLinked >= 1)? 1 : 0);

    /* Note: This check does not require that Father != NULL */
    // Display_Errors is set by the caller...
    if ((sex_linked  || y_linked) && IS_MALE(*PedEntry)) {
        /* if child is male,  must be a homozygote */
        if (!R(A1,A2)) {
            SECTION_ERR(check_inheritance);
            warnvf("Ped %s: Male %s is a heterozygote at %c-linked locus %s.\n",
		   name,
		   ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
		   ((y_linked)? 'Y' : 'X'),
		   LTop->Locus[locus].LocusName);
            return 0;
        }
    }

    if (Father == NULL) { return 1; }
    /* PTALLELE1(locus) = father->lentry->Data[locus].Allele_1 etc. */

    get_2alleles(Father->Marker, locus, &PA[0], &PA[1]);
    get_2alleles(Mother->Marker, locus, &MA[0], &MA[1]);

    /* Sex linked allele inheritance check */
    if (!(sex_linked || y_linked) || (sex_linked && IS_FEMALE(*PedEntry))) {
        /* sexLinked = zero, if autosomal,
           female child in either case */

        if ((R(A1, PA[0]) || R(A1, PA[1]))  &&
            (R(A2, MA[0]) || R(A2, MA[1]))) {
            return 1;
        }
        if ((R(A1, MA[0]) || R(A1, MA[1]))  &&
            (R(A2, PA[0]) || R(A2, PA[1]))) {
            return 1;
        }

        SECTION_ERR(check_inheritance);
        warnvf("Ped %s: Entry %s has non-mendelian genotype at locus %s (chr %d).\n",
	       name,
	       ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
	       LTop->Locus[locus].LocusName,
	       LTop->Marker[locus].chromosome);
        warnvf("\t%s <- %s X %s\n",
	       ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
	       ((uniqueids == 1)? Father->LEntry->UniqueID : Father->LEntry->OrigID),
	       ((uniqueids == 1)? Mother->LEntry->UniqueID : Mother->LEntry->OrigID));
        // This will write the numeric allele always...
        warnvf("\t%d/%d <- %d/%d X %d/%d\n",
	       A1, A2, PA[0], PA[1], MA[0], MA[1]);

        /* Father->LEntry->OrigID, Mother->LEntry->OrigID, PedEntry->LEntry->OrigID); */

        return 0;
    }

    /* These two checks are for hemizygotes */
    HEMI2HOMO(A1, A2);

    // Display_Errors is set by the caller...
    if (sex_linked && IS_MALE(*PedEntry)){
        /* if child is male,  both alleles must match the mother's X-chromosome chromosome */
        /* We have already checked if this person is a homozygote */
        /* This returns 1 if any allele is untyped */
        if (R(A1,MA[0]) || R(A1, MA[1])) {
            return 1;
        }

        SECTION_ERR(check_inheritance);
        warnvf("Ped %s: Male %s has non-mendelian genotype at X-linked locus %s\n",
	       name,
	       ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
	       LTop->Locus[locus].LocusName);
        warnvf("\t%s <- %s X %s\n",
	       ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
	       ((uniqueids == 1)? Father->LEntry->UniqueID : Father->LEntry->OrigID),
	       ((uniqueids == 1)? Mother->LEntry->UniqueID : Mother->LEntry->OrigID));
        // This will write the numeric allele always...
        warnvf("\t%d/ <- %d/ X %d/%d\n",
	       A1, PA[0], MA[0], MA[1]);

        return 0;
    }

    /* y-linked check for males */
    if (y_linked && IS_MALE(*PedEntry)){
        /* if child is male,  must match father's genotype */
        if (R(A1, PA[0]) || R(A1, PA[1])) {
            /* This takes care of father's missing allele */
            return 1;
        }

        SECTION_ERR(check_inheritance);
        warnvf("Ped %s: Male %s has non-mendelian genotype at locus %s\n",
	       name, PedEntry->LEntry->OrigID,
	       LTop->Locus[locus].LocusName);
        warnvf("\t%s <- %s\n",
	       PedEntry->LEntry->OrigID,
	       Father->LEntry->OrigID);
        // This will write the numeric allele always...
        warnvf("%d/ <- %d/\n", A1, PA[0]);

        return 0;
    }

    /* put in this bit of unreachable code and see if the compiler
       complains */
    return 1;
}

/* check the PA_sibs and MA_sibs */
/* locus is the index of locus being checked in locus_top */
SECTION_ERR_EXTERN(check_sibship_alleles);
static int check_sibship_alleles(ped_tree *PedTree, ped_rec **Sibs,
                                 const int locus,
				 linkage_locus_top *LTop, int *checked_sibship,
				 const int uniqueids)

{

    int alleles[4], all1, all2;
    int entry, off, haserr, inc, allelecnt, a;
    int max_allelecnt, sib_cnt;
//  ped_rec *Entry, *Sib, **Sibs;
    ped_rec *Entry, *Sib;
    int sex_linked = (((LTop->SexLinked == 2 &&
                        LTop->Marker[locus].chromosome == SEX_CHROMOSOME) ||
                       (LTop->SexLinked == 1)) ? 1 : 0);
    int y_linked = (((LTop->SexLinked == 2 &&
                      LTop->Marker[locus].chromosome == MALE_CHROMOSOME)) ? 1 : 0);

    max_allelecnt = ((sex_linked)? 3 : 4);

//  Sibs = CALLOC((size_t) PedTree->EntryCnt, ped_rec *);
    haserr = 0;
    /* set all the checked flags to 0 */
    for(entry=0; entry < PedTree->EntryCnt; entry++) {
        PedTree->Entry[entry].sibship_checked=0;
    }

    /* Now check the full-sibs of each entry that has not
       been checked */
    // Display_Errors is set by the caller...
    for(entry=0; entry < PedTree->EntryCnt; entry++) {
        Entry = &(PedTree->Entry[entry]);

        int P1, P2, M1, M2;
        if (Entry->Father == NULL) {
            Entry->sibship_checked=1;
            continue;
        }
        get_2alleles(Entry->Father->Marker, locus, &P1, &P2);
        get_2alleles(Entry->Mother->Marker, locus, &M1, &M2);
        if ((!Entry->sibship_checked) &&
            (P1 == 0 || P2 == 0 || M1 == 0 || M2 == 0)) {

            sib_cnt=0;   inc=0;

            Sib=Entry;	/* Check the full-sibs */
            while(Sib != NULL) {
                Sibs[sib_cnt] = Sib;
                sib_cnt++;
                if (Sib->PA_Sib != NULL && Sib->MA_Sib != NULL &&
                   Sib->PA_Sib == Sib->MA_Sib) {
                    Sib = Sib->PA_Sib;
                } else {
                    Sib = NULL;
                }
            }

            if (sib_cnt > 1) {
                /* check only if there are more than 2 offspring */
                *checked_sibship=1; /* set the checked sibship flag */
                /* first set all 4 alleles to 0, and allelecnt to 0 */
                allelecnt=0;
                for(a=0; a<4; a++) { alleles[a]=0; }

                for(off=0; off < sib_cnt; off++) {
                    Sib=Sibs[off];
                    get_2alleles(Sib->Marker, locus, &all1, &all2);
                    if (all1 != 0 && all2 != 0) {
                        inc = ((!sex_linked || Sib->Sex == 2)?
                               autosomal_or_female_xlinked(all1, all2, alleles,
                                                           &allelecnt, max_allelecnt) :
                               male_xlinked(all1, alleles, &allelecnt));
                        if (inc) break;
                    }
                }
                if (inc) {
                    SECTION_ERR(check_sibship_alleles);
                    sprintf(err_msg,
                            "Ped %s: Too many distinct alleles in sibship at %slocus %s:",
                            PedTree->Name,
                            ((sex_linked)? "x-linked " : ""),
                            LTop->Locus[locus].LocusName);
                    errorf(err_msg);

                    /* only one sibling needs to be marked as inconsistent for
                       all other siblings to be inconsistent as well */
                    for(off=0; off<sib_cnt; off++) {
                        int a1, a2;
                        Sib=Sibs[off];
                        get_2alleles(Sib->Marker, locus, &a1, &a2);
                        if (!sex_linked) {
                            sprintf(err_msg,
                                    "     Ped %s, entry %s has genotype [%d/%d].", PedTree->Name,
                                    ((uniqueids == 1)? Sib->LEntry->UniqueID : Sib->LEntry->OrigID),
                                    a1, a2);
                        } else {
                            sprintf(err_msg,
                                    "     Ped %s, entry %s (%s) has genotype [%d/%d].", PedTree->Name,
                                    ((uniqueids == 1)? Sib->LEntry->UniqueID : Sib->LEntry->OrigID),
                                    ((Sib->Sex == 1)? "male" : "female"),
                                    a1, a2);
                        }
                        errorf(err_msg);
                    }
                    haserr=1;
                }

                /* added y-chromosome genotype check NM - Aug 2008 */
                if (y_linked) {
                    inc=0;
                    alleles[0]=0;
                    for (off=0; off < sib_cnt; off++) {
                        int a1, a2;
                        Sib=Sibs[off];
                        get_2alleles(Sib->Marker, locus, &a1, &a2);
                        if (Sib->Sex == 2) {
                            Sib->sibship_checked =1;
                            continue;
                        }
                        if (alleles[0] > 0) {
                            if (a1 == alleles[0]) continue;
                            else {
                                inc = 1;
                                Sib->sibship_checked =1;
                                break;
                            }
                        } else {
                            alleles[0] = a1;
                            Sib->sibship_checked = 1;
                        }
                    }
                    if (inc) {
                        SECTION_ERR(check_sibship_alleles);
                        sprintf(err_msg,
                                "Ped %s: Too many distinct Y-alleles in sibship at locus %s:",
                                PedTree->Name,
                                LTop->Locus[locus].LocusName);
                        errorf(err_msg);

                        for(off=0; off<sib_cnt; off++) {
                            int a1, a2;
                            if (Sibs[off]->Sex == 2) continue;

                            Sib=Sibs[off];
                            get_2alleles(Sib->Marker, locus, &a1, &a2);
                            sprintf(err_msg,
                                    "     Ped %s, entry %s has Y-allele %d.", PedTree->Name,
                                    Sib->LEntry->OrigID, a1);
                            errorf(err_msg);
                        }
                        haserr=1;
                    }
                }
            }
        }
        if (inc) break;
    }
//  free(Sibs);
    return haserr;
}


/* connected_traverse()
 *
 * A slave function to check_ped_connected() which recursively
 * traverses the pedigree updating the connected[] array.
 */
static void connected_traverse(ped_rec *Entry, const int current_ped)
{
    int child;

    if (Entry->Father != NULL) {
        if (!(Entry->Father->connected)) {
            Entry->Father->connected = current_ped;
            connected_traverse(Entry->Father, current_ped); // Tail recursion!
        }
    }
    if (Entry->Mother != NULL) {
        if (!(Entry->Mother->connected)) {
            Entry->Mother->connected = current_ped;
            connected_traverse(Entry->Mother, current_ped);
        }
    }

    for (child = 0; child < Entry->OffspringCnt; child++)
        if (!(Entry->Offspring[child]->connected)) {
            Entry->Offspring[child]->connected = current_ped;
            connected_traverse(Entry->Offspring[child], current_ped);
        }

    if (!(Entry->connected)) {
        Entry->connected = current_ped;
    }

    /*  printf("%d %d\n", Entry->LEntry->ID, Entry->connected); */
}


/* check_ped_connected()
 *
 * Traverse the tree and check if the pedigree is connected
 * (so that there is no member which cannot be reached
 * from any oher member by traversing the tree).
 *
 * Return the number of members found to be disconnected
 * from the first member.
 */
static int check_ped_connected(ped_tree *Ped)
{
    int entry, conn;
    int found_disconnect=1, current_ped=1;
    int root=0;
    char names[71*2];

    if (Ped->EntryCnt == 1) return 0;

    /* maximum number of unconnected peds */

    root=0;
    while (found_disconnect) {
        connected_traverse(&(Ped->Entry[root]), current_ped);
        found_disconnect=0;
        for (entry = 0; entry < Ped->EntryCnt; entry++) {
            if (!Ped->Entry[entry].connected) {
                current_ped++;
                found_disconnect=1;
                root=entry;
                break;
            }
        }
    }

    if (current_ped > 1) {
        errorvf("Ped %s: Found %d disconnected sub-pedigrees:\n",
                Ped->Name, current_ped);
        for (conn = 0; conn < current_ped; conn++) {
            sprintf(names,"       ");
            entry =0 ;
            // This could easily overflow 'names[]' since "char UniqueID[NAMELEN]"...
            while(strlen(names) < 70 && entry < Ped->EntryCnt) {
                if (Ped->Entry[entry].connected == (conn+1)) {
                    grow(names, "  %s",
                         Ped->Entry[entry].LEntry->UniqueID); //full name; keep only; post loop
                }
                entry++;
            }

            warnvf("Sub-pedigree %d including (but not limited to) individuals\n%s\n",
		   conn+1, names);
        }
    }

    return (current_ped - 1);
}

int cmp_ids(const void *aa, const void *bb)
{
    const uniq_id_type *a, *b;
    int val;

    a = (const uniq_id_type *) aa;
    b = (const uniq_id_type *) bb;

    val = strcmp(a->id, b->id);
    val = val ? val : a->ped - b->ped;
    return val;

}

int check_unique_ids(linkage_ped_top *Top)
{
    int i, ped, per, uniq;

    uniq_id_type *uniq_ids = CALLOC((size_t) Top->IndivCnt, uniq_id_type);
    int  nonuniq=0;
    int (*compare_ids)(const void *, const void *);
    int id_count = 1;

    uniq=0;
    compare_ids = cmp_ids;
    SECTION_ERR_INIT(check_unique_ids);
    for (ped=0; ped < Top->PedCnt; ped++) {
        for (per=0; per < Top->Ped[ped].EntryCnt; per++) {
            if ((Top->Ped[ped].Entry[per].loopbreakers != NULL) ||
                (Top->Ped[ped].Entry[per].OrigProband > 1)) {
/*
                mssgvf("ped/per %d/%d, ID: %s, loopbreakers %p, Proband %d\n",
                       (int) ped, (int) per,
                       Top->Ped[ped].Entry[per].UniqueID,
                       Top->Ped[ped].Entry[per].loopbreakers,
                       Top->Ped[ped].Entry[per].OrigProband);
*/
                continue;
            }
            /* Possible error: casting size_t to an int */
            uniq_ids[uniq].ped = (int) ped;
            uniq_ids[uniq].per = (int) per;
            uniq_ids[uniq].id = strdup(Top->Ped[ped].Entry[per].UniqueID);
            uniq++;
        }
    }

    /* This will group identical IDs together */
    qsort(((void *) uniq_ids), (size_t) uniq, sizeof(uniq_id_type), compare_ids);

    for (i=0; i < uniq - 1; i++) {
        if (!strcmp(uniq_ids[i].id, uniq_ids[i+1].id)) {
            id_count++;
        } else {
            if (id_count > 1) {
                SECTION_ERR(check_unique_ids);
                warnvf("\"ID\" value %s duplicated in %d pedigrees.\n", uniq_ids[i].id, id_count);
                nonuniq=1;
            }
            id_count = 1;
        }
    }
    if (!strcmp(uniq_ids[i-1].id, uniq_ids[i].id)) {
        if (id_count > 1) {
            SECTION_ERR(check_unique_ids);
            warnvf("\"ID\" value %s duplicated in %d pedigrees.\n", uniq_ids[i].id, id_count);
            nonuniq=1;
        }
    }
    SECTION_ERR_FINI(check_unique_ids);
    for (i=0; i < Top->IndivCnt; i++) {
        free(uniq_ids[i].id);
    }
    free(uniq_ids);
    return nonuniq;

}
/* check_ped_relations()
   only check the pedigree structure, IDs, parents etc. */

/**
   @brief Checks on the pedigree relationships

   @param PedTree    A pointer to a ped_tree structure (see pedtree.h)
   @param PedStatus  A vector of ped_status structures
   @return           0 if all OK, 2 on error
 */
int check_ped_relations(ped_tree *PedTree, ped_status *PedStatus)

{
    int entry;
    ped_rec *PedEntry;
    char *entryID;

    /* if a real file pointer is specified, then also print to stderr */
    /*   if (PedStatus == NULL) return 0; */

    /*  clear_ped_status(PedStatus); */
    if (PedTree->Name == NULL) {
        // In this particular case 'Name' is used in many places following this code
        // without a further check for if being == NULL. So, if allowed to continue and
        // an error occurred Mega2 would segfault.
        errorf("While checking pedigree relationships, a name was found to be undefined");
        EXIT(DATA_TYPE_ERROR);
    }

    if (PedTree->Proband == NULL) {
        warnvf("Pedigree %s: proband is undefined\n", PedTree->Name);
    }

    PedStatus->entry_unconnected += check_ped_connected(PedTree);

    for (entry = 0; entry < PedTree->EntryCnt; entry++) {
        PedEntry = &(PedTree->Entry[entry]);

        /* 1. check ID */
        if (PedEntry->ID == UNDEF) {
            warnvf("Pedigree %s: entry ID (%d) is undefined\n",
		   PedTree->Name, entry + 1);
            PedStatus->bad_ID++;
            sprintf(PedEntry->LEntry->OrigID, "%d", entry + 1); /* temporarily */
        }

        entryID = &(PedEntry->LEntry->OrigID[0]);

        /* 2, 3, 4, check parents */
        if (PedEntry->Father != NULL) {
            if (!IS_MALE((*(PedEntry->Father)))) {
                errorvf("Pedigree %s: father %s of %s is not male!\n",
                        PedTree->Name, PedEntry->Father->LEntry->OrigID, entryID);
                PedStatus->father_invalid++;
            }
        }
        if (PedEntry->Mother != NULL) {
            if (!IS_FEMALE((*(PedEntry->Mother)))) {
                errorvf("Pedigree %s: mother %s  of %s is not female!\n",
                        PedTree->Name, PedEntry->Mother->LEntry->OrigID, entryID);
                PedStatus->mother_invalid++;
            }
        }
        if (PedEntry->Mother == NULL && PedEntry->Father != NULL) {
            errorvf("Pedigree %s: Missing record for mother of %s\n", PedTree->Name, entryID);
            PedStatus->parents_not_defined++;
        } else if (PedEntry->Mother != NULL && PedEntry->Father == NULL) {
            errorvf("Pedigree %s: Missing record for father of %s\n", PedTree->Name, entryID);
            PedStatus->parents_not_defined++;
        }

        /* 6. check sex */
        if ((!IS_MALE(*(PedEntry))) && (!IS_FEMALE(*(PedEntry)))) {
            errorvf("Pedigree %s: member %s has unknown sex (%d)\n",
		    PedTree->Name, entryID, PedEntry->Sex);
            PedStatus->bad_sex++;
        }
    }

    // Later the code will EXIT when the return status == 2.
    if (PedStatus->bad_ID || PedStatus->parents_not_defined ||
        PedStatus->mother_invalid || PedStatus->father_invalid ||
        PedStatus->incomplete_entry || PedStatus->bad_sex ) {
        return 2;
    }

    return 0;
}

int check_out_of_bounds(ped_tree *PedTree, ped_status *PedStatus,
                        linkage_locus_top *LTop, int ped_num, int locus,
                        int uniqueids,
                        FILE **reset_fp, bool *first, int reset)
{
    int entry, iserr;
    ped_rec *PedEntry;
    int abortf = 0;
    int all1, all2;

    int AlleleCnt = LTop->Locus[locus].AlleleCnt;

    SECTION_ERR_EXTERN(check_oob);

    for (entry = 0; entry < PedTree->EntryCnt; entry++) {
        PedEntry = &(PedTree->Entry[entry]);
        iserr = 0;

        get_2alleles(PedEntry->Marker, locus, &all1, &all2);

        if (all1 > AlleleCnt) {
            SECTION_ERR(check_oob);
            warnvf("Ped %s: Entry %s allele1 %d out of bounds at locus %s.\n",
                   PedTree->Name,
                   ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
                   all1,
                   LTop->Locus[locus].LocusName);
            iserr = 1;
        }

        if (all2 > AlleleCnt) {
            SECTION_ERR(check_oob);
            warnvf("Ped %s: Entry %s allele2 %d out of bounds at locus %s.\n",
                   PedTree->Name,
                   ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
                   all2,
                   LTop->Locus[locus].LocusName);
            iserr=1;
        }

        if (iserr) {
            PedStatus->exceed_allcnt++;
            abortf=1;

            if (reset)
                set_2alleles(PedEntry->Marker, locus,
                             &LTop->Locus[locus], 0, 0);

            if (*reset_fp == NULL) {
                if (*first) {
                    *reset_fp = fopen(Mega2ResetRun, "w");
                    summary_time_stamp(mega2_input_files, *reset_fp, "");
                    *first = 0;
                } else
                    *reset_fp = fopen(Mega2ResetRun, "a");
                fprintf(*reset_fp, "Genotypes with out-of-bounds alleles %s:\n",
                        reset ? "RESET" : "ALLOWED");
                fprintf(*reset_fp, "Pedigree   Person    Marker\n");
            }
            fprintf(*reset_fp, "%s   %s   %s\n",
                    PedTree->Name,
                    ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
                    LTop->Locus[locus].LocusName);
        }
    }
    return abortf;
}

int check_half_type(ped_tree *PedTree, ped_status *PedStatus,
                    linkage_locus_top *LTop, int ped_num, int locus,
                    int uniqueids,
                    FILE **reset_fp, bool *first, int reset)
{
    int entry;
    ped_rec *PedEntry;
    int abortf = 0;
    int all1, all2;

    SECTION_ERR_EXTERN(check_half_type);

    for (entry = 0; entry < PedTree->EntryCnt; entry++) {
        PedEntry = &(PedTree->Entry[entry]);
//      iserr = 0;

        get_2alleles(PedEntry->Marker, locus, &all1, &all2);

        /* semi-typed */
        if ((all1 == 0 && all2 != 0) ||
            (all1 != 0 && all2 == 0)) {
            abortf=1;
            PedStatus->halftyped++;
            HalfTypedReset=1;
            SECTION_ERR(check_half_type);
            warnvf("Ped %s: Entry %s is half-typed at locus %s.\n",
                   PedTree->Name,
                   ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
                   LTop->Locus[locus].LocusName);

            if (reset)
                set_2alleles(PedEntry->Marker, locus,
                             &LTop->Locus[locus], 0, 0);

            if (*reset_fp == NULL) {
                if (*first) {
                    *reset_fp = fopen(Mega2ResetRun, "w");
                    summary_time_stamp(mega2_input_files, *reset_fp, "");
                    *first = 0;
                } else
                    *reset_fp = fopen(Mega2ResetRun, "a");
                fprintf(*reset_fp, "Half-typed genotypes %s:\n",
                        reset ? "RESET" : "ALLOWED");
                fprintf(*reset_fp, "Pedigree   Person    Marker\n");
            }
            fprintf(*reset_fp, "%s   %s   %s\n",
                    PedTree->Name,
                    ((uniqueids == 1)? PedEntry->LEntry->UniqueID : PedEntry->LEntry->OrigID),
                    LTop->Locus[locus].LocusName);
        }
    }
    return abortf;
}

int check_invalid_fam(ped_tree *PedTree, ped_status *PedStatus,
                      linkage_locus_top *LTop, int ped_num, int locus,
                      int uniqueids,
                      ped_rec **Sibs)
{
    int entry;
    ped_rec *PedEntry;
    int err = 0;

/*
    PedTree->EntryCnt == 1
    means there is no family structure so we should return 0.
    BUT check_inheritance also checks for anomalies on make with X-linked
*/

    for (entry = 0; entry < PedTree->EntryCnt; entry++) {
        PedEntry = &(PedTree->Entry[entry]);
        /* out of bounds */
        int all1, all2;
        get_2alleles(PedEntry->Marker, locus, &all1, &all2);

        if (!check_inheritance(PedEntry, LTop, locus, all1, all2,
                               PedTree->Name, uniqueids)) {
            PedStatus->genotype_invalid++;
            err = 1;
        }
    }

    if (PedTree->EntryCnt == 1) return err;

   /* Check sibships to see that they have only 4 distinct alleles,
       only need to check sibships of size greater than 2 */

    /* send in the index of locus in locus_top */
    if (check_sibship_alleles(PedTree, Sibs, locus,
                              LTop, &(PedStatus->checked_sibship),
                              uniqueids)) {
        PedStatus->genotype_invalid++;
        err = 1;
    }
    return err;
}

/* NM: locus errors :
   Fatal (2):  No name, <= 0 alleles, empty allele list,
   Non-fatal(1): Strange frequency values, do not sum to 1.0,
   less than 2 alleles, more than 2 allele, if PLINK
*/
/* check_locus()
 *
 * Check the locus and its alleles.
 *
 * Return non-zero if there is an error.
 * (We'll soon have a LocusStatus structure like PedStatus.)
 */
int check_locus(linkage_locus_rec *Locus,
		analysis_type analysis,
		int *plink_locus_num)

{
    int allele;
    double fsum = 0.0;
    int retval = 0;
    char messg[200];

    SECTION_ERR_EXTERN(check_locus);
    SECTION_ERR_EXTERN(biallele);

    if (Locus->LocusName == NULL) {
        SECTION_ERR(check_locus);
        errorf("Name is UNDEFINED");
        strcpy(Locus->LocusName, "UNDEF");
        retval= 2;
    }

    if (Locus->AlleleCnt <= 0) {
        SECTION_ERR(check_locus);
        sprintf(messg, "Locus %s has invalid number of alleles %d",
                Locus->LocusName, Locus->AlleleCnt);
        errorf(messg);
        retval = 2;
    }

    if ((Locus->Allele == NULL) && (Locus->AlleleCnt != 0)) {
        SECTION_ERR(check_locus);
        sprintf(messg, "Locus %s has corrupt or incomplete record",
                Locus->LocusName);
        errorf(messg);
        retval = 2;
    }

    /* check the alleles */
    if (Locus->Allele == NULL) {
        SECTION_ERR(check_locus);
        sprintf(messg, "Locus %s has NULL allele array",
                Locus->LocusName);
        errorf(messg);
        retval= 2;
    } else /* if (analysis != TO_PLINK) */  {
        /*  check allele frequencies (only if not PLINK) */

        for (allele = 0; allele < Locus->AlleleCnt; allele++) {
            // I can't remember which analysis programs didn't like it when some alleles had a
            // frequency of zero, but now that we are analyzing more and more rare variants,
            // we'll much more often encounter non-polymorphic markers where allele 2 has a
            // true frequency of zero.
            // It should still warn about negative allele frequencies (but I'm not
            // sure how that would occur - maybe an error in an input file?).

            if ((Locus->Allele[allele].Frequency < 0.0) || (Locus->Allele[allele].Frequency > 1.0)) {
                SECTION_ERR(check_locus);
                warnvf(//"Allele %d of locus %s has a negative or zero frequency %5.4f\n",
                       "Allele %s of locus %s has a frequency %5.4f which is not in the range 0.0 - 1.0 (inclusive)\n",
                       format_allele(Locus, allele+1),
                       Locus->LocusName, Locus->Allele[allele].Frequency);
                retval = 1;
            }

            fsum += Locus->Allele[allele].Frequency;
        }

        if (fsum == 0) {
            SECTION_ERR(check_locus);
            errorvf("Sum of allele frequencies of locus %s is zero\n",
                    Locus->LocusName);
            retval = 1;

        } else if (((fsum - 1.0) > 0.001001) || ((1.0 - fsum) > 0.001001)) {
            SECTION_ERR(check_locus);

            sprintf(messg, "Sum of allele frequencies of locus %s",
                    Locus->LocusName);
            errorf(messg);
            sprintf(messg, "is not within the allowed deviation of 0.001.");
            errorf(messg);
            sprintf(messg, "The sum is %6.5f.", fsum);
            errorf(messg);
            retval = 1;
        }
    }

    if (Locus->AlleleCnt < 2) {
        SECTION_ERR(check_locus);

        sprintf(messg, "Locus %s has less than 2 alleles", Locus->LocusName);
        warnf(messg);
        retval = 1;
    } else if (Locus->AlleleCnt > 2 /* && (analysis == TO_PLINK || analysis == IQLS) */) {
        SECTION_ERR(biallele);
        warnvf("Locus %s has more than 2 alleles; unusable with PLINK and IQLS\n", Locus->LocusName);
        (*plink_locus_num)--;
        retval = 1;
    }

    return retval;
}

static int autosomal_or_female_xlinked(int all1, int all2,
				       int *alleles, int *allelecnt,
				       int max_alleles)
{
    int a;
    int found1 = -1, found2 = -1;

    for (a=0; a < *allelecnt; a++) {
        if (all1 == alleles[a]) {
            found1 = a;
            break;
        }
    }
    if (found1 < 0) {
        if ((*allelecnt) < max_alleles) {
            alleles[(*allelecnt)++] = all1;
            found1 = (*allelecnt) - 1;
        } else {
            return 1;
        }
    }

    for (a=0; a < *allelecnt; a++) {
        if (all2 == alleles[a] && a != found1) {
            found2 = a;
            break;
        }
    }
    if (found2 < 0) {
        if ((*allelecnt) < max_alleles) {
            alleles[(*allelecnt)++] = all2;
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}

static int male_xlinked(int all, int *alleles, int *allelecnt)
{
    int a;
    int found = -1;

    for (a=0; a < *allelecnt; a++) {
        if (all == alleles[a]) {
            found = a;
            break;
        }
    }
    if (found < 0) {
        if ((*allelecnt) < 3) {
            alleles[(*allelecnt)++] = all;
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}

void mito_transmission_report(ped_top *PTop,
                              int *num_hetero,
                              int *num_non_maternal)
{
    int *locus1;

    *num_hetero=*num_non_maternal = 0;

    SECTION_ERR_INIT(check_mito);
    for (locus1 = &(PedTreeMito[0]); *locus1 != -99; locus1++) {
        int lloc = *locus1;
        int ped;
        for (ped = 0; ped < PTop->PedCnt; ped++) {
            int entry;
            for (entry = 0; entry < PTop->PedTree[ped].EntryCnt; entry++) {
                int A1, A2;
                get_2alleles(PTop->PedTree[ped].Entry[entry].Marker, lloc, &A1, &A2);
                if (! R(A1,A2)) {
                    (*num_hetero)++;
                    sprintf(err_msg,
                            "Ped %s: %s is a heterozygote at mitochondrial locus %s.",
                            PTop->PedTree[ped].Name,
                            ((InputFileFormat == ANNOTATED)?
                             PTop->PedTree[ped].Entry[entry].LEntry->UniqueID :
                             PTop->PedTree[ped].Entry[entry].LEntry->OrigID),
                            PTop->LocusTop->Locus[*locus1].LocusName);

                    SECTION_ERR(check_mito);
                    warnf(err_msg);
                }
                if (PTop->PedTree[ped].Entry[entry].Mother != NULL) {
                    int MA1, MA2;
                    get_2alleles(PTop->PedTree[ped].Entry[entry].Mother->Marker, lloc, &MA1, &MA2);
                    /* See if child matches mother's genotype */
                    if ((R(A1, MA1) && R(A2, MA2))  || (R(A2, MA1) && R(A1, MA2))) {
                        ;
                    } else {
                        (*num_non_maternal)++;
                        SECTION_ERR(check_mito);
                        warnvf("Ped %s: %s does not match maternal alleles at mitochondrial locus %s.\n",
                               PTop->PedTree[ped].Name,
                               ((InputFileFormat == ANNOTATED)?
                                PTop->PedTree[ped].Entry[entry].LEntry->UniqueID :
                                PTop->PedTree[ped].Entry[entry].LEntry->OrigID),
                               PTop->LocusTop->Locus[*locus1].LocusName);
                        // This will write the numeric allele always...
                        warnvf("\t%d/%d -> %d/%d\n",
                               MA1, MA2, A1, A2);
                        warnvf("\t%s  -> %s\n",
                               ((InputFileFormat == ANNOTATED)?
                                PTop->PedTree[ped].Entry[entry].Mother->LEntry->UniqueID :
                                PTop->PedTree[ped].Entry[entry].Mother->LEntry->OrigID),
                               ((InputFileFormat == ANNOTATED)?
                                PTop->PedTree[ped].Entry[entry].LEntry->UniqueID :
                                PTop->PedTree[ped].Entry[entry].LEntry->OrigID));
                    }
                }
            }
        }
    }
    SECTION_ERR_FINI(check_mito);
}
