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

#ifndef PEDTREE_H
#define PEDTREE_H


/*
  Certain data is labeled as Ordinal in this header
  (and perhaps in other headers). An ordinal is used
  here to specify a strictly positive integer (except
  in odd cases such as those requiring error handling).
  The ordinal may be reduced by 1 to create an array
  index.
*/

#define MALE_ID   1
#define FEMALE_ID 2
#define IS_MALE(e)    ((e).Sex == MALE_ID)
#define IS_FEMALE(e)  ((e).Sex == FEMALE_ID)
#define IS_FOUNDER(e) (((e).Mother == NULL) && ((e).Father == NULL))

#define AUTOSOMAL_KEYWORD "Autosomal"
#define X_LINKED_KEYWORD  "X-linked"

/*
  Note: Status != 2 doesn't imply that Status == 1
  The first two should be interpreted as:
  AFFECTED_STATUS   -> affected but nothing else assumed
  UNAFFECTED_STATUS -> not affected but nothing else assumed
*/
#define AFFECTED_STATUS   2
#define UNAFFECTED_STATUS 1
#define UNKNOWN_STATUS    0

#define IS_AFFECTED(e) ((e).Status == AFFECTED_STATUS)
#define IS_UNAFFECTED(e) ((e).Status == UNAFFECTED_STATUS)


/* These two are used by is_typed() and remove_untyped_affecteds()
   to set how to determine which affecteds are typed. */
/* NOTE: the mode should only be anded with TYPED_AT_ANY_LOCUS
   ((mode & TYPED_AT_ALL_LOCI) is always false). */
#define TYPED_AT_ALL_LOCI  0  /* bit 0 */
#define TYPED_AT_ANY_LOCUS 1  /* bit 0 */
#define APM_MODE           2  /* bit 1 */

/* structure containing the data we have on an allele */
typedef struct _allele_rec {
    double Frequency;
} allele_rec;


/* structure containing the information on a locus and
   it alleles' substructures */
typedef struct _locus_rec {
    char       *Name;
    int         AlleleCnt;
    int         chromosome;
    int         linkage_loc_num;
    linkage_locus_rec         *linkage_loc_rec;
    allele_rec *Allele;       /* will be Allele[] */
} locus_rec;


/* the top of the locus information, containing the
   locus structures */
typedef struct _locus_top {
    int        LocusCnt;
    locus_rec *Locus;     /* will be Locus[] */
    // SexLinked is defined as follows...
    // ==0 the markers found in the map file were on MT or unknown chromosomes.
    // ==1 says that at least one of the markers found in the map file was on the X or Y chromosome
    // ==2 says the same as ==1 but in addition that at least one of the markers was on the XY or an autosome
    int SexLinked;        /* non-zero if yes */
    linkage_locus_top *LTop;
} locus_top;


/* Structure containing the information on an individual */
/* Important note: The ID is an ordinal and must be equal
   to the index in the pedigree's Entry[] array (increased
   by one). Some functions, such as those copying this
   structure, reset the ID to a value (called UNDEF) which
   means it is undefined. This is because the calling function
   is more qualified to set it than the copying function.
   An UNDEF'd ID is used internally as a sign of an incomplete
   (or otherwise possibly invalid) member record. */
typedef struct _ped_rec {
    int              ID;     /* index in Entry[] */
    int Status;    /* 2 == affected, 1 == not affected, 0 == unknown  */
    /* ^This comment is not correct */
    int Orig_status; /* Original status */
    struct _ped_rec *Father;     /* pointer to father in Entry[] */
    struct _ped_rec *Mother;     /* pointer to mother in Entry[] */
    struct _ped_rec *PA_Sib;        /* pointer to next sib in Entry[] */
    struct _ped_rec *MA_Sib;        /* pointer to next sib in Entry[] */
    int              OffspringCnt;
    struct _ped_rec **Offspring; /* will be *Offspring[]     */
    int              Sex;        /* use macros above         */
    int             *Allele_1;   /* will be Allele_1[locus] */ /* Ordinal */
    int             *Allele_2;   /* will be Allele_2[locus] */ /* Ordinal */
    void            *TmpData;
    int             Father_aff;  /* 1 if father is affected, 0 otherwise */
    int             Mother_aff;  /* 1 if mother is affected, 0 otherwise */
    int             Off_aff;     /* Count of affected offspring */
    int             Dau_nrm;     /* Count of normal daughters */
    int             Son_nrm;     /* Count of normal sons */
    int             Dau_aff;     /* Count of affected daughters */
    int             Son_aff;     /* Count of affected sons */
    int             N_son;       /* Number of sons */
    int             N_dau;       /* Number of daughters */
    int             IsTyped;     /* lwa, just a flag to indicate if genotyped */
    int             connected;   /* flag for checking unconnected components  */
    int             sibship_checked; /* checked sibship for mendelian errors */
    linkage_ped_rec *LEntry;     /* pointer to corresponding linkage record */
} ped_rec;


/* structure containing the actual pedigree */
/* member Entry is array of ped_rec                         */
/* so we can treat entries as members of a tree or an array */
typedef struct _ped_tree {
    char     *Name;
    ped_rec *Entry;       /* will be Entry[]    */
    int      EntryCnt;
    int      AffectedCnt;
    ped_rec **Affected;   /* will be *Affected[] */
    ped_rec *Proband;
    void    *TmpData;
    int      Naff;        /* Total number of affecteds */
    int      Ntyped;      /* Number of people typed at at least one locus */
    int      Ntyped_all;  /* Number of people typed at ALL loci */
    int      Naff_typed;  /* Number of affecteds typed at at least one locus */
    int      Naff_typed_all;  /* Number of affecteds typed at ALL loci */
    int      Nsibships;       /* Number of sibships with at least one affected */
    int      Nsibs_typed;    /* Number of genotyped affected sibpairs */
    int      Nsibpairs;       /* Number of affected sibpairs */
    int      Nhalfsibpairs;   /* Number of affected half-sib pairs */
    int      Nhalfsibs_typed;   /* Number of genotyped affected half-sib pairs */
    int      Navuncularpairs; /* Number of affected avuncular pairs */
    int      Navunc_typed; /* Number of genotyped affected avuncular pairs */
    int      Nparchildpairs; /* Number of affected parent/child pairs */
    int      Nparchild_typed; /* Number of genotyped affected parent/child pairs */
    int      min_off_aff;     /* Min # of affecteds in a sibship */
    int      max_off_aff;     /* Max # of affecteds in a sibship */
    int      rel[33];         /* Count of relatives with kinship k/32 */
} ped_tree;


/* structure containing the pedigrees */
typedef struct _ped_top {
    int        PedCnt;
    ped_tree  *PedTree;   /* will be PedTree[] */
    locus_top *LocusTop;
} ped_top;

/* from check.h */

#define R(a, b) ((!a) || (!b) || (a == b))
#define HEMI2HOMO(a, b) if (a == 0 && b > 0) {a=b;}	\
    else if (a > 0 && b == 0) {b=a;}

#define ENTRY_DEFINED(e) (   (e.ID  != UNDEF)           \
                             && (e.Allele_1 != NULL)    \
                             && (e.Allele_2 != NULL))

/* from apmfile.h */

#define SL_FILE    0
#define ML_FILE   -1
#define MULT_FILE -2
#define IS_SL_FILE(t)   ((t) == SL_FILE)
#define IS_ML_FILE(t)   ((t) == ML_FILE)
#define IS_MULT_FILE(t) ((t) == MULT_FILE)


#endif
