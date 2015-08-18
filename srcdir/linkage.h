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

#ifndef LINKAGE_H
#define LINKAGE_H

#include <stdio.h>
#include "typedefs.h"

#define MALE_ID   1
#define FEMALE_ID 2
#define UNKNOWN_SEX 0
#define IS_LFOUNDER(e) ((e).Mother == 0 && (e).Father == 0)

typedef struct _loci_markers{
    int NumberInSubOrder;
    int total_files;
    int **MarkerSubsetMatrix;
} loci_markers;

typedef struct _allele_freq {
    int **genotype_matrix;
    int *Allele_Bin;
    double *allele_freq;
    double *input_freq;
    int HeteroGenoCount;
    int GenoCount;
    int TotalPeople;
    int HalfTyped;
    int TotalAlleles;
    int TotalKnownAlleles;
    int PhenoCount;
    int *Male_Allele_Bin;
    int *Female_Allele_Bin;
    double *male_allele_freq;
    double *female_allele_freq;
    int FemaleCnt;
    int MaleCnt;
    int FemaleHetzCnt;
    int MalesGenotyped;
    int FemalesGenotyped;
    int TotalMaleAlleles;
    int TotalFemaleAlleles;
    double FemaleHetzObs;
    double Hetz;
} allele_freq_struct;


// Locus->Class = TRAIT; Locus->Type = QUANT // a quantitative trait, a number (e.g., pulse rate)
// Locus->Class = TRAIT; Locus->Type = AFFECTION // affection status, a boolean (e.g., yes they have the disease)
// Locus->Classs = MARKER; Locus->Type = (BINARY | NUMBERED | XLINKED | YLINKED)
// Locus->Class = CLASS_UNSET; Locus->Type = TYPE_UNSET // uninitialized
// NOTE: BINARY is historical and not used any more.
typedef enum {
    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
} linkage_locus_type;

typedef enum {
    CLASS_UNSET, TRAIT, COVARIATE, MARKER
} linkage_locus_class;

typedef struct _linkage_allele_rec {
    double Frequency;
/*  char name[ALL_LEN+1];*/
    const char *name;
    int index; /* initially 1 .. num_alleles in ascending order */
} linkage_allele_rec;

#define allelecmp(a1, a2) (a1 != a2)

typedef struct _allele_prop {
    void *prop;
    int idx;
    char name[ALL_LEN+1];
} allele_prop;

#define allele_prop_allele(p) (((allele_prop *)(p))->name)
/*
#define allele2allele_prop_prop(p,t) ((t)(&( ((allele_prop *)(((char *)(p)) - ((int)allele_prop_allele(0))))->prop )))
#define allele2allele_prop_prop(p) ( ((allele_prop *)(((char *)(p)) - ((long)allele_prop_allele(0))))->prop )
*/
#define allele2allele_prop_prop(p) ( ((allele_prop *)(((char *)(p)) - ((char *)allele_prop_allele(0))))->prop )
#define allele2allele_prop_idx(p)  ( ((allele_prop *)(((char *)(p)) - ((char *)allele_prop_allele(0))))->idx )

/* added separate lists for average, male and female penetrances
   for linkage-format input, only the SexPen and Penetrances
   will be filled in. In annotated mode, only the last 3 will be
   filled in, then we can make the first two point to one the
   last three depending on chromosome or analysis mode. */

typedef struct _linkage_affection_class {
    unsigned char  MaleDef,  FemaleDef,  AutoDef;
    double *MalePen, *FemalePen, *AutoPen;
} linkage_affection_class;


typedef struct _linkage_quant_data {      /* Quantitative Trait */
    int ClassCnt;      /* should be at least 1 */
    double **Mean;
    double **Variance;
    double Multiplier;
} linkage_quant_data;


typedef struct _linkage_numbered_data {   /* Numbered Alelles   */
    int Recoded; /* will be 0 for not recoded, 1 for recoded */
    int NumAlleles; /* if recoded, # alleles for freq calculation */
    int SelectOpt; /* Type of individuals selected for freq calc */
    int EstimateFrequencies; /* 0 if not estimated, 1 if allele frequencies were estimated */
} linkage_numbered_data;


typedef struct _linkage_binary_data {     /* Binary Factor      */
    int FactorCnt;
    char **Factor;      /* just a boolean array, will be Factor[factor][allele] */
} linkage_binary_factor;


typedef struct _linkage_affection_data {  /* Affection Status   */
    int ClassCnt;    /* should be at least 1 */
    int PenCnt;      /* number of penetrances */
    linkage_affection_class *Class;   /* will be Class[] */
    int NumLabels;
    int *Labels;
} linkage_affection_data;

/* new */

typedef union _pheno_data {
    struct _linkage_quant_data Quant;
    struct _linkage_affection_data Affection;
} pheno_data;

typedef struct _pheno_rec {
    char *Name;
    pheno_data Props;
    int col_num; /* for annotated files only for now */ //X
} pheno_rec;

typedef union _marker_data {
    struct _linkage_numbered_data Numbered;
    struct _linkage_binary_data Binary;
} marker_data;

typedef struct _marker_rec {
    char *Name;
    marker_data Props;
    // UNKNOWN_POSITION is used in the following fields pos[ition|_male|_female], & error_prob
    // when a non-numeric or negative value is read, a value is missing, or un-initialized.
    double pos_avg, pos_male, pos_female, error_prob;
    int chromosome;
    int col_num;
} marker_rec;

/* end new */

typedef union _linkage_locus_data {
    struct _linkage_quant_data Quant;
    struct _linkage_numbered_data Numbered;
    struct _linkage_binary_data Binary;
    struct _linkage_affection_data Affection;
} linkage_locus_data;

typedef struct _linkage_locus_rec {
    char *Name;
    int AlleleCnt;
    linkage_allele_rec *Allele;     /* will be Allele[] */
    linkage_locus_type Type;
    linkage_locus_class Class;
    pheno_rec *Pheno;
    marker_rec *Marker;
    int number;
    int col_num; /* for annotated files only for now */ //X
} linkage_locus_rec;

/* Recombination mode constants */
/*  sex difference */
#define NO_SEX_DIFF       0
#define CONSTANT_SEX_DIFF 1
#define VARIABLE_SEX_DIFF 2
/* interference */
#define NO_INTERFERENCE   0
#define NO_MAPPING_FUN    1
#define WITH_MAPPING_FUN  2

/* added by nandita */
#define MLINK             5
#define LINKMAP           4

typedef struct mlink_params_ {
    double start_theta, increment, stop_theta;
    int num_evals;
} mlink_params;

typedef struct linkmap_params_ {
    int trait_marker,  num_evals;
    double stop_theta;
    double *recomb_frac;
} linkmap_params;

typedef struct slink_params_ {
    double *trait_positions;
} slink_params;

typedef union run_data {
    linkmap_params linkmap;
    mlink_params mlink;
    slink_params slink;
} run_data_type;

typedef struct _plink_info_type {
    int plinkf;

    int SNP_major;

    FILE *bed_filep;

    // CPK: In PLINK an allele is only one character. So we store alleles as follows...
    // alleles[n] == allele1
    // alleles[n+1] == allele2
    // where 'n' corresponds the the genotype as found in the .BIM file.
    char *alleles;

    // This is used to get the count of alleles in the .bim file.
    // It is the length of 'alleles' if it is filled
    int allele_count;
} plink_info_type;

typedef struct _linkage_locus_top {
    int LocusCnt;
    int PhenoCnt;
    int MarkerCnt;
    int NumPedigreeCols;
    linkage_locus_rec *Locus;       /* will be Locus[] */
    pheno_rec *Pheno;
    marker_rec *Marker;
    int RiskLocus, RiskAllele;
    int SexLinked, Program;    /* program is here because that's where it is */
    int MutLocus, Haplotype;
    double MutMale, MutFemale;
    char map_distance_type;  /* haldane(h) or kosambi(k) */
    /* the following is recombination information */
    int SexDiff, Interference;
    double *MaleRecomb;    /* will contain info if no sex difference also */
    double *FemaleRecomb;  /* only if there sex difference */
    double FlankRecomb;    /* used if interference == 1 */
    double FMRatio;       /* used if sex diff == 1 */
    run_data_type Run;  /* for mlink and linkmap */
    record_type PedRecDataType; /* recode, post-make etc. */
} linkage_locus_top;


typedef union _linkage_pedrec_data {
    struct {
        int Status, Class;
    } Affection;
    double Quant;
    struct {
        int Allele_1, Allele_2;
    } Alleles;
    struct {
        const char *Allele_1, *Allele_2;
    } RAlleles;

} linkage_pedrec_data;

/* new */

typedef union _pheno_pedrec_data {
    struct {
        int Status, Class;
    } Affection;
    double Quant;
} pheno_pedrec_data;

typedef struct _Alleles {
    int Allele_1, Allele_2;
} Alleles_int;

typedef struct _RAlleles {
    const char *Allele_1, *Allele_2;
} Alleles_str;

typedef union _marker_pedrec_data {
    Alleles_int Alleles;
    Alleles_str RAlleles;
} marker_pedrec_data;

typedef struct _marker_pedrec_char {
    unsigned char Allele_1, Allele_2;
} marker_pedrec_char;

/*
 * We're keeping Allele's independent of locus type
 * because they will probably be useful in all cases.
 */
typedef struct _linkage_ped_rec {
    char UniqueID[NAMELEN], OrigID[NAMELEN], FamName[NAMELEN];
    char PerPre[MAX_NAMELEN];
    int ID;              /* linkage ID */
    int Father, Mother;  /* linkage ID */
    int First_Offspring; /* linkage ID */
    int Next_PA_Sib, Next_MA_Sib;   /* linkage ID */
    int Sex;         /* 1 == male, 2 == female  */
    int OrigProband; /* Original proband field */
    pheno_pedrec_data   *Pheno;
    void                *Marker;
    void *TmpData;
    int genocnt;
    int Orig_status;  /* Original status */
    int Ngeno; /* Number of genotypes */
    int IsTyped; /* original linkage record status */
    double PercentTyped;
    int *loopbreakers;
    /* the next two fields are used to store original pednum and person ID
       when a pedigree is broken up into nuclear families */
    int ext_ped_num, ext_per_num;
    int MZTwin, DZTwin, Group; /* new fields added  6-27-05 */
    int rec_num; /* for storing line number of input file */

    /* This field stores the ids of loop-breakers corresponding to this entry */
} linkage_ped_rec;


/* modified by nandita (2/27/2000) to take a variable number of
   entries (at most 10 = max marriages) */
/* Keeps a record of which members are identical */
typedef struct _linkage_loop_rec {
    int Num;      /* 1 to N (equal to proband field - 1) */
    int Entry[MAXMARRIAGES];   /* ID's */
    int numEntry;
} linkage_loop_rec;

/* List is defined in list.c */

typedef struct _linkage_ped_tree {
    int Num;  /* May not be ordinal, but usually is. */
    int EntryCnt;
    linkage_ped_rec *Entry;         /* will be Entry[] */
    int IsTyped;
    char Name[MAX_NAMELEN];
    char PedPre[MAX_NAMELEN];
    int OriginalID;
    int origped; /* only to be used for conversion to nukepeds */
    int Proband;
    list *Loops; /* will be NULL if no loops */
} linkage_ped_tree;

typedef struct pre_makeped_record_ {
    char uniqueid[MAX_NAMELEN];
    int ped, indiv, father, mother, gender;
    int rec_num;
    pheno_pedrec_data   *pheno;
    void                *marker;
    int genocnt;

} pre_makeped_record;


typedef struct _person_node_ {
    char uniqueid[MAX_NAMELEN],  origid[MAX_NAMELEN], famname[MAX_NAMELEN];
    char perpre[MAX_NAMELEN];
    int node_id;
    int indiv, father, mother, gender;
    int proband;
    pheno_pedrec_data   *pheno;
    void                *marker;
    int genocnt;
    int from_marriage_node_id; /* offspring of marriage ? */
    int num_to_marriages, to_marriage_node_id[MAXMARRIAGES]; /* participates in marriages ? */
    int degree_genotyped;
    int removed;
    int loop_breaker_id;
    int first_pa_sib, first_ma_sib, first_off;
    int MZTwin, DZTwin, Group; /* new fields added  6-27-05 */

} person_node_type;

typedef struct _lcollapse_pass {
    int old, new_linkage;
} id_pair_type;


#endif /* LINKAGE_H */
