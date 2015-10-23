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

#ifndef MRECODE_H
#define MRECODE_H

#include "typedefs.h"

#define DOM_G11  0.05
#define DOM_G12  0.90
#define DOM_G22  0.90
#define REC_G11  0.05
#define REC_G12  0.05
#define REC_G22  0.90

/* structure to store liability class information for each class */
typedef struct liability_class_ {
//    char name[LIABILITY_LEN+1]; // previous to 4.5.8
    int num;
    double *autosomal_pen;
    double *male_pen;
    double *female_pen;
} liability_class;

/* structure to store frequency information for each allele */

typedef struct allele_freq_ {
/*  char name[ALL_LEN+1];*/
    const char *AlleleName;
    int index; /* initially 1 .. num_alleles in ascending order */
    double freq;
    int count;
    int founder_count, random_count, unique_count, everyone_count;

} allele_freq_type;

/* List of classes for each affection status locus */

typedef struct class_list_ {
    struct class_list_ *next;
    liability_class l_class;
} class_list_type;

/* list of alleles for each marker */

typedef struct allele_list_ {
    struct allele_list_ *next;
    allele_freq_type allele_freq;
} allele_list_type;

/* to store num_alleles and num_classes */
typedef struct _mdata_ {
    int num_alleles;
    int num_classes;
} recode_marker_data;


/* The structure below can be used for markers, quant and aff */
/* NM: 7-15-2008, added new field recode_allele */

typedef struct marker_ {
    allele_list_type *first_allele;
    recode_marker_data data;
    int estimate_frequencies, recode_alleles;
    int num_total_alleles;
    int num_people;
/* These are cumulative counts, e.g.
   num_random = num_founders + num randomly selected individuals */
    int num_founders, num_random, num_unique, num_everyone;
    int num_half_typed;
    int ht_founders, ht_random, ht_unique, ht_everyone;
} marker_type;

//r fields ? && above
typedef struct pheno_ {
    allele_list_type *first_allele;
    class_list_type *class_list;
    recode_marker_data data;
    int estimate_frequencies;
} pheno_type;

typedef struct {
    const char *all1;
    const char *all2;
    int num;
} allelecnt;

#endif
