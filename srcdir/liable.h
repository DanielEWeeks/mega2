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

#ifndef LIABLE_H
#define LIABLE_H

#define MARKERLIABLE(n)                MarkerLiability[n]
#define TOTAL(n)                       MarkerLiability[n].Total
#define GTOTAL(n)                      MarkerLiability[n].GTotal
#define ALLELETOTAL(n,l)               MarkerLiability[n].AlleleTotal[l]
#define GENOTOTAL(n,l)                 MarkerLiability[n].GenoTotal[l]
#define ALLELELIABLE(n,l)              MARKERLIABLE(n).AlleleLiability[l]
#define LIABLESTATUS(n,l,j,s)          ALLELELIABLE(n,l).Liabls[j].Status[s]
#define STATUSPERCENT(n,l,j,s)         ALLELELIABLE(n,l).Liabls[j].StatusPercent[s]
#define GENOLIABLE(n,l)                MARKERLIABLE(n).GenoLiability[l]
#define GLIABLESTATUS(n,l,j,s)         GENOLIABLE(n,l).Liabls[j].Status[s]
#define GSTATUSPERCENT(n,l,j,s)        GENOLIABLE(n,l).Liabls[j].StatusPercent[s]
#define MARKERTOTAL(n,l,s)             MARKERLIABLE(n).MarkerTotal[l][s]
#define MARKERGTOTAL(n,l,s)            MARKERLIABLE(n).MarkerGTotal[l][s]


typedef struct _liabl_type {
    int Class;
    int Status[3];
    float StatusPercent[3];
} liabl_type;

typedef struct _allele_liabl_type {
    int Allele;
    liabl_type *Liabls;
} allele_liabl_type;

typedef struct _geno_liabl_type {
    int Allele1, Allele2;
    liabl_type *Liabls;
} geno_liabl_type;

typedef struct _marker_liabl_type {
    int MarkerLocus;
    int num_alleles;
    allele_liabl_type *AlleleLiability;
    int *AlleleTotal; /* one total for each allele */
    int **MarkerTotal; /* column total for each liability, each status */
    int Total; /* Total number of alleles */
    int num_genos;
    geno_liabl_type *GenoLiability;
    int *GenoTotal; /* one total for each genotype */
    int **MarkerGTotal; /* column total for each liability, each status */
    int GTotal; /* Total number of genotypes */
} marker_liabl_type;

typedef struct _liable_allele_dist{
    int locus_AFFECTION;
    int num_liabls;
    int *loci_NUMBERED;
    marker_liabl_type *MarkerLiability;
} liable_allele_dist;

#endif
