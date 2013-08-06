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


/* linkage_ext.h defines some extensions on the existing linkage data structure
 */


#ifndef LINKAGE_EXTERNAL_H
#define LINKAGE_EXTERNAL_H


/* This is a matrix of frequencies, where a non-zero
   frequency implies that a specific allele is present
   within a specific group.
   AlleleFreq[#groups][#alleles]
   This should not be very sparse in reality.
*/

typedef struct _ext_linkage_locus_rec {
    char **RAlleles;
    double **AlleleFreq;     /* will be AlleleFreq[][] */
    double *positions; /* one position per map */
    double *pos_male, *pos_female;
    int *AlleleCnt; /* AlleleCnts per group */
} ext_linkage_locus_rec;

typedef struct _ext_linkage_locus_top {
    int LocusCnt;
    int GroupCnt; /* Maximum number of groups in pedigree data */
    ext_linkage_locus_rec *EXLocus;

    char *map_functions; /* will be map_function[], one for each map */
    char **MapNames; /* Names of each map taken from Map file header */
    int MapCnt; /* For different locus maps */
//cpk    sex_map_types **SexMaps; /* record whether sex-maps were provided */
    int **SexMaps; /* record whether sex-maps were provided */
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    int **valid_map_p; // a predicate that tells if the map is valid
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
} ext_linkage_locus_top;

#endif
