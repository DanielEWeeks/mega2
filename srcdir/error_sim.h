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

#ifndef ERROR_SIM_H
#define ERROR_SIM_H

#include "typedefs.h"

/* definitions for error_sim.c and error_models.c */

typedef struct _gen {
    int all1, all2;
} genotype;


extern int ***GenotypeIndexes;
/* for indexing and reverse-indexing genotypes
   for all markers */

extern double UniformErrProb;
extern double SW2ErrProb[5];
extern int MarkerErrorProb; /* flag to specify if map file contains errors */

typedef enum _sw2_err_class {
    E1, E2, E3_homo, E3_hetero, E4, E5_homo, E5_hetero, Uniform, Marker
} sw2_error_class;

typedef struct _err_sim_out {
    int ped, entry;
    genotype original, new_genotype;
    sw2_error_class err_type;
} err_sim_out;


#endif
