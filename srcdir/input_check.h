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


#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H

#define MAX_PED_ERRORS 9
typedef struct _uniq_ids {
    char *id;
    int ped, per;
} uniq_id_type;

typedef struct _loc_list {
    int ped, locus, written;
    struct _loc_list *next;
} loc_list;

typedef struct _person_loc_list {
    int ped, locus, person;
    struct _person_loc_list *next;
} person_loc_list;

typedef struct _ped_status {
    int father_invalid, mother_invalid;
    int parents_not_defined;
    int incomplete_entry;
    int offspring_out_of_range;
    int genotype_invalid, halftyped, exceed_allcnt;
    int bad_ID;
    int bad_sex;
    int entry_unconnected;
    int unknown;
    int checked_sibship;
    loc_list *invalid_genotypes; /* for storing loci that are non-mendelian */
    person_loc_list *half_types; /* for storing half-typed genotypes*/
    person_loc_list *allele_outof_bounds; /* for storing out of bound alleles */
} ped_status;


typedef struct _file_info
{
    file_format infl_type;
    int line_count;
} file_info;

#endif
