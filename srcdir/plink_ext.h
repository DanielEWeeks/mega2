/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

#ifndef PLINK_EXT_H
#define PLINK_EXT_H

enum PLINK_FORMAT {
    not_plink_format = 0,
    binary_PED_format = 1,
    PED_format = 2
};

#ifndef _WIN
//cl: C requires that a structure of union has at least one member
typedef struct _plink_stats {
} plink_stats;
#endif /* _WIN */

// This structure is initialized in annotated_ped_file.c: PLINK_clr
typedef 
struct PLINK_s {
  int xcf;
  enum PLINK_FORMAT plink;
  int no_fid;
  int no_parents;
  int no_pheno;
  int map3;
  int cM;
  int missing_pheno; // default is 1
  int geneticMapType;
  int traitType;
  double pheno_value; // default is -9
  int phenotypeCoding;

  // satistics...
    int markers_total;
    int markers_included;
    int individuals;
    int phenotyped_individuals;
    int cases;
    int controls;
    int missing;
    int males;
    int females;
    int unspecified_sex;
    int founders;
    int non_founders;

  char trait[32]; // default is "default"
} PLINK_t;

extern PLINK_t PLINK;

extern void PLINK_clr(enum PLINK_FORMAT plink);

extern int PLINK_args(char *str, int xcf);

extern void PLINK_usage(int xcf);

extern void PLINK_str(char *str, int len);

#endif
