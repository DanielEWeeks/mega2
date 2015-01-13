/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
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

*/

#ifndef ANALYSIS_EXT_H
#define ANALYSIS_EXT_H

#include "common.h"

#include "write_fbat_ext.h"
#include "write_plink_ext.h"
#include "write_eigenstrat_ext.h"
#include "write_pangaea_ext.h"
#include "write_beagle_ext.h"
#include "write_structure_ext.h"
#include "write_pseq_ext.h"

extern void  prog_name_to_num(char *prog_name, analysis_type *analysis);

extern CLASS_HAPLOTYPE          *HAPLOTYPE;
extern CLASS_LOCATION           *LOCATION;
extern CLASS_NONPARAMETRIC      *NONPARAMETRIC;
extern CLASS_IBD_EST            *IBD_EST;
extern CLASS_MISTYPING          *MISTYPING;
extern CLASS_MENDEL             *TO_MENDEL;
extern CLASS_ASPEX              *TO_ASPEX;
extern CLASS_CREATE_SUMMARY     *CREATE_SUMMARY;
extern CLASS_GENOTYPING_SUMMARY *GENOTYPING_SUMMARY;
extern CLASS_GENEHUNTER         *TO_GeneHunter;
extern CLASS_APM                *TO_APM;
extern CLASS_APM_MULT           *TO_APM_MULT;
extern CLASS_NUKE               *TO_NUKE;
extern CLASS_LIABLE_FREQ        *TO_LIABLE_FREQ;
extern CLASS_ALLELE_FREQ        *TO_ALLELE_FREQ;
extern CLASS_SLINK              *TO_SLINK;
extern CLASS_SPLINK             *TO_SPLINK;
extern CLASS_LOD2               *TO_LOD2;
extern CLASS_GENEHUNTERPLUS     *TO_GeneHunterPlus;
extern CLASS_SIMULATE           *TO_SIMULATE;
extern CLASS_SAGE               *TO_SAGE;
extern CLASS_TDTMAX             *TO_TDTMAX;
extern CLASS_SOLAR              *TO_SOLAR;
extern CLASS_HWETEST            *TO_HWETEST;
extern CLASS_VITESSE            *TO_VITESSE;
extern CLASS_LINKAGE            *TO_LINKAGE;
extern CLASS_ALLEGRO            *TO_Allegro;
extern CLASS_GHMLB              *TO_GHMLB;
extern CLASS_SAGE4              *TO_SAGE4;
extern CLASS_PREMAKEPED         *TO_PREMAKEPED;
extern CLASS_MERLIN             *TO_MERLIN;
extern CLASS_PREST              *TO_PREST;
extern CLASS_PAP                *TO_PAP;
extern CLASS_MERLINONLY         *TO_MERLINONLY;
extern CLASS_QUANT_SUMMARY      *QUANT_SUMMARY;
extern CLASS_LOKI               *TO_LOKI;
extern CLASS_MENDEL4            *TO_MENDEL4;
extern CLASS_SUP                *TO_SUP;
extern CLASS_PLINK              *TO_PLINK;
extern CLASS_MENDEL7_CSV        *TO_MENDEL7_CSV;
extern CLASS_CRANEFOOT          *CRANEFOOT;
extern CLASS_MEGA2ANNOT         *MEGA2ANNOT;
extern CLASS_IQLS               *IQLS;

extern CLASS_SIMWALK2           *TO_SIMWALK2;
extern CLASS_SUMMARY            *TO_SUMMARY;

extern CLASS_FBAT               *FBAT;
extern CLASS_PANGAEA            *PANGAEA;
extern CLASS_BEAGLE             *BEAGLE;
extern CLASS_EIGENSTRAT         *EIGENSTRAT;
extern CLASS_STRUCTURE          *STRUCTURE;
extern CLASS_PSEQ               *TO_PSEQ;

/*
extern CLASS_XXXXX              *XXXXX;
*/

#endif
