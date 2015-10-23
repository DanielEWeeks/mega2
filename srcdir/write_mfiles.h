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

/* moved macros to here */

#ifndef WRITE_MFILES_H
#define WRITE_MFILES_H

#include "typedefs.h"

#define MENDEL_AFFECTION fprintf(fp, "%-8s",  strtail(Locus->LocusName, MENDEL_MAX_LOCUS_NAME_LEN));	\
    if ((sex_linked  == 1) ||                                           \
        (sex_linked == 2 && Locus->Marker->chromosome == SEX_CHROMOSOME)) {     \
	fprintf(fp, "X-LINKED");                                        \
    }                                                                   \
    else {                                                              \
        fprintf(fp, "AUTOSOME");                                        \
    }                                                                   \
    fprintf(fp, "%2d", Locus->AlleleCnt);                               \
    if (simwalk2==0)                                                    \
	fprintf(fp, "%2d  <= Locus,Type,# Alleles,# Phenotypes \n",     \
		2*Locus->Pheno->Props.Affection.ClassCnt);                            \
    else                                                                \
	fprintf(fp, " 2  <= Locus,Type,# Alleles,# Phenotypes \n");     \
    for (allele = 0; allele < Locus->AlleleCnt; allele++)  {            \
	fprintf(fp, "%7d %8.6f      <= Allele, Frequency\n", allele + 1, \
		Locus->Allele[allele].Frequency);                       \
    }                                                                   \
    if (simwalk2==0) {                                                  \
	for (tmpi = 1; tmpi <= 2; tmpi++)                               \
            for (tmpi2 = 1; tmpi2 <= Locus->Pheno->Props.Affection.ClassCnt; tmpi2++)  { \
                if (Locus->Pheno->Props.Affection.ClassCnt > 1) {                     \
                    if (analysis == TO_MENDEL4) {                       \
                        sprintf(trait_phen, "%1d_%d", tmpi, tmpi2);     \
                        fprintf(fp, "%7s ", trait_phen);                \
                    }                                                   \
                    else {                                              \
                        fprintf(fp, "%2d%3d   ", tmpi, tmpi2);          \
                    }                                                   \
                }                                                       \
                else                                                    \
                    fprintf(fp, "%2d      ", tmpi);                     \
                numgen = 0;                                             \
                if (tmpi == 1)	{                                       \
                    for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) \
                        if (1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen] > 0.0) \
                            numgen++;                                   \
                }                                                       \
                else {                                                  \
                    for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) \
                        if (Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen] > 0.0) \
                            numgen++;                                   \
                }                                                       \
                fprintf(fp, "%2d # Genotypes", numgen);                 \
                for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++)   \
                    if (tmpi == 1)                                      \
                        fprintf(fp,  "%5.2f",                           \
                                1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen]); \
                    else                                                \
                        fprintf(fp, "%5.2f",                            \
                                Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen]); \
                fprintf(fp, " <= Penetrances\n");                       \
                if (tmpi == 2)  {                                       \
                    if (Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[0] > 0.0) \
                        fprintf(fp, "  1/ 1\n");                        \
                    if (Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[1] > 0.0) \
                        fprintf(fp, "  1/ 2\n");                        \
                    if (Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[2] > 0.0) \
                        fprintf(fp, "  2/ 2\n");                        \
                }                                                       \
                else  {                                                 \
                    if (1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[0] > 0.0) \
                        fprintf(fp, "  1/ 1\n");                        \
                    if (1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[1] > 0.0) \
                        fprintf(fp, "  1/ 2\n");                        \
                    if (1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[2] > 0.0) \
                        fprintf(fp, "  2/ 2\n");                        \
                }                                                       \
            }                                                           \
    }                                                                   \
    else {                                                              \
	fprintf(fp, " 1       3 # Genotypes\n");                        \
	fprintf(fp, "  1/ 1\n");                                        \
	fprintf(fp, "  1/ 2\n");                                        \
	fprintf(fp, "  2/ 2\n");                                        \
	fprintf(fp, " 2       3 # Genotypes\n");                        \
	fprintf(fp, "  1/ 1\n");                                        \
	fprintf(fp, "  1/ 2\n");                                        \
	fprintf(fp, "  2/ 2\n");                                        \
    }


#ifdef csh_lines
#undef csh_lines
#endif

#define csh_lines        fprintf(filep, "if ( -e RERUN-0%d.BCH ) then\n", i); \
    fprintf(filep, "   rm -f BATCH2.DAT\n");                            \
    fprintf(filep, "   cp RERUN-0%d.BCH BATCH2.DAT\n", i);              \
    fprintf(filep, "   rm RERUN-0%d.BCH \n", i)

#ifdef else_lines
#undef else_lines
#endif

#define else_lines 	fprintf(filep, "else \n");              \
    fprintf(filep, "  echo No more runs were suggested\n");     \
    fprintf(filep, "  echo \"************\"\n");                \
    fprintf(filep, "  exit 0\n");                               \
    fprintf(filep, "endif\n")

#ifdef errf_lines
#undef errf_lines
#endif

#define errf_lines  fprintf(filep, "if ( -e ERROR-0%d.TXT ||  -e ERROR-00.TXT ) then \n", i); \
    fprintf(filep, "   echo Please check ERROR files.\n");              \
    fprintf(filep, "   exit 1\nendif\n");                               \
    fprintf(filep, "echo \"************\"\n")


#ifdef TRAIT_PHEN1
#undef TRAIT_PHEN1
#endif

#define TRAIT_PHEN1 if (Locus->Pheno->Props.Affection.ClassCnt > 1) { \
        sprintf(trait_phen, "1_%d", tmpi2+1);           \
    }                                                   \
    else {                                              \
        strcpy(trait_phen, "1");                        \
    }

#ifdef TRAIT_PHEN2
#undef TRAIT_PHEN2
#endif

#define TRAIT_PHEN2 if (Locus->Pheno->Props.Affection.ClassCnt > 1) { \
        sprintf(trait_phen, "2_%d", tmpi2+1);           \
    }                                                   \
    else {                                              \
        strcpy(trait_phen, "2");                        \
    }


#define MENDEL_QUANT

#ifdef TAB_TO_GENO_COL
#undef TAB_TO_GENO_COL
#endif

#define TAB_TO_GENO_COL	    if (loc > 0 && (loc % 3) == 0) {    \
        fprintf(filep, "\n");                                   \
        for(i=0; i < 3*pwid+3; i++) {                           \
            fprintf(filep, " ");                                \
        }                                                       \
    }                                                           \

#endif
