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

#ifndef MENU_VALUE_MISSING_H
#define MENU_VALUE_MISSING_H

Missing_Value missing_value;

Missing_Value  missing_values[] = {

    { "Allegro format",             TO_Allegro,
      /*
       * Reads pre-makeped LINKAGE format; expects a single affection status locus: 
         0 is unknown, 1 is unaffected and 2 is affected.
       */
      /*quant*/       "",  Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },

    { "Mega2 format",               MEGA2ANNOT,
      /*
       *
       */
      /*quant*/       "NA", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "NA", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      "",   Missing_Value::Any, Missing_Value::Fixed },

    { "ASPEX format",               TO_ASPEX,
      /*
       *
       */
      /*quant*/       "",  Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "",  Missing_Value::Any, Missing_Value::Fixed },

    { "Beagle format",              BEAGLE,
      /*
       *
       */
      /*quant*/       "?", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "?", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      "",  Missing_Value::Any, Missing_Value::Fixed },

    { "CRANEFOOT format",           CRANEFOOT,
      /*
       *
       */
      /*quant*/       "unknown", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "\\t11",   Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      "",        Missing_Value::Any, Missing_Value::Fixed },

    { "Eigenstrat format",          EIGENSTRAT,
      /*
       * For a quantitative variable: 'The value -100.0 signifies "missing data".'
       */
      /*quant*/       "-100.0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "-100",   Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "",        Missing_Value::Any, Missing_Value::Fixed },

    { "FBAT format",                FBAT,
      /*
       *
       */
      /*quant*/       " - ", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0",   Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0",   Missing_Value::Any, Missing_Value::Fixed },

    { "GeneHunter format",          TO_GeneHunter,
      /*
       * A 0 in any of the disease phenotype or marker genotype positions indicates missing data. 
         A - in the phenotype/covariate data indicates missing data.  
         NB: 0 is a real value that a phenotype may take on and DOES NOT represent missing data.
       */
      /*quant*/       "-", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },

    { "GeneHunter-Plus format",     TO_GeneHunterPlus,
      /*
       * Reads pre-makeped LINKAGE format; expects a single affection status locus: 
         1=UNAFFECTED, 2=AFFECTED.  A 0 in any of the disease phenotype or marker genotype
         positions indicates missing data.
       */
      /*quant*/       "-", Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },

    { "Homogeneity analyses",       TO_LOD2,
/*?*/ /*
       * LINKAGE format
       */
      /*quant*/       "-99.99", Missing_Value::Num, Missing_Value::Varies,
      /*affect*/      "0",      Missing_Value::Num, Missing_Value::Varies,
      /*allele*/      "0",      Missing_Value::Num, Missing_Value::Fixed },

    { "Test loci for HWE",          TO_HWETEST,
/*?*/ /*
       *
       */
      /*quant*/       "", Missing_Value::Any, Missing_Value::Varies,
      /*Affect*/      "", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "IQLS/Idcoefs format",        IQLS,
      /*
       * Analyzes a single affection status trait:  
         0=unknown, 1=unaffected, 2=affected; NN for missing genotype.
       */
      /*quant*/       "",   Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "0" , Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "NN", Missing_Value::Any, Missing_Value::Fixed },

    { "Linkage format",             TO_LINKAGE,
      /*
       *
       */
      /*quant*/       "0.0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "0",   Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0",   Missing_Value::Num, Missing_Value::Fixed },

    { "Loki format",                TO_LOKI,
/*?*/ /*
       * 'Versions after 2.4.5 have a default where `0' is the missing code
          for all categorical variables (pedigree, genotype and other discrete
          variables).  This will be overridden if any Missing command is found.'
	  Quant: Custom missing values defined via the MISSING command, so 'NA' can be used:    X any ch?
	  Affect: Defined via the MISSING command: X any ch?
       */
      /*quant*/       "x", Missing_Value::Any, Missing_Value::Varies,
      /*affect*/      "x", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "Mendel format",              TO_MENDEL7_CSV,
/*?*/ /*
       *
       */
      /*quant*/       " ", Missing_Value::Any, Missing_Value::Varies,
      /*affect*/      " ", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "",  Missing_Value::Any, Missing_Value::Fixed },

    { "Merlin format",              TO_MERLINONLY,
      /*
       * Do not give the user a choice - use 'x' for the missing value throughout.
         Quant: Quantitative traits are encoded as numeric values with X denoting missing values; any? Nch?
         Affect: X or 0 for missing  Any? Nch?
         Allele: either a 0, an X or an N can be used
       */
      /*quant*/      "x", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/     "x", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/     "x", Missing_Value::Any, Missing_Value::Fixed },

    { "Merlin/SimWalk2-NPL format", TO_MERLIN,
      /*
       * Do not give the user a choice - use 'x' for the missing value in the Merlin-format output files
         and blank for the missing value in the SimWalk2-format output files.
         Quant: x in Merlin files; blank in SimWalk2 files  any? Nch?
         Affect: x in Merlin files, blank in SimWalk2 files.
         Allele: x in Merlin files, blank in SimWalk2 files.
       */
      /*quant*/      "x", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/     "x", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/     "x", Missing_Value::Any, Missing_Value::Fixed },

    { "MLBQTL format",              TO_GHMLB,
      /*
       * Based on GeneHunter
       */
      /*quant*/       "-99.0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "0",     Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0",     Missing_Value::Num, Missing_Value::Fixed },

    { "PANGAEA MORGAN format",      PANGAEA,
      /*
       * The 'input pedigree record traits' command is used when real (non-integer)
       data for traits are included. 'A real value with integer part 999 indicates a missing value.'
       Allele: Marker data are specified for each marker locus as a pair of integer alleles, and 
               '0' indicates a missing value.
       */
      /*quant*/       "999.0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "0",     Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0",     Missing_Value::Any, Missing_Value::Fixed },

    { "Create nuclear families",    TO_NUKE,
/*?*/ /*
       *
       */
      /*quant*/       "0", Missing_Value::Num, Missing_Value::Varies,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "Old SAGE format",            TO_SAGE,
/*?*/ /*
       *
       */
      /*quant*/       "", Missing_Value::Any, Missing_Value::Varies,
      /*affect*/      "", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "PAP format",                 TO_PAP,
      /*
       * "Optional records in header.dat specify a missing value code (default -9999)". 
         "To indicate an unknown phenotype, the phenotype equals -9999 or the missing value
          code specified for the variable in header.dat."
         Quant:
         Affect:
         Allele:
       */
      /*quant*/       "-9999", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0",     Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      "0",     Missing_Value::Any, Missing_Value::Fixed },

    { "PLINK format",               TO_PLINK,
      /*
       * The missing phenotype value for quantitative traits is, by default, -9 
         (this can also be used for disease traits as well as 0)
       */
      /*quant*/       "-9", Missing_Value::Num, Missing_Value::Varies,
      /*affect*/      "-9", Missing_Value::Num, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "PSEQ format",                TO_PSEQ,
      /*
       * Missing value is defined in the header lines of the phenotype file, and can be set to strings, so use 'NA'.
       */
      /*quant*/       "NA", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0",  Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "Pre-makeped format",         TO_PREMAKEPED,
      /*
       *
       */
      /*quant*/       "0", Missing_Value::Num, Missing_Value::Varies,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Varies,
      /*allele*/      "",  Missing_Value::Any, Missing_Value::Fixed },

    { "PREST format",               TO_PREST,
      /*
       * In the ped file, expects a single affection status locus in column 6, 
         coded 0=unknown, 1=unaffected, 2=affected.  However, as this column is not used,
         entries in this column could just as well be all zeros.
       */
      /*quant*/       "",  Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },

    { "SAGE format",                TO_SAGE4,
      /*
       * Missing value code defined by the 'missing' attribute of the trait declaration.  Can be a string, so use 'NA'.
         different missing codes may be used for different variables.
       */
      /*quant*/       "NA", Missing_Value::Any, Missing_Value::Varies,
      /*affect*/      "NA", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "SIMULATE format",            TO_SIMULATE,
      /*
       * If the first locus is an affection status locus, it is carried through unchanged; 
         otherwise affection and trait loci are simulated and are expected to be coded as
         '1' = simulate a value; '0' indicates do not simulate a value. 
       */
/*?*/ /*quant*/       "", Missing_Value::Any, Missing_Value::Varies,
      /*affect*/      "", Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "SimWalk2 format",            TO_SIMWALK2,
      /*
       * Missing values for any field are represented by blanks
       */
      /*quant*/       " ", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      " ", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      " ", Missing_Value::Any, Missing_Value::Fixed },

    { "SLINK format",               TO_SLINK,
      /*
       *
       */
      /*quant*/       "0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Any, Missing_Value::Fixed },

    { "SOLAR format",               TO_SOLAR,
      /*
       * SOLAR uses comma-delimited files, so uses blanks for all missing values.
       */
      /*quant*/       " ", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      " ", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      " ", Missing_Value::Any, Missing_Value::Fixed },

    { "SPLINK format",              TO_SPLINK,
      /*
       *
       */
/*?*/ /*quant*/       "",    Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "0",   Missing_Value::Any, Missing_Value::Varies,
      /*allele*/      "0/0", Missing_Value::Any, Missing_Value::Fixed },

    { "Structure format",           STRUCTURE,
      /*
       * Phenotype information is not actually used in structure, but optionally can be in 
         the input file as an integer if the "PHENOTYPE" parameter is set to "1".
       */
      /*quant*/       "",   Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "-9", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/      "-9", Missing_Value::Any, Missing_Value::Fixed },

    { "Create summary files",       TO_SUMMARY,
      /*
       * This option does not generate output files.
       */
      /*quant*/       "", Missing_Value::Any, Missing_Value::NoQnt,
      /*affect*/      "", Missing_Value::Any, Missing_Value::NoQnt,
      /*allele*/      "", Missing_Value::Any, Missing_Value::NoQnt },

    { "SUP format",                 TO_SUP,
      /*
       *
       */
      /*quant*/       "",  Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },

    { "Vintage Mendel format",      TO_MENDEL,
      /*
       *
       */
      /*quant*/      " ", Missing_Value::Any, Missing_Value::Fixed,
      /*affect*/     " ", Missing_Value::Any, Missing_Value::Fixed,
      /*allele*/     " ", Missing_Value::Any, Missing_Value::Fixed },

    { "Vitesse format",             TO_VITESSE,
      /*
       *
       */
      /*quant*/       "0", Missing_Value::Num, Missing_Value::Fixed,
      /*affect*/      "0", Missing_Value::Num, Missing_Value::Fixed,
      /*allele*/      "0", Missing_Value::Num, Missing_Value::Fixed },



    { "APM format [DISABLED]",      TO_APM,
      /*
       *
       */
      /*quant*/       "", Missing_Value::Numeric, Missing_Value::Change,
      /*affect*/      "", Missing_Value::Numeric, Missing_Value::Change,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "APM-MULT format [DISABLED]", TO_APM_MULT,
      /*
       *
       */
      /*quant*/       "", Missing_Value::Numeric, Missing_Value::Change,
      /*affect*/      "", Missing_Value::Numeric, Missing_Value::Change,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

    { "TDTMax analyses [DISABLED]", TO_TDTMAX,
      /*
       *
       */
      /*quant*/       "", Missing_Value::Numeric, Missing_Value::Change,
      /*affect*/      "", Missing_Value::Numeric, Missing_Value::Change,
      /*allele*/      "", Missing_Value::Any, Missing_Value::Fixed },

};

#endif
