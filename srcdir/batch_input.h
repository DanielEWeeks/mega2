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

#ifndef BATCH_INPUT_H
#define BATCH_INPUT_H


/* allocate memory in chunks of 50 for mult_names */
#define OPT_CHUNK 50
#define KEYWORD_LEN 50
// batch file comment character...
#define COMMENT_CHAR '#'

/* To add a new output option:

   Define a new batch input integer variable - this should be the name of
   your new output option, prefixed by the string "BI_".  Assign it the
   next highest integer.
*/

/* These are the types of values options can take:
   YORN = yes or no (boolean)
   INT = integer
   CHAR = char
   FLOAT = float
   STRING  = character string
   INT_LIST =  list of integers
   NAME_LIST = list of strings separated by white-space
   FLOAT_LIST = list of floatng point numbers separated by white-space
*/
typedef enum _item_value_type {
    YORN, INT, FLOAT, CHAR, STRING, INT_LIST, NAME_LIST, FLOAT_LIST
} item_value_type;

typedef union batch_item_value_ {
    char *name;
    char copt;
    int option;
    int *mult_opts;
    char **mult_names;
    double fvalue;
    double *mult_fvalues;
} batch_item_value;

typedef struct _batch_item_ {
    int item_number, item_read;
    char keyword[KEYWORD_LEN];
    item_value_type value_type;
    batch_item_value value;
} batch_item_type;


/**
 These are the indices to the array "Mega2BatchItems".

 This is the way that batchfile items were maintained before the object
 system was developed for them. While the object system is aware of them
 for backward compatability, nothing should be added to this list and the
 object system used from here on out.
 */
#define Input_Pedigree_File 0
#define Input_Locus_File 1
#define Input_Map_File 2
#define Input_Omit_File 3
#define Input_Untyped_Ped_Option 4
#define Analysis_Option 5
#define Analysis_Sub_Option 6
#define Chromosome_Single 7
#define Chromosomes_Multiple_Num 8
#define Chromosomes_Multiple 9
#define Loci_Selected_Num 10
#define Loci_Selected 11
#define Trait_Single 12
#define Traits_Num 13
#define Traits_Loop_Over 14
#define Traits_Combine 15
#define Trait_Subdirs 16
#define Value_Missing_Quant 17
#define Value_Missing_Quant_On_Input 17
#define Value_Affecteds 18
#define Error_Loci 19
#define Error_Except_Loci 20
#define Error_Loci_Num 21
#define Error_Model 22
#define Error_Probabilities 23
#define Input_Do_Error_Sim 24
#define Default_Outfile_Names 25
#define Default_Reset_Invalid 26
#define Default_Other_Values 27
#define Default_Ignore_Nonfatal 28
#define Default_Rplot_Options 29
#define Xlinked_Analysis_Mode 30
#define Covariates_Num 31
#define Covariates_Selected 32
#define Output_Path 33
#define Count_Genotypes 34
#define Rplot_Statistics 35
#define Count_Halftyped 36
#define AlleleFreq_SquaredDev 37
#define Count_HWE_genotypes 38
#define Input_Frequency_File 39
#define Input_Penetrance_File 40
#define Output_Map_Num 41
#define Value_Missing_Allele 42
#define PLINK_Args 43
#define Input_Phenotype_File 44
#define Input_Aux_File 45
#define Value_Genetic_Distance_Index 46
#define Value_Base_Pair_Position_Index 47
#define Value_Genetic_Distance_SexTypeMap 48
#define Value_Missing_Quant_On_Output 49
#define Loop_Over_Chromosomes 50
#define Structure$PopDataPheno 51
#define Value_Marker_Compression 52
#define Input_Format_Type 53
#define Input_Path 54
#define Input_PLINK_Map_File 55
#define VCF_Args 56
#define VCF_Marker_Alternative_INFO_Key 57
#define Value_Missing_Affect_On_Input 58
#define Value_Missing_Affect_On_Output 59
#define Output_File_Stem 60
// Number of keywords which will appear in the batch file.
// The above list is zero based, so this would be the last number + 1
// You must update this number when you add an additional keyword.
#define NUM_KEYS 61

// THIS VARIABLE IS USED EVERYWHERE...
// It is defined in batch_input.cpp or by the new object system for
// backward compatability.
extern batch_item_type *Mega2BatchItems;


#define ITEM_READ(n)     (Mega2BatchItems[n].item_read == 1)

#ifndef DEFAULT_OUTFILES
#define DEFAULT_OUTFILES (Mega2BatchItems[Default_Outfile_Names].value.copt == 'y' || \
                          Mega2BatchItems[Default_Outfile_Names].value.copt == 'Y')
#endif

#ifndef DEFAULT_OPTIONS
#define DEFAULT_OPTIONS (Mega2BatchItems[Default_Other_Values].value.copt == 'y' || \
                         Mega2BatchItems[Default_Other_Values].value.copt == 'Y')
#endif

#ifndef DEFAULT_ROPTS
#define DEFAULT_ROPTS (Mega2BatchItems[Default_Rplot_Options].value.copt == 'y' || \
                       Mega2BatchItems[Default_Rplot_Options].value.copt == 'Y')
#endif


extern int batchINPUTFILES, batchANALYSIS, batchREORDER;
extern int batchTRAIT, batchQVALUE, batchAFFVALUE, batchERROR;

#endif /* BATCH_INPUT_H */
