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


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "common.h"
#include "typedefs.h"

#include "batch_input.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "genetic_utils_ext.h"
#include "reorder_loci_ext.h"
#include "user_input_ext.h"
#include "plink_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc warnf
      genetic_utils_ext.h:  write_time
       reorder_loci_ext.h:  parse_string
         user_input_ext.h:  invalid_analysis
              plink_ext.h:  PLINK_args
              utils_ext.h:  EXIT draw_line get_line
*/

batch_item_type *Mega2BatchItems;

int batchINPUTFILES, batchANALYSIS, batchREORDER;
int batchTRAIT, batchAFFVALUE, batchERROR;

/* end of exports */
/*    Current keywords

      0    Input_Pedigree_File
      1.    Input_Locus_File
      2.    Input_Map_File
      3.    Input_Omit_File
      4    Input_Untyped_Ped_Option
      5    Analysis_Option
      6    Analysis_Sub_Option
      7    Chromosome_Single
      8    Chromosomes_Multiple_Num
      9    Chromosomes_Multiple
      10    Loci_Selected_Num
      11    Loci_Selected
      12    Trait_Single
      13    REMTraits_Num
      14    Traits_Loop_Over
      15    Traits_Combine
      16    Trait_Subdirs
      17    Value_Missing_Quant_On_Input (formerly Value_Missing_Quant which is deprecated)
      18    Value_Affecteds
      19    Error_Loci
      20    Error_Except_Loci
      21    Error_Loci_Num
      22    Error_Model
      23    Error_Probabilities
      24    Input_Do_Error_Sim
      25    Default_Outfile_Names
      26    Default_Reset_Invalid
      27    Default_Other_Values
      28    Default_Ignore_Nonfatal
      29    Default_Rplot_Options
      30    Xlinked_Analysis_Mode
      31    Covariates_Num
      32    Covariates_Selected
      33    Output_Path
      34    Count_Genotypes
      35    Rplot_Statistics
      36    Count_Halftyped
      37    AlleleFreq_SquaredDev
      38    Count_HWE_genotypes
      39    Input_Frequency_File
      40    Input_Penetrance_File
      41    Output_Map_Num
      42    Value_Missing_Allele
      43    PLINK
      44    Input_Phenotype_File
      45    Input_Binary_File
      46    Value_Genetic_Distance_Index
      47    Value_Base_Pair_Position_Index
      48    Value_Genetic_Distance_SexTypeMap
      49    Value_Missing_Quant_On_Output
      50    Loop_Over_Chromosomes
*/

static char keywords[NUM_KEYS][KEYWORD_LEN] = {
    "Input_Pedigree_File",
    "Input_Locus_File",
    "Input_Map_File",
    "Input_Omit_File",
    "Input_Untyped_Ped_Option",
    "Analysis_Option",
    "Analysis_Sub_Option",
    "Chromosome_Single",
    "Chromosomes_Multiple_Num",
    "Chromosomes_Multiple",
    "Loci_Selected_Num",
    "Loci_Selected",
    "Trait_Single",
    "REMTraits_Num",
    "Traits_Loop_Over",
    "Traits_Combine",
    "Trait_Subdirs",
    "Value_Missing_Quant_On_Input", // "Value_Missing_Quant" is an alias
    "Value_Affecteds",
    "Error_Loci",
    "Error_Except_Loci",
    "Error_Loci_Num",
    "Error_Model",
    "Error_Probabilities",
    "Input_Do_Error_Sim",
    "Default_Outfile_Names",
    "Default_Reset_Invalid",
    "Default_Other_Values",
    "Default_Ignore_Nonfatal",
    "Default_Rplot_Options",
    "Xlinked_Analysis_Mode",
    "REMCovariates_Num",
    "Covariates_Selected",
    "Output_Path",
    "Count_Genotypes",
    "Rplot_Statistics",
    "Count_Halftyped",
    "AlleleFreq_SquaredDev",
    "Count_HWE_genotypes",
    "Input_Frequency_File",
    "Input_Penetrance_File",
    "REMOutput_Map_Num",
    "Value_Missing_Allele_Aff",
    "PLINK",
    "Input_Phenotype_File",
    "Input_Binary_File",
    "Value_Genetic_Distance_Index",
    "Value_Base_Pair_Position_Index",
    "Value_Genetic_Distance_SexTypeMap",
    "Value_Missing_Quant_On_Output",
    "Loop_Over_Chromosomes",
    "Structure.PopDataPheno"
};

void batch_file_doc(FILE *batchfp)
{
    int key;

    fprintf(batchfp, "%c Lines beginning with %c are comments.\n",
            COMMENT_CHAR, COMMENT_CHAR);
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
    fprintf(batchfp, "%c Currently implemented keywords:\n",
            COMMENT_CHAR);
    for (key =0; key < NUM_KEYS; key++) {
        if (strncasecmp(keywords[key], "REM", (size_t) 3)) {
            fprintf(batchfp, "%c   %d) %s \n", COMMENT_CHAR,
                    key+1, keywords[key]);
        }
    }
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
    fprintf(batchfp, "%c Restrictions on usage:\n", COMMENT_CHAR);
    fprintf(batchfp,
            "%c   Use either Chromosome_Single or \n",
            COMMENT_CHAR);
    fprintf(batchfp,
            "%c     Chromosome_Single and Loci_Selected or\n",
            COMMENT_CHAR);
    fprintf(batchfp,
            "%c     Chromosomes_Multiple and Chromsomes_Multiple_Num\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
    fprintf(batchfp, "%c  Use either Trait_Single or\n", COMMENT_CHAR);
    fprintf(batchfp, "%c             Traits_Loop_Over or\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c             Traits_Combine.\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c  Keyword    Trait_Subdirs is allowed only if\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c  keyword    Traits_Loop_Over is defined.\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
    fprintf(batchfp, "%c     Error_Loci and Error_Loci_Num or\n",
            COMMENT_CHAR);
    fprintf(batchfp,
            "%c     Error_Except_Loci and Error_Loci_Num.\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
    fprintf(batchfp, "%c Default settings :\n", COMMENT_CHAR);
    fprintf(batchfp, "%c Default_Reset_Invalid: \n", COMMENT_CHAR);
    fprintf(batchfp, "%c   \"yes\"= set inconsistent genotypes to 0 and continue.\n", COMMENT_CHAR);
    fprintf(batchfp, "%c   \"no\"= continue without setting inconsistent genotypes to 0.\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c Default_Ignore_Nonfatal : \n", COMMENT_CHAR);
    fprintf(batchfp, "%c   Don't pause for other non-fatal errors in input:\n", COMMENT_CHAR);
    fprintf(batchfp, "%c Default_Other_Values: \n", COMMENT_CHAR);
    fprintf(batchfp, "%c   Use Mega2's default values instead of asking user\n", COMMENT_CHAR);
    fprintf(batchfp, "%c   inside analysis-option menus. \n", COMMENT_CHAR);
    fprintf(batchfp, "%c Xlinked_Analysis_Mode:\n", COMMENT_CHAR);
    fprintf(batchfp, "%c   Set x-linked based on chromosome number(human).\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c   treat all markers as autosomal,\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c   or treat all markers as x-linked.\n", COMMENT_CHAR);
    fprintf(batchfp, "%c Count_Genotypes and Count_Halftyped\n", COMMENT_CHAR);
    fprintf(batchfp, "%c   These correspond to the Select Individuals menu\n",
            COMMENT_CHAR);
    fprintf(batchfp,
            "%c   within the allele-recoding step.\n",
            COMMENT_CHAR);
    fprintf(batchfp, "%c \n", COMMENT_CHAR);
}

void create_batchfile(void)
{
    char syscmd[100];

    sprintf(syscmd, "/bin/cp /dev/null %s\n", Mega2Batch);
    System(syscmd);
    BatchFileCreated=1;
}


/* The type of keyword is used only for writing out batch items */
void batchfile_init_Mega2BatchItems(void)
{
    int i;

    Mega2BatchItems = CALLOC((size_t)NUM_KEYS, batch_item_type);

    for (i=0;  i < NUM_KEYS; i++) {
        Mega2BatchItems[i].item_number=i;
        Mega2BatchItems[i].item_read=0;
        strcpy(Mega2BatchItems[i].keyword, keywords[i]);
        switch(i) {
        case /* 0 */ Input_Pedigree_File:
        case /* 1 */ Input_Locus_File:
        case /* 2 */ Input_Map_File:
        case /* 3 */ Input_Omit_File:
        case /* 5 */ Analysis_Option:
        case /* 6 */ Analysis_Sub_Option:
        case /* 11 */ Loci_Selected:
        case /* 14 */ Traits_Loop_Over:
        case /* 15 */ Traits_Combine:
        case /* 32 */ Covariates_Selected:
        case /* 33 */ Output_Path:
        case /* 35 */ Rplot_Statistics:
        case /* 39 */ Input_Frequency_File:
        case /* 40 */ Input_Penetrance_File:
        case /* 42 */ Value_Missing_Allele:
        case /* 43 */ PLINK_Args:
        case /* 44 */ Input_Phenotype_File:
        case /* 45 */ Input_Binary_File:
        case /* 51 */ Structure$PopDataPheno:
            Mega2BatchItems[i].value.name = CALLOC((size_t)FILENAME_LENGTH, char);
            strcpy(Mega2BatchItems[i].value.name, "");
            Mega2BatchItems[i].value_type = STRING;
            break;
        case /* 49 */ Value_Missing_Quant_On_Output:
            Mega2BatchItems[i].value.name = CALLOC((size_t)FILENAME_LENGTH, char);
            strcpy(Mega2BatchItems[i].value.name, "0");
            Mega2BatchItems[i].value_type = STRING;
            break;
        case /* 9 */ Chromosomes_Multiple:
        case /* 19 */ Error_Loci:
        case /* 20 */ Error_Except_Loci:
            Mega2BatchItems[i].value.mult_opts=NULL;
            Mega2BatchItems[i].value_type = INT_LIST;
            break;
        case /* 4 */ Input_Untyped_Ped_Option:
            Mega2BatchItems[i].value.option=-1;
            Mega2BatchItems[i].value_type = INT;
            break;
        case /* 46 */ Value_Genetic_Distance_Index:
        case /* 47 */ Value_Base_Pair_Position_Index:
        case /* 48 */ Value_Genetic_Distance_SexTypeMap:
            Mega2BatchItems[i].value.option=-1;
            Mega2BatchItems[i].value_type = INT;
            break;
        case /* 7 */ Chromosome_Single:
        case /* 41 */ Output_Map_Num:
            Mega2BatchItems[i].value.option=1;
            Mega2BatchItems[i].value_type = INT;
            break;
        case /* 8 */ Chromosomes_Multiple_Num:
        case /* 10 */ Loci_Selected_Num:
        case /* 12 */ Trait_Single:
        case /* 13 */ Traits_Num:
        case /* 21 */ Error_Loci_Num:
        case /* 31 */ Covariates_Num:
            Mega2BatchItems[i].value.option=0;
            Mega2BatchItems[i].value_type = INT;
            break;
        case /* 24 */ Input_Do_Error_Sim:
        case /* 25 */ Default_Outfile_Names:
        case /* 27 */ Default_Other_Values:
        case /* 28 */ Default_Ignore_Nonfatal:
        case /* 29 */ Default_Rplot_Options:
        case /* 36 */ Count_Halftyped:
        case /* 50 */ Loop_Over_Chromosomes:
            Mega2BatchItems[i].value.copt='n';
            Mega2BatchItems[i].value_type = YORN;
            break;
        case /* 26 */ Default_Reset_Invalid:
            /* Reset invalid genotypes */
            Mega2BatchItems[i].value.copt='y';
            Mega2BatchItems[i].value_type = YORN;
            break;
        case /* 16 */ Trait_Subdirs:
            /* trait numbers and trait directory names */
        case /* 18 */ Value_Affecteds:
            /* affection labels */
            Mega2BatchItems[i].value.mult_names=NULL;
            Mega2BatchItems[i].value_type = NAME_LIST;
            break;
        case /* 17 */ Value_Missing_Quant_On_Input:
        case /* 37 */ AlleleFreq_SquaredDev:
            Mega2BatchItems[i].value.fvalue=0.0;
            Mega2BatchItems[i].value_type = FLOAT;
            break;
        case /* 22 */ Error_Model:
            Mega2BatchItems[i].value.copt='X';
            Mega2BatchItems[i].value_type = CHAR;
            break;
        case /* 23 */ Error_Probabilities:
            Mega2BatchItems[i].value.mult_fvalues=NULL;
            Mega2BatchItems[i].value_type = FLOAT_LIST;
            break;
        case /* 30 */ Xlinked_Analysis_Mode:
        case /* 34 */ Count_Genotypes:
        case /* 38 */ Count_HWE_genotypes:
            Mega2BatchItems[i].value.option=2;
            Mega2BatchItems[i].value_type = INT;
            break;
        default:
            strcpy(Mega2BatchItems[i].value.name, "");
            Mega2BatchItems[i].value_type = STRING;
            break;
        }
    }
}


/* Read the batch file and set the keywords,
   this is the slow version,  but avoids restrictions on the
   order of items inside the batch file.
   The keywords are designed so that we have only a few categories
   that can be tested right away.
   The batch is read twice, so that some keywords which depend
   on other keywords are not forced  into sequential ordering.
*/

#ifdef READ_ITER1
#undef READ_ITER1
#endif

#define READ_ITER1(it) ((it == -1)? 0 : (Mega2BatchItems[it].item_read == 1) && (iter == 1))


static void missing_item_goto_menu(int batch_item, const char *menu_name)
{
    if (batch_item >= 0) {
        warnvf("Missing necessary batch item %s, going to %s.\n",
                Mega2BatchItems[batch_item].keyword, menu_name);
    } else {
        warnvf("Going to %s.\n", menu_name);
    }
}

void check_batch_items(void)
{
    int read_a_section;
    /* Input files */

    if (ITEM_READ(/* 0 */ Input_Pedigree_File) && 
        (ITEM_READ(/* 1 */ Input_Locus_File) || Mega2BatchItems[PLINK_Args].item_read) &&
        ITEM_READ(/* 2 */ Input_Map_File) && ITEM_READ(/* 4 */ Input_Untyped_Ped_Option)) {
        batchINPUTFILES=1;
        mssgf("Input filenames and missing value indicator read in from batch file.");
    } else {
        if (!ITEM_READ(/* 0 */ Input_Pedigree_File)) {
            missing_item_goto_menu(0, "Input menu");
        } else if (!ITEM_READ(/* 1 */ Input_Locus_File) && !Mega2BatchItems[PLINK_Args].item_read){
            missing_item_goto_menu(1, "Input menu");
        } else if (!ITEM_READ(/* 2 */ Input_Map_File)) {
            missing_item_goto_menu(2, "Input menu");
        } else if (!ITEM_READ(/* 4 */ Input_Untyped_Ped_Option)) {
            missing_item_goto_menu(4, "Input menu");
        }
        batchINPUTFILES=0;
    }

    if (ITEM_READ(/* 5 */ Analysis_Option)) {
        batchANALYSIS = 1;
        mssgf("Analysis option read in from batch file.");
    } else {
        missing_item_goto_menu(5, "Analysis option menu");
        batchANALYSIS = 0;
    }

    if (ITEM_READ(/* 7 */ Chromosome_Single) || (ITEM_READ(/* 8 */ Chromosomes_Multiple_Num) && ITEM_READ(/* 9 */ Chromosomes_Multiple))) {
        batchREORDER = 1;
        mssgf("Markers, chromosome(s) and read in from batch file.");
    } else {
        warnf("Locus selections not specified in batch file.");
        missing_item_goto_menu(-1, "Reorder menu");
        batchREORDER = 0;
    }

    if (ITEM_READ(/* 12 */ Trait_Single) || ITEM_READ(/* 14 */ Traits_Loop_Over) || ITEM_READ(/* 15 */ Traits_Combine)) {
        batchTRAIT = 1;
        mssgf("Trait selection(s) read in from batch file.");
    } else {
        warnf("Trait selections not specified in batch file.");
        missing_item_goto_menu(-1, "Trait selection menu");
        batchTRAIT = 0;
    }

    batchAFFVALUE = ((ITEM_READ(/* 18 */ Value_Affecteds))? 1: 0);

    if ((ITEM_READ(/* 19 */ Error_Loci) || ITEM_READ(/* 20 */ Error_Except_Loci)) && ITEM_READ(/* 21 */ Error_Loci_Num)) {
        batchERROR = 1;
        mssgf("Error parameters read in from batch file.");
    } else  {
        batchERROR = 0;
    }

    read_a_section =
      batchINPUTFILES + batchANALYSIS + batchREORDER + batchTRAIT +
      (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].item_read ? 1 : 0) +
      batchAFFVALUE + batchERROR;

    if (read_a_section) {
        draw_line();
    }
}

#ifndef NEW_BATCH


/* item dependent on multiple items, item dependent on any one
   of multiple items being defined */
static void missing_items(int batch_item, int list_item)
{
    errorf("Not enough items to read:");
    errorvf("Item %s specifies %d items, \n",
            Mega2BatchItems[batch_item].keyword,
            Mega2BatchItems[batch_item].value.option);
    errorvf("Item %s has less than %d items.\n",
            Mega2BatchItems[list_item].keyword,
            Mega2BatchItems[batch_item].value.option);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void missing_dependency_keyword(int dependent_item, int dependent_on_item)
{
    errorvf("Keyword %s is defined but keyword %s is missing from batch file.\n",
            Mega2BatchItems[dependent_item].keyword,
            Mega2BatchItems[dependent_on_item].keyword);
    errorvf("Keyword %s depends on keyword %s.\n",
            Mega2BatchItems[dependent_item].keyword,
            Mega2BatchItems[dependent_on_item].keyword);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void missing_mult_dependency_keyword(int dependent_item, int dependent_on_items[],
					    int num_dependent_items, int and_or)
{
    int i;
    /* and=1, or=0 */
    if (and_or) {
        sprintf(err_msg,
                "Keyword %s requires all these keywords to be defined:",
                Mega2BatchItems[dependent_item].keyword);
    } else {
        sprintf(err_msg,
                "Keyword %s depends on one of these keywords to be defined:",
                Mega2BatchItems[dependent_item].keyword);
    }

    for (i=0; i < num_dependent_items; i++) {
        strcat(err_msg, Mega2BatchItems[dependent_on_items[i]].keyword);
    }
    errorf(err_msg);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void mult_decl(int item)
{
    warnvf("Found multiple occurrences of keyword %s,\n",
            Mega2BatchItems[item].keyword);
    warnf("Using the first definition, and ignoring all later definitions.");
}

static void malformed_batch_line(char *keyword)
{
    errorvf("Malformed line for keyword %s.\n", keyword);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

/* global variables */
static char analysis_name[50]="", sub_analysis_name[100]="";

static void set_batch_items(char *batch_file_name, int iter, analysis_type *analysis)
{
    char nextline[FILENAME_LENGTH];
    char keyword[KEYWORD_LEN];
    char value[FILENAME_LENGTH];
    char *intstr;

    int j, it = -1, version4;
    FILE *fp;

    version4 = 0;
    /* First the the verseion number */
    fp = fopen(batch_file_name, "r");
    if (fp == (FILE *)NULL) {
        errorvf("could not open %s for reading!\n", batch_file_name);
        EXIT(FILE_READ_ERROR);
    }

    (void)fgets(nextline, FILENAME_LENGTH - 1, fp);
    if (!strncasecmp(nextline, "#Version4.4", (size_t) 11)) {
        version4 = 1;
    }
    fclose(fp);

    fp = fopen(batch_file_name, "r");
    while (!feof(fp)) {
        get_line(fp, nextline);
        // if the line is not the empty string, though this would match a line with one space on it...
        // better would be to search for the leack of a comment character in the first column and
        // an equal sign somewhere in the line....
        if (strcmp(nextline, "") != 0) {
            // get everything before the '=' sign...
            // NOTE: there must be NO white space before the start of the line or between the keyword of the '='?!
            sscanf(strtok(nextline, "="), "%s", keyword);
            // get everything after the '=' sign...
            intstr = strtok(NULL, "\n");
            if (intstr == NULL) {
                malformed_batch_line(keyword);
            }
            strcpy(value, intstr);
            if (!strncmp(keyword, "Input", (size_t) 5)) {
                /* file names section */
                it = -1;
                for (j = /* 0 */ Input_Pedigree_File ; j <= /* 5 - 1 */ Input_Untyped_Ped_Option ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (it == -1) {
                    if (!strcmp(keyword, keywords[24])) {
                        it = /* 24 */ Input_Do_Error_Sim ;
                    } else if (!strcmp(keyword, keywords[39])) {
                        it = /* 39 */ Input_Frequency_File ;
                    } else if (!strcmp(keyword, keywords[40])) {
                        it = /* 40 */ Input_Penetrance_File ;
                    } else if (!strcmp(keyword, keywords[44])) {
                        it = /* 44 */ Input_Phenotype_File ;
                    } else if (!strcmp(keyword, keywords[45])) {
                        it = /* 45 */ Input_Binary_File ;
                    }
                } else {
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                }
                if ((it >= /* 0 */ Input_Pedigree_File  && it < /* 4 */ Input_Untyped_Ped_Option )
                    || (it == /* 39 */ Input_Frequency_File )
                    || (it == /* 40 */ Input_Penetrance_File )
                    || (it == /* 44 */ Input_Phenotype_File )
                    || (it == /* 45 */ Input_Binary_File )) {
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%s", Mega2BatchItems[it].value.name);
                } else if (it == /* 4 */ Input_Untyped_Ped_Option ) {
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 24 */ Input_Do_Error_Sim ) {
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%c", &(Mega2BatchItems[it].value.copt));
                }
            }

            /*  Analysis option section */
            else if (!strncmp(keyword, "Analysis", (size_t) 8)) {
                it = -1;
                for (j = /* 5 */ Analysis_Option ; j <= /* 6 */ Analysis_Sub_Option ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }

                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }
                Mega2BatchItems[it].item_read = 1;
                if (iter == 1) {
                    if (it == /* 5 */ Analysis_Option ) {
                        sscanf(value, "%s", analysis_name);
                    } else {
                        //sscanf(value, "%s", sub_analysis_name);
                        // Want the remainder of the line not just the first substring...
                        strcpy(sub_analysis_name, value);
                    }
                }

                if (iter == 2) {
                    if (it == /* 5 */ Analysis_Option ) {
                        // This is where the 'analysis' variable get's assiged too...
                        prog_name_to_num(analysis_name, analysis);
                        if (*analysis == NULL || !(*analysis)->is_enabled()) {
                            unknown_prog(analysis_name);
                        }
                        /* Also check if this analysis has a sub-option and if sub-option was
                           provided in the batch file */
                    }
#if 0
                    if ((*analysis)->has_sub_options() && !ITEM_READ(/* 6 */ Analysis_Sub_Option )) {
                            // For backwards compatibility this is not required with PLINK so don't warn about it...
                        if (*analysis == TO_PLINK || *analysis == FBAT) {
                                // Since it is not required for backwards compatibility we fudge it here...
                                (*analysis)->_suboption = 1; // 'lgen
                        } else {
                            errorvf("%s sub-program needs to be specified.\n",
                                    analysis_name);
			    errorvf("Required keyword %s missing from batch file.\n",
				    Mega2BatchItems[Analysis_Sub_Option].keyword);
			    EXIT(BATCH_FILE_ITEM_ERROR);
                        }
                        strcpy(Mega2BatchItems[/* 5 */ Analysis_Option].value.name, analysis_name);
                    } else {
                        if ((*analysis)->has_sub_options()) {
                            (*analysis)->sub_prog_name_to_sub_option(sub_analysis_name, analysis);
                        }
                    }
#endif /* 0 */
                    strcpy(Mega2BatchItems[/* 5 */ Analysis_Option].value.name, analysis_name);
                }
                // The Analysis_Sub_Option is not guaranteed to have been processed till the third pass
                // if it has been specified before the Analysis_Option in the batch file.
                if (iter == 3 && it == /* 5 */ Analysis_Option ) {
                    if ((*analysis)->has_sub_options()) {
                        if (! ITEM_READ(/* 6 */ Analysis_Sub_Option )) 
                            warnvf("%s: using empty (i.e. default) sub-program name.\n",
                                    analysis_name);
                        (*analysis)->sub_prog_name_to_sub_option(sub_analysis_name, analysis);
                    }
#ifdef DEBUGOPT
                    printf("option number %d, name %s.\n",
                           Mega2BatchItems[/* 5 */ Analysis_Option].value.option, analysis_name);
                    if (HAS_SUB_OPTION(Mega2BatchItems[/* 5 */ Analysis_Option].value.option)) {
                        printf("sub-option number %d, name %s.\n",
                               Mega2BatchItems[/* 6 */ Analysis_Sub_Option].value.option, sub_analysis_name);
                    }
                    exit(0);
#endif /* DEBUGOPT */
                }
            }
            /* Chromosome section  */
            else if (!strncmp(keyword, "Chromosome", (size_t) 10)) {
                it = -1;
                for (j = /* 7 */ Chromosome_Single ; j <= /* 9 */ Chromosomes_Multiple ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }
                if (it == /* 7 */ Chromosome_Single ) { /*  single chromosome */
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 8 */ Chromosomes_Multiple_Num ) { /*  Number of chromosomes */
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }

                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 9 */ Chromosomes_Multiple ) { /*  list  of chromosomes */
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    if (!Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].item_read) {
                        /* if iter == 1 do nothing */
                        if (iter == 2) {
                            /* declare error */
                            missing_dependency_keyword(9, 8);
                        }
                    } else {
		        int ch;
#ifdef CFREE
                        if (iter == 1) continue;
#endif /* CFREE */
                        if (Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].value.option <= 0) {
                            invalid_value_field(8);
                        }
                        Mega2BatchItems[it].value.mult_opts
                            =CALLOC((size_t)Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].value.option, int);
                        Mega2BatchItems[it].item_read = 1;
                        intstr = strtok(value, " ");
                        if (intstr == NULL) {
                            missing_items( /*8 */ Chromosomes_Multiple_Num , /* 9 */ Chromosomes_Multiple );
                        }
                        if ((ch=STR_CHR(intstr)) != -1) {
                            Mega2BatchItems[it].value.mult_opts[0] = ch;
                        } else {
                            errorf("For batch file item 'Chromosomes_Multiple'.");
                            errorvf("A chromosome must be a positive integer less then %d or one of the\n", lastautosome+4);
                            errorf("following: X, Y, XY, MT, U.");
                            EXIT(BATCH_FILE_ITEM_ERROR);
                            //errorf("A chromosome must be a positive integer or one of the following: X, Y, XY, MT, U.");
                        }
                        for (j = 1; j < Mega2BatchItems[/* 8 */ Chromosomes_Multiple_Num].value.option; j++) {
                            intstr = strtok(NULL, " ");
                            if (intstr == NULL) {
                                missing_items( /*8 */ Chromosomes_Multiple_Num , /* 9 */ Chromosomes_Multiple );
                            }
                            if ((ch=STR_CHR(intstr)) != -1) {
                                Mega2BatchItems[it].value.mult_opts[j] = ch;
                            } else {
                                errorf("For batch file item 'Chromosomes_Multiple'.");
                                errorvf("A chromosome must be a positive integer less then %d or one of the\n", lastautosome+4);
                                errorf("following: X, Y, XY, MT, U.");
                                EXIT(BATCH_FILE_ITEM_ERROR);
                                //errorf("A chromosome must be a positive integer or one of the following: X, Y, XY, MT, U.");
                            }
                        }
                    }
                }
            }
            /* Marker selection section */
            else if (!strncmp(keyword, "Loci", (size_t) 4)) {
                it = -1;
                for (j =  /*10 */ Loci_Selected_Num ; j <= /* 11 */ Loci_Selected ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }
                if (it == /* 10 */ Loci_Selected_Num ) { /* Number of markers */
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 11 */ Loci_Selected ) { /* Marker  numbers */
                    strcpy(Mega2BatchItems[it].value.name, value);
                    Mega2BatchItems[it].item_read = 1;
                }
            }
            /* Trait selection section */
            else if ((!strncmp(keyword, "Trait", (size_t) 5)) ||
		     (!strncmp(keyword, "Covar", (size_t) 5))) {
                it = -1;
                for (j = /* 12 */ Trait_Single ; j <= /* 16 */ Trait_Subdirs ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (it == -1) {
                    if (!strcmp(keyword, keywords[31])) {
                        it = /* 31 */ Covariates_Num ;
                    } else if (!strcmp(keyword, keywords[32])) {
                        it = /* 32 */ Covariates_Selected ;
                    }
                }
                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }

                if (it == /* 12 */ Trait_Single  || it == /* 13 */ Traits_Num  || it == /* 31 */ Covariates_Num ) {
                    /* Single trait number, number of traits or number of covars */
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 14 */ Traits_Loop_Over  || it == /* 15 */ Traits_Combine  || it == /* 32 */ Covariates_Selected ) {
                    /* Remove dependency on the number items, set these directly
                       from the string using parse_string */
                    Mega2BatchItems[it].item_read = 1;
                    strcpy(Mega2BatchItems[it].value.name, value);
                } else if (it == /* 16 */ Trait_Subdirs ) {
                    int num_dirs = 0;
                    /* Check if the value has the keyword Use trait names */
                    if (!strcasecmp(value, "use trait names")) {
                        Mega2BatchItems[/* 16 */ Trait_Subdirs].item_read = 1;
                        Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names
                            = CALLOC((size_t)1, char *);
                        Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[0]
                            = CALLOC((size_t)16, char);
                        strcpy(Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[0],
                               "use trait names");
                    } else {
                        /* User defined Trait sub directories,
                           num_traits must  be set using the trait-selection item */
                        if (!Mega2BatchItems[/* 14 */ Traits_Loop_Over].item_read && iter == 2) {
                            /* error */
                            missing_dependency_keyword(16, 14);
                        } else if (!Mega2BatchItems[/* 14 */ Traits_Loop_Over].item_read) {
                            continue;
                        } else {
#ifdef CFREE
                            if (iter == 1) continue;
#endif
                            parse_string(Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name, NULL,
                                         &num_dirs, (int) LARGE, 0);
                            Mega2BatchItems[/* 16 */ Trait_Subdirs].item_read = 1;
                            Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names
                                = CALLOC((size_t)num_dirs, char *);
                            Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[0]
                                = CALLOC((size_t)FILENAME_LENGTH, char);
                            strcpy(Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[0],
                                   strtok(value, " "));
                            for (j = 1; j < num_dirs; j++) {
                                Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[j]
                                    = CALLOC((size_t)FILENAME_LENGTH, char);
                                strcpy(Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[j],
                                       strtok(NULL, " "));
                            }
                        }
                    }
                }
            } else if (!strncmp(keyword, "Value", (size_t) 5)) {
                if (!strcmp(keyword, "Value_Missing_Quant_On_Input") ||
                    !strcmp(keyword, "Value_Missing_Quant") // deprecated but not yet forgotten...
		    ) {
                    // "Value_Missing_Quant_On_Input" is an alias for "Value_Missing_Quant" which has been deprecated
                    it = /* 17 */ Value_Missing_Quant_On_Input ;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].item_read = 1;
                    // QMISSING is the internal numeric value of NA.
                    if (strcasecmp(value, "NA") == 0) {
                        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue = QMISSING;
                    } else {
                        // Interesting, but MissingQuant was not assigned here in the past!
                        // NOTE: the use of 'strtok' here assumes that there is no comment following the value
                        MissingQuant = atof(strtok(value, " "));
                        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue = MissingQuant;
                    }
                } else if (!strcmp(keyword, "Value_Missing_Quant_On_Output")) {
                    // This is the missing value code that is used on output analysises those targets
                    // support that support optional missing value codes.
                    if (READ_ITER1(/* 49 */ Value_Missing_Quant_On_Output)) {
                        // if it's the first pass and has been declared more than once, complain...
                        mult_decl(/* 49 */ Value_Missing_Quant_On_Output);
                        continue;
                    }

                    // The 'analysis' variable can only be guaranted to hold data until the third pass.
                    // This is because the sub-analysis option may not have been read till the end of the second pass.
                    if (iter == 3) {
                        // We need to determine if the option will support this batch item.
                        if ((*analysis)->output_quant_can_define_missing_value()) {
                            Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read = 1;
                            strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name, value);
                        } else {
                            if (ITEM_READ(/* 6 */ Analysis_Sub_Option))
                                // NOTE: ProgName does not exist at this point.
                                warnvf("'Analysis_Option' = '%s' with 'Analysis_Sub_Option' = '%s' does not allow the definition of 'Value_Missing_Quant_On_Output'.\n",
                                       analysis_name, sub_analysis_name);
                            else
                                warnvf("'Analysis_Option' = '%s' does not allow the definition of 'Value_Missing_Quant_On_Output'.\n",
                                       analysis_name);
                        }
                        // We need to determine if this option must be numeric...
                        if ((*analysis)->output_quant_must_be_numeric()) {
                            char *end;
                            // The converted value is thrown away, because we just want to know if it is valid...
                            (void) strtod((const char *)value, &end);
                            // See user_input.cpp:set_missing_quant_input() for an explaination of this test...
                            if (strlen(value) == 0 || strlen(end) != 0 || errno == ERANGE) {
                                // Conversion of the entire string was not successful or some other error...
                                errorvf("'Analysis_Option' = '%s' while specifying 'Value_Missing_Quant_On_Output'.\n",
                                        (*analysis)->_name);
                                errorvf("Analysis type requires quantitative values to be numeric.\n");
                                EXIT(BATCH_FILE_ITEM_ERROR);
                            }
                            // NOTE: It is not possible at this point to check to see if the
                            // missing quantitative value defined here is in the input quantitative
                            // phenotype data. This will be done later in the routine
                            // user_input.cpp:set_missing_quant_output()
                        } // if ((*analysis)->output_quant_must_be_numeric()) {
                    } // if (iter == 3) {
                } else if (!strcmp(keyword, "Value_Affecteds")) {
                    int num_trs = 0;
                    char *affdata_str;
                    it = /* 18 */ Value_Affecteds ;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    if (ITEM_READ(/* 12 */ Trait_Single) || ITEM_READ(/* 14 */ Traits_Loop_Over) || ITEM_READ( /*15 */ Traits_Combine)) {
                        if (!ITEM_READ(/* 18 */ Value_Affecteds)) {
                            Mega2BatchItems[/* 18 */ Value_Affecteds].item_read = 1;
                            if (ITEM_READ(/* 12 */ Trait_Single)) {
                                /* single trait */
                                Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names
                                    = CALLOC((size_t)1, char *);
                                Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[0]
                                    = CALLOC((size_t)FILENAME_LENGTH, char);
                                sscanf(value, "%s",
                                       Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[0]);
                            } else if (ITEM_READ(/* 14 */ Traits_Loop_Over) || ITEM_READ(/* 15 */ Traits_Combine)) {
                                /* Traits Loop Over */
                                if (ITEM_READ(/* 14 */ Traits_Loop_Over)) {
                                    parse_string(
                                        Mega2BatchItems[/* 14 */ Traits_Loop_Over].value.name,
                                        NULL, &num_trs, (int) LARGE, 0);
                                } else if (ITEM_READ(/* 15 */ Traits_Combine)) {
                                    parse_string(
                                        Mega2BatchItems[/* 15 */ Traits_Combine].value.name,
                                        NULL, &num_trs, (int) LARGE, 0);
                                }

                                Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names
                                    = CALLOC((size_t)num_trs, char *);
                                affdata_str = NULL;
                                affdata_str = strtok(value, " \n");
                                j = 0;
                                while (affdata_str != NULL) {
                                    Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[j]
                                        = strdup(affdata_str);
                                    affdata_str = strtok(NULL, " \n");
                                    j++;
                                }
                            }
                        }
                    } else {
                        if (iter == /* 2 */ Input_Map_File ) {
                            int opt_array[3] = { 12, 14, 15 };
                            /* error */
                            missing_mult_dependency_keyword(18, opt_array, 3, 0);
                        }
                    }
                }
#ifdef USER_UNKNOWN
                else if (!strcmp(keyword, "Value_Missing_Allele_Aff")) {
                    it = /* 42 */ Value_Missing_Allele;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[/* 42 */ Value_Missing_Allele].item_read=1;
                    strcpy(Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name, value);
                }
#endif
                else if (!strcmp(keyword, "Value_Genetic_Distance_Index")) {
                    it = /* 46 */ Value_Genetic_Distance_Index;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[it].item_read=1;
                    sscanf(value, "%d", &Mega2BatchItems[it].value.option);
                }
                else if (!strcmp(keyword, "Value_Base_Pair_Position_Index")) {
                    it = /* 47 */ Value_Base_Pair_Position_Index;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[it].item_read=1;
                    sscanf(value, "%d", &Mega2BatchItems[it].value.option);
                }
                else if (!strcmp(keyword, "Value_Genetic_Distance_SexTypeMap")) {
                    it = /* 48 */ Value_Genetic_Distance_SexTypeMap;
                    if (READ_ITER1(it)) {
                        mult_decl(it);
                        continue;
                    }
                    Mega2BatchItems[it].item_read=1;
                    sscanf(value, "%d", &Mega2BatchItems[it].value.option);
                }
            } else if (!strncmp(keyword, "Error", (size_t) 5)) {
                it = -1;
                for (j = /* 19 */ Error_Loci ; j <= /* 23 */ Error_Probabilities ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }

                if (it == /* 19 */ Error_Loci  || it == /* 20 */ Error_Except_Loci ) {
                    /* Error loci or error_except_loci,
                       both need error_loci_num (item 21) to be
                       specified.
                    */
                    if (!Mega2BatchItems[/* 21 */ Error_Loci_Num].item_read) {
                        if (iter == 2) {
                            missing_dependency_keyword(it, 21);
                        }
                    } else {
                        if (Mega2BatchItems[/* 21 */ Error_Loci_Num].value.option <= 0) {
                            invalid_value_field(21);
                        }
                        Mega2BatchItems[it].item_read = 1;
                        Mega2BatchItems[it].value.mult_opts
                            = CALLOC((size_t)Mega2BatchItems[/* 21 */ Error_Loci_Num].value.option, int);
                        intstr = strtok(value, " ");
                        if (intstr == NULL) {
                            missing_items(/* 21 */ Error_Loci_Num, it);
                        }
                        Mega2BatchItems[it].value.mult_opts[0] = atoi(intstr);
                        for (j = 1; j < Mega2BatchItems[/* 21 */ Error_Loci_Num].value.option; j++) {
                            intstr = strtok(NULL, " ");
                            if (intstr == NULL) {
                                missing_items(/* 21 */ Error_Loci_Num, it);
                            }
                            Mega2BatchItems[it].value.mult_opts[j] = atoi(
                                intstr);
                        }
                    }
                } else if (it == /* 21 */ Error_Loci_Num ) {
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (it == /* 22 */ Error_Model ) {
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%c", &(Mega2BatchItems[it].value.copt));

                } else if (it == /* 23 */ Error_Probabilities ) {
                    Mega2BatchItems[/* 23 */ Error_Probabilities].item_read = 1;
                    if (!Mega2BatchItems[/* 22 */ Error_Model].item_read) {
                        if (iter == 2) {
                            missing_dependency_keyword(it, 22);
                        }
                    }

                    switch (Mega2BatchItems[/* 22 */ Error_Model].value.copt) {
                    case 'U':
                        /* read a single value and ignore the rest */
                        Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues
                            = CALLOC((size_t)1, double);
                        intstr = strtok(value, " ");
                        if (intstr == NULL) {
                            errorvf("Item %s has to have 5 values.\n",
                                    Mega2BatchItems[/* 22 */ Error_Model].keyword);
                            EXIT(BATCH_FILE_ITEM_ERROR);
                        }
                        Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[0]
                            = atof(intstr);
                        break;

                    case 'S':
                        /* read 5 values */
                        Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues
                            = CALLOC((size_t)5, double);
                        intstr = strtok(value, " ");
                        if (intstr == NULL) {
                            errorvf("Item %s has to have 5 values.\n",
                                    Mega2BatchItems[/* 23 */ Error_Probabilities].keyword);
                            EXIT(BATCH_FILE_ITEM_ERROR);
                        }
                        Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[0]
                            = atof(intstr);
                        for (j = 1; j < 5; j++) {
                            intstr = strtok(NULL, " ");
                            if (intstr == NULL) {
                                errorvf("Item %s has to have 5 values.\n",
                                        Mega2BatchItems[/* 23 */ Error_Probabilities].keyword);
                                EXIT(BATCH_FILE_ITEM_ERROR);
                            }
                            Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[j] = atof(
                                intstr);
                        }
                        break;

                    case 'M':
                        /* Ignore everything */
                        break;
                    default:
                        invalid_value_field(/*22 */ Error_Model);
                        break;
                    }
                }
            } else if (!strncmp(keyword, "Defau", (size_t) 5)) {
                it = -1;
                for (j = /* 25 */ Default_Outfile_Names ; j <= /* 29 */ Default_Rplot_Options ; j++) {
                    if (!strcmp(keyword, keywords[j])) {
                        it = j;
                        break;
                    }
                }
                if (READ_ITER1(it)) {
                    mult_decl(it);
                    continue;
                }

                Mega2BatchItems[it].item_read = 1;
                sscanf(value, "%c", &(Mega2BatchItems[it].value.copt));
                if (Mega2BatchItems[it].value.copt != 'y'
                    && Mega2BatchItems[it].value.copt != 'Y'
                    && Mega2BatchItems[it].value.copt != 'n'
                    && Mega2BatchItems[it].value.copt != 'N') {
                    invalid_value_field(it);
                    break;
                }
            } else if (!strncmp(keyword, "Xlinked", (size_t) 7)) {
                it = /* 30 */ Xlinked_Analysis_Mode ;
                Mega2BatchItems[it].item_read = 1;
                sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
            } else if (!strncmp(keyword, "Output_Path", (size_t) 11)) {
                it = /* 33 */ Output_Path ;
                Mega2BatchItems[it].item_read = 1;
                sscanf(value, "%s", Mega2BatchItems[it].value.name);
            } else if (!strncmp(keyword, "Count", (size_t) 5)) {
                if (!strcmp(keyword, "Count_Genotypes")) {
                    it = /* 34 */ Count_Genotypes ;
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                } else if (!strcmp(keyword, "Count_Halftyped")) {
                    it = /* 36 */ Count_Halftyped ;
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%c", &(Mega2BatchItems[it].value.copt));
                    if (Mega2BatchItems[it].value.copt != 'y'
                        && Mega2BatchItems[it].value.copt != 'Y'
                        && Mega2BatchItems[it].value.copt != 'n'
                        && Mega2BatchItems[it].value.copt != 'N') {
                        invalid_value_field(it);
                        break;
                    }
                } else {
                    it = /* 38 */ Count_HWE_genotypes ;
                    Mega2BatchItems[it].item_read = 1;
                    Mega2BatchItems[it].item_read = 1;
                    sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
                }
            } else if (!strncmp(keyword, "Rplot", (size_t) 5)) {
                /* for now only one item */
                it = /* 35 */ Rplot_Statistics ;
                Mega2BatchItems[it].item_read = 1;
                strcpy(Mega2BatchItems[it].value.name, value);
            } else if (!strncmp(keyword, "Allel", (size_t) 5)) {
                it = /* 37 */ AlleleFreq_SquaredDev ;
                Mega2BatchItems[it].item_read = 1;
                Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue = atof(strtok(value, " "));
            } else if (!strncmp(keyword, "Output_Map_Num", (size_t) 10)) {
                it = /* 41 */ Output_Map_Num;
                Mega2BatchItems[it].item_read = 1;
                sscanf(value, "%d", &(Mega2BatchItems[it].value.option));
            } else if (!strncasecmp(keyword, "PLINK", (size_t) 5)) {
                // Shouldn't there be a check here to see if the Analysis_Option == Plink?
                if (iter >= 2) continue;
                it = /* 43 */ PLINK_Args;
                Mega2BatchItems[it].item_read = 1;
/*              sscanf(value, "%s", Mega2BatchItems[it].value.name);*/
                strcpy(Mega2BatchItems[it].value.name, value);
                PLINK_args(value);
            } else if (!strcasecmp(keyword, "Loop_Over_Chromosomes")) {
                it = /* 50 */ Loop_Over_Chromosomes;
                Mega2BatchItems[it].item_read = 1;
                sscanf(value, "%c", &(Mega2BatchItems[it].value.copt));
                if (Mega2BatchItems[it].value.copt != 'y'
                    && Mega2BatchItems[it].value.copt != 'Y'
                    && Mega2BatchItems[it].value.copt != 'n'
                    && Mega2BatchItems[it].value.copt != 'N') {
                    invalid_value_field(it);
                    break;
                }
                // The 'analysis' variable can only be guaranted to hold data until the third pass.
                // This is because the sub-analysis option may not have been read till the end of the second pass.
                if (iter == 3) {
                    if (!(*analysis)->Loop_Over_Chromosomes_implemented()) {
                        warnvf("The batch file item 'Loop_Over_Chromosomes' is not implemented for analysis option '%s'.\n", analysis_name);
                    }
                }
            } else if (!strcasecmp(keyword, "Structure.PopDataPheno")) {
                if (iter == 3) {
                    if (*analysis == STRUCTURE) {
                        it = /* 51 */ Structure$PopDataPheno;
                        Mega2BatchItems[it].item_read = 1;
                        sscanf(value, "%s", Mega2BatchItems[it].value.name);
                    } else {
                        warnvf("The batch file item 'Structure.PopDataPheno' is not implemented for analysis option '%s'.\n", analysis_name);
                    }
                }
            }
            
            if (iter == 1 && it == -1) {
                warn_unknown(keyword);
            }
        } // if (strcmp(nextline, "") != 0) {
    } //  while (!feof(fp)) {

    // If the Value_Missing_Quant_On_Output is either not specified or was not permitted
    // to have been specified in the batch file (e.g., output_quant_can_define_missing_value()
    // returns false), then check to see if there is a default value for that analysis mode (e.g.,
    // output_quant_default_value() != NULL). If so, then make it look as if
    // Value_Missing_Quant_On_Output was read as that default.
    if (iter == 3) {
        if (!ITEM_READ(Value_Missing_Quant_On_Output) &&
            (*analysis)->output_quant_default_value() != (const char *)NULL) {
            Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read = 1;
            strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name,
                   (*analysis)->output_quant_default_value());
            
            // Since the user did not specify the default value, this is an internal check.
            if ((*analysis)->output_quant_must_be_numeric()) {
                char *end;
                const char *value = (*analysis)->output_quant_default_value();
                // The converted value is thrown away, because we just want to know if it is valid...
                (void) strtod(value, &end);
                if (strlen(value) == 0 || strlen(end) != 0 || errno == ERANGE) {
                    // Conversion of the entire string was not successful or some other error...
                    errorvf("INTERNAL ERROR:'Analysis_Option' = '%s' the method\n",
                            (*analysis)->_name);
                    errorf(" 'output_quant_default_value()' must return a numeric string.");
                    EXIT(DATA_INCONSISTENCY);
                }
            }
        }
    }

    fclose(fp);
}

void batchfile_process(char *batch_file_name, analysis_type *analysis)
{
	// this is where the batch items are read from the batch file, and checked for consistency...
    set_batch_items(batch_file_name, 1, analysis);
    set_batch_items(batch_file_name, 2, analysis);
	// Ahh... because analysis doesn't get defined until iteration 2 and so any logic that depends
	// on it's existance won't work till after that!!!!
    set_batch_items(batch_file_name, 3, analysis);

    check_batch_items();
}


//
// Write the item associated with the Mega2BatchItem to the batch file.
// It will be written based on the format (batch_item_value) for which it was defined (item_value_type).
// see batch_input.h
void batchf(int item)
{
    int index, *opt;
    double *fopt;
    char **str;
    FILE *batchfp;
    int first_time = ((access(Mega2Batch, F_OK) != 0)? 1: 0);

    /* Mega2Batch does not exist, first call to batchf */

    batch_item_type batch_item = Mega2BatchItems[item];
    if (first_time) {
        create_batchfile();
    }

    batchfp = fopen(Mega2Batch, "a");
    if (batchfp == NULL) {
        errorvf("Could not open Mega2 batch file (no write permission?).\n");
        EXIT(FILE_WRITE_ERROR);
    }

    if (first_time) {
        fprintf(batchfp, "%cVersion4.4\n", COMMENT_CHAR);
        fprintf(batchfp, "%c          ", COMMENT_CHAR);
        write_time(batchfp);
        batch_file_doc(batchfp);
    }

    fprintf(batchfp, "%s=", batch_item.keyword);
    switch(batch_item.value_type) {
    case INT:
        fprintf(batchfp, "%d", batch_item.value.option);
        break;
    case YORN:
        index = (batch_item.value.copt == 'y' ||
		 batch_item.value.copt == 'Y')? 1 : 0;
        fprintf(batchfp, "%s", yorn[index]);
        break;
    case STRING:
        fprintf(batchfp, "%s", batch_item.value.name);
        break;
    case FLOAT:
        if (fabs(batch_item.value.fvalue - QMISSING) <= EPSILON) {
            fprintf(batchfp, "NA");
        } else {
            fprintf(batchfp, "%f", batch_item.value.fvalue);
        }
        break;
    case CHAR:
        fprintf(batchfp, "%c", batch_item.value.copt);
        break;
    case INT_LIST:
        opt = batch_item.value.mult_opts;
        while (*opt != -99) {
            fprintf(batchfp, "%d ", *opt);
            opt++;
        }
        break;
    case FLOAT_LIST:
        fopt = batch_item.value.mult_fvalues;
        while ((int)(*fopt) != -9999) {
            fprintf(batchfp, "%g ", *fopt);
            fopt++;
        }
        break;

    case NAME_LIST:
        str = &(batch_item.value.mult_names[0]);
        while(*str  != NULL) {
            fprintf(batchfp, "%s ", *str);
            str++;
        }
        break;

    default:
        break;
    }
    fprintf(batchfp, "\n");
    fclose(batchfp);
}

#endif /* NEW_BATCH */

#ifdef CFREE

void Free_batch_items(void) {
    item_value_type vt;
    int i, j, num_trs;

    if (ITEM_READ(Traits_Loop_Over)) {
        if (LoopOverTrait)
            num_trs = parse_string(Mega2BatchItems[Traits_Loop_Over].value.name,
                                   NULL, &num_trs, (int)LARGE, 0);
        else
            num_trs = 1;
    } else if (ITEM_READ(Traits_Combine)) {
        num_trs = parse_string(Mega2BatchItems[Traits_Combine].value.name,
                               NULL, &num_trs, (int)LARGE, 0);
    } else {
        num_trs = 1;
    }
    for (i=0; i < NUM_KEYS; i++) {
        vt = Mega2BatchItems[i].value_type;
        if (vt == STRING) {
            free(Mega2BatchItems[i].value.name);
            continue;
        } else if (!ITEM_READ(i))
            continue;

        if (vt == INT_LIST) {
            free(Mega2BatchItems[i].value.mult_opts);
        } else if (vt == FLOAT_LIST)
            free(Mega2BatchItems[i].value.mult_fvalues);
        else if (vt == NAME_LIST) {
            for (j = 0; j < num_trs; j++) {
                /* mult_names may be copied into the output_path array */
                if (output_paths[j+1] != Mega2BatchItems[i].value.mult_names[j]) 
                    free(Mega2BatchItems[i].value.mult_names[j]);
            }
            free(Mega2BatchItems[i].value.mult_names);
        }
    }
    free(Mega2BatchItems);
}

#endif /* CFREE */
