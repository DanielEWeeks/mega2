/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "common.h"
#include "typedefs.h"
#include "analysis.h"

#include "batch_input.h"

#include "input_ops.hh"
#include "batch_input_ext.h"
#include "compress_ext.h"
#include "error_messages_ext.h"
#include "menu_value_missing_ext.h"
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

using namespace std;

#include <iostream>

#include "str_utils.hh"

typedef map<Cstr, batch_item_type *> Mapsb;
typedef map<Cstr, batch_item_type *>::const_iterator Mapsbp;
typedef std::list<batch_item_type *> Listb;
typedef std::list<batch_item_type *>::const_iterator Listbp;

Mapsb BatchItemMap;
Listb BatchItemList;

batch_item_type *Mega2BatchItems;

int batchINPUTFILES, batchANALYSIS, batchREORDER;
int batchTRAIT, batchAFFVALUE, batchERROR;

extern void prog_name_to_num(char *prog_name, analysis_type *analysis);


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
      45    Input_Aux_File
      46    Value_Genetic_Distance_Index
      47    Value_Base_Pair_Position_Index
      48    Value_Genetic_Distance_SexTypeMap
      49    Value_Missing_Quant_On_Output
      50    Loop_Over_Chromosomes
      51    Structure$PopDataPheno
      52    Value_Marker_Compression
      53    Input_Format_Type
      54    Input_Path
      55    Input_PLINK_Map_File
      56    VCF_Args
      57    VCF_INFO_Marker_Alternative_Key
      58    Value_Missing_Affect_On_Input
      59    Value_Missing_Affect_On_Output
      60    Output_File_Stem
      61    Value_Imputed_Threshold
*/

// NOTE: the default value is always represented as a string and converted if necessary.
static keyw_t keywords[] = {
    {"Input_Pedigree_File",                   STRING,     ""},
    {"Input_Locus_File",                      STRING,     ""},
    {"Input_Map_File",                        STRING,     ""},
    {"Input_Omit_File",                       STRING,     ""},
    {"Input_Untyped_Ped_Option",              INT,      "-1"},
    {"Analysis_Option",                       STRING,     ""},
    {"Analysis_Sub_Option",                   STRING,     ""},
    {"Chromosome_Single",                     INT,       "1"},
    {"Chromosomes_Multiple_Num",              INT,       "0"},
    {"Chromosomes_Multiple",                  CHRM_LIST,  ""},
    {"Loci_Selected_Num",                     INT,       "0"},
    {"Loci_Selected",                         LINE,       ""},
    {"Trait_Single",                          INT,       "0"},
    {"REMTraits_Num",                         INT,       "0"},
    {"Traits_Loop_Over",                      LINE,       ""},
    {"Traits_Combine",                        LINE,       ""},
    {"Trait_Subdirs",                         NAME_LIST,  ""},
    {"Value_Missing_Quant_On_Input",          FLOAT,   "0.0"}, 
    {"Value_Affecteds",                       NAME_LIST,  ""},
    {"Error_Loci",                            INT_LIST,   ""},
    {"Error_Except_Loci",                     INT_LIST,   ""},
    {"Error_Loci_Num",                        INT,       "0"},
    {"Error_Model",                           CHAR,      "X"},
    {"Error_Probabilities",                   FLOAT_LIST, ""},
    {"Input_Do_Error_Sim",                    YORN,      "n"},
    {"Default_Outfile_Names",                 YORN,      "n"},
    {"Default_Reset_Invalid",                 YORN,      "y"},
    {"Default_Other_Values",                  YORN,      "n"},
    {"Default_Ignore_Nonfatal",               YORN,      "n"},
    {"Default_Rplot_Options",                 YORN,      "n"},
    {"Xlinked_Analysis_Mode",                 INT,       "2"},
    {"REMCovariates_Num",                     INT,       "0"},
    {"Covariates_Selected",                   LINE,       ""},
    {"Output_Path",                           STRING,     ""},
    {"Count_Genotypes",                       INT,       "2"},
    {"Rplot_Statistics",                      LINE,       ""},
    {"Count_Halftyped",                       YORN,      "n"},
    {"AlleleFreq_SquaredDev",                 FLOAT,   "0.0"},
    {"Count_HWE_genotypes",                   INT,       "2"},
    {"Input_Frequency_File",                  STRING,     ""},
    {"Input_Penetrance_File",                 STRING,     ""},
    {"REMOutput_Map_Num",                     INT,       "1"},
    {"Value_Missing_Allele_Aff",              LINE,       ""},
    {"PLINK_Args",                            LINE,       ""},
    {"Input_Phenotype_File",                  STRING,     ""},
    {"Input_Aux_File",                        STRING,     ""},
    {"Value_Genetic_Distance_Index",          INT,      "-1"},
    {"Value_Base_Pair_Position_Index",        INT,      "-1"},
    {"Value_Genetic_Distance_SexTypeMap",     INT,      "-1"},
    {"Value_Missing_Quant_On_Output",         LINE,    "0.0"},
    {"Loop_Over_Chromosomes",                 YORN,      "n"},
    {"Structure.PopDataPheno",                STRING,     ""},
    {"Value_Marker_Compression",              INT,       "1"}, // MARKER_SCHEME_BITS
    {"Input_Format_Type",                     INT,       "0"},
    {"Input_Path",                            STRING,     ""},
    {"Input_PLINK_Map_File",                  STRING,     ""},
    {"VCF_Args",                              LINE,       ""},
    {"VCF_Marker_Alternative_INFO_Key",       LINE,       ""},
    {"Value_Missing_Affect_On_Input",         LINE,      "0"},
    {"Value_Missing_Affect_On_Output",        LINE,      "0"},
    {"Output_File_Stem",                      STRING,     ""},

    {"Default_Reset_Halftype",                YORN,      "y"},
    {"Default_Reset_Mendelerr",               YORN,      "y"},
    {"Default_Reset_Alleleerr",               YORN,      "y"},
    {"Default_Set_Uniq",                      YORN,      "y"},

    {"Imputed_Oxford_Single_Chr",             STRING,   "--"},
    {"Imputed_Info_Metric_Threshold",         FLOAT,   "0.3"},
    {"Imputed_Hard_Call_Threshold",           FLOAT,   "0.9"},
    {"Imputed_Genotype_Missing_Fraction",     FLOAT,   "0.1"},
    {"Imputed_Missing_Codes",                 NAME_LIST,  ""},
    {"Imputed_Allow_Duplicate_Markers",       YORN,      "n"},
    {"Imputed_Allow_Indels",                  YORN,      "n"},
    {"Imputed_RSID_Separator",                STRING,    ":"},
    {"Input_Imputed_Info_File",               STRING,     ""},

    {"ID_pedigree",                           INT,       "0"},
    {"ID_person",                             INT,       "0"},

    {"Shapeit_recomb_rdir",                   STRING,     ""},
    {"Shapeit_recomb_rfile",                  STRING,     "?"},
    {"Shapeit_file_stem",                     STRING,     ""},

    {"file_name_stem",                        STRING,     ""},
    {"additional_program_args",               LINE,       ""},
    {"RoadTrips_male_prevalence",             FLOAT,      "0.123"},
    {"RoadTrips_female_prevalence",           FLOAT,      "0.123"},

    {"minimac_reference_haplotype_file",      STRING,     "?"},
    {"batch_cpu_count",                       INT,        "1"},

    {"shapeit_reference_map_file",            STRING,     "?"},
    {"shapeit_reference_haplotype_file",      STRING,     "?"},
    {"shapeit_reference_legend_file",         STRING,     "?"},
    {"shapeit_reference_sample_file",         STRING,     ""},
    {"shapeit_haps_file_selected",            INT,        "0"},
    {"shapeit_reference_map_directory",       STRING,     "."},
    {"shapeit_reference_panel_directory",     STRING,     "."},
    {"minimac_reference_panel_directory",    STRING,     "."},

    {"DBfile_name",                           STRING,     "dbmega2.db"},

};

int NUM_KEYS = sizeof(keywords)  / sizeof (keyw_t);

pair<Cstr,Cstr> keyword_aliases[] = {
    make_pair("Value_Missing_Quant_On_Input", "Value_Missing_Quant"),
    make_pair("VCF_Args", "VCF_ARGS"),
    make_pair("PLINK_Args", "PLINK"),
    make_pair("Input_Aux_File", "Input_Binary_File"),
    make_pair("AlleleFreq_SquaredDev", "AlleleFreq_SquaredError"),
    make_pair("Count_Halftyped", "Count_Halftypes"),
    make_pair("REMOutput_Map_Num", "Output_Map_Num"),
    make_pair("Imputed_Allow_Duplicate_Markers", "Imputed_Allow_Duplicates"),
    make_pair("","")  //sentinel
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
        if (keywords[key].keyword.compare(0, (size_t) 3, "REM")) {
            fprintf(batchfp, "%c   %d) %s \n", COMMENT_CHAR,
                    key+1, C(keywords[key].keyword));
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
    extern int env;

//  Mega2BatchItems = CALLOC((size_t)NUM_KEYS, batch_item_type); // Only on OSX
    Mega2BatchItems = new batch_item_type[NUM_KEYS];
    for (i=0;  i < NUM_KEYS; i++) {
        batch_item_type *bi    = &Mega2BatchItems[i];
        item_value_type& type  = keywords[i].type;
        Str& deflt             = keywords[i].default_string;
        bi->keyword     = keywords[i].keyword;
        bi->value_type  = type;
        bi->item_number = i;
        bi->line_number = 0;
        bi->items_read  = 0;
        bi->flag        = Add2BatchItemList;

        Cstr& name      = bi->keyword;
        BatchItemMap[name]  = bi;

        if (env & 1) {
            char *xp = getenv(C(name));
            if (xp != NULL) {
                warnvf("parameter %s\n\treplacing default value (\"%s\") with value from environment (\"%s\")\n",
                       C(bi->keyword), C(deflt), xp);
                deflt = string(xp);
            }
        }
        switch(type) {
        case STRING:
        case LINE:
            bi->value.name         = CALLOC((size_t)FILENAME_LENGTH, char);
            strcpy(bi->value.name, deflt.c_str());
            break;

        default:
            bi->value.name         = CALLOC((size_t)FILENAME_LENGTH, char);
            strcpy(bi->value.name, deflt.c_str());
            break;

        case INT:
            bi->value.option       = atoi(deflt.c_str());
            break;

        case FLOAT:
            bi->value.fvalue       = atof(deflt.c_str());
            break;

        case CHAR:
            bi->value.copt         = deflt[0];
            break;

        case YORN:
            bi->value.copt         = (unsigned char)tolower(deflt[0]);
            break;

        case INT_LIST:
        case CHRM_LIST:
            bi->value.mult_opts     = NULL;
            break;

        case NAME_LIST:
            bi->value.mult_names   = NULL;
            break;

        case FLOAT_LIST:
            bi->value.mult_fvalues = NULL;
            break;
        }
    }

    // We won't list these options to BatchItemList since they need special treatment
    int err = 0;
    batch_item_type *bi;
    Cstr clear[] = {"Analysis_Option", "Analysis_Sub_Option", "Input_Format_Type"};
    for (i = 0; i < 3; i++) {
        if (map_get(BatchItemMap, clear[i], bi))
            bi->flag = Clear;
        else {
            errorvf("Internal Error!  %s not found %s\n:", C(clear[i]));
            err++;
        }
    }
    if (err) {
        EXIT(BATCH_FILE_ITEM_ERROR);
    }
    
    // Aliases
    for (i = 0; keyword_aliases[i].first != ""; i++) {
        if (map_get(BatchItemMap, keyword_aliases[i].first, bi))
	    BatchItemMap[keyword_aliases[i].second] = bi;
        else {
            errorvf("Internal Error!  %s not found %s\n:", C(keyword_aliases[i].first));
            err++;
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

static void missing_item_goto_menu(int batch_item, const char *menu_name)
{
    if (batch_item >= 0) {
        warnvf("Missing necessary batch item %s, going to %s.\n",
               C(Mega2BatchItems[batch_item].keyword), menu_name);
    } else {
        warnvf("Going to %s.\n", menu_name);
    }
}

void check_batch_items(void)
{
    int read_a_section;
    /* Input files */

    if (database_read)
        batchINPUTFILES=1;

    else if (ITEM_READ(/* 0 */ Input_Pedigree_File) &&
// #1
        (ITEM_READ(/* 1 */ Input_Locus_File) || Input->req_locus_file == 0) &&
// #2
        (ITEM_READ(/* 2 */ Input_Map_File)   || Input->req_map_file == 0)  &&
// #3
        ITEM_READ(/* 4 */ Input_Untyped_Ped_Option)) {
        // NOTE: You don't need a Map_File if you are working with a VCF file because it contains one map (VCF.p)...
        batchINPUTFILES=1;
        mssgf("Input filenames and missing value indicator read in from batch file.");
    } else {
        if (!ITEM_READ(/* 0 */ Input_Pedigree_File)) {
            missing_item_goto_menu(0, "Input menu");
        } else if (!ITEM_READ(/* 1 */ Input_Locus_File) && 
                   (Mega2BatchItems[Input_Format_Type].items_read && Mega2BatchItems[Input_Format_Type].value.option < in_format_binary_PED)) {
            missing_item_goto_menu(1, "Input menu");
        } else if (!ITEM_READ(/* 2 */ Input_Map_File)) {
            missing_item_goto_menu(2, "Input menu");
//pp is this necessary
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
        mssgf("Chromosome(s) and markers read in from batch file.");
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
      (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].items_read ? 1 : 0) +
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
            C(Mega2BatchItems[batch_item].keyword),
            Mega2BatchItems[batch_item].value.option);
    errorvf("Item %s has less than %d items.\n",
            C(Mega2BatchItems[list_item].keyword),
            Mega2BatchItems[batch_item].value.option);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void missing_dependency_keyword(int dependent_item, int dependent_on_item)
{
    errorvf("Keyword %s is defined but keyword %s is missing from batch file.\n",
            C(Mega2BatchItems[dependent_item].keyword),
            C(Mega2BatchItems[dependent_on_item].keyword));
    errorvf("Keyword %s depends on keyword %s.\n",
            C(Mega2BatchItems[dependent_item].keyword),
            C(Mega2BatchItems[dependent_on_item].keyword));
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
                C(Mega2BatchItems[dependent_item].keyword));
    } else {
        sprintf(err_msg,
                "Keyword %s depends on one of these keywords to be defined:",
                C(Mega2BatchItems[dependent_item].keyword));
    }

    for (i=0; i < num_dependent_items; i++) {
        strcat(err_msg, C(Mega2BatchItems[dependent_on_items[i]].keyword));
    }
    errorf(err_msg);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void mult_decl(int item)
{
    errorvf("Found multiple occurrences of keyword %s,\n",
            C(Mega2BatchItems[item].keyword));
//  errorf("Using the first definition, and ignoring all later definitions.");
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void malformed_batch_line(const char *keyword)
{
    errorvf("Malformed line for keyword %s.\n", keyword);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

static void check_dependencies(analysis_type *analysis)
{
    batch_item_type& bi8 = Mega2BatchItems[Chromosomes_Multiple_Num];
    batch_item_type& bi9 = Mega2BatchItems[Chromosomes_Multiple];
    if (bi9.items_read) {
        if (!bi8.items_read)
            missing_dependency_keyword(9, 8);
        else if (bi8.value.option <= 0 || bi8.value.option != bi9.items_read-1)
            invalid_value_field(8);
    }

    batch_item_type& bi14 = Mega2BatchItems[Traits_Loop_Over];
    batch_item_type& bi16 = Mega2BatchItems[Trait_Subdirs];
    if (bi16.items_read) {
        if (!bi14.items_read)
            missing_dependency_keyword(16, 14);
    }
    batch_item_type& bi12 = Mega2BatchItems[Trait_Single];
    batch_item_type& bi15 = Mega2BatchItems[Traits_Combine];
    batch_item_type& bi18 = Mega2BatchItems[Value_Affecteds];
    if (bi18.items_read) {
        int num_trs = 1;
        if (bi12.items_read)
            num_trs = 1;
        else if (bi14.items_read)
            parse_string(bi14.value.name, NULL, &num_trs, (int) LARGE, 0);
        else if (bi15.items_read)
            parse_string(bi15.value.name, NULL, &num_trs, (int) LARGE, 0);
        else {
            int opt_array[3] = { 12, 14, 15 };
            /* error */
            missing_mult_dependency_keyword(18, opt_array, 3, 0);
        }
        if (num_trs != bi18.items_read-1) {
            invalid_value_field(18);
        }
    }

    batch_item_type& bi19 = Mega2BatchItems[Error_Loci];
    batch_item_type& bi20 = Mega2BatchItems[Error_Except_Loci];
    batch_item_type& bi21 = Mega2BatchItems[Error_Loci_Num];
    if (bi19.items_read) {
        if (!bi21.items_read)
            missing_dependency_keyword(19, 21);
        else if (bi21.value.option <= 0)
            invalid_value_field(21);

        if (bi21.value.option != bi19.items_read-1) {
            missing_items(/* 21 */ Error_Loci_Num, 19);
        }
    } else if (bi20.items_read) {
        if (!bi21.items_read)
            missing_dependency_keyword(20, 21);
        else if (bi21.value.option <= 0)
            invalid_value_field(21);

        if (bi21.value.option != bi20.items_read-1) {
            missing_items(/* 21 */ Error_Loci_Num, 20);
        }
    }

    batch_item_type& bi22 = Mega2BatchItems[Error_Model];
    batch_item_type& bi23 = Mega2BatchItems[Error_Probabilities];
    if (bi23.items_read) {
        if (!bi22.items_read) 
            missing_dependency_keyword(23, 22);

        switch (bi22.value.copt) {
        case 'U':
            if (bi23.items_read != 1 + 1) {
                errorvf("Item %s has to have 1 value.\n", C(bi22.keyword));
		EXIT(BATCH_FILE_ITEM_ERROR);
            }
            break;
        case 'S':
            if (bi23.items_read != 5 + 1) {
                errorvf("Item %s has to have 5 values.\n", C(bi22.keyword));
		EXIT(BATCH_FILE_ITEM_ERROR);
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

    batch_item_type& bi50 = Mega2BatchItems[Loop_Over_Chromosomes];
    if (bi50.items_read) {
        if (!(*analysis)->Loop_Over_Chromosomes_implemented()) {
//            warnvf("The batch file item 'Loop_Over_Chromosomes' is not implemented for analysis option '%s/%s'.\n",
//                   (*analysis)->_name, (*analysis)->_subname);
            warnvf("The batch file item 'Loop_Over_Chromosomes' is not implemented for analysis option '%s'.\n",
                   (*analysis)->_name);
        }
    }

}

/* global variables */
static char analysis_name[50]="", sub_analysis_name[100]="";

void RawBatchValueSet(char *value, batch_item_type *bi)
{
    char *intstr;
    int capacity = 8;
    int j;

    bi->items_read = 1;

    switch(bi->value_type) {
    case STRING:
        sscanf(value, "%s", bi->value.name);
        break;

    case LINE:
        strcpy(bi->value.name, value);
        break;

    default:
        sscanf(value, "%s", bi->value.name);
        break;

    case INT:
        sscanf(value, "%d", &(bi->value.option));
        break;

//err: why not sscan %f
    case FLOAT:
        bi->value.fvalue = atof(strtok(value, " "));
        break;

    case CHAR:
        sscanf(value, "%c", &(bi->value.copt));
        break;

    case YORN:
        sscanf(value, "%c", &(bi->value.copt));
        if (bi->value.copt != 'y' && bi->value.copt != 'Y' &&
            bi->value.copt != 'n' && bi->value.copt != 'N')
            invalid_value_field(bi->item_number);
        break;

    case INT_LIST:
//      asm("int $3");
        warnvf("INT_LIST: %s\n", value);
        bi->value.mult_opts = CALLOC((size_t)capacity, int);
        intstr = NULL;
        intstr = strtok(value, " \n");
        for (j = 0; intstr != NULL; j++) {
            if (j >= capacity) {
                capacity *= 2;
                bi->value.mult_opts = REALLOC(bi->value.mult_opts, (size_t)capacity, int);
            }
            bi->value.mult_opts[j] = atoi(intstr);
            intstr = strtok(NULL, " \n");
        }
        bi->value.mult_opts = REALLOC(bi->value.mult_opts, (size_t)j, int);
        bi->items_read += j;
        break;

    case NAME_LIST:
//      tested
        capacity = 8;
        bi->value.mult_names = CALLOC((size_t)capacity, char *);
        intstr = NULL;
        intstr = strtok(value, " \n");
        for (j = 0; intstr != NULL; j++) {
            if (j >= capacity) {
                capacity *= 2;
                bi->value.mult_names = REALLOC(bi->value.mult_names, (size_t)capacity, char *);
            }
            bi->value.mult_names[j] = CALLOC((size_t)FILENAME_LENGTH, char);
            strcpy(bi->value.mult_names[j], intstr);
            intstr = strtok(NULL, " \n");
        }
        bi->value.mult_names = REALLOC(bi->value.mult_names, (size_t)j, char *);
        bi->items_read += j;
        break;

    case FLOAT_LIST:
//      asm("int $3");
        warnvf("FLOAT_LIST: %s\n", value);
        bi->value.mult_fvalues = CALLOC((size_t)capacity, double);
        intstr = NULL;
        intstr = strtok(value, " \n");
        for (j = 0; intstr != NULL; j++) {
            if (j >= capacity) {
                capacity *= 2;
                bi->value.mult_fvalues = REALLOC(bi->value.mult_fvalues, (size_t)capacity, double);
            }
            bi->value.mult_fvalues[j] = atof(intstr);
            intstr = strtok(NULL, " \n");
        }
        bi->value.mult_fvalues = REALLOC(bi->value.mult_fvalues, (size_t)j, double);
        bi->items_read += j;
        break;

    case CHRM_LIST:
//      tested
        bi->value.mult_opts = CALLOC((size_t)capacity, int);
        int ch;
        intstr = NULL;
        intstr = strtok(value, " \n");
        for (j = 0; intstr != NULL; j++) {
            if (j >= capacity) {
                capacity *= 2;
                bi->value.mult_opts = REALLOC(bi->value.mult_opts, (size_t)capacity, int);
            }
            if ((ch=STR_CHR(intstr)) != -1) {
                bi->value.mult_opts[j] = ch;
            } else {
                errorf("For batch file item 'Chromosomes_Multiple'.");
                errorvf("A chromosome must be a positive integer less then %d or one of the\n", lastautosome+4);
                errorf("following: X, Y, XY, MT, U.");
                EXIT(BATCH_FILE_ITEM_ERROR);
            }
            intstr = strtok(NULL, " \n");
        }
        bi->value.mult_opts = REALLOC(bi->value.mult_opts, (size_t)j, int);
        bi->items_read += j;
        break;
    }
}

class Mega2BatchItemException : public std::exception {
public:
    Str exc;
    Mega2BatchItemException(Cstr& s): exc(s) {};
    virtual
//  ~Mega2BatchItemException() _NOEXCEPT throw () {};
    ~Mega2BatchItemException() throw () {};
    const char* what() const throw() { return exc.c_str(); };
};

void RawBatchValueSet(char *value, int i)
{
    if (i >= NUM_KEYS) {
        warnvf("Mega2BatchItemSet: index %d too large\n:", i);
        return;
    }
    batch_item_type *bi = &Mega2BatchItems[i];
    RawBatchValueSet(value, bi);
}

void RawBatchValueSet(char *value, Cstr& key)
{
    batch_item_type *bi = NULL;
    if (map_get(BatchItemMap, key, bi))
	RawBatchValueSet(value, bi);
    else {
        warnvf("Mega2BatchItemSet: key %s not found\n:", C(key));
//      throw std::runtime_error(Str("Mega2BatchItemSet: bad key"));
        throw Mega2BatchItemException("Mega2BatchItemSet: bad key");
    }

}

void BatchValueSet(Veci& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    bi->items_read = 1;

    int capacity = vec.size()+1;
    int j;

    switch(bi->value_type) {

    case CHRM_LIST:
    case INT_LIST:
        bi->value.mult_opts = CALLOC((size_t)capacity, int);
        for (j = 0; j < capacity - 1; j++) {
            bi->value.mult_opts[j] = vec[j];
        }
        bi->value.mult_opts[j] = 0;
        bi->items_read += j;
        break;
    default:
        break;
    }
}

void BatchValueSet(Vecs& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    bi->items_read = 1;

    int capacity = vec.size()+1;
    int j;

    switch(bi->value_type) {

    case NAME_LIST:
//      tested
        bi->value.mult_names = CALLOC((size_t)capacity, char *);
        for (j = 0; j < capacity - 1; j++) {
            bi->value.mult_names[j] = CALLOC((size_t)vec[j].size()+1, char);
            strcpy(bi->value.mult_names[j], C(vec[j]));
        }
        bi->value.mult_names[j] = CALLOC((size_t)1, char);
        bi->value.mult_names[j] = 0;
        bi->items_read += j;
        break;
    default:
        break;
    }
}

void BatchValueSet(Vecd& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    bi->items_read = 1;

    int capacity = vec.size()+1;
    int j;

    switch(bi->value_type) {

    case FLOAT_LIST:
        bi->value.mult_fvalues = CALLOC((size_t)capacity, double);
        for (j = 0; j < capacity - 1; j++) {
            bi->value.mult_fvalues[j] = vec[j];
        }
        bi->value.mult_fvalues[j] = 0.0;
        bi->items_read += j;
        break;
    default:
        break;
    }
}

/*
void BatchValueSet(Veci& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    bi->items_read = 1;

    int capacity = vec.size()+1;
    int j;

    switch(bi->value_type) {

    case CHRM_LIST:
//      tested
        bi->value.mult_opts = CALLOC((size_t)capacity, int);
        int ch;
        for (j = 0; j < capacity - 1; j++) {
            bi->value.mult_names[j] = CALLOC((size_t)vec[j].size()+1, char);
            sprintf(bi->value.mult_names[j], "%s", vec[j]);
        }
        bi->value.mult_names[j] = CALLOC((size_t)1, char);
        bi->value.mult_names[j] = 0;
        bi->items_read += j;
        break;
    default:
        break;
    }
}
*/

batch_item_type *BatchItemGet(int i)
{
    if (i >= NUM_KEYS) {
        warnvf("Mega2BatchItemGet: index %d too large\n:", i);
        throw Mega2BatchItemException("Mega2BatchItemGet: index too large");
//      return (batch_item_type *) NULL;
    }
    return &Mega2BatchItems[i];
}

batch_item_type *BatchItemGet(Cstr& key)
{
    batch_item_type *bi = NULL;
    if (map_get(BatchItemMap, key, bi))
        return bi;
    else {
        warnvf("Mega2BatchItemGet: key \"%s\" not found\n:", C(key));
        throw Mega2BatchItemException("Mega2BatchItemGet: bad key");
//      return (batch_item_type *) NULL;
    }
}

void BatchValueGet(Vecs& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    vec.clear();
    if (bi->items_read <= 0) return;

    int capacity = bi->items_read - 1;
    int j = 0;

    switch(bi->value_type) {

    case NAME_LIST:
//      tested
        vec.insert(vec.begin(), &bi->value.mult_names[j], &bi->value.mult_names[capacity]);
        break;
    default:
        break;
    }
}

void BatchValueGet(Vecd& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    vec.clear();
    if (bi->items_read <= 0) return;

    int capacity = bi->items_read - 1;
    int j = 0;

    switch(bi->value_type) {

    case FLOAT_LIST:
        vec.insert(vec.begin(), &bi->value.mult_fvalues[j], &bi->value.mult_fvalues[capacity]);
        break;
    default:
        break;
    }
}

void BatchValueGet(Veci& vec, Cstr &item)   // for lists
{
    batch_item_type *bi = BatchItemGet(item);
    vec.clear();
    if (bi->items_read <= 0) return;

    int capacity = bi->items_read - 1;
    int j = 0;

    switch(bi->value_type) {

    case CHRM_LIST:
    case INT_LIST:
        vec.insert(vec.begin(), &bi->value.mult_opts[j], &bi->value.mult_opts[capacity]);
        break;
    default:
        break;
    }
}

static void parse_batch_file(char *batch_file_name, analysis_type *analysis)
{
    char nextline[FILENAME_LENGTH];
    char value[FILENAME_LENGTH];
    int errtok = 0, errline = 0, cnt = 0;
    extern int debug;
    extern int env;

    /* First the the verseion number */
    FILE *fp = fopen(batch_file_name, "r");
    if (fp == (FILE *)NULL) {
        errorvf("could not open %s for reading!\n", batch_file_name);
        EXIT(FILE_READ_ERROR);
    }

#ifdef later
//  int version4 = 0;
    (void)fgets(nextline, FILENAME_LENGTH - 1, fp);
    if (!strncasecmp(nextline, "#Version4.4", (size_t) 11)) {
//      version4 = 1;
    }
    rewind(fp);
#endif

    int line_n = 1;
    int lz;
    Token token;
    Str lhs, rhs;
    Mapsbp lookup;
    batch_item_type *bi = (batch_item_type *) NULL;

    while (!feof(fp)) {
        get_line(fp, nextline);
        line_n++;
        if (*nextline == 0) continue;
        lz = strlen(nextline) - 1;
        if (nextline[lz] == '\n') nextline[lz] = 0;
        if (*nextline == 0) continue;

        // if the line is not the empty string, though this would match a line with one space on it...
        // better would be to search for the leack of a comment character in the first column and
        // an equal sign somewhere in the line....
        token.set(nextline);
        if (token.assign(lhs, rhs)) {
//
            if (debug) warnvf("lhs: |%s| rhs: |%s|\n", C(lhs), C(rhs));
            if (map_get(BatchItemMap, lhs, bi)) {
                if (bi->items_read > 0) {
                    errorvf("Duplicate %s option at lines %d and %d in batch file.\n",
                            C(lhs), bi->line_number, line_n);
                    mult_decl(bi->item_number);
                }
                bi->line_number = line_n;
                bi->items_read = 1;
                bi->value_str = rhs;
                if (bi->flag == Add2BatchItemList)
                    BatchItemList.push_back(bi);
                cnt++;
            } else {
                errorvf("Unexpected batch file option: %s\n", C(lhs));
                errtok++;
            }
        } else {
            errorvf("Odd format batch file line %s\n:", nextline);
            errline++;
        }
    }
    fclose(fp);

    int err = 0;

    // First Select Input_Format_Type
    if (map_get(BatchItemMap, "Input_Format_Type", bi)) {
        if (bi->items_read)
            sscanf(bi->value_str.c_str(), "%d", &(bi->value.option));
        else
            bi->value.option = 0;
        Input_Format = (INPUT_FORMAT_t)bi->value.option;
        if (debug) warnvf("Option: %s %s\n", "Input_Format_Type", bi->value_str.c_str());
    } else {
        errorvf("%s option not found.\n", "Input_Format_Type");
        err++;
    }

    // Then Select Analysis Mode
    if (map_get(BatchItemMap, "Analysis_Option", bi)) {
//xx        if (bi->items_read == 0 && database_dump == 0)
        if (bi->items_read == 0)
        {
            errorvf("%s option not set\n", "Analysis_Option");
            err++;
        } else {
//xx
/*
            if (bi->items_read == 1)
                strcpy(analysis_name, bi->value_str.c_str());
            else
                strcpy(analysis_name, "dump");
*/
            strcpy(analysis_name, bi->value_str.c_str());
            strcpy(bi->value.name, analysis_name);

            // This is where the 'analysis' variable get's assiged too...
            prog_name_to_num(analysis_name, analysis);
            if (*analysis == NULL || !(*analysis)->is_enabled()) {
                err++;
                unknown_prog(analysis_name);
            }
        }
        if (debug) warnvf("Option: %s %s\n", "Analysis_Option", bi->value_str.c_str());
    } else {
        errorvf("%s option not found.\n", "Analysis_Option");
        err++;
    }
    /* Also check if this analysis has a sub-option and if sub-option was
       provided in the batch file */

    if ((*analysis) != NULL && (*analysis)->has_sub_options()) {
        if (map_get(BatchItemMap, "Analysis_Sub_Option", bi)) {
            if (bi->items_read == 0) {
                warnvf("Trying empty sub-program name (no name specified).\n",
                       analysis_name);
            } else {
                strcpy(sub_analysis_name, bi->value_str.c_str());
                strcpy(bi->value.name, sub_analysis_name);
            }
            (*analysis)->sub_prog_name_to_sub_option(sub_analysis_name, analysis);
            if (*analysis == NULL) {
                errorvf("%s sub-program needs to be specified.\n", analysis_name);
                errorvf("Required keyword %s missing from batch file.\n",
                        C(Mega2BatchItems[Analysis_Sub_Option].keyword));
                err++;
            }
            if (debug) warnvf("Option: %s %s\n", "Analysis_Sub_Option", bi->value_str.c_str());
        } else {
            errorvf("%s option not found.\n", "Analysis_Sub_Option");
            err++;
        }
    }

    // gross hack: needs better plan
    pair<Cstr,char *> hacks[] = {
	make_pair("Value_Missing_Quant_On_Input", Str_Missing_Quant_On_Input),
	make_pair("Value_Missing_Quant_On_Output", Str_Missing_Quant_On_Output),
	make_pair("Value_Missing_Affect_On_Input", Str_Missing_Affect_on_Input),
	make_pair("Value_Missing_Affect_On_Output", Str_Missing_Affect_On_Output),
    };
    for (int i = 0; i < 4; i++) {
        if (map_get(BatchItemMap, hacks[i].first, bi)) {
	    if (bi->items_read > 0) {
		strcpy(hacks[i].second, bi->value_str.c_str());
                if (i < 2)
                    bi->value_str = "0.0";  // reset to default value
                else
                    bi->value_str = "0";    // reset to default value
	    }
	} else {
            errorvf("%s option not found %s\n:", C(hacks[i].first));
		err++;
        }
    }
    if (err) EXIT(BATCH_FILE_ITEM_ERROR);

    // Finally process all the other batch arguments
    for (Listbp BatchItemListp = BatchItemList.begin(); BatchItemListp != BatchItemList.end(); BatchItemListp++) {
        batch_item_type *bibatch = (*BatchItemListp);
        Str& keyword = bibatch->keyword;
        Str& bivalue = bibatch->value_str;
        strcpy(value, bibatch->value_str.c_str());
        if (debug) msgvf("key %s, val %s\n", C(keyword), value);

        if (env & 2) {
            char *xp = getenv(C(keyword));
            if (xp != NULL) {
                warnvf("parameter %s\n\treplacing batch file value (\"%s\") with value from environment (\"%s\")\n",
                       C(keyword), value, xp);
                strcpy(value, xp);
            }
        }
//err: does this do anything
        if (bivalue == "" && !keyword.compare(0, 14, "value_missing_") ) {
            *value = 0;
        } else if (bivalue == "") {
            malformed_batch_line(C(keyword));
        }

        RawBatchValueSet(value, bibatch);
    }

    if (errtok + errline) {
        errorvf("\nBatch file Status: %d arguments, %d bad lines, %d unknown options\n", cnt, errline, errtok);
        EXIT(BATCH_FILE_ITEM_ERROR);
    }
}

void batchfile_process(char *batch_file_name, analysis_type *analysis)
{
	// this is where the batch items are read from the batch file, and checked for consistency...
    parse_batch_file(batch_file_name, analysis);
    check_dependencies(analysis);

    int input_set = 0;
    int xcf = 0;
    if (Mega2BatchItems[/* 53 */ Input_Format_Type].items_read) {
        Input_Format = (INPUT_FORMAT_t) Mega2BatchItems[/* 53 */ Input_Format_Type].value.option;

        input_set = 1;
        if (Input_Format == in_format_binary_VCF || Input_Format == in_format_compressed_VCF ||
            Input_Format == in_format_VCF)
            xcf = 1;

    }
    if (Mega2BatchItems[/* 43 */ PLINK_Args].items_read) {
        PLINK_args(Mega2BatchItems[PLINK_Args].value.name, xcf);
        if (!Input_Format) {
            input_set = 1;
            if (PLINK.plink == binary_PED_format)
                Input_Format = in_format_binary_PED;
            else if (PLINK.plink == PED_format)
                Input_Format = in_format_PED;
        }
    }
    if (!input_set)
        Input_Format = in_format_traditional;

    extern Input_Base *createinput(INPUT_FORMAT in_format);
    Input = createinput(Input_Format);

    check_batch_items();
}

//
// Write the item associated with the Mega2BatchItem to the batch file.
// It will be written based on the format (batch_item_value) for which it was defined (item_value_type).
// see batch_input.h
// NOTE: this always outputs an item to BATCH FILE
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
        fprintf(batchfp, "%c          %s",
                COMMENT_CHAR, write_time());
        batch_file_doc(batchfp);
    }

    fprintf(batchfp, "%s=", C(batch_item.keyword));
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
    case LINE:
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
    case CHRM_LIST:
//err: i = batch_item.items_read
        opt = batch_item.value.mult_opts;
        while (*opt != -99) {
            fprintf(batchfp, "%d ", *opt);
            opt++;
        }
        break;
    case FLOAT_LIST:
//err:
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

// NOTE: these two MAY output an item to BATCH FILE
void batchf(Cstr& keyword)
{
    batch_item_type *bi = NULL;
    if (map_get(BatchItemMap, keyword, bi)) {
        if (bi->items_read)
            batchf(bi->item_number);
    } else {
        warnvf("batchf: key %s not found %s\n:", C(keyword));
        return;
    }
}

void batchf(batch_item_type *bi) {
    if (bi->items_read)
        batchf(bi->item_number);
}

#endif /* NEW_BATCH */

void Free_batch_items(void) 
{
    item_value_type vt;
    int i, j, num_trs;

//err: use num_trs = items_read - 1
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
        } else if (vt == CHRM_LIST) {
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
    delete[] Mega2BatchItems;
}
