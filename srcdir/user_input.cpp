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

#include "batch_input_ext.h"
#include "cw_routines_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "output_file_names_ext.h"
#include "plink_ext.h"
#include "utils_ext.h"
/*
        batch_input_ext.h:  batchf
        cw_routines_ext.h:  strsep
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  makedir prog_name sub_prog_name
        grow_string_ext.h:  grow
              plink_ext.h:  PLINK_args
              utils_ext.h:  EXIT chomp draw_line log_line
*/



/* prototypes */
int             ReOrderMenu(int num_chromo, int *chromsomes, int *selection);
int             analysis_menu1(analysis_type  *analysis);

int             menu1(file_format *infl_type,
		      char **pedfl_name, char **locusfl_name,
		      char **mapfl_name, char **omitfl_name,
		      char **freqfl_name, char **penfl_name,
                      char **bedfl_name, char **phefl_name,
		      int *Untyped_ped_opt, int *Error_sim_opt,
		      char **output_path, double *freq_mismatch_thresh);

void            define_affection_labels(linkage_ped_top *Top,
					analysis_type analysis);
void           set_missing_quant_input(linkage_ped_top *Top, analysis_type analysis);

int            gh_cov_selection(int num_select, int *trait_order,
				int *covariates,
				linkage_locus_top *LTop);
char           *untyped_ped_messg(int opt, char *messg);

void           get_genetic_distance_index(ext_linkage_locus_top *EXLTop);
void           get_base_pair_position_index(ext_linkage_locus_top *EXLTop);

static  void define_labels(linkage_ped_top *Top, int tr, char *affdata_str);

static void     untyped_ped_menu(int *opt);

/*===========prototypes============= */

extern struct PLINK PLINK;


// see R_output.h & R_output.c write_plink_map_file...
int genetic_distance_index;
//cpk genetic_distance_map_type genetic_distance_sex_type_map; // see common.h
int genetic_distance_sex_type_map; // see common.h
int base_pair_position_index;

/*===========================================================================*/
/* menu for user loci selection for the user's selection/reordering option
   an attempt to make the entry crash proof is made here */
int             ReOrderMenu(int num_chromo, int *chromosomes, int *selection)
{

    char select[10];
    int i;
    char stars[3];
    int myselect = -1;

    for (i = 0; i < 3; i++) stars[i] = ' ';
    stars[*selection - 1] = '*';

    draw_line();
    if (num_chromo == 0) {
        while (myselect < 0 || myselect > 2) {
            printf("Locus Reordering Menu\n");
            printf("0) Done with this menu - please proceed.\n");
            printf("%c1) Select all loci.\n", stars[0]);
            printf("%c2) Select by locus number.\n", stars[1]);
            printf("Enter 0-2 > ");
            fcmap(stdin, "%s", select);
            newline;
            draw_line();
            sscanf(select, "%d", &myselect);
        }
        *selection = myselect;
    } else {
        char prchr[4];
        while (myselect < 0 || myselect > 3) {
            printf("Locus Reordering Menu\n");
            printf("0) Done with this menu - please proceed.\n");
            CHR_STR(chromosomes[0], prchr);
            printf("%c1) Select all loci in map order on chromosome: %s\n",
                   stars[0], prchr);
            printf("%c2) Select by locus number.\n", stars[1]);
            printf("%c3) Select loci on multiple chromosomes[", stars[2]);
#if 1
            for (i = 0; i < num_chromo; i++)
            {
                if (Ignore_Unmapped(chromosomes[i])) continue;

                CHR_STR(chromosomes[i], prchr);
                printf("%s ", prchr);
            }
#else
            for (i = 0; i < NumChromo; i++)
            {
                if (Ignore_Unmapped(i)) continue;

                CHR_STR(global_chromo_entries[i], prchr);
                printf("%s ", prchr);
            }
#endif
            printf("]\n");
            printf("Enter 0 - 3 > ");
            fcmap(stdin, "%s", select);
            newline;
            draw_line();
            sscanf(select, "%d", &myselect);
        }
    }
    *selection = myselect;
    return *selection;
}

int invalid_analysis(int opt)
{ 
    if (opt == 6 || opt == 7 || opt == 15) {
        switch(opt) {
        case 6:
            errorf("APM option is disabled");
            break;
        case 7:
            errorf("APM-Mult option is disabled");
            break;
        case 15:
            errorf("TDT-Max program is disabled.");
            break;
        }
        return 1;
    } else
        return 0;
}

extern analysis_types analysis_list[];
extern int count_analysis_list;

int             analysis_menu1(analysis_type  *analysis)
{
    char            choice[10], sub_prog[20];
    int             i, choice_=0;
    int             num_per_screen;
    char            option_col1[35], option_col2[35];
/* When adding a new output option:

   The function analysis_menu1() displays the current list of output
   choices within Mega2.
   a) Add the full name for your program to the constant array
   option_name[][]. Keep it short, but descriptive
   enough. The length of this is currently set to 29 characters,
   so that the list can be displayed in 2 columns.
*/

    //
    // The selection from the suboption generally is used to create a new option.
    // In some instances this does not make sense, like in the instance where you
    // are just changing the type of file that you are writing, and everything else
    // about the option is the same.

    if (batchANALYSIS) {
/*
   analysis was set already in batch_input.c: set_batch_items()
        choice_ = Mega2BatchItems[/ * 5 * / Analysis_Option].value.option;
        choice_--;
 */
    } else {
        num_per_screen = HALF(count_analysis_list);

        while (choice_ < 1 || choice_ > (count_analysis_list) || !analysis_list[choice_-1].analysis->is_enabled())  {
            printf("                 ANALYSIS MENU \n");
            draw_line();
            for (i = 0; i < num_per_screen; i++) {
                sprintf(option_col1, "%2d %s", i+1, analysis_list[i].option_name);
                if ((num_per_screen +i) < (count_analysis_list)) {
                    sprintf(option_col2, "%2d %s",
                            num_per_screen+i+1, analysis_list[num_per_screen+i].option_name);
                } else {
                    strcpy(option_col2, "");
                }
                printf("%-34s %-34s\n", option_col1, option_col2);
            }

            printf("\nSelect an option between 1-%2d > ", count_analysis_list);
            fcmap(stdin, "%s", choice); newline;
            choice_ = 0;
            sscanf(choice, "%d", &choice_);
            if (choice_ >= 1 && choice_ <= (count_analysis_list)) {
                if (analysis_list[choice_-1].analysis->is_enabled())
                    printf("    You have selected: %s\n", analysis_list[choice_-1].option_name);
                else
                    printf("Disabled choice \"%s\". Please select again.\n",
                           analysis_list[choice_-1].option_name);
            } else
                printf("Choice out of range. Please select again from options 1-%2d.\n",
                    count_analysis_list);
        }
        *analysis = analysis_list[choice_-1].analysis;
        /* note the analysis may be changed */
        if ((*analysis)->has_sub_options())
            (*analysis)->interactive_sub_prog_name_to_sub_option(analysis);
    }

    AllowUnmapped = (ALLOW_NO_CHR(*analysis) ? 1 : 0);

    (*analysis)->prog_name(ProgName);
    mssgvf("Analysis option: %s.\n", ProgName);

    if (InputMode == INTERACTIVE_INPUTMODE) {
        // Write the Analysis_Option, and Analysis_Sub_Option to the batch file ...
        strcpy(Mega2BatchItems[/* 5 */ Analysis_Option].value.name, ProgName);
        batchf(Analysis_Option);
        if  ((*analysis)->_suboption > 0) {
            /* We already know the sub-option as well, so we can map and store both */
            (*analysis)->sub_prog_name((*analysis)->_suboption, &(sub_prog[0]));
            strcpy(Mega2BatchItems[/* 6 */ Analysis_Sub_Option].value.name, sub_prog);
	    // Write out the Mega2BatchItem given. The value written will be the '.value.option'...
            batchf(Analysis_Sub_Option);
            grow(ProgName, "-%s", sub_prog);
        }
    } else {
        if ((*analysis)->is_sub_option() || (*analysis)->has_sub_options()) {

            (*analysis)->sub_prog_name((*analysis)->_suboption, &(sub_prog[0]));

            grow(ProgName, "-%s", sub_prog);
            mssgvf("Analysis sub-option: %s.\n", sub_prog);
        }
    }

    return 1;
}   /* end of analysis_menu */

static void     thresh_string(double val, char *str)

{
    if (val >= LARGE) {
        sprintf(str, "No limit");
    } else {
        sprintf(str, "%g", val);
    }
    return;

}

static void gen_files(int plinkf, char *extension, char **locusfl, char **pedfl, char **mapfl, char **omitfl, char **freqfl, char **penfl, char **bedfl, char **phefl)
{
    if (plinkf) {
        sprintf(*locusfl, "%s.datain", extension);

        if (access(*locusfl, F_OK)) {
            sprintf(*locusfl, "%s.names", extension);
        }

        if (PLINK.plink == binary_PED_format) {
            sprintf(*pedfl,  "%s.fam", extension);
            sprintf(*mapfl,  "%s.bim", extension);
            sprintf(*bedfl,  "%s.bed", extension);
        } else if (PLINK.plink == PED_format) {
            sprintf(*pedfl,  "%s.ped", extension);
            sprintf(*mapfl,  "%s.map", extension);
        }
        sprintf(*phefl, "%s.phe", extension);
        sprintf(*omitfl, "%s.omit", extension);
        sprintf(*freqfl, "%s.frequency", extension);
        sprintf(*penfl,  "%s.penetrance", extension);
    } else {
        sprintf(*locusfl, "datain.%s", extension);

        if (access(*locusfl, F_OK)) {
            sprintf(*locusfl, "names.%s", extension);
        }

        sprintf(*pedfl,  "pedin.%s", extension);
        sprintf(*mapfl,  "map.%s", extension);
        sprintf(*omitfl, "omit.%s", extension);
        sprintf(*freqfl, "frequency.%s", extension);
        sprintf(*penfl,  "penetrance.%s", extension);
    }
}


static void missing_optional_keyword(int batch_item, const char *message)
{
#ifndef HIDESTATUS
    mssgvf("Keyword %s not in batch file, %s.\n",
            Mega2BatchItems[batch_item].keyword, message);
#endif
}


static void missing_mandatory_keyword(int batch_item)
{
    errorvf("Required keyword %s missing from batch file.\n",
            Mega2BatchItems[batch_item].keyword);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

int             menu1(file_format *infl_type,
		      char **pedfl_name, char **locusfl_name,
		      char **mapfl_name, char **omitfl_name,
		      char **freqfl_name, char **penfl_name,
                      char **bedfl_name, char **phefl_name,
		      int *Untyped_ped_opt, int *Error_sim_opt,
		      char **output_path, double *freq_mismatch_thresh)
{
    int             i, j, choice_ = -1;
    char            cchoice[10], extension_name[FILENAME_LENGTH], float_str[11];
    char            messg[100];
    char            oldmapname[FILENAME_LENGTH], oldomitname[FILENAME_LENGTH];
    char            oldpedname[FILENAME_LENGTH], oldlocusname[FILENAME_LENGTH];
    char            oldfreqname[FILENAME_LENGTH], oldpenname[FILENAME_LENGTH];
    char            oldbedname[FILENAME_LENGTH], oldphenname[FILENAME_LENGTH];
    int             exit_loop=0;
    boolean         LocusDataFileSpecified = false, PedigreeDataFileSpecified = false;
    boolean         MapFileNameSpecified = false, FreqFileSpecified = false;
    boolean         PenFileSpecified = false,  OmitFileSpecified = false;
    boolean         BinaryFileSpecified = false, PhenotypeFileSpecified = false;
    const int      file_batch_items[] = {Input_Pedigree_File,
					 Input_Locus_File,
					 Input_Map_File,
					 Input_Omit_File,
					 Input_Frequency_File,
					 Input_Penetrance_File,
					 Input_Phenotype_File,
					 Input_Binary_File};
    int            file_batch_items_size = sizeof(file_batch_items) / sizeof(int);
    int            ext_i=1, loc_i=2, ped_i=3, map_i=4, omit_i = 5, freq_i=6, pen_i=7;
    int            out_i=8, err_i=9, untyp_i=10, thresh_i=11, miss_i=12;
    int            plink_i = 14, plink_j = 15, plink_k = 16;
    int            plink_args_i = 17, plink_phe_i = 18, plink_bed_i = 19;
    int            plinkf = 0, idx, choiceA[20]; /* idx should be 1+ largest <>_i value (above)*/
    char           PLINKArgs[FILENAME_LENGTH] = "";

    /* specifying chromosome or extension overrides previous specifications
       unless user specified each file separately */
    *Untyped_ped_opt=2; /* Exclude any pedigree with 1 or less untyped people */
    *Error_sim_opt = 0;
    *infl_type = LINKAGE;
    *freq_mismatch_thresh = LARGE;

    if (*output_path == NULL)
        *output_path = CALLOC((size_t) FILENAME_LENGTH, char);
    if (*locusfl_name == NULL)
        *locusfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*pedfl_name == NULL)
        *pedfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*mapfl_name == NULL)
        *mapfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*omitfl_name == NULL)
        *omitfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*freqfl_name == NULL)
        *freqfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*penfl_name == NULL)
        *penfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*bedfl_name == NULL)
        *bedfl_name =  CALLOC((size_t) FILENAME_LENGTH, char);
    if (*phefl_name == NULL)
        *phefl_name =  CALLOC((size_t) FILENAME_LENGTH, char);

    if (batchINPUTFILES) {
        /* necessary files ped, locus and map
           or ped, name and map
        */
        if (Mega2BatchItems[/* 43 */ PLINK_Args].item_read == 1) {
            plinkf = 1;
        }
        for(j=0; j < file_batch_items_size; j++) {
            i = file_batch_items[j];
            if (Mega2BatchItems[i].item_read==1) {
                if (access(Mega2BatchItems[i].value.name, F_OK) == 0) {
                    switch(i) {
                    case Input_Pedigree_File:
                        strcpy(*pedfl_name,Mega2BatchItems[i].value.name);
                        break;
                    case Input_Locus_File:
                        strcpy(*locusfl_name,Mega2BatchItems[i].value.name);
                        break;
                    case Input_Map_File:
                        strcpy(*mapfl_name, Mega2BatchItems[i].value.name);
                        break;
                    case Input_Omit_File:
                        strcpy(*omitfl_name, Mega2BatchItems[i].value.name);
                        break;
                    case Input_Frequency_File:
                        strcpy(*freqfl_name, Mega2BatchItems[i].value.name);
                        break;
                    case Input_Penetrance_File:
                        strcpy(*penfl_name, Mega2BatchItems[i].value.name);
                        break;
                    case Input_Phenotype_File:
                        strcpy(*phefl_name, Mega2BatchItems[i].value.name);
                        break;
                    case Input_Binary_File:
                        if (PLINK.plink == binary_PED_format)
                            strcpy(*bedfl_name, Mega2BatchItems[i].value.name);
                        break;
                    }
                } else {
                    errorvf("Could not find file or path %s named by keyword %s.\n",
                            Mega2BatchItems[i].value.name,
                            Mega2BatchItems[i].keyword);
                    EXIT(FILE_NOT_FOUND);
                }
            } else {
                // These were not read from the batchfile...
                if (i < 3) {
                    if (i == Input_Locus_File && plinkf) continue;
                    missing_mandatory_keyword(i);
                } else if (i == Input_Binary_File && PLINK.plink == binary_PED_format) {
                    missing_mandatory_keyword(i);
                } else {
                    switch(i) {
                    case /* 3 */ Input_Omit_File:
                        missing_optional_keyword(Input_Omit_File, "omit file assumed to be unspecified");
                        if (*omitfl_name != NULL) {
                            free(*omitfl_name);
                            *omitfl_name = NULL;
                        }
                        break;

                    case /* 39 */ Input_Frequency_File:
                        missing_optional_keyword(Input_Frequency_File, "frequency file assumed to be unspecified");
                        if (*freqfl_name != NULL) {
                            free(*freqfl_name);
                            *freqfl_name = NULL;
                        }
                        break;

                    case /* 40 */ Input_Penetrance_File:
                        missing_optional_keyword(Input_Penetrance_File, "penetrance file assumed to be unspecified");
                        if (*penfl_name != NULL) {
                            free(*penfl_name);
                            *penfl_name = NULL;
                        }
                        break;

                    case /* 44 */ Input_Phenotype_File:
                        if (plinkf)
                            missing_optional_keyword(Input_Phenotype_File, "phenotype file assumed to be unspecified");
                        if (*phefl_name != NULL) {
                            free(*phefl_name);
                            *phefl_name = NULL;
                        }
                        break;

                    case /* 45 */ Input_Binary_File:
                        if (*bedfl_name != NULL) {
                            free(*bedfl_name);
                            *bedfl_name = NULL;
                        }
                        break;
                    }
                }
            }
        }

        /* untyped ped option */

        if (Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].item_read  == 1) {
            if (Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option >= 0) {
                *Untyped_ped_opt=Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option;
            } else {
                invalid_value_field(Input_Untyped_Ped_Option);
            }
        } else {
            missing_mandatory_keyword(Input_Untyped_Ped_Option);
        }

        /* Error sim option */

        if (Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].item_read  == 1) {
            if (tolower((unsigned char)Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt) == 'y' ||
                tolower((unsigned char)Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt) == 'n') {
                *Error_sim_opt= Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt;
            } else {
                invalid_value_field(Input_Do_Error_Sim);
            }
        }
        *Error_sim_opt =
            ((tolower(*Error_sim_opt) == 'y')? 1:0);

        /*     if (Mega2BatchItems[/ * 7 * / Chromosome_Single].item_read == 1) { */
        /*       *numchr = Mega2BatchItems[/ * 7 * / Chromosome_Single].value.option; */
        /*     } */

        if (Mega2BatchItems[/* 33 */ Output_Path].item_read == 1) {
            if (access(Mega2BatchItems[/* 33 */ Output_Path].value.name, W_OK) == 0 &&
                is_dir(Mega2BatchItems[/* 33 */ Output_Path].value.name)) {
               strcpy(*output_path, Mega2BatchItems[/* 33 */ Output_Path].value.name);
            } else {
                errorvf("file path %s named by keyword %s is not a writable directory.\n",
                        Mega2BatchItems[/* 33 */ Output_Path].value.name,
                        Mega2BatchItems[/* 33 */ Output_Path].keyword);
                EXIT(FILE_NOT_FOUND);
            }

        } else {
            missing_optional_keyword(Output_Path,  "using default '.' (current directory)");
            strcpy(*output_path, ".");
        }

        /* Threshold for comparing input and observed allele frequencies */

        if (Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].item_read == 1) {
            if (Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue < 0) {
                invalid_value_field(AlleleFreq_SquaredDev);
            }

            *freq_mismatch_thresh = Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue;
        } else {
            missing_optional_keyword(AlleleFreq_SquaredDev, "using default = 'no limit'");
        }
#ifdef USER_UNKNOWN
        if (Mega2BatchItems[/* 42 */ Value_Missing_Allele].item_read == 1) {
            strcpy(REC_UNKNOWN, Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name);
        }
#endif
        if (plinkf) {
            free(*locusfl_name);
            *locusfl_name = NULL;
        }

        return plinkf;
    }

    if (plinkf)
        strcpy(extension_name, "plink");
    else
        strcpy(extension_name, "01");

    gen_files(plinkf, extension_name, locusfl_name, pedfl_name, mapfl_name, omitfl_name, freqfl_name, penfl_name, bedfl_name, phefl_name);

    sprintf(*output_path, ".");

    while (!exit_loop) {
        printf("              Mega2 %s %s input menu:\n", Mega2Version, plinkf ? "PLINK" : "Traditional");
        draw_line();
        printf(" 0) Done with this menu - please proceed\n");
        idx=1;
        if (plinkf) {
            printf("%2d) Enter PLINK parameters:                 %s\n", idx, PLINKArgs);
            choiceA[idx++] = plink_args_i;
        }
        if (plinkf)
            printf("%2d) Input file stem:                        %s\n", idx,
                   extension_name);
        else
            printf("%2d) Input file extension:                   %s\n", idx,
                   extension_name);
        choiceA[idx++] = ext_i;

        if (! plinkf) {
            printf("%2d) Locus datafile:      (%s)        %s\n", idx,
                   "Mega2 loc",
                   (access(*locusfl_name, F_OK)? "_": *locusfl_name));
            choiceA[idx++] = loc_i;
        }

        printf("%2d) Pedigree datafile:   (%s)        %s\n", idx,
               (PLINK.plink == not_plink_format ? "Mega2 ped" : (PLINK.plink == binary_PED_format ? "PLINK fam" : "PLINK ped")),
               (access(*pedfl_name, F_OK)? "_": *pedfl_name));
        choiceA[idx++] = ped_i;

        printf("%2d) Map datafile:        (%s)        %s\n", idx,
               (PLINK.plink == not_plink_format ? "Mega2 map" : (PLINK.plink == binary_PED_format ? "PLINK bim" : "PLINK map")),
               (access(*mapfl_name, F_OK)? "_": *mapfl_name));
        choiceA[idx++] = map_i;

        if (plinkf) {
            if (PLINK.plink == binary_PED_format) {
                printf("%2d) Binary datafile:     (%s)        %s\n", idx,
                       "PLINK bed",
                       (access(*bedfl_name, F_OK)? "_": *bedfl_name));
                choiceA[idx++] = plink_bed_i;
            }

            printf("%2d) Phenotype datafile:  {%s}        %s\n", idx,
                   "PLINK phe",
                   (access(*phefl_name, F_OK)? "_": *phefl_name));
            choiceA[idx++] = plink_phe_i;
        }

        printf("%2d) Omit datafile:       {Mega2 omit}       %s\n", idx,
               (access(*omitfl_name, F_OK)? "_": *omitfl_name));
        choiceA[idx++] = omit_i;

        printf("%2d) Frequency datafile:  {annotated fmt}    %s\n", idx,
               (access(*freqfl_name, F_OK)? "_": *freqfl_name));
        choiceA[idx++] = freq_i;

        printf("%2d) Penetrance datafile: {annotated fmt}    %s\n", idx,
               (access(*penfl_name, F_OK)? "_": *penfl_name));
        choiceA[idx++] = pen_i;

        printf("%2d) Directory for writing output:           %s\n", idx,
               ((!strcmp(*output_path, "."))?
                "[ Current directory ]" : *output_path));
        choiceA[idx++] = out_i;

        printf("%2d) Simulate genotyping errors:             [%s]\n", idx,
               yorn[*Error_sim_opt]);
        choiceA[idx++] = err_i;

        printf("%2d) %s (based on fully_typed genotypes)\n", idx,
               untyped_ped_messg(*Untyped_ped_opt, messg));
        thresh_string(*freq_mismatch_thresh, float_str);
        choiceA[idx++] = untyp_i;

        printf("%2d) Warn if allele frequency error measure exceeds: %s\n", idx,
               float_str);
        choiceA[idx++] = thresh_i;

        if (plinkf) {
            printf("%2d) Switch to standard input menu\n", idx);
            choiceA[idx++] = plink_i;
            if (PLINK.plink == PED_format) {
                printf("%2d) Switch to PLINK input menu (binary format)\n", idx);
                choiceA[idx++] = plink_j;
            } else if (PLINK.plink == binary_PED_format) {
                printf("%2d) Switch to PLINK input menu (ped format)\n", idx);
                choiceA[idx++] = plink_k;
            }
        } else {
            printf("%2d) Switch to PLINK input menu (binary format)\n", idx);
            choiceA[idx++] = plink_j;
            printf("%2d) Switch to PLINK input menu (ped format)\n", idx);
            choiceA[idx++] = plink_k;
        }

#ifdef USER_UNKNOWN
        printf("%2d) Unknown allele & affection value(annotated only): %s\n",
               miss_i, REC_UNKNOWN);
        printf("Select from options 0-%d (%d to toggle) > ", idx-1, err_i);
#else
        printf("Select from options 0-%d > ", idx-1);
#endif
        while (1) {
            fcmap(stdin, "%s", cchoice); newline;
            if (!strcmp(cchoice, "q")) {
                printf("Exiting Mega2.\n");
                EXIT(1);
            }
            choice_ = -1;
            sscanf(cchoice, "%d", &choice_);
            if (choice_ < idx) {
                if (choice_ > 0) choice_ = choiceA[choice_];
                break;
            } else
                printf("Select from options 0-%d > ", idx-1);
        }

        if (choice_ ==  0) {
            exit_loop=1;
            if (! plinkf) {
                if (access(*locusfl_name, F_OK) != 0)   {
                    printf("ERROR: You must specify a locus datafile.\n");
                    exit_loop=0;
                }
                // Since we are not using PLINK format we don't need the .bed file...
                if (*bedfl_name != (char *)NULL ) {
                    free(*bedfl_name);
                    *bedfl_name = NULL;
                }
            }
            if (access(*pedfl_name, F_OK) != 0) {
                printf("ERROR: You must specify a pedigree datafile.\n");
                exit_loop=0;
            }
            if (access(*mapfl_name, F_OK) != 0) {
                printf("ERROR: You did not specify a map datafile.\n");
                exit_loop=0;
            }
            if (plinkf && PLINK.plink == not_plink_format) {
                printf("ERROR: You did not specify any PLINK parameters.\n");
                exit_loop=0;
            }
            if (PLINK.plink == binary_PED_format && access(*bedfl_name, F_OK) != 0) {
                printf("ERROR: You did not specify a PLINK binary datafile.\n");
                exit_loop=0;
            }
            if (PLINK.plink && PLINK.no_pheno == 0 && PLINK.trait[0] == 0) {
                printf("ERROR: You did not specify the pedigree file trait name.\n");
                printf("ERROR: Please choose option 1 and do so now.\n");
                exit_loop=0;
            }
            if (exit_loop == 0)
                draw_line();
            else {
                if (access(*omitfl_name, F_OK) != 0) {
                    free(*omitfl_name);
                    *omitfl_name = NULL;
                } else {
                    printf("NOTE: Marker untyping will take place according to the omit file.\n");
                    printf(
                        "      Please check the '%s' log file after MEGA2 is finished.\n",
                        Mega2Log);
                }

                if (access(*freqfl_name, F_OK) != 0) {
                    free(*freqfl_name);
                    *freqfl_name = NULL;
                }
                if (access(*penfl_name, F_OK) != 0) {
                    free(*penfl_name);
                    *penfl_name = NULL;
                }
                if (access(*phefl_name, F_OK) != 0) {
                    free(*phefl_name);
                    *phefl_name = NULL;
                }
                if (access(*bedfl_name, F_OK) != 0) {
                    free(*bedfl_name);
                    *bedfl_name = NULL;
                }
            }
        } else  if (choice_ == ext_i) {
            if (plinkf)
                printf("Please enter file stem > ");
            else
                printf("Please enter file extension > ");
            fcmap(stdin, "%s", extension_name); newline;
            if (strcasecmp(extension_name, "clear")==0) {
                strcpy(extension_name, "-");
            } else if (plinkf) {
                if (!LocusDataFileSpecified) {
                    sprintf(*locusfl_name, "%s.datain", extension_name);
                }
                if (!PedigreeDataFileSpecified) {
                    if (PLINK.plink == binary_PED_format)
                        sprintf(*pedfl_name,  "%s.fam", extension_name);
                    else if (PLINK.plink == PED_format)
                        sprintf(*pedfl_name,  "%s.ped", extension_name);
                }
                if (!MapFileNameSpecified) {
                    if (PLINK.plink == binary_PED_format)
                        sprintf(*mapfl_name,  "%s.bim", extension_name);
                    else if (PLINK.plink == PED_format)
                        sprintf(*mapfl_name,  "%s.map", extension_name);
                }
                if (PLINK.plink == binary_PED_format && !BinaryFileSpecified)
                    sprintf(*bedfl_name, "%s.bed", extension_name);
                if (!PhenotypeFileSpecified)
                    sprintf(*phefl_name, "%s.phe", extension_name);
                if (!OmitFileSpecified)
                    sprintf(*omitfl_name, "%s.omit", extension_name);
                if (!FreqFileSpecified)
                    sprintf(*freqfl_name, "%s.frequency", extension_name);
                if (!PenFileSpecified)
                    sprintf(*penfl_name,  "%s.penetrance", extension_name);
            } else {
                if (!LocusDataFileSpecified) {
                    sprintf(*locusfl_name, "datain.%s", extension_name);
                }
                if (access(*locusfl_name, F_OK)) {
                    sprintf(*locusfl_name, "names.%s", extension_name);
                }
                if (!PedigreeDataFileSpecified)
                    sprintf(*pedfl_name, "pedin.%s", extension_name);
                if (!MapFileNameSpecified)
                    sprintf(*mapfl_name, "map.%s", extension_name);
                if (!OmitFileSpecified)
                    sprintf(*omitfl_name, "omit.%s", extension_name);
                if (!FreqFileSpecified)
                    sprintf(*freqfl_name, "frequency.%s", extension_name);
                if (!PenFileSpecified)
                    sprintf(*penfl_name, "penetrance.%s", extension_name);
            }
        } else if (choice_ == loc_i) {
            draw_line();
            printf("Please enter locus file name ('clear' to clear) > ");
            strcpy(oldlocusname, *locusfl_name);
            fcmap(stdin, "%s", *locusfl_name); newline;

            if (strcasecmp(*locusfl_name, "clear") == 0) {
                strcpy(*locusfl_name, "-");
            } else {
                if (access(*locusfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *locusfl_name);
                    strcpy(*locusfl_name, oldlocusname);
                    LocusDataFileSpecified = false;
                } else {
                    LocusDataFileSpecified = true;
                }
            }
        } else if (choice_ == ped_i) {
            draw_line();
            printf("Please enter pedigree file name ('clear' to clear) > ");
            strcpy(oldpedname, *pedfl_name);
            fcmap(stdin, "%s", *pedfl_name); newline;

            if (strcasecmp(*pedfl_name, "clear") == 0) {
                strcpy(*pedfl_name, "-");
            } else {
                if (access(*pedfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *pedfl_name);
                    strcpy(*pedfl_name, oldpedname);
                    PedigreeDataFileSpecified = false;
                } else {
                    PedigreeDataFileSpecified = true;
                }
            }
        } else if (choice_ == map_i) {
            draw_line();
            printf("Please enter map file name ('clear' to clear) > ");
            strcpy(oldmapname, *mapfl_name);
            fcmap(stdin, "%s", *mapfl_name); newline;
            if (strcasecmp(*mapfl_name, "clear") == 0) {
                strcpy(*mapfl_name, "-");
            } else {
                if (access(*mapfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *mapfl_name);
                    MapFileNameSpecified = false;
                    strcpy(*mapfl_name, oldmapname);
                }
                else
                    MapFileNameSpecified = true;
            }
        } else if (PLINK.plink == binary_PED_format && choice_ == plink_bed_i) {
            draw_line();
            printf("Please enter binary file name ('clear' to clear) > ");
            strcpy(oldbedname, *bedfl_name);
            fcmap(stdin, "%s", *bedfl_name); newline;
            if (strcasecmp(*bedfl_name, "clear") == 0) {
                strcpy(*bedfl_name, "-");
            } else {
                if (access(*bedfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *bedfl_name);
                    BinaryFileSpecified = false;
                    strcpy(*bedfl_name, oldbedname);
                } else {
                    BinaryFileSpecified = true;
                }
            }
        } else if (choice_ == plink_phe_i) {
            draw_line();
            printf("Please enter phenotype file name ('clear' to clear) > ");
            strcpy(oldphenname, *phefl_name);
            fcmap(stdin, "%s", *phefl_name); newline;
            if (strcasecmp(*phefl_name, "clear") == 0) {
                strcpy(*phefl_name, "-");
            } else {
                if (access(*phefl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *phefl_name);
                    PhenotypeFileSpecified = false;
                    strcpy(*phefl_name, oldphenname);
                } else {
                    PhenotypeFileSpecified = true;
                }
            }
        } else if (choice_ == omit_i) {
            draw_line();
            printf("Please enter omit file name ('clear' to clear) > ");
            strcpy(oldomitname, *omitfl_name);
            fcmap(stdin, "%s", *omitfl_name); newline;
            if (strcasecmp(*omitfl_name, "clear") == 0) {
                strcpy(*omitfl_name, "-");
            } else {
                if (access(*omitfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *omitfl_name);
                    OmitFileSpecified = false;
                    strcpy(*omitfl_name, oldomitname);
                } else {
                    OmitFileSpecified = true;
                }
            }
        } else if (choice_ == freq_i) {   /* The 'Freq' datafile */
            draw_line();
            printf("Please enter frequency file name ('clear' to clear) > ");
            strcpy(oldfreqname, *freqfl_name);
            fcmap(stdin, "%s", *freqfl_name); newline;
            if (strcasecmp(*freqfl_name, "clear") == 0) {
                strcpy(*freqfl_name, "-");
            } else {
                if (access(*freqfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *freqfl_name);
                    strcpy(*freqfl_name, oldfreqname);
                    FreqFileSpecified = false;
                }
                else
                    FreqFileSpecified = true;
            }
        } else if (choice_ == pen_i) {   /* The penetrance datafile */
            draw_line();
            printf("Please enter penetrance file name ('clear' to clear) > ");
            strcpy(oldpenname, *penfl_name);
            fcmap(stdin, "%s", *penfl_name); newline;
            if (strcasecmp(*penfl_name, "clear") == 0) {
                strcpy(*penfl_name, "-");
            } else {
                if (access(*penfl_name, F_OK)) {
                    printf("WARNING: Could not find file %s\n", *penfl_name);
                    strcpy(*penfl_name, oldpenname);
                    PenFileSpecified = false;
                }
                else
                    PenFileSpecified = true;
            }
        } else if (choice_ == out_i) {   /* The output directory */
            draw_line();
            printf("Please enter output directory name > ");
            fcmap(stdin, "%s", *output_path); newline;

            if (access(*output_path, F_OK)) {
                char y[100];
                printf("WARNING: Could not find directory %s\n", *output_path);
                printf("Create this new directory \n");
                printf("   %s\n", *output_path);
                printf("   Yes or No (y/n)[default no] > ");
                fflush(stdout);
                (void)fgets(y, 9, stdin); newline; fflush(stdin);
                if (y[0] == 'Y' || y[0] == 'y') {
                    makedir(*output_path);
                } else {
                    strcpy(*output_path, ".");
                }
            } else if (! is_dir(*output_path)) {
                    printf("WARNING: %s is not a directory.\n", *output_path);
                    printf("Please specify a new or valid directory.\n");
                    strcpy(*output_path, ".");
            } else if (access(*output_path, W_OK)) {
                    printf("WARNING: %s is not a writable directory.\n", *output_path);
                    printf("Please specify a new or valid directory.\n");
            }

        } else if (choice_ == err_i) {
            *Error_sim_opt = TOGGLE(*Error_sim_opt);
        } else if (choice_ == untyp_i) {
            draw_line();
            untyped_ped_menu(Untyped_ped_opt);
        } else if (choice_ == thresh_i) {
            draw_line();
            printf("Please enter threshold value > ");
            fcmap(stdin, "%g", freq_mismatch_thresh); newline;
        } else if (choice_ == miss_i) {
            draw_line();
            printf("Please enter missing value indicator > ");
            fcmap(stdin, "%s", REC_UNKNOWN); newline;
        } else if (choice_ == plink_i || choice_ == plink_j || choice_ == plink_k) {
            if (choice_ == plink_i) {
                plinkf = 0;
                PLINK_clr(not_plink_format);
                strcpy(extension_name, "01");
            } else if (choice_ == plink_j) {
                plinkf = 1;
                PLINK_clr(binary_PED_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(extension_name, "plink");
            } else if (choice_ == plink_k) {
                plinkf = 1;
                PLINK_clr(PED_format);
/*             free(*bedfl_name); *bedfl_name = NULL; */
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(extension_name, "plink");
            }
            gen_files(plinkf, extension_name, locusfl_name, pedfl_name, mapfl_name, omitfl_name, freqfl_name, penfl_name, bedfl_name, phefl_name);
           LocusDataFileSpecified = PedigreeDataFileSpecified = 0;
           MapFileNameSpecified = OmitFileSpecified = 0;
           FreqFileSpecified = PenFileSpecified = 0;
           BinaryFileSpecified = PhenotypeFileSpecified = 0;
        } else if (choice_ == plink_args_i) {
            PLINK_usage();
            while (1) {
                fflush(stdout);
                (void)fgets(PLINKArgs, sizeof(PLINKArgs)-1, stdin); newline;
                i = (int)strlen(PLINKArgs);
                if (PLINKArgs[i-1] == '\n') PLINKArgs[i-1] = 0;
                if (PLINKArgs[i-1] == '\r') PLINKArgs[i-1] = 0;
                if (PLINK_args(PLINKArgs)) break;
                printf("Enter     UPDATED PLINK parameters:  %s\n", PLINKArgs);
            }
        } else {
            printf("Invalid option %s, select from options 0-%d.\n", cchoice, idx-1);
        }

        if (choice_) draw_line();
    }

    if (InputMode == INTERACTIVE_INPUTMODE) {
        strcpy(Mega2BatchItems[/* 0 */ Input_Pedigree_File].value.name, *pedfl_name);
        batchf(Input_Pedigree_File);
        if (! plinkf) {
            strcpy(Mega2BatchItems[/* 1 */ Input_Locus_File].value.name, *locusfl_name);
            batchf(Input_Locus_File);
        }
        strcpy(Mega2BatchItems[/* 2 */ Input_Map_File].value.name, *mapfl_name);
        batchf(Input_Map_File);
        if (PLINK.plink == binary_PED_format) {
            strcpy(Mega2BatchItems[/* 45 */ Input_Binary_File].value.name, *bedfl_name);
            batchf(Input_Binary_File);
        }
        if (*omitfl_name  != NULL) {
            strcpy(Mega2BatchItems[/* 3 */ Input_Omit_File].value.name, *omitfl_name);
            batchf(Input_Omit_File);
        }
        if (*freqfl_name  != NULL) {
            strcpy(Mega2BatchItems[/* 39 */ Input_Frequency_File].value.name, *freqfl_name);
            batchf(Input_Frequency_File);
        }
        if (*penfl_name  != NULL) {
            strcpy(Mega2BatchItems[/* 40 */ Input_Penetrance_File].value.name, *penfl_name);
            batchf(Input_Penetrance_File);
        }
        if (*phefl_name  != NULL) {
            strcpy(Mega2BatchItems[/* 44 */ Input_Phenotype_File].value.name, *phefl_name);
            batchf(Input_Phenotype_File);
        }
        Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option= *Untyped_ped_opt;
        batchf(Input_Untyped_Ped_Option);
        Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt = yorn[*Error_sim_opt][0];
        batchf(Input_Do_Error_Sim);
        strcpy(Mega2BatchItems[/* 33 */ Output_Path].value.name, *output_path);
        batchf(Output_Path);
        Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue=*freq_mismatch_thresh;
        batchf(AlleleFreq_SquaredDev);
#ifdef USER_UNKNOWN
        strcpy(Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name, REC_UNKNOWN);
        batchf(Value_Missing_Allele);
#endif
        strcpy(Mega2BatchItems[/* 43 */ PLINK_Args].value.name, PLINKArgs);
        if (plinkf) batchf(PLINK_Args);
    }

    if (plinkf) {
        free(*locusfl_name);
        *locusfl_name = NULL;
    }

    return plinkf;
}



/* static void default_labels(char *msg, int *liability, int *status, int num_classes) */

/* { */
/*   int i; */
/*   if (msg != NULL) */
/*     printf("%s", msg); */
/*   for (i=0; i< num_classes; i++) { */
/*     liability[i]=i+1; */
/*     status[i]=2; */
/*   } */
/*   return; */
/* } */

/* This is the routine to parse the affection label strings, and create the
   affected phenotype list Labels. We can have at most 2*num_classes labels.
   Set NumLabels = num_classes*2, and set undefined labels to -1?
*/

static int check_affdata_str(char *affdata_str, int num_classes)
{
    int no_error=1;
    char *status_str, *aff_copy, *aff_ptr;
    int liability, status;

    /* first see if there is a hyphen */

    aff_copy = strdup(affdata_str);
    aff_ptr = strtok(aff_copy, ",\n");

    while (1) {
        if (aff_ptr == NULL) {
            errorf("Empty affection label.\n");
            no_error=0;
            break;
        }
        status_str = strsep(&aff_ptr, "-");
        if (status_str == NULL || aff_ptr == NULL) {
            sprintf(err_msg, "Format of affection label %s is incorrect.",
                    affdata_str);
            errorf(err_msg);
            no_error=0;
        } else {
            if (strcmp(status_str, "*")) {
                status = atoi(status_str);
                if (status == 0) {
                    errorf("Cannot set unknown status to affected.");
                    no_error=0;
                } else if (status != 1 && status != 2) {
                    sprintf(err_msg, "Invalid status value %s.", status_str);
                    errorf(err_msg);
                    no_error=0;
                }
            }
            if (strcmp(aff_ptr, "*")) {
                liability = atoi(aff_ptr);
                if (liability < 1 || liability > num_classes) {
                    sprintf(err_msg, "Invalid liability class %s.", aff_ptr);
                    errorf(err_msg);
                    no_error=0;
                }
            }

        }
        aff_ptr = strtok(NULL, ",\n");
        if (aff_ptr == NULL) {
            break;
        }

    }
    free(aff_copy);
    return no_error;
}

static int *affected_labels(int num_classes, const char *affdata_str,
			    int *num_labels)

{
    char *aff_ptr, *status_str, *class_str;
    int l, liability, status;

    char *aff_copy;

    int *labels=CALLOC((size_t) num_classes*2, int);

    for (l=0; l < num_classes; l++) {
        labels[2*l]=-1;
        labels[2*l+1]=-1;
    }

    *num_labels = 2*num_classes;

    /* create the copy to use strsep on */
    aff_copy = strdup(affdata_str);
    aff_ptr = strtok(aff_copy, ",\n");

    while(1) {
        /* There is at least one label */
        status_str = strsep(&aff_ptr, "-");
        class_str = aff_ptr;
        if (!strcmp(status_str, "*") && !strcmp(class_str, "*"))  {
            /* ignore the rest of the affdata_str string */
            for (l=0; l < num_classes; l++) {
                labels[2*l] = (l+1) + liability_multiplier;
                labels[2*l+1] = (l+1) + liability_multiplier*2;
            }
            break;
        } else if (!strcmp(status_str, "*")) {
            /* single liability class, status = 1 or 2 */
            liability = atoi(class_str);
            labels[2*(liability-1)] = liability  + liability_multiplier;
            labels[2*liability - 1] = liability + liability_multiplier*2;
        } else if (!strcmp(class_str, "*")) {
            /* single status value, all liability classes */
            status = atoi(status_str);
            for (l=0; l < num_classes; l++) {
                labels[2*l] = (l+1)  + liability_multiplier*status;
            }
        } else {
            /* single liability, single status */
            status = atoi(status_str);
            liability = atoi(class_str);
            labels[2*(liability-1) + status - 1] = liability + liability_multiplier * status ;
        }

        aff_ptr = strtok(NULL, ",\n");
        if (aff_ptr == NULL) {
            break;
        }
    }

    free(aff_copy);
    return(labels);

}

/* This routine now simply displays counts and sets the input string,
   called only in the input mode.
   Removed check for ClassCnt, as this will be called only for traits with
   multiple liability */

static void define_labels(linkage_ped_top *Top, int tr, char *affdata_str)

{

    int index, i, ped, entry;
    int num_classes=1;
    int *found, stat, classidx, num_found=0;

    linkage_ped_rec *Entry;
    char *aff_copy;

    /* First, find out the maximum number of
       classes across all affection loci
    */

    num_classes = Top->LocusTop->Pheno[tr].Props.Affection.ClassCnt;

    found=CALLOC((size_t)(3*num_classes), int);

    if (Top->pedfile_type == POSTMAKEPED_PFT) {
        for (ped= 0; ped < Top->PedCnt; ped++) {
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry=&(Top->Ped[ped].Entry[entry]);
                index=3*(Entry->Pheno[tr].Affection.Class - 1) +  Entry->Pheno[tr].Affection.Status;
                found[index]++;
            }
        }
    } else {
        for (ped= 0; ped < Top->PedCnt; ped++) {
            for (entry = 0; entry < Top->PTop[ped].num_persons; entry++)  {
                index=
                    3*(Top->PTop[ped].persons[entry].pheno[tr].Affection.Class - 1) +
                    Top->PTop[ped].persons[entry].pheno[tr].Affection.Status;
                found[index]++;
            }
        }
    }

    num_found=0;  draw_line();
    printf("Found these status-class pairs:\n");
    for (i = 0; i < 3*num_classes; i++)  {
        if (found[i] > 0)   {
            stat = i % 3;
            classidx = (int) (i / 3) + 1;
            num_found++;
            printf("\t\t%1d-%-2d (%d entries)", stat, classidx, found[i]);
            if ((num_found % 2) == 0) printf("\n");
        }
        if ((num_found % 2) != 0) printf("\n");
    }

    /* Then enter the string */

    printf("Now enter the \"trait status\"-\"liability class\" pairs \n");
    printf("you wish to be considered as \"affected\" by %s (e.g. 2-1).\n\n",
           ProgName);

    printf("You may enter a * in either field as a wildcard \n");
    printf("e.g., 2-* means status 2 and all classes). \n");
    printf("Separate pairs with commas (e.g. 2-1,2-3) > ");
    fflush(stdout);
    (void)fgets(affdata_str, MAX_NAMELEN-1, stdin); newline;
    aff_copy = strdup(affdata_str);
    if (!check_affdata_str(aff_copy, num_classes)) {
        /* put back to the default */
        strcpy(affdata_str, "2-*");
    }
    free(aff_copy);

    return ;

}


/* get affected labels from the user for all trait_loci */
/* routines to set affected phenotype identifiers */

/* 1) default_labels() has been commented out

   2) affected_labels() <- computes the status*LM + liability values
   that are to be considered affected.
   This is slightly different, the number of labels is always 2*num_classes,
   initialized to -1. Then only the ones that are affected are set to the
   correct value. For a single class locus, there are 2  labels
   200001, 200002 (201, 202 for simwalk2).

   3) define_labels()
   Prompts the user for affection strings for each multiple class trait.

   4) check_affdata_str() checks each label to make sure that it is c-c.

*/
void define_affection_labels(linkage_ped_top *Top, analysis_type analysis)

{
    int tr, *ml_traits, i, max_classes = 0;
    int num_mult_tr=0;
    char *affdata_str, *trait_name;

/*
  If your option cannot handle affection status variables with
  liability classes, then Mega2 needs to be able to decide which
  class-status combinations are to be considered affected. So, add
  KEYWORD to the list of options inside the first case block of the
  function define_affection_labels().
*/
    if (!analysis->allow_affection_liability_class()) {
        Labels = NULL;
        return;
    }

    liability_multiplier = LIABILITYMULTIPLIER;
    ml_traits = CALLOC((size_t) num_traits, int);

    for(i=0; i < num_traits; i++) {
        SKIP_TRI(i)
            if (Top->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
                max_classes =
                    ((Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt >
                      max_classes)?
                     Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt :
                     max_classes);
                if (Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt > 1) {
                    ml_traits[num_mult_tr] = global_trait_entries[i];
                    num_mult_tr++;
                }
            }
    }

    if (SIMWALK2(analysis)) {
        if (max_classes < 100) {
            /* set the liability multiplier correctly */
            liability_multiplier = 100;
        }
    }

    /* For the single class traits, set the default labels */

    for(i=0; i < num_traits; i++) {
        SKIP_TRI(i)
            if (Top->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION) {
                if (Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt == 1) {
                    /* This simply defines a single label 201 or 200001 for single class */
                    Labels = affected_labels(1, "2-1", &NumLabels);
                    Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.Labels =
                        Labels;
                    Top->LocusTop->Pheno[global_trait_entries[i]].Props.Affection.NumLabels =
                        NumLabels;
                }
            }
    }

    if (batchAFFVALUE && Mega2BatchItems[/* 18 */ Value_Affecteds].item_read) {
        /* check each individual string */
        for (i=0; i < num_mult_tr; i++) {
            int found = -1;
            if (Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i] == NULL) {
                /* Not enough labels */
                errorvf("Not enough labels, read in %d, needed %d labels.\n", i, num_mult_tr);
                EXIT(BATCH_FILE_ITEM_ERROR);
            }

            affdata_str=strdup(Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i]);
            trait_name = strsep(&affdata_str, ":");

            /* find the trait number */
            for (tr = 0; tr < num_mult_tr; tr ++) {
                if (!strcmp(Top->LocusTop->Locus[ml_traits[tr]].Name,
                            trait_name)) {
                    found=tr;
                    break;
                }
            }
            if (found == -1) {
                /* This trait is wrong */
                errorvf("Trait %s is not an affection trait with multiple liability classes.\n",
                        trait_name);
                EXIT(BATCH_FILE_ITEM_ERROR);
            } else {
                if (check_affdata_str(affdata_str,
                                      Top->LocusTop->Pheno[ml_traits[found]].Props.Affection.ClassCnt)) {
                    Labels =
                        affected_labels(Top->LocusTop->Pheno[ml_traits[found]].Props.Affection.ClassCnt,
                                        affdata_str, &NumLabels);
                    Top->LocusTop->Pheno[ml_traits[found]].Props.Affection.Labels = Labels;
                    Top->LocusTop->Pheno[ml_traits[found]].Props.Affection.NumLabels = NumLabels;
                } else {
                    /* This string is wrong */
		    errorvf("Found problems with affection labels.\n");
                    EXIT(BATCH_FILE_ITEM_ERROR);
                }
            }
            /* commented out because strsep, advances the pointer */
            /*       free(affdata_str); */
        }
    } else {

        /* Now print the menu, one item for each trait that has multiple liability classes
         */

        int select = -1, done = 0;
        char cselect[10];
        affdata_str = CALLOC((size_t) FILENAME_LENGTH, char);

        /* set the defaults */

        if (num_mult_tr == 0) return;

        Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names = CALLOC((size_t) num_mult_tr+1, char*);

        for (i=0; i < num_mult_tr; i++) {
            Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i] = CALLOC((size_t) 100, char);
            strcpy(Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i], "2-*");
        }

        Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[num_mult_tr] = NULL;

        while(!done) {
            printf("Affection label menu:\n");
            printf("0) Done with this menu - please proceed\n");
            for(i=0; i < num_mult_tr; i++) {
                printf(" %d) %s [%s]\n", i+1, Top->LocusTop->Locus[ml_traits[i]].Name,
                       Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i]);
            }
            printf("Enter 0 or a trait item between 1-%d >",
                   num_mult_tr);
            fcmap(stdin, "%s", cselect); newline;
            sscanf(cselect, "%d", &select);
            if (select == 0) done = 1;
            else {
                define_labels(Top, ml_traits[select-1], affdata_str);
                /* check if this string is okay */
                chomp(affdata_str);
                strcpy(Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[select-1], affdata_str);
            }
        }

        /* Now set the liability and status labels */
        for(i=0; i < num_mult_tr; i++) {
            Labels = affected_labels(Top->LocusTop->Pheno[ml_traits[i]].Props.Affection.ClassCnt,
                                     Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i],
                                     &NumLabels);

            Top->LocusTop->Pheno[ml_traits[i]].Props.Affection.Labels = Labels;
            Top->LocusTop->Pheno[ml_traits[i]].Props.Affection.NumLabels = NumLabels;
            /* Now prepend each affdata-string with the trait name */
            affdata_str = strdup(Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i]);
            sprintf(Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i], "%s:%s",
                    Top->LocusTop->Locus[ml_traits[i]].Name,
                    affdata_str);
            free(affdata_str);
        }
    }

    mssgf("Status-liability combinations to consider affected:");
    for(i=0; i < num_mult_tr; i++) {
        sprintf(err_msg, "%s", Mega2BatchItems[/* 18 */ Value_Affecteds].value.mult_names[i]);
        mssgf(err_msg);
    }
    log_line(mssgf);
    if (InputMode == INTERACTIVE_INPUTMODE) {
        batchf(Value_Affecteds);
    }
}

/**
   Search the quantitative phenotype data for a "match" with the value specified.
   Return true (1) if found, else false (0) if not.
*/
static int check_quant_phenotype_data_has_value(linkage_ped_top *Top, double value) {
    int i, i1, ped, entry;

    // It's OK if the value matches that of the value used internally to mark missing values.
    // This routine is checking against the "data" and not our missing value marker.
    // MissingQuant get's set in batch_input.cpp:set_batch_items() from the batch item keyword
    // Value_Missing_Quant_On_Input
    if (fabs(MissingQuant - value) <= EPSILON) return 0;
    
    for (i1 = 0; i1 < num_traits ; i1++) {
        i = global_trait_entries[i1];
        if (i == -1) continue;
        if (Top->LocusTop->Locus[i].Type == QUANT) {
            for(ped = 0; ped < Top->PedCnt; ped++) {
                if (Top->pedfile_type == POSTMAKEPED_PFT) {
                    for(entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                        if (fabs(Top->Ped[ped].Entry[entry].Pheno[i].Quant - value) <= EPSILON) return 1;
                    }
                } else {
                    for(entry = 0; entry < Top->PTop[ped].num_persons; entry++) {
                        if (fabs(Top->PTop[ped].persons[entry].pheno[i].Quant - value) <= EPSILON) return 1;
                    }
                }
            }
        }
    }
    return 0;
}

/**
 This routine should be called from within reorder_loci.cpp:ReOrderLoci() when
 write_quant_stats() returns true indicating that a missing quantitative trait
 has been found in the trait markers that have been selected by the user. This
 routine ASSUMES that the user has selected a QTL(s) AND that one or more of
 the QTL(s) is/are missing a quant (e.g., contains the MissingQuant value).
 
 If 'Value_Missing_Quant_On_Output' has not been specified, and the analysis mode
 allows the user to define a missing value (e.g., output_quant_can_define_missing_value() == true)
 then the routine will present the user with an interactive menu that will allow
 them to set the missing quantitative value on output. If the analysis mode requires
 a numeric missing value the code will both check for a valid number, and then
 search all selected quantitative traits for that value and disallow it if found.

 If 'Value_Missing_Quant_On_Output' has been previously been specified, then the
 routine will perform a search of the selected quantitative values if the analysis
 type requires that the missing output quant be numeric. This search is done at this
 time, because it is the first time that we know exactly which traits have been
 selected by the user, and so know what we should search.
 
 There are several methods kept in the anlaysis option class that are used by this
 routine. For a complete definition of them, please see their virtual function
 definitions of CLASS_ANALYSIS in analysis.h. They are: output_quant_can_define_missing_value(),
 output_quant_must_be_numeric(), and output_quant_default_value().
 */
void set_missing_quant_output(linkage_ped_top *Top, analysis_type analysis)
{
    int select=-1;
    double value;
    char missingq[200], *end;
    
    // Some analysis options have a 'hard wired' missing quantitative value.
    // For these options that value will be used later in the program so there
    // is no point in asking the user now or checking the value of the batch
    // file item Value_Missing_Quant_On_Output.
    if (!analysis->output_quant_can_define_missing_value()) return;
    
    if (!Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read) {
        
        if (InputMode == BATCH_FILE_INPUTMODE) {
            errorvf("The batch file item 'Value_Missing_Quant_On_Output' must be specified because\n");
            errorvf("one or more QTL has been found to be missing in the input data.\n");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
        
        // The value returned by analysis->output_quant_default_value() is based on the analysis mode
	// in the class definition (see the 'analysis.h' file). By default it is NULL, but that shouldn't
	// matter because we don't get here unless analysis->output_quant_can_define_missing_value() == true
	// but we check it anyway just in case the analysis mode was misconfigured.
        strcpy(missingq, (analysis->output_quant_default_value() != NULL ? analysis->output_quant_default_value() : ""));
        while (select != 0) {
            draw_line();
            printf("Missing output QTL value definition menu:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) Assign missing quantitative value on ouput to : %s\n", missingq);
            printf("Enter 0 or 1 > ");
            fcmap(stdin, "%d", &select); newline;
            if (select == 1) {
                printf("Enter missing quant value > ");
                fcmap(stdin, "%s", missingq); newline;
                if (analysis->output_quant_must_be_numeric()) {
                    // Since the missing quantitative value must be numeric,
                    // attempt to turn the string that was entered to a floating point value...
                    value = strtod((const char *)missingq, &end);
                    // Did a successful numeric conversion take place?
                    if (strlen(missingq) == 0 || strlen(end) != 0 || errno == ERANGE) {
                        select = -1;
                        printf("INPUT ERROR: Analysis type '%s' requires quantitative values to be numeric.\n",
                               analysis->_name);
                        printf("INPUT ERROR: Converting specified value '%s' to a floating point number.\n", missingq);
                        printf("INPUT ERROR: The entire string is not a valid floating point number.\n");
                        strcpy(missingq, analysis->output_quant_default_value());
                    } else if (check_quant_phenotype_data_has_value(Top, value)) {
                        // Since we have a valid numeric value, we look to see if it is in the input now...
                        select = -1;
                        printf("INPUT ERROR: Value specified '%s' was found in the input quantitative data.\n", missingq);
                        strcpy(missingq, analysis->output_quant_default_value());
                    }
                }
            }
        }
        
        Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read = 1;
        strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name, missingq);
        if (InputMode == INTERACTIVE_INPUTMODE) batchf(Value_Missing_Quant_On_Output);
        
    } else if (analysis->output_quant_must_be_numeric()) {
        // In 'InputMode == BATCH_FILE_INPUTMODE' we were not able to make this check till now
        // since before this point we were not able to see that the value specified in the
        // batch file item was in the input data...
        char *end;
        double value;
        // We will not check for an error in the number conversion since we assume that this
        // was done on batch file processing, or we would not have gotten this far in the code.
        value = strtod((const char *)Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name, &end);
        if (check_quant_phenotype_data_has_value(Top, value)) {
            // The value was found in the QTL data selected by the user...
            if (analysis->output_quant_can_define_missing_value() == false &&
                analysis->output_quant_default_value() != NULL &&
                strcmp(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name,
                       analysis->output_quant_default_value()) == 0) {
                    // Here the user was not permitted to specify a missing value, but still
                    // their data contained the default.
                    warnvf("The default QTL value on output '%lf' was found in the\n", value);
                    warnf("selected input quantitative trait data.");
                } else {
                    // Here the user was permitted to specify a missing value, and their data
                    // contained it.
                    errorvf("'Analysis_Option' = '%s' while specifying 'Value_Missing_Quant_On_Output'.\n",
                            analysis->_name);
                    errorvf("Value specified '%lf' was found in the selected input quantitative trait data.\n", value);
                    EXIT(BATCH_FILE_ITEM_ERROR);
                }
        }
    }
    mssgvf("NOTE: The Missing QTL value on output will be assigned as '%s'.\n",
           Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
}

/**
 This routine will set the MissingQuant value for the input data.
 
 This routine should be called from within reorder_loci.cpp:ReOrderLoci() before write_quant_stats()
 so that missing values in the input can be identified. We always need to know what the missing quant
 value is on input, because we need to look for it...
 
 The missing quant can come from one of three places:
 1) for PLINK input if --missing_phenotype is specified it comes from there (a warning is
 given if Value_Missing_Quant_On_Input is also specified),
 2) if the batch file item Value_Missing_Quant_On_Input is specified (not PLINK input) it comes from there,
 3) finally we need to ask the user.
 
 NOTE: in the past, we NEVER asked the user if we were processing a PLINK file and --missing_phenotype
 was not specified, we forced them to take the default (-9). Since the user may have used a different
 value for missing phenotype, but did not specify it with the --missing_phenotype parameter, this was
 not good.
 */
void set_missing_quant_input(linkage_ped_top *Top, const analysis_type analysis)
{
    int i, i1, ped, entry;
    
    // If we will never encounter a quant, then there is no point in caring about the value for a missing one.
    if ((!HasQuant && num_covariates == 0) || !QTL_ALLOW) return;

    //
    // The missing quant can come from one of three places:

    if (PLINK.plink && PLINK.missing_pheno) {
        // 1) for PLINK input if --missing_phenotype is specified it comes from there (a warning is
        // given if Value_Missing_Quant_On_Input is also specified),
#ifndef HIDESTATUS
        if (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].item_read) {
            warnvf("Missing QTL phenotype value specified by Value_Missing_Quant_On_Input will be ignored.\n");
        }
        mssgvf("Missing QTL phenotype value read from --missing_phenotype (%f).\n", PLINK.pheno_value);
#endif
        MissingQuant = PLINK.pheno_value;
        if (analysis == TO_SAGE) {
#ifndef HIDESTATUS
            mssgf("Decimal places in will be ignored since SAGE");
            mssgf("accepts only integers as missing values.");
#endif
            MissingQuant = (MissingQuant > 0 ? floor(MissingQuant) : ceil(MissingQuant));
        }
        
    } else if (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].item_read) {
        // 2) if the batch file item Value_Missing_Quant_On_Input is specified (not PLINK input) it comes from there,
        //
        // NOTE: batch_input.cpp:set_batch_items() will look for the string 'NA' in this field, and if
        // found substitute the value QMISSING. So the 'NA' to QMISSING conversion has already taken
        // place at this point in the code.
        MissingQuant = Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue;
        if (analysis == TO_SAGE) {
#ifndef HIDESTATUS
            mssgf("Missing QTL phenotype value read in from batch file.");
            mssgf("Decimal places in will be ignored since SAGE");
            mssgf("accepts only integers as missing values.");
#endif
            MissingQuant = (MissingQuant > 0 ? floor(MissingQuant) : ceil(MissingQuant));
        }
        
    } else {
        // 3) finally since it was not otherwize specified, we need to ask the user.
        char missingq[200];
        int select=-1;
        strcpy(missingq, PLINK.plink ? "-9" : "NA");
        MissingQuant = PLINK.plink ? -9 : QMISSING;   // NOTE: QMISSING is the internal numeric representation of "NA".
        printf("Selected loci contain one or more QTLs or covariates.\n");
        // Here, the user is forced to enter either 'NA' or a number representable as a double...
        while (select != 0) {
            draw_line();
            if (analysis == TO_SAGE) {
                printf("Decimal places will be ignored since SAGE\n");
                printf("accepts only integers as missing values.\n");
            }
            printf("Missing input quantitative value definition menu:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) Assign missing value for input: %s\n", missingq);
            printf("Enter 0 or 1 > ");
            fcmap(stdin, "%d", &select); newline;
            if (select == 1) {
                char *end;
                double missingq_value;
                printf("Enter missing quant value > ");
                fcmap(stdin, "%s", missingq); newline;
                
                if (strcasecmp(missingq, "NA") == 0) {
                    MissingQuant =  QMISSING;
                } else {
                    // At this point, the missing quant value must be numeric...
                    
                    // strtod(const char *restrict nptr, char **restrict endptr);
                    // 1) if strlen(nptr) == 0 then there is nothing to convert.
                    // 2) If endptr is not NULL, a pointer to the character after the last character used
                    // in the conversion is stored in the location referenced by endptr. So strlen(endptr)
                    // should be zero if the entire string was converted.
                    // 3) If no conversion is performed, zero is returned and the value of nptr is stored
                    // in the location referenced by endptr.
                    // 4) If the correct value would cause underflow or overflow, ERANGE is stored in errno.
                    // NOTE: Numbers in exponential notation are valid here (e.g., 1.2e5).
                    missingq_value = strtod((const char *)missingq, &end);
                    if (strlen(missingq) == 0 || strlen(end) != 0 || errno == ERANGE) {
                        // Conversion of the entire string was not successful or some other error...
                        printf("ERROR: The string that you entered '%s' was not valid?\n", missingq);
                        printf("ERROR: It is either not 'NA' or a number representable as a double.\n");
                        strcpy(missingq, PLINK.plink ? "-9" : "NA");
                        select = -1;
                        continue;
                    }
                    MissingQuant = missingq_value;
                    if (analysis == TO_SAGE) {
                        MissingQuant = (MissingQuant > 0 ? floor(MissingQuant) : ceil(MissingQuant));
                    }
                }
            } // if (select == 1) {
        } // while (select != 0) {
        
        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].item_read = 1;
        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue = MissingQuant;
        if (InputMode == INTERACTIVE_INPUTMODE) batchf(Value_Missing_Quant_On_Input);
    }
    
    // At this point MissingQuant is either QMISSING (e.g., "NA") or some other number...
    
    if (Top == NULL) return;
    
    //
    // Since MissingQuant has just been firmly established, this is the point in the code
    // where we convert any "NA"s (QMISSINGs) that we have read to the value of MissingQuant.
    //
    // NOTE: It is pointless to do this if 'MissingQuant == QMISSING' as would be the case if we
    // were reading Mega2 annotated files with NA used as the missing value.
    if (fabs(MissingQuant - QMISSING) > EPSILON) {
        
        if (HasQuant) {
            for (i1 = 0; i1 < num_traits ; i1++) {
                SKIP_TRI(i1)
                i = global_trait_entries[i1];
                if (Top->LocusTop->Locus[i].Type == QUANT) {
                    for(ped = 0; ped < Top->PedCnt; ped++) {
                        if (Top->pedfile_type == POSTMAKEPED_PFT) {
                            for(entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                                if (fabs(Top->Ped[ped].Entry[entry].Pheno[i].Quant - QMISSING) <= EPSILON) {
                                    /* Set this to the proper missing quant value */
                                    Top->Ped[ped].Entry[entry].Pheno[i].Quant = MissingQuant;
                                }
                            }
                        } else {
                            for(entry = 0; entry < Top->PTop[ped].num_persons; entry++) {
                                if (fabs(Top->PTop[ped].persons[entry].pheno[i].Quant - QMISSING) <= EPSILON) {
                                    /* Set this to the proper missing quant value */
                                    Top->PTop[ped].persons[entry].pheno[i].Quant = MissingQuant;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        if (num_covariates > 0) {
            for (i1 = 0; i1 < num_covariates ; i1++) {
                i = covariates[i1];
                for(ped = 0; ped < Top->PedCnt; ped++) {
                    if (Top->pedfile_type == POSTMAKEPED_PFT) {
                        for(entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                            if (fabs(Top->Ped[ped].Entry[entry].Pheno[i].Quant - QMISSING) <= EPSILON) {
                                /* Set this to the proper missing quant value */
                                Top->Ped[ped].Entry[entry].Pheno[i].Quant = MissingQuant;
                            }
                        }
                    } else {
                        for(entry = 0; entry < Top->PTop[ped].num_persons; entry++) {
                            if (fabs(Top->PTop[ped].persons[entry].pheno[i].Quant - QMISSING) <= EPSILON) {
                                /* Set this to the proper missing quant value */
                                Top->PTop[ped].persons[entry].pheno[i].Quant = MissingQuant;
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (analysis == TO_GHMLB) {
        /*    prog_name(analysis, progname); */
        warnf("***************");
        warnvf("Quantitative phenotypes should be normalized for %s.\n", ProgName);
        warnf("***************\n");
    }
}

char      *untyped_ped_messg(int opt, char *untyped_mssg)

{
    int opt1;


    opt1 = ((opt >= 4)? 4: opt);

    switch(opt1) {
    case 1:
        strcpy(untyped_mssg, "Omit any completely untyped pedigrees");
        break;

    case 2:
        strcpy(untyped_mssg, "Include all pedigrees whether typed or not");
        break;

    case 3:
        strcpy(untyped_mssg,
               "Exclude pedigrees not fully typed at one or more markers");
        break;

    case 4:
        ((opt == 4)?
         sprintf(untyped_mssg,
                 "Exclude any pedigree with 0 marker-typed people") :
         sprintf(untyped_mssg,
                 "Exclude any pedigree with %d or less marker-typed people",
                 opt - 4));
        break;
    default:
        /*    for taking case of negative numbers */
        printf("Unknown option %d.\n", opt);
        break;
    }
    return &(untyped_mssg[0]);
}

static void     untyped_ped_menu(int *opt)

{
    int new_opt=-1;
    char opt_[10], select[5] = {"    "};
    char messg[100];

    if (*opt <= 4) {
        select[*opt - 1]='*';
    } else {
        select[3]='*';
    }
    printf("Untyped pedigree exclusion options:\n");
    while(new_opt != 0) {
        new_opt = ((*opt < 4)? 5 : *opt);
        printf("0) Done with this menu - please proceed.\n");
        printf("%c1) %s\n", select[0], untyped_ped_messg(1, messg));
        printf("%c2) %s\n", select[1], untyped_ped_messg(2, messg));
        printf("%c3) %s\n", select[2], untyped_ped_messg(3, messg));
        printf("%c4) %s\n", select[3], untyped_ped_messg(new_opt, messg));
        printf("Select from options 0 - 4 > ");
        fcmap(stdin, "%s", opt_); newline;
        new_opt = atoi(opt_);
        switch(new_opt) {
        case 0:
            break;
        case 1:
            strcpy(select, "*   ");
            *opt=new_opt;
            break;
        case 2:
            strcpy(select, " *  ");
            *opt=new_opt;
            break;
        case 3:
            strcpy(select, "  * ");
            *opt=new_opt;
            break;
        case 4:
            strcpy(select, "   *");
            *opt = -1;
            while(*opt < 0) {
                printf("Enter minimum number of marker-typed people in a pedigree\n");
                printf("for it to be included in the output file > ");
                fcmap(stdin, "%s", opt_); newline;
                *opt = atoi(opt_);
                if (*opt < 0) {
                    printf("Negative number not allowed, please enter again.\n");
                }
                else *opt -= 1 ;
            }
            *opt += 4;

            break;
        default:
            warn_unknown(opt_);
            break;
        }
    }
    return;
}


int  gh_cov_selection(int num_traits, int *traits,
		      int *covariates, linkage_locus_top *LTop)

{

    int t, n, numq, num_cov=0, done=0, *quants, *selected;
    char cselect[200]="", *str_p, cnum[10];


    numq=0;
    for(t=0; t<num_traits; t++) {
        if (LTop->Locus[traits[t]].Type == QUANT) numq++;
    }

    if (numq == 0) return 0;

    selected = CALLOC((size_t) numq, int);
    quants = CALLOC((size_t) numq, int);
    for(t=0; t<numq; t++) selected[t]=0;

    numq=0;
    for(t=0; t<num_traits; t++) {
        if (LTop->Locus[traits[t]].Type == QUANT) {
            quants[numq]=t; /* index of element in trait */
            numq++;
        }
    }

    draw_line();
    printf("Covariate selection menu:\n");
    printf("Traits denoted with an asterisk will be ");
    printf("processed as covariates.\n");
    printf("Enter only 'e' if no covariates are desired.\n");
    while(!done) {
        draw_line();
        printf("Trait loci list:\n");
        for (t=0; t<numq; t++) {
            printf("%c %d) %s\n", (selected[t]? '*': ' '),
                   t+1, LTop->Locus[traits[quants[t]]].Name);
        }
        draw_line();
        printf("\nEnter only 'e' if no covariates are desired.\n");
        printf("Enter trait numbers ('e' to terminate) > ");
        fflush(stdout);
        (void)fgets(cselect, 199, stdin);
        newline;
        str_p=&(cselect[0]);
        while(*str_p != '\0') {
            if (num_cov >= num_traits) {
                done=1; break;
            } else if (*str_p == 'e' || *str_p == 'E') {
                done=1; break;
            } else if (*str_p == '\n') {
                break;
            } else if (isspace((int) *str_p)) str_p++;
            else if (isdigit((int) *str_p)) {
                n=0;
                while(isdigit((int) *str_p)) {
                    cnum[n]=*str_p; n++; str_p++;
                    if (n>6) {
                        errorvf("Locus number too big!\n");
			EXIT(DATA_TYPE_ERROR);
                    }
                }
                cnum[n]='\0';
                if (atoi(cnum) >= 1 && atoi(cnum) <= num_traits) {
                    selected[atoi(cnum)-1]=1;
                    num_cov++;
                } else {
                    printf("Invalid selection %s, please re-select.\n",
                           cnum);
                    for (t=0; t<num_traits; t++) {
                        selected[t]=0;
                    }
                    strcpy(cnum, ""); strcpy(cselect, "");
                    str_p = &(cselect[0]);
                    break;
                }
            } else {
                printf("Non-numeric characters in input, please re-select.\n");
                for (t=0; t<num_traits; t++) {
                    selected[t]=0;
                }
                strcpy(cnum, ""); strcpy(cselect, "");
                str_p = &(cselect[0]);
                break;
            }
        }
    }

    if (num_cov > 0) {
        for(t=0; t<numq; t++) {
            if (selected[t]==1) {
                covariates[quants[t]]=1;
            }
        }

        sprintf(err_msg, "QTL/covariate list:");
        mssgf(err_msg);
        for (t=0; t<numq; t++) {
            sprintf(err_msg, "%d) %s %s", t+1,
                    LTop->Locus[traits[quants[t]]].Name,
                    (covariates[quants[t]]? "[covariate]":""));
            mssgf(err_msg);
        }
        draw_line();
    }
    free(selected);
    free(quants);

    return num_cov;
}

static char *ind_id_choice_messg(int opt, char *ind_id_messg)

{
    switch(opt) {
    case 1:
        sprintf(ind_id_messg, "Individual id");
        break;

    case 2:
        sprintf(ind_id_messg, "Per: field");
        break;

    case 3:
        sprintf(ind_id_messg, "ID: field");
        break;

    case 4:
        sprintf(ind_id_messg, "Unique id e.g 1_2, 1=ped, 2=ind");
        break;

    case 5:
        sprintf(ind_id_messg, "Renumber consecutively in pedigree");
        break;

    default:
        break;
    }

    return &(ind_id_messg[0]);
}


int individual_id_item(int item_number, analysis_type analysis,
                       int current_opt_val, int align_right_pos,
                       int get_disp_log, int has_orig, int has_uniq)
{
    /* This displays a string as a menu item
       inside the file name menu as item number <item_number> for
       certain analysis types. Its value is set by <current_opt_val>
       which is 0,1, or,2(default 1), and the string is
       right-justified at  <align_right_pos > 35> (the item message will
       be a constant width = 39). */

    char ind_id_messg[45];
    char opt_[10];
    int new_opt = -1;
    int i;


    if (get_disp_log == 3) {
        /* log */
        sprintf(err_msg, "Person id in output pedigree file = %s",
                ind_id_choice_messg(current_opt_val, &(ind_id_messg[0])));
        mssgf(err_msg);
        log_line(mssgf);
        return 0;
    }

    if (get_disp_log==2) {
        printf(" %d) Person id in output pedigree file:        %s\n",
               item_number,
               ind_id_choice_messg(current_opt_val, &(ind_id_messg[0])));
        return 0;
    }


    if (get_disp_log == 1) {

        /* print the menu */
        printf("Output person id selection menu:\n");
        while (new_opt != 0) {
            printf("0) Done with this menu - please proceed.\n");
            i=1;
            printf("%c%d) %s\n",
                   ((current_opt_val == 1) ? '*' : ' '),
                   i,
                   ind_id_choice_messg(1, &(ind_id_messg[0])));
            i++;
            if (has_orig) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(2, &(ind_id_messg[0])));
                i++;
            }
            if (has_uniq) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 3) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(3, &(ind_id_messg[0])));
            } else {
                printf("%c%d) %s\n",
                       ((current_opt_val == 4) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(4, &(ind_id_messg[0])));
            }
            i++;

            printf("%c%d) %s\n",
                   ((current_opt_val == 5) ? '*' : ' '),
                   i,
                   ind_id_choice_messg(5, &(ind_id_messg[0])));

            printf("Select from options 0 - %d > ", i);
            fcmap(stdin, "%s", opt_);
            new_opt = atoi(opt_);

            switch(new_opt) {
            case 0:
                break;
            case 1:
                /* Always individual id */
                current_opt_val = new_opt;
                break;
            case 2:
                /* Could be per: id: or regenerated */
                if (has_orig) {
                    current_opt_val = 2;
                } else {
                    current_opt_val = ((has_uniq)? 3 : 4);
                }
                break;
            case 3:
                /* could be unique or generated or renumbered*/
                if (has_orig) {
                    current_opt_val = ((has_uniq)? 3: 4);
                    break;
                } else {
                    current_opt_val = 5;
                }
                break;
            case 4:
                if (has_orig) {
                    current_opt_val = 5;
                    break;
                }
                /* otherwise default behaviour */
            default:
                warn_unknown(opt_);
                break;
            }
        }

    }
    return current_opt_val;
}

static char *ped_id_choice_messg(int opt, char *id_messg)

{

    switch(opt) {
    case 1:
        if (pedfile_type == POSTMAKEPED_PFT) {
            sprintf(id_messg, "Pedigree number");
            break;
        }

    case 2:
        if (pedfile_type == POSTMAKEPED_PFT) {
            sprintf(id_messg, "Ped: field");
            break;
        } else {
            sprintf(id_messg, "Premakeped pedigree number.");
        }
        break;

    case 3:
        sprintf(id_messg, "Renumbered consecutively");
        break;

    case 4:
        /* This is only for nuclear pedigrees */
        sprintf(id_messg, "With extensions e.g. 1_2 etc.");
        break;

    case 5:
        /* This is a different type of an extension which is
           purely numeric */
        sprintf(id_messg, "With multipliers e.g 1002 etc.");
        break;

    default:
        break;
    }

    return &(id_messg[0]);
}

int pedigree_id_item(int item_number, analysis_type analysis,
                     int current_opt_val, int align_right_pos,
                     int get_disp_log, int has_orig)

{

    /* This displays a string as a menu item
       inside the file name menu as item number <item_number> for
       certain analysis types. Its value is set by <current_opt_val>
       which is 0,1, or,2(default 1), and the string is
       right-justified at  <align_right_pos > 35> (the item message will
       be a constant width = 39). */

    char id_messg[45];
    char opt_[10];
    int new_opt = -1;
    int i;

    if (get_disp_log == 3) {
        /* log */
        sprintf(err_msg, "Pedigree id in output pedigree file = %s",
                ped_id_choice_messg(current_opt_val, &(id_messg[0])));
        mssgf(err_msg);
        log_line(mssgf);
        return 0;
    }

    if (get_disp_log == 2) {
        printf(" %d) Pedigree id in output pedigree file:      %s\n",
               item_number,
               ped_id_choice_messg(current_opt_val, &(id_messg[0])));
        return 0;
    }

    if (get_disp_log == 1) {

        printf("Output pedigree id selection menu:\n");
        while(new_opt != 0) {
            printf("0) Done with this menu - please proceed\n");
            i=1;
            if (analysis == TO_NUKE) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 5) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(5, &(id_messg[0])));
                i++;
            } else if (pedfile_type == POSTMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 1) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(1, &(id_messg[0])));
                i++;
            }

            if (analysis == TO_NUKE) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 4) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(4, &(id_messg[0])));
                i++;
            } else if (has_orig && pedfile_type == POSTMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(2, &(id_messg[0])));
                i++;
            } else if (pedfile_type == PREMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(2, &(id_messg[0])));
                i++;
            }
            printf("%c%d) %s\n",
                   ((current_opt_val == 3) ? '*' : ' '),
                   i,
                   ped_id_choice_messg(3, &(id_messg[0])));

            printf("Select from options 0 - %d > ", i);
            fcmap(stdin, "%s", opt_);
            new_opt = atoi(opt_);

            switch(new_opt) {
            case 0: break;
            case 1:
                if (analysis == TO_NUKE) {
                    current_opt_val = 5;
                } else if (pedfile_type == PREMAKEPED_PFT) {
                    current_opt_val = 2;
                }
                break;

            case 2:
                if (analysis == TO_NUKE) {
                    current_opt_val = 4;
                } else if (has_orig && pedfile_type == POSTMAKEPED_PFT) {
                    current_opt_val = 2;
                } else {
                    current_opt_val = 3;
                }
                break;

            case 3:
                if (has_orig) {
                    current_opt_val = 3;
                    break;
                }
            default:
                warn_unknown(opt_);
                break;
            }
        }
    }
    if (current_opt_val == 1 && pedfile_type == PREMAKEPED_PFT) {
        current_opt_val = 2;
    }
    return current_opt_val;

}

void ped_ind_defaults(int unique, analysis_type analysis)
{

    /* set default values based on analysis */

    /*
      Add keyword to ped_ind_defaults(). This function decides which of
      the three individual ids, and two pedigree ids will be selected as
      the default output ids. This is critical if the target analysis has
      restrictions on ids (e.g., linkage-format only handles numerical values).

      Pedigree ID choices:
      1=Pedigree field for post-makeped format
      2=Pedigree field for pre-makeped
      3=Consecutive numerical ids for nuclear pedigrees

      Individual ID choices:
      1= Person field for post- or pre-makeped files
      2 = Unique Ids
      5 = Consecutively numbered person ids
      3,4 = 	not used, although these are handled.
    */

    analysis->ped_ind_defaults(unique);

}


void test_modified(int choice)
{
    if (FirstIterMenu == 1 && choice > 0) {
        FirstIterMenu = 0;
    }
}



// Cranefoot, and smmary needs order but not position. In these cases we could use physical also.
// For a map there is order and relative/absolute position.
// Programs that require genetic maps are only interested in relative positions
// Something that uses phusical maps wants the ablolute position.
// Something that wants thetas only cares about relative positoin.

// Issues:
// * Restriction issues: (mega2 might catch these; make sure; those would lead to early exists from mega2)
// 1) Trait type restrictions
// 1a) if you need quantitative traits, and you don't give one mega2 should catch this.
// 1b) no markers required for SOLAR (h^2 estimation), or CRANEFOOT.
// 2) How many traits (mega2 could fail if you have too many; does mega2 stop the user if they pick too many or not enought)
// 3) whether you have to split it into chromosome specific marker files (e.g. MERLIN).
// 4) looping over traits makes that you have only one trait per file for a set of output files. Could create subdirectories.
// 5) some things won't run if no markers are specified.
// * the use of a physical map to disambiguate same position genetic maps (e.g., one or more marker that has the same cM position)
// sometime the genetic positions are not precise enough, you can sometimes use the phuysical to disambiguate.
// Now it assumes that the order in the map file is in the right order, so it adds a little increment.

// NOTE: SEE Summary requires q genetic map but mifght be able to get away with a physical


// For a list of valid analysis types, see common.h:enum analysis_type;
// 
// For 'AnalysisOpt' in 'RequiresGeneticMap', do not give the user the 'None' option when
// allowing them to choose the Genetic Distance Map.
// NOTE: Currently this should be all but PLINK abnd IQLS.
// For 'AnalysisOpt' in 'PhysicalMapOnly', do not present the user with the Genetic Distance Map menu.
analysis_type PhysicalMapOnly[] = {
  IQLS
};
#define SIZEOF_PHUYSICALMAPONLY (sizeof(PhysicalMapOnly) / sizeof(analysis_type))

// Any program that supports linkage should be able to support sex-specific maps)...
// In the spread sheet look at "Output Format Type" Linkage.
// where is SimWalk2???
analysis_type SexSpecificMapSupport[] = {
  IQLS, TO_MENDEL, TO_MENDEL4, TO_LINKAGE, TO_PREMAKEPED, TO_NUKE, TO_VITESSE, TO_SLINK, TO_GeneHunter, TO_GeneHunterPlus
};


static void print_information_based_on_genetic_distance_sex_type_map(genetic_distance_map_type gdsm) {
    switch (gdsm) {
        case FEMALE_GDMT:
            printf("With the current option, only analysis on the X-chromosome\n will be possible.\n");
            break;
        case SEX_AVERAGED_GDMT:
            printf("With the current option, any X-chromosome processing will be done\n with a sex-averaged genetic map.\n");
            break;
        case SEX_SPECIFIC_GDMT:
            printf("With the currently option, any processing requiring a genetic map\n will be done with sex-specific maps only\n");
            break;
        case UNKNOWN_GDMT:
        case NONE_AVAILABLE_GDMT:
            // This should not occur, and if it does we won't be printing anything for it.
            break;
    }
}

/*
Genetic Distance Map Selection [cpk]

Determine which genetic distance index from EXLTop that we should use for subsequent
processing. This is necessary since the input map file may specify several possible maps.

This routine is used both for checking existing values (from a batch file?), or
allowing the user to specify the map in interactive input mode. Possibly the thing
to do is to separate these, but they may be set anywhere, and doing this catches any
problems that the setting code may have overlooked.

>>For interactive input mode, the algorithm is roughly as follows:

========================================
Genetic map selection menu
0) Done with this menu - please proceed.
* 1) MapF: female           Kosambi
  2) Map:  sex-averaged     Kosambi
  3) Map:  sex-specific     Kosambi
  4) MapA: sex-averaged     Haldane
  5) MapQ: sex-averaged     Haldane
  6) MapY: sex-averaged     Kosambi
  7) MapY: female           Kosambi
  8) None
Select from options 0 - 8 >
========================================

RULES:
1) If ave, present only sex-averaged
2) If ave + (m), present only sex-averaged
3) If ave + (f), present as sex-averated on one line, and female on the next line.
4) if ave + (m,f), present as sex-averaged on one line, and as sex-specific on the next line.
5) if a map cannot be selected (e.g. male only map) then it shouldn't be listed in the selection menu as an option.
6) if the user has chosen an analysis option that requires a genetic map, then they should
   not be presented with the option 'None' in the Genetic map selection menu.
7) if the map is empty (e.g., all 0's) then do not present it to the user.

INFORMATIVE MESSAGES:
After making a choice, the user give informative messages (NOTE) as appropriate...
1) None: No genetic map will be used. Only analysis requiring only a physical map will be possible.
2) Male map only: This map cannot be used for analysis. No genetic map will be used. Only analysis requiring only a physical map will be possible.
3) Female map only: Only analysis on the X chromosome will be possible
4) sex-averaged map: X-chromosome processing will be done with a sex-averaged map
5) Make & Female map only: processing requiring a genetic map will be done with sex specific maps only

NOTE:
For the moment we always require a genetic map since these things are entangled into the guts of MEGA2.
*/
void get_genetic_distance_index(ext_linkage_locus_top *EXLTop) {
    int i, l, gds = 0, gds2, *gdi = NULL, option, option_selected, max_map_name_len = 0;
    int requires_genetic_map_p = 0, valid_map_p = 0, valid_sex_map_p = 0;
    genetic_distance_map_type *gdsm = NULL;
    char format_str[100];
    
    //
    // At this point we are either in INTERACTIVE_INPUTMODE, or BATCH_INPUTMODE...
    
    requires_genetic_map_p = AnalysisOpt->allow_no_genetic_map() ? 0 : 1;
    // always require a genetic map for now since genetic maps are entangled in the code...
    requires_genetic_map_p = 1;
    
    //
    // Here we get a list of the valid maps from the structure that holds the parsing of
    // the map file. This is needed to both make sure an existing value specifies a real
    // map in the map file, and to present the user a list of maps for them to choose
    // from in interactive input mode.
    
    // NOTE: For annotated map files, indexes are assigned to maps in a map file by reading
    // it from left to right in order of the first occurrance of a map-name. It is an error
    // for the map-name to be associated with more than one map-type. It is an error for the
    // map-name to have a duplicate sex-designator.
    //
    // Annotated maps are defined as follows:
    // map-designator == map-name . map-type. sex-designator
    // where: map-type is one of <h, k p>; sex-designator is one of <a, m, f>
    
    // gdi holds the indicies of valid genetic distance maps.
    gdi = CALLOC((size_t)(3*EXLTop->MapCnt), int);
    // gdsm tells whether the map is: sex-averaged (a), sex-specific (m,f), or female (f)
    gdsm = CALLOC((size_t)(3*EXLTop->MapCnt), genetic_distance_map_type);
    // gds is the index into the gdsm (map) array which specifies a valid genetic map (Index & SexType).
    // For example, if the (annotated) map file looks like this...
    // Chromosome      Map.k.a Name    Map.k.m Map.k.f BP.p
    // 3 genetic maps will be possible: Map.k (sex-averaged), Map.k (sex-specific), Map.k (female)
    // and so after exiting this for loop gds will be 3.
    
    // For each available map name determine the valid maps and types...
    for (i = 0, gds = 0; i < EXLTop->MapCnt; i++) {
        // The only map functions that specify a genetic distance are haldane and kosambi
        // and then only if they contain a valid map...
        if (EXLTop->map_functions[i] == 'h' || EXLTop->map_functions[i] == 'k') {
            if (EXLTop->SexMaps[i][SEX_AVERAGED_MAP] != 0
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
		&& EXLTop->valid_map_p[i][SEX_AVERAGED_MAP]
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
		) {
                // there is a sex averaged map and it is valid...
                gdi[gds] = i; gdsm[gds++] = SEX_AVERAGED_GDMT;
                l = (int)strlen(EXLTop->MapNames[i]);
                max_map_name_len = (l > max_map_name_len ? l : max_map_name_len);
            }
            // you can have a sex-specific and a female map...
            // male only maps are not permitted...
            if (EXLTop->SexMaps[i][FEMALE_SEX_MAP] != 0
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
		&& EXLTop->valid_map_p[i][FEMALE_SEX_MAP]
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
		) {
                if (EXLTop->SexMaps[i][MALE_SEX_MAP] != 0
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
		    && EXLTop->valid_map_p[i][MALE_SEX_MAP]
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
		    ) {
                    // with a male and female map, we can have a sex-specific map...
                    gdi[gds] = i; gdsm[gds++] = SEX_SPECIFIC_GDMT;
                    l = (int)strlen(EXLTop->MapNames[i]);
                    max_map_name_len = (l > max_map_name_len ? l : max_map_name_len);
                }
                
                // there is a female map, and it is valid...
                gdi[gds] = i; gdsm[gds++] = FEMALE_GDMT;
                l = (int)strlen(EXLTop->MapNames[i]);
                max_map_name_len = (l > max_map_name_len ? l : max_map_name_len);
            }
            // if only male, then the map will be ignored here...
            if (EXLTop->SexMaps[i][SEX_AVERAGED_MAP] == 0 &&
                EXLTop->SexMaps[i][FEMALE_SEX_MAP] == 0 &&
                EXLTop->SexMaps[i][MALE_SEX_MAP] != 0
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
		&& EXLTop->valid_map_p[i][MALE_SEX_MAP]
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
		) {
                mssgvf("Male-only maps are not supported, as they cannot be used for analyis.\n");
            }
        }
    }
    // If no genetic map is available...
    if (gds == 0 && requires_genetic_map_p == 1) {
        // yet a map is required...
        errorvf("For the analysis type specified, a genetic map is required,\n       but none are available in the map file.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    
    // Check to see if a map has been specified.
    // The only way that this would happen to this point is that it was read from the batch file, or
    // set someplace where the EXLTop stucture is being created (e.g., read_files.c:make_EXLTop_from_LTop).
    // Here we check for 'None' and a specific map value (e.g., 'unknown == -1 is not checked here
    // since it is the default value when mega2 is started up).
    if (genetic_distance_index == -2) {
        // 'None' was specified for a map in the batch file...
        if (requires_genetic_map_p == 1) {
            errorvf("For the analysis type specified, a genetic map is required.\n");
            errorvf("However, Value_Genetic_Distance_Index was specified as -2 (None).\n");
            EXIT(INPUT_DATA_ERROR);
        } else {
            // So, 'None' is OK...
            free(gdi); free(gdsm);
            return;
        }
    } else if (genetic_distance_index >= 0) {
        // A map has been specifid, so we make sure that it is valid by searching for it
        // in the list created above and set the appropriate booleans if encountered...
        for (i=0; i<gds; i++)
            if (genetic_distance_index == gdi[i]) {
                valid_map_p = 1;
                if (genetic_distance_sex_type_map == gdsm[i]) {
                    valid_sex_map_p = 1;
                    break;
                }
            }
        
        if (!valid_map_p) {
            // didn't find the map...
           errorvf("The Value_Genetic_Distance_Index specified (%d) does not reference a genetic map in the map file.\n",
                    genetic_distance_index);
            EXIT(INPUT_DATA_ERROR);
        } else if (!valid_sex_map_p) {
            // didn't find the sex_type...
            errorvf("The Value_Genetic_Distance_SexTypeMap specified (%d) in the Map file does not reference a valid map type for the Value_Genetic_Distance_Index specified (%d).\n",
                    genetic_distance_sex_type_map, genetic_distance_index);
            EXIT(INPUT_DATA_ERROR);
        }
        
        // So, the batch file referenced a valid map!
        free(gdi); free(gdsm);
        return;
    } else if (genetic_distance_index == -1 &&
               ITEM_READ(Value_Genetic_Distance_Index)) {
        // This should not be set in the batch file, since we cannot change it...
        errorvf("A Value_Genetic_Distance_Index must not be specified as -1 in the batch file.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    
    if (genetic_distance_sex_type_map != -1) {
        // At this point no map has been specified (unknown), so no map type should be specified...
        errorvf("You cannot specify a Value_Genetic_Distance_SexTypeMap (%d) without specifying a genetic map (Value_Genetic_Distance_Index).\n",
                genetic_distance_sex_type_map);
        EXIT(INPUT_DATA_ERROR);
    }
    
    //
    // At this point 'should' be in either INTERACTIVE_INPUTMODE, or BATCH_INPUTMODE and no map has been specified...
    
    if (requires_genetic_map_p == 0 && gds == 0) {
#ifndef HIDESTATUS
        mssgvf("NOTE: no genetic map is available, and none is required based on the analysis type chosen.\n");
#endif /* HIDESTATUS */
        genetic_distance_index = -2; // choose 'None' for the user...
        genetic_distance_sex_type_map = NONE_AVAILABLE_GDMT;
        
        // We should only write to the batch file if we believe that this value has not already been
        // written. Should probably check, but at this point is should not have been.
        Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option = genetic_distance_index;
        batchf(Value_Genetic_Distance_Index);
        
        Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
        batchf(Value_Genetic_Distance_SexTypeMap);
        
        free(gdi); free(gdsm);
        return;
    }
    
    if (requires_genetic_map_p == 1 &&
        (gds == 1 ||
         (InputMode == BATCH_FILE_INPUTMODE && gds > 1))) {
            // A genetic map is required;
            // there is only one of them, or
            // we are in batch input mode, multiple maps are available, and we can't ask the user which one they want
            // so we always choose the first map and warn them.
#ifndef HIDESTATUS
            if (gds == 1) {
                mssgf("NOTE: Only one genetic map is available for use, so it will be used.\n");
            } else {
                mssgf("NOTE: No genetic map was specified in the batch file, so the first map will be used.\n");
            }
#endif /* HIDESTATUS */
            genetic_distance_index = gdi[0];
            genetic_distance_sex_type_map = gdsm[0];
            
            Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option = genetic_distance_index;
            batchf(Value_Genetic_Distance_Index);
            
            Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
            batchf(Value_Genetic_Distance_SexTypeMap);
            
            free(gdi); free(gdsm);
            return;
        }
    
    // We have done as much automatic 'fixes' that we can if in batch input mode, so we error if things
    // are not completely specified and need to be.
    if (InputMode == BATCH_FILE_INPUTMODE) {
        if (requires_genetic_map_p == 1) {
            errorvf("The Value_Genetic_Distance_Index and Value_Genetic_Distance_SexTypeMap must specified in the batch file.\n",
                    genetic_distance_sex_type_map);
            EXIT(INPUT_DATA_ERROR);
        } else {
            // BATCH_FILE_INPUTMODE and no map file required.
            // Here Value_Genetic_Distance_Index == -1 so we do not need to write this to the batch file...
            free(gdi); free(gdsm);
            return;
        }
    }
    
    //
    // At this point we should be in INTERACTIVE_INPUTMODE only...
    
    // create a format string that will accomodate the longest map name...
    sprintf(format_str, "%%c%%d) %%%ds: %%13s %%s\n", max_map_name_len);
    
    option_selected = 1; // pick the first option as the default
    
    // setup for counting the 'None' option if the analysis type selected by the user allows it...
    gds2 = gds + (requires_genetic_map_p == 0 ? 1 : 0) ;
    
    // Build the menu for the user...
    do {
        draw_line();
        // make the selected option the current option, if it is valid...
        if (option_selected < 1 || option_selected > gds2)
            printf("? Invalid option selected. It must be between 0 - %d\n", gds2);
        else {
            option = option_selected;
        }
        
        // Ask the user which option to use...
        printf("\nGenetic map selection menu:\n");
        printf("0) Done with this menu - please proceed\n");
        for (i = 0; i < gds; i++) {
            //printf("%c%d) %s: %s\t%s\n",
            printf(format_str,
                   ((i+1) == option) ? '*' : ' ', i+1,
                   EXLTop->MapNames[gdi[i]],
                   genetic_distance_map_type_string[gdsm[i]],
                   (EXLTop->map_functions[gdi[i]] == 'h') ? "Haldane" : "Kosambi"
                   );
        }
        // if the analysis does not require a map, allow the user to select 'None'...
        if (requires_genetic_map_p == 0)
            printf("%c%d) None\n",
                   ((i+1) == option) ? '*' : ' ', i+1);
        
        print_information_based_on_genetic_distance_sex_type_map(gdsm[option-1]);
        printf("Select from options 0 - %d > ", gds2);
        
        // Get the response from the user....
        fcmap(stdin, "%d", &option_selected); newline;
        
        // is the user done with this menu?
    } while (option_selected != 0);
    
    // NOTE: from above... the last option 'option == gds2' is the 'None' option,
    // but only when a genetic map is not required...
    if (requires_genetic_map_p == 0 && option == gds2) {
        // No genetic map was required, and the user selected the 'None' option...
        // NOTE that above we only give the user the 'None' option if a map file is not required
        // (e.g., 'requires_genetic_map_p == 0').
        genetic_distance_index = -2;
        genetic_distance_sex_type_map = NONE_AVAILABLE_GDMT;
        mssgvf("No genetic map will be used. Only analyses requiring\n only a physical map will be possible.\n");
    } else {
        // The user selected this map...
        genetic_distance_index = gdi[option-1];
        genetic_distance_sex_type_map = gdsm[option-1];
    }
    
    // Store the relevant distance and sex map according to the user's selection...
    Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option = genetic_distance_index;
    batchf(Value_Genetic_Distance_Index);
    
    Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
    batchf(Value_Genetic_Distance_SexTypeMap);
    
    free(gdi);
    free(gdsm);
}

/*
Physical Distance Map Selection [cpk]

Determine which physical distance index from EXLTop that we should use for subsequent
processing. This is necessary since the input map file may specify several possible maps.

The algorithm is roughly as follows:

========================================
Physical map selection menu
0) Done with this menu - please proceed.
* 1) MapP
  2) MapL
  3) MapB
  4) None
Select from options 0 - 4 >
========================================

RULES:
1) if the user has chosen an analysis option that requires a physical map (such as IQLS), then the option 'None' should not appear in the menu at all.
2) if the map is empty (e.g., all 0's) then do not present it to the user.
3) always select the first map by default.
 
NOTE:
 IQLS acquires physical position from (int)LTop->Marker[markers[locus]].pos_avg where
 LTop->map_distance_type == 'p'
*/

// For a list of valid analysis types, see common.h:enum analysis_type;
//
// For 'AnalysisOpt' in 'RequiresPhysicalMap', do not give the user the 'None' option when
// allowing them to choose the physical map.
void get_base_pair_position_index(ext_linkage_locus_top *EXLTop) {
    int i, requires_physical_map_p = 0, bpps = 0, bpps2, *bppi = NULL;
    int option, option_selected, valid_map_p = 0;
    
    requires_physical_map_p = AnalysisOpt->require_physical_map() ? 1 : 0;
    // Determine how many maps are available....
    // bppi holds the indicies of valid maps.
    bppi = CALLOC((size_t) EXLTop->MapCnt, int);
    for (i = 0, bpps = 0; i < EXLTop->MapCnt; i++)
        // only take physical maps that are not empty...
        if (EXLTop->map_functions[i] == 'p'
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
	    && EXLTop->valid_map_p[i][SEX_AVERAGED_MAP]
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
	    )
            bppi[bpps++] = i;
    
    // If no physical map is available, but one is required...
    if (bpps == 0 && requires_physical_map_p == 1) {
        errorvf("For the analysis type specified, a physical map is required,\n       but none are available in the Map file.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    
    // If a map has already been specified, make sure that is is a valid one...
    if (base_pair_position_index == -2) {
        // 'None' was specified for a map in the batch file...
        if (requires_physical_map_p == 1) {
            errorvf("For the analysis type specified, a physical map is required.\n");
            errorvf("However, Value_Base_Pair_Position_Index was specified as -2 (None).\n");
            EXIT(INPUT_DATA_ERROR);
        } else {
            // So, 'None' is OK...
            free(bppi);
            return;
        }
    } else if (base_pair_position_index >= 0) {
        // A map has been specifid, so we make sure that it is valid by searching for it
        // in the list created above and set the appropriate booleans if encountered...
        for (i=0; i<bpps; i++)
            if (base_pair_position_index == bppi[i]) {
                valid_map_p = 1;
                break;
            }
        
        if (!valid_map_p) {
            errorvf("The Value_Base_Pair_Position_Index specified (%d) does not reference a physical map in the map file.\n",
                    base_pair_position_index);
            EXIT(INPUT_DATA_ERROR);
        }
        
        // So, the batch file referenced a valid map!
        free(bppi);
        return;
    } else if (base_pair_position_index == -1 &&
               ITEM_READ(Value_Base_Pair_Position_Index)) {
        // This should not be set in the batch file, since we cannot change it...
        errorvf("A Value_Base_Pair_Position_Index must not be specified as -1 in the batch file.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    
    
    //
    // At this point 'should' be in either INTERACTIVE_INPUTMODE, or BATCH_INPUTMODE and no map has been specified...
    
    if (bpps == 0 && requires_physical_map_p == 0) {
#ifndef HIDESTATUS
        mssgf("NOTE: no physical map was available, and none is required based on the analysis type chosen.\n");
#endif /* HIDESTATUS */
        base_pair_position_index = -2; // choose 'None' for the user...
        
        Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option = base_pair_position_index;
        batchf(Value_Base_Pair_Position_Index);
        
        free(bppi);	
        return;
    }
    
    // If there is only one physical map then choose it, whether in batch mode or interactive mode,
    // whether required or not!
    if (bpps == 1) {
#ifndef HIDESTATUS
            if (bpps == 1) {
                mssgf("NOTE: Only one physical map is available for use, so it will be used.\n");
            }
#endif /* HIDESTATUS */
            
            base_pair_position_index = bppi[0];
            
            Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option = base_pair_position_index;
            batchf(Value_Base_Pair_Position_Index);
            
            free(bppi);
            return;
        }

    //
    // At this point we should be in INTERACTIVE_INPUTMODE, or in BATCH_INPUTMODE and there are more
    // than one maps and none has been specified in the batch file...
    
    // setup for counting the 'None' option if the analysis type selected by the user allows it...
    bpps2 = bpps + (requires_physical_map_p == 0 ? 1 : 0);

    // pick the first option as the default, of no physical map is required, then choose that as the default...
    //option_selected = (requires_physical_map_p == 0 ? bpps2 : 1);
    option_selected = 1; // pick the first option as the default

    // Build the manu for the user...
    do {
        draw_line();
        // make the selected option the current option, if it is valid...
        if (option_selected < 1 || option_selected > bpps2)
            printf("? Invalid option selected. It must be between 0 - %d\n", bpps2);
        else
            option = option_selected;
        // Ask the user which option to use...
        printf("Physical map selection menu:\n");
        printf("0) Done with this menu - please proceed\n");
        for (i = 0; i < bpps; i++) {
            printf("%c%d) %s\n",
                   ((i+1) == option) ? '*' : ' ', i+1,
                   EXLTop->MapNames[bppi[i]]
                   );
        }
        // if the analysis does not require a map, allow the user to select 'None'...
        if (requires_physical_map_p == 0)
            printf("%c%d) None\n",
                   ((i+1) == option) ? '*' : ' ', i+1);
        
        printf("Select from options 0 - %d > ", bpps2);
        
        // Get the response from the user....
        fcmap(stdin, "%d", &option_selected); newline;
        // is the user done with this menu?
    } while (option_selected != 0);
    
    // NOTE: from above... the last option 'option == gds2' is the 'None' option...
    // but only when a physical map is not required...
    if (requires_physical_map_p == 0 && option == bpps2) {
        // No physical map was required, and the user selected the 'None' option...
        base_pair_position_index = -2;
        mssgf("No physical map will be used. Only analysis requiring only\n a genetic map will be possible.");
    } else {
        base_pair_position_index = bppi[option-1];
    }
    
    if (InputMode == INTERACTIVE_INPUTMODE) {
        Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option = base_pair_position_index;
        batchf(Value_Base_Pair_Position_Index);
    }
    
    free(bppi);
}
