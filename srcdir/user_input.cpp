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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "common.h"
#include "typedefs.h"

#include "input_ops.hh"
#include "input.hh"

#include "batch_input_ext.h"
#include "cw_routines_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "reorder_loci_ext.h"
#include "grow_string_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "plink_ext.h"
#include "utils_ext.h"
#include "user_input_ext.h"
#include "vcftools/mega2_vcftools_interface.h"

#include "class_old.h"
#include "database_dump_ext.h"

#ifdef _WIN
#define R_OK 4
#define access(str,type) _access(str,type)
#endif


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

void batchfdb(int i) {
//    if (! database_dump && ! database_read) 
    if (database_dump)
        batchf(i);
}

/* prototypes */
int             ReOrderMenu(int num_chromo, int *chromsomes, int *selection);
int             analysis_menu1(analysis_type  *analysis);

void  menu1(file_format *infl_type,
            char **pedfl_name, char **locusfl_name,
            char **mapfl_name, char **pmapfl_name,
            char **input_path, char **omitfl_name,
            char **freqfl_name, char **penfl_name,
            char **auxfl_name, char **phefl_name,
            int *Untyped_ped_opt, int *Error_sim_opt,
            char **output_path, char **db_name,
            double *freq_mismatch_thresh);

void            define_affection_labels(linkage_ped_top *Top,
					analysis_type analysis);
void           set_missing_quant_input(linkage_ped_top *Top, analysis_type analysis);

int            gh_cov_selection(int num_select, int *trait_order,
				int *local_covariates,
				linkage_locus_top *LTop);
char           *untyped_ped_messg(int opt, char *messg);

void           get_genetic_distance_index(ext_linkage_locus_top *EXLTop);
void           get_base_pair_position_index(ext_linkage_locus_top *EXLTop);

static  void define_labels(linkage_ped_top *Top, int tr, char *affdata_str);

static void     untyped_ped_menu(int *opt);

/*===========prototypes============= */

INPUT_FORMAT_t Input_Format = in_format_mega2;
Input_Base *Input;

const char *INPUT_FORMAT_STR[] = {
     "Mega2 format with header",
     "Linkage format",
     "Linkage with Mega2 names file",
     "PLINK binary PED format (bed)",
     "PLINK PED format (ped)",
     "BCF format (bcf)",
     "VCF compressed format (vcf.gz)",
     "VCF format (vcf)",
     "IMPUTE2 GEN format (gen/impute2)",
     "IMPUTE2 BGEN format (bgen)",
//   "IMPUTE2 BGEN format2 (bgen)",
};
const char *INPUT_FORMAT_STR100 = "Traditional (4.6.1) format";

double Imputed_Info_Metric_Threshold = 0.30;

// see R_output.h & R_output.c write_plink_map_file...
int genetic_distance_index;
//cpk genetic_distance_map_type genetic_distance_sex_type_map; // see common.h
int genetic_distance_sex_type_map; // see common.h
int base_pair_position_index;

/*===========================================================================*/
/*
 *      Requires corresponding addition/change to
 *              input.hh: add to enum INPUT_FORMAT {                   | around line 44
 ?              batch_input.cpp: check_batch_items                     | around line 427
 *              user_input.cpp: INPUT_FORMAT_STR[]                     | around line  107
 ?              user_input.cpp: "if (choice_ == file_format_i)" case   | around line 1187
 *              user_input.cpp: "if (Input_Format == in_format_xxx)" case | around line 963
 */
Input_Base *createinput(INPUT_FORMAT in_format) {
    switch(in_format) {
    case in_format_mega2:
        return new Input_Mega2(in_format);
        break;
    case in_format_linkage:
        return new Input_Linkage(in_format);
        break;
    case in_format_extended_linkage:
        return new Input_Extended_Linkage(in_format);
        break;
    case in_format_binary_PED:
        return new Input_PED_Binary(in_format);
        break;
    case in_format_PED:
        return new Input_PED(in_format);
        break;
    case in_format_binary_VCF:
        return new Input_VCF_Binary(in_format);
        break;
    case in_format_compressed_VCF:
        return new Input_VCF_Compressed(in_format);
        break;
    case in_format_VCF:
        return new Input_VCF(in_format);
        break;
    case in_format_imputed:
        return new Input_Impute(in_format);
        break;
    case in_format_bgen:
        return new Input_BGEN(in_format);
        break;
    case in_format_bgen2:
        return new Input_BGEN2(in_format);
        break;
    case in_format_traditional:
        return new Input_Traditional(in_format);
        break;
    default:
        printf("Invalid input format selected. #%d\n", in_format+1);
//extern void         Exit(int arg, const char *file, const int line, const char *err);
        EXIT(DATA_INCONSISTENCY);
        return new Input_Traditional(in_format);
    }
}

/*===========================================================================*/
/* menu for user loci selection for the user's selection/reordering option
   an attempt to make the entry crash proof is made here */
int ReOrderMenu(int num_chromo, int *chromosomes, int *selection)
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

void menu0()
{
    int line_len = 40;
    char selectstr[16];
    int select = -1, set = -1;

    BatchValueGet(set, "Input_Database_Mode");
    if (InputMode == INTERACTIVE_INPUTMODE) {
        while (select != 0) {
            printf("              Mega2 %s database mode menu:\n", Mega2Version);
            draw_line();
            printf("0) Done with this menu - please proceed\n");
            int idx=1;

            printf("%c%d", idx == set ? '*' : ' ', idx);
            printf(") %-*s\n", line_len, "Select Mega2 database create mode");
            idx++;

            printf("%c%d", idx == set ? '*' : ' ', idx);
            printf(") %-*s\n", line_len, "Select Mega2 database read mode");
            idx++;

            printf("%c%d", idx == set ? '*' : ' ', idx);
            printf(") %-*s\n", line_len, "Select Mega2 database create \"then use\" mode");
            idx++;

            printf("%c%d", idx == set ? '*' : ' ', idx);
            printf(") %-*s\n", line_len, "Select Mega2 \"NO database\" legacy mode");
            idx++;

            printf("Select from options 0-%d> ", idx-1);

            fcmap(stdin, "%s", selectstr);
            newline;
            draw_line();
            sscanf(selectstr, "%d", &select);
            if (select && select < 5)
                set = select;
        }
        BatchValueSet(set, "Input_Database_Mode");
        batchf("Input_Database_Mode");
    }

    if (set == 1) {
        database_dump = 1;
        database_read = 0;
    } else if (set == 2) {
        database_dump = 0;
        database_read = 1;
    } else if (set == 3) {
        database_dump = 1;
        database_read = 1;
    } else if (set == 3) {
        database_dump = 0;
        database_read = 0;
    }

//  printf("dump %d, read %d\n", database_dump, database_read);

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

extern Missing_Value missing_value;
extern Missing_Value missing_values[];
extern int count_missing_values;

static int sort_analysis(const void *a, const void *b)
{
    return strcasecmp(analysis_list[*(const int *)a].option_name,
                      analysis_list[*(const int *)b].option_name);
}

int analysis_menu1(analysis_type  *analysis)
{
    char            choice[50], sub_prog[50];
    int             i, choice_=-1, toggle = 0;
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
        int *sort_array = CALLOC(count_analysis_list, int);
        for (i = 0; i < count_analysis_list; i++)
            sort_array[i] = i;

        for (;;) {
            printf("                 ANALYSIS MENU \n");
            draw_line();
            if (toggle == 0)
                printf("%17s %2d Sort by Analysis Name\n", "", 0);
            else
                printf("%17s %2d Sort by Analysis Number\n", "", 0);

            for (i = 0; i < num_per_screen; i++) {
                sprintf(option_col1, "%2d %s", i+1, analysis_list[sort_array[i]].option_name);
                if ((num_per_screen +i) < (count_analysis_list)) {
                    sprintf(option_col2, "%2d %s",
                            num_per_screen+i+1, analysis_list[sort_array[num_per_screen+i]].option_name);
                } else {
                    strcpy(option_col2, "");
                }
                printf("%-34s %-34s\n", option_col1, option_col2);
            }

            printf("\nSelect an option between 0-%2d > ", count_analysis_list);
            fcmap(stdin, "%s", choice); newline;
            choice_ = -1;
            sscanf(choice, "%d", &choice_);
            if (choice_ == 0) {
                if (toggle == 0) {
                    qsort((void *) sort_array, (size_t) count_analysis_list, sizeof (int), sort_analysis);
                } else {
                    for (i = 0; i < count_analysis_list; i++)
                        sort_array[i] = i;
                }
                toggle = !toggle;
            } else if (choice_ >= 1 && choice_ <= (count_analysis_list)) {
                choice_ = sort_array[choice_ - 1];
                if (analysis_list[choice_].analysis->is_enabled()) {
                    printf("    You have selected: %s\n", analysis_list[choice_].option_name);
                    break;
                } else
                    printf("Disabled choice \"%s\". Please select again.\n",
                           analysis_list[choice_].option_name);
            } else
                printf("Choice out of range. Please select again from options 0-%2d.\n",
                    count_analysis_list);
        }
        *analysis = analysis_list[choice_].analysis;
        /* note the analysis may be changed */
        if ((*analysis)->has_sub_options())
            (*analysis)->interactive_sub_prog_name_to_sub_option(analysis);
    }

    const char *key = (*analysis)->_missing_value_key;
    if (*key == 0)
        for (i = 0; i < count_missing_values; i++) {
            if (! strcmp((*analysis)->_name, (missing_values[i].analysis)->_name)) {
                missing_value = missing_values[i];
                break;
            }
	}
    else
        for (i = 0; i < count_missing_values; i++) {
            if (! strcmp(key, (missing_values[i].analysis)->_subname)) {
                missing_value = missing_values[i];
                break;
            }
	}
    if (i == count_missing_values) {
        errorvf("Internal error: analysis type %s/%s was not found in missing_values table.\n",
                (*analysis)->_name, (*analysis)->_subname);
        EXIT(SYSTEM_ERROR);
    }

    AllowUnmapped = (ALLOW_NO_CHR(*analysis) ? 1 : 0);

    (*analysis)->prog_name(ProgName);
    mssgvf("Analysis Class: %s.\n", ProgName);

    int subo = 0;
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
            mssgvf("Sub-option: %s.\n", sub_prog);
            subo++;
        }
    }
    if (subo) mssgvf("Analysis: %s.\n", ProgName);

    return 1;
}   /* end of analysis_menu */

static void thresh_string(double val, char *str)
{
    if (val >= LARGE) {
        sprintf(str, "No limit");
    } else {
        sprintf(str, "%g", val);
    }
    return;

}

static void missing_optional_keyword(batch_item_type *bi, const char *message)
{
#ifndef HIDESTATUS
    mssgvf("Keyword %s not in batch file, %s.\n", C(bi->keyword), message);
#endif
}
static void missing_optional_keyword(int batch_item, const char *message)
{
    missing_optional_keyword(&Mega2BatchItems[batch_item], message);
}

static void missing_mandatory_keyword(batch_item_type *bi)
{
    errorvf("Required keyword %s missing from batch file.\n",
            C(bi->keyword));
    EXIT(BATCH_FILE_ITEM_ERROR);
}
static void missing_mandatory_keyword(int batch_item)
{
    missing_mandatory_keyword(&Mega2BatchItems[batch_item]);
}

typedef struct fln {
    const char *title;
    int batch;
    Str item;
    const char *type;
    const char *typefx;
    const char *stat;
    const char *sfx1;
    const char *sfx2;
          char *name;
          char **nameref;
    int on;
    int specified;
} fln_t;

fln_t pedO  = {"Pedigree file:",   Input_Pedigree_File} ,  *pedo  = &pedO;
fln_t locO  = {"Locus file:",      Input_Locus_File},      *loco  = &locO;
fln_t mapO  = {"Map file:",        Input_Map_File},        *mapo  = &mapO;
fln_t pmapO = {"PLINK map file:",  Input_PLINK_Map_File},  *pmapo = &pmapO;
fln_t omitO = {"Omit file:",       Input_Omit_File},       *omito = &omitO;
fln_t freqO = {"Frequency file:",  Input_Frequency_File},  *freqo = &freqO;
fln_t penO  = {"Penetrance file:", Input_Penetrance_File}, *peno  = &penO;
fln_t auxO  = {"Aux file:",        Input_Aux_File},        *auxo  = &auxO;
fln_t pheO  = {"Phenotype file:",  Input_Phenotype_File},  *pheo  = &pheO;
fln_t infO  = {"Imputed Info file:", -1, "Input_Imputed_Info_File"};
fln_t *info  = &infO;
fln_t refO  = {"Reference Allele File:",  -1, "Reference_Allele_File"},  *refo  = &refO;

fln_t *fln_array[] = {pedo, loco, mapo, pmapo, omito, freqo, peno, auxo, pheo, refo, info, 0};

static void fln_init(fln_t *fln, const char *type, const char *typefx, const char *stat, const char *sfx1)
{
    if (fln->on) return;

    fln->on     = 1;

    fln->type   = type;
    fln->typefx = typefx;
    fln->stat   = stat;
    fln->sfx1   = sfx1;
    fln->sfx2   = 0;
}

static void fln_init(fln_t *fln, const char *type, const char *typefx, const char *stat, const char *sfx1, const char *sfx2)
{
    if (fln->on) return;

    fln->on     = 1;

    fln->type   = type;
    fln->typefx = typefx;
    fln->stat   = stat;
    fln->sfx1   = sfx1;
    fln->sfx2   = sfx2;
}

#define MAP_REQ 1
static void fln_init_mega2(int map_req) {
    if (map_req < 0) 
        ;
    else if (map_req)
        fln_init(mapo,  "Mega2", "map", "[required]", "map");
    else
        fln_init(mapo,  "Mega2", "map", "[optional]", "map");
    fln_init(omito, "Mega2", "omit", "[optional]", "omit");
    fln_init(freqo, "Mega2", "freq", "[optional]", "freq", "frequency");
    fln_init(peno,  "Mega2", "pen", "[optional]", "pen", "penetrance");
    fln_init(refo,  "Alleles", "", "[optional]", "ref", "reference");
}

#define PMAP_REQ 1
static void fln_init_plink(int map) {
    fln_init(pedo,  "PLINK", "ped", "[required]", "ped");
    if (map)
        fln_init(pmapo, "PLINK", "map", "[required]", "map");
    fln_init(pheo,  "PLINK", "phe", "[optional]", "phe");
}

static void fln_off() {
    fln_t **fln = fln_array;
    while (*fln) {
        (*fln)->on = 0;
        (*fln)->specified = 0;
        (*fln)->name[0] = 0;
        fln++;
    }
}

static void fln_extension_off()  {
    fln_t **fln = fln_array;
    while (*fln)
        (*fln++)->specified = 0;
}

static void fln_alloc(char **name) {
    if (*name == NULL)
        *name = CALLOC((size_t) FILENAME_LENGTH, char);
}

static void fln_alloc(char **name, fln_t *fln) {
    if (*name == NULL) {
        *name = CALLOC((size_t) FILENAME_LENGTH, char);
        fln->name = *name;
        fln->nameref = name;
    }
}

static void fln_free(fln_t *fln) {
    free(fln->name);
    fln->on = 0;
    *(fln->nameref) = NULL;
}

static void fln_free_not_present(fln_t *fln) {
    if (fln->specified == 0 || access(fln->name, F_OK) != 0) {
        free(fln->name);
        fln->on = 0;
        *(fln->nameref) = NULL;
    }
}

int  fln_stem = 0;
char fln_tmp[FILENAME_LENGTH];
char extension_name[FILENAME_LENGTH];
int  fln_col1, fln_col2, fln_col4, fln_col1234;

static void fln_col_sizes() {
    fln_col1 = fln_col2 = fln_col4 = fln_col1234 = 0;
    int l;
    fln_t **fln = fln_array, *f;
    while (*fln) {
        f = *fln++;
        if (f->on == 0) continue;

        l = (int)strlen(f->title);
        if (l > fln_col1) fln_col1 = l;
        l = (int)strlen(f->type) + (int)strlen(f->typefx);
        if (l > fln_col2) fln_col2 = l;
        l = (int)strlen(f->stat);
        if (l > fln_col4) fln_col4 = l;
    }

    fln_col1234 = fln_col1 + 2 + fln_col2 + 4 + fln_col4 + 1;
}

static int fln_print(fln_t *fln, int idx, int _i) {
    if (! fln->on) return 0;

    if (fln_stem)
        sprintf(fln_tmp, "(%s .%s)", fln->type, fln->typefx);
    else
        sprintf(fln_tmp, "(%s %s.)", fln->type, fln->typefx);

    if (fln->specified == 0) {
        if (fln_stem)
            sprintf(fln->name, "%s.%s", extension_name, fln->sfx1);
        else
            sprintf(fln->name, "%s.%s", fln->sfx1, extension_name);

        if (access(fln->name, F_OK) && fln->sfx2) { // try again
            if (fln_stem)
                sprintf(fln->name, "%s.%s", extension_name, fln->sfx2);
            else
                sprintf(fln->name, "%s.%s", fln->sfx2, extension_name);
            if (access(fln->name, F_OK)) { // done trying
                if (fln_stem)
                    sprintf(fln->name, "%s.%s", extension_name, fln->sfx1);
                else
                    sprintf(fln->name, "%s.%s", fln->sfx1, extension_name);
            } else 
          fln->specified = 1;
        } else
          fln->specified = 1;
    } // else
//        fln->specified = 1;

    printf("%2d) %-*s%-*s%-*s%s\n",
           idx, fln_col1+2, fln->title, fln_col2+4, fln_tmp, fln_col4+1, fln->stat,
           (access(fln->name, F_OK) ? "_" : fln->name));
    return _i;
}

static void fln_get(fln_t *fln, const char *title) {
    draw_line();
    printf("Please enter %s file name ('clear' to clear) > ", title);
    strcpy(fln_tmp, fln->name);
    fcmap(stdin, "%s", fln->name); newline;
    if (strcasecmp(fln->name, "clear") == 0) {
        strcpy(fln->name, "-");
    } else {
        if (access(fln->name, F_OK)) {
            printf("WARNING: Could not find file %s\n", fln->name);
            strcpy(fln->name, fln_tmp);
            fln->specified = false;
        } else
            fln->specified = true;
    }
}

static void fln_batchf(fln_t *fln) {
    if (*fln->name != 0) {
	if (fln->batch == -1) {
            strcpy(BatchItemGet(fln->item)->value.name, fln->name);
            batchf(fln->item);
        } else {
	    strcpy(Mega2BatchItems[/* ? */ fln->batch].value.name, fln->name);
	    batchf(fln->batch);
	}
    }
}

static void menu1_batch_set_files(file_format *infl_type,
                                  char **pedfl_name, char **locusfl_name,
                                  char **mapfl_name, char **pmapfl_name,
                                  char **input_path, char **omitfl_name,
                                  char **freqfl_name, char **penfl_name,
                                  char **auxfl_name, char **phefl_name,
                                  int plinkf, int xcf)
{
    int i;

    batch_item_type *bi;
    fln_t **fln = fln_array;
    while (*fln) {
        i = (*fln)->batch;
        if (i == -1) {
            bi = BatchItemGet((*fln)->item);
        } else {
            bi = BatchItemGet(i);
        }
        if (bi->items_read) {
            if (access(bi->value.name, F_OK) == 0) {
                strcpy((*fln)->name, bi->value.name);
            } else {
                errorvf("Could not find file or path %s named by keyword %s.\n",
                        bi->value.name,
                        C(bi->keyword));
                EXIT(FILE_NOT_FOUND);
            }
        } else {
            if (i == Input_Pedigree_File || 
                ((i == Input_Locus_File) && Input->req_locus_file) ||
                ((i == Input_Aux_File) && Input->req_aux_file) ) {
                missing_mandatory_keyword(i);
            } else {
                strcpy(fln_tmp, (*fln)->title);
                size_t len = strlen(fln_tmp);
                fln_tmp[len-1] = ' ';
                strcat(&fln_tmp[len], "assumed to be unspecified");
                missing_optional_keyword(bi, fln_tmp);
            }
            fln_free(*fln);
        }
        fln++;
    }

    if (Mega2BatchItems[/* 54 */ Input_Path].items_read) {
        if (access(Mega2BatchItems[/* 54 */ Input_Path].value.name, R_OK) == 0 &&
            is_dir(Mega2BatchItems[/* 54 */ Input_Path].value.name)) {
            strcpy(*input_path, Mega2BatchItems[/* 54 */ Input_Path].value.name);
        } else {
            errorvf("file path %s named by keyword %s is not a readable directory.\n",
                    Mega2BatchItems[/* 54 */ Input_Path].value.name,
                    C(Mega2BatchItems[/* 33 */ Input_Path].keyword));
            EXIT(FILE_NOT_FOUND);
        }
    } else {
        missing_optional_keyword(Input_Path,  "using default '.' (current directory)");
        strcpy(*input_path, ".");
    }
}

static void menu1_batch_set_outfiles(char **output_path, char **db_name)
{
    if (Mega2BatchItems[/* 33 */ Output_Path].items_read) {
        if (access(Mega2BatchItems[/* 33 */ Output_Path].value.name, W_OK) == 0 &&
            is_dir(Mega2BatchItems[/* 33 */ Output_Path].value.name)) {
           strcpy(*output_path, Mega2BatchItems[/* 33 */ Output_Path].value.name);
        } else {
            errorvf("file path %s named by keyword %s is not a writable directory.\n",
                    Mega2BatchItems[/* 33 */ Output_Path].value.name,
                    C(Mega2BatchItems[/* 33 */ Output_Path].keyword));
            EXIT(FILE_NOT_FOUND);
        }

    } else {
        missing_optional_keyword(Output_Path,  "using default '.' (current directory)");
        strcpy(*output_path, ".");
    }

    char *cp = *db_name;
    if (*cp == 0)
        BatchValueGet(cp, "DBfile_name");
}

static void menu1_batch_set_misc(int *Untyped_ped_opt, int *Error_sim_opt,
                                 double *freq_mismatch_thresh)
{

    /* untyped ped option */
    if (Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].items_read) {
        if (Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option >= 0) {
            *Untyped_ped_opt=Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option;
        } else {
            invalid_value_field(Input_Untyped_Ped_Option);
        }
    } else {
        missing_mandatory_keyword(Input_Untyped_Ped_Option);
    }

    /* Error sim option */
    if (Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].items_read) {
        if (tolower((unsigned char)Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt) == 'y' ||
            tolower((unsigned char)Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt) == 'n') {
            *Error_sim_opt= Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt;
        } else {
            invalid_value_field(Input_Do_Error_Sim);
        }
    }
    *Error_sim_opt = ((tolower(*Error_sim_opt) == 'y') ? 1:0);

    /* Threshold for comparing input and observed allele frequencies */
    if (Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].items_read) {
        if (Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue < 0) {
            invalid_value_field(AlleleFreq_SquaredDev);
        }
        *freq_mismatch_thresh = Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue;
    } else {
        missing_optional_keyword(AlleleFreq_SquaredDev, "using default = 'no limit'");
    }
#ifdef USER_UNKNOWN
    if (Mega2BatchItems[/* 42 */ Value_Missing_Allele].items_read) {
        strcpy(REC_UNKNOWN, Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name);
    }
#endif
}

static int menu1_show_misc(int *Untyped_ped_opt, int *Error_sim_opt,
                           double *freq_mismatch_thresh,
                           int err_i, int untyp_i, int _thresh_i, int compress_i,
                           int *choiceA, int idx, int line_len)
{
    char messg[100], float_str[11];

    printf("%2d) %-*s[ %s]\n", idx, line_len,
           "Simulate genotyping errors:", yorn[*Error_sim_opt]);
    choiceA[idx++] = err_i;

    printf("%2d) %-*s\n", idx, line_len,
           untyped_ped_messg(*Untyped_ped_opt, messg));
    choiceA[idx++] = untyp_i;

    if (_thresh_i) {
        thresh_string(*freq_mismatch_thresh, float_str);
        printf("%2d) %-*s%s\n", idx, line_len,
               "Allele frequency error measure threshold:", float_str);
        choiceA[idx++] = _thresh_i;
    }

    if (compress_i) {
        printf("%2d) %-*s%s\n", idx, line_len,
               "Maximum number of alleles per marker:",
               MARKER_SCHEME == 1 ? "2 alleles" : (MARKER_SCHEME == 2 ? "255 alleles" : (MARKER_SCHEME == 3 ? "256 or more alleles" : "???")) );
        choiceA[idx++] = compress_i;
    }

    return idx;
}

static int menu1_set_misc(int *Untyped_ped_opt, int *Error_sim_opt,
                          double *freq_mismatch_thresh,
                          int err_i, int untyp_i, int thresh_i, int compress_i,
                          int choice_)
{
    int ret = 1;

    if (choice_ == err_i) {
        *Error_sim_opt = TOGGLE(*Error_sim_opt);
    } else if (choice_ == untyp_i) {
        draw_line();
        untyped_ped_menu(Untyped_ped_opt);
    } else if (choice_ == thresh_i) {
        draw_line();
        printf("Please enter threshold value > ");
        fcmap(stdin, "%g", freq_mismatch_thresh); newline;
    } else if (choice_ == compress_i) {
        int ans;
        while (1) {
            fflush(stdout);
            draw_line();
            printf("              Mega2 %s Maximum alleles per marker menu:\n", Mega2Version);
            draw_line();
            printf("%1d) 2 alleles (2 bits/marker)\n", 1);
            printf("%1d) 255 alleles (2 bytes/marker)\n", 2);
            printf("%1d) 256 or more alleles (16 bytes/marker)\n", 3);
            printf("Select from options 1-3 > ");

            fcmap(stdin, "%d", &ans); newline;
            if (ans <= 3 && ans >= 1) {
                MARKER_SCHEME = ans;
                break;
            } else
                printf("MARKER_SCHEME allowed values are 1, 2 or 3\n");
        }
    } else
        ret = 0;
    return ret;
}

static void menu1_batch_save_misc(int *Untyped_ped_opt, int *Error_sim_opt,
                                 double *freq_mismatch_thresh)
{
    Mega2BatchItems[/* 4 */ Input_Untyped_Ped_Option].value.option = *Untyped_ped_opt;
    batchf(Input_Untyped_Ped_Option);

    Mega2BatchItems[/* 24 */ Input_Do_Error_Sim].value.copt = yorn[*Error_sim_opt][0];
    batchf(Input_Do_Error_Sim);

    Mega2BatchItems[/* 37 */ AlleleFreq_SquaredDev].value.fvalue = *freq_mismatch_thresh;
    batchf(AlleleFreq_SquaredDev);
}

void menu1(file_format *infl_type,
           char **pedfl_name, char **locusfl_name,
           char **mapfl_name, char **pmapfl_name,
           char **input_path, char **omitfl_name,
           char **freqfl_name, char **penfl_name,
           char **auxfl_name, char **phefl_name,
           int *Untyped_ped_opt, int *Error_sim_opt,
           char **output_path, char **db_name,
           double *freq_mismatch_thresh,
           char **reffl_name, int *strand_flip_opt)
{
    int            i, choice_ = -1;
    char           cchoice[10];
    int            exit_loop=0;

    int            ext_i=1, loc_i=2, ped_i=3, map_i=4, omit_i = 5, freq_i=6, pen_i=7, ref_i = 8;
    int            out_i=9, err_i=10, untyp_i=11, thresh_i=12, miss_i=13, _thresh_i;
    int            plink_args_i = 15, plink_phe_i = 16, plink_bed_i = 17;
    int            compress_i = 18, file_format_i = 19, vcf_args_i = 20;
    int            vcf_mak_i = 21, site_vcf_i = 22, site_bcf_i = 23, site_vcf_gz_i = 24, _aux_i = 0;
    int            db_file_i = 25, in_dir_i = 26, pmap_i = 27;
    int	           imputed_i = 28, inf_i = 29/*, flip_i = 30*/;
    int            idx, choiceA[31]; /* idx should be 1+ largest <>_i value (above)*/

    int            plinkf = 0, xcf = 0;

    char           PLINKArgs[FILENAME_LENGTH] = "";
    char           VCFArgs[FILENAME_LENGTH] = "";
    char           VCFMarkerAlternativeKey[FILENAME_LENGTH] = "";
    int            reset = 1, reset_extension = 1;

    char           *infofl_array = NULL;
    char          **infofl_name = &infofl_array;

    /* specifying chromosome or extension overrides previous specifications
       unless user specified each file separately */
    *infl_type = LINKAGE;

    *Untyped_ped_opt=2; /* Exclude any pedigree with 1 or less untyped people */
    *Error_sim_opt = 0;
    *freq_mismatch_thresh = LARGE;

    fln_alloc(output_path);
    fln_alloc(input_path);

    fln_alloc(pedfl_name,   pedo);
    fln_alloc(locusfl_name, loco);
    fln_alloc(mapfl_name,   mapo);
    fln_alloc(pmapfl_name,  pmapo);
    fln_alloc(omitfl_name,  omito);
    fln_alloc(freqfl_name,  freqo);
    fln_alloc(penfl_name,   peno);
    fln_alloc(auxfl_name,   auxo);
    fln_alloc(phefl_name,   pheo);
    fln_alloc(infofl_name,  info);
    fln_alloc(reffl_name,   refo);
    if (batchINPUTFILES) {

        Input = createinput(Input_Format);

        if (Input_Format == in_format_PED)
            plinkf = 1;
        else if (Input_Format == in_format_binary_PED) {
            strcpy(&mega2_input_file_type[BED][0],  "PLINK Bed file");
            plinkf = 1;
        } else if (Input_Format == in_format_binary_VCF) {
            strcpy(&mega2_input_file_type[BED][0],  "Binary Variant file");
            xcf = 1;
        } else if (Input_Format == in_format_compressed_VCF) {
            strcpy(&mega2_input_file_type[BED][0],  "Compressed Variant file");
            xcf = 1;
        } else if (Input_Format == in_format_VCF) {
            strcpy(&mega2_input_file_type[BED][0],  "Variant file");
            xcf = 1;
        } else if (Input_Format == in_format_imputed) {
            strcpy(&mega2_input_file_type[BED][0],  "IMPUTE2 GEN file");
        } else if (Input_Format == in_format_bgen || Input_Format == in_format_bgen2) {
            strcpy(&mega2_input_file_type[BED][0],  "IMPUTE2 BGEN file");
        }


        menu1_batch_set_files(infl_type, pedfl_name, locusfl_name,
                              mapfl_name, pmapfl_name, input_path, omitfl_name,
                              freqfl_name, penfl_name, auxfl_name, phefl_name,
                              plinkf, xcf);

        menu1_batch_set_outfiles(output_path, db_name);

        menu1_batch_set_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh);

        return;
    }

    Input = createinput(Input_Format);  // may be changed later but need something for now

    sprintf(*output_path, ".");
    sprintf(*input_path, ".");

    int line_len;
    while (!exit_loop) {
        if (reset) {
            fln_off();
            reset = 0;
            plinkf = 0;
            xcf = 0;
            if (Input_Format == in_format_traditional || Input_Format == in_format_mega2) {
                strcpy(extension_name, "01");

                fln_init(loco, "Mega2", "names", "[required]", "names");
                fln_init(pedo, "Mega2", "pedin", "[required]", "pedin");
                fln_init_mega2(MAP_REQ);

            } else if (Input_Format == in_format_linkage) {
                strcpy(extension_name, "01");

                fln_init(loco, "Linkage", "datain", "[required]", "datain");
                fln_init(pedo, "Linkage", "pedin", "[required]", "pedin");
                fln_init(mapo, "Mega2 simple", "map", "[required]", "map");
                fln_init(omito, "Mega2", "omit", "[optional]", "omit");

            } else if (Input_Format == in_format_extended_linkage) {
                strcpy(extension_name, "01");

                fln_init(loco, "Mega2", "hdrless names", "[required]", "names");
                fln_init(pedo, "Linkage", "pedin", "[required]", "pedin");
//              fln_init_mega2( MAP_REQ);
                fln_init(mapo, "Mega2 simple", "map", "[required]", "map");
                fln_init(omito, "Mega2", "omit", "[optional]", "omit");

            } else if (Input_Format == in_format_binary_PED) {
                plinkf = 1;
                PLINK_clr(binary_PED_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(extension_name, "plink");

                fln_init(pedo, "PLINK", "fam", "[required]", "fam");
                fln_init(pmapo, "PLINK", "bim", "[required]", "bim");
                fln_init(auxo, "PLINK", "bed", "[required]", "bed");
                auxo->title = "Binary data file:";
                strcpy(&mega2_input_file_type[BED][0],  "PLINK Bed file");
                _aux_i = plink_bed_i;
                fln_init_plink(  PMAP_REQ);
//              fln_init_mega2(! MAP_REQ);
                fln_init_mega2(-1);

            } else if (Input_Format == in_format_PED) {
                plinkf = 1;
                PLINK_clr(PED_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(extension_name, "plink");

                fln_init_plink(  PMAP_REQ);
//              fln_init_mega2(! MAP_REQ);
                fln_init_mega2(-1);

            } else if (Input_Format == in_format_binary_VCF) {
                xcf = 1;
                PLINK_clr(not_plink_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(VCFArgs, "--remove-indels");
                strcpy(extension_name, "study");

                fln_init(pedo, "PLINK", "fam", "[required]", "fam");
                fln_init_plink(0);
                fln_init(auxo, "binary", "bcf", "[required]", "bcf");
                auxo->title = "Variant file:";
                strcpy(&mega2_input_file_type[BED][0],  "Binary Variant file");
                _aux_i = site_bcf_i;

                fln_init_mega2(! MAP_REQ);

            } else if (Input_Format == in_format_compressed_VCF) {
                xcf = 1;
                PLINK_clr(not_plink_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(VCFArgs, "--remove-indels");
                strcpy(extension_name, "study");

                fln_init(pedo, "PLINK", "fam", "[required]", "fam");
                fln_init(auxo, "compressed", "vcf.gz", "[required]", "vcf.gz");
                auxo->title = "Variant file:";
                strcpy(&mega2_input_file_type[BED][0],  "Compressed Variant file");
                _aux_i = site_vcf_gz_i;
                fln_init_plink(! PMAP_REQ);
                fln_init_mega2(! MAP_REQ);

            } else if (Input_Format == in_format_VCF) {
                xcf = 1;
                PLINK_clr(not_plink_format);
                PLINK_str(PLINKArgs, FILENAME_LENGTH);
                strcpy(VCFArgs, "--remove-indels");
                strcpy(extension_name, "study");

                fln_init(pedo, "PLINK", "fam", "[required]", "fam");
                fln_init(auxo, "text", "vcf", "[required]", "vcf");
                auxo->title = "Variant file:";
                strcpy(&mega2_input_file_type[BED][0],  "Variant file");
                _aux_i = site_vcf_i;
                fln_init_plink(! PMAP_REQ);
                fln_init_mega2(! MAP_REQ);

            } else if (Input_Format == in_format_imputed) {
                strcpy(extension_name, "impute");

                fln_init(pedo, "IMPUTE2", "sample", "[required]", "sample");
                pedo->title = "Sample file:";
                fln_init(auxo, "IMPUTE2", "gen", "[required]", "gen", "impute2");
                auxo->title = "IMPUTE2 file:";
                strcpy(&mega2_input_file_type[BED][0],  "IMPUTE2 GEN file");
                _aux_i = imputed_i;
                fln_init(info, "IMPUTE2", "gen_info", "[optional]", "gen_info", "impute2_info");
//              fln_init_plink(! PMAP_REQ);
                fln_init_mega2(! MAP_REQ);

            } else if (Input_Format == in_format_bgen || Input_Format == in_format_bgen2) {
                strcpy(extension_name, "bgen");

                fln_init(pedo, "IMPUTE2", "sample", "[required]", "sample");
                pedo->title = "Sample file:";
                fln_init(auxo, "IMPUTE2", "bgen", "[required]", "bgen");
                auxo->title = "IMPUTE2 bgen file:";
                strcpy(&mega2_input_file_type[BED][0],  "IMPUTE2 BGEN file");
                _aux_i = imputed_i;
                fln_init(info, "IMPUTE2", "gen_info", "[optional]", "gen_info", "impute2_info");
//              fln_init_plink(! PMAP_REQ);
                fln_init_mega2(! MAP_REQ);
            }
            reset_extension = 1;
        }
        if (reset_extension) {
            reset_extension = 0;
            fln_extension_off();

            if (strcasecmp(extension_name, "clear")==0) {
                strcpy(extension_name, "-");
            }
        }

        fln_col_sizes();
        line_len = fln_col1234;

        printf("              Mega2 %s file input menu:\n", Mega2Version);
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        idx=1;
        printf("%2d) %-*s%s\n", idx, line_len,
               "Select input file format:", (Input_Format != in_format_traditional ) ?
               INPUT_FORMAT_STR[Input_Format] : INPUT_FORMAT_STR100);
        choiceA[idx++] = file_format_i;

        if (plinkf || xcf) {
            if (plinkf) {
                printf("%2d) %-*s%s\n", idx, line_len, "Enter PLINK parameters:", PLINKArgs);
                choiceA[idx++] = plink_args_i;
            } else if (xcf) {
                char tmp[2000];
                printf("%2d) %-*s%s\n", idx, line_len, "Enter VCF parameters:", VCFArgs);
                choiceA[idx++] = vcf_args_i;
                
                if (strcmp(VCFMarkerAlternativeKey,"") == 0) sprintf(tmp, "%s", "ID field");
                else sprintf(tmp, "%s sub-field of the INFO field", VCFMarkerAlternativeKey);
                printf("%2d) %-*s%s\n", idx, line_len, "Read marker names from the:", tmp);
                choiceA[idx++] = vcf_mak_i;
                
                printf("%2d) %-*s%s\n", idx, line_len, "Enter PLINK parameters:", PLINKArgs);
                choiceA[idx++] = plink_args_i;
            }
        }

        if (Input->req_stem_flag) {
            printf("%2d) %-*s%s\n", idx, line_len, "Input file stem:", extension_name);
            fln_stem = 1;

        } else {
            printf("%2d) %-*s%s\n", idx, line_len, "Input file suffix:", extension_name);
            fln_stem = 0;
        }
        choiceA[idx++] = ext_i;

        Input->GetOps()->do_menu_display(idx, line_len, choiceA);

        choiceA[idx] = fln_print(auxo, idx, _aux_i);
        if (choiceA[idx]) idx++;

        if (Input_Format == in_format_imputed) {
            choiceA[idx] = fln_print(info, idx, inf_i);
            if (choiceA[idx]) idx++;
	}

	// BUG? _aux_i is only initialized in certain cases...
        // No. auxo is off iff _aux_i is not initialized and fln_print just returns 0 
        //  w/o doing anything

        choiceA[idx] = fln_print(loco, idx, loc_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(pedo, idx, ped_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(pmapo, idx, pmap_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(pheo, idx, plink_phe_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(mapo, idx, map_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(omito, idx, omit_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(freqo, idx, freq_i);
        if (choiceA[idx]) idx++;

        choiceA[idx] = fln_print(peno, idx, pen_i);
        if (choiceA[idx]) idx++;



        printf("%2d) %-*s%s\n", idx, line_len,
               "Output Directory:",
               ((!strcmp(*output_path, "."))?"[ Current directory ]" : *output_path));
        choiceA[idx] = out_i;
        idx++;

        if (database_dump ^ database_read) {
            extern int db_exists_db();
            int db_exists = db_exists_db();

            printf("%2d) %-*s%s%14s\n", idx, line_len,
                   "SQLite3 Database file:", *db_name,
                   db_exists ? (database_dump ? "[overwrite]" : "[exists]") : "[create]");
            choiceA[idx] = db_file_i;
            idx++;
        }

        if(database_dump)
            choiceA[idx] = fln_print(refo, idx, ref_i);

        if (choiceA[idx]) idx++;

/*
 * to be turned on some day
        if (xcf) {
            printf("%2d) %-*s%s\n", idx, line_len,
                   "Variant File(s) Directory:",
                   (!strcmp(*input_path, ".")?"[ Current directory ]" : *input_path));
            choiceA[idx++] = in_dir_i;
        }
*/

        _thresh_i = (access(*freqfl_name, R_OK) == 0) ?  thresh_i : 0;
        idx = menu1_show_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh,
                              err_i, untyp_i, _thresh_i, compress_i,
                              choiceA, idx, line_len);

        printf("Select from options 0-%d > ", idx-1);

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
            if (Input->req_locus_file)
                {
                if (access(*locusfl_name, F_OK) != 0)   {
                    printf("ERROR: You must specify a locus file.\n");
                    exit_loop=0;
                }
            }
            if (access(*pedfl_name, F_OK) != 0) {
                printf("ERROR: You must specify a pedigree file.\n");
                exit_loop=0;
            }
            if (Input->req_map_file)
            {
                // There is a map file that is extracted from the VCF file
                // and presented to the user as 'VCF.p'. So, the map file is optional.
                if (access(*mapfl_name, F_OK) != 0 && access(*pmapfl_name, F_OK) != 0) {
                    printf("ERROR: You did not specify a map or PLINK map file.\n");
                    exit_loop=0;
                }
            }
            if (plinkf) {
                if (*PLINKArgs == 0) {
                    printf("ERROR: You did not specify any PLINK parameters.\n");
                    exit_loop=0;
                } else if (exit_loop) {
                    Mega2BatchItems[PLINK_Args].items_read = 1;
                    printf("NOTE: You have indicated the trait in the fam/ped file is named \"%s\" and\n",
                           PLINK.trait);
                    printf("      is %s trait.  You also have specified the --missing-phenotype is %g.\n",
                           PLINK.traitType ? "a quantitative" : "an affective",  PLINK.pheno_value);
                    printf("NOTE: If this is not what you intended, go back and edit menu line 2.\n");
                }

            }
            if (Input_Format == in_format_binary_PED && access(*auxfl_name, F_OK) != 0) {
                printf("ERROR: You did not specify a PLINK binary data file.\n");
                exit_loop=0;
            }
            if (xcf) {
                if (access(*auxfl_name, F_OK) != 0) {
                    printf("ERROR: You did not specify a Variant file.\n");
                    exit_loop=0;
                }
                //
                // So that we can use the already available VCFtools command line
                // checking. The decision was made to 'splice' the file option into
                // the args just like the user would invoke VCFtools.
                char *VCFArgs_w_file = CALLOC(FILENAME_LENGTH, char);
                strcpy(VCFArgs_w_file, VCFArgs);
                if (Input_Format == in_format_binary_VCF) {
                    strcat(VCFArgs_w_file, " --bcf ");
                } else if (Input_Format == in_format_compressed_VCF) {
                    strcat(VCFArgs_w_file, " --gzvcf ");
                } else if (Input_Format == in_format_VCF) {
                    strcat(VCFArgs_w_file, " --vcf ");
                }
                strcat(VCFArgs_w_file, *auxfl_name);
                
                if (VCFtools_process_cmd_line_w_file(VCFArgs_w_file) != -1) {
                    printf("ERROR: VCF argument list is invalid.  Please fix it.\n");
                    exit_loop=0;
                } else {
                    Mega2BatchItems[VCF_Args].items_read = 1;
                }
                free(VCFArgs_w_file);
            }
            if (PLINK.plink && PLINK.no_pheno == 0 && PLINK.trait[0] == 0) {
                printf("ERROR: You did not specify the pedigree file trait name.\n");
                printf("ERROR: Please choose option 2 and do so now.\n");
                exit_loop=0;
            }

            if (exit_loop == 0)
                draw_line();
            else {
                if (omito->specified == 0)
                    fln_free(omito);
                else if (access(*omitfl_name, F_OK) != 0)
                    fln_free(omito);
                else {
                    printf("NOTE: Marker untyping will take place according to the omit file.\n");
                    printf("      Please check the '%s' log file after MEGA2 is finished.\n",
                           Mega2Log);
                }

                fln_free_not_present(freqo);
                fln_free_not_present(peno);
                fln_free_not_present(auxo);
                fln_free_not_present(pheo);
                fln_free_not_present(mapo);
                fln_free_not_present(pmapo);
                fln_free_not_present(refo);
                fln_free_not_present(info);
            }

        } else if (choice_ == file_format_i) {
            int ans;
            while (1) {
                fflush(stdout);
                draw_line();
                printf("              Mega2 %s input file type menu:\n", Mega2Version);
                draw_line();

                int cnt = sizeof (INPUT_FORMAT_STR) / sizeof (char *);
                for (ans = 0; ans < cnt; ans++)
                    printf("%1d) %s\n", ans+1, INPUT_FORMAT_STR[ans]);
                printf("Select from options 1-%d > ", cnt);

                fcmap(stdin, "%d", &ans); newline;
                if (ans <= cnt && ans >= 1) {
                    Input_Format = (INPUT_FORMAT_t) (ans - 1);
                    if (Input != 0) delete Input;
		    Input = createinput(Input_Format);
                    break;
                } else
                    printf("allowed values are 1 - %d.\n", cnt);
            }
            reset = 1;

        } else  if (choice_ == ext_i) {
            if (fln_stem)
                printf("Please enter file stem > ");
            else
                printf("Please enter file extension > ");
            fcmap(stdin, "%s", extension_name); newline;
            reset_extension = 1;

        } else if (Input->GetOps()->do_menu_parse(choice_) ) {
            ; // work handled in do_menu_parse iff it returns 1

        } else if (choice_ == inf_i) {
            fln_get(info, "imputed info");

        } else if (choice_ == loc_i) {
            fln_get(loco, "locus");

        } else if (choice_ == ped_i) {
            fln_get(pedo, "pedigree");

        } else if (choice_ == map_i) {
            fln_get(mapo, "map");

        } else if (choice_ == pmap_i) {
            fln_get(pmapo, "PLINK map");

        } else if (choice_ == omit_i) {
            fln_get(omito, "omit");

        } else if (choice_ == freq_i) {   /* The 'Freq' file */
            fln_get(freqo, "frequency");

        } else if (choice_ == pen_i) {   /* The penetrance file */
            fln_get(peno, "penetrance");

        } else if (choice_ == plink_bed_i) {
            fln_get(auxo, "binary");

        } else if (choice_ == site_bcf_i) {
            fln_get(auxo, "binary");

        } else if (choice_ == site_vcf_gz_i) {
            fln_get(auxo, "compressed");

        } else if (choice_ == site_vcf_i) {
            fln_get(auxo, "text");

        } else if (choice_ == imputed_i) {
            fln_get(auxo, "text");

        } else if (choice_ == plink_phe_i) {
            fln_get(pheo, "phenotype");

        } else if (choice_ == ref_i) {
            printf("You can use an external reference panel to get a set of reference alleles.\n");
            printf("This process is described in the section called 'External Reference Allele Panel in the Database'\n");
            printf("in the Mega2 documentation.\n\n");
            printf("Reference panels are 3 column files of CHR POS REF that are then gzipped.\n");
            printf("They can be constructed by hand or using a shell script included with Mega2\n");
            printf("called GetRefAlleles.sh.  Additionally we provide a reference of 1000 genomes\n");
            printf("most recent build at ____________. \n\n");
            fln_get(refo, "reference allele");
            BatchValueSet(refo->name,"Reference_Allele_File");
            newline;
            draw_line();
            char buildname[255];

            Str rfile = refo->name;

            if(rfile.find("B37")!=std::string::npos || rfile.find("b37")!=std::string::npos)
                strcpy(buildname,"B37");
            else if(rfile.find("HG37")!=std::string::npos || rfile.find("hg37")!=std::string::npos)
                strcpy(buildname,"HG37");
            else if(rfile.find("B38")!=std::string::npos || rfile.find("b38")!=std::string::npos)
                strcpy(buildname,"B38");
            else if(rfile.find("HG38")!=std::string::npos || rfile.find("hg38")!=std::string::npos)
                strcpy(buildname,"HG38");
            else if(rfile.find("B19")!=std::string::npos || rfile.find("b19")!=std::string::npos)
                strcpy(buildname,"B19");
            else if(rfile.find("HG19")!=std::string::npos || rfile.find("hg19")!=std::string::npos)
                strcpy(buildname,"HG19");
            else {
                printf("Enter genome build for Reference File and Dataset > ");
                fcmap(stdin, "%s", buildname);
                newline;
            }
            printf("Genome build has been set to %s\n",buildname);
            Str hgbuild = buildname;

            BatchValueSet(hgbuild,"human_genome_build");
            batchf(BatchItemGet("human_genome_build"));



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
                IgnoreValue(fgets(y, 9, stdin)); newline; fflush(stdin);
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
        } else if (choice_ == db_file_i) {   /* The database file */
            draw_line();
            printf("Please enter SQLite3 database file name > ");
            fcmap(stdin, "%s", *db_name); newline;
            BatchValueSet(*db_name, "DBfile_name");

        } else if (choice_ == in_dir_i) {   /* The input directory */
            draw_line();
            printf("Please enter input variants directory name > ");
            fcmap(stdin, "%s", *input_path); newline;

            if (access(*input_path, F_OK)) {
                printf("WARNING: Could not find directory %s\n", *input_path);
                strcpy(*input_path, ".");
            } else if (! is_dir(*input_path)) {
                printf("WARNING: %s is not a directory.\n", *input_path);
                printf("Please specify a new or valid directory.\n");
                strcpy(*input_path, ".");
            } else if (access(*input_path, R_OK)) {
                printf("WARNING: %s is not a readable directory.\n", *input_path);
                printf("Please specify a new or valid directory.\n");
                strcpy(*input_path, ".");
            }

        } else if (choice_ == miss_i) {
            draw_line();
            printf("Please enter missing value indicator > ");
            fcmap(stdin, "%s", REC_UNKNOWN); newline;

        } else if (choice_ == plink_args_i) {
            PLINK_usage(xcf);
            while (1) {
                fflush(stdout);
                IgnoreValue(fgets(PLINKArgs, sizeof(PLINKArgs)-1, stdin)); newline;
                i = (int)strlen(PLINKArgs);
                if (PLINKArgs[i-1] == '\n') PLINKArgs[i-1] = 0;
                if (PLINKArgs[i-1] == '\r') PLINKArgs[i-1] = 0;
                if (PLINK_args(PLINKArgs, xcf)) break;
                printf("Enter     UPDATED PLINK parameters:  %s\n", PLINKArgs);
            }

        } else if (choice_ == vcf_args_i) {

	    char new_VCFArgs[FILENAME_LENGTH];
	    VCFtools_printf_supported_cmd_line_options();
            while (1) {
                printf("\nCurrent VCF parameters:  %s\n", VCFArgs);
                // 3) When doing the vcftools filtering, it wasn't clear to me at first which set
                // of IDs/positions one should use.  Presumably you must use the IDs/positions as
                // given in the VCF file itself.  So we have to filter individuals by VCF sample IDs,
                //and positions by VCF positions.
                printf("Note that filtering must be done using position information and sample IDs from\n");
                printf("the input VCF file. \n");
                printf("Enter new VCF parameters: ");
                strcpy(new_VCFArgs, VCFArgs);
                fflush(stdout);
                IgnoreValue(fgets(new_VCFArgs, sizeof(new_VCFArgs)-1, stdin)); newline;
                i = (int)strlen(new_VCFArgs);
                if (new_VCFArgs[i-1] == '\n') new_VCFArgs[i-1] = 0;
                if (new_VCFArgs[i-1] == '\r') new_VCFArgs[i-1] = 0;
                if (strlen(new_VCFArgs) == 0) break;
                if (VCFtools_process_cmd_line_wo_file(new_VCFArgs) == -1) {
                    strcpy(VCFArgs, new_VCFArgs);
                    break;
                }
            }

        } else if (choice_ == vcf_mak_i) {
            char select[100];
            int ans;
            strcpy(VCFMarkerAlternativeKey, "");
            while (1) {
                fflush(stdout);
                draw_line();
                printf("Read the marker names from:\n");
                printf("1) the ID field\n");
                printf("2) an INFO sub-field\n");
                printf("Select from options 1-2 > ");
                fcmap(stdin, "%s", select);
                sscanf(select, "%d", &ans);
                if (ans < 1 || ans > 2) {
                    printf("Please enter a 1, or 2.\n");
                    continue;
                } else if (ans == 2) {
                    printf("\n Please input the name of the INFO sub-field 'key'\n that contains the marker names: ");
                    fcmap(stdin, "%s", VCFMarkerAlternativeKey); newline;
                }
                break;
            }
        } else if (menu1_set_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh,
                                  err_i, untyp_i, thresh_i, compress_i, choice_)) {
                       // above function looks for match and does action
        } else {
            printf("Invalid option %s, select from options 0-%d.\n", cchoice, idx-1);
        }

        if (choice_) draw_line();
    }

    if (InputMode == INTERACTIVE_INPUTMODE) {

        Mega2BatchItems[/* 53 */ Input_Format_Type].value.option = Input_Format;
        batchf(Input_Format_Type);

        fln_t **fln = fln_array;
        while (*fln) {
            if ((*fln)->on)
                fln_batchf(*fln);
            fln++;
        }

        strcpy(Mega2BatchItems[/* 33 */ Output_Path].value.name, *output_path);
        batchf(Output_Path);
        strcpy(Mega2BatchItems[/* 54 */ Input_Path].value.name, *input_path);
        batchf(Input_Path);

#ifdef USER_UNKNOWN
        strcpy(Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name, REC_UNKNOWN);
        batchf(Value_Missing_Allele);
#endif
        // Handle the batch file item only if the user input something...
        // Note: these PLINKArgs and friends are stack variables so we don't have to check for != NULL...
        if (strlen(PLINKArgs) > 0) {
            strcpy(Mega2BatchItems[/* 43 */ PLINK_Args].value.name, PLINKArgs);
            Mega2BatchItems[/* 43 */ PLINK_Args].items_read = 1;
            batchf(PLINK_Args);
        }
        if (strlen(VCFArgs) > 0) {
            strcpy(Mega2BatchItems[/* 56 */ VCF_Args].value.name, VCFArgs);
            if (xcf) batchf(VCF_Args);
        }
        if (strlen(VCFMarkerAlternativeKey) > 0) {
            strcpy(Mega2BatchItems[/* 57 */ VCF_Marker_Alternative_INFO_Key].value.name, VCFMarkerAlternativeKey);
            if (xcf) batchf(VCF_Marker_Alternative_INFO_Key);
        }

        Input->GetOps()->do_menu2batch();

        menu1_batch_save_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh);

        Mega2BatchItems[/* 52 */ Value_Marker_Compression].value.option = MARKER_SCHEME;
        batchf(Value_Marker_Compression);

    }

    if (! Input->req_locus_file) fln_free(loco);
}

void menu1a(int *Untyped_ped_opt, int *Error_sim_opt,
            char **output_path, char **db_name,
            double *freq_mismatch_thresh, int *strand_flip_opt)
{
    int            choice_ = -1;
    char           cchoice[10];
    int            exit_loop=0;

    int            db_i=7, out_i=8, err_i=9, untyp_i=10, thresh_i=11, miss_i=12, _thresh_i;
    int            compress_i = 17;
    int            flip_i = 30;
    int            idx, choiceA[31]; /* idx should be 1+ largest <>_i value (above)*/
    extern int     db_exists_db();
    int            db_exists = 0;
    extern char    DBfile[255];
    char          *fn = DBfile;

    *Untyped_ped_opt=2; /* Exclude any pedigree with 1 or less untyped people */
    *Error_sim_opt = 0;
    *freq_mismatch_thresh = LARGE;

    fln_alloc(output_path);

    if (batchINPUTFILES) {

        menu1_batch_set_outfiles(output_path, db_name);

        if (*fn == 0)
            BatchValueIfSet(fn,   "DBfile_name");

        menu1_batch_set_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh);

        return;
    }

    sprintf(*output_path, ".");
    if (*fn == 0)
        BatchValueGet(fn,   "DBfile_name");

    *strand_flip_opt=0;

/*
    printf("\n");
    printf("Mega2 will read the input from the database specified in menu item 2 below.\n");
    printf("To instead read from a set of input files and create a new database, hit 'q'\n");
    printf("to exit Mega2 and then restart Mega2 with the --DBdump command line argument.\n");
    printf("\n");
*/
    int line_len = 45;
    while (!exit_loop) {

        printf("              Mega2 %s Database input menu:\n", Mega2Version);
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        idx=1;

        printf("%2d) %-*s%s\n", idx, line_len,
               "Output Directory:",
               ((!strcmp(*output_path, "."))?"[ Current directory ]" : *output_path));
        choiceA[idx] = out_i;
        idx++;

        db_exists = db_exists_db();
        printf("%2d) %-*s%s%10s\n", idx, line_len, "Database filename:", fn, 
               db_exists ? "[exists]" : "[create]");
        choiceA[idx] = db_i;
        idx++;

        _thresh_i = 1 ?  thresh_i : 0;
        idx = menu1_show_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh,
                              err_i, untyp_i, _thresh_i, 0 /*compress_i*/,
                              choiceA, idx, line_len);

        printf("%2d) %-*s[ %s]\n", idx, line_len,
               "Simulate genotyping errors:", yorn[*strand_flip_opt]);
        choiceA[idx++] = flip_i;

        printf(" q) %-*s\n", line_len, "Exit Mega2.");

        printf("Select from options 0-%d or q> ", idx-1);

        while (1) {
            fcmap(stdin, "%s", cchoice); newline;
            if (!strcmp(cchoice, "q")) {
                close_logs();
                exit(0);
            }
            choice_ = -1;
            sscanf(cchoice, "%d", &choice_);
            if (choice_ < idx) {
                if (choice_ > 0) choice_ = choiceA[choice_];
                break;
            } else
                printf("Select from options 0-%d > ", idx-1);
        }

        if (choice_ == 0) {
            exit_loop = 1;

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
                IgnoreValue(fgets(y, 9, stdin)); newline; fflush(stdin);
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

        } else if (choice_ == db_i) {   /* The database file */
            draw_line();
            printf("Please enter SQLite database filename > ");
            fcmap(stdin, "%s", fn); newline;
            
/*
            if (access(*output_path, F_OK)) {
                char y[100];
            } else if (access(*output_path, W_OK)) {
                    printf("WARNING: %s is not a writable directory.\n", *output_path);
                    printf("Please specify a new or valid directory.\n");
            }
*/
        } else if (choice_ == miss_i) {
            draw_line();
            printf("Please enter missing value indicator > ");
            fcmap(stdin, "%s", REC_UNKNOWN); newline;

        } else if (choice_ == flip_i){
            *strand_flip_opt = (*strand_flip_opt + 1) % 2;
        }

        else if (menu1_set_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh,
                                  err_i, untyp_i, thresh_i, compress_i, choice_)) {
                       // above function looks for match and does action
        } else {
            printf("Invalid option %s, select from options 0-%d.\n", cchoice, idx-1);
        }

        if (choice_) draw_line();
    }

    if (InputMode == INTERACTIVE_INPUTMODE) {

        BatchValueSet(fn, "DBfile_name");

        strcpy(Mega2BatchItems[/* 33 */ Output_Path].value.name, *output_path);
        batchf(Output_Path);

#ifdef USER_UNKNOWN
        strcpy(Mega2BatchItems[/* 42 */ Value_Missing_Allele].value.name, REC_UNKNOWN);
        batchf(Value_Missing_Allele);
#endif
        // Handle the batch file item only if the user input something...
        // Note: these PLINKArgs and friends are stack variables so we don't have to check for != NULL...

        menu1_batch_save_misc(Untyped_ped_opt, Error_sim_opt, freq_mismatch_thresh);
    }

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
    IgnoreValue(fgets(affdata_str, MAX_NAMELEN-1, stdin)); newline;
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

    if (batchAFFVALUE && Mega2BatchItems[/* 18 */ Value_Affecteds].items_read) {
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
                if (!strcmp(Top->LocusTop->Pheno[ml_traits[tr]].TraitName,
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
                printf(" %d) %s [%s]\n", i+1, Top->LocusTop->Pheno[ml_traits[i]].TraitName,
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
                    Top->LocusTop->Pheno[ml_traits[i]].TraitName,
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
int check_quant_phenotype_data_has_value(linkage_ped_top *Top, double value) {
    int i, i1, ped, entry;

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
#if 0 /* defunct code ... see menu_missing_value.cpp */
    int select=-1;
    double value;
    char missingq[200], *end;
    
    // Some analysis options have a 'hard wired' missing quantitative value.
    // For these options that value will be used later in the program so there
    // is no point in asking the user now or checking the value of the batch
    // file item Value_Missing_Quant_On_Output.
    if (!analysis->output_quant_can_define_missing_value()) {
        if (analysis->output_quant_default_value() != NULL) {
            const char *missingq = analysis->output_quant_default_value();
            char *end;
            double value;
            if (analysis->output_quant_must_be_numeric()) {
                value = strtod((const char *)missingq, &end);
                    // Did a successful numeric conversion take place?
                if (strlen(missingq) == 0 || strlen(end) != 0 || errno == ERANGE) {
                    errorvf("Analysis type '%s' requires quantitative values to be numeric.\n",
                            analysis->_name);
                    errorvf("But converting specified value '%s' to a floating point number failed.\n", missingq);
                    EXIT(DATA_TYPE_ERROR);
                } else if (check_quant_phenotype_data_has_value(Top, value)) {
                    errorvf("The default missing quantitative value is among the values found in the set of quantitative values");
                    EXIT(DATA_TYPE_ERROR);
                }
            }
            strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name, missingq);
            Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read = 1;
        }
        return;
    }

    if (!Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read) {
        
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
        
        Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read = 1;
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
                    errorvf("Select an output missing quantitative trait value that is not in your data.\n",
                            analysis->_name);
                    errorvf("The output missing value '%lf' exists in the quantitative trait data.\n", value);
                    EXIT(BATCH_FILE_ITEM_ERROR);
                }
        }
    }

    mssgvf("NOTE: The Missing QTL value on output will be assigned as '%s'.\n",
           Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
#endif
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

    // At this point MissingQuant is either QMISSING (e.g., "NA") or some other number...
    // Now determined by menu_value_missing.cpp:Value_Missing_get()
    
    if (Top == NULL) return;
    
    //
    // Since MissingQuant has just been firmly established, this is the point in the code
    // where we convert any "NA"s (QMISSINGs) that we have read to the value of MissingQuant.
    //
    // NOTE: It is pointless to do this if 'MissingQuant == QMISSING' as would be the case if we
    // were reading Mega2 files with NA used as the missing value.
    if (fabs(MissingQuant - QMISSING) > EPSILON) {
        if (HasQuant) {
            for (i1 = 0; i1 < num_traits ; i1++) {
                SKIP_TRI(i1)
                i = global_trait_entries[i1];
                if (Top->LocusTop->Locus[i].Type == QUANT) {
                    for(ped = 0; ped < Top->PedCnt; ped++) {
                        if (Top->pedfile_type == POSTMAKEPED_PFT) {
                            linkage_ped_tree *Ped;
                            if (database_dump) 
                                Ped = Top->PedBroken; // larger set of person id's 
                            else
                                Ped = Top->Ped;
                            for(entry = 0; entry < Ped[ped].EntryCnt; entry++) {
                                if (fabs(Ped[ped].Entry[entry].Pheno[i].Quant - QMISSING) <= EPSILON) {
                                    /* Set this to the proper missing quant value */
                                    Ped[ped].Entry[entry].Pheno[i].Quant = MissingQuant;
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

char *untyped_ped_messg(int opt, char *untyped_mssg)
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
                 "Exclude any pedigree with 0 genotyped people") :
         sprintf(untyped_mssg,
                 "Exclude any pedigree with %d or less genotyped people",
                 opt - 4));
        break;
    default:
        /*    for taking case of negative numbers */
        printf("Unknown option %d.\n", opt);
        break;
    }
    return &(untyped_mssg[0]);
}

static void untyped_ped_menu(int *opt)
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

int gh_cov_selection(int num_select, int *traits,
                     int *local_covariates, linkage_locus_top *LTop)
{

    int t, n, numq, num_cov=0, done=0, *quants, *selected;
    char cselect[200]="", *str_p, cnum[10];


    numq=0;
    for(t=0; t<num_select; t++) {
        if (LTop->Locus[traits[t]].Type == QUANT) numq++;
    }

    if (numq == 0) return 0;

    selected = CALLOC((size_t) numq, int);
    quants = CALLOC((size_t) numq, int);
    for(t=0; t<numq; t++) selected[t]=0;

    numq=0;
    for(t=0; t<num_select; t++) {
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
                   t+1, LTop->Pheno[traits[quants[t]]].TraitName);
        }
        draw_line();
        printf("\nEnter only 'e' if no covariates are desired.\n");
        printf("Enter trait numbers ('e' to terminate) > ");
        fflush(stdout);
        IgnoreValue(fgets(cselect, 199, stdin));
        newline;
        str_p=&(cselect[0]);
        while(*str_p != '\0') {
            if (num_cov >= num_select) {
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
                if (atoi(cnum) >= 1 && atoi(cnum) <= num_select) {
                    selected[atoi(cnum)-1]=1;
                    num_cov++;
                } else {
                    printf("Invalid selection %s, please re-select.\n",
                           cnum);
                    for (t=0; t<num_select; t++) {
                        selected[t]=0;
                    }
                    strcpy(cnum, ""); strcpy(cselect, "");
                    str_p = &(cselect[0]);
                    break;
                }
            } else {
                printf("Non-numeric characters in input, please re-select.\n");
                for (t=0; t<num_select; t++) {
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
                local_covariates[quants[t]]=1;
            }
        }

        sprintf(err_msg, "QTL/covariate list:");
        mssgf(err_msg);
        for (t=0; t<numq; t++) {
            sprintf(err_msg, "%d) %s %s", t+1,
                    LTop->Pheno[traits[quants[t]]].TraitName,
                    (local_covariates[quants[t]]? "[covariate]":""));
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

    case 6:
        sprintf(ind_id_messg, "Original pre makeped person");
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
    int choice[10];
    int i;


    if (get_disp_log == 3) {
        /* log */
        sprintf(err_msg, "Person id in output pedigree file = %s",
                ind_id_choice_messg(current_opt_val, &(ind_id_messg[0])));
        mssgf(err_msg);
        log_line(mssgf);
        if (InputMode == INTERACTIVE_INPUTMODE) {
            batchf("ID_person");
        }
        return 0;
    }

    if (get_disp_log==2) {
        printf(" %d) %sPerson id in output pedigree file:   %s\n",
               item_number, (item_number == 9) ? " " : "",
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
            choice[i++] = 1;

            if (has_orig) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(2, &(ind_id_messg[0])));
                choice[i++] = 2;
            }

            if (has_uniq) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 3) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(3, &(ind_id_messg[0])));
                choice[i] = 3;
            } else {
                printf("%c%d) %s\n",
                       ((current_opt_val == 4) ? '*' : ' '),
                       i,
                       ind_id_choice_messg(4, &(ind_id_messg[0])));
                choice[i] = 4;
            }
            i++;

            printf("%c%d) %s\n",
                   ((current_opt_val == 5) ? '*' : ' '),
                   i,
                   ind_id_choice_messg(5, &(ind_id_messg[0])));
            choice[i++] = 5;

            printf("%c%d) %s\n",
                   ((current_opt_val == 6) ? '*' : ' '),
                   i,
                   ind_id_choice_messg(6, &(ind_id_messg[0])));
            choice[i] = 6;

            printf("Select from options 0 - %d > ", i);
            fcmap(stdin, "%s", opt_);
            new_opt = atoi(opt_);

            if (new_opt < 0 || new_opt > i) {
                warn_unknown(opt_);
                continue;
            } else if (new_opt > 0)
                current_opt_val = choice[new_opt];
#if 0
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
            case 5:
                current_opt_val = 6;
                break;
            }
#endif
        }

    }
    BatchValueSet(current_opt_val, "ID_person");
    return current_opt_val;
}

static char *ped_id_choice_messg(int opt, char *id_messg)
{
    switch(opt) {
    case 1:
        if (basefile_type == POSTMAKEPED_PFT) {
            sprintf(id_messg, "Pedigree number");
            break;
        }

    case 2:
        if (basefile_type == POSTMAKEPED_PFT) {
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

    case 6:
        sprintf(id_messg, "Original pre makeped pedigree");
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
    int choice[10];
    int new_opt = -1;
    int i;

    if (get_disp_log == 3) {
        /* log */
        sprintf(err_msg, "Pedigree id in output pedigree file = %s",
                ped_id_choice_messg(current_opt_val, &(id_messg[0])));
        mssgf(err_msg);
        log_line(mssgf);
        if (InputMode == INTERACTIVE_INPUTMODE) {
            batchf("ID_pedigree");
        }
        return 0;
    }

    if (get_disp_log == 2) {
        printf(" %d) Pedigree id in output pedigree file: %s\n",
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
                choice[i++] = 5;

            } else if (pedfile_type == POSTMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 1) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(1, &(id_messg[0])));
                choice[i++] = 1;
            }

            if (analysis == TO_NUKE) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 4) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(4, &(id_messg[0])));
                choice[i++] = 4;

            } else if (has_orig && pedfile_type == POSTMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(2, &(id_messg[0])));
                choice[i++] = 2;

            } else if (pedfile_type == PREMAKEPED_PFT) {
                printf("%c%d) %s\n",
                       ((current_opt_val == 2) ? '*' : ' '),
                       i,
                       ped_id_choice_messg(2, &(id_messg[0])));
                choice[i++] = 2;
            }
            printf("%c%d) %s\n",
                   ((current_opt_val == 3) ? '*' : ' '),
                   i,
                   ped_id_choice_messg(3, &(id_messg[0])));
            choice[i++] = 3;

            printf("%c%d) %s\n",
                   ((current_opt_val == 6) ? '*' : ' '),
                   i,
                   ped_id_choice_messg(6, &(id_messg[0])));
            choice[i] = 6;

            printf("Select from options 0 - %d > ", i);
            fcmap(stdin, "%s", opt_);
            new_opt = atoi(opt_);

            if (new_opt < 0 || new_opt > i) {
                warn_unknown(opt_);
                continue;
            } else if (new_opt > 0)
                current_opt_val = choice[new_opt];
#if 0
            switch(new_opt) {
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
            case 4:
                current_opt_val = 6;
                break;
            }
#endif
        }
    }
    if (current_opt_val == 1 && pedfile_type == PREMAKEPED_PFT) {
        current_opt_val = 2;
    }
    BatchValueSet(current_opt_val, "ID_pedigree");
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

    BatchValueIfSet(OrigIds[1],  "ID_pedigree");
    BatchValueIfSet(OrigIds[0],  "ID_person");

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
//  requires_genetic_map_p = 1;
    
    //
    // Here we get a list of the valid maps from the structure that holds the parsing of
    // the map file. This is needed to both make sure an existing value specifies a real
    // map in the map file, and to present the user a list of maps for them to choose
    // from in interactive input mode.
    
    // NOTE: For Mega2 map files, indexes are assigned to maps in a map file by reading
    // it from left to right in order of the first occurrance of a map-name. It is an error
    // for the map-name to be associated with more than one map-type. It is an error for the
    // map-name to have a duplicate sex-designator.
    //
    // Mega2 maps are defined as follows:
    // map-designator == map-name . map-type. sex-designator
    // where: map-type is one of <h, k p>; sex-designator is one of <a, m, f>
    
    // gdi holds the indicies of valid genetic distance maps.
    gdi = CALLOC((size_t)(3*EXLTop->MapCnt), int);
    // gdsm tells whether the map is: sex-averaged (a), sex-specific (m,f), or female (f)
    gdsm = CALLOC((size_t)(3*EXLTop->MapCnt), genetic_distance_map_type);
    // gds is the index into the gdsm (map) array which specifies a valid genetic map (Index & SexType).
    // For example, if the Mega2 map file looks like this...
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
        errorvf("For the analysis type specified, a genetic map is required,\nbut none are available in the map file.\n");
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
        batchfdb(Value_Genetic_Distance_Index);
        
        Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
        batchfdb(Value_Genetic_Distance_SexTypeMap);
        
        free(gdi); free(gdsm);
        return;
    }
    
    if (// requires_genetic_map_p == 1 &&
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
            batchfdb(Value_Genetic_Distance_Index);
            
            Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
            batchfdb(Value_Genetic_Distance_SexTypeMap);
            
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
    sprintf(format_str, "%%d) %%%ds: %%13s %%s\n", max_map_name_len);
    option_selected = 1; // pick the first option as the default
    
    // setup for counting the 'None' option if the analysis type selected by the user allows it...
    gds2 = gds + (requires_genetic_map_p == 0 ? 1 : 0) ;
    
    // Build the menu for the user...
    do {
        draw_line();
        // make the selected option the current option, if it is valid...
        if (option_selected < 1 || option_selected > gds2) {
            printf("? Invalid option selected. It must be between 0 - %d\n", gds2);
            option = 1; // as it says above ... pick the first option as the default
        } else {
            option = option_selected;
        }
        
        // Ask the user which option to use...
        printf("\nGenetic map selection menu:\n");
        printf("0) Done with this menu - please proceed\n");
        for (i = 0; i < gds; i++) {
            //printf("%c%d) %s: %s\t%s\n",
            printf("%c", (((i+1) == option) ? '*' : ' '));
            printf(format_str,
                   i+1,
                   EXLTop->MapNames[gdi[i]],
                   genetic_distance_map_type_string[gdsm[i]],
                   ((EXLTop->map_functions[gdi[i]] == 'h') ? "Haldane" : "Kosambi")
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
    batchfdb(Value_Genetic_Distance_Index);
    
    Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
    batchfdb(Value_Genetic_Distance_SexTypeMap);
    
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
void  check_both_pos_index() {

    if (base_pair_position_index == -2 && AnalysisOpt->require_physical_map() ) {
        errorvf("For the analysis type specified, a physical map is required.\n");
        errorvf("However, none is in the database (value specified as -2).\n");
        EXIT(INPUT_DATA_ERROR);
    }

    int requires_genetic_map_p = AnalysisOpt->allow_no_genetic_map() ? 0 : 1;
    // always require a genetic map for now since genetic maps are entangled in the code...
//  requires_genetic_map_p = 1;

    if (genetic_distance_index == -2 && requires_genetic_map_p) {
        errorvf("For the analysis type specified, a genetic map is required.\n");
        errorvf("However, none is in the database (value specified as -2).\n");
        EXIT(INPUT_DATA_ERROR);
    }
}

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
        errorvf("For the analysis type specified, a physical map is required,\n       but none are available.\n");
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
            warnvf("No physical map was specified for the database.  So analysis as\n");
            warnvf("PLINK, PSEQ, Eigenstrat, SHAPEIT, IQLS and others will not be allowed.\n");
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
        batchfdb(Value_Base_Pair_Position_Index);
        
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
            batchfdb(Value_Base_Pair_Position_Index);
            
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
        if (option_selected < 1 || option_selected > bpps2) {
            printf("? Invalid option selected. It must be between 0 - %d\n", bpps2);
            option = 1; // as it says above ... pick the first option as the default
        } else
            option = option_selected;
        // Ask the user which option to use...
        printf("Physical map selection menu:\n");
        printf("0) Done with this menu - please proceed\n");
        for (i = 0; i < bpps; i++) {
            printf("%c", ((i+1) == option) ? '*' : ' ');        // msvc needs these two lines
            printf("%d) %s\n", i+1, EXLTop->MapNames[bppi[i]]   // separated
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
        batchfdb(Value_Base_Pair_Position_Index);
    }
    
    free(bppi);
}
