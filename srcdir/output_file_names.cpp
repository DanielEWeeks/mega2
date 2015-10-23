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

/* These routines set the default names for the output files, and
   handle chromosome number changing in the file names

   Instructions re adding a new output option:

   Routines inside this file set the default output names of options, and
   also decide which "output_file_name" items are initialized/modified
   via option-specific output file names menus.

   1) Add KEYWORD to the case statement in default_outfile_names(), and
   statements to initialize output file names (starting from 0
   typically).

   2) When multiple chromosomes are selected, for each chromosome, the
   suffixes of output files can be changed simultaneously with a single
   call to replace_chr_number(). Add KEYWORD to the case statement in
   replace_chr_number, and add calls to change_output_chr(), one for each
   file whose suffix is to be changed.

   Sometimes, output files are a mixture of chromosome-specific and
   combined files (e.g., the shell script may be a combined file for all
   chromosomes, or there might be a single overall phenotype file). This
   has to be taken into consideration, when you add calls to
   change_output_chr().

   3) Routine set_output_paths(): this sets trait-specific sub-directory
   names if trait-specific output option was chosen via the
   trait-selection menu. Usually, there is nothing to be done inside this
   function, unless the option is "special" e.g. trait-specific analysis
   is meaningless for PREST, HWE etc. The first block of the case lists
   these options, where only a single set of output files is created
   inside the output folder. If your option does not need/handle traits,
   add it to this list.

   For all other options, trait-specific sub directory names are
   initialized. SAGE 3.0 is special, since some of its files are
   trait-specific, irrespective of the user's selection.

*/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "output_file_names_ext.h"
#include "utils_ext.h"
/*
        batch_input_ext.h:  batchf
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  makedir
  output_file_names_ext.h:  CHR_STR
              utils_ext.h:  draw_line
*/


/*--------------exports---------------------*/
void          change_output_chr(char *OP_files, int numchr);

void          set_output_paths(analysis_type analysis, linkage_ped_top *Top,
			       char *output_path);
void          create_mssg(analysis_type analysis);

void          print_outfile_mssg(void);
/*--------------end of prototypes-------------*/

void default_outfile_names(const analysis_type  analysis,
                           int *numchr, char *file_names[],
                           const char *logdir)
{
    char            num[3];
#ifdef LOGDIR
    char            tmp[FILENAME_LENGTH];
    int             i;
#endif

    /* This is approximately the contents of the file_names[] array
       file_names[0] = outfl_name (pedigree file)
       file_names[1] = loutfl_name (locus file);
       file_names[2] = batfl_name (simwalk and mendel batch file);
       file_names[3] = cshfl_name (cshell file);
       file_names[4] = aspin_name (aspex tcl file);
       file_names[5] = aspdat_name (aspex pedigree file);
       file_names[6] = poutfl_name (mendel penetrance file);
       file_names[7] = sibfl_name;
       file_names[8] = goutfl_name (genotype file, SOLAR?);
       file_names[9] = moutfl_name;
       file_names[10] = cbatfl_name;
       file_names[11] = cntcshfl_name;
       file_names[12] = sumfl_name (summary file);
       file_names[15] = error_genos (genotypes table);
       file_names[16] = error_sum (summary table);

       HWE files are not defined here

       NOTE: The naming scheme assumed here is such that
       each name consists of 3 components separated by a period:
       comp 1: file_name
       comp 2: extension
       comp 3: rest
    */

    // If we are dealing with human chromosomes, use the symbolic names...
    if (*numchr >= 0) CHR_STR(*numchr, num);
 
    // This is not really the right place to call this.
    // 'num' get's over written later...
    analysis->file_names(file_names, num);

    if (ErrorSimOpt == 1) {
        sprintf(file_names[15], "error_genos.%s", num);
        sprintf(file_names[16], "error_sum.%s", num);
    }
#ifdef LOGDIR
    if (strcmp(Mega2OutputPath, ".") != 0)  return;
    if (strcmp(logdir, ".") == 0)  return;
    for (i = 0; i < 17; i++) {
        strcpy(tmp, file_names[i]);
        sprintf(file_names[i], "%s/%s", logdir, tmp);
    }
#endif
}


//
// Insert the chromosome number into OP_files (in place).
//
// OP_files can be of the following forms:
// 1) <file_name> . <extension> . <rest>
// 2) <file_name> . <extension>
// Typically the <extension> is the character string '%s'.
// 
// numchr == -9 returns only the <file_name> removing <extension> and <rest> if they exist.
// numchr == 0 replaces <extension> with 'all', keeping <extension> and <rest> if they exist.
// numchr > 0 replaces <extension> with numchr, keeping <extension> and <rest> if they exist.
void change_output_chr(char *OP_files, const int numchr)
{
    char chr_str[4]={""};
    int i, pos[2], numdots=0, OP_files_len = (int)strlen(OP_files);
    char *file_name, *rest;

    pos[0]=pos[1]= -1;

    if (numchr == 0) {
        strcpy(chr_str, "all");
    } else if (numchr > 0) {
        CHR_STR(numchr, chr_str);
    }

    for(i = OP_files_len-1 ; i >= 0; i--) {
        if (numdots == 2) break;
        if (OP_files[i] == '.'){
            numdots++;
            switch(numdots) {
            case 1:
                // Pretend that an ending ".gz" isn't really an extension, but a part of the overall extension...
                if (i == OP_files_len-3 && strcasecmp(".gz", &OP_files[i]) == 0) { numdots--; continue; }
                pos[0]=i;
                break;
            case 2:
                /* shift the first dot */
                pos[1] = pos[0];
                pos[0] = i;
                break;
            default:
                break;
            }
        }
    }

    file_name = CALLOC(strlen(OP_files)+1, char);
    rest = CALLOC(strlen(OP_files)+1, char);

    if (numdots == 0) {
        /* only one component */
        strcpy(file_name, OP_files);
    } else if (numdots >= 1) {
        /* two components */
        for(i=0; i < pos[0]; i++) {
            file_name[i]=OP_files[i];
        }
        file_name[pos[0]]='\0';
        if (numdots == 2) {
            int r=0;
            for(i=pos[1]+1; i < (int) strlen(OP_files); i++) {
                rest[r]=OP_files[i]; r++;
            }
            rest[r]='\0';
        } else {
            strcpy(rest,"");
        }
    }

    if (numchr == -9) {
        sprintf(OP_files, "%s", file_name);
    } else {
        if (numdots < 2) {
            sprintf(OP_files, "%s.%s", file_name, chr_str);
        } else {
            sprintf(OP_files, "%s.%s.%s", file_name, chr_str, rest);
        }
    }

    free(file_name);
    free(rest);
}

/*---------------------------------------------------------*/

void           set_output_paths(analysis_type analysis,
				linkage_ped_top *Top)


{
    int i, j, choice = -1;
    char cchoice[10];
    char dirname[2*FILENAME_LENGTH];
    int ntraits=0, *trait_loc_num;
    int *trait_skip;

    /* LOCATION, NPL, IBD, HAPLO
       ASPEX, GH, GHP, APM-MULT
       SPLINK, LOD2,  SIMULATE, SAGE,
       TDTMAX (affection), VITESSE */

/* Instructions re adding a new output option:

   3) Routine set_output_paths(): this sets trait-specific sub-directory
   names if trait-specific output option was chosen via the
   trait-selection menu. Usually, there is nothing to be done inside this
   function, unless the option is "special" e.g. trait-specific analysis
   is meaningless for PREST, HWE etc. The first block of the case lists
   these options, where only a single set of output files is created
   inside the output folder. If your option does not need/handle traits,
   add it to this list.

   For all other options, trait-specific sub directory names are
   initialized. SAGE 3.0 is special, since some of its files are
   trait-specific, irrespective of the user's selection.

*/

    if (LoopOverTrait == 0) {
        ntraits = num_traits - 1;
    } else {
        ntraits = num_traits;
    }

    trait_loc_num = &(global_trait_entries[0]);

    if ((ntraits <= 1) ||
        ((LoopOverTrait == 0) &&
         (analysis != TO_SAGE || (analysis == TO_SAGE && !HasAff))) ||
        analysis->forbid_trait_directories() ) {
        output_paths = (char **) CALLOC_PTR((size_t) 2, char *);
        output_paths[0] = CALLOC((size_t) FILENAME_LENGTH, char);
        trait_paths = (char **) CALLOC_PTR((size_t) 2, char *);
        trait_paths[0] = CALLOC((size_t) FILENAME_LENGTH , char);
        strcpy(trait_paths[0],".");
        strcpy(output_paths[0], Mega2OutputPath);

        /* specified output directory */
        output_paths[1] = NULL;
    }

    if (analysis->forbid_trait_directories()) {
        LoopOverTrait=0;
        return;
    }
//      if (analysis == TO_SAGE)
//          if (LoopOverTrait == 0 && !HasAff) {
//              break;
//          }
        
    if (!((analysis == TO_SAGE) && (LoopOverTrait == 0 && !HasAff)) &&
        (ntraits > 1  && (LoopOverTrait == 1 || (analysis == TO_SAGE && HasAff)))) {
        /* re-arrange the output_paths data_structure */
        output_paths = CALLOC((size_t) ntraits+1, char *);
        trait_paths  = CALLOC((size_t) ntraits+1, char *);
        trait_skip   = CALLOC((size_t) ntraits+1, int);

        /* first set the default names */
        output_paths[0] = CALLOC((size_t) FILENAME_LENGTH, char);
        trait_paths[0] = CALLOC((size_t) FILENAME_LENGTH, char);
        sprintf(output_paths[0], "%s", Mega2OutputPath);
        strcpy(trait_paths[0], "");
        if (InputMode != INTERACTIVE_INPUTMODE && Mega2BatchItems[/* 16 */ Trait_Subdirs].items_read) {
            j=1;
            for (i=1; i <= num_traits; i++) {
                SKIP_TRII(trait_loc_num[i-1]);
                if (analysis->skip_trait(Top->LocusTop, trait_loc_num[i-1])) continue;
                if (strcmp(Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[0],
                           "use trait names")) {
                    /* I do not like this alaising for only this case of Trait_Subdirs. */
                    output_paths[j] = Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[j-1];
                } else {
                    output_paths[j] =
                        CALLOC(strlen(Top->LocusTop->Pheno[trait_loc_num[i-1]].TraitName)+1,
                               char);
                    strcpy(output_paths[j],
                           Top->LocusTop->Pheno[trait_loc_num[i-1]].TraitName);
                }
                j++;
            }
        } else {
            j=1;
            for (i=1; i <= ntraits; i++) {
                SKIP_TRII(trait_loc_num[i-1]);
                output_paths[j] = CALLOC((size_t) 2*FILENAME_LENGTH, char);
                if (analysis->skip_trait(Top->LocusTop, trait_loc_num[i-1])) continue;
                sprintf(output_paths[j], "%s",
                        Top->LocusTop->Pheno[trait_loc_num[i-1]].TraitName);
                j++;
            }
            /* first print the n path_names */
            printf("Output directory selection menu\n");
            printf("Directories will be created within %s.\n",
                   (!strcmp(Mega2OutputPath, ".") ? "Current directory" : Mega2OutputPath));
            while (choice != 0) {
                draw_line();
                printf("0) Done with this menu - please proceed\n");
                j=1;
                for (i=1; i <= ntraits; i++) {
                    SKIP_TRII(trait_loc_num[i-1]);
                    if (analysis->skip_trait(Top->LocusTop, trait_loc_num[i-1])) continue;
                    if (strcmp(Mega2OutputPath, ".")) {
                        sprintf(dirname, "%s/%s", Mega2OutputPath, output_paths[j]);
                    } else {
                        strcpy(dirname, output_paths[j]);
                    }

                    printf(" %d) Output directory for %-15s:    %s[%s]\n",
                           j,
                           Top->LocusTop->Pheno[trait_loc_num[i-1]].TraitName,
                           output_paths[j],
                           ((access(dirname, F_OK)? "new":"existing")));
                    trait_skip[j] = i-1;
                    j++;
                }
                printf("Select from 1-%d to enter a new directory name > ",
                       j-1);
                fcmap(stdin, "%s", cchoice); newline;
                sscanf(cchoice, "%d", &choice);
                if (choice == 0) break;

                if (choice > 0 && choice <= j-1) {
                    printf("Enter new directory name for trait %s > ",
                           Top->LocusTop->Pheno[trait_loc_num[trait_skip[choice]]].TraitName);
                    fcmap(stdin, "%s", output_paths[choice]);
                }
                else
                    warn_unknown(cchoice);
                draw_line();
            }
            if (InputMode == INTERACTIVE_INPUTMODE) {
                /* write the trait sub-directory names  into batch file */
                Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names = CALLOC((size_t) num_traits+1, char *);
                j = 1;
                for(i=1; i <= num_traits; i++) {
                    if (analysis->skip_trait(Top->LocusTop, trait_loc_num[i-1])) continue;
                    Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[j-1] = output_paths[j];
                    j++;
                }
                Mega2BatchItems[/* 16 */ Trait_Subdirs].value.mult_names[j] = NULL;
                batchf(/* 16 */ Trait_Subdirs);
            }
        }
    }

    if ((ntraits > 1) &&
        ((LoopOverTrait == 1) || (analysis == TO_SAGE && HasAff))) {
        j = 1;
        for (i=1; i<=ntraits; i++) {
            SKIP_TRII(trait_loc_num[i-1]);
            if (analysis->skip_trait(Top->LocusTop, trait_loc_num[i-1])) continue;
            if (strcmp(Mega2OutputPath, ".")) {
                trait_paths[j] = CALLOC((size_t) FILENAME_LENGTH, char);
                strcpy(trait_paths[j], output_paths[j]);
                sprintf(dirname, "%s/%s", Mega2OutputPath, output_paths[j]);
            } else {
                /* don't allocate these are the same */
                trait_paths[j] = &(output_paths[j][0]);
                strcpy(dirname, output_paths[j]);
            }
            if ((access(dirname, F_OK)) != 0) {
                makedir(dirname);
                sprintf(err_msg, "Created directory %s", output_paths[j]);
                mssgf(err_msg);
            }
            strcpy(output_paths[j], dirname);
            j++;
        }
    }
    return;
}

void Free_output_paths(analysis_type analysis) {
    int i;
    int ntraits;
    if (LoopOverTrait == 0) {
        ntraits = num_traits - 1;
    } else {
        ntraits = num_traits;
    }
    
    if (ntraits > 1  && (LoopOverTrait == 1 || (analysis == TO_SAGE && HasAff))) {
        for (i=1; i <= ntraits; i++) {
            /* output_path[it] may be aliased to 
               Mega2BatchItems[/ * 16 * / Trait_Subdirs].value.mult_names[it+1]
               Check for it!! in Free_batch_items()
             */
            free(output_paths[i]);
        }
    }
    free(output_paths[0]);
    free(output_paths);

    if (ntraits > 1  && (LoopOverTrait == 1 || (analysis == TO_SAGE && HasAff))) {
        for (i=1; i <= ntraits; i++)
            if (trait_paths[i] != output_paths[i])
                free(trait_paths[i]);
    }
    free(trait_paths[0]);
    free(trait_paths);
}

char *file_status(char *fl_name, char *status)

{
    /* set the status to "[overwrite]", "[new]", or "[each trait]"
       depending on whether file is being written for one trait or many */
    char file_name[2*FILENAME_LENGTH+1];

    sprintf(file_name, "%s/%s", output_paths[0], fl_name);
    if (LoopOverTrait==0 || num_traits <= 1) {
        /* one file */
        if (access(file_name, F_OK) == 0)
            sprintf(status, "[overwrite]");
        else
            sprintf(status, "[new]");
    }
    else
        sprintf(status, "[each trait]");

    return &(status[0]);
}

void create_mssg(analysis_type analysis)
{
    /*  prog_name(analysis, ProgName); */
    if (LoopOverTrait == 1 && num_traits > 1) {
        if (analysis->loops())
            mssgvf("Mega2 created the following file(s) for %s\n", ProgName);
        else
            mssgvf("For each trait, Mega2 created the following file(s) for %s\n", ProgName);
    } else {
        mssgvf("Mega2 created the following file(s) for %s:\n", ProgName);
    }
}

void shorten_path(char *old_path, char *new_path)

{
    char *path_end, *p;

    path_end=strrchr(old_path, '/');
    if (path_end == NULL) {
        strcpy(new_path, old_path);
    } else {
        if ((path_end+1) == NULL) {
            strcpy(new_path, old_path);
        } else {
            p = (path_end + 1);
            strcpy(new_path, p);
        }
    }
    return;
}

void          print_outfile_mssg(void)

{
    if (strcmp(Mega2OutputPath, ".")) {
        printf("Output files will be written in directory %s.\n", Mega2OutputPath);
    } else {
        printf("Output files will be written in current directory.\n");
    }
}

//
// The chromosome number is used unless the variables 'human_*' (globals) are set, in which
// case the symbolic value for a human chromosome is used. The 'human_*' variables are only set in
// annotated_ped_file.c:read_common_map_file() when a chromosome is first encountered.
// I'm assuming that this has something to do with extending mega2 to be used in some other species.
void CHR_STR(int numchr, char chr_str[])
{
    if (numchr == UNKNOWN_CHROMO) {
        if (human_unknown) {
            strcpy(chr_str, "U");
        } else {
            sprintf(chr_str, "%d", numchr);
        }

    } else if (numchr ==  SEX_CHROMOSOME) {
        if (human_x) {
            strcpy(chr_str, "X");
        } else {
            sprintf(chr_str, "%d", numchr);
        }

    } else if (numchr ==  PSEUDO_X) {
        if (human_xy) {
            strcpy(chr_str, "XY");
        } else {
            sprintf(chr_str, "%d", numchr);
        }

    } else if (numchr ==  MALE_CHROMOSOME) {
        if (human_y) {
            strcpy(chr_str, "Y");
        } else {
            sprintf(chr_str, "%d", numchr);
        }

    } else if (numchr ==  MITO_CHROMOSOME) {
        if (human_mt) {
            strcpy(chr_str, "MT");
        } else {
            sprintf(chr_str, "%d", numchr);
        }
    } else {

        // the autosomes 'human_auto' have no symbolic value in humans...
        if (numchr <= 9 ) {
           sprintf(chr_str, "0%1d", numchr);
        } else if (numchr <= 99 ) {
	   sprintf(chr_str, "%2d", numchr);
        } else if (numchr == MISSING_CHROMO ) {
           errorvf("INTERNAL: Processing missing chromosome.\n");
           EXIT(OUTOF_BOUNDS_ERROR);
        } else {
           errorvf("chromosome number (%d) greater than 99.\n", numchr);
           EXIT(OUTOF_BOUNDS_ERROR);
        } 

    }
}
