/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_premakeped_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  connect_loops get_loci_on_chromosome
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
   write_premakeped_ext.h:  save_premakeped_peds
*/


#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

#ifdef FL_STAT
#undef FL_STAT
#endif

#define DEFAULT_LOKI_iterations 1000
#define DEFAULT_LOKI_startoutput 1
#define DEFAULT_LOKI_freq 1
#define DEFAULT_LOKI_outputfile "Loki_log"
#define DEFAULT_LOKI_start_variance 1.0
#define DEFAULT_LOKI_start_mean 0.0

#define FL_STAT(x) ((main_chromocnt > 1)? &(null_str[0]) :      \
                    file_status(x, fl_stat))

typedef struct loki_params_ {
    int iterations, startoutput, frequency;
    char outputfile[FILENAME_LENGTH];
    double start_variance, start_mean;
} loki_param_type;

static void loki_run_options(loki_param_type *loki_params, int numchr)

{

    int choice=-1, dummy;
    double fdum;
    char cchoice[10], logfile[FILENAME_LENGTH], stem[5];

    loki_params->iterations =  DEFAULT_LOKI_iterations;
    loki_params->startoutput = DEFAULT_LOKI_startoutput;
    loki_params->frequency =   DEFAULT_LOKI_freq;
    strcpy(loki_params->outputfile, DEFAULT_LOKI_outputfile);
    loki_params->start_variance = DEFAULT_LOKI_start_variance;
    loki_params->start_mean=DEFAULT_LOKI_start_mean;

    if (main_chromocnt > 1) {
        strcpy(stem, "stem");
    } else {
        strcpy(stem, "");
        change_output_chr(loki_params->outputfile, numchr);
    }

    if (DEFAULT_OPTIONS) {
        mssgf("LOKI run parameters set to defaults");
        mssgf("as requested in the batch file.");
        return;
    }
    while (choice != 0) {
        printf("Loki run parameter selection menu: \n");
        printf("0) Done with this menu - please proceed.\n");
        printf(" 1) Number of iterations:        %d\n",
               loki_params->iterations);
        printf(" 2) Start output at iteration:   %d\n",
               loki_params->startoutput);
        printf(" 3) Output frequency:            %d iterations\n",
               loki_params->frequency);
        printf(" 4) Log file name %s:               %s\n", stem,
               loki_params->outputfile);
        printf(" 5) Start residual variance:     %f\n",
               loki_params->start_variance);
        printf(" 6) Start mean:                  %f\n",
               loki_params->start_mean);
        printf("Select from options 0-6 > ");
        fflush(stdout);
        IgnoreValue(fgets(cchoice, 9, stdin));
        sscanf(cchoice, "%d", &(choice)); newline;
        dummy=0; fdum=0.0;
        strcpy(logfile, "");
        switch(choice) {
        case 0:
            break;
        case 1:
            printf("Enter number of iterations > ");
            fcmap(stdin, "%d", dummy); newline;
            ((dummy > 0)? loki_params->iterations = dummy :
             printf("Invalid value entered for iteration.\n"));
            break;

        case 2:
            printf("Enter first output iteration number > ");
            fcmap(stdin, "%d", dummy); newline;
            ((dummy > 0)? loki_params->startoutput = dummy :
             printf("Invalid value entered for starting output iteration number.\n"));
            break;

        case 3:
            printf("Enter frequency of output > ");
            fcmap(stdin, "%d", dummy); newline;
            ((dummy > 0)? loki_params->frequency = dummy :
             printf("Invalid value entered for output frequency.\n"));
            break;

        case 4:
            printf("Enter log file name > ");
            fcmap(stdin, "%s", logfile); newline;
            if (strcmp(logfile, "")) {
                strcpy(loki_params->outputfile, logfile);
            } else {
                printf("Invalid value entered for log file name.\n");
            }
            break;

        case 5:
            printf("Enter starting residual variance > ");
            fcmap(stdin, "%g", &fdum ); newline;
            ((fdum >= 0.0)? loki_params->start_variance = fdum :
             printf("Invalid value entered for starting varianec.\n"));
            break;

        case 6:
            printf("Enter starting residual variance > ");
            fcmap(stdin, "%g", &(loki_params->start_mean)); newline;
            break;
        default:
            printf("Unknown option  %s\n", cchoice);
            break;
        }
    }
    draw_line();
    return;

}



static void loki_file_names(char *file_names[], int orig_ids[2],
			    analysis_type analysis, int has_orig,
			    int has_uniq)
{
    int choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12], null_str[3] = "";

//  orig_ids[0] = 1;
//  orig_ids[1] = 1;


    if (DEFAULT_OUTFILES) {
        mssgf("LOKI output file names set to defaults.");
        return;
    }

    ((num_traits > 1 && LoopOverTrait)?
     strcpy(stem, "in each subdir") : strcpy(stem, ""));

    if (main_chromocnt > 1) {
        strcpy(stem, "stem");
        // Here -9 is a flag to return only the file name portion...
        analysis->replace_chr_number(file_names, -9);
    } else {
        strcpy(stem, "");
    }

    while (choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("  Loki file name menu\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Pedigree file name %s:                  %-15s\t%s\n",
               stem, file_names[0], FL_STAT(file_names[0]));
        printf(" 2) Frequency file name %s:                 %-15s\t%s\n",
               stem, file_names[1], FL_STAT(file_names[1]));
        printf(" 3) Map file name %s:                       %-15s\t%s\n",
               stem, file_names[2], FL_STAT(file_names[2]));
        printf(" 4) Locus control file name %s:             %-15s\t%s\n",
               stem, file_names[3], FL_STAT(file_names[3]));
        printf(" 5) Link control file name %s:              %-15s\t%s\n",
               stem, file_names[4], FL_STAT(file_names[4]));
        printf(" 6) Overall control file name %s:           %-15s\t%s\n",
               stem, file_names[5], FL_STAT(file_names[5]));
        printf(" 7) Overall parameter file name %s:         %-15s\t%s\n",
               stem, file_names[6], FL_STAT(file_names[6]));
        printf(" 8) Shell script file name %s:              %-15s\t%s\n",
               stem, file_names[7], FL_STAT(file_names[7]));

        individual_id_item(9, analysis, orig_ids[0], 39, 2,
                           0, 0);
        pedigree_id_item(10, analysis, orig_ids[1], 32, 2, 0);

        printf("Select options 0-10 > ");

        choice = -1;

        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        case 1:
            printf("New pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);
            newline;
            break;
        case 2:
            printf("New frequency file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]);
            newline;
            break;
        case 3:
            printf("New map file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]);
            newline;
            break;
        case 4:
            printf("New locus control file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]);
            newline;
            break;
        case 5:
            printf("New link control file name %s > ", stem);
            fcmap(stdin, "%s", file_names[4]);
            newline;
            break;
        case 6:
            printf("New overall control file name %s > ", stem);
            fcmap(stdin, "%s", file_names[5]);
            newline;
            break;
        case 7:
            printf("New overall parameter file name %s > ", stem);
            fcmap(stdin, "%s", file_names[6]);
            newline;
            break;
        case 8:
            printf("New shell script file name %s > ", stem);
            fcmap(stdin, "%s", file_names[7]);
            newline;
            break;
        case 9:
            orig_ids[0] =
                individual_id_item(0, analysis, orig_ids[0],
                                   35, 1, has_orig, has_uniq);
            break;
        case 10:
            orig_ids[1] =
                pedigree_id_item(0, analysis, orig_ids[1],
                                 35, 1, has_orig);
            break;
        default:
            warn_unknown(cchoice);
            break;
        }
    }

    individual_id_item(0, analysis, orig_ids[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, analysis, orig_ids[1], 0, 3, has_orig);
    if (main_chromocnt > 1) {
        change_output_chr(file_names[7], global_chromo_entries[0]);
        strcat(file_names[7], ".sh");
    }
}

static void write_LOKI_map(char *map_file, int numchr,
			   linkage_locus_top *LTop)

{

    int nloop, tr, locus, num_affec=num_traits;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;

    if (genetic_distance_index < 0 || genetic_distance_sex_type_map == NONE_AVAILABLE_GDMT) {
      warnvf("Unable to write map file(s) '%s' since no usable genetic map information was available.\n", map_file);
      return;
    }

    NLOOP;
    for(tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], map_file);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open map file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (LTop->map_distance_type == 'k') {
            fprintf(fp, "MAP FUNCTION KOSAMBI\n");
        } else {
            fprintf(fp, "MAP FUNCTION HALDANE\n");

        }

        for (locus = 0; locus < NumChrLoci; locus++) {
            if (LTop->Locus[ChrLoci[locus]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[locus]].Type == BINARY) {
	        // In the current implementation of mega2 there are three
	        // different types of sex map types....
		if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
		  fprintf(fp, "Position  %-15s %6.3f, %10.7f\n",
			  LTop->Locus[ChrLoci[locus]].LocusName,
			  LTop->Marker[ChrLoci[locus]].pos_male,
			  LTop->Marker[ChrLoci[locus]].pos_female);
		} else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
		  fprintf(fp, "Position  %-15s %10.7f\n",
			  LTop->Locus[ChrLoci[locus]].LocusName,
			  LTop->Marker[ChrLoci[locus]].pos_female);
		} else if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
		  fprintf(fp, "Position  %-15s %10.7f\n",
			  LTop->Locus[ChrLoci[locus]].LocusName,
			  LTop->Marker[ChrLoci[locus]].pos_avg);
		} else {
		  // This should never happen...
		  errorvf("Internal error writing map file.\n");
		  EXIT(SYSTEM_ERROR);
		}
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
}

static void write_LOKI_freq(char *freq_file, int numchr,
			    linkage_locus_top *LTop)

{

    int nloop, tr, locus, num_affec=num_traits ;
    int all;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], freq_file);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open frequency file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        for (locus = 0; locus < NumChrLoci; locus++) {
            if (LTop->Locus[ChrLoci[locus]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[locus]].Type == BINARY) {
                fprintf(fp, "Frequency  %s",
                        LTop->Locus[ChrLoci[locus]].LocusName);
                for (all=0; all < LTop->Locus[ChrLoci[locus]].AlleleCnt; all++) {
                    fprintf(fp, " %s,%7.6f",
			    format_allele(&LTop->Locus[ChrLoci[locus]], all+1),
			    LTop->Locus[ChrLoci[locus]].Allele[all].Frequency);
                }
                fprintf(fp, "\n");
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
}


static void write_LOKI_marker(char *marker_file, int numchr,
			      linkage_locus_top *LTop,
			      char *pedigree_file)
{
    int nloop, tr, *trp, loc1, locus, num_affec=num_traits ;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;
    int mrk, lastmrk;

    NLOOP;

    if (num_traits > 0) {
        trp = &(global_trait_entries[0]);
    } else {
        trp= NULL;
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], marker_file);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open marker file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }

        mrk=0;
        for(locus=0; locus < NumChrLoci; locus++) {
            if (LTop->Locus[ChrLoci[locus]].Type == NUMBERED) {
                mrk++;
            }
        }

        if (mrk > 0) {
            fprintf(fp, "ARRAY hap1_%d(%d),hap2_%d(%d)\n",
                    numchr, mrk, numchr, mrk);
        }

        /* To write the pedigree data columns in compact notatoon for markers
           as well as allow interspersing of trait loci,
           when we see a marker locus, increment marker count,
           when we see a trait locus, write out marker loci then write out the
           trait locus.
        */

        fprintf(fp, "# Columns in pedigree data\n");
        fprintf(fp, "file \"%s\",family,id,father,mother,sx",
                pedigree_file);

        if (LoopOverTrait == 1) {
            fprintf(fp, ",%s", LTop->Pheno[*trp].TraitName);
        }

        mrk=0; lastmrk=1;

        for (loc1 = 0; loc1 < NumChrLoci; loc1++) {
            locus = ChrLoci[loc1];
            switch(LTop->Locus[locus].Type) {
            case AFFECTION:
            case QUANT:
                if (LoopOverTrait == 0) {
                    /* See if we have collected an array of markers */
                    if (mrk > lastmrk) {
                        /* write in array form */
                        fprintf(fp, ",(hap1_%d(j), hap2_%d(j), j=%d,%d)",
                                numchr, numchr, lastmrk, mrk);
                        /* update the beginning of the next array */
                    } else {
                        if (mrk == lastmrk) {
                            fprintf(fp, ",hap1_%d(%d), hap2_%d(%d)",
                                    numchr, mrk, numchr, mrk);
                        }
                    }
                    lastmrk=mrk+1;
                    fprintf(fp, ",%s", LTop->Locus[locus].LocusName);
                }
                break;

            case BINARY:
            case NUMBERED:
                mrk++;
                break;

            default:
                break;
            }
        }

        /* Dump markers if any are left */
        if (mrk > lastmrk) {
            /* write in array form */
            fprintf(fp, ",(hap1_%d(j), hap2_%d(j), j=%d,%d)",
                    numchr, numchr, lastmrk, mrk);
            /* update the beginning of the next array */
            lastmrk = mrk+1;
        } else {
            if (mrk == lastmrk) {
                fprintf(fp, ",hap1_%d(%d), hap2_%d(%d)",
                        numchr, mrk, numchr, mrk);
            }
            lastmrk=mrk+1;
        }

        mrk=0;
        fprintf(fp, "\n");
        for (locus = 0; locus < NumChrLoci; locus++) {
            if (LTop->Locus[ChrLoci[locus]].Type == NUMBERED) {
                mrk++;
                fprintf(fp,
                        "marker locus %s [hap1_%d(%d), hap2_%d(%d)]\n",
                        LTop->Locus[ChrLoci[locus]].LocusName,
                        numchr, mrk, numchr, mrk);
            }
        }
        fclose(fp);
        if (nloop == 1) break;
        trp++;
    }
}

/* chromosome-specific linkage file */

static void write_LOKI_link(char *chromo_file, int numchr,
			    linkage_locus_top *LTop)

{

    int nloop, tr, locus, num_affec=num_traits ;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], chromo_file);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open link file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(fp, "link \'chr%d\'", numchr);
        for (locus = 0; locus < NumChrLoci; locus++) {
            if (LTop->Locus[ChrLoci[locus]].Type == BINARY ||
                LTop->Locus[ChrLoci[locus]].Type == NUMBERED) {
                fprintf(fp, ",\n %s", LTop->Locus[ChrLoci[locus]].LocusName);
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }

    return;

}
/*--------------end of write_LOKI_link-----------------*/

/* can be used for global control and parameter files */

static void write_LOKI_global_file(char *cfl_name, int numchr)

{

    int nloop, tr, num_affec=num_traits;
    char cfl[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (tr==0 && nloop > 1) continue;
        change_output_chr(cfl_name, 0);
        sprintf(cfl, "%s/%s", output_paths[tr], cfl_name);
        fp = ((numchr > 0)? fopen(cfl, "a") : fopen(cfl, "w"));
        if (fp == NULL) {
            errorvf("Could not open global control file %s.\n", cfl);
            EXIT(FILE_WRITE_ERROR);
        }

        if (numchr == 0) {
            script_time_stamp(fp);
        } else {
            change_output_chr(cfl_name, numchr);
            fprintf(fp, "include \"%s\"\n", cfl_name);
        }
        fclose(fp);
        if (nloop == 1) break;

    }

}

/*-------------end of write_LOKI_global_file-----------------*/

/**
   It appears that there are a number of keywords (case does not matter) that cannot
   appear as the name of a field (e.g., phenotype). If these keywords do appear as a trait
   name the LOKI 'prep' program will generate an error indicating the keyword, and so
   the shell script created by Mega2 will not run.

   These keywords include:
*/
const char *LOKI_keywords[] = {
  "AFFECTED", "ARRAY", "CENSORED", "CONSTANT", "DISCRETE", "DO", "FILE", "FILTER", "GROUP",
  "INTEGER", "LINK", "LOG", "MARKER", "MISSING", "MODEL", "MULTIPLE", "OUTPUT", "PEDIGREE",
  "RANDOM", "REAL", "SET", "SEX", "TRAIT", "USE", "WHILE"
};
#define LOKI_KEYWORDS_SIZE                (sizeof(LOKI_keywords) / sizeof(const char *))

/**
   Is the given string a LOKI keyword?
   @param[in] str Potential LOKI keyword string
   @return 1 (Yes/true), 0 (No/false)
*/
static const int is_LOKI_keyword(const char *str) {
  for (size_t i=0; i<LOKI_KEYWORDS_SIZE; i++)
    if (strcasecmp(str, LOKI_keywords[i]) == 0) return 1;
  return 0;
}

/**
   Warn the user if the trait name is a LOKI keyword...
   @param[in] traitName Potential LOKI keyword string
*/
static void warn_if_LOKI_keyword(const char *traitName) {
  if (is_LOKI_keyword(traitName))
    warnvf("Trait name '%s' is a LOKI keyword.\n", traitName);
}

static void write_LOKI_control(char *file_names[], int numchr,
                               linkage_locus_top *LTop,
                               char *loki_log_file,
                               int global, int first_time)
{
    int nloop, *trp, tr, locus, num_affec=num_traits;
    char mfl[2*FILENAME_LENGTH];
    FILE *fp;
    
    NLOOP;
    
    if (num_traits > 0) {
        trp = &(global_trait_entries[0]);
    } else {
        trp= NULL;
    }
    
    for (tr=0; tr <= nloop; tr++) {
        if ((nloop > 1) && (tr == 0)) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], file_names[5]);
        
        fp = fopen(mfl, (global && (first_time == 0)) ? "a" : "w");
        if (fp == NULL) {
            errorvf("Could not open control file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        
        if (first_time) {
            script_time_stamp(fp);
            fprintf(fp, "# The output file\n");
            
            change_output_chr(loki_log_file, numchr);
            fprintf(fp, "log \"%s\"  \n", loki_log_file);
            fprintf(fp, "# Set the default missing value for pedigree and genotype variables.\n");
            fprintf(fp, "MISSING [\"PG\"] \"0\" \n");
            fprintf(fp, "SEX sx 1,2\n");
            
            if (LoopOverTrait == 1 && trp != NULL) {
                switch (LTop->Locus[*trp].Type) {
                    case AFFECTION:
                        // This really is not necessary because the MISSING [PG] statement covers this...
                        //
                        // NOTE: affection status (traits) themselves are written by 'write_files.cpp:write_affection_data()'
                        // which is called by 'write_premakeped.cpp:save_premakeped_peds()' which is called by
                        // 'create_LOKI_files()'. 'write_affection_data()' always uses '0' for a missing affection status.
                        fprintf(fp, "MISSING \"0\", %s\n",
                                LTop->Pheno[*trp].TraitName);
                        warn_if_LOKI_keyword(LTop->Pheno[*trp].TraitName);
                        break;
                    case QUANT:
                        // Only write the MISSING command if the trait has a missing value (e.g., .item.read true).
                        //
                        // NOTE: quantitative traits themselves are written by 'write_files.cpp:write_quantitative_data()'
                        // which is called by 'write_premakeped.cpp:save_premakeped_peds()' which is called by
                        // 'create_LOKI_files()'. The batch item 'Value_Missing_Quant_On_Output' is used by
                        //'write_quantitative_data()' when the quant == MissingQuant.
                        if (Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read) {
                            fprintf(fp, "MISSING \"%s\", %s\n",
                                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name,
                                    LTop->Pheno[*trp].TraitName);
                        }
                        warn_if_LOKI_keyword(LTop->Pheno[*trp].TraitName);
                        break;
                        
                    default:
                        break;
                }
            } else {
                for (locus = 0; locus < num_traits; locus++) {
                    if (global_trait_entries[locus] < 0) continue;
                    switch(LTop->Locus[global_trait_entries[locus]].Type) {
                        case AFFECTION:
                            // This really is not necessary because the MISSING [PG] statement covers this...
                            fprintf(fp, "MISSING \"0\", %s\n",
                                    LTop->Pheno[global_trait_entries[locus]].TraitName);
                            warn_if_LOKI_keyword(LTop->Pheno[global_trait_entries[locus]].TraitName);
                            break;
                            
                        case QUANT:
                            // Only write the MISSING command if one or more QTLs are missing (e.g., .item.read true).
                            // Even so, we can't be sure that it is this one that is missing.
                            if (Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read) {
                                fprintf(fp, "MISSING \"%s\", %s\n",
                                        Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name,
                                        LTop->Pheno[global_trait_entries[locus]].TraitName);
                            }
                            warn_if_LOKI_keyword(LTop->Pheno[global_trait_entries[locus]].TraitName);
                            break;
                            
                        default:
                            break;
                    }
                }
            }
            
            fprintf(fp, "PEDIGREE family,id,father,mother\n");
            fprintf(fp, "TRAIT LOCUS QTL\n");
            fprintf(fp, "# Insert your own model here. \n");
            if (LoopOverTrait == 1) {
                fprintf(fp, "model %s=QTL\n", LTop->Pheno[*trp].TraitName);
            } else {
                if (num_traits > 1) {
                    // THIS LOOKS LIKE A BUG TO ME....
                    // Why are we assuming that global_trait_entries[1] >= 0 ????
                    if (global_trait_entries[0] < 0) {
                        /* write only the first QTL */
                        fprintf(fp, "model %s=QTL\n",
                                LTop->Pheno[global_trait_entries[1]].TraitName);
                    } else {
                        fprintf(fp, "model %s=QTL\n",
                                LTop->Pheno[global_trait_entries[0]].TraitName);
                    }
                }
            }
        }
        /* file_names[3] = Loki_locus, file_names[4] = Loki_link */
        if (!global || !first_time) {
            // Define the fields in the .ped file...
            fprintf(fp, "include \"%s\"\n", file_names[3]);
            fprintf(fp, "include \"%s\"\n", file_names[4]);
        }
        fclose(fp);
        if (nloop == 1) break;
        trp++;
    }
}

static void write_LOKI_param(char *file_names[],
			     int numchr,
			     linkage_locus_top *LTop,
			     loki_param_type loki_params)
{

    int nloop, tr, locus, num_affec=num_traits ;
    char mfl[2*FILENAME_LENGTH];
    double averaged_max=0.0, male_max=0.0, female_max=0.0;
    FILE *fp;

    for(locus=0; locus < NumChrLoci; locus++) {
        if (LTop->Locus[ChrLoci[locus]].Type == NUMBERED ||
           LTop->Locus[ChrLoci[locus]].Type == BINARY) {
            if (LTop->Marker[ChrLoci[locus]].pos_avg > averaged_max) {
                averaged_max = LTop->Marker[ChrLoci[locus]].pos_avg;
            }
            if (LTop->Marker[ChrLoci[locus]].pos_male > male_max) {
                male_max = LTop->Marker[ChrLoci[locus]].pos_male;
            }
            if (LTop->Marker[ChrLoci[locus]].pos_female > female_max) {
                female_max = LTop->Marker[ChrLoci[locus]].pos_female;
            }
        }
    }

    NLOOP;
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], file_names[6]);
        if ((fp = fopen(mfl, "w")) == NULL) {
            errorvf("Could not open link file %s.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        script_time_stamp(fp);
        fprintf(fp, "Iterations %d\n", loki_params.iterations);
        fprintf(fp, "Start Output %d, %d\n",
                loki_params.startoutput, loki_params.startoutput);
        fprintf(fp, "Output frequency %d, %d\n",
                loki_params.frequency, loki_params.frequency);
        fprintf(fp, "Start Residual Variance 1.0\n");
        fprintf(fp, "Start Mean 0.0\n");
        // In the current implementation of mega2 there are three
        // different types of sex map types....
        if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
            fprintf(fp, "Male Map \'Chromosome %d\' 0.0, %10.7f\n", numchr, male_max);
            fprintf(fp, "Female Map \'Chromosome %d\' 0.0, %10.7f\n", numchr, female_max);
            fprintf(fp, "Total Male Map 3600.0\n");
            fprintf(fp, "Total Female Map 3600.0\n");
        } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
            fprintf(fp, "Female Map \'Chromosome %d\' 0.0, %10.7f\n", numchr, female_max);
            fprintf(fp, "Total Female Map 3600.0\n");
        } else if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            fprintf(fp, "Map \'Chromosome %d\' 0.0, %10.7f\n", numchr, averaged_max);
            fprintf(fp, "Total Map 3600.0\n");
        } else {
            // This should never happen...
            errorvf("Internal error writing '%s' file.\n", file_names[6]);
            EXIT(SYSTEM_ERROR);
        }
        
        /* Not writing the total map */
        fprintf(fp, "include \"%s\"\n", file_names[1]);
        fprintf(fp, "include \"%s\"\n", file_names[2]);
        fclose(fp);
        if (nloop == 1) break;
    }
}

static void write_LOKI_script(char *output_files[],
			      int numchr,
			      int time_stamp)
{

    /* control file = file_names[5]
       param_file = file_names[6]
    */

    int nloop, tr, num_affec=num_traits;
    FILE *filep;
    char mfl[2*FILENAME_LENGTH], *fl;
    char *syscmd;
    char append[2];

    NLOOP;

    if (time_stamp) {
        strcpy(append, "w");
    } else {
        strcpy(append, "a");
    }

    fl = &(output_files[7][0]);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(mfl, "%s/%s", output_paths[tr], fl);
        if ((filep = fopen(mfl, append)) == NULL) {
            errorvf("Unable to open %s for writing.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
	fprintf(filep, "#!/bin/csh -f\n");
	/* print some identification */
	fprintf(filep, "# C-shell file name: %s\n", fl);
        if (time_stamp) {
            script_time_stamp(filep);
        }

        fprintf(filep, "# Chromosome number:    %d\n", numchr);
#ifdef RUNSHELL_SETUP
        fprintf_env_checkset_csh(filep, "_PREP", "prep");
        fprintf(filep, "$_PREP %s\n", output_files[5]);
	fprintf_status_check_csh(filep, "prep", 1);

        fprintf_env_checkset_csh(filep, "_LOKI", "loki");
        fprintf(filep, "$_LOKI %s\n", output_files[6]);
	fprintf_status_check_csh(filep, "loki", 1);
#else /* RUNSHELL_SETUP */
        fprintf(filep, "prep %s\n", output_files[5]);
        fprintf(filep, "loki %s\n", output_files[6]);
#endif /* RUNSHELL_SETUP */
        fclose(filep);
        if (time_stamp) {
            syscmd = CALLOC((strlen(mfl) + 10), char);
            sprintf(syscmd, "chmod +x %s", mfl);
            System(syscmd);
        }
        if (nloop == 1) break;
    }
}


void create_LOKI_files(linkage_ped_top **LPedTop, int *numchr,
		       file_format *infl_type, char *file_names[],
		       int untyped_ped_opt, analysis_type analysis)
{

    int i;
    int xlinked;
    linkage_locus_top *LTop = Top->LocusTop;

    loki_param_type loki_params;

    linkage_ped_tree *LPed;

    loki_file_names(file_names, OrigIds, analysis,
                    Top->OrigIds, Top->UniqueIds);
    loki_run_options(&(loki_params), *numchr);

    if (*infl_type == LINKAGE) {
        for (i = 0; i < Top->PedCnt; i++)  {
            LPed = &(Top->Ped[i]);
            if (LPed->Loops != NULL)   {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)    {
                    errorf("Fatal error! Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    for (i=0; i < main_chromocnt; i++) {
        *numchr=global_chromo_entries[i];
        xlinked=(((LTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
                  (LTop->SexLinked == 1))? 1 : 0);
        if (xlinked) {
            warnf("Files not created for X-Chromosome.");
            continue;
        }

        get_loci_on_chromosome(*numchr);
        if (main_chromocnt > 1) {
            if (i==0) {
                change_output_chr(file_names[5], 0);
                write_LOKI_control(file_names, *numchr, Top->LocusTop,
                                   &(loki_params.outputfile[0]), 1, 1);
                write_LOKI_global_file(file_names[6], 0);
            }
        }
        analysis->replace_chr_number(file_names, *numchr);
        omit_peds(untyped_ped_opt, Top);
        if (i==0)   {
            create_mssg(analysis);
        }

        save_premakeped_peds(file_names[0], Top, TO_LOKI, 0); // write_premakeped.cpp
        write_LOKI_freq(file_names[1], *numchr, Top->LocusTop);
        write_LOKI_map(file_names[2], *numchr, Top->LocusTop);
        write_LOKI_marker(file_names[3], *numchr, Top->LocusTop, file_names[0]);
        write_LOKI_link(file_names[4], *numchr, Top->LocusTop);
        write_LOKI_control(file_names, *numchr, Top->LocusTop,
                           &(loki_params.outputfile[0]), 0, 1);
        write_LOKI_param(file_names, *numchr, Top->LocusTop, loki_params);
        write_LOKI_script(file_names, *numchr, 1);
        sprintf(err_msg, "     Pedigree file:               %s",
                file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg, "     Frequency file:              %s",
                file_names[1]);
        mssgf(err_msg);
        sprintf(err_msg, "     Map file:                    %s",
                file_names[2]);
        mssgf(err_msg);
        sprintf(err_msg, "     Marker control file:         %s",
                file_names[3]);
        mssgf(err_msg);
        sprintf(err_msg, "     Link control file:           %s",
                file_names[4]);
        mssgf(err_msg);
        sprintf(err_msg, "     Overall Control file:        %s",
                file_names[5]);
        mssgf(err_msg);
        sprintf(err_msg, "     Overall Parameter file:      %s",
                file_names[6]);
        mssgf(err_msg);

        sprintf(err_msg, "     C-shell script:              %s",
                file_names[7]);
        mssgf(err_msg);

        if (main_chromocnt > 1) {
            change_output_chr(file_names[5], 0);
            write_LOKI_control(file_names, *numchr, Top->LocusTop,
                               &(loki_params.outputfile[0]), 1, 0);
            write_LOKI_global_file(file_names[6], *numchr);
        }
        if ((i == (main_chromocnt - 1)) && (main_chromocnt > 1)) {
            change_output_chr(file_names[5], 0);
            sprintf(err_msg, "     Combined Control file for all chromosomes:   %s",
                    file_names[5]);
            mssgf(err_msg);
            change_output_chr(file_names[6], 0);
            sprintf(err_msg, "     Combined Parameter file for all chromosomes: %s",
                    file_names[6]);
            mssgf(err_msg);
        }

    }
    draw_line();
}
