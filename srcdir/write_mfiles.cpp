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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "common.h"
#include "typedefs.h"
#include "R_output.h"

#include "write_mfiles.h"

#include "batch_input_ext.h"
#include "R_output_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_mendel7_files_ext.h"
#include "write_files_ext.h"
#define WRITE_MFILES_EXT_H
#include "write_mfiles_ext.h"

/*
           R_output_ext.h:  append_R_commands sw2_R_setup
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file haldane_theta init_file kosambi_theta move_file
            linkage_ext.h:  connect_loops get_loci_on_chromosome switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line log_line script_time_stamp strtail
write_mendel7_files_ext.h:  csv_mendel7_pen_file csv_save_mendel_peds csv_write_mendel_locus_file csv_write_mendel_map
*/


/*---------------static functions-----------------------*/
static void mendel_file_names(int glob_files,
                              char *file_names[], int *orig_id_opt,
                              int has_aff, int has_orig,
                              int has_uniq,
                              analysis_type *analysis, int *write_pen_file);
static void simwalk2_file_names(char *file_names[], int *orig_id_opt,
				int *r_setup, int has_orig, int has_uniq);
/* static int  put_quant_locus_last(); */
static void write_cshell(int numchr, analysis_type analysis,
			 char *file_names[], int global, int xlinked,
			 int qtdt);
static void write_batch(char *loutfl_name, char *outfl_name,
			char *batfl_name, linkage_ped_top *LLTop,
			int write_trait);
static void write_batch2(char *batfl_name, int numchr,
			 linkage_ped_top *LLTop,
			 analysis_type analysis);
static void write_batch3(char *loutfl_name, char *outfl_name,
			 char *batfl_name, linkage_ped_top *LLTop);
static void mwrite_quantitative_data(FILE *filep, int locusnm, linkage_locus_rec *locus,
				     linkage_ped_rec *entry);
static void lwrite_affection_data(FILE *filep, int locusnm, linkage_locus_rec *locus,
				  linkage_ped_rec *entry);
static void swrite_affection_data(FILE *filep, int locusnm,
				  linkage_locus_rec *locus,
				  linkage_ped_rec *entry);
static void mwrite_affection_data(FILE *filep, int locusnm,
				  linkage_locus_rec *locus,
				  linkage_ped_rec *entry);
static void mwrite_binary_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry);
static void mwrite_numbered_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry);
static void save_mendel_peds(char *file_name, linkage_ped_top *Top,
			     analysis_type analysis);

static int  write_mendel_pen_file(char *poutfl_name,
                                  linkage_locus_top * LTop,
                                  analysis_type analysis);
static void  init_mendel_files(char *file_names[],
			       int glob_files,
			       analysis_type analysis,
			       linkage_locus_top *Ltop,
			       int write_pen_file);

/*-----------------exports----------------*/

void  write_mendel_locus_file(char *file_name,
			      linkage_locus_top *LTop,
			      analysis_type analysis, int sex_linked);

int  create_mendel_file(linkage_ped_top **LPedTreeTop, char *mapfl_name,
			ped_top *PedTreeTop, file_format *infl_type,
			file_format *outfl_type,
			int *numchr, analysis_type *analysis,
			char *file_names[], int untyped_ped_opt);
void create_mendel_formats(int famwid, int perwid, char *fformat,
			   char *pformat);
void write_mendel_map(char *mapfl, linkage_locus_top *LTop,
		      analysis_type analysis);

int get_aff_status(linkage_locus_rec *Locus, analysis_type analysis,
		   linkage_ped_rec *Entry, int locus1);

/*-------------end of prototypes----------*/

/* Change the ChrLoci structure */
static int  put_quant_locus_last(linkage_ped_top *Top)

{
    struct ql {
        int qlnum;
        struct ql *next;
    };

    struct ql *quantloc=NULL, *qfirst=NULL, *qnext=NULL;
    int loc, loc1, *order, nonq=0, qlnum=0;

    /* find out the number of quant locus and their positions,
       then move these to the end using linkage_set_locus_order, only
       if LoopOverTrait = 0 */

    if (LoopOverTrait == 1) {
        for (loc = 0; loc < NumChrLoci; loc++) {
            if (Top->LocusTop->Locus[ChrLoci[loc]].Type == QUANT) {
                qlnum++;
            }
        }
        return qlnum;
    }
    order = CALLOC((size_t) NumChrLoci, int);
    for (loc1 = 0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        if (Top->LocusTop->Locus[loc].Type == QUANT) {
            if (qlnum == 0) {
                quantloc = CALLOC((size_t) 1, struct ql);
                qfirst = quantloc;
            } else {
                quantloc->next = CALLOC((size_t) 1, struct ql);
                quantloc = quantloc->next;
            }
            quantloc->qlnum = loc;
            qlnum++;
        } else {
            order[nonq] = loc;
            nonq++;
        }
    }

    if (qlnum > 0) {
        for (loc=0, quantloc=qfirst; loc < qlnum;
             loc++, quantloc=quantloc->next) {
            order[nonq+loc] = quantloc->qlnum;
        }
    }
    /* increment all indexes by 1 to feed into linkage_set_locus_order */
    for (loc = 0; loc < NumChrLoci; loc++) {
        ChrLoci[loc] = order[loc];
    }
    if (qlnum > 0) {
        for (loc=0, quantloc=qfirst; loc < qlnum; loc++) {
            qnext = quantloc->next;
            free(quantloc);
            quantloc = qnext;
        }
    }
    free(order);
    return qlnum;

}

/*---------------------end of put_quant_locus_last--------------*/

static void write_cshell(int numchr, analysis_type analysis,
			 char *file_names[], int global,
			 int xlinked, int qtdt)

{
    char           *batfl_name, *outfl_name, *moutfl_name;
    char           *loutfl_name, *cshfl_name, *poutfl_name;
    FILE           *filep;
    char            *syscmd, chr_str[3];
    char            fl_name[2*FILENAME_LENGTH];
    char            global_shell[FILENAME_LENGTH];
    int             nloop, tr, i, num_affec=num_traits;
    char            merlin_prog[7];
#ifdef RUNSHELL_SETUP
    char            merlin_prog_env[10];
    strcpy(merlin_prog_env, (xlinked ? "_MINX" : "_MERLIN"));
#endif /* RUNSHELL_SETUP */

    // MINX (MERLIN in X) is an X-specific version of Merlin.
    strcpy(merlin_prog, (xlinked ? "minx" : "merlin"));

    batfl_name=file_names[2]; outfl_name=file_names[0];
    loutfl_name=file_names[1]; cshfl_name=file_names[3];
    poutfl_name=file_names[6];  moutfl_name=file_names[4];
    if (strcmp(output_paths[0], ".")) {
        sprintf(global_shell, "%s/%s", output_paths[0], cshfl_name);
    } else {
        strcpy(global_shell, cshfl_name);
    }
    change_output_chr(global_shell, 0);

    NLOOP;
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        sprintf(fl_name, "%s/%s", output_paths[tr], cshfl_name);
        delete_file(fl_name);
        if ((filep = fopen(fl_name, "a")) == NULL) {
            errorvf("Could not open %s/%s for writing.\n", output_paths[tr], cshfl_name);
            EXIT(FILE_WRITE_ERROR);
	}
        /* added by Nandita- */
        fprintf(filep, "#!/bin/csh -f\n");
        /* print some identification */
        fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                NOPATH
#else /* HIDEPATH */
                cshfl_name
#endif /* HIDEPATH */
            );
        script_time_stamp(filep);
        fprintf(filep, "# Chromosome number:    %d\n", numchr);
#ifdef RUNSHELL_SETUP
        // This handles the environment variable setup to allow the checking
        // functions in 'batch_run' to work correctly...
        fprintf_env_checkset_csh(filep, merlin_prog_env, merlin_prog);
        fprintf_env_checkset_csh(filep, "_MERLIN2SW2", "merlin2sw2.pl");
        fprintf_env_checkset_csh(filep, "_SIMWALK2", "simwalk2");
#endif /* RUNSHELL_SETUP */
        fprintf(filep, "echo \n");
        if (analysis == TO_MERLIN) {
            CHR_STR(numchr, chr_str);
            analysis->replace_chr_number(file_names, numchr);
            fprintf(filep, "echo Running Merlin to generate exact npl scores.\n");
            
            // http://www.sph.umich.edu/csg/abecasis/Merlin/reference.html
            // -d datafile
            //    Selects input data file, in linkage or QTDT format.
            // -p pedfile
            //    Selects pedigree file, with genotype, phenotype and family structure information.
            // -m mapfile
            //    File indicating chromosome and centimorgan position for each marker.
            //    USE WITH QTDT FORMAT INPUT FILES. Recombination fractions will be derived
            //    from marker positions using the Haldane mapping function.
            // --simwalk2
            //    Perform a smart linkage analysis in conjuction with Simwalk2. MERLIN tackles
            //    the small pedigrees, Simwalk2 does the larger ones
            //    (See notes in analysis.h for building SimWalk2)
            // --minutes:n
            //    Do not attempt to analyse families where calculations for the forward portion
            //    of the Markov-Chain require more than n minutes.
            //    NOTE: a 1 min cutoff is used, those that take longer will be handled by simwalk2.
            if (qtdt) {
                // Pedigree and data files in QTDT format. See: http://www.well.ox.ac.uk/asthma/QTDT
                // QTDT format requires a FIXED WIDTH 3-column map file containing: CHROMOSOME MARKER POSITION
                fprintf(filep,
#ifdef RUNSHELL_SETUP
                        "$%s -p %s -d %s -m %s --simwalk2 --minutes 1\n",
                        merlin_prog_env,
#else /* RUNSHELL_SETUP */
                        "%s -p %s -d %s -m %s --simwalk2 --minutes 1\n",
                        merlin_prog,
#endif /* RUNSHELL_SETUP */
                        file_names[10], file_names[11], file_names[12]);
            } else {
                // Pedigree and data files in Linkage format. See: http://linkage.rockefeller.edu/
                fprintf(filep,
#ifdef RUNSHELL_SETUP
                        "$%s -p %s -d %s --simwalk2 --minutes 1\n",
                        merlin_prog_env,
#else /* RUNSHELL_SETUP */
                        "%s -p %s -d %s --simwalk2 --minutes 1\n",
                        merlin_prog,
#endif /* RUNSHELL_SETUP */
                        file_names[10], file_names[11]);
            }
	    fprintf_status_check_csh(filep, "merlin", 0);

            if (!xlinked) {
                // merlin2sw2.pl: Formats merlin's NPL output for reading into SimWalk2
                // input: npl statistics file (merlin.npl), and order file (file_names[13])
                // exit status -1 on error; 0 OK.
                fprintf(filep,
                        "%s merlin.npl %s merlin_out.%s\n",
                        file_names[14], file_names[13], chr_str);
#ifdef RUNSHELL_SETUP
		fprintf_status_check_csh(filep, "merlin2sw2", 0);
            } else {
                // NOTE: SimWalk2 does not support x-linked analysis, so just say it ran OK.
                fprintf(filep, "echo '0' > merlin2sw2_status\n");
#endif /* RUNSHELL_SETUP */
            }
        }
        if (analysis == HAPLOTYPE)
            fprintf(filep, "echo Commencing haplotyping run for chromosome %2d using SimWalk2\n", numchr);
        if (analysis == LOCATION)
            fprintf(filep, "echo Commencing location score run for chromosome %2d using SimWalk2\n", numchr);
        if (analysis == NONPARAMETRIC ||
            (analysis ==  TO_MERLIN && xlinked == 0))
            fprintf(filep, "echo Commencing nonparametric analysis run for chromosome %2d using SimWalk2\n", numchr);
        if (analysis == IBD_EST)
            fprintf(filep, "echo Commencing IBD estimation run for chromosome %2d using SimWalk2\n", numchr);
        if (analysis == MISTYPING)
            fprintf(filep, "echo Commencing mistyping analysis run for chromosome %2d using SimWalk2\n", numchr);
        fprintf(filep, "echo \n");
        fprintf(filep, "unalias rm\n");
        fprintf(filep, "unalias cp\n");
        if (analysis == HAPLOTYPE) {
            fprintf(filep, "foreach i (0 1 2 3 4 5 6 7)\n");
            fprintf(filep, "\t if (-e ERROR-0$i.TXT) then\n");
            fprintf(filep, "\t\trm -f ERROR-0$i.TXT\n");
            fprintf(filep, "\tendif\n");
            fprintf(filep, "end\n");
            fprintf(filep, "if (-e RERUN-01.BCH) then\n");
            fprintf(filep, "\trm -f RERUN-01.BCH\n");
            fprintf(filep, "endif\n");
        }
/*     if (!xlinked || analysis != TO_MERLIN) { */
        fprintf(filep, "rm -f PEDIGREE.DAT\n");
        fprintf(filep, "rm -f LOCUS.DAT\n");
        fprintf(filep, "rm -f MAP.DAT\n");
        fprintf(filep, "rm -f PEN.DAT\n");
        fprintf(filep, "rm -f BATCH2.DAT\n");
        fprintf(filep, "cp %s PEDIGREE.DAT\n", outfl_name);
        fprintf(filep, "cp %s LOCUS.DAT\n", loutfl_name);
        fprintf(filep, "cp %s MAP.DAT\n", moutfl_name);
        if (num_affec > 1 || (LoopOverTrait == 1 && num_affec > 0)) {
            fprintf(filep, "cp %s PEN.DAT\n", poutfl_name);
        }
        // This is the input control data file which contains the instruction parameters...
        fprintf(filep, "cp %s BATCH2.DAT\n", batfl_name);
        // See notes in analysis.h for instructions on building SimWalk2...
#ifdef RUNSHELL_SETUP
        fprintf(filep, "$_SIMWALK2\n");
        // NOTE: SimWalk2 can be run multiple times (see below). However, we will only
        // check for an error on the first run.
	fprintf_status_check_csh(filep, "simwalk2", 0);
#else /* RUNSHELL_SETUP */
        fprintf(filep, "simwalk2\n");
#endif /* RUNSHELL_SETUP */
/*     } */
        if (analysis == HAPLOTYPE) {
            i=1;
            errf_lines;
            fprintf(filep,
                    "echo SimWalk2 finished the first complete haplotype run.\n");
            /* Rerun loops */
            switch(i) {
            case 1:
                csh_lines;
                fprintf(filep, "echo Simwalk2 suggested a second haplotype run\n");
                fprintf(filep, "echo starting the second haplotype run\n    simwalk2\n");
                else_lines;
                i++;
                errf_lines;
                fprintf(filep, "echo SimWalk2 finished the second complete haplotype run.\n");
            case 2:
                csh_lines;
                fprintf(filep, "echo Simwalk2 suggested a third haplotype run\n");
                fprintf(filep, "echo starting the third haplotype run\n     simwalk2\n");
                else_lines;
                i++;
                errf_lines;
                fprintf(filep, "echo SimWalk2 finished the third complete haplotype run.\n");
            case 3:
                csh_lines;
                fprintf(filep, "echo Simwalk2 suggested a fourth haplotype run\n");
                fprintf(filep, "echo starting the fourth haplotype run\n    simwalk2\n");
                else_lines;
                i++;
                errf_lines;
                fprintf(filep, "echo SimWalk2 finished the fourth complete haplotype run.\n");
            case 4:
                csh_lines;
                fprintf(filep, "echo Simwalk2 suggested a fifth haplotype run\n");
                fprintf(filep, "echo starting the fifth haplotype run\n    simwalk2\n");
                else_lines;
                i++;
                errf_lines;
                fprintf(filep, "echo SimWalk2 finished the fifth complete haplotype run.\n");
            case 5:
                csh_lines;
                fprintf(filep, "echo Simwalk2 suggested a sixth haplotype run\n");
                fprintf(filep, "echo NOT starting the sixth haplotype run\n");
                fprintf(filep, "echo \"*****************\" \n");
                fprintf(filep, "exit 0\n");
                else_lines;
                break;
            default:
                break;
            }
        }
        if (xlinked && analysis == TO_MERLIN) {
            fprintf(filep,
                    "echo SimWalk2 does not support x-linked analysis\n");
        }
        fclose(filep);
        syscmd=CALLOC(strlen(fl_name)+11, char);
        sprintf(syscmd, "chmod +x %s", fl_name);
        System(syscmd);
        free(syscmd);
        if (global) {
            /* write the invokation command into the top-level file */
            if ((filep = fopen(global_shell, "a")) == NULL) {
                errorvf("Could not open %s for writing.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "echo SimWalk2 analysis for chromosome %d", numchr);
            if (tr > 0) {
                fprintf(filep, " ,trait %d ...", tr);
            }
            fprintf(filep, "\n");
            fprintf(filep, "echo This may take a while.\n");
            if (nloop > 1) {
                fprintf(filep, "cd %s\n", trait_paths[tr]);
            }
            fprintf(filep, "./%s\n", cshfl_name);
            if (nloop > 1) {
                fprintf(filep, "cd ..\n");
            }
            fclose(filep);
        }

        if (nloop == 1) break;
    }

    /* local shell scripts are being created */
    sprintf(err_msg,
            "        C-Shell script:       %s", cshfl_name);
    mssgf(err_msg);
}

#undef csh_lines
#undef else_lines
#undef errf_lines


/*----------------- End of write_cshell------------ */

/* Set up a batch file for userm15, in search mode */
static void write_batch(char *loutfl_name, char *outfl_name,
			char *fl_name, linkage_ped_top *LLTop,
			int write_trait)
{
    linkage_locus_top *LTop;
    int            *trp, cont=0, tr, i, num_affec=num_traits;
    FILE           *filep = NULL; //compiler: can not tell filep is assigned
    char           batfl_name[2*FILENAME_LENGTH];
    int first_time = 1;
    /* implementing the loop over traits a little differently,
       if LoopOverTraits is 1, then current file will be closed
       and the next one opened, also the last marker will say yes
       until the last trait is encountered.
    */

    trp = &(global_trait_entries[0]);
    LTop = LLTop->LocusTop;

    for (tr=0; tr < num_affec; tr++) {
        if (LoopOverTrait == 1 && num_affec > 1) {
            sprintf(batfl_name, "%s/%s", output_paths[tr+1], fl_name);
            if ((filep = fopen(batfl_name, "w")) == NULL) {
                errorvf("Could not open M15 batch file %s.\n", batfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
        } else {
            if (*trp == -1) {
                trp++;
                continue;
            }
/*       if (LLTop->LocusTop->Locus[*trp].Type != AFFECTION) { */
/* 	continue; */
/*       } */
            if (first_time) {
                first_time=0;
                /* loop being executed the first time, open file
                   in current directory */
                sprintf(batfl_name, "%s/%s", output_paths[0], fl_name);
                if ((filep = fopen(batfl_name, "w")) == NULL) {
                    errorvf("Could not open M15 batch file %s.\n", batfl_name);
                    EXIT(FILE_WRITE_ERROR);
                }
            }
            /* if tr > 0, file is already open */
        }

        fprintf(filep, "15\n");  /* Switch from grid to search */
        fprintf(filep, "2\n");
        fprintf(filep, "%s\n",loutfl_name);  /* locus file name */
        fprintf(filep, "3\n");
        fprintf(filep, "%s\n",outfl_name);  /* pedigree file name */
        fprintf(filep, "8\n");   /* Extra array */
        fprintf(filep, "%d\n", ((LoopOverTrait == 1)? 1: *trp+1));
        /* USERM15 needs to be told which one is the trait */
        fprintf(filep, "\n");
        /*    } */

        for (i = 0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Type != AFFECTION &&
                LTop->Locus[ChrLoci[i]].Type != QUANT) {
                if (cont==1) {
                    fprintf(filep, "yes\n");  /* yes wrote a previous run*/
                }
                fprintf(filep, "9\n");   /* Option 9 */
                /* Have the loci names been truncated by now? */
                fprintf(filep, "%s\n", strtail(LTop->Locus[*trp].Name, MENDEL_MAX_LOCUS_NAME_LEN));
                /* Name of trait locus */
                fprintf(filep, "%s\n", strtail(LTop->Locus[ChrLoci[i]].Name, MENDEL_MAX_LOCUS_NAME_LEN));
                /* Name of marker locus */
                fprintf(filep, "\n");
                fprintf(filep, "21\n");
                cont = 1;
            }
        }
        if (LoopOverTrait == 1) {
            fprintf(filep, "no\n");  /* no more runs in this file */
            fclose(filep);
            cont = 0;
        }
        trp++;
    }

    if (LoopOverTrait == 0) {
        fprintf(filep, "no\n");  /* no more runs in this file */
        fclose(filep);

    }

    sprintf(err_msg, "        M15 batch file:       %s", fl_name);
    mssgf(err_msg);
    return;
}

/*----------------- end of write_batch------------------- */

static  void  write_batch2(char *fl_name, int numchr,
			   linkage_ped_top *LLTop,
			   analysis_type analysis)
{

    FILE           *filep;
    linkage_locus_top *LTop;
    int             nloop, num_affec=num_traits, tr, i,j;
/*   double           difff, diffm; */
#ifndef HIDEDATE
    time_t          now;
#endif
    char            chr_str[3], comment_char[5], batfl_name[2*FILENAME_LENGTH];
    int num_markers, *markers;

    CHR_STR(numchr, chr_str);
    strcpy(comment_char, "!");
    LTop = LLTop->LocusTop;

    if (num_traits > 0) {
        if (LoopOverTrait == 0) {
            num_markers = NumChrLoci - num_traits + 1;
        } else {
            num_markers = NumChrLoci - num_traits;
        }
    } else {
        num_markers = NumChrLoci;
    }

    markers = CALLOC((size_t) num_markers, int);
    j=0;
    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY) {
            markers[j] = ChrLoci[i]; j++;
        }
    }

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(batfl_name, "%s/%s", output_paths[tr], fl_name);
        if ((filep = fopen(batfl_name, "w")) == NULL) {
            errorvf("Could not open file %s for writing.\n", batfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "\n");
        fprintf(filep, "000001                 %s batch item number\n", comment_char);

        if (analysis == HAPLOTYPE)
            fprintf(filep, "1                      %s", comment_char);

        else if (analysis == LOCATION)
            fprintf(filep, "2                      %s", comment_char);

        else if (analysis == NONPARAMETRIC || analysis == TO_MERLIN)
            fprintf(filep, "3                      %s", comment_char);

        else if (analysis == IBD_EST)
            fprintf(filep, "4                      %s", comment_char);

        else if (analysis == MISTYPING)
            fprintf(filep, "5                      %s", comment_char);

        fprintf(filep,
                " Analysis: 1=Haplo, 2=LOD, 3=NPL, 4=IBD, 5=Mistyping\n");
        fprintf(filep, "\n");
        fprintf(filep, "000002                 %s batch item number\n", comment_char);
        ((analysis == HAPLOTYPE)?
         fprintf(filep, "1                      %s integer label for this run of the program\n", comment_char):
         fprintf(filep, "%-2d                     %s integer label for this run of the program\n",
                 numchr, comment_char));
        fprintf(filep, "\n");
        fprintf(filep, "000003                 %s batch item number\n", comment_char);
#ifndef HIDEDATE
        now = time(NULL);
#endif
        fprintf(filep, "Chr %2d; Infiles: %s", numchr,
#ifdef HIDEDATE
                NOcTIME
#else
                ctime(&now)
#endif
);
        fprintf(filep, "\n");
        fprintf(filep, "000004                 %s batch item number\n", comment_char);
        fprintf(filep, "N                      %s continuation of a previous job?\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000005                 %s batch item number\n", comment_char);
        fprintf(filep, "N                      %s run silently?\n", comment_char);
        fprintf(filep, "\n");
/*      fprintf(filep, "000009                 %s batch item number\n", comment_char); */
/*      fprintf(filep, "MAP.DAT                %s name of map file\n", comment_char); */
/*      fprintf(filep, "\n"); */
        fprintf(filep, "000010                 %s batch item number\n", comment_char);
        fprintf(filep, "LOCUS.DAT              %s name of locus file\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000011                 %s batch item number\n", comment_char);
        fprintf(filep, "PEDIGREE.DAT           %s name of pedigree file\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000012                 %s batch item number\n", comment_char);
        fprintf(filep, "F                      %s symbol for female (case insensitive)\n", comment_char);
        fprintf(filep, "M                      %s symbol for   male (case insensitive)\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000013                 %s batch item number\n", comment_char);
        if ((num_affec >=1 && LoopOverTrait == 1) ||
            (num_affec >= 2 && LoopOverTrait == 0))  {
            fprintf(filep, "Y                      %s is trait listed in locus and pedigree files?\n", comment_char);
        } else  {
            fprintf(filep, "N                      %s is trait listed in locus and pedigree files?\n", comment_char);
        }
        fprintf(filep, "1.00                   %s the a-priori proportion of linked pedigrees\n", comment_char);
        fprintf(filep, "\n");

/*     fprintf(filep, "000014                 %s batch item number\n", comment_char); */
/*     fprintf(filep, "%6d                 %s number of marker loci\n",  */
/* 	    num_markers, comment_char); */
/*     for (i = 0; i < num_markers; i++) { */
/*       fprintf(filep, "%6d                 %s genomic position for marker in position %d\n",  */
/* 	      i + 1, comment_char, i + 1);  */
/*     } */
/*     fprintf(filep, "\n"); */
/*     fprintf(filep, "000015                 %s batch item number (only REQUIRED batch item)\n", comment_char); */
/*     fprintf(filep, "%6d                 %s number of marker loci\n",  */
/* 	    num_markers, comment_char); */

/*     /\* Check if marker distances are monotonic  */
/*        Assume that the num_affec trait loci are first in order ... *\/ */
/*     for (i = 0; i < num_markers-1; i++)  { */
/*       if (analysis == LOCATION) { */
/* 	/\* use the sex-averaged map *\/ */
/* 	diffm = difff =  */
/* 	  LTop->Marker[markers[i+1]].pos_avg - LTop->Marker[markers[i]].pos_avg; */
/*       } */
/*       else { */
/* 	diffm = LTop->Marker[markers[i+1]].pos_male -  */
/* 	  LTop->Marker[markers[i]].pos_male; */
/* 	difff = LTop->Marker[markers[i+1]].pos_female - */
/* 	  LTop->Marker[markers[i]].pos_female; */
/*       } */
/*       if (diffm < 0.0) {	diffm=0.0999; } */
/*       if (difff < 0.0) {	difff=0.0999; } */
/*       else { */
/* 	if (LTop->map_distance_type == 'k') { */
/* 	  diffm=kosambi_theta(diffm); */
/* 	  difff=kosambi_theta(difff); */
/* 	} */
/* 	else { */
/* 	  diffm=haldane_theta(diffm); */
/* 	  difff=haldane_theta(difff); */
/* 	} */
/*       } */
/*       fprintf(filep,  */
/* 	      " %7.5f %7.5f       %s recomb freq in interval %d for female & male\n",  */
/* 		difff, diffm, comment_char, i+1); */
/*       j=i;  */
/*     } */
/*     fprintf(filep, "\n"); */
        fprintf(filep, "000016                 %s batch item number\n", comment_char);
        fprintf(filep, "2                      %s label indicating an affected individual\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000017                 %s batch item number\n", comment_char);
        if ((num_affec > 0 && LoopOverTrait ==1) ||
           (num_affec > 1 && LoopOverTrait == 0)) {
            fprintf(filep, "0                      %s affected status in:0=>trait locus; n=>n-th qant var\n",
                    comment_char);
        } else {
            fprintf(filep, "-1                     %s affected status in:0=>trait locus; n=>n-th quant var\n", comment_char);
        }
        fprintf(filep, "\n");
        fprintf(filep, "000018                 %s batch item number\n", comment_char);
        if ((num_affec > 0 && LoopOverTrait ==1) ||
            (num_affec > 1 && LoopOverTrait == 0))
            fprintf(filep, "1                      ");
        else
            fprintf(filep, "0                      ");
        fprintf(filep, "%s number of quantitative variables in pedigree\n",
                comment_char);
        fprintf(filep, "\n");

        if ((num_affec > 0 && LoopOverTrait ==1) ||
            (num_affec > 1 && LoopOverTrait == 0)) {
            fprintf(filep, "000019                 %s batch item number\n",
                    comment_char);
            fprintf(filep, "PEN.DAT                %s Name of penetrance file; blank impies no file\n", comment_char);
            fprintf(filep, "\n");
        }
        fprintf(filep, "000020                 %s batch item number\n", comment_char);
        fprintf(filep, "1000                   %s number of pedigrees sampled for each input pedigree\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000021                 %s batch item number {for simulated annealing}\n", comment_char);
        fprintf(filep, "1                      %s number of parallel runs to perform\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000022                 %s batch item number {for simulated annealing}\n", comment_char);
        fprintf(filep, "100.0                  %s temperature at start of simulated annealing\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000023                 %s batch item number {for simulated annealing}\n", comment_char);
        fprintf(filep, "0.99                   %s factor by which the temperature changes\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000024                 %s batch item number {for simulated annealing}\n", comment_char);
        fprintf(filep, "800                    %s number of temperature changes\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000025                 %s batch item number {for simulated annealing}\n", comment_char);
        fprintf(filep, "0                      %s number of pre-simulated annealing steps\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "00026                  %s batch item number\n", comment_char);
        fprintf(filep, "64                     %s maximum number of attempts to find a legal state\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000029                 %s batch item number\n", comment_char);
        fprintf(filep, "1000                   %s min. # of steps between sampled peds in sim. annealing\n", comment_char);
        fprintf(filep, "1000                   %s min. # of steps between sampled peds in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000030                 %s batch item number\n", comment_char);
        fprintf(filep, "10                     %s factor for number of steps in sim. annealing\n", comment_char);
        fprintf(filep, "10                     %s factor for number of steps in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000031                 %s batch item number\n", comment_char);
        fprintf(filep, "2.00                   %s mean # of transitions/step in sim. annealing\n", comment_char);
        fprintf(filep, "2.00                   %s mean # of transitions/step in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000032                 %s batch item number\n", comment_char);
        fprintf(filep, "10                     %s weight for untyped people in sim. annealing\n", comment_char);
        fprintf(filep, "10                     %s weight for untyped people in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000033                 %s batch item number\n", comment_char);
        fprintf(filep, "0.50                   %s weight for neighbor-pivot in sim. annealing\n", comment_char);
        fprintf(filep, "0.50                   %s weight for neighbor-pivot in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000034                 %s batch item number\n", comment_char);
        fprintf(filep, "0.25                   %s freq of transition type 2 in sim. annealing\n", comment_char);
        fprintf(filep, "0.25                   %s freq of transition type 2 in random walk\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000035                 %s batch item number\n", comment_char);
        fprintf(filep, "0                      %s number of steps between 'free' runs\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000036                 %s batch item number\n", comment_char);
        fprintf(filep, "27713                  %s random seed: positive integer < 30000\n", comment_char);
        fprintf(filep, " 2321                  %s random seed: positive integer < 30000\n", comment_char);
        fprintf(filep, "18777                  %s random seed: positive integer < 30000\n", comment_char);
        fprintf(filep, "\n");
/*      fprintf(filep, "000037                 %s batch item number\n", comment_char); */
/*      fprintf(filep, "0                      %s error model used, 0=>uniform, 1=>empirical\n", comment_char); */
/*      fprintf(filep, "0.025                  %s overall error rate\n", comment_char); */
/*      fprintf(filep, "\n"); */
        if (analysis == TO_MERLIN) {
            fprintf(filep, "000038                 %s batch item for exact npl scores\n",
                    comment_char);
            fprintf(filep, "merlin_out.%s          %s exact npl score file name\n",
                    chr_str, comment_char);
            fprintf(filep, "\n");
        }
        fprintf(filep, "000040                 %s batch item number\n", comment_char);
        fprintf(filep, "N                      %s create the INPED-nn.mmm files?\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000041                 %s batch item number {only for location scores}\n", comment_char);
        fprintf(filep, "Y                      %s create the SCORE-nn.mmm files?\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000043                 %s batch item number {only for haplotyping}\n", comment_char);
        fprintf(filep, "Y                      %s create the haplotype export files?\n", comment_char);
        fprintf(filep, "\n");
        /*
          Menu item #44(1) has changed slightly:  call it HAPLO_THRESHOLD.
          Menu item #44(2) is new (but optional): call it RECOMB_THRESHOLD.
          A pedigree is now included in RERUN-nn.PED if there are at least
          HAPLO_THRESHOLD number of haplotypes with RECOMB_THRESHOLD number of
          recombinants in each OR if there are any haplotypes with more than
          RECOMB_THRESHOLD number of recombinants. Default value is 2 for both.
        */
        fprintf(filep, "000044                 %s batch item number {only for haplotyping}\n", comment_char);
        fprintf(filep, "0                      %s Haplotype threshold for rerunning a pedigree\n", comment_char);
        fprintf(filep, "0                      %s Recombination threshold for rerunning a pedigree\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000045                 %s batch item number {only for haplotyping}\n", comment_char);
        fprintf(filep, "0.85                   %s Step count threshold for rerunning a pedigree\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000046                 %s batch item number {only for haplotyping}\n", comment_char);
        fprintf(filep, "N                      %s haplotypes displayed vertically?\n", comment_char);
        fprintf(filep, "N                      %s haplotypes displayed using founder-source label?\n", comment_char);
        fprintf(filep, "N                      %s allele specified even when no information available?\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000047                 %s batch item number {only for sampling}\n", comment_char);
        fprintf(filep, "N                      %s output the sampled pedigrees in LINKAGE-format?\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000048                 %s batch item number {only for NPL statistics}\n", comment_char);
        fprintf(filep, "10000                  %s number of unconditional simulations for p-value\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000049                 %s batch item number {only for IBD analysis}\n", comment_char);
        fprintf(filep, "4                      %s number of points analyzed between markers\n", comment_char);
        fprintf(filep, "\n");
        fprintf(filep, "000050                 %s batch item number signifying end of data\n", comment_char);
        fclose(filep);
        if (nloop == 1) break;
    }
    sprintf(err_msg, "        Simwalk2 batch file:  %s", fl_name);
    mssgf(err_msg);
    return;
}
/*-----------------End of write_batch2-------------------- */


// Missing quantitative data is not written to the output file...
static void mwrite_quantitative_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON)
        fprintf(filep, "        ");
    else
        fprintf(filep, "%7.3f ", entry->Pheno[locusnm].Quant);
}

static void mwrite_affection_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    if (entry->Pheno[locusnm].Affection.Status > 0)
        fprintf(filep, "%1d", entry->Pheno[locusnm].Affection.Status);
    else
        fprintf(filep, " ");
    if (locus->Pheno->Props.Affection.ClassCnt > 1)  {
        if (entry->Pheno[locusnm].Affection.Status > 0)
            fprintf(filep, "%7d", entry->Pheno[locusnm].Affection.Class);
        else
            fprintf(filep, "       "); /*  entry->Pheno[locusnm].Affection.Class); */
    }
    else
        fprintf(filep, "       ");
}

static void            m5write_affection_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
/*   FILE           *filep;
     int             locusnm;
     linkage_locus_rec *locus;
     linkage_ped_rec *entry; */
{
    char trait_phen[8];

    if (entry->Pheno[locusnm].Affection.Status == 0) {
        fprintf(filep, "        ");
    } else {
        if (locus->Pheno->Props.Affection.ClassCnt > 1)  {
            sprintf(trait_phen, "%d_%d ", entry->Pheno[locusnm].Affection.Status,
                    entry->Pheno[locusnm].Affection.Class);
            fprintf(filep, "%7s ", trait_phen);
        } else {
            fprintf(filep, "%7d ", entry->Pheno[locusnm].Affection.Status);
        }
    }
}

static  void            swrite_affection_data(FILE *filep,
					      int locusnm,
					      linkage_locus_rec *locus,
					      linkage_ped_rec *entry)

{
    int affected;


    affected = aff_status_entry(entry->Pheno[locusnm].Affection.Status, entry->Pheno[locusnm].Affection.Class,
                                locus);
    entry->Orig_status=affected;
    if (affected == 0) {
        fprintf(filep, "    ");
    } else {
        fprintf(filep, "%4d", affected);
    }
    fprintf(filep, "    "); /*  entry->Pheno[locusnm].Affection.Class); */
}

static void lwrite_affection_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int aff_status;

    aff_status=entry->Pheno[locusnm].Affection.Status;
    entry->Orig_status=aff_status;
    if (aff_status == 0)
        fprintf(filep, "    ");
    else
        fprintf(filep, "%4d", aff_status);
    fprintf(filep, "    ");
}

static void mwrite_binary_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int a1, a2;
    int allele;
    /* Need to check for completely untyped person */
    for (allele = 0; allele < locus->AlleleCnt; allele++)
        {
            get_2alleles(entry->Marker, locusnm, &a1, &a2);
            if ((allele == a1 - 1)
                || (allele == a2 - 1))
                fprintf(filep, "1");
            else
                fprintf(filep, "0");
            for (allele = locus->AlleleCnt; allele < 7; allele++)
                fprintf(filep, " ");
        }
}


static void mwrite_numbered_data(FILE *filep, const int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int a1, a2;
    get_2alleles(entry->Marker, locusnm, &a1, &a2);
    if (a1 == 0)
        fputs("        ", filep);
    else
        fprintf(filep, "%3s/%3s ", format_allele(locus, a1), format_allele(locus, a2));
}

/* This is a hack, copied from create_formats,
   in order to maximise the allowable width (7), don't print a space
*/
void create_mendel_formats(int famwid, int perwid,
			   char *fformat, char *pformat)

{

    if (fformat != NULL) {
        if (OrigIds[1] == 2 || OrigIds[1] == 4) {
            sprintf(fformat, "%%%ds", famwid);
        } else {
            sprintf(fformat, "%%%dd", famwid);
        }
    }

    if (pformat != NULL) {
        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
            sprintf(pformat, "%%%ds", perwid);
        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
            sprintf(pformat, "%%%ds", perwid);
        } else {
            sprintf(pformat, "%%%dd", perwid);
        }
    }
    return;

}

/*
 * Save all the pedigrees to a file in MENDEL format. Return the
 * number of members saved (total).
 * Since we are allowing the original entry ids to be written as
 * strings, these have to be checked and shortened if they are more
 * than 7 characters long. If we are writing the oridinal ids,
 * then this check is not being done as presumably, we have far
 * fewer than 10**8 entries.
 * Pedigree names are also being checked for the name length.
 *        NM, Aug 25, 2000

 * NM, using a pointer to the correct position in ID strings,
 * instead of modifying them.
 */


static void  save_mendel_peds(char *outfl_name, linkage_ped_top *Top,
			      analysis_type analysis)
{
#define LTOP Top->LocusTop

    linkage_locus_rec *Locus;
    int      nloop, tr, i, ped, entry, locus1, *trp, loc;
    int      num_affec = num_traits, naff=num_affection, *aff_status = NULL, *paff_status;
    linkage_ped_rec *Entry, *TEntry;
    int      fwid, pwid;
    char     ofl[2*FILENAME_LENGTH],  fformat[6], pformat[6];
    FILE     *filep;
    int      loc1, j;

    /* We already have the number of affection status loci in a global */
    /* set the number of affection loci to decide how many quantitative
       variables need to be written at the end */

    linkage_ped_tree *TailPed = CALLOC((size_t) Top->PedCnt, linkage_ped_tree);
    for (ped = 0; ped < Top->PedCnt; ped++)  {
        strcpy(TailPed[ped].Name, strtail(Top->Ped[ped].Name, MENDEL_MAX_PEDIGREE_NAME_LEN));

        TailPed[ped].EntryCnt=Top->Ped[ped].EntryCnt;
        TailPed[ped].Entry = CALLOC((size_t) TailPed[ped].EntryCnt, linkage_ped_rec);
        for (entry = 0; entry < TailPed[ped].EntryCnt; entry++)  {
            TEntry = &(TailPed[ped].Entry[entry]);
             Entry = &(Top->Ped[ped].Entry[entry]);

            strcpy(TEntry->OrigID,   strtail(Entry->OrigID,   MENDEL_MAX_PEDIGREE_NAME_LEN));
            strcpy(TEntry->UniqueID, strtail(Entry->UniqueID, MENDEL_MAX_PEDIGREE_NAME_LEN));
//          strcpy(TEntry->PerPre,   strtail(Entry->PerPre,   MENDEL_MAX_PEDIGREE_NAME_LEN));
            TEntry->ID     = Entry->ID;
            TEntry->Father = Entry->Father;
            TEntry->Mother = Entry->Mother;
        }
    }

    if (analysis != TO_MENDEL4) {
        if (naff > 0 && LoopOverTrait == 0) {
            aff_status = CALLOC((size_t) naff, int);
        } else {
            aff_status = CALLOC((size_t) 1, int);
        }
    }

    NLOOP;

    trp = (num_affec > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;

        sprintf(ofl, "%s/%s", output_paths[tr], outfl_name);
        if ((filep = fopen(ofl, "a")) == NULL) {
            errorvf("Could not open file %s for writing.\n", ofl);
            EXIT(FILE_WRITE_ERROR);
        }

        field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);

        if (fwid >= 7) fwid=8;
        else fwid = fwid+1;
        if (pwid >= 8) pwid=8;
        else pwid = pwid+1;

        create_mendel_formats(fwid, pwid, fformat, pformat);

        fprintf(filep, "(I5,1x,A%d)\n", fwid);
        /*  If we want a line of length 70, then we can have 6 phenotype fields */
        fprintf(filep, "(3A%d,1X,2A1,(T%d,3(A8),:))\n", pwid, 3*pwid+4);

        for (ped = 0; ped < Top->PedCnt; ped++)  {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            /* Write total number of persons and the pedigree id */
            fprintf(filep, "%5d ", Top->Ped[ped].EntryCnt);
            prID_ped(filep, ped, fformat, &TailPed[ped], "");
            fprintf(filep, "\n");

            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++)  {
                Entry = &(Top->Ped[ped].Entry[entry]);
                /* write the entry number then parents */

                prID_fam(filep, pformat, &TailPed[ped].Entry[entry], TailPed[ped].Entry, "");
                if (Entry->Father == 0) {
                    for (i=0; i < 2*pwid; i++) {
                        fprintf(filep, " ");
                    }
                }

                /* write the sex */
                if (Entry->Sex == MALE_ID)	   fprintf(filep, " M ");
                if (Entry->Sex == FEMALE_ID)      fprintf(filep, " F ");
                if (Entry->Sex != MALE_ID && Entry->Sex != FEMALE_ID)
                    fprintf(filep, " 0 ");

                /* first write the trait locus only if LoopOverTrait is 1 */
                loc=0;
                if (LoopOverTrait == 1 && num_affec > 0) {
                    locus1 = *trp;
                    if (Top->LocusTop->Locus[locus1].Type == AFFECTION) {
                        if (analysis == LOCATION)
                            lwrite_affection_data(filep, locus1,
                                                  &(Top->LocusTop->Locus[locus1]),
                                                  Entry);

                        else if (analysis == HAPLOTYPE || analysis == NONPARAMETRIC ||
                                 analysis == TO_MERLIN || analysis == IBD_EST ||
                                 analysis == MISTYPING)
                            swrite_affection_data(filep, locus1,
                                                  &(Top->LocusTop->Locus[locus1]), Entry);

                        else if (analysis == TO_MENDEL)
                            mwrite_affection_data(filep, locus1,
                                                  &(Top->LocusTop->Locus[locus1]), Entry);

                        else if (analysis == TO_MENDEL4)
                            m5write_affection_data(filep, locus1,
                                                   &(Top->LocusTop->Locus[locus1]), Entry);

                        if (analysis != TO_MENDEL4) {
                            aff_status[0] =
                                get_aff_status(&(Top->LocusTop->Locus[locus1]),
                                               analysis, Entry, locus1);
                        }
                        loc=1;
                    }
                }

                /* now write the genotypes and QTLs */

//xx            if (naff > 0) paff_status=&(aff_status[0]);
                paff_status=&(aff_status[0]);

                for (loc1 = 0; loc1 < NumChrLoci; loc1++)   {
                    locus1 = ChrLoci[loc1];
                    switch (Top->LocusTop->Locus[locus1].Type)  {
                    case QUANT:
                        if (LoopOverTrait==0) {
                            TAB_TO_GENO_COL
                                mwrite_quantitative_data(filep, locus1,
                                                         &(Top->LocusTop->Locus[locus1]), Entry);
                            loc++;
                        }
                        break;
                    case NUMBERED:
                        TAB_TO_GENO_COL
                            mwrite_numbered_data(filep, locus1,
                                                 &(Top->LocusTop->Locus[locus1]), Entry);
                        loc++;
                        break;
                    case BINARY:
                        TAB_TO_GENO_COL
                            mwrite_binary_data(filep, locus1, &(Top->LocusTop->Locus[locus1]),
                                               Entry);
                        loc++;
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            TAB_TO_GENO_COL

                                if (analysis == LOCATION)
                                    lwrite_affection_data(filep, locus1,
                                                          &(Top->LocusTop->Locus[locus1]),
                                                          Entry);

                                else if (analysis == HAPLOTYPE || analysis == NONPARAMETRIC ||
                                         analysis == TO_MERLIN || analysis == IBD_EST ||
                                         analysis == MISTYPING)
                                    swrite_affection_data(filep, locus1,
                                                          &(Top->LocusTop->Locus[locus1]),
                                                          Entry);

                                else if (analysis == TO_MENDEL)
                                    mwrite_affection_data(filep, locus1,
                                                          &(Top->LocusTop->Locus[locus1]),
                                                          Entry);

                                else if (analysis == TO_MENDEL4)
                                    m5write_affection_data(filep, locus1,
                                                           &(Top->LocusTop->Locus[locus1]),
                                                           Entry);

                            if (analysis != TO_MENDEL4) {
                                *paff_status =
                                    get_aff_status(&(Top->LocusTop->Locus[locus1]),
                                                   analysis, Entry, locus1);
                                paff_status++;
                            }
                            loc++;
                        }
                        break;
                    default:
                        fprintf(stderr,
                                "unknown locus type (locus %d) while writing output\n",
                                locus1 + 1);
                    } /* switch(locus type) */
                } /* End of loop locus1 to Top->LocusTop->LocusCnt */

                /* write actual quantitative phenotypes as well as the aff-related
                   quants, quant only present for the Mendel options anyway */
                if (LoopOverTrait == 1) {
                    /* write the quant for all options as usual,
                       do not write the aff-quant for Mendel 4 */

                    locus1 = *trp;
                    Locus = &(Top->LocusTop->Locus[*trp]);

                    if (Locus->Type == QUANT) {
                        /* If the trait was quantitative, write the phenotype */
                        TAB_TO_GENO_COL
                            mwrite_quantitative_data(filep, locus1,
                                                     &(Top->LocusTop->Locus[locus1]), Entry);
                    } else {
                        if (analysis != TO_MENDEL4) {
                            TAB_TO_GENO_COL
                                if (aff_status[0] == 0) {
                                    fprintf(filep, "        ");
                                } else {
                                    fprintf(filep,"%8d", aff_status[0]);
                                }
                        }
                    }
                } else {
                    /* All options besides Mendel4 and combine traits,
                       real quant phenos have already been written in the genotype
                       loop above, just have to write the aff-quant phenos for non-Mendel4
                    */
                    if (analysis != TO_MENDEL4) {
                        for (j=0; j < naff; j++) {
                            TAB_TO_GENO_COL
                                if (aff_status[j] == 0) {
                                    fprintf(filep, "        ");
                                } else {
                                    fprintf(filep,"%8d", aff_status[j]);
                                }
                            loc++;
                        }
                    }
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;     trp++;
    }
    sprintf(err_msg, "        Pedigree file:        %s", outfl_name);
    mssgf(err_msg);
    free(aff_status);

    for (ped = 0; ped < Top->PedCnt; ped++)  {
        free(TailPed[ped].Entry);
    }
    free(TailPed);

    return;
}

/*------------------ End of save_mendel_peds------------ */



/*
 * Write the locus data file.
 */
void            write_mendel_locus_file(char *file_name,
					linkage_locus_top *LTop,
					analysis_type analysis,
					int sex_linked)
{

    int      nloop, tr, locus1, allele, tmpi, tmpi2, *trp;
    register linkage_locus_rec *Locus;
    char     loutfl_name[2*FILENAME_LENGTH];
    int      numgen, ipen, num_affec = num_traits;
    int      simwalk2;
    FILE     *fp;
    char trait_phen[8];

    simwalk2  = (SIMWALK2(analysis)? 1 : 0);
    NLOOP;

//  if (num_affec > 0) trp = &(global_trait_entries[0]);
    trp = global_trait_entries;
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(loutfl_name, "%s/%s", output_paths[tr], file_name);
        if ((fp = fopen(loutfl_name, "a")) == NULL) {
            errorvf("Could not open locus file %s.\n", loutfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        /* first write the <tr>th affection locus */
        if (LoopOverTrait == 1 && num_affec > 0) {
            /* this pointer is advanced only if there are multiple
               traits being looped over. */
            Locus = &(LTop->Locus[*trp]);   trp++;
            if (Locus->Type == AFFECTION && analysis != TO_HWETEST) {
                MENDEL_AFFECTION
            } else if (Locus->Type == QUANT && analysis == TO_MENDEL4) {
                MENDEL_QUANT
            }
        }
        /* LOCUS LIST MOD */
        /* Now write the numbered and binary loci */
        for (locus1 = 0; locus1 < NumChrLoci; locus1++) {
            Locus = &(LTop->Locus[ChrLoci[locus1]]);
            if (Locus->Type == NUMBERED || Locus->Type == BINARY) {
                if (Locus->Name != NULL) {
                    fprintf(fp, "%-8s", strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN));
                }
                else
                    fprintf(fp, "%8d ", locus1 + 1);
                if (sex_linked == 1) {
                    fprintf(fp, "X-LINKED");
                } else if (sex_linked == 2 &&
                         Locus->Marker->chromosome == SEX_CHROMOSOME) {
                    fprintf(fp, "X-LINKED");
                } else {
                    fprintf(fp, "AUTOSOME");
                }
                fprintf(fp, "%2d", Locus->AlleleCnt);
            }
            switch (Locus->Type)  {
            case BINARY:
                fprintf(fp, "%d\n", Locus->Marker->Props.Binary.FactorCnt);
                for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++)  {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(fp, " %d", (int) Locus->Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', fp);
                }
                break;
            case NUMBERED:
                fprintf(fp, " 0%3d %8f <= # Alleles,# Phenotypes. %5.2f cM\n",
                        Locus->Marker->chromosome, (Locus->Marker->pos_avg/100.0), Locus->Marker->pos_avg);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(fp, "%7s %8f\n",
			    format_allele(Locus, allele+1),
			    Locus->Allele[allele].Frequency);
                }
                break;
            case AFFECTION:
                if (LoopOverTrait == 0 && analysis != TO_HWETEST) {
                    MENDEL_AFFECTION
                        }
                break;
            case QUANT:
                /* do nothing */
                break;
            default:
                errorvf("Unknown locus type %d\n", (int) Locus->Type);
                EXIT(DATA_TYPE_ERROR);
                break;
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
    if (analysis != TO_HWETEST) {
        sprintf(err_msg, "        Locus file:           %s", file_name);
        mssgf(err_msg);
    }
    return;
}
/*------------------ End of write_mendel_locus_file------------ */


/* Write out a pen.dat file in the proper format for use by USERM15 and
   by SIMWALK2 */

static int  write_mendel_pen_file(char *fl_name, linkage_locus_top * LTop,
                                  analysis_type analysis)
{
    int       tr, tmpi2, num_affec=num_traits;
    register  linkage_locus_rec *Locus;
    int       ipen, i, first_time=1;
    FILE      *filep = NULL; //compiler can't figure ...
    char      poutfl_name[2*FILENAME_LENGTH];
    int       naff=0;


    /* TRP SETTING MOD */
    for (tr=0; tr < num_affec; tr++) {
        if (global_trait_entries[tr] < 0) {
            continue;
        } else if (LTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
            naff++;
        }
    }

    for (tr=0; tr < num_affec; tr++) {
        if (LoopOverTrait == 1 && num_affec > 1) {
            sprintf(poutfl_name, "%s/%s", output_paths[tr+1], fl_name);
            if ((filep = fopen(poutfl_name, "w")) == NULL) {
                errorvf("Could not open penetrance file %s.\n", poutfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            /* first the affection locus */
            fprintf(filep,
                    "    1                            <= Number of Affection Status locus\n");
        } else {
            if (global_trait_entries[tr] < 0) continue;
            if (first_time) {
                /* loop being executed the first time, open file
                   in current directory */
                sprintf(poutfl_name, "%s/%s", output_paths[0], fl_name);
                first_time=0;
                if ((filep = fopen(poutfl_name, "w")) == NULL) {
                    errorvf("Could not open penetrance file %s.\n", poutfl_name);
                    EXIT(FILE_WRITE_ERROR);
                }
                fprintf(filep,
                        "%5d                            <= Number of Affection Status locus\n",
                        naff);
            }
            /* if tr > 0, file is already open */
        }

        Locus = &(LTop->Locus[global_trait_entries[tr]]);

        if (Locus->Type != AFFECTION) {
            continue;
        }

        ((Locus->Name != NULL) ?
         fprintf(filep, "%-8s", strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN)) :
         fprintf(filep, "%8d", global_trait_entries[tr] + 1));
        fprintf(filep, "%8d                    <= Trait name, Number of liability classes\n",
                Locus->Pheno->Props.Affection.ClassCnt);
        for (tmpi2 = 1; tmpi2 <= Locus->Pheno->Props.Affection.ClassCnt; tmpi2++) {
            for (i = 1; i <= 2; i++) {
                fprintf(filep, "%8d", ((Locus->Pheno->Props.Affection.ClassCnt > 1) ?
                                       (liability_multiplier*i + tmpi2) : i));
                if (i == 1)    {
                    /* Use 1-Pen */
                    for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++)
                        fprintf(filep, "%8.5f",
                                (1.0 - Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen]));
                } else {
                    /* Use Pen */
                    for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++)
                        fprintf(filep, "%8.5f", (Locus->Pheno->Props.Affection.Class[tmpi2 - 1].AutoPen[ipen]));
                }
                fprintf(filep, " <= Penetrances : Affection = %d & Liability = %d\n",
                        i, tmpi2);
            }
        }
        if (LoopOverTrait == 1) {
            fclose(filep);
        }
    }
    if (LoopOverTrait == 0) {
        fclose(filep);
    }
    sprintf(err_msg, "        Penetrance file:      %s", fl_name);
    mssgf(err_msg);
    return 1;
}

/*------------------ End of write_mendel_pen_file------------ */

/* Mendel 5 penetrance file */

/* NOTE: Mega2's Mendel options assume that there are three
   penetrances because there are two alleles, but the linkage
   locus file uses the standard marker locus-type number of
   penetrances N*(N+1)/2.
   Maintaining the assumtion of 2 alleles, 3 genotypes for now
*/

static int write_mendel5_pen_file(char *fl_name,
				  linkage_locus_top *LTop)
{

    int   *trp, tr, tmpi2, ipen, num_affec=num_traits;
    register linkage_locus_rec *Locus;
    FILE *filep = NULL; //compiler: too hard to tell its assigned.
    char poutfl_name[2*FILENAME_LENGTH];
    int naff=0, first_time=1;
    char trait_phen[8];

    for(tr=0; tr < num_affec; tr++) {
        if (global_trait_entries[tr] < 0) {
            continue;
        } else if (LTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
            naff++;
        }
    }

    /* TRP SETTING MOD */
    if (naff == 0) return 0;

    trp=&(global_trait_entries[0]);
    for (tr=0; tr < num_affec; tr++) {
        if (*trp == -1) {
            trp++;
            continue;
        }
        Locus = &(LTop->Locus[*trp]); trp++;
        if (Locus->Type != AFFECTION) {
            continue;
        }

        if (LoopOverTrait == 1 && num_affec > 1) {
            sprintf(poutfl_name, "%s/%s", output_paths[tr+1], fl_name);
            if ((filep = fopen(poutfl_name, "w")) == NULL) {
                errorvf("Could not open penetrance file %s.\n", poutfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            /* write single affection locus */

            fprintf(filep, "%8sPROB    %8s        %2d\n",
                    strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN),
                    strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN),
                    2*Locus->Pheno->Props.Affection.ClassCnt);
            for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.ClassCnt; tmpi2++) {
                /* For unaffected status */
                TRAIT_PHEN1;
                fprintf(filep, "%7s ", trait_phen);

                for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) {
                    fprintf(filep, "%8.5f",
                            1 - Locus->Pheno->Props.Affection.Class[tmpi2].AutoPen[ipen]);
                }
                fprintf(filep, "\n");

                /* For affected status */
                TRAIT_PHEN2;
                fprintf(filep, "%7s ", trait_phen);

                for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) {
                    fprintf(filep, "%8.5f",
                            Locus->Pheno->Props.Affection.Class[tmpi2].AutoPen[ipen]);
                }
                fprintf(filep, "\n");
            }
        } else {
            if (first_time) {
                /* loop being executed the first time, open file
                   in current directory */
                sprintf(poutfl_name, "%s/%s", output_paths[0], fl_name);
                if ((filep = fopen(poutfl_name, "w")) == NULL) {
                    errorvf("Could not open penetrance file %s.\n", poutfl_name);
                    EXIT(FILE_WRITE_ERROR);
                }
                first_time=0;
            }
            /* if tr > 0, file is already open */
            /* print the current locus info */
            fprintf(filep, "%8sPROB    %8s        %2d\n",
                    strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN),
                    strtail(Locus->Name, MENDEL_MAX_LOCUS_NAME_LEN),
                    2*Locus->Pheno->Props.Affection.ClassCnt);
            for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.ClassCnt; tmpi2++) {
                /* For unaffected status */
                TRAIT_PHEN1;
                fprintf(filep, "%7s ", trait_phen);

                for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) {
                    fprintf(filep, "%8.5f",
                            1 - Locus->Pheno->Props.Affection.Class[tmpi2].AutoPen[ipen]);
                }
                fprintf(filep, "\n");

                /* For affected status */

                TRAIT_PHEN2;
                fprintf(filep, "%7s ", trait_phen);

                for (ipen = 0; ipen < Locus->Pheno->Props.Affection.PenCnt; ipen++) {
                    fprintf(filep, "%8.5f",
                            Locus->Pheno->Props.Affection.Class[tmpi2].AutoPen[ipen]);
                }
                fprintf(filep, "\n");
            }
        }
        if (LoopOverTrait == 1) {
            fclose(filep);
        }

    }
    if (LoopOverTrait == 0) {
        fclose(filep);
    }
    sprintf(err_msg, "        Penetrance file:      %s", fl_name);
    mssgf(err_msg);
    return 1;
}

/* routine to create the variable file for Mendel5
 */

static void write_mendel5_variable_file(char *varfl_name,
					linkage_locus_top *LTop)

{
    /* MOD TR */
    char vfl_name[2*FILENAME_LENGTH];
    int            i, *tr;
    FILE           *filep = NULL; //compiler: too hard to tell its assigned
    int            numq=0;
    int created=0;

    tr=&(global_trait_entries[0]);
    numq=0;

    for(i=0; i < num_traits; i++) {
        if (LTop->Locus[*tr].Type != QUANT) {
            continue;
        }
        created=1;
        if (LoopOverTrait == 0 || num_traits == 1) {
            if (i==0) {
                sprintf(vfl_name, "%s/%s", output_paths[0], varfl_name);
                filep=fopen(vfl_name, "w");
            }
        } else {
            sprintf(vfl_name, "%s/%s", output_paths[numq+1], varfl_name);
            filep=fopen(vfl_name, "w");
        }
        if (*tr == -1) {
            continue;
        }

        fprintf(filep, "%8s               \n",
                strtail(LTop->Locus[*tr].Name, MENDEL_MAX_LOCUS_NAME_LEN));
        tr++;
        if (LoopOverTrait == 1 && num_traits > 1) {
            fclose(filep);
            numq++;
        }

    }

    if (LoopOverTrait == 0 || num_traits == 1) {
        fclose(filep);
    }

    if (created) {
        sprintf(err_msg, "        Variable file:        %s",
                varfl_name);
        mssgf(err_msg);
    }
}


static void write_batch3(char *loutfl_name, char *outfl_name,
			 char *batfl_name, linkage_ped_top *LLTop)
{
    linkage_locus_top *LTop;
    int            *trp, tr, i, j, nloop;
    FILE           *filep;
    char           bat13fl_name[2*FILENAME_LENGTH];
    int            numql=0, num_affec=0;

    LTop = LLTop->LocusTop;

    for (tr=0; tr < num_traits; tr++) {
        if (global_trait_entries[tr] < 0) {
            continue;
        } else if (LTop->Locus[global_trait_entries[tr]].Type == QUANT) {
            numql++;
        } else if (LTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
            num_affec++;
        }
    }
    /* num_affec needs to be set before calling macro nloop */
    NLOOP;

    if (num_affec > 0) {
        trp=&(global_trait_entries[0]);
    } else {
        trp = NULL;
    }
    for (tr=0; tr <= num_affec; tr++) {
        if (tr==0 && nloop > 1) continue;

        if (trp != NULL && *trp == -1) {
            trp++;
        }
        sprintf(bat13fl_name, "%s/%s", output_paths[tr], batfl_name);

        if ((filep = fopen(bat13fl_name, "w")) == NULL) {
            errorvf("Could not open batch file %s.\n", bat13fl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        /* Set up a batch file for userm13 */
        fprintf(filep, "2\n");
        fprintf(filep, "%s\n",loutfl_name);  /* locus file name */
        fprintf(filep, "3\n");
        fprintf(filep, "%s\n",outfl_name);  /* pedigree file name */
        /* Number of quantitative variables is one if there is a trait
           locus; zero otherwise */
        fprintf(filep, "7\n");   /* Number of quantitative variables*/
        if (num_affec > 0) {
            if (LoopOverTrait == 1)	fprintf(filep, "%d\n", numql+1);
            else                      fprintf(filep, "%d\n", numql+num_affec);
        }
        else
            fprintf(filep, "%d\n", numql);   /* numql+0 if there is no trait */
        j=0;
        for (i = 0; i < NumChrLoci; i++) {
            if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
                LTop->Locus[ChrLoci[i]].Type == BINARY)  {
                if (j > 0) {
                    fprintf(filep, "yes\n");  /* yes if another run */
                }
                fprintf(filep, "9\n");   /* Option 9 */
                /* Have the loci names been truncated by now? */
                fprintf(filep, "%s\n", strtail(LTop->Locus[ChrLoci[i]].Name, MENDEL_MAX_LOCUS_NAME_LEN));  /* Name of marker locus */
                fprintf(filep, "\n");
                fprintf(filep, "17\n");  /* Set number of parameters */
                fprintf(filep, "%d\n",LTop->Locus[ChrLoci[i]].AlleleCnt); /* to the # of alleles */
                fprintf(filep, "21\n");
                j++;
            }
        }
        fprintf(filep, "no\n");  /* no more runs */
        fclose(filep);
        if (nloop == 1) break;
    }
    sprintf(err_msg, "        M13 batch file:       %s", batfl_name);
    mssgf(err_msg);
}

/*----------------------end of write_batch3-------------- */

static void write_mendel_map_using_sex_specific(char *mapfile, linkage_locus_top *LTop,
						analysis_type analysis)
{

    int i, j, nloop, tr, *trp = 0, num_affec=num_traits;
    char mapfl_name[2*FILENAME_LENGTH];
    double difff = INVALID_POS_DIFF, diffm = INVALID_POS_DIFF; //compiler: this makes diff f/m defined for truly odd cases

    FILE *filep;

    int num_markers, *markers;

    if (!HasMarkers) return;

    if (LoopOverTrait == 0) {
        num_markers = NumChrLoci + 1;
    } else {
        num_markers = NumChrLoci;
    }

    markers = CALLOC((size_t) num_markers, int);
    j=0;
    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY ||
            LTop->Locus[ChrLoci[i]].Type == XLINKED ||
            LTop->Locus[ChrLoci[i]].Type == YLINKED ||
            LoopOverTrait==0) {
            markers[j] = ChrLoci[i]; j++;
        }
    }
    num_markers=j;

    NLOOP;

    if (analysis != TO_HWETEST) {
        if (num_affec > 0) trp = &(global_trait_entries[0]);
    }
    /* NLOOP is always 1 for HWE_TEST*/

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        if (analysis != TO_HWETEST && num_affec > 0 && trp && *trp == -1) {
            trp++;
        }
        sprintf(mapfl_name, "%s/%s", output_paths[tr], mapfile);

        if (num_markers <= 1) {
            delete_file(mapfl_name); continue;
        }
        if ((filep = fopen(mapfl_name, "a")) == NULL) {
            errorvf("Could not open map file %s.\n", mapfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_affec > 0) {
            if (LoopOverTrait == 1 && trp && LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "%-8s\n", strtail(LTop->Locus[*trp].Name, MENDEL_MAX_LOCUS_NAME_LEN));
                fprintf(filep, "        0.50000 0.50000 ! %s\n",
                        (LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
            }
        }

        for (i = 1; i < num_markers; i++)  {
            switch(LTop->Locus[markers[i]].Type) {
            case NUMBERED:
            case BINARY:
	    case XLINKED:
	    case YLINKED:
                switch(LTop->Locus[markers[i-1]].Type) {
                case NUMBERED:
                case BINARY:
		case XLINKED:
		case YLINKED:
                    if (LTop->Marker[markers[i]].chromosome == LTop->Marker[markers[i-1]].chromosome) {
                        diffm = LTop->Marker[markers[i]].pos_male - LTop->Marker[markers[i-1]].pos_male;
                        difff = LTop->Marker[markers[i]].pos_female - LTop->Marker[markers[i-1]].pos_female;
                    } else {
                        diffm = difff = INVALID_POS_DIFF;
                    }
                    break;
                case AFFECTION :
                    if (LoopOverTrait == 0 && analysis != TO_HWETEST) {
                        difff = diffm = INVALID_POS_DIFF;
                    }
                    break;
                default:
                    break;
                }
                break;
            case AFFECTION :
                if (LoopOverTrait == 0 && analysis != TO_HWETEST) {
                    difff = diffm = INVALID_POS_DIFF;
                }
                break;
            default:
                break;
            }

            if (diffm <= INVALID_POS_DIFF) {
                diffm = difff = 0.5;
            } else {
                if (diffm < 0.0) {
                    sprintf(err_msg,
                            "%s (%7.4g) and %s (%7.4g) are not ordered by increasing male map distance!",
                            LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_male,
                            LTop->Marker[markers[i+1]].Name, LTop->Marker[markers[i+1]].pos_male);
                    warnf(err_msg);
                    warnf("Setting the distance between these to 0.001 cM.");
                    diffm=0.001;
                }
                if (difff < 0.0) {
                    sprintf(err_msg,
                            "%s (%7.4g) and %s (%7.4g) are not ordered by increasing female map distance!",
                            LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_female,
                            LTop->Marker[markers[i+1]].Name, LTop->Marker[markers[i+1]].pos_female);
                    warnf(err_msg);
                    warnf("Setting the distance between these to 0.001 cM.");
                    difff=0.001;
                }

                if (LTop->map_distance_type == 'k') {
                    diffm=kosambi_theta(diffm);
                    difff=kosambi_theta(difff);
                } else {
                    diffm=haldane_theta(diffm);
                    difff=haldane_theta(difff);
                }
            }

            fprintf(filep, "%-8s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL_MAX_LOCUS_NAME_LEN));
            fprintf(filep, "        %7.5f %7.5f ! %s\n",
                    difff, diffm,
                    (LTop->map_distance_type  == 'h' ? "Haldane" : "Kosambi"));
        }

        fprintf(filep, "%-8s\n",
                strtail(LTop->Locus[markers[num_markers-1]].Name, MENDEL_MAX_LOCUS_NAME_LEN));
        fclose(filep);
        if (nloop == 1) break;
        if (num_affec > 0) trp++;
    }
    if (analysis != TO_HWETEST) {
        sprintf(err_msg, "        Map file:             %s", mapfile);
        mssgf(err_msg);
    }
    free(markers);
}

static void write_mendel_map_using_sex_averaged(char *mapfile, linkage_locus_top *LTop,
						analysis_type analysis)

{

    int i, j, nloop, tr, *trp = 0, num_affec=num_traits;
    char mapfl_name[2*FILENAME_LENGTH];
    double diff = INVALID_POS_DIFF; //compiler: this makes diff defined for truly odd cases

    FILE *filep;

    int num_markers, *markers;

    if (!HasMarkers) return;

    if (LoopOverTrait == 0) {
        num_markers = NumChrLoci + 1;
    } else {
        num_markers = NumChrLoci;
    }

    markers = CALLOC((size_t) num_markers, int);
    j=0;
    for (i=0; i < NumChrLoci; i++) {
        if (LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY || LoopOverTrait==0) {
            markers[j] = ChrLoci[i]; j++;
        }
    }
    num_markers=j;

    NLOOP;

    if (analysis != TO_HWETEST) {
        if (num_affec > 0) trp = &(global_trait_entries[0]);
    }
    /* NLOOP is always 1 for HWE_TEST*/

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        if (analysis != TO_HWETEST && num_affec > 0 && trp && *trp == -1) {
            trp++;
        }
        sprintf(mapfl_name, "%s/%s", output_paths[tr], mapfile);

        if (num_markers <= 1) {
            delete_file(mapfl_name); continue;
        }
        if ((filep = fopen(mapfl_name, "a")) == NULL) {
	    errorvf("Could not open map file %s.\n", mapfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        if (num_affec > 0) {
            if (LoopOverTrait == 1 && trp && LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "%-8s\n", strtail(LTop->Locus[*trp].Name, MENDEL_MAX_LOCUS_NAME_LEN));
                fprintf(filep, "        0.50000\n");
            }
        }

        for (i = 1; i < num_markers; i++)  {
            switch(LTop->Locus[markers[i]].Type) {
            case NUMBERED:
            case BINARY:
                switch(LTop->Locus[markers[i-1]].Type) {
                case NUMBERED:
                case BINARY:
                    if (LTop->Marker[markers[i]].chromosome == LTop->Marker[markers[i-1]].chromosome) {
                        diff = LTop->Marker[markers[i]].pos_avg - LTop->Marker[markers[i-1]].pos_avg;
                    } else {
                        diff = INVALID_POS_DIFF;
                    }
                    break;
                case AFFECTION :
                    if (LoopOverTrait == 0 && analysis != TO_HWETEST) {
                        diff = INVALID_POS_DIFF;
                    }
                    break;
                default:
                    break;
                }
                break;
            case AFFECTION :
                if (LoopOverTrait == 0 && analysis != TO_HWETEST) {
                    diff = INVALID_POS_DIFF;
                }
                break;
            default:
                break;
            }

            if (diff <= INVALID_POS_DIFF) {
                diff = 0.5;
            } else {
                if (diff < 0.0) {
                    sprintf(err_msg,
                            "%s (%7.4g) and %s (%7.4g) are not ordered by increasing map distance!",
                            LTop->Marker[markers[i]].Name, LTop->Marker[markers[i]].pos_avg,
                            LTop->Marker[markers[i+1]].Name, LTop->Marker[markers[i+1]].pos_avg);
                    warnf(err_msg);
                    diff=0.001;
                    warnf("Setting the distance between these to 0.001 cM.");
                }

                if (LTop->map_distance_type == 'k') {
                    diff=kosambi_theta(diff);
                } else {
                    diff=haldane_theta(diff);
                }
            }

            fprintf(filep, "%-8s\n", strtail(LTop->Locus[markers[i-1]].Name, MENDEL_MAX_LOCUS_NAME_LEN));
            fprintf(filep, "        %7.5f\n", diff);
        }

        fprintf(filep, "%-8s\n", strtail(LTop->Locus[markers[num_markers-1]].Name, MENDEL_MAX_LOCUS_NAME_LEN));
        fclose(filep);
        if (nloop == 1) break;
        if (num_affec > 0) trp++;
    }
    if (analysis != TO_HWETEST) {
        sprintf(err_msg, "        Map file:             %s", mapfile);
        mssgf(err_msg);
    }
    free(markers);
}

void write_mendel_map(char *mapfl, linkage_locus_top *LTop, analysis_type analysis)
{
    if (genetic_distance_index > -1) {
        if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
  	    write_mendel_map_using_sex_averaged(mapfl, LTop, analysis);
        } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
            write_mendel_map_using_sex_specific(mapfl, LTop, analysis);
        } else {
            errorvf("Either a sex averaged or sex specific map must be provided to write the map file.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    } else {
        // write an error message if no genetic map was specified (currently one must be specified).
    }
}

/* Mendel control file,
   NM - Dec 8th, 2009, only write the first affection trait.
   breaks out of the trait look in Combine Mode after writing
   the first trait.
*/

static void write_mendel4_control(char *file_names[], int numchr,
				  linkage_locus_top *LTop,
				  analysis_type analysis,
				  int has_mult_liab)

{
    int nloop, tr, trp, num_affec=num_traits;
    char cfl[2*FILENAME_LENGTH], chr_str[3];
    int i;
    FILE *fp;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        if (LoopOverTrait == 1 && tr > 0) {
            trp = tr - 1;
        } else {
            trp = tr;
        }

        sprintf(cfl, "%s/%s", output_paths[tr], file_names[3]);
        if ((fp = fopen(cfl, "a")) == NULL) {
            errorvf("could not open control file %s.\n", file_names[3]);
            EXIT(FILE_WRITE_ERROR);
        }
        if (numchr > 0) {
            CHR_STR(numchr, chr_str);
            fprintf(fp, "OUTPUT_FILE = Mendel.%s.out\n", chr_str);
        } else {
            fprintf(fp, "OUTPUT_FILE = Mendel.All.out\n");
        }
        if (analysis == TO_MENDEL7_CSV) {
            fprintf(fp, "DEFINITION_FILE = %s\n", file_names[1]);
        } else {
            fprintf(fp, "LOCUS_FILE = %s\n", file_names[1]);
            if (HasQuant) {
                fprintf(fp, "VARIABLE_FILE=%s\n", file_names[7]);
            }
        }
        fprintf(fp, "MAP_FILE = %s\n", file_names[4]);
        fprintf(fp, "PEDIGREE_FILE = %s\n", file_names[0]);

        if (LoopOverTrait == 0 && HasAff && has_mult_liab) {
            fprintf(fp, "PENETRANCE_FILE = %s\n", file_names[6]);
        }

        fprintf(fp, "SUMMARY_FILE=Summary.%s.out\n", chr_str);
        fprintf(fp, "ANALYSIS_OPTION = Allele_frequencies\n");
        if (analysis == TO_MENDEL7_CSV) {
            fprintf(fp, "DEFAULT_LIST_READ = True\n");
        }

        /* write the penetrances */
        if (LoopOverTrait == 0) {
            for (i=0; i < num_traits; i++) {
                if (global_trait_entries[i] < 0 ||
		    LTop->Locus[global_trait_entries[i]].Type == QUANT) continue;

                fprintf(fp, "AFFECTED = 2\n");
                fprintf(fp, "AFFECTED_LOCUS_OR_FACTOR = %s\n",
                        ((analysis== TO_MENDEL7_CSV)?
                         strtail(LTop->Locus[global_trait_entries[i]].Name, MENDEL7_MAX_LOCUS_NAME_LEN) :
                         strtail(LTop->Locus[global_trait_entries[i]].Name, MENDEL_MAX_LOCUS_NAME_LEN)));

                if (LTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt == 1) {
                    fprintf(fp, "PENETRANCE = %5.4f :: 1/1\n",
                            LTop->Pheno[global_trait_entries[i]].Props.Affection.Class[0].AutoPen[0]);
                    fprintf(fp, "PENETRANCE = %5.4f :: 1/2\n",
                            LTop->Pheno[global_trait_entries[i]].Props.Affection.Class[0].AutoPen[1]);
                    fprintf(fp, "PENETRANCE = %5.4f :: 2/2\n",
                            LTop->Pheno[global_trait_entries[i]].Props.Affection.Class[0].AutoPen[2]);
                }
                break;
            }
        } else {
            if (LTop->Locus[global_trait_entries[trp]].Type == AFFECTION) {
                fprintf(fp, "AFFECTED = 2\n");
                fprintf(fp, "AFFECTED_LOCUS_OR_FACTOR = %s\n",
                        ((analysis== TO_MENDEL7_CSV)?
                         strtail(LTop->Locus[global_trait_entries[trp]].Name, MENDEL7_MAX_LOCUS_NAME_LEN) :
                         strtail(LTop->Locus[global_trait_entries[trp]].Name, MENDEL_MAX_LOCUS_NAME_LEN)));

                if (LTop->Pheno[global_trait_entries[trp]].Props.Affection.ClassCnt > 1) {
                    fprintf(fp, "PENETRANCE_FILE = %s\n", file_names[6]);
                } else {
                    fprintf(fp, "PENETRANCE = %5.4f :: 1/1\n",
                            LTop->Pheno[global_trait_entries[trp]].Props.Affection.Class[0].AutoPen[0]);
                    fprintf(fp, "PENETRANCE = %5.4f :: 1/2\n",
                            LTop->Pheno[global_trait_entries[trp]].Props.Affection.Class[0].AutoPen[1]);
                    fprintf(fp, "PENETRANCE = %5.4f :: 2/2\n",
                            LTop->Pheno[global_trait_entries[trp]].Props.Affection.Class[0].AutoPen[2]);
                }
            }
        }
        fclose(fp);
        if (nloop == 1) break;
    }
    sprintf(err_msg, "        Control file:         %s",
            file_names[3]);
    mssgf(err_msg);
    return;

}

/* remove Mendel4 format */

static void mendel_file_names(int glob_files,
			      char *file_names[], int *orig_id_opt,
			      int has_aff, int has_orig,
			      int has_uniq, analysis_type *analysis,
			      int *write_pen_file)

{
    int i=1, choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12];
    int var_it, map_it, pen_it, bat_it, m13_it, pedid_it, perid_it;
    int csv_it, ped_it, loc_it, csv_format;

    ped_it = loc_it = map_it = pen_it = m13_it =
        csv_it = bat_it = pedid_it = perid_it= var_it= -1;


    ((num_traits > 1 && LoopOverTrait)?
     strcpy(stem, "in each subdir") : strcpy(stem, ""));

    csv_format = 1;
    if (main_chromocnt > 1 && glob_files != 1) strcpy(stem, "stem");

    if (glob_files) {
        change_output_chr(file_names[0], 0);
        change_output_chr(file_names[1], 0);
        change_output_chr(file_names[3], 0);
        if (*analysis == TO_MENDEL4 || *analysis == TO_MENDEL7_CSV) {
            change_output_chr(file_names[4], 0);
/*       if (HasQuant) { */
/* 	change_output_chr(file_names[7], 0); */
/*       } */
        } else {
            change_output_chr(file_names[2], 0);
        }
        change_output_chr(file_names[6], 0);
    }

    if (!glob_files && main_chromocnt > 1) {
        change_output_chr(file_names[0], -9);
        change_output_chr(file_names[1], -9);
        change_output_chr(file_names[3], -9);
        if (*analysis == TO_MENDEL4 || *analysis == TO_MENDEL7_CSV) {
            change_output_chr(file_names[4], -9);
            if (HasQuant) {
                change_output_chr(file_names[7], -9);
            }
        } else {
            change_output_chr(file_names[2], -9);
        }
        change_output_chr(file_names[6], -9);
    }

    if (DEFAULT_OUTFILES) {
        /* log the id option */
        individual_id_item(0, TO_MENDEL, orig_id_opt[0], 0, 3, has_orig, has_uniq);
        pedigree_id_item(0, TO_MENDEL, orig_id_opt[1], 0, 3, has_orig);
        if (*analysis == TO_MENDEL7_CSV && csv_format == 0) {
            *analysis = TO_MENDEL4;
        }
        return;
    }

    /* added by Nandita - 7/28 get alternate output file names from the user*/
    while (choice != 0) {
        draw_line(); i=1;
        print_outfile_mssg();
        if (*analysis == TO_MENDEL7_CSV) {
            printf("  MENDEL 7+ file name menu\n");
        } else {
            printf("  MENDEL 3 file name menu\n");
        }

        draw_line();
        if (csv_format == 0) {
            warnf("Mendel 8 may not work on non CSV-format files.");
        }
        printf("0) Done with this menu - please proceed\n");
/*     if (*analysis == TO_MENDEL7_CSV) { */
/*       printf(" %d) Write CSV-format output files:     %s\n",  */
/* 	     i, yorn[csv_format]); csv_it=i; */
/*       i++; */

/*     } */

        if (main_chromocnt == 1 || glob_files) {
            printf(" %d) Locus file name:                   %-15s\t%s\n", i,
                   file_names[1],
                   file_status(file_names[1], fl_stat));
            loc_it = i; i++;
            printf(" %d) Pedigree filename:                 %-15s\t%s\n", i,
                   file_names[0],
                   file_status(file_names[0], fl_stat));
            ped_it = i; i++;
            if (*analysis == TO_MENDEL7_CSV) {
                printf(" %d) Map filename:                      %-15s\t%s\n",
                       i, file_names[4],
                       file_status(file_names[4], fl_stat));
                map_it=i; i++;
            }
/*       if (*analysis == TO_MENDEL7_CSV && HasQuant) { */
/* 	printf(" %d) Variable filename:                 %-15s\t%s\n", */
/* 	       i, file_names[7],  */
/* 	       file_status(file_names[7], fl_stat));  */
/* 	var_it=i; i++; */
/*       } */
            if (has_aff &&
                ((csv_format == 0) || (*write_pen_file == 1))) {
                printf(" %d) Penetrance file name:              %-15s\t%s\n", i,
                       file_names[6], file_status(file_names[6], fl_stat));
                pen_it = i; i++;
            }
        } else {
            printf(" %d) Locus file name %s:             %-15s\n", i, stem,
                   file_names[1]);
            loc_it = i; i++;
            printf(" %d) Pedigree file name %s:          %-15s\n", i, stem,
                   file_names[0]);
            ped_it = i; i++;
            if (*analysis == TO_MENDEL7_CSV) {
                printf(" %d) Map file name %s:               %-15s\n",
                       i, stem, file_names[4]);
                map_it = i; i++;
            }

            if (*analysis == TO_MENDEL7_CSV && !csv_format) {
                printf(" %d) Variable file name %s:          %-15s\n",
                       i, stem, file_names[7]);
                var_it = i; i++;
            }

            if ((has_aff > 0) &&
               ((csv_format == 0) || (*write_pen_file == 1))) {
                printf(" %d) Penetrance file name %s:        %-15s\n", i, stem,
                       file_names[6]);
                pen_it = i; i++;
            }

        }

        if (main_chromocnt == 1 || glob_files) {
            if (*analysis == TO_MENDEL && (num_traits >= 2 || LoopOverTrait == 1)) {
                printf(" %d) Batch file name:                   %-15s\t%s\n",
                       i, file_names[2],
                       file_status(file_names[2], fl_stat));
                bat_it = i; i++;
            }
            if (*analysis == TO_MENDEL) {
                printf(" %d) M13 batch file name:               %-15s\t%s\n", i,
                       file_names[3],
                       file_status(file_names[3], fl_stat));
            } else {
                printf(" %d) Control file name:                 %-15s\t%s\n", i,
                       file_names[3],
                       file_status(file_names[3], fl_stat));
            }
            m13_it = i; i++;
        } else {
            if (*analysis == TO_MENDEL && (num_traits >= 2 || LoopOverTrait == 1)) {
                printf(" %d) Batch file name %s:             %-15s\n", i, stem,
                       file_names[2]);
                bat_it = i; i++;
            }
            if (*analysis == TO_MENDEL) {
                printf(" %d) M13 batch file name %s:         %-15s\n", i, stem,
                       file_names[3]);
            } else {
                printf(" %d) Control file name %s:           %-15s\n",
                       i, stem,
                       file_names[3]);
            }
            m13_it = i; i++;
        }

        individual_id_item(i, TO_MENDEL, orig_id_opt[0], 43, 2, 0, 0);
        perid_it = i;  i++;
        pedigree_id_item(i, TO_MENDEL, orig_id_opt[1], 43, 2, 0);
        pedid_it = i;

        printf("Select options 0-%d to enter new file names/options > ", i);
        choice = -1;
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        default:
            if (csv_it > 0 && choice == csv_it) {
                csv_format=TOGGLE(csv_format);
            } else if (loc_it > 0 && choice == loc_it) {
                printf("New locus file name %s > ", stem);
                fcmap(stdin, "%s", file_names[1]);
                printf("\n");
            } else if (ped_it > 0 && choice == ped_it) {
                printf("New pedigree file name %s > ", stem);
                fcmap(stdin, "%s", file_names[0]);
                printf("\n");
            } else if (map_it > 0 && choice == map_it) {
                printf("New map file name %s > ", stem);
                fcmap(stdin, "%s", file_names[4]);
                printf("\n");
            } else if (var_it > 0 && choice == var_it) {
                printf("New variable file name %s > ", stem);
                fcmap(stdin, "%s", file_names[7]);
                printf("\n");
            } else if (pen_it > 0 && choice == pen_it) {
                printf("New penetrance file name %s > ", stem);
                fcmap(stdin, "%s", file_names[6]);
                printf("\n");
            } else if (bat_it > 0 && choice == bat_it) {
                printf("Name of batch file name %s > ", stem);
                fcmap(stdin, "%s", file_names[2]);
                printf("\n");
            } else if (m13_it > 0 && choice == m13_it) {
                (*analysis == TO_MENDEL)?
                    printf("Name of m13 batch file name %s > ", stem) :
                    printf("Name of control file name %s > ", stem);
                fcmap(stdin, "%s", file_names[3]);
                printf("\n");
            } else if (perid_it > 0 && choice == perid_it) {
                orig_id_opt[0] =
                    individual_id_item(0, TO_MENDEL, orig_id_opt[0], 35, 1,
                                       has_orig, has_uniq);
            } else if (pedid_it > 0 && choice == pedid_it) {
                orig_id_opt[1] =
                    pedigree_id_item(0, TO_MENDEL,
                                     orig_id_opt[1], 35, 1, has_orig);
            } else {
                warn_unknown(cchoice);
            }
            break;
        }
    }

    if (*analysis == TO_MENDEL7_CSV) {
        *analysis = ((*analysis == TO_MENDEL7_CSV && csv_format == 1)?
                     (analysis_type)TO_MENDEL7_CSV : (analysis_type)TO_MENDEL4);
    }

    /* log the id option */
    individual_id_item(0, TO_MENDEL, orig_id_opt[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, TO_MENDEL, orig_id_opt[1], 0, 3, has_orig);

    return;
}

/*------------------ End of mendel_file_names------------ */

static void  glob_mendel_files(int *glob_files)

{
    int select=-1;
    char cselect[10];

    if (main_chromocnt <= 1) {
        *glob_files=0;
        return;
    }

//SL
    if (batchANALYSIS && Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read) {
        *glob_files = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
        return;
    }

    *glob_files=0;
    if (DEFAULT_OPTIONS) {
        mssgf("Separate output files will be created for each chromosome (DEFAULT)");
        mssgf("as requested in the batch file.");
        return;
    }
    while(select != 0) {
        draw_line();
        printf("Mendel output options\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        printf(" 1) Separate files per chromosome\t[%s]\n",
               ((*glob_files==1)? "no" : "yes"));
        printf("Enter 1 to toggle, 0 to proceed > ");
        select=-1;
        fcmap(stdin, "%s", cselect); newline;
        sscanf(cselect, "%d", &select);
        switch(select) {
        case 0:
            break;
        case 1:
            *glob_files = ((*glob_files)? 0 : 1);
            break;
        default:
            warn_unknown(cselect);
            break;
        }
    }

    Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt = *glob_files == 0 ? 'y': 'n';

    if (InputMode == INTERACTIVE_INPUTMODE)
        batchf(Loop_Over_Chromosomes);

    return;

}
/*------------------ End of mendel_glob_files------------ */

static void simwalk2_file_names(char *file_names[], int *orig_id_opt,
				int *r_setup, int has_orig, int has_uniq)

{
    int i, choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12];
    int is_npl = *r_setup;
    /* r_setup is 1 if NPL analysis, so store this value to decide
       if this option should be displayed in the menu. From then on,
       the *r_setup value can be toggled to indicate if plots
       should be drawn, and is used in the rest of the analysis.
    */


    ((num_traits > 1 && LoopOverTrait)?
     strcpy(stem, "in each subdir"): strcpy(stem,""));

    if (main_chromocnt > 1) strcpy(stem, "stem");

    if (main_chromocnt > 1) {
        change_output_chr(file_names[0], -9);
        change_output_chr(file_names[1], -9);
        change_output_chr(file_names[6], -9);
        change_output_chr(file_names[3], -9);
        change_output_chr(file_names[2], -9);
        change_output_chr(file_names[4], -9);
    }

    if (DEFAULT_OUTFILES) {
        /* log the id option */
        if (orig_id_opt != NULL) {
            individual_id_item(0, TO_MENDEL, orig_id_opt[0],
                               0, 3, has_orig, has_uniq);
            pedigree_id_item(0, TO_MENDEL, orig_id_opt[1],
                             0, 3, has_orig);
        }
        return;
    }

    while (choice != 0) {
        draw_line(); i=1;
        print_outfile_mssg();
        printf("  SimWalk2 file name menu\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt == 1) {
            printf(" %d) Locus file name :                 %-15s\t%s\n", i, file_names[1],
                   file_status(file_names[1], fl_stat)); i++;
            printf(" %d) Pedigree filename :               %-15s\t%s\n", i, file_names[0],
                   file_status(file_names[0], fl_stat)); i++;
            printf(" %d) Map filename :                    %-15s\t%s\n", i, file_names[4],
                   file_status(file_names[0], fl_stat));
            if ((LoopOverTrait == 1 && num_traits > 0) || num_traits > 1) {
                i++;
                printf(" %d) Penetrance file name :            %-15s\t%s\n", i, file_names[6],
                       file_status(file_names[6], fl_stat));
            }
            i++;
            printf(" %d) Batch file name :                 %-15s\t%s\n", i, file_names[2],
                   file_status(file_names[2], fl_stat)); i++;
            printf(" %d) C-shell file name :               %-15s\t%s\n", i, file_names[3],
                   file_status(file_names[3], fl_stat));
        } else {
            printf(" %d) Locus file name %s :            %s\n", i, stem, file_names[1]);
            i++;
            printf(" %d) Pedigree filename %s :          %s\n", i, stem, file_names[0]);
            i++;
            printf(" %d) Map filename %s :               %s\n", i, stem, file_names[4]);
            if (num_traits > 0) {
                i++;
                printf(" %d) Penetrance file name %s :       %s\n", i, stem, file_names[6]);
            }
            i++;
            printf(" %d) Batch file name %s :            %s\n", i, stem, file_names[2]);
            i++;
            printf(" %d) C-shell file name %s :          %s\n", i, stem, file_names[3]);
        }

        if (orig_id_opt != NULL) {
            i++;
            individual_id_item(i, TO_MENDEL, orig_id_opt[0], 43, 2, 0, 0);
            i++;
            pedigree_id_item(i, TO_MENDEL, orig_id_opt[1], 43, 2, 0);
        }

        if (is_npl) {
            i++;
            printf(" %d) Generate R-plots for SimWalk2-NPL:       \t[%s]\n",
                   i, yorn[*r_setup]);
        }
        newline;
        printf("Select options 0-%d to enter new values",
               ((is_npl)? i - 1 : i));
        ((is_npl)?
         printf(",\n                              %d to toggle > ", i):
         printf(" > "));

        fcmap(stdin, "%s", cchoice); newline;
        choice = -1;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);
        if ((num_traits < 1 || (LoopOverTrait == 0 && num_traits < 2))
            && choice >= 4) choice++;
        switch (choice) {
        case 0:
            break;
        case 1:
            printf("New locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]);
            printf("\n");
            break;
        case 2:
            printf("New pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);
            printf("\n");
            break;
        case 3:
            printf("New map file name %s > ", stem);
            fcmap(stdin, "%s", file_names[4]);
            printf("\n");
            break;
        case 4:
            printf("New penetrance file name %s > ", stem);
            fcmap(stdin, "%s", file_names[6]);
            printf("\n");
            break;
        case 5:
            printf("Name of Batch file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]);
            printf("\n");
            break;
        case 6:
            printf("Name of C-Shell file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]);
            printf("\n");
            break;
        case 7:
            if (orig_id_opt != NULL) {
                orig_id_opt[0] =
                    individual_id_item(0, TO_MENDEL, orig_id_opt[0], 35, 1,
                                       has_orig, has_uniq);
                break;
            }
        case 8:
            if (orig_id_opt != NULL) {
                orig_id_opt[1] = pedigree_id_item(0, TO_MENDEL,
                                                  orig_id_opt[1], 35, 1,
                                                  has_orig);
                break;
            }
        case 9:
            if (is_npl) {
                *r_setup = TOGGLE(*r_setup);
                break;
            }

        default:
            warn_unknown(cchoice);
            break;
        }
    }
    if (orig_id_opt != NULL) {
        individual_id_item(0, TO_MENDEL, orig_id_opt[0], 0, 3, has_orig, has_uniq);
        pedigree_id_item(0, TO_MENDEL, orig_id_opt[1], 0, 3, has_orig);
    }
    return;
}

/*------------------ End of simwalk2_file_names------------ */

/* Note on having multiple traits and looping over each trait:
   The output_paths will be set up such that paths[0] contains ".",
   the current directory, and the sub-directory names are stored
   in paths[1 .. num_affec].
   Thus, if we are looping over multiple traits, we only execute
   the loop for i=1 to i=num_affec, otherwise we execute the loop
   once for i=0. Penetrance file is created if num_affec is at least
   one.
*/

static void init_mendel_files(char *file_names[], int glob_files,
			      analysis_type analysis,
			      linkage_locus_top *LTop,
			      int write_pen_file)

{
    int i, j, *trp, nloop, num_affec = num_traits;
    char fl_name[2*FILENAME_LENGTH];
    FILE *fp;

    NLOOP;

    if ((LoopOverTrait == 0 && num_affec > 1) ||  (LoopOverTrait == 1)) {
        trp = &(global_trait_entries[0]);
    } else {
        trp = NULL;
    }
    for(i=0; i<= nloop; i++) {
        if (nloop > 1 && i==0) continue;

        if (!glob_files  && main_chromocnt > 1) {
            for (j=0; j < main_chromocnt; j++) {
                /* pedigree file */
                if ((SIMWALK2(analysis)) &&
                    global_chromo_entries[j] == SEX_CHROMOSOME
                    && LTop->SexLinked) {
                    continue;
                }
                change_output_chr(file_names[0], global_chromo_entries[j]);
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[0]);
                init_file(fl_name);
                /* locus file */
                change_output_chr(file_names[1], global_chromo_entries[j]);
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[1]);
                init_file(fl_name);

/* 	if (analysis == TO_MENDEL4 || analysis == TO_MENDEL7_CSV) { */
                if (analysis != TO_MENDEL) {
                    change_output_chr(file_names[4], global_chromo_entries[j]);
                    sprintf(fl_name, "%s/%s", output_paths[i], file_names[4]);
                    init_file(fl_name);
                }

                if (trp != NULL) {
                    /* pen file */
                    change_output_chr(file_names[6], global_chromo_entries[j]);
                    sprintf(fl_name, "%s/%s", output_paths[i], file_names[6]);
                    if (*trp == -1) {
                        trp++;
                    }

                    if (LoopOverTrait == 0) {
                         init_file(fl_name);
                    } else if (LTop->Locus[*trp].Type == AFFECTION) {
                        if (analysis != TO_MENDEL7_CSV ||  LTop->Pheno[*trp].Props.Affection.ClassCnt > 1) {
                            init_file(fl_name);
                        } 
                    }
                }
                if (SIMWALK2(analysis) || (analysis == TO_MENDEL && num_affec > 0)) {
                    /* simwalk2 batch file or Mendel m15 batch file*/
                    change_output_chr(file_names[2], global_chromo_entries[j]);
                    sprintf(fl_name, "%s/%s", output_paths[i], file_names[2]);
                    init_file(fl_name);
                }
                /* m13 batch file or simwalk2 shell file*/
                change_output_chr(file_names[3], global_chromo_entries[j]);
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[3]);
                init_file(fl_name);
            }
        } else {
            /* pedigree file */
            sprintf(fl_name, "%s/%s", output_paths[i], file_names[0]);
            fp=fopen(fl_name, "w"); fclose(fp);
            /* locus file */
            sprintf(fl_name, "%s/%s", output_paths[i], file_names[1]);
            init_file(fl_name);
            /* map file */
/*       if (analysis == TO_MENDEL4 || analysis == TO_MENDEL7_CSV) { */
            if (analysis != TO_MENDEL) {
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[4]);
                init_file(fl_name);
            }
            if (trp != NULL) {
                /* pen file */
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[6]);
                if (*trp == -1) {
                    trp++;
                }
                if (LoopOverTrait == 0) {
                    if (write_pen_file == 1) {
                        init_file(fl_name);
                    } 
                } else if (LTop->Locus[*trp].Type == AFFECTION) {
                    if (analysis != TO_MENDEL7_CSV ||  LTop->Pheno[*trp].Props.Affection.ClassCnt > 1) {
                        init_file(fl_name);
                    } 
                }
            }

            if (SIMWALK2(analysis) || (analysis == TO_MENDEL && HasAff > 0)) {
                /* simwalk2 BATCH2 or Mendel m15*/
                sprintf(fl_name, "%s/%s", output_paths[i], file_names[2]);
                init_file(fl_name);
            }
            /* m13 batch file or simwalk2 shell script */
            sprintf(fl_name, "%s/%s", output_paths[i], file_names[3]);
            init_file(fl_name);
        }
        if (nloop==1) break;
        trp++;
    }
    return;

}

static  int count_aff_classes(linkage_locus_top *LTop)

{
    int i, write_pen_file = 0;

    for (i=0; i < num_traits; i++) {
        if (global_trait_entries[i] < 0) continue;
        if (LTop->Locus[global_trait_entries[i]].Type == AFFECTION &&
            LTop->Pheno[global_trait_entries[i]].Props.Affection.ClassCnt > 1) {
            write_pen_file = 1;
            break;
        }
    }

    return write_pen_file;
}
/*------------------ End of init_mendel_files------------ */

/* #ifdef LPedTreeTop  */
/* #undef LPedTreeTop */
/* #endif */
/* #define LPedTreeTop (*Top) */

int  create_mendel_file(linkage_ped_top **Top, char *mapfl_name,
			ped_top *PedTreeTop,
			file_format *infl_type,
			file_format *outfl_type, int *numchr,
			analysis_type *analysis,
			char *file_names[], int untyped_ped_opt)
{
    linkage_ped_tree *LPed;
    int             glob_files, sex_linked;
    int             i, j, first_time, simwalk2=0;
    int             R_setup, global=0;
    char            syscmd[2*FILENAME_LENGTH + 10];
    int             qtdt, write_pen_file=0; /* checke donly for csv */
    linkage_ped_top *LPedTreeTop = *Top;

/*    qtdt = ((LPedTreeTop->LocusTop->PedRecDataType == Premakeped || */
/*  	   LPedTreeTop->LocusTop->PedRecDataType == Postmakeped)?  */
/*  	  0 : 1); */
    qtdt=1;
    R_setup = (*analysis == NONPARAMETRIC || *analysis == TO_MERLIN);

    /* assuming only 3 status classes, 0, 1, 2 */
    /*    if (*analysis == HAPLOTYPE || *analysis == LOCATION || */
    /*        *analysis == NONPARAMETRIC || *analysis == IBD_EST || */
    /*        *analysis == MISTYPING) */

    if (SIMWALK2(*analysis)) simwalk2=1;

    if (simwalk2) {
        glob_files =0;
        if (*analysis == TO_MERLIN) {
            simwalk2_file_names(file_names, NULL, &R_setup,
                                LPedTreeTop->OrigIds, LPedTreeTop->UniqueIds);
        } else {
            simwalk2_file_names(file_names, &(OrigIds[0]), &R_setup,
                                LPedTreeTop->OrigIds, LPedTreeTop->UniqueIds);
        }
    } else {
        glob_mendel_files(&glob_files);
//SL
        if (LPedTreeTop->LocusTop->SexLinked == 2 && glob_files == 1) {
            if (batchANALYSIS)
                if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
                    errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
                else
                    errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

            EXIT(DATA_INCONSISTENCY);
        }
        write_pen_file = count_aff_classes(LPedTreeTop->LocusTop);
        mendel_file_names(glob_files, file_names, &(OrigIds[0]), HasAff,
                          LPedTreeTop->OrigIds, LPedTreeTop->UniqueIds,
                          analysis, &write_pen_file);
    }

    if (R_setup) {
        R_setup = sw2_R_setup(*analysis, *numchr, *Top);
    }

    if ((LoopOverTrait == 0 && num_traits < 2) || num_traits == 0) {
        if (*analysis != LOCATION && *analysis != NONPARAMETRIC) {
            printf("WARNING: there are NO affection status loci present.\n");
        } else {
            /* LOCATION, NONPARAMETRIC */
            errorvf("INTERNAL: NO affection status loci present.\n");
            EXIT(SYSTEM_ERROR);
        }
    }

    if (simwalk2 && main_chromocnt > 1) {
        /* just set it with a chromosome number*/
        change_output_chr(file_names[3], global_chromo_entries[0]);
        strcat(file_names[3], ".sh");
    }

    init_mendel_files(file_names, glob_files, *analysis, LPedTreeTop->LocusTop, write_pen_file);

    if (*infl_type == LINKAGE && *analysis != TO_MENDEL) {
        for (i = 0; i < LPedTreeTop->PedCnt; i++)  {
            LPed = &(LPedTreeTop->Ped[i]);
            if (LPed->Loops != NULL)   {
                sprintf(err_msg, "Connecting loops in pedigree %d...", LPed->Num);
                mssgf(err_msg);
                if (connect_loops(LPed, LPedTreeTop) < 0)    {
                    errorf("Fatal error! Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }


    first_time=1;

    if ((SIMWALK2(*analysis)) &&
       ((main_chromocnt > 1) || (num_traits > 2) ||
        (LoopOverTrait == 1 && num_traits > 1))) {
        FILE *filep;
        char fl[2*FILENAME_LENGTH];

        change_output_chr(file_names[3], 0);
        if (!strcmp(output_paths[0], ".")) {
            strcpy(fl, file_names[3]);
        } else {
            sprintf(fl, "%s/%s", output_paths[0], file_names[3]);
        }
        if ((filep = fopen(fl, "w")) == NULL) {
            errorvf("Could not open %s for writing\n", fl);
            EXIT(FILE_WRITE_ERROR);
        }
        script_time_stamp(filep);
        fclose(filep);
        if (main_chromocnt == 1) {
            change_output_chr(file_names[3], *numchr);
        }

    }

    if (main_chromocnt > 1 || (num_traits > 1 && LoopOverTrait == 1)) {
        global=1;
    } else {
        global=0;
    }
    LPedTreeTop = *Top;

    for (j=0; j < main_chromocnt; j++) {
        if (!glob_files && main_chromocnt > 1) {
            *numchr = global_chromo_entries[j];
            get_loci_on_chromosome(*numchr);
            change_output_chr(file_names[0], *numchr); /* peds */
            change_output_chr(file_names[1], *numchr); /* locus */
            change_output_chr(file_names[6], *numchr); /* penetrance */
            change_output_chr(file_names[2], *numchr); /* batch.chr or batch2.chr */
            change_output_chr(file_names[3], *numchr); /* batchm13, cshell file
                                                          or control file */
            if (*analysis != TO_MENDEL) {
                change_output_chr(file_names[4], *numchr); /* mendel4 map*/
            }
        }

        if (glob_files) {
            get_loci_on_chromosome(0);
            sex_linked = LPedTreeTop->LocusTop->SexLinked;
        } else {
            get_loci_on_chromosome(*numchr);
            sex_linked =
                (((*numchr == SEX_CHROMOSOME &&
                   LPedTreeTop->LocusTop->SexLinked == 2) ||
                  (LPedTreeTop->LocusTop->SexLinked == 1))? 1 : 0);
        }

        if (sex_linked) {
            if (SIMWALK2(*analysis)) {
                warnf("Files not created for X-Chromosome.");
                continue;
            }
        }

        if (*analysis == TO_MENDEL || *analysis == TO_MENDEL4 ||
            *analysis == TO_MENDEL7_CSV ) {
            put_quant_locus_last(LPedTreeTop);
        }

        omit_peds(untyped_ped_opt, LPedTreeTop);

        if (j==0) create_mssg(*analysis);
        /* 1. locus and pedigree files */
        if (*analysis == TO_MENDEL7_CSV) {
            csv_write_mendel_locus_file(file_names[1],
                                        LPedTreeTop->LocusTop,
                                        sex_linked);
            csv_save_mendel_peds(file_names[0], LPedTreeTop);
        } else {
            write_mendel_locus_file(file_names[1], LPedTreeTop->LocusTop,
                                    *analysis,
                                    sex_linked);
            save_mendel_peds(file_names[0], LPedTreeTop, *analysis);
        }

        /* 2. Map file */
        if (*analysis == TO_MENDEL4 || SIMWALK2(*analysis)) {
            write_mendel_map(file_names[4], LPedTreeTop->LocusTop,
                             *analysis);
        } else if (*analysis == TO_MENDEL7_CSV) {
            csv_write_mendel_map(file_names[4], LPedTreeTop->LocusTop,
                                 *analysis, sex_linked);

        }
        /* 3. Penetrance and variable files for MENDEL5 - non-csv */
        if (*analysis == TO_MENDEL4) {
            if (HasAff) {
                write_mendel5_pen_file(file_names[6], LPedTreeTop->LocusTop);
            }
            if (HasQuant) {
                write_mendel5_variable_file(file_names[7],
                                            LPedTreeTop->LocusTop);
            }
        }

        /* 4. Penetrance file */
        if (*analysis == TO_MENDEL7_CSV) {
            if (write_pen_file) {
                csv_mendel7_pen_file(file_names[6], LPedTreeTop, sex_linked);
            }
        } else {
            if (HasAff) {
                write_mendel_pen_file(file_names[6], LPedTreeTop->LocusTop,
                                      *analysis);
            }
        }

        if (*analysis == TO_MENDEL7_CSV) {
            if (!glob_files || main_chromocnt > 1) {
                write_mendel4_control(file_names, *numchr,
                                      LPedTreeTop->LocusTop, *analysis,
                                      write_pen_file);
            } else {
                write_mendel4_control(file_names, 0,
                                      LPedTreeTop->LocusTop, *analysis,
                                      write_pen_file);
            }
        }

        /* 4. Control files */
        if (simwalk2) {
            write_batch2(file_names[2], *numchr, LPedTreeTop, *analysis);
            write_cshell(*numchr, *analysis, file_names, global,
                         sex_linked, qtdt);
        }

        if (*analysis == TO_MENDEL) {
            if (HasAff) {
                write_batch(file_names[1], file_names[0], file_names[2], LPedTreeTop,
                            first_time);
            }
            write_batch3(file_names[1], file_names[0], file_names[3],
                         LPedTreeTop);
            first_time=0;
        }

        log_line(mssgf);
/*     if (main_chromocnt > 1 && !glob_files) { */
/*       /\*      free_all_including_lpedtop(&LPedTreeTop); *\/ */
/*       free_chr_loci(); */
/*     } */
        if (glob_files) break;
    }

    if (simwalk2 && (global || R_setup)) {
        mssgf("Additional files:");
    }
    if (R_setup && LPedTreeTop->LocusTop->SexLinked != 1) {
        if (global) {
            change_output_chr(file_names[3], 0);
        }
        append_R_commands(file_names[3], Rshell_file);
        sprintf(err_msg, "        R-script file:        %s", Rscript_file);
        mssgf(err_msg);
        sprintf(err_msg, "        Shell file to run R:  %s", Rshell_file);
        mssgf(err_msg);
        if (base_pair_position_index >= 0) {
            // Needed to produce UCSC Custom Tracks which is only done if a physical map file exists.
            sprintf(err_msg, "        Physical map for R:   %s", Rmap_file);
            mssgf(err_msg);
        }
    }

    if (global && simwalk2) {
        change_output_chr(file_names[3], 0);
        if (strcmp(output_paths[0], ".")) {
            sprintf(syscmd, "chmod +x %s/%s\n", output_paths[0], file_names[3]);
        } else {
            sprintf(syscmd, "chmod +x %s\n", file_names[3]);
        }
        System(syscmd);
        sprintf(err_msg, "        Combined shell script:       %s", file_names[3]);
        mssgf(err_msg);
/*     if (main_chromocnt > 1) { */
/*       mssgf("        for all chromosomes "); */
/*     } */
/*     if (num_traits > 1 && LoopOverTrait == 1) { */
/*       mssgf("        for all traits "); */
/*     } */
    }

    if (R_setup && LPedTreeTop->LocusTop->SexLinked != 1) {
        log_line(mssgf);
        sprintf(err_msg,
                "Run %s in %s for all analysis and creation of R-plots",
                file_names[3],
#ifdef HIDEPATH
                NOPATH
#else
                ((strcmp(output_paths[0], ".")?
                  output_paths[0] : "current directory"))
#endif
            );
        mssgf(err_msg);
    }

    if (liability_multiplier > 100) {
        if (simwalk2) {
            warnf("*********");
            warnf("Found trait loci with more than 99 liability classes.");
            warnf("SimWalk2 will not accept the penetrance files for these loci.");
            warnf("*********");
        }
    }
    if (simwalk2 && LPedTreeTop->LocusTop->SexLinked==1) {
        warnf("*********");
        warnf("Input data contains x-linked markers.");
        warnf("SimWalk2 does not support x-linked analysis.");
        warnf("*********");
    }

/*   if (renamed_peds) { */
/*     sprintf(err_msg, "Pedigree names shortened to %d characters,", */
/* 	    ((analysis == TO_MENDEL_CSV)? 16: 7)); */
/*     warnf(err_msg); */
/*     warnf("please check to see that pedigree names are unique."); */
/*   } */

/*   if (renamed_loci) {  */
/*     sprintf(err_msg, "Marker names shortened to %d characters,", */
/* 	    ((analysis == TO_MENDEL_CSV)? 16: 7)); */
/*     warnf(err_msg); */
/*   } */

    return 1;
}

#undef LPedTreeTop

/*--------------- End of create_mendel_file-----------------*/



int get_aff_status(linkage_locus_rec *Locus, analysis_type analysis,
		   linkage_ped_rec *Entry, int locus1)
{

    if (analysis == LOCATION || analysis == HAPLOTYPE ||
        analysis == NONPARAMETRIC || analysis == IBD_EST ||
        analysis == MISTYPING || analysis == TO_MERLIN) {

        if (Entry->Pheno[locus1].Affection.Status == 0)  return 0;
        else if (Locus->Pheno->Props.Affection.ClassCnt == 1)  return Entry->Pheno[locus1].Affection.Status;
        else {
            return(liability_multiplier*Entry->Pheno[locus1].Affection.Status +
                   Entry->Pheno[locus1].Affection.Class);
        }
    } else if (analysis == TO_MENDEL || analysis == TO_MENDEL4 ||
               analysis == TO_MERLINONLY) {

        if (Entry->Pheno[locus1].Affection.Status == 0) return 0;
        else {
            if (Locus->Pheno->Props.Affection.ClassCnt > 1)
                return(liability_multiplier*Entry->Pheno[locus1].Affection.Status +
                       Entry->Pheno[locus1].Affection.Class);
            else
                return(Entry->Pheno[locus1].Affection.Status);
        }
    }
    return 0;
}
