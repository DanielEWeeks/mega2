/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"

#include "loop.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

/*
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_loci_on_chromosomes get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp summary_time_stamp
*/


static void inner_file_names(char **file_names, const char *num, const char *stem /* = "structure" */)
{
    sprintf(file_names[0], "%s.%s.data", stem, num);
    sprintf(file_names[3], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.all.sh", stem);
    sprintf(file_names[5], "%s.%s.mainparams", stem, num);
    sprintf(file_names[6], "%s.%s.results", stem, num);
    sprintf(file_names[7], "%s.%s.extraparams", stem, num);
    sprintf(file_names[8], "%s.all.data", stem);
}

void CLASS_STRUCTURE::file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num, "structure");
}

void CLASS_STRUCTURE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[6], numchr);
}

/**
 @brief Structure

Home Page...
http://pritch.bsd.uchicago.edu/structure.html

2.3.4 documentation...
http://pritch.bsd.uchicago.edu/structure_software/release_versions/v2.3.4/structure_doc.pdf

Structure is used for inferring population structure using genotype data. Mega2 will produce files
that are used by the "computational part" (sometimes referred to the structure kernel) of the
program (stricture) that is written in C.

The structure kernel reads two "control" files and one data file. Examples of both of these control files
and a data file are included as a part of the structure download. One of the control files "extraparams"
contains parameters that describes algorithm variable settings (e.g. how the program runs) and so is not
of concern for our purposes in data file conversion. The "mainparams" control file contins parameters
that describe the data file used as input to stricture, and so is of concern here. Mega2 will produce a
data file which stores data for each individual on one line (rather than the default two lines per individual).
Mega2 allows you to choose to associate a quantitative phenotype to a population identifier through the
use of the batch file parameter Structure.PopDataPheno.

Section (2.1 "Components of the data file") specifies the order of rows listed in the data file.
As a point of confirmation, the data file is read in datain.c:ReadData() which processes lines in
the data file in a specified order, according to settings in the "mainparams" file.
The first line to be output is the marker names (mainparams: #define MARKERNAMES 1).

For our purpose "Inter-Marker Distances" is the first row that will appear in the data file is that
of the inter-marker distances. These will be genetic distances in cM. Markers from different chromosomes
are separated by the value of '-1'.

The variables in the "mainparams" control file that Mega2 will modify are:
NUMINDS the number of individuals in the structure input data file produced by Mega2
NUMLOCI the number of loci in the data file
ONEROWPERIND the data for each individual will be stored in a single line (= 1) where each locus is
in two consecutive columns (diploid organizms).
LABEL the file will contain individual labels (= 1)
POPDATA the user may chooses to associate a quantitative phenotype to a population identifier (= 1)
PHENOTYPE is set (= 1) if the user specifies a phenotype.
MAPDISTANCES The user may choose to include map distnaces in the output (required in order to run
the linkage model. see extraparams)
 */

static int mapdistances; // flag telling if MAPDISTANCES should be set in 'mainparams'
static int popdata; // flag telling if POPDATA should be set in 'mainparams'
static int phenotype; // flag telling if (the optional) phenotype is present
static int *markers_per_chromo; // number of markers encountered per file or chromosome
static int *numinds; // number of individuals who has at least one marker selected per file or chromosome

static void write_INFILE(linkage_ped_top *Top, char *file_names[],
                         const int pwid, const int fwid)
{
    // #define MARKERNAMES      1  // (B) data file contains row of marker names
    struct save_marker_names: public loop::outer, loop::loci {
        typedef char *str;
        str *file_names;
        int markers_i;
        
        save_marker_names(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::loci(Top) { markers_i = -1; }
        void make_file() {
            run_loop(file_names[0]);
        }
        void chr_start() {
            markers_i++;
        }
        void inner() {
            pr_printf("%15s ", _tle->Name);
            markers_per_chromo[markers_i]++;
        }
        void file_trailer() {
            pr_nl();
        }
    } *smn = new save_marker_names(Top);
    
    smn->file_names = file_names;
    smn->load_formats(fwid, pwid, -1);
    smn->iterate();
    delete smn;

    mapdistances = 0;
    if (genetic_distance_index >= 0) {
        // The saving of 'LoopOverChrm' is a feature of the looping code. It was never ment to handle
        // looping over chromosomes when the files did not change. Here we fake it into appending into
        // the same file because we need this behavior.
        int LoopOverChrm_save = LoopOverChrm;
        LoopOverChrm = 1;

        mapdistances = 1;
        // #define MAPDISTANCES     1  // (B) data file contains row of map distances between loci
        struct save_genetic_distance_markers: public loop::outer, loop::loci {
            typedef char *str;
            str *file_names;
            double last_genetic_distance;
            int new_chromo, warnp, LoopOverChrm_save;
            
            save_genetic_distance_markers(linkage_ped_top *Top, int LoopOverChrm_save) : person_locus_entry(Top), loop::outer(Top), loop::loci(Top) { this->LoopOverChrm_save = LoopOverChrm_save; }
            void make_file() {
                char *fn = file_names[(LoopOverChrm_save == 1 ? 0 : 8)];
                run_loop(*_opath, fn, "a");
            }
            void chr_start() {
                new_chromo = 1;
            }
            void inner() {
                // The markers must be in map order within linkage groups. When consecutive markers are
                // from different linkage groups (e.g., different chromosomes), this should be indicated
                // by the value -1. The first marker is also assigned the value -1. All other distances
                // are non-negative. This row contains L real numbers.
                if (new_chromo == 1) {
                    pr_printf("-1 ");
                    new_chromo = 0;
                    last_genetic_distance = get_genetic_distance(&warnp);
                } else {
                    double genetic_distance = get_genetic_distance(&warnp);
                    double diff = genetic_distance - last_genetic_distance;
#if defined(_WIN) || defined(MINGW)
                    pr_printf("%.6f ", diff);
#else
                    pr_printf("%.6lf ", diff);
#endif
                    last_genetic_distance = genetic_distance;
                }
            }
            void file_trailer() {
                // Print the new line only if these are different files...
                if (LoopOverChrm_save == 1) pr_nl();
            }
        } *sgdM = new save_genetic_distance_markers(Top, LoopOverChrm_save);
        
        sgdM->file_names = file_names;
        sgdM->load_formats(fwid, pwid, -1);
        sgdM->iterate();
        delete sgdM;

        if (LoopOverChrm_save == 0) {
            char fn[2*FILENAME_LENGTH];
            FILE *filep;
            sprintf(fn, "%s/%s", output_paths[0], file_names[8]);
            filep = fopen(fn, "a");
            fprintf(filep, "\n");
            fclose(filep);
        }

        LoopOverChrm = LoopOverChrm_save;
    }

    phenotype = popdata = 0;
    // 2.3 Individual/genotype data
    // Each row of individual data contains the following elements. These form columns in the data file.
    struct save_pers: public loop::outer, loop::ped_per_loci {
        typedef char *str;
        str *file_names;
        int PopDataPheno_i, personHasMarkers, chr_i;
        
        save_pers(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::ped_per_loci(Top) { chr_i = -1; }
        void make_file() {
            int i;
            // The user can specify an integer (quantitative phenotype) designating a user-defined
            // population from which the individual was obtained.
            PopDataPheno_i = -1;
            if (Mega2BatchItems[Structure$PopDataPheno].item_read) {
                char *PopDataPheno = Mega2BatchItems[Structure$PopDataPheno].value.name;
                for (i=0; i<_LTop->LocusCnt; i++) {
                    if (_LTop->Locus[i].Type == QUANT && strcasecmp(PopDataPheno, _LTop->Locus[i].Name) == 0) {
                        PopDataPheno_i = i;
                        popdata = 1;
                    }
                }
                if (PopDataPheno_i == -1) {
                    warnvf("Batch file item 'Structure.PopDataPheno' was specified as '%s',\n", PopDataPheno);
                    warnf("but no such quantitative phenotype exists.");
                }
            }

            if (_tte != (linkage_locus_rec *)NULL) phenotype = 1;

            numinds[++chr_i] = _Top->IndivCnt; // initially consider all of the individuals

            run_loop(*_opath, file_names[0], "a");
        }
        void per_start() {
            int m;
            personHasMarkers = 0;
            // Since Structure does not do anything useful with completely untyped people (and does not
            // take family structure into account), when generating the Structure data file, please
            // exclude completely untyped people who have no known genotypes at the selected markers at all.
            for (m=0; m < NumChrLoci; m++) {
                int _locus = ChrLoci[m];
                linkage_locus_rec *_tle = &(_LTop->Locus[_locus]);
                if (_tle->Class == MARKER) {
                    int _allele1, _allele2;
                    get_2alleles(_tpe->Marker, _locus, &_allele1, &_allele2);
                    if (_allele1 || _allele2) { personHasMarkers = 1; break; }
                }
            }
            // exceude individuals who do not have a marker on any of the loci selected...
            if (personHasMarkers == 0) { numinds[chr_i]--; return; }
            
            // #define LABEL     1     // (B) Input file contains individual labels
            // Label (Optional; string) A string of integers or characters used to designate each
            // individual in the sample.
            pr_printf("%s ", _tpe->UniqueID);
            
            // #define POPDATA   1     // (B) Input file contains a population identifier
            // Only the integer portion is valid here. We already know that it's a QUANT from 'make_file()'
            // above so there is no need to check here, just use it...
            if (PopDataPheno_i != -1) {
                double popData = _tpe->Pheno[PopDataPheno_i].Quant;
                if (popData < 0) {
                    errorvf("Processing a user-defined population-of-origin data (POPDATA) value\n");
                    errorvf("using quantitative phenotype information '%s' specified by the\n",
                            Mega2BatchItems[Structure$PopDataPheno].value.name);
                    errorvf("Batch Item 'Structure.PopDataPheno'.\n");
                    if (popData == QMISSING)
                        errorvf("Found a missing value indicator for individual '%s'.\n",
                                _tpe->UniqueID);
                    else
                        errorvf("Found a negative value of '%lf' for individual '%s'.\n",
                                popData, _tpe->UniqueID);
                    EXIT(DATA_TYPE_ERROR);
                }
                pr_printf("%10.0f ", popData);
            }
            
            // #define PHENOTYPE 1     // (B) Input file contains phenotype information
            // Phenotype (Optional; integer) An integer designating the value of a phenotype of
            // interest, for each individual.
            if (phenotype == 1) pr_pheno();
        }
        void inner() {
            if (personHasMarkers == 0) return;
            // Genotype Data (Required; integer) Each allele at a given locus should be coded by a unique integer
            // The integer printed here is associated with the allele name printed by 'pr_marker_alleles()'
            pr_marker();
        }
        void per_end() {
            if (personHasMarkers == 0) return;
            pr_nl();
        }
        void file_trailer() {
            mssgvf("        STRUCTURE Data File:          %s/%s\n", *_opath, file_names[0]);
            if (numinds[chr_i] != _Top->IndivCnt) {
                warnvf("%d People were filtered out because they were not genotyped at the selected markers.\n",
                       _Top->IndivCnt - numinds[chr_i]);
            }
        }
    } *sI = new save_pers(Top);
    
    sI->file_names = file_names;
    sI->load_formats(fwid, pwid, -1);
    sI->iterate();
    delete sI;
}

static void write_mainparams(linkage_ped_top *Top, char *file_names[],
                             const int pwid, const int fwid)
{
    struct save_mainparams: public loop::outer, loop::null {
        typedef char *str;
        str *file_names;
        int markers_i, inds_i;
        
        save_mainparams(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::null(Top) { markers_i = inds_i = -1; }
        void make_file() {
            mssgvf("        STRUCTURE mainparams:         %s/%s\n", *_opath, file_names[5]);
            run_loop(file_names[5]);
        }
        void inner () {
            //pr_printf("#define MAXPOPS          %d\n", _Top->PedCnt);
            pr_printf("#define MAXPOPS          3\n");
            // In the documentation for Structure, we see that typical burnin times are in the range
            // of 10,000 to 100,000, so I think the defaults are too small.
            //pr_printf("#define BURNIN           2000\n");
            //pr_printf("#define NUMREPS          2000\n");
            // Let's instead use the larger defaults that were used in the example file distributed with Structure
            pr_printf("#define BURNIN           10000\n");
            pr_printf("#define NUMREPS          20000\n");
                                          
            pr_printf("#define INFILE           %s\n", file_names[0]);
            pr_printf("#define OUTFILE          %s\n", file_names[6]);
            
            pr_printf("#define NUMINDS          %d\n", numinds[++inds_i]);
            pr_printf("#define NUMLOCI          %d\n", markers_per_chromo[++markers_i]);
            pr_printf("#define PLOIDY           2\n");
            pr_printf("#define MISSING          0\n");
            // for diploid data, the two alleles for each locus are in consecutive order in
            // the same row, rather than being arranged in the same column, in two consecutive rows.
            pr_printf("#define ONEROWPERIND     1\n");
            
            pr_printf("#define LABEL            1\n");
            pr_printf("#define POPDATA          %s\n", (popdata != 0 ? "1" : "0"));
            pr_printf("#define POPFLAG          0\n");
            pr_printf("#define LOCDATA          0\n");
            pr_printf("#define PHENOTYPE        %s\n", (phenotype != 0 ? "1" : "0"));
            pr_printf("#define EXTRACOLS        0\n");
            
            pr_printf("#define MARKERNAMES      1\n");
            pr_printf("#define RECESSIVEALLELES 0\n");
            pr_printf("#define MAPDISTANCES     %s\n", (mapdistances == 1 ? "1" : "0"));
            
            pr_printf("#define PHASED           0\n");
            pr_printf("#define PHASEINFO        0\n");
            pr_printf("#define MARKOVPHASE      0\n");
            pr_printf("#define NOTAMBIGUOUS     -999\n");

        }
    } *sm = new save_mainparams(Top);
    
    sm->file_names = file_names;
    sm->load_formats(fwid, pwid, -1);
    sm->iterate();
    delete sm;
}

#if 0
static void write_extraparams(linkage_ped_top *Top, char *file_names[],
			      const int pwid, const int fwid)
{
    struct save_extraparams: public loop::outer, loop::null {
        typedef char *str;
        str *file_names;
        int markers_i;
        
        save_extraparams(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::null(Top) { markers_i = -1; }
        void make_file() {
            mssgvf("        STRUCTURE extraparams:         %s/%s\n", *_opath, file_names[7]);
            run_loop(file_names[7]);
        }
        void inner () {
            // Data taken from the file 'structure2.3.1_console/extraparams'...
            pr_printf("PROGRAM OPTIONS\n");
            pr_printf("#define NOADMIX     0\n");
            pr_printf("#define LINKAGE     0\n");
            pr_printf("#define USEPOPINFO  0\n");
            pr_printf("#define LOCPRIOR    0\n");
            pr_printf("#define FREQSCORR   1\n");
            pr_printf("#define ONEFST      0\n");
            pr_printf("#define INFERALPHA  1\n");
            pr_printf("#define POPALPHAS   0\n");
            pr_printf("#define ALPHA     1.0\n");
            pr_printf("#define INFERLAMBDA 0\n");
            pr_printf("#define POPALPHAS   0\n");
            pr_printf("#define ALPHA     1.0\n");
            pr_printf("#define INFERLAMBDA 0\n");
            pr_printf("#define POPSPECIFICLAMBDA 0\n");
            pr_printf("#define LAMBDA    1.0\n");
            pr_printf("PRIORS\n");
            pr_printf("#define FPRIORMEAN 0.01\n");
            pr_printf("#define FPRIORSD   0.05\n");
            pr_printf("#define UNIFPRIORALPHA 1\n");
            pr_printf("#define ALPHAMAX     10.0\n");
            pr_printf("#define ALPHAPRIORA   1.0\n");
            pr_printf("#define ALPHAPRIORB   2.0\n");
            pr_printf("#define LOG10RMIN     -4.0\n");
            pr_printf("#define LOG10RMAX      1.0\n");
            pr_printf("#define LOG10RPROPSD   0.1\n");
            pr_printf("#define LOG10RSTART   -2.0\n");
            pr_printf("USING PRIOR POPULATION INFO (USEPOPINFO)\n");
            pr_printf("#define GENSBACK    2\n");
            pr_printf("#define MIGRPRIOR 0.01\n");
            pr_printf("#define PFROMPOPFLAGONLY 0\n");
            pr_printf("LOCPRIOR MODEL FOR USING LOCATION INFORMATION\n");
            pr_printf("#define LOCISPOP      1\n");
            pr_printf("#define LOCPRIORINIT  1.0\n");
            pr_printf("#define MAXLOCPRIOR  20.0\n");
            pr_printf("OUTPUT OPTIONS\n");
            pr_printf("#define PRINTNET     1\n");
            pr_printf("#define PRINTLAMBDA  1\n");
            pr_printf("#define PRINTQSUM    1\n");
            pr_printf("#define SITEBYSITE   0\n");
            pr_printf("#define PRINTQHAT    0\n");
            pr_printf("#define UPDATEFREQ   100\n");
            pr_printf("#define PRINTLIKES   0\n");
            pr_printf("#define INTERMEDSAVE 0\n");
            pr_printf("#define ECHODATA     1\n");
            // ...
        }
    } *se = new save_extraparams(Top);
    
    se->file_names = file_names;
    se->load_formats(fwid, pwid, -1);
    se->iterate();
    delete se;
}
#endif /* 0 */

/**
 @brief Output an initial shell control file

 This file will run Structure setting up the environment that is necessary.

 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] analysis   A class instance associated with the current analysis.
 @param[in] file_names List of file neames used by this analysis mode.
 @return void
*/
static void write_sh(linkage_ped_top *Top,
                    analysis_type *analysis,
                    char *file_names[])
{
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1);
    
    struct all_sh: public sh_util {
        all_sh(linkage_ped_top *Top) : sh_util(Top) {}
        virtual void file_post() {
            chmod_X_file(path_);
        }
    } *sh = 0;
    
    if (top_shell) {
        sh = new all_sh(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }
    
    struct STRUCTURE_sh_script: public loop::outer, all_sh {
        typedef char *str;
        str *file_names;
        all_sh *sh;
        analysis_type *analysis;
        
        STRUCTURE_sh_script(linkage_ped_top *Top, analysis_type *analysis) : person_locus_entry(Top), loop::outer(Top), all_sh(Top) {
            this->analysis = analysis;
        }
        void make_file() {
            mssgvf("        STRUCTURE shell file:         %s/%s\n", *_opath, file_names[3]);
            run_loop(file_names[3]);
        }
        void file_header() {
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
	    // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
            fprintf_env_checkset_csh(_filep, "_STRUCTURE", "structure");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];
            pr_printf("# Process %s format data...\n", (*analysis)->_name);
            pr_printf("# For further information on Structure please see:\n");
            pr_printf("# http://pritch.bsd.uchicago.edu/structure.html\n");
            // If there isn't one then create an empty one.
            pr_printf("echo 'You may need to customize the extraparams file for your specific needs.'\n");
            pr_printf("touch extraparams\n");
#ifdef RUNSHELL_SETUP
            sprintf(cmd, "$_STRUCTURE -m %s\n", file_names[5]);
#else /* RUNSHELL_SETUP */
            pr_printf("structure -m %s\n", file_names[5]);
#endif /* RUNSHELL_SETUP */
            sh_run("STRUCTURE", cmd);
            fprintf_status_check_csh(_filep, "STRUCTURE", 1);
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new STRUCTURE_sh_script(Top, analysis);
    
    xp->file_names = file_names;
    xp->sh         = sh;
    xp->iterate();
    
    delete xp;
    
    if (top_shell) {
        mssgvf("        STRUCTURE top shell file:     %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("                  the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}

static void user_queries(char *file_names[],
                         int *combine_chromo,
                         analysis_type anal)
{
    int i, choice = -1, icc=-1, istem=-1;
    
    while (choice != 0) {
        print_outfile_mssg();
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        i=1;
        
        if (main_chromocnt > 1) {
            printf(" %d) Combine chromosomes?                      %s\n",
                   i, yorn[*combine_chromo]);
            icc=i++;
        }

        printf(" %d) Change file names stem?                   %s\n",
               i, anal->file_name_stem);
        istem=i++;
        
        printf("Enter options 0-%d > ", i-1);
        fcmap(stdin, "%d", &choice); printf("\n");
        
        if (choice == 0) {
            ; // OK...
        } else if (choice == icc) {
            *combine_chromo = TOGGLE(*combine_chromo);
        } else if (choice == istem) {
            printf("Enter new stem for the output file names > ");
            fcmap(stdin, "%s", anal->file_name_stem);    newline;
            // It doesn't matter what the parameter 'num' in the method file_names() is. It will get
            // changed to the appropriate thing later in the code. The method should be rewritten
            // globally without num and a place holder inserted instead.
            inner_file_names(file_names, (char *)"xx", anal->file_name_stem);
        } else {
            printf("Unknown option %d\n", choice);
        }
    }
}

void CLASS_STRUCTURE::create_output_file(linkage_ped_top *LPedTreeTop,
                                         analysis_type *analysis,
                                         char *file_names[],
                                         int untyped_ped_opt,
                                         int *numchr,
                                         linkage_ped_top **Top2)
{
    linkage_ped_top *Top = LPedTreeTop;
    int pwid, fwid, mwid;
    int combine_chromo = 0;
    
    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    combine_chromo = (main_chromocnt > 1) ? 1 : 0;
    
    if (main_chromocnt > 1) {
        // This is the batch file item that controls whether you wish to comnine the
        // chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
        // This is a derective from the user which will override the default...
        // NOTE that for interactive input, reorder_loci.cpp:ReOrderLoci() will ask the user
        // if they want to write chromosomes to one file or one file for each cromosome IF the user selects
        // multiple chromosomes in the preceeding menu. The answer will be stored in "Loop_Over_Chromosomes".
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].item_read)
            combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
    }

    // NOTE: comnine_chromo is passed in because the user can change it via a menu.
    if (InputMode == INTERACTIVE_INPUTMODE || (! DEFAULT_OUTFILES))
        user_queries(file_names, &combine_chromo, *analysis);

    LoopOverChrm = ! combine_chromo;
    
    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);

    // either the 1 chosen chromosome # or .all
    if (main_chromocnt > 1)
            change_output_chr(file_names[0], 0);
    else
            change_output_chr(file_names[0], global_chromo_entries[0]);
    
    // computed when writing 'INFILE', referenced in 'mainparams'...
    markers_per_chromo = (int *)calloc((size_t)main_chromocnt, sizeof(int));
    numinds = (int *)calloc((size_t)main_chromocnt, sizeof(int));

    write_INFILE(Top, file_names, pwid, fwid);
    write_mainparams(Top, file_names, pwid, fwid);
    //write_extraparams(Top, file_names, pwid, fwid);

    free(markers_per_chromo);
    free(numinds);

    write_sh(Top, analysis, file_names);
}
