/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
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
#include <ctime>

#include "common.h"
#include "typedefs.h"

#include "loop.h"
#include "sh_util.h"

#include "fcmap_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"

#include "write_merlin_ext.h"

void create_MERLIN_files(linkage_ped_top **LPedTop,
                       char *file_names[],
                       const int untyped_ped_opt,
                       const int output_format);

static void write_MERLIN_peds(linkage_ped_top *Top, char *file_names[],
                           const int pwid, const int fwid,
                           const bool has_x);

static void write_MERLIN_map(linkage_ped_top *Top, char *file_names[]);

static void write_MERLIN_data(linkage_ped_top *Top, char *file_names[]);

static void write_MERLIN_freq(linkage_ped_top *Top, char *file_names[]);

static void write_MERLIN_sh(linkage_ped_top *Top, char *file_names[]);

static int merlin_option_check(char *option,
                               merlin_opt_type *merlin_opt,
                               char *file_names[],
                               int *create_merlin_model);

static void init_merlin_opts(merlin_opt_type *merlin_opt);

static void merlin_options(char *option,
                           merlin_opt_type *merlin_opt,
                           char *file_names[],
                           int *create_liability_model);


static void inner_file_names(char **file_names, const char *num, const char *stem = "merlin");




static void write_MERLIN_peds(linkage_ped_top *Top, char **file_names,
                              const int pwid, const int fwid,
                              const bool has_x)
{

    vlpCLASS(merlin_ped, chr, ped_per_loci) {
        vlpCTOR(merlin_ped, chr, ped_per_loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Merlin Phenotype File:   %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }

        //Format of pedigree file: ID, Person, Father, Mother, Sex
        void per_start() {
            pr_id();
            //this appears to be the id again
            //pr_per();
            pr_father();
            pr_mother();
            pr_sex();
            //trait affection status?
            pr_aff();
        }

        void per_end() {
            pr_nl();
        }

        void inner() {
            pr_marker();
        }

    } *merlin_peds = new merlin_ped(Top);

    merlin_peds->file_names = file_names;

    merlin_peds->load_formats(fwid, pwid, -1);

    merlin_peds->iterate();

    delete merlin_peds;
}


static void write_MERLIN_map(linkage_ped_top *Top, char *file_names[])
{

    vlpCLASS(merlin_map,chr,loci) {
        vlpCTOR(merlin_map, chr, loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Merlin Map File:   %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "w");
        }

        void file_header(){
            pr_printf("CHR\tMARKER\tPOSITION");
            pr_nl();
        }

        //format of Merlin map: Chromosome, Marker, Position
        void inner() {
            int chr = _tlocusp->Marker->chromosome;
            str marker = _tlocusp->Marker->MarkerName;
            double pos = _tlocusp->Marker->pos_avg;

            pr_printf("%d\t",chr);
            pr_printf("%s\t",marker);
            pr_printf(" %f",pos);
            pr_nl();

        }
    } *merlin_maps = new merlin_map(Top);

    merlin_maps->file_names  = file_names;

    merlin_maps->iterate();

    delete merlin_maps;

}

//this should be the the Locus File according to the documentation
static void write_MERLIN_data(linkage_ped_top *Top, char *file_names[])
{
    /*
     * <contents of basic2.dat>
     * From the Merlin Documentation, the .dat format is as follows:
     * A  some_disease
     * T  some_trait
     * M  some_marker
     * M  another_marker
     *
     * A affection status
     * T quantitative trait
     *
     */


    vlpCLASS(merlin_dat,both,loci) {
        vlpCTOR(merlin_dat,both,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Merlin Data File:   %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
        }

        //this is for some reason causing a seg fault but I believe is the correct way to get this file looking the way it should

        void trait_start(){
            str trait = _ttraitp->LocusName;
            //str trait = _ttraitp->Pheno->TraitName;
            pr_printf("A ");
            pr_printf(trait);
            pr_nl();



//            int disease = _tlocusp->Type;
//            pr_printf("A %d", disease);
//            pr_nl();
        }



        void inner() {

            //need to figure out a way to split up these loops
            //it looks like in the merlin files we've been outputting we only want trait and marker
            //i.e. A and M can stay
            //then it should be a list like
            //a  trait
            //m m1
            //m m2 etc.

            //moved to trait_start()
//            str trait = _ttraitp->Pheno->TraitName;
//            pr_printf("T ");
//            pr_printf(trait);
//            pr_nl();

            str marker = _tlocusp->Marker->MarkerName;
            pr_printf("M ");
            pr_printf(marker);
            pr_nl();

        }

    } *merlin_dats = new merlin_dat(Top);

    merlin_dats->file_names = file_names;

    merlin_dats->iterate();

    delete merlin_dats;
}

//Create a freq file? Mega2 currently outputs one
//Need to find where it is in the Merlin documentation
static void write_MERLIN_freq(linkage_ped_top *Top, char *file_names[])
{
    vlpCLASS(merlin_freq,chr,loci) {
        vlpCTOR(merlin_freq, chr, loci) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Merlin Frequency file:   %s/%s\n", *_opath, file_names[3]);
            data_loop(*_opath, file_names[3], "w");
        }
        /*
         * From the Merlin Documentation, the format is as follows:
         *<contents of basic2.freq>
         * M some_marker
         * F 0.1
         * F 0.2
         * M another_marker
         * F 0.6
         * F 0.4
         *
         * It appears this was incorrect and the original merlin output uses the extended frequency format
         * that goes:
         * M some_marker
         * A size_of_base_pairs frequency
         */

        //loci start is not defined for loci...
        void loci_start(){
            str marker = _tlocusp->Marker->MarkerName;
            pr_printf("M ");
            pr_printf(marker);
            pr_nl();
        }

        void inner (){
            str marker = _tlocusp->Marker->MarkerName;
            pr_printf("M ");
            pr_printf(marker);
            pr_nl();

            double freq = _tlocusp->Allele->Frequency;
            //const char * name = _tlocusp->Allele->AlleleName;
            pr_printf("A %s %f\n", _tlocusp->Allele->AlleleName, freq);
        }

    } *merlin_freqs = new merlin_freq(Top);

    merlin_freqs->file_names = file_names;

    merlin_freqs->iterate();

    delete merlin_freqs;

}

//this will generate a c-shell script
static void write_MERLIN_sh(linkage_ped_top *Top, char *file_names[])
    {
        vlpCLASS(merlin_sh, both, sh_exec) {
            vlpCTOR(merlin_sh, both, sh_exec) { }

            typedef char *str;
            str *file_names;
            sh_exec *sh;
            bool has_x;

            void file_loop() {
                mssgvf("        Merlin shell file:       %s/%s\n", *_opath, file_names[4]);
                data_loop(*_opath, file_names[4], "w");
            }

            void file_header() {
                pr_printf("#!/bin/csh -f\n");
                pr_printf("#----------------------------------------------\n");
                //this should be the curren shell scripts file name
                pr_printf("# C-shell file name: %s\n",file_names[4]);
                //how do I get a reference to the current version?
                pr_printf("#   Mega2 version 4.8.2\n");
                //current date/time
                //yyyy-mm-dd-hh-mm
                time_t rawtime1;
                struct tm *timeinfo1;
                char buffer1 [50];
                time (&rawtime1);
                timeinfo1 = localtime (&rawtime1);
                strftime(buffer1,50,"%Y-%m-%d-%h-%M",timeinfo1);
                pr_printf("#   Run date:                %s\n",buffer1);
                //current date/time again, but this time formatted differently!
                //{day of week abreviated}{space}{month abreviated}{space}{date}{space}{hh:mm}{space}{yyyy}
                time_t rawtime2;
                struct tm *timeinfo2;
                char buffer2 [100];
                time (&rawtime2);
                timeinfo2 = localtime (&rawtime2);
                strftime(buffer2,100, "%c",timeinfo2);
                pr_printf("#   This script created on   %s\n",buffer2);
                pr_printf("#   Input file names:\n");
                //perhaps?
                pr_printf("#       Pedigree file:              %s\n",Input->input_files.pedfl[0]);
                pr_printf("#          Locus file:              %s\n",Input->input_files.locusfl[0]);
                pr_printf("#            Map file:              %s\n",Input->input_files.locusfl[0]);
                //these need to come in from inputs
                pr_printf("#    Untyped pedigree option  Include all pedigrees whether typed or not\n");
                pr_printf("#----------------------------------------------\n");
                pr_printf("# Chromosome number:    %d\n",_numchr);

            }


            void inner() {
                pr_printf("if ( ! $?_MERLIN) then\n");
                pr_printf("\tset _MERLIN='merlin'\n");
                pr_printf("\tendif\n");
                pr_printf("if ( ! $?_MERLIN2SW2) then\n");
                pr_printf("\tset _MERLIN2SW2='merlin2sw2.pl'\n");
                pr_printf("\tendif\n");
                pr_printf("if ( ! $?_SIMWALK2) then\n");
                pr_printf("\tset _SIMWALK2='simwalk2'\n");
                pr_printf("\tendif\n");
                pr_printf("if (-e merlin_out.%d) then\n",_numchr);
                pr_printf("\t/bin/rm merlin_out.%d\n",_numchr);
                pr_printf("\tendif\n");
                pr_printf("if (-e merlin_table.%d) then\n",_numchr);
                pr_printf("\t/bin/rm merlin_table.%d\n",_numchr);
                pr_printf("\tendif\n");
                pr_printf("  if (-e merlin_out.%d.pdf) then\n",_numchr);
                pr_printf("\t/bin/rm merlin_out.%d.pdf\n",_numchr);
                pr_printf("\tendif\n");

                //actual merlin command
                pr_printf("echo Running merlin on %s and  %s\n",file_names[0],file_names[2]);
                pr_printf("$_MERLIN -p %s -d %s -m %s -f %s --npl --pairs --tabulate  --markerNames >> merlin_out.%d\n"
                        ,file_names[0],file_names[2],file_names[1],file_names[3],_numchr);
                //?
                pr_printf("set merlin_status=$status\n");
                pr_printf("echo $merlin_status > merlin_status\n");
                pr_printf("if ($merlin_status != 0) then\n");
                pr_printf("  echo \"ERROR: Run of 'merlin' failed with status code: $merlin_status\"\n");
                pr_printf("endif\n");
                pr_printf("cat merlin_out.%d\n",_numchr);
                pr_printf("echo ===================== >> merlin_out.%d\n",_numchr);
                pr_printf("cat %s >> merlin_out.%d\n",file_names[1],_numchr);
                pr_printf("echo Created output file merlin_out.%d\n",_numchr);
            }

            void file_post() {
                chmod_X_file(path_);
            }

        } *merlin_shs = new merlin_sh(Top);

        merlin_shs->file_names = file_names;

        merlin_shs->iterate();

        delete merlin_shs;

}



void CLASS_NEWMERLIN::create_output_file(
        linkage_ped_top *LPedTreeTop,
        analysis_type *analysis,
        char *file_names[],
        int untyped_ped_opt,
        int *numchr,
        linkage_ped_top **Top2)
{
    int pwid, fwid, mwid;
    int combine_chromo=0;
    char prefix[100];
    linkage_ped_top *Top = LPedTreeTop;
    bool use_map;

    use_map = true;


    //add in code to call the merlin options menu...
    merlin_opt_type merlin_opt;
    init_merlin_opts(&merlin_opt);

    /* Input merlin options */
    if (analysis == TO_MERLINONLY) {
        merlin_options(&(opts[0]), &(merlin_opt), file_names, &model_file);

        if (!DEFAULT_OPTIONS) {
            if (strcmp(opts, "")) {
                sprintf(err_msg, "User options to Merlin: %s.", opts);
                mssgf(err_msg);
            } else {
                if (HasAff) {
                    sprintf(opts, "--npl --pairs");
                    mssgf("Using default Merlin options --npl and --pairs");
                    merlin_opt.npl=merlin_opt.pairs=1;
                    merlin_opt.markernames=0;
                } else if (HasQuant) {
                    sprintf(opts, "--qtl");
                    merlin_opt.qtl=1;
                    mssgf("Using default Merlin option --qtl.");
                } else {
                    sprintf(opts, "--ibd");
                    mssgf("Using default Merlin option --ibd.");
                }
            }
            log_line(mssgf);
        }
    }

    for (i=0; i < num_traits; i++) {
        SKIP_TRI(i)
        if ((*LPedTop)->LocusTop->Locus[global_trait_entries[i]].Type == AFFECTION ||
            (*LPedTop)->LocusTop->Locus[global_trait_entries[i]].Type == QUANT) {
            num_affec++;
        }
    }

    if ((analysis == TO_MERLINONLY) &&
        (merlin_opt.npl==1 ||  merlin_opt.pairs==1 ||
         merlin_opt.qtl==1 ||  merlin_opt.vc==1 ||
         merlin_opt.parametric == 1)) {
        R_setup = 1;
    }




    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    //
    // Set the default value...
    // write everything to one file (unless only one chromosome has been selected).

    combine_chromo = main_chromocnt > 1;

    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    LoopOverChrm = ! combine_chromo;

    //adding in looping over traits for the data file
    LoopOverTrait = num_traits != 1;

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    write_MERLIN_map(Top, file_names);

    write_MERLIN_data(Top, file_names);

    write_MERLIN_freq(Top,file_names);

    write_MERLIN_peds(Top, file_names, pwid, fwid,
                      use_map);

    write_MERLIN_sh(Top, file_names);


}



/*
 * The following functions merlin_option_check, init_merlin_options, and merlin_options all came from write_premakeped.cpp
 * THe hope is that they can be refitted to create the correct option menu functionality for this "newmerlin"
 * */

static int merlin_option_check(char *option,
                               merlin_opt_type *merlin_opt,
                               char *file_names[],
                               int *create_merlin_model)

{
    char opt_copy[FILENAME_LENGTH];
    char opt_output[FILENAME_LENGTH];
    char *nextopt, prevopt[100];
    float num;
    int valid, i;
    char y[10];


    char knownopts[NUM_MERLIN_OPTS][20] =
            {
                    "--error", "--information", "--likelihood", "--model",
                    "--ibd", "--kinship", "--matrices", "--extended", "--select",
                    "--npl", "--pairs", "--qtl", "--deviates", "--exp",
                    "--vc", "--useCovariates", "--ascertainment", "--unlinked",
                    "--infer", "--assoc", "--fastAssoc", "--filter", "--custom",
                    "--best", "--sample", "--all", "--founders", "--horizontal",
                    "--zero", "--one", "--two", "--three", "--singlepoint",
                    "--steps", "--maxStep", "--minStep", "--grid", "--start", "--stop",
                    "--clusters", "--distance", "--rsq", "--cfreq",
                    "--bits", "--megabytes", "--minutes",
                    "--trim", "--noCoupleBits", "--swap", "--smallSwap",
                    "--quiet", "--markerNames", "--frequencies", "--perFamily", "--pdf",
                    "--tabulate", "--prefix",
                    "--simulate", "--reruns", "--save", "--trait"
            };

    init_merlin_opts(merlin_opt);


    strcpy(opt_copy, option);
    nextopt=strtok(opt_copy, " ");
    strcpy(opt_output, "");

    while (nextopt != NULL) {
        if (!strncmp(nextopt, "--",(size_t) 2)) {
            valid=0;
            for (i=0; i < NUM_MERLIN_OPTS; i++) {
                if (!strcmp(nextopt, knownopts[i])) {
                    valid=1; break;
                }
            }
            if (!valid) { sprintf(err_msg, "Unknown option %s.\n", nextopt); warnf(err_msg);}
        }
        nextopt = strtok(NULL, " ");
    }

    /* Now re-initialize the string */
    strcpy(opt_copy, option);
    nextopt=strtok(opt_copy, " ");

    while(nextopt != NULL) {
        strcat(opt_output, nextopt);
        strcat(opt_output, " ");
        /* Check that these options are followed by another arguments */
        if (!strcmp(nextopt, "--unlinked") || !strcmp(nextopt, "--steps") ||
            !strcmp(nextopt, "--maxStep") || !strcmp(nextopt, "--minStep") ||
            !strcmp(nextopt, "--grid") || !strcmp(nextopt, "--start") ||
            !strcmp(nextopt, "--stop") || !strcmp(nextopt, "--bits") ||
            !strcmp(nextopt, "--rsq") || !strcmp(nextopt, "--distance") ||
            !strcmp(nextopt, "--reruns") ||
            !strcmp(nextopt, "--megabytes") || !strcmp(nextopt, "--minutes")) {

            strcpy(prevopt, nextopt);
            nextopt = strtok(NULL, " ");
            if (nextopt != NULL) {
                if (sscanf(nextopt, "%f", &num) != 1) {
                    printf("ERROR: Option %s needs a numeric argument.\n", prevopt);
                    return 0;
                }
            } else {
                printf("ERROR: Option %s needs a numeric argument.\n", prevopt);
                return 0;
            }
            strcat(opt_output, nextopt); strcat(opt_output, " ");
        } else if (!strcmp(nextopt, "--model") || !strcmp(nextopt, "--custom")) {
            strcpy(prevopt, nextopt);
            /* Before reading in the next option we need to figure out whether
               this is a constant or another option */
            nextopt = strtok(NULL, " ");
            if (nextopt == NULL || ! strncmp(nextopt, "--", (size_t) 2)) {
                printf("ERROR: Option %s needs a file name argument.\n", prevopt);
                return 0;
            }

            if (!strcmp(prevopt, "--model")) {
                shorten_path(nextopt, file_names[14]);
                if (access(file_names[14], F_OK) != 0) {
                    if (*create_merlin_model) {
                        printf("File %s will be created from traits automatically.\n",
                               file_names[14]);
                    } else {
                        printf("WARNING: File %s not found, Merlin may not run correctly.\n",
                               file_names[14]);

                        printf("Create model file(s) %s in output folder? [y/n](default n)",
                               file_names[14]);
                        fflush(stdout);
                        IgnoreValue(fgets(y, 9, stdin)); newline; fflush(stdin);
                        if (y[0] == 'Y' || y[0] == 'y') {
                            *create_merlin_model=1;
                        } else {
                            printf("WARNING: No model file, Merlin may not run correctly.\n");
                        }
                    }
                }
                if (strcmp(file_names[14], "None")) {
                    strcat(opt_output, file_names[14]); strcat(opt_output, " ");
                }
                merlin_opt->parametric = 1;
            } else {
                strcat(opt_output, nextopt); strcat(opt_output, " ");
                if (access(nextopt, F_OK) != 0) {
                    printf("WARNING: file %s provided for %s does not exist or is unreadable.\n",
                           prevopt, nextopt);
                }
            }
        } else if (!strcmp(nextopt, "--prefix")) {
            nextopt = strtok(NULL, " ");
            if (nextopt == NULL || !strncmp(nextopt, "--", (size_t) 2)) {
                printf("ERROR: Option %s needs a string name argument.\n", prevopt);
                return 0;
            }
            strcat(opt_output, nextopt); strcat(opt_output, " ");
            warnf("Mega2 generated shell script recognizes only the default prefix \"merlin\"");
            warnf("Merlin shell script will fail with alternate prefixes.");
            merlin_opt->prefix=1;
            /* otherwise process this string as an option */
        } else {
            if (!strcmp(nextopt, "--npl")) { merlin_opt->npl=1; }
            if (!strcmp(nextopt, "--pairs")) { merlin_opt->pairs=1; }
            if (!strcmp(nextopt, "--qtl")) { merlin_opt->qtl=1; }
            if (!strcmp(nextopt, "--vc")) { merlin_opt->vc=1; }
            if (!strcmp(nextopt, "--markerNames")) { merlin_opt->markernames=1; }
            if (!strcmp(nextopt, "--pdf")) { merlin_opt->pdf=1; }
            if (!strcmp(nextopt, "--tabulate")) { merlin_opt->tabulate=1; }
        }
        nextopt = strtok(NULL, " ");
    }

    strcpy(option, opt_output);

    return 1;

}

static void init_merlin_opts(merlin_opt_type *merlin_opt)
{
    merlin_opt->npl=
    merlin_opt->qtl =
    merlin_opt->pairs =
    merlin_opt->vc =
    merlin_opt->parametric =
    merlin_opt->pdf =
    merlin_opt->markernames =
    merlin_opt->prefix =
    merlin_opt->tabulate = 0;
}

static void merlin_options(char *option,
                           merlin_opt_type *merlin_opt,
                           char *file_names[],
                           int *create_liability_model)

{
    int done_ = 0;
    char cdone_[10];


    if (HasAff) {
        strcpy(option, "--npl --pairs --tabulate");
    } else if (HasQuant) {
        strcpy(option, "--qtl --tabulate");
    } else {
        strcpy(option, "--ibd");
    }

    *create_liability_model = 0;

    if (DEFAULT_OPTIONS) {
        /* This is called only to set the correct flags in merlin_opt structure */
        merlin_option_check(option, merlin_opt, file_names, create_liability_model);
    } else {
        done_ = -1;
        while(done_) {
            /* New merlin run parameters menu */
            draw_line();
            printf("Merlin run options menu:\n");
            printf("0) Done with this menu, please proceed.\n");
            printf(" 1) Select Merlin option [%s].\n", option);
            if (strcmp(file_names[14], "None")) {
                printf(" 2) %s Merlin model file %s from affection trait parameters [%s].\n",
                       ((access(file_names[14], F_OK) == 0)?
                        "Recreate" : "Create"),
                       file_names[14],
                       yorn[*create_liability_model]);
            } else {
                printf(" 2) Create Merlin model file from affection trait parameters [%s].\n",
                       yorn[*create_liability_model]);
            }
            if (*create_liability_model && strcmp(file_names[14], "None")) {
                printf("    To use this model file, use the Merlin analysis option:\n");
                printf("       --model %s\n", file_names[14]);
            }
            printf("Select 0-2 (2 to toggle) > ");
            fcmap(stdin, "%s", cdone_);
            newline;
/*       done_ = atoi(cdone_); */
            done_ = (int)strtol(cdone_, (char **)NULL, 10);
/*       printf("%d\n", done_); */
/*       sleep(1); */
            switch(done_) {
                case 0:
                    if (!merlin_option_check(option, merlin_opt, file_names, create_liability_model)) {
                        printf("Error: bad option string, discarding options.\n");
                        strcpy(option, "");
                    }
                    break;
                case 1:
                    draw_line();
                    printf("Merlin option selection:\n");
                    printf("Enter one or more of the options listed below, ");
                    printf("separated by spaces, \n");
                    printf("   e.g. --steps 2 --bits 32 \n");
                    printf("Some options must be followed by a number as indicated\n");

                    printf("For no options press <RETURN> or <Enter>\n\n");

                    printf("         General : --error, --information, --likelihood, --model\n");
                    printf("      IBD States : --ibd, --kinship, --matrices, --extended, --select\n");
                    printf("     NPL Linkage : --npl, --pairs, --qtl, --deviates, --exp\n");
                    printf("     VC Linkage  : --vc, --useCovariates, --ascertainment\n");
                    printf("     Association : --infer, --assoc, --fastAssoc, --filter, --custom\n");
                    printf("     Haplotyping : --best, --sample, --all, --founders, --horizontal\n");
                    printf("   Recombination : --zero, --one, --two, --three, --singlepoint\n");
                    printf("       Positions : --steps, --maxStep, --minStep, --grid, --start, --stop\n");
                    printf(" Marker Clusters : --clusters, --distance, --rsq, --cfreq\n");
                    printf("          Limits : --bits, --megabytes, --minutes\n");
                    printf("     Performance : --trim, --noCoupleBits, --swap, --cache\n");
                    printf("          Output : --quiet, --markerNames, --frequencies, --perFamily, --pdf,\n");
                    printf("                   --prefix\n");
                    printf("      Simulation : --simulate, --reruns, --save\n");

                    printf("These options can be followed by one additional argument each:\n");
                    printf(" model, custom, steps, maxStep, minStep, grid, start, stop, \n");
                    printf(" clusters, megabytes, minutes, cache, prefix\n");
                    printf("\nEnter options as text (at most 200 characters) > ");
                    strcpy(option, "");
                    fflush(stdout);
                    IgnoreValue(fgets(option, FILENAME_LENGTH-1, stdin)); newline;
                    option[strlen(option)-1]='\0';
                    if (!merlin_option_check(option, merlin_opt, file_names, create_liability_model)) {
                        printf("Error: bad option string, discarding options.\n");
                        strcpy(option, "");
                    }
                    break;
                case 2:
                    printf("Set value to \"yes\" to automatically generate model file\n");
                    printf("using affection status trait locus parameters from input files.\n");
                    *create_liability_model = TOGGLE(*create_liability_model);
                    break;
                default:
                    warn_unknown(cdone_);
                    break;
            }
        }
    }
}




void CLASS_NEWMERLIN::get_file_names(char *file_names[], char *prefix,
                                int has_orig, int has_uniq, int *combine_chromo)
{
    int i, choice;
    char fl_stat[12];
    int igl, ipre, iphen, ish, ioui, ioup, isum, isumf;
    analysis_type analysis = this;

    strcpy(prefix, "fbat");

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        choice = 0;
    } else {
        choice = -1;
    }
    if (main_chromocnt > 1) {
        // This is the batch file item that controls whether you wish to comnine the
        // chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
        // This is a derective from the user which will override the default...
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
            *combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
        else
            *combine_chromo=0;
    }

    if (main_chromocnt > 1 && *combine_chromo) {
        // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
        analysis->replace_chr_number(file_names, 0);
    } else {
        analysis->replace_chr_number(file_names, global_chromo_entries[0]);
    }

    /* output file name menu */
    igl = ipre = iphen = ish = ioui = ioup = isum = isumf = -1;
    // If the default output files are used we do not go here.
    // Otherwise, enter with choice -- -1.
    while (choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("Output file names menu:\n");
        printf("0) Done with this menu - please proceed\n");
        i=1;

        if (main_chromocnt > 1) {
            printf(" %d) Combine chromosomes?                      %s\n",
                   i, yorn[*combine_chromo]);
            igl=i++;
        }

        if (num_traits > 2 && LoopOverTrait == 0) {
            printf(" %d) Phenotype file name                       %-15s\t%s\n",
                   i, file_names[2],
                   file_status(file_names[2], fl_stat));
            iphen=i++;
        }

        printf(" %d) File name stem:                           %-15s\n", i, prefix);

        ipre=i++;

        printf(" %d) Shell file name:                          %-15s\t%s\n",
               i, file_names[3],
               ((main_chromocnt <= 1 || *combine_chromo == 1) ?
                file_status(file_names[3], fl_stat) : ""));
        ish=i++;

        individual_id_item(i, analysis, OrigIds[0], 43, 2,0, 0);
        ioui=i++;

        pedigree_id_item(i, analysis, OrigIds[1], 43, 2, 0);
        ioup=i;

        printf("Enter options 0-%d > ", i);
        fcmap(stdin, "%d", &choice); printf("\n");
        test_modified(choice);

        if (choice < 0) {
            printf("Unknown option %d\n", choice);
        } else if (choice == 0) {
            ;

        } else if (choice == igl) {
            *combine_chromo = TOGGLE(*combine_chromo);
            if (main_chromocnt > 1 && *combine_chromo) {
                // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
                analysis->replace_chr_number(file_names, 0);
            } else {
                analysis->replace_chr_number(file_names, global_chromo_entries[0]);
            }

        } else if (choice == ipre) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);    newline;
            inner_file_names(file_names, "", prefix);
            if (main_chromocnt > 1 && *combine_chromo)
                analysis->replace_chr_number(file_names, 0);
            else
                analysis->replace_chr_number(file_names, global_chromo_entries[0]);

        } else if (choice == iphen) {
            printf("Enter new phenotype file name > ");
            fcmap(stdin, "%s", file_names[2]);    newline;

        } else if (choice == ish) {
            printf("Enter new shell script name %s > ", file_names[3]);
            fcmap(stdin, "%s", file_names[3]);    newline;

        } else if (choice == ioui) {
            OrigIds[0] = individual_id_item(0, analysis, OrigIds[0], 35, 1, has_orig, has_uniq);
            individual_id_item(0, analysis, OrigIds[0], 0, 3, has_orig, has_uniq);

        } else if (choice == ioup) {
            OrigIds[1] = pedigree_id_item(0, analysis, OrigIds[1], 35, 1, has_orig);
            pedigree_id_item(0, analysis, OrigIds[1], 0, 3, has_orig);

        } else {
            printf("Unknown option %d\n", choice);
        }
    }
}

    static void inner_file_names(char **file_names, const char *num, const char *stem) {

        //this should reflect the original naming convention of merlin in mega2
        sprintf(file_names[0], "%s_ped.%s",stem,  num);
        sprintf(file_names[1], "%s_map.%s",stem, num);
        sprintf(file_names[2], "%s_data.%s",stem, num);
        sprintf(file_names[3], "%s_freq.%s",stem,  num);
        sprintf(file_names[4], "%s.%s.sh",stem, num);
    }

void CLASS_NEWMERLIN::gen_file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_NEWMERLIN::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
}
