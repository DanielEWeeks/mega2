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
            pr_per();
            pr_father();
            pr_mother();
            pr_sex();
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

        //format of Merlin map: Chromosome, Marker, Position
        void inner() {

            int chr = _tlocusp->Marker->chromosome;
            str marker = _tlocusp->Marker->MarkerName;
            double pos = _tlocusp->Marker->pos_avg;

            pr_printf("%d ",chr);
            pr_printf(marker);
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


    vlpCLASS(merlin_dat,chr,loci) {
        vlpCTOR(merlin_dat,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Merlin Data File:   %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
        }
        void inner() {

            int disease = _tlocusp->Type;
            pr_printf("A %d", disease);
            pr_nl();

            //this seems correct for trait

            str trait = _ttraitp->Pheno->TraitName;
            pr_printf("T ");
            pr_printf(trait);
            pr_nl();

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
         */
        void inner (){
            str marker = _tlocusp->Marker->MarkerName;
            pr_printf("M ");
            pr_printf(marker);
            pr_nl();

            double freq = _tlocusp->Allele->Frequency;
            pr_printf("F %d\n", freq);
        }

    } *merlin_freqs = new merlin_freq(Top);

    merlin_freqs->file_names = file_names;

    merlin_freqs->iterate();

    delete merlin_freqs;

}


//this is what the original merlin make shell file looked like, may need to add merlin options looks like

///* A single overall shell script is created for all traits and all chromosomes
//Then, for each trait and each chromosome, we also create a separate
//shell-script.
//*/
//static void merlin_shell_script(int numchr, int global,
//                                char *output_files[],
//                                char *opts,
//                                linkage_locus_top *LTop,
//                                int qtdt,
//                                merlin_opt_type merlin_opt)
//
//

//this will generate a c-shell script
static void write_MERLIN_sh(linkage_ped_top *Top, char *file_names[]) {

    //this is the original shell script code
    int nloop, tr, *trp, num_affec=num_traits;
    FILE *filep;
    char mfl[2*FILENAME_LENGTH];
    char *syscmd;
    char chr_str[3];
    char merlinoutfile[FILENAME_LENGTH];
    char merlinpdffile[FILENAME_LENGTH];
    char merlintablefile[FILENAME_LENGTH];
    char global_shell[FILENAME_LENGTH];
    char merlin_prog[7];
    int xlinked;
    int tabulate=0;
#ifdef RUNSHELL_SETUP
    char            merlin_prog_env[8];
#endif /* RUNSHELL_SETUP */

    xlinked=(((LTop->SexLinked == 2 && numchr == SEX_CHROMOSOME) ||
              (LTop->SexLinked == 1))? 1 : 0);

#ifdef RUNSHELL_SETUP
    strcpy(merlin_prog_env, (xlinked ? "_MINX" : "_MERLIN"));
#endif /* RUNSHELL_SETUP */
    strcpy(merlin_prog, (xlinked? "minx" : "merlin"));

    if (strcmp(output_paths[0], ".")) {
        sprintf(global_shell, "%s/%s", output_paths[0], output_files[13]);
    } else {
        strcpy(global_shell, output_files[13]);
    }
    change_output_chr(global_shell, 0);

    NLOOP;

    CHR_STR(numchr, chr_str);
    MERLINOUT(chr_str, merlinoutfile);
    MERLINPDF(chr_str, merlinpdffile);
    MERLINTABLE(chr_str, merlintablefile);

    trp = (num_traits > 0) ? &(global_trait_entries[0]) : NULL;

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && global == 1) {
            /*  create a top level shell script */
            if ((filep = fopen(global_shell, "w")) == NULL) {
                errorvf("could not open shell file %s.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "#!/bin/csh -f\n");
            /* print some identification */
            fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    global_shell
#endif
            );
            script_time_stamp(filep);
            fclose(filep);
            syscmd = CALLOC((strlen(global_shell) + 10), char);
            sprintf(syscmd, "chmod +x %s", global_shell);
            System(syscmd);
            free(syscmd);
        }

        if (nloop > 1 && tr == 0) {
            continue;
        }

        sprintf(mfl, "%s/%s", output_paths[tr], Outfile_Names[13]);
        if ((filep = fopen(mfl, "w")) == NULL) {
            errorvf("Unable to open %s for writing.\n", mfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "#!/bin/csh -f\n");
        /* print some identification */
        fprintf(filep, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                NOPATH
#else
                Outfile_Names[13]
#endif
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
        fprintf(filep, " if (-e %s) then\n", merlinoutfile);
        fprintf(filep, "    /bin/rm %s\n", merlinoutfile);
        fprintf(filep, "endif\n");
        fprintf(filep, " if (-e %s) then\n", merlintablefile);
        fprintf(filep, "    /bin/rm %s\n", merlintablefile);
        fprintf(filep, "endif\n");
        fprintf(filep, " if (-e %s) then\n", merlinpdffile);
        fprintf(filep, "    /bin/rm %s\n", merlinpdffile);
        fprintf(filep, "endif\n");

        fprintf(filep,
                "echo Running merlin on %s and  %s \n",
                output_files[10], output_files[11]);
#ifdef RUNSHELL_SETUP
        fprintf(filep, "$%s -p %s -d %s ", merlin_prog_env, output_files[10], output_files[11]);
#else /* RUNSHELL_SETUP */
        fprintf(filep, "%s -p %s -d %s ", merlin_prog, output_files[10], output_files[11]);
#endif /* RUNSHELL_SETUP */

        if (qtdt) {
            fprintf(filep, "-m %s -f %s ", output_files[12], output_files[9]);
        }

        if (opts == NULL) {
            if (trp==NULL) {
                fprintf(filep, "--ibd >> %s\n", merlinoutfile);
            } else if (LTop->Locus[*trp].Type == QUANT) {
                fprintf(filep, " --qtl --tabulate >> %s\n", merlinoutfile);
                tabulate=1;
                trp++;
            } else if (LTop->Locus[*trp].Type == AFFECTION) {
                fprintf(filep, "--npl --pairs --tabulate >> %s\n", merlinoutfile);
                tabulate=1;
                trp++;
            }
        } else {
            /* Add tabulate if one of the following options are selected */
            tabulate = merlin_opt.tabulate;

            if ((merlin_opt.npl == 1 || merlin_opt.pairs == 1 || merlin_opt.qtl == 1 ||
                 merlin_opt.parametric == 1 || merlin_opt.vc == 1)) {
                /* tabulate should always be 1 */
                if (merlin_opt.tabulate != 1) {
                    fprintf(filep, "%s --tabulate >> %s\n", opts, merlinoutfile);
                } else {
                    fprintf(filep, "%s >> %s\n", opts, merlinoutfile);
                }
                tabulate = 1;
            } else {
                fprintf(filep, "%s >> %s\n", opts, merlinoutfile);
            }
        }
        fprintf_status_check_csh(filep, "merlin", 0);

        fprintf(filep, "cat %s\n", merlinoutfile);
        fprintf(filep, "echo ===================== >> %s\n",
                merlinoutfile);
        fprintf(filep, "cat %s >> %s\n", output_files[12],
                merlinoutfile);
        fprintf(filep, "echo Created output file %s\n",
                merlinoutfile);

        if (tabulate == 1) {
            if (merlin_opt.pairs == 1 || merlin_opt.npl == 1 ||
                merlin_opt.qtl == 1) {
                fprintf(filep, "if (-e merlin-nonparametric.tbl) then\n");
                fprintf(filep, "   mv merlin-nonparametric.tbl %s\n", merlintablefile);
                fprintf(filep, "   echo Renamed merlin-nonparametric.tbl to %s\n",
                        merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");
            } else if (merlin_opt.parametric == 1) {
                fprintf(filep, "if (-e merlin-parametric.tbl) then\n");
                fprintf(filep, "   mv merlin-parametric.tbl %s\n", merlintablefile);
                fprintf(filep, "   echo Renamed merlin-parametric.tbl to %s\n", merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");
            } else if (merlin_opt.vc == 1) {
                fprintf(filep, "if (-e merlin-vc-chr%2s.tbl) then\n", chr_str);
                fprintf(filep, "   mv merlin-vc-chr%2s.tbl %s\n", chr_str, merlintablefile);
                fprintf(filep, "   echo Renamed merlin-vc-chr%2s.tbl to %s\n", chr_str, merlintablefile);
                fprintf(filep, "else \n");
                fprintf(filep, "   echo Computation failed!\n");
                fprintf(filep, "   exit(-1)\n");
                fprintf(filep, "endif\n");

            } else {
                sprintf(err_msg, "Tabulate option not available for selected options %s.\n",
                        opts);
                errorf(err_msg);
            }
        }
        if (merlin_opt.pdf == 1) {
            fprintf(filep, "if (-e merlin.pdf) then \n");
            fprintf(filep, "   mv merlin.pdf %s\n", merlinpdffile);
            fprintf(filep, "   echo Renamed merlin.pdf to %s\n", merlinpdffile);
            fprintf(filep, "else \n");
            fprintf(filep, "   echo Failed to create merlin pdf file!\n");
            fprintf(filep, "endif\n");

        }
        fclose(filep);
        syscmd = CALLOC((strlen(mfl) + 10), char);
        sprintf(syscmd, "chmod +x %s", mfl);
        System(syscmd);
        free(syscmd);

        if (global) {
            if ((filep = fopen(global_shell, "a")) == NULL) {
                errorvf("Unable to open %s for writing.\n", global_shell);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(filep, "echo Running Merlin shell script %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    mfl
#endif
            );
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd %s\n", trait_paths[tr]);
            }
            fprintf(filep, "./%s\n", Outfile_Names[13]);
            if (strcmp(trait_paths[tr], ".")) {
                fprintf(filep, "cd ..\n");
            }
            fclose(filep);
        }

        if (nloop == 1) break;

    }


//    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
//                    strcmp(output_paths[0], ".");
//
//    dataloop::sh_exec *sh = 0;
//    if (top_shell) {
//        sh = new dataloop::sh_exec(Top);
//        sh->filep_open(output_paths[0], file_names[4], "w");
//        sh->sh_main();
//    }
//    vlpCLASS(merlin_sh, both, sh_exec) {
//        vlpCTOR(merlin_sh, both, sh_exec) { }
//
//        typedef char *str;
//        str *file_names;
//        sh_exec *sh;
//        bool has_x;
//
//        void file_loop() {
//            mssgvf("        Merlin shell file:       %s/%s\n", *_opath, file_names[3]);
//            data_loop(*_opath, file_names[3], "w");
//        }
//
//        void file_header() {
//            if (sh)
//                sh->sh_sh(this);
//            sh_shell_type();
//            sh_id();
//            script_time_stamp(_filep);
//            #ifdef RUNSHELL_SETUP
//            // This handles the environment variable setup to allow the checking
//            // functions in 'batch_run' to work correctly...
//            fprintf_env_checkset_csh(_filep, "_FBAT", "fbat");
//            #endif /* RUNSHELL_SETUP */
//        }
//
//
//        void inner () {
//            char out_fl[2*FILENAME_LENGTH];
//
//            if (_numchr > 0) {
//                char strchr[4];
//                CHR_STR(_numchr, strchr);
//                pr_printf("echo Running FBAT on chromosome %s markers\n", strchr + (strchr[0] == '0' ? 1 : 0 ));
//            }
//
//            sprintf(out_fl, "%s.log", this->file_names[3]);
//            sh_rm(out_fl);
//
//#ifdef RUNSHELL_SETUP
//            sh_inline("FBAT", "cat - %s%s.cmd.txt <<EOF | $_FBAT >& /dev/null",
//                      ((LoopOverTrait && num_traits > 1) ? "../" : ""),
//                      file_names[6]);
//#else /* RUNSHELL_SETUP */
//            sh_inline("FBAT", "cat - %s%s.cmd.txt <<EOF | fbat >& /dev/null",
//                      ((LoopOverTrait && num_traits > 1) ? "../" : ""),
//                      file_names[6]);
//#endif /* RUNSHELL_SETUP */
//
//            pr_printf("log %s.log\n", file_names[3]);
////          pr_printf("log off\n");
//            if (has_x) {
//                pr_printf("load map %s%s\n", ((LoopOverTrait && num_traits > 1) ? "../" : ""), file_names[1]);
//                pr_printf("load -x ped %s\n", file_names[0]);
//            } else {
//                pr_printf("load ped %s\n", file_names[0]);
//            }
//            if (num_traits > 1) {
//                if (LoopOverTrait)
//                    pr_printf("load phe ../%s\n", file_names[2]);
//                else
//                    pr_printf("load phe %s\n", file_names[2]);
//            }
//
//            pr_printf("EOF\n");
//
//            pr_printf("echo The raw FBAT output can be found in %s/%s.log\n",
//                      *_opath, file_names[3]);
//
//            fprintf_status_check_csh(_filep, "fbat", 0);
//
//            if (has_x) {
//                pr_nl();
//                if (_numchr >= 0)
//                    CHR_STR(_numchr, out_fl);
//                else
//                    strcpy(out_fl, "all");
//                pr_printf("Rscript --vanilla %s%s.R %s%s %s %s %s\n",
//                          ((LoopOverTrait && num_traits > 1) ? "../" : ""),
//                          file_names[6],  //fbat.R
//                          ((LoopOverTrait && num_traits > 1) ? "../" : ""),
//                          file_names[1],  // map
//                          _LTop->Pheno[_trait].TraitName,  // trait
//                          file_names[6],              // stem
//                          out_fl);                    // chr
//                pr_printf("echo A table of the FBAT results can be found in %s/%s.%s.tbl\n",
//                          *_opath, file_names[6], out_fl);
//                pr_printf("echo A graph of the FBAT results can be found in %s/R%s.%s.pdf\n",
//                          *_opath, file_names[6], out_fl);
//            }
//        }
//        void file_post() {
//            chmod_X_file(path_);
//        }
//
//    } *merlin_shs = new merlin_sh(Top);
//
//    merlin_shs->file_names = file_names;
//
//    merlin_shs->sh = sh;
//
//    merlin_shs->iterate();
//
//    delete merlin_shs;
//
//    if (top_shell) {
//        mssgvf("        Merlin top shell file:   %s/%s\n", output_paths[0], file_names[4]);
//        mssgvf("             the above shell runs all shells\n");
//        sh->filep_close();
//        delete sh;
//    }


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


    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    //
    // Set the default value...
    // write everything to one file (unless only one chromosome has been selected).
    combine_chromo = main_chromocnt > 1;

    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    LoopOverChrm = ! combine_chromo;

    write_MERLIN_map(Top, file_names);

    write_MERLIN_data(Top, file_names);

    write_MERLIN_freq(Top,file_names);

    write_MERLIN_sh(Top, file_names);

    write_MERLIN_peds(Top, file_names, pwid, fwid,
                      use_map);
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
        sprintf(file_names[4], "%s.all.sh",stem, stem);
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
