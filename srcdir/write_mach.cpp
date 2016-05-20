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

#include "write_mach_ext.h"

void create_MACH_files(linkage_ped_top **LPedTop, char *file_names[], const int untyped_ped_opt, const int output_format);

static void write_MACH_peds(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MACH_data(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MACH_sh(linkage_ped_top *Top, char *file_names[]);

//general function to output options for mach and then parse them.
static void mach_option_menu();


static void inner_file_names(char **file_names, const char *num, const char *stem = "mach");


// Output like this
//  FAM1001   ID1234  0   0   M  A A   A C   C C
//  FAM1002   ID1234  0   0   F  A C   C C   G G
//  or
//  FAM1001   ID1234  0   0   M   1 1   1 2   2 2
//  FAM1002   ID5678  0   0   F   1 2   2 2   3 3
static void write_MACH_peds(linkage_ped_top *Top, char **file_names, const int pwid, const int fwid)
{

    vlpCLASS(mach_ped, chr, ped_per_loci) {
        vlpCTOR(mach_ped, chr, ped_per_loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Mach Pedigree File:   %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }

        void per_start() {
            //pr_id prints both family id and person id looks like
            pr_id();
            pr_father();
            pr_mother();
            pr_sex_l();
        }

        void per_end() {
            pr_nl();
        }

        void inner() {
            pr_marker();
        }

    } *mach_peds = new mach_ped(Top);

    mach_peds->file_names = file_names;

    mach_peds->load_formats(fwid, pwid, -1);

    mach_peds->iterate();

    delete mach_peds;
}


static void write_MACH_data(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
{

/*<Example of a simple data file>
 *  M marker1
 *  M marker2
 *  ...
 *  <End of simple data file>
 *  */

    vlpCLASS(mach_dat,chr,loci) {
        vlpCTOR(mach_dat,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Mach Data File:   %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "w");
        }

        void inner() {
            pr_printf("M ");
            //seems to not pull up the marker name and print it out
            //pr_marker_name();

            //this works though
            char *markername = _tlocusp->Marker->MarkerName;
            pr_printf(markername);
            pr_nl();
        }

    } *mach_dats = new mach_dat(Top);

    mach_dats->file_names = file_names;

    mach_dats->load_formats(fwid, pwid, -1);

    mach_dats->iterate();

    delete mach_dats;
}

/*Shell should look like the following in order to produce prephased VCF according to minimac documentation
 *
 * mach1 -d Gwas.chr20.Unphased.dat \
      -p Gwas.chr20.Unphased.ped \
      --rounds 20 \
      --states 200 \
      --phase \
      --interim 5 \
      --sample 5 \
      --prefix Gwas.Chr20.Phased.Output
      */
static void write_MACH_sh(linkage_ped_top *Top, char *file_names[]) {

    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
                    strcmp(output_paths[0], ".");

    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names[3], "w");
        sh->sh_main();
    }

    vlpCLASS(mach_sh, both, sh_exec) {
        vlpCTOR(mach_sh, both, sh_exec) { }

        void file_loop() {
            mssgvf("        MaCH shell file:       %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
        }

        typedef char *str;
        str *file_names;
        sh_exec *sh;
        bool has_x;

        void file_header() {
            //if (sh)
            //    sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
            // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
	    fprintf_env_checkset_csh(_filep, "_MACH", "mach");
#endif /* RUNSHELL_SETUP */
        }
        void inner() {
            pr_printf ("mach1 -d %s -p %s --rounds 20 --states 200 --phase --interim 5 --sample 5 --prefix Chr%d.Phased.Output", file_names[1], file_names[0],_numchr);
        }

    } *mach_shs = new mach_sh(Top);

    mach_shs->file_names = file_names;

    mach_shs->iterate();

    delete mach_shs;
}

void CLASS_MACH::create_output_file(
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


    combine_chromo = main_chromocnt > 1;

    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    LoopOverChrm = ! combine_chromo;

    //adding in looping over traits for the data file
    LoopOverTrait = num_traits != 1;

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    write_MACH_data(Top, file_names, pwid, fwid);

    write_MACH_peds(Top, file_names, pwid, fwid);

    write_MACH_sh(Top, file_names);


}

//this will create and parse the options
void mach_option_menu(){

}


//there should be a new way to do this.
void CLASS_MACH::get_file_names(char *file_names[], char *prefix,
                                     int has_orig, int has_uniq, int *combine_chromo)
{
    int i, choice;
    char fl_stat[12];
    int igl, ipre, iphen, ish, ioui, ioup, isum, isumf;
    analysis_type analysis = this;

    strcpy(prefix, "mach");

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

    //this should reflect the original naming convention of mach in mega2
    sprintf(file_names[0], "%s_ped.%s", stem, num);
    sprintf(file_names[1], "%s_data.%s", stem, num);
    sprintf(file_names[2], "%s.%s.sh", stem, num);
    sprintf(file_names[3], "%s.all.sh", stem);
}

void CLASS_MACH::gen_file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_MACH::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
}

