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
#include "sh_util.h"
#include "loop_templates.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "write_fbat_ext.h"
#include "fbat.R.h"

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


void create_FBAT_files(linkage_ped_top **LPedTop,
                       char *file_names[],
                       const int untyped_ped_opt,
                       const int output_format);

static void save_FBAT_pheno(linkage_ped_top *Top, char *file_names[],
                            const int pwid, const int fwid);

static void save_FBAT_peds(linkage_ped_top *Top, char *file_names[],
                           const int pwid, const int fwid,
                           const bool has_x);

static void write_FBAT_sh(linkage_ped_top *Top, char *file_names[], bool has_x);

static void write_FBAT_Rhdr(linkage_ped_top *Top, char *file_names[]);

static double get_gp(ext_linkage_locus_top *EXLTop, int LType, int chr, char *snp, int m);

static void write_FBAT_map(linkage_ped_top *Top, char *file_names[]);

static void inner_file_names(char **file_names, const char *num, const char *stem = "fbat");


static void save_FBAT_pheno(linkage_ped_top *Top, char *file_names[],
                            const int pwid, const int fwid)
{
/*
    struct save_pheno_loop: public fileloop::once {
        typedef char *str;
        str *file_names;

        save_pheno_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::once(Top) { }
        void make_file() {
            msgvf("        FBAT phenotype file:   %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2]);
        }
    } *fl = new save_pheno_loop(Top);
*/

    FLPonce *floop = new FLPonce(Top, file_names[2], "w");
    floop->file_type = "        FBAT phenotype file:   ";

    struct save_pheno: public dataloop::ped_per_trait {

        save_pheno(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per_trait(Top, fl) { }
        void file_header() {
            int tr;
            int trait;
            /* print the file header */
            for (tr = 0; tr < num_traits; tr++) {
                trait = global_trait_entries[tr];
                // If this is a marker, skip it...
                if (trait == -1) continue;
                // not quant, skip it ...
//              if (_LTop->Locus[trait].Type != QUANT) continue;

                pr_printf("%s ", _LTop->Pheno[trait].TraitName);
            }
            pr_nl();
        }
        void per_start() { pr_id(); }
        void per_end()   { pr_nl(); }
        void inner() {
            /* trait locus */
            pr_pheno();
        }
    } *sp = new save_pheno(Top, floop);
    sp->load_formats(fwid, pwid, -1);

    floop->iterate();

    delete sp;
    delete floop;
}

static void save_FBAT_peds(linkage_ped_top *Top, char *file_names[],
                           const int pwid, const int fwid,
                           const bool has_x)
{
/*
    struct save_peds_loop: public fileloop::both {
        typedef char *str;
        str *file_names;

        save_peds_loop(linkage_ped_top *Top) : fileloop::both(Top) { }
        void make_file() {
            mssgvf("        FBAT pedigree file:    %s/%s\n", *_opath, file_names[0]);
            if (_ftte == NULL) {
                mssgvf("        FBAT pedigree file:    needs a trait value to be defined\nExiting\n");
                EXIT(DATA_INCONSISTENCY);
            }
            dataloop->data_loop(*_opath, file_names[0]);
        }
    } *fl = new save_peds_loop(Top);
    fl->file_names    = file_names;
    fl->_trait_affect = true;
*/
//HERE
    FLPboth *floop = new FLPboth(Top, file_names[0], "w");
    floop->file_type     = "        FBAT pedigree file:    ";
    floop->_trait_affect = true;
    floop->_needs_trait  = true;

    struct save_peds: public dataloop::ped_per_loci {
        bool has_x;

//HERE
        save_peds(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per_loci(Top, fl) { }
/*
        save_peds(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per_loci(Top, fl) {
            fileloop = fl;
            fl->dataloop = this;
        }
*/
        void file_header() {
            int locus;
            int m;
            linkage_locus_rec *lpe;
            if (has_x) return;
            for (m=0; m < NumChrLoci; m++) {
                locus = ChrLoci[m];
                lpe = &(_LTop->Locus[locus]);
                if (lpe->Class == MARKER)
                    pr_printf("%s ", _Top->LocusTop->Locus[locus].LocusName);
            }
            pr_nl();
        }
        void per_start() {
            pr_id();
            pr_parent();
            pr_sex();
            pr_aff();
        }
        void per_end() {
            pr_nl();
        }
        void inner() {
            pr_marker();
        }
    } *sp = new save_peds(Top, floop);
    sp->has_x              = has_x;
    sp->_loci_allele_limit = 40;
    sp->load_formats(fwid, pwid, -1);

    floop->iterate();

//HERE
    delete sp;
    delete floop;
}

/**
   @param output_format determines how the files are written:
*/
static void write_FBAT_sh(linkage_ped_top *Top, char *file_names[], bool has_x)
{
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
        strcmp(output_paths[0], ".");

/*
    struct all_sh: public sh_util {
        all_sh(linkage_ped_top *Top) : sh_util(Top) {}
        all_sh(linkage_ped_top *Top, fileloop::fileloop_data *fl) : sh_util(Top, fl) { }
        virtual ~all_sh() {}
        virtual void file_post() {
            chmod_X_file(path_);
        }
    }  *sh = 0;
*/
    DTshell *sh = 0;
    if (top_shell) {
        sh = new DTshell(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }

/*
    struct FBAT_sh_script_loop: public fileloop::both {
        typedef char *str;
        str *file_names;

        FBAT_sh_script_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            mssgvf("        FBAT shell file:       %s/%s\n", *_opath, file_names[3]);
            run_loop(*_opath, file_names[3]);
        }
    } *fl = new FBAT_sh_script_loop(Top);
*/

    FLPboth *floop = new FLPboth(Top, file_names[3], "w");
    floop->file_type = "        FBAT shell file:       ";
    floop->_trait_affect = true;

    struct FBAT_sh_script: public DTshell {
        typedef char *str;
        str *file_names;
        DTshell *sh;
        bool has_x;

        FBAT_sh_script(linkage_ped_top *Top, fileloop::fileloop_data *fl) : DTshell(Top, fl) { }
        void file_header() {
            if (sh)
                sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
	    // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
	    fprintf_env_checkset_csh(_filep, "_FBAT", "fbat");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            char out_fl[2*FILENAME_LENGTH];

            if (_numchr > 0) {
                char strchr[4];
                CHR_STR(_numchr, strchr);
                pr_printf("echo Running FBAT on chromosome %s markers\n", strchr + (strchr[0] == '0' ? 1 : 0 ));
            }

            sprintf(out_fl, "%s.log", this->file_names[3]);
            sh_rm(out_fl);

#ifdef RUNSHELL_SETUP
            sh_inline("FBAT", "cat - %s%s.cmd.txt <<EOF | $_FBAT >& /dev/null",
                      ((LoopOverTrait && num_traits > 1) ? "../" : ""),
                      file_names[6]);
#else /* RUNSHELL_SETUP */
            sh_inline("FBAT", "cat - %s%s.cmd.txt <<EOF | fbat >& /dev/null",
                      ((LoopOverTrait && num_traits > 1) ? "../" : ""),
                      file_names[6]);
#endif /* RUNSHELL_SETUP */

            pr_printf("log %s.log\n", file_names[3]);
//          pr_printf("log off\n");
            if (has_x) {
                pr_printf("load map %s%s\n", ((LoopOverTrait && num_traits > 1) ? "../" : ""), file_names[1]);
                pr_printf("load -x ped %s\n", file_names[0]);
            } else {
                pr_printf("load ped %s\n", file_names[0]);
            }
            if (num_traits > 1) {
                if (LoopOverTrait)
                    pr_printf("load phe ../%s\n", file_names[2]);
                else
                    pr_printf("load phe %s\n", file_names[2]);
            }

            pr_printf("EOF\n");

            pr_printf("echo The raw FBAT output can be found in %s/%s.log\n",
                      *_fileloop->_opath, file_names[3]);

            fprintf_status_check_csh(_filep, "fbat", 0);

            if (has_x) {
                pr_nl();
                if (_numchr >= 0)
                    CHR_STR(_numchr, out_fl);
                else
                    strcpy(out_fl, "all");

                pr_printf("Rscript --vanilla %s%s.R %s%s %s %s %s\n",
                          ((LoopOverTrait && num_traits > 1) ? "../" : ""),
                          file_names[6],  //fbat.R
                          ((LoopOverTrait && num_traits > 1) ? "../" : ""),
                          file_names[1],  // map
                          _LTop->Pheno[_trait].TraitName,  // trait
                          file_names[6],              // stem
                          out_fl);                    // chr
                pr_printf("echo A table of the FBAT results can be found in %s/%s.%s.tbl\n",
                          *_fileloop->_opath, file_names[6], out_fl);
                pr_printf("echo A graph of the FBAT results can be found in %s/R%s.%s.pdf\n",
                          *_fileloop->_opath, file_names[6], out_fl);
            }
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new FBAT_sh_script(Top, floop);
    xp->file_names = file_names;
    xp->sh         = sh;
    xp->has_x      = has_x;

    floop->iterate();

    delete xp;
    delete floop;

    if (top_shell) {
        mssgvf("        FBAT top shell file:   %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("             the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}

static void write_FBAT_Rhdr(linkage_ped_top *Top, char *file_names[])
{
/*
    struct FBAT_Rhdr_loop: public fileloop::both {
        typedef char *str;
        str *file_names;

        FBAT_Rhdr_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            mssgvf("        FBAT R hdr file:       %s/%s\n", *_opath, file_names[7]);
            data_loop(*_opath, file_names[7]);
        }
    } *fl = new FBAT_Rhdr_loop(Top);
*/

    FLPboth *floop = new FLPboth(Top, file_names[7], "w");
    floop->file_type = "        FBAT R hdr file:       ";
    floop->_trait_affect = true;

    struct FBAT_Rhdr: public dataloop::null {
        char strchr[4];

        FBAT_Rhdr(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner () {

            pr_printf("# reference line, y-axis minimum and y-axis maximum\n");
            pr_printf("yline <- 4.00\n");
            pr_printf(" ymin <- 0.00\n");
            pr_printf(" ymax <- 6.00\n");
            pr_printf(" # Enforce y-axis bounds?\n");
            pr_printf("yfix=FALSE\n");
            pr_printf("# Plot subtitle\n");
            if (_numchr >= 0)
                CHR_STR(_numchr, strchr);
            else
                strcpy(strchr, "all");
            pr_printf("title <- \"Chromosome %s\"\n", strchr + (strchr[0] == '0' ? 1 : 0 ));
            pr_printf("# Score units\n");
            pr_printf("ylabel <- \"-log10(P)\"\n");
            pr_printf("# Automatically adjusted by nplplot.multi.R\n");
            pr_printf("# draw.lgnd <- FALSE\n");
            pr_printf("# xlabl <- ""\n");
            pr_printf("# lgndx <- NULL\n");
            pr_printf("# lgndy <- NULL\n");
            pr_printf("# Font scaling for legend\n");
            pr_printf("cex.legend <- 0.7000\n");
            pr_printf("# Font scaling for axis\n");
            pr_printf("cex.axis <- 0.9000\n");
            pr_printf("# Tick length for marker labels at the top of plots\n");
            pr_printf("tcl <- -0.3000 \n");
            pr_printf("# NO colors\n");
            pr_printf("bw <- TRUE \n");
            pr_printf("ltypes <- 0 \n");
            pr_printf("ptypes <- 1 \n");
            pr_printf("# Remove NAs before plotting each curve\n");
            pr_printf("na.rm <- TRUE \n");
            pr_printf("#---inserted by Mega2----- \n");

        }
    } *xp = new FBAT_Rhdr(Top, floop);

    floop->iterate();

    delete xp;
    delete floop;
}

int display_error_sagm = 0, display_error_ssgm = 0, display_error_fgm = 0;
static double get_gp(ext_linkage_locus_top *EXLTop, int LType, int chr, char *snp, int m)
{
    double genetic_distance = 0.0;

    if (genetic_distance_index >= 0) {
        // The routines that build the EXLTop structure should be making sure that
        // sex maps can only be of the form <a>, <m,f>, <f>.
        if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            // we were told to only use the average map, or only an average map was specified...
            genetic_distance = EXLTop->EXLocus[m].positions[genetic_distance_index];
            // But tell the user about it if used for a sex linked chromosome...
//          if ((LType == XLINKED || LType == YLINKED) && display_error_sagm == 0)
            if ((chr == SEX_CHROMOSOME || chr == MALE_CHROMOSOME) && display_error_sagm == 0) {
                warnf("Since only a sex-averaged genetic map was used, these positions have been used");
                warnvf("for markers (%d, %s) on the X and Y chromosomes.\n", chr, snp);
                display_error_sagm++;
            }
        } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
            // Male, Female are available, so pick the right one given the current locus type...
//          if (LType == XLINKED) // X or SEX_CHROMOSOME
            if (chr == SEX_CHROMOSOME) // X or SEX_CHROMOSOME  
                genetic_distance = EXLTop->EXLocus[m].pos_female[genetic_distance_index];
//          else if (LType == YLINKED) // Y or MALE_CHROMOSOME
            else if (chr == MALE_CHROMOSOME) // Y or MALE_CHROMOSOME
                genetic_distance = EXLTop->EXLocus[m].pos_male[genetic_distance_index];
            else if (display_error_ssgm == 0) {
                // autozomes will be set to missing is a result of the decision made to allow the user
                // to select a SEX_SPECIFIC_GDMT
                warnvf("Since a sex specific genetic map was used, only analysis on the X & Y chromosome is possible.\n");
                display_error_ssgm++;
            }
        } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
//          if (LType == XLINKED) // X or SEX_CHROMOSOME
            if (chr == SEX_CHROMOSOME) // X or SEX_CHROMOSOME
                genetic_distance = EXLTop->EXLocus[m].pos_female[genetic_distance_index];
            else if (display_error_fgm == 0) {
                warnvf("Since only a female genetic map was used, only analysis on the X chromosome is possible.\n");
                display_error_fgm++;
            }
        }
    } // if (genetic_distance_index >= 0) {

    return genetic_distance;
}

/**
   @param output_format determines how the files are written:
*/
static void write_FBAT_map(linkage_ped_top *Top, char *file_names[])
{
/*
    struct FBAT_map_loop: public fileloop::chr {
        typedef char *str;
        str *file_names;

        FBAT_map_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::chr(Top) { }
        void make_file() {
            mssgvf("        FBAT map file:         %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1]);
        }
    } *fl = new FBAT_map_loop(Top);
*/

    FLPchr *floop = new FLPchr(Top, file_names[1], "w");
    floop->file_type = "        FBAT map file:         ";

    struct FBAT_map: public dataloop::loci {

        FBAT_map(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) { }
//         marker_name   chr#   genetic_pos   physical_pos   sex_link
        void inner () {
            int chr = _tle->Marker->chromosome;
            if (chr == UNKNOWN_CHROMO || chr >= MITO_CHROMOSOME) chr = 0;
            else if (chr == PSEUDO_X) chr = SEX_CHROMOSOME;

            double pp = base_pair_position_index >= 0 ?_EXLTop->EXLocus[_locus].positions[base_pair_position_index] : 0;
            if (pp < 0) pp = 0;

            double gp = get_gp(_EXLTop, _tle->Type, _tle->Marker->chromosome, _tle->LocusName, _locus);
            pr_printf("%s %d %.6f %.0f %d\n", _tle->LocusName, chr, gp, pp, _tle->Marker->chromosome == SEX_CHROMOSOME);
        }
    } *xp = new FBAT_map(Top, floop);
    xp->_loci_allele_limit = 40;

    floop->iterate();

    delete xp;
    delete floop;
}

void CLASS_FBAT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2)
{
    int pwid, fwid, mwid;
    int combine_chromo = 0;
    char prefix[100];
    linkage_ped_top *Top = LPedTreeTop;
    bool has_x = false;
    bool use_map;

/*
 * always use map
    if ((*analysis)->_suboption == 2)
        use_map = false;
    else
 */
        use_map = true;

    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    //
    // Set the default value...
    // write everything to one file (unless only one chromosome has been selected).
    combine_chromo = main_chromocnt > 1;

    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    LoopOverChrm = ! combine_chromo;

    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);

    //Tilt combine_chromo = main_chromocnt > 1;

    sh_util *Rpgm = NULL;
    Rpgm = new sh_util(Top);
    Rpgm->filep_open(output_paths[0], "fbat.R", "w");
    mssgvf("        FBAT R code file:      %s/%s\n", output_paths[0], "fbat.R");
    Rpgm->pr_printf(fbatR);
    delete Rpgm;
    Rpgm = NULL;

    sh_util *Fpgm = new sh_util(Top);
    char outfl[FILENAME_LENGTH];
    sprintf(outfl, "%s.cmd.txt", file_names[6]);
    Fpgm->filep_open(output_paths[0], outfl, "w");
    mssgvf("        FBAT commands file:    %s/%s\n", output_paths[0], outfl);
//  Fpgm->pr_printf("log on\n");
    Fpgm->pr_printf("minsize 2\n");
    Fpgm->pr_printf("fbat\n");
    Fpgm->pr_printf("quit\n");
    delete Fpgm;
    Fpgm = NULL;
    
    for (int i = 0; i < main_chromocnt; i++) {
        if (global_chromo_entries[i] == SEX_CHROMOSOME) {
            has_x = true;
            break;
        }
    }

    if (num_traits > 1) {
        save_FBAT_pheno(Top, file_names, pwid, fwid);
    }

    write_FBAT_map(Top, file_names);

    save_FBAT_peds(Top, file_names, pwid, fwid,
                   use_map || (combine_chromo && has_x));

    write_FBAT_Rhdr(Top, file_names);

    write_FBAT_sh(Top, file_names, use_map || (combine_chromo && has_x));

}

void CLASS_FBAT::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 0:
    case 1:  strcpy(subprog, "map");                    break;
    case 2:  strcpy(subprog, "no map");                 break;
    default:                                             break;
    }
}

void CLASS_FBAT::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int selection = 1;
    int selected  = 1;
    char select[10];

    if (batchANALYSIS) {
        if (Mega2BatchItems[/* 6 */ Analysis_Sub_Option].items_read) {
            selection = (*analysis)->_suboption;
        } else {
            selection = 1;
        }
    } else {
        while (selected != 0) {
            draw_line();
            printf("Selection Menu: FBAT output file options\n");
            printf("0) Done with this menu - please proceed\n");
            printf("%c1) generate FBAT map file\n", selection == 1 ? '*' : ' ');
            printf("%c2) do not generate FBAT map file\n", selection == 2 ? '*' : ' ');
            printf("Enter selection: 0 - 2 > ");
            fcmap(stdin,"%s", select); newline;
            sscanf(select, "%d", &selected);
            if (selected < 0 || selected > 2) warn_unknown(select);
            else if (selected) selection = selected;
        }
    }
    (*analysis)->_suboption = selection;
}

void CLASS_FBAT::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'm': // map
        (*analysis)->_suboption = 1; break;
    case 'n': // no map
        (*analysis)->_suboption = 2; break;
    case 0: // missing
    default:
        (*analysis)->_suboption = 1; break;
    }
}

void CLASS_FBAT::get_file_names(char *file_names[], char *prefix,
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

    sprintf(file_names[0], "%s.%s.ped", stem, num);
    sprintf(file_names[1], "%s.%s.map", stem, num);
    sprintf(file_names[2], "%s.phe", stem);
    sprintf(file_names[3], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.all.sh", stem);
    sprintf(file_names[5], "%s.%s.", stem, num);
    sprintf(file_names[6], "%s", stem);
    sprintf(file_names[7], "R%s.%s.hdr", stem, num);
}

void CLASS_FBAT::file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_FBAT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[7], numchr);
}
