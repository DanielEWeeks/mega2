/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "write_roadtrips_ext.h"

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

static void inner_file_names(char **file_names, const char *num, const char *stem = "roadtrips");


static void save_ROADTRIPS_phes(linkage_ped_top *Top, const char *prefix,
                                const int pwid, const int fwid)
{
    vlpCLASS(roadtrips_phes,trait,ped_per) {
     vlpCTOR(roadtrips_phes,trait,ped_per) { }

        void file_loop() {
            mssgvf("        ROADTRIPS phenotype file:     %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_id();
            pr_parent();
            pr_sex();
            pr_pheno();
            pr_nl();
        }
    } *sp = new roadtrips_phes(Top);

    sp->setfln(prefix, ".phes");
    
    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;
}

static void save_ROADTRIPS_peds(linkage_ped_top *Top, const char *prefix,
                                const int pwid, const int fwid)
{
    vlpCLASS(roadtrips_peds,trait,ped_per) {
     vlpCTOR(roadtrips_peds,trait,ped_per) { }

        void file_loop() {
            mssgvf("        ROADTRIPS pedigree file:      %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            linkage_ped_rec *ntpe;
            int nper;
            for (ntpe = _tpersonp, nper = _per; nper < _Top->Ped[_ped].EntryCnt; nper++, ntpe++) {
                pr_fam();
                pr_per();
                pr_per(ntpe);
                pr_printf(" 0\n");
            }
        }
    } *sp = new roadtrips_peds(Top);

    sp->setfln(prefix, ".peds");
    
    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;
}

static void save_ROADTRIPS_tpeds(linkage_ped_top *Top, const char *prefix,
                                const int pwid, const int fwid)
{
    vlpCLASS(roadtrips_tpeds,chr,loci_ped_per) {
     vlpCTOR(roadtrips_tpeds,chr,loci_ped_per) { }

        void file_loop() {
            mssgvf("        ROADTRIPS genotype file:      %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void loci_start() {
            pr_printf("%d %s ", _tlocusp->Marker->chromosome, _tlocusp->LocusName);
            pr_genetic_distance(0, 0);
            pr_physical_distance(0);
            pr_printf("  ");
        }
        void inner() {
            pr_printf("%d %d   ", _allele1, _allele2);
        }
        void loci_end() {
            pr_nl();
        }
    } *sp = new roadtrips_tpeds(Top);

    sp->setfln(prefix, ".XX.tpeds");
    
    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();
    delete sp;
}

static void save_ROADTRIPS_kins(linkage_ped_top *Top, const char *prefix,
                                const int pwid, const int fwid)
{
    vlpCLASS(roadtrips_kins,trait,ped_per) {
     vlpCTOR(roadtrips_kins,trait,ped_per) { }

        void file_loop() {
            mssgvf("        ROADTRIPS family file:        %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_id();
            pr_parent();
            pr_nl();
        }
    } *sp = new roadtrips_kins(Top);

    sp->setfln(prefix, ".fams");
    
    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;

    vlpCLASS(roadtrips_lsts,trait,ped_per) {
     vlpCTOR(roadtrips_lsts,trait,ped_per) { }

        void file_loop() {
            mssgvf("        ROADTRIPS family list file:   %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_id();
            pr_nl();
        }
    } *lp = new roadtrips_lsts(Top);

    lp->setfln(prefix, ".lst");
    
    lp->load_formats(fwid, pwid, -1);

    lp->_trait_affect = true;
    lp->iterate();

    delete lp;
}

static void save_ROADTRIPS_prvl(linkage_ped_top *Top, const char *prefix,
                                const int pwid, const int fwid,
                                double male_prevalence, double female_prevalence)
{
    vlpCLASS(roadtrips_prvl,once,null) {
     vlpCTOR(roadtrips_prvl,once,null) { }
     double male_prevalence, female_prevalence;

        void file_loop() {
            mssgvf("        ROADTRIPS prevalence file:    %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_printf("%.4f %.4f\n", male_prevalence, female_prevalence);
        }

    } *sp = new roadtrips_prvl(Top);

    sp->male_prevalence = male_prevalence;
    sp->female_prevalence = female_prevalence;
    sp->setfln(prefix, ".prvl");
    
    sp->load_formats(fwid, pwid, -1);
    sp->iterate();
    delete sp;
}

void CLASS_ROADTRIPS::create_output_file(
    linkage_ped_top *LPedTreeTop,
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
    // if 'combine_chromo == 1' all information goes into one file with the '.all' suffix.
    combine_chromo = 0;

    get_file_names(file_names, Top->OrigIds, Top->UniqueIds, &combine_chromo);
//.    combine_chromo = 0;
    LoopOverChrm   = ! combine_chromo;

    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);

    save_ROADTRIPS_phes(Top, file_name_stem, pwid, fwid);
    save_ROADTRIPS_peds(Top, file_name_stem, pwid, fwid);
    save_ROADTRIPS_tpeds(Top, file_name_stem, pwid, fwid);
    save_ROADTRIPS_kins(Top, file_name_stem, pwid, fwid);
    save_ROADTRIPS_prvl(Top, file_name_stem, pwid, fwid, male_prevalence, female_prevalence);
    write_key_file(Mega2KeysRun, Top);

}

/**
   @param output_format determines how the files are written:
*/
void CLASS_ROADTRIPS::create_sh_file(linkage_ped_top *Top, char *file_names[], const int numchr)
{

    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
        strcmp(output_paths[0], ".");

    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);

        sh->setfln(file_name_stem, ".top.sh");

        sh->filep_open(output_paths[0], sh->_fln, "w");
        sh->sh_main();
    }

    vlpCLASS(ROADTRIPS_sh_script,both,sh_exec) {
     vlpCTOR(ROADTRIPS_sh_script,both,sh_exec) {
            strcpy(pfx, (LoopOverTrait && num_traits > 1) ? "../" : "");
            do_pedigree = 1;
        }
        dataloop::sh_exec *sh;
        const char *pgm;
        const char *prefix;
        const char *additional_program_args;
        char pfx[4];
        int do_pedigree;

        void file_loop() {
            mssgvf("        ROADTRIPS shell file:         %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void file_header() {
            if (sh)
                sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
	    // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
	    fprintf_env_checkset_csh(_filep, "_ROADTRIPS", "roadtrips");
#endif /* RUNSHELL_SETUP */
        }

        void inner () {
            char cmd[2*FILENAME_LENGTH];
            char target[FILENAME_LENGTH];
            
            if (do_pedigree) {
                const char *lpgm = "KinInbcoef";
/*
 *              GEN kinship matrix
 */

                sprintf(cmd, "%s/%s", "$ROADTRIPS", "KinInbcoef");
                sh_find_pgm("ROADTRIPS", cmd, "KinInbcoef");
                sprintf(cmd, "$%s_program ", "KinInbcoef");

                mkfln(target, prefix, ".fams");
                sh_need_data(lpgm, target);
                sprintf(cmd, "%s %s", cmd, target);

                mkfln(target, prefix, ".lst");
                sh_need_data(lpgm, target);
                sprintf(cmd, "%s %s", cmd, target);

                mkfln(target, prefix, ".kins");
                sprintf(cmd, "%s %s", cmd, target);
                pr_nl();

                mkfln(target, prefix, ".log");
                sh_rm(target);
                sh_run(lpgm, cmd, target);
                fprintf_status_check_csh(_filep, "KinInbcoef", 0);
                mkfln(target, prefix, ".kins");
                sh_need_data(lpgm, target);
                pr_nl();
                pr_nl();
            }
/*
 *          Run ROADTRIPS
 */
            sprintf(cmd, "%s/%s", "$ROADTRIPS", "roadtrips");
            sh_find_pgm("ROADTRIPS", cmd, "roadtrips");
            sprintf(cmd, "$%s_program ", "roadtrips");

            mkfln(target, prefix, ".phes");
            sh_need_data(pgm, target);
            sprintf(cmd, "%s -p %s", cmd, target);

            if (do_pedigree) 
                mkfln(target, prefix, ".kins");
            else
                mkfln(target, prefix, ".peds");

            sh_need_data(pgm, target);
            sprintf(cmd, "%s -k %s", cmd, target);

            mkfln(target, pfx, prefix, ".XX.tpeds");
            sh_need_data(pgm, target);
            sprintf(cmd, "%s -g %s", cmd, target);

            mkfln(target, pfx, prefix, ".prvl");
            sh_need_data(pgm, target);
            sprintf(cmd, "%s -r %s", cmd, target);
            pr_nl();

            sprintf(cmd, "%s %s ", cmd, additional_program_args);
            pr_nl();

            mkfln(target, prefix, ".XX.log");
            sh_rm(target);
            sh_run(pgm, cmd, target);
            fprintf_status_check_csh(_filep, "roadtrips", 0);

            sh_show(pgm, target);
        }
        void sh_find_pgm(const char *NAME, const char *fullpath, const char *path) {
            pr_printf("if ( $?%s  ) then\n", NAME);
            pr_printf("  set %s_def=1\n", NAME);
            pr_printf("else\n");
            pr_printf("  set %s_def=0\n", NAME);
            pr_printf("endif\n");
            pr_printf("echo\n");
            pr_printf("if ( \"`type -t %s`\" == \"file\" ) then\n", path);
            pr_printf("  echo set %s_program=`type -p %s`\n", path, path);
            pr_printf("  set %s_program=`type -p %s`\n", path, path);
            pr_printf("else if ( $%s_def && -x \"%s\" ) then\n", NAME, fullpath);
            pr_printf("  echo set %s_program=%s\n", path, fullpath);
            pr_printf("  set %s_program=%s\n", path, fullpath);
            pr_printf("else\n");
            pr_printf("  echo The %s executable was not found - \n", fullpath);
            pr_printf("  echo please set your %s environment variable properly so %s can be found.\n", NAME, path);
            pr_printf("  echo\n");
            pr_printf("    if (! $%s_def) then\n", NAME);
            pr_printf("      echo %s is not defined.\n", NAME);
            pr_printf("    else\n");
            pr_printf("      echo %s is set to \"$%s\".\n", NAME, NAME);
            pr_printf("    endif\n");
            pr_printf("  echo\n");
            pr_printf("  echo If using Bash and ksh you would use something like this:\n");
            pr_printf("  echo export %s=dir_to_%s\n", NAME, path);
            pr_printf("  echo\n");
            pr_printf("  echo If using csh you would use something like this:\n");
            pr_printf("  echo setenv %s dir_to_%s\n", NAME, path);
            pr_printf("  echo\n");
//          pr_printf("  echo \"Be sure to run 'make %s' to build %s in the %s\"\n",
//                      pgm, pgm, path);
//          pr_printf("  echo sub directory of %s.\n", NAME);
            pr_printf("  echo\n");
            pr_printf("  echo \"For further details, please see '%s' section of the Mega2 documentation.\"\n", NAME);
            pr_printf("  exit 0\n");
            pr_printf("endif\n");
            pr_nl();
        }
    } *xp = new ROADTRIPS_sh_script(Top);

    xp->setfln(file_name_stem, ".XX.sh");
    xp->prefix = file_name_stem;
    xp->additional_program_args = C(additional_program_args);

    xp->sh         = sh;
    xp->pgm        = "ROADTRIPS";
    xp->iterate();
    delete xp;

    if (top_shell) {
        mssgvf("        ROADTRIPS top shell file:     %s/%s\n", output_paths[0], sh->_fln);
        mssgvf("                the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}

void CLASS_ROADTRIPS::get_file_names(char *file_names[], int has_orig, int has_uniq,
                                     int *combine_chromo)
{
    int i, choice, nl;
    int iarg, igl, ipre, iphen, ish, ioui, ioup, isum, isumf, iprevm, iprevf;
    int Outfile_Names = 0;
    analysis_type analysis = this;
    char selection[MAX_NAMELEN];
    char *selectionp = selection;

    batch_in();

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        choice = 0;
    } else {
        choice = -1;
    }
    *combine_chromo = 1;
    int tmp = ! *combine_chromo;
    BatchValueSet(tmp, "Loop_Over_Chromosomes");

    /* output file name menu */
    iarg = igl = ipre = iphen = ish = ioui = ioup = isum = isumf = iprevm = iprevf = -1;
    // If the default output files are used we do not go here.
    // Otherwise, enter with choice -- -1.
    while (choice != 0) {
//        draw_line();
        print_outfile_mssg();
        printf(" %s menu:\n", (_subname && (*_subname != 0)) ? _subname : _name);
        printf("0) Done with this menu - please proceed\n");
        i=1;
/*
        if (main_chromocnt > 1) {
            printf(" %d) Combine chromosomes?                      %s\n",
                   i, yorn[*combine_chromo]);
            igl=i++;
        }
*/
        printf(" %d) Filename stem:                       %-15s\n", i, file_name_stem);
        ipre=i++;

        printf(" %d) Additional ROADTRIPS arguments:      %-15s\n", i, 
               (additional_program_args.size() > 0 ? C(additional_program_args) : 
                "<none specified>"));
        iarg=i++;

        printf(" %d) Male prevalence fraction:            %-15.4f\n", i, male_prevalence);
        iprevm=i++;

        printf(" %d) Female prevalence fraction:          %-15.4f\n", i, female_prevalence);
        iprevf=i++;

        individual_id_item(i, analysis, OrigIds[0], 48, 2, 0, 0);
        ioui=i++;

        pedigree_id_item(i, analysis, OrigIds[1], 48, 2, 0);
        ioup=i;

        printf("Enter options 0-%d > ", i);
        fcmap(stdin, "%d", &choice); printf("\n");
        test_modified(choice);

        if (choice < 0) {
            printf("Unknown option %d\n", choice);
        } else if (choice == 0) {
            ;

//      } else if (choice == igl) {
//          *combine_chromo = TOGGLE(*combine_chromo);

        } else if (choice == ipre) {
            printf("Enter output filename stem > ");
            fcmap(stdin, "%s", selection);    newline;
            free(file_name_stem);
            file_name_stem = strdup(selection);
            BatchValueSet(file_name_stem, "file_name_stem");
            Outfile_Names++;
            selection[0] = 'n';
            BatchValueSet(selection[0], "Default_Outfile_Names");

        } else if (choice == iarg) {
            printf("Enter additional arguments for ROADTRIP > ");
            IgnoreValue(fgets(selection, MAX_NAMELEN-1, stdin));
            nl = strlen(selection);
            if (selection[nl-1] == '\n') selection[nl-1] = 0;
            BatchValueSet(selectionp, "additional_program_args");
            additional_program_args = selectionp;

        } else if (choice == iprevm) {
            printf("Enter male prevalence fraction > ");
            fcmap(stdin, "%g", male_prevalence);    newline;
            BatchValueSet(male_prevalence, "RoadTrips_male_prevalence");

        } else if (choice == iprevf) {
            printf("Enter female prevalence fraction > ");
            fcmap(stdin, "%g", female_prevalence);    newline;
            BatchValueSet(female_prevalence, "RoadTrips_female_prevalence");

        } else if (choice == ioui) {
            OrigIds[0] = individual_id_item(0, analysis, OrigIds[0], 35, 1, has_orig, has_uniq);
            individual_id_item(0, analysis, OrigIds[0], 0, 3, has_orig, has_uniq);

        } else if (choice == ioup) {
            OrigIds[1] = pedigree_id_item(0, analysis, OrigIds[1], 35, 1, has_orig);
            pedigree_id_item(0, analysis, OrigIds[1], 0, 3, has_orig);

        } else {
            printf("Unknown option %d\n", choice);
        }
        draw_line();
    }

    if (Outfile_Names == 0) {
        selection[0] = 'y';
        BatchValueSet(selection[0], "Default_Outfile_Names");
    }

    if (! batchANALYSIS) batch_out();
    batch_show();
}


/*
static keyw_t keywords[] = {
    {"RoadTrips_male_prevalence",                      STRING,     ""},
    {"RoadTrips_female_prevalence",                    STRING,     ""},
};
*/

void CLASS_ROADTRIPS::batch_out()
{
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "additional_program_args",
                       "RoadTrips_male_prevalence",
                       "RoadTrips_female_prevalence",
                       "Loop_Over_Chromosomes",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_ROADTRIPS::batch_in()
{
    char *fn = this->file_name_stem;

    BatchValueIfSet(                fn,   "file_name_stem");
    
    BatchValueGet(additional_program_args, "additional_program_args");
    BatchValueGet(male_prevalence,         "RoadTrips_male_prevalence");
    BatchValueGet(female_prevalence,       "RoadTrips_female_prevalence");

    BatchValueGet(LoopOverChrm,            "Loop_Over_Chromosomes");
}

void CLASS_ROADTRIPS::batch_show()
{
    msgvf("\n");
    if (! DEFAULT_OUTFILES) {
        msgvf("Output file stem:                         %s\n",    C(file_name_stem));
    }
    msgvf("Additional ROADTRIPS program args:       %s\n",
          (additional_program_args.size() > 0 ? C(additional_program_args) : 
            "<none specified>"));

    msgvf("ROADTRIPS Male Prevalence:                %.4f\n",   male_prevalence);
    msgvf("ROADTRIPS Female Prevalence:              %.4f\n",   female_prevalence);

    msgvf("\n");
 }

static void inner_file_names(char **file_names, const char *num, const char *stem /* = "roadtrips" */) {

    sprintf(file_names[0], "%s", stem);
    sprintf(file_names[1], "%s.%s", stem, num);
}

void CLASS_ROADTRIPS::gen_file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_ROADTRIPS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[1], numchr);
}
