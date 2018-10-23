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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA .

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

===========================================================================
*/

// for centos ...
#define __STDC_LIMIT_MACROS 1

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"
#include "types.hh"
#include "loop.h"
#include "sh_util.h"
#include "batch_input.h"
#include "utils_ext.h"

#include "fcmap_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "write_files_ext.h"

#include "write_MQLS_ext.h"

void CLASS_MQLS::create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    linkage_ped_top *Top = LPedTreeTop;

    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all information goes into one file with the '.all' suffix.
    int combine_chromo = 0;
    LoopOverChrm   = ! combine_chromo;

    this->prevalencefilename = strdup("prevalence.txt");

    if ( InputMode == INTERACTIVE_INPUTMODE ) {
        option_menu(file_names,file_name_stem, &combine_chromo, Top);
    }
    else {
        batch_in();
    }

    hasXdata=false;
    for (int i = 0; i < main_chromocnt; i++) {
        if(global_chromo_entries[i]==23)
            hasXdata=true;
    }

    LoopOverTrait = 0;

    omit_peds(untyped_ped_opt, Top);
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    inner_file_names(file_names,"",file_name_stem, &combine_chromo);

    printf("Mega2 created the following file(s) for MQLS-XM/KinInbcoef Format:\n");

    write_MQLS_genofile(Top, file_name_stem, file_names, pwid, fwid);
    write_MQLS_phenofile(Top, file_name_stem, file_names, pwid, fwid);
    write_KinInbcoef_ped(Top, file_name_stem, file_names, pwid, fwid);
    if(hasXdata)
        write_KinInbcoefX_ped(Top, file_name_stem, file_names, pwid, fwid);
    write_MQLS_listfile(Top, file_name_stem, file_names, pwid, fwid);
    write_IQLS_shell_script(Top, file_name_stem, file_names);

}

void CLASS_MQLS::option_menu(char **file_names, char *prefix, int *combine_chromo, linkage_ped_top *Top) {
    int i, done, choice, istem, iprev;
    i = 1;
    istem = 0;
    iprev = 0;
    choice = -1;
    done = 0;


    while(choice != 0){
        draw_line();
        printf("MQLS-XM/KinInbcoef Analysis Menu:\n");
        printf(" 0) Done with this menu - please proceed\n");

        printf(" %d) Filenames stem:                                   \"%s\"\n", i++, this->file_name_stem);
        istem=i;
        printf(" %d) Prevalence file name:                             \"%s\"\n", i++, this->prevalencefilename);
        iprev=i;

        printf("Enter options 0-%d > ", i-1);
        fcmap(stdin, "%d", &choice); printf("\n");
        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }
        else if(choice == done){
            BatchValueSet(file_name_stem,"file_name_stem");
            BatchValueSet(this->prevalencefilename, "prevalence_file");
            //force seperated chromosomes
            *combine_chromo = 1;
            int tmp = (! *combine_chromo) ? 'y' : 'n';
            BatchValueSet(tmp, "Loop_Over_Chromosomes");

        }
        else if(choice == istem) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", this->file_name_stem);
            newline;
        }
        else if(choice == iprev){
            printf("Enter prevalence file name > ");
            fcmap(stdin, "%s", this->prevalencefilename);
            newline;
        }
        else {
            printf("Unknown option %d\n", choice);
        }
    }

}

//MQLS-XM genofile, it is effectively a tped file
void CLASS_MQLS::write_MQLS_genofile(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    vlpCLASS(mqls_tpeds,chr,loci_ped_per) {
        vlpCTOR(mqls_tpeds,chr,loci_ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        MQLS genotype file:                  %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
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
    } *sp = new mqls_tpeds(Top);

    sp->file_names = file_names;

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();
    delete sp;
}

//this will write a PLINK ped file which is used by MQLS-XM
//MQLS_XM refers to this as the phenofile
//the OrigID flags should set the values for fam and per appropriately
void CLASS_MQLS::write_MQLS_phenofile(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    vlpCLASS(MQLS_peds,trait,ped_per) {
        vlpCTOR(MQLS_peds,trait,ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        MQLS phenotype file:                 %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_father();
            pr_mother();
            pr_sex();
            pr_pheno();
            pr_nl();
        }
    } *sp = new MQLS_peds(Top);

    sp->file_names = file_names;

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();
    delete sp;
}

//creates a 4 column pedigree with no sex values or phenotypes
//this file is used by KinInbcoef
void CLASS_MQLS::write_KinInbcoef_ped(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    vlpCLASS(MQLS_peds2,trait,ped_per) {
        vlpCTOR(MQLS_peds2,trait,ped_per) { }
        typedef char *str;
        str *file_names;
        void file_loop() {
            mssgvf("        KinInbcoef pedigree file:            %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_father();
            pr_mother();
            pr_nl();
        }
    } *sp = new MQLS_peds2(Top);

    sp->file_names = file_names;

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();

    delete sp;
}

//creates a 5 column pedigree with sex values but no phenotypes
//this file is used by KinInbcoefX
void CLASS_MQLS::write_KinInbcoefX_ped(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    vlpCLASS(MQLS_peds3,trait,ped_per) {
        vlpCTOR(MQLS_peds3,trait,ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        KinInbcoefX pedigree file:           %s/%s\n", *_opath, file_names[3]);
            data_loop(*_opath, file_names[3], "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_father();
            pr_mother();
            pr_sex();
            pr_nl();
        }
    } *sp = new MQLS_peds3(Top);

    sp->file_names = file_names;

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();

    delete sp;
}

//this creates a list of all pairs within a family
//there is a loop over a pedigree within the ped per loop so each person loops over all others they haven't been listed with
void CLASS_MQLS::write_MQLS_listfile(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    vlpCLASS(MQLS_list, trait, ped_per){
        vlpCTOR(MQLS_list, trait,ped_per){ }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        MQLS list file:                      %s/%s\n", *_opath, file_names[4]);
            data_loop(*_opath, file_names[4], "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_nl();
            //linkage_ped_rec *ntpe;
            //int nper;
            //for (ntpe = _tpersonp, nper = _per; nper < _Top->Ped[_ped].EntryCnt; nper++, ntpe++) {
            //    pr_per();
            //    pr_per(ntpe);
            //    pr_nl();
            //}
        }
    } *sp = new MQLS_list(Top);

    sp->file_names = file_names;

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();

    delete sp;
}

void CLASS_MQLS::write_IQLS_shell_script(linkage_ped_top *Top, const char *prefix, char **file_names) {
    int top_shell = 1;

    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names[5], "w");
        sh->sh_main();
    }

    vlpCLASS(mqls_sh, both, sh_exec) {
        vlpCTOR(mqls_sh, both, sh_exec) { }

        void file_loop() {
            mssgvf("        MQLS-XM/KinInbcoef Shell File:       %s/%s\n", *_opath, file_names[6]);
            data_loop(*_opath, file_names[6], "w");
        }

        typedef char *str;
        str *file_names;
        sh_exec *sh;
        bool has_x;

        void file_header() {
            if (sh)
                sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
            // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
	    fprintf_env_checkset_csh(_filep, "_MQLS", "mqls");
#endif /* RUNSHELL_SETUP */
        }
        void inner() {
            char cmd1[2*FILENAME_LENGTH];
            char cmd2[2*FILENAME_LENGTH];
            char cmd3[2*FILENAME_LENGTH];

            //check our applications are available by running sh_find_pgm
            //only Kininbcoef or KininbcoefX is valid in one file so we don't need to check or use both
            if(_numchr != 23) {
                sprintf(cmd1, "%s/%s", "KININBCOEF", "KinInbcoef");
                sh_find_pgm("KININBCOEF", cmd1, "KinInbcoef");
                sprintf(cmd1, "$%s_program ", "KinInbcoef");
            }

            else {
                sprintf(cmd2, "%s/%s", "KININBCOEFX", "KinInbcoefX");
                sh_find_pgm("KININBCOEFX", cmd2, "KinInbcoefX");
                sprintf(cmd2, "$%s_program ", "KinInbcoefX");
            }

            //always use MQLS-XM
            sprintf(cmd3, "%s/%s", "MQLS_XM", "MQLS");
            sh_find_pgm("MQLS_XM", cmd3, "MQLS");
            sprintf(cmd3, "$%s_program ", "MQLS");

            pr_nl();
            if(_numchr != 23) {
                pr_printf("#use KinInbcoef on output to create to calculate autosomal breeding coeficients\n");
                // Try ./KinInbcoef pedtest listtest out
                //we want to add a 0 if under 10
                if(_numchr < 10)
                pr_printf("%s %s %s %s0%d%s\n", cmd1, file_names[2]/*.kininbcoef*/, file_names[4]/*.list*/,
                          "kininbcoef.", _numchr, ".out");
                else
                    pr_printf("%s %s %s %s%d%s\n", cmd1, file_names[2]/*.kininbcoef*/, file_names[4]/*.list*/,
                              "kininbcoef.", _numchr, ".out");
            }
            else {
                pr_printf("#use KinInbcoefX on output to create to calculate X-Chromosome breeding coeficients\n");
                // ./KinInbcoefX pedtestX listtestX out error
                pr_printf("%s %s %s %sX%s %s\n", cmd2, file_names[3]/*.kininbcoefx*/, file_names[4]/*.list*/,
                          "kininbcoefx.",".out", "kininbcoefx.error");
            }

            pr_nl();
            pr_printf("#use MQLS-XM to perform single-SNP, case-control association testing on the autosomal chromosomes and the X-chromosome \n");
            // ./MQLS-XM -g genofile -p phenofile -k kinfile -r prevalence -x -u -m -h
            //kinfile is the output from KinInbcoef/X
            //if chromsome 23 we need the -x flag otherwise no
            if(_numchr == 23)
                pr_printf("%s -g %s -p %s -k %sX%s -r %s -x -u -m -h\n",cmd3, file_names[0]/*.gen*/, file_names[1]/*.fam*/,"kininbcoef.", ".out", file_names[7]/*prevalencefilename*/);
            //add a zero if under 10
            else if(_numchr < 10)
                pr_printf("%s -g %s -p %s -k %s0%d%s -r %s -u -m -h\n",cmd3, file_names[0]/*.gen*/, file_names[1]/*.fam*/,"kininbcoef.", _numchr, ".out", file_names[7]/*prevalencefilename*/);
            else
                pr_printf("%s -g %s -p %s -k %s%d%s -r %s -u -m -h\n",cmd3, file_names[0]/*.gen*/, file_names[1]/*.fam*/,"kininbcoef.", _numchr, ".out", file_names[7]/*prevalencefilename*/);

            pr_printf("exit 0\n");

            pr_nl();
        }

        //finds the program to run dynamically and gives an error if it can't be found.
        void sh_find_pgm(const char *NAME, const char *fullpath, const char *path) {
            pr_printf("if ( $?%s  ) then\n", NAME);
            pr_printf("  set %s_def=\"%s\"\n", NAME, fullpath);
            pr_printf("else\n");
            pr_printf("  set %s_def=0\n", NAME);
            pr_printf("endif\n");
            pr_printf("echo\n");
            pr_printf("if ( \"`type -t %s`\" == \"file\" ) then\n", path);
            pr_printf("  echo set %s_program=`type -p %s`\n", path, path);
            pr_printf("  set %s_program=`type -p %s`\n", path, path);
            pr_printf("else if ( \"$%s_def\" != \"0\" && -x \"$%s_def\" ) then\n", NAME, NAME);
            pr_printf("  echo set %s_program=\"$%s_def\"\n", path, NAME);
            pr_printf("  set %s_program=\"$%s_def\"\n", path, NAME);
            pr_printf("else\n");
            pr_printf("  echo The %s executable was not found - \n", fullpath);
            pr_printf("  echo please set your %s environment variable properly so %s can be found.\n", NAME, path);
            pr_printf("  echo\n");
            pr_printf("    if (\"$%s_def\" == \"0\") then\n", NAME);
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
    } *mqls_shs = new mqls_sh(Top);

    mqls_shs->file_names = file_names;

    mqls_shs->sh  = sh;
    mqls_shs->iterate();

    if (top_shell) {
        mssgvf("        MQLS-XM/Kininbcoef Top Shell File:   %s/%s\n", output_paths[0], file_names[5]);
        mssgvf("        The above shell runs all shells.\n");

        sh->filep_close();
        delete sh;
    }

    delete mqls_shs;
}

void CLASS_MQLS::batch_in() {
    char c;
    char *fn = this->file_name_stem;
    BatchValueIfSet(fn,   "file_name_stem");
    BatchValueGet(c,"Loop_Over_Chromosomes");
    BatchValueGet(prevalencefilename, "prevalence_file");
    LoopOverChrm = (c == 'y' || c == 'Y');
}

void CLASS_MQLS::batch_out() {
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "prevalence_file",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_MQLS::batch_show() {
    msgvf("\n");
    msgvf("Output file stem:                         %s\n",    C(file_name_stem));
    msgvf("Prevalence file:                          %s\n",    C(prevalencefilename));
    msgvf("\n");
}

void CLASS_MQLS::inner_file_names(char **file_names, const char *num, const char *stem /* = "MQLS" */, int *combine_chromo) {
    sprintf(file_names[0], "%s.%s.gen", stem, num);
    sprintf(file_names[1], "%s.fam", stem);
    sprintf(file_names[2], "%s.kininbcoef", stem);
    sprintf(file_names[3], "%s.kininbcoefX", stem);
    sprintf(file_names[4], "%s.txt", stem);
    sprintf(file_names[5], "%s.top.sh", stem);
    sprintf(file_names[6], "%s.%s.sh", stem, num);
    sprintf(file_names[7], "%s", prevalencefilename);
}

void CLASS_MQLS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[6], numchr);
}


