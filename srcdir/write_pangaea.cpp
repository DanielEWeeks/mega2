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
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "write_pangaea_ext.h"

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

static void inner_file_names(char **file_names, const char *num, const char *stem = "pangaea");


static void save_PANGAEA_peds(linkage_ped_top *Top, char *file_names[],
                              const int pwid, const int fwid, const int subopt)
{
/*
    struct save_peds_loop: public fileloop::trait {
        save_peds_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA pedigree file:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0]);
        }
    } *floop = new save_peds_loop(Top);
*/

    FLPtrait *floop = new FLPtrait(Top, file_names[0], "w");
    floop->file_type = "        PANGAEA pedigree file:      ";

    struct save_peds: public dataloop::ped_per {
        int subopt;
        int skip;

        save_peds(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }
        void file_header() {
            int nquant = 0, nint = 0, col = 5;

            if (_Top->IndivCnt > 20000)
                pr_printf("allow pedigree size %d\n", _Top->IndivCnt);
            pr_printf("input pedigree size %d\n", _Top->IndivCnt);
            pr_printf("input pedigree record father mother\n");
            pr_printf("input pedigree record gender present\n");
            if (_tte == 0) {
            } else if (_tte->Type == AFFECTION && _tte->Pheno->Props.Affection.ClassCnt > 1) {
                pr_printf("# file columns %d,%d (integer pair) correspond to trait name %s\n",
                          col, col+1,  _tte->LocusName);
                col += 2;
                nint += 2;
            } else if (_tte->Type == AFFECTION) {
                pr_printf("# file column %d (integer) corresponds to trait name %s\n",
                          col,  _tte->LocusName);
                col++;
                nint++;
            } else if (_tte->Type == QUANT) {
                pr_printf("# file column %d (real) corresponds to trait name %s\n",
                          nquant+1+nint+1,  _tte->LocusName);
                nquant++;
            }
            pr_printf("input pedigree record names 3 integers %d", nint + 1 /* sex */);
            if (nquant)
                pr_printf(" reals %d", nquant);
            pr_nl();

            pr_printf("****************************************\n");
        }
        void ped_start() {
            if (subopt == 6) {
/*
                if (_tp->EntryCnt > 50 || _ped > 100) {
                    warnvf("Skipping family %s; too many families for analysis\n",_tp->Name);
                    skip = 1;
                } else
*/
                    skip = 0;
            } else
                skip = 0;
        }
        void inner() {
            if (skip) return;
            pr_per();
            pr_parent();
            pr_sex();
            if (_tte == 0) {
            } else if (_tte->Type == AFFECTION) {
/*
                int ase;
                ase = aff_status_entry(_tpe->Data[_trait].Affection.Status,
                                       _tpe->Data[_trait].Affection.Class,
                                       &(_Top->LocusTop->Locus[_trait]));
                pr_printf("%1d   ", ase);
*/
                pr_printf("%1d   ", _tpe->Pheno[_trait].Affection.Status);
                if (_tte->Pheno->Props.Affection.ClassCnt > 1)
                    pr_printf("%1d   ", _tpe->Pheno[_trait].Affection.Class);
            } else if (_tte->Type == QUANT) {
                if (fabs(_tpe->Pheno[_trait].Quant - MissingQuant) <= EPSILON) {
                    pr_printf("%10.5f ",  999.0);
                } else {
                    pr_printf("%10.5f ", _tpe->Pheno[_trait].Quant);
                }
            } else {
                pr_printf("%d ", 0);
            }
            pr_nl();
        }
    } *sp = new save_peds(Top, floop);
    sp->subopt     = subopt;
    sp->load_formats(fwid, pwid, -1);

    floop->iterate();
    delete sp;
    delete floop;

/*
    struct save_peds_all_loop: public fileloop::once {
        save_peds_all_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::once(Top) { }
        void make_file() {
            mssgvf("        PANGAEA pedigree file:      %s/%s\n", *_opath, file_names[9]);
            data_loop(*_opath, file_names[9]);
        }
    } *floop = new save_peds_all_loop(Top);
*/

    FLPonce *floop1 = new FLPonce(Top, file_names[9], "w");
    floop1->file_type = "        PANGAEA pedigree file:      ";

    struct save_peds_all: public dataloop::ped_per_trait {

        save_peds_all(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per_trait(Top, fl) { }
        void file_header() {
            int i, j;
            int nquant = 0, nint = 0;

            if (_Top->IndivCnt > 20000)
                pr_printf("allow pedigree size %d\n", _Top->IndivCnt);
            pr_printf("input pedigree size %d\n", _Top->IndivCnt);
            pr_printf("input pedigree record father mother\n");
            pr_printf("input pedigree record gender present\n");
            int col = 5;
            for (i = 0; i < num_traits; i++) {
                j = global_trait_entries[i];
                if (j == -1)  continue;
                linkage_locus_rec *loc = &_LTop->Locus[j];
                if (loc->Type == AFFECTION && loc->Pheno->Props.Affection.ClassCnt > 1) {
                    pr_printf("# file columns %d,%d (integer pair) correspond to trait name %s\n",
                              col, col+1,  _LTop->Locus[j].LocusName);
                    col += 2;
                    nint += 2;
                } else if (loc->Type == AFFECTION) {
                    pr_printf("# file column %d (integer) corresponds to trait name %s\n",
                              col,  _LTop->Locus[j].LocusName);
                    col++;
                    nint++;
                } else if (loc->Type == QUANT) {
                    pr_printf("# file column %d (real) corresponds to trait name %s\n",
                              nquant+1+nint+1,  _LTop->Locus[j].LocusName);
                    nquant++;
                }
            }
            pr_printf("input pedigree record names 3 integers %d", nint + 1 /* sex */);
            if (nquant)
                pr_printf(" reals %d", nquant);
            pr_nl();

            pr_printf("****************************************\n");
        }
        void per_start() {
            pr_per();
            pr_parent();
            pr_sex();
        }
        void per_end() {
            pr_nl();
        }
        void inner() {
            if (_tte == 0) {
            } else if (_tte->Type == AFFECTION) {
/*
                int ase;
                ase = aff_status_entry(_tpe->Data[_trait].Affection.Status,
                                       _tpe->Data[_trait].Affection.Class,
                                       &(_Top->LocusTop->Locus[_trait]));
                pr_printf("%1d   ", ase);
*/
                pr_printf("%1d   ", _tpe->Pheno[_trait].Affection.Status);
                if (_tte->Pheno->Props.Affection.ClassCnt > 1)
                    pr_printf("%1d   ", _tpe->Pheno[_trait].Affection.Class);
            } else if (_tte->Type == QUANT) {
                if (fabs(_tpe->Pheno[_trait].Quant - MissingQuant) <= EPSILON) {
                    pr_printf("%10.5f ",  999.0);
                } else {
                    pr_printf("%10.5f ", _tpe->Pheno[_trait].Quant);
                }
            } else {
                pr_printf("%d ", 0);
            }
        }
        void sort_quant_last() // first names; then integer values; then real values
        {
            int i = 0, j = 0;
            int *straits = CALLOC(num_traits, int);
            for (i = 0; i < num_traits; i++)
                if (_LTop->Locus[i].Type != QUANT)
                    straits[j++] = global_trait_entries[i];
            for (i = 0; i < num_traits; i++)
                if (_LTop->Locus[i].Type == QUANT)
                    straits[j++] = global_trait_entries[i];
            for (i = 0; i < num_traits; i++)
                global_trait_entries[i] = straits[i];
            free(straits);
        }
    } ;
/*
   *sp1 = new save_peds_all(Top, floop1);

    sp1->load_formats(fwid, pwid, -1);
    sp1->sort_quant_last();
    floop1->iterate();
    delete sp1;
    delete floop1;
*/

}

/**
   @param output_format determines how the files are written:
*/
const char *pgm_path[] = {"",
                          "PedComp/pedcheck",    //1
                          "PedComp/kin",         //2
                          "PedComp/translink",   //3
                          "Lodscore/lm_linkage", //4
                          "Lodscore/lm_bayes",   //5
                          "Autozyg/lm_ibdtests", //6
                          "Autozyg/lm_ibdtests", //7
};

static void write_PANGAEA_sh(linkage_ped_top *Top, char *file_names[], char *pgm, int subopt)
{
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
        strcmp(output_paths[0], ".");
/*
    struct all_sh: public sh_util {
        all_sh(linkage_ped_top *Top) : sh_util(Top) {}
//      virtual ~all_sh() {}
        virtual void file_post() {
            chmod_X_file(path_);
        }
    } *sh = 0;
*/
    DTshell *sh = 0;
    if (top_shell) {
        sh = new DTshell(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }

/*
    struct PANGAEA_sh_script_loop: public fileloop::both {
        PANGAEA_sh_script_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            mssgvf("        PANGAEA shell file:         %s/%s\n", *_opath, file_names[3]);
            run_loop(*_opath, file_names[3]);
        }
    } *floop = new PANGAEA_sh_script_loop(Top);
*/

    FLPboth *floop = new FLPboth(Top, file_names[3], "w");
    floop->file_type = "        PANGAEA shell file:         ";

    struct PANGAEA_sh_script: public DTshell {
        typedef char *str;
        str *file_names;
        DTshell *sh;
        char *pgm;
        int subopt;
        char pfx[4];

        PANGAEA_sh_script(linkage_ped_top *Top, fileloop::fileloop_data *fl) : DTshell(Top, fl) {
            strcpy(pfx, (LoopOverTrait && num_traits > 1) ? "../" : "");
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
	    fprintf_env_checkset_csh(_filep, "_PANGAEA", "pangaea");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];
            char target[FILENAME_LENGTH];

            sprintf(target, "%spar_chr_trt", file_names[5]);
            sh_need_data(pgm, target);

            sprintf(cmd, "cat %spar_chr_trt %s.par_trt",
                    file_names[5], pgm);

            sprintf(target, "%spar", file_names[5]);
            sh_cat(cmd, target);

            sprintf(cmd, "%s/%s", "$MORGAN_RELEASE", pgm_path[subopt]);
            sh_find_pgm(cmd, pgm_path[subopt]);

            sprintf(cmd+strlen(cmd), " %s", target);
            strcat(target, ".log");
            sh_rm(target);
            sh_run(pgm, cmd, target);
            sh_show(pgm, target);

            fprintf_status_check_csh(_filep, "pangaea", 0);

            if (subopt == 3) { // translink
                sprintf(target, "datafile%s", file_names[11]);
                sh_save_output(pgm, "datafile.dat", target);
                sprintf(target, "pedfile%s", file_names[11]);
                sh_save_output(pgm, "pedfile.dat", target);
            }
        }
        void sh_find_pgm(const char *fullpath, const char *path) {
            pr_printf("if ( $?MORGAN_RELEASE  ) then\n");
            pr_printf("  set morgan_def=1\n");
            pr_printf("else\n");
            pr_printf("  set morgan_def=0\n");
            pr_printf("  set MORGAN_RELEASE=MORGAN_RELEASE\n");
            pr_printf("endif\n");
            pr_printf("echo\n");
            pr_printf("if (! $morgan_def  || ! -x \"%s\" ) then\n", fullpath);
            pr_printf("  echo The %s executable was not found - \n", fullpath);
            pr_printf("  echo please set your MORGAN_RELEASE environment variable properly so %s can be found.\n", pgm);
            pr_printf("  echo\n");
            pr_printf("    if (! $morgan_def) then\n");
            pr_printf("      echo MORGAN_RELEASE is not defined.\n");
            pr_printf("    else\n");
            pr_printf("      echo MORGAN_RELEASE is $MORGAN_RELEASE.\n");
            pr_printf("    endif\n");
            pr_printf("  echo\n");
            pr_printf("  echo If using Bash and ksh you would use something like this:\n");
            pr_printf("  echo export MORGAN_RELEASE=path_to/MORGAN_V311_RELEASE\n");
            pr_printf("  echo\n");
            pr_printf("  echo If using csh you would use something like this:\n");
            pr_printf("  echo setenv MORGAN_RELEASE path_to/MORGAN_V311_RELEASE\n");
            pr_printf("  echo\n");
            pr_printf("  echo \"Be sure to run 'make %s' to build %s in the %s\"\n",
                      pgm, pgm, path);
            pr_printf("  echo sub directory of MORGAN_RELEASE.\n");
            pr_printf("  echo\n");
            pr_printf("  echo \"For further details, please see 'PANGAEA/MORGAN' section of the Mega2 documentation.\"\n");
            pr_printf("  exit 0\n");
            pr_printf("endif\n");
            pr_nl();
        }
    } *xp = new PANGAEA_sh_script(Top, floop);
    xp->file_names = file_names;
    xp->sh         = sh;
    xp->pgm        = pgm;
    xp->subopt     = subopt;

    floop->iterate();

    delete xp;
    delete floop;

    if (top_shell) {
        mssgvf("        PANGAEA top shell file:     %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("                the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}

typedef FLPboth par_var_loop;

/*
    struct par_var_loop: public fileloop::both {
        par_var_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
    void make_file() {
        mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_names[7]);
        data_loop(*_opath, file_names[7]);
    }
    } *floop = new par_var_loop(Top);
*/

struct par_var: public dataloop::null {
    typedef char *str;
    str *file_names;
    char pfx[4];

    par_var(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) {
        strcpy(pfx, (LoopOverTrait && num_traits > 1) ? "../" : "");
    }

    void inner () {
        pr_printf("input pedigree file '%s'\n", file_names[0]);
        pr_printf("output overwrite pedigree file '%s.out'\n", file_names[0]);
        pr_printf("# chromosome data does not depend on traits so is in the above directory\n");
        pr_printf("input marker data file '%s%s'\n", pfx, file_names[1]);
        pr_nl();
        if (_tte != 0) {
            if (_tte->Type == AFFECTION && _tte->Pheno->Props.Affection.ClassCnt > 1) {
                pr_printf("input extra file '%s.liability.extra'\n", _tte->LocusName);
                pr_printf("input pedigree record trait 1 integer pairs 2 3\n");
                pr_printf("set trait 1 data discrete with liability\n");
            } else if (_tte->Type == AFFECTION) {
                pr_printf("input pedigree record trait 1 integer 2\n");
                pr_printf("set trait 1 data discrete\n");
            } else if (_tte->Type == QUANT) {
                pr_printf("input pedigree record trait 1 real 1\n");
                pr_printf("set trait 1 data quant\n");
            }
        }
    }
};

struct liability_traits_loop: public fileloop::trait {
    liability_traits_loop(linkage_ped_top *Top) : fileloop::trait(Top) { }

    void make_file() {
        if (_ftte == 0 || _ftte->Type != AFFECTION || _ftte->Pheno->Props.Affection.ClassCnt <= 1)
            return;
        char outfl[FILENAME_LENGTH];
        sprintf(outfl, "%s.liability.extra", _ftte->LocusName);
        mssgvf("        PANGAEA liability file:     %s/%s\n", *_opath, outfl);
        _dataloop->data_loop(*_opath, outfl);
    }
};

struct liability_traits: public dataloop::null {
    liability_traits(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }

    void inner() {
        int i = 0;
        int cnt  = _tte->Pheno->Props.Affection.ClassCnt;
        linkage_affection_class *pen = _tte->Pheno->Props.Affection.Class; //.{Male/Female/Auto}Pen[0]@3 MaleDef

        for (i = 1; i <= cnt; i++) {
            pr_printf("%d %.6f %.6f %.6f\n",
                      i,
                      pen->AutoPen[0],
                      pen->AutoPen[1],
                      pen->AutoPen[2]);
            pen++;
        }
    }
};

/**
   @param output_format determines how the files are written:
*/
static void write_PANGAEA_par_template(linkage_ped_top *Top, char *file_names[], char *pgm)
{
    par_var_loop *floop1 = new par_var_loop(Top, file_names[7], "w");
    floop1->file_type    = "        PANGAEA chr var  par file:  ";

    par_var *xp1 = new par_var(Top, floop1);
    xp1->file_names = file_names;
    floop1->iterate();
    delete xp1;
    delete floop1;

/*
    struct template_par_user_loop: public fileloop::trait {
        template_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new template_par_user_loop(Top);
*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct template_par_user: public dataloop::null {

        template_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner() {
            pr_printf("set printlevel 5\n");
            pr_printf("output pedigree chronological\n");
            pr_printf("output pedigree record father mother\n");

            pr_printf("set printlevel 4\n");
            pr_printf("select all markers\n");
            pr_printf("select trait 1\n");
            pr_printf("set trait 1 data discrete\n");
            pr_printf("set traits 1 for tlocs 1 incomplete penetrances .05 .90 .90\n");

            pr_printf("set traits 1 tlocs 1\n");
            pr_printf("set tloc 1 allele frequency .09 .91\n");
            pr_printf("map tloc 1 marker 00 recom frac .0005\n");

            pr_printf("use single meiosis sampler\n");
            pr_printf("set MC iterations 1\n");
            pr_printf("set proband gametes 202 0 202 1\n");
        }
    } *xp2 = new template_par_user(Top, floop2);
    floop2->iterate(); 
    delete xp2;
    delete floop2;
}

static void write_PANGAEA_par_pedcheck(linkage_ped_top *Top, char *file_names[], char *pgm)
{

/*
    struct pedcheck_xx_par_var_loop: public fileloop::both {
        pedcheck_xx_par_var_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_names[7]);
            data_loop(*_opath, file_names[7]);
        }
    } *floop = new pedcheck_xx_par_var_loop(Top);
*/

    FLPboth *floop1 = new FLPboth(Top, file_names[7], "w");
    floop1->file_type = "        PANGAEA chr var  par file:  ";

    struct pedcheck_xx_par_var: public dataloop::null {
        typedef char *str;
        str *file_names;

        pedcheck_xx_par_var(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner () {
            pr_printf("input pedigree file '%s'\n", file_names[0]);
            pr_printf("output overwrite pedigree file '%s.out'\n", file_names[0]);
        }
    } *xp1 = new pedcheck_xx_par_var(Top, floop1);
    xp1->file_names     = file_names;
    floop1->iterate();
    delete xp1;
    delete floop1;

/*
    struct pedcheck_par_user_loop: public fileloop::trait {
        pedcheck_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new pedcheck_par_user_loop(Top);
*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct pedcheck_par_user: public dataloop::null {
        pedcheck_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner() {
            pr_printf("set printlevel 5\n");
            pr_printf("output pedigree chronological\n");
            pr_printf("output pedigree record father mother\n");

        }
    } *xp2 = new pedcheck_par_user(Top, floop2);
    floop2->iterate(); 
    delete xp2;
    delete floop2;
}

static void write_PANGAEA_par_kin(linkage_ped_top *Top, char *file_names[], char *pgm)
{
/*
    struct kin_xx_par_varA_loop: public fileloop::both {
        kin_xx_par_varA_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            i = 0;
            mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_names[7]);
            data_loop(*_opath, file_names[7]);
        }
    } *floop = new kin_xx_par_varA_loop(Top);
*/

    FLPboth *floop1A = new FLPboth(Top, file_names[7], "w");
    floop1A->file_type = "        PANGAEA chr var  par file:  ";

    struct kin_xx_par_varA: public dataloop::ped_per {
        typedef char *str;
        str *file_names;
        int i;

        kin_xx_par_varA(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }
        void file_header () {
            i = 0;
            pr_printf("input pedigree file '%s'\n", file_names[0]);
            pr_printf("output overwrite pedigree file '%s.out'\n", file_names[0]);
            pr_nl();
        }
        void ped_start() {
            ++i;
            if (_tp->EntryCnt == 1) pr_printf("# ");
            pr_printf("compute component %d kinship coeff ", i);
        }
        void inner() {
            int nper;
            for (nper = _per+1; nper < _Top->Ped[_ped].EntryCnt; nper++) {
                pr_per();
                pr_per(&_tp->Entry[nper]);
            }
        }
        void ped_end() {
            pr_nl();
        }
        void file_trailer() {
            pr_nl();
        }
    } *xp1A = new kin_xx_par_varA(Top, floop1A);
    xp1A->file_names = file_names;
    xp1A->load_formats(5, 5, -1);   //?

    floop1A->iterate();

    delete xp1A;
    delete floop1A;

/*
    struct kin_xx_par_varB_loop: public fileloop::both {
        kin_xx_par_varB_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            i = 0;
            mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_names[7]);
            data_loop(*_opath, file_names[7], "a");
        }
    } *floop = new kin_xx_par_varB_loop(Top);
*/

    FLPboth *floop1B = new FLPboth(Top, file_names[7], "a");
    floop1B->file_type = "        PANGAEA chr var  par file:  ";

    struct kin_xx_par_varB: public dataloop::ped_per {
        int i;

        kin_xx_par_varB(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }
        void file_header() {
            i = 0;
        }
        void ped_start() {
            pr_printf("compute component %d inbreeding coeff ", ++i);
        }
        void inner() {
            pr_per();
        }
        void ped_end() {
            pr_nl();
        }
        void file_trailer() {
            pr_nl();
        }
    } *xp1B = new kin_xx_par_varB(Top, floop1B);
    xp1B->load_formats(5, 5, -1);   //?

    floop1B->iterate();

    delete xp1B;
    delete floop1B;

/*
    struct kin_xx_par_varC_loop: public fileloop::both {
        kin_xx_par_varC_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::both(Top) { }
        void make_file() {
            i = 0;
            mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_names[7]);
            data_loop(*_opath, file_names[7], "a");
        }
    } *floop = new kin_xx_par_varC_loop(Top);
*/

    FLPboth *floop1C = new FLPboth(Top, file_names[7], "a");
    floop1C->file_type = "        PANGAEA chr var  par file:  ";

    struct kin_xx_par_varC: public dataloop::ped_per {
        int i;

        kin_xx_par_varC(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }
        void file_header() {
            i = 0;
        }
        void ped_start() {
            pr_printf("compute component %d two-locus inbreeding coeff ", ++i);
        }
        void inner() {
            pr_per();
        }
        void ped_end() {
            pr_nl();
        }
        void file_trailer() {
            pr_nl();
        }
    } *xp1C = new kin_xx_par_varC(Top, floop1C);
    xp1C->load_formats(5, 5, -1);   //?

    floop1C->iterate();

    delete xp1C;
    delete floop1C;

/*
    struct kin_par_user_loop: public fileloop::trait {
        kin_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new kin_par_user_loop(Top);
*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct kin_par_user: public dataloop::null {

        kin_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner() {
            pr_printf("set printlevel 5\n");
            pr_printf("# output pedigree chronological\n");

        }
    } *xp2 = new kin_par_user(Top, floop2);

    floop2->iterate(); 

    delete xp2;
    delete floop2;
}

static void write_PANGAEA_par_translink(linkage_ped_top *Top, char *file_names[], char *pgm)
{

    struct translink_xx_par_var_loop: public FLPboth {
        translink_xx_par_var_loop(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
            FLPboth(Top, f_name, f_mode) {}

//HERE
        void make_file() {
            if (_ftte == 0) {
            } else if (_ftte->Type == AFFECTION && _ftte->Pheno->Props.Affection.ClassCnt > 1) {
                warnvf("%s: Trait with Liability Class is not supported.  Trait will be ignored.\n", _ftte->LocusName);
            } else if (_ftte->Type == QUANT) {
                warnvf("%s: QUANT type is not supported.  Trait will be ignored.\n", _ftte->LocusName);
            }
            mssgvf("        PANGAEA chr var  par file:  %s/%s\n", *_opath, file_name);
            _dataloop->data_loop(*_opath, file_name);
        }
    } *floop1 = new translink_xx_par_var_loop(Top, file_names[7], "w");
//  FLPboth *floop = new FLPboth(Top);

    struct translink_xx_par_var: public dataloop::null {
        typedef char *str;
        str *file_names;
        char pfx[4];

        translink_xx_par_var(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) {
            strcpy(pfx, (LoopOverTrait && num_traits > 1) ? "../" : "");
        }
        void inner () {
            if (_tte == 0) {
                return;
            } else if (_tte->Type == AFFECTION && _tte->Pheno->Props.Affection.ClassCnt > 1) {
                return; // leave empty file
            } else if (_tte->Type == QUANT) {
                return;
            }
            pr_printf("input pedigree file '%s'\n", file_names[0]);
            pr_printf("output overwrite pedigree file '%s.out'\n", file_names[0]);
            pr_printf("# chromosome data does not depend on traits so is in the above directory\n");
            pr_printf("input marker data file '%s%s'\n", pfx, file_names[1]);
            if (_tte != 0) {
                if (_tte->Type == AFFECTION && _tte->Pheno->Props.Affection.ClassCnt > 1) {
                    pr_printf("# liability classes are not supported for translink.\n");
                    pr_printf("# input extra file '%s.liability.extra'\n", _tte->LocusName);
                    pr_printf("input pedigree record trait 1 integer pairs 2 3\n");
                    pr_printf("set trait 1 data discrete with liability\n");
                } else if (_tte->Type == AFFECTION) {
                    pr_printf("input pedigree record trait 1 integer 2\n");
                    pr_printf("set trait 1 data discrete\n");
                } else if (_tte->Type == QUANT) {
                    pr_printf("input pedigree record trait 1 real 1\n");
                    pr_printf("set trait 1 data quant\n");
                }
            }
        }
    } *xp1 = new translink_xx_par_var(Top, floop1);

    xp1->file_names     = file_names;

    floop1->iterate();

    delete xp1;
    delete floop1;

/*
    struct translink_par_user_loop: public fileloop::trait {
        translink_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new translink_par_user_loop(Top);
*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct translink_par_user: public dataloop::null {
        typedef char *str;
        str *file_names;
        char pfx[4];

        translink_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) {
            strcpy(pfx, (LoopOverTrait && num_traits > 1) ? "../" : "");
        }
        void inner() {
            pr_printf("input extra file '%s%s.extra'\n", pfx, file_names[6]);
            pr_printf("set printlevel 4\n");
            pr_printf("select all markers\n");

            pr_printf("select trait 1\n");
            pr_printf("set traits 1 for tlocs 1 incomplete penetrances .05 .90 .90\n");

            pr_printf("set traits 1 tlocs 1\n");
            pr_printf("set tloc 1 allele frequency .09 .91\n");
            pr_printf("map tloc 1 marker 00 recom frac .0005\n");

            pr_printf("use single meiosis sampler\n");
            pr_printf("set MC iterations 1\n");
            pr_printf("set proband gametes 202 0 202 1\n");
        }
    } *xp2 = new translink_par_user(Top, floop2);

    xp2->file_names     = file_names;

    floop2->iterate();

    delete xp2;
    delete floop2;

    sh_util *Xpgm = new sh_util(Top);
    char outfl[FILENAME_LENGTH];
    sprintf(outfl, "%s.extra", file_names[6]);
    Xpgm->filep_open(output_paths[0], outfl, "w");
    mssgvf("        translink extra file:       %s/%s\n", output_paths[0], outfl);
    Xpgm->pr_printf(" 1 1 1 10\n\
# first item: numnam; if 1, output using MORGAN index as name.\n\
#                     if 0, output using MORGAN names\n\
# second item: siblinks; if 1, first-kid, next_patsib, next_matsib links\n\
#                 will be set from Marlist and output ; if 0, not\n\
# third  item: dexit; effective only if siblinks=1;\n\
#                     if 1, output by Marlist order; \n\
#                     if 0, use MORGAN index order\n\
# fourth item:   the number of lod score evaluations\n\
#     The info on the interval\n\
#           is now taken from the  \"map tloc\" statement \n\
#     In the case of first or last interval, the initial 0.5 or 0,\n\
#          and the final value is 0.0 or 0.5\n\
");
    delete Xpgm;
    Xpgm = NULL;
}

static void write_PANGAEA_par_lod(linkage_ped_top *Top, char *file_names[], char *pgm, int subopt)
{
    par_var_loop *floop1 = new par_var_loop(Top, file_names[7], "w");
    floop1->file_type    = "        PANGAEA chr var  par file:  ";

    par_var *xp1 = new par_var(Top, floop1);
    xp1->file_names = file_names;
    floop1->iterate();

    delete xp1;
    delete floop1;

/*
    struct lod_par_user_loop: public fileloop::trait {
        lod_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new lod_par_user_loop(Top);
, file_names[7], "w"*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct lod_par_user: public dataloop::null {
        int subopt;

        lod_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner() {
            pr_printf("#For default seeds comment out the following lines\n");
            pr_printf("#To save seeds remove comments\n");
            pr_printf("#input seed file         'seedfil'\n");
            pr_printf("#output seed file        'seedfil'\n");
            pr_printf("set printlevel 3\n");

            pr_printf("select all markers\n");
            pr_printf("select trait 1\n");

            pr_printf("set traits 1 for tlocs 11 incomplete penetrances .05 .90 .90\n");
            pr_printf("set trait 1 tloc 11\n");
            pr_printf("map test tloc 11 all interval proportions .1 .3 .5 .7 .9\n");
            pr_printf("map test tloc 11 external recomb fracts .05 .2 .3 .4 .45\n");
            pr_printf("set tloc 11 allele freqs  .95 .05\n");

            if (subopt == 4 || subopt == 5) { // lm_linkage and lm_bayes common
                pr_printf("use single meiosis sampler\n");
                pr_printf("sample by scan\n");
                pr_printf("set L-sampler probability 0.2\n");
                pr_printf("set MC iterations 3000\n");
                pr_printf("set burn-in iterations 150\n");
                pr_printf("check progress MC iterations 500\n");
            } 
            if (subopt == 4) { //lm_bayes
                pr_printf("use locus-by-locus sampling for setup\n");
            } else if (subopt == 5) { //lm_bayes
                pr_printf("use sequential imputation for setup\n");
                pr_printf("use 1000 sequential imputation realizations for setup\n");
                pr_printf("set sequential imputation proposals every 100 iterations\n");
                pr_printf("set burn-in iterations 500\n");
                pr_printf("set pseudo-prior iterations 3000\n");
                pr_printf("compute scores every 10 iterations\n\n");
            }

        }
    } *xp2 = new lod_par_user(Top, floop2);
    xp2->subopt         = subopt;

    floop2->iterate(); 
    delete xp2;
    delete floop2;

    liability_traits_loop *floop3 = new liability_traits_loop(Top);

    liability_traits *xp3 = new liability_traits(Top, floop3);
    floop3->iterate();

    delete xp3;
    delete floop3;
}

static void write_PANGAEA_par_ibd_tests(linkage_ped_top *Top, char *file_names[], char *pgm, int subopt)
{
    par_var_loop *floop1 = new par_var_loop(Top, file_names[7], "w");
    floop1->file_type    = "        PANGAEA chr var  par file:  ";

    par_var *xp1 = new par_var(Top, floop1);
    xp1->file_names = file_names;

    floop1->iterate();
    delete xp1;
    delete floop1;

/*
    struct lod_par_user_loop: public fileloop::trait {
        lod_par_user_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::trait(Top) { }
        void make_file() {
            mssgvf("        PANGAEA user par file:      %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10]);
        }
    } *floop = new lod_par_user_loop(Top);
*/

    FLPtrait *floop2 = new FLPtrait(Top, file_names[10], "w");
    floop2->file_type = "        PANGAEA user par file:      ";

    struct lod_par_user: public dataloop::null {
        int subopt;

        lod_par_user(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::null(Top, fl) { }
        void inner() {
            pr_printf("#For default seeds comment out the following lines\n");
            pr_printf("#To save seeds remove comments\n");
            pr_printf("#output permutat seeds only\n");
            pr_printf("#output sampler seeds only\n");
            pr_printf("#input seed file         'seedfil'\n");
            pr_printf("#output seed file        'seedfil'\n");
            pr_printf("set printlevel 4\n\n");

            pr_printf("select all markers\n");
            pr_printf("select trait 1\n\n");

            pr_printf("sample by scan\n");
            pr_printf("set L-sampler probability 0.5\n");
            pr_printf("set MC iterations 3000\n");
            pr_printf("#set burn-in iterations 150\n");
            pr_printf("#check progress MC iterations 500\n");
            pr_printf("compute scores every 100 iterations\n");
            if (subopt == 6) { //lm_ibdtests
                pr_printf("compute ibd statistics\n");
                pr_printf("set ibd measures Spairs Srobdom\n");
                pr_printf("set ibd tests norm permu\n");
                pr_printf("set ibd permutations 999\n");
                if (_Top->PedCnt > 100) {
                    pr_printf("\n");
                    pr_printf(" comp COMMENT OUT THIS LINE AND THE NEXT IF YOU WANT TO RUN LM_IBDTESTS\n");
                    pr_printf(" comp BUT BE WARNED WITH THESE OPTIONS LM_IBDTEST USES MANY GIGABYTES OF MEMORY\n");
                }
            } else if (subopt == 7) { //lm_ibdtests_lr (lr == likelyhood ratio)
                pr_printf("compute likelihood-ratio stats\n");
                pr_printf("set like measures delta lambda-p\n");
                pr_printf("set like lambda-p model grid 6 9\n");
            }
        }
    } *xp2 = new lod_par_user(Top, floop2);
    xp2->subopt         = subopt;

    floop2->iterate(); 

    delete xp2;
    delete floop2;

    liability_traits_loop *floop3 = new liability_traits_loop(Top);

    liability_traits *xp3 = new liability_traits(Top, floop3);
    floop3->iterate();

    delete xp3;
    delete floop3;
}

static double get_gp(ext_linkage_locus_top *EXLTop, int LType, int chr, char *snp, int m, int choice)
{
    double genetic_distance = 0.0;

    if (genetic_distance_index >= 0) {
        if (choice == 1) {
            // we were told to only use the average map, or only an average map was specified...
            genetic_distance = EXLTop->EXLocus[m].positions[genetic_distance_index];
        } else if (choice == 2) {
            genetic_distance = EXLTop->EXLocus[m].pos_female[genetic_distance_index];
        } else if (choice == 3) {
            genetic_distance = EXLTop->EXLocus[m].pos_male[genetic_distance_index];
        }
    }

    return genetic_distance;
}

/**
   @param output_format determines how the files are written:
*/
static void write_PANGAEA_map(linkage_ped_top *Top, char *file_names[], int subopt)
{

/*
    struct PANGAEA_map_names_loop: public fileloop::chr {
        PANGAEA_map_names_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::chr(Top) { }
        void make_file() {
            mssgvf("        PANGAEA names map file:     %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1]);
        }
    } *floop = new PANGAEA_map_names_loop(Top);
*/

    FLPchr *floop = new FLPchr(Top, file_names[1], "w");
    floop->file_type = "        PANGAEA names map file:     ";

    struct PANGAEA_map_names: public dataloop::loci {
        int token;

        PANGAEA_map_names(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) {  }
        void file_header() {
            pr_printf("set marker names ");
            token = 0;
        }
        void inner() {
            if (token++ >= 100) { pr_nl(); token = 1; }
            pr_printf("%s ", _tle->LocusName);
        }
        void file_trailer() {
            pr_nl();
        }

    } *xp = new PANGAEA_map_names(Top, floop);

    floop->iterate(); 

    delete xp;
    delete floop;

/*
    struct PANGAEA_map_dist_loop: public fileloop::chr {
        PANGAEA_map_dist_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::chr(Top) { }
        void make_file() {
            ogp = 0 - delta;
            mssgvf("        PANGAEA dist map file:      %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "a");
        }
    } *floop = new PANGAEA_map_dist_loop(Top);
*/

    FLPchr *floop1 = new FLPchr(Top, file_names[1], "a");
    floop1->file_type = "        PANGAEA dist map file:      ";

    struct PANGAEA_map_dist: public dataloop::loci {
        double ogp;
        char *name;
        int choice;
        int token;
        double delta;

        PANGAEA_map_dist(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) {
            delta = .000001;
        }
        void file_header() {
            pr_printf("map %s marker Kosambi positions ", name);
            token = 0;
            ogp = 0 - delta;
        }
        void inner () {
            int chr = _tle->Marker->chromosome;
            if (chr == UNKNOWN_CHROMO || chr >= MITO_CHROMOSOME) chr = 0;

            double gp = get_gp(_EXLTop, _tle->Type, _tle->Marker->chromosome, _tle->LocusName, _locus, choice);
//          gp = kosambi_to_haldane(gp);
//          pr_printf(" <%.6f, %.6f> ", gp, ogp);
            if (ogp > gp - delta) gp = ogp + delta;
            ogp = gp;
            if (token++ >= 100) { pr_nl(); token = 1; }
            pr_printf("%.6f ", gp);
        }
        void file_trailer() {
            pr_nl();
        }
    };

    if (genetic_distance_sex_type_map == SEX_AVERAGED_MAP) {
        PANGAEA_map_dist *xp1 = new PANGAEA_map_dist(Top, floop1);
        xp1->name       = (char *)"";
        xp1->choice     = 1;
        floop1->iterate();
        delete xp1;
        delete floop1;
    } else {
        PANGAEA_map_dist *xp1 = new PANGAEA_map_dist(Top, floop1);
        xp1->name       = (char *)"gender f";
        xp1->choice     = 2;
        floop1->iterate();
        delete xp1;
        delete floop1;

        floop1 = new FLPchr(Top, file_names[1], "a");
        floop1->file_type = "        PANGAEA dist map file:      ";

        xp1 = new PANGAEA_map_dist(Top, floop1);
        xp1->name       = (char *)"gender m";
        xp1->choice     = 3;
        floop1->iterate();
        delete xp1;
        delete floop1;
    }
/*
    struct PANGAEA_map_freq_loop: public fileloop::chr {
        PANGAEA_map_freq_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::chr(Top) { }
        void make_file() {
            mssgvf("        PANGAEA freq map file:      %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "a");
        }
    } *floop = new PANGAEA_map_freq_loop(Top);
*/

    FLPchr *floop2 = new FLPchr(Top, file_names[1], "a");
    floop2->file_type = "        PANGAEA freq map file:      ";

    struct PANGAEA_map_freq: public dataloop::loci {
        int i;
        int token;

        PANGAEA_map_freq(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) { i = 0; }
//HERE        void chr_start() { i = 0; }
        void file_header() { i = 0; }
        void inner() {
            int j;
            pr_printf("set marker %d allele freqs ", ++i);
            linkage_allele_rec *allele_info = _tle->Allele;
            for (j = 0; j < _tle->AlleleCnt; j++) {
                pr_printf("%.4f ", allele_info[j].Frequency);
            }
            pr_nl();
        }
    } *xp2 = new PANGAEA_map_freq(Top, floop2);

    floop2->iterate(); 

    delete xp2;
    delete floop2;

/*
    struct PANGAEA_map_data_loop: public fileloop::chr {
        PANGAEA_map_data_loop(linkage_ped_top *Top, fileloop::fileloop_data *dl) : fileloop::chr(Top) { }
        void make_file() {
            mssgvf("        PANGAEA data map file:      %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "a");
        }
    } *floop = new PANGAEA_map_data_loop(Top);
*/

    FLPchr *floop3 = new FLPchr(Top, file_names[1], "a");
    floop3->file_type = "        PANGAEA data map file:      ";

    struct PANGAEA_map_data: public dataloop::ped_per_loci {
        int token;
        int subopt;
        int skip;

        PANGAEA_map_data(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per_loci(Top, fl) { }
        void file_header() { 
            int i, j = NumChrLoci;
            for (i = 0; i < num_traits; i++) { // num_traits as a count includes space for -1 and -99
                if (global_trait_entries[i] != -1 && global_trait_entries[i] != -99) j--;
            }
            pr_printf("set MARKERS %d data\n", j);  // true count of markers
            token = 0;
        }
        void ped_start() {
            if (subopt == 6) {
/*
                if (_tp->EntryCnt > 50 || _ped > 100) {
                    skip = 1;
                } else
*/
                    skip = 0;
            } else
                skip = 0;
        }
        void per_start() {
            if (skip) return;
            pr_per();
        }
        void inner () {
            if (skip) return;
            if (token++ >= 100) { pr_nl(); token = 1; }
            pr_marker();
        }
        void per_end() {
            if (skip) return;
            pr_nl();
        }
    } *xp3 = new PANGAEA_map_data(Top, floop3);
    xp3->subopt         = subopt;
    xp3->load_formats(5, 5, -1);   //?

    floop3->iterate();

    delete xp3;
    delete floop3;
}

void CLASS_PANGAEA::create_output_file(
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
    char prefix[100];

    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    combine_chromo = 0;

    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    combine_chromo = 0;
    LoopOverChrm   = 1;

    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);


    save_PANGAEA_peds(Top, file_names, pwid, fwid, _suboption);

    if (_suboption > 2) {
        // need trait ...
        if (! have_trait_b(Top->LocusTop, /*affect=*/true, /*quant=*/false)) {
            char sub_prog[32];
            sub_prog_name(_suboption, sub_prog);
            errorvf("Morgan option %s requires that a trait be specified.\n", sub_prog);
            EXIT(DATA_INCONSISTENCY);
        }
        write_PANGAEA_map(Top, file_names, _suboption);
    }

    sub_prog_name(_suboption, prefix);
    sprintf(file_names[3], "%s..sh", prefix);
    sprintf(file_names[4], "%s.all.sh", prefix);
    sprintf(file_names[5], "%s..", prefix);
    sprintf(file_names[6], "%s", prefix);
    sprintf(file_names[7], "%s..par_chr_trt", prefix);
    sprintf(file_names[8], "%s..par_user", prefix);
    sprintf(file_names[9], "%s.par_var", prefix);
    sprintf(file_names[10], "%s.par_trt", prefix);
    sprintf(file_names[11], "..dat");

    write_PANGAEA_sh(Top, file_names, prefix, _suboption);

    switch (_suboption) {
    case 0:
    case 1:  write_PANGAEA_par_pedcheck(Top, file_names, prefix);   break;
    case 2:  write_PANGAEA_par_kin(Top, file_names, prefix);   break;
    case 3:  write_PANGAEA_par_translink(Top, file_names, prefix);   break;
    case 4:  write_PANGAEA_par_lod(Top, file_names, prefix, _suboption);   break;
    case 5:  write_PANGAEA_par_lod(Top, file_names, prefix, _suboption);   break;
    case 6:  write_PANGAEA_par_ibd_tests(Top, file_names, prefix, _suboption);   break;
    case 7:  write_PANGAEA_par_ibd_tests(Top, file_names, prefix, _suboption);   break;
    default:
        break;

/* so compiler won't complain, "use" template */
    case -99:  write_PANGAEA_par_template(Top, file_names, prefix);   break;
    }
}

void CLASS_PANGAEA::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 0:
    case 1:  strcpy(subprog, "pedcheck");              break;
    case 2:  strcpy(subprog, "kin");                   break;
    case 3:  strcpy(subprog, "translink");             break;
    case 4:  strcpy(subprog, "lm_linkage");            break;
    case 5:  strcpy(subprog, "lm_bayes");              break;
    case 6:  strcpy(subprog, "lm_ibdtests");           break;
    case 7:  strcpy(subprog, "lm_ibdtests_lr");        break;
    default:                                           break;
    }
}

void CLASS_PANGAEA::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
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
            printf("Selection Menu: PANGAEA output file options\n");
            printf("0) Done with this menu - please proceed\n");
            printf("%c1) generate PANGAEA pedcheck Output files\n",
                   selection == 1 ? '*' : ' ');
            printf("%c2) generate PANGAEA kinship Output files\n",
                   selection == 2 ? '*' : ' ');
            printf("%c3) generate PANGAEA translink Output files\n",
                   selection == 3 ? '*' : ' ');
            printf("%c4) generate PANGAEA LOD lm_linkage Output files\n",
                   selection == 4 ? '*' : ' ');
            printf("%c5) generate PANGAEA LOD lm_bayes Output files\n",
                   selection == 5 ? '*' : ' ');
            printf("%c6) generate PANGAEA lm_ibdtests Output files\n",
                   selection == 6 ? '*' : ' ');
            printf("%c7) generate PANGAEA lm_ibdtests_lr Output files\n",
                   selection == 7 ? '*' : ' ');
            printf("Enter selection: 0 - 7 > ");
            fcmap(stdin,"%s", select); newline;
            sscanf(select, "%d", &selected);
            if (selected < 0 || selected > 7) warn_unknown(select);
            else if (selected) selection = selected;
        }
    }
    (*analysis)->_suboption = selection;
}

void CLASS_PANGAEA::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'p': // pedcheck
        (*analysis)->_suboption = 1; break;
    case 'k': // kinship
        (*analysis)->_suboption = 2; break;
    case 't': // translink
        (*analysis)->_suboption = 3; break;
    case 0: // missing
//      (*analysis)->_suboption = 1; break;
        interactive_sub_prog_name_to_sub_option(analysis); break;

    default:
        switch(tolower((unsigned char)sub_prog_name[3])) {
        case 'l': // lm_linkage
            (*analysis)->_suboption = 4; break;
        case 'b': // lm_bayes
            (*analysis)->_suboption = 5; break;
        case 'i': // lm_ibdtests
        default:
            switch(strlen(sub_prog_name)) {
            case 11:
                (*analysis)->_suboption = 6; break;
            case 14:
                (*analysis)->_suboption = 7; break;
            }
        }
        break;
    }
}

void CLASS_PANGAEA::get_file_names(char *file_names[], char *prefix,
                                   int has_orig, int has_uniq, int *combine_chromo)
{
    int i, choice;
    int igl, ipre, iphen, ish, ioui, ioup, isum, isumf;
    analysis_type analysis = this;

    strcpy(prefix, "pangaea");

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

        printf(" %d) File name stem:                           %-15s\n", i, prefix);
        ipre=i++;

        individual_id_item(i, analysis, OrigIds[0], 43, 2, 0, 0);
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

static void inner_file_names(char **file_names, const char *num, const char *stem /* = "pangaea" */) {

    sprintf(file_names[0], "%s.ped", stem);
    sprintf(file_names[1], "%s.%s.map", stem, num);
    sprintf(file_names[3], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.all.sh", stem);
    sprintf(file_names[5], "%s.%s.", stem, num);
    sprintf(file_names[6], "%s", stem);
    sprintf(file_names[7], "%s.%s.par_chr_trt", stem, num);
    sprintf(file_names[8], "%s.%s.par_user", stem, num);
    sprintf(file_names[9], "%s.all.ped", stem);
    sprintf(file_names[10], "%s.par_trt", stem);
    sprintf(file_names[11], ".%s.dat", num);
}

void CLASS_PANGAEA::file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_PANGAEA::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[7], numchr);
    change_output_chr(file_names[8], numchr);
    change_output_chr(file_names[11], numchr);
}
