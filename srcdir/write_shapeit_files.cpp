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

#include <set>
#include <sstream>
#include <string>

#include "common.h"
#include "typedefs.h"

#include "loop.h"
#include "sh_util.h"

#include "batch_input_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "str_utils.hh"
#include "write_shapeit_ext.h"
#include "write_shapeit_args.h"

// I don't believe shapeit cares about affection status

/*
     create_summary_ext.h:  aff_status_entry marker_typing_summary
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_loci_on_chromosomes get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp summary_time_stamp
*/

/*

Implementation Notes...


For PLINK documentation see: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

This describes column 6 of the .PED/.FAM file:
A PED file must have 1 and only 1 phenotype in the sixth column. The phenotype can be either a
quantitative trait or an affection status column: PLINK will automatically detect which type
(i.e. based on whether a value other than 0, 1, 2 or the missing genotype code is observed).

*/

void CLASS_SHAPEIT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, "shapeit", analysis);
}


void CLASS_SHAPEIT::user_queries(char **file_names_array,
                                 int *combine_chromo, int *create_summary)
{
    int i, choice = -1, istem = -1;
    int idir = -1, ifile = -1;
    char selection[100], *sp = selection;

    *combine_chromo = 0;

    while (choice != 0) {
        print_outfile_mssg();
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        i=1;

        printf(" %d) Specify genetic recombination map directory?       \"%s\"\n",
               i, BatchItemGet("Shapeit_recomb_rdir")->value.name);
        idir=i++;

        printf(" %d) Specify genetic recombination map file name?       \"%s\"\n",
               i, BatchItemGet("Shapeit_recomb_rfile")->value.name);
        ifile=i++;

        printf(" %d) Change file names stem?                            \"%s\"\n",
               i, BatchItemGet("Shapeit_file_stem")->value.name);
//               i, this->file_name_stem);
        istem=i++;

// specify recomb map

        printf("Enter options 0-%d > ", i-1);
        fcmap(stdin, "%d", &choice); printf("\n");

        if (choice == 0) {
            ; // OK...

        } else if (choice == idir) {
            printf("Enter directory for genetic map with recombination rate > ");
            IgnoreValue(fgets(selection, sizeof(selection)-1, stdin)); newline;
            i = strlen(selection) - 1;
            if (selection[i] == '\n')
                selection[i] = 0;
            BatchValueSet(sp, "Shapeit_recomb_rdir");
            rdir = sp;

        } else if (choice == ifile) {
            while (1) {
                printf("Enter recombination map file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");
                IgnoreValue(fgets(selection, sizeof(selection)-1, stdin)); newline;
                i = strlen(selection) - 1;
                if (selection[i] == '\n')
                    selection[i] = 0;
                BatchValueSet(sp, "Shapeit_recomb_rfile");

                Cstr file(sp);
                Vecs filesplit;
                split(filesplit, file, "?");
                if (filesplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }
                rpre  = filesplit[0];
                rpost = filesplit[1];
                break;
            }
        } else if (choice == istem) {
            char *fn = this->file_name_stem;
            printf("Enter new stem for the output file names > ");
            fcmap(stdin, "%s", this->file_name_stem);    newline;
            // It doesn't matter what the parameter 'num' in the method file_names() is. It will get
            // changed to the appropriate thing later in the code. The method should be rewritten
            // globally without num and a place holder inserted instead.
            this->file_names(file_names_array, (char *)"xx");
            BatchValueSet(fn, "Shapeit_file_stem");
            file_stem = fn;

        } else {
            printf("Unknown option %d\n", choice);
        }
    }
}

/*
static keyw_t keywords[] = {
    {"Shapeit_recomb_rdir",                            STRING, ""},
    {"Shapeit_recomb_rfile",                           STRING, ""} ,
    {"Shapeit_file_stem",                              STRING, ""}
};
*/

void CLASS_SHAPEIT::batch_out()
{
    Cstr Values[] = { "Shapeit_recomb_rdir",
                      "Shapeit_recomb_rfile",
                      "Shapeit_file_stem",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_SHAPEIT::batch_in()
{
    Str file;
    Vecs filesplit;

    BatchValueGet(this->rdir,  "Shapeit_recomb_rdir");
    BatchValueGet(       file, "Shapeit_recomb_rfile");
    BatchValueGet(this->file_stem, "Shapeit_file_stem");

    split(filesplit, file, "?");
    if (filesplit.size() != 2) {
        printf("Please include one and only one ? in the file name\n");
        return;
    }
    rpre  = filesplit[0];
    rpost = filesplit[1];
}

void CLASS_SHAPEIT::batch_show()
{
    msgvf("\n");
    msgvf("Shapeit recombination data directory:     %s\n",    C(this->rdir));
    msgvf("Shapeit recombination file:               %s?%s\n", C(this->rpre), C(this->rpost));
    msgvf("Shapeit data file stem:                   %s\n",    C(this->file_stem));
    msgvf("\n");
 }


void CLASS_SHAPEIT::create_sh_file(linkage_ped_top *Top,
                                char *file_names_array[],
                                const int numchr)
{
/*
p Outfile_Names[0]  "2015-11-17-10-44/shapeit.05.fam"
p Outfile_Names[1]  "2015-11-17-10-44/shapeit.05.bim"
p Outfile_Names[2]  "2015-11-17-10-44/shapeit.phe"
p Outfile_Names[3]  "2015-11-17-10-44/shapeit.05.bed"
p Outfile_Names[4]  "2015-11-17-10-44/shapeit.all.sh"
p Outfile_Names[5]  "2015-11-17-10-44/shapeit_geno_summary.05"
p Outfile_Names[6]  "2015-11-17-10-44/shapeit.05.fam"
p Outfile_Names[7]  "2015-11-17-10-44/shapeit.05"
p Outfile_Names[8]  "2015-11-17-10-44/shapeit.05.sh"
p Outfile_Names[9]  "2015-11-17-10-44/shapeit.05.ref"
p Outfile_Names[10] "2015-11-17-10-44/"

0. make pedigree be unique ID (ORIGID[1] option ??)
2. indicate where output is generated.
4. which types of shapeit runs to make
5. --duohmm
*/

//    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
//        strcmp(output_paths[0], ".");
    int top_shell = 1;

    char argfile[2*FILENAME_LENGTH];
    sprintf(argfile, "%s.args", file_names_array[4]);
    dataloop::sh_exec arg(Top);
    arg.filep_open(output_paths[0], argfile, "w");
    arg.pr_puts(SHAPEIT_ARGS);
//  arg.pr_puts("set MoreArgs=--duohmm\n");
    arg.filep_close();

    
    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names_array[4], "w");
        sh->sh_main();
    }

    vlpCLASS(SHAPEIT_sh_script,both,sh_exec) {
     vlpCTOR(SHAPEIT_sh_script,both,sh_exec) { }
        typedef char *str;
        str *file_names_intrnl;
        dataloop::sh_exec *sh;
        CLASS_SHAPEIT *clss;
        
        void file_loop() {
            mssgvf("      SHAPEIT shell file:          %s/%s\n", *_opath, file_names_intrnl[8]);
            data_loop(*_opath, file_names_intrnl[8], "w");
        }
        void file_header() {
//            asm("int $3");
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
            pr_nl();

            // This handles the environment variable setup to allow the checking
            // functions in 'batch_run' to work correctly...
            fprintf_env_checkset_csh(_filep, "_SHAPEIT", "shapeit");
            pr_nl();
            pr_printf("alias usage 'echo \"Usage: %s \"\\\n", file_names_intrnl[8]);
            pr_printf("  exit'\n");
            pr_nl();
            pr_printf("if ($1 == '?' || $1 == 'help' || $#argv > 2) then\n");
            pr_printf("  usage\n");
            pr_printf("endif\n");
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];

            pr_nl();
            sprintf(cmd, "source %s.args", file_names_intrnl[4]);
            sh_echo(cmd);
            pr_nl();

            pr_printf("echo\n");
            sprintf(cmd, "$_SHAPEIT --input-bed %s %s %s --input-map %s/%s%d%s --output-max %s.haps %s.sample %s\n",
                    file_names_intrnl[3], file_names_intrnl[1], file_names_intrnl[0],
                    C(clss->rdir), C(clss->rpre), _numchr, C(clss->rpost),
                    file_names_intrnl[7], file_names_intrnl[7],
                    "$MoreArgs");

            sh_run("SHAPEIT", cmd);
            pr_nl();
            fprintf_status_check_csh(_filep, "SHAPEIT", 1);

            // can't do this because we get the status from the 'head' that we pipe the data to...
            //fprintf_status_check_csh(_filep, "SHAPEIT", 1);
        }
    } *xp = new SHAPEIT_sh_script(Top);
    
    xp->file_names_intrnl = file_names_array;
    xp->clss              = this;
    xp->sh                = sh;
    
    xp->iterate();
    delete xp;

    if (top_shell) {

        sh->filep_close();
        delete sh;
    }
}
