/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

void CLASS_SHAPEIT::save_pedsix_file(linkage_ped_top *Top,
                                     const int pwid,
                                     const int fwid)
{
    vlpCLASS(plink_pedsix,trait,ped_per) {
     vlpCTOR(plink_pedsix,trait,ped_per) { }

        void file_loop() {
            mssgvf("        PLINK pedigree file:       %s/%s\n", *_opath, Outfile_Names[0]);  //fam
            data_loop(*_opath, Outfile_Names[0], "w");
        }
        void inner() {
#if 0
            if (PLINK_OUT.no_fid != 1) pr_fam();
            pr_per();
            if (PLINK_OUT.no_parents != 1) pr_parent(); // father-id & mother-id
            if (PLINK_OUT.no_sex != 1) pr_sex();
            if (PLINK_OUT.no_pheno != 1) pr_pheno();
#endif
            pr_uid();
            pr_parent();
            pr_sex();
            pr_pheno();

            pr_nl();
        }
    } *sp = new plink_pedsix(Top);

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();

    delete sp;
}

void CLASS_SHAPEIT::user_queries(char **file_names_array,
                                 int *combine_chromo, int *create_summary)
{
    int i, choice = -1, istem = -1, idef = -1;
    int idir = -1, ifile = -1;
    char selection[100], *sp = selection;

    *combine_chromo = 0;

    print_outfile_mssg();
    while (choice != 0) {
        printf("             SHAPEIT input menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        i=1;

        batch_in();

        if (_suboption == 1) {
            printf(" %d) Specify genetic recombination map directory?       \"%s\"\n",
                   i, C(rdir));
            idir=i++;

            printf(" %d) Specify genetic recombination map file name?       \"%s\"\n",
                   i, BatchItemGet("Shapeit_recomb_rfile")->value.name);
            ifile=i++;
        }
	printf(" %d) Use default filenames?                             %s\n",
	       i, yorn[DEFAULT_OUTFILES]);
	idef=i++;

        if (! DEFAULT_OUTFILES) {
            printf(" %d) Change filenames stem?                            \"%s\"\n",
                   i, this->file_name_stem);
        istem=i++;
        }

// specify recomb map

        printf("Enter options 0-%d > ", i-1);
        fcmap(stdin, "%d", &choice); printf("\n");

        if (choice == 0) {
	    if (_suboption == 1 && rpre == "" && rpost == "") {
		printf("You must specify the recombination rate file (and if necessary the path)\n");
		choice = i;
	    }

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

        } else if (choice == idef) {
	    printf("Please type \"yes\" or \"no\" > ");
            fcmap(stdin, "%s", selection);    newline;
	    BatchValueSet(selection[0], "Default_Outfile_Names");

        } else if (choice == istem) {
            char *fn = this->file_name_stem;
            printf("Enter new stem for the output file names > ");
            fcmap(stdin, "%s", this->file_name_stem);    newline;
            BatchValueSet(fn, "Shapeit_file_stem");

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
    char *fn = this->file_name_stem;
    Str file;
    Vecs filesplit;

    BatchValueIfSet(           fn,   "Shapeit_file_stem");

    if (_suboption == 1) {
        BatchValueGet(this->rdir,  "Shapeit_recomb_rdir");
        BatchValueGet(       file, "Shapeit_recomb_rfile");

        if (file == "") return;
        split(filesplit, file, "?");
        if (filesplit.size() != 2) {
            printf("Please include one and only one ? in the file name \"%s\"\n",
                   C(file));
            return;
        }
        rpre  = filesplit[0];
        rpost = filesplit[1];
    }
}

void CLASS_SHAPEIT::batch_show()
{
    msgvf("\n");
    if (_suboption == 1) {
        msgvf("Shapeit recombination data directory:     %s\n",    C(this->rdir));
        msgvf("Shapeit recombination file:               %s?%s\n", C(this->rpre), C(this->rpost));
    }
    if (! DEFAULT_OUTFILES) {
        msgvf("Shapeit data file stem:                   %s\n",    C(this->file_name_stem));
    }
    msgvf("\n");
 }


void CLASS_SHAPEIT::create_sh_file(linkage_ped_top *Top,
                                char *file_names[],
                                const int numchr)
{
    char prefix[100];

    sub_prog_name(_suboption, prefix);
    sprintf(file_names[4], "%s_%s.all.sh", file_name_stem, prefix);
    sprintf(file_names[8], "%s_%s.%02d.sh", file_name_stem, prefix, numchr);
    add_sumdir(file_names[4]);
    add_sumdir(file_names[8]);

    switch (_suboption) {
    case 0:
    case 1:  create_sh_file_phased(Top, file_names, numchr);   break;
    case 2:  create_sh_file_check (Top, file_names, numchr);   break;
    default: break;
    }
}

void CLASS_SHAPEIT::create_sh_file_phased(linkage_ped_top *Top,
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
*/

//    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
//        strcmp(output_paths[0], ".");
    int top_shell = 1;

 
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
            pr_printf("if ($1 == '?' || $1 == 'help') then\n");
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
            sprintf(cmd, "$_SHAPEIT --input-bed %s %s %s --input-map ",
                    file_names_intrnl[3], file_names_intrnl[1], file_names_intrnl[0]);
	    if (clss->rdir != "")
		sprintf(cmd, "%s %s/", cmd, C(clss->rdir));
            sprintf(cmd, "%s%s%d%s --output-max %s.haps %s.sample %s\n",
                    cmd, C(clss->rpre), _numchr, C(clss->rpost),
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
        mssgvf("      SHAPEIT top shell file:      %s/%s\n", output_paths[0], file_names_array[4]);
        mssgvf("              the above shell runs all shells\n");

        sh->filep_close();
        delete sh;
    }

    char argfile[2*FILENAME_LENGTH];
    sprintf(argfile, "%s.args", file_names_array[4]);
    dataloop::sh_exec arg(Top);
    arg.filep_open(output_paths[0], argfile, "w");
    arg.pr_puts(SHAPEIT_ARGS);
    arg.filep_close();
    mssgvf("      SHAPEIT shell argument file: %s/%s.args\n", output_paths[0], file_names_array[4]);
}

void CLASS_SHAPEIT::create_sh_file_check(linkage_ped_top *Top,
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

2. indicate where output is generated.
*/

//    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
//        strcmp(output_paths[0], ".");
    int top_shell = 1;

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
            pr_printf("if ($1 == '?' || $1 == 'help') then\n");
            pr_printf("  usage\n");
            pr_printf("endif\n");
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];

            pr_printf("echo\n");
            sprintf(cmd, "$_SHAPEIT -check --input-bed %s %s %s\n",
                    file_names_intrnl[3], file_names_intrnl[1], file_names_intrnl[0]);

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
        mssgvf("      SHAPEIT top shell file:      %s/%s\n", output_paths[0], file_names_array[4]);
        mssgvf("              the above shell runs all shells\n");

        sh->filep_close();
        delete sh;
    }
}

void CLASS_SHAPEIT::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 0:
    case 1:  strcpy(subprog, "phased");              break;
    case 2:  strcpy(subprog, "check");               break;
    default:                                         break;
    }
}

void CLASS_SHAPEIT::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
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
            printf("Selection Menu: Shapeit shell file options\n");
            printf("0) Done with this menu - please proceed\n");
            printf("%c1) generate Shapeit phased shell files\n",
                   selection == 1 ? '*' : ' ');
            printf("%c2) generate Shapeit check shell files\n",
                   selection == 2 ? '*' : ' ');
            printf("Enter selection: 0 - 2 > ");
            fcmap(stdin,"%s", select); newline;
            sscanf(select, "%d", &selected);
            if (selected < 0 || selected > 2) warn_unknown(select);
            else if (selected) selection = selected;
        }
    }
    (*analysis)->_suboption = selection;
}

void CLASS_SHAPEIT::sub_prog_name_to_sub_option(char *subprog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)subprog_name[0])) {
    case 'p': // phased
        (*analysis)->_suboption = 1; break;
    case 'c': // check
        (*analysis)->_suboption = 2; break;
    default:
        break;
    }
}
