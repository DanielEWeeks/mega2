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
#include "loop_templates.h"


#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "plink_core_ext.h"
#include "write_plink_ext.h"

#include "write_shapeit_ext.h"

// This will be the first affection status that we find in the list of traits.
// What if there are no affection status in the list of traits?
char *fam_file_phenotype_nameX = NULL;

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

    // The missing phenotype value for quantitative traits is, by default, -9.
    // Here we use the value of 'MissingQuant' which should be derived from the batch
    // file item "Value_Missing_Quant_On_Input".
        create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, "shapeit", analysis);
}

void CLASS_SHAPEIT::create_sh_file(linkage_ped_top *Top,
                                char *file_names[],
                                const int numchr)
{
    DTshell *sh = 0;

    FLPboth *floop = new FLPboth(Top, file_names[8], "w");
    floop->file_type = "      SHAPEIT shell file:          ";

    struct SHAPEIT_sh_script: public DTshell {
        SHAPEIT_sh_script(linkage_ped_top *Top, fileloop::fileloop_data *fl) : DTshell(Top, fl) { }

        typedef char *str;
        str *file_names;
        DTshell *sh;
        
        void file_header() {
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
            
            // This handles the environment variable setup to allow the checking
            // functions in 'batch_run' to work correctly...
            fprintf_env_checkset_csh(_filep, "_SHAPEIT", "pseq");
            pr_printf("\n");
            pr_printf("alias usage 'echo \"Usage: %s [ SHAPEIT_PROJ [ SHAPEIT_RESDIR ] ]\"\\\n", file_names[8]);
            pr_printf("  echo \" SHAPEIT_PROJ    the project name\"\\\n");
            pr_printf("  echo \" SHAPEIT_RESDIR  the resource directory\"\\\n");
            pr_printf("  exit'\n");
            pr_printf("\n");
            pr_printf("if ($1 == '?' || $1 == 'help' || $#argv > 2) then\n");
            pr_printf("  usage\n");
            pr_printf("endif\n");
        }
        void inner () {

            char cmd[2*FILENAME_LENGTH];
            char out_fl[2*FILENAME_LENGTH];
            
            pr_printf("\n");

            if (_numchr > 0) {
                pr_printf("echo Running SHAPEIT on chromosome %d markers\n", _numchr);
                pr_printf("echo\n");
            }

            if (num_traits > 1) {
                if (LoopOverTrait) {
                    sprintf(out_fl, "../%s", file_names[3]); //bed
                    sh_ln(out_fl, file_names[3]);
                    
                    sprintf(out_fl, "../%s", file_names[1]); //map&bim
                    sh_ln(out_fl, file_names[1]);
                }
            }
            
            if (strcmp(file_names[0], file_names[6]))        //fam
                sh_ln(file_names[0], file_names[6]);
            
            pr_printf("\n");
            pr_printf("# Assign a project name...\n");
            pr_printf("if ($#argv > 0) then\n");
            pr_printf("  set SHAPEIT_PROJ=$argv[1]\n");
            pr_printf("  echo Using the user specified project name of \\\"$SHAPEIT_PROJ\\\".\n");
            pr_printf("else\n");
            pr_printf("  set SHAPEIT_PROJ=%s\n", file_names[7]);
            pr_printf("  echo Using the default project name of \\\"$SHAPEIT_PROJ\\\".\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("# Assign a resource directory...\n");
            pr_printf("if ($#argv > 1) then\n");
            pr_printf("  set SHAPEIT_RESDIR=$argv[2]\n");
            pr_printf("  echo Using the user specified resource directory of \\\"$SHAPEIT_RESDIR\\\".\n");
            pr_printf("else\n");
            pr_printf("  set SHAPEIT_RESDIR=${SHAPEIT_PROJ}_res\n");
            pr_printf("  echo Using the default resource directory of \\\"$SHAPEIT_RESDIR\\\".\n");
            pr_printf("endif\n");

	    // TODO: I still need to figure out how to tell what the former resource directory
	    // for a project it so that the user doesn't use something different, or warn them if
	    // they do.
            pr_printf("\n");
            pr_printf("# It is an error if the resource directory does not exist...\n");
            pr_printf("if (! -d $SHAPEIT_RESDIR) then\n");
            pr_printf("  echo\n");
            pr_printf("  echo The resource directory \\\"$SHAPEIT_RESDIR\\\" does not exist.\n");
            pr_printf("  echo\n");
            pr_printf("  usage\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("# Do not create a new project if it already exists...\n");
            pr_printf("if (! -f $SHAPEIT_PROJ) then\n");
            pr_printf("echo\n");
            pr_printf("echo ... Creating a new project ...\n");
            sprintf(cmd, "$_SHAPEIT $SHAPEIT_PROJ new-project --resources $SHAPEIT_RESDIR\n");
            sh_run("SHAPEIT", cmd);
            fprintf_status_check_csh(_filep, "SHAPEIT", 1);
            pr_printf("else\n");
            pr_printf("  echo Using existing project \\\"${SHAPEIT_PROJ}\\\".\n");
            pr_printf("endif\n");
            
            pr_printf("\n");
            pr_printf("mkdir -p ${SHAPEIT_PROJ}_out\n");
            pr_printf("if (-f ${SHAPEIT_PROJ}_out/%s.bed) then\n",file_names[7]);
            pr_printf("  echo\n");
            pr_printf("  echo ERROR: Attemping to move your trio of PLINK files into the \\\"${SHAPEIT_PROJ}_out\\\" SHAPEIT project folder.\n");
            pr_printf("  echo ERROR: The SHAPEIT project folder \\\"${SHAPEIT_PROJ}_out\\\" already contains PLINK files of the same name.\n");
            pr_printf("  exit\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("cp %s.bed %s.bim %s.fam ${SHAPEIT_PROJ}_out\n", file_names[7], file_names[7], file_names[7]);
            pr_printf("echo\n");
            pr_printf("echo Your trio of PLINK files has been moved into the \\\"${SHAPEIT_PROJ}_out\\\" SHAPEIT project folder.\n");
            pr_printf("echo Do not move or alter these files for the duration of your project.\n");

            pr_printf("\n");
            pr_printf("echo\n");
            pr_printf("echo ... Loading plink binary files ...\n");
            sprintf(cmd, "$_SHAPEIT $SHAPEIT_PROJ load-plink --file ${SHAPEIT_PROJ}_out/%s --phenotype %s --id $SHAPEIT_PROJ --check-reference\n",
                    file_names[7], fam_file_phenotype_nameX);
            sh_run("SHAPEIT", cmd);
            fprintf_status_check_csh(_filep, "SHAPEIT", 1);

            
            if (num_traits > 2) {
                pr_printf("\n");
                pr_printf("echo ... Load additional pheotypes ...\n");
                sprintf(cmd, "$_SHAPEIT $SHAPEIT_PROJ load-pheno --file %s\n", file_names[2]);
                sh_run("SHAPEIT", cmd);
                fprintf_status_check_csh(_filep, "SHAPEIT", 1);
            }
            
            pr_printf("\n");
            pr_printf("echo\n");
            pr_printf("echo ... Listing some individuals in the project/file ...\n");
            sprintf(cmd, "$_SHAPEIT $SHAPEIT_PROJ i-view | head\n");
            sh_run("SHAPEIT", cmd);

            // can't do this because we get the status from the 'head' that we pipe the data to...
            //fprintf_status_check_csh(_filep, "SHAPEIT", 1);
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new SHAPEIT_sh_script(Top, floop);

    xp->file_names = file_names;
    xp->sh         = sh;

    floop->iterate();

    delete xp;
    delete floop;
}
