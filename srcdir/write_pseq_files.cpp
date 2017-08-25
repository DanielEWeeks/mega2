/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

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
#include "write_pseq_ext.h"

// This will be the first affection status that we find in the list of traits.
// What if there are no affection status in the list of traits?
char *fam_file_phenotype_name = NULL;

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

For PSEQ documentation see: http://atgu.mgh.harvard.edu/plinkseq/input.shtml#plink

By default, the phenotype encoded in column 6 of the FAM file is assumed to be a dichotomous
case/control phenotype, and is labelled 'phe1'. To give a different AlleleName, add: --phenotype t2d

Discussion...

It is also clear that there is no explicit way telling PSEQ if there is a quantitative phenotype
(from above PSEQ documenation it is assumed to be affection status (e.g., case/control)) in the
.FAM file (column 6). However the code in 'plinkseq/lib/bed.cpp:BEDReader::read_fam()' figures
out if the trait is quantitative or case/control by looking at the values. If all values are in
the set {-9,0,1,2} then the trait is defined as case/control any value outside of that set causes
the trait to be defined as quantitative.

In reading over the PSEQ code, it does not seem that it can take a five column .FAM file (e.g.,
PLINK command line option --no-pheno), if there is such a thing. The PLINK documentation states
the the .FAM file should be the first six columns of the .PED file. None the less, in PSEQ mode
you must give one phenotype.


Possible source of error in running PDEQ...

If the PSEQ database is not deleted between runs of the script there is an extra column geneated,
one with the label given with the '--phenotype' command line argument and one with a label of 'phe1'.

*/

static void save_PSEQ_pheno(const char *phenofl_name, linkage_ped_top *Top,
                            const int pwid, const int fwid)
{
    vlpCLASS(pseq_pheno,once,ped_per_trait) {
     vlpCTOR(pseq_pheno,once,ped_per_trait) { }
        
        bool use_fid, use_iid, use_joint;
        // The '_trait' value of the first trait which will be found in the .FAM file.
        // This trait will not go into this file (the PSEQ pheno file).
        int skip_trait;
        
        void file_loop() {
            // Here we loop over the pedigrees and individuals to determine what to use
            // as an "ID" for the PSEQ phenotype file. This should be the same as what PSEQ
            // uses for the "ID". This code is found in...
            // lib/bed.cpp:int BEDReader::read_fam()
            //
            // Similar methods are used below...
            
            std::ostringstream convert;   // stream used for the conversion
            std::set<std::string> sfid;
            std::set<std::string> siid;
            size_t ni = 0;
            
            // for each pedigree (consulting the linkage_ped_top structure)...
            for (_ped=0; _ped < _Top->PedCnt; _ped++) {
                // It there is some indication as to why this pedigree should not be included, don't...
                if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
                _tpedtreep = &(_Top->Ped[_ped]);
//yy
                if (OrigIds[1] == 2 || OrigIds[1] == 4) sfid.insert(_tpedtreep->Name);
                else if (OrigIds[1] == 3) {convert << _ped+1; sfid.insert(convert.str()); }
                else if (OrigIds[1] == 6) sfid.insert(_tpedtreep->PedPre);
                else {convert << _tpedtreep->Num; sfid.insert(convert.str()); }
                
                for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
                    _tpersonp = &(_tpedtreep->Entry[_per]);
                    if (OrigIds[0] == 1 || OrigIds[0] == 2)  siid.insert(_tpersonp->OrigID);
                    else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) siid.insert(_tpersonp->UniqueID);
                    else if (OrigIds[0] == 6) siid.insert(_tpersonp->PerPre);
                    else {convert << _tpersonp->ID; siid.insert(convert.str()); }
                    ni++;
                }
            }
            
            use_fid = ni == sfid.size();
            use_iid = ni == siid.size();
            use_joint = !use_fid && !use_iid;
            
            msgvf("         PSEQ phenotype file:      %s/%s\n", *_opath, Outfile_Names[2]);
            data_loop(*_opath, Outfile_Names[2], "w");
        }
        void file_header() {
            int tr, first;
            // The PSEQ header takes one line per phenotype listed in the file, followed
            // by a header line that has the "ID" followed by the phenotype names.
            // NOTE: PLINK/SEQ is really picky about having or not havings spaces/tabs in the right places.
            
            // The first lines (one line for each trait) in the file contain trait meta data...
            for (tr = 0, first = 1; tr < num_traits; tr++) {
                if (global_trait_entries[tr] == -1) continue; // If this is a marker, skip it...
                if (first == 1) {
                    skip_trait = global_trait_entries[tr];
                    //fam_file_phenotype_name = _Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName;
                    first = 0;
                    continue;
                }
                switch (_Top->LocusTop->Locus[global_trait_entries[tr]].Type) {
                    case QUANT:
                        pr_printf("##%s,Float,%s,\"Quantitative trait description\"\n",
                                  _Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName,
                                  Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].items_read ?
                                  Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name : "-9"
                                  );
                        break;
                    case AFFECTION:
                        // Since we have no way of specifying the value of the missing affection
                        // status, we will use the customary default of '0'...
                        // The assumption here is that: Case == 2, Control == 1, and Missing == 0
                        pr_printf("##%s,Integer,0,\"Affection status description\"\n",
                                  _Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                        break;
                    default:
                        // This is an internal error....
                        break;
                }
            }
            
            // The last line of the header containing the ID followed by tab delimited trait names...
            pr_printf("#ID");
            for (tr = 0, first = 1; tr < num_traits; tr++) {
                if (global_trait_entries[tr] == -1) continue; // If this is a marker, skip it...
                if (first == 1) { first = 0; continue; }
                pr_printf("\t%s", _Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName);
            }
            pr_nl();
        }
        void per_start() {
            // The "ID" comes in one of three forms depending on the input data:
            // A) if each family is unique then the family ID is used (use_id == 1),
            // B) if the family ids are not unique, but the individual ids are, then the individual id is used (use_id == 2),
            // C) finally if neither famiy, nor individual ids are unique a concatination of the two is used (use_id == 3).
            // with an interveinign underscore.
            // Consult the PSEQ code found in this method:
            // lib/bed.cpp:int BEDReader::read_fam( )
            if (use_joint || use_fid) pr_fam();
            if (use_joint) pr_printf("_");
            if (use_joint || use_iid) pr_per();
        }
        void inner()     { if (skip_trait != _trait) { pr_printf("\t"); pr_pheno(); } }
        void per_end()   { pr_nl(); }
    } *sp = new pseq_pheno(Top);
    
    // So that we can make up a string like "FID_IID".
    // In addition PSEQ requires tab delimited data, it seems to get confused if sperious spaces are introduced.
    sp->load_formats_no_space(-1);
    
    sp->iterate();
    
    delete sp;
}

void CLASS_PSEQ::save_pheno_file(linkage_ped_top *Top,
                                 const int pwid,
                                 const int fwid)
{
    int tr;
    // Even if there is no PSEQ pheno file generated, we still need to know the name of the
    // phenotype that will go into the .FAM file for the PSEQ --phenotype command line argument.
    for (tr = 0; tr < num_traits; tr++) {
        if (global_trait_entries[tr] == -1) continue; // If this is a marker, skip it...
        fam_file_phenotype_name = Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName;
        break;
    }
    // The first trait will be placed in the .FAM file, all others in the PSEQ pheno file.
    // Remember to account for the '-1' in this list, so there are N-1 traits actually
    // listed in global_trait_entries[]
    if (num_traits > 2) save_PSEQ_pheno(Outfile_Names[2], Top, pwid, fwid);
}

void CLASS_PSEQ::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    // The missing phenotype value for quantitative traits is, by default, -9.
    // Here we use the value of 'MissingQuant' which should be derived from the batch
    // file item "Value_Missing_Quant_On_Input".
    //if (have_trait_b(LPedTreeTop->LocusTop) != 0 && PLINK_OUT.no_pheno == 1) {
    if (have_trait_b(LPedTreeTop->LocusTop) != 0) {
        // See #defines for PLINK_SUB_OPTION_*_INT for the value of _suboption
        // NOTE: Analysis_Sub_Option MUST BE 'SNP major binary', the use of 'Individual major binary'                       
        // will generate the pseq message, and data will not be properly loaded...                                          
        // plinkseq warning: problems detected loading BED file, ./plink.all     
        create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, "pseq", analysis);
    } else {
        errorf("There are no traits available for column 6 of the PLINK phenotype file.");
        EXIT(EARLY_TERMINATION);
    }
}

void CLASS_PSEQ::create_sh_file(linkage_ped_top *Top,
                                char *file_names[],
                                const int numchr)
{
    dataloop::sh_exec *sh = 0;

    vlpCLASS(PSEQ_sh_script,both,sh_exec) {
     vlpCTOR(PSEQ_sh_script,both,sh_exec) { }
        typedef char *str;
        str *file_names;
        dataloop::sh_exec *sh;
        
        void file_loop() {
            mssgvf("         PSEQ shell file:          %s/%s\n", *_opath, file_names[8]);
            data_loop(*_opath, file_names[8], "w");
        }
        void file_header() {
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
            
            // This handles the environment variable setup to allow the checking
            // functions in 'batch_run' to work correctly...
            fprintf_env_checkset_csh(_filep, "_PSEQ", "pseq");
            pr_printf("\n");
            pr_printf("alias usage 'echo \"Usage: %s [ PSEQ_PROJ [ PSEQ_RESDIR ] ]\"\\\n", file_names[8]);
            pr_printf("  echo \" PSEQ_PROJ    the project name\"\\\n");
            pr_printf("  echo \" PSEQ_RESDIR  the resource directory\"\\\n");
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
                pr_printf("echo Running PSEQ on chromosome %d markers\n", _numchr);
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
            pr_printf("  set PSEQ_PROJ=$argv[1]\n");
            pr_printf("  echo Using the user specified project name of \\\"$PSEQ_PROJ\\\".\n");
            pr_printf("else\n");
            pr_printf("  set PSEQ_PROJ=%s\n", file_names[7]);
            pr_printf("  echo Using the default project name of \\\"$PSEQ_PROJ\\\".\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("# Assign a resource directory...\n");
            pr_printf("if ($#argv > 1) then\n");
            pr_printf("  set PSEQ_RESDIR=$argv[2]\n");
            pr_printf("  echo Using the user specified resource directory of \\\"$PSEQ_RESDIR\\\".\n");
            pr_printf("else\n");
            pr_printf("  set PSEQ_RESDIR=${PSEQ_PROJ}_res\n");
            pr_printf("  echo Using the default resource directory of \\\"$PSEQ_RESDIR\\\".\n");
            pr_printf("endif\n");

	    // TODO: I still need to figure out how to tell what the former resource directory
	    // for a project it so that the user doesn't use something different, or warn them if
	    // they do.
            pr_printf("\n");
            pr_printf("# It is an error if the resource directory does not exist...\n");
            pr_printf("if (! -d $PSEQ_RESDIR) then\n");
            pr_printf("  echo\n");
            pr_printf("  echo The resource directory \\\"$PSEQ_RESDIR\\\" does not exist.\n");
            pr_printf("  echo\n");
            pr_printf("  usage\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("# Do not create a new project if it already exists...\n");
            pr_printf("if (! -f $PSEQ_PROJ) then\n");
            pr_printf("echo\n");
            pr_printf("echo ... Creating a new project ...\n");
            sprintf(cmd, "$_PSEQ $PSEQ_PROJ new-project --resources $PSEQ_RESDIR\n");
            sh_run("PSEQ", cmd);
            fprintf_status_check_csh(_filep, "PSEQ", 1);
            pr_printf("else\n");
            pr_printf("  echo Using existing project \\\"${PSEQ_PROJ}\\\".\n");
            pr_printf("endif\n");
            
            pr_printf("\n");
            pr_printf("mkdir -p ${PSEQ_PROJ}_out\n");
            pr_printf("if (-f ${PSEQ_PROJ}_out/%s.bed) then\n",file_names[7]);
            pr_printf("  echo\n");
            pr_printf("  echo ERROR: Attemping to move your trio of PLINK files into the \\\"${PSEQ_PROJ}_out\\\" PSEQ project folder.\n");
            pr_printf("  echo ERROR: The PSEQ project folder \\\"${PSEQ_PROJ}_out\\\" already contains PLINK files of the same name.\n");
            pr_printf("  exit\n");
            pr_printf("endif\n");

            pr_printf("\n");
            pr_printf("cp %s.bed %s.bim %s.fam ${PSEQ_PROJ}_out\n", file_names[7], file_names[7], file_names[7]);
            pr_printf("echo\n");
            pr_printf("echo Your trio of PLINK files has been moved into the \\\"${PSEQ_PROJ}_out\\\" PSEQ project folder.\n");
            pr_printf("echo Do not move or alter these files for the duration of your project.\n");

            pr_printf("\n");
            pr_printf("echo\n");
            pr_printf("echo ... Loading plink binary files ...\n");
            sprintf(cmd, "$_PSEQ $PSEQ_PROJ load-plink --file ${PSEQ_PROJ}_out/%s --phenotype %s --id $PSEQ_PROJ --check-reference\n",
                    file_names[7], fam_file_phenotype_name);
            sh_run("PSEQ", cmd);
            fprintf_status_check_csh(_filep, "PSEQ", 1);

            
            if (num_traits > 2) {
                pr_printf("\n");
                pr_printf("echo ... Load additional pheotypes ...\n");
                sprintf(cmd, "$_PSEQ $PSEQ_PROJ load-pheno --file %s\n", file_names[2]);
                sh_run("PSEQ", cmd);
                fprintf_status_check_csh(_filep, "PSEQ", 1);
            }
            
            pr_printf("\n");
            pr_printf("echo\n");
            pr_printf("echo ... Listing some individuals in the project/file ...\n");
            sprintf(cmd, "$_PSEQ $PSEQ_PROJ i-view | head\n");
            sh_run("PSEQ", cmd);

            // can't do this because we get the status from the 'head' that we pipe the data to...
            //fprintf_status_check_csh(_filep, "PSEQ", 1);
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new PSEQ_sh_script(Top);
    
    xp->file_names = file_names;
    xp->sh         = sh;
    
    xp->iterate();
    
    delete xp;
}
