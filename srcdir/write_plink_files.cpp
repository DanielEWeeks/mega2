/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"

#include "loop.h"

#include "plink_core_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

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



void CLASS_PLINK::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (strcasecmp(sub_prog_name, PLINK_SUB_OPTION_LGEN) == 0 || strlen(sub_prog_name) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_LGEN_INT; // the default
    else if (strcasecmp(sub_prog_name, PLINK_SUB_OPTION_SNP_MAJOR) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_SNP_MAJOR_INT;
    else if (strcasecmp(sub_prog_name, PLINK_SUB_OPTION_INDIVIDUAL_MAJOR) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT;
    else if (strcasecmp(sub_prog_name, PLINK_SUB_OPTION_PED) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_PED_INT;
    else {
        errorvf("Valid values for PLINK 'Analysis_Sub_Option' are:\n");
        errorvf("'%s', '%s', '%s', and '%s'.\n",
                PLINK_SUB_OPTION_LGEN,
                PLINK_SUB_OPTION_PED,
                PLINK_SUB_OPTION_SNP_MAJOR,
                PLINK_SUB_OPTION_INDIVIDUAL_MAJOR);
        EXIT(BATCH_FILE_ITEM_ERROR);
    }
}

void CLASS_PLINK::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
        case 0: // missing. For backward compatability when this was the only PLINK output type produced by Mega2.
        case PLINK_SUB_OPTION_LGEN_INT:
            strcpy(subprog, PLINK_SUB_OPTION_LGEN);
            break;
        case PLINK_SUB_OPTION_SNP_MAJOR_INT:
            strcpy(subprog, PLINK_SUB_OPTION_SNP_MAJOR);
            break;
        case PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT:
            strcpy(subprog, PLINK_SUB_OPTION_INDIVIDUAL_MAJOR);
            break;
        case PLINK_SUB_OPTION_PED_INT:
            strcpy(subprog, PLINK_SUB_OPTION_PED);
            break;
        default:
            errorvf("INTERNAL: Unknown PLINK suboption value.\n");
            EXIT(SYSTEM_ERROR);
            break;
    }
}

void CLASS_PLINK::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int selection;
    char select[10];
    
    if (batchANALYSIS) {
        if (Mega2BatchItems[/* 6 */ Analysis_Sub_Option].item_read) {
            selection = (*analysis)->_suboption;
        } else {
#ifndef HIDESTATUS
            mssgvf("For backward compatablity, the output file type used for PLINK analysis\n");
            mssgvf("will be the Long Format (LGEN) since no 'Analysis_Sub_Option' was specified\n");
            mssgvf("in the MEGA2.BATCH file. Valid values for PLINK 'Analysis_Sub_Option' are:\n");
            mssgvf("'%s', '%s', '%s', and '%s'.\n",
                   PLINK_SUB_OPTION_LGEN,
                   PLINK_SUB_OPTION_PED,
                   PLINK_SUB_OPTION_SNP_MAJOR,
                   PLINK_SUB_OPTION_INDIVIDUAL_MAJOR);
#endif
            selection = PLINK_SUB_OPTION_LGEN_INT;
        }
    } else { // Interactive...
        selection = 0;
        while(!(selection == PLINK_SUB_OPTION_LGEN_INT ||
                selection == PLINK_SUB_OPTION_SNP_MAJOR_INT ||
                selection == PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT ||
                selection == PLINK_SUB_OPTION_PED_INT)) {
            draw_line();
            printf("Selection Menu: PLINK output file options\n");
            printf("%d) %s long-format\n", PLINK_SUB_OPTION_LGEN_INT, PLINK_SUB_OPTION_LGEN);
            printf("%d) %s (PLINK's default for .bed files)\n", PLINK_SUB_OPTION_SNP_MAJOR_INT, PLINK_SUB_OPTION_SNP_MAJOR);
            printf("%d) %s\n", PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT, PLINK_SUB_OPTION_INDIVIDUAL_MAJOR);
            printf("%d) %s\n", PLINK_SUB_OPTION_PED_INT, PLINK_SUB_OPTION_PED);
            printf("Enter selection: 1 - 4 > ");
            fcmap(stdin,"%s", select); newline;
            selection=0;
            sscanf(select, "%d", &selection);
            if (!(selection == PLINK_SUB_OPTION_LGEN_INT ||
                  selection == PLINK_SUB_OPTION_SNP_MAJOR_INT ||
                  selection == PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT ||
                  selection == PLINK_SUB_OPTION_PED_INT))
                warn_unknown(select);
        }
    }
    (*analysis)->_suboption = selection;
}

//
// This processes only the first six entries of the .PED file, or the complete .FAM file.
// The genotype information that is included at the end fo the .PED file
// is not included, this information is written to the .LGEN file. In binary format
// the genotype information is written to the .BIM (allele names), and .BED (values).

void CLASS_PLINK::save_pedsix_file(linkage_ped_top *Top,
				   const int pwid,
				   const int fwid)
{
    struct plink_pedsix: public loop::trait, loop::ped_per {
        plink_pedsix(linkage_ped_top *Top) : person_locus_entry(Top), loop::trait(Top), loop::ped_per(Top) { }
        void make_file() {
            mssgvf("        PLINK pedigree file:       %s/%s\n", *_opath, ::file_names[0]);  //fam
            run_loop(::file_names[0]);
        }
        void inner() {
#if 0
            if (PLINK_OUT.no_fid != 1) pr_fam();
            pr_per();
            if (PLINK_OUT.no_parents != 1) pr_parent(); // father-id & mother-id
            if (PLINK_OUT.no_sex != 1) pr_sex();
            if (PLINK_OUT.no_pheno != 1) pr_pheno();
#endif
            pr_id();
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

void CLASS_PLINK::save_ped_file(linkage_ped_top *Top,
				const int pwid,
				const int fwid,
				const int mwid)
{
    struct plink_ped: public loop::chr, loop::ped_per_loci {

        plink_ped(linkage_ped_top *Top) : person_locus_entry(Top), loop::chr(Top), loop::ped_per_loci(Top) {}
        void make_file() {
            mssgvf("        PLINK ped file:            %s/%s\n", *_opath, ::file_names[0]);
            run_loop(::file_names[0]);
        }
        void per_start() {
#if 0
            if (PLINK_OUT.no_fid != 1) pr_fam();
            pr_per();
            if (PLINK_OUT.no_parents != 1) pr_parent(); // father-id & mother-id
            if (PLINK_OUT.no_sex != 1) pr_sex();
            if (PLINK_OUT.no_pheno != 1) pr_pheno();
#endif
            pr_id();
            pr_parent();
            pr_sex();
            pr_pheno();
        }
        void inner() {
            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                // How we represent the allele in the output...
                // display the allele name according to the value associated with this person...
                pr_printf("%2s %2s ",
                          recode_name(_allele1, "0", "?"),
                          recode_name(_allele2, "0", "?"));
            } else {
                // an allele with no name....
                pr_printf("%2d %2d ", _allele1, _allele2);
            }
        }
        void per_end() {
            pr_nl();
        }
    } *sp = new plink_ped(Top);

    sp->_loci_allele_limit = 2;

    sp->load_formats(fwid, pwid, mwid);
    sp->iterate();

    delete sp;
}

//
// For a description of the PLINK binary (.bed) file, please see...
// http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
//
// mode_flag determines the .BED mode (SNP major == 1; Individual major == 2).

//
// This processes only the first six entries of the .PED file, or the complete .FAM file.
// The genotype information that is included at the end fo the .PED file
// is not included, this information is written to the .LGEN file. In binary format
// the genotype information is written to the .BIM (allele names), and .BED (values).


void CLASS_PLINK::save_bed_file(const char *bedfl_name,
                                linkage_ped_top *Top,
                                const int binary_mode_flag)
{
    if (binary_mode_flag == 1) {
        struct plink_snp_major: public loop::chr, loop::loci_ped_per, public plink_binary {

            plink_snp_major(linkage_ped_top *Top) : person_locus_entry(Top), loop::chr(Top), loci_ped_per(Top) { }
           ~plink_snp_major() {}
            void make_file() {
                mssgvf("        PLINK binary file snp:     %s/%s\n", *_opath, ::file_names[3]);
                run_loop(::file_names[3]);
            }
            void file_header() {
                plink_binary::file_header(_filep);

                // PLINK .bed file v1.00 SNP-major mode (lists all individuals for a SNP).
                fputc(0x01, _filep);
            }
            void loci_start() { SNP_count = 0; SNP_data = 0; }
              void inner() { plink_binary::inner(_filep, _allele1, _allele2); }
            void loci_end() { if ((SNP_count & 0x7) != 0) fputc(SNP_data, _filep); }
        } *xp = new plink_snp_major(Top);

        xp->iterate();
        delete xp;

    } else if (binary_mode_flag == 2) {
        struct plink_indiv_major: public loop::chr, loop::ped_per_loci, public plink_binary {
            plink_indiv_major(linkage_ped_top *Top) : person_locus_entry(Top), loop::chr(Top), ped_per_loci(Top) { }
           ~plink_indiv_major() {}
            void make_file() {
                mssgvf("        PLINK binary file indiv:   %s/%s\n", *_opath, ::file_names[3]);
                run_loop(::file_names[3]);
            }
            void file_header() {
                plink_binary::file_header(_filep);

                // PLINK .bed file v1.00 Individual-major mode (lists all SNPs for the individual).
                fputc(0x00, _filep);
            }
            void per_start() { SNP_count = 0; SNP_data = 0; }
            void inner() { plink_binary::inner(_filep, _allele1, _allele2); }
            void per_end() { if ((SNP_count & 0x7) != 0) fputc(SNP_data, _filep); }
        } *xp = new plink_indiv_major(Top);

        xp->iterate();
        delete xp;

    } else {
        errorvf("INTERNAL: binary_mode_flag unknown value?\n");
        EXIT(SYSTEM_ERROR);
    }
}

void CLASS_PLINK::create_output_file(
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
//  why not let me have no traits ...
    if (have_trait_b(LPedTreeTop->LocusTop) != 0 || 1) {
        // See #defines for PLINK_SUB_OPTION_*_INT for the value of _suboption
        create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, _suboption-1, "plink", analysis);
    } else {
        errorf("There are no traits available for column 6 of the PLINK phenotype file.");
        EXIT(EARLY_TERMINATION);
    }
}

void CLASS_PLINK::file_names_w_stem(char **file_names, char *num, const char *stem)
{
    if (_suboption == PLINK_SUB_OPTION_PED_INT)
        sprintf(file_names[0], "%s.%s.ped", stem, num);
    else // for PLINK_SUB_OPTION_SNP_MAJOR, or PLINK_SUB_OPTION_INDIVIDUAL_MAJOR.
        sprintf(file_names[0], "%s.%s.fam", stem, num);
    
    if (_suboption == PLINK_SUB_OPTION_LGEN_INT) {
        sprintf(file_names[1], "%s.%s.map", stem, num);
        sprintf(file_names[3], "%s.%s.lgen", stem, num);

    } else if (_suboption == PLINK_SUB_OPTION_PED_INT) {
        sprintf(file_names[1], "%s.%s.map", stem, num);

    } else { // for PLINK_SUB_OPTION_SNP_MAJOR, or PLINK_SUB_OPTION_INDIVIDUAL_MAJOR.
        sprintf(file_names[1], "%s.%s.bim", stem, num);
        sprintf(file_names[3], "%s.%s.bed", stem, num);
    }
    sprintf(file_names[2], "%s.phe", stem);
    //  sprintf(file_names[4], "%s.%s.sh", stem, num);
    //  sprintf(file_names[5], "%s_geno_summary.%s", stem, num);
    //  sprintf(file_names[6], "%s.%s.fam", stem, num);
    //
    sprintf(file_names[4], "%s.all.sh", stem);
    sprintf(file_names[5], "%s_geno_summary.%s", stem, num);
    sprintf(file_names[6], "%s.%s.fam", stem, num);
    //
    sprintf(file_names[7], "%s.%s", stem, num);
    sprintf(file_names[8], "%s.%s.sh", stem, num);
}

/**
 @param output_format determines how the files are written:
 .LGEN mode (== 0),
 .BED mode: (SNP major == 1; Individual major == 2).
 */

void CLASS_PLINK::create_sh_file(linkage_ped_top *Top,
                                 char *file_names[],
                                 const int numchr)
{
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1) ||
        strcmp(output_paths[0], ".");
    int output_format = _suboption - 1;

    struct all_sh: public sh_util {
        all_sh(linkage_ped_top *Top) : sh_util(Top) {}
        virtual void file_post() {
            chmod_X_file(path_);
        }
    } *sh = 0;
    
    if (top_shell) {
        sh = new all_sh(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }
    
    struct PLINK_sh_script: public loop::outer, all_sh {
        typedef char *str;
        str *file_names;
        all_sh *sh;
        int outputFormat;
        
        PLINK_sh_script(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), all_sh(Top) { }
        void make_file() {
            mssgvf("        PLINK shell file:          %s/%s\n", *_opath, file_names[8]);
            run_loop(file_names[8]);
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
            fprintf_env_checkset_csh(_filep, "_PLINK", "plink");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];
            char out_fl[2*FILENAME_LENGTH];
            
            if (_numchr > 0)
                pr_printf("echo Running PLINK on chromosome %d markers\n", _numchr);
            
            sprintf(out_fl, "%s.assoc", file_names[7]);
            sh_rm(out_fl);
            
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
            
            // From the PLINK code (input.cpp) it appears that the '--missing-phenotype' flag can
            // be associated with an affection status (what the PLINK code refers to as a "binary trait")
            // as well as a quantitative trait. For a binary trait the value "is also" (in addition to '0'
            // in the !par:coding01 case) considered a missing value. For quantitative traits it "is" the
            // missing value.
            //
            // So, it is "not wrong" to output the '--missing-phenotype' flag for a affection status written
            // by Mega2, because we only use (0, 1, 2). The missing phenotype flag value will simply not appear
            // in the input that PLINK reads.
            
#ifdef RUNSHELL_SETUP
            sprintf(cmd, "$_PLINK --noweb --%cfile %s --missing-phenotype %s --assoc --out %s\n",
                    (outputFormat == 0 ? 'l' : 'b'), file_names[7],
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read ?
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name : "-9",
                    file_names[7]);
#else /* RUNSHELL_SETUP */
            sprintf(cmd, "plink --%cfile %s --missing-phenotype %s --assoc --out %s\n",
                    (outputFormat == 0 ? 'l' : 'b'), file_names[7],
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].item_read ?
                    Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name : "-9" ,
                    file_names[7]);
#endif /* RUNSHELL_SETUP */
            sh_run("PLINK", cmd);
            fprintf_status_check_csh(_filep, "PLINK", 1);
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new PLINK_sh_script(Top);
    
    xp->file_names = file_names;
    xp->sh         = sh;
    xp->outputFormat = output_format;
    
    xp->iterate();
    
    delete xp;
    
    if (top_shell) {
        mssgvf("        PLINK top shell file:      %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("              the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}
