/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

#include "common.h"
#include "typedefs.h"

#include "loop.h"
#include "sh_util.h"

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

#include "write_eigenstrat_ext.h"

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
  EIGENSTRAT output file created by Mega2.

  The EIGENSTRAT is a family of programs available through the EIGENSOFT software site
  which is located at this link:
  http://www.hsph.harvard.edu/alkes-price/software/

  There appear to be two verions that are currently in use EIGENSOFT 4.2 (which supports
  multi-threading), and EIGENSOFT 5.0.1 (which supports additional features and bug fixes).
  Our testing has been done with the 4.2 version.

  The EIGENSTRAT control progrema are: smartpca.perl, and smarteigenstrat.perl. 
  There is a description of the programs and their options in the file 'EIGENSTRAT/README'.
  These programs run other programs that area located in the bin directory of the release
  (e.g., EIG4.2/bin). FOR THIS REASON THE BIN DIRECTORY MUST BE IN YOUR CURRENT PATH!
 */

void CLASS_EIGENSTRAT::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (strcasecmp(sub_prog_name, EIGENSTRAT_SUB_OPTION_SNP_MAJOR) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_SNP_MAJOR_INT;
    else if (strcasecmp(sub_prog_name, EIGENSTRAT_SUB_OPTION_PED) == 0)
        (*analysis)->_suboption = PLINK_SUB_OPTION_PED_INT;
    else {
        errorvf("Valid values for Eigenstrat 'Analysis_Sub_Option' are:\n");
        errorvf("'%s', and '%s'.\n",
                EIGENSTRAT_SUB_OPTION_PED,
                EIGENSTRAT_SUB_OPTION_SNP_MAJOR);
        EXIT(BATCH_FILE_ITEM_ERROR);
    }
}

void CLASS_EIGENSTRAT::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
        case PLINK_SUB_OPTION_SNP_MAJOR_INT:
            strcpy(subprog, EIGENSTRAT_SUB_OPTION_SNP_MAJOR);
            break;
        case PLINK_SUB_OPTION_PED_INT:
            strcpy(subprog, EIGENSTRAT_SUB_OPTION_PED);
            break;
        default:
            errorvf("INTERNAL: Unknown EIGENSTRAT suboption value.\n");
            EXIT(SYSTEM_ERROR);
            break;
    }
}

void CLASS_EIGENSTRAT::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int selection = 0;
    char select[10];
    
    if (batchANALYSIS) {
        if (Mega2BatchItems[/* 6 */ Analysis_Sub_Option].items_read) {
            selection = (*analysis)->_suboption;
        } else {
            errorvf("Valid values for Eigenstrat 'Analysis_Sub_Option' are:\n");
            errorvf("'%s', and '%s'.\n",
                    EIGENSTRAT_SUB_OPTION_PED,
                    EIGENSTRAT_SUB_OPTION_SNP_MAJOR);
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
    } else { // Interactive...
        selection = 0;
        while(!(selection == PLINK_SUB_OPTION_SNP_MAJOR_INT ||
                selection == PLINK_SUB_OPTION_PED_INT)) {
            draw_line();
            printf("Selection Menu: output file options\n");
            printf("%d) %s\n", PLINK_SUB_OPTION_SNP_MAJOR_INT, EIGENSTRAT_SUB_OPTION_SNP_MAJOR);
            printf("%d) %s\n", PLINK_SUB_OPTION_PED_INT, EIGENSTRAT_SUB_OPTION_PED);
            printf("Enter selection: %d or %d > ", PLINK_SUB_OPTION_SNP_MAJOR_INT, PLINK_SUB_OPTION_PED_INT);
            fcmap(stdin,"%s", select); newline;
            selection=0;
            sscanf(select, "%d", &selection);
            if (!(selection == PLINK_SUB_OPTION_SNP_MAJOR_INT ||
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

void CLASS_EIGENSTRAT::save_pedsix_file(linkage_ped_top *Top,
                                        const int pwid,
                                        const int fwid)
{
    vlpCLASS(eigenstrat_pedsix,trait,ped_per) {
     vlpCTOR(eigenstrat_pedsix,trait,ped_per) { }
        int missing_affection_status;

        void file_loop() {
            //if (_ttraitp == (linkage_locus_rec *)NULL) return;
            mssgvf("        EIGENSTRAT pedigree file:  %s/%s\n", *_opath, Outfile_Names[0]);  //fam
            missing_affection_status = 0;
            data_loop(*_opath, Outfile_Names[0], "w");
        }
        void inner() {
            if (!has_pheno()) {
                missing_affection_status++;
                // Do not write out data associated with individuals that have unknown affection status
                return;
            }
            pr_id();
            pr_parent();
            pr_sex();
            pr_pheno(1);
            pr_nl();
        }
        void file_post() {
            if (missing_affection_status > 0) {
                warnvf("%d entries with missing affection status removed from the pedigree file.\n",
		       missing_affection_status);
                warnf("This is because smarteigenstrat does not accept individuals with unknown affection status.");
            }
        }
    } *sp = new eigenstrat_pedsix(Top);
    
    sp->load_formats(fwid, pwid, -1);
    
    sp->iterate();
    
    delete sp;
}

void CLASS_EIGENSTRAT::save_ped_file(linkage_ped_top *Top,
                                     const int pwid,
                                     const int fwid,
                                     const int mwid)
{
    vlpCLASS(eigenstrat_ped,both,ped_per_loci) {
     vlpCTOR(eigenstrat_ped,both,ped_per_loci) { }
        int missing_affection_status;
        int process_per;

        void file_loop() {
            //if (_ttraitp == (linkage_locus_rec *)NULL) return;
            mssgvf("        EIGENSTRAT pedigree file:  %s/%s\n", *_opath, Outfile_Names[0]);
            missing_affection_status = 0;
            data_loop(*_opath, Outfile_Names[0], "w");
        }
        void per_start() {
            process_per = 1;
            if (!has_pheno()) {
                missing_affection_status++;
                process_per = 0;
                // Do not write out data associated with individuals that have unknown affection status
                return;
            }
            pr_id();
            pr_parent();
            pr_sex();
            pr_pheno();
        }
        void inner() {
            if (process_per == 0) return;
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
            if (process_per == 0) return;
            pr_nl();
        }
        void file_post() {
            if (missing_affection_status > 0) {
                warnvf("%d entries with missing affection status removed from the pedigree file.\n",
		       missing_affection_status);
                warnf("This is because smarteigenstrat does not accept individuals with unknown affection status.");
            }
        }
    } *sp = new eigenstrat_ped(Top);
    
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

// Eigenstrat only uses a PLINK SNP Major binary file.
// We also must not write out data associated with individuals that have unknown affection status

void CLASS_EIGENSTRAT::save_bed_file(const char *bedfl_name,
                                     linkage_ped_top *Top,
                                     const int binary_mode_flag)
{
    if (binary_mode_flag == 1) {
        struct eigenstrat_snp_major: public fileloop::both, public dataloop::loci_ped_per, public plink_binary {
         vlpCTOR(eigenstrat_snp_major,both,loci_ped_per) { }

            ~eigenstrat_snp_major() {}
            void file_loop() {
                //if (_ttraitp == (linkage_locus_rec *)NULL) return;
                mssgvf("        EIGENSTRAT binary file snp: %s/%s\n", *_opath, Outfile_Names[3]);
                data_loop(*_opath, Outfile_Names[3], write_binary);
            }
            void file_header() {
                plink_binary::file_header(_filep);
                
                // PLINK .bed file v1.00 SNP-major mode (lists all individuals for a SNP).
                fputc(0x01, _filep);
                SNP_bufsz = ( _Top->IndivCnt + 3 ) >> 2;
                SNP_buf = CALLOC((size_t) SNP_bufsz, unsigned char);
            }
            void file_trailer() { free(SNP_buf); }
            void loci_start() { SNP_count = 0; SNP_data = 0; SNP_cp = SNP_buf; }
            void inner() {
                // NOTE: For Eigenstrat non-quantitative (affection status) phenotype,
                // values must be either 'Case' or 'Control'. Eigenstrat seems to have no
                // concept of 'missing' for these.
                if (!has_pheno()) return;
                plink_binary::inner(_filep, _allele1, _allele2);
            }
//          void loci_end() { if ((SNP_count & 0x7) != 0) fputc(SNP_data, _filep); }
            void loci_end() { if ((SNP_count & 0x7) != 0) *SNP_cp++ = (unsigned char)SNP_data;
                              fwrite(SNP_buf, 1, SNP_bufsz, _filep); }
        } *xp = new eigenstrat_snp_major(Top);
        
        xp->iterate();
        delete xp;
        
    } else if (binary_mode_flag == 2) {
        errorvf("Eigenstrat does not use individual major bed files?\n");
        EXIT(SYSTEM_ERROR);
        
    } else {
        errorvf("INTERNAL: binary_mode_flag unknown value?\n");
        EXIT(SYSTEM_ERROR);
    }
}

void CLASS_EIGENSTRAT::create_output_file(linkage_ped_top *LPedTreeTop,
                                          analysis_type *analysis,
                                          char *file_names[],
                                          int untyped_ped_opt,
                                          int *numchr,
                                          linkage_ped_top **Top2) {
    // Eigenstrat was designed for the correction of trait values based on ancestry
    // (population stratification). This implies that you need a trait.
    // http://www.nature.com/ng/journal/v38/n8/full/ng1847.html
    if (have_trait_b(LPedTreeTop->LocusTop) != 0) {
        // See #defines for PLINK_SUB_OPTION_*_INT for the value of _suboption
        create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, _suboption-1, "eigenstrat", analysis);
    } else {
        errorvf("Eigenstrat was designed for the correction of trait values based on ancestry\n");
        errorvf("(population stratification).\n");
        errorvf("This implies that you need a trait, but none was specified.\n");
        EXIT(EARLY_TERMINATION);
    }
}


void CLASS_EIGENSTRAT::file_names_w_stem(char **file_names, char *num, const char *stem,
                                         const int suboption)
{
    if (suboption == PLINK_SUB_OPTION_PED_INT) {
        sprintf(file_names[1], "%s.%s.map", stem, num); // use .map file format not .snp
        sprintf(file_names[0], "%s.%s.ped", stem, num); // use .ped file format not .geno
    } else { // for PLINK_SUB_OPTION_SNP_MAJOR
        sprintf(file_names[0], "%s.%s.pedind", stem, num);
        //sprintf(file_names[0], "%s.%s.fam", stem, num);
        sprintf(file_names[1], "%s.%s.pedsnp", stem, num);
        //sprintf(file_names[1], "%s.%s.bim", stem, num);
        sprintf(file_names[3], "%s.%s.bed", stem, num);
    }
    sprintf(file_names[2], "%s.phe", stem);
    
    sprintf(file_names[4], "%s.top.sh", stem);
    sprintf(file_names[5], "%s_geno_summary.%s", stem, num);
    sprintf(file_names[6], "%s.%s.pedind", stem, num);
    //sprintf(file_names[6], "%s.%s.fam", stem, num);
    
    sprintf(file_names[7], "%s.%s", stem, num);
    sprintf(file_names[8], "%s.%s.sh", stem, num);
    
    sprintf(file_names[9], "%s.%s.pca", stem, num);
    sprintf(file_names[10], "%s.%s.plot", stem, num);
    sprintf(file_names[11], "%s.%s.eval", stem, num);
    sprintf(file_names[12], "%s.%s.p_log", stem, num);
    sprintf(file_names[13], "%s.%s.chisq", stem, num);
    sprintf(file_names[14], "%s.%s.e_log", stem, num);
    
    // NOTE: 15, and 16 are used as error files. See output_file_names.cpp:default_outfile_names()
}


void CLASS_EIGENSTRAT::replace_chr_number(char *file_names[], int numchr) {
    
    CLASS_PLINK_CORE::replace_chr_number(file_names, numchr);
    
    // Files specific to Eigenstrat....
    change_output_chr(file_names[9], numchr);
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[12], numchr);
    change_output_chr(file_names[13], numchr);
    change_output_chr(file_names[14], numchr);
}


void CLASS_EIGENSTRAT::create_sh_file(linkage_ped_top *Top,
                                      char *file_names[],
                                      const int numchr) {
    
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1);
    
    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }
    
    vlpCLASS(EIGENSTRAT_sh_script,both,sh_exec) {
     vlpCTOR(EIGENSTRAT_sh_script,both,sh_exec) { }
        typedef char *str;
        str *file_names;
        dataloop::sh_exec *sh;
        int subOption;
        
        void file_loop() {
            //if (_ttraitp == (linkage_locus_rec *)NULL) return;
            mssgvf("        EIGENSTRAT shell file:     %s/%s\n", *_opath, file_names[8]);
            data_loop(*_opath, file_names[8], "w");
        }
        void file_header() {
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
            pr_printf("\n");
            pr_printf("echo 'IMPORTANT NOTICE'\n");
            pr_printf("echo 'You MUST put the EIGENSTRAT bin directory in your PATH BEFORE RUNNING THIS SCRIPT.'\n");
            pr_printf("echo\n");
            pr_printf("echo 'This will allow the files run by the perl scripts (smartpcl.perl, and'\n");
            pr_printf("echo 'smarteigenstrat.perl) to be found (as suggested by the documentation found in'\n");
            pr_printf("echo 'EIGx.y/EIGENSTRAT/README).'\n");
            pr_printf("echo \"To do this let us assume that the EIGENSTRAT bin directory is located at '/usr/local/src/EIG4.2/bin'.\"\n");
            pr_printf("echo 'If using the csh use the command line:'\n");
            pr_printf("echo 'set path = ($path /usr/local/src/EIG4.2/bin)'\n");
            pr_printf("echo 'If using bash/ksh/sh use the command line:'\n");
            pr_printf("echo 'export PATH=$PATH:/usr/local/src/EIG4.2/bin'\n");
            pr_printf("echo\n");
            pr_printf("\n");
#ifdef TEST
            fprintf_env_checkset_csh(_filep, "_PLINK", "plink");
#endif /* TEST */
        }
        void inner () {
            char out_fl[2*FILENAME_LENGTH];
#ifdef TEST
            char cmd[2*FILENAME_LENGTH];
#endif /* TEST */
            
            if (_numchr > 0)
                pr_printf("echo Running EIGENSTRAT on chromosome %d markers\n", _numchr);
            
            if (num_traits > 1) {
                if (LoopOverTrait) {
                    sprintf(out_fl, "../%s", file_names[3]); //bed
                    sh_ln(out_fl, file_names[3]);
                    
                    sprintf(out_fl, "../%s", file_names[1]); //map&bim
                    sh_ln(out_fl, file_names[1]);
                }
            }

            // For Eigenstrat only two of the PLINK suboptions are valid...
            if (subOption == PLINK_SUB_OPTION_PED_INT) {
                pr_printf("smartpca.perl -i %s -a %s -b %s -k 2 -o %s -p %s -e %s -l %s -m 5 -t 2 -s 6.0\n",
                          file_names[0], // .ped (.geno)
                          file_names[1], // .map (.snp)
                          file_names[0], // .ped (.ind)
                          
                          file_names[9], // .pca
                          file_names[10], // .plot
                          file_names[11], // .eval
                          file_names[12] // .p_log
                    );
                pr_printf("\n");
                // NOTE: Phenotype information is contained in '-b file.ind' and is dependent
                // on the '-q' flag.
                // Affection data if '-q NO' (default) encoded as "Case"/"Control" labels
                // Quantitative phenotypes if '-q YES' implies real numbers with -100.0
                // signifying missing data.
                pr_printf("smarteigenstrat.perl -i %s -a %s -b %s -p %s -k 1 -o %s -l %s%s\n",
                          file_names[0], // .ped
                          file_names[1], // .map
                          file_names[0], // .ped
                          
                          file_names[9], // .pca
                          file_names[13], // .chisq
                          file_names[14], // .e_log
                          (_ttraitp->Type == QUANT ? " -q YES" : "")
                    );
                pr_printf("\n");
#ifdef TEST
                // It seems that Eigenstrat is happy to take broken PLINK files and not complain under
                // some circumstances. Here, we run the output produced by Mega2 through plink and check
                // the log file for errors.
                pr_printf("echo\n");
                pr_printf("echo 'As a test run the Mega2 output through PLINK to see if it complains.'\n");
                sprintf(cmd, "$_PLINK --noweb --file %s --out %s.plink_check --test-all\n",
                        file_names[7], // file root, expects a .ped & .map file
                        file_names[7]  // for log file
                        );
                sh_run("PLINK", cmd);
                fprintf_status_check_csh(_filep, "PLINK", 1);
#endif /* TEST */
            } else { // for PLINK_SUB_OPTION_SNP_MAJOR
                pr_printf("smartpca.perl -i %s -a %s -b %s -k 2 -o %s -p %s -e %s -l %s -m 5 -t 2 -s 6.0\n",
                          file_names[3], // .bed
                          file_names[1], // .bim (.pedsnp)
                          file_names[0], // .fam (.pedind)
                          
                          file_names[9], // .pca
                          file_names[10], // .plot
                          file_names[11], // .eval
                          file_names[12] // .p_log
                    );
                pr_printf("\n");
                pr_printf("smarteigenstrat.perl -i %s -a %s -b %s -p %s -k 1 -o %s -l %s%s\n",
                          file_names[3], // .bed
                          file_names[1], // .bim
                          file_names[0], // .fam
                          
                          file_names[9], // .pca
                          file_names[13], // .chisq
                          file_names[14], // .e_log
                          (_ttraitp->Type == QUANT ? " -q YES" : "")
                    );
                pr_printf("\n");
#ifdef TEST
                pr_printf("echo\n");
                pr_printf("echo 'As a test run the Mega2 output through PLINK to see if it complains.'\n");
                pr_printf("awk '{if ($6 == \"Case\") {$6 = \"2\"} else if ($6 == \"Control\") {$6 = \"1\"}; print $1, $2, $3, $4, $5, $6;}' %s > %s.fam\n",
                          file_names[0], // .fam
                          file_names[0]
			  );
                // Use either --bfile {root} OR --bed {name} --bim {name} --fam {name}
                sprintf(cmd, "$_PLINK --noweb --bed %s --bim %s --fam %s.fam --out %s.plink_check --test-all\n",
                        file_names[3], // .bed
                        file_names[1], // .bim (.pedsnp)
                        file_names[0], // .fam (.pedind)
                        file_names[7]  // for log file
                        );
                sh_run("PLINK", cmd);
                fprintf_status_check_csh(_filep, "PLINK", 1);
#endif /* TEST */
            }
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new EIGENSTRAT_sh_script(Top);
    
    xp->file_names = file_names;
    xp->sh         = sh;
    xp->subOption = _suboption;
    
    xp->iterate();
    
    delete xp;
    
    if (top_shell) {
        mssgvf("        EIGENSTRAT top shell file: %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("                   the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}
