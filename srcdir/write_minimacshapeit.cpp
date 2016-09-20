/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
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
#include <ctime>

#include "common.h"
#include "typedefs.h"

#include "loop.h"
#include "sh_util.h"
#include "batch_input.h"
#include "utils_ext.h"

#include "fcmap_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"

#include "write_minimacshapeit_ext.h"


//static void write_MINIMAC_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MINIMAC_sh(linkage_ped_top *Top, char *file_names[]);

static void inner_file_names(char **file_names, const char *num, const char *stem = "minimac");

//We need some general reference varablies, these are necessary to get things into the shell output.
extern int g_cpus;
Str s_haplotype_file;
Str legend_file;
Str sample_file;
Str map_file;
Str m_haplotype_file;
int haps_sample_selected;
Str genetic_map_directory;
Str shapeit_panel_directory;
Str minimac_panel_directory;


//Main function, sets properties, calls PLINK file creation, calls shell creation
void CLASS_MINIMAC::create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    char prefix[100];
    linkage_ped_top *Top = LPedTreeTop;


    if ( InputMode == INTERACTIVE_INPUTMODE ) {
        minimac_option_menu(file_names,prefix);
    }
    else {
        batch_in();
        inner_file_names(file_names, "", file_name_stem);
    }

    int combine_chromo = 0;
    LoopOverChrm  = 1;

    LoopOverTrait = 0;
    num_traits = 0;

    //We need to check for only ACTG as alleles.
    const char *allele_name;
    for (int locus = Top->LocusTop->PhenoCnt; locus < Top->LocusTop->LocusCnt; locus++){
        for (int allele = 0; allele < Top->LocusTop->Locus[locus].AlleleCnt; allele++){
            allele_name = Top->LocusTop->Locus[locus].Allele[allele].AlleleName;
            //printf("%s\n",allele_name);
            if (allele_name == 0 || !((strcmp(allele_name, "A") == 0) || (strcmp(allele_name, "C") == 0) || (strcmp(allele_name, "G") == 0) || (strcmp(allele_name, "T") == 0) || (strcmp(allele_name, "0") == 0) || (strcmp(allele_name, "dummy") == 0))) {
                char error[255];
                strcpy(error, "The SHAPEIT Minimac3 pipeline requires alleles to be labeled as A,C,T,G.\nInvalid allele label: ");
                if (allele_name)
                    strcat(error, allele_name);
                else {
                    int l = strlen(error);
                    sprintf(&error[l], "%d <numeric>", Top->LocusTop->Locus[locus].Allele[allele].index);
                }
                errorf(error);
                EXIT(DATA_TYPE_ERROR);
            }
        }
    }


    omit_peds(untyped_ped_opt, Top);

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    (*analysis)->_suboption = PLINK_SUB_OPTION_SNP_MAJOR_INT;

    create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, file_name_stem, analysis);

    //Hindsight 20/20 this is not necessary as it was used for MaCH2VCF, I'll keep the function in, in case I realize we do need it but for now I don't think we need a snps file
    //write_MINIMAC_snps(Top, file_names, pwid, fwid);

    write_MINIMAC_sh(Top, file_names);
}

// Format of file:
// CHR:PhysicalMapDistance
// This file is used by Minimac3
//static void write_MINIMAC_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
//{
//    vlpCLASS(minimac_snp,chr,loci) {
//        vlpCTOR(minimac_snp,chr,loci) { }
//
//        typedef char *str;
//        str *file_names;
//
//        void file_loop() {
//            msgvf("        Minimac3 SNP File:         %s/%s\n", *_opath, file_names[10]);
//            data_loop(*_opath, file_names[10], "w");
//        }
//
//        void inner() {
//            pr_printf("%d:",_numchr);
//            pr_physical_distance(NULL);
//            pr_nl();
//        }
//
//    } *minimac_snps = new minimac_snp(Top);
//
//    minimac_snps->file_names = file_names;
//
//    minimac_snps->load_formats(fwid, pwid, -1);
//
//    minimac_snps->iterate();
//
//    delete minimac_snps;
//}

//Outputs the Shell file
//Looks for Shapeit/Minimac3 if it can find them it runs them
//First it runs shapeit check, then phase then convert
//Finally it runs Minimac3
static void write_MINIMAC_sh(linkage_ped_top *Top, char *file_names[]) {

    int top_shell = 1;

    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }

    vlpCLASS(minimac_sh, both, sh_exec) {
        vlpCTOR(minimac_sh, both, sh_exec) { }

        void file_loop() {
            mssgvf("        Minimac3 Shell File:       %s/%s\n", *_opath, file_names[8]);
            data_loop(*_opath, file_names[8], "w");
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
	    fprintf_env_checkset_csh(_filep, "_MINIMAC", "minimac");
#endif /* RUNSHELL_SETUP */
        }
        void inner() {
            char cmd1[2*FILENAME_LENGTH];
            char cmd2[2*FILENAME_LENGTH];

            sprintf(cmd1, "%s/%s", "SHAPEIT", "shapeit");
            sh_find_pgm("SHAPEIT", cmd1, "shapeit");
            sprintf(cmd1, "$%s_program ", "shapeit");

            //we can make the split between minimac3 and minimac3-omp here instead
            if (g_cpus == 1) {
                sprintf(cmd2, "%s/%s", "MINIMAC3", "Minimac3");
                sh_find_pgm("MINIMAC3", cmd2, "Minimac3");
                sprintf(cmd2, "$%s_program ", "Minimac3");
            }

            if (g_cpus > 1) {
                sprintf(cmd2, "%s/%s", "MINIMAC3OMP", "Minimac3-omp");
                sh_find_pgm("MINIMAC3OMP", cmd2, "Minimac3-omp");
                sprintf(cmd2, "$%s_program ", "Minimac3-omp");
            }

            Vecs s_hapsplit;
            split(s_hapsplit, s_haplotype_file,"?");
            Vecs legsplit;
            split(legsplit, legend_file,"?");
            Vecs mapsplit;
            split(mapsplit, map_file,"?");
            Vecs m_hapsplit;
            split(m_hapsplit, m_haplotype_file,"?");


            //need a split for if the user has selected the HAPS, SAMPLE, and LEGEND files for shapeit or not
            //something for users to note, without the haps, sample file for shapeit, there is no prebuilt exclusion (so I removed it from that corresponding script)
            //this section will run print out the shapeit check and phase with arguments depending on whether or not we're using HAPS etc. for shapeit
            //If it's trying to use a file it checks to make sure it can insert the chromosome and if not errors
            
            if(haps_sample_selected) {
                if (s_hapsplit.size() == 2 && legsplit.size() == 2 && mapsplit.size() == 2) {
                    pr_nl();
                    pr_printf("#use mega2 plink formatted output to run shapeit checks\n");
                    pr_printf("%s -check --input-bed %s %s %s --input-map %s/%s%d%s --input-ref %s/%s%d%s %s/%s%d%s %s/%s --output-log Chr%d.checks",
                              cmd1, file_names[3], file_names[1], file_names[0],
                              genetic_map_directory.c_str(), mapsplit[0].c_str(), _numchr, mapsplit[1].c_str(),
                              shapeit_panel_directory.c_str(), s_hapsplit[0].c_str(), _numchr, s_hapsplit[1].c_str(),
                              shapeit_panel_directory.c_str(), legsplit[0].c_str(), _numchr, legsplit[1].c_str(),
                              shapeit_panel_directory.c_str(), sample_file.c_str(), _numchr);
                    pr_nl();
                    pr_nl();

                    pr_printf("#use mega2 plink formatted output to run shapeit phase mode\n");
                    pr_printf("%s --input-bed %s %s %s --input-ref %s/%s%d%s %s/%s%d%s %s/%s --exclude-snp Chr%d.checks.snp.strand.exclude -O Chr%d.Phased.Output --thread %d",
                              cmd1, file_names[3], file_names[1], file_names[0],
                              shapeit_panel_directory.c_str(), s_hapsplit[0].c_str(), _numchr, s_hapsplit[1].c_str(),
                              shapeit_panel_directory.c_str(), legsplit[0].c_str(), _numchr, legsplit[1].c_str(),
                              shapeit_panel_directory.c_str(), sample_file.c_str(), _numchr, _numchr, g_cpus);
                    pr_nl();
                    pr_nl();
                }
                else
                    errorf("Error in input of Shapeit reference files: haplotype, legend, and map.\nMake sure there to include one (and only one) '?' in the file name to be replaced with the chromosome number.\n");
            }
            else
            {
                if (mapsplit.size() == 2) {
                    pr_nl();
                    pr_printf("#use mega2 plink formatted output to run shapeit checks\n");
                    pr_printf("%s -check --input-bed %s %s %s --input-map %s/%s%d%s  --output-log Chr%d.checks", cmd1, file_names[3], file_names[1], file_names[0], genetic_map_directory.c_str(), mapsplit[0].c_str(), _numchr, mapsplit[1].c_str(), _numchr);
                    pr_nl();
                    pr_nl();

                    pr_printf("#use mega2 plink formatted output to run shapeit phase mode\n");
                    pr_printf("%s --input-bed %s %s %s -O Chr%d.Phased.Output --thread %d", cmd1, file_names[3], file_names[1], file_names[0], _numchr, g_cpus);
                    pr_nl();
                    pr_nl();
                }
                else
                    errorf("Error in input of Shapeit reference files: map.\nMake sure there to include one (and only one) '?' in the filename to be replaced with the chromosome number.\n");

            }

            //finally we have to convert to VCF for a usable input for Minimac3
            pr_printf ("#use mega2 plink formatted output to run shapeit phase mode\n");
            pr_printf ("%s -convert --input-haps Chr%d.Phased.Output --output-vcf Chr%d.Phased.Output.vcf",cmd1,_numchr,_numchr);
            pr_nl();
            pr_nl();

            //And last of all we run Minimac3 imputation
            if (m_hapsplit.size() == 2) {
                if (g_cpus == 1)
                    pr_printf("%s --refHaps %s/%s%d%s --haps Chr%d.Phased.Output.vcf --prefix Chr%d.Imputed.Output --chr %d\n",
                            cmd2, minimac_panel_directory.c_str(), m_hapsplit[0].c_str(), _numchr, m_hapsplit[1].c_str(), _numchr, _numchr, _numchr);
                if (g_cpus > 1)
                    pr_printf("%s --refHaps %s/%s%d%s --haps Chr%d.Phased.Output.vcf --prefix Chr%d.Imputed.Output --chr %d --cpus %d\n",
                            cmd2, minimac_panel_directory.c_str(), m_hapsplit[0].c_str(), _numchr, m_hapsplit[1].c_str(), _numchr, _numchr, _numchr, g_cpus);
            }
            else
                errorf("Error in input or Minimac3 reference files: haplotype.\n Make sure there to include one (and only one) '?' in the filename to be replaced with the chromosome number.\n");

            pr_nl();

        }

        //finds the program to run dynamically and gives an error if it can't be found.
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
            pr_printf("  echo\n");
            pr_printf("  echo \"For further details, please see SHAPEIT/Minimac3 section of the Mega2 documentation.\"\n");
            pr_printf("  exit 0\n");
            pr_printf("endif\n");
            pr_nl();
        }

    } *minimac_shs = new minimac_sh(Top);

    minimac_shs->file_names = file_names;

    minimac_shs->sh         = sh;
    minimac_shs->iterate();

    if (top_shell) {
        mssgvf("        Minimac3 Top Shell File:   %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("        The above shell runs all shells.\n");

        sh->filep_close();
        delete sh;
    }

    delete minimac_shs;
}

//Since this is a PLINK output we needed to call this function while outputting shell files
void CLASS_MINIMAC::create_sh_file(linkage_ped_top *Top, char **file_names, const int numchr) {

    char prefix[100];

    sub_prog_name(_suboption, prefix);

    switch (_suboption) {
        case 0:
        case 1:  write_MINIMAC_sh(Top, file_names);   break;
        default: break;
    }
}

//Override of user_queries method, don't want to query the user about combining chromosomes, genotyping summaries, or filename stem (again)
//Since I wrote a different option menu function before realizing plink core had this method to override I realize now that it's easier to have a blank method in it's place and hide it's behavior
void CLASS_MINIMAC::user_queries(char **file_names_array, int *combine_chromo, int *create_summary) {
    *combine_chromo = 0;
    LoopOverChrm = 1;
    //do nothing since we have another menu.
}


/*
 * Renders the following menu:
 *
Shapeit/Minimac3 Analysis Menu:
0) Done with this menu - please proceed
1) File name stem:                                             minimac
2) Number of CPUS for Shapeit Prephasing/Minimac3 Imputation:  1
3) Choose Minimac3 reference panel directory:                  .
4) Choose Minimac3 reference sample file (in VCF/M3VCF):       ?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz
5) Choose Shapeit genetic recombination map directory:         .
6) Choose Shapeit genetic recombination map file:              genetic_map_chr?_combined_b37.txt
7) Use reference panel in HAPS/SAMPLE format for shapeit?      Yes
8) Choose Shapeit reference panel directory:                   .
9) Choose Shapeit reference haplotype file:                    1000GP_Phase3_chr?.hap.gz
10) Choose Shapeit reference legend file:                      1000GP_Phase3_chr?.legend.gz
11) Choose Shapeit reference sample file:                      1000GP_Phase3.sample
Enter selection: 0 - 11 >
 Initially the 8-11 options are hidden until you select to use the HAPS SAMPLE format
 */
void CLASS_MINIMAC::minimac_option_menu (char *file_names[], char *prefix){
    int done, s_hap, cpus, choice, s_hap_renamed, legend_renamed, stem, legend, sample, map, map_renamed, m_hap, m_hap_renamed, ref_toggle, map_dir, s_ref_dir, m_ref_dir;
    done = 0;
    stem = 1;
    map_dir = 5;
    map = 6;
    m_ref_dir = 3;
    s_hap = 9;
    legend = 10;
    sample = 11;
    s_ref_dir = 8;
    m_hap = 4;
    cpus = 2;
    ref_toggle = 7;
    choice = -1;
    map_renamed = 0;
    s_hap_renamed = 0;
    legend_renamed = 0;
    m_hap_renamed = 0;
    haps_sample_selected = 0;

    strcpy(prefix, file_name_stem);

    char map_file_input[255];
    char s_hap_file_input[255];
    char legend_file_input[255];
    char sample_file_input[255];
    char m_hap_file_input[255];
    char map_dir_input[255];
    char s_dir_input[255];
    char m_dir_input[255];
    char response[255];

    Vecs mapsplit;
    Vecs s_hapsplit;
    Vecs legendsplit;
    Vecs m_hapsplit;

    //check environment variables or set to default.
    if(getenv("minimac_reference_panel_directory")!=NULL)
        minimac_directory_name = getenv("minimac_reference_panel_directory");
    else
        minimac_directory_name = ".";

    if(getenv("shapeit_reference_panel_directory")!=NULL)
        shapeit_directory_name = getenv("shapeit_reference_panel_directory");
    else
        shapeit_directory_name = ".";

    if(getenv("shapeit_reference_map_directory")!=NULL)
        map_directory_name = getenv("shapeit_reference_map_directory");
    else
        map_directory_name = ".";

    if(getenv("minimac_reference_haplotype_file")!=NULL) {
        minimac_reference_haplotype_file = getenv("minimac_reference_haplotype_file");
        split(m_hapsplit, minimac_reference_haplotype_file, "?");
        m_haplotype_pre = m_hapsplit[0];
        m_haplotype_post = m_hapsplit[1];
    }
    else{
        m_haplotype_pre = "1000GP_Phase3_chr";
        m_haplotype_post = ".hap.gz";
    }

    if(getenv("shapeit_reference_haplotype_file")!=NULL) {
        shapeit_reference_haplotype_file = getenv("shapeit_reference_haplotype_file");
        split(s_hapsplit, shapeit_reference_haplotype_file, "?");
        s_haplotype_pre = s_hapsplit[0];
        s_haplotype_post = s_hapsplit[1];
    }
    else{
        s_haplotype_pre = "1000GP_Phase3_chr";
        s_haplotype_post = ".hap.gz";
    }

    if(getenv("shapeit_reference_map_file")!=NULL) {
        map_file = getenv("shapeit_reference_map_file");
        split(mapsplit, map_file, "?");
        map_pre = mapsplit[0];
        map_post = mapsplit[1];
    }
    else{
        map_pre = "genetic_map_chr";
        map_post = "_combined_b37.txt";
    }

    if(getenv("shapeit_reference_legend_file")!=NULL) {
        legend_file = getenv("shapeit_reference_legend_file");
        split(legendsplit, legend_file, "?");
        legend_pre = legendsplit[0];
        legend_post = legendsplit[1];
    }
    else{
        legend_pre = "1000GP_Phase3_chr";
        legend_post = ".legend.gz";
    }

    if(getenv("shapeit_reference_sample_file")!=NULL)
        reference_sample_file = getenv("shapeit_reference_sample_file");
    else
        reference_sample_file = "1000GP_Phase3.sample";


    while (choice != 0) {
        draw_line();
        printf("Shapeit/Minimac3 Analysis Menu:\n");
        printf("%d) Done with this menu - please proceed\n",                              done);
        printf("%d) File name stem:                                             %-15s\n", stem, prefix);
        printf("%d) Number of CPUS for Shapeit Prephasing/Minimac3 Imputation:  %d\n",    cpus, g_cpus);
        printf("%d) Choose Minimac3 reference panel directory:                  %s\n",    m_ref_dir, minimac_directory_name.c_str());
        printf("%d) Choose Minimac3 reference sample file (in VCF/M3VCF):       %s?%s\n", m_hap, m_haplotype_pre.c_str(),m_haplotype_post.c_str());
        printf("%d) Choose Shapeit genetic recombination map directory:         %s\n",    map_dir, map_directory_name.c_str());
        printf("%d) Choose Shapeit genetic recombination map file:              %s?%s\n", map, map_pre.c_str(), map_post.c_str());
        if(!haps_sample_selected)
            printf("%d) Use reference panel in HAPS/SAMPLE format for shapeit?      No   \n", ref_toggle);
        else{
            printf("%d) Use reference panel in HAPS/SAMPLE format for shapeit?      Yes  \n", ref_toggle);
            printf("%d) Choose Shapeit reference panel directory:                   %s\n",    s_ref_dir, shapeit_directory_name.c_str());
            printf("%d) Choose Shapeit reference haplotype file:                    %s?%s\n", s_hap, s_haplotype_pre.c_str(), s_haplotype_post.c_str());
            printf("%d) Choose Shapeit reference legend file:                      %s?%s\n", legend, legend_pre.c_str(), legend_post.c_str());
            printf("%d) Choose Shapeit reference sample file:                      %s\n",    sample, reference_sample_file.c_str());
        }



        if(!haps_sample_selected)
            printf("Enter selection: 0 - %d > ",7);
        else
            printf("Enter selection: 0 - %d > ",11);

        fcmap(stdin,"%d", &choice); newline;

        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }

        else if ( choice == done ) {
            if (!map_renamed) {
                if (getenv("SHAPEIT_MAP") != NULL) {
                    strcpy(map_file_input,getenv("SHAPEIT_MAP"));
                    reference_map_file = getenv("SHAPEIT_MAP");
                }
                else {
                    strcpy(map_file_input, "genetic_map_chr?_combined_b37.txt");
                    reference_map_file = "genetic_map_chr?_combined_b37.txt";
                }
            }
            if (!s_hap_renamed){
                if(getenv("SHAPEIT_REF_HAP")!=NULL) {
                    strcpy(s_hap_file_input, getenv("SHAPEIT_REF_HAP"));
                    shapeit_reference_haplotype_file = getenv("SHAPEIT_REF_HAP");
                }
                else{
                    strcpy(s_hap_file_input, "1000GP_Phase3_chr?.hap.gz");
                    shapeit_reference_haplotype_file = "1000GP_Phase3_chr?.hap.gz";
                }

            }


            if (!legend_renamed){
                if(getenv("SHAPEIT_LEGEND")!=NULL) {
                    strcpy(legend_file_input, getenv("SHAPEIT_LEGEND"));
                    reference_legend_file = getenv("SHAPEIT_LEGEND");
                }
                else{
                    strcpy(legend_file_input, "1000GP_Phase3_chr?.legend.gz");
                    reference_legend_file = "1000GP_Phase3_chr?.legend.gz";
                }
            }

            if (!m_hap_renamed){
                if(getenv("MINIMAC_REF_HAP")!=NULL){
                    strcpy(m_hap_file_input, getenv("MINIMAC_REF_HAP"));
                    minimac_reference_haplotype_file = getenv("MINIMAC_REF_HAP");
                }
                else{
                    strcpy(m_hap_file_input, "?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz");
                    minimac_reference_haplotype_file = "?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz";
                }
            }

            map_file = reference_map_file;
            s_haplotype_file = shapeit_reference_haplotype_file;
            legend_file = reference_legend_file;
            sample_file = reference_sample_file;
            m_haplotype_file = minimac_reference_haplotype_file;

            genetic_map_directory = map_directory_name;
            shapeit_panel_directory = shapeit_directory_name;
            minimac_panel_directory = minimac_directory_name;


            BatchValueSet (reference_map_file, "shapeit_reference_map_file");
            BatchValueSet (shapeit_reference_haplotype_file, "shapeit_reference_haplotype_file");
            BatchValueSet (reference_legend_file, "shapeit_reference_legend_file");
            BatchValueSet (reference_sample_file, "shapeit_reference_sample_file");
            BatchValueSet (haps_sample_selected, "shapeit_haps_file_selected");
            BatchValueSet (minimac_reference_haplotype_file, "minimac_reference_haplotype_file");
            BatchValueSet (map_directory_name, "shapeit_reference_map_directory");
            BatchValueSet (shapeit_directory_name, "shapeit_reference_panel_directory");
            BatchValueSet (minimac_directory_name, "minimac_reference_panel_directory");


            BatchValueSet (g_cpus, "batch_cpu_count");

            free(file_name_stem);
            file_name_stem = strdup(prefix);
            BatchValueSet(file_name_stem, "file_name_stem");
        }

        else if ( choice == stem ) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);
            newline;
            inner_file_names(file_names, "", prefix);
        }

        else if( choice == map) {
            while (1) {
                printf("Enter genetic recombination map file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");

                fcmap(stdin, "%s", &map_file_input);
                newline;
                split(mapsplit, map_file_input, "?");

                if (mapsplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }

                map_pre = mapsplit[0];
                map_post = mapsplit[1];
                reference_map_file = map_file_input;
                map_renamed =1;
                break;
            }
        }

        else if ( choice == s_hap ) {
            while (1) {
                printf("Enter Shapeit reference haplotype file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");

                fcmap(stdin, "%s", &s_hap_file_input);
                newline;
                split(s_hapsplit, s_hap_file_input, "?");

                if (s_hapsplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }

                s_haplotype_pre = s_hapsplit[0];
                s_haplotype_post = s_hapsplit[1];
                shapeit_reference_haplotype_file = s_hap_file_input;
                s_hap_renamed =1;
                break;
            }
        }

        else if (choice == legend) {
            while (1) {
                printf("Enter Shapeit Reference legend file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");

                fcmap(stdin, "%s", &legend_file_input);
                newline;
                split(legendsplit, legend_file_input, "?");

                if (legendsplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }

                legend_pre = legendsplit[0];
                legend_post = legendsplit[1];
                reference_legend_file = legend_file_input;
                legend_renamed =1;
                break;
            }
        }

        else if (choice == sample) {
            while (1) {
                printf("Enter Shapeit Reference legend file name >\n");

                fcmap(stdin, "%s", &sample_file_input);
                newline;
                reference_sample_file = sample_file_input;
                break;
            }
        }

        else if ( choice == m_hap ) {
            while (1) {
                printf("Enter Minimac3 reference haplotype file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");

                fcmap(stdin, "%s", &m_hap_file_input);
                newline;
                split(m_hapsplit, m_hap_file_input, "?");

                if (m_hapsplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }

                m_haplotype_pre = m_hapsplit[0];
                m_haplotype_post = m_hapsplit[1];
                minimac_reference_haplotype_file = m_hap_file_input;
                m_hap_renamed =1;
                break;
            }
        }

        else if ( choice == cpus ) {
            while (1) {
                printf("Number of CPUS for Minimac3 imputation > ");
                fcmap(stdin, "%d", &g_cpus);
                newline;

                if (g_cpus < 1){
                    printf("Please enter a valid number of cpus\n");
                    continue;
                }
                break;
            }
        }

        else if (choice == ref_toggle)
        {
            if (!haps_sample_selected)
                haps_sample_selected = 1;

            else
                haps_sample_selected = 0;
        }

        else if (choice == map_dir){
            while (1) {
                bool no = 0;
                printf("Enter Shapeit genetic recombination map directory>\n");
                fcmap(stdin, "%s", &map_dir_input);
                if(!is_dir(map_dir_input)){
                    while (1){
                        printf("This directory does not exist on the current machine, are you sure you want to continue with this value? (Y/N)");
                        fcmap(stdin, "%s", &response);
                        if(response[0] =='Y' || response[0] == 'y')
                            break;
                        else if (response[0] =='N' || response[0] == 'n'){
                            no = 1;
                            break;
                        }
                        else{
                            printf("Unknown response\n");
                            continue;
                        }
                    }
                }
                newline;
                if(!no)
                    map_directory_name = map_dir_input;
                break;
            }
        }
        else if (choice == s_ref_dir){
            while (1) {
                bool no = 0;
                printf("Enter Shapeit reference panel directory>\n");
                fcmap(stdin, "%s", &s_dir_input);
                if(!is_dir(s_dir_input)) {
                    while (1){
                        printf("This directory does not exist on the current machine, are you sure you want to continue with this value? (Y/N)");
                        fcmap(stdin, "%s", &response);
                        if(response[0] =='Y' || response[0] == 'y')
                            break;
                        else if (response[0] =='N' || response[0] == 'n'){
                            no = 1;
                            break;
                        }
                        else{
                            printf("Unknown response\n");
                            continue;
                        }
                    }
                }

                newline;
                if(!no)
                    shapeit_directory_name = s_dir_input;
                break;
            }
        }
        else if( choice == m_ref_dir){
            while (1) {
                bool no = 0;
                printf("Enter Minimac3 reference panel directory>\n");
                fcmap(stdin, "%s", &m_dir_input);
                if(!is_dir(m_dir_input)) {
                    while (1){
                        printf("This directory does not exist on the current machine, are you sure you want to continue with this value? (Y/N)");
                        fcmap(stdin, "%s", &response);
                        if(response[0] =='Y' || response[0] == 'y')
                            break;
                        else if (response[0] =='N' || response[0] == 'n'){
                            no = 1;
                            break;
                        }
                        else {
                            printf("Unknown response\n");
                            continue;
                        }
                    }
                }

                newline;
                if(!no)
                    minimac_directory_name = m_dir_input;
                break;
            }
        }


        else {
            printf("Unknown option %d\n", choice);
            newline;

        }
    }
}

//change... at this point just the shell file names
static void inner_file_names(char **file_names, const char *num, const char *stem) {
    //sprintf(file_names[10], "%s_snps.%s", stem,num);
    sprintf(file_names[8], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.top.sh", stem);
}

//need to call inner file names for the names we change
//need to call file_names_w_stem for PLINK output files
void CLASS_MINIMAC::gen_file_names(char **file_names, char *num){
    file_names_w_stem(file_names, num, file_name_stem, PLINK_SUB_OPTION_SNP_MAJOR_INT);
    inner_file_names(file_names, num);
}

//Replace the chromosomes for the PLINK files, bed, bim, fam
void CLASS_MINIMAC::replace_chr_number(char *file_names[], int numchr) {
    CLASS_PLINK_CORE::replace_chr_number(file_names, numchr);

    //change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[8], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);

}

/*
 * Using some old Batch items from the MaCH/Minimac format (batch, minimac reference
 * As well as new reference items for References fro prephasing with Shapeit
 * The seleccted flag indicates if the user wanted to use those references
 * the Map is required
 * then we have the Sample haplotype file indicators
 */
void CLASS_MINIMAC::batch_out()
{
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "shapeit_reference_map_file",
                       "shapeit_reference_haplotype_file",
                       "shapeit_reference_legend_file",
                       "shapeit_reference_sample_file",
                       "shapeit_haps_file_selected",
                       "minimac_reference_haplotype_file",
                       "batch_cpu_count",
                       "shapeit_reference_map_directory",
                       "shapeit_reference_panel_directory",
                       "minimac_reference_panel_directory",

    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_MINIMAC::batch_in()
{
    char *fn = this->file_name_stem;

    BatchValueIfSet(fn,   "file_name_stem");
    BatchValueGet(reference_map_file, "shapeit_reference_map_file");
    BatchValueGet(shapeit_reference_haplotype_file, "shapeit_reference_haplotype_file");
    BatchValueGet(reference_legend_file, "shapeit_reference_legend_file");
    BatchValueGet(reference_sample_file, "shapeit_reference_sample_file");
    BatchValueGet(haps_sample_selected, "shapeit_haps_file_selected");
    BatchValueGet(minimac_reference_haplotype_file,"minimac_reference_haplotype_file");
    BatchValueGet(g_cpus, "batch_cpu_count");
    BatchValueGet(map_directory_name, "shapeit_reference_map_directory");
    BatchValueGet(minimac_directory_name, "minimac_reference_panel_directory");
    BatchValueGet(shapeit_directory_name, "shapeit_reference_panel_directory");

    map_file = reference_map_file;
    sample_file = reference_sample_file;
    legend_file = reference_legend_file;
    s_haplotype_file = shapeit_reference_haplotype_file;
    m_haplotype_file = minimac_reference_haplotype_file;
    genetic_map_directory = map_directory_name;
    minimac_panel_directory = minimac_directory_name;
    shapeit_panel_directory = shapeit_directory_name;
}

//overriding these functions but we only have one "suboption"
void CLASS_MINIMAC::sub_prog_name(int sub_opt, char *subprog) {
    strcpy(subprog, "Check And Phase");
}

void CLASS_MINIMAC::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int selection = 1;
    (*analysis)->_suboption = selection;
}

void CLASS_MINIMAC::sub_prog_name_to_sub_option(char *subprog_name, analysis_type *analysis) {
    (*analysis)->_suboption = 1;
}

