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

#include "fcmap_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"

#include "write_minimacshapeit_ext.h"


static void write_MINIMAC_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MINIMAC_sh(linkage_ped_top *Top, char *file_names[]);

static void inner_file_names(char **file_names, const char *num, const char *stem = "minimac");

extern int g_cpus;
Str s_haplotype_file;
Str legend_file;
Str sample_file;
Str map_file;
Str m_haplotype_file;

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
    LoopOverChrm  = ! combine_chromo;

    LoopOverTrait = 0;
    num_traits = 0;

    //We need to check for only ACTG as alleles.
    const char *allele_name;
    for (int locus = Top->LocusTop->PhenoCnt; locus < Top->LocusTop->LocusCnt; locus++){
        for (int allele = 0; allele < Top->LocusTop->Locus[locus].AlleleCnt; allele++){
            allele_name = Top->LocusTop->Locus[locus].Allele[allele].AlleleName;
            //printf("%s\n",allele_name);
            if ( ! ((strcmp(allele_name,"A") == 0) || (strcmp(allele_name,"C") == 0) || (strcmp(allele_name,"G") == 0)|| (strcmp(allele_name,"T") == 0) || (strcmp(allele_name,"0") == 0) || (strcmp(allele_name,"dummy") == 0))) {
                char error[255];
                strcpy(error, "The MaCH Minimac3 pipeline requires alleles to be labeled as A,C,T,G.\nInvalid allele label: ");
                strcat(error, allele_name);
                errorf(error);
                EXIT(DATA_TYPE_ERROR);
            }
        }
    }


    omit_peds(untyped_ped_opt, Top);

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    (*analysis)->_suboption = PLINK_SUB_OPTION_SNP_MAJOR_INT;
    create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, file_name_stem, analysis);

    printf("Mega2 created the following file(s) for SHAPEIT/MINIMAC3:\n");
    write_MINIMAC_snps(Top, file_names, pwid, fwid);

    write_MINIMAC_sh(Top, file_names);
}

// Format of file:
// CHR:PhysicalMapDistance
static void write_MINIMAC_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
{
    vlpCLASS(minimac_snp,chr,loci) {
        vlpCTOR(minimac_snp,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("        Minimac3 SNP File:         %s/%s\n", *_opath, file_names[10]);
            data_loop(*_opath, file_names[10], "w");
        }

        void inner() {
            pr_printf("%d:",_numchr);
            pr_physical_distance(NULL);
            pr_nl();
        }

    } *minimac_snps = new minimac_snp(Top);

    minimac_snps->file_names = file_names;

    minimac_snps->load_formats(fwid, pwid, -1);

    minimac_snps->iterate();

    delete minimac_snps;
}

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
            mssgvf("        Minimac Shell File:        %s/%s\n", *_opath, file_names[8]);
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
            //char cmd2[2*FILENAME_LENGTH];
            char cmd3[2*FILENAME_LENGTH];

            //get our program names from sh_find_pgm instead of having them static and hoping they're on the path
            sprintf(cmd1, "%s/%s", "SHAPEIT", "shapeit");
            sh_find_pgm("SHAPEIT", cmd1, "shapeit");
            sprintf(cmd1, "$%s_program ", "shapeit");

            //we can make the split between minimac3 and minimac3-omp here instead
            if (g_cpus == 1) {
                sprintf(cmd3, "%s/%s", "MINIMAC3", "minimac3");
                sh_find_pgm("MINIMAC3", cmd3, "minimac3");
                sprintf(cmd3, "$%s_program ", "minimac3");
            }

            if (g_cpus > 1) {
                sprintf(cmd3, "%s/%s", "MINIMAC3OMP", "minimac3-omp");
                sh_find_pgm("MINIMAC3OMP", cmd3, "minimac3-omp");
                sprintf(cmd3, "$%s_program ", "minimac3-omp");
            }

            Vecs s_hapsplit;
            split(s_hapsplit, s_haplotype_file,"?");
            Vecs legsplit;
            split(legsplit, legend_file,"?");
            Vecs mapsplit;
            split(mapsplit, map_file,"?");
            Vecs m_hapsplit;
            split(m_hapsplit, m_haplotype_file,"?");


            //first we run check to get snps to exclude
            pr_nl();
            pr_printf ("#use mega2 plink formatted output to run shapeit checks\n");
            pr_printf ("%s -check --input-bed %s %s %s --input-map %s%d%s --input-ref %s%d%s %s%d%s %s --output-log Chr%d.checks"
                    ,cmd1,file_names[3],file_names[1],file_names[0],mapsplit[0].c_str(),_numchr,mapsplit[1].c_str(),s_hapsplit[0].c_str(),_numchr,s_hapsplit[1].c_str(),legsplit[0].c_str(),_numchr,legsplit[1].c_str(), sample_file.c_str(),_numchr);
            pr_nl();
            pr_nl();


            //next we run the shapeit phasing mode
            pr_printf ("#use mega2 plink formatted output to run shapeit phase mode\n");
            pr_printf ("%s --input-bed %s %s %s --input-ref %s%d%s %s%d%s %s --exclude-snp Chr%d.checks.snp.strand.exclude -O Chr%d.Phased.Output --thread %d"
                    ,cmd1,file_names[3],file_names[1],file_names[0],s_hapsplit[0].c_str(),_numchr,s_hapsplit[1].c_str(),legsplit[0].c_str(),_numchr,legsplit[1].c_str(), sample_file.c_str(),_numchr,_numchr,g_cpus);
            pr_nl();
            pr_nl();

            //finally we have to convert to VCF for a usable input for Minimac3
            pr_printf ("#use mega2 plink formatted output to run shapeit phase mode\n");
            pr_printf ("%s -convert --input-haps Chr%d.Phased.Output --output-vcf Chr%d.Phased.Output.vcf"
                    ,cmd1,_numchr,_numchr);
            pr_nl();
            pr_nl();

            //And last of all we run Minimac3 imputation
            if (g_cpus == 1)
                pr_printf ("%s --refHaps %s%d%s --haps Chr%d.Phased.Output.vcf --prefix Chr%d.Imputed.Output --chr %d\n",cmd3, m_hapsplit[0].c_str(),_numchr,m_hapsplit[1].c_str(),_numchr,_numchr,_numchr);
            if (g_cpus > 1)
                pr_printf ("%s --refHaps %s%d%s --haps Chr%d.Phased.Output.vcf --prefix Chr%d.Imputed.Output --chr %d --cpus %d\n",cmd3, m_hapsplit[0].c_str(),_numchr,m_hapsplit[1].c_str() ,_numchr,_numchr,_numchr,g_cpus);

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
            pr_printf("  echo \"For further details, please see MaCH/Minimac3 section of the Mega2 documentation.\"\n");
            pr_printf("  exit 0\n");
            pr_printf("endif\n");
            pr_nl();
        }

    } *minimac_shs = new minimac_sh(Top);

    minimac_shs->file_names = file_names;

    minimac_shs->sh         = sh;
    minimac_shs->iterate();

    if (top_shell) {
        mssgvf("        MaCH Top Shell File:       %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("        The above shell runs all shells.\n");

        sh->filep_close();
        delete sh;
    }

    delete minimac_shs;
}

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
    combine_chromo = 0;
    //do nothing since we have another menu.
}


//need to switch some of the inputs for shapeit
void CLASS_MINIMAC::minimac_option_menu (char *file_names[], char *prefix){
    int done, s_hap, cpus, choice, s_hap_renamed, legend_renamed, stem, legend, sample, map, map_renamed, m_hap, m_hap_renamed, ref_toggle, toggled;
    done = 0;
    stem = 1;
    map = 4;
    s_hap = 6;
    legend = 7;
    sample = 8;
    m_hap = 3;
    cpus = 2;
    ref_toggle = 5;
    choice = -1;
    map_renamed = 0;
    s_hap_renamed = 0;
    legend_renamed = 0;
    m_hap_renamed = 0;
    toggled = 0;

    strcpy(prefix, file_name_stem);

    char map_file_input[255];
    char s_hap_file_input[255];
    char legend_file_input[255];
    char sample_file_input[255];
    char m_hap_file_input[255];
    Vecs mapsplit;
    Vecs s_hapsplit;
    Vecs legendsplit;
    Vecs m_hapsplit;
    map_pre = "genetic_map_chr";
    map_post = "_combined_b37.txt";
    s_haplotype_pre = "1000GP_Phase3_chr";
    s_haplotype_post = ".hap.gz";
    legend_pre = "1000GP_Phase3_chr";
    legend_post = ".legend.gz";
    reference_sample_file = "1000GP_Phase3.sample";
    m_haplotype_pre = "";
    m_haplotype_post = ".1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz";

    while (choice != 0) {
        draw_line();
        printf("Shapeit/Minimac3 Analysis Menu:\n");
        printf("%d) Done with this menu - please proceed\n",                              done);
        printf("%d) File name stem:                                             %-15s\n", stem, prefix);
        printf("%d) Number of CPUS for Shapeit Prephasing/Minimac3 Imputation:  %d\n",    cpus, g_cpus);
        printf("%d) Choose Minimac3 reference sample file:                      %s?%s\n", m_hap, m_haplotype_pre.c_str(),m_haplotype_post.c_str());
        printf("%d) Choose Shapeit map file:                                    %s?%s\n", map, map_pre.c_str(), map_post.c_str());
        if(!toggled)
            printf("%d) Use reference panel in HAPS/SAMPLE format for shapeit?      No   \n", ref_toggle);
        else{
            printf("%d) Use reference panel in HAPS/SAMPLE format for shapeit?      Yes  \n", ref_toggle);
            printf("%d) Choose Shapeit reference haplotype file:                    %s?%s\n", s_hap, s_haplotype_pre.c_str(), s_haplotype_post.c_str());
            printf("%d) Choose Shapeit reference legend file:                       %s?%s\n", legend, legend_pre.c_str(), legend_post.c_str());
            printf("%d) Choose Shapeit reference sample file:                       %s\n",    sample, reference_sample_file.c_str());
        }




        printf("Enter selection: 0 - %d > ",8);

        fcmap(stdin,"%d", &choice); newline;

        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }

        else if ( choice == done ) {
            if (!map_renamed){
                strcpy(map_file_input, "genetic_map_chr?_combined_b37.txt");
                reference_map_file = "genetic_map_chr?_combined_b37.txt";
            }
            if (!s_hap_renamed){
                strcpy(s_hap_file_input, "1000GP_Phase3_chr?.hap.gz");
                shapeit_reference_haplotype_file = "1000GP_Phase3_chr?.hap.gz";
            }


            if (!legend_renamed){
                strcpy(legend_file_input, "1000GP_Phase3_chr?.legend.gz");
                reference_legend_file = "1000GP_Phase3_chr?.legend.gz";
            }

            if (!m_hap_renamed){
                strcpy(m_hap_file_input, "?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz");
                minimac_reference_haplotype_file = "?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz";
            }
            map_file = reference_map_file;
            s_haplotype_file = shapeit_reference_haplotype_file;
            legend_file = reference_legend_file;
            sample_file = reference_sample_file;
            m_haplotype_file = minimac_reference_haplotype_file;

            BatchValueSet (shapeit_reference_haplotype_file, "shapeit_reference_haplotype_file");
            BatchValueSet (reference_legend_file, "shapeit_reference_legend_file");
            BatchValueSet (reference_sample_file, "shapeit_reference_sample_file");
            BatchValueSet (minimac_reference_haplotype_file, "minimac_reference_haplotype_file");

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
                printf("Enter map file name >\n");
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
            if (toggled == 0)
                toggled = 1;
            else
                toggled = 0;
        }

        else {
            printf("Unknown option %d\n", choice);
            newline;

        }
    }
}

static void inner_file_names(char **file_names, const char *num, const char *stem) {
    sprintf(file_names[10], "%s_snps.%s", stem,num);
    sprintf(file_names[8], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.top.sh", stem);
}

void CLASS_MINIMAC::gen_file_names(char **file_names, char *num){
    file_names_w_stem(file_names, num, file_name_stem,
                      PLINK_SUB_OPTION_SNP_MAJOR_INT);
    inner_file_names(file_names, num);
}

void CLASS_MINIMAC::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[8], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_MINIMAC::batch_out()
{
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "shapeit_reference_haplotype_file",
                       "shapeit_reference_legend_file",
                       "shapeit_reference_sample_file",
                       "minimac_reference_haplotype_file",
                       "batch_cpu_count",
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
    BatchValueIfSet(shapeit_reference_haplotype_file, "shapeit_reference_haplotype_file");
    BatchValueIfSet(reference_legend_file, "shapeit_reference_legend_file");
    BatchValueIfSet(reference_sample_file, "shapeit_reference_sample_file");
    BatchValueIfSet(minimac_reference_haplotype_file,"minimac_reference_haplotype_file");
    BatchValueGet(g_cpus, "batch_cpu_count");

    sample_file = reference_sample_file;
    legend_file = reference_legend_file;
    s_haplotype_file = shapeit_reference_haplotype_file;
    m_haplotype_file = minimac_reference_haplotype_file;
}

//overriding these functions but we only have one "suboption"
void CLASS_MINIMAC::sub_prog_name(int sub_opt, char *subprog) {
    strcpy(subprog, "CheckAndPhase");
}

void CLASS_MINIMAC::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int selection = 1;
    (*analysis)->_suboption = selection;
}

void CLASS_MINIMAC::sub_prog_name_to_sub_option(char *subprog_name, analysis_type *analysis) {
    (*analysis)->_suboption = 1;
}

