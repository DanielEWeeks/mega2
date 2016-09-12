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

#include "write_mach_ext.h"

void create_MACH_files(linkage_ped_top **LPedTop, char *file_names[], const int untyped_ped_opt, const int output_format);

static void write_MACH_peds(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MACH_data(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MACH_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid);

static void write_MACH_sh(linkage_ped_top *Top, char *file_names[]);

static void inner_file_names(char **file_names, const char *num, const char *stem = "mach");

int g_cpus =1;

extern int allele_count;
extern int ALLELE_ARRAY;
extern allele_prop **Allele_Array;
extern allele_prop *canonical_allele_internal(const char *v);

Str reference_directory = ".";


//  Output will like this:
//  FAM1001   ID1234  0   0   M  A A   A C   C C
//  FAM1002   ID1234  0   0   F  A C   C C   G G
//  Input should only have ACTG as alleles so the output should as well.
static void write_MACH_peds(linkage_ped_top *Top, char **file_names, const int pwid, const int fwid)
{

    vlpCLASS(mach_ped, chr, ped_per_loci) {
        vlpCTOR(mach_ped, chr, ped_per_loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("     MaCH Pedigree File:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }

        void per_start() {
            //pr_id prints both family id and person id looks like
            pr_id();
            pr_father();
            pr_mother();
            pr_sex_l();
        }

        void per_end() {
            pr_nl();
        }

        void inner() {
            pr_marker();
        }

    } *mach_peds = new mach_ped(Top);

    mach_peds->file_names = file_names;

    mach_peds->load_formats(fwid, pwid, -1);

    mach_peds->iterate();

    delete mach_peds;
}


/*<Example of a simple data file>
 *  M marker1
 *  M marker2
 *  ...
 *  <End of simple data file>
 *  */
static void write_MACH_data(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
{
    vlpCLASS(mach_dat,chr,loci) {
        vlpCTOR(mach_dat,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("     MaCH Data File:          %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "w");
        }

        void inner() {
            pr_printf("M ");

            //this works though
            char *markername = _tlocusp->Marker->MarkerName;
            pr_printf(markername);
            pr_nl();
        }

    } *mach_dats = new mach_dat(Top);

    mach_dats->file_names = file_names;

    mach_dats->load_formats(fwid, pwid, -1);

    mach_dats->iterate();

    delete mach_dats;
}


// Format of file:
// CHR:PhysicalMapDistance
// Needed to create this file for mach2VCF as an input
static void write_MACH_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
{
    vlpCLASS(mach_snp,chr,loci) {
        vlpCTOR(mach_snp,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("     MaCH SNP File:           %s/%s\n", *_opath, file_names[4]);
            data_loop(*_opath, file_names[4], "w");
        }

        void inner() {
            pr_printf("%d:",_numchr);
            pr_physical_distance(NULL);
            pr_nl();
        }

    } *mach_snps = new mach_snp(Top);

    mach_snps->file_names = file_names;

    mach_snps->load_formats(fwid, pwid, -1);

    mach_snps->iterate();

    delete mach_snps;
}


/*Shell should look like the following in order to produce prephased file according to minimac documentation
 *
 * mach1 -d Gwas.chr20.Unphased.dat \
      -p Gwas.chr20.Unphased.ped \
      --rounds 20 \
      --states 200 \
      --phase \
      --interim 5 \
      --sample 5 \
      --prefix Gwas.Chr20.Phased.Output

      It looks like we'll the shell to create a pipeline to go MaCH->VCF->minimac3.

      For this full pipeline users will need mach1, mach2vcf, and minimac3/minimac3-omp installed on their machines and in their path

      To get the vcf from the mach output we need a command like:
      mach2VCF --haps Gwas.Chr20.Phased.Output.hap \
         --snps Gwas.Chr20.Phased.Output.snps \
         --prefix Gwas.Chr20.Phased.Output.VCF.format

      Finally with a phased vcf we need to set up a script for minimac3 or minimac3-omp
      something like this:
      Minimac3 --refHaps ReferencePanel.Chr20.1000Genomes.vcf \
                --haps Gwas.Chr20.Phased.Output.VCF.format.vcf \
                --prefix Gwas.Chr20.Imputed.Output

      or alternatively it Minimac3 with Minimac3-omp and add a flag of --cpus # if more than one cpu is selected.
      */


static void write_MACH_sh(linkage_ped_top *Top, char *file_names[]) {

    int top_shell = 1;

    dataloop::sh_exec *sh = 0;
    if (top_shell) {
        sh = new dataloop::sh_exec(Top);
        sh->filep_open(output_paths[0], file_names[3], "w");
        sh->sh_main();
    }

    vlpCLASS(mach_sh, both, sh_exec) {
        vlpCTOR(mach_sh, both, sh_exec) { }

        void file_loop() {
            mssgvf("     MaCH Shell File:         %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
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
	    fprintf_env_checkset_csh(_filep, "_MACH", "mach");
#endif /* RUNSHELL_SETUP */
        }
        void inner() {
            char cmd1[2*FILENAME_LENGTH];
            char cmd2[2*FILENAME_LENGTH];
            char cmd3[2*FILENAME_LENGTH];

            //get our program names from sh_find_pgm instead of having them static and hoping they're on the path
            sprintf(cmd1, "%s/%s", "MACH1", "mach1");
            sh_find_pgm("MACH1", cmd1, "mach1");
            sprintf(cmd1, "$%s_program ", "mach1");

            sprintf(cmd2, "%s/%s", "MACH2VCF", "Mach2VCF");
            sh_find_pgm("MACH2VCF", cmd2, "Mach2VCF");
            sprintf(cmd2, "$%s_program ", "Mach2VCF");

            //we can make the split between minimac3 and minimac3-omp here instead
            if (g_cpus == 1) {
                sprintf(cmd3, "%s/%s", "MINIMAC3", "Minimac3");
                sh_find_pgm("MINIMAC3", cmd3, "Minimac3");
                sprintf(cmd3, "$%s_program ", "Minimac3");
            }

            if (g_cpus > 1) {
                sprintf(cmd3, "%s/%s", "MINIMAC3OMP", "Minimac3-omp");
                sh_find_pgm("MINIMAC3OMP", cmd3, "Minimac3-omp");
                sprintf(cmd3, "$%s_program ", "Minimac3-omp");
            }

            pr_nl();
            //mach1 run, might need additional parameters, specifically do we need hapmap and snps files or are those only used with mach2vcf
            pr_printf ("#use mega2 output in mach1 to prephasedata\n");
            pr_printf ("%s -d %s -p %s --rounds 20 --states 200 --phase --interim 5 --sample 5 --prefix Chr%d.Phased.Output", cmd1 , file_names[1], file_names[0],_numchr);
            pr_nl();
            pr_nl();

            //based on what I've read, we want the out put file from mach1 to be the haps file for mach2VCF, I'm assuming the user needs to input a snpfile using the menu
            pr_printf ("#use mach2VCF to create a VCF output of prephased data\n");
            pr_printf ("%s --haps Chr%d.Phased.Output.gz --snps %s --prefix Chr%d.Phased.Output.VCF.format",cmd2, _numchr,file_names[4],_numchr);
            pr_nl();
            pr_nl();


            Vecs hapsplit;
            split(hapsplit, file_names[5], "?");

            if(hapsplit.size() ==2) {
                if (g_cpus == 1)
                    pr_printf(
                            "%s --refHaps %s/%s%d%s --haps Chr%d.Phased.Output.VCF.format.vcf.gz --prefix Chr%d.Imputed.Output --chr %d\n",
                            cmd3, reference_directory.c_str(), hapsplit[0].c_str(), _numchr, hapsplit[1].c_str(), _numchr, _numchr, _numchr);

                if (g_cpus > 1)
                    pr_printf(
                            "%s --refHaps %s/%s%d%s --haps Chr%d.Phased.Output.VCF.format.vcf.gz --prefix Chr%d.Imputed.Output --chr %d --cpus %d\n",
                            cmd3, reference_directory.c_str(), hapsplit[0].c_str(), _numchr, hapsplit[1].c_str(), _numchr, _numchr, _numchr, g_cpus);
            }
            else
                errorf("Error in input of Minimac3 reference files: haplotype.\n Make sure there to include one (and only one) '?' in the filename to be replaced with the chromosome number.\n");

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

    } *mach_shs = new mach_sh(Top);

    mach_shs->file_names = file_names;

    mach_shs->sh         = sh;
    mach_shs->iterate();

    if (top_shell) {
        mssgvf("     MaCH Top Shell File:     %s/%s\n", output_paths[0], file_names[3]);
        mssgvf("     The above shell runs all shells.\n");

        sh->filep_close();
        delete sh;
    }

    delete mach_shs;
}

void CLASS_MACH::create_output_file(
        linkage_ped_top *LPedTreeTop,
        analysis_type *analysis,
        char *file_names[],
        int untyped_ped_opt,
        int *numchr,
        linkage_ped_top **Top2)
{
    int pwid, fwid, mwid;
    int combine_chromo = 0;
    char prefix[100];
    linkage_ped_top *Top = LPedTreeTop;


    //AnalyInputMode will be new variabl
    if ( InputMode == INTERACTIVE_INPUTMODE ) {
        //get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
        mach_option_menu(file_names,prefix);
        //for some reason getting two sets of batch outputs with this on?
        //batch_out();
    }
    else {
        batch_in();
        sprintf( file_names[5], "%s", mach_reference_haplotype_file.c_str());
        inner_file_names(file_names, "", file_name_stem);
    }

    combine_chromo = 0;
    LoopOverChrm  = ! combine_chromo;


    //There are no traits for MaCH/Minimac3
    LoopOverTrait = 0;
    num_traits = 0;

    //need to make sure we only have A,C,T,G for allele markers
    //allele_prop *current_allele;
    const char *allele_name;


    //so instead of looking at the allele array we have to look at every allele for every marker
    //the reason for this is the possibility that a phenotype can come up as a marker and be coded as a number
    //while mach doesn't do anything with the trait we still don't want to error for it
    for (int locus = Top->LocusTop->PhenoCnt; locus < Top->LocusTop->LocusCnt; locus++){
        for (int allele = 0; allele < Top->LocusTop->Locus[locus].AlleleCnt; allele++) {
            if (Top->LocusTop->Locus[locus].Allele[allele].AlleleName != NULL) {
            allele_name = Top->LocusTop->Locus[locus].Allele[allele].AlleleName;
                if(allele_name == 0){
                    char error[255];
                    strcpy(error, "The SHAPEIT Minimac3 pipeline requires alleles to be labeled as A,C,T,G.\nInvalid allele label: ");
                    strcat(error, allele_name);
                    errorf(error);
                    EXIT(DATA_TYPE_ERROR);
                }
                else {
                    if (!((strcmp(allele_name, "A") == 0) || (strcmp(allele_name, "C") == 0) || (strcmp(allele_name, "G") == 0) || (strcmp(allele_name, "T") == 0) || (strcmp(allele_name, "0") == 0) || (strcmp(allele_name, "dummy") == 0))) {
                        char error[255];
                        strcpy(error, "The MaCH Minimac3 pipeline requires alleles to be labeled as A,C,T,G.\nInvalid allele label: ");
                        strcat(error, allele_name);
                        errorf(error);
                        EXIT(DATA_TYPE_ERROR);
                    }
                }
            }
            else{
                errorf("The MaCH Minimac3 pipeline requires alleles to be labeled as A,C,T,G. Entry without an Allele Name was found.  Check to make sure you aren't using numeric Alleles as input.\n");
                EXIT(DATA_TYPE_ERROR);
            }
        }
    }

    // Omit pedigrees under certain circumstances, this was added to resolve a bug
    omit_peds(untyped_ped_opt, Top);

    printf("Mega2 created the following file(s) for MaCH/Minimac3:\n");

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    write_MACH_data(Top, file_names, pwid, fwid);

    write_MACH_peds(Top, file_names, pwid, fwid);

    write_MACH_snps(Top, file_names, pwid, fwid);

    write_MACH_sh(Top, file_names);


}

//this will create and parse the options
void CLASS_MACH::mach_option_menu (char *file_names[], char *prefix){
    int done, hap, cpus, choice, renamed, stem,dir;
    done = 0;
    stem = 1;
    dir = 2;
    hap = 3;
    cpus = 4;
    choice = -1;
    renamed = 0;



    strcpy(prefix, file_name_stem);

    char hap_file_input[255];
    char dir_input[255];
    char response[255];
    Vecs hapsplit;
    haplotype_post = ".1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz";
    Str directory = ".";

    while (choice != 0) {
        draw_line();
        printf("MaCH/Minimac3 Analysis Menu:\n");
        printf("%d) Done with this menu - please proceed\n",done);
        printf("%d) File name stem:                                     %-15s\n", stem, prefix);
        printf("%d) Choose reference haplotype directory:               %s\n", dir, directory.c_str());
        printf("%d) Choose reference haplotype file:                    %s?%s\n", hap, haplotype_pre.c_str(), haplotype_post.c_str());
        printf("%d) Number of CPUS for Minimac3 Imputation:             %d\n",    cpus, g_cpus);
        printf("Enter selection: 0 - %d > ",4);

        fcmap(stdin,"%d", &choice); newline;

        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }

        else if ( choice == done ) {
            if (!renamed)
                strcpy (hap_file_input, "?.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz");

            reference_directory = directory;

            BatchValueSet (g_cpus, "batch_cpu_count");

            free(file_name_stem);
            file_name_stem = strdup(prefix);
            BatchValueSet(file_name_stem, "file_name_stem");
            BatchValueSet(directory,"minimac_reference_panel_directory");
        }

        else if ( choice == stem ) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);
            newline;
            inner_file_names(file_names, "", prefix);
        }

        else if ( choice == hap ) {
            while (1) {
                printf("Enter reference haplotype file name >\n");
                printf(" Reserve space for the chromosome number with a ? > ");

                fcmap(stdin, "%s", &hap_file_input);
                newline;
                split(hapsplit, hap_file_input, "?");

                if (hapsplit.size() != 2) {
                    printf("Please include one and only one ? in the file name\n");
                    continue;
                }

                haplotype_pre = hapsplit[0];
                haplotype_post = hapsplit[1];
                renamed =1;
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

        else if ( choice == dir) {
            while (1) {
                bool no = 0;
                printf("Enter reference haplotype directory>\n");
                fcmap(stdin, "%s", &dir_input);
                if(!is_dir(dir_input)){
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
                    directory = dir_input;
                break;
            }
        }

        else {
            printf("Unknown option %d\n", choice);
        }
    }

    sprintf(file_names[5], "%s", hap_file_input);
    mach_reference_haplotype_file = hap_file_input;
    BatchValueSet (mach_reference_haplotype_file, "minimac_reference_haplotype_file");

}

static void inner_file_names(char **file_names, const char *num, const char *stem) {
    sprintf(file_names[0], "%s_ped.%s", stem, num);
    sprintf(file_names[1], "%s_data.%s", stem, num);
    sprintf(file_names[2], "%s.%s.sh", stem, num);
    sprintf(file_names[3], "%s.top.sh", stem);
    sprintf(file_names[4], "%s_snps.%s", stem,num);
}

void CLASS_MACH::gen_file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_MACH::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_MACH::batch_out()
{
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "minimac_reference_haplotype_file",
                       "batch_cpu_count",
                       "minimac_reference_panel_directory",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_MACH::batch_in()
{
    char *fn = this->file_name_stem;

    BatchValueIfSet(                fn,   "file_name_stem");
    BatchValueGet(mach_reference_haplotype_file, "minimac_reference_haplotype_file");
    BatchValueGet(reference_directory,"minimac_reference_panel_directory");
    BatchValueGet(g_cpus, "batch_cpu_count");
}
