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

void CLASS_MINIMAC::create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    linkage_ped_top *Top = LPedTreeTop;

    int combine_chromo = 0;
    LoopOverChrm  = ! combine_chromo;

    LoopOverTrait = 0;
    num_traits = 0;

    omit_peds(untyped_ped_opt, Top);

    create_PLINK_files(&LPedTreeTop, file_names, UntypedPedOpt, PLINK_SUB_OPTION_SNP_MAJOR_INT-1, "minimac", analysis);
    save_bed_file();


    for(int i = 0; i<16; i++){
        printf("%s\n",file_names[i]);
    }

    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    printf("Mega2 created the following file(s) for SHAPEIT/Minimac3:\n");
    write_MINIMAC_snps(Top, file_names, pwid, fwid);

    write_MINIMAC_sh(Top, file_names);
}

// Format of file:
// CHR:PhysicalMapDistance
// Needed to create this file for mach2VCF as an input
static void write_MINIMAC_snps(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid)
{
    vlpCLASS(minimac_snp,chr,loci) {
        vlpCTOR(minimac_snp,chr,loci) { }

        typedef char *str;
        str *file_names;

        void file_loop() {
            msgvf("     Minimac3 SNP File:           %s/%s\n", *_opath, file_names[10]);
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

}



//this will create and parse the options
void CLASS_MINIMAC::minimac_option_menu (char *file_names[], char *prefix){

}

static void inner_file_names(char **file_names, const char *num, const char *stem) {
    sprintf(file_names[10], "%s_snps.%s", stem,num);
}

void CLASS_MINIMAC::gen_file_names(char **file_names, char *num){
    inner_file_names(file_names, num);
}

void CLASS_MINIMAC::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[10], numchr);
}

void CLASS_MINIMAC::batch_out()
{

}

void CLASS_MINIMAC::batch_in()
{

}
