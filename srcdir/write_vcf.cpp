/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#include "write_vcf_ext.h"


void CLASS_VCF::create_output_files(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    linkage_ped_top *Top = LPedTreeTop;


    if ( InputMode == INTERACTIVE_INPUTMODE ) {

    }
    else {
        batch_in();
        //inner_file_names(file_names, "", file_name_stem);
    }

    //I believe everything goes into one VCF file
    int combine_chromo = 1;
    LoopOverChrm  = ! combine_chromo;

    //We'll want a phenotype file, not sure if we'll want to loop over traits or have them split out
    //LoopOverTrait = 0;

    omit_peds(untyped_ped_opt, Top);
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);


    write_VCF_file(Top,file_names,pwid,fwid);
    write_VCF_ped(Top, file_names,pwid,fwid);
    write_VCF_pheno(Top,file_names,pwid,fwid);
    write_VCF_sh(Top, file_names);
}

void CLASS_VCF::write_VCF_file(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_peds,both,ped_per) {
        vlpCTOR(vcf_peds,both,ped_per) { }

        void file_loop() {
            mssgvf("        VCF format file:      %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }

        //goal for inner loop:
        //#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
        //20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
        void inner() {
            pr_printf("%d\t", _numchr);
            pr_physical_distance(0);
            pr_marker_name();
            //pr_ref()
            //pr_alt()
            //pr_qual()
            //pr_filter
            //pr_info
            //pr_format
            //etc...

        }
    } *sp = new vcf_peds(Top);

    sp->setfln(prefix, ".vcf");

    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;
}

//this will write the pedigree in the PLINK fam file format
void CLASS_VCF::write_VCF_ped(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_peds,trait,ped_per) {
        vlpCTOR(vcf_peds,trait,ped_per) { }

        void file_loop() {
            mssgvf("        VCF pedigree file:      %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_father();
            pr_mother();
            pr_sex();
            pr_pheno();

        }
    } *sp = new vcf_peds(Top);

    sp->setfln(prefix, ".peds");

    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;
}

// this will write the phenotypic data in
void CLASS_VCF::write_VCF_pheno(linkage_ped_top *Top, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_phenos,trait,ped_per) {
        vlpCTOR(vcf_phenos,trait,ped_per) { }

        void file_loop() {
            mssgvf("        VCF phenotype file:     %s/%s\n", *_opath, _fln);
            data_loop(*_opath, _fln, "w");
        }
        void inner() {
            pr_id();
            pr_parent();
            pr_pheno();
            //what is sampleid?
            //pr_sampleid
            pr_nl();
        }
    } *sp = new vcf_phenos(Top);

    sp->setfln(prefix, ".phes");

    sp->load_formats(fwid, pwid, -1);

    sp->_trait_affect = true;
    sp->iterate();

    delete sp;
}

//this will write the output of the VCF file
void CLASS_VCF::write_VCF_sh(linkage_ped_top *Top, char *file_names[]) {

}

void CLASS_VCF::batch_in(){

}

void CLASS_VCF::batch_out(){

}