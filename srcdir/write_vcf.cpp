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
#include "write_files_ext.h"

#include "write_vcf_ext.h"


void CLASS_VCF::create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    linkage_ped_top *Top = LPedTreeTop;

//    if ( InputMode == INTERACTIVE_INPUTMODE ) {
//
//    }
//    else {
//        batch_in();
//        //inner_file_names(file_names, "", file_name_stem);
//    }

    //I believe everything goes into one VCF file
    int combine_chromo = 1;
    LoopOverChrm  = ! combine_chromo;

    LoopOverTrait = 0;

    //get_file_names(file_names, Top->OrigIds, Top->UniqueIds, &combine_chromo);

    //We'll want a phenotype file, not sure if we'll want to loop over traits or have them split out
    //LoopOverTrait = 0;

    omit_peds(untyped_ped_opt, Top);
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    file_name_stem = strdup("vcf");
    inner_file_names(file_names,"",file_name_stem);


    printf("Mega2 created the following file(s) for VCF Format:\n");
    write_VCF_file(Top, file_name_stem,file_names ,pwid, fwid);
    write_VCF_ped(Top, file_name_stem, file_names ,pwid, fwid);
    //we only want a phenotype file if we have more than one trait, the first trait is always put into the pedigree fam file by convention
    if(num_traits>1)
        write_VCF_pheno(Top, file_name_stem, file_names ,pwid, fwid);
    write_VCF_sh(Top, file_name_stem, file_names);
}


//CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1_1 1_2 2_1
//20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
void CLASS_VCF::write_VCF_file(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_vcfs_header,chr,ped_per) {
        vlpCTOR(vcf_vcfs_header, chr, ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF format file:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }

        void file_header() {
            pr_printf("##fileformat=VCFv4.2\n");
            pr_printf("##filedate=%s\n", __DATE__);
            pr_printf("##source=MEGA2\n");
            pr_printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n");
            pr_printf("##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
            pr_printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
            pr_printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            pr_printf("##FILTER=<ID=PASS,Description=\"Passed variant FILTERs\">\n");
//            pr_printf("%-7s","#CHROM");
//            pr_printf("%-10s","POS");
//            pr_printf("%-15s","ID");
//            pr_printf("%-6s","REF");
//            pr_printf("%-6s","ALT");
//            pr_printf("%-5s","QUAL");
//            pr_printf("%-7s","FILTER");
//            pr_printf("%-65s","INFO");
//            pr_printf("%-7s","FORMAT");
            pr_printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
        }

        void inner(){
            pr_fam();
            pr_printf("_");
            pr_per();
            pr_printf(" ");
        }

//        void filep_close(){
//            pr_nl();
//        }

    } *hlp = new vcf_vcfs_header(Top);

    hlp->setfln(prefix, ".vcf");
    hlp->file_names = file_names;

    hlp->load_formats_no_space(-1);
    //lp->load_formats(fwid, pwid, -1);

    hlp->iterate();



    vlpCLASS(vcf_vcfs,chr,loci_ped_per) {
        vlpCTOR(vcf_vcfs,chr,loci_ped_per) { }
        typedef char *str;
        str *file_names;
        int firstrun;
        linkage_locus_top *LTop ;

        void file_loop() {
            //mssgvf("        VCF format file:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "a");
        }
        //here we can put the VCF header data
        void file_header() {
            //want this here since the new line at filep_close() didn't seem to do the trick
            pr_nl();
        }

        void inner() {
            str a1, a2;
            if (_allele1 == 1)
                a1 = strdup("0");
            else if (_allele1 == 2)
                a1 = strdup("1");
            else
                a1 = strdup(".");

            if (_allele2 == 1)
                a2 = strdup("0");
            else if (_allele2 == 2)
                a2 = strdup("1");
            else
                a2 = strdup(".");

            pr_printf("%s/%s\t", a1, a2);
            //pr_printf("%s/%s\t", format_allele(_tlocusp, _allele1), format_allele(_tlocusp, _allele1));
            //pr_printf("%s/%s ", format_allele(_tlocusp, _allele1), format_allele(_tlocusp, _allele1));
            //pr_printf("%d/%d    ",_allele1,_allele2);

        }



        void loci_start() {
            pr_printf("%d\t", _numchr);
            pr_physical_distance(0);
            pr_printf("\t");
            pr_printf("%s\t",_tlocusp->LocusName);
            str a1, a2;
            if(!strcmp(_tlocusp->Allele[_allele1].AlleleName,"dummy"))
                a1 = strdup(".");
            else
                a1 = strdup(_tlocusp->Allele[_allele1].AlleleName);

            if(!strcmp(_tlocusp->Allele[_allele2].AlleleName,"dummy"))
                a2 = strdup(".");
            else
                a2 = strdup(_tlocusp->Allele[_allele2].AlleleName);
            pr_printf("%s\t%s\t",a1,a2);
            pr_printf(".\t");
            pr_printf("PASS\t");
            //info line: genetic distance in centimorgans, a (triple float avg, male, female), Allele frequency (float), genotype count (triple of ints?), and # of samples for data (int)
            pr_printf("CM=%.2f,%.2f,%.2f;AF=%.2f;GC=%s,%s,%s;NS=%d;\t",_tlocusp->Marker->pos_avg,_tlocusp->Marker->pos_male,_tlocusp->Marker->pos_female,_tlocusp->Allele[_allele2].Frequency,"count1","count2","count3",0);
            pr_printf("GT\t");


//            pr_printf("%-7d", _numchr);
//            pr_physical_distance(0);
//            pr_printf("%-1s","");
//            pr_printf("%-15s",_tlocusp->LocusName);
//            pr_printf("%-6s%-6s",_tlocusp->Allele[_allele1].AlleleName,_tlocusp->Allele[_allele2].AlleleName);
//            pr_printf("%-5s",".");
//            pr_printf("%-7s","PASS");
//            //info line, we want Allele frequency, genotype count, and # of samples
//            pr_printf("CM=%.2f,%.2f,%.2f;AF=%.2f;GC=%s,%s,%s;NS=%d%-6s ",_tlocusp->Marker->pos_avg,_tlocusp->Marker->pos_male,_tlocusp->Marker->pos_female,_tlocusp->Allele[_allele2].Frequency,"count1","count2","count3",_tlocusp->Marker->col_num,";");
//            pr_printf("%-7s","GT");
        }

        void loci_end(){
            pr_nl();
        }
    } *lp = new vcf_vcfs(Top);

    lp->setfln(prefix, ".vcf");
    lp->file_names = file_names;

    lp->load_formats_no_space(-1);
    //lp->load_formats(fwid, pwid, -1);

    lp->iterate();

    delete lp;
    delete hlp;
}

//this will write the pedigree in the PLINK fam file format
void CLASS_VCF::write_VCF_ped(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_peds,trait,ped_per) {
        vlpCTOR(vcf_peds,trait,ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF pedigree file:      %s/%s\n", *_opath, file_names[1]);
            data_loop(*_opath, file_names[1], "w");
        }
        void inner() {
            pr_fam();
            pr_per();
            pr_father();
            pr_mother();
            pr_sex();
            pr_pheno();
            pr_nl();

        }
    } *lp = new vcf_peds(Top);

    lp->setfln(prefix, ".fam");
    lp->file_names = file_names;

    lp->load_formats(fwid, pwid, -1);

    lp->iterate();

    delete lp;
}

// this will write the phenotypic data in
void CLASS_VCF::write_VCF_pheno(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_phenos,trait,ped_per) {
        vlpCTOR(vcf_phenos,trait,ped_per) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF phenotype file:     %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "w");
        }
        void file_header() {
            pr_printf("FID\tIID\tA1\tSAMPLEID\n");
//            pr_printf("FID");
//            pr_printf("IID");
//            pr_printf("A1");
//            pr_printf("SAMPLEID");
//            pr_nl();
        }
        void inner() {
            pr_fam();
            pr_printf("\t");
            pr_per();
            pr_printf("\t");
            pr_pheno();
            pr_printf("\t");
            pr_fam();
            pr_printf("_");
            pr_per();
            pr_nl();
        }
    } *lp = new vcf_phenos(Top);

    lp->setfln(prefix, ".phe");
    lp->file_names = file_names;

    //lp->load_formats(fwid, pwid, -1);
    lp->load_formats_no_space(-1);

    lp->iterate();

    delete lp;
}

//this will write the output of the VCF file
void CLASS_VCF::write_VCF_sh(linkage_ped_top *Top, const char *prefix, char *file_names[]) {

}

void CLASS_VCF::batch_in(){

}

void CLASS_VCF::batch_out(){

}

void CLASS_VCF::inner_file_names(char **file_names, const char *num, const char *stem) {

     sprintf(file_names[0], "%s", stem);
     sprintf(file_names[1], "%s", stem);
     sprintf(file_names[2], "%s", stem);
}

void CLASS_VCF::gen_file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num, file_name_stem);
}

void CLASS_VCF::replace_chr_number(char *file_names[], int numchr) {
    //change_output_chr(file_names[1], numchr);
}