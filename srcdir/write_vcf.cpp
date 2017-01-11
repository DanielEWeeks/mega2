/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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
#include "vcftools/mega2_vcftools_interface.h"
#include "vcftools/bcf_file.h"
#include "vcftools/vcf_file.h"
#include "vcftools/parameters.h"
#include "zlib-1.2.8/zlib.h"

#include "dblite.hh"
#include "dbmisc.hh"
#include "dbrefallele.h"

extern DBlite MasterDB;

extern int  db_exists_db();
extern void db_open_db();
extern void db_fini_all();
extern char DBfile[255];

Str hg_build;
int combinechromovcf;
int outfiletype;
Str ref_choice;
vector<string> references;


void CLASS_VCF::create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2) {
    int pwid, fwid, mwid;
    linkage_ped_top *Top = LPedTreeTop;

    int combine_chromo = main_chromocnt > 1;
    //1 VCF
    //2 BCF
    //3 VCF.gz
    outfiletype = 1;

    db_open_db();

    if ( InputMode == INTERACTIVE_INPUTMODE ) {
        option_menu(file_names,file_name_stem, &combine_chromo, Top);
        LoopOverChrm  = ! combine_chromo;
    }
    else {
        batch_in();
        combine_chromo = ! LoopOverChrm;
    }

    printf("%d",outfiletype);
    LoopOverTrait = 0;

    omit_peds(untyped_ped_opt, Top);
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    inner_file_names(file_names,"",file_name_stem, &combine_chromo);


    printf("Mega2 created the following file(s) for VCF Format:\n");
    write_VCF_file(Top, file_name_stem,file_names ,pwid, fwid);
    write_VCF_ped(Top, file_name_stem, file_names ,pwid, fwid);


    //we only want a phenotype file if we have more than one trait, the first trait is always put into the pedigree fam file by convention
    //if(num_traits>1)
    //I think we want this all the time since it contains sample id
    write_VCF_pheno(Top, file_name_stem, file_names ,pwid, fwid);
    //write_VCF_sh(Top, file_name_stem, file_names);

    write_VCF_map(Top, file_name_stem, file_names, pwid, fwid);
    write_VCF_freq(Top, file_name_stem, file_names, pwid, fwid);
    write_VCF_pen(Top, file_name_stem, file_names, pwid, fwid);

    if(outfiletype == 2) {
        printf("\nMega2 is using VCF tools to convert to BCF format:\n");
        convert_vcf_bcf(Top, file_name_stem, file_names, pwid, fwid);
    }

    if(outfiletype == 3) {
        printf("Mega2 is using zlib to convert to VCF.gz format:\n");
        convert_vcf_vcfgz(Top, file_name_stem, file_names, pwid, fwid);
    }
    printf("\n");
}


//CHROM POS ID REF ALT QUAL FILTER INFO FORMAT 1_1 1_2 2_1
//20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
void CLASS_VCF::write_VCF_file(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    //needed a loop to calculate the length for the contig flag

    vlpCLASS(vcf_vcfs_header_start,chr,loci) {
        vlpCTOR(vcf_vcfs_header_start, chr, loci) { }
        typedef char *str;
        str *file_names;
        bool first;
        int dummychr;
        int dummylocus;
        double diff;

        void file_loop() {
            mssgvf("        VCF format file:        %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }
        void file_header(){
            pr_printf("##fileformat=VCFv4.1\n");
            pr_printf("##filedate=%s\n",__TIMESTAMP__);
            pr_printf("##source=MEGA2\n");
            if(base_pair_position_index > 0)
                pr_printf("##INFO=<ID=CM,Number=3,Type=Float,Description=\"Genetic Distance in centimorgans (avg, male, female)\">\n");
            pr_printf("##INFO=<ID=RF,Number=1,Type=Float,Description=\"Allele Frequency of reference allele\">\n");
            pr_printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency of alternate allele(s)\">\n");
            pr_printf("##INFO=<ID=NO,Number=0,Type=Flag,Description=\"No external reference allele panel match to this position. Major Allele used instead.\">");
            //don't know these for now
            //pr_printf("##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
            //pr_printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
            pr_printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            pr_printf("##FILTER=<ID=PASS,Description=\"Passed variant FILTERs\">\n");
            first = true;
        }

        //calculate length and print out number chromosomes
        //need different behavior for combine chromosome vs not (show all chromosome rows vs just one)
        //to get the length we loop over the ChrLoci array for the currently saved chromosome to determine the longest one, we need to make sure nothing is a trait so we check for markers
        //If the chromosome is different from the saved one we assume we need to print a new row
        void inner(){
            if(combinechromovcf){
                if(first) {
                    if (_tlocusp->Type == BINARY || _tlocusp->Type == NUMBERED)
                        diff = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                    dummychr = _tlocusp->Marker->chromosome;
                    first = false;
                }
                else {
                    if(dummychr == _tlocusp->Marker->chromosome) {
                        if (_tlocusp->Type == BINARY || _tlocusp->Type == NUMBERED) {
                            if (_EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1 > diff) {
                                diff = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                            }
                        }
                    }
                    else{
                        pr_printf("##contig=<ID=%d,length=%.0lf,assembly=%s>\n", dummychr, diff, hg_build.c_str());
                        first = true;
                    }
                }
            }
            else {
                if (first) {
                    if (_tlocusp->Type == BINARY || _tlocusp->Type == NUMBERED)
                        diff = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                    dummychr = _tlocusp->Marker->chromosome;
                    first = false;
                }
                else {
                    if(dummychr == _tlocusp->Marker->chromosome) {
                        if (_tlocusp->Type == BINARY || _tlocusp->Type == NUMBERED) {
                            if (_EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1 > diff) {
                                diff = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                            }
                        }
                    }
                }
            }
        }
        void file_trailer(){
            pr_printf("##contig=<ID=%d,length=%.0lf,assembly=%s>\n", dummychr, diff, hg_build.c_str());
        }


    } *hlps = new vcf_vcfs_header_start(Top);

    hlps->file_names = file_names;

    hlps->load_formats_no_space(-1);

    hlps->iterate();

    //then we need a loop to make the header for the individuals 1_1 1_2 etc.
    vlpCLASS(vcf_vcfs_header,chr,ped_per) {
        vlpCTOR(vcf_vcfs_header, chr, ped_per) { }
        typedef char *str;
        str *file_names;
        bool first;

        void file_loop() {
            //mssgvf("        VCF format file:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "a");
        }

        void file_header() {
            pr_printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
        }


        void inner() {
            pr_fam();
            pr_printf("_");
            pr_per();
            pr_printf("\t");

        }

    } *hlp = new vcf_vcfs_header(Top);

    hlp->file_names = file_names;

    hlp->load_formats_no_space(-1);

    hlp->iterate();


    //finally a large loop for the data
    vlpCLASS(vcf_vcfs,chr,loci_ped_per) {
        vlpCTOR(vcf_vcfs,chr,loci_ped_per) { }
        typedef char *str;
        str *file_names;
        //need to save the ref for comparison
        std::string ref;
        //need to save a list of alts for comparison
        std::vector<std::string> alt;
        std::string altstring;
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
            if(_allele1 - 1 == -1)
                pr_printf(".");
            else
                pr_printf("%d",_allele1-1);

            pr_printf("/");

            if(_allele2 - 1 == -1)
                pr_printf(".");
            else
                pr_printf("%d",_allele2-1);

            pr_printf("\t");
        }



        void loci_start() {
            pr_printf("%d\t", _tlocusp->Marker->chromosome);
            pr_physical_distance(0);
            pr_printf("\t");
            pr_printf("%s\t",_tlocusp->LocusName);

            //Here we want to loop over all Alleles, first value is the ref allele, all other comma seperated Alt alelles.
            std::string a1;
            std::string a2;

            //extremum means minimum or maximum
            double extremum_frequency = 0;
            int extremum_allele = 0;
            int reference_exists = 0;

            if(ref_choice == "Major Allele"){
                extremum_frequency = 0;
                for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                    if (_tlocusp->Allele[allele].Frequency > extremum_frequency) {
                        extremum_frequency = _tlocusp->Allele[allele].Frequency;
                        extremum_allele = allele;
                    }
                }
            }
            //change to 1.0 in case so we find the minimum
            if(ref_choice == "Minor Allele"){
                extremum_frequency = 1.0;
                for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                    if (_tlocusp->Allele[allele].Frequency < extremum_frequency) {
                        extremum_frequency = _tlocusp->Allele[allele].Frequency;
                        extremum_allele = allele;
                    }
                }
            }

            if(ref_choice == "External Reference"){
                //db_open_db();
                MasterDB.begin();

                DBstmt *select;
                char select_string[100];
                if(base_pair_position_index >0)
                    sprintf(select_string,"SELECT chr, pos, marker, ref FROM ref_allele_table WHERE pos = %.0f AND chr = %d;",_EXLTop->EXLocus[_locus].positions[base_pair_position_index],_tlocusp->Marker->chromosome);
                select = MasterDB.prep( select_string);
                int ret = select && select->abort();
                while (ret){
                    int chromosome = 0;
                    int position = 0;
                    int marker = 0;
                    char *reference;

                    ret = select->step();
                    if (ret == SQLITE_ROW) {
                        select->column(0, chromosome);
                        select->column(1, position);
                        select->column(2, marker);
                        select->column(3, reference);

                        //printf("%d:%d:%d:%s\n",chromosome,position,marker,reference);
                        for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                            if(!strcmp(_tlocusp->Allele[allele].AlleleName,reference)){
                                extremum_allele = allele;
                                //printf("%d/%s\n",extremum_allele,_tlocusp->Allele[extremum_allele].AlleleName);
                                reference_exists = 1;
                            }
                        }
                    }
                    //check for empty return then select major allele
                    else{
                        extremum_frequency = 0;
                        for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                            if (_tlocusp->Allele[allele].Frequency > extremum_frequency) {
                                extremum_frequency = _tlocusp->Allele[allele].Frequency;
                                extremum_allele = allele;
                                reference_exists = 0;
                            }
                        }

                    }
                }

                MasterDB.commit();
            }


            for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                if (allele==extremum_allele){
                    if(!strcmp(_tlocusp->Allele[allele].AlleleName,"dummy"))
                        a1 = ".";
                    else
                        a1 = _tlocusp->Allele[allele].AlleleName;
                    ref = a1;
                }
                else {
                    if (a2.empty()) {
                        if (!strcmp(_tlocusp->Allele[allele].AlleleName, "dummy")) {
                            alt.push_back(".");
                            a2 = a2 + ".";
                        }
                        else {
                            alt.push_back(_tlocusp->Allele[allele].AlleleName);
                            a2 = a2 + _tlocusp->Allele[allele].AlleleName;
                        }
                    }
                    else{
                        if (!strcmp(_tlocusp->Allele[allele].AlleleName, "dummy")) {
                            alt.push_back(".");
                            a2 = a2 + ",.";
                        }
                        else {
                            alt.push_back(_tlocusp->Allele[allele].AlleleName);
                            a2 = a2 + ","+ _tlocusp->Allele[allele].AlleleName;
                        }}

                }

            }

            pr_printf("%s\t%s\t",a1.c_str(),a2.c_str());
            pr_printf(".\t");
            pr_printf("PASS\t");
            if(base_pair_position_index > 0)
               pr_printf("CM=%.2f,%.2f,%.2f;",_tlocusp->Marker->pos_avg,_tlocusp->Marker->pos_male,_tlocusp->Marker->pos_female);
            //double alternate_frequency = 0;
            pr_printf("RF=%.6f;",_tlocusp->Allele[extremum_allele].Frequency);
            pr_printf("AF=");
            for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                if(allele== extremum_allele)
                    continue;
                else
                    pr_printf("%.6f,",_tlocusp->Allele[allele].Frequency);
            }
            pr_printf(";");
            if(reference_exists)
                pr_printf("NO;");
            //pr_printf("AF=%.6f;",alternate_frequency);
            //pr_printf("GC=%s,%s,%s;","count1","count2","count3");
            //pr_printf("NS=%d;",0);
            pr_printf("\t");
            pr_printf("GT\t");
        }

        void loci_end(){
            pr_nl();
        }
    } *lp = new vcf_vcfs(Top);

    lp->file_names = file_names;

    lp->load_formats_no_space(-1);

    lp->iterate();

    delete lp;
    delete hlp;
    delete hlps;
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

    lp->file_names = file_names;

    lp->load_formats(fwid, pwid, -1);

    lp->iterate();

    delete lp;
}

// this will write the phenotypic data in
void CLASS_VCF::write_VCF_pheno(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    vlpCLASS(vcf_phenos_header,trait,ped_per_trait) {
        vlpCTOR(vcf_phenos_header, trait, ped_per_trait) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            data_loop(*_opath, file_names[2], "w");
        }

        void file_header() {
            int tr;
            pr_printf("FID\tIID\t");
                for (tr=0; tr < num_traits; tr++) {
                    if(global_trait_entries[tr] < 0)
                        continue;
                    //if(global_trait_entries[tr] != _trait)
                        pr_printf("%s\t",_LTop->Pheno[global_trait_entries[tr]].TraitName);
            }
            pr_printf("SAMPLEID\n");
        }
    } *header = new vcf_phenos_header(Top);

    header->file_names = file_names;
    header->load_formats_no_space(-1);
    header->iterate();

    vlpCLASS(vcf_phenos,trait,ped_per_trait) {
        vlpCTOR(vcf_phenos,trait,ped_per_trait) { }
        typedef char *str;
        str *file_names;
        linkage_locus_rec *dummy;

        void file_loop() {
            mssgvf("        VCF phenotype file:     %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "a");
        }
        void trait_start(){
            dummy = _ttraitp;
        }
        void per_start(){
            pr_fam();
            pr_printf("\t");
            pr_per();
            pr_printf("\t");
        }
        void per_end(){
            pr_fam();
            pr_printf("_");
            pr_per();
            pr_nl();
        }
        void inner() {
            //if(_ttraitp != dummy) {
                pr_pheno();
                pr_printf("\t");
            //}
        }
    } *lp = new vcf_phenos(Top);

    lp->file_names = file_names;

    //lp->load_formats(fwid, pwid, -1);
    lp->load_formats_no_space(-1);

    lp->iterate();

    delete lp;
    delete header;
}

void CLASS_VCF::write_VCF_map(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid) {
    vlpCLASS(VCF_map,chr,loci) {
        vlpCTOR(VCF_map,chr,loci) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF map file:           %s/%s\n", *_opath, file_names[3]);
            data_loop(*_opath, file_names[3], "w");
        }

        void file_header(){
            pr_printf("Chromosome\t");
            pr_printf("Name\t");
            if (genetic_distance_index >= 0) {
                if (_LTop->map_distance_type == 'k' || _LTop->map_distance_type == 'h') {
                    if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT)
                        pr_printf("Map.%c.a\t",_LTop->map_distance_type);
                    else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT || genetic_distance_sex_type_map == FEMALE_GDMT)
                        pr_printf("Map.%c.f\tMap.%c.m\t",_LTop->map_distance_type,_LTop->map_distance_type);
                }
            }
            if (base_pair_position_index >= 0)
                pr_printf("%s.p\t",_EXLTop->MapNames[base_pair_position_index]);
            pr_nl();
        }

        void inner() {
            pr_printf("%d\t", _tlocusp->Marker->chromosome);
            pr_printf("%s\t",_tlocusp->LocusName);
            if(genetic_distance_index >= 0){
                if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT)
                    pr_printf("%10.6f\t", _EXLTop->EXLocus[_locus].positions[genetic_distance_index]);
                else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT)
                    pr_printf("%10.6f\t%10.6f\t", _EXLTop->EXLocus[_locus].pos_female[genetic_distance_index], _EXLTop->EXLocus[_locus].pos_male[genetic_distance_index]);
            }
            if (base_pair_position_index >= 0)
                pr_printf("%.0lf\t", _EXLTop->EXLocus[_locus].positions[base_pair_position_index]);
            pr_nl();
        }

    } *xp = new VCF_map(Top);
    xp->file_names = file_names;
    xp->iterate();
    delete xp;
}

void CLASS_VCF::write_VCF_freq(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid) {
    vlpCLASS(VCF_freq,chr,loci) {
        vlpCTOR(VCF_freq,chr,loci) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF freq file:          %s/%s\n", *_opath, file_names[4]);
            data_loop(*_opath, file_names[4], "w");
        }

        void file_header(){
            pr_printf("Name\tAllele\tFrequency\n");
            for (int tr=0; tr < num_traits; tr++) {
                if(global_trait_entries[tr] < 0)
                    continue;
                for(int al = 0; al < _LTop->Locus[global_trait_entries[tr]].AlleleCnt; al++) {
                    pr_printf("%s\t%d\t%.4f\n", _LTop->Pheno[global_trait_entries[tr]].TraitName, al+1, _LTop->Locus[global_trait_entries[tr]].Allele[al].Frequency);
                }
            }
        }
        void inner() {
            for(int i = 0; i < _tlocusp->AlleleCnt; i++) {
                if(strcmp(_tlocusp->Allele[i].AlleleName,"dummy") != 0)
                    pr_printf("%s\t%s\t%.4f\n", _tlocusp->LocusName, _tlocusp->Allele[i].AlleleName, _tlocusp->Allele[i].Frequency);
            }
        }
    } *xp = new VCF_freq(Top);
    xp->file_names = file_names;
    xp->iterate();
    delete xp;
}

void CLASS_VCF::write_VCF_pen(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid) {
    vlpCLASS(VCF_pen,trait,null) {
        vlpCTOR(VCF_pen,trait,null) { }
        typedef char *str;
        str *file_names;

        void file_loop() {
            mssgvf("        VCF pen file:           %s/%s\n", *_opath, file_names[5]);
            data_loop(*_opath, file_names[5], "w");
        }

        void file_header(){
            pr_printf("Name\tClass\tPen.11\tPen.12\tPen.22\tType\n");
            for (int tr=0; tr < num_traits; tr++) {
                if (global_trait_entries[tr] < 0)
                    continue;
                if(_LTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
                    for(int cl=0; cl <_LTop->Pheno[global_trait_entries[tr]].Props.Affection.ClassCnt; cl++) {
                        if(_LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].AutoPen != NULL)
                            for(int ch =0; ch < main_chromocnt; ch++) {
                                if (global_chromo_entries[ch] != 23) {
                                    pr_printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%s\n",
                                              _LTop->Pheno[global_trait_entries[tr]].TraitName, cl + 1,
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].AutoPen[0],
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].AutoPen[1],
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].AutoPen[2],
                                              "autosomal");
                                    break;
                                }
                            }
                        for(int ch =0; ch < main_chromocnt; ch++) {
                            if (global_chromo_entries[ch] == 23) {
                                if (_LTop->Pheno[tr].Props.Affection.Class[cl].FemalePen != NULL)
                                    pr_printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%s\n",
                                              _LTop->Pheno[global_trait_entries[tr]].TraitName, cl + 1,
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].FemalePen[0],
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].FemalePen[1],
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].FemalePen[2],
                                              "female");
                                if(_LTop->Pheno[tr].Props.Affection.Class[cl].MalePen != NULL)
                                    pr_printf("%s\t%d\t%.4f\t0.0000\t%.4f\t%s\n",
                                              _LTop->Pheno[global_trait_entries[tr]].TraitName,
                                              cl+1, _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].MalePen[0],
                                              _LTop->Pheno[global_trait_entries[tr]].Props.Affection.Class[cl].MalePen[1],
                                              "male");
                            }
                        }
                    }
                }
            }

        }

        void trait_start() {
            //pr_printf("%s\t%d\t%.4f\t%.4f\t%.4f\t%s\n", _ttraitp->LocusName,1,0.00,0.00,0.00,"type");
        }


    } *xp = new VCF_pen(Top);
    xp->file_names = file_names;
    xp->iterate();
    delete xp;

}

//use VCFTools to turn our VCF output into a BCF file
void CLASS_VCF::convert_vcf_bcf(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid){
    char vcftools[10] = "vcftools";
    char vcfflag[10] = "--vcf";
    char recode[15] = "--recode-bcf";
    char outflag[10] = "--out";
    char outname[10] = "out";
    char *argv[] = {vcftools, vcfflag, file_names[0], recode,outflag,outname, NULL};
    int argc = sizeof(argv) / sizeof(char*) - 1;

    parameters params(argc,argv);

    params.read_parameters();

    params.vcf_filename=file_names[0];
    params.vcf_compressed = false;

    params.recode_all_INFO = true;
    params.recode_bcf = true;
    params.output_prefix = file_names[0];
    params.recode_bcf_to_stream = false;

    params.print_params();

    variant_file *vcf;
    vcf = new vcf_file(params.vcf_filename,params.vcf_compressed,params.chrs_to_keep,params.chrs_to_exclude,params.force_write_index);
    vcf->print_bcf(params.output_prefix,params.recode_INFO_to_keep,params.recode_all_INFO,params.recode_bcf_to_stream);
}


void CLASS_VCF::convert_vcf_vcfgz(linkage_ped_top *Top, const char *prefix, char **file_names, const int pwid, const int fwid) {
    FILE *infile = fopen(file_names[0], "rb");
    char outfilename[255];
    strcpy(outfilename,file_names[0]);
    strcat(outfilename,".gz");

    gzFile outfile = gzopen(outfilename, "wb");
    //if (!infile || !outfile) return -1;

    char inbuffer[128];
    int num_read = 0;
    unsigned long total_read = 0;

    while ((num_read = fread(inbuffer, 1, sizeof(inbuffer), infile)) > 0) {
        total_read += num_read;
        gzwrite(outfile, inbuffer, num_read);
    }

    fclose(infile);
    gzclose(outfile);

    printf("Read %ld bytes, Wrote %ld bytes,Compression factor %4.2f%%\n",total_read, file_size(outfilename), (1.0-file_size(outfilename)*1.0/total_read)*100.0);

}

unsigned long CLASS_VCF::file_size(char *filename)
{
    FILE *pFile = fopen(filename, "rb");
    fseek (pFile, 0, SEEK_END);
    unsigned long size = ftell(pFile);
    fclose (pFile);
    return size;
}


void CLASS_VCF::option_menu (char *file_names[], char *prefix, int *combine_chromo, linkage_ped_top *Top) {
    int choice, choice2, choice3, done, stem, build, chromo, fileout, ref,reftableexists;
    done = 0;
    stem = 1;
    build = 2;
    chromo = 5;
    fileout = 4;
    ref = 3;
    reftableexists = 0;

    //check for reference table:
    //db_open_db();
    MasterDB.begin();
    DBstmt *select;
    char select_string[255] = "SELECT name FROM sqlite_master WHERE type='table' AND name='ref_allele_table';";

    select = MasterDB.prep( select_string);
    int ret = select && select->abort();
    while (ret){
        ret = select->step();
        if (ret == SQLITE_ROW)
            reftableexists = 1;
        else
            break;
    }
    MasterDB.commit();

    //for testing remove afterwards.
    reftableexists = 0;

    //set more variables
    strcpy(prefix,file_name_stem);
    char buildname[255] = "B37";
    choice = -1;
    std::string refchoice;
    if(!reftableexists)
        refchoice = "Major Allele";
    else
        refchoice = "External Reference";


    // actual menu loop
    while (choice != 0) {
        draw_line();
        printf("VCF Analysis Menu:\n");
        printf("%d) Done with this menu - please proceed\n", done);
        printf("%d) File name stem:                                             %-15s\n", stem, prefix);
        printf("%d) Human Genome Build                                          %s\n", build, buildname);
        printf("%d) Reference Alleles                                           %s\n", ref, refchoice.c_str());
        if(outfiletype == 1)
            printf("%d) VCF/BCF/VCF.gz:                                             VCF\n", fileout);
        else if(outfiletype == 2)
            printf("%d) VCF/BCF/VCF.gz:                                             BCF\n", fileout);
        else if(outfiletype == 3)
            printf("%d) VCF/BCF/VCF.gz:                                             VCF.gz\n", fileout);
        if(main_chromocnt > 1) {
            if (*combine_chromo)
                printf("%d) Combine Chromosomes                                         Yes\n", chromo);
            else
                printf("%d) Combine Chromosomes                                         No\n", chromo);
            printf("Enter selection: 0 - %d > ", 5);
        }
        printf("Enter selection: 0 - %d > ", 4);
        fcmap(stdin,"%d", &choice); newline;

        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }

        else if ( choice == done ) {
            strcpy(file_name_stem,prefix);
            hg_build = buildname;
            BatchValueSet(hg_build,"human_genome_build");
            BatchValueSet(outfiletype,"VCF_output_file");

            ref_choice = refchoice;
        }

        else if ( choice == stem ) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);
            newline;
        }

        else if ( choice == build ) {
            printf("Enter human genome build > ");
            fcmap(stdin, "%s", buildname);
            newline;
        }

        else if(choice == chromo){
            if(*combine_chromo)
                *combine_chromo = 0;
            else
                *combine_chromo = 1;
            int tmp = (! combine_chromo) ? 'y' : 'n';
            BatchValueSet(tmp, "Loop_Over_Chromosomes");
        }

        else if(choice == ref){
            choice3 = -1;
            while(choice3 != 1 || choice3 != 2 || choice3 != 3) {
                printf("Reference Allele Menu\n");
                draw_line();
                printf("1) Use Major Allele Frequency\n");
                printf("2) Use Minor Allele Frequency\n");
                printf("3) Use external reference panel in database (see documentation)\n");
                printf("Enter selection: 1 - 3 > ");
                fcmap(stdin, "%d", &choice3);
                newline;
                if (choice3 == 1) {
                    refchoice = "Major Allele";
                    break;
                }
                else if (choice3 == 2){
                    refchoice = "Minor Allele";
                    break;
                }
                else if (choice3 == 3) {
                    if (reftableexists) {
                        //if the ref table is already there we're all good
                        refchoice = "External Reference";
                        break;
                    }
                    else {
                        //otherwise we want to load it up
                        //keep this here in case someone wants to add to an existing database.
                        //just now the way the file is handled is different
                        Reference_Allele_Table *reference_allele_table = new Reference_Allele_Table();
                        reference_allele_table->read_ref_allele_file(Top, reference_allele_table->get_filename());
                        reftableexists = 1;
                        break;
                    }
                }
                else
                    printf("Unknown option %d\n", choice3);
            }
        }

        else if(choice == fileout){
            choice2 = -1;
            while(choice2 != 1 || choice2 != 2 || choice2 != 3) {
                draw_line();
                printf("Choose output format:\n");
                printf("1) VCF output\n");
                printf("2) BCF output\n");
                printf("3) VCF.gz output\n");
                printf("Enter selection: 1 - 3 > ");
                fcmap(stdin, "%d", &choice2);
                newline;
                if (choice2 == 1 || choice2 == 2 || choice2 == 3) {
                    outfiletype = choice2;
                    break;
                }
                else
                    printf("Unknown option %d\n", choice2);
            }
        }

        else {
            printf("Unknown option %d\n", choice);
            newline;
        }

    }
}

//this will write the output of the VCF file
void CLASS_VCF::write_VCF_sh(linkage_ped_top *Top, const char *prefix, char *file_names[]) {
}


void CLASS_VCF::batch_in(){
    char c;
    char *fn = this->file_name_stem;
    BatchValueIfSet(fn,   "file_name_stem");
    BatchValueGet(hg_build, "human_genome_build");
    BatchValueGet(c,   "Loop_Over_Chromosomes");
    LoopOverChrm = c == 'y' || c == 'Y';
    BatchValueGet(outfiletype,"VCF_output_file");
}

void CLASS_VCF::batch_out(){
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "human_genome_build",
                       "Loop_Over_Chromosomes",
                       "VCF_output_file",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_VCF::inner_file_names(char **file_names, const char *num, const char *stem, int *combine_chromo) {

     sprintf(file_names[0], "%s.%s.vcf", stem, num);
     if(main_chromocnt == 1 || *combine_chromo)
        sprintf(file_names[1], "%s.%s.fam", stem, num);
    else
         sprintf(file_names[1], "%s.fam", stem);
    if(main_chromocnt == 1 || *combine_chromo)
        sprintf(file_names[2], "%s.%s.phe", stem, num);
    else
        sprintf(file_names[2], "%s.phe", stem);
     sprintf(file_names[3], "%s.%s.map", stem, num);
     sprintf(file_names[4], "%s.%s.freq", stem, num);
    if(main_chromocnt == 1 || *combine_chromo)
        sprintf(file_names[5], "%s.%s.pen", stem, num);
    else
        sprintf(file_names[5], "%s.pen", stem);

    combinechromovcf = *combine_chromo;
}


void CLASS_VCF::gen_file_names(char **file_names, char *num)
{
    //inner_file_names(file_names, num, file_name_stem);
}

void CLASS_VCF::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    if(main_chromocnt == 1 || combinechromovcf)
        change_output_chr(file_names[1], numchr);
    if(main_chromocnt == 1 || combinechromovcf)
        change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    if(main_chromocnt == 1 || combinechromovcf)
        change_output_chr(file_names[5], numchr);
}
