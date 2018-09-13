/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

// for centos ...
#define __STDC_LIMIT_MACROS 1

#include "common.h"
#include "typedefs.h"
#include "types.hh"

#include "loop.h"
#include "sh_util.h"
#include "batch_input.h"
#include "utils_ext.h"

#include "fcmap_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "write_files_ext.h"
#include "read_files_ext.h"

#include "write_vcf_ext.h"
#include "zlib-1.2.8/zlib.h"

#include <ctime>
#include "dbrefallele.h"
#include "mega2_bcftools_interface.h"

extern "C" {
    #include "bcftools.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/vcf.h"
    //#include "lib/bcftools-1.6/filter.h"
    //#include "lib/bcftools-1.6/vcfview.h"
    //#include "lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h"
}

extern DBlite MasterDB;

extern int  db_exists_db();
extern void db_open_db();
extern void db_fini_all();
extern char DBfile[255];

Str hg_build;
int combinechromovcf;
int outfiletype;
Str ref_choice;

bcf_hdr_t *bcfheader[50];
int contigchrs[50];
float contigdist[50];
kstring_t *kstringbcf;
bcf1_t *bcfline;
htsFile *htsfileout;
int bcf_hdr_cnt;


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

    //printf("%d",outfiletype);
    LoopOverTrait = 0;

    omit_peds(untyped_ped_opt, Top);
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);

    inner_file_names(file_names,"",file_name_stem, &combine_chromo);

    printf("Mega2 created the following file(s) for VCF Format:\n");

    mssgvf("VCF file created using allele ordering setting: %s\n",ref_choice.c_str());
    if(_strand_flips)
        mssgvf("VCF output is aligned using the reference panel:\n   %s\n", mega2_input_files[REFfl]);

    write_VCF_file(Top, file_name_stem, file_names, pwid, fwid);
    write_VCF_ped(Top, file_name_stem, file_names, pwid, fwid);


    //we only want a phenotype file if we have more than one trait, the first trait is always put into the pedigree fam file by convention
    //if(num_traits>1)
    //I think we want this all the time since it contains sample id
    write_VCF_pheno(Top, file_name_stem, file_names, pwid, fwid);
    //write_VCF_sh(Top, file_name_stem, file_names);

    write_VCF_map(Top, file_name_stem, file_names, pwid, fwid);
    write_VCF_freq(Top, file_name_stem, file_names, pwid, fwid);
    write_VCF_pen(Top, file_name_stem, file_names, pwid, fwid);

//    if(outfiletype == 2) {
//        printf("\nMega2 is using VCF tools to convert to BCF format:\n");
//        convert_vcf_bcf(Top, file_name_stem, file_names, pwid, fwid);
//    }
//
//    if(outfiletype == 3) {
//        printf("Mega2 is using zlib to convert to VCF.gz format:\n");
//        convert_vcf_vcfgz(Top, file_name_stem, file_names, pwid, fwid);
//    }
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
        int dummychr;
        int dummylocus;
        double diff;
        int first;

        void file_loop() {
            mssgvf("        VCF format file:        %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "w");
        }

        void chr_start(){
            if (first != 1) {
                first = 0;
                bcf_hdr_cnt = 0;
            }
        }
        void loci_start() {
            if (first == 0)
                first = 1;
            else
                bcf_hdr_cnt++;

            contigchrs[bcf_hdr_cnt] = 0;
            contigdist[bcf_hdr_cnt] = 0;

            bcfheader[bcf_hdr_cnt] = bcf_hdr_init("w");

            bcf_hdr_append(bcfheader[bcf_hdr_cnt], "##fileformat=VCFv4.3\n");
            //pr_printf("##fileformat=VCFv4.1\n");
#ifdef HIDEDATE
            bcf_hdr_printf(bcfheader[bcf_hdr_cnt],"##filedate=%d%02d%02d\n", 1997, 8, 29);
            //pr_printf("##filedate=%d%02d%02d\n", 1997, 8, 29);
#else
            time_t now = time(0);
            tm *ltm = localtime(&now);
            bcf_hdr_printf(bcfheader[bcf_hdr_cnt], "##filedate=%d%02d%02d\n", 1900 + ltm->tm_year, 1 + ltm->tm_mon,
                           ltm->tm_mday);
            //pr_printf("##filedate=%d%02d%02d\n", 1900 + ltm->tm_year,1 + ltm->tm_mon,ltm->tm_mday );
#endif
            bcf_hdr_append(bcfheader[bcf_hdr_cnt], "##source=MEGA2\n");
            //pr_printf("##source=MEGA2\n");
            if (base_pair_position_index >= 0)
                bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                               "##INFO=<ID=CM,Number=3,Type=Float,Description=\"Genetic Distance in centimorgans (avg, male, female)\">\n");
            //pr_printf("##INFO=<ID=CM,Number=3,Type=Float,Description=\"Genetic Distance in centimorgans (avg, male, female)\">\n");
            bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                           "##INFO=<ID=RF,Number=1,Type=Float,Description=\"Allele Frequency of reference allele\">\n");
            //pr_printf("##INFO=<ID=RF,Number=1,Type=Float,Description=\"Allele Frequency of reference allele\">\n");
            bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                           "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency of alternate allele(s)\">\n");
            //pr_printf("##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency of alternate allele(s)\">\n");
            //add conditional
            if (_strand_flips) {
                bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                               "##INFO=<ID=NO,Number=0,Type=Flag,Description=\"No external reference allele panel match to this position. Major Allele was used instead.\">\n");
                //pr_printf("##INFO=<ID=NO,Number=0,Type=Flag,Description=\"No external reference allele panel match to this position. Major Allele was used instead.\">\n");
                bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                               "##INFO=<ID=AMBIG,Number=2,Type=String,Description=\"Reference panel has a match for this position (REF, ALT), but it is ambiguous.\">\n");
                //pr_printf("##INFO=<ID=AMBIG,Number=2,Type=String,Description=\"Reference panel has a match for this position (REF, ALT), but it is ambiguous.\">\n");
                bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                               "##INFO=<ID=ORIG,Number=2,Type=String,Description=\"REF and ALT values (REF,ALT) from original dataset flipped according to T/G <-> A/C.\">\n");
                //pr_printf("##INFO=<ID=ORIG,Number=2,Type=String,Description=\"REF and ALT values (REF,ALT) from original dataset flipped according to T/G <-> A/C.\">\n");
                bcf_hdr_append(bcfheader[bcf_hdr_cnt],
                               "##INFO=<ID=FLIP,Number=0,Type=Flag,Description=\"REF and ALT values (REF, ALT) from original dataset if REF and ALT alleles were switched in reference.\">\n");
                //pr_printf("##INFO=<ID=FLIP,Number=0,Type=Flag,Description=\"REF and ALT values (REF, ALT) from original dataset if REF and ALT alleles were switched in reference.\">\n");
            }
            //don't know these for now
            //pr_printf("##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
            //pr_printf("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
            bcf_hdr_append(bcfheader[bcf_hdr_cnt], "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            //pr_printf("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bcf_hdr_append(bcfheader[bcf_hdr_cnt], "##FILTER=<ID=PASS,Description=\"Passed variant FILTERs\">\n");
            //pr_printf("##FILTER=<ID=PASS,Description=\"Passed variant FILTERs\">\n");
        }

        //calculate length and print out number chromosomes
        //need different behavior for combine chromosome vs not (show all chromosome rows vs just one)
        //to get the length we loop over the ChrLoci array for the currently saved chromosome to determine the longest one, we need to make sure nothing is a trait so we check for markers
        //If the chromosome is different from the saved one we assume we need to print a new row
        void inner() {
            if (_tlocusp->Type == BINARY || _tlocusp->Type == NUMBERED) {
                if (contigchrs[bcf_hdr_cnt]  == 0)
                    contigchrs[bcf_hdr_cnt]  = _tlocusp->Marker->chromosome;
                if (contigdist[bcf_hdr_cnt] == 0)
                    contigdist[bcf_hdr_cnt] = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                if(_tlocusp->Marker->chromosome != contigchrs[bcf_hdr_cnt] ) {
                    if(combinechromovcf)
                        bcf_hdr_printf(bcfheader[bcf_hdr_cnt],"##contig=<ID=%d,length=%.0lf,assembly=%s>\n", contigchrs[bcf_hdr_cnt], contigdist[bcf_hdr_cnt], hg_build.c_str());
                    contigchrs[bcf_hdr_cnt]  = _tlocusp->Marker->chromosome;
                    contigdist[bcf_hdr_cnt] = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
                }
                if (_EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1 > diff)
                    contigdist[bcf_hdr_cnt] = _EXLTop->EXLocus[_locus].positions[base_pair_position_index] + 1;
            }
        }
        void loci_end(){
            bcf_hdr_printf(bcfheader[bcf_hdr_cnt],"##contig=<ID=%d,length=%.0lf,assembly=%s>\n", contigchrs[bcf_hdr_cnt], contigdist[bcf_hdr_cnt], hg_build.c_str());

            //make our own ped/per loop so that it plays nicely with bcftools header writer.
            //when it was a sepreate loop structure the contig wasn't saved properly and it wasn't writing the sample names
            //to each of the files, we didn't need the full functionality of the loop structure, just a list of samples
            for (int ped=0; ped < _Top->PedCnt; ped++) {
                linkage_ped_tree *tpedtreep = &(_Top->Ped[ped]);
                for (int per = 0; per < _Top->Ped[ped].EntryCnt; per++) {
                    linkage_ped_rec  *tpersonp = &(tpedtreep->Entry[per]);
                    char pedigree[16];
                    char person[16];

                    if (OrigIds[1] == 2) {
                        sprintf(pedigree, "%s", tpedtreep->Name);
                    } else if (OrigIds[1] == 3) {
                        sprintf(pedigree, "%d", ped + 1);
                    } else if (OrigIds[1] == 4) {
                        sprintf(pedigree, "%s", tpedtreep->Name);
                    } else if (OrigIds[1] == 6) {
                        sprintf(pedigree, "%s", tpedtreep->PedPre);
                    } else {
                        sprintf(pedigree, "%s", tpedtreep->PedPre);
                    }

                    if (OrigIds[0] == 1) {
                        sprintf(person, "%s", tpersonp->OrigID);
                    } else if (OrigIds[0] == 2) {
                        sprintf(person, "%s", tpersonp->OrigID);
                    } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                        sprintf(person, "%s", tpersonp->UniqueID);
                    } else if (OrigIds[0] == 6) {
                        sprintf(person, "%s", tpersonp->PerPre);
                    } else {
                        sprintf(person, "%s", tpersonp->PerPre);
                    }

                    char sample[33];
                    sprintf(sample, "%s_%s", tpedtreep->PedPre, tpersonp->PerPre);

                    bcf_hdr_add_sample(bcfheader[bcf_hdr_cnt], sample);
                }
            }
        }
    } *hlps = new vcf_vcfs_header_start(Top);

    hlps->file_names = file_names;

    hlps->load_formats_no_space(-1);

    hlps->iterate();

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
        bool first;
        char * dummycanon;
        char * canonA;
        char * canonC;
        char * canonT;
        char * canonG;

        HMapis references;
        HMapis alternates;
        HMapii strand_flips;
        HMapii major_minor_flips;
        HMapii dummies;
        HMapis oldrefs;
        HMapis oldalts;
        int lastchr;
        int extremum_allele;
        int bcf_first;

        void file_loop() {
            //mssgvf("        VCF format file:      %s/%s\n", *_opath, file_names[0]);
            data_loop(*_opath, file_names[0], "a");
        }

        void chr_start() {
            dummycanon = canonical_allele("dummy");
            canonA = canonical_allele("A");
            canonC = canonical_allele("C");
            canonG = canonical_allele("G");
            canonT = canonical_allele("T");

            if(bcf_first != 1) {
                bcf_hdr_cnt = 0;
                bcf_first = 0;
            }

            if (bcf_first == 0)
                bcf_first = 1;
            else
                bcf_hdr_cnt++;

            //this had to get moved since path_ isn't defined in chr_start()
            //file_names[0] doesn't include the directory so causes an error in output.
            //keeping this commented here for now in case there's an issue
//            if(outfiletype == 1)
//                htsfileout = hts_open(file_names[0], "w");
//            else if(outfiletype == 2)
//                htsfileout = hts_open(file_names[0], "wb");
//            else if(outfiletype == 3)
//                htsfileout = hts_open(file_names[0], "wg");
//
//            bcf_hdr_write(htsfileout,bcfheader[bcf_hdr_cnt]);


        }

        //here we can put the VCF header data
        void file_header() {
            if(_strand_flips)
                first = true;
            else
                first = false;

            lastchr = global_chromo_entries[0];

            if(outfiletype == 1)
                htsfileout = hts_open(path_, "w");
            else if(outfiletype == 2)
                htsfileout = hts_open(path_, "wb");
            else if(outfiletype == 3)
                htsfileout = hts_open(path_, "wg");

            bcf_hdr_write(htsfileout,bcfheader[bcf_hdr_cnt]);


        }

        void inner() {
            if(ref_choice == "Original_Order" ||  (_strand_flips && extremum_allele == 0)) {
                ksprintf(kstringbcf,"%s","\t");
                //pr_printf("\t");

                if (_allele1 == 0)
                    ksprintf(kstringbcf,"%s",".");
                    //pr_printf(".");
                else
                    ksprintf(kstringbcf,"%d",_allele1 - 1);
                    //pr_printf("%d", _allele1 - 1);

                ksprintf(kstringbcf,"%s","/");
                //pr_printf("/");

                if (_allele2 == 0)
                    ksprintf(kstringbcf,"%s",".");
                    //pr_printf(".");
                else
                    ksprintf(kstringbcf,"%d",_allele2 - 1);
                    //pr_printf("%d", _allele2 - 1);
            }
            else {
                int allele1;
                int allele2;


                if (_allele1 == 0)
                    allele1 = 0;
                else if (_allele1 - 1 == extremum_allele)
                    allele1 = 1;
                else if (_allele1 - 1 < extremum_allele)
                    allele1 = _allele1 + 1;
                else
                    allele1 = _allele1;

                if (_allele2 == 0)
                    allele2 = 0;
                else if (_allele2 - 1 == extremum_allele)
                    allele2 = 1;
                else if (_allele2 - 1 < extremum_allele)
                    allele2 = _allele2 + 1;
                else
                    allele2 = _allele2;

                ksprintf(kstringbcf,"%s","\t");
                //pr_printf("\t");


                if (allele1 == 0)
                    ksprintf(kstringbcf,"%s",".");
                    //pr_printf(".");
                else
                    ksprintf(kstringbcf,"%d", allele1 - 1);
                    //pr_printf("%d", allele1 - 1);

                ksprintf(kstringbcf,"%s","/");
                //pr_printf("/");

                if (allele2 == 0)
                    ksprintf(kstringbcf,"%s",".");
                    //pr_printf(".");
                else
                    ksprintf(kstringbcf,"%d", allele2 - 1);
                    //pr_printf("%d", allele2 - 1);
            }
        }

        void loci_start() {
            kstringbcf = new kstring_t();

            //check the chromosome, if it's changed we'll want to select the next set
            if(lastchr != _tlocusp->Marker->chromosome) {
                first = true;
                lastchr = _tlocusp->Marker->chromosome;
            }
            //here we grab our reference values and put them into a map
            if(first && _strand_flips){
                //clear the map if it has values (saves a bit of time)
                if(!references.empty())
                    references.clear();
                if(!alternates.empty())
                    alternates.clear();

                //db connection
                MasterDB.begin();
                DBstmt *select;
                char select_string[255];
                sprintf(select_string, "SELECT pos, ref, alt FROM ref_allele_table WHERE chr = %d;", _tlocusp->Marker->chromosome);
                select = MasterDB.prep(select_string);
                int ret = select && select->abort();

                //get position and reference and put them into our map
                while (ret) {
                    int position = 0;
                    char *reference;
                    char *alternate;
                    ret = select->step();
                    if (ret == SQLITE_ROW) {
                        select->column(0, position);
                        select->column(1, reference);
                        select->column(2, alternate);
                        references[position] = canonical_allele(reference);
                        alternates[position] = canonical_allele(alternate);
                    } else
                        break;
                }

                MasterDB.commit();
                delete select;

                DBstmt *select2;
                char select_string2[255];
                sprintf(select_string2, "SELECT marker, strand, major_minor, dummy, oldref, oldalt FROM ref_allele_flips");
                select2 = MasterDB.prep(select_string2);
                int ret2 = select2 && select2->abort();

                while (ret2) {
                    int locus = 0;
                    int strand_flip = 0;
                    int major_minor = 0;
                    int dummy = 0;
                    char *oldalt;
                    char *oldref;
                    ret2 = select2->step();
                    if (ret2 == SQLITE_ROW) {
                        select2->column(0, locus);
                        select2->column(1, strand_flip);
                        select2->column(2, major_minor);
                        select2->column(3, dummy);
                        select2->column(4, oldref);
                        select2->column(5, oldalt);
                        strand_flips[locus] = strand_flip;
                        major_minor_flips[locus] = major_minor;
                        dummies[locus] = dummy;
                        oldrefs[locus] = oldref;
                        oldalts[locus] = oldalt;
                    } else
                        break;
                }


                delete select2;
                first = false;
            }

            ksprintf(kstringbcf,"%d\t",_tlocusp->Marker->chromosome);
            //pr_printf("%d\t", _tlocusp->Marker->chromosome);
            //pr_physical_distance has an extra space
            //so this is basically the same code from in that function with the space removed
            double base_pair_position = 0.0;
            if (base_pair_position_index >= 0) {
                base_pair_position = _EXLTop->EXLocus[_locus].positions[base_pair_position_index];
            }
            ksprintf(kstringbcf,"%.0lf\t",base_pair_position);
            //pr_printf("%.0lf\t", base_pair_position);

            //pr_physical_distance(0);
            //pr_printf("\t");
            ksprintf(kstringbcf,"%s\t",_tlocusp->LocusName);
            //pr_printf("%s\t",_tlocusp->LocusName);

            //Here we want to loop over all Alleles, first value is the ref allele, all other comma separated Alt alleles.
            std::string a1;
            std::string a2;

            //extremum means minimum or maximum
            double extremum_frequency = 0;
            extremum_allele = 0;
            int reference_exists = 0;

            if(ref_choice == "Original_Order"){
                extremum_allele = 0;
            }

            if(ref_choice == "Major_Allele"){
                extremum_frequency = 0;
                for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                    if (_tlocusp->Allele[allele].Frequency > extremum_frequency) {
                        extremum_frequency = _tlocusp->Allele[allele].Frequency;
                        extremum_allele = allele;
                    }
                }
            }
            //change to 1.0 in case so we find the minimum
            if(ref_choice == "Minor_Allele"){
                extremum_frequency = 1.0;
                for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                    if (_tlocusp->Allele[allele].Frequency < extremum_frequency) {
                        extremum_frequency = _tlocusp->Allele[allele].Frequency;
                        extremum_allele = allele;
                    }
                }
            }

            std::string auxillary_ref;
            std::string auxillary_alt;
            if(_strand_flips) {
                //get references
                std::string ref_return;
                std::string alt_return;
                int lookup = (int) _EXLTop->EXLocus[_locus].positions[base_pair_position_index];
                //if we have the reference in the data set it's the extremum
                if (map_get(references, lookup, ref_return)) {
                    if (strcmp(ref_return.c_str(), "*") != 0) {
                        extremum_allele = 0;
                        reference_exists = 1;
                        map_get(alternates, lookup, alt_return);
                        auxillary_ref = ref_return;
                        auxillary_alt = alt_return;
                    }
                    else
                        reference_exists = 0;
                }
                else
                        reference_exists = 0;
                //otherwise calculate the extremum
                if(reference_exists != 1) {
                    extremum_frequency = 0;
                    for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                        if (_tlocusp->Allele[allele].Frequency > extremum_frequency) {
                            extremum_frequency = _tlocusp->Allele[allele].Frequency;
                            extremum_allele = allele;
                        }
                    }
                }
            }

            for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                if (allele == extremum_allele){
                    if(_tlocusp->Allele[allele].AlleleName == dummycanon)
                        a1 = ".";
                    else
                        a1 = _tlocusp->Allele[allele].AlleleName;
                    ref = a1;
                }
                else {
                    if (a2.empty()) {
                        if (_tlocusp->Allele[allele].AlleleName == dummycanon) {
                            alt.push_back(".");
                            a2 = a2 + ".";
                        }
                        else {
                            alt.push_back(_tlocusp->Allele[allele].AlleleName);
                            a2 = a2 + _tlocusp->Allele[allele].AlleleName;

                        }
                    }
                    else {
                        if (_tlocusp->Allele[allele].AlleleName == dummycanon) {
                            alt.push_back(".");
                            a2 = a2 + ",.";
                        } else {
                            alt.push_back(_tlocusp->Allele[allele].AlleleName);
                            a2 = a2 + ",";
                            a2 = a2 + _tlocusp->Allele[allele].AlleleName;
                        }
                    }
                }
            }

            ksprintf(kstringbcf,"%s\t%s\t.\tPASS\t",a1.c_str(),a2.c_str());
            //pr_printf("%s\t%s\t",a1.c_str(),a2.c_str());
            //pr_printf(".\t");
            //pr_printf("PASS\t");
            if(base_pair_position_index >= 0){
                ksprintf(kstringbcf,"%s","CM=");
                //pr_printf("CM=");
                if(_tlocusp->Marker->pos_avg != UNKNOWN_POSITION)
                    ksprintf(kstringbcf,"%.2f,",_tlocusp->Marker->pos_avg);
                    //pr_printf("%.2f,",_tlocusp->Marker->pos_avg);
                else
                    ksprintf(kstringbcf,"%s",".,");
                    //pr_printf(".,");
                if(_tlocusp->Marker->pos_male != UNKNOWN_POSITION)
                    ksprintf(kstringbcf,"%.2f,",_tlocusp->Marker->pos_male);
                    //pr_printf("%.2f,",_tlocusp->Marker->pos_male);
                else
                    ksprintf(kstringbcf,"%s",".,");
                    //pr_printf(".,");
                if(_tlocusp->Marker->pos_female != UNKNOWN_POSITION)
                    ksprintf(kstringbcf,"%.2f;",_tlocusp->Marker->pos_female);
                    //pr_printf("%.2f;",_tlocusp->Marker->pos_female);
                else
                    ksprintf(kstringbcf,"%s",".;");
                    //pr_printf(".;");

            }
            //pr_printf("CM=%.2f,%.2f,%.2f;",_tlocusp->Marker->pos_avg,_tlocusp->Marker->pos_male,_tlocusp->Marker->pos_female);
            //double alternate_frequency = 0;
            if(extremum_allele == -1)
                ksprintf(kstringbcf,"RF=%.6f;",0.0);
                //pr_printf("RF=%.6f;",0);
            else
                ksprintf(kstringbcf,"RF=%.6f;",_tlocusp->Allele[extremum_allele].Frequency);
                //pr_printf("RF=%.6f;",_tlocusp->Allele[extremum_allele].Frequency);
            ksprintf(kstringbcf,"%s","AF=");
            //pr_printf("AF=");
            for (int allele = 0; allele < _tlocusp->AlleleCnt; allele++) {
                if(allele == extremum_allele)
                    continue;
                else if((allele == 1 && extremum_allele == 0) || allele == 0)
                    ksprintf(kstringbcf,"%.6f",_tlocusp->Allele[allele].Frequency);
                    //pr_printf("%.6f",_tlocusp->Allele[allele].Frequency);
                else
                    ksprintf(kstringbcf,",%.6f",_tlocusp->Allele[allele].Frequency);
                    //pr_printf(",%.6f",_tlocusp->Allele[allele].Frequency);
            }

            ksprintf(kstringbcf,"%s",";");
            //pr_printf(";");


            if(_strand_flips) {
                std::string oldref = oldrefs[_tlocusp->locus_link];
                std::string oldalt = oldalts[_tlocusp->locus_link];
                if (reference_exists == 0 && strand_flips[_tlocusp->locus_link] == 0)
                    ksprintf(kstringbcf,"%s","NO;");
                    //pr_printf("NO;");
                else if(oldref == a2 && oldalt == a1)
                    ksprintf(kstringbcf,"%s","FLIP;");
                    //pr_printf("FLIP;");
                else if(strand_flips[_tlocusp->locus_link]) {
                    if(oldref == "dummy")
                        ksprintf(kstringbcf,"ORIG=%s,.;",oldref.c_str());
                        //pr_printf("ORIG=%s,.;",oldref.c_str());
                    else if(oldalt == "dummy")
                        ksprintf(kstringbcf,"ORIG=.,%s;",oldref.c_str());
                        //pr_printf("ORIG=.,%s;",oldalt.c_str() );
                    else
                        ksprintf(kstringbcf,"ORIG=%s,%s;",oldref.c_str(),oldalt.c_str());
                        //pr_printf("ORIG=%s,%s;",oldref.c_str(),oldalt.c_str());
                }
                else if(major_minor_flips[_tlocusp->locus_link] == 0 && (oldref != auxillary_ref || oldalt != auxillary_alt))
                    ksprintf(kstringbcf,"AMBIG=%s,%s;",auxillary_ref.c_str(),auxillary_alt.c_str());
                    //pr_printf("AMBIG=%s,%s;", auxillary_ref.c_str(),auxillary_alt.c_str());
            }
            //pr_printf("AF=%.6f;",alternate_frequency);
            //pr_printf("GC=%s,%s,%s;","count1","count2","count3");
            //pr_printf("NS=%d;",0);
            ksprintf(kstringbcf,"%s","\tGT");
            //pr_printf("\t");
            //pr_printf("GT");
        }

        void loci_end(){
            bcfline = bcf_init();
            vcf_parse(kstringbcf,bcfheader[bcf_hdr_cnt],bcfline);
            bcf_write(htsfileout,bcfheader[bcf_hdr_cnt],bcfline);
            ks_release(kstringbcf);
            bcf_empty(bcfline);
            //pr_nl();
        }

        void file_trailer() {
            vcf_close(htsfileout);
        }
    } *lp = new vcf_vcfs(Top);

    lp->file_names = file_names;

    lp->load_formats_no_space(-1);

    lp->iterate();

    delete lp;
    //delete hlp;
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
                    if (global_trait_entries[tr] < 0)
                        continue;

                    pr_printf("%s\t", _LTop->Pheno[global_trait_entries[tr]].TraitName);
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

        void file_loop() {
            mssgvf("        VCF phenotype file:     %s/%s\n", *_opath, file_names[2]);
            data_loop(*_opath, file_names[2], "a");
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
            pr_pheno();
            pr_printf("\t");
        }
    } *lp = new vcf_phenos(Top);

    lp->file_names = file_names;

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
                        pr_printf("Map.%c.a",_LTop->map_distance_type);
                    else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT || genetic_distance_sex_type_map == FEMALE_GDMT)
                        pr_printf("Map.%c.f\tMap.%c.m",_LTop->map_distance_type,_LTop->map_distance_type);
                }
            }
            if (base_pair_position_index >= 0)
                pr_printf("\t%s.p",_EXLTop->MapNames[base_pair_position_index]);
            pr_nl();
        }

        void inner() {
            pr_printf("%d\t", _tlocusp->Marker->chromosome);
            pr_printf("%s\t",_tlocusp->LocusName);
            if(genetic_distance_index >= 0){
                if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT)
                    pr_printf("%10.6f", _EXLTop->EXLocus[_locus].positions[genetic_distance_index]);
                else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT)
                    pr_printf("%10.6f\t%10.6f", _EXLTop->EXLocus[_locus].pos_female[genetic_distance_index], _EXLTop->EXLocus[_locus].pos_male[genetic_distance_index]);
            }
            if (base_pair_position_index >= 0)
                pr_printf("\t%.0lf", _EXLTop->EXLocus[_locus].positions[base_pair_position_index]);
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
        char * dummycanon;

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
                    pr_printf("%s\t%d\t%.6f\n", _LTop->Pheno[global_trait_entries[tr]].TraitName, al+1, _LTop->Locus[global_trait_entries[tr]].Allele[al].Frequency);
                }
            }
            dummycanon = canonical_allele(("dummy"));
        }
        void inner() {
            for(int i = 0; i < _tlocusp->AlleleCnt; i++) {
                if(_tlocusp->Allele[i].AlleleName != dummycanon)
                    pr_printf("%s\t%s\t%.6f\n", _tlocusp->LocusName, _tlocusp->Allele[i].AlleleName, _tlocusp->Allele[i].Frequency);
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


void CLASS_VCF::convert_vcf_bcf(char *filename){
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
    std::vector<std::string> args;
    args.push_back("bcftools");
    args.push_back("--output-type");
    args.push_back("b");
    args.push_back("--output-file");
    const char * extension = ".bcf";
    char * filename2 = strdup(filename);
    args.push_back(strcat(filename,extension));
    args.push_back(filename2);

    size_t argc = 6;
    char **argv;
    argv = (char **) malloc(argc * sizeof(char *));
    for (size_t ii = 0; ii < argc; ii += 1) {
        argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
        argv[ii] = &args[ii][0];
    }

    mbi->mega2_main_vcfview(argc,argv);

}


void CLASS_VCF::convert_vcf_vcfgz(char *filename) {
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
    std::vector<std::string> args;
    args.push_back("bcftools");
    args.push_back("--output-type");
    args.push_back("z");
    args.push_back("--output-file");
    const char * extension = ".gz";
    char * filename2 = strdup(filename);
    args.push_back(strcat(filename,extension));
    args.push_back(filename2);

    size_t argc = 6;
    char **argv;
    argv = (char **) malloc(argc * sizeof(char *));
    for (size_t ii = 0; ii < argc; ii += 1) {
        argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
        argv[ii] = &args[ii][0];
    }

    mbi->mega2_main_vcfview(argc,argv);
}

void CLASS_VCF::option_menu (char *file_names[], char *prefix, int *combine_chromo, linkage_ped_top *Top) {
    int choice, choice2, choice3, done, stem, build, chromo, fileout, ref,reftableexists,change_build_allowed, indmenu, pedmenu;
    reftableexists = 0;
    change_build_allowed = 1;

    std::string refchoice = "Major_Allele";
    strcpy(prefix, file_name_stem);
    char buildname[255] = "B37";
    choice = -1;

    //don't want any database activities unless we're in db mode
    if(database_read) {
        //check for reference table:
        //db_open_db();
        MasterDB.begin();
        DBstmt *check;
        char select_string[255] = "SELECT name FROM sqlite_master WHERE type='table' AND name='ref_allele_table';";

        check = MasterDB.prep(select_string);
        int retcheck = check && check->abort();
        while (retcheck) {
            retcheck = check->step();
            if (retcheck == SQLITE_ROW)
                reftableexists = 1;
            else
                reftableexists = 0;

            break;
        }

        MasterDB.commit();
        delete check;

        //if we find a table we want to check it has values
        //only do this if the reftable exists
        if (reftableexists) {
            DBstmt *select;
            char select_string2[255] = "SELECT COUNT(pos) FROM 'ref_allele_table'";

            int locuscnt = 0;
            select = MasterDB.prep(select_string2);
            int ret = select && select->abort();
            while (ret) {
                ret = select->step();
                if (ret == SQLITE_ROW) {
                    select->column(0, locuscnt);
                    //printf("%d,%d",locuscnt,Top->LocusTop->LocusCnt - 1);
                    if (locuscnt > 0)
                        reftableexists = 1;
                    else
                        reftableexists = 0;

                    break;
                }
            }

            delete select;
        }

        if (!reftableexists)
            refchoice = "Original_Order";
        else {
            refchoice = "Use Mega2 Allele DB Table";
            std::string rfile = mega2_input_files[REFfl];
            if (rfile.find("B37") != std::string::npos || rfile.find("b37") != std::string::npos)
                strcpy(buildname, "B37");
            if (rfile.find("HG37") != std::string::npos || rfile.find("hg37") != std::string::npos)
                strcpy(buildname, "HG37");
            if (rfile.find("B38") != std::string::npos || rfile.find("b38") != std::string::npos)
                strcpy(buildname, "B38");
            if (rfile.find("HG38") != std::string::npos || rfile.find("hg38") != std::string::npos)
                strcpy(buildname, "HG38");
            if (rfile.find("B19") != std::string::npos || rfile.find("b19") != std::string::npos)
                strcpy(buildname, "B19");
            if (rfile.find("HG19") != std::string::npos || rfile.find("hg19") != std::string::npos)
                strcpy(buildname, "HG19");
            change_build_allowed = 0;
        }
        if(!_strand_flips)
            refchoice = "Original_Order";

        if(_strand_flips)
            refchoice = "Use Mega2 Allele DB Table";
    }

    OrigIds[0] = 6;
    OrigIds[1] = 6;

    // actual menu loop
    while (choice != 0) {
        //change this up to be more dynamic, incrementing menu count and assigning the variables here.
        draw_line();
        int menu_count = 0;
        done = menu_count;
        printf("VCF Analysis Menu:\n");
        printf("%d) Done with this menu - please proceed\n", done);
        printf(" %d) File name stem:                                  %-15s\n", ++menu_count, prefix);
        stem = menu_count;

        printf(" %d) Human Genome Build                               %s\n", ++menu_count, buildname);
        build = menu_count;

        if(!_strand_flips)
            printf(" %d) Allele Ordering                                  %s\n", ++menu_count, refchoice.c_str());
        ref = menu_count;

        menu_count++;
        if(outfiletype == 1)
            printf(" %d) Choose Format:                                   VCF\n", menu_count);
        else if(outfiletype == 2)
            printf(" %d) Choose Format:                                   BCF\n", menu_count);
        else if(outfiletype == 3)
            printf(" %d) Choose Format:                                   VCF.gz\n", menu_count);
        fileout = menu_count;

        if(main_chromocnt > 1) {
            if (*combine_chromo)
                printf(" %d) Combine Chromosomes                              Yes\n", ++menu_count);
            else
                printf( " %d) Combine Chromosomes                              No\n", ++menu_count);

            chromo = menu_count;
        }

        else {
            //make sure this is set to remove a compiler warning.
            chromo = -1;
        }


        indmenu = ++menu_count;
        individual_id_item(indmenu, VCF, OrigIds[0], 43, 2, 0, 0);
        pedmenu = ++menu_count;
        pedigree_id_item(pedmenu, VCF, OrigIds[1], 43, 2, 0);


        printf("Enter selection: 0 - %d > ", menu_count);

        fcmap(stdin,"%d", &choice); newline;

        if ( choice < done ) {
            printf("Unknown option %d\n", choice);
        }

        else if ( choice == done ) {
            strcpy(file_name_stem,prefix);
            hg_build = buildname;
            ref_choice = refchoice;
            BatchValueSet(hg_build,"human_genome_build");
            BatchValueSet(outfiletype,"VCF_output_file_type");
            BatchValueSet(file_name_stem,"file_name_stem");
            BatchValueSet(ref_choice,"VCF_Allele_Order");
        }

        else if ( choice == stem ) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);
            newline;
        }

        else if ( choice == build ) {
            if(change_build_allowed) {
                printf("Enter human genome build > ");
                fcmap(stdin, "%s", buildname);
                newline;
            }
            else
                printf("Cannot change build when it's been read in from reference panel.\n");
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
            while(choice3 != 1 || choice3 != 2 || choice3 != 3 || choice3 != 4) {
                printf("Allele Ordering Menu\n");
                draw_line();
                printf("1) Use Original Allele Order from Input Data\n");
                printf("2) Use Major Allele Frequency\n");
                printf("3) Use Minor Allele Frequency\n");
                if(database_read) {
                    if (reftableexists && _strand_flips) {
                        printf("4) Use Mega2 Allele DB Table\n");
                        printf("Enter selection: 1 - 4 > ");
                    }
                    else
                        printf("Enter selection: 1 - 3 > ");
                }
                else
                    printf("Enter selection: 1 - 3 > ");


                fcmap(stdin, "%d", &choice3);
                newline;
                if (choice3 == 1) {
                    refchoice = "Original_Order";
                    break;
                }
                if (choice3 == 2) {
                    refchoice = "Major_Allele";
                    break;
                }
                else if (choice3 == 3){
                    refchoice = "Minor_Allele";
                    break;
                }
//                else if (choice3 == 4 && database_read) {
//                    if (reftableexists) {
//                        //if the ref table is already there we're all good
//                        printf("Using the existing external reference table.\n");
//                        refchoice = "Use Mega2 Allele DB Table";
//                        break;
//                    }
//                    else {
//                        printf("Reference Allele Table does not exist. You must create a new database with a reference panel to use this option.");
//                        //otherwise we want to load it up
//                        //keep this here in case someone wants to add to an existing database.
//                        //just now the way the file is handled is different
//                        //Reference_Allele_Table *reference_allele_table = new Reference_Allele_Table();
//                        //reference_allele_table->read_ref_allele_file(Top, reference_allele_table->get_filename(), false, 0);
//                        //reftableexists = 1;
//                        //refchoice = "Use Mega2 Allele DB Table";
//                        //break;
//                    }
//                }
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
        else if(choice == indmenu) {
            OrigIds[0] = individual_id_item(0, VCF, OrigIds[0], 35, 1,
                                       Top->OrigIds, Top->UniqueIds);
        }
        else if(choice == pedmenu) {
            OrigIds[1] = pedigree_id_item(0, VCF, OrigIds[1], 35, 1, Top->OrigIds);
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
    BatchValueGet(outfiletype,"VCF_output_file_type");
    BatchValueGet(ref_choice,"VCF_Allele_Order");
}

void CLASS_VCF::batch_out(){
    extern void batchf(batch_item_type *bi);

    Cstr Values[] =  { "file_name_stem",
                       "human_genome_build",
                       "Loop_Over_Chromosomes",
                       "VCF_output_file_type",
                       "VCF_Allele_Order"
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}

void CLASS_VCF::inner_file_names(char **file_names, const char *num, const char *stem, int *combine_chromo) {

    if(outfiletype == 1)
        sprintf(file_names[0], "%s.%s.vcf", stem, num);
    else if( outfiletype == 2)
        sprintf(file_names[0], "%s.%s.bcf", stem, num);
    else if (outfiletype == 3)
        sprintf(file_names[0], "%s.%s.vcf.gz", stem, num);
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
