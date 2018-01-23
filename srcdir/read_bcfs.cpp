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


/*
 *  read lines reseting stringstream/vector each line
 */

#include <stdlib.h>

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "common.h"
#include "tod.hh"
#include "error_messages_ext.h"

extern "C" {
    #include "bcftools.h"
    #include "lib/bcftools-1.6/filter.h"
    #include "lib/bcftools-1.6/vcfview.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h"
}

extern void           Exit(int arg, const char *file, const int line, const char *err);
#include "typedefs.h"
#include "fcmap_ext.h"
#include "batch_input.h"
#include "batch_input_ext.h"
#include "mrecode.h"
#include "mrecode_ext.h"
#include "reorder_loci_ext.h"
#include "plink_ext.h"
#include "annotated_ped_file_ext.h"
#include "annotated_ped_file.h"
#include "read_files_ext.h"
#include "utils_ext.h"

#include "str_utils.hh"
#include "input.hh"
#include "input_ops.hh"

#include "read_bcfs.h"
#include "mega2_bcftools_interface.h"

#ifdef _WIN
#define R_OK 4
#define access(str,type) _access(str,type)
#endif

using namespace std;

/*
 * Displays menu, adds additional options in menu 1 after input suffix
 * this gets called in user_input.cpp for whatever input opt is defined
 */
void ReadBCFs::do_menu_display(int &idx, int line_len, int choiceA[]) {
    printf("%2d) %-*s%s\n", idx, line_len, "BCF File Directory:",
           ((!strcmp(BatchItemGet("BCFs_File_Directory")->value.name, "."))?"[ Current directory ]" :
            BatchItemGet("BCFs_File_Directory")->value.name));
    choiceA[idx] = site_bcfs_dir_i;
    idx++;

    printf("%2d) %-*s%s\n", idx, line_len, "BCF File Template:", BatchItemGet("BCFs_File_Template")->value.name);
    choiceA[idx] = site_bcfs_template_i;
    idx++;
}
/*
 * Parses the menu displayed by ReadBCFs
 * We are getting a directory and a template for BCF files
 */
int ReadBCFs::do_menu_parse(int choice_) {
    int ret = 0;
    char bcfs_path_array[255], bcfs_template_array[255];
    char *bcfs_path = &bcfs_path_array[0];
    char *bcfs_template = &bcfs_template_array[0];
    if(choice_ == site_bcfs_dir_i) {
        while (1) {
            draw_line();
            printf("Please enter BCF directory name > ");
            fcmap(stdin, "%s", bcfs_path);
            newline;

            if (access(bcfs_path, F_OK)) {
                printf("WARNING: Could not find directory %s\n", bcfs_path);
                continue;
            } else if (!is_dir(bcfs_path)) {
                printf("WARNING: %s is not a directory.\n", bcfs_path);
                printf("Please specify a new or valid directory.\n");
                strcpy(bcfs_path, ".");
                continue;
            } else if (access(bcfs_path, W_OK)) {
                printf("WARNING: %s is not a writable directory.\n", bcfs_path);
                printf("Please specify a new or valid directory.\n");
                continue;
            }
            else {
                BatchValueSet(bcfs_path, "BCFs_File_Directory");
                break;
            }
        }
        ret = 1;
    } else if(choice_ == site_bcfs_template_i) {
        while (1) {
            draw_line();
            printf("To enter a template please enter a value of the form:\n");
            printf("[data?.bcf]\nWhere the wildecard '?' will replace the CHR number for all chromosomes.\n");
            printf("Please enter BCF file template format > ");
            fcmap(stdin, "%s", bcfs_template);
            newline;

            Vecs bcfsplit;
            split(bcfsplit, bcfs_template, "?");

            if (bcfsplit.size() != 2) {
                printf("Please include one and only one ? in the template name\n");
                continue;
            } else {
                BatchValueSet(bcfs_template, "BCFs_File_Template");
                break;
            }
        }
        ret = 1;
    }
    return ret;
}

/*
 * Interpret batch variables
 */
void ReadBCFs::do_menu2batch() {
    Cstr Values[] = {
            "BCFs_File_Directory",
            "BCFs_File_Template",
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read)
            batchf(bip);
    }
}
/*
 * Save batch values to internal variables
 */
void ReadBCFs::do_batch2local(){
    BatchValueGet(this->BCF_path, "BCFs_File_Directory");
    BatchValueGet(this->BCF_template, "BCFs_File_Template");
    show_settings();
}

/*
 * Prints the settings
 */
void ReadBCFs::show_settings() {
    msgvf("\n");
    msgvf("BCF File Directory:                         %s\n", C(this->BCF_path));
    msgvf("BCF File Template:                          %s\n", C(this->BCF_template));
}

/*
 * Called in Annotated_Ped_File.cpp
 * This is where we will call handlers to do various other tasks
 */
void ReadBCFs::do_init(Input_Base *inp)
{
    this->input = inp;
    this->pedfile = *inp->input_files.pedfl;
    this->phefile = *inp->input_files.phefl;

    check_bcf_files();
    this->num_samples = build_samples();
    //for(int i = 0; i< this->num_samples; i++)
    //    printf("%s\n",this->samples[i].c_str());
    build_markers();
    //checkindelsndups();

}

/*
 * here we will try to check that the files are available and read them
 */
void ReadBCFs::read_BCFs( linkage_locus_top *LPedTreeTop )
{
    //linkage_locus_top *LTop = LPedTreeTop;

    vector<string> temp;
    temp.push_back("bcftools");
    //temp.push_back("view");
    int argc = 1;
    int count = 0;

    Str directory = this->BCF_path;
    Str file_template = this->BCF_template;
    Vecs filesplit;

    split(filesplit,file_template,"?");

    for(int chr = 1; chr < 23; chr++){
        ifstream ifs;
        char file[255];
        if(chr >= 1  && chr < 10)
            sprintf(file, "./%s/%s0%d%s",directory.c_str(),filesplit[0].c_str(),chr,filesplit[1].c_str());
        else
            sprintf(file, "./%s/%s%d%s",directory.c_str(),filesplit[0].c_str(),chr,filesplit[1].c_str());
        ifs.open(file);
        if (! ifs.is_open() )
            errorvf("read_BCFs: Can not open \"%s\" file\n", file);
        else {
            mssgvf("read_BCFs: Found file \"%s\"\n", file);
            temp.push_back(file);
            count++;
        }
    }

    char** argv;
    argv = (char**)malloc(argc * sizeof(char*));
    for (size_t i = 0; i < argc; i += 1)
        argv[i] = (char*)malloc(255 * sizeof(char));

    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
    argv[0] = &temp[0][0];
    argc = 2;

    for( int i = 1; i < count +1; i++ ) {
        argv[1] = &temp[i][0];
        //printf("%d, %s %s\n",argc, argv[0],argv[1]);
        mbi->mega2_main_vcfview(argc, argv);
    }

    //MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
    //mbi->mega2_main_vcfview(argc, argv);

   //for (size_t i = 0; i < argc; i += 1)
    //    free(argv[i]);
    //free(argv);
}

void ReadBCFs::check_bcf_files() {

    vector<string> files;
    int count = 0;

    Str directory = this->BCF_path;
    Str file_template = this->BCF_template;
    Vecs filesplit;

    split(filesplit,file_template,"?");

    for(int chr = 1; chr < 23; chr++){
        ifstream ifs;
        char file[255];
        if(chr >= 1  && chr < 10)
            sprintf(file, "./%s/%s0%d%s",directory.c_str(),filesplit[0].c_str(),chr,filesplit[1].c_str());
        else
            sprintf(file, "./%s/%s%d%s",directory.c_str(),filesplit[0].c_str(),chr,filesplit[1].c_str());
        ifs.open(file);
        if (! ifs.is_open() )
            errorvf("read_BCFs: Can not open \"%s\" file\n", file);
        else {
            mssgvf("read_BCFs: Found file \"%s\"\n", file);
            files.push_back(file);
            count++;
        }
    }

    this->filelist = files;
    this->filecount = count;
    //for (auto i = fileslist.begin(); i != fileslist.end(); ++i)
    //    std::cout << *i << "\n";

}

int ReadBCFs::build_samples() {

    vector<string> files = this->filelist;
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

    int argc = 2;

    vector<string> temp;
    temp.push_back("bcftools");
    temp.push_back(&files[0][0]);

    char** argv;
    argv = (char**)malloc(argc * sizeof(char*));
    for (size_t i = 0; i < argc; i += 1) {
        argv[i] = (char*)malloc(255 * sizeof(char));
        argv[i] = &temp[i][0];
    }

    args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
    bcfargs = mbi->get_args(argc, argv);

    bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

    for(int i = 0; i < hdr->n[2]; i++){
        //printf("%s\n",hdr->samples[i]);
        this->samples.push_back(hdr->samples[i]);
    }

    return hdr->n[2];
}

void ReadBCFs::build_markers() {

    vector<string> files = this->filelist;
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

    int argc = 2;

    vector<string> args;
    args.push_back("bcftools");

    char** argv;
    argv = (char**)malloc(argc * sizeof(char*));
    for (size_t i = 0; i < argc; i += 1) {
        argv[i] = (char*)malloc(255 * sizeof(char));
        argv[i] = &args[i][0];
    }

    int total_markers = 0;
    for(int i = 0; i < this->filecount; i++) {
        args.push_back(files[i]);
        //for (auto i = temp.begin(); i != temp.end(); ++i)
        //    std::cout << *i << "\n";

        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t i = 0; i < argc; i += 1) {
            argv[i] = (char *) malloc(255 * sizeof(char));
            argv[i] = &args[i][0];
        }

        args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        int count = 0;

        while ( bcf_sr_next_line(bcfargs->files) ) {
            bcf1_t *line = bcfargs->files->readers[0].buffer[0];
            //the 1 references BCF_DT_CTG in vcf.h
            //this seems to be the actual way to get the chromosome
            //strcpy(markers[total_markers].chr ,hdr->id[1][line->rid].key);
            //printf("%s %s %d %s %s\n", line->d.id, hdr->id[1][line->rid].key, line->pos, line->d.allele[0],
            //       line->d.allele[1]);

            Vecs alleles(2);
            alleles.clear();
            char *al1 = canonical_allele(line->d.allele[0]);
            char *al2 = canonical_allele(line->d.allele[1]);
            //printf("%s %s",canonical_allele(line->d.allele[0]), canonical_allele(line->d.allele[0]));
            alleles.push_back(al1);
            alleles.push_back(al2);

            this->markers.push_back(new BCFMarker(line->d.id, hdr->id[BCF_DT_CTG][line->rid].key, line->pos, alleles));

            //printf("%s,%s\n",al1,al2);
            count++;
        }

        //printout to check that the marker array was populated.
        //BCFMarker *bp = NULL;
        //for (MarkerVectorP bpp = markers.cbegin(); bpp != markers.cend(); bpp++){
        //    bp     = *bpp;
        //    printf("%s %s %d %s %s\n", C(bp->name) , C(bp->chr), bp->pos, C(bp->alleles[0]), C(bp->alleles[1]));
        //}

        total_markers += count;
        printf("%d\n",total_markers);
        args.pop_back();
    }

    this->marker_count = total_markers;
}

linkage_ped_top *ReadBCFs::do_ped(linkage_locus_top *LTop)
{
//    Tod tod_pf("read ped file");
//    FILE *filep = this->pedfile ? fopen(this->pedfile, "r") : NULL;
//    if (filep == NULL) {
//        errorvf("could not open %s for reading!\n", pedfile);
//        EXIT(FILE_READ_ERROR);
//    }
//
//    linkage_ped_top *Top = read_common_ped_file(filep, this->pedfile, plink_info, LTop, file_desc,
//                                                phecols, num_groups, groups,
//                                                num_ped_records, has_extra_ids);
//    tod_pf();
//
//    return Top;

//    extern vector<Vecc> VecAlleles;
//    int num_peds = 0;
//
//    annotated_ped_rec *persons = build_impute2_ped(LTop, &num_peds);
//
//    build_genotypes(LTop, persons, VecAlleles);
//
//    linkage_ped_top *Top;
//    Top = mk_ped_top(persons, this->people.size(), LTop, num_peds,
//            /*untyped*/ 0, /*totaltyped*/ this->people_filtered,
//            /*groups*/ NULL, 0, 0,
//            /*num_err*/0, 1);
//
//    return Top;
//    return NULL;
    return NULL;

}



void ReadBCFs::do_map(std::vector<m2_map>& additional_maps)
{
    m2_map bcf_map;

    build_bcf_map(bcf_map);

    additional_maps.push_back(bcf_map);
}

void ReadBCFs::build_bcf_map(m2_map &bcf_map) {
    std::string alternative_key = std::string(Mega2BatchItems[/* 57 */ VCF_Marker_Alternative_INFO_Key].value.name);

    BCFMarker *bp = NULL;
    m2_map_entry map_entry;
    for (MarkerVectorP bpp = markers.cbegin(); bpp != markers.cend(); bpp++){
        bp     = *bpp;
        map_entry.set_chr(bp->chr);
        map_entry.set_POS(bp->pos);
        map_entry.set_REF(bp->alleles[0]);
        map_entry.set_marker_name(bp->name);
        bcf_map.push_back_entry(map_entry);

        int chr = STR_CHR(C(bp->chr));
        if (chr == SEX_CHROMOSOME) {
            human_x++;
        } else if (chr == MALE_CHROMOSOME) {
            human_y++;
        } else if (chr == PSEUDO_X) {
            human_xy++;
        } else if (chr == MITO_CHROMOSOME) {
            human_mt++;
        } else if (chr == UNKNOWN_CHROMO) {
            human_unknown++;
            NumUnmapped++;
        } else {
            human_auto++;
        }
    }
}


/*
 * We want to build "names" and genotypes in one pass to only load these BCF files once
 */

linkage_locus_top *ReadBCFs::do_names(const char *&names_fn)
{
    linkage_locus_top *LTop = build_BCFs_names();
    return LTop;
}

void ReadBCFs::do_phe_names(char *phe_file, char ***phe_names, int **phe_types, int phe_cols) {
    this->phefile = phe_file;
    this->phecols = phe_cols;
    this->phenames = phe_names;
    this->phetypes = phe_types;
}

linkage_locus_top *ReadBCFs::build_BCFs_names()
{
    int num_pheno   = this->phecols;
    int num_markers = this->marker_count;
    int num_all     = num_pheno + num_markers;

    char **names    = CALLOC(num_all, char *);
    char  *types    = CALLOC(num_all, char);

    int count;
    count = 0;

    for(int i = 0; i < num_pheno; i++){
        //names[i] = CALLOC(phenames[i].size()+1, char);
        //strcpy(names[i], this->phenames[i]);
        //types[i] = ( ((*typep) == "B") || ((*typep) == "D") ) ? 'A' : 'T';
        count++;
    }

    //int i;
    //i = 0;
//    Vecsp typep = sample_file_hdr2b.cbegin()+5;
//    for (Vecsp phep = sample_file_hdr1b.cbegin()+5; i < num_pheno; i++, phep++, typep++) {
//        names[i] = CALLOC((*phep).size()+1, char);
//        strcpy(names[i], (*phep).c_str());
//        types[i] = ( ((*typep) == "B") || ((*typep) == "D") ) ? 'A' : 'T';
//    }

    BCFMarker *bp;
    for (MarkerVectorP vecp = markers.cbegin(); vecp != markers.cend(); vecp++) {
        bp = *vecp;
        if (count >= num_all) {
            errorvf("Internal Error: build_BCFs_names() count of skipped markers too large\n");
            EXIT(DATA_INCONSISTENCY);
        }
        names[count] = CALLOC(bp->name.size()+1, char);
        strcpy(names[count], bp->name.c_str());
        types[count] = 'M';
        count++;
    }

    return read_common_marker_data(num_pheno + num_markers, num_markers, names, types,
            /*annotated*/ 1, /*penetrances_read*/ 0, 0, NULL);
}

void ReadBCFs::do_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons) {
    vector<string> files = this->filelist;
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

    int argc = 2;

    vector<string> args;
    args.push_back("bcftools");

    char** argv;
    argv = (char**)malloc(argc * sizeof(char*));
    for (size_t i = 0; i < argc; i += 1) {
        argv[i] = (char*)malloc(255 * sizeof(char));
        argv[i] = &args[i][0];
    }

    int total_markers = 0;
    for(int i = 0; i < this->filecount; i++) {
        args.push_back(files[i]);

        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t i = 0; i < argc; i += 1) {
            argv[i] = (char *) malloc(255 * sizeof(char));
            argv[i] = &args[i][0];
        }

        args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        int count = 0;

        while ( bcf_sr_next_line(bcfargs->files) ) {
            bcf1_t *line = bcfargs->files->readers[0].buffer[0];
            //bcf_unpack(line, BCF_UN_ALL);

            //vcf_format gets a kstring without having to use bcf_write
            kstring_t str = { 0, 0, NULL };
            vcf_format(hdr,line, &str);
            

            if(str.l > 0) {
                Vecs linesplit;
                const char *delimiter = "\t";
                char outline[str.m];

                strcpy(outline, str.s);
                split(linesplit, outline, delimiter);
                int person = 0;
                for(int s = 9; s < sizeof(linesplit); s++){
                    //printf("linesplit at %d: %s\n",s, linesplit[s].c_str());
                    Vecs personsplit;
                    const char *delimter2 = "/";

                    split(personsplit,linesplit[s].c_str(),delimter2);
                    if(strcmp(personsplit[0].c_str(),".") !=0) {

                        int index1 = stoi(personsplit[0]);
                        int index2 = stoi(personsplit[1]);
                        //printf("allele values of s: %s %s\n", personsplit[0].c_str(),
                        //       personsplit[1].c_str());
                        set_2Ralleles(persons[person].marker,count,&LTop->Locus[person],canonical_allele(line->d.allele[index1]),canonical_allele(line->d.allele[index2]));
                        //printf("%s %s-",canonical_allele(line->d.allele[index1]),canonical_allele(line->d.allele[index2]));

                    }
                    else{
                        set_2Ralleles(persons[person].marker,count,&LTop->Locus[person],canonical_allele("."),canonical_allele("."));
                        //printf("allele values of s: . .\n");
                    }
                    person++;
                }
            }
            free(str.s);
            count++;
        }
        args.pop_back();
    }
}
