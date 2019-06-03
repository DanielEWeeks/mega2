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


/*
 *  read lines reseting stringstream/vector each line
 */

// for centos ...
#define __STDC_LIMIT_MACROS 1

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
#include "phe_lookup_ext.h"

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
#include "phe_lookup_ext.h"

#include "str_utils.hh"

#include "read_bcfs.h"
#include "mega2_bcftools_interface.h"

#ifdef _WIN
#define R_OK 4
#define access(str,type) _access(str,type)
#endif

#define cbegin() begin()
#define cend()   end()

SECTION_LOG_INIT(dup_marker_relabeled);

using namespace std;

/*
 * Displays menu, adds additional options in menu 1 after input suffix
 * this gets called in user_input.cpp for whatever input opt is defined
 */
void ReadBCFs::do_menu_display(int &idx, int line_len, int choiceA[]) {
    printf("%2d) %-*s%s\n", idx, line_len, "BCFtools Parameters", BatchItemGet("BCF_Args")->value.name);
    choiceA[idx] = site_bcfs_args_i;
    idx++;

    printf("%2d) %-*s%s%s\n", idx, line_len-11, "Data file or Manifest file:","[required] ", BatchItemGet("BCFs_File")->value.name);
    choiceA[idx] = site_bcfs_file_i;
    idx++;
}
/*
 * Parses the menu displayed by ReadBCFs
 * We are getting a directory and a template for BCF files
 */
int ReadBCFs::do_menu_parse(int choice_) {
    int ret = 0;
    char bcfs_file_array[FILENAME_LENGTH];
    char *bcfs_file = &bcfs_file_array[0];
    char bcf_args[FILENAME_LENGTH] = "";
    vector<string> files = this->filelist;
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

    if(choice_ == site_bcfs_file_i) {
        while (1){
            draw_line();
            if(Input->input_format == in_format_bcfs)
                printf("Enter a file in BCF, VCF.gz, or VCF format or enter a file that contains\na list of files in these formats:\n");
            fcmap(stdin, "%s", bcfs_file);
            newline;

            if (access(bcfs_file, F_OK)) {
                printf("WARNING: Could not find file %s\n", bcfs_file);
                continue;
            }
            else {
                BatchValueSet(bcfs_file, "BCFs_File");
                ret = 1;
                break;
            }
        }
    }

    else if(choice_ == site_bcfs_args_i) {
        draw_line();
        printf("\nCurrent BCF parameters:  %s\n", BatchItemGet("BCF_Args")->value.name);
        printf("   For more additional information on BCFTools flags\n");
        printf("see documentation at samtools.github.io/bcftools/bcftools.\n");
        printf("Mega2's BCFTools option is based off of BCFTools view\n");
        printf("and allows a subset of BCFTools view flags.\n");
        printf("\nValid BCFTools view options for Mega2 include:\n");
        printf("--known                     --novel \n");
        printf("--phased                    --exclude-phased \n");
        printf("--uncalled                  --exclude-uncalled\n");
        printf("--min-ac [INT]              --max-ac [INT] \n");
        printf("--min-alleles [INT]         --max-alleles[INT]\n");
        printf("--min-af [FLOAT]            --max-af [FLOAT]\n");
        printf("--exclude [EXPRESSION]      --include [EXPRESSION]\n");
        printf("--types [LIST]              --exclude-types [LIST]\n");
        printf("--regions [chr:to-from]     --regions-file [FILE]\n");
        printf("--targets [chr:to-from]     --targets-file [FILE]\n");
        printf("--samples [LIST]            --samples-file [FILE]\n");
        printf("--force-samples             --apply-filters [LIST]\n");
        printf("and their single letter counterparts.\n\n");
        printf("To remove all current BCFTools options enter \"clear\"\n");

        while (1) {
            printf("Please enter BCFTools arguments > \n");

            //this code was in user input for vcf arguments
            //it seems fcmap terminates on whitespace but we want a whole line
            fflush(stdout);
            IgnoreValue(fgets(bcf_args, sizeof(bcf_args)-1, stdin)); newline;
            int l = (int)strlen(bcf_args);
            if (bcf_args[l-1] == '\n') bcf_args[l-1] = 0;
            if (bcf_args[l-1] == '\r') bcf_args[l-1] = 0;

            // We want to test out the arguments we get
            //to do this we construct an argc and argv
            unsigned int argc = 2;
            vector<string> args;
            args.push_back("bcftools");

            Vecs argssplit;
            string extraargs = strdup(bcf_args);
            if(extraargs.find("clear") == 0) {
                mssgvf("Clearing BCFTools options.\n");
                this->BCF_args = "";
                BatchValueSet(this->BCF_args, "BCF_Args");
                break;
            }
            if(extraargs.find("-") != 0) {
                mssgvf("BCFTools options must begin with \"-\", \"--\", or be \"clear\".\n");
                continue;
            }

            if(!extraargs.empty()) {
                split(argssplit, extraargs, " ");
                for (unsigned int a = 0; a < argssplit.size(); a++) {
                    args.push_back(argssplit[a]);
                    argc++;
                }
            }
            //we need some sort of file to end our testargs
            args.push_back("dummy.bcf");

            char **argv;
            argv = (char **) malloc(argc * sizeof(char *));
            for (size_t ii = 0; ii < argc; ii += 1) {
                argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
                strcpy(argv[ii], &args[ii][0]);
            }

            int test = mbi->test_args(argc, argv);
            for (size_t ii = 0; ii < argc; ii += 1) {
                free(argv[ii]);
            }
            free(argv);

            args.clear();
            std::vector<std::string>().swap(args);

            if(test == 1) {
                this->BCF_args = strdup(bcf_args);
                BatchValueSet(this->BCF_args, "BCF_Args");
                ret = 1;
                break;
            }
            else {
                continue;
            }

        }
    }

    free(mbi);
    return ret;
}

/*
 * Interpret batch variables
 */
void ReadBCFs::do_menu2batch() {
    Cstr Values[] = {
            "BCF_Args",
            "BCFs_File",
            "BCF_Sample_Style",
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
    BatchValueGet(this->BCF_args, "BCF_Args");
    BatchValueGet(this->BCF_file, "BCFs_File");
    BatchValueGet(this->sample_ind, "BCF_Sample_Style");

    this->inpfile = *Input->input_files.bedfl;

    show_settings();
}

/*
 * Prints the settings
 */
void ReadBCFs::show_settings() {
    msgvf("\n");
    msgvf("Input File:                                 %s\n", ((this->inpfile) ? C(this->inpfile) : C(this->BCF_file)));
    msgvf("BCF Arguments:                              %s\n", C(this->BCF_args));
    msgvf("Sample Index:                               %d\n", this->sample_ind);
}


void ReadBCFs::init_filters() {
    int it = PLINK_Args;

    if (Mega2BatchItems[it].items_read)
        PLINK_args(Mega2BatchItems[it].value.name, 1);
}

/*
 * Called in Annotated_Ped_File.cpp
 * This is where we will call handlers to do various other tasks
 */
void ReadBCFs::do_init(Input_Base *inp)
{
    init_filters();

    this->input = inp;
    this->pedfile = *inp->input_files.pedfl;
    this->phefile = *inp->input_files.phefl;

    check_bcf_files();
    build_markers_and_samples();
    check_dups();
}


void ReadBCFs::check_bcf_files() {
    vector<string> files;
    int count = 0;
    ifstream ifs;

    char *ofile = *Input->input_files.bedfl;
    if (ofile) { // no ?
        ifs.open(ofile);
        if (!ifs.is_open())
            errorvf("read_BCFs: Can not open \"%s\" file\n", ofile);
        else {
            mssgvf("read_BCFs: Found file \"%s\"\n", ofile);
            files.push_back(ofile);
            count++;
        }
    }
        //single bcf/vcf/vcf.gz etc.
    else if ( this->BCF_file.substr(this->BCF_file.find_last_of(".") + 1) == "bcf"
              || this->BCF_file.substr(this->BCF_file.find_last_of(".") + 1) == "vcf"
              || this->BCF_file.substr(this->BCF_file.find_last_of(".") + 1) == "gz"){
      ifs.open(C(this->BCF_file));
        if (! ifs.is_open() )
                warnvf("read_BCFs: Can not open \"%s\" file\n", this->BCF_file.c_str());
            else {
                mssgvf("read_BCFs: Found file \"%s\"\n", this->BCF_file.c_str());
                files.push_back(this->BCF_file);
                count++;
                ifs.close();
            }
    }
        //manifest file split by lines with escape characters
    else {
      ifs.open(C(this->BCF_file));
        Str line;
        while (getline(ifs, line)){
            line = rtrim(line);
            line = ltrim(line);
            if (line.size() == 0 || "#" == line.substr(0, 1)) {
                mssgvf("read_BCFs: ignoring \"%s\"\n",line.c_str());
                continue;
            }
            else {
                ifstream innerifs;
                innerifs.open(C(line));
                if (! innerifs.is_open() )
                    warnvf("read_BCFs: Can not open \"%s\" file\n", line.c_str());
                else {
                    mssgvf("read_BCFs: Found file \"%s\"\n", line.c_str());
                    files.push_back(line);
                    count++;
                    innerifs.close();
                }
            }
        }
        ifs.close();
    }

    this->filelist = files;
    this->filecount = count;
}

void ReadBCFs::build_markers_and_samples() {
    vector<string> files = this->filelist;

    unsigned int argc = 3;
    vector<string> args;
    args.push_back("bcftools");

    Vecs argssplit;
    string extraargs = this->BCF_args;
    if(!extraargs.empty()) {
        split(argssplit, extraargs, " ");
        for (unsigned int a = 0; a < argssplit.size(); a++) {
            args.push_back(argssplit[a]);
            argc++;
        }
    }
    //only read header to get samples
    args.push_back("-h");


    for(int i = 0; i < this->filecount; i++) {
        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for file: %s\n", files[i].c_str());
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            strcpy(argv[ii], &args[ii][0]);
            if(ii == 1)
                printf("view ");
            printf("%s ", argv[ii]);
        }
        printf("\n");

        args_t *bcfargs = (args_t *) calloc(1, sizeof(args_t));
        MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
        bcfargs = mbi->get_args(argc, argv);

        for (size_t ii = 0; ii < argc; ii += 1) {
            free(argv[ii]);
        }
        free(argv);
        free(mbi);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);
        //consistancy between the number of samples (only way to tell if there's less in subsequent files
        if(i == 0)
            this->num_samples = hdr->n[2];
        else if( this->num_samples != hdr->n[2])
            errorvf("Different number of sample columns between first file in manifest and file %s\n",files[i].c_str());

        int err = 0;
        for (int j = 0; j < this->num_samples; j++) {
            Str name = hdr->samples[j];
            //first file just add
            if (i == 0) {
                this->samples.push_back(name);
                Pairsi pair(name, j);
                std::pair<HMapsip, bool> ans = hdrMap.insert(pair);
                if (ans.second == false) {
                    errorvf("name %s is repeated.  It was sample %d and now reappears as sample %d\n",
                            C(name), ans.first->second, j);
                    err++;
                }
            }
            //other files use hdrMap to check the names are the same
            else {
                if (hdrMap.find(name) == hdrMap.end())
                    errorvf("Sample in file %s not found in first file from manifest\n", files[i].c_str());
            }
        }
        if (err)
            EXIT(DATA_INCONSISTENCY);

        args.pop_back();

        destroy_data_vcfview(bcfargs);
        free(bcfargs);

        std::cout << "Samples for " << files[i] << " completed\n";// << std::ctime(&time)
    }

    args.pop_back();
    //no header since we have the samples
    args.push_back("-H");
    //drop genotypes to not read them while getting allele labels
    args.push_back("-G");
    argc++;



    int total_markers = 0;
    for(int i = 0; i < this->filecount; i++) {

        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for file: %s\n",files[i].c_str());
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            strcpy(argv[ii], &args[ii][0]);
            if(ii == 1)
                printf("view ");
            printf("%s ",argv[ii]);
        }
        printf("\n");

        args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
        MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();
        bcfargs = mbi->get_args(argc, argv);

        for (size_t ii = 0; ii < argc; ii += 1) {
            free(argv[ii]);
        }
        free(argv);
        free(mbi);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        int count = 0;
        char sname[50];
        if (! hdr->id[BCF_DT_CTG]) {
            errorf("Probably missing ##contig entries.\n"
                   "Possibly problems with marker position.\n"
                   "ABORTING!\n");
            EXIT(DATA_INCONSISTENCY);
        }
        while (bcf_sr_next_line(bcfargs->files) ) {
            bcf1_t * line = bcfargs->files->readers[0].buffer[0];
            if ( subset_vcf(bcfargs, line) ) {
                Vecs alleles;

                for (int al = 0; al < line->n_allele; al++) {
                    alleles.push_back(canonical_allele(line->d.allele[al]));
                }

                Str name = line->d.id;
                const char *posp = hdr->id[BCF_DT_CTG][line->rid].key;
                if (strncmp(posp, "chr", 3) == 0) posp += 3;
//rvb: and another
                if (name == ".") {
                    sprintf(sname, "chr%s_%d", posp, line->pos + 1);
                    name = sname;
                }

                BCFMarker * line_marker = new BCFMarker(name, posp, line->pos + 1, alleles);// pos + 1 matches the VCF line pos field.
                this->markers.push_back(line_marker);

                alleles.clear();
                std::vector<std::string>().swap(alleles);
                count++;
            }
            if(line) bcf_clear(line);
        }
        int err = bcfargs->files->errnum;
        if ( err ) {
            fprintf(stderr,"Error: %s\n", bcf_sr_strerror(bcfargs->files->errnum));
            errorf("BCFTools encountered a problem with the given command.\n"
                   "Mega2 was unable to process allele labels and is exiting.\n"
                   "Please verify your BCFTools command is valid.\n");
            EXIT(DATA_INCONSISTENCY);
        }


        total_markers += count;
        args.pop_back();

        std::cout << "Alleles for " << files[i] << " completed\n";

        destroy_data_vcfview(bcfargs);
        free(bcfargs);

    }
    args.clear();
    std::vector<std::string>().swap(args);

    this->marker_count = total_markers;
}



void ReadBCFs::check_dups() {
    Hmapsi markerMap;

    SECTION_LOG_INIT(dup_marker_relabeled);
    for(int i = 0; i < this->marker_count; i++){
        Str name = this->markers[i]->name;
        if (markerMap.find(name) == markerMap.end())
            markerMap[name] = 1;
        else {
            int value = markerMap[name] + 1;
            char newname[FILENAME_LENGTH] ;
            sprintf(newname, "%s_%d", name.c_str(), value);
            markerMap[newname] = 1;
            markerMap[name] = value;
            this->markers[i]->name = newname;
            SECTION_LOG(dup_marker_relabeled);
            mssgvf("Duplicate marker name found for %s relabeling %s:%d as %s\n", name.c_str(),
                   this->markers[i]->chr.c_str(), this->markers[i]->pos, newname);
        }
    }

    SECTION_LOG_FINI(dup_marker_relabeled);

    markerMap.clear();
    Hmapsi().swap(markerMap);
}

void ReadBCFs::do_map(std::vector<m2_map>& additional_maps)
{
    m2_map bcf_map("VCF (created from input)", 'p');

    extern int just_gen_batch_file;

    if (! just_gen_batch_file)
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
    //we need to call this as a getopt function
    linkage_locus_top *LTop = build_BCFs_names();

    //deprecated
    //phenotype_file_SAMPLEID_entry_checks();

    return LTop;
}

void ReadBCFs::do_phe_names(char *phe_file, char **phe_names, int *phe_types, int phe_cols) {
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

    int count = 0;

    for(int i = 0; i < num_pheno; i++){
        names[i] = strdup(phenames[i]);
        types[i] = ( phetypes[i] == 0 ) ? 'A' : 'T';
        count++;
    }

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

//static void phenotype_file_SAMPLEID_entry_checks()
//{
    //this was removed as it didn't work the way I anticipated
    //sampleID checking against the .phe file
    // and PED_PER /PER now happens in do_genotypes.
//}

void ReadBCFs::do_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons,
                            std::vector<Vecc> &VecAlleles, int num_ped_recs) {
    vector <string> files = this->filelist;
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

    unsigned int argc = 2;

    vector <string> args;
    args.push_back("bcftools");

    Vecs argssplit;
    string extraargs = this->BCF_args;
    if(!extraargs.empty()) {
        split(argssplit, extraargs, " ");
        for (unsigned int a = 0; a < argssplit.size(); a++) {
            args.push_back(argssplit[a]);
            argc++;
        }
    }

    int mrkindex = LTop->PhenoCnt;
    int zeroed_genotypes_nonx = 0;
    int zeroed_genotypes_x = 0;

    hdrInd = new int[num_ped_recs];

    if(InputMode == INTERACTIVE_INPUTMODE) {
        extern int _bcf_sample_index;
        sample_ind = _bcf_sample_index;
    }


    SECTION_LOG_INIT(sampleid_mismatch);

    if(sample_ind != 0) {
        mssgvf("\nChecking SAMPLEID consistency within files.\n");
        if(sample_ind == 1)
            mssgvf("Since 'BCF_Sample_Style' is 1, the VCF header is checked against the .phe file SampleID Column.\n");
        else if(sample_ind == 2)
            mssgvf("Since 'BCF_Sample_Style' is 2, the VCF header is checked against '<pedID>_<perID>'.\n");
        else if(sample_ind == 3)
            mssgvf("Since 'BCF_Sample_Style' is 3, the VCF header is checked against '<perID>'.\n");
        else {
            errorvf(" Since 'BCF_Sample_Style' is not defined properly as 1, 2 or 3, we cannot proceed.\n"
                    " Please rerun in interactive mode to create a new batch file with this value set,\n"
                    " or add the value to the mega2 batch file:\n"
                    " BCF_Sample_Style=1 the VCF header is checked against the .phe file SampleID Column\n"
                    " BCF_Sample_Style=2 the VCF header is checked against '<pedID>_<perID>\n"
                    " BCF_Sample_Style=3 the VCF header is checked against '<<perID>\n");
            EXIT(DATA_INCONSISTENCY);
        }
        draw_line();

        for(int per = 0; per < num_ped_recs; per++) {
            bool done = false;
            char ped_per[FILENAME_LENGTH];
            sprintf(ped_per,"%s_%s",persons[per].PedID,persons[per].PerID);
            //is SAMPLEIDS map null? no phe file column
            if (sample_ind == 1) {
                if(SAMPLEIDS != NULL) {
                    Pairss pp(string(persons[per].PedID),string(persons[per].PerID));
                    Str sample;
                    int idx;
#if 1
                    if (map_get(*PEDPERIDS, pp, sample)) {
                        if (map_get(hdrMap, sample, idx)) {
                            hdrInd[per] = idx;
                            done = true;
                        } else {
                            SECTION_LOG(sampleid_mismatch);
                            warnvf("sample %s [%s (#%d)] Not found in BCF file\n",
                                   C(sample), ped_per, per);
                        }
                    } else {
                        SECTION_LOG(sampleid_mismatch);
                        warnvf("%s (#%d) Not found in SAMPLEID map\n",
                               ped_per, per);
                    }
#else
            /*Initial way I wrote this code which was a tremendously inefficient use of the map structures */

                //check SAMPLEIDS map from phenotype file
//                    for(sample_map_type::iterator sampleit = SAMPLEIDS->begin(); sampleit != SAMPLEIDS->end(); sampleit++) {
//                        for(HMapsip hdrit = hdrMap.begin(); hdrit != hdrMap.end(); hdrit++){
//                            //check if the sample is in both the SAMPLEIDS and the headers from the BCF FILE
//                            if(sampleit->first == hdrit->first) {
//                                //check to see that the values we got from the SAMPLIDS map match the persons table
//                               if(sampleit->second.first == persons[per].PedID && sampleit->second.second ==persons[per].PerID) {
//                                    //printf("%s %s %d \n", sampleit->first.c_str(), hdrit->first.c_str(), hdrMap.find(sampleit->first)->second);
//                                   hdrInd[per] = hdrMap.find(sampleit->first)->second;
//                                    done = true;
//                                    break;
//                                }
//                           }
//                        }
//                    }
#endif
                }
            } else if(sample_ind == 2) {
                //check BCF FILE header name against PED_PER
                if(!(hdrMap.find(ped_per) == hdrMap.end()) && !done){
                    hdrInd[per] = hdrMap.find(ped_per)->second;
                    done = true;
                }
            } else if(sample_ind == 3) {
                //check against PER
                if(!(hdrMap.find(persons[per].PerID)== hdrMap.end()) && !done){
                    hdrInd[per] = hdrMap.find(persons[per].PerID)->second;
                    done = true;
                }
            }
             if(done){
                //need a catch for the sample map if it's found so we don't put an warning
                //but we don't need to do anything
            }
            //finally warning and zero out the value, we either got this value in the BCF FILE with  no match
            //or we got it in the pedigree and couldn't find the match in the BCF FILE
            else {
                SECTION_LOG(sampleid_mismatch);
                mssgvf("Cannot find match in VCF file within provided pedigree for pedigree: %s person: %s\n",persons[per].PedID, persons[per].PerID);
                hdrInd[per] = -1;
            }
        }
    } else {
        // map 1:1 (for no map ped)
        for(int per = 0; per < num_ped_recs; per++) {
            hdrInd[per] = per;
        }
    }

    SECTION_LOG_FINI(sampleid_mismatch);


    for (int i = 0; i < this->filecount; i++) {
        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for file: %s\n",files[i].c_str());
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            strcpy(argv[ii], &args[ii][0]);
            if(ii == 1)
                printf("view ");
            printf("%s ",argv[ii]);
        }
        printf("\n");

        args_t *bcfargs = (args_t *) calloc(1, sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        for (size_t ii = 0; ii < argc; ii += 1) {
            free(argv[ii]);
        }
        free(argv);
        free(mbi);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        while (bcf_sr_next_line(bcfargs->files)) {
            bcf1_t *line = bcfargs->files->readers[0].buffer[0];
                if ( subset_vcf(bcfargs, line) && mrkindex < LTop->LocusCnt) {
                int m, n, p;

                //uses convert object so we need our own void * instead
                void *dat = NULL;

                m = 0;
                n = bcf_get_genotypes(hdr, line, &dat, &m);

                if (n <= 0) {
                    error("Error parsing GT tag at %s:%d\n", bcf_seqname(hdr, line), line->pos + 1);
                }

                Vecc canons;
                extern void set_2Ralleles_2bits(int marker, linkage_locus_rec *locus, const char *all1, const char *all2);

                for (int al = 0; al < line->n_allele; al++) {
                    canons.push_back(canonical_allele(line->d.allele[al]));
                }
                VecAlleles.push_back(canons);
                set_2Ralleles_2bits(mrkindex, &LTop->Locus[mrkindex], canons[0], canons[1]);

                const char * chromosome = hdr->id[BCF_DT_CTG][line->rid].key;

                //should give number of allele options per marker
                n /= num_samples;
                //this used to be p, but now we want p to traverse all ped records from external ped files
                //int vcf_col = 0;
                for (p = 0; p < num_ped_recs; p++) {
                    if(hdrInd[p] == -1) {
                        set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                        continue;
                    }
                    int32_t *ptr = (int32_t *) dat + hdrInd[p] * n;
                    int j;
                    for (j = 0; j < n; j++)
                        if (ptr[j] == bcf_int32_vector_end) break;


                    int is_x_chr = 0;
                    if(((strcasecmp(chromosome,"X") == 0) ||
                        (strcasecmp(chromosome,"CHRX") == 0) ||
                        (strcmp(chromosome,"23") == 0)) && (persons[p].Sex == 1))
                        is_x_chr = 1;

                    // diploid
                    if (j == 2) {
                        // ./.
                        if (bcf_gt_is_missing(ptr[0]) && bcf_gt_is_missing(ptr[1]))
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                        // ./ALT
                        else if (bcf_gt_is_missing(ptr[0]) && bcf_gt_allele(ptr[1]) == 1) {
                            if (MARKER_SCHEME > 1)
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, canons[1]);
                            else {
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                                if(is_x_chr)
                                    zeroed_genotypes_x++;
                                else
                                    zeroed_genotypes_nonx++;
                            }
                        }
                        // ./REF
                        else if (bcf_gt_is_missing(ptr[0])) {
                            if (MARKER_SCHEME > 1)
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, canons[0]);
                            else {
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                                if(is_x_chr)
                                    zeroed_genotypes_x++;
                                else
                                    zeroed_genotypes_nonx++;
                            }
                        }
                        // ALT/.
                        else if (bcf_gt_is_missing(ptr[1]) && bcf_gt_allele(ptr[0]) == 1) {
                            if(MARKER_SCHEME > 1)
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[1], zero);
                            else {
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                                if(is_x_chr)
                                    zeroed_genotypes_x++;
                                else
                                    zeroed_genotypes_nonx++;
                            }
                        }
                        // REF/.
                        else if (bcf_gt_is_missing(ptr[1])) {
                            if(MARKER_SCHEME > 1)
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], zero);
                            else {
                                set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                                if(is_x_chr)
                                    zeroed_genotypes_x++;
                                else
                                    zeroed_genotypes_nonx++;
                            }
                        }
                        // ALT/REF
                        else if (bcf_gt_allele(ptr[0]) != bcf_gt_allele(ptr[1]))
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[1]);
                        // ALT/ALT
                        else if (bcf_gt_allele(ptr[0]) == 1)
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[1], canons[1]);
                        // REF/REF
                        else
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[0]);
                        // haploid
                    } else if (j == 1) {
                        // single missing genotypes
                        if (bcf_gt_is_missing(ptr[0]))
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                        // single ALT
                        else if (bcf_gt_allele(ptr[0]) == 1)
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[1], canons[1]);
                        // single REF
                        else
                            set_2Ralleles(persons[p].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[0]);
                    } else error("FIXME: not ready for ploidy %d\n", j);

                }

                canons.clear();
                Vecc().swap(canons);
                free(dat);

                mrkindex++;
            }
            if(line) bcf_clear(line);
        }
        int err = bcfargs->files->errnum;
        if ( err ) {
            fprintf(stderr,"Error: %s\n", bcf_sr_strerror(bcfargs->files->errnum));
            errorf("BCFTools encountered a problem with the given command.\n"
                   "Mega2 was unable to process allele labels and is exiting.\n"
                   "Please verify your BCFTools command is valid.\n");
        }

        args.pop_back();

        std::cout << "Genotypes for " << files[i] << " completed\n";

        if(zeroed_genotypes_nonx + zeroed_genotypes_x > 0)
            warnvf("We encountered %d half-typed genotypes that were set to missing.\n"
                   "If you would like to read these in, set the maximum number of alleles per marker in the 'File Input Menu' to more than 2.\n"
                   "Also, there were %d half-typed genotypes in males on chromosome X that were treated as valid.\n",zeroed_genotypes_nonx,zeroed_genotypes_x);

        args.clear();
        std::vector<std::string>().swap(args);
        destroy_data_vcfview(bcfargs);
        free(bcfargs);
    }
}

linkage_ped_top* ReadBCFs::do_ped(linkage_locus_top *LTop)   {
    extern vector<Vecc> VecAlleles;

    //mssgvf("\nAs a pedigree (.fam) file was not provided, we have assumed everyone is unrelated.\n"
    //       "If you have pedigree information that you did not include please rerun providing a pedigree file.\n"
    //       "All sex values have been set to male as a default.\n");
    annotated_ped_rec *persons = build_bcf_ped(LTop);

    do_genotypes(LTop, persons, VecAlleles, num_samples);
    linkage_ped_top *Top;

    Top = mk_ped_top(persons, this->num_samples, LTop, this->num_samples,
            /*untyped*/ 0, /*totaltyped*/ 0,
            /*groups*/ NULL, 0, 0,
            /*num_err*/0, 1);

    return Top;
}

annotated_ped_rec * ReadBCFs::build_bcf_ped(linkage_locus_top *LTop) {
    annotated_ped_rec *persons;
    int num_ped_records = this->num_samples;
    if ((persons = CALLOC((size_t) num_ped_records, annotated_ped_rec)) == NULL) {
        errorf("Build_BCF_Ped: Could not allocate enough memory, exiting.");
        EXIT(MEMORY_ALLOC_ERROR);
    }

    for (int i = 0; i < this->num_samples; i++) {
        annotated_ped_rec *entry = &persons[i];

        entry->pheno = (LTop->PhenoCnt > 0) ?
                       CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;

        entry->marker = (LTop->MarkerCnt > 0) ?
                        marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;

        entry->rec_num = i + 1;  // really p+1 now

        entry->Pedigree = canonicalColName(this->samples[i].c_str());
        entry->ID = canonicalColName(this->samples[i].c_str());

        PLINK.individuals += 1;
        entry->per_index = i;
        entry->ped_index = i + 1;

        entry->Father = canonicalColName("0");
        entry->Mother = canonicalColName("0");

        PLINK.unspecified_sex += 1;
        entry->Sex = '1';


        entry->LinkPerID = i;
        entry->LinkPedID = i;

        entry->PedID = entry->Pedigree;
        entry->PerID = entry->ID;

    }

    return persons;
}

void ReadBCFs::do_gc() {

    hdrMap.clear();
    Hmapsi().swap(hdrMap);

    samples.clear();
    vector<string>().swap(samples);

    filelist.clear();
    vector<string>().swap(filelist);

    for (MarkerVectorP bpp = markers.cbegin(); bpp != markers.cend(); bpp++) {
        delete *bpp;
    }
    markers.clear();
    MarkerVector().swap(markers);
}
