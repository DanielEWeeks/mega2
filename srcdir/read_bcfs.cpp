/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2018 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

using namespace std;

/*
 * Displays menu, adds additional options in menu 1 after input suffix
 * this gets called in user_input.cpp for whatever input opt is defined
 */
void ReadBCFs::do_menu_display(int &idx, int line_len, int choiceA[]) {
    printf("%2d) %-*s%s\n", idx, line_len, "BCFtools Parameters", BatchItemGet("BCF_Args")->value.name);
    choiceA[idx] = site_bcfs_args_i;
    idx++;

    printf("%2d) %-*s%s%s\n", idx, line_len-11, "Variant File:","[required] ", BatchItemGet("BCFs_File")->value.name);
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
            printf("Enter a BCF file or a file that contains a list of BCF files to be read:\n");
            fcmap(stdin, "%s", bcfs_file);
            newline;

            if (access(bcfs_file, F_OK)) {
                printf("WARNING: Could not find file %s\n", bcfs_file);
                continue;
            }
            else {
                BatchValueSet(bcfs_file, "BCFs_File");
                break;
            }
        }
    }

    else if(choice_ == site_bcfs_args_i) {
        draw_line();
        printf("\nCurrent BCF parameters:  %s\n", BatchItemGet("BCF_Args")->value.name);
        printf("  For more additional information on BCFTools flags\n");
        printf("see the documentation at samtools.github.io/bcftools/bcftools\n");
        printf("Valid options in Mega2 include:\n");
        printf("--known  --novel --phased --exclude-phased --uncalled --exclude-uncalled\n");
        printf("--min-ac --max-ac --min-alleles --max-alleles[INT]\n");
        printf("--min-af --max-af [FLOAT]");
        printf("--exclude --include [EXPRESSION]\n");
        printf("--regions [chr:to-from] --regions-file [FILE]\n");
        printf("--apply-filters [LIST]\n");

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
                mssgvf("BCFTools options must begin with \"--\"\n");
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
                argv[ii] = &args[ii][0];
            }

            int test = mbi->test_args(argc, argv);

            if(test == 1) {
                this->BCF_args = strdup(bcf_args);
                BatchValueSet(this->BCF_args, "BCF_Args");
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
            "BCFs_File"
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
              || this->BCF_file.substr(this->BCF_file.find_last_of(".") + 1) == "vcf.gz"){
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
    MEGA2_BCFTOOLS_INTERFACE *mbi = new MEGA2_BCFTOOLS_INTERFACE();

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
        auto start = std::chrono::system_clock::now();

        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for chromosome %d:\n",i+1);
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            argv[ii] = &args[ii][0];
            printf("%s ",argv[ii]);
        }
        printf("\n");

        args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);
        this->samples.push_back(hdr->samples[i]);
        this->num_samples = hdr->n[2];
        args.pop_back();

        auto end = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(end);
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Samples for chromosome " << i + 1 << " completed\n" << std::ctime(&time)
                  << "Duration: " << elapsed_seconds.count() << "\n";
    }

    args.pop_back();
    //no header since we have the samples
    args.push_back("-H");
    //drop genotypes to not read them while getting allele labels
    args.push_back("-G");
    argc++;



    int total_markers = 0;
    for(int i = 0; i < this->filecount; i++) {
        auto start = std::chrono::system_clock::now();

        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for chromosome %d:\n",i+1);
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            argv[ii] = &args[ii][0];
            printf("%s ",argv[ii]);
        }
        printf("\n");

        args_t *bcfargs  = (args_t*) calloc(1,sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        int count = 0;

        while ( bcf_sr_next_line(bcfargs->files) ) {
            bcf1_t *line = bcfargs->files->readers[0].buffer[0];
            if ( subset_vcf(bcfargs, line) ) {
                Vecs alleles(line->n_allele);
                alleles.clear();

                for (int al = 0; al < line->n_allele; al++) {
                    alleles.push_back(canonical_allele(line->d.allele[al]));
                }

                Str name = line->d.id;
//rvb: and another
                if (name == ".") {
                    char pos[50];
                    sprintf(pos, "chr%s_%d", hdr->id[BCF_DT_CTG][line->rid].key, line->pos + 1);
                    name = pos;
                }
                this->markers.push_back(new BCFMarker(name, hdr->id[BCF_DT_CTG][line->rid].key, line->pos + 1,
                                                      alleles)); // pos + 1 matches the VCF line pos field.

                count++;
            }
        }

        total_markers += count;
        args.pop_back();

        auto end = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(end);
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Alleles for chromosome " << i + 1 << " completed\n" << std::ctime(&time)
                  << "Duration: " << elapsed_seconds.count() << "\n";
    }

    this->marker_count = total_markers;
}

void ReadBCFs::do_map(std::vector<m2_map>& additional_maps)
{
    m2_map bcf_map("VCF (created from input)", 'p');

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

static void phenotype_file_SAMPLEID_entry_checks();

linkage_locus_top *ReadBCFs::do_names(const char *&names_fn)
{
    linkage_locus_top *LTop = build_BCFs_names();

    phenotype_file_SAMPLEID_entry_checks();

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

/**
 This routine processes the SAMPLEID information in the phenotype file, and the sample information
 in the VCF file and will give the user a message if any of the following occurs:
 1) a SAMPLEID in the Phenotype File is not found (or found multiple times) as a VCF file
    sample label (found in the VCF file header line)
 2) a SAMPLEID entry has been excluded through the vcftools filtering mechanism
 3) a sample label in the VCF file is not found in the SAMPLEID column of the phenotype file.
 */
extern sample_map_type *SAMPLEIDS;
static void phenotype_file_SAMPLEID_entry_checks()
{
    if (SAMPLEIDS == NULL) {
        mssgf("The SAMPLEID column was not included in the phenotype file.");
        return;
    }
    mssgf("The SAMPLEID column was included in the phenotype file.");
    
#if 0
    // Iterate over the SAMPLEID information from the phenotype file...
    for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
        string map_sample = it->first;
        int found = 0;
        
        // Look for the SAMPLEID (map_sample) in the VCF file...
        for (int ui=0; ui < vf->N_total_indv(); ui++) {
            string vcf_file_sample = vf->indv[(size_t)ui];
            if (vcf_file_sample == map_sample) {
                // It was found...
                if (vf->include_indv[(size_t)ui] == false && found == 0) {
                    mssgvf("SAMPLEID '%s' has been excluded by the vcftools filters.\n", map_sample.c_str());
                }
                found++; // SAMPLEID found in the VCF file...
            }
        }
        if (found == 0) {
            mssgvf("SAMPLEID '%s' was not found in the VCF file.\n", map_sample.c_str());
        } else if (found > 1) {
            mssgvf("SAMPLEID '%s' was found multiple times in the VCF file.\n", map_sample.c_str());
        } // else if (fount == 1) all is well!
    } // for (sample_map_type::iterator it ...
    
    // Look through all of the samples in the VCF file....
    for (int ui=0; ui < vf->N_total_indv(); ui++) {
        if (vf->include_indv[(size_t)ui] == false) continue;
        string vcf_file_sample = vf->indv[(size_t)ui];
        int found = 0;
        // Determine if it matches any sample in the SAMPLEID column of the phenotype file...
        for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
            string map_sample = it->first;
            if (vcf_file_sample == map_sample) {
                found++;
                break;
            }
        }
        if (found == 0) {
            mssgvf("The VCF file sample '%s' was not found in the SAMPLEID column\n", vcf_file_sample.c_str());
            mssgf("included in the phenotype file.");
        }
    }
#endif
}

void ReadBCFs::do_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons,
                            std::vector<Vecc> &VecAlleles) {
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

    for (int i = 0; i < this->filecount; i++) {
        auto start = std::chrono::system_clock::now();

        args.push_back(files[i]);

        printf("\nRunning the following bcftools command for chromosome %d:\n",i+1);
        char **argv;
        argv = (char **) malloc(argc * sizeof(char *));
        for (size_t ii = 0; ii < argc; ii += 1) {
            argv[ii] = (char *) malloc(FILENAME_LENGTH * sizeof(char));
            printf("%s ",argv[ii]);
        }
        printf("\n");

        args_t *bcfargs = (args_t *) calloc(1, sizeof(args_t));
        bcfargs = mbi->get_args(argc, argv);

        bcf_hdr_t *hdr = bcfargs->hnull ? bcfargs->hnull : (bcfargs->hsub ? bcfargs->hsub : bcfargs->hdr);

        while (bcf_sr_next_line(bcfargs->files)) {
            bcf1_t *line = bcfargs->files->readers[0].buffer[0];

            if ( subset_vcf(bcfargs, line) && mrkindex < LTop->LocusCnt) {
                int m, n, i;

                //uses convert object so we need our own void * instead
                void *dat = NULL;

                m = 0;
                n = bcf_get_genotypes(hdr, line, &dat, &m);
                //probably not necessary
                //convert->ndat = m * sizeof(int32_t);

                if (n <= 0) {
                    error("Error parsing GT tag at %s:%d\n", bcf_seqname(hdr, line), line->pos + 1);
                }



                Vecc canons;

                for (int al = 0; al < line->n_allele; al++) {
                    canons.push_back(canonical_allele(line->d.allele[al]));
                }
                VecAlleles.push_back(canons);

                //should give number of allele options per marker
                n /= num_samples;
                for (i = 0; i < num_samples; i++) {
                    int32_t *ptr = (int32_t *) dat + i * n;
                    int j;
                    for (j = 0; j < n; j++)
                        if (ptr[j] == bcf_int32_vector_end) break;

                    // diploid
                    if (j == 2) {
                        if (bcf_gt_is_missing(ptr[0]))
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                            //kputs(" 0.33 0.33 0.33", str);
                        else if (bcf_gt_allele(ptr[0]) != bcf_gt_allele(ptr[1]))
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[1]);
                            //kputs(" 0 1 0", str);       // HET
                        else if (bcf_gt_allele(ptr[0]) == 1)
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], canons[1], canons[1]);
                            //kputs(" 0 0 1", str);       // ALT HOM, first ALT allele
                        else
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[0]);
                        //kputs(" 1 0 0", str);       // REF HOM or something else than first ALT
                        // haploid
                    } else if (j == 1) {
                        if (bcf_gt_is_missing(ptr[0]))
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], zero, zero);
                            //kputs(" 0.5 0.0 0.5", str);
                        else if (bcf_gt_allele(ptr[0]) == 1)
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], canons[1], canons[1]);
                            //kputs(" 0 0 1", str);       // first ALT allele
                        else
                            set_2Ralleles(persons[i].marker, mrkindex, &LTop->Locus[mrkindex], canons[0], canons[0]);
                        //kputs(" 1 0 0", str);       // REF or something else than first ALT
                    } else error("FIXME: not ready for ploidy %d\n", j);

                }
                mrkindex++;
            }
        }
        args.pop_back();

        auto end = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(end);
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "Genotypes for chromosome " << i + 1 << " completed\n" << std::ctime(&time)
                  << "Duration: " << elapsed_seconds.count() << "\n";
    }
}


