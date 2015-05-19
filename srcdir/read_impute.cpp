/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

#include "str_utils.hh"
#include "input.hh"

// g++ does not like these const operators on "vectordb" classes
#define cbegin() begin()
#define cend()   end()

using namespace std;

void ReadImputed::do_menu_pr(int &idx, int line_len, int choiceA[])
{
    printf("%2d) %-*s%s\n", idx, line_len, "Enter imputation chromosome [required]:",
           Mega2BatchItemGet("Imputed_Chromosome")->items_read ? 
           Mega2BatchItemGet("Imputed_Chromosome")->value.name : "--");
    choiceA[idx++] = imputed_chromosome_i;

    printf("%2d) %-*s%.4f\n", idx, line_len, "Enter imputation info metric threshold:", 
           Mega2BatchItemGet("Imputed_Info_Metric_Threshold")->value.fvalue);
    choiceA[idx++] = imputed_info_metric_threshold_i;

    printf("%2d) %-*s%.4f\n", idx, line_len, "Enter probability threshold:", 
           Mega2BatchItemGet("Imputed_Probability_Threshold")->value.fvalue);
    choiceA[idx++] = imputed_probability_threshold_i;

}

int ReadImputed::do_menu_parse(int choice_)
{
    int ret = 0;
    if (choice_ == imputed_chromosome_i) {
        char select[100], *sel = &select[0];
        int chrm;
        while (1) {
            printf("Please specify a chromosome for unspecified markers: ");
            fcmap(stdin, "%s", select); newline;
            if ((chrm = STR_CHR(select)) == -1) {
                printf("%s is not a valid chromosome. Valid chromosomes are numbers 1-26, X, Y, XY, or MT\n", select);
            } else break;
        }
        BatchItemSet(sel, "Imputed_Chromosome");
        ret = 1;
    } else if (choice_ == imputed_info_metric_threshold_i) {
	double ansd;
	while (1) {
	    printf("Please enter threshold for acceptable imputed value ");
	    fcmap(stdin, "%g", &ansd); newline;
	    if (ansd < 0.0 || ansd > 1.0) {
		printf("threshold must be between 0.0 and 1.0\n");
	    } else break;
	}

	BatchItemSet(ansd, "Imputed_Info_Metric_Threshold");
	ret = 1;
    } else if (choice_ == imputed_probability_threshold_i) {
	double ansd;
	while (1) {
	    printf("Please enter threshold for acceptable probabilities ");
	    fcmap(stdin, "%g", &ansd); newline;
	    if (ansd < 0.0 || ansd > 1.0) {
		printf("threshold must be between 0.0 and 1.0\n");
	    } else break;
	}

	BatchItemSet(ansd, "Imputed_Probability_Threshold");
	ret = 1;
    }
    return ret;
}

void ReadImputed::do_menu2batch()
{
    Cstr Values[] = { "Imputed_Chromosome",
                      "Imputed_Info_Metric_Threshold", 
                      "Imputed_Probability_Threshold" };
                      
    for(int i = 0; i < 2; i++) {
        batch_item_type *bip = Mega2BatchItemGet(Values[i]);
        if (bip->items_read) 
            batchf(bip);
    }
}

void ReadImputed::do_batch2local()
{
    BatchItemGet(this->default_chrm, "Imputed_Chromosome");
    BatchItemGet(this->info_threshold, "Imputed_Info_Metric_Threshold");
    BatchItemGet(this->probability_threshold, "Imputed_Probability_Threshold");
}


void ReadImputed::do_init(Input_Impute *inp)
{
    Str t = "";
    this->input = inp;
    this->files(*inp->input_files.bedfl, t, *inp->input_files.pedfl);

///    err_fn = mega2_input_files[6];

    read_imputed_file();

    if (this->read_info) read_info_file();

    read_sample_file();
}

linkage_locus_top *ReadImputed::do_names()
{
    linkage_locus_top *LTop = build_impute2_names();

    return LTop;
}

void ReadImputed::do_map(std::vector<m2_map>& additional_maps)
{
    m2_map impute_map("IMPUTE", 'p');

    build_impute2_map(impute_map);

    additional_maps.push_back(impute_map);
}

linkage_ped_top *ReadImputed::do_ped(linkage_locus_top *LTop)
{
    int num_peds = 0;

    annotated_ped_rec *persons = build_impute2_ped(LTop, &num_peds);

    build_impute2_genotypes(LTop, persons);

    linkage_ped_top *Top;
    Top = mk_ped_top(persons, this->people.size(), LTop, num_peds,
                     /*totaltyped*/ this->people_filtered,
                     /*groups*/ NULL, 0, 0, 
                     /*num_err*/0, 1);

    return Top;

}

void ReadImputed::read_imputed_file ()
{
//  asm("int $3");
    int dbg = 0;

    ifstream ifs;
    ifstream infs;
    ifs.open(impute_file);
    if (! ifs.is_open() ) {
        errorvf("read_imputed_file: Can not open \"%s\" file\n", impute_file);
        EXIT(1);
    }

    if (!info_file.empty()) {
        infs.open(info_file.c_str());
        if (! infs.is_open() ) {
            warnvf("read_imputed_file: Info file explicitly specified but can not be opened: \"%s\"\n", C(info_file));
            EXIT(1);
        }
        read_info = true;
    } else {
        info_file = string(impute_file) + "_info";
        infs.open(info_file.c_str());
        if (! infs.is_open() ) {
            warnvf("read_imputed_file: Implicitly guessed Info file can not be opened: \"%s\"\n", C(info_file));
            read_info = false;
        } else {
            infs.close();
            read_info = true;
        }
    }

    Str hmm, rsid, pos, A, B;
    VecsDB fields;
    int idx = 0;

    Str  name;
    Str  chrm;
    SUPPRESS_MSSG_NESTED_INIT(bad_line_msg);
    while (! ifs.eof() ) {

        ifs >> hmm;
        ifs >> rsid;
        ifs >> pos;
        ifs >> A;
        ifs >> B;
        ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        if (ifs.eof()) break;
        idx++;

        if (dbg) {
            cout << idx << ": ";
            cout << hmm << " ";
        }
        if (inMap(hmm, input->G.chrm_set)) {
            chrm = hmm;
            fields.clear();
            split(fields, rsid, ":", 3);
            name = fields[0];  // it seems to be this way
        } else if (hmm == "---") {
            fields.clear();
            split(fields, rsid, ":", 3);

            if (dbg) {
                cout << rsid << " ";
                cout << "#" << fields.size() << " ";
            }
            if (fields.size() > 1) {
                if (inMap(fields[0], input->G.chrm_set)) {
                    chrm = fields[0];
                    name = "chr" + fields[0] + "_" + fields[1];
                } else if (fields[0].compare(0, 2, "rs") == 0) {
                    chrm  =  default_chrm;
                    name = fields[0];
                } else {
                    SUPPRESS_MSSG_NESTED(bad_line_msg);
                    warnvf("impute2 file: bad line(%d): %s %s %s\nrs_id field first item unexpected\n", 
                           C(idx), C(hmm), C(rsid), C(pos));
                    continue;
                }
                if (dbg) {
                    cout << fields[0] << " ";
                    cout << fields[1] << " ";
                    cout << fields[2] << " ";
                    cout << fields[3] << " ";
                }
                if (fields[1] != pos) {
                    SUPPRESS_MSSG_NESTED(bad_line_msg);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id pos field (%s) does not match position column (%s)\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[1]), C(pos));
                    }
                if (fields[2] != A) {
                    SUPPRESS_MSSG_NESTED(bad_line_msg);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id A allele field (%s) does not match A column (%s)\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[2]), C(A));
                    }
                if (fields[3] != B) {
                    SUPPRESS_MSSG_NESTED(bad_line_msg);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id B allele field (%s) does not match B column (%s)\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[3]), C(B));
                    }
            }
        } else {
            SUPPRESS_MSSG_NESTED(bad_line_msg);
            warnvf("impute file: bad line(%d): %s %s %s %s %s\n         first field is not --- or chromosome\n",
                   C(idx), C(hmm), C(rsid), C(pos), C(A), C(B));
        }
        if (dbg) {
            cout << pos << " ";
            cout << A << " ";
            cout << B << " ";
            cout << endl;
        }
        markers.push_back(new ImpMarker(name, chrm, pos, A, B, read_info));
    }

    markers_all = markers.size();
    markers_filtered = markers_all;
    SUPPRESS_MSSG_NESTED_FINI(bad_line_msg);

}

Str ReadImputed::info_file_hdr = "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type";

void ReadImputed::read_info_file ()
{
//  asm("int $3");
    int dbg = 0;

    ifstream ifs;

    ifs.open(info_file.c_str());
    if (! ifs.is_open() ) {
        errorvf("read_info_file: Can not open \"%s\" file\n", C(info_file));
        EXIT(1);
    }

    Str line;
    getline(ifs, line);
    if (line.compare(0, info_file_hdr.size(), info_file_hdr) != 0) {
        errorvf("Bad header for \"%s\" file.\n", C(info_file));
        errorvf(" expecting: %s\n", C(info_file_hdr));
        errorvf(" found: %s\n", C(line));
        EXIT(1);
    }

    Vecs        fields;
    Vecmarkerpp mpp;
    size_t line_n = 1;
    int    skip_count = 0;
    double info;

    SUPPRESS_MSSG_NESTED_INIT(info_threshold_msg);
    SUPPRESS_MSSG_NESTED_INIT(indel_msg);
    SUPPRESS_MSSG_NESTED_INIT(dup_msg);
    ImpMarker *mp;
    for (mpp = markers.cbegin(); ! ifs.eof(); mpp++) {
        mp = *mpp;
        getline(ifs, line);
        if (ifs.fail()) break;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        line_n++;
        fields.clear();
        split(fields, line);
        if (fields[2] != mp->pos) {
            errorvf("Files \"%s\" and \"%s\" do not list markers in the identical order starting at line %d: %s.\n", 
                    impute_file, C(info_file), line_n, C(line));
                EXIT(1);
        }
        info = atof(fields[6].c_str());
        mp->info = info;
        mp->certainty = atof(fields[7].c_str());

        if (info < info_threshold) {
            skip_count++;
            mp->skip = true;

            SUPPRESS_MSSG_NESTED(info_threshold_msg);
            warnvf("Marker: %s %s %s %s %s info (%.4f) < threshold (%.4f)\n", 
                   C(fields[0]), C(fields[1]), C(fields[2]), C(fields[3]), C(fields[4]),
                   info, info_threshold);
        } else if ( mp->A.find(mp->B)!=string::npos || mp->B.find(mp->A)!=string::npos ) {
            SUPPRESS_MSSG_NESTED(indel_msg);
            warnvf("Markers: %s (bp %s) indel alleles[(%s) (%s)] ignored.\n",
                   C(mp->name), C(mp->pos), C(mp->A), C(mp->B));
            skip_count++;
            mp->skip = true;
        } else if (mpp == markers.cbegin())
            ;
        else if ((*(mpp-1))->pos == mp->pos) {
            SUPPRESS_MSSG_NESTED(dup_msg);
            warnvf("Markers: %s (bp %s) repeated with different alleles [(%s/%s) (%s/%s)] ignored.\n",
                   C(mp->name), C(mp->pos), C((*(mpp-1))->A), C((*(mpp-1))->B), C(mp->A), C(mp->B));
            skip_count++;
            mp->skip = true;
        }
    }

    markers_filtered = markers_all - skip_count;
    SUPPRESS_MSSG_NESTED_FINI(info_threshold_msg);
    SUPPRESS_MSSG_NESTED_FINI(indel_msg);
    SUPPRESS_MSSG_NESTED_FINI(dup_msg);

    if (mpp != markers.cend() || line_n != markers.size() + 1 /*hdr*/) {
        errorvf("Files \"%s\" and \"%s\" are different lengths: %d vs %d\n",
                impute_file, C(info_file), markers.size(), line_n);
        EXIT(1);
    }
    ifs.close();

    warnvf("%d of %d markers that are filtered by  the info_metric_threshold or indels.\n", 
           skip_count, markers.size());
}

Str ReadImputed::sample_file_hdr = "ID_1 ID_2 missing";

void ReadImputed::read_sample_file ()
{
//  asm("int $3");
    int dbg = 0;

    ifstream ifs;

    ifs.open(sample_file.c_str());
    if (! ifs.is_open() ) {
        errorvf("read_sample_file: Can not open \"%s\" file\n", C(sample_file));
        EXIT(1);
    }

    Str  line;
    int  line_n = 2;
    VecsDB fields;
    VecsDB mappedfields;

//line 1
    getline(ifs, line);
    Vecs lineV, hdrV;
    split(lineV, line);
    split(hdrV,  sample_file_hdr);
    if (lineV[0] != hdrV[0] || lineV[1] != hdrV[1] || lineV[2] != hdrV[2]) {
        errorvf("Bad header for \"%s\" file.\n", C(sample_file));
        errorvf(" expecting: %s\n", C(sample_file_hdr));
        errorvf(" found: %s\n", C(line));
        EXIT(1);
    }
    split(sample_file_hdr1a, line);
    split(fields, line);
//  column_type.reserve(fields.size());

    Mapsi  col2idx;
    Mapsip col2idxp;
    Mapii  idx2fixed;
    Mapiip idx2fixedp;
    int    i;
    int    ret;
    i = 0;
    for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
        col2idx[*fp]    = i;
//      column_type.push_back(PHENO);
    }
#if 0
    Cstr special[] = {"ID_1", "ID_2", "father", "mother", "sex"};
    for (i = 0; i < 5; i++) {
        if (map_get(col2idx, special[i], ret)) {
            column_type[ret] = RESERVED;
            reserved.push_back(ret);
        } else {
            reserved.push_back(-1);
        }
    }

//line 2
    getline(ifs, line);
    if (ifs.fail()) {
        errorvf("Bad second header line for \"%s\" file.\n", C(sample_file));
        EXIT(1);
    }
    if (dbg) {
        cout << line_n << ": " << line << endl;
    }
    split(sample_file_hdr2a, line);

//line rest
    for (; ! ifs.eof(); ) {
        getline(ifs, line);
        if (ifs.fail()) break;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        line_n++;
        fields.clear();
        split(fields, line);
        mappedfields.clear();
        for (i = 0; i < 5; i++) {
            if (reserved[i] < 0)
                mappedfields.push_back("0");
            else
                mappedfields.push_back(fields[reserved[i]]);
        }

        i = 0;
        for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
            if (column_type[i] != RESERVED)
                mappedfields.push_back(*fp);
        }
        people.push_back(mappedfields);
    }
#else
    Cstr special[] = {"ID_1", "ID_2", "father", "mother", "sex", "missing"};
    for (i = 0; i < 6; i++) {
        if (map_get(col2idx, special[i], ret))
            idx2fixed[ret] = i;
    }

//line 2
    sample_file_hdr1b.clear();
    Cstr fill1[] = {"", "", "0", "0", "0"};
    sample_file_hdr1b.insert(sample_file_hdr1b.end(), fill1, fill1+5);
    i = 0;
    for (VecspDB fp = sample_file_hdr1a.cbegin(); fp != sample_file_hdr1a.cend(); i++, fp++) {
        if (map_get(idx2fixed, i, ret)) {
            if (ret < 5)
                sample_file_hdr1b[ret] = *fp;
        } else
            sample_file_hdr1b.push_back(*fp);
    }

    getline(ifs, line);
    if (ifs.fail()) {
        errorvf("Bad second header line for \"%s\" file.\n", C(sample_file));
        EXIT(1);
    }
    if (dbg) {
        cout << line_n << ": " << line << endl;
    }
    split(sample_file_hdr2a, line);
    Cstr fill2[] = {"", "", "", "", ""};
    sample_file_hdr2b.insert(sample_file_hdr2b.end(), fill2, fill2+5);
    i = 0;
    for (VecspDB fp = sample_file_hdr2a.cbegin(); fp != sample_file_hdr2a.cend(); i++, fp++) {
        if (map_get(idx2fixed, i, ret)) {
            if (ret < 5)
                sample_file_hdr2b[ret] = *fp;
        } else
            sample_file_hdr2b.push_back(*fp);
    }

//line rest
    for (; ! ifs.eof(); ) {
        getline(ifs, line);
        if (ifs.fail()) break;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        line_n++;
        fields.clear();
        split(fields, line);
        
        mappedfields.clear();
        Cstr fill[] = {"", "", "0", "0", "0"};
        mappedfields.insert(mappedfields.end(), fill, fill+5);
        i = 0;
        int ret = 0;
        for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
            if (map_get(idx2fixed, i, ret)) {
                if (ret < 5)
                    mappedfields[ret] = *fp;
            } else {
                mappedfields.push_back(*fp);
            }
        }
        people.push_back(mappedfields);
    }
#endif

    people_all = people.size();
    people_filtered = people_all - 0;
}

linkage_locus_top *
ReadImputed::build_impute2_names()
{
//  asm("int $3");
//  int dbg = 0;

    int num_pheno   = sample_file_hdr1b.size() - 5;
    int num_markers = markers_filtered;
    int num_all     = num_pheno + num_markers;

    char **names    = CALLOC(num_all, char *);
    char  *types    = CALLOC(num_all, char);
    int i;
    i = 0;
    Vecsp typep = sample_file_hdr2b.cbegin()+5;
    for (Vecsp phep = sample_file_hdr1b.cbegin()+5; i < num_pheno; i++, phep++, typep++) {
        names[i] = CALLOC((*phep).size()+1, char);
        strcpy(names[i], (*phep).c_str());
        types[i] = (*typep) == "P" ? 'T' : 'A';
    }

    ImpMarker *mp;
    for (Vecmarkerpp marp = markers.cbegin(); marp != markers.cend(); marp++) {
        mp = *marp;
        if (! mp->skip ) {
            if (i >= num_all) {
                errorvf("Internal Error: build_impute2_names() count of skipped markers too large\n");
                EXIT(1);
            }
            names[i] = CALLOC(mp->name.size()+1, char);
            strcpy(names[i], mp->name.c_str());
            types[i] = 'M';
            i++;
        }
    }

    return read_common_marker_data(num_pheno + num_markers, num_markers, names, types,
                                  /*annotated*/ 0, /*penetrances_read*/ 0, 0, NULL);
}

void ReadImputed::build_impute2_map(m2_map& map)
{
//  asm("int $3");
//  int dbg = 0;

    human_x = human_y = human_xy = human_unknown = human_mt = human_auto = 0;

    int num_pheno   = sample_file_hdr1b.size() - 5;
    int num_markers = markers_filtered;
    int num_all     = num_pheno + num_markers;

    int i;
    i = 0;

    double pos;
    ImpMarker *mp;
    m2_map_entry map_entry;
    for (Vecmarkerpp marp = markers.cbegin(); marp != markers.cend(); marp++) {
        mp = *marp;
        if (! mp->skip ) {
            if (i >= num_all) {
                errorvf("Internal Error: build_impute2_map() count of skipped markers too large\n");
                EXIT(1);
            }
            i++;
            map_entry.set_chr(mp->chr);
            pos = atof(C(mp->pos));
            map_entry.set_POS(pos);
            map_entry.set_REF(mp->A);
            map_entry.set_marker_name(mp->name);
            map.push_back_entry(map_entry);

            int chr = STR_CHR(C(mp->chr));
            if (chr == SEX_CHROMOSOME) {
                human_x++;
    //          if (LTop->Locus[mrk_num].Type != XLINKED) {
    //              LTop->Locus[mrk_num].Type = XLINKED;
    //          }
            } else if (chr == MALE_CHROMOSOME) {
                human_y++;
    //          if (LTop->Locus[mrk_num].Type != YLINKED) {
    //              LTop->Locus[mrk_num].Type = YLINKED;
    //          }
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
}

annotated_ped_rec *
ReadImputed::build_impute2_ped(linkage_locus_top *LTop, int *num_peds) 
{
//    asm("int $3");
    int dbg = 0;

    int p = 0;
    int p_idx;
    annotated_ped_rec *persons;
    int num_ped_records = people.size();
    if ((persons = CALLOC((size_t)num_ped_records, annotated_ped_rec)) == NULL) {
        errorf("Build_Impute2_Ped: Could not allocate enough memory, exiting.");
        EXIT(MEMORY_ALLOC_ERROR);
    }
    Mapsi ped_names;
    int curr_per_index = 0;
    int curr_ped_index = 0;
    int num_errors     = 0;

    SUPPRESS_MSSG_NESTED_INIT(illegal_affect_msg);
    SUPPRESS_MSSG_NESTED_INIT(illegal_quant_msg);
    for (Vecvecsp peop = people.cbegin(); peop != people.cend(); peop++) {
        const VecsDB& pp = (*peop);

        annotated_ped_rec *entry = &persons[p++]; // incr p for next cycle

        entry->pheno  = (LTop->PhenoCnt > 0) ? 
            CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;

        // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
        entry->marker = (LTop->MarkerCnt > 0) ? 
            marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;

        entry->rec_num = p;  // really p+1 now

        if (dbg)
            printf("%s %s %s %s %s\n  ", C(pp[0]), C(pp[1]), C(pp[2]), C(pp[3]), C(pp[4]));

// 0 ID_1 -> Pedigree
        strcpy(entry->Pedigree, C(pp[0]));
        if (map_get(ped_names, pp[0], p_idx)) {
            entry->ped_index = p_idx;
            curr_per_index++;
        } else {
            curr_ped_index++;
            entry->ped_index = curr_ped_index;
            ped_names[pp[0]] = curr_ped_index;
            curr_per_index = 1;
        }
        entry->per_index = curr_per_index;

// 1 ID_2 -> ID
        strcpy(entry->ID, C(pp[1]));
        PLINK.individuals += 1;
#ifdef IDS
        printf("Ped/per: %s/%s; pedi/peri %d/%d\n",
               entry->Pedigree, entry->ID, curr_ped_index, curr_per_index);
#endif

// 2,3 father,mother
        strcpy(entry->Father,  C(pp[2]));
        strcpy(entry->Mother,  C(pp[3]));

// 4 sex
        entry->Sex = pp[4][0];
        if (tolower((unsigned char)entry->Sex) == 'm' || entry->Sex == '1') {
            PLINK.males += 1;
        } else if (tolower((unsigned char)entry->Sex) == 'f' || entry->Sex == '2') {
            PLINK.females += 1;
        } else {
            PLINK.unspecified_sex += 1;
            entry->Sex = '0' ;
        }

	entry->LinkPerID = curr_per_index;
        entry->LinkPedID = curr_ped_index;

	strcpy(entry->PedID, entry->Pedigree);
	strcpy(entry->PerID, entry->ID);

#ifdef IDS
	printf("Ped/per: %s/%s; Lped/Lperi %d/%d; %d; %s\n",
	       entry->PedID, entry->PerID, entry->LinkPedID, entry->LinkPerID, curr_per_index,
	       entry->ID);
#endif


// 5+ pheno
        int ret;
        int i;
        i = 0;
        Vecsp namep = sample_file_hdr1b.cbegin()+5;
        Vecsp typep = sample_file_hdr2b.cbegin()+5;
        for (VecspDB vp = pp.cbegin()+5; vp != pp.cend(); vp++, i++) {
            if (dbg)
                printf("%s %s %s\n", C(namep[i]), C(typep[i]), C(*vp));
            switch (typep[i] == "P" ? 'T' : 'A') {
            case 'A':
                ret = plink_annot_string_aff_phen(entry->rec_num+HDR, &(LTop->Pheno[i]),
                                                  &entry->pheno[i], C(*vp));
                if (!ret) {
                    SUPPRESS_MSSG_NESTED(illegal_affect_msg);
                }
                if (entry->pheno[i].Affection.Status == UNDEF) {
                    SUPPRESS_MSSG_NESTED(illegal_affect_msg);
                    errorvf("File %s, Entry %d : Invalid status at locus %s\n",
                            C(sample_file), entry->rec_num+HDR, LTop->Locus[i].Name);
                    num_errors++;
                }
                if (LTop->Pheno[i].Props.Affection.ClassCnt != 1) {
                    if (entry->pheno[i].Affection.Class == UNDEF) {
                        SUPPRESS_MSSG_NESTED(illegal_affect_msg);
                        errorvf("File \"%s\", Entry %d : Invalid liability class at locus %s\n\n",
                                C(sample_file), entry->rec_num+HDR, LTop->Locus[i].Name);
                        num_errors++;
                    }
                }
                break;
            case 'T':
                ret = plink_annot_string_quant_phen(entry->rec_num+HDR, &(LTop->Pheno[i]),
                                                    &entry->pheno[i], C(*vp));
                if (!ret) {
                    SUPPRESS_MSSG_NESTED(illegal_quant_msg);
                }
                // Here we are not checking for it being undefined, just "invalid".
                // There should be a better way...
                if (entry->pheno[i].Quant == QUNDEF) {
                    SUPPRESS_MSSG_NESTED(illegal_quant_msg);
                    errorvf("File \"%s\", Entry %d : Invalid quantitative phenotype at locus %s\n",
                            C(sample_file), entry->rec_num+HDR, LTop->Locus[i].Name);
                    num_errors++;
                }
                break;
            }
        }
    }
    SUPPRESS_MSSG_NESTED_FINI(illegal_affect_msg);
    SUPPRESS_MSSG_NESTED_FINI(illegal_quant_msg);
    printf("num errors: %d\n", num_errors);

#ifdef SHOWSTATUS
    int fudge;
    if (MARKER_SCHEME == MARKER_SCHEME_BITS)
        fudge = LTop->LocusCnt * (sizeof (Alleles_int) + sizeof (Alleles_str));
    else
        fudge = 0;
    msgvf("ALLOC SPACE: ALL pedrec: %d MB (%d x %d)\n",
          (p * marker_size(LTop->MarkerCnt) + fudge) / 1024 / 1024,
          p,  marker_size(LTop->MarkerCnt));
#endif

    *num_peds = curr_ped_index;
    return persons;
}

void ReadImputed::build_impute2_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons)
{
//  asm("int $3");
    int dbg = 0;

    ifstream ifs;
    string line;

    string hmm, rsid, pos, A, B;
//    vector<double> nums;
    double nums[3];

    ifs.open(impute_file);
    if (! ifs.is_open() ) {
        errorvf("read_imputed_genotype_file: Can not open \"%s\" file\n",
                impute_file);
        EXIT(1);
    }

    int line_n = 0;
    Token token(3);
//  Str   genotype[3];
    Vecmarkerpp mpp;
    ImpMarker   *mp;

    int   mrk_idx = sample_file_hdr2b.size() - 5 -1;
    linkage_locus_rec *locus;

    char *callele1;
    char *callele2;
    char *callele0  = canonical_allele(C("0"));

    SUPPRESS_MSSG_NESTED_INIT(skip_msg);
    SUPPRESS_MSSG_NESTED_INIT(prob_msg);

    Tod tod_gen("impute genotypes");
    for (mpp = markers.cbegin(); ! ifs.eof(); mpp++) {
        getline(ifs, line);
        if (ifs.eof()) break;
        line_n++;

        token.set(line);

//	ifs >> hmm; // +A+B+ nums() all read via >> @ 14.62 sec
//                               vector<double> nums  12.26
//                               double nums[3]       11.59
        token.more(hmm);
        token.more(rsid);
        token.more(pos);
        mp = *mpp;
        if (pos != mp->pos) {
            errorvf("internal error: impute_file (\"%s\") second pass does not match first pass at line %d\n",
                    impute_file, line_n);
                EXIT(1);
        }
        if (mp->skip) {
//          SUPPRESS_MSSG_NESTED(skip_msg);
//          warnvf("Marker: skipped %s %s %s info (%.4f) < threshold (%.4f)\n", 
//                 C(mp->name), C(mp->chr), C(mp->pos), mp->info, info_threshold);
            continue;
        }
        mrk_idx++;  // hence the -1 above
        locus = &LTop->Locus[mrk_idx];
// ?? totaltyped++

        token.more(A);
        token.more(B);
        callele1    = canonical_allele(C(A));
        callele2    = canonical_allele(C(B));
//      genotype[0] = A + A;
//      genotype[1] = A + B;
//      genotype[2] = B + B;
        
        int i;
        int sam = 0;
        double maxx;
        int p = 0;

	for( ; p < people_filtered; p++) {
            annotated_ped_rec *entry = &persons[p];
            sam++;
            token.getDC(nums);

            maxx = 0.0;

            if (dbg) {
                cout << "#";
                cout << line_n << ": ";
                cout << mp->name << " ";
                cout << nums[0] << " " << nums[1] << " " << nums[2] << "\n";
            }

            if (nums[0] > nums[1]) {
                maxx = nums[0];
                i = 0;
            } else {
                maxx = nums[1];
                i = 1;
            }
            if (nums[2] > maxx) {
                maxx = nums[2];
                i = 3;
            }

            if (maxx == 0)
                i = 2;
            else if (maxx < probability_threshold) {
                SUPPRESS_MSSG_NESTED(prob_msg);
                warnvf("Marker: %s %s %s max probability (%.4f) < probability threshold (%.4f)\n", 
                       C(mp->name), C(mp->chr), C(mp->pos), maxx, probability_threshold);
                i = 2;
            }

            switch (i) {
            case 0:
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele1);
                break;
            case 1:
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele2);
                break;
            case 2:
                set_2Ralleles(entry->marker, mrk_idx, locus, callele0, callele0);
                break;
            case 3:
                set_2Ralleles(entry->marker, mrk_idx, locus, callele2, callele2);
                break;
            }

        }
/*
        token.getD(nums);  // last one on line
        if (nums.size() != 0) {
            errorvf("internal error: impute_file (\"%s\") second pass does not match first pass at line %d for number of genotypes\n",
                    impute_file, line_n);
        }
*/
    }
    tod_gen();
    SUPPRESS_MSSG_NESTED_FINI(skip_msg);
    SUPPRESS_MSSG_NESTED_FINI(prob_msg);
}
