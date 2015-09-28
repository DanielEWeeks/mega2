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
#include "input_ops.hh"

#include "read_impute.hh"

// g++ does not like these const operators on "vectordb" classes
#define cbegin() begin()
#define cend()   end()

using namespace std;

void ReadImputed::do_menu_display(int &idx, int line_len, int choiceA[])
{
    printf("%2d) %-*s%s\n", idx, line_len, "Oxford-single-chr [if required]:",
           BatchItemGet("Imputed_Oxford_Single_Chr")->value.name);
    choiceA[idx++] = imputed_oxford_single_chr_i;

    printf("%2d) %-*s%.4f\n", idx, line_len, "Imputation info metric threshold:", 
           BatchItemGet("Imputed_Info_Metric_Threshold")->value.fvalue);
    choiceA[idx++] = imputed_info_metric_threshold_i;

    printf("%2d) %-*s%.4f\n", idx, line_len, "\"Impute2 style\" hard call uncertainty:", 
           BatchItemGet("Imputed_Hard_Call_Uncertainty")->value.fvalue);
    choiceA[idx++] = imputed_hard_call_uncertainty_i;

    printf("%2d) %-*s%.4f\n", idx, line_len, "Genotype missing fraction:", 
           BatchItemGet("Imputed_Genotype_Missing_Fraction")->value.fvalue);
    choiceA[idx++] = imputed_genotype_missing_fraction_i;

    printf("%2d) %-*s", idx, line_len, "Missing trait value codes: (space separated)");
    char **lst = BatchItemGet("Imputed_Missing_Codes")->value.mult_names;
    if (lst) {
        for (; *lst; ) printf("%s ", *lst++);
    }
    printf("\n");
    choiceA[idx++] = imputed_missing_code_i;

    printf("%2d) %-*s[ %s ]\n", idx, line_len, "Allow indels (toggle):",
           BatchItemGet("Imputed_Allow_Indels")->value.copt == 'y' ? "yes" : "no");
    choiceA[idx++] = imputed_allow_indels_i;

    printf("%2d) %-*s[ %s ]\n", idx, line_len, "Allow duplicate markers (toggle):",
           BatchItemGet("Imputed_Allow_Duplicates")->value.copt == 'y' ? "yes" : "no");
    choiceA[idx++] = imputed_allow_dups_i;

}

int ReadImputed::do_menu_parse(int choice_)
{
    int ret = 0;
    if (choice_ == imputed_oxford_single_chr_i) {
        char select[100], *sel = &select[0];
        int chrm;
        while (1) {
            printf("Please specify a chromosome for unspecified markers: ");
            fcmap(stdin, "%s", select); newline;
            if ((chrm = STR_CHR(select)) == -1) {
                printf("%s is not a valid chromosome. Valid chromosomes are numbers 1-26, X, Y, XY, or MT\n", select);
            } else break;
        }

        BatchValueSet(sel, "Imputed_Oxford_Single_Chr");
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

	BatchValueSet(ansd, "Imputed_Info_Metric_Threshold");
	ret = 1;
    } else if (choice_ == imputed_hard_call_uncertainty_i) {
	double ansd;
	while (1) {
	    printf("Please enter hard call fraction ");
	    fcmap(stdin, "%g", &ansd); newline;
	    if (ansd < 0.0 || ansd > 1.0) {
		printf("uncertainty must be between 0.0 and 1.0\n");
	    } else break;
	}

	BatchValueSet(ansd, "Imputed_Hard_Call_Uncertainty");
	ret = 1;
    } else if (choice_ == imputed_genotype_missing_fraction_i) {
	double ansd;
	while (1) {
	    printf("Please enter genotype missing fraction ");
	    fcmap(stdin, "%g", &ansd); newline;
	    if (ansd < 0.0 || ansd > 1.0) {
		printf("uncertainty must be between 0.0 and 1.0\n");
	    } else break;
	}

	BatchValueSet(ansd, "Imputed_Genotype_Missing_Fraction");
	ret = 1;
    } else if (choice_ == imputed_missing_code_i) {
	char selection[100];
        if (1) {
	    printf("Please enter space separated list of missing trait value codes: ");
            IgnoreValue(fgets(selection, sizeof(selection)-1, stdin)); newline;
	}
        Str MissingCodes = string(selection);
        Vecs vec;
        split(vec, MissingCodes);
        join(vec, Input->MissingCodes, " ");

        Input->MissingCodesSet.clear();
        Input->MissingCodesSet.insert(vec.cbegin(), vec.cend());

	BatchValueSet(vec, "Imputed_Missing_Codes");
	ret = 1;
    } else if (choice_ == imputed_allow_indels_i) {
        batch_item_type *bip = BatchItemGet("Imputed_Allow_Indels");
        bip->value.copt = bip->value.copt == 'y' ? 'n' : 'y';

        bip->items_read = 1;
	ret = 1;
    } else if (choice_ == imputed_allow_dups_i) {
        batch_item_type *bip = BatchItemGet("Imputed_Allow_Duplicates");
        bip->value.copt = bip->value.copt == 'y' ? 'n' : 'y';

        bip->items_read = 1;
	ret = 1;
    }
    return ret;
}

void ReadImputed::do_menu2batch()
{
    Cstr Values[] = { "Imputed_Oxford_Single_Chr",
                      "Imputed_Info_Metric_Threshold", 
                      "Imputed_Hard_Call_Uncertainty", 
                      "Imputed_Genotype_Missing_Fraction", 
                      "Imputed_Missing_Codes",
                      "Imputed_Allow_Indels",
                      "Imputed_Allow_Duplicates"
    };

    for(size_t i = 0; i < ((sizeof Values) / sizeof (Cstr)); i++) {
        batch_item_type *bip = BatchItemGet(Values[i]);
        if (bip->items_read) 
            batchf(bip);
    }
}

void ReadImputed::do_batch2local()
{
    BatchValueGet(this->oxford_single_chr, "Imputed_Oxford_Single_Chr");
    BatchValueGet(this->info_threshold, "Imputed_Info_Metric_Threshold");
    BatchValueGet(this->hard_call_uncertainty, "Imputed_Hard_Call_Uncertainty");
    BatchValueGet(this->genotype_missing_fraction, "Imputed_Genotype_Missing_Fraction");
    BatchValueGet(this->allow_indels, "Imputed_Allow_Indels");
    BatchValueGet(this->allow_dups, "Imputed_Allow_Duplicates");
    BatchValueGet(this->info_file, "Input_Imputed_Info_File");

    Vecs vec;
    BatchValueGet(vec, "Imputed_Missing_Codes");
    Input->MissingCodesSet.clear();
    Input->MissingCodesSet.insert(vec.cbegin(), vec.cend());
}

void ReadImputed::show_settings()
{
    msgvf("\n");
    msgvf("Impute2 Analysis Oxford Single Chr:         %s\n", C(this->oxford_single_chr));
    msgvf("Impute2 Analysis Info Metric Threshold:     %.3f\n", this->info_threshold);
    msgvf("Impute2 Analysis Hard Call Uncertainty:     %.3f\n", this->hard_call_uncertainty);
    msgvf("Impute2 Analysis Genotype Missing Fraction: %.3f\n", this->genotype_missing_fraction);
    msgvf("Impute2 Analysis Allow Duplicate Markers:   %c\n",   this->allow_dups);
    msgvf("Impute2 Analysis Allow Indels:              %c\n",   this->allow_indels);

    Vecs vec;
    BatchValueGet(vec, "Imputed_Missing_Codes");
    join(vec, Input->MissingCodes, " ");
    msgvf("Impute2 Analysis Missing Trait Value Codes: %s\n", C(Input->MissingCodes));

    msgvf("Impute2 Analysis Imputed Sample File:       %s\n", C(this->sample_file));
    msgvf("Impute2 Analysis Imputed File:              %s\n",  this->impute_file);
    msgvf("Impute2 Analysis Imputed Info File:         %s\n", C(this->info_file));

}

void ReadImputed::do_init(Input_Base *inp)
{
    this->input = inp;
    this->files(*inp->input_files.bedfl, *inp->input_files.pedfl);
    show_settings();

    read_info = false;  // this should be default ... but not for win mvc
    if (!info_file.empty()) {
        ifstream infs;
        infs.open(info_file.c_str());
        if (! infs.is_open() ) {
            errorvf("read_imputed_file: Info file can not be opened: \"%s\"\n", C(info_file));
            EXIT(FILE_NOT_FOUND);
        }
        infs.close();
        read_info = true;
    }

    read_input_file();

    if (read_info)
        read_info_file();

    check_indelsNdups();

    read_sample_file();
}

linkage_locus_top *ReadImputed::do_names(const char *&names_fn)
{
    linkage_locus_top *LTop = build_impute2_names();
    names_fn = impute_file;

    return LTop;
}

void ReadImputed::do_map(std::vector<m2_map>& additional_maps)
{
//  m2_map 
    impute_map = m2_map("IMPUTE", 'p');

    build_impute2_map(impute_map);

    additional_maps.push_back(impute_map);
}

linkage_ped_top *ReadImputed::do_ped(linkage_locus_top *LTop)
{
    int num_peds = 0;

    annotated_ped_rec *persons = build_impute2_ped(LTop, &num_peds);

    build_genotypes(LTop, persons);

    linkage_ped_top *Top;
    Top = mk_ped_top(persons, this->people.size(), LTop, num_peds,
                     /*untyped*/ 0, /*totaltyped*/ this->people_filtered,
                     /*groups*/ NULL, 0, 0, 
                     /*num_err*/0, 1);

    return Top;

}

void ReadImputed::do_gc()
{
    for (int m = 0; m < markers_all; m++) {
        delete markers[m];
    }
    markers.clear();
    Vecmarkerp().swap(markers);
//  msgvf("read_impute markers %d %d\n", markers.size(), markers.capacity());

    for (int p = 0; p < people_all; p++) {
        people[p].clear();
        VecsDB().swap(people[p]);
//      msgvf("read_impute people[i] %d %d\n", people[p].size(), people[p].capacity());
    }
    people.clear();
    Vecvecs().swap(people);
//  msgvf("read_impute people %d %d\n", people.size(), people.capacity());

    impute_map.gc();
}

void ReadImputed::read_input_file ()
{
    int dbg = 0;

    ifstream ifs;
    ifs.open(impute_file);
    if (! ifs.is_open() ) {
        errorvf("read_imputed_file: Can not open \"%s\" file\n", impute_file);
        EXIT(FILE_NOT_FOUND);
    }

    Str hmm, rsid, pos, A, B;
    VecsDB fields;
    int line_n = 0;
    int nochr = 0;

    Str  name;
    Str  chrm;
    SECTION_ERR_INIT(bad_line);
    SECTION_ERR_INIT(halftyped);
    while (! ifs.eof() ) {

        if (Ncol == 0) {
            Str l;
            Vecs ls;
            getline(ifs, l);
            split(ls, l);
            fflush(stdout);
            Ncol = ( (ls.size() - 5) % 3 == 0) ? 5 : 6;
            if (Ncol == 5) {
                hmm  = ls[0];
                rsid = ls[1];
                pos  = ls[2];
                A    = ls[3];
                B    = ls[4];
            } else {
                hmm  = ls[0];
                chrm = ls[1];
                rsid = ls[2];
                pos  = ls[3];
                A    = ls[4];
                B    = ls[5];
            }
        } else if (Ncol == 5) {
            ifs >> hmm;
            ifs >> rsid;
            ifs >> pos;
            ifs >> A;
            ifs >> B;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else if (Ncol == 6) {
            ifs >> hmm;
            ifs >> chrm;
            ifs >> rsid;
            ifs >> pos;
            ifs >> A;
            ifs >> B;
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        if (ifs.eof()) break;
        line_n++;

        if (dbg) {
            cout << line_n << ": ";
            cout << hmm << " ";
        }
        if (hmm.compare(0, 3, "chr") == 0) {
            hmm.erase(0, 3);
        }
        if (inMap(hmm, input->G.chrm_set)) {
            chrm = hmm;
            fields.clear();
            split(fields, rsid, ":", 3);
            if (rsid == "." || rsid == "NA" || rsid == "na")
                name = "chr" + chrm + "_" + pos;
            else
                name = fields[0];  // it seems to be this way
        } else /* if (hmm == "---") */ {   // first column can be -9 ... anything
            fields.clear();
            split(fields, rsid, ":", 3);

            if (dbg) {
                cout << rsid << " ";
                cout << "#" << fields.size() << " ";
            }
            if (Ncol == 6) {
//              chrm = oxford_single_chr;  // already set
                name = "chr" + chrm + "_" + pos;
            } else if (rsid == "." || rsid == "NA" || rsid == "na") {
                chrm = oxford_single_chr;
                name = "chr" + chrm + "_" + pos;
            } else if (fields.size() > 0) {
                if (fields[0].compare(0, 3, "chr") == 0) {
                    fields[0].erase(0, 3);
                }
                if (inMap(fields[0], input->G.chrm_set)) {
                    chrm = fields[0];
                    name = "chr" + fields[0] + "_" + pos;
                } else if (fields[0].compare(0, 2, "rs") == 0 && oxford_single_chr != "--") {
                    chrm  =  oxford_single_chr;
                    name = fields[0];
                } else if (oxford_single_chr != "--") {
                    chrm  =  oxford_single_chr;
                    name = fields[0];
                } else {
                    SECTION_ERR(bad_line);
                    errorvf("impute2 file: bad line(%d): %s %s %s\nrs_id field first item unexpected: no chromosome was specified and can not be inferred\n", 
                           C(line_n), C(hmm), C(rsid), C(pos));
                    nochr++;
                    continue;
                }
                if (dbg) {
                    cout << fields[0] << " ";
                    cout << fields[1] << " ";
                    cout << fields[2] << " ";
                    cout << fields[3] << " ";
                }
                if (fields.size() > 1 && fields[1] != pos) {
                    SECTION_ERR(bad_line);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id pos field (%s) does not match position column (%s)\n",
                           C(line_n), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[1]), C(pos));
                }
                if (check_format && fields.size() > 2 && fields[2] != A) {
                    SECTION_ERR(bad_line);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id A allele field (%s) does not match A column (%s)\n",
                           C(line_n), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[2]), C(A));
                }
                if (check_format && fields.size() > 3 && fields[3] != B) {
                    SECTION_ERR(bad_line);
                    warnvf("impute2 file: bad line(%d): %s %s %s %s %s\n         rs_id B allele field (%s) does not match B column (%s)\n",
                           C(line_n), C(hmm), C(rsid), C(pos), C(A), C(B), C(fields[3]), C(B));
                }
                if (A != B && (A == "0" || B == "0")) {
                    SECTION_ERR(halftyped);
                    warnvf("impute2 file: odd line(%d): %s %s %s %s %s\n         heterozygote samples will appear halftyped because only one allele is 0\n",
                           C(line_n), C(hmm), C(rsid), C(pos), C(A), C(B));

                }
            }
            if (check_format && hmm != "---") {
                SECTION_ERR(bad_line);
                warnvf("impute file: bad line(%d): %s %s %s %s %s\n         first field is not --- or chromosome\n",
                       C(line_n), C(hmm), C(rsid), C(pos), C(A), C(B));
            }

        }
        if (dbg) {
            cout << pos << " ";
            cout << A << " ";
            cout << B << " ";
            cout << endl;
        }
        markers.push_back(new ImpMarker(name, chrm, pos, A, B, read_info));
    }

    ifs.close();

    SECTION_ERR_FINI(bad_line);
    SECTION_ERR_FINI(halftyped);
    if (nochr) {
        errorvf("read_imputed_file: can not determine chromosomes for %d markers\n",
                nochr);
        EXIT(DATA_INCONSISTENCY);
    }

    markers_all = markers.size();
    markers_filtered = markers_all;
}

void ReadImputed::check_indelsNdups()
{
// some checks: indels and dups
    ImpMarker *mp = NULL, *mpprev = NULL;
    int line_n = 0;
    int skip_count = 0;

    SECTION_ERR_INIT(indel);
    for (Vecmarkerpp mpp = markers.cbegin(); mpp != markers.cend(); mpp++) {
        mp     = *mpp;

        line_n++;

        if (mp->skip) continue;

        if ( ( mp->A.size() > 1 && mp->A.compare(0, 5, "dummy" )) || 
             ( mp->B.size() > 1 && mp->B.compare(0, 5, "dummy" )) ) {
            if (allow_indels == 'y') {
                SECTION_ERR(indel);
                warnvf("Marker: %s (bp %s) indel alleles[%s/%s] allowed [line %d].\n",
                       C(mp->name), C(mp->pos), C(mp->A), C(mp->B), line_n);
            } else {
                SECTION_ERR(indel);
                warnvf("Marker: %s (bp %s) indel alleles[%s/%s] ignored [line %d].\n",
                       C(mp->name), C(mp->pos), C(mp->A), C(mp->B), line_n);
                skip_count++;
                mp->skip = true;
            }
        }
    } 
    SECTION_ERR_FINI(indel);
    markers_filtered = markers_filtered - skip_count;
    if (allow_indels == 'y')
        mssgvf("Markers remaining %d; NO Filtering by indels\n", markers_filtered);
    else
        mssgvf("Markers remaining %d; Filtered by indels %d\n", markers_filtered, skip_count);

    int same_name = 0;
    int same_run  = 0;
    mp = NULL;
    line_n = 0;
    skip_count = 0;
    SECTION_ERR_INIT(duplicate);
    SECTION_ERR_INIT(duplicate_map);
    for (Vecmarkerpp mpp = markers.cbegin(); mpp != markers.cend(); mpp++) {
        mpprev = mp;
        mp     = *mpp;

        line_n++;

        if (mp->skip) continue;

        same_run = 0;

        if (mpprev == NULL)
            ;
        else if (mpprev->pos == mp->pos) {

            SECTION_ERR(duplicate);
            if (mp->name == mpprev->name) {
                same_name++;
                same_run++;
                warnvf("Markers: %s ", C(mp->name));
            } else {
                warnvf("Markers: %s,%s ", C(mpprev->name), C(mp->name));
            }
            if (allow_dups == 'y') {
                errvf("(bp %s) repeated with different alleles [%s/%s %s/%s] both allowed [line %d].\n",
                       C(mp->pos), C(mpprev->A), C(mpprev->B), C(mp->A), C(mp->B), line_n);
                if (same_run) {
                    char N[10];
                    sprintf(N, "@%d", same_run+1);

                    SECTION_ERR(duplicate_map);
                    warnvf("Remapping %s to ", C(mp->name));
                    mp->name.append(N);
                    errvf("%s\n", C(mp->name));
                }
            } else {
                errvf("(bp %s) repeated with different alleles [%s/%s %s/%s] ignored repeat [line %d].\n",
                       C(mp->pos), C(mpprev->A), C(mpprev->B), C(mp->A), C(mp->B), line_n);
                skip_count++;
                mp->skip = true;
            }
        }
    }
    SECTION_ERR_FINI(duplicate);
    SECTION_ERR_FINI(duplicate_map);
    markers_filtered = markers_filtered - skip_count;

    if (allow_dups == 'y')
        mssgvf("Markers remaining %d; NO Filtering by duplicate markers\n", markers_filtered);
    else
        mssgvf("Markers remaining %d; Filtered by duplicate markers %d\n",
               markers_filtered, skip_count);
    if (same_name)
        mssgvf("\t%d duplicate markers had same name.\n\n", same_name);
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
        EXIT(FILE_NOT_FOUND);
    }

    Str line;
    getline(ifs, line);
    Vecs l, h;
    split(l, line);
    split(h, info_file_hdr);
//  if (l[0] != h[0] || l[1] != h[1] || l[2] != h[2] || l[3] != h[3] || 
//      l[4] != h[4] || l[5] != h[5] || l[6] != h[6] || l[7] != h[7] || 
//      l[8] != h[8])
//  if (line.compare(0, info_file_hdr.size(), info_file_hdr) != 0)
    if (l[0] != h[0] || l[1] != h[1] || l[2] != h[2])
    {
        errorvf("Bad header for \"%s\" file.\n", C(info_file));
        errorvf(" expecting: %s %s %s\n", C(h[0]), C(h[1]), C(h[2]));
        errorvf(" found: %s %s %s\n", C(l[0]), C(l[1]), C(l[2]));
        EXIT(INPUT_DATA_ERROR);
    }

    Mapsi col2idx;
    int i, info_col = 6, certainty_col = 7;
    i = 0;
    for (VecspDB fp = l.cbegin(); fp != l.cend(); i++, fp++)
        col2idx[*fp] = i;
    if (!map_get(col2idx, "info", info_col)) {
        errorvf("Bad header for \"%s\" file.\n", C(info_file));
        errorvf("Can not find \"info\" column in %s.\n", C(line));
        EXIT(INPUT_DATA_ERROR);
    }
    if (!map_get(col2idx, "certainty", certainty_col)) {
        errorvf("Bad header for \"%s\" file.\n", C(info_file));
        errorvf("Can not find \"certainty\" column in %s.\n", C(line));
        EXIT(INPUT_DATA_ERROR);
    }

    Vecs        fields;
    Vecmarkerpp mpp;
    size_t line_n = 0;
    int    skip_count = 0;
    double info;

    SECTION_ERR_INIT(info_threshold);
    ImpMarker *mp;
    for (mpp = markers.cbegin(); ! ifs.eof(); mpp++) {
        mp = *mpp;
        getline(ifs, line);
        if (ifs.fail()) break;
        line_n++;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        fields.clear();
        split(fields, line);
        if (fields[2] != mp->pos) {
            errorvf("Files \"%s\" and \"%s\" do not list markers in the identical order starting at line %d: %s.\n", 
                    impute_file, C(info_file), line_n, C(line));
                EXIT(DATA_INCONSISTENCY);
        }

        info = atof(fields[info_col].c_str());
        mp->info = info;
        mp->certainty = atof(fields[certainty_col].c_str());

        if (! mp->skip) {
            if (info < info_threshold) {
                skip_count++;
                mp->skip = true;

                SECTION_ERR(info_threshold);
                warnvf("Marker: %s %s %s %s %s info (%.4f) < threshold (%.4f) [line %d]\n", 
                       C(fields[0]), C(fields[1]), C(fields[2]), C(fields[3]), C(fields[4]),
                       info, info_threshold, line_n);
            }
        }
    }

    markers_filtered = markers_filtered - skip_count;
    SECTION_ERR_FINI(info_threshold);

    if (mpp != markers.cend() || line_n != markers.size()) {
        errorvf("Files \"%s\" and \"%s\" are different lengths: %d vs %d\n",
                impute_file, C(info_file), markers.size(), line_n);
        EXIT(DATA_INCONSISTENCY);
    }
    ifs.close();

    mssgvf("Markers remaining %d; Filtered by info threshold %d\n", markers_filtered, skip_count);
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
        EXIT(FILE_NOT_FOUND);
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
        EXIT(INPUT_DATA_ERROR);
    }
    split(sample_file_hdr1a, line);
    split(fields, line);
//  column_type.reserve(fields.size());

    Mapsi  col2idx;
    Mapii  idx2rsvd;
    int    i;
    int    ret;
    i = 0;
    for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
        col2idx[*fp]    = i;
//      column_type.push_back(PHENO);
    }

    Cstr rsvd[] = {"ID_1", "ID_2", "father", "mother", "sex", "missing"};
    for (i = 0; i < 6; i++) {
        if (map_get(col2idx, rsvd[i], ret))
            idx2rsvd[ret] = i;
    }

//line 2
    sample_file_hdr1b.clear();
    Str missin;
    Cstr fill1[] = {"", "", "0", "0", "0"};
    sample_file_hdr1b.insert(sample_file_hdr1b.end(), fill1, fill1+5);
    i = 0;
    for (VecspDB fp = sample_file_hdr1a.cbegin(); fp != sample_file_hdr1a.cend(); i++, fp++) {
        if (map_get(idx2rsvd, i, ret)) {
            if (ret < 5)
                sample_file_hdr1b[ret] = *fp;
            else
                missin = *fp;
        } else
            sample_file_hdr1b.push_back(*fp);
    }
    sample_file_hdr1b.push_back(missin);

    getline(ifs, line);
    if (ifs.fail()) {
        errorvf("Bad second header line for \"%s\" file.\n", C(sample_file));
        EXIT(INPUT_DATA_ERROR);
    }
    if (dbg) {
        cout << line_n << ": " << line << endl;
    }
    split(sample_file_hdr2a, line);
    Cstr fill2[] = {"", "", "", "", ""};
    sample_file_hdr2b.insert(sample_file_hdr2b.end(), fill2, fill2+5);
    i = 0;
    for (VecspDB fp = sample_file_hdr2a.cbegin(); fp != sample_file_hdr2a.cend(); i++, fp++) {
        if (map_get(idx2rsvd, i, ret)) {
            if (ret < 5)
                sample_file_hdr2b[ret] = *fp;
            else
                missin = *fp;
        } else
            sample_file_hdr2b.push_back(*fp);
    }
    sample_file_hdr2b.push_back(missin);


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
        ret = 0;
        for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
            if (map_get(idx2rsvd, i, ret)) {
                if (ret < 5)
                    mappedfields[ret] = *fp;
                else
                    missin = *fp;
            } else {
                mappedfields.push_back(*fp);
            }
        }
        mappedfields.push_back(missin);
        people.push_back(mappedfields);
    }

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
        types[i] = ( ((*typep) == "B") || ((*typep) == "D") ) ? 'A' : 'T';
    }

    ImpMarker *mp;
    for (Vecmarkerpp marp = markers.cbegin(); marp != markers.cend(); marp++) {
        mp = *marp;
        if (! mp->skip ) {
            if (i >= num_all) {
                errorvf("Internal Error: build_impute2_names() count of skipped markers too large\n");
                EXIT(DATA_INCONSISTENCY);
            }
            names[i] = CALLOC(mp->name.size()+1, char);
            strcpy(names[i], mp->name.c_str());
            types[i] = 'M';
            i++;
        }
    }

    return read_common_marker_data(num_pheno + num_markers, num_markers, names, types,
                                  /*annotated*/ 1, /*penetrances_read*/ 0, 0, NULL);
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
                EXIT(DATA_INCONSISTENCY);
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
//              if (LTop->Locus[mrk_num].Type != XLINKED) {
//                  LTop->Locus[mrk_num].Type = XLINKED;
//              }
            } else if (chr == MALE_CHROMOSOME) {
                human_y++;
//              if (LTop->Locus[mrk_num].Type != YLINKED) {
//                  LTop->Locus[mrk_num].Type = YLINKED;
//              }
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
//  asm("int $3");
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

    SECTION_ERR_EXTERN(illegal_affect);
    SECTION_ERR_EXTERN(illegal_quant);
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

// 1 ID_2 -> ID
        strcpy(entry->ID, C(pp[1]));
        PLINK.individuals += 1;
        entry->per_index = curr_per_index;
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
            entry->Sex = '1';
            PLINK.males += 1;
        } else if (tolower((unsigned char)entry->Sex) == 'f' || entry->Sex == '2') {
            entry->Sex = '2';
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
	    switch (LTop->Locus[i].Type) {
	    case AFFECTION:
                ret = plink_annot_string_aff_phen(entry->rec_num+HDR, &(LTop->Pheno[i]),
                                                  &entry->pheno[i], C(*vp));
                num_errors += ret;

                if (entry->pheno[i].Affection.Status == UNDEF) {
                    SECTION_ERR(illegal_affect);
                    errorvf("File %s, Entry %d : Invalid status at locus %s\n",
                            C(sample_file), entry->rec_num+HDR, LTop->Locus[i].Name);
                    num_errors++;
                }
                if (LTop->Pheno[i].Props.Affection.ClassCnt != 1) {
                    if (entry->pheno[i].Affection.Class == UNDEF) {
                        SECTION_ERR(illegal_affect);
                        errorvf("File \"%s\", Entry %d : Invalid liability class at locus %s\n\n",
                                C(sample_file), entry->rec_num+HDR, LTop->Locus[i].Name);
                        num_errors++;
                    }
                }
                break;
	    case QUANT:
	    default:
                ret = plink_annot_string_quant_phen(entry->rec_num+HDR, &(LTop->Pheno[i]),
                                                    &entry->pheno[i], C(*vp));
                num_errors += ret;
                break;
            }
        }
    }
    SECTION_ERR_FINI(illegal_affect);
    SECTION_ERR_FINI(illegal_quant);
    printf("number phenotype errors: %d\n", num_errors);

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

//need to share ifs and token
void ReadImputedGenotypeReadHelper::genotypes_init()
{
    ifs.open(C(impute_file));
    if (! ifs.is_open() ) {
        errorvf("read_imputed_genotype_file: Can not open \"%s\" file\n",
                C(impute_file));
        EXIT(FILE_NOT_FOUND);
    }
    line_n = 0;
}

boolean ReadImputedGenotypeReadHelper::genotypes_marker_hdr(int mrk_idx, string& hmm, string& chrm, string& rsid,
                                                            string& pos, const char* &A, const char* &B)
{
    getline(ifs, line);
    line_n++;
//?    tod_im_line("impute read line  ");
    if (ifs.eof()) return false;
    token.set(line);
//?    tod_im_line("impute set line   ");

//	ifs >> hmm; // +A+B+ nums() all read via >> @ 14.62 sec
//                               vector<double> nums  12.26
//                               double nums[3]       11.59
    if (rip->Ncol == 5) {
        token.more(hmm);
        token.more(rsid);
        token.more(pos);
        token.more(A);
        token.more(B);
    } else {
        token.more(hmm);
        token.more(chrm);
        token.more(rsid);
        token.more(pos);
        token.more(A);
        token.more(B);
    }

    return true;
}

void ReadImputedGenotypeReadHelper::genotypes_sample_prob(Token::d3& nums)
{
    token.getDC(nums, 3);
}

void ReadImputedGenotypeReadHelper::genotypes_end()
{
    ifs.close();
}

void ReadImputed::build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons)
{
    ReadImputedGenotypeReadHelper  gh;
    gh.impute_file = impute_file;
    gh.rip = this;

    build_internal_genotypes(LTop, persons, gh);
}

void ReadImputed::build_internal_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons,
                                           GenotypeReadHelper &gh)
{
    int dbg = 0;

    string hmm, chrm, rsid, pos;
    const char *A, *B;
    Token::d3 nums;

    gh.genotypes_init();
    ImpMarker   *mp;

    int   mrk_idx = sample_file_hdr2b.size() - 5 -1;
    linkage_locus_rec *locus;

    char *callele1;
    char *callele2;
    char *callele0  = canonical_allele(C("0"));

    Tod tod_gen("impute genotypes");
    Tod tod_per(30);
    Tod tod_im_line(30);
    Tod tod_im_line_cpy(30);
    SECTION_ERR_INIT(genotype_missing_fraction);
    SECTION_ERR(genotype_missing_fraction);
    warnvf("%8s %8s %8s %8s  %s\n         %8s %8s %8s %8s  %s\n         %8s %8s %8s %8s  %s\n",
           "untyped", "uncertain", "good", "geno-", "Marker + chr:pos",
           "marker", "hard", "hard", "typing", "",
           "0/0", "call", "call", "rate", "");
    for (Vecmarkerpp mpp = markers.begin(); mpp != markers.end(); mpp++) {
        tod_im_line.reset();
        tod_im_line_cpy.reset();

        gh.genotypes_marker_hdr(mrk_idx, hmm, chrm, rsid, pos, A, B);
        tod_im_line("impute read line  ");

        mp = *mpp;

        if (pos != mp->pos) {
            errorvf("internal error: impute_bgen_file (\"%s\") second pass does not match first pass at marker #%d %s\n",
                    impute_file, mrk_idx, C(rsid));
                EXIT(DATA_INCONSISTENCY);
        }
        if (mp->skip) {
            gh.genotypes_skip();
            continue;
        }
        mrk_idx++;  // hence the -1 above
        locus = &LTop->Locus[mrk_idx];

        callele1    = canonical_allele(A);
        callele2    = canonical_allele(B);

        int i;
        int sam = 0;
        int zero = 0, uncertain = 0, good = 0;

        tod_per.reset();
        annotated_ped_rec *entry = persons;
//      if (mrk_idx == 237)
//          warnvf("Marker 237 == %s, scale = \n", locus->Name);
	for(int p = 0 ; p < people_filtered; p++, entry++) {
            sam++;

            gh.genotypes_sample_prob(nums);
//          if (mrk_idx == 237)
//              warnvf("sam %d, %.4f %.4f %.4f\n", sam, nums[0], nums[1], nums[2]);
            if (nums[0] > nums[1]) {
                i = 0;
            } else {
                i = 1;
            }
            if (nums[2] > nums[i]) {
                i = 2;
            }

            if (nums[i] == 0) {
                i = 3;
                zero++;
            } else if (nums[i] <= 1 - hard_call_uncertainty) {
                i = 3;
                uncertain++;
            } else {
                good++;
            }

            if (dbg) {
                cout << "#";
//              cout << line_n << "@";
                cout << sam << ": ";
                cout << mp->name << " ";
                cout << nums[0] << " " << nums[1] << " " << nums[2] << "  " << i << "\n";
            }

            if (i == 0)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele1);
            else if (i == 1)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele2);
            else if (i == 2)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele2, callele2);
            else if (i == 3)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele0, callele0);
        }
        tod_per("impute person loop");
        if (good < genotype_missing_fraction * people_filtered) {
/* some day
??          mp->skip = true;
            LTop->Marker[mrk_idx].chromosome = MISSING_CHROMO;
            locus->??
            skip_count++;
*/
            SECTION_ERR(genotype_missing_fraction);
            warnvf("%8d %8d %8d %8.3f  %s %s:%s\n",
                   zero, uncertain, good, ((double)good)/people_filtered,
                   C(mp->name), C(mp->chr), C(mp->pos));
        }
    }
    tod_gen();
    SECTION_ERR_FINI(genotype_missing_fraction);

    gh.genotypes_end();
}
