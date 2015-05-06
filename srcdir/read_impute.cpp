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
#include "error_messages_ext.h"
extern void           Exit(int arg, const char *file, const int line, const char *err);

#include "str_utils.hh"
#include "read_impute.hh"

// g++ does not like these const operators on "vectordb" classes
#define cbegin() begin()
#define cend()   end()

using namespace std;

int dbg = 1;

void pr_str(const string& x) {
    cout << x << " ";
}

void ReadImputed::read_imputed_file ()
{

    ifstream ifs;
    ifstream infs;
    int dbg = 0;

//  asm("int $3");
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
        if (inMap(hmm, G.chrm_set)) {
            chrm = hmm;
            name = rsid;  // it seems to be this way
        } else if (hmm == "---") {
            fields.clear();
            split(fields, rsid, ":", 3);

            if (dbg) {
                cout << rsid << " ";
                cout << "#" << fields.size() << " ";
            }
            if (fields.size() > 1) {
                if (inMap(fields[0], G.chrm_set)) {
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
        markers.push_back(new Marker(name, chrm, pos, A, B, read_info));
    }
    SUPPRESS_MSSG_NESTED_FINI(bad_line_msg);

}

Str ReadImputed::info_file_hdr = "snp_id rs_id position a0 a1 exp_freq_a1 info certainty type";

void ReadImputed::read_info_file () {
    ifstream ifs;
    int dbg = 0;

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
    Vecmarkerpp mp;
    int    line_n = 1;
    int    skip_count = 0;
    double info;

    SUPPRESS_MSSG_NESTED_INIT(info_threshold_msg);
    for (mp = markers.cbegin(); ! ifs.eof(); mp++) {
        getline(ifs, line);
        if (ifs.fail()) break;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        line_n++;
        fields.clear();
        split(fields, line);
        if (fields[2] != (*mp)->pos) {
            errorvf("Files \"%s\" and \"%s\" do not list markers in the identical order starting at line %d: %s.\n", 
                    impute_file, C(info_file), line_n, C(line));
                EXIT(1);
        }
        info = atof(fields[6].c_str());
        (*mp)->info = info;
        (*mp)->certainty = atof(fields[7].c_str());

        if (info < info_threshold) {
            skip_count++;
            (*mp)->skip = true;

            SUPPRESS_MSSG_NESTED(info_threshold_msg);
            warnvf("Marker: %s %s %s %s %s info (%.4f) < threshold (%.4f)\n", 
                   C(fields[0]), C(fields[1]), C(fields[2]), C(fields[3]), C(fields[4]),
                   info, info_threshold);
        }
            
    }
    SUPPRESS_MSSG_NESTED_FINI(info_threshold_msg);

    if (mp != markers.cend() || line_n != markers.size() + 1 /*hdr*/) {
        errorvf("Files \"%s\" and \"%s\" are different lengths: %d vs %d\n",
                impute_file, C(info_file), markers.size(), line_n);
        EXIT(1);
    }
    ifs.close();

    warnvf("%d of %d markers that are below the info_metric_threshold will be skipped.\n", 
           skip_count, markers.size());
}

Str ReadImputed::sample_file_hdr = "ID_1 ID_2 missing";

void ReadImputed::read_sample_file () {
    ifstream ifs;
    int dbg = 0;

    asm("int $3");
    ifs.open(sample_file.c_str());
    if (! ifs.is_open() ) {
        errorvf("read_sample_file: Can not open \"%s\" file\n", C(sample_file));
        EXIT(1);
    }

    Str  line;
    int  line_n = 2;
    VecsDB fields;
    Vecs mappedfields;

//line 1
    getline(ifs, line);
    if (line.compare(0, sample_file_hdr.size(), sample_file_hdr) != 0) {
        errorvf("Bad header for \"%s\" file.\n", C(sample_file));
        errorvf(" expecting: %s\n", C(sample_file_hdr));
        errorvf(" found: %s\n", C(line));
        EXIT(1);
    }
    split(sample_file_hdr1a, line);
    split(fields, line);

    Mapsi  col2idx;
    Mapsip col2idxp;
    Mapii  idx2fixed;
    Mapiip idx2fixedp;
    int    i;
    int    ret;
    i = 0;
    for (VecspDB fp = fields.cbegin(); fp != fields.cend(); i++, fp++) {
        col2idx[*fp] = i;
    }

    Cstr special[] = {"ID_1", "ID_2", "father", "mother", "sex"};
    for (i = 0; i < 5; i++) {
        if (map_get(col2idx, special[i], ret))
            idx2fixed[ret] = i;
    }

//line 2
    sample_file_hdr1b.clear();
    Cstr fill1[] = {"", "", "0", "0", "0"};
    sample_file_hdr1b.insert(sample_file_hdr1b.end(), fill1, fill1+5);
    i = 0;
    for (VecspDB fp = sample_file_hdr1a.cbegin(); fp != sample_file_hdr1a.cend(); i++, fp++) {
        if (map_get(idx2fixed, i, ret))
            sample_file_hdr1b[ret] = *fp;
        else
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
        if (map_get(idx2fixed, i, ret))
            sample_file_hdr2b[ret] = *fp;
        else
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
                mappedfields[ret] = *fp;
            } else {
                mappedfields.push_back(*fp);
            }
        }
        people.push_back(mappedfields);
    }
}

void ReadImputed::read_genotypes_file () {

    ifstream ifs;
    string line;

    string hmm, rsid, pos, A, B;
    vector<double> nums;

    string field;
    int idx = 0;

    asm("int $3");
    ifs.open(impute_file);
    if (! ifs.is_open() ) {
        errorvf("read_imputed_genotype_file: Can not open \"%s\" file\n",
                impute_file);
        EXIT(1);
    }
    Token token(3);
    Vecs  triple_genotype;
    Str   genotype[3];
    Vecmarkerpp mpp;
    Marker     *mp;
    SUPPRESS_MSSG_NESTED_INIT(skip_msg);

    for (mpp = markers.cbegin(); ! ifs.eof(); mpp++) {
        getline(ifs, line);
        if (ifs.eof()) break;
        idx++;

        token.set(line);

        token.more(hmm);
        token.more(rsid);
        token.more(pos);
        mp = *mpp;
        if (pos != mp->pos) {
            errorvf("internal error: impute_file (\"%s\") second pass does not match first pass at line %d\n",
                    impute_file, idx);
                EXIT(1);
        }
        if (mp->skip) {
            SUPPRESS_MSSG_NESTED(skip_msg);
            warnvf("Marker: skipped %s %s %s info (%.4f) < threshold (%.4f)\n", 
                   C(mp->name), C(mp->chr), C(mp->pos), mp->info, info_threshold);
        }
        token.more(A);
        token.more(B);
        genotype[0] = A + A;
        genotype[1] = A + B;
        genotype[2] = B + B;

        triple_genotype.clear();
        token.getD(nums);

        int i;
        int sam = 1;
        int triple = 0;
        double maxx;
        do {
            sam++;
            token.getD(nums);
            if (nums[0] > 0.0 && nums[1] > 0.0 && nums[2] > 0.0) {
                triple++;
                cout << "#";
                cout << idx << ": ";
                cout << mp->name << " ";
                cout << nums[0] << " " << nums[1] << " " << nums[2] << "\n";
                if (nums[0] > nums[1]) {
                    maxx = nums[0];
                    i = 0;
                } else {
                    maxx = nums[1];
                    i = 1;
                }
                if (nums[2] > maxx) {
                    maxx = nums[2];
                    i = 2;
                }
                triple_genotype.push_back(genotype[i]);
            }
        } while (nums.size() != 0);

        if (triple) {
            cout << "#";
            cout << idx << ": ";
            cout << hmm << " ";

            cout << rsid << " (";
            vector<string> fields;
            split(fields, rsid, ":");

            fields[0] = mp->name;
            cout << fields[0] << " ";
            cout << fields[1] << " ";
            cout << fields[2] << " ";
            cout << fields[3] << ") ";

            cout << pos << " ";

            cout << A << " ";
            cout << B << " ";

            cout << nums[0] << " ";
            cout << nums[1] << " ";
            cout << nums[2] << " ";

            for_each (triple_genotype.begin(), triple_genotype.end(), pr_str);

            cout << "sample: " << sam;
            cout << " triple: " << triple;

            cout << endl;
        }
    }
    SUPPRESS_MSSG_NESTED_FINI(skip_msg);
}
