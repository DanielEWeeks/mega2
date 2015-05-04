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
                    warnvf("impute file: bad line(%d) rs_id field first item: %s %s %s\n", 
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
                    warnvf("impute2 file: bad line(%d) rs_id field pos: %s %s %s %s %s\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B));
                    }
                if (fields[2] != A) {
                    warnvf("impute2 file: bad line(%d) rs_id field A: %s %s %s %s %s\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B));
                    }
                if (fields[3] != B) {
                    warnvf("impute2 file: bad line(%d) rs_id field B: %s %s %s %s %s\n",
                           C(idx), C(hmm), C(rsid), C(pos), C(A), C(B));
                    }
            }
        } else {
            warnvf("impute file: bad line(%d) first field not --- or chromosome: %s %s %s %s %s\n",
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

    Vecs fields;
    int  line_n = 1;
    Vecmarkerpp mp;

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
        (*mp)->info = atof(fields[6].c_str());
        (*mp)->certainty = atof(fields[7].c_str());
    }
    if (mp != markers.cend() || line_n != markers.size() + 1 /*hdr*/) {
        errorvf("Files \"%s\" and \"%s\" are different lengths: %d vs %d\n",
                impute_file, C(info_file), markers.size(), line_n);
        EXIT(1);
    }
    ifs.close();
}

void ReadImputed::read_imputed_file_genotype (const char *imp) {

    ifstream ifs;
    string line;

    string hmm, rsid, pos, A, B;
    vector<double> nums;

    string field;
    int idx = 0;

    asm("int $3");
    ifs.open(imp);
    if (! ifs.is_open() ) {
        cout << "Can not open file: " << "ps.20.impute" << endl;
        errorvf("read_imputed_genotype_file: Can not open \"%s\" file\n", imp);
        EXIT(1);
    }
    Token token(3);
    vector<string> triple_genotype;
    string genotype[3];

    while (! ifs.eof() ) {
        getline(ifs, line);
        if (ifs.eof()) break;
        idx++;

        token.set(line);

        token.more(hmm);
        token.more(rsid);
        token.more(pos);
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
            cout << idx << ": ";

            cout << hmm << " ";

            cout << rsid << " ";
            vector<string> fields;
            split(fields, rsid, ":");

            if (G.chrm_set.find(fields[0]) != G.chrm_set.end())
            {
                fields[0] = "chr" + fields[0] + "_" + fields[1];
            } else {
                // pass
            }
            cout << fields[0] << " ";
            cout << fields[1] << " ";
            cout << fields[2] << " ";
            cout << fields[3] << " ";

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
}

Str ReadImputed::sample_file_hdr = "ID_1 ID_2 missing";

void ReadImputed::read_sample_file () {
    ifstream ifs;
    int dbg = 0;

    ifs.open(sample_file.c_str());
    if (! ifs.is_open() ) {
        errorvf("read_sample_file: Can not open \"%s\" file\n", C(sample_file));
        EXIT(1);
    }

    Str line;
    getline(ifs, line);
    if (line.compare(0, sample_file_hdr.size(), sample_file_hdr) != 0) {
        errorvf("Bad header for \"%s\" file.\n", C(sample_file));
        errorvf(" expecting: %s\n", C(sample_file_hdr));
        errorvf(" found: %s\n", C(line));
        EXIT(1);
    }

    Vecs fields;
    int  line_n = 1;
    Vecmarkerpp mp;

    for (mp = markers.cbegin(); ! ifs.eof(); mp++) {
        getline(ifs, line);
        if (ifs.fail()) break;
        if (dbg) {
            cout << line_n << ": " << line << endl;
        }
        line_n++;
        fields.clear();
        split(fields, line);
        
    }
}
