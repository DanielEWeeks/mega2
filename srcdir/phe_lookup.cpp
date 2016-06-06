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
#include <stdlib.h>

#include "common.h"
#include "typedefs.h"
#include "utils_ext.h"
#include "genetic_utils_ext.h"
#include "error_messages_ext.h"

#include "plink_ext.h"
#include "phe_lookup_ext.h"

#include <vector>

#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

class simpleobj
{
public:
    simpleobj();

//    simpleobj(simpleobj& cpy);
    simpleobj& operator=(const simpleobj& cpy);
};

simpleobj::simpleobj() {}
//simpleobj::simpleobj(simpleobj& cpy) {}
simpleobj& simpleobj::operator=(const simpleobj& cpy) {
    if (this != &cpy) *this = cpy;
    return *this;
}



class pedper
{
public:
    pedper(char *ped, char *per);
    pedper(string ped, string per);
   ~pedper();

    struct LTpedper
    {
        bool operator()(const pedper s1, const pedper s2) const {
            bool ans;
            ans = (s1.ped > s2.ped) ? false : (s1.ped < s2.ped) || s1.per < s2.per;
            return ans;
        }
    };

//private:
    string ped;
    string per;
};

pedper::pedper(char *ped, char *per) {
    this->ped = ped;
    this->per = per;
}

pedper::pedper(string ped, string per) {
    this->ped = ped;
    this->per = per;
}

pedper::~pedper() {
}



class trait
{
public:
    trait(string name);
   ~trait();
//private:
    string name;
    string type;
    long   quant;
    long   values[4];

};

trait::trait(string name)
{
    this->name      = name;
    this->quant     = 0;
    this->values[0] = 0;
    this->values[1] = 0;
    this->values[2] = 0;
    this->values[3] = 0;
}

trait::~trait()
{
}



class phe
{
public:
    phe();
    phe(const char *path);
   ~phe();
    int  make();
    void hdr(char **names, int *types);
    void debug_dump_pydict();
    void test(const string path);
    int  search(const char *ped, const char *per, vector<string> &v);
    
    void free();
private:
    string path;
    vector<trait> traits;
    int debug;
    typedef map<pedper,vector<string>,pedper::LTpedper>  dict;
    typedef map<pedper,vector<string>,pedper::LTpedper>::iterator dictp;
    dict PHES;
};

phe::phe() {
}

phe::phe(const char *path) {
    this->path = string(path);
    PHES = dict();
    debug = 0;
}

phe::~phe() {
}


//
// For mapping the VCF Tools sample into the pedigree,person found in this file.
sample_map_type *SAMPLEIDS = (sample_map_type *)NULL;

int phe::make()
{
    string ped, per, val;
    string rest;
    const char *c_ptr;
    char  *endptr;
    double quant;
    int i = 0, sampleid_column_number = 0;
    map<string,int> samples;

    SAMPLEIDS = (sample_map_type *)NULL;

//    printf("phe::make\n"); fflush(stdout);
    ifstream ifs(path.c_str(), ifstream::in);

    if (true) {
        ifs >> ped;
        ifs >> per;
        if (ped != "FID" && per != "IID") {
            printf("YOU MUST PROVIDE A HEADER WITH THE PHENOTYPE FILE.\n");
            printf("And the first two IDs must be FID and IID.\n");
            ifs.close();
            return -1;
        }
        getline(ifs, rest);
        if (debug)
            printf("ped %s, per %s, rest %s\n", ped.c_str(), per.c_str(), rest.c_str());

        istringstream fields(rest, istringstream::in);

        while (fields.good()) {
            fields >> val;
	    // The comlumn name "SAMPLEID" is used to map the sample id into the <FID, IID>
	    if (val == "SAMPLEID") {
	      SAMPLEIDS = new sample_map_type();
	      sampleid_column_number = i;
	    } else {
	      traits.push_back(trait(val));
	    }
	    i++;
        }
    }

    while (true) {
        ifs >> ped;
        if (!ifs.good()) break;
        ifs >> per;
        getline(ifs, rest);
        if (debug)
            printf("ped %s, per %s, rest %s\n", ped.c_str(), per.c_str(), rest.c_str());

        vector<string> &v = PHES[pedper(ped, per)] = vector<string>();
        vector<trait>::iterator tp = traits.begin();

        istringstream fields(rest, istringstream::in);
        i = 0;
        while (fields.good()) {
            fields >> val;
            if (fields.fail()) continue;

	    // If we found the "SAMPLEID" column in the header file, and we are processing the
	    // SAMPLEID data column, then make a mapping between the sample name (val) and the <ped,per>
	    // Note that the user can choose not to fill in this column or fill it with the string "NA".
            // In this case there is no mapping between the sample ID and the person.
	    if (SAMPLEIDS != NULL && sampleid_column_number == i) {
              if (val != "NA") {
	        SAMPLEIDS->insert(pair<string,pair<string,string> >(val, pair<string,string>(ped, per)));
                samples[val]++; // count the number of uses of this string
	        i++;
              }
	      continue;
	    } else {
	      v.push_back(val);
	    }

	    if (val == "NA" || val == "0")
	      tp->values[0]++;
	    else if (val == "1")
	      tp->values[1]++;
	    else if (val == "2")
	      tp->values[2]++;
	    else if (PLINK.missing_pheno) {
	      c_ptr = val.c_str();
	      quant = strtod(c_ptr, &endptr);
	      if (endptr == c_ptr || (*endptr)) { /* conversion failed */
                tp->values[3]++;
                tp->quant++;
	      } else if (quant == PLINK.pheno_value)
                tp->values[0]++;
	      else {
                tp->values[3]++;
                tp->quant++;
	      }
	    } else {
	      tp->values[3]++;
	      tp->quant++;
	    }
	    tp++;
	    i++;
        } // while (fields.good()) {

//      printf("%ld# 1: %s; 2: %s; 3: %s; 4: %s; 5: %s\n", v->size(), (*v)[0].c_str(), (*v)[1].c_str(), (*v)[2].c_str(), (*v)[3].c_str(), (*v)[4].c_str());

    } // while (true) {
    ifs.close();

    // Check for multiple occurranaces of a SAMPLEID string and give an error if found.
    for (map<string,int>::iterator it=samples.begin(); it!=samples.end(); ++it)
        if (it->second > 1) {
            errorvf("The phenotype file SAMPLEID column entry '%s' occurs more than once.\n", it->first.c_str());
            EXIT(DATA_TYPE_ERROR);
        }

    return (int)traits.size();
}

#ifdef _WIN
#define strdup(str) _strdup(str)
#endif
void phe::hdr(char **names, int *types)
{
    char **p = names;
    int *t   = types;

    vector<trait>::iterator nn = traits.begin(), ee = traits.end();
    for (; nn != ee; nn++) {
        *p++ = strdup(nn->name.c_str());
        *t++ = nn->quant ? 1 : 0;
        if (debug) printf("name %s: quant %d.\n", p[-1], t[-1]);
    }
}

int phe::search(const char *ped, const char *per, vector<string> &v)
{
    int success = 0;

    dictp ans = PHES.find(pedper(ped, per));

    if (ans != PHES.end()) {
        v = ans->second;
        if (debug) {
            size_t lim = v.size();
            printf(
#ifdef MS_PRINTF
                "ped %s, per %s; %Iu# %s %s %s %s %s\n",
#else
                "ped %s, per %s; %zu# %s %s %s %s %s\n",
#endif
                ped, per, lim, 
                lim > 0 ? v[0].c_str() : "",
                lim > 1 ? v[1].c_str() : "",
                lim > 2 ? v[2].c_str() : "",
                lim > 3 ? v[3].c_str() : "",
                lim > 4 ? v[4].c_str() : "");
        }
        success = 1;
    } else {
        fprintf(stderr, "WARNING: pedigree %s, person %s; phenotype lookup failed\n", ped, per);
        success = 0;
    }

    return success;
}

void phe::test(const string filepath)
{
    vector<string> v;
    string rest;
    string ped, per, val;
    int i;

//    printf("phe::test\n"); fflush(stdout);
    ifstream ifs(filepath.c_str(), ifstream::in);
    while (true) {
        ifs >> ped;
        if (!ifs.good()) break;
        ifs >> per;
        getline(ifs, rest);

        if (!search(ped.c_str(), per.c_str(), v))
            continue;

        istringstream fields(rest, istringstream::in);
        i = 0;
        while (fields.good()) {
            fields >> val;
            printf("%d: s/f %s/%s; ", i, v[i].c_str(), val.c_str());
            if (v[i] != val)  printf("***");
            i++;
        }
        printf("\n"); fflush(stdout);

    }
    ifs.close();
}

void phe::debug_dump_pydict()
{
//    printf("phe::dump\n"); fflush(stdout);
    vector<trait>::iterator nn = traits.begin(), ee = traits.end();
    for (; nn != ee; nn++) {
        printf("name %s; v[0] %ld, v[1] %ld, v[2] %ld, v[3] %ld\n",
               nn->name.c_str(),
               nn->values[0], nn->values[1],
               nn->values[2], nn->values[3]);
    }

    dictp n = PHES.begin(), e = PHES.end();
    for (; n != e; n++) {
        pedper pp        = n->first;
        vector<string> v = n->second;
        printf(
#ifdef MS_PRINTF
            "ped %s, per %s; %Iu# %s %s %s %s %s\n",
#else
            "ped %s, per %s; %zu# %s %s %s %s %s\n",
#endif
            pp.ped.c_str(), pp.per.c_str(), v.size(), v[0].c_str(), v[1].c_str(), v[2].c_str(), v[3].c_str(), v[4].c_str());
        fflush(stdout);
    }
}

void phe::free()
{
}

phe PHE;

int phemake(const char *s)
{
    PHE = phe(s);
    return PHE.make();
}

void phehdr(char **names, int *types)
{
    PHE.hdr(names, types);
////    PHE.debug_dump_pydict();
}

void phefree()
{
    PHE.free();
}

int phesearch(const char *ped, const char *per, char **out)
{
    char **p = out;
    vector<string> v;
    int dbg = 0;
    int i, l;
    int ret = PHE.search(ped, per, v);

    if (ret == 0) {
        return 0;
    }
    if (dbg)
        printf(
#ifdef MS_PRINTF
            "phesearch: %s %s %d #%Iu# ",
#else
            "phesearch: %s %s %d #%zu# ",
#endif
            ped, per, ret, v.size());
    for (i = 0, l = (int)v.size(); i < l; i++) {
        if (dbg)
            printf("%d      : %s ", i, v[i].c_str());
        *p++ = strdup(v[i].c_str());
    }
    if (dbg) printf("\n");
    return l;
}

void cdump()
{
    PHE.debug_dump_pydict();
}


#ifdef notdef

void ctest()
{
    PHE.test("/Users/rbaron/tx/mega2test/samoan_iddm_ann/plink.phe");
}

main()
{
    phemake("/Users/rbaron/tx/mega2test/samoan_iddm_ann/plink.phe");
    cdump();
    ctest();
    phefree();
}

#endif
