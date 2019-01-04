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

#ifndef READ_IMPUTE_HH
#define READ_IMPUTE_HH

#include <fstream>
#include <string>
#include <vector>
#include <queue>

#include "input_ops.hh"
#include "str_utils.hh"
#include "annotated_ped_file.h"

class ImpMarker {
public:
    ImpMarker(Cstr& name, Cstr& chr, Cstr& pos, std::vector<std::string>& alleles, bool read_info): name(name), chr(chr), pos(pos), alleles(alleles) {
        if (! read_info)
            info = 1.0;
        else
            info = 0.0;
        certainty = 0.0;
        skip = false;
    }

public:
    Str name;
    Cstr chr;
    Cstr pos;
    std::vector<std::string> alleles;
    double info;
    double certainty;
    bool skip;
};

typedef vectordb<ImpMarker *> Vecmarkerp;
typedef vectordb<ImpMarker *>::const_iterator Vecmarkerpp;

typedef vectordb<VecsDB> Vecvecs;
typedef vectordb<VecsDB>::const_iterator Vecvecsp;

enum Column_Type {RESERVED, MISSING, PHENO, GENO};
typedef vectordb<Column_Type> Vecct;
typedef vectordb<Column_Type>::const_iterator Vecctp;

struct ProbID {
    ProbID(): dd(0), ii(0), jj(0) {};
    ProbID(double d, int i, int j): dd(d), ii(i), jj(j) {};
    double dd;
    int    ii;
    int    jj;
};
struct ProbIDComp {
    bool operator()(ProbID& a, ProbID& b) {return a.dd < b.dd;}
};

typedef std::priority_queue<ProbID, std::vector<ProbID>, ProbIDComp> ProbQ;


class GenotypeReadHelper
{
public:
    virtual void genotypes_init() {}
    virtual bool genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid,
                                         std::string& pos, std::vector<std::string>& alleles) { return true; }
    virtual void genotypes_skip() {}
    virtual bool genotypes_sample_prob(ProbQ& Q) {return true;}
    virtual void genotypes_end() {}
    virtual bool genotypes_eol(int person) {return false;}
};

class ReadImputed : public Input_Ops
{
public:
    ReadImputed() :  check_format(0), Ncol(0), HDR(2) { };  // two line header
    ~ReadImputed() { (void) markers[0]; };

    virtual bool  use_getops() {return true;}

    virtual void do_menu_display(int &idx, int line_len, int choiceA[]);
    virtual int  do_menu_parse(int choice);
    virtual void do_menu2batch();
    virtual void do_batch2local();

    void show_settings();

    virtual void do_init(Input_Base *inp);
    virtual linkage_locus_top *do_names(const char *&names_fn);
    virtual void do_map(std::vector<m2_map>& additional_maps);
    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop);
    virtual void do_gc();

    void files(char *imp, Cstr& sam) {
        impute_file = imp;
        sample_file = sam;
    }

    virtual void read_input_file();
    virtual void build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons, std::vector<Vecc> &veca);

    void build_internal_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons, GenotypeReadHelper &gh, std::vector<Vecc> &VA);

    void check_indelsNdups();
    void read_info_file();
    void read_sample_file();
    linkage_locus_top *build_impute2_names();
    void build_impute2_map(m2_map& impute_map);
    annotated_ped_rec *build_impute2_ped(linkage_locus_top *LTop, int *num_peds);


public:
    Vecmarkerp    markers;
    int           markers_all;
    int           markers_filtered;
    Vecvecs       people;
    int           people_all;
    int           people_filtered;
    bool          read_info;

    Str           oxford_single_chr;
    double        info_threshold;
    double        hard_call_threshold;
    double        genotype_missing_fraction;
    char          allow_indels;
    char          allow_dups;
    Str           rsid_sep;
    Input_Base   *input;
    int           check_format;

public:
    Vecs          sample_file_hdr2b;
    Str           info_file;
    int           Ncol;


protected:
    const char   *impute_file;
    Str           sample_file;
    m2_map        impute_map;

static
    Str     info_file_hdr;
static
    Str     sample_file_hdr;
    int     HDR;
    Vecs    sample_file_hdr1a;
    Vecs    sample_file_hdr1b;
    Vecs    sample_file_hdr2a;

static const
    int imputed_oxford_single_chr_i         = 8000001;
static const
    int imputed_info_metric_threshold_i     = 8000002;
static const
    int imputed_hard_call_threshold_i       = 8000003;
static const
    int imputed_genotype_missing_fraction_i = 8000004;
static const
    int imputed_missing_code_i              = 8000005;
static const
    int imputed_allow_indels_i              = 8000006;
static const
    int imputed_allow_dups_i                = 8000007;
static const
    int imputed_rsid_sep_i                  = 8000008;

};

class ReadImputedGenotypeReadHelper : public GenotypeReadHelper
{
public:
    virtual void genotypes_init();
    virtual bool genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid, 
                                         std::string& pos, std::vector<std::string>& alleles);
    virtual bool genotypes_sample_prob(ProbQ& Q);
    virtual void genotypes_end();
    virtual bool genotypes_eol(int person);

public:
    std::string impute_file;
    std::ifstream ifs;
    Token token;
    ReadImputed *rip;
    std::string line;
    int line_n;
};

#endif
