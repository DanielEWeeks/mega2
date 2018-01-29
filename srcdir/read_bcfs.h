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

#ifndef READ_BCFS_H
#define READ_BCFS_H

#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <string>

#include "input_ops.hh"
#include "str_utils.hh"
#include "annotated_ped_file.h"

using namespace std;
class BCFMarker {
public:    BCFMarker(Cstr& name, Cstr& chr, int pos, std::vector<std::string>& alleles) {
        this->name = name;
        this->chr = chr;
        this->pos = pos;
        this->alleles = alleles;
    }

public:
    Str name;
    Str chr;
    int pos;
    std::vector<std::string> alleles;
};

typedef vectordb<BCFMarker *> MarkerVector;
typedef vectordb<BCFMarker *>::const_iterator MarkerVectorP;

class ReadBCFs : public Input_Ops
{
public:
    ReadBCFs() { };
    ~ReadBCFs() { };

    virtual bool  use_getops() {return true;}

    virtual void do_menu_display(int &idx, int line_len, int choiceA[]);
    virtual int  do_menu_parse(int choice);
    virtual void do_menu2batch();
    virtual void do_batch2local();

    virtual void do_init(Input_Base *inp);
    void show_settings();

    void read_BCFs(linkage_locus_top *LPedTreeTop);
    linkage_locus_top * do_names(const char *&names_fn);
    void  do_phe_names (char *phe_file, char ***phe_names, int **phe_types, int phe_cols);
    linkage_locus_top * build_BCFs_names();
    void check_bcf_files();
    int build_samples();
    void build_markers();

    virtual void do_map(std::vector<m2_map>& additional_maps);
    void build_bcf_map(m2_map& bcf_map);

    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop);
    void do_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons);

public:
    Str     BCF_path;
    Str     BCF_template;
    Input_Base   *input;
    MarkerVector markers;
    int marker_count;


protected:
    static const
    int site_bcfs_dir_i         = 8000001;
    static const
    int site_bcfs_template_i     = 8000002;

    vector<string> filelist;
    int filecount;

    char * pedfile;

    vector<string> samples;
    int num_samples;

    char * phefile;
    int phecols;
    char ***phenames;
    int **phetypes;

};


#endif

