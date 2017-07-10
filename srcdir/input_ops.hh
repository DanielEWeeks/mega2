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

#ifndef INPUT_OPS_HH
#define INPUT_OPS_HH

#include "str_utils.hh"

class Input_Base;
extern Input_Base *Input;

typedef
enum INPUT_FORMAT {
    in_format_mega2 = 0,
    in_format_linkage = 1,
    in_format_extended_linkage = 2,
    in_format_binary_PED = 3,
    in_format_PED = 4,
    in_format_binary_VCF = 5,
    in_format_bcfs = 6,
    in_format_compressed_VCF = 7,
    in_format_VCF = 8,
    in_format_imputed = 9,
//  in_format_bgen  = 9,
//  in_format_bgen2 = 10,
    in_format_bgen2 = 10,
    in_format_bgen  = 11,

    in_format_traditional = 100,
} INPUT_FORMAT_t;

class Input_Files {
public:
    char **pedfl;
    char **locusfl;
    char **mapfl;
    char **pmapfl;
    char **omitfl;
    char **freqfl;
    char **penfl;
    char **bedfl;
    char **phefl;
    char **reffl;
    char **bcfsdir;
    char **bcfstemp;

    Input_Files() :
        pedfl(&mega2_input_files[PEDIGREE]),
        locusfl(&mega2_input_files[LOCUS]),
        mapfl(&mega2_input_files[MAP]),
        pmapfl(&mega2_input_files[PMAP]),
        omitfl(&mega2_input_files[OMIT]),
        freqfl(&mega2_input_files[FREQ]),
        penfl(&mega2_input_files[PEN]),
        bedfl(&mega2_input_files[BED]),
        phefl(&mega2_input_files[PHEfl]),
        reffl(&mega2_input_files[REFfl]),
        bcfsdir(&mega2_input_files[BCFSDIR]),
        bcfstemp(&mega2_input_files[BCFSTEMP]) {  };

};

class Input_Ops
{
public:
    virtual bool use_getops() { return false; }

    virtual void do_menu_display(int &idx, int line_len, int choiceA[])   { }
    virtual int  do_menu_parse(int choice) { return 0; }   // 0 indicates no match ; but false (above) means this is not called.
    virtual void do_menu2batch()  { }
    virtual void do_batch2local() { }

    virtual void do_init(Input_Base *inp)  { }
    virtual linkage_locus_top *do_names(const char *&names_fn) { return (linkage_locus_top *)0; }
    virtual void do_map(std::vector<m2_map>& additional_maps)  { }
    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop)   { return (linkage_ped_top *) 0; }
    virtual void do_gc() { }
};

class Input_Base {
public:
    Input_Base(INPUT_FORMAT_t i): input_format(i), plink(0), xcf(0), impute(0), req_aux_file(0), req_locus_file(1), req_map_file(1), req_stem_flag(0)  {}
    virtual ~Input_Base() {}
    virtual Input_Ops *GetOps() {return (Input_Ops *) 0;}

public:
    INPUT_FORMAT_t input_format;
    Input_Files input_files;
    Globals G;
    Str     MissingCodes;
    Sets    MissingCodesSet;

    int plink;
    int xcf;
    int impute;
    int req_aux_file;
    int req_locus_file;
    int req_map_file;
    int req_stem_flag;
};

#endif
