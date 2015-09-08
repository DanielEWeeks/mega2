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

#ifndef INPUT_HH
#define INPUT_HH

class Input_Impute;
class Input_BGEN;

// #include "types.hh" // common.h includes types.hh
#include "str_utils.hh"
#include "read_impute.hh"
#include "read_impute_bgen.hh"

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
    in_format_compressed_VCF = 6,
    in_format_VCF = 7,
    in_format_imputed = 8,
    in_format_bgen    = 9,
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

    Input_Files() :
        pedfl(&mega2_input_files[PEDIGREE]),
        locusfl(&mega2_input_files[LOCUS]),
        mapfl(&mega2_input_files[MAP]),
        pmapfl(&mega2_input_files[PMAP]),
        omitfl(&mega2_input_files[OMIT]),
        freqfl(&mega2_input_files[FREQ]),
        penfl(&mega2_input_files[PEN]),
        bedfl(&mega2_input_files[BED]),
        phefl(&mega2_input_files[PHEfl]) {  };

};

class Input_Base {
public:
    Input_Base(INPUT_FORMAT_t i): input_format(i), plink(0), xcf(0), req_aux_file(0), req_locus_file(1), req_map_file(1), req_stem_flag(0)  {}
    virtual ~Input_Base() {};

    virtual boolean has_menu_display() { return false; }
    virtual boolean has_menu_parse()   { return false; }
    virtual boolean has_menu2batch()   { return false; }
    virtual boolean has_batch2local()  { return false; }

    virtual void do_menu_display(int &idx, int line_len, int choiceA[])   { }
    virtual int  do_menu_parse(int choice) { return 0; }   // 0 indicates no match ; but false (above) means this is not called.
    virtual void do_menu2batch() { }
    virtual void do_batch2local() { }


    virtual boolean has_init()  { return false; }
    virtual boolean has_names() { return false; }
    virtual boolean has_map()   { return false; }
    virtual boolean has_ped()   { return false; }
    virtual boolean has_gc()    { return false; }

//  these functions are define the corresponding function above returns true;
    virtual void do_init()  { }
    virtual linkage_locus_top *do_names(const char *&names_fn) { return (linkage_locus_top *)0; }
    virtual void do_map(std::vector<m2_map>& additional_maps) { }
    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop) { return (linkage_ped_top *) 0; }
    virtual void do_gc() { }

public:
    INPUT_FORMAT_t input_format;
    Input_Files input_files;
    Globals G;
    Str     MissingCodes;
    Sets    MissingCodesSet;

    int plink;
    int xcf;
    int req_aux_file;
    int req_locus_file;
    int req_map_file;
    int req_stem_flag;
};

class Input_Old : public Input_Base {
public:
    Input_Old(INPUT_FORMAT_t i): Input_Base(i) {}
    virtual ~Input_Old() {};
};

////////////////

class Input_Mega2 : public Input_Old {
public:
    Input_Mega2(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Mega2() {};
};

class Input_Linkage : public Input_Old {
public:
    Input_Linkage(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Linkage() {};
};

class Input_Extended_Linkage : public Input_Old {
public:
    Input_Extended_Linkage(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Extended_Linkage() {};
};

////////////////

class Input_PLINK_Common : public Input_Old {
public:
    Input_PLINK_Common(INPUT_FORMAT_t i): Input_Old(i) {
        plink          = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
    }
    virtual ~Input_PLINK_Common() {};
};

class Input_PED_Binary : public Input_PLINK_Common {
public:
    Input_PED_Binary(INPUT_FORMAT_t i): Input_PLINK_Common(i) {
        req_aux_file   = 1;
    }
    virtual ~Input_PED_Binary() {};
};

class Input_PED : public Input_PLINK_Common {
public:
    Input_PED(INPUT_FORMAT_t i): Input_PLINK_Common(i) {}
    virtual ~Input_PED() {};
};

////////////////

class Input_VCF_Common : public Input_Old {
public:
    Input_VCF_Common(INPUT_FORMAT_t i): Input_Old(i) {
        xcf = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
}
    virtual ~Input_VCF_Common() {};
};

class Input_VCF_Binary : public Input_VCF_Common {
public:
    Input_VCF_Binary(INPUT_FORMAT_t i): Input_VCF_Common(i) {}
    virtual ~Input_VCF_Binary() {};
};

class Input_VCF_Compressed : public Input_VCF_Common {
public:
    Input_VCF_Compressed(INPUT_FORMAT_t i): Input_VCF_Common(i) {}
    virtual ~Input_VCF_Compressed() {};
};

class Input_VCF : public Input_VCF_Common {
public:
    Input_VCF(INPUT_FORMAT_t i): Input_VCF_Common(i) {}
    virtual ~Input_VCF() {};
};

////////////////

class Input_Impute;

class Input_Impute : public Input_Base {
public:
    Input_Impute(INPUT_FORMAT_t i): Input_Base(i) {
        req_aux_file   = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
    }
    virtual ~Input_Impute() {};

    virtual boolean has_menu_display() { return true; }
    virtual boolean has_menu_parse()   { return true; }
    virtual boolean has_menu2batch()   { return true; }
    virtual boolean has_batch2local()  { return true; }

    virtual void do_menu_display(int &idx, int line_len, int choiceA[])    { Obj.do_menu_display(idx, line_len, choiceA); }
    virtual int  do_menu_parse(int choice) { return Obj.do_menu_parse(choice); }
    virtual void do_menu2batch() { Obj.do_menu2batch(); }
    virtual void do_batch2local() { Obj.do_batch2local(); }


    virtual boolean has_init()  { return true; }
    virtual boolean has_names() { return true; }
    virtual boolean has_map()   { return true; }
    virtual boolean has_ped()   { return true; }
    virtual boolean has_gc()    { return true; }

    virtual void do_init()  { Obj.do_init(this); }
    virtual linkage_locus_top *do_names(const char *&names_fn) { return Obj.do_names(names_fn); }
    virtual void do_map(std::vector<m2_map>& additional_maps) { Obj.do_map(additional_maps); }
    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop) { return Obj.do_ped(LTop); }
    virtual void do_gc() { Obj.do_gc(); }

public:

    ReadImputed Obj;

};

class Input_BGEN : public Input_Impute {
public:
    Input_BGEN(INPUT_FORMAT_t i) : Input_Impute(i) {};
   ~Input_BGEN() {};

public:

    virtual void do_init()  { bgen.do_init(this); }
    virtual linkage_ped_top *do_ped(linkage_locus_top *LTop) { return bgen.do_ped(LTop); }

public:

    ReadBgen bgen;

};

class Input_Traditional : public Input_Old {
public:
    Input_Traditional(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Traditional() {};
};

/*
 *      Requires corresponding addition/change to
 *              input.hh: add to enum INPUT_FORMAT {                   | around line 44
 *              bstch_input.cpp: check_batch_items                     | around line 427
 *              user_input.cpp: INPUT_FORMAT_STR[]                     | around line  107
 *              user_input.cpp: "if (choice_ == file_format_i)" case   | around line 1187
 *              user_input.cpp: "if (Input_Format == in_format_xxx)" case | around line 963
 *              mega2.cpp: "if (Input_Format == in_format_bgen)        | around line 1087
 */
class InputCreate {
public:
    static
    Input_Base *createinput(INPUT_FORMAT in_format) {
        switch(in_format) {
        case in_format_mega2:
	default:
            return new Input_Mega2(in_format);
            break;
        case in_format_linkage:
            return new Input_Linkage(in_format);
            break;
        case in_format_extended_linkage:
            return new Input_Extended_Linkage(in_format);
            break;
        case in_format_binary_PED:
            return new Input_PED_Binary(in_format);
            break;
        case in_format_PED:
            return new Input_PED(in_format);
            break;
        case in_format_binary_VCF:
            return new Input_VCF_Binary(in_format);
            break;
        case in_format_compressed_VCF:
            return new Input_VCF_Compressed(in_format);
            break;
        case in_format_VCF:
            return new Input_VCF(in_format);
            break;
        case in_format_imputed:
            return new Input_Impute(in_format);
            break;
        case in_format_bgen:
            return new Input_BGEN(in_format);
            break;
        case in_format_traditional:
            return new Input_Traditional(in_format);
            break;
        }
    }
};


#endif
