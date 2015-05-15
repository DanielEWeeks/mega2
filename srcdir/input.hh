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

// #include "types.hh" // common.h includes types.hh
#include "str_utils.hh"
#include "read_impute.hh"

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

class Input {
public:
    Input(INPUT_FORMAT_t i): input_format(i) {}
    virtual ~Input() {};

public:
    INPUT_FORMAT_t input_format;
    Input_Files input_files;
    Globals G;
};

class Input_Old : public Input {
public:
    Input_Old(INPUT_FORMAT_t i): Input(i) {}
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
    Input_PLINK_Common(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_PLINK_Common() {};
};

class Input_PED_Binary : public Input_PLINK_Common {
public:
    Input_PED_Binary(INPUT_FORMAT_t i): Input_PLINK_Common(i) {}
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
    Input_VCF_Common(INPUT_FORMAT_t i): Input_Old(i) {}
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

class Input_Impute : public Input {
public:
    Input_Impute(INPUT_FORMAT_t i): Input(i) {}
    virtual ~Input_Impute() {};

public:
    ReadImputed Obj;
};

class Input_Traditional : public Input_Old {
public:
    Input_Traditional(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Traditional() {};
};

class InputCreate {
public:
    static
    Input *createinput(INPUT_FORMAT in_format) {
        switch(in_format) {
        case in_format_mega2:
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
        case in_format_traditional:
            return new Input_Traditional(in_format);
            break;
        }
    }
};


#endif
