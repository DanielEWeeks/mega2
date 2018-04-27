/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

#ifndef INPUT_HH
#define INPUT_HH

#include "input_ops.hh"

#include "read_impute.hh"
#include "read_impute_bgen.hh"
#include "read_impute_bgen2.hh"
#include "read_bcfs.h"


class Input_Old : public Input_Base {
public:
    Input_Old(INPUT_FORMAT_t i): Input_Base(i) {}
    virtual ~Input_Old() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

////////////////

class Input_Mega2 : public Input_Old {
public:
    Input_Mega2(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Mega2() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_Linkage : public Input_Old {
public:
    Input_Linkage(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Linkage() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_Extended_Linkage : public Input_Old {
public:
    Input_Extended_Linkage(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Extended_Linkage() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
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
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_PED_Binary : public Input_PLINK_Common {
public:
    Input_PED_Binary(INPUT_FORMAT_t i): Input_PLINK_Common(i) {
        req_aux_file   = 1;
    }
    virtual ~Input_PED_Binary() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_PED : public Input_PLINK_Common {
public:
    Input_PED(INPUT_FORMAT_t i): Input_PLINK_Common(i) {}
    virtual ~Input_PED() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
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
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_VCF_Binary : public Input_VCF_Common {
public:
    Input_VCF_Binary(INPUT_FORMAT_t i): Input_VCF_Common(i) {}
    virtual ~Input_VCF_Binary() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_VCF_Compressed : public Input_VCF_Common {
public:
    Input_VCF_Compressed(INPUT_FORMAT_t i): Input_VCF_Common(i) {}
    virtual ~Input_VCF_Compressed() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

class Input_VCF : public Input_VCF_Common {
public:
    Input_VCF(INPUT_FORMAT_t i): Input_VCF_Common(i) {
        xcf = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
    }
    virtual ~Input_VCF() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    Input_Ops Ops;
};

////////////////

class Input_Impute : public Input_Base {
public:
    Input_Impute(INPUT_FORMAT_t i): Input_Base(i) {
        req_aux_file   = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
        impute         = 1;
    }
    virtual ~Input_Impute() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    ReadImputed Ops;
};

class Input_BGEN2 : public Input_Base {
public:
    Input_BGEN2(INPUT_FORMAT_t i) : Input_Base(i) {
        req_aux_file   = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
        impute         = 1;
    };
   ~Input_BGEN2() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    ReadBgen2 Ops;
};

class Input_BGEN : public Input_Base {
public:
    Input_BGEN(INPUT_FORMAT_t i) : Input_Base(i) {
        req_aux_file   = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag  = 1;
        impute         = 1;
    };
   ~Input_BGEN() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    ReadBgen Ops;
};


class Input_BCFs: public Input_Base{
public:
    Input_BCFs(INPUT_FORMAT_t i): Input_Base(i) {
        extern int SORT_HETEROZYGOTE;
        xcf = 1;
        req_locus_file = 0;
        req_map_file   = 0;
        req_stem_flag = 1;
        SORT_HETEROZYGOTE = 0;
    }
    virtual ~Input_BCFs() {};
    virtual Input_Ops *GetOps() {return &Ops;};
    ReadBCFs Ops;
};

class Input_BCFS: public Input_BCFs {
public:
    Input_BCFS(INPUT_FORMAT_t i): Input_BCFs(i) {};
};

class Input_GZCFS: public Input_BCFs {
public:
    Input_GZCFS(INPUT_FORMAT_t i): Input_BCFs(i) {};
};

class Input_VCFS: public Input_BCFs {
public:
    Input_VCFS(INPUT_FORMAT_t i): Input_BCFs(i) {};
};

class Input_Traditional : public Input_Old {
public:
    Input_Traditional(INPUT_FORMAT_t i): Input_Old(i) {}
    virtual ~Input_Traditional() {};
    virtual Input_Ops *GetOps() {return &Ops;};
//  virtual Input_Ops &GetOps() {return Ops;};
    Input_Ops Ops;
};

#endif
