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


#ifndef WRITE_VCF_EXT_H
#define WRITE_MINIWRITE_VCF_EXT_HMAC_EXT_H

#include "analysis.h"

class CLASS_VCF: public CLASS_ANALYSIS {
public:
    CLASS_VCF() : CLASS_ANALYSIS() {
        _name = "VCF";
        file_name_stem = strdup("vcf");
    }

    ~CLASS_VCF() { }

    virtual bool allow_trait_combination()  { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
//  virtual bool allow_no_map()  { return true; }

    virtual bool loops() { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    //virtual const char* output_quant_default_value() { return "-9"; }
    //virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true;  }
    //virtual const char* output_affect_default_value() { return "-9"; }
    //virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool require_physical_map() { return true; }


    virtual void ped_ind_defaults(int unique) {
        OrigIds[0] = 6; /* Keeps Per ID that was input*/
        OrigIds[1] = 6; /* Keeps Ped ID that was input*/
    }

    void create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2);

    void write_VCF_file(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_ped(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_pheno(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_map(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_freq(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_pen(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void convert_vcf_bcf(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void convert_vcf_vcfgz(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    void write_VCF_sh(linkage_ped_top *Top, const char *prefix, char *file_names[]);

    void option_menu (char *file_names[], char *prefix, int *combine_chromo);

    virtual void batch_in();

    virtual void batch_out();

    void gen_file_names(char **file_names, char *num);

    void inner_file_names(char **file_names, const char *num, const char *stem, int *combine_chromo);

    void replace_chr_number(char *file_names[], int numchr);


};

extern CLASS_VCF            *VCF;

#endif

