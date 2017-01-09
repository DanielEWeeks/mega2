/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

*/

#ifndef WRITE_PLINK_FILES_EXT_H
#define WRITE_PLINK_FILES_EXT_H

#include "plink_core_ext.h"

// See plink_core_ext.h for the corresponding numbers associated with these strings...
#define PLINK_SUB_OPTION_LGEN                      "lgen"
#define PLINK_SUB_OPTION_SNP_MAJOR                 "SNP major binary"
#define PLINK_SUB_OPTION_INDIVIDUAL_MAJOR          "Individual major binary"
#define PLINK_SUB_OPTION_PED                       "ped"


class CLASS_PLINK: public CLASS_PLINK_CORE {
public:
    CLASS_PLINK() : CLASS_PLINK_CORE() {
        _name = "PLINK";
	file_name_stem = strdup("plink");
    }
   ~CLASS_PLINK() {}

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "-9"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true;  }
    virtual const char* output_affect_default_value() { return "-9"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    void replace_chr_number(char *file_names[], int numchr);

    void sub_prog_name(int sub_opt, char *subprog);
    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);
    void save_pedsix_file(linkage_ped_top *Top,
			   const int pwid,
			   const int fwid);
    void save_ped_file(linkage_ped_top *Top,
		       const int pwid,
		       const int fwid,
		       const int mwid);
    void save_bed_file(const char *bedfl_name,
		       linkage_ped_top *Top,
		       const int binary_mode_flag);
    virtual void create_sh_file(linkage_ped_top *Top,
			char *file_names[],
			const int numchr);
    void file_names_w_stem(char **file_names, char *num, const char *stem,
                           const int suboption);
};

extern CLASS_PLINK              *TO_PLINK;

#endif /* WRITE_PLINK_FILES_EXT_H */
