/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 2012-2018, University of Pittsburgh. All Rights Reserved.

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

*/

#ifndef WRITE_EIGENSTRAT_FILES_EXT_H
#define WRITE_EIGENSTRAT_FILES_EXT_H

#include "plink_core_ext.h"

// See plink_core_ext.h for the corresponding numbers associated with these strings...
#define EIGENSTRAT_SUB_OPTION_SNP_MAJOR                 "PACKEDPED"
#define EIGENSTRAT_SUB_OPTION_PED                       "PED"


class CLASS_EIGENSTRAT: public CLASS_PLINK_CORE {
public:
    CLASS_EIGENSTRAT() : CLASS_PLINK_CORE() {
        _name = "Eigenstrat";
	file_name_stem = strdup("eigenstrat");
    }
   ~CLASS_EIGENSTRAT() {}

    virtual const char* output_quant_default_value() { return "-100.0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "-100"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    void replace_chr_number(char *file_names[], int numchr);

    void sub_prog_name(int sub_opt, char *subprog);
    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    void create_output_file(linkage_ped_top *LPedTreeTop,
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
    void create_sh_file(linkage_ped_top *Top,
			char *file_names[],
			const int numchr);
    void file_names_w_stem(char **file_names, char *num, const char *stem,
                           const int suboption);
};

extern CLASS_EIGENSTRAT         *EIGENSTRAT;

#endif /* WRITE_EIGENSTRAT_FILES_EXT_H */
