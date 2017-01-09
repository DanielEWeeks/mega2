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

#ifndef WRITE_SHAPEIT_FILES_EXT_H
#define WRITE_SHAPEIT_FILES_EXT_H

#include "types.hh"
#include "write_plink_ext.h"


class CLASS_SHAPEIT: public CLASS_PLINK {
public:
    CLASS_SHAPEIT() : CLASS_PLINK() {
        _name = "SHAPEIT";
	file_name_stem = strdup("shapeit");
    }
   ~CLASS_SHAPEIT() {}

    bool has_sub_options()    { return true; }
    virtual void sub_prog_name(int sub_opt, char *subprog);
    virtual void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);
    virtual void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);

    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

/*  this is necessary because we are inheriting from CLASS_PLINK not CLASS_ANALYSIS */
    virtual const char* output_quant_default_value() { return "-9"; }

    virtual const char* output_affect_default_value() { return "-9"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num) {
	file_names_w_stem(file_names, num, file_name_stem,
                          PLINK_SUB_OPTION_SNP_MAJOR_INT);
    }

    void save_pedsix_file(linkage_ped_top *Top,
                          const int pwid,
                          const int fwid);

    void save_pheno_file(linkage_ped_top *Top,
			 const int pwid, const int fwid) {}; // Nothing here for shapeit

    void create_sh_file(linkage_ped_top *Top,
			char *file_names_array[],
			const int numchr);

    void create_sh_file_phased(linkage_ped_top *Top,
                               char *file_names_array[],
                               const int numchr);

    void create_sh_file_check(linkage_ped_top *Top,
                              char *file_names_array[],
                              const int numchr);

    void user_queries(char **file_names_array,
                      int *combine_chromo, int *create_summary);

    virtual void batch_in();

    virtual void batch_out();

    virtual void batch_show();

public:
    Str rdir;
    Str rpre;
    Str rpost;
};

extern CLASS_SHAPEIT *SHAPEIT;

#endif /* WRITE_SHAPEIT_FILES_EXT_H */
