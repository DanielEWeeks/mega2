/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
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

#ifndef WRITE_PANGAEA_EXT_H
#define WRITE_PANGAEA_EXT_H

#include "analysis.h"

class CLASS_PANGAEA: public CLASS_ANALYSIS {
public:
    CLASS_PANGAEA() : CLASS_ANALYSIS() {
        _name = "PANGAEA";
    }
   ~CLASS_PANGAEA() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_genetic_map()  { return (_suboption <= 2) ? false : true; }
    virtual bool allow_no_map()     { return (_suboption <= 2) ? false : true; }
    virtual bool allow_trait_combination()  { return true; }
//  virtual bool forbid_trait_directories()  { return true; }
    virtual bool has_sub_options()  { return true; }
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

    virtual const char* output_quant_default_value() { return "999.0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow()        { return true; }

    virtual void ped_ind_defaults(int unique)  {
        OrigIds[0] = 4; /* uniqueIds */
        OrigIds[1] = 2; /* Ped num */
    }

    virtual void sub_prog_name(int sub_opt, char *subprog);

    virtual void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    virtual void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void get_file_names(char *file_names[], char *prefix,
                        int has_orig, int has_uniq, int *combine_chromo);
    void gen_file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

extern CLASS_PANGAEA            *PANGAEA;

#endif
