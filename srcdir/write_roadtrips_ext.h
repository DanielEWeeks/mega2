/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 2012-2019, University of Pittsburgh. All Rights Reserved.

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

#ifndef WRITE_ROADTRIPS_EXT_H
#define WRITE_ROADTRIPS_EXT_H

#include "analysis.h"

class CLASS_ROADTRIPS: public CLASS_ANALYSIS {
public:
    CLASS_ROADTRIPS() : CLASS_ANALYSIS() {
        _name = "Roadtrips";
        file_name_stem = strdup("roadtrips");
    }
   ~CLASS_ROADTRIPS() {}

//  virtual bool allow_trait_combination()  { return true; }

    virtual bool loops()  { return true; }

    virtual const char* output_quant_default_value() { return "999.0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool forbid_sex_linked_loci() {return true;}

//  virtual bool qtl_allow()        { return true; }
    virtual bool qtl_disallow()        { return true; }

    virtual void ped_ind_defaults(int unique)  {
        OrigIds[0] = 5; /* uniqueIds */
        OrigIds[1] = 3; /* Ped num */
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void create_sh_file(linkage_ped_top *Top, char *file_names[], const int numchr);

    void get_file_names(char *file_names[], int has_orig, int has_uniq,
                        int *combine_chromo);
    void gen_file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);

    virtual void batch_in();

    virtual void batch_out();

    virtual void batch_show();

public:

    double male_prevalence;
    double female_prevalence;
    Str    additional_program_args;
};

extern CLASS_ROADTRIPS            *ROADTRIPS;

#endif
