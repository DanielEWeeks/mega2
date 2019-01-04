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

#ifndef WRITE_STRUCTURE_EXT_H
#define WRITE_STRUCTURE_EXT_H

#include "analysis.h"

class CLASS_STRUCTURE: public CLASS_ANALYSIS {
public:
    CLASS_STRUCTURE() : CLASS_ANALYSIS() {
        _name = "Structure";
        _suboption = -1;
        file_name_stem = strdup("structure");
    }
   ~CLASS_STRUCTURE() {}

//  not currently needed
//  virtual bool has_sub_options()    { return true; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
    virtual bool allow_no_map()     { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    virtual const char* output_quant_default_value() { return "?"; }
    virtual const char* output_affect_default_value() { return "?"; }
    virtual bool qtl_allow()        { return true; }

    virtual void ped_ind_defaults(int unique)  {
        OrigIds[0] = 2; /* uniqueIds */
        OrigIds[1] = 2; /* Ped num */
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                linkage_ped_top **Top2);

//    void get_file_names(char *file_names[], char *prefix,
//                        int has_orig, int has_uniq, int *combine_chromo);
    void user_queries(char **file_names,
                      int *combine_chromo, int *create_summary);

    void gen_file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);

};

extern CLASS_STRUCTURE          *STRUCTURE;

#endif
