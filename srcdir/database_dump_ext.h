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

#ifndef DATABASE_DUMP_EXT_H
#define DATABASE_DUMP_EXT_H

#include "analysis.h"

class CLASS_DUMP: public CLASS_ANALYSIS {
public:
    CLASS_DUMP() : CLASS_ANALYSIS() {
        _name = "Dump";
        file_name_stem = strdup("dump");
    }
   ~CLASS_DUMP() {}

    virtual bool allele_data_use_name_if_available() { return true; }
    virtual bool allow_covariates()  { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    virtual bool qtl_allow()     { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Always unique ids for persons */
        /* Set the output pedigree and per field choices */
        OrigIds[0] = 5;
        OrigIds[1] = 1;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2) {};


    void create_sh_file(linkage_ped_top *Top, char *file_names[], const int numchr) {};

    void get_file_names(char *file_names[], int has_orig, int has_uniq,
                        int *combine_chromo) {};
    void gen_file_names(char **file_names, char *num) {};
    void replace_chr_number(char *file_names[], int numchr) {};

};

extern CLASS_DUMP            *DUMP;

#endif
