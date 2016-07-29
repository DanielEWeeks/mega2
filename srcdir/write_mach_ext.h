/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
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

#ifndef WRITE_MACH_EXT_H
#define WRITE_MACH_EXT_H

#include "analysis.h"
#include "write_plink_ext.h"


class CLASS_MACH: public CLASS_ANALYSIS {
public:
    CLASS_MACH() : CLASS_ANALYSIS() {
        _name = "MACH";
        file_name_stem = strdup("mach");
    }
    ~CLASS_MACH() {}

//	virtual bool allow_affection_liability_class()  { return true; }
//  virtual bool allow_no_aff_trait()     { return true; }
//	virtual bool allow_no_chr()     { return true; }
    virtual bool allow_no_trait()  { return true; }
/*
 *  turn off suboptions until there is a good use for them
    virtual bool has_sub_options()  { return true; }
    virtual bool is_sub_option()  { return true; }
 */


    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

    //gives an error if user selects X,Y,etc.
    virtual bool forbid_sex_linked_loci()  { return true; }

    // missing values are always replaced with an '-'


//  virtual bool output_quant_can_define_missing_value() { return false; }
//	virtual const char* output_quant_default_value() { return " - "; }
//	virtual const char* output_affect_default_value() { return "0"; }
//	virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool require_physical_map()  { return false; }

    virtual bool allele_data_use_name_if_available() { return true; }

/*	virtual bool skip_trait(linkage_locus_top *LocusTop, int trait)  {
		return (LocusTop->Locus[trait].Type != AFFECTION);
	}
 */

    virtual void ped_ind_defaults(int unique)  {
        OrigIds[0] = 6; /* Keeps Per ID that was input*/
        OrigIds[1] = 6; /* Keeps Ped ID that was input*/
    }


    void create_output_file(linkage_ped_top *LPedTreeTop,
                            analysis_type *analysis,
                            char *file_names[],
                            int untyped_ped_opt,
                            int *numchr, linkage_ped_top **Top2);

    void mach_option_menu (char *file_names[],char *prefix);

    void gen_file_names(char **file_names, char *num);

    void replace_chr_number(char *file_names[], int numchr);

    virtual void batch_in();

    virtual void batch_out();

    public:
        Str mach_reference_haplotype_file;
        Str haplotype_pre;
        Str haplotype_post;
};


extern CLASS_MACH            *MACH;

#endif
