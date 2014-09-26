/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2014 Robert Baron, Charles P. Kollar,
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

#ifndef WRITE_FBAT_EXT_H
#define WRITE_FBAT_EXT_H

class CLASS_FBAT: public CLASS_ANALYSIS {
public:
    CLASS_FBAT() : CLASS_ANALYSIS() {
        _name = "FBAT";
    }
   ~CLASS_FBAT() {}

    virtual bool allow_affection_liability_class()  { return true; }
//  virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_chr()     { return true; }
//  virtual bool allow_no_trait()  { return true; }
/*
 *  turn off suboptions until there is a good use for them
    virtual bool has_sub_options()  { return true; }
    virtual bool is_sub_option()  { return true; }
 */
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    // missing values are always replaced with an '-'

//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual const char* output_quant_default_value() { return " - "; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool require_physical_map()  { return false; }
/*
    You can select them ... they just can not be a dir ... just .phe file
*/
    virtual bool skip_trait(linkage_locus_top *LocusTop, int trait)  { 
        return (LocusTop->Locus[trait].Type != AFFECTION);
    }

    virtual void ped_ind_defaults(int unique)  {
        /* options that require ids, mostly linkage options and
           gh-like options*/
        OrigIds[0] = 1; /* Renumbered */
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
    void file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

#endif
