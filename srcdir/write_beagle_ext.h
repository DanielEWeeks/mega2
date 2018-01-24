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

#ifndef WRITE_BEAGLE_EXT_H
#define WRITE_BEAGLE_EXT_H

#include "analysis.h"

class CLASS_BEAGLE: public CLASS_ANALYSIS {
public:
    CLASS_BEAGLE() : CLASS_ANALYSIS() {
        _name = "Beagle";
    }
   ~CLASS_BEAGLE() {}

    virtual bool has_sub_options()    { return true; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
    virtual bool allow_no_map()     { return true; }
    virtual bool allow_trait_combination()  { return true; }
    //virtual bool forbid_trait_directories()  { return true; }
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    virtual const char* output_quant_default_value() { return "?"; }
    virtual const char* output_affect_default_value() { return "?"; }
    virtual bool qtl_allow()        { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

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

    void get_file_names(char *file_names[], char *prefix,
                        int has_orig, int has_uniq, int *combine_chromo);
    void gen_file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);

    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);
};

#define BEAGLE_UNPHASED_UNRELATED_SUBOPTION         1
class CLASS_BEAGLE_UNPHASED_UNRELATED: public CLASS_BEAGLE {
public:
    CLASS_BEAGLE_UNPHASED_UNRELATED() : CLASS_BEAGLE() {
        _subname = "Unphased_Unrelated";
        _suboption = BEAGLE_UNPHASED_UNRELATED_SUBOPTION;
    }
    // Here the parents are not important, so we want all of the data tocome through...
    virtual bool allow_missing_parent_in_linkage_input()  { return true; }
};

#define BEAGLE_UNPHASED_TRIO_SUBOPTION               2
class CLASS_BEAGLE_UNPHASED_TRIO: public CLASS_BEAGLE {
public:
    CLASS_BEAGLE_UNPHASED_TRIO() : CLASS_BEAGLE() {
        _subname = "Unphased_Trio";
        _suboption = BEAGLE_UNPHASED_TRIO_SUBOPTION;
    }
    // Trio requires both parents to exist, and imposes additional requirements...
};

#define BEAGLE_UNPHASED_PAIR_SUBOPTION               3
class CLASS_BEAGLE_UNPHASED_PAIR: public CLASS_BEAGLE {
public:
    CLASS_BEAGLE_UNPHASED_PAIR() : CLASS_BEAGLE() {
        _subname = "Unphased_Pair";
        _suboption = BEAGLE_UNPHASED_PAIR_SUBOPTION;
    }
    // So that the input will allow <mother, child> and <father, child> pairs
    // to come through with Linkage input...
    virtual bool allow_missing_parent_in_linkage_input()  { return true; }
};

extern CLASS_BEAGLE             *BEAGLE;
/**
   These classes define the sub-options associated with the different Genotype file
   that Beagle can produce. See the section on 'Beagle Genotype files'.
 */
extern CLASS_BEAGLE_UNPHASED_UNRELATED *BEAGLE_UNPHASED_UNRELATED;
extern CLASS_BEAGLE_UNPHASED_TRIO *BEAGLE_UNPHASED_TRIO;
extern CLASS_BEAGLE_UNPHASED_PAIR *BEAGLE_UNPHASED_PAIR;

#endif
