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

#ifndef WRITE_MINIMAC_EXT_H
#define WRITE_MINIMAC_EXT_H

#include "analysis.h"
#include "write_plink_ext.h"


class CLASS_MINIMAC: public CLASS_PLINK {
public:
    CLASS_MINIMAC() : CLASS_PLINK() {
        _name = "MINIMAC";
        file_name_stem = strdup("minimac");
    }

    ~CLASS_MINIMAC() { }


    virtual bool allow_no_trait() { return true; }

    virtual bool loops() { return true; }

    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

    /*  this is necessary because we are inheriting from CLASS_PLINK not CLASS_ANALYSIS */
    virtual const char* output_quant_default_value() { return "-9"; }

    virtual const char* output_affect_default_value() { return "-9"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool no_missing_menu() { return true; }

    virtual bool no_trait_covariate_menu() { return true; }

    virtual bool forbid_sex_linked_loci() { return true; }

    virtual bool require_physical_map() { return false; }

    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique) {
        OrigIds[0] = 6; /* Keeps Per ID that was input*/
        OrigIds[1] = 6; /* Keeps Ped ID that was input*/
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
                            analysis_type *analysis,
                            char *file_names[],
                            int untyped_ped_opt,
                            int *numchr, linkage_ped_top **Top2);

    void user_queries(char **file_names_array, int *combine_chromo, int *create_summary);

    void minimac_option_menu(char *file_names[], char *prefix);

    void gen_file_names(char **file_names, char *num);
    //{
//        file_names_w_stem(file_names, num, file_name_stem,
//                          PLINK_SUB_OPTION_SNP_MAJOR_INT);
    //};

    void replace_chr_number(char *file_names[], int numchr);

    virtual void batch_in();

    virtual void batch_out();

    bool has_sub_options()    { return true; }
    void sub_prog_name(int sub_opt, char *subprog);
    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *subprog_name, analysis_type *analysis);

    void create_sh_file(linkage_ped_top *Top,
                        char *file_names_array[],
                        const int numchr);

    //void write_MINIMAC_sh(linkage_ped_top *Top, char *file_names[]);

public:
    Str reference_map_file;
    Str map_pre;
    Str map_post;
    Str shapeit_reference_haplotype_file;
    Str s_haplotype_pre;
    Str s_haplotype_post;
    Str reference_legend_file;
    Str legend_pre;
    Str legend_post;
    Str reference_sample_file;
    Str minimac_reference_haplotype_file;
    Str m_haplotype_pre;
    Str m_haplotype_post;
};


extern CLASS_MINIMAC            *MINIMAC;

#endif
