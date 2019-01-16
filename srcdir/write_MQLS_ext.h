/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

===========================================================================
*/

#ifndef WRITE_MQLS_EXT_H
#define WRITE_MQLS_EXT_H

#include "analysis.h"

class CLASS_MQLS: public CLASS_ANALYSIS {
public:
    CLASS_MQLS() : CLASS_ANALYSIS() {
        _name = "MQLS-XM/KinInbcoef";
        file_name_stem = strdup("mqls");
    }

    ~CLASS_MQLS() {}

    virtual bool loops()  { return true; }

    virtual const char* output_quant_default_value() { return "999.0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

//  virtual bool qtl_allow()        { return true; }
    virtual bool qtl_disallow()        { return true; }

    virtual void ped_ind_defaults(int unique)  {
        OrigIds[0] = 2; /* Unique person name */
        OrigIds[1] = 3; /* consecutive numbers for pedigree */
    }

    void option_menu (char *file_names[], char *prefix, int *combine_chromo, linkage_ped_top *Top);

    void genetic_distance_menu(linkage_ped_top *Top);

    //called to create all files
    void create_output_file(linkage_ped_top *LPedTreeTop, analysis_type *analysis, char *file_names[], int untyped_ped_opt, int *numchr, linkage_ped_top **Top2);

    //creates PLINK tped
    void write_MQLS_genofile(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    //ped file with some condtions on numbering
    void write_MQLS_phenofile(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    //ped file with 4 columns
    void write_KinInbcoef_ped(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    //ped file with 5 columns
    void write_KinInbcoefX_ped(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    //creates a list file of all permutations of family members
    void write_MQLS_listfile(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid);

    //creates a two value prevalence file
    void write_MQLS_prevalence_file(linkage_ped_top *Top, const char *prefix, char *file_names[], const int pwid, const int fwid, double mprev, double fprev);

    //creates the Kininbcoef/X shells and the MQLS shells
    //these shells create the kinfile for MQLS
    void write_MQLS_shell_script(linkage_ped_top *Top, const char *prefix, char *file_names[]);

    void inner_file_names(char **file_names, const char *num, const char *stem /* = "MQLS" */, int *combine_chromo);

    void replace_chr_number(char *file_names[], int numchr);

    virtual void batch_in();

    virtual void batch_out();

    virtual void batch_show();

    char * prevalencefilename;
    double maleprevalence;
    double femaleprevalence;
    bool hasXdata;
    Str additional_arguments;

};

extern CLASS_MQLS            *MQLS;

#endif
