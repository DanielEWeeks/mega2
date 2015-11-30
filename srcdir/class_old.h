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

#ifndef CLASS_OLD_H
#define CLASS_OLD_H

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

#include "analysis.h"
#include "write_plink_ext.h"

// documentation shows this as disabled.
class CLASS_APM: public CLASS_ANALYSIS {
public:
    CLASS_APM() : CLASS_ANALYSIS() {
        _name = "APM";
    }
   ~CLASS_APM() {}

    virtual bool is_enabled() { return false; }

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

// documentation shows this as disabled.
class CLASS_APM_MULT: public CLASS_ANALYSIS {
public:
    CLASS_APM_MULT() : CLASS_ANALYSIS() {
        _name = "APM-Mult";
    }
   ~CLASS_APM_MULT() {}

    virtual bool is_enabled() { return false; }

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_ASPEX: public CLASS_ANALYSIS {
public:
    CLASS_ASPEX() : CLASS_ANALYSIS() {
        _name = "Aspex";
    }
   ~CLASS_ASPEX() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool has_sub_options()    { return true; }
    virtual bool nuke_opts()    { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_disallow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void sub_prog_name(int sub_opt, char *subprog);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_ALLEGRO: public CLASS_ANALYSIS {
public:
    CLASS_ALLEGRO() : CLASS_ANALYSIS() {
        _name = "Allegro";
    }
   ~CLASS_ALLEGRO() {}

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_CRANEFOOT: public CLASS_ANALYSIS {
public:
    CLASS_CRANEFOOT() : CLASS_ANALYSIS() {
        _name = "CRANEFOOT";
    }
   ~CLASS_CRANEFOOT() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_chr()  { return true; }
    virtual bool allow_no_map()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }

//  virtual const char* output_quant_default_value() { return "unknown"; }
//  virtual const char* output_affect_default_value() { return "\t11"; }

    virtual bool qtl_allow()     { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_GHMLB: public CLASS_ANALYSIS {
public:
    CLASS_GHMLB() : CLASS_ANALYSIS() {
        _name = "MLBQTL";
    }
   ~CLASS_GHMLB() {}

    virtual bool allow_no_trait()  { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool nuke_opts()     { return true; }

    // MLB_3.0_Source/gh.help
    // A 0 in any of the disease phenotype or marker genotype positions
    // (as in the the genotypes for the third marker above) indicates
    // missing data. See the file linkped.pre as an example.
    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Person ids and consecutive numbering for pedigrees */
        OrigIds[0] = 5;
        OrigIds[1] = 3;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_GENEHUNTER: public CLASS_ANALYSIS {
public:
    CLASS_GENEHUNTER() : CLASS_ANALYSIS() {
        _name = "GeneHunter";
    }
   ~CLASS_GENEHUNTER() {}

    virtual bool allow_trait_combination()  { return true; }
    virtual bool allow_covariates()  { return true; }
    virtual bool forbid_sex_linked_loci()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

//
// NOTE: Needs an sex average average map.
class CLASS_GENEHUNTERPLUS: public CLASS_ANALYSIS {
public:
    CLASS_GENEHUNTERPLUS() : CLASS_ANALYSIS() {
        _name = "Genehunter-Plus";
    }
   ~CLASS_GENEHUNTERPLUS() {}

    virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_trait()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

/**
   Simwalk is a fortran program that is created by by linking with "vintage" Mendel
   (v 3.3.5). Download vintage Mendel, and any of the SimWalk executables; you will
   be using the fortran source which is in the "CODE" directory. Download from:
   <http://watson.hgen.pitt.edu/register>, or <http://www.genetics.ucla.edu/software>.

   You can build it as follows:

   #!/bin/sh
   cd /usr/local/src
   # unpack Mendel and SimWalk...
   tar xf SimWalk291_mac_g4-g5.tar
   tar xf Mendel-335_unix.tar
   cd SimWalk291/CODE
   # create a link to the Mendel library, and compile...
   ln -s ../../Mendel-335_unix/MENDEL.FOR mendel.f
   g77 -o simwalk2 -O simwalk2.f mendel.f 
   cp simwalk2 /usr/local/bin

   OS X FORTRAN COMPILER NOTES:
   To get g77 working on OS X Lion (10.7) you will need to install the the latest version
   of the command line tools from Apple 'xcode452cltools10_76938212a.dmg', and then install
   the gnu fortran code from sourceforge 'gfortran-lion.tar.gz' according to their
   instructions (e.g., as user root, cd /; tar zxf gfortran-lion.tar.gz). The source forge
   package will install files in '/usr/local'. Installation of the command line tools 'alone'
   will install a 'gfortran' but it is broken. To get a working version you need to install
   the sourceforge version over top of it.
*/
class CLASS_SIMWALK2: public CLASS_ANALYSIS {
public:
    CLASS_SIMWALK2() : CLASS_ANALYSIS() {
        _name = "SimWalk2";
    }
   ~CLASS_SIMWALK2() {}

    virtual bool has_sub_options()    { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual const char* output_quant_default_value() { return ""; }
    virtual const char* output_affect_default_value() { return ""; }

    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);
};

class CLASS_HAPLOTYPE: public CLASS_SIMWALK2 {
public:
    CLASS_HAPLOTYPE() : CLASS_SIMWALK2() {
        _subname = "Haplotype analysis";
      _suboption = 1;
    }
   ~CLASS_HAPLOTYPE() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool is_sub_option()  { return true; }
//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual bool qtl_disallow()  { return true; }
    virtual bool simwalk2()      { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_LOCATION: public CLASS_SIMWALK2 {
public:
    CLASS_LOCATION() : CLASS_SIMWALK2() {
        _subname = "Parametric linkage";
      _suboption = 2;
    }
   ~CLASS_LOCATION() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool is_sub_option()  { return true; }
//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual bool qtl_disallow() { return true; }
    virtual bool simwalk2()     { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_NONPARAMETRIC: public CLASS_SIMWALK2 {
public:
    CLASS_NONPARAMETRIC() : CLASS_SIMWALK2() {
        _subname = "Nonparametric linkage";
      _suboption = 3;
    }
   ~CLASS_NONPARAMETRIC() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool is_sub_option()  { return true; }
//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual bool qtl_disallow()  { return true; }
    virtual bool simwalk2()      { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_IBD_EST: public CLASS_SIMWALK2 {
public:
    CLASS_IBD_EST() : CLASS_SIMWALK2() {
        _subname = "IBD estimation";
      _suboption = 4;
    }
   ~CLASS_IBD_EST() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool is_sub_option()  { return true; }
//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual bool qtl_disallow()  { return true; }
    virtual bool simwalk2()      { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MISTYPING: public CLASS_SIMWALK2 {
public:
    CLASS_MISTYPING() : CLASS_SIMWALK2() {
        _subname = "Mistyping detection";
      _suboption = 5;
    }
   ~CLASS_MISTYPING() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool is_sub_option()  { return true; }

//  virtual bool output_quant_can_define_missing_value() { return false; }

    virtual bool qtl_disallow()  { return true; }
    virtual bool simwalk2()      { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_HWETEST: public CLASS_ANALYSIS {
public:
    CLASS_HWETEST() : CLASS_ANALYSIS() {
        _name = "Hardy-Weinberg";
    }
   ~CLASS_HWETEST() {}

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool has_sub_options()    { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* for Mendel option only */
        OrigIds[0] = 5;
        OrigIds[1] = 3;
    }

    void sub_prog_name(int sub_opt, char *subprog);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_IQLS: public CLASS_ANALYSIS {
public:
    CLASS_IQLS() : CLASS_ANALYSIS() {
        _name = "IQLS/Idcoefs";
    }
   ~CLASS_IQLS() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool loops()  { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
    virtual bool allow_no_trait()  { return true; }
    // IQLS always writes one file per chromosome. see write_IQLS.cpp:create_IQLS_files()
    // So we can't allow the user to turn it off....
    //virtual bool Loop_Over_Chromosomes_implemented() { return false; }
    virtual bool qtl_disallow() { return true; }
    virtual bool require_physical_map()  { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Always unique ids for persons */
        /* Set the output pedigree and per field choices */
        OrigIds[0] = 3;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_LINKAGE: public CLASS_ANALYSIS {
public:
    CLASS_LINKAGE() : CLASS_ANALYSIS() {
        _name = "Linkage";
    }
   ~CLASS_LINKAGE() {}

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }

    // For quantitative variables a missing or unobserved value is coded as a 0.0
    // NOTE: considers an observation with absolute value less than 0.0000001 to indicate a missing value.
    virtual const char* output_quant_default_value() { return "0.0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow()    { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_LOD2: public CLASS_ANALYSIS {
public:
    CLASS_LOD2() : CLASS_ANALYSIS() {
        _name = "Homogeneity";
    }
   ~CLASS_LOD2() {}

    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "-99.99"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "-99"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_LOKI: public CLASS_ANALYSIS {
public:
    CLASS_LOKI() : CLASS_ANALYSIS() {
        _name = "Loki";
    }
   ~CLASS_LOKI() {}

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool forbid_sex_linked_loci()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "x"; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "x"; }

    virtual bool qtl_allow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
//      OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
        OrigIds[1] = 1;  //code expects this.
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MEGA2ANNOT: public CLASS_ANALYSIS {
public:
    CLASS_MEGA2ANNOT() : CLASS_ANALYSIS() {
        _name = "Mega2";
    }
   ~CLASS_MEGA2ANNOT() {}

    virtual bool allow_covariates()  { return true; }
    virtual bool allow_no_chr()  { return true; }
    virtual bool allow_no_map()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    virtual bool maintain_broken_loops()  { return true; }

    virtual const char* output_quant_default_value() { return "NA"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual const char* output_affect_default_value() { return "NA"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow()     { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Always unique ids for persons */
        /* Set the output pedigree and per field choices */
        OrigIds[0] = 3;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MENDEL: public CLASS_ANALYSIS {
public:
    CLASS_MENDEL() : CLASS_ANALYSIS() {
        _name = "Mendel";
    }
   ~CLASS_MENDEL() {}

    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual const char* output_quant_default_value() { return ""; }
    virtual const char* output_affect_default_value() { return ""; }

    virtual bool qtl_allow()     { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MENDEL4: public CLASS_ANALYSIS {
public:
    CLASS_MENDEL4() : CLASS_ANALYSIS() {
        _name = "Mendel7+";
    }
   ~CLASS_MENDEL4() {}

    virtual bool allow_sex_map() { return true; }

//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual const char* output_quant_default_value() { return ""; }
    virtual const char* output_affect_default_value() { return ""; }

    virtual bool qtl_allow()     { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MENDEL7_CSV: public CLASS_ANALYSIS {
public:
    CLASS_MENDEL7_CSV() : CLASS_ANALYSIS() {
        _name = "Mendel7+";
    }
   ~CLASS_MENDEL7_CSV() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }

//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual const char* output_quant_default_value() { return ""; }
    virtual const char* output_affect_default_value() { return ""; }

    virtual bool qtl_allow()     { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MERLIN: public CLASS_ANALYSIS {
public:
    CLASS_MERLIN() : CLASS_ANALYSIS() {
        _name = "Merlin/Simwalk2";
    }
   ~CLASS_MERLIN() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait()     { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool allow_sex_map() { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    // missing values are always replaced with an 'x'
//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual const char* output_quant_default_value() { return "x"; }
    virtual const char* output_affect_default_value() { return "x"; }

    virtual bool qtl_disallow()  { return true; }
    virtual bool simwalk2()      { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_MERLINONLY: public CLASS_ANALYSIS {
public:
    CLASS_MERLINONLY() : CLASS_ANALYSIS() {
        _name = "Merlin";
    }
   ~CLASS_MERLINONLY() {}

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }
    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_covariates()  { return true; }
    virtual bool allow_sex_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }

//  virtual bool output_quant_can_define_missing_value() { return false; }
    virtual const char* output_quant_default_value() { return "x"; }
    virtual const char* output_affect_default_value() { return "x"; }

    virtual bool qtl_allow()     { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_NUKE: public CLASS_ANALYSIS {
public:
    CLASS_NUKE() : CLASS_ANALYSIS() {
        _name = "Nuclear families";
    }
   ~CLASS_NUKE() {}

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool nuke_opts()    { return true; }

    // There is no shell script or control file created for this analysis
    // mode for which to output the quant missing value too...
    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow()    { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Person ids and consecutive numbering for pedigrees */
        OrigIds[0] = 5;
        OrigIds[1] = 3;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_PAP: public CLASS_ANALYSIS {
public:
    CLASS_PAP() : CLASS_ANALYSIS() {
        _name = "PAP";
    }
   ~CLASS_PAP() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }

    // http://hasstedt.genetics.utah.edu/pap/papmanual70.pdf
    // II.3. Marker/Trait Descriptions: header.dat
    // To indicate an unknown phenotype, the phenotype equals -9999 or the missing value code
    // specified for the variable in header.dat (which we no longer do)
//  virtual bool output_quant_can_define_missing_value() { return false; }

    virtual bool qtl_allow()    { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Always unique ids for persons */
        /* Set the output pedigree and per field choices */
        OrigIds[0] = 3;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_PREMAKEPED: public CLASS_ANALYSIS {
public:
    CLASS_PREMAKEPED() : CLASS_ANALYSIS() {
        _name = "Premakeped";
    }
   ~CLASS_PREMAKEPED() {}

    virtual bool allow_no_chr()     { return true; }
    virtual bool allow_trait_combination()  { return true; }

    // TO_MERLIN, TO_MERLINONLY, and TO_LOKI also use premakeped code
    // but use a default missing output quant.
    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_PREST: public CLASS_ANALYSIS {
public:
    CLASS_PREST() : CLASS_ANALYSIS() {
        _name = "Prest";
    }
   ~CLASS_PREST() {}

    virtual bool forbid_sex_linked_loci()  { return true; }
    virtual bool forbid_trait_directories()  { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SAGE: public CLASS_ANALYSIS {
public:
    CLASS_SAGE() : CLASS_ANALYSIS() {
        _name = "S.A.G.E.3.0";
    }
   ~CLASS_SAGE() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_trait_combination()  { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    // I have not been able to find the manuals for SAGE v3 online.
    // So, I have decided not to break things by trying to incorporate the
    // missing output quant batch file item in this version of SAGE.
//  virtual bool output_quant_can_define_missing_value() { return false; }

    virtual bool qtl_allow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SAGE4: public CLASS_ANALYSIS {
public:
    CLASS_SAGE4() : CLASS_ANALYSIS() {
        _name = "S.A.G.E.4.0";
    }
   ~CLASS_SAGE4() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_trait_combination()  { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "NA"; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "NA"; }

    virtual bool qtl_allow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SIMULATE: public CLASS_ANALYSIS {
public:
    CLASS_SIMULATE() : CLASS_ANALYSIS() {
        _name = "Simulate";
    }
   ~CLASS_SIMULATE() {}

    virtual bool allow_trait_combination()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};


class CLASS_SLINK: public CLASS_ANALYSIS {
public:
    CLASS_SLINK() : CLASS_ANALYSIS() {
        _name = "Slink";
    }
   ~CLASS_SLINK() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SOLAR: public CLASS_ANALYSIS {
public:
    CLASS_SOLAR() : CLASS_ANALYSIS() {
        _name = "SOLAR";
    }
   ~CLASS_SOLAR() {}

    virtual bool allow_covariates()  { return true; }
    virtual bool allow_trait_combination()  { return true; }

    virtual const char* output_quant_default_value() { return ""; }
    virtual const char* output_affect_default_value() { return ""; }

    virtual bool qtl_allow() { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SPLINK: public CLASS_ANALYSIS {
public:
    CLASS_SPLINK() : CLASS_ANALYSIS() {
        _name = "Splink";
    }
   ~CLASS_SPLINK() {}

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }

    // Accepts character alleles...
    virtual bool allele_data_use_name_if_available() { return true; }

    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }
    virtual bool nuke_opts()    { return true; }
    virtual bool qtl_disallow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Person ids and consecutive numbering for pedigrees */
        OrigIds[0] = 5;
        OrigIds[1] = 3;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SUMMARY: public CLASS_ANALYSIS {
public:
    CLASS_SUMMARY() : CLASS_ANALYSIS() {
        _name = "Summary";
    }
   ~CLASS_SUMMARY() {}

    virtual bool has_sub_options()    { return true; }

    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);
};

class CLASS_CREATE_SUMMARY: public CLASS_SUMMARY {
public:
    CLASS_CREATE_SUMMARY() : CLASS_SUMMARY() {
        _subname = "Segregation count";
      _suboption = 1;
    }
   ~CLASS_CREATE_SUMMARY() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool is_sub_option()  { return true; }
    virtual bool qtl_disallow() { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_ALLELE_FREQ: public CLASS_SUMMARY {
public:
    CLASS_ALLELE_FREQ() : CLASS_SUMMARY() {
        _subname = "Allele frequency";
      _suboption = 2;
    }
   ~CLASS_ALLELE_FREQ() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool is_sub_option()  { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_LIABLE_FREQ: public CLASS_SUMMARY {
public:
    CLASS_LIABLE_FREQ() : CLASS_SUMMARY() {
	   _subname = "Liability"; // "Counts within groups";
         _suboption = 3;
    }
   ~CLASS_LIABLE_FREQ() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool is_sub_option()  { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_GENOTYPING_SUMMARY: public CLASS_SUMMARY {
public:
    CLASS_GENOTYPING_SUMMARY() : CLASS_SUMMARY() {
        _subname = "Genotyping success summary";
      _suboption = 4;
    }
   ~CLASS_GENOTYPING_SUMMARY() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool is_sub_option()  { return true; }

    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
           set individual default to person field */
        OrigIds[0] = 1;
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_QUANT_SUMMARY: public CLASS_SUMMARY {
public:
    CLASS_QUANT_SUMMARY() : CLASS_SUMMARY() {
        _subname = "Quantitative Phenotypes summary";
      _suboption = 5;
    }
   ~CLASS_QUANT_SUMMARY() {}

    virtual bool has_sub_options()    { return false; }

    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool forbid_trait_directories()  { return true; }
    virtual bool is_sub_option()  { return true; }
    virtual bool qtl_allow()    { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_SUP: public CLASS_ANALYSIS {
public:
    CLASS_SUP() : CLASS_ANALYSIS() {
        _name = "SUP";
    }
   ~CLASS_SUP() {}

    virtual bool allow_no_trait()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }
    virtual bool qtl_allow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

// documentation shows this as disabled.
class CLASS_TDTMAX: public CLASS_ANALYSIS {
public:
    CLASS_TDTMAX() : CLASS_ANALYSIS() {
        _name = "TDT-max";
    }
   ~CLASS_TDTMAX() {}

    virtual bool is_enabled() { return false; }

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool maintain_broken_loops()  { return true; }
    virtual bool nuke_opts()    { return true; }
    virtual bool qtl_disallow() { return true; }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

class CLASS_VITESSE: public CLASS_ANALYSIS {
public:
    CLASS_VITESSE() : CLASS_ANALYSIS() {
        _name = "Vitesse";
    }
   ~CLASS_VITESSE() {}

    virtual bool allow_no_aff_trait() { return true; }
    virtual bool allow_no_trait()  { return true; }
    virtual bool break_loops()  { return true; }
    virtual bool has_sub_options()    { return true; }
    virtual bool maintain_broken_loops()  { return true; }

    virtual bool output_quant_can_define_missing_value() { return true; }
    virtual const char* output_quant_default_value() { return "0"; }
    virtual bool output_quant_must_be_numeric() { return true; }

    virtual bool output_affect_can_define_missing_value() { return true; }
    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual bool qtl_allow() { return true; }

    void sub_prog_name(int sub_opt, char *subprog);
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void gen_file_names(char *file_names[], char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

//extern void  prog_name_to_num(char *prog_name, analysis_type *analysis);

extern CLASS_HAPLOTYPE          *HAPLOTYPE;
extern CLASS_LOCATION           *LOCATION;
extern CLASS_NONPARAMETRIC      *NONPARAMETRIC;
extern CLASS_IBD_EST            *IBD_EST;
extern CLASS_MISTYPING          *MISTYPING;
extern CLASS_MENDEL             *TO_MENDEL;
extern CLASS_ASPEX              *TO_ASPEX;
extern CLASS_CREATE_SUMMARY     *CREATE_SUMMARY;
extern CLASS_GENOTYPING_SUMMARY *GENOTYPING_SUMMARY;
extern CLASS_GENEHUNTER         *TO_GeneHunter;
extern CLASS_APM                *TO_APM;
extern CLASS_APM_MULT           *TO_APM_MULT;
extern CLASS_NUKE               *TO_NUKE;
extern CLASS_LIABLE_FREQ        *TO_LIABLE_FREQ;
extern CLASS_ALLELE_FREQ        *TO_ALLELE_FREQ;
extern CLASS_SLINK              *TO_SLINK;
extern CLASS_SPLINK             *TO_SPLINK;
extern CLASS_LOD2               *TO_LOD2;
extern CLASS_GENEHUNTERPLUS     *TO_GeneHunterPlus;
extern CLASS_SIMULATE           *TO_SIMULATE;
extern CLASS_SAGE               *TO_SAGE;
extern CLASS_TDTMAX             *TO_TDTMAX;
extern CLASS_SOLAR              *TO_SOLAR;
extern CLASS_HWETEST            *TO_HWETEST;
extern CLASS_VITESSE            *TO_VITESSE;
extern CLASS_LINKAGE            *TO_LINKAGE;
extern CLASS_ALLEGRO            *TO_Allegro;
extern CLASS_GHMLB              *TO_GHMLB;
extern CLASS_SAGE4              *TO_SAGE4;
extern CLASS_PREMAKEPED         *TO_PREMAKEPED;
extern CLASS_MERLIN             *TO_MERLIN;
extern CLASS_PREST              *TO_PREST;
extern CLASS_PAP                *TO_PAP;
extern CLASS_MERLINONLY         *TO_MERLINONLY;
extern CLASS_QUANT_SUMMARY      *QUANT_SUMMARY;
extern CLASS_LOKI               *TO_LOKI;
extern CLASS_MENDEL4            *TO_MENDEL4;
extern CLASS_SUP                *TO_SUP;
extern CLASS_MENDEL7_CSV        *TO_MENDEL7_CSV;
extern CLASS_CRANEFOOT          *CRANEFOOT;
extern CLASS_MEGA2ANNOT         *MEGA2ANNOT;
extern CLASS_IQLS               *IQLS;

extern CLASS_SIMWALK2           *TO_SIMWALK2;
extern CLASS_SUMMARY            *TO_SUMMARY;

#include "write_plink_ext.h"
extern CLASS_PLINK              *TO_PLINK;

#endif
