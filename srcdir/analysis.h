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

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

class   CLASS_ANALYSIS;
typedef CLASS_ANALYSIS *analysis_type;

class CLASS_ANALYSIS {
public:
    const char  *_name;
    const char  *_subname;
    int          _suboption;
    char file_name_stem[100];

protected:
    file_format infl_type;
    file_format outfl_type;
    ped_top    *PedTreeTop;
//not yet
    char *mapfl_name;
public:
    CLASS_ANALYSIS() {

        _name      = "";
        _subname   = "";
        _suboption = -1;

        infl_type  = LINKAGE;
        outfl_type = UNKNOWN;
        PedTreeTop = NULL;
    }
    virtual ~CLASS_ANALYSIS() {}

/*
  If your option cannot handle affection status variables with
  liability classes, then Mega2 needs to be able to decide which
  class-status combinations are to be considered affected. So, add
  KEYWORD to the list of options inside the first case block of the
  function define_affection_labels().

  NOTE: The default behavior as defined here (return false) says that
  the analysis mode cannot handle affection status phenotypes deirectly.
  This can be a source of bugs in outputing affection status.
*/
    virtual bool allow_affection_liability_class()  { return false; }

/*
 * If QTLs can be specially designated as
 * covariates for analysis by your program, return true.
*/
    virtual bool allow_covariates()  { return false; }

/*
 */
    virtual bool allow_no_aff_trait()     { return false; }

/*
  Should the analysis option allow selection of markers whose
  chromosome number is unknown (special menu item in reorder_loci)?
 */
    virtual bool allow_no_chr()     { return false; }

/*
  pretty much all analysis requires a genetic map.
 */
    virtual bool allow_no_genetic_map()  { return false; }

    //
    // Used to allow the "you must have both parents defined" Linkage rule
    // to be violated by a specific Analysis type (see Beagle).
    // See makeped1.cpp:check_and_reassign_parents()
    virtual bool allow_missing_parent_in_linkage_input()  { return false; }

/*
   Should analysis option allow selection of markers
   whose map position is unknown (even though chromosome is known)?
*/
    virtual bool allow_no_map()     { return false; }

/*
  This function returns FALSE(0) if the analysis option requires at
  least one trait, TRUE (1) otherwise.
 */
    virtual bool allow_no_trait()  { return false; }

/*
  Should analysis option allow selection of markers on X-chromosome?
*/
    virtual bool allow_sex_map()    { return false; }

/*
 * This function return TRUE if traits are
 * allowed to be combined in the output as a single set of files, and
 * traits do not need to come first in the locus list. 
 */
    virtual bool allow_trait_combination()  { return false; }

/*
 * pre-makeped pedigrees can be converted to post-makeped,
 * if your analysis option requires this conversion.
 * 
 * Convert pre-makeped to post-makeped, breaking loops */
    virtual bool break_loops()  { return false; }

/*
*/
    virtual bool forbid_sex_linked_loci()  { return false; }

/*
*/
    virtual bool forbid_trait_directories()  { return false; }

    //
    // Returns true if the analysis is enabled, false if the analysis is DISABLED...
    virtual bool is_enabled() { return true; }

/*
  Whether this analysis consists of more than one
   type of analyses. This entails creating a menu to select
   sub-option, and using both analysis and sub-analysis to decide what
   output files should be created. Some output file names are
   sub-option specific, as are output file formats.
*/
    virtual bool has_sub_options()  { return false; }

/*

*/
    virtual bool is_sub_option()  { return false; }

/*
  analysis uses new loop classes
*/
    virtual bool loops()  { return false; }

/*
  True if the 'Loop_Over_Chromosomes' batch file item has been
  implemented for this analysis option.
*/
    virtual bool Loop_Over_Chromosomes_implemented() { return false; }

/*
 * If your option requires broken loops to be
 * maintained when post-makeped pedigrees are provided, then make the
 * function return true. Otherwise, post-makeped
 * pedigrees will have their loops identified and reconnected.
 *
 * Maintain broken-loops when post-makeped pedigrees are provided
 * ALL these analyses get the ID/Father/Mother in LinkPerID space
 * and must deal with it if necessary, viz. MEGA2ANNOT: (set_link_IDS)
 */
    virtual bool maintain_broken_loops()  { return false; }

/*
  whether pedigrees should be broken up into nuclear components.
*/
    virtual bool nuke_opts()        { return false; }

    //
    // NOTE: A distinction is made between the analysis types that take numeric and
    // non-numeric hard coded values. Those that use a letter (e.g., Merlin, Mendel7,
    // FBAT) are OK because the letter will not conflict with a numeric quantity.
    // Those that use a number (e.g., SUP, LINKAGE) generally have their code recompiled
    // to maket the default more useful (e.g., '0.0' used by linkage is a really bad
    // choice because it is likely to show up in data). So, for these analysis modes
    // we say that the user can define the output quant even though it is hard coded.

    //
    // Some analysis modes take a 'hard coded' default missing quantitative
    // value, and some allow the user to specify a value. For analysis that
    // take a hard coded value, warn the user if they give one beacuse it will
    // be ignored.
    virtual bool output_quant_can_define_missing_value() { return false; }
    virtual bool output_affect_can_define_missing_value() { return false; }
    //
    // This value is presented to the user as a default value in the "Missing output QTL
    // value definition menu" found in user_input.cpp:set_missing_quant_output(). If you
    // override 'output_quant_can_define_missing_value()' to return true then you should
    // also override 'output_quant_default_value()' to return a non-null character string.
    //
    // This value is also given to the batch file item 'Value_Missing_Quant_On_Output'
    // when 'output_quant_can_define_missing_value() == false'. In this manner
    // 'Value_Missing_Quant_On_Output' will seem to have "been read" as either what the user
    // set it to in the batch file if they were permitted too
    // (e.g., 'output_quant_can_define_missing_value() == true'), or the value of
    // 'output_quant_default_value()' if it is not NULL.
    //
    // HINT: In the default situation (e.g., 'output_quant_can_define_missing_value()' returns false,
    // and 'output_quant_default_value()' returns NULL), ITEM_READ(Value_Missing_Quant_On_Output)
    // will return false.
    virtual const char* output_quant_default_value() { return (const char *)NULL; }
    virtual const char* output_affect_default_value() { return (const char *)NULL; }
    //
    // This is used in batch_input.cpp:set_batch_items(), and user_input.cpp:set_missing_quant_output()
    // to determine if the string given for the output quantitative value must be numeric.
    //
    // HINT: If this is set as 'true' in an analysis sub-class then the value returned
    // by the method 'output_quant_default_value()' should be numeric as well...
    virtual bool output_quant_must_be_numeric() { return false; }
    virtual bool output_affect_must_be_numeric() { return false; }

    // True if the target analysis supports character alleles...
    // This implies that the named representation of an allele will be used if it is available.
    // If false alleles will always be written in their 'number'ed representation even if
    // there is a 'name'ed (character) representation available.
    // Use by write_files.cpp::write_numbered_data() and entry.cpp::person_locus_entry::pr_marker.
    virtual bool allele_data_use_name_if_available() { return false; }

/*
  whether QTLS are allowed.
*/
    virtual bool qtl_allow()        { return false; }

/*
  whether QTLS are allowed.
  This list is not the negation of the allow_qtl, since
  some of the options do not take trait loci at all.
*/
    virtual bool qtl_disallow()     { return false; }

/*
 */
    virtual bool require_physical_map()  { return false; }

/*
*/
    virtual bool simwalk2()         { return false; }

/*
  should we list/process specified trait number
*/
    virtual bool skip_trait(linkage_locus_top *LocusTop, int trait)  { return false; }

/*
  Add keyword to ped_ind_defaults(). This function decides which of
  the three individual ids, and two pedigree ids will be selected as
  the default output ids. This is critical if the target analysis has
  restrictions on ids (e.g., linkage-format only handles numerical values).

  Pedigree ID choices:
  1=Pedigree field for post-makeped format
  2=Pedigree field for pre-makeped
  3=Consecutive numerical ids for nuclear pedigrees

  Individual ID choices:
  1= Person field for post- or pre-makeped files
  2 = Unique Ids
  5 = Consecutively numbered person ids
  3,4 = 	not used, although these are handled.
*/
    virtual void ped_ind_defaults(int unique)  {
        /* options that require ids, mostly linkage options and
           gh-like options*/
        OrigIds[0] = 5; /* Renumbered */
        OrigIds[1] = 1; /* Ped num */
    }

/*
  copy out program name (_name) member, if set.
 */
    virtual void prog_name(char *prog) { if (*_name) strcpy(prog, _name); }

/*
  copy out sub program name (_subname) member, if set.
*/
    virtual void sub_prog_name(int sub_opt, char *subprog) {
        if (*_subname) strcpy(subprog, _subname); }

    // This method behaves in two different ways. In the older analysises where a sub-option
    // generated a new sub-class, this method will return that sub-class.
    // The new way of doing things is to have this method fill in the _suboption from
    // the analysis sub-option name 'sub_prog_name' as shown below.
    virtual void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
        (*analysis)->_suboption = 0;
        //analysis = sub_class_of_this_analysis; // a side effect used in the older code (see comments above)
    }

/*
   menu driven request for sub program name; and set id
*/
    virtual void interactive_sub_prog_name_to_sub_option(analysis_type *analysis) {  }

/*
   do the work to produce the output files for the analysis
 */
    virtual void create_output_file(linkage_ped_top *LPedTreeTop,
                                    CLASS_ANALYSIS **analysis,
                                    char *file_names[],
                                    int untyped_ped_opt,
                                    int *numchr,
                                    linkage_ped_top **Top2) {}
    // In the past, the method 'create_output_file' would create all of the output files.
    // With some widely used output formats (e.g., PLINK) it is usefull to build the
    // data files, and create shell files based on the way that the analysis uses the data files.
    // The method 'create_sh_files' is currently only used for PLINK and Eigensoft.
    // For completeness sake, it should be migrated to others.
    virtual void create_sh_file(linkage_ped_top *Top,
				char *file_names[],
				const int numchr) {}

/*
   initialize the file_names array with the names of outputs for this analysis
 */
    virtual void file_names(char **file_names, char *num) {}

/*
    change those file names that have the chromosome number embeded.
 */
    virtual void replace_chr_number(char *file_names[], int numchr) {}
};

struct analysis_types {
public:
    char       option_name[30];
    CLASS_ANALYSIS   *analysis;
};

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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }

    void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
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

    void file_names(char **file_names, char *num);
    void replace_chr_number(char *file_names[], int numchr);
};

#endif
