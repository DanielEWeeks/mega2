/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 2012-2017, University of Pittsburgh. All Rights Reserved.

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
    const char  *_missing_value_key;
    int          _suboption;
    char        *file_name_stem;

protected:
    file_format infl_type;
    file_format outfl_type;
    ped_top    *PedTreeTop;
//not yet
    char *mapfl_name;
public:
    CLASS_ANALYSIS() : file_name_stem(NULL) {

        _name      = "";
        _subname   = "";
        _missing_value_key = "";
        _suboption = -1;

        infl_type  = LINKAGE;
        outfl_type = UNKNOWN;
        PedTreeTop = NULL;
        mapfl_name = NULL;
    }
    virtual ~CLASS_ANALYSIS() {
        if (file_name_stem) {
            free(file_name_stem);
            file_name_stem = NULL;
        }
    }

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
//  virtual bool allow_no_aff_trait()     { return false; }
    virtual bool requires_aff_trait()     { return false; }

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
//  virtual bool allow_no_trait()  { return false; }
    virtual bool require_traits()  { return false; }

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
     * New function to hide nonautosomal Chromosomes.
     * The Reorder Loci Menu will only display Chr 1-22 if set to true.
     */
    virtual bool only_display_autosomes() { return false;}

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

/*
 * Do not request the missing menu for Mega2.
 */
    virtual bool no_missing_menu() { return false; }
/*
 * Do not request the trait and covariate menu for Mega2.
 */
    virtual bool no_trait_covariate_menu() { return false; }



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
 */
    virtual void user_queries(char **file_names_array,
                              int *combine_chromo, int *create_summary) {}

/*
 */
    virtual void batch_in() {}

/*
 */
    virtual void batch_out() {}

/*
 */
    virtual void batch_show() {}

/*
   initialize the file_names array with the names of outputs for this analysis
 */
    virtual void gen_file_names(char *file_names[], char *num) {}

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

class Missing_Value {
/* some day
    enum class numeric {OFF=0, Num,   Any};
   but scoped enum are a c++11 feature and possibly ... not there
*/

public:
    enum numeric {Numeric=0, Num=1,   Any=2};
    enum change  {Change=0,  Fixed=1, Varies=2, NoQnt=3};

public:
    std::string       Description;
    CLASS_ANALYSIS   *analysis;

    const char      *quant_str;
    enum numeric     quant_num;
    enum change      quant_change;

    const char      *affect_str;
    enum numeric     affect_num;
    enum change      affect_change;

    const char      *allele_str;
    enum numeric     allele_num;
    enum change      allele_change;
};

#endif
