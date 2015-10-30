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

#ifndef PLINK_CORE_EXT_H
#define PLINK_CORE_EXT_H

#include "analysis.h"

//
// The names associated with these selections are defined in the subclass files
// (e.g., write_plink_ext.h, and write_eigenstat_ext.h).
// However, since the PLINK core code uses numbers they must be defned consistently here.
#define PLINK_SUB_OPTION_LGEN_INT                  1
#define PLINK_SUB_OPTION_SNP_MAJOR_INT             2
#define PLINK_SUB_OPTION_INDIVIDUAL_MAJOR_INT      3
#define PLINK_SUB_OPTION_PED_INT                   4


class   CLASS_PLINK_CORE;
typedef CLASS_PLINK_CORE *plink_analysis_type;


class CLASS_PLINK_CORE: public CLASS_ANALYSIS {

public:
    CLASS_PLINK_CORE() : CLASS_ANALYSIS() {
    }
    ~CLASS_PLINK_CORE() {}
    
    virtual bool allow_affection_liability_class()  { return true; }
    virtual bool allow_no_chr() { return true; }
    virtual bool allow_no_genetic_map()  { return true; }
    virtual bool allow_no_map() { return true; }
    virtual bool allow_trait_combination()  { return true; }
    virtual bool has_sub_options()    { return true; }
    virtual bool loops()  { return true; }
    virtual bool Loop_Over_Chromosomes_implemented() { return true; }
    virtual bool qtl_allow()    { return true; }
    virtual bool require_physical_map()  { return true; }
    
    virtual void ped_ind_defaults(int unique)  {
        /* Set the pedigree default to pedigree field
         set individual default to person field */
//      OrigIds[0] = 1; // OrigID
//      OrigIds[0] = 3; // UniqueID
        OrigIds[0] = 3; // OrigID 
        OrigIds[1] = (pedfile_type == POSTMAKEPED_PFT)? 1 : 2;
    }
    
    void replace_chr_number(char *file_names[], int numchr);
    
    //
    // These [pure virtual] funnctions must be implemented in every sub-class....
    
    void sub_prog_name(int sub_opt, char *subprog) = 0;
    void interactive_sub_prog_name_to_sub_option(analysis_type *analysis) = 0;
    void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) = 0;
    
    // NOTE: it is important for all
    //void sub_prog_name(int sub_opt, char *subprog);
    //void interactive_sub_prog_name_to_sub_option(analysis_type *analysis);
    //void sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis);
    virtual void create_output_file(linkage_ped_top *LPedTreeTop,
                            analysis_type *analysis,
                            char *file_names[],
                            int untyped_ped_opt,
                            int *numchr,
                            linkage_ped_top **Top2) = 0;
    virtual void save_pheno_file(linkage_ped_top *Top,
                                  const int pwid,
                                  const int fwid);
    virtual void save_pedsix_file(linkage_ped_top *Top,
                                  const int pwid,
                                  const int fwid) = 0;
    virtual void save_ped_file(linkage_ped_top *Top,
                               const int pwid,
                               const int fwid,
			       const int mwid) = 0;
    virtual void save_bed_file(const char *bedfl_name,
			       linkage_ped_top *Top,
			       const int binary_mode_flag) = 0;
    virtual void create_sh_file(linkage_ped_top *Top,
                        char *file_names[],
                        const int numchr) = 0;

    // Allows the stem to be changed as well. Usually called by 'file_names()' with
    // the default, or with a "custom stem" during user interactive mode.
    virtual void file_names_w_stem(char **file_names,
                                   char *num,
                                   const char *stem) = 0;
    // This is the method that must be defined for the CLASS_ANALYSIS.
    // It is just a stub that calls the real method with the default stem.
    void file_names(char **file_names, char *num) {
      file_names_w_stem(file_names, num, file_name_stem);
    }

    void user_queries(char **file_names, int *combine_chromo, int *create_summary, plink_analysis_type analysis);
};


class plink_binary {
public:
    int  SNP_count;      // always start on a new byte
    int  SNP_data;       // byte to write to the binary file
    int  SNP_bufsz;
    unsigned char *SNP_buf;
    unsigned char *SNP_cp;

    void file_header(FILE *filep) {
        // Write the magic numbers for a PLINK V1.00 .bed file
        fputc(0x6c, filep);
        fputc(0x1b, filep);
    }
    void inner (FILE *filep, int al1, int al2) {
        // 4 two bit nibbles are placed in a byte (8 bits),
        // right to left. So we need to left shift the nibble
        // before we or it in.
        if (al1 == 1 && al2 == 1) {
            // 00  Homozygote "1"/"1"
            SNP_data |= 0 << SNP_count;
        } else if (al1 != al2) {
            // 01  Heterozygote ( "1"/"2" or "2"/"1" )
            SNP_data |= 2 << SNP_count; // read backwards
        } else if (al1 == 2 && al2 == 2) {
            // 11  Homozygote "2"/"2"
            SNP_data |= 3 << SNP_count;
        } else {
            // 10  Missing genotype (in this case allele1 == 0)
            SNP_data |= 1 << SNP_count; // read backwards
        }
        SNP_count += 2;

        if ((SNP_count & 0x7) == 0) {
            *SNP_cp++ = (unsigned char) SNP_data;
            /*   fputc(SNP_data, filep) */;
            SNP_data = SNP_count = 0;
        }
    }
};

extern void create_PLINK_files(linkage_ped_top **LPedTop,
			       char *file_names[],
			       const int untyped_ped_opt,
			       const int output_format,
			       const char *prefix,
			       analysis_type *analysis);

extern void write_R_PLINK_map_file(linkage_ped_top *LPTop, char *map_file_name);

#endif /* PLINK_CORE_EXT_H */
