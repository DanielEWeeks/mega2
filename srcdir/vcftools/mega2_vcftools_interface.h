/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

===========================================================================
*/

/*
  Mega2 interface to VCFtools.

  VCFtools is a program package used for processing VCF files. Further information
  on VCFtools can be found at the URL...
  http://vcftools.sourceforge.net
  We would like to thank those who developed this software and we include it in Mega2
  as a part of our continuing service to the Genetics community.
*/

#ifndef MEGA2_VCFTOOLS_INTERFACE_H
#define MEGA2_VCFTOOLS_INTERFACE_H

#include <string>
#include <vector>

#include "../annotated_ped_file.h"

//
// The classes 'm2_map' and 'm2_map_entry' define a simple three column map (for now).
// In the future this can be customized/extended to add other types of maps (e.g., plink).

class m2_map_entry {
private:
    int chr;
    double POS;
    string marker_name;
    string REF; // Reference Allele
    
public:
    m2_map_entry() : chr(UNKNOWN_CHROMO), POS(0), marker_name(""), REF("") {};
    m2_map_entry(string CHROM, double POS, string marker_name, string REF) {
        set_chr(CHROM); set_POS(POS); set_marker_name(marker_name); set_REF(REF);
    };
    m2_map_entry(int chr, double POS, string marker_name, string REF) {
        this->chr = chr; this->POS = POS; this->marker_name = marker_name; this->REF = REF;
    };
    ~m2_map_entry() {};
    
    const int get_chr() { return chr; };
    const linkage_locus_type get_locus_type (); // e.g. XLINKED, YLINKED, or NUMBERED (see linkage.h)
    const string get_chr_type_string (); // e.g. 'X', 'Y', or, 'M'
    // If the CHROM is specified as a string we must convert it to
    // an integer, because that is what Mega2 want's to process it as. Normal CHROM strings
    // are things like '1' 'chr1'. If it's something like 'X', 'Y', or 'chrX' or 'chrY' then
    // we need to convert it to a number based on the organism.
    void set_chr(const string CHROM);
    
    const double get_POS() { return POS; };
    void set_POS(const double POS) { this->POS = POS; };
    
    const string get_marker_name() { return marker_name; };
    void set_marker_name(const string marker_name) { this->marker_name = marker_name; };
    
    const string get_REF() { return REF; };
    void set_REF(const string REF) { this->REF = REF; };
};

class m2_map {
private:
    string full_name; // name.function
    vector <m2_map_entry> entries;
    
public:
    // where function is: 'h', 'k', or 'p'
    m2_map(const string name, const char function);
    m2_map() { full_name = "invalid.x"; }
    ~m2_map() {};
    
    const string get_name() { return full_name; };
    
    const char get_function() { return full_name.at(full_name.length() -1 ); };
    
    const m2_map_entry get_entry(const unsigned int i) { return entries[i]; };
    void push_back_entry(m2_map_entry e) { entries.push_back(e); };
    
    const size_t size() { return entries.size(); };
};


#ifndef INTERNAL_MEGA2_VCFTOOLS_INTERFACE

extern void VCFtools_printf_supported_cmd_line_options();
extern int VCFtools_process_cmd_line_w_file(const char *VCFArgs);
extern int VCFtools_process_cmd_line_wo_file(const char *VCFArgs);
extern int VCFtools_process_cmd_line_if_necessary_inclusive(const char *VCFArgs);

extern void VCFtools_process_file_meta_information_and_header();

//
// Both 'VCFtools_process_entries()' and 'VCFtools_get_map()' may be exected after
// 'VCFtools_process_file_meta_information_and_header'.
extern void VCFtools_process_entries(annotated_ped_rec persons[],
                                     const unsigned int persons_n,
                                     linkage_locus_top *LTop);
//
// Here we strip the .MAP file from the side of the VCF file for those entries
// that have passed the VCFtools filtering criteria. The 'file' is created as a
// vector of vcf_map_entry(s). The entries in the vector are read as they appear
// in the VCF file. No tests are made to see that the CHROM and POS within them
// are increasing as is stated in the VCF File standard.
// http://www.1000genomes.org/node/101
extern m2_map VCFtools_get_map(const std::string info_id_alternative_key,
                               const std::string unknown_marker_prefix);

// Can be called at any time.
extern void VCFtools_close();

#endif /* INTERNAL_MEGA2_VCFTOOLS_INTERFACE */


#endif /* MEGA2_VCFTOOLS_INTERFACE_H */
