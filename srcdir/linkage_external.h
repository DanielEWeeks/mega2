/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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


/* linkage_ext.h defines some extensions on the existing linkage data structure
 */


#ifndef LINKAGE_EXTERNAL_H
#define LINKAGE_EXTERNAL_H


/* This is a matrix of frequencies, where a non-zero
   frequency implies that a specific allele is present
   within a specific group.
   AlleleFreq[#groups][#alleles]
   This should not be very sparse in reality.
*/

typedef struct _ext_linkage_locus_rec {
    double *positions; /* one position per map */
    double *pos_male, *pos_female;
} ext_linkage_locus_rec;

typedef struct _ext_linkage_locus_top {
    int LocusCnt;
    int GroupCnt; /* Maximum number of groups in pedigree data */
    ext_linkage_locus_rec *EXLocus;

    char *map_functions; /* will be map_function[], one for each map */
    char **MapNames; /* Names of each map taken from Map file header */
    int MapCnt; /* For different locus maps */
//cpk    sex_map_types **SexMaps; /* record whether sex-maps were provided */
    int **SexMaps; /* record whether sex-maps were provided */
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    int **valid_map_p; // a predicate that tells if the map is valid
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
} ext_linkage_locus_top;

//
// The classes 'm2_map' and 'm2_map_entry' define a simple three column map (for now).
// In the future this can be customized/extended to add other types of maps (e.g., plink).

class m2_map_entry {
private:
    int chr;
    double POS;
    std::string marker_name;
    std::string REF; // Reference Allele
    
public:
    m2_map_entry() : chr(UNKNOWN_CHROMO), POS(0), marker_name(""), REF("") {};
    m2_map_entry(std::string CHROM, double POS, std::string marker_name, std::string REF) {
        set_chr(CHROM); set_POS(POS); set_marker_name(marker_name); set_REF(REF);
    };
    m2_map_entry(int chr, double POS, std::string marker_name, std::string REF) {
        this->chr = chr; this->POS = POS; this->marker_name = marker_name; this->REF = REF;
    };
    ~m2_map_entry() {};
    
    const int get_chr() { return chr; };
    const linkage_locus_type get_locus_type (); // e.g. XLINKED, YLINKED, or NUMBERED (see linkage.h)
    const std::string get_chr_type_string (); // e.g. 'X', 'Y', or, 'M'
    // If the CHROM is specified as a string we must convert it to
    // an integer, because that is what Mega2 want's to process it as. Normal CHROM strings
    // are things like '1' 'chr1'. If it's something like 'X', 'Y', or 'chrX' or 'chrY' then
    // we need to convert it to a number based on the organism.
    void set_chr(const std::string CHROM);
    
    const double get_POS() { return POS; };
    void set_POS(const double POSx) { this->POS = POSx; };
    
    const std::string get_marker_name() { return marker_name; };
    void set_marker_name(const std::string marker_nameX) { this->marker_name = marker_nameX; };
    
    const std::string get_REF() { return REF; };
    void set_REF(const std::string REFx) { this->REF = REFx; };
};

class m2_map {
private:
    std::string full_name; // name.function
    std::vector <m2_map_entry> entries;
    
public:
    // where function is: 'h', 'k', or 'p'
    m2_map(const std::string name, const char function);
    m2_map() { full_name = "invalid.x"; }
   ~m2_map() {};
    
    const std::string get_name() { return full_name; };
    
    const char get_function() { return full_name.at(full_name.length() -1 ); };
    
    const m2_map_entry get_entry(const unsigned int i) { return entries[i]; };
    void push_back_entry(m2_map_entry e) { entries.push_back(e); };

    const size_t size() { return entries.size(); };

    void gc();
};

#endif
