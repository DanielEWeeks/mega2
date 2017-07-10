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
                                     linkage_locus_top *LTop,
                                     std::string info_id_alternative_key,
                                     std::string unknown_marker_prefix,
                                     std::vector<Vecc> &VCFalleles);
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
