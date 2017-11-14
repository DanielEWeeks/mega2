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
  Mega2 interface to BCFTOOLS.

*/

#ifndef MEGA2_BCFTOOLS_INTERFACE_H
#define MEGA2_BCFTOOLS_INTERFACE_H

#include <string>
#include <vector>
#include "lib/bcftools-1.6/filter.h"
#include <lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h>

#ifndef INTERNAL_MEGA2_BCFTOOLS_INTERFACE


//typedef struct _args_t
//{
//    filter_t *filter;
//    char *filter_str;
//    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
//
//    bcf_srs_t *files;
//    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
//    char **argv, *format, *sample_names, *subset_fname, *targets_list, *regions_list;
//    int argc, clevel, n_threads, output_type, print_header, update_info, header_only, n_samples, *imap, calc_ac;
//    int trim_alts, sites_only, known, novel, min_alleles, max_alleles, private_vars, uncalled, phased;
//    int min_ac, min_ac_type, max_ac, max_ac_type, min_af_type, max_af_type, gt_type;
//    int *ac, mac;
//    float min_af, max_af;
//    char *fn_ref, *fn_out, **samples;
//    int sample_is_file, force_samples;
//    char *include_types, *exclude_types;
//    int include, exclude;
//    int record_cmd_line;
//    htsFile *out;
//}
//        args_t;


class MEGA2_BCFTOOLS_INTERFACE {

public:
    MEGA2_BCFTOOLS_INTERFACE( ) {

    }

    ~MEGA2_BCFTOOLS_INTERFACE() { }

    int mega2_main_vcfview(int argc, char *argv[]);


    //void init_data(args_t *args);

   // void destroy_data(args_t *args);

    //int bcf_all_phased(const bcf_hdr_t *header, bcf1_t *line);

    //int subset_vcf(args_t *args, bcf1_t *line);

    //void set_allele_type (int *atype, char *atype_string);

    //void usage(args_t *args);

};

extern MEGA2_BCFTOOLS_INTERFACE *mega2_bcftools_interface;


#endif /* INTERNAL_MEGA2_BCFTOOLS_INTERFACE */


#endif /* MEGA2_BCFTOOLS_INTERFACE_H */

