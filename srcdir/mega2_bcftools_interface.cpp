/*
Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2018 Robert Baron, Charles P. Kollar,
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

#include <string.h>
#include <vector>
#include <errno.h>
#include <iostream>
#include <fstream>

// Mega2 includes....
#include "common.h"
#include "typedefs.h"
#include "tod.hh"
#include "utils_ext.h"
#include "genetic_utils_ext.h"
#include "error_messages_ext.h"
#include "makeped.h"
#include "annotated_ped_file.h"
#include "linkage.h"
#include "read_files_ext.h"
#include "compress_ext.h"
#include "marker_lookup_ext.h"
#include "phe_lookup_ext.h"
#include "reorder_loci_ext.h"

// BCFTOOLS view includes..
extern "C" {
    #include "bcftools.h"
    #include "lib/bcftools-1.6/filter.h"
    #include "lib/bcftools-1.6/vcfview.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/khash_str2int.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/hts.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/vcf.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h"
    #include "lib/bcftools-1.6/htslib-1.6/htslib/vcfutils.h"
}
#include <getopt.h>

#define INTERNAL_MEGA2_VCFTOOLS_INTERFACE
#include "mega2_bcftools_interface.h"

// This is what Mega2 uses...
using namespace std;

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define ALLELE_NONREF 1
#define ALLELE_MINOR 2
#define ALLELE_ALT1 3
#define ALLELE_MAJOR 4
#define ALLELE_NONMAJOR 5

#define GT_NEED_HOM 1
#define GT_NEED_HET 2
#define GT_NO_HOM   3
#define GT_NO_HET   4
#define GT_NEED_MISSING 5
#define GT_NO_MISSING 6

int MEGA2_BCFTOOLS_INTERFACE::mega2_main_vcfview(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->min_ac = args->max_ac = args->min_af = args->max_af = -1;
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] =
            {
                    {"genotype",required_argument,NULL,'g'},
                    {"compression-level",required_argument,NULL,'l'},
                    {"threads",required_argument,NULL,9},
                    {"header-only",no_argument,NULL,'h'},
                    {"no-header",no_argument,NULL,'H'},
                    {"exclude",required_argument,NULL,'e'},
                    {"include",required_argument,NULL,'i'},
                    {"trim-alt-alleles",no_argument,NULL,'a'},
                    {"no-update",no_argument,NULL,'I'},
                    {"drop-genotypes",no_argument,NULL,'G'},
                    {"private",no_argument,NULL,'x'},
                    {"exclude-private",no_argument,NULL,'X'},
                    {"uncalled",no_argument,NULL,'u'},
                    {"exclude-uncalled",no_argument,NULL,'U'},
                    {"apply-filters",required_argument,NULL,'f'},
                    {"known",no_argument,NULL,'k'},
                    {"novel",no_argument,NULL,'n'},
                    {"min-alleles",required_argument,NULL,'m'},
                    {"max-alleles",required_argument,NULL,'M'},
                    {"samples",required_argument,NULL,'s'},
                    {"samples-file",required_argument,NULL,'S'},
                    {"force-samples",no_argument,NULL,1},
                    {"output-type",required_argument,NULL,'O'},
                    {"output-file",required_argument,NULL,'o'},
                    {"types",required_argument,NULL,'v'},
                    {"exclude-types",required_argument,NULL,'V'},
                    {"targets",required_argument,NULL,'t'},
                    {"targets-file",required_argument,NULL,'T'},
                    {"regions",required_argument,NULL,'r'},
                    {"regions-file",required_argument,NULL,'R'},
                    {"min-ac",required_argument,NULL,'c'},
                    {"max-ac",required_argument,NULL,'C'},
                    {"min-af",required_argument,NULL,'q'},
                    {"max-af",required_argument,NULL,'Q'},
                    {"phased",no_argument,NULL,'p'},
                    {"exclude-phased",no_argument,NULL,'P'},
                    {"no-version",no_argument,NULL,8},
                    {NULL,0,NULL,0}
            };
    char *tmp;
    while ((c = getopt_long(argc, argv, "l:t:T:r:R:o:O:s:S:Gf:knv:V:m:M:auUhHc:C:Ii:e:xXpPq:Q:g:",loptions,NULL)) >= 0)
    {
        char allele_type[8] = "nref";
        switch (c)
        {
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'l':
                args->clevel = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --compression-level %s\n", optarg);
                args->output_type |= FT_GZ;
                break;
            case 'o': args->fn_out = optarg; break;
            case 'H': args->print_header = 0; break;
            case 'h': args->header_only = 1; break;

            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;

            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case  1 : args->force_samples = 1; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;

            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm':
                args->min_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-alleles %s\n", optarg);
                break;
            case 'M':
                args->max_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --max-alleles %s\n", optarg);
                break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c':
            {
                args->min_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->min_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->min_ac)!=1 )
                    error("Error: Could not parse --min-ac %s\n", optarg);
                set_allele_type(&args->min_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'C':
            {
                args->max_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->max_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->max_ac)!=1 )
                    error("Error: Could not parse --max-ac %s\n", optarg);
                set_allele_type(&args->max_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'q':
            {
                args->min_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->min_af, allele_type)!=2 && sscanf(optarg,"%f",&args->min_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->min_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'Q':
            {
                args->max_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->max_af, allele_type)!=2 && sscanf(optarg,"%f",&args->max_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->max_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }

            case 'x': args->private_vars |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'X': args->private_vars |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'u': args->uncalled |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'U': args->uncalled |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'p': args->phased |= FLT_INCLUDE; break; // phased
            case 'P': args->phased |= FLT_EXCLUDE; break; // exclude-phased
            case 'g':
            {
                if ( !strcasecmp(optarg,"hom") ) args->gt_type = GT_NEED_HOM;
                else if ( !strcasecmp(optarg,"het") ) args->gt_type = GT_NEED_HET;
                else if ( !strcasecmp(optarg,"miss") ) args->gt_type = GT_NEED_MISSING;
                else if ( !strcasecmp(optarg,"^hom") ) args->gt_type = GT_NO_HOM;
                else if ( !strcasecmp(optarg,"^het") ) args->gt_type = GT_NO_HET;
                else if ( !strcasecmp(optarg,"^miss") ) args->gt_type = GT_NO_MISSING;
                else error("The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"%s\".\n", optarg);
                break;
            }
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( args->private_vars > FLT_EXCLUDE ) error("Only one of -x or -X can be given.\n");
    if ( args->uncalled > FLT_EXCLUDE ) error("Only one of -u or -U can be given.\n");
    if ( args->phased > FLT_EXCLUDE ) error("Only one of -p or -P can be given.\n");

    if ( args->sample_names && args->update_info) args->calc_ac = 1;

    char *fname = NULL;
    fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    else if ( optind+1 < argc )
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
            error("Failed to read the regions: %s\n", tmp.s);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }

    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

    init_data_vcfview(args);
    bcf_hdr_t *out_hdr = args->hnull ? args->hnull : (args->hsub ? args->hsub : args->hdr);
    if (args->print_header)
        bcf_hdr_write(args->out, out_hdr);
    else if ( args->output_type & FT_BCF )
        error("BCF output requires header, cannot proceed with -H\n");

    int ret = 0;
    if (!args->header_only)
    {
        while ( bcf_sr_next_line(args->files) )
        {
            bcf1_t *line = args->files->readers[0].buffer[0];
            if ( line->errcode && out_hdr!=args->hdr ) error("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
            if ( subset_vcf(args, line) ) {
                //bcf_unpack(line, BCF_UN_FMT);
                bcf_write1(args->out, out_hdr, line);
            }

        }
        ret = args->files->errnum;
        if ( ret ) fprintf(stderr,"Error: %s\n", bcf_sr_strerror(args->files->errnum));
    }

    //hts_close(args->out);
    //destroy_data_vcfview(args);
    //bcf_sr_destroy(args->files);
    //free(args);
    return ret;
}


//return an args object instead
args_t * MEGA2_BCFTOOLS_INTERFACE::get_args(int argc, char *argv[]){
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->min_ac = args->max_ac = args->min_af = args->max_af = -1;
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] =
            {
                    {"threads",required_argument,NULL,9},
                    {"header-only",no_argument,NULL,'h'},
                    {"no-header",no_argument,NULL,'H'},
                    {"exclude",required_argument,NULL,'e'},
                    {"include",required_argument,NULL,'i'},
                    {"drop-genotypes",no_argument,NULL,'G'},
                    {"uncalled",no_argument,NULL,'u'},
                    {"exclude-uncalled",no_argument,NULL,'U'},
                    {"apply-filters",required_argument,NULL,'f'},
                    {"known",no_argument,NULL,'k'},
                    {"novel",no_argument,NULL,'n'},
                    {"min-alleles",required_argument,NULL,'m'},
                    {"max-alleles",required_argument,NULL,'M'},
                    {"samples",required_argument,NULL,'s'},
                    {"samples-file",required_argument,NULL,'S'},
                    {"force-samples",no_argument,NULL,1},
                    {"types",required_argument,NULL,'v'},
                    {"exclude-types",required_argument,NULL,'V'},
                    {"targets",required_argument,NULL,'t'},
                    {"targets-file",required_argument,NULL,'T'},
                    {"regions",required_argument,NULL,'r'},
                    {"regions-file",required_argument,NULL,'R'},
                    {"min-ac",required_argument,NULL,'c'},
                    {"max-ac",required_argument,NULL,'C'},
                    {"min-af",required_argument,NULL,'q'},
                    {"max-af",required_argument,NULL,'Q'},
                    {"phased",no_argument,NULL,'p'},
                    {"exclude-phased",no_argument,NULL,'P'},
                    {NULL,0,NULL,0}
            };
    char *tmp;
    while ((c = getopt_long(argc, argv, "l:t:T:r:R:o:O:s:S:Gf:knv:V:m:M:auUhHc:C:Ii:e:xXpPq:Q:g:",loptions,NULL)) >= 0)
    {
        char allele_type[8] = "nref";
        switch (c)
        {
            case 'H': args->print_header = 0; break;
            case 'h': args->header_only = 1; break;

            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;

            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case  1 : args->force_samples = 1; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;

            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm':
                args->min_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-alleles %s\n", optarg);
                break;
            case 'M':
                args->max_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --max-alleles %s\n", optarg);
                break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c':
            {
                args->min_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->min_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->min_ac)!=1 )
                    error("Error: Could not parse --min-ac %s\n", optarg);
                set_allele_type(&args->min_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'C':
            {
                args->max_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->max_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->max_ac)!=1 )
                    error("Error: Could not parse --max-ac %s\n", optarg);
                set_allele_type(&args->max_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'q':
            {
                args->min_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->min_af, allele_type)!=2 && sscanf(optarg,"%f",&args->min_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->min_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'Q':
            {
                args->max_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->max_af, allele_type)!=2 && sscanf(optarg,"%f",&args->max_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->max_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }

            case 'u': args->uncalled |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'U': args->uncalled |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'p': args->phased |= FLT_INCLUDE; break; // phased
            case 'P': args->phased |= FLT_EXCLUDE; break; // exclude-phased
            case 'g':
            {
                if ( !strcasecmp(optarg,"hom") ) args->gt_type = GT_NEED_HOM;
                else if ( !strcasecmp(optarg,"het") ) args->gt_type = GT_NEED_HET;
                else if ( !strcasecmp(optarg,"miss") ) args->gt_type = GT_NEED_MISSING;
                else if ( !strcasecmp(optarg,"^hom") ) args->gt_type = GT_NO_HOM;
                else if ( !strcasecmp(optarg,"^het") ) args->gt_type = GT_NO_HET;
                else if ( !strcasecmp(optarg,"^miss") ) args->gt_type = GT_NO_MISSING;
                else error("The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"%s\".\n", optarg);
                break;
            }
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            default : error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( args->private_vars > FLT_EXCLUDE ) error("Only one of -x or -X can be given.\n");
    if ( args->uncalled > FLT_EXCLUDE ) error("Only one of -u or -U can be given.\n");
    if ( args->phased > FLT_EXCLUDE ) error("Only one of -p or -P can be given.\n");

    if ( args->sample_names && args->update_info) args->calc_ac = 1;

    char *fname = NULL;
    fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    else if ( optind+1 < argc )
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
            error("Failed to read the regions: %s\n", tmp.s);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }

    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

    init_data_vcfview(args);

    //need to reset optind, since bcftools doesn't
    optind = 1;
    return args;

}

int MEGA2_BCFTOOLS_INTERFACE::test_args(int argc, char *argv[]){
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->min_ac = args->max_ac = args->min_af = args->max_af = -1;
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] =
            {
                    {"threads",required_argument,NULL,9},
                    {"exclude",required_argument,NULL,'e'},
                    {"include",required_argument,NULL,'i'},
                    {"uncalled",no_argument,NULL,'u'},
                    {"exclude-uncalled",no_argument,NULL,'U'},
                    {"apply-filters",required_argument,NULL,'f'},
                    {"known",no_argument,NULL,'k'},
                    {"novel",no_argument,NULL,'n'},
                    {"min-alleles",required_argument,NULL,'m'},
                    {"max-alleles",required_argument,NULL,'M'},
                    {"samples",required_argument,NULL,'s'},
                    {"samples-file",required_argument,NULL,'S'},
                    {"force-samples",no_argument,NULL,1},
                    {"types",required_argument,NULL,'v'},
                    {"exclude-types",required_argument,NULL,'V'},
                    {"targets",required_argument,NULL,'t'},
                    {"targets-file",required_argument,NULL,'T'},
                    {"regions",required_argument,NULL,'r'},
                    {"regions-file",required_argument,NULL,'R'},
                    {"min-ac",required_argument,NULL,'c'},
                    {"max-ac",required_argument,NULL,'C'},
                    {"min-af",required_argument,NULL,'q'},
                    {"max-af",required_argument,NULL,'Q'},
                    {"phased",no_argument,NULL,'p'},
                    {"exclude-phased",no_argument,NULL,'P'},
                    {NULL,0,NULL,0}
            };
    int ret = 1;
    char *tmp;
    while ((c = getopt_long(argc, argv, "l:t:T:r:R:o:O:s:S:Gf:knv:V:m:M:auUhHc:C:Ii:e:xXpPq:Q:g:",loptions,NULL)) >= 0)
    {
        char allele_type[8] = "nref";
        switch (c)
        {
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;

            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case  1 : args->force_samples = 1; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;

            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm':
                args->min_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-alleles %s\n", optarg);
                break;
            case 'M':
                args->max_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --max-alleles %s\n", optarg);
                break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c':
            {
                args->min_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->min_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->min_ac)!=1 )
                    error("Error: Could not parse --min-ac %s\n", optarg);
                set_allele_type(&args->min_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'C':
            {
                args->max_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->max_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->max_ac)!=1 )
                    error("Error: Could not parse --max-ac %s\n", optarg);
                set_allele_type(&args->max_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'q':
            {
                args->min_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->min_af, allele_type)!=2 && sscanf(optarg,"%f",&args->min_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->min_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'Q':
            {
                args->max_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->max_af, allele_type)!=2 && sscanf(optarg,"%f",&args->max_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->max_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }

            case 'u': args->uncalled |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'U': args->uncalled |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'p': args->phased |= FLT_INCLUDE; break; // phased
            case 'P': args->phased |= FLT_EXCLUDE; break; // exclude-phased
            case 'g':
            {
                if ( !strcasecmp(optarg,"hom") ) args->gt_type = GT_NEED_HOM;
                else if ( !strcasecmp(optarg,"het") ) args->gt_type = GT_NEED_HET;
                else if ( !strcasecmp(optarg,"miss") ) args->gt_type = GT_NEED_MISSING;
                else if ( !strcasecmp(optarg,"^hom") ) args->gt_type = GT_NO_HOM;
                else if ( !strcasecmp(optarg,"^het") ) args->gt_type = GT_NO_HET;
                else if ( !strcasecmp(optarg,"^miss") ) args->gt_type = GT_NO_MISSING;
                else error("The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"%s\".\n", optarg);
                break;
            }
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            default: {ret = 0; break;}
        }
    }
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( args->private_vars > FLT_EXCLUDE ) error("Only one of -x or -X can be given.\n");
    if ( args->uncalled > FLT_EXCLUDE ) error("Only one of -u or -U can be given.\n");
    if ( args->phased > FLT_EXCLUDE ) error("Only one of -p or -P can be given.\n");

    if ( args->sample_names && args->update_info) args->calc_ac = 1;

    char *fname = NULL;
    fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    else if ( optind+1 < argc )
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
            error("Failed to read the regions: %s\n", tmp.s);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }

    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");

    //optind needs to be reset to 1 run bcftools again after testing that this is a correct string
    optind = 1;
    destroy_data_vcfview(args);
    bcf_sr_destroy(args->files);
    free(args);
    return ret;

}

