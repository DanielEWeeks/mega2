#ifndef __VCFVIEW_H__
#define __VCFVIEW_H__

#include <htslib/synced_bcf_reader.h>

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_list, *regions_list;
    int argc, clevel, n_threads, output_type, print_header, update_info, header_only, n_samples, *imap, calc_ac;
    int trim_alts, sites_only, known, novel, min_alleles, max_alleles, private_vars, uncalled, phased;
    int min_ac, min_ac_type, max_ac, max_ac_type, min_af_type, max_af_type, gt_type;
    int *ac, mac;
    float min_af, max_af;
    char *fn_ref, *fn_out, **samples;
    int sample_is_file, force_samples;
    char *include_types, *exclude_types;
    int include, exclude;
    int record_cmd_line;
    htsFile *out;
} args_t;

void init_data_vcfview(args_t *args);

void destroy_data_vcfview(args_t *args);

int bcf_all_phased(const bcf_hdr_t *header, bcf1_t *line);

int subset_vcf(args_t *args, bcf1_t *line);

void set_allele_type (int *atype, char *atype_string);

void usage(args_t *args);

#endif