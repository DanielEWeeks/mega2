/*
Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

  Our goal in integrating VCftools with Mega2 is to make minimal changes to VCFtools.

  NOTE: The following changes were made to the vcftools source. Additional changes
  are noted below.
  variant_file.h: comment out //#include "gamma.h"
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

// VCFtools includes...
// NOTE: that vcftools was modified to use the 'vcftools' namespace because
// of symbol clashes with Mega2 (in particular 'entity').
#include "output_log.h"
#include "parameters.h"
#include "variant_file.h"
#include "bcf_file.h"
#include "vcf_file.h"
#define INTERNAL_MEGA2_VCFTOOLS_INTERFACE
#include "mega2_vcftools_interface.h"

// This is what Mega2 uses...
using namespace std;

////////////////////////////////////////////////////////////////////////////////

//
// The 'output_log' class is used everywhere in VCFtools.
// Here we graft on the Mega2 logging mechanism to the methods used by output_log.
// These methods take the place of anything that was in the file 'output_log.cpp' from the VCFtools distribution.

output_log LOG;

// warning: private field 'output_to_screen' is not used
//output_log::output_log() : output_to_screen(true) { }
output_log::output_log() { }

void output_log::open(const string &filename_prefix ) { }

void output_log::close() { }

void output_log::set_screen_output(bool do_screen_output) { }

void output_log::printLOG(string s) {
    mssgf(s.c_str());
}

void output_log::error(string err_msg, int error_code) {
    errorf(err_msg.c_str());
    EXIT(error_code); // We need to factor this into a Mega2 error code somehow
}

void output_log::error(string err_msg, double value1, double value2, int error_code) {
    stringstream ss;
    ss << "Value1=" << value1 << " Value2=" << value2 << endl;
    errorf(ss.str().c_str());
    EXIT(error_code);
}

void output_log::warning(string err_msg) {
    warnf(err_msg.c_str());
}

void output_log::one_off_warning(string err_msg) {
    static set<string> previous_warnings;
    if (previous_warnings.find(err_msg) == previous_warnings.end())
    {
        // if the warning was not found before....
        warnf(err_msg.c_str());
        previous_warnings.insert(err_msg);
    }
}

string output_log::int2str(int n) {
    std::ostringstream s2( std::stringstream::out );
    s2 << n;
    return s2.str();
}

string output_log::longint2str(long int n) {
    std::ostringstream s2( std::stringstream::out );
    s2 << n;
    return s2.str();
}

string output_log::dbl2str(double n, int prc) {
    std::ostringstream s2;
    if ( prc > 0 )
        s2.precision(prc);
    s2 << n;
    return s2.str();
}

string output_log::dbl2str_fixed(double n, int prc) {
    std::ostringstream s2;
    s2 << setiosflags( ios::fixed );
    if ( prc > 0 )
        s2.precision(prc);
    s2 << n;
    return s2.str();
}

////////////////////////////////////////////////////////////////////////////////



//
// vcftools: Usage and Options
// http://vcftools.sourceforge.net/options.html#basic

//
// Here we define the data structures for parsing the VCFtools command line.
// We do this because we want to do a better job at telling the user what went
// wrong with the options string, we also need to list and accept what for us are
// valid options 'supported_vcftools_options[]' since we only process the
// options related to filtering.

enum argType {
    NONE,
    CONSTANT,
    FILENAME
};

typedef struct _command_line_option {
#ifdef _WIN
    char *str;       // The option string e.g., '--foo'
    enum argType at; // e.g., 0: '--foo'; 1: '--foo value'; 2: -foo filename
#else
    const char *str;       // The option string e.g., '--foo'
    const enum argType at; // e.g., 0: '--foo'; 1: '--foo value'; 2: -foo filename
#endif
} command_line_option;

//
// Mega2 only allows for the use of the VCFtools 'Basic Options', and 'Filters'.

static command_line_option supported_vcftools_options[] = {
    // Basic Options...
    // The three file name options are inserted into the option line by Mega2.
    {"--vcf", FILENAME},
    {"--gzvcf", FILENAME},
    {"--bcf", FILENAME},
    // Site Filter Options...
    {"--chr", CONSTANT}, // <chromosome>
    {"--not-chr", CONSTANT}, // <chromosome>
    {"--from-bp", CONSTANT}, // <integer>
    {"--to-bp", CONSTANT}, // <integer>
    {"--snp", CONSTANT}, // <string>
    {"--snps", FILENAME},
    {"--exclude", FILENAME},
    {"--positions", FILENAME},
    {"--exclude-positions", FILENAME},
    {"--keep-only-indels", NONE},
    {"--remove-indels", NONE},
    {"--bed", FILENAME},
    {"--exclude-bed", FILENAME},
    {"--remove-filtered-all", NONE},
    {"--remove-filtered", CONSTANT}, // <string>
    {"--keep-filtered", CONSTANT}, // <string>
    {"--remove-INFO", CONSTANT}, // <string>
    {"--keep-INFO", CONSTANT}, // <string>
    {"--minQ", CONSTANT}, // <float>
    {"--min-meanDP", CONSTANT}, // <float>
    {"--max-meanDP", CONSTANT}, // <float>
    {"--maf", CONSTANT}, // <float>
    {"--max-maf", CONSTANT}, // <float>
    {"--non-ref_af", CONSTANT}, // <float>
    {"--max-non-ref-af", CONSTANT}, // <float>
    {"--mac", CONSTANT}, // <int>
    {"--max-mac", CONSTANT}, // <int>
    {"--non-ref-ac", CONSTANT}, // <float>
    {"--max-non-ref-ac", CONSTANT}, // <float>
    {"--hwe", CONSTANT}, // <float>
    {"--geno", CONSTANT}, // <float>
    {"--max-missing-count", CONSTANT}, // <int>
    {"--min-alleles", CONSTANT}, // <int>
    {"--max-alleles", CONSTANT}, // <int>
    {"--thin", CONSTANT}, // <int>
    {"--mask", FILENAME},
    {"--invert-mask", FILENAME},
    {"--mask-min", CONSTANT}, // <int>
    // Individual Filters...
    {"--indv", CONSTANT}, // <string>
    {"--keep", FILENAME},
    {"--remove-indv", CONSTANT},
    {"--remove", FILENAME},
    {"--min-indv-meanDP", CONSTANT}, // <float>
    {"--max-indv-meanDP", CONSTANT}, // <float>
    {"--mind", CONSTANT}, // <float>
    {"--phased", NONE},
    {"--max-indv", CONSTANT}, // <int>
    // Genotype Filters...
    {"--remove-filtered-geno-all", NONE},
    {"--remove-filtered-geno", CONSTANT}, // <string>
    {"--minGQ", CONSTANT}, // <float>
    {"--minDP", CONSTANT}, // <float>
    {"--maxDP", CONSTANT} // <float>
};
#define SUPPORTED_VCFTOOLS_OPTIONS_COUNT (int)(sizeof(supported_vcftools_options) / sizeof(command_line_option))

//
// Mega2 only supports the filtering command line options for now.
// Anything in this list will not be supported by Mega2.
// We have this list so that we can recognise the string as a valid vcftools
// command line option as opposed to something totally unsupported by vcftools.

static const command_line_option unsupported_vcftools_options[] = {
    // Basic Options...
    {"--out", CONSTANT}, // <string> file prefix for generating output.
    // Output Statistics...
    {"--freq", NONE},
    {"--counts", NONE},
    {"--freq2", NONE},
    {"--counts2", NONE},
    {"--depth", NONE},
    {"--site-depth", NONE},
    {"--site-mean-depth", NONE},
    {"--geno-depth", NONE},
    {"--site-quality", NONE},
    {"--het", NONE},
    {"--hardy", NONE},
    {"--missing", NONE},
    {"--hap-r2", NONE},
    {"--geno-r2", NONE},
    {"--geno-chisq", NONE},
    {"--ld-window", CONSTANT},
    {"--ld-window-bp", CONSTANT},
    {"--min-r2", CONSTANT},
    {"--SNPdensity", CONSTANT},
    {"--TsTv", CONSTANT},
    {"--TsTv-by-count", NONE},
    {"--TsTv-by-qual", NONE},
    {"--FILTER-summary", NONE},
    {"--filtered-sites", NONE},
    {"--singletons", NONE},
    {"--site-pl", NONE},
    {"--window-pl", CONSTANT},
    {"--window-pl-step", CONSTANT},
    {"--hist-indel-len", NONE},
    {"--TajimaD", CONSTANT},
    {"--hapmap-fst-pop", FILENAME},
    {"--weir-fst-pop", FILENAME},
    {"--fst-window-size", CONSTANT},
    {"--fst-window-step", CONSTANT},
    // Output in Other Formats...
    {"--012", NONE},
    {"--IMPUTE", NONE},
    {"--ldhat", NONE},
    {"--ldhat-geno", NONE},
    {"--BEAGLE-GL", NONE},
    {"--BEAGLE-PL", NONE},
    {"--plink", NONE},
    {"--plink-tped", NONE},
    {"--recode", NONE},
    {"--recode-bcf", NONE},
    {"--recode-to-stream", NONE},
    {"--recode-bcf-to-stream", NONE},
    // Miscellaneous...
    {"--extract-FORMAT-info", CONSTANT},
    {"--get-INFO", CONSTANT},
    {"--force-index-write", NONE},
    // VCF File Comparison Options...
    {"--diff", FILENAME},
    {"--gzdiff", FILENAME},
    {"--diff-bcf", FILENAME},
    {"--diff-site-discordance", NONE},
    {"--diff-indv-discordance", NONE},
    {"--diff-indv-map", FILENAME},
    {"--diff-discordance-matrix", NONE},
    {"--diff-switch-error", NONE},
    // Options still in development...
    {"--LROH", NONE},
    {"--relalatedness", NONE}
};
#define UNSUPPORTED_VCFTOOLS_OPTIONS_COUNT (int)(sizeof(unsupported_vcftools_options) / sizeof(command_line_option))

//
// Print for the user a list of the command line options that are supported by Mega2
// in such a manner that the line length never exceeds 'max_len'.
void VCFtools_printf_supported_cmd_line_options()
{
    static const int max_len = 80;
    int j, i, len;
    static const char *arg = "v";
    
    // 'j=3' to skip the first three options which can only be filled in by Mega2,
    // never the user...
    // 'i' is a running count of the leangth of the line...
    for (j=3,i=0; j<SUPPORTED_VCFTOOLS_OPTIONS_COUNT; j++) {
        // NOTE: if the length of any of these parameters are ever >= max_len
        // then this will result in an infinite loop...
        len = (int)strlen(supported_vcftools_options[j].str) + 1;
        len += (supported_vcftools_options[j].at != NONE ? strlen(arg)+1 : 0);
        if ((i += len) >= max_len) {
            printf("\n"); // line is never >= max_len characters
            i = len;
        }
        printf("%s ", supported_vcftools_options[j].str);
        if (supported_vcftools_options[j].at != NONE) printf("%s ", arg);
    }
    printf("\n");
}

//
// VCFtools class instances that are used by Mega2...
static parameters *params = (parameters *)NULL;
static variant_file *vf = (variant_file *)NULL;
// Each data line in the VCF file will be parsed into an 'entry', but since we
// don't refer back to them we use the same one for all lines in the file.
static entry *e = (entry *)NULL;

//
// Close the file and destroy the objects used to process the file.
// Keep 'params' around in case we need to process the file again in the same way.
// NOTE: this does not delete the m2_map this is created from the file.
void VCFtools_close() {
    delete params;
    params = (parameters *)NULL;
    delete vf;
    vf = (variant_file *)NULL;
    delete e;
    e = (entry *)NULL;
}

//
// These functions take the place of vcftools.cpp:main() and serves as the initialization of VCFtools for Mega2.

/**
 This function is called in user_input.cpp:menu1(). It will create and fill 'params' which
 is used for further processing from the 'VCFArgs' command line arguments given. Only VCFtools
 command line options supported by Mega2 are permitted. It tries to provide some additional
 checking of the arguments that VCFtools does not.
 
 @return -1 for OK, otherwise a VCFtools error code.
 */
static int VCFtools_process_cmd_line(const char *VCFArgs,
                                     const int supported_vcftools_options_i)
{
    int argc = 1, i, j, found;
    char **argv = (char **)malloc(sizeof(char *));
    char *cmd_line = strdup(VCFArgs);
    
    // Strip any newline character from the end of the string...
    if (cmd_line[strlen(cmd_line)-1] == '\n') cmd_line[strlen(cmd_line)-1] = '\0';
    if (cmd_line[strlen(cmd_line)-1] == '\r') cmd_line[strlen(cmd_line)-1] = '\0';
    
    // The VCFtools code assumes the name of the program is the first argument on the
    // command line, as it would be if the program were run from a shell, so put it there...
    argv[0] = strdup("VCFtools");
    
    // Break the VCFtools command line into argc, and argv used by the parameters constructor...
    // NOTE: Potential leak of memory pointed to by cmd_line...
    char *p = strtok((char *)cmd_line, " ");
    while (p) {
        if ((argv = (char **)realloc(argv, sizeof(char *) * (size_t)++argc)) == (char **)NULL)  {
            errorf("Could not allocate enough memory, exiting.");
            EXIT(MEMORY_ALLOC_ERROR);
        }
        argv[argc-1] = p;
        p = strtok (NULL, " ");
    }
    
    // Search the command line (following argv[0]) for supported options...
    for (i=1; i<argc; i++)
        // All command line options begin with '--'...
        if (strncmp(argv[i], "--", 2) == 0) {
            // Try to locate the option in the list that we support...
            for (j=supported_vcftools_options_i, found=0; j<SUPPORTED_VCFTOOLS_OPTIONS_COUNT; j++) {
                if (strcasecmp(argv[i], supported_vcftools_options[j].str) == 0) {
                    found = 1;
                    break;
                }
            }
            if (found == 1) {
                // Good, the option found!
                // If the option does not have an argument, don't check it...
                if (supported_vcftools_options[j].at == NONE) continue;
                
                // Check for the existence of what appears to be a valid argument.
                // Valid is define as: not being at the end of the argument list, or
                // not beginning with '--'.
                if (i+1 >= argc || strncmp(argv[i+1], "--", 2) == 0) {
                    errorvf("The argument required for option '%s' appears to be missing.\n",
                            supported_vcftools_options[j].str);
                    free(argv[0]); free(argv); free(cmd_line);
                    return 0; // an error occurred...
                }
                
                // Depending on the argument type, do some checking...
                if (supported_vcftools_options[j].at == FILENAME) {
                    FILE *file;
                    // This is more portable than using POSIX 'stat()'...
                    if ((file = fopen(argv[i+1], "rb")) == (FILE *)NULL) {
                        // Tell the user about the error...
                        if (errno == ENOENT) {
                            errorvf("The file '%s' associated with option '%s' does not exist.\n",
                                    argv[i+1], supported_vcftools_options[j].str);
                        } else {
                            errorvf("Unknown error accessing file '%s' associated with option '%s'.\n",
                                    argv[i+1], supported_vcftools_options[j].str);
                        }
                        free(argv[0]); free(argv); free(cmd_line);
                        return 0; // an error occurred...
                    } else {
                        // Found it, but it could be a directory as well as a file.
                        // We don't check further to see if it is a directory.
                        fclose (file); // no error occurred, keep processing...
                    }
                }
            } else { // if (found == 1) {
                     // Here the command line option was not found so we are just trying to
                     // figure out what type of error message to give the user...
                for (j=0, found=0; j<UNSUPPORTED_VCFTOOLS_OPTIONS_COUNT; j++) {
                    if (strcasecmp(argv[i], unsupported_vcftools_options[j].str) == 0) {
                        found = 1;
                        break;
                    }
                }
                if (found == 1)
                    errorvf("Mega2 does not support the VCFtools command line option '%s'.\n", argv[i]);
                else
                    errorvf("Invalid VCFtools command line option '%s'.\n", argv[i]);
                free(argv[0]); free(argv); free(cmd_line);
                return 0; // an error code rather than just removing the offender and continuing...
            }
        } // if (strncmp(argv[i], "--", 2) == 0) {
    
    params = new parameters(argc, argv);
//  free(argv[0]); free(argv); free(cmd_line);
    
    // Needed to modify:
    // parameters.cpp:parameters::error() to do a throw rather than an exit().
    //
    // NOTE: parameters.cpp:check_parameters() will throw an error if the vcf_filename
    // (e.g., '--vcf file', or '--bcf file', or '--gzvcf file') option is not specified
    // at this point in time, so it needs to be specified by Mega2 by appending to VCFArgs.
    try {
        params->read_parameters();
    } catch (int code) {
        // VCFtools determined a problem while reading/processing the command line arguments.
        // Clean up, abort, and give the error code returned by VCFtools.
        delete params;
        params = NULL;
        return code; // an error was found
    }
    
    // For now we do not deal with indels...
    //printf("keep the exons, forcing remove the indels...");
    //params->keep_only_indels = false; params->remove_indels = true;
    
    // ask vcftools to print the command options that it has successfully parsed...
    //params->print_params();
    
    return -1; // no error
}

//
// Use this one to check the parameters that the user is permitted to select...
int VCFtools_process_cmd_line_wo_file(const char *VCFArgs)
{
    // Skip over the three --file types...
    return VCFtools_process_cmd_line(VCFArgs, 3);
}

//
// This one includes the user selected parameters as well as the file...
// Mega2 will add the file to the command string. The user is not permitted
// to do this.
int VCFtools_process_cmd_line_w_file(const char *VCFArgs)
{
    return VCFtools_process_cmd_line(VCFArgs, 0);
}

/**
 If the user takes the default from the Mega2 menu (e.g., continue), then
 the default VCFArgs will not have been process. If that is the case, then
 we process the arguments that we have on hand (which should be the default).
 It's important not to use this in the menu system where the user may be changing
 command line parameters. It is a last check done after menu processing.
 
 @return -1 for OK, otherwise a VCFtools error code.
 */
int VCFtools_process_cmd_line_if_necessary_inclusive(const char *VCFArgs)
{
    int ret = 0;
    
    // Should only be 'NULL' if the user took the default...
    if (params == (parameters *)NULL) {
        if ((ret = VCFtools_process_cmd_line_w_file(VCFArgs)) != -1) return ret;
    }
    
    // no error since params had already been successfully processed...
    return -1;
}

/**
 Process the file up to the header line (e.g., line beginning with a single '#' character).
 This is done using the information collected in the 'params' structure, so the function
 'VCFtools_process_cmd_line()' must have been called before calling this.
 */
void VCFtools_process_file_meta_information_and_header()
{
    string mssg;
    
    if (params == (parameters *)NULL) {
        // This should never happen, as the calling code should have made sure that this was done.
        errorf("INTERNAL: The VCFtools command line arguments were not parsed.");
        EXIT(SYSTEM_ERROR);
    } else if (vf != (variant_file *)NULL || e != (entry *)NULL) {
        errorf("INTERNAL: The VCF file meta information and header has already been parsed.");
        EXIT(SYSTEM_ERROR);
    }
    
    if (!params->bcf_format)
        vf = new vcf_file(params->vcf_filename, params->vcf_compressed, params->chrs_to_keep, params->chrs_to_exclude, params->force_write_index);
    else
        vf = new bcf_file(params->vcf_filename, params->chrs_to_keep, params->chrs_to_exclude, params->force_write_index, params->gatk);
    
    Tod vcfap("VCF apply filters");
    vf->apply_filters(*params);
    vcfap();

    mssgf("After application of VCFtools filtering:");
    mssg = "Kept " + output_log::int2str(vf->N_kept_individuals()) +
        " out of " + output_log::int2str(vf->N_total_indv()) + " Individuals";
    mssgf(mssg.c_str());
    if (vf->N_kept_sites() == 0) {
        errorf("Filtering constraints have excluded all locus/sites.");
        errorf("No locus/sites left for analysis in the VCF file.");
        EXIT(DATA_TYPE_ERROR);
    }
    mssg = "Kept " + output_log::int2str(vf->N_kept_sites()) +
        " out of a possible " + output_log::int2str(vf->N_total_sites()) + " Loci/Sites";
    mssgf(mssg.c_str());
    mssgf("");
    
    // Create an empty VCF entry that can hold information for the 'N' samples that
    // we will process using the filters ...
    e = vf->get_entry_object((unsigned int)vf->N_total_indv());
}

/**

 The phenotype file allows for an additional SAMPLEID column when processing VCF files.
 phe_lookup.cpp::phe::make() will detect the SAMPLEID column and create a map in
 phe_lookup_ext.h::SAMPLEID. This column maps the 'sample' found in the VCF file to
 the Mega2 <ped, per>.

 */

/**
 This routine processes the SAMPLEID information in the phenotype file, and the sample information
 in the VCF file and will give the user a message if any of the following occurs:
 1) a SAMPLEID in the Phenotype File is not found (or found multiple times) as a VCF file
    sample label (found in the VCF file header line)
 2) a SAMPLEID entry has been excluded through the vcftools filtering mechanism
 3) a sample label in the VCF file is not found in the SAMPLEID column of the phenotype file.
 */
static void phenotype_file_SAMPLEID_entry_checks()
{
    if (SAMPLEIDS == NULL) {
        mssgf("The SAMPLEID column was not included in the phenotype file.");
        return;
    }
    mssgf("The SAMPLEID column was included in the phenotype file.");
    
    // Iterate over the SAMPLEID information from the phenotype file...
    for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
        string map_sample = it->first;
        int found = 0;
        
        // Look for the SAMPLEID (map_sample) in the VCF file...
        for (int ui=0; ui < vf->N_total_indv(); ui++) {
            string vcf_file_sample = vf->indv[(size_t)ui];
            if (vcf_file_sample == map_sample) {
                // It was found...
                if (vf->include_indv[(size_t)ui] == false && found == 0) {
                    mssgvf("SAMPLEID '%s' has been excluded by the vcftools filters.\n", map_sample.c_str());
                }
                found++; // SAMPLEID found in the VCF file...
            }
        }
        if (found == 0) {
            mssgvf("SAMPLEID '%s' was not found in the VCF file.\n", map_sample.c_str());
        } else if (found > 1) {
            mssgvf("SAMPLEID '%s' was found multiple times in the VCF file.\n", map_sample.c_str());
        } // else if (fount == 1) all is well!
    } // for (sample_map_type::iterator it ...
    
    // Look through all of the samples in the VCF file....
    for (int ui=0; ui < vf->N_total_indv(); ui++) {
        if (vf->include_indv[(size_t)ui] == false) continue;
        string vcf_file_sample = vf->indv[(size_t)ui];
        int found = 0;
        // Determine if it matches any sample in the SAMPLEID column of the phenotype file...
        for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
            string map_sample = it->first;
            if (vcf_file_sample == map_sample) {
                found++;
                break;
            }
        }
        if (found == 0) {
            mssgvf("The VCF file sample '%s' was not found in the SAMPLEID column\n", vcf_file_sample.c_str());
            mssgf("included in the phenotype file.");
        }
    }
}

//
// Return a vector that contains the index into vf->indv[] for each person in persons[].
// -1 is used to indicate that that no match was found for a persons[] in the VCFfile.
static const vector<int> build_person_indv_v(const annotated_ped_rec persons[],
                                             const unsigned int person_n)
{
    vector<int> person_indv_v(person_n); // created on the stack
    SECTION_ERR_INIT(untyped);
    
    // For each person in Mega2....
    for (unsigned int pi=0; pi<person_n; pi++) {
        string map_sample; // what we match on in the sample
        const char *ped = persons[pi].Pedigree;
        const char *per = persons[pi].ID;
        
        // Initially assume that there is no data in the VCF file for this <ped,per>...
        person_indv_v[pi] = -1;
        
        // Determine if we match on the SAMPLEID value or just the per.
        //
        // If the SAMPLEID column exists in the phenotype file, then only the
        // [ VCF ID -> SAMPLEID -> (ped,per)] mapping should be used.
        // If the SAMPLEID column does not exist in the fam file, then only the
        // [VCF ID -> per ] mapping should be used.
        if (SAMPLEIDS == NULL) {
            // Search for a per mapping...
            map_sample = string(per);
        } else {
            pair<string,string> key = pair<string,string>(string(ped),string(per));
            // Look for a <ped,per> mapping from the SAMPLEID value of the phenotype file...
            for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
                if (key == it->second) {
                    map_sample = it->first; // value of the SAMPLEID column for this <ped,per>
                    break;
                }
            }
        }
        
        // Scan the VCF file heater sample strings for a match. If it is found then record the index
        // of the vf->indv[] where it occurred. We also look to detect no matches and multiple matches.
        vector<int> vcf_sample_matches;
        for (int ii=0; ii < vf->N_total_indv(); ii++) {
            string vcf_sample = vf->indv[(size_t)ii];
            if (vf->include_indv[(size_t)ii] == false) continue;
            if (map_sample == vcf_sample)
                vcf_sample_matches.push_back(ii);
        }
        
        if (vcf_sample_matches.size() == 0) {
            // NOTE: The actual assignment of the missing genotype happens in VCFtools_process_next_entry() when
            // persons_indv_v[pi] == -1 where the individual genotype is made 0/0.
            SECTION_ERR(untyped);
            warnvf("Pedigree %s person %s not in vcf file so will be untyped.\n", ped, per);
        } else if (vcf_sample_matches.size() == 1) {
            // On a good day, we should find just one match....
            person_indv_v[pi] = vcf_sample_matches[0];
            // But, if there is no SAMPLEID column at all in the phenotype file and if any ID in the
            // VCF file is non-unique in the fam file, then stop with an error.
            if (SAMPLEIDS == NULL) {
                // So here after having found a match in the VCF file, we look for duplicate 'per's in the fam file...
                unsigned int per_found = 0;
                for (unsigned int pi2=0; pi2<person_n; pi2++) {
                    const char *per2 = persons[pi2].ID;
                    if (strcasecmp(per, per2) == 0) per_found++;
                    if (per_found > 1 ) {
                        errorvf("Person '%s' is not unique in the fam file.\n", per);
                        errorf("To map person IDs to VCF samples, the person IDs must be unique.");
                        errorf("Or add a SAMPLEID column to the phenotype file.");
                        EXIT(DATA_TYPE_ERROR);
                    }
                }
            }
        } else { // vcf_sample_matches.size() > 1
            if (SAMPLEIDS == NULL)
                errorvf("Person '%s' maps to multiple VCF file entries.\n", map_sample.c_str());
            else
	        errorvf("SAMPLEID '%s' maps to multiple VCF file entries.\n", map_sample.c_str());
            EXIT(DATA_TYPE_ERROR);
        }
    }
    SECTION_ERR_FINI(untyped);

    return person_indv_v;
}

static string get_marker_name(const string info_id_alternative_key,
                              const string unknown_marker_prefix);

/**
   Process the next entry in the VCF file after doing some sanity checking.

   @return true(1)/false(0) depending on whether entry is processed...
 */
static int VCFtools_process_next_entry(const unsigned int entry_i,
                                       const annotated_ped_rec persons[],
                                       const unsigned int person_n,
                                       const vector<int> &persons_indv_v,
                                       linkage_locus_top *LTop,
                                       string info_id_alternative_key, 
                                       string unknown_marker_prefix)
{
    vector<string> Alleles;
    vector<char *> alleles;
    char phase;
    pair<int, int> genotype;
    vector<char> variant_line;
    char *allele1, *allele2;
    char *allele0 = canonical_allele("0");
    string ID;
    int Locus_i;

    // Only entries (markers) that have passed the VCFtools filtering will be built into the names file
    // and thus will have Mega2 Loci...
    if (vf->include_entry[entry_i] == false) return 0; // This entry was not processed

    // read, load, and parse the VCF file data entry (marker and associated samples)...
//    asm("int $3");
    vf->get_entry(entry_i, variant_line);
    e->reset(variant_line);
    e->parse_basic_entry(true,true,true); // one marker's worth of data...
    e->get_alleles_vector(Alleles);
  
    for (vector<string>::iterator I = Alleles.begin();
         I != Alleles.end();) {
        alleles.push_back(canonical_allele((*I++).c_str()));
    }

    // Get the index in LTop->Locus[] for the marker/locus associated with the entry ID string...
    ID = get_marker_name(info_id_alternative_key, unknown_marker_prefix);

    if (!search_marker(ID.c_str(), &Locus_i)) {
        errorvf("INTERNAL: Unable to find ID '%s' from VCF file in Mega2 Locus.\n", ID.c_str());
        EXIT(SYSTEM_ERROR);
    }
    // Consistency check: Make sure that the marker names (Locus and VCF file marker entry) really do match....
    if (strcmp(LTop->Locus[Locus_i].Name,ID.c_str()) != 0) {
        errorvf("INTERNAL: Mega2 Locus name '%s' does not match ID '%s' from VCF file.\n",
                LTop->Locus[Locus_i].Name, ID.c_str());
        EXIT(SYSTEM_ERROR);
    }
    // Consistency check: Make sure that the type of the Locus is really a marker...
    linkage_locus_rec *LLR = &LTop->Locus[Locus_i];
    linkage_locus_type lltype = LLR->Type;
    linkage_locus_class llclass = LLR->Class;
    if (llclass != MARKER || !(lltype == NUMBERED || lltype == XLINKED || lltype == YLINKED)) {
        errorvf("INTERNAL: Mega2 Locus name '%s' is not a Marker.\n", LTop->Locus[Locus_i].Name);
        EXIT(SYSTEM_ERROR);
    }
    
    // CHECK THIS:
    // The interactively built batch file for VCF to PLINK (.bed) did not run in batch mode.
    // Asks what PLINK output format you want.
    
    // We are iterating over the persons[] vector. Because the .fam file is what drives the lookup.
    // Cases:
    // 1) If an individual in the fam file is not found in the VCF file then they get assigned a missing genotype
    //   (tell the user that the individual has been assigned the missing genotype once per individual).
    // 2) If an individual in the VCF file is excluded from the .fam file then they are not considered for processing
    //   (tell the user that the entry in the VCF file will be ignored).
    for (unsigned int pi=0; pi<person_n; pi++) {
        
        if (persons_indv_v[pi] < 0) {
            // No match was found for the persons[] entry in the VCF file, OR
            // a match was found in the VCF file BUT the person was excluded because of VCFtools filtering.

            // Don't give a warning here because we give one earlier in build_person_indv_v()

            // In this case, make the individual 0/0...
            set_2Ralleles(persons[pi].marker, Locus_i, &LTop->Locus[Locus_i], allele0, allele0);
        } else {
            // Using code patterned after variant_file_format_convert.cpp::variant_file::output_as_plink()

            // A match was found for the persons[] entry in the VCF file, AND
            // the person was not excluded because of VCFtools filtering...
            const unsigned int ui = (const unsigned int)persons_indv_v[pi]; // it must be  >= 0
            genotype = make_pair(-1,-1);
            // NOTE: '|' is used to denote phased, and '/' to denote unphased.
            phase = '/';
            
            if (vf->include_genotype[entry_i][ui] == true) {
                e->parse_genotype_entry(ui, true);
                e->get_indv_GENOTYPE_ids(ui, genotype);
                phase = e->get_indv_PHASE(ui);
            }
            
            // I have to find the appropriate slot in Mega2 for the Locus (ID) associated with the
            // line ID of the entry.
            // Then I need to take all of the REF and ALT alleles that are found in the VCFfile entry
            // and canonicalize the strings.
            
            allele1 = (genotype.first == -1 ?
                       allele0 :
                       alleles[(size_t)genotype.first]);
            // NOTE: Male X-chr, Y-chr etc double allele1 in Mega2 ??? (check this)
            // NOTE: Male X-chr, Y-chr is represented in a VCF file as 'allele1'
            // and not 'allele1|.' but vcf_entry::set_indv_GENOTYPE_and_PHASE() stores
            // it as if it were 'allele1|.'.
            allele2 = (genotype.second != -1 ?
                       alleles[(size_t)genotype.second] :
                       (phase == '/' ? allele0 : allele1));
            
            // annotated_ped_rec *, int, const char *, const char *
            set_2Ralleles(persons[pi].marker, Locus_i, LLR, allele1, allele2);
        } // } else {
    } // for (unsigned int pi=0; pi<person_n; pi++) {
    
    return 1; // This entry was processed
}

//
// Process all of the entries in the VCF file. These are the lines in the file that
// occur after the header.
// This routine should only be called after calling 'VCFtools_process_file_meta_information_and_header()'
// as it assumes the existence of the 'params', and 'vf' objects.
//
// It is a feature that the .fam file can be longer than the number of persons (samples associated with)
// in the VCF file. This means that there are people in the family structure that never gave DNA so they
// won't be genotyped.
void VCFtools_process_entries(annotated_ped_rec persons[],
                              const unsigned int person_n,
                              linkage_locus_top *LTop,
                              string info_id_alternative_key,
                              string unknown_marker_prefix)
{
    if (params == (parameters *)NULL) {
        errorf("INTERNAL: The VCFtools command line arguments were not parsed.");
        EXIT(SYSTEM_ERROR);
    } else if (vf == (variant_file *)NULL || e == (entry *)NULL) {
        errorf("INTERNAL: The VCF file meta information and header was not parsed.");
        EXIT(SYSTEM_ERROR);
    }

    phenotype_file_SAMPLEID_entry_checks();

    // A vector that gives the index of the corresponding value in vf->indv[] for an entry in persons[].
    const vector<int> &person_indv_v = build_person_indv_v(persons, person_n);
    
    // Allow rewinding to a point just after the header...
    streampos file_pos = vf->get_filepos();

    // Loop over all of the lines (entries) in the VCF file...
    // These lines follow the header file which should have been read by the routine
    // VCFtools_process_file_meta_information_and_header()
    Tod vcfall("VCF all entries");
    Tod vcfpne(50);
    for (unsigned int entry_i = 0; entry_i < (unsigned int)vf->N_total_sites(); entry_i++) {
        vcfpne.reset();
        VCFtools_process_next_entry(entry_i, persons, person_n, person_indv_v, LTop,
                                    info_id_alternative_key, unknown_marker_prefix);
        vcfpne("VCF next entry");
    }
    vcfall();
    vf->set_filepos(file_pos);
}

#if 0

/**
  Determine if the pedigree/person is excluded according to the VCFtools filtering criterial.
*/
const bool VCFtools_persons_indv_included(const char *ped, const char *per) {
    const string ped_s = string(ped);
    const string per_s = string(per);
    string id;
    if (SAMPLEIDS != NULL) {
        // Attempt to find the ID from the <ped,per>....
        for (sample_map_type::iterator it = SAMPLEIDS->begin(); it != SAMPLEIDS->end(); it++) {
            pair<string,string> p = it->second; // pair<string,string>(ped, per);
            string pped = p.first;
            string pper = p.second;
            if (ped_s == pped && per_s == pper) {
                id = it->first; // found it...
                break;
            }
        }
        if (id.size() > 0) {
            // found an entry for the <ped,per>, so look for it's inclusion/exclusion entry and return it...
            for (unsigned int ui=0; ui < vf->N_total_indv(); ui++)
                if (id == vf->indv[ui]) return vf->include_indv[ui];
            // We should never get here...
            //return false;
        }
    }
    // Since there was no SAMPLEIDS maping, just look or the per...
    for (unsigned int ui=0; ui < vf->N_total_indv(); ui++)
        if (id == vf->indv[ui]) return vf->include_indv[ui];
    
    return false;
}

/**
 The phenotype file allows for an additional SAMPLEID column when processing VCF files.
 phe_lookup.cpp::phe::make() will detect the SAMPLEID column and create a map in
 phe_lookup_ext.h::SAMPLEID. This column maps the 'sample' found in the VCF file to the Mega2 <ped, per>.
 If no such column is found in the phenotype file then assume that the 'sample' found in the VCF
 file is the Mega2 person id (per). phe_lookup.cpp::phe::make() will insert available mappings
 (SAMPLEID column values) from VCF file 'sample ID' into <'ped', 'per'> that it finds into the map.
 In some instances there will be no mapping, if for example the individual <ped, per> is a
 founder and was never genotyped.
 
 @return -1 if no SAMPLEID column exists, 0 if no mapping found, 1 if mapping found.
 */
static int sample_to_ped_per(const string sample, char **ped, char **per)
{
    // There was no SAMPLEIDS column in the phenotype file...
    if (SAMPLEIDS == NULL) return -1;
    
    // create an iterator to search for the 'sample' in the 'sample_map_type'...
    sample_map_type::iterator it = (*SAMPLEIDS).find(sample);
    // test to see if the 'sample' was not found in the map...
    if (it == SAMPLEIDS->end()) return 0; // not there
    
    // The sample was found...
    pair<string,string> item = it->second;
    *ped = strdup(item.first.c_str());
    *per = strdup(item.second.c_str());
    return 1; // found it
}

/**
   Create a vector of indicies into the persons[] array that correspond to the
   individuals associated with the samples from the VCF file. This vector will
   be used when processing 'entries' in the VCF file to get to the associated
   persons[i] for all entries.

   @return vector<int> indicies into persons[] in the order that they are read from the VCF file.
 */
static const vector<int> gen_mega2_ui(const annotated_ped_rec persons[])
{
    vector<int> mega2_ui_v(vf->N_total_indv()); // created on the stack
    
    // Here we search each sample ID in the VCF file for a corresponding <ped, per> in
    // the Mega2 data structure.
    for (unsigned int ui=0; ui < vf->N_total_indv(); ui++) {
        char *PedID, *PerID;
        int mega2_ui = -1;
        string sample = vf->indv[ui];
        int i;
        
        // Algorithm for matching VCF File sample to Mega2 person/individual:
        // 1) Find the <ped,per> using 'sample_to_ped_per()' and look for it in persons (PedID, PerID).
        // If you don't find it then throw an error.
        // 2) If it is not returned with 'sample_to_ped_per()' then, loop over persons[].ID to find one
        // that matches the sample ID from the VCF file.
        int ret = sample_to_ped_per(sample, &PedID, &PerID);
        if (ret == -1) {
            // There is no SAMPLEID column on in the phenotype file, so attempt to
            // match the sample name with the persons[] id.
            for (i=0; i < vf->N_total_indv(); i++)
                if (strcasecmp(sample.c_str(), persons[i].ID) == 0) {
                    mega2_ui = i;
                    break;
                }
        } else if (ret == 1) {
            // There is a SAMPLEID column in the phenotype file and we found a match, so
            // find the index into the persons[] vector associated with this sample.
            for (i=0; i < vf->N_total_indv(); i++)
                if (strcasecmp(PedID, persons[i].PedID) == 0 && strcasecmp(PerID, persons[i].PerID) == 0) {
                    mega2_ui = i;
                    break;
                }
        }
        // NOTE: That if mega2_ui == -1 at this point, that we were unable to match a sample ID with
        // an individual in Mega2.
        mega2_ui_v[ui] = mega2_ui;
        free(PedID); free(PerID);
    }
    
    // pass by value, so a copy is made of the vector on the stack by the compiler and that copy is returned
    // to the caller. This gets us out of the business of dealing with 'new' and pointers in the code.
    return mega2_ui_v;
}

#endif /* 0 */

/**
 Marker Selection Algorithm
 --------------------------
 
 The marker can be specified as the string associated with the 'ID' field of the VCF file, or as the
 value associated with a key in the 'INFO' field. If the user specifies a key (not empty) then use only
 the INFO field, otherwise use only the 'ID' field. One or the other is used, never switching between them.
 
 If from this the string value is a '.' (unknown) then construct a marker name from a prefix ('chr') appending it to the
 CHROM and POS in the following way (chrCHROM_POS). This is not guaranteed to be unique since the VCF file specification
 allows multiple records to have the same POS. In the construction of this marker name use 'X', 'Y', 'XY', and 'MT' rather
 than the associated numbers for the chromosomes. This is done in accordance to the species option (Mega2 manual 20.1)
 that Mega2 is told to use. There is a table for supported species in the code.
 
 The routine assumes that 'params' has been setup using 'VCFtools_process_cmd_line_if_necessary_inclusive()'
 */

//
// Convert the VCF file CHROM to something that Mega2 might like...
static const string get_Mega2_CHROM_name()
{
    stringstream ss(e->get_CHROM());
    int chrom_num;
    if ((ss >> chrom_num).good()) {
        // Here CHROM was successfully converted to an int, so we convert it to the 'char *'
        // representation as appropriate for the organism...
        char chrom_name[1000];
        chrom_num_to_name(chrom_num, chrom_name); // from genetic_utils.cpp
        return string(chrom_name);
    } else {
        // We could not convert CHROM to an integer. So, just use whatever the string is...
        return e->get_CHROM();
    }
}

//
// Determine what is used for the marker name (see "Marker Selection Algorithm")...
static string get_marker_name(const string info_id_alternative_key,
                              const string unknown_marker_prefix)
{
    string id = (info_id_alternative_key.length() == 0 ? e->get_ID() : e->get_INFO_value(info_id_alternative_key));
    // get_INFO_value() will return "?" if it cannot find the key
    // You must remember to '-DNDEBUG' while compiling to turn off the assert macros in the VCFtools code
    if (id == "?") {
        errorvf("The requested '%s' sub-field of the VCF INFO field does not exist.\n", info_id_alternative_key.c_str());
        EXIT(DATA_TYPE_ERROR);
    } else if (id == ".") {
        // Create a marker name from the prefix, processed CHROM, and POS...
        ostringstream id_s;
        id_s << unknown_marker_prefix << get_Mega2_CHROM_name() << "_" << e->get_POS();
        id = id_s.str();
    }
    return string(id);
}

//
// Here we strip the .MAP file from the 'left side' of the VCF file for those entries
// ""that have passed the VCFtools filtering criteria"". The 'file' is created as a
// m2_map that contains one m2_map_entry for each 'line' in the VCF file.. The m2_map_entries
// are stored in a vector which orders elements as they appear in the VCF file.
// No tests are made to see that the CHROM and POS within them are increasing as is stated
// in the VCF File standard.
// http://www.1000genomes.org/node/101
m2_map VCFtools_get_map(const std::string info_id_alternative_key,
                        const std::string unknown_marker_prefix)
{
    vector<char> variant_line;
    m2_map map("VCF", 'p');

    if (params == (parameters *)NULL) {
        errorf("INTERNAL: The VCFtools command line arguments were not parsed.");
        EXIT(SYSTEM_ERROR);
    } else if (vf == (variant_file *)NULL || e == (entry *)NULL) {
        errorf("INTERNAL: The VCF file meta information and header was not parsed.");
        EXIT(SYSTEM_ERROR);
    }
    
    // Allow rewinding to a point just after the header...
    streampos file_pos = vf->get_filepos();
    
    human_x = human_y = human_xy = human_unknown = human_mt = human_auto = 0;

    Tod vcfme("vcf map entry: one line", 50);
    Tod vcffetch(50);
    Tod vcfstore(50);
    Tod vcfl1(50);
    for (unsigned int entry_i = 0; entry_i < (unsigned int)vf->N_total_sites(); entry_i++) {
        // If it didn't pass the filtering criteria, don't include the marker in the map...
        if (vf->include_entry[entry_i] == false) continue;
	vcfme.reset();
        vcffetch.reset();

        vcfl1.reset();
        vf->get_entry(entry_i, variant_line); // 256
        vcfl1("l1");
        e->reset(variant_line);               // 31
        e->parse_basic_entry(true);           //  5
        vcffetch("vcf fetch");                // 299

        vcfstore.reset();
        m2_map_entry map_entry;
        map_entry.set_chr(e->get_CHROM());
        map_entry.set_POS(e->get_POS());
        map_entry.set_marker_name(get_marker_name(info_id_alternative_key, unknown_marker_prefix));
	map_entry.set_REF(e->get_REF());
        map.push_back_entry(map_entry);
        vcfstore("vcf store");                   // 2?

        int chr = STR_CHR((e->get_CHROM()).c_str());
        if (chr == SEX_CHROMOSOME) {
            human_x++;
//          if (LTop->Locus[mrk_num].Type != XLINKED) {
//              LTop->Locus[mrk_num].Type = XLINKED;
//          }
        } else if (chr == MALE_CHROMOSOME) {
            human_y++;
//          if (LTop->Locus[mrk_num].Type != YLINKED) {
//              LTop->Locus[mrk_num].Type = YLINKED;
//          }
        } else if (chr == PSEUDO_X) {
            human_xy++;
        } else if (chr == MITO_CHROMOSOME) {
            human_mt++;
        } else if (chr == UNKNOWN_CHROMO) {
            human_unknown++;
            NumUnmapped++;
        } else {
            human_auto++;
        }
        vcfme();
    }
    
    vf->set_filepos(file_pos); // rewind the stream...
    
    return map;
}
