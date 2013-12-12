/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2013 Robert Baron, Charles P. Kollar,
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
  Mega2 interface to VCFTools.

  The goal is to make minimal changes to VCFTools.

  The following changes were made to the vcftools source. Additional changes
  are noted below.

  variant_file.h: comment out //#include "gamma.h"
*/

#include <string.h>
#include <errno.h>
#include <iostream>
#include <fstream>

// Mega2 includes....
#include "common.h"
#include "typedefs.h"
#include "utils_ext.h"
#include "error_messages_ext.h"
#include "makeped.h"
#include "annotated_ped_file.h"
#include "linkage.h"
#include "read_files_ext.h"
#include "compress_ext.h"

// VCFTOOLS includes...
// NOTE: that vcftools was modified to use the 'vcftools' namespace because
// of symbol clashes with Mega2 (in particular 'entity').
#include "output_log.h"
#include "parameters.h"
#include "variant_file.h"
#include "bcf_file.h"
#include "vcf_file.h"

// This is what Mega2 uses...
using namespace std;

//
// The 'output_log' class is used eveywhere in VCFTools.
// Here we graft on the Mega2 logging mechanism to the methods used by output_log.
// These methods take the place of anything that was in the file 'output_log.cpp' from the VCFTools distribution.

output_log LOG;

output_log::output_log() : output_to_screen(true) { }

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


//
// This is the flag by which VCFTools is interfaced with Mega2. In general a vcf
// file is processed (for the moment) as a PLINK SNP Major Mode file where...
// ==0, PLINK Binary mode; ==1, VCF mode.
// in VCF mode...
// 1) PLINKArgs contains the VCFTools commandline arguments. These behave the same
// as with VCFTools with the excepiton that the keywords --vcf, --bcf, and --gzvcf
// don't take an argument filename since the filename is passed in through the
// Mega2 interface.
// The .vcf/.bcf/.vcf.gz file name is given as the .bed file.
// 2) The .bim file is actually a .map file since the last two columns of
// a .bim (extended map) file are not needed when processing a .vcf file. Only the
// four columns of the .map file (chr, snp, pd, gd) are needed.
int VCF=1; 

/*
A Mega2 batch file would look like this:

# In PLINK binary mode this would be a .bim file (two extra columns for allele names).
Input_Map_File=ped.map
# In VCF mode this is either a .vcf/.bcf/.vcf.gz file with the 'PLINK' argumentset
# appropriately, e.g. absent or--vcf for a .vcf file; --bcf for a .bcffile; -gzvcf for a .vcf.gz file.
Input_Binary_File=bed.vcf
PLINK=--remove-indels
*/

//
// vcftools: Usage and Options
// http://vcftools.sourceforge.net/options.html#basic

enum argType {
  NONE,
  CONSTANT,
  FILENAME
};

typedef struct _command_line_option {
    const char *str;       // The option string e.g., '--foo'
    const enum argType at; // e.g., 0: '--foo'; 1: '--foo value'; 2: -foo filename
} command_line_option;

//
// Mega2 only allows for the use of the VCFTools 'Basic Options', and 'Filters'.

static const command_line_option supported_vcftools_options[] = {
  // Basic Options (all but '--out' are used)...
#ifdef VCFTOOLS_FILENAME_REQUIRED
  {"--vcf", FILENAME},
  {"--gzvcf", FILENAME},
  {"--bcf", FILENAME},
#else /* VCFTOOLS_FILENAME_REQUIRED */
  {"--vcf", NONE},
  {"--gzvcf", NONE},
  {"--bcf", NONE},
#endif /* VCFTOOLS_FILENAME_REQUIRED */
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
  {"--remove-indv", FILENAME},
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
#define SUPPORTED_VCFTOOLS_OPTIONS_COUNT (sizeof(supported_vcftools_options) / sizeof(command_line_option))

//
// Mega2 only supports the filtering command line options for now.
// Anything in this list will not be supported by mega2.
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
#define UNSUPPORTED_VCFTOOLS_OPTIONS_COUNT (sizeof(unsupported_vcftools_options) / sizeof(command_line_option))

//
// Print for the user a list of the command line options that are supported by Mega2.
void VCFTools_printf_supported_cmd_line_options()
{
  int j, i, len;
  static const char *arg = "v";

  for (j=0,i=0; j<SUPPORTED_VCFTOOLS_OPTIONS_COUNT; j++) {
    // NOTE: if the length of any of these parameters are ever >= 80
    // then this will result in an infinite loop...
    len = (int)strlen(supported_vcftools_options[j].str) + 1;
    len += (supported_vcftools_options[j].at != NONE ? strlen(arg)+1 : 0);
    if ((i += len) >= 80) {
      printf("\n"); // line is never >= 80 characters
      i = len;
    }
    printf("%s ", supported_vcftools_options[j].str);
    if (supported_vcftools_options[j].at != NONE) printf("%s ", arg);
  }
  printf("\n");
}

//
// VCFTools class insances that are used by Mega2...
static parameters *params = NULL;
static variant_file *vf = NULL;
// The name 'entry' defined in vcftools source 'entry.h/cpp' conflicts with a symbol in Mega2. 
static variant_file_entry *e = NULL;
static int entry_i = 0;

//
// Thse functions take the place of vcftools.cpp:main() and serves as the initialization of VCFTools for Mega2.

//
// This function is called in user_input.cpp:menu1(). It will create and fill 'params' which
// is used for further processing from the 'VCFArgs' command line arguments given. Only VCFTools
// command line options supported by Mega2 are permitted. It tries to provide some additional
// checking of the arguments that VCFTools does not.
//
// @return -1 for OK, otherwise a VCFTools error code.
int VCFTools_process_cmd_line(const char *VCFArgs)
{
    int argc = 1, i, j, found;
    char **argv = (char **)malloc(sizeof(char *));
    char *cmd_line = strdup(VCFArgs);
    
    // Strip any newline character from the end of the string...
    if (cmd_line[strlen(cmd_line)-1] == '\n') cmd_line[strlen(cmd_line)-1] = '\0';
    if (cmd_line[strlen(cmd_line)-1] == '\r') cmd_line[strlen(cmd_line)-1] = '\0';
    
    // The VCFTools code assumes the name of the program is the first argument on the
    // command line, as it would be if the program were run from a shell, so put it there...
    argv[0] = strdup("VCFTools");
    
    // Break the VCFTools command line into argc, and argv used by the parameters constructor...
    char *p = strtok((char *)cmd_line, " ");
    while (p) {
      if ((argv = (char **)realloc(argv, sizeof(char *) * ++argc)) == (char **)NULL)  exit (-1); /* memory allocation failed */
        argv[argc-1] = p;
        p = strtok (NULL, " ");
    }
    
    // Search the command line (following argv[0]) for supported options...
    for (i=1; i<argc; i++)
        // All command line options begin with '--'...
        if (strncmp(argv[i], "--", 2) == 0) {
            // Try to locate the option in the list that we support...
            for (j=0, found=0; j<SUPPORTED_VCFTOOLS_OPTIONS_COUNT; j++) {
                if (strcasecmp(argv[i], supported_vcftools_options[j].str) == 0) {
                    found = 1;
                    break;
                }
            }
            if (found == 1) {
                // Good, the option found!
                // If the option does not have an argument, don't check it...
                if (supported_vcftools_options[j].at == NONE) continue;
                
                // Check for the exinstance of what appears to be a valid argument.
                // Valid is define as: not being at the end of the argument list, or
                // not beginning with '--'.
                if (i+1 >= argc || strncmp(argv[i+1], "--", 2) == 0) {
                    errorvf("The argument required for option '%s' appears to be missing.\n",
                           supported_vcftools_options[j].str);
                    free(argv);
                    return 0; // an error occurred...
                }
                
                // Depending on the agument type, do some checking...
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
                        free(argv);
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
                    errorvf("Mega2 does not support the VCFTools command line option '%s'.\n", argv[i]);
                else
                    errorvf("Unknown VCFTools command line option '%s'.\n", argv[i]);
                free(argv);
                return 0; // an error code rather than just removing the offender and continuing...
            }
        } // if (strncmp(argv[i], "--", 2) == 0) {
    
    params = new parameters(argc, argv);
//    free(argv); argv = NULL; argc = 0;
    
    // Needed to modify:
    // parameters.cpp:parameters::error() to do a throw rather than an exit().
    //
    // NOTE: parameters.cpp:check_parameters() will throw an error if the vcf_filename
    // (e.g., '--vcf file', or '--bcf file', or '--gzvcf file') option is not specified
    // at this point in time, so it needs to be specified by the user or appended to VCFArgs.
    try {
        params->read_parameters();
    } catch (int code) {
        // VCFTools determined a problem while reading/processing the command line arguments.
        // Clean up, abort, and give the error code returned by VCFTools.
        delete params;
        params = NULL;
        return code; // an error was found
    }
    
    // For now we do not deal with indels...
    //printf("keep the exons, forcing remove the indels...");
    //params->keep_only_indels = false; params->remove_indels = true;
    
    //params->print_params();
    
    return -1; // no error
}

//
// If the user takes the default from the Mega2 menu (e.g., continue), then
// the default VCFArgs will not have been process. If that is the case, then
// we process the arguments that we have on hand (which should be the default).
// 
// @return -1 for OK, otherwise a VCFTools error code.
int VCFTools_process_cmd_line_if_necessary(const char *VCFArgs, const char *filename)
{
    int ret = 0;
    
    // Should only be 'NULL' if the user took the default...
    if (params == (parameters *)NULL)
      if ((ret = VCFTools_process_cmd_line(VCFArgs)) != -1) return ret;

    params->vcf_filename = filename;

    // no error since params had already been successfully processed...
    return -1;
}

//
// Process the file up to the header line (e.g., line beginning with a single '#' character).
void VCFTools_process_file_meta_information_and_header()
{
    string mssg;
    
    if (params == (parameters *)NULL) {
        // This should never happen, as the calling code should have made sure that this was done.
        errorf("INTERNAL: The VCFTools command line arguments were not parsed.");
        EXIT(SYSTEM_ERROR);
    }

    if (!params->bcf_format)
        vf = new vcf_file(params->vcf_filename, params->vcf_compressed, params->chrs_to_keep, params->chrs_to_exclude, params->force_write_index);
    else
        vf = new bcf_file(params->vcf_filename, params->chrs_to_keep, params->chrs_to_exclude, params->force_write_index, params->gatk);
    
    vf->apply_filters(*params);
    
    unsigned int N_indv = vf->N_kept_individuals();
    unsigned int N_sites = vf->N_kept_sites();
    unsigned int N_total_indv = vf->N_total_indv();
    unsigned int N_total_sites = vf->N_total_sites();
    mssgf("After application of VCFTools filtering:");
    mssg = "Kept " + output_log::int2str(N_indv) + " out of " + output_log::int2str(N_total_indv) + " Individuals";
    mssgf(mssg.c_str());
    mssg = "Kept " + output_log::int2str(N_sites) + " out of a possible " + output_log::int2str(N_total_sites) + " Sites";
    mssgf(mssg.c_str());
    if (N_sites == 0) {
        errorf("No sites left for analysis in the VCF file.");
        exit(-1);
    }

    // Create an empty VCF entry that can hold information for the 'N' individuals that
    // we are considering from the filters ...
    e = vf->get_entry_object(N_indv);
    entry_i = 0;
}

int VCFTools_process_next_entry(annotated_ped_rec persons[],
                                const int person_n,
                                linkage_locus_rec Locus[],
                                const int Locus_i)
{
    vector<string> alleles;
    char phase;
    pair<int, int> genotype;
    vector<char> variant_line;
    string allele1, allele2;
    
    if (params == (parameters *)NULL) {
        errorf("INTERNAL: The VCFTools command line arguments were not parsed.");
        EXIT(SYSTEM_ERROR);
    } else if (vf == (variant_file *)NULL || e == (variant_file_entry *)NULL) {
        errorf("INTERNAL: The VCF file was not parsed.");
        EXIT(SYSTEM_ERROR);
    } else if (vf->N_indv != person_n) {
        errorf("INTERNAL: The VCF file person count error.");
        EXIT(SYSTEM_ERROR);
    }
    
    if (vf->include_entry[entry_i] == false) {
        // Mark the locus so that it will be ignore by subsequent Mega2 processing
        // since it was excluded by the VCFTools rules...
        Locus[Locus_i].Class = CLASS_UNSET; Locus[Locus_i].Type = TYPE_UNSET;
        entry_i++;
        return 0; // didn't process this entry
    }
    
    // read, load, and parse the VCF file data entry...
    vf->get_entry(entry_i, variant_line);
    e->reset(variant_line);
    e->parse_basic_entry(true);
    e->get_alleles_vector(alleles);
    
    // For each individual in the VCF file...
    for (unsigned int ui=0; ui<vf->N_indv; ui++) {
        char *canonical_allele1, *canonical_allele2;
        
        if (vf->include_indv[ui] == false) {
            // make the individual 0/0 in this case.
            canonical_allele1 = canonical_allele("0");
            set_2Ralleles(persons[ui].marker, Locus_i, canonical_allele1, canonical_allele1);
        } else {
            genotype = make_pair(-1,-1);
            // NOTE: '|' is used to denote phased, and '/' to denote unphased.
            phase = '/';
            if (vf->include_genotype[entry_i][ui] == true) {
                e->parse_genotype_entry(ui, true);
                e->get_indv_GENOTYPE_ids(ui, genotype);
                phase = e->get_indv_PHASE(ui);
            }
            
            allele1 = (genotype.first == -1 ?
                       "0" :
                       alleles[genotype.first]);
            // NOTE: Male X-chr, Y-chr etc double allele1 in Mega2 ??? (check this)
            // NOTE: Male X-chr, Y-chr is represented in a VCF file as 'allele1'
            // and not 'allele1|.' but vcf_entry::set_indv_GENOTYPE_and_PHASE() stores
            // it as if it were 'allele1|.'.
            allele2 = (genotype.second != -1 ?
                       alleles[genotype.second] :
                       (phase == '/' ? "0" : allele1));
            
            canonical_allele1 = canonical_allele(allele1.c_str());
            canonical_allele2 = canonical_allele(allele2.c_str());
            // annotated_ped_rec *, int, const char *, const char *
            set_2Ralleles(persons[ui].marker, Locus_i, canonical_allele1, canonical_allele2);
        }
    }
    entry_i++;
    return 1; // process this entry
}
