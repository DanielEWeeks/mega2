/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

/* mega.c - original MEGA modified to MEGA1 by Lee Almasy and
   Daniel E. Weeks.
   Further modified by William Mulvihill. 12/97
   Distributed as-is. We make no warranties or guarantees. The only
   promise that we make is that we try to make this code as reliable
   and accurate as may be reasonably expected.

   Based upon chapm (written by Mark Schroeder).

   ============================= mega2 =================================

   This program converts LINKAGE and APM format files and APM
   files to different APM formats. The basic scheme is:

   o Read the input file(s) into data structures.

   o Convert the data structures to APM (pedtree) structures,
   if necessary. This may require further information
   from the user. Reconnect loops in LINKAGE files.

   o Prompt for the loci to save and their order. Delete
   unused loci and change their order in the internal
   structures.

   o Remove affecteds from the affecteds list that are not
   typed at any locus that we are saving.

   o Check the pedigrees and delete those that have fewer
   than two typed affecteds for the loci we are saving.

   o Finally, save the file in the specified APM format.

   The Mega1's file conversion options will create files to be used as
   input for the various genetic linkage simulation programs. In addition,
   Mega1 can create summary files and has the ability to create nuclear
   families.

   The specific file conversion options of mega1:
   o Haplotype - creation of input files for SIMWALK
   o Random-walk location score analysis
   - create input files for SIMWALK location score analysis
   o Convert to Mendel-format
   o Convert to ASPEX-format
   o Create summary files
   - affected relative count
   - segregation analysis [unknown, unaffected, unaffected]
   - allele distribution tables
   o Convert to GeneHunter-format
   o Convert to APM-format
   o Convert to APM MULT multiple locus version
   o Create nuclear families
   o Create SLINK-format files
   o Create SPLINK-format files

   Mega1 provides the user with the ability to re-arrange the loci order and/or
   select a subset of loci for the output files.

   ======================================================================== */

/* When adding a new output option:

   The case statement in the main body of the program decides which
   top-level function to call for each option. Add your KEYWORD and the
   appropriate call there.

   Typically each analysis option takes as input a POINTER to the
   linkage-pedigree data-structure (LPedTreeTop), output files names list
   (file_names), and the global untyped pedigree option (UntypedPedOpt).

   Post-processing is generally the same across all options, or handled
   by one of the macros. SAGE needed some special handling upon
   completion of the Mega2 run. If your option needs any special
   post-processing, add code after the statement setting MEGA2 status to
   TERM_MEGA2.

*/

#include <stdio.h>
#include <string.h>
#ifdef SOLARIS
#define _STDC_C99
#endif
#include <math.h>
#include <ctype.h>
#include <errno.h>
#ifdef _WIN
#include <process.h>
#endif

#include "common.h"
#include "typedefs.h"
#include "tod.hh"
#include "R_output.h"

#include "annotated_ped_file.h"
#include "annotated_ped_file_ext.h"
#include "aspex_ext.h"
#include "batch_input_ext.h"
#include "cranefoot_ext.h"
#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "error_sim_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "hwe_user_input_ext.h"
#include "input_check_ext.h"
#include "makenucs_ext.h"
#include "makeped_ext.h"
#include "marker_lookup_ext.h"
#include "mega2annot_ext.h"
#include "menu_value_missing_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "pedtree_ext.h"
#include "plink_ext.h"
#include "read_files_ext.h"
#include "reorder_loci_ext.h"
#include "scripts_ext.h"
#include "slink_ext.h"
#include "splink_ext.h"
#include "user_input_ext.h"
#include "vcftools/mega2_vcftools_interface.h"
#include "utils_ext.h"
#include "write_IQLS_ext.h"
#include "write_SUP_files_ext.h"
#include "write_apm_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
#include "write_loki_ext.h"
#include "write_mfiles_ext.h"
#include "write_pap_ext.h"
#include "write_premakeped_ext.h"
#include "write_prest_ext.h"
#include "write_sage_files_ext.h"
#include "write_simulate_files_ext.h"
#include "write_solar_files_ext.h"
#include "write_vitesse_ext.h"
#include "version.h"

#include "class_old.h"
#include "write_shapeit_ext.h"
#include "database_dump_ext.h"

/*
 annotated_ped_file_ext.h:  check_annotated_file_format read_annotated_files Free_annotated_files
              aspex_ext.h:  create_aspex_files
        batch_input_ext.h:  batchfile_init_Mega2BatchItems batchfile_process batchf
          cranefoot_ext.h:  create_CRANEFOOT_files
     create_summary_ext.h:  create_summary_file
     error_messages_ext.h:  errorf mssgf my_calloc open_logs time_stamp_logs warnf
          error_sim_ext.h:  simulate_errors
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  LogFileNames backup_file getRunDate mega2_version_check
     hwe_user_input_ext.h:  hwe_user_input
        input_check_ext.h:  full_check
           makenucs_ext.h:  create_nuclear_families
            makeped_ext.h:  makeped
      marker_lookup_ext.h:  add_allele
         mega2annot_ext.h:  create_mega2annot_files
  output_file_names_ext.h:  default_outfile_names set_output_paths shorten_path
    output_routines_ext.h:  write_key_file write_nuc_ped_key_file
            pedtree_ext.h:  convert_to_pedtree free_all_including_ped_top
         read_files_ext.h:  check_dos_file check_empty check_size_dos free_lpedtop_pedinfo read_linkage2
       reorder_loci_ext.h:  ReOrderLoci
            scripts_ext.h:  create_LOD2_format_files create_TDTMAX
              slink_ext.h:  create_SLINK_format_files
             splink_ext.h:  create_SPLINK
         user_input_ext.h:  analysis_menu1 menu1 ped_ind_defaults
              utils_ext.h:  EXIT goodbye hello log_line mega2_opts mklogdir seed_random
         write_IQLS_ext.h:  create_IQLS_files
    write_SUP_files_ext.h:  create_SUP_files
          write_apm_ext.h:  create_APMULT create_APM_file
        write_files_ext.h:  create_linkage_files write_ped_stats
      write_ghfiles_ext.h:  create_gh_file
         write_loki_ext.h:  create_LOKI_files
       write_mfiles_ext.h:  create_mendel_file
          write_pap_ext.h:  create_pap_files
   write_premakeped_ext.h:  create_merlin_files create_premakeped_files
        write_prest_ext.h:  create_prest_files
   write_sage_files_ext.h:  create_SAGE4_file create_SAGE_file
write_simulate_files_ext.h:  create_SIMULATE_format_files
  write_solar_files_ext.h:  create_SOLAR_files
      write_vitesse_ext.h:  create_vitesse_files
*/


/*Global vars*/

/* int             quiet; */
/* global variables that control output behaviour*/
// strings associated with common.h: genetic_distance_map_type
int             MARKER_SCHEME;
int             marker_scheme_mega2_opts = 0;
const char *genetic_distance_map_type_string[3] = {
  "sex-averaged", "sex-specific", "female"
};
int             LoopOverTrait; /* flag for looping over traits */
int             LoopOverChrm; /* flag for looping over chromosomes */
int             ManualReorder; /* flag denoting whether option 2 of reorder menu was used */
char            ProgName[100]; /* descriptive option name */
analysis_type   AnalysisOpt; /* analysis option */
char            Mega2Log[FILENAME_LENGTH], Mega2LogRun[FILENAME_LENGTH]; /* Log file name */
char            Mega2Err[FILENAME_LENGTH], Mega2ErrRun[FILENAME_LENGTH]; /* Error file name */
char            Mega2Sim[FILENAME_LENGTH], Mega2SimRun[FILENAME_LENGTH]; /* Error sim log file name */
char            Mega2Batch[FILENAME_LENGTH]; /* Batch file name */
char            Mega2Keys[FILENAME_LENGTH], Mega2KeysRun[FILENAME_LENGTH]; /* ID-ped-per map */
char            Mega2Recode[FILENAME_LENGTH], Mega2RecodeRun[FILENAME_LENGTH]; /* Recoded markers */
char            Mega2Reset[FILENAME_LENGTH], Mega2ResetRun[FILENAME_LENGTH]; /* Reset genotypes */
const char      *NOPATH; /* what to print for file name iff in HIDEPATH mode */
const char      *NOcTIME;
char            Mega2Version[20];
int             Mega2Ver, Mega2Rev, Mega2Patch; /* Version */
char            Mega2WebVersion[50]; /* web-site file name */
char            err_msg[2*FILENAME_LENGTH]; /* array for error, warning and log messages */
char            **output_paths; /* trait directories prepended with output_dir */
char            **trait_paths;  /* list of directory names for each trait */
char            InputPath[FILENAME_LENGTH];
char            yorn[2][4]; /* two element list "yes", "no" */
char            nory[2][4]; /* two element list "no ", "yes" */
int             *Labels; /* affected labels (status-liability pairs) */
int             NumLabels; /* Number of Labels */
int             UntypedPedOpt; /* Which peds to omit based on typing */
int             liability_multiplier;
int             ErrorSimOpt; /* Flag for whether to introduce errors */
int             OrigIds[2];
/* which id field to use from the linkage ped file,
   ele1 - for ind id, ele2, for ped id */
char           *Outfile_Names[NUM_OUTFILES];
int            CreateRunFolder;  /* create time-stamped folder for each run? This is the default. Turned off with the -nosave option */
int            check_web_ver;
int            BatchFileCreated; /* 1 if batchfile is created during run */
char           *REC_UNKNOWN;
int            MissingIsNa;

/*  outputfile_ext_type OutFileExtType; */
/*  char           **OutFileExts; */

/* globals to store progress */
int             TimeStampWritten[2];
int             Mega2Status;
char           *Mega2OutputPath;
int             FirstIterMenu;
/* globals that describe input data */
InputModeType   InputMode; // see common.h
InputModeType   AnalyInputMode = NOEXEC_INPUTMODE; // see common.h
file_format     InputFileFormat; /* Annotated or linkage */
char            mega2_path[256]; /* path to mega2 executable */
char            *mega2_input_files[NUMBER_OF_MEGA2_INPUT_FILES]; /* ped, loc, map, freq, pen, and omit files */
char            mega2_input_file_type[NUMBER_OF_MEGA2_INPUT_FILES][24];
int             pedfile_type; /* whether input-file is pre-makeped or not */
int             basefile_type; /* type for original pedfile type */
int             HasLoops;
int             HasAff;
int             HasMarkers;
int             HasQuant;
int             LinkageLocusFile; /* added Dec 2007 */

int             MaxChromo; /* Largest chrom number */
int             NumChromo; /* Number of chromosomes in locus data */
int             *PedTreeMito; /* mitochondrial marker indexes when creating pedtree */
/* variables for storing markers without chromosome numbers */
int             NumUnmapped, *unmapped_markers, UnmappedSelected;
int             AllowUnmapped;

/* reordered_marker_loci = list of marker loci in user-specified or map order
   after locus reordering,
   num_reordered = length of reordered_marker_loci

   global_chromo_entries = chromosome numbers
   chromo_loci_count = number of loci on chromosome (cumulative)
   main_chromocnt = #chromosomes =
   length of global_chromo_entries and chromo_loci_count
*/
int             main_chromocnt, *global_chromo_entries, *chromo_loci_count;
int             num_reordered, *reordered_marker_loci;
int             num_traits, num_affection, num_quantitative;

/*  global _trait_entries = list of trait loci indexes
    may contain an entry -1 to denote position of markers
    may also contain an entry -99
    num_traits = number of traits selected = length of global_trait_entries
    ditto num_covariates and covariates
*/
int             *global_trait_entries;
int             num_covariates, *covariates;

int             seed1, seed2, seed3; /* random seeds */
int             FirstTime; /* First time pedigrees are constructed */
double          MissingQuant;
#ifdef _WIN
double          MissingOutQuant, QuantOutLow, QuantOutHi;
#else
double          MissingOutQuant, QuantOutLow = NAN, QuantOutHi = NAN;
#endif
int             MissingOutQuantSet;
int             NukedMultiplier;

int             SelectIndividuals;
int             HalfTypedReset; /* 1 if user resets half-typed genotypes */
int             NonMendelianReset; /* 1 if user resets inconsistent genotypes */
char            RunDate[FILENAME_LENGTH];
char            bindir[40];
double          FreqMismatchThreshold;

int             HasMapFileBeenRead;
int             HasFreqFileBeenRead;
int             HasPenFileBeenRead;
int             HasPedFileBeenRead;
int             HasLocFileBeenRead;
/* Global variables that control output data */

int             *UntypedPeds; /* vector of flags to indicate whether to leave out a pedigree */
int             NumTypedPeds; /* just a count of positive flags */
int             *ChrLoci; /* Contains selected marker loci on selected
			     chromosomes
			     For LoopOverTrait=0, list may have trait_loci as well */
int             NumChrLoci; /* Number of ChrLoci */
int             NumChrSite; /* Number of markers */
int             human_x, human_xy, human_y, human_mt, human_unknown, human_auto;
int             default_output_filenames;
/* Global flags. */
int             abortflag;
int             SetMarkerPosToSpecial;
int             force_numeric_alleles = 0;

int lastautosome    = 0;
int pseudoautosome  = 0;
int mitoautosome    = 0;

int SEX_CHROMOSOME  = -2;
int MALE_CHROMOSOME = -2;
int PSEUDO_X        = -2;
int MITO_CHROMOSOME = -2;

char SEX_CHROMOSOME_STR[4]  = ".x.";
char MALE_CHROMOSOME_STR[4] = ".y.";
char PSEUDO_X_STR[4]        = "x.y";
char MITO_CHROMOSOME_STR[4] = "...";
#if defined(_WIN) || defined(MINGW)
const char write_binary[] = "wb";
const char read_binary[]  = "rb";
#else
const char write_binary[] = "w";
const char read_binary[]  = "r";
#endif

SubjectOrganismType subject_organism;

/* global that control processing of the data */

int debug; /* for all the little things */
int found_path; /* for the makeped part */
int expire_year, expire_date, expire_mo; /* validation */
char months[12][10];
char awk_str[5];

int zero_female_ylinked = 0;

extern file_format check_locus_file_format(FILE *input_files);

/* set all global variables to NULL so that they can
   be freed safely */

static void    init_globals(char *argv0)
{
    int i;
    char *path_end, *p;
    char mega2_path1[FILENAME_LENGTH]=" ";

    /* set the mega2 path */
#ifdef _WIN
    ((long *)&QuantOutHi)[0] = 0xffffffff;
    ((long *)&QuantOutHi)[1] = 0x7fffffff;
    QuantOutLow = QuantOutHi;
#endif
    check_web_ver = 1;
    MARKER_SCHEME = MARKER_SCHEME_BITS;
    IgnoreValue(getcwd(InputPath, (size_t) FILENAME_LENGTH));
    strcpy(mega2_path1, argv0);
    path_end=strrchr(mega2_path1, '/');
    if (path_end != NULL) {
        p = &(mega2_path1[0]);
        i=0;
        while(p != path_end) {
            mega2_path[i]=mega2_path1[i];
            p++;
            i++;
        }
        mega2_path[i]='/';
        mega2_path[i+1]='\0';
    }
    else {
        strcpy(mega2_path,"");
    }

    for (i=0; i<NUMBER_OF_MEGA2_INPUT_FILES; i++) {
        mega2_input_files[i] = NULL;
    }
    global_chromo_entries=NULL;
    global_trait_entries=NULL;
    covariates=NULL;
    output_paths=NULL;
    Labels=NULL;
    NumLabels=0;
    HasQuant=HasAff=HasMarkers=0;
    Mega2OutputPath = NULL;
#if defined(SOLARIS)
    strcpy(awk_str, "nawk");
#elif defined(DARWIN_OS)
    strcpy(awk_str, "awk");
#else
    strcpy(awk_str, "gawk");
#endif
    UntypedPedOpt = -1;
    UntypedPeds = NULL;
    HalfTypedReset=-1;
    NonMendelianReset=-1;
    NOcTIME = "Fri Aug 29 02:14:00 1997\n";
    getRunDate();
    NukedMultiplier=1000;

    /* initialize seed */
    seed_random();
    seed2 = 23723;
    seed3 = 24724;

    strcpy(yorn[0], "no ");
    strcpy(yorn[1], "yes");
    strcpy(nory[0], "yes");
    strcpy(nory[1], "no ");
    CreateRunFolder=1;
    strcpy(Mega2Batch, "none");
    BatchFileCreated=0;
    Display_Errors=Display_Messages=1;
    ChrLoci= NULL;
    MissingIsNa=0;
    LinkageLocusFile=0;
    AllowUnmapped=0;
#ifdef MINGW
    strcpy(bindir, "c:/MinGW/msys/1.0/bin");
#else
    strcpy(bindir, "/bin");
#endif

    strcpy(Rplot_function,"nplplot.multi");
    strcpy(months[0], "January");
    strcpy(months[1], "February");
    strcpy(months[2], "March");
    strcpy(months[3], "April");
    strcpy(months[4], "May");
    strcpy(months[5], "June");
    strcpy(months[6], "July");
    strcpy(months[7], "August");
    strcpy(months[8], "September");
    strcpy(months[9], "October");
    strcpy(months[10], "November");
    strcpy(months[11], "December");
    debug = 0;
    Mega2Ver = MEGA2VER;
    Mega2Rev = MEGA2REV;
    Mega2Patch = MEGA2PATCH;
    /* This is a temporary fix,
       web comparison should still work since
       it compares the integer values. */
#ifdef HIDESTATUS
    strcpy(Mega2Version, "X.Y.Z");
#else
    sprintf(Mega2Version, "%d.%d.%d", Mega2Ver, Mega2Rev, Mega2Patch);
#endif
#ifdef HIDEPATH
    NOPATH = "...";
#else
    NOPATH = NULL;
#endif
    default_output_filenames=1;
    FirstIterMenu = 1;
    ChrLoci = NULL;
}

static void free_globals(void)

{
    int ii;

/*  superceded by Free_batch_items() */

    /* global_trait_entries */
    if (global_trait_entries != NULL)
        free(global_trait_entries);

    /* Free labels */
    if (NumLabels > 0) {
        free(Labels);
    }

    if (Mega2OutputPath != NULL)
        free(Mega2OutputPath);
/*
    for (ii=0; ii < NUMBER_OF_MEGA2_INPUT_FILES; ii++) {
	if (mega2_input_files[ii] != NULL)
	    free(mega2_input_files[ii]);
    }
*/
    for (ii=0; ii < NUM_OUTFILES; ii++) {
        /*    if (Outfile_Names[ii] != NULL)
              check not necessary, as they are always allocated */
        free(Outfile_Names[ii]);
    }
/*  superceded by Free_batch_items() */
    return;
}

/**
   Set the chromosome numbers that will be used by MEGA2.

   @param autosome the number of autosomes that the organism has
   @param xy  if == 0 then undefined; else if == 1 then autosome + 3; else the number given
   @param mito  if == 0 then undefined; if == 1 then autosome + 4; else the number given
 */
void set_chromosomes(const int autosome, const int xy, const int mito)
{
    lastautosome = autosome;

    SEX_CHROMOSOME  = autosome + 1;
    MALE_CHROMOSOME = autosome + 2;

    sprintf(SEX_CHROMOSOME_STR,  "%d", SEX_CHROMOSOME);
    sprintf(MALE_CHROMOSOME_STR, "%d", MALE_CHROMOSOME);

    if (xy) {
        PSEUDO_X        = xy == 1 ? autosome + 3 : xy;
        sprintf(PSEUDO_X_STR,        "%d", PSEUDO_X);
        pseudoautosome = xy;
    }
    if (mito) {
        MITO_CHROMOSOME = mito == 1 ? autosome + 4 : mito;
        sprintf(MITO_CHROMOSOME_STR, "%d", MITO_CHROMOSOME);
        mitoautosome = mito;
    }
}

void set_chromosomes_horse() { subject_organism = HORSE_SO; set_chromosomes(31, 0, 0); }
void set_chromosomes_sheep() { subject_organism = SHEEP_SO; set_chromosomes(26, 0, 0); }
void set_chromosomes_dog()   { subject_organism = DOG_SO; set_chromosomes(38, 1, 0); }
void set_chromosomes_mouse() { subject_organism = MOUSE_SO; set_chromosomes(19, 0, 0); }
void set_chromosomes_cow()   { subject_organism = COW_SO; set_chromosomes(29, 0, 0); }
void set_chromosomes_human() { subject_organism = HUMAN_SO; set_chromosomes(22, 1, 1); }

int exceeded_max_morgan_value(const double position) {
    switch (subject_organism) {
    case HUMAN_SO:
        return position > 5;
    default:
        return 0;
  }
}


/* -------------------------- here is main()----------------------------*/

extern void init_analysis();

extern int  db_exists_db();
extern void db_open_db();
extern void db_init_all();
extern void db_fini_all();
extern char DBfile[255];

int             main(int argc, char **argv, char **env)
{
    ped_top        *PedTreeTop = NULL;
    linkage_ped_top *LPedTreeTop = NULL;
    linkage_ped_top *Top2 = NULL;
    file_format     infl_type;
    analysis_type   analysis;
    char           *pedfl_name=NULL, *locusfl_name = NULL;
    int             it, ii;
    int             numchr;
    char           *omitfl_name = NULL;
    char           *freqfl_name = NULL;
    char           *mapfl_name  = NULL;
    char           *pmapfl_name = NULL;
    char           *penfl_name  = NULL;
    char           *bedfl_name  = NULL;
    char           *phefl_name  = NULL;
    char           *input_path  = NULL;
    char const     *logdir;
    int            num_cols=0;
/*  int            check_web_ver = 1; */

    FILE          *fp;
    int            ferr = 0; /* error opening one or more files*/

    plink_info_type *plink_info = (plink_info_type *)NULL;
    char          *zero = canonical_allele("0");
    
    Tod tod_all("total elapsed time");
    Tod tod_init("mega2 init");

    Mega2Status = INIT;
    /* Set the default modes/options/initializations */
    init_globals(argv[0]);
    TimeStampWritten[0]=0; TimeStampWritten[1]=0;

    mega2_opts(argc, argv);

//  ask IS READ BY DEFAULT
/*
    if (database_off)
        database_read = database_dump = 0;
    else if (database_dump || database_read) ;
    else if (db_exists_db())
        database_read = 1;
    else
        database_dump = 1;
 */
    init_analysis();
    // Initialize these just in case we are not getting the data from a batch file...

    genetic_distance_index = -1;
    base_pair_position_index = -1;
    genetic_distance_sex_type_map = UNKNOWN_GDMT;

#ifdef TEST
#undef EXPIRE
#endif

#ifdef EXPIRE
    check_expiration();
#endif
    logdir = mklogdir();
    LogFileNames();
    open_logs();
    hello(stdout);

    /* delete error, log and batch files after making copies,
       the user will be warned that copies were made when
       these files are created again (utils.c). */

    tod_init();
    Tod tod_batch("batch file init");
    batchfile_init_Mega2BatchItems();
    Mega2BatchItems[/* 52 */ Value_Marker_Compression].value.option = MARKER_SCHEME; // set on cmd line

    // Try to open the batch file if it's name is anything other than 'none'.
    // This signifies BATCH_FILE_INPUTMODE...
    if (*Mega2Batch != 0 && strcmp(Mega2Batch, "none")) {
        if (access(Mega2Batch, F_OK) != 0) {
            printf("ERROR: Cannot find batch file %s!\n", Mega2Batch);
            printf("       Terminating Mega2.\n");
            EXIT(FILE_NOT_FOUND);
        } else {
            if (check_empty(Mega2Batch)) {
                printf("ERROR: Batch file %s appears to be empty!\n", Mega2Batch);
                printf("       Terminating Mega2.\n");
                EXIT(FILE_READ_ERROR);
            }
            fp = fopen(Mega2Batch, "r");
            if (check_dos_file(fp) != 0) {
                printf("ERROR: Batch file %s is not a UNIX format file!\n",
                       Mega2Batch);
                printf("       Terminating Mega2.\n");
                EXIT(FILE_READ_ERROR);
            }
            fclose(fp);
        }

#ifdef HIDEFILE
        sprintf(err_msg, "Running Mega2 in batch mode.\n");
#else
	// DO NOT CHANGE: This line is parsed by 'utils.c:#define LOG2HTML'...
        sprintf(err_msg, "Running Mega2 in batch mode from %s.", Mega2Batch);
#endif
        mssgf(err_msg);

        batchfile_process(Mega2Batch, &analysis);

	// initialize globals form the batch file input if any...
        genetic_distance_index = Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option;
        base_pair_position_index = Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option;
        genetic_distance_sex_type_map = Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option;
//      MARKER_SCHEME = Mega2BatchItems[/* 52 */ Value_Marker_Compression].value.option;
        BatchValueIfSet(MARKER_SCHEME, "Value_Marker_Compression");
        if (marker_scheme_mega2_opts)
            MARKER_SCHEME = marker_scheme_mega2_opts;
        if (MARKER_SCHEME > MARKER_SCHEME_PTR || MARKER_SCHEME < MARKER_SCHEME_BITS) {
            printf("MARKER_SCHEME is %d.  Allowed values are 1, 2 or 3\n", MARKER_SCHEME);
            EXIT(OUTOF_BOUNDS_ERROR);
        }
//  printf("MARKER_SCHEME = %d\n", MARKER_SCHEME);

        InputMode=BATCH_FILE_INPUTMODE;
    } else {
        InputMode=INTERACTIVE_INPUTMODE;

        // Since we are in interactive mode, we need to setup for writing a batch file...
        strcpy(Mega2Batch, "MEGA2.BATCH");
        backup_file(&(Mega2Batch[0]));
    }

    if (AnalyInputMode == NOEXEC_INPUTMODE)
        AnalyInputMode = InputMode;

    if (database_off)
        database_read = database_dump = 0;
    else if (database_dump || database_read) ;
    else menu0();
    if (InputMode == BATCH_FILE_INPUTMODE)  // may have changed database_read/dump value
        check_batch_items();
    tod_batch();
    // determine if we should go out to the web and check to see if the user is running the latest release of MEGA2...
    if ((InputMode == BATCH_FILE_INPUTMODE && access(Mega2Batch, F_OK) == 0) ||
	!check_web_ver) {
        ;
    } else {
        mega2_version_check();
    }

    main_chromocnt = 1; numchr = 1;

    for (ii=0; ii<NUM_OUTFILES; ii++) {
        Outfile_Names[ii]=CALLOC((size_t)FILENAME_LENGTH, char);
    }

    /* init local variables */
    infl_type = InputFileFormat = LINKAGE;

    plink_info = CALLOC((size_t)1, plink_info_type);
/*
 Display the input menu, and return when the user enters option '0'... 
 
 ==========================================================
 Mega2 X.Y.Z Traditional input menu:
 ==========================================================
 0) Done with this menu - please proceed
 1) Input file extension:                   01
 2) Locus datafile:      (Mega2 loc)        names.01
 3) Pedigree datafile:   (Mega2 ped)        pedin.01
 4) Map datafile:        (Mega2 map)        map.01
 5) Omit datafile:       {Mega2 omit}       _
 6) Frequency datafile:  {annotated fmt}    _
 7) Penetrance datafile: {annotated fmt}    _
 8) Directory for writing output:           [ Current directory ]
 9) Simulate genotyping errors:             [no ]
 10) Include all pedigrees whether typed or not (based on fully_typed genotypes)
 11) Warn if allele frequency error measure exceeds: No limit
 12) Switch to PLINK input menu (binary format)
 13) Switch to PLINK input menu (ped format)
 Select from options 0-13 >
 */

    if (database_dump || ! database_read) {
        mega2_input_file_type[BED][0] = 0;  // just to be safe
        Tod tod_menu1("menu1");
        char *dbf = DBfile;
        menu1(&infl_type, &pedfl_name, &locusfl_name,
              &mapfl_name,  &pmapfl_name, &input_path, &omitfl_name,
              &freqfl_name, &penfl_name, &bedfl_name, &phefl_name,
              &UntypedPedOpt, &ErrorSimOpt, &Mega2OutputPath, &dbf,
              &FreqMismatchThreshold);
        tod_menu1();

        Input_Files& inf = Input->input_files;  // Input is set in menu1 as soon as possible.

        plink_info->plinkf = (Input_Format == in_format_binary_PED) ? binary_PED_format : 
                               (Input_Format == in_format_PED) ? PED_format : not_plink_format;

        *inf.pedfl   = pedfl_name;
        *inf.locusfl = locusfl_name;
        *inf.mapfl   = mapfl_name;
        *inf.pmapfl  = pmapfl_name;
        *inf.omitfl  = omitfl_name;
        *inf.freqfl  = freqfl_name;
        *inf.penfl   = penfl_name;
        *inf.bedfl   = bedfl_name;
        *inf.phefl   = phefl_name;

        strcpy(&mega2_input_file_type[PEDIGREE][0],  "Pedigree file");
        strcpy(&mega2_input_file_type[LOCUS][0],  "Locus file");
        strcpy(&mega2_input_file_type[MAP][0],  "Map file");
        strcpy(&mega2_input_file_type[PMAP][0],  "PLINK Map file");
        strcpy(&mega2_input_file_type[OMIT][0],  "Omit file");
        strcpy(&mega2_input_file_type[FREQ][0],  "Frequency file");
        strcpy(&mega2_input_file_type[PEN][0],  "Penetrance file");
        if (mega2_input_file_type[BED][0] == 0)
            strcpy(&mega2_input_file_type[BED][0],  "Aux file");
        strcpy(&mega2_input_file_type[PHEfl][0],  "PLINK Phenotype file");

        Mega2Status = FILE_NAMES_READ;
        time_stamp_logs();

        // Check to see if all of the files that were specified are available...
        for (ii = 0; ii < NUMBER_OF_MEGA2_INPUT_FILES; ii++) {
            if (mega2_input_files[ii]) {
                if ((fp=fopen(mega2_input_files[ii], "r")) == NULL) {
                    errorvf("Could not open %s (\"%s\") for reading!\n",
                            mega2_input_file_type[ii], mega2_input_files[ii]);
                    ferr += 1;
                } else
                    fclose(fp);
            }
        }
        if (ferr) EXIT(FILE_READ_ERROR);

        check_size_dos();
        log_line(mssgf);
    } else {
        Tod tod_menu1a("menu1a");
        char *dbf = DBfile;
        menu1a(&UntypedPedOpt, &ErrorSimOpt, &Mega2OutputPath, &dbf,
               &FreqMismatchThreshold);
        tod_menu1a();
    }

//    log_line(mssgf);
    /*
     Display the Analysis Menu, and return when the user enters a valid analysis number...
     
     ===========================================================
     ANALYSIS MENU 
     ==========================================================
      1 SimWalk2 format                 18 Linkage format                 
      2 Vintage Mendel format           19 Test loci for HWE              
      3 ASPEX format                    20 Allegro format                 
      4 GeneHunter-Plus format          21 MLBQTL format                  
      5 GeneHunter format               22 SAGE format                    
      6 APM format [DISBALED]           23 Pre-makeped format             
      7 APM-MULT format [DISABLED]      24 Merlin/SimWalk2-NPL format     
      8 Create nuclear families         25 PREST format                   
      9 SLINK format                    26 PAP format                     
     10 SPLINK format                   27 Merlin format                  
     11 Homogeneity analyses            28 Loki format                    
     12 SIMULATE format                 29 Mendel format                  
     13 Create summary files            30 SUP format                     
     14 Old SAGE format                 31 PLINK format                   
     15 TDTMax analyses [DISABLED]      32 CRANEFOOT format               
     16 SOLAR format                    33 Mega2 format         
     17 Vitesse format                  34 IQLS/Idcoefs format            
     
     Select an option between 1-34 > 
     */
    if (database_dump && ! database_read) {

        extern Missing_Value missing_value;
        extern Missing_Value missing_values[];
        extern int count_missing_values;

        analysis = DUMP;
        Mega2Status = ANALYSIS_NAME_READ;
        AnalysisOpt = analysis;
        if (InputMode == INTERACTIVE_INPUTMODE) {
            // Write the Analysis_Option
            strcpy(Mega2BatchItems[/* 5 */ Analysis_Option].value.name, analysis->_name);
            batchf(Analysis_Option);
        }

        analysis->prog_name(ProgName);
        mssgvf("Analysis Class: %s.\n", ProgName);
        for (int i = 0; i < count_missing_values; i++) {
            if (! strcmp(analysis->_name, (missing_values[i].analysis)->_name)) {
                missing_value = missing_values[i];
                break;
            }
        }

    } else {

        analysis_menu1(&analysis);
        Mega2Status = ANALYSIS_NAME_READ;
        AnalysisOpt = analysis;
#ifdef DARWIN_OS
        if (SIMWALK2(analysis)) {
            char *xmap = mega2_input_files[2] != 0 ? mega2_input_files[MAP] : mega2_input_files[PMAP];
            if ( !(xmap && strcasecmp(xmap, "map.dat")) &&
                !(strcasecmp(Mega2OutputPath, InputPath))) {
                char reply;
                warnf("Darwin file names are case-insensitive.");
                sprintf(err_msg, "Simwalk2 shell script will overwrite %s.",
                        mega2_input_files[2]);
                warnf(err_msg);
                printf("Terminate mega2 ? [y/n] (default 'y') > ");
                fcmap(stdin, "%c", &reply);
                if (tolower(reply) != 'n') {
                    EXIT(EARLY_TERMINATION);
                }
                fflush(stdin);
            }
        }
#endif
    }

    Value_Missing_get(&analysis);  // needed before file read

#ifndef HIDESTATUS
    int guess;
#endif
    if (database_dump || ! database_read) {
        int af;
        if (Input_Format == in_format_traditional) {
            af = check_annotated_file_format(mega2_input_files);
            if (af == 0) {
                if (locusfl_name && *locusfl_name != 0) {
                    FILE *tfp = fopen(locusfl_name, "r");
                    if (check_locus_file_format(tfp) == NAMES)
                        Input_Format = in_format_extended_linkage;
                    else
                        Input_Format = in_format_linkage;
                    fclose(tfp);
                } else {
                    errorvf("The input format is linkage and the locus file is missing.\n");
                    EXIT(EARLY_TERMINATION);
                }
            } else if ((af == 3) || (af == 4))
                Input_Format = in_format_mega2;
            msgvf("Input Format Deduced as: %s\n", INPUT_FORMAT_STR[Input_Format]);
    #ifndef HIDESTATUS
            guess = 1;
    #endif
        } else {
            msgvf("Input Format: %s\n", INPUT_FORMAT_STR[Input_Format]);
    #ifndef HIDESTATUS
            guess = 0;
    #endif
        }
    }

    extern linkage_ped_top SQLop;
    Tod tod_files("read all files");
    if (!database_dump && database_read) {
        extern void dbmega2_import(linkage_ped_top *Top);
        extern void check_both_pos_index();

        add_allele("NA", zero);
        REC_UNKNOWN = zero;

        Tod dbimport("db_import");

        db_open_db();
        if (InputMode == INTERACTIVE_INPUTMODE && BatchValueRead("DBfile_name"))
            batchf("DBfile_name");
        db_init_all();
        dbmega2_import(&SQLop);
        db_fini_all();

        dbimport();

        check_both_pos_index();

        LPedTreeTop = &SQLop;

        if (NumUnmapped > 0)
            unmapped_markers = CALLOC((size_t) NumUnmapped, int);

        global_chromo_entries=CALLOC((size_t)MaxChromo + 1, int);
        chromo_loci_count = CALLOC((size_t)MaxChromo + 1, int);

        NumChromo = get_chromosome_list(LPedTreeTop->LocusTop, global_chromo_entries,
                                        chromo_loci_count, 0, analysis);

        get_trait_list(LPedTreeTop->LocusTop, 1);

        log_line(mssgf);

//      Mega2OutputPath = strdup((char *)".");

    } else if (Input_Format == in_format_mega2) {
#ifndef HIDESTATUS
        if (mega2_input_files[3] == NULL) {
            msgvf("Pedigree, names and map file %s Mega2 format.\n",
                  guess ? "appear to be in" : "specified as");
        } else {
            msgvf("Pedigree, names, map and omit file %s Mega2 format.\n",
                  guess ? "appear to be in" : "specified as");
        }
        mssgf("Input files will be read in as Mega2 format files.");
#endif
        add_allele("NA", zero);
        REC_UNKNOWN = zero;

        LPedTreeTop = read_annotated_files(pedfl_name,  locusfl_name,  mapfl_name,
                                           freqfl_name, penfl_name,    omitfl_name,
                                           bedfl_name,  phefl_name,
                                           UntypedPedOpt, analysis, plink_info);
    } else if (Input_Format == in_format_binary_PED || Input_Format == in_format_PED) {
        it = PLINK_Args;
        if (Mega2BatchItems[it].items_read == 0) {
            errorvf("PLINK arguments not specified.\n");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }

#ifndef HIDESTATUS
        msgvf("Pedigree and map files %s PLINK format.\n",
              guess ? "appear to be in" : "specified as");
        mssgf("omit, penetrance, and frequency files are always in Mega2 format.");
        mssgf("Input files will be read in as PLINK or Mega2 format files as appropriate.");
#endif
        add_allele("NA", zero);
        REC_UNKNOWN = zero;

        LPedTreeTop = read_annotated_files(pedfl_name,  pmapfl_name,
                                           ((mapfl_name && *mapfl_name != 0) ? mapfl_name : pmapfl_name),
                                           freqfl_name, penfl_name,    omitfl_name,
                                           bedfl_name,  phefl_name,
                                           UntypedPedOpt, analysis, plink_info);
    } else if (Input_Format == in_format_binary_VCF || Input_Format == in_format_compressed_VCF ||
               Input_Format == in_format_VCF) {

        it = PLINK_Args;
        if (Mega2BatchItems[it].items_read)
            PLINK_args(Mega2BatchItems[it].value.name, 1);

        it = VCF_Args;
        if (Mega2BatchItems[it].items_read == 0) {
            errorvf("VCF arguments not specified.\n");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
	char *cp = CALLOC(FILENAME_LENGTH, char);
	strcpy(cp, Mega2BatchItems[it].value.name);
	if (Input_Format == in_format_binary_VCF) {
	  strcat(cp, " --bcf ");
	} else if (Input_Format == in_format_compressed_VCF) {
	  strcat(cp, " --gzvcf ");
	} else if (Input_Format == in_format_VCF) {
	  strcat(cp, " --vcf ");
	}
	strcat(cp, bedfl_name);

	if (VCFtools_process_cmd_line_if_necessary_inclusive(cp) != -1) {
	  errorvf("VCF arguments can not be processed.\n");
	  EXIT(BATCH_FILE_ITEM_ERROR);
	}
        free(cp);
        
        mssgf("\nProcessing VCF file meta information and header.");
        VCFtools_process_file_meta_information_and_header();
        
        
#ifndef HIDESTATUS
        mssgf("Pedigree and bim files appear to be in PLINK format.");
        mssgf("omit, penetrance, and frequency files are always in Mega2 format.");
#endif
        add_allele("NA", zero);
        REC_UNKNOWN = zero;
        
        LPedTreeTop = read_annotated_files(pedfl_name,  pmapfl_name, mapfl_name,
                                           freqfl_name, penfl_name, omitfl_name,
                                           bedfl_name,  phefl_name,
                                           UntypedPedOpt, analysis, plink_info);
    } else if (Input_Format == in_format_linkage || Input_Format == in_format_extended_linkage) {

#ifndef HIDESTATUS
        if (Input_Format == in_format_linkage) {
            msgvf("Pedigree, names and map file %s LINKAGE format.\n",
                  guess ? "appear to be in" : "specified as");
            mssgf("Input files will be read in as LINKAGE format files.");
        } else if (Input_Format == in_format_extended_linkage) {
            msgvf("Pedigree and map file %s LINKAGE format.\n",
                  guess ? "appear to be in" : "specified as");
            mssgf("Names (aka locus) file appears to be in Mega2 format w/o header.");
            mssgf("Input files will be read appropriately.");
        }

        if (penfl_name != NULL) {
            warnvf("Penetrance file %s will not be read in (only available in Mega2 format.\n", penfl_name);
        }
        if (freqfl_name != NULL) {
            warnvf("Frequency file %s will not be read in (only available in Mega2 format.\n", freqfl_name);
        }
#endif
        REC_UNKNOWN = zero;
        LPedTreeTop = read_linkage2(pedfl_name, locusfl_name,
                                    omitfl_name, mapfl_name,
                                    UntypedPedOpt, analysis);
    } else if (Input_Format < 8 ||
               Input->GetOps()->use_getops() ||
               Input_Format == in_format_traditional) { 

        Input->GetOps()->do_batch2local();

        add_allele("NA", zero);
        REC_UNKNOWN = zero;
        
        LPedTreeTop = read_annotated_files(pedfl_name,  pmapfl_name, mapfl_name,
                                           freqfl_name, penfl_name, omitfl_name,
                                           bedfl_name,  phefl_name,
                                           UntypedPedOpt, analysis, plink_info);
    } else {
        errorf("Input files appear to be in mixed Mega2 format and LINKAGE format.");
        errorf("Please use only Mega2 files or only LINKAGE files.");
        errorf("Unsuccessful in reading input files - aborting mega2!\n");
        EXIT(INPUT_DATA_ERROR);
    }
    if (LPedTreeTop == NULL) {
        errorvf("Unsuccessful in reading input files - aborting mega2!\n");
        EXIT(INPUT_DATA_ERROR);
    }

    tod_files();

    if (LPedTreeTop->pedfile_type == POSTMAKEPED_PFT) {
        infl_type = LINKAGE;
    } else {
        infl_type = PREMAKEPED;
    }

    Mega2Status = INPUT_FILES_READ;
    LPedTreeTop->analysis = analysis;

    if (database_dump || ! database_read) {
        Input->GetOps()->do_gc();
    } else {
        LPedTreeTop->pedfile_type = POSTMAKEPED_PFT;
        if (UntypedPeds == NULL) {
            UntypedPeds = CALLOC((size_t) LPedTreeTop->PedCnt, int);
        }
    }

    if (database_dump || ! database_read) {
        Tod tod_makeped("makeped");
        makeped(LPedTreeTop, analysis);  // if --db, might connect loops based on analysis
        tod_makeped();
    }
    if (database_dump) {
        LPedTreeTop->Ped = LPedTreeTop->PedRaw;  // largest set of persons
    } else {
        if ( (basefile_type == POSTMAKEPED_PFT && analysis->maintain_broken_loops()) ||
             (basefile_type != POSTMAKEPED_PFT && analysis->break_loops()) ) 
            LPedTreeTop->Ped = LPedTreeTop->PedBroken;
        else
            LPedTreeTop->Ped = LPedTreeTop->PedRaw;
    }

    LPedTreeTop->IndivCnt = 0;
    for (int ped = 0; ped < LPedTreeTop->PedCnt; ped++) {
        LPedTreeTop->IndivCnt += LPedTreeTop->Ped[ped].EntryCnt;
    }


    if (database_dump || ! database_read) {
        // Since the input mapfile may specify more than one genetic distance,
        // or base pair_position. Determine the appropriate one to use.
        // These get_* routines are found in the file user_input.c.
        //
        // fills genetic_distance_index and genetic_distance_sex_type_map; as a side effect...
        // see user_input.c:analysis_type RequiresGeneticlMap[]
        //
        // If these variables are -1 then they have not been specified by the batch file...
        if (LPedTreeTop->EXLTop == NULL)  {
            genetic_distance_index = -3;
            base_pair_position_index = -3;
        }

        get_genetic_distance_index(LPedTreeTop->EXLTop);
        // see user_input.c:analysis_type RequiresPhysicalMap[]
        get_base_pair_position_index(LPedTreeTop->EXLTop);
    }

    Tod tod_reorder("ReOrderLoci");
    /* Reorder the loci */
    LPedTreeTop = ReOrderLoci(LPedTreeTop, &numchr, &analysis);

    LPedTreeTop->LocusTop->SexLinked = x_linked_check(main_chromocnt, global_chromo_entries, analysis);
    /* Now set the default names of the output files */
    if (main_chromocnt >= 1) {
        // When the user selects more than one chromosome, this is the vector that holds their numbers.
        // Here we are choosing the first entry...
        numchr = global_chromo_entries[0];
    } else {
        numchr = 0;
    }
    tod_reorder();

    /* Convert to PedTree for checking purposes,
       don't need to assign affecteds */
    if (database_dump || ! database_read) {
        FirstTime=1;

        Tod tod_pedtree("convert_to_pedtree");
        PedTreeTop=convert_to_pedtree(LPedTreeTop, 0);
        /*   printf("Mega2Status = %d\n", Mega2Status); sleep(2);  */
        tod_pedtree();

    //NrvB: PedTreeTop->Locus ONLY contains markers!!
        Tod tod_fcheck("full_check");
        full_check(PedTreeTop, LPedTreeTop, analysis);
        /*   printf("Mega2Status = %d\n", Mega2Status); sleep(2);  */
        tod_fcheck();

        free_all_including_ped_top(PedTreeTop, NULL, NULL);
        /*   printf("Mega2Status = %d\n", Mega2Status); sleep(2);  */
        PedTreeTop = NULL;
        FirstTime=0;

    } else {
        extern int set_uniq_check(linkage_ped_top *LPedTop, analysis_type analysis);
        if (set_uniq_check(LPedTreeTop, analysis)) {
            extern void create_unique_ids(linkage_ped_top *Top, analysis_type analysis);
            create_unique_ids(LPedTreeTop, analysis);
        }

        extern void allelecnt_check(linkage_ped_top *Top, analysis_type analysis);
        if ((analysis == TO_PLINK || analysis == IQLS) /* && plink_locus_num < LTop->LocusCnt */) {
            allelecnt_check(LPedTreeTop, analysis);
        }

        extern void dbgenotype_import_genotype(linkage_ped_top *Top);
        Tod import_genotype("import_genotype");
        dbgenotype_import_genotype(LPedTreeTop);
        import_genotype();
    }

    if (true || database_dump || ! database_read) {
        Tod tod_stat1("write_ped_stat [again]");
        /* Output ped stats one more time */
        if (analysis != QUANT_SUMMARY) {
            log_line(mssgf);
            mssgf("Pedigree statistics after selecting chromosomes and marker loci:");
            write_ped_stats(LPedTreeTop);
        }
        /* Insert error simulation steps here */
        if (ErrorSimOpt) {
            simulate_errors(LPedTreeTop, numchr, Outfile_Names);
        }
        tod_stat1();
    }

    {
        /*  Mega2Status=TRAIT_SELECTED_M2S; */
        default_outfile_names(analysis, &(global_chromo_entries[0]), Outfile_Names, logdir);

        /*   printf("Mega2Status = %d\n", Mega2Status); sleep(2);  */
        set_output_paths(analysis, LPedTreeTop);

        /*   printf("Mega2Status = %d\n", Mega2Status); sleep(2);  */

        /* Now to analysis-specific options */
        ped_ind_defaults(LPedTreeTop->UniqueIds, analysis);
    }

    Mega2Status = INSIDE_ANALYSIS;
    /*  printf("num-traits = %d\n", num_traits); sleep(2); */

    // EXPORTS
    if (database_dump) {
        Tod dbexport("db_export");

        extern void dbmega2_export(linkage_ped_top *Top);
        db_open_db();
        if (InputMode == INTERACTIVE_INPUTMODE && BatchValueRead("DBfile_name"))
            batchf("DBfile_name");
        db_init_all();
        dbmega2_export(LPedTreeTop);
        db_fini_all();

        dbexport();
    }

    if (database_dump && database_read) {
        int i, j;
        char *name = argv[0];
        char **argvn = CALLOC((size_t) argc+6, char *);
        argvn[0] = argv[0];
        argvn[1] = (char *)"--dbread";
        argvn[2] = (char *)"--dbfile";
        argvn[3] = (char *)DBfile;
        j = 4;

        if (InputMode == BATCH_FILE_INPUTMODE)
            argvn[j++] = (char *)"--batch_file";
        else if (InputMode == INTERACTIVE_INPUTMODE)
            argvn[j++] = (char *)"--interactive";

        argvn[j++] = (char *)"--run_date";
        argvn[j++] = RunDate;

        for (i = 1; i < argc; i++) {
            if (strcasecmp(argv[i], "--dbdump") && strcasecmp(argv[i], "--dbread"))
                argvn[j++] = argv[i];
        }

        if (InputMode == INTERACTIVE_INPUTMODE)
            argvn[j++] = Mega2Batch;

        argvn[j++] = 0;
        int eans = 0;
        fflush(stdin);
        fflush(stdout);
        fflush(stderr);
        close_logs();
#if defined(_WIN) || defined(MINGW)
        eans = _spawnvpe(P_WAIT, name, argvn, env);
        exit(eans);
#else
        eans = execvp(name, argvn);
        printf("exec failed: eans = %d, errno %d\n", eans, errno);
        fflush(stdout);
#endif
    }

    InputMode = AnalyInputMode;  //What was it before the exec
    // Create the data files, and then the shell scripts...
    Tod tod_out("create_output_files");
    analysis->create_output_file(LPedTreeTop, 
                                 &analysis, Outfile_Names,
                                 UntypedPedOpt, &numchr, &Top2);
    tod_out();
    Tod tod_sh("create_shell_file");
    analysis->create_sh_file(LPedTreeTop, Outfile_Names, numchr);
    tod_sh();

    if (FirstIterMenu == 1 && InputMode == INTERACTIVE_INPUTMODE) {
	analysis->batch_out();
        if (analysis != SHAPEIT) {
            Mega2BatchItems[/* 25 */ Default_Outfile_Names].value.copt = 'y';
        }
	batchf(/* 25 */ Default_Outfile_Names);
    }
    Mega2Status = TERM_MEGA2;

    Tod tod_nuke_key("nuclear key");
    if (NUKE_OPTS) {
        if (Top2 != NULL) {
            write_nuc_ped_key_file(Mega2KeysRun, LPedTreeTop, Top2);
            free_lpedtop_pedinfo(Top2);
        }
    } else {
        write_key_file(Mega2KeysRun, LPedTreeTop);
    }
    tod_nuke_key();


    Tod tod_fin("epilogue");
#ifndef HIDEPATH
    if (database_read)
        msgvf("SQLite3 database \"%s\" was processed to generate this output.\n",
              DBfile);
#endif
    if (output_paths == 0) {
    } else if (!strcmp(output_paths[0], ".") &&
        ((LoopOverTrait == 1 && num_traits > 1) ||
         (LoopOverTrait == 0 && num_traits > 2 && (analysis == TO_SAGE && HasAff)))) {
        sprintf(err_msg,
                "Output is in %s, trait-specific files in subdirectories:",
#ifdef HIDEPATH
                NOPATH
#else
                output_paths[0]
#endif
            );
        mssgf(err_msg);
    } else if (!strcmp(output_paths[0], ".") &&
             ((LoopOverTrait == 1 && num_traits > 1) ||
              (LoopOverTrait == 0 && num_traits > 2 && (analysis == TO_SAGE && HasAff)))) {
        sprintf(err_msg, "Trait-specific output files are in subdirectories:");
        mssgf(err_msg);
    } else if (strcmp(output_paths[0], ".") &&
             (LoopOverTrait == 0 || num_traits == 1 || analysis != TO_SAGE)) {

        sprintf(err_msg, "Output is in %s",
#ifdef HIDEPATH
                NOPATH
#else
                output_paths[0]
#endif
            );
        mssgf(err_msg);
    }


    if ((LoopOverTrait == 1 && num_traits > 1) ||
       (LoopOverTrait == 0 && num_traits > 2 && (analysis == TO_SAGE && HasAff))) {
        char fl[FILENAME_LENGTH];
        sprintf(err_msg, "   ");
        ii = 1;
        for(it=1; it <= num_traits; it++) {
            int trait = global_trait_entries[it-1];
            if (trait < 0) continue;
            if (analysis->skip_trait(LPedTreeTop->LocusTop, trait)) continue;

            if (strcmp(output_paths[0], ".")) {
                shorten_path(output_paths[ii], fl);
            } else {
                strcpy(fl, output_paths[ii]);
            }
            /* Possible error: cast of a size_t to an int */
            num_cols += 2 + (int) strlen(fl);
            if (num_cols > 70 ) {
                mssgf(err_msg);
                sprintf(err_msg, "   ");
                num_cols=0;
            }
            strcat(err_msg, fl); strcat(err_msg, " ");
            ii++;
        }
        if (num_cols > 0) {
            mssgf(err_msg);
        }
    }

    log_line(mssgf);

    goodbye(0);
    Free_batch_items();
    Free_output_paths(analysis);
    Free_reorder_loci();
    Free_annotated_files(LPedTreeTop);
    free_globals();
    /* free local */

    if (pedfl_name != NULL)
        free(pedfl_name);

    if (locusfl_name != NULL)
        free(locusfl_name);

    if (mapfl_name != NULL)
        free(mapfl_name);

    if (omitfl_name != NULL)
        free(omitfl_name);

    if (freqfl_name != NULL)
        free(freqfl_name);

    if (penfl_name != NULL)
        free(penfl_name);

    if (bedfl_name != NULL)
        free(bedfl_name);

    if (phefl_name != NULL)
        free(phefl_name);

    tod_fin();

    tod_all();

    /*  free_all_including_lpedtop(&LPedTreeTop); */

#ifdef __INTEL__
    printf("Mega2 completed successfully, press any key to terminate Mega2 > ");
    fgetc(stdin); newline;
#endif

    return SUCCESS;
}  /* End of program! */
