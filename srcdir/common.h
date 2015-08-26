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

#ifndef COMMON_H
#define COMMON_H common_h

#define USE_PROTOS prototypes
#define ANSI ansi

#include "types.hh"

#ifdef _WIN
#include <io.h>
#include <direct.h>
#ifdef _WIN
#pragma warning(disable:4127)
#endif

#define F_OK 00
#define F_WR 02
#define F_RO 04
#define F_RW 06

/* apparently directories in windows are "special" */
#define W_OK 0

#define strcasecmp(s1,s2) _stricmp(s1,s2)
#define strncasecmp(s1,s2,len) _strnicmp(s1,s2,len)

#define mkdir(path,mode) _mkdir(path)
#define strdup(str) _strdup(str)
#define getcwd(buf,len) _getcwd(buf,len)
#define access(str,type) _access(str,type)
#define unlink(str) _unlink(str)

extern int snprintf(char *buf, int cnt, const char *fmt, ...);

/* warning C4701: potentially uninitialized local variable ... used */
#pragma warning(disable: 4100 4701)
#endif


#ifdef BORLAND
#include <alloc.h>
#endif


// UNKNOWN_POSITION is used for a position when a non-numeric or negative value is read, a missing value, or on initialization.
#define UNKNOWN_POSITION           -99.99

// INVALID_POS_DIFF would be used when a difference in positions would not be valid under certain assumptions.
// used in write_mendel7_files.c and write_mfiles.c
#define INVALID_POS_DIFF           -999.0

/* QUNDEF is a value of QMISSING which will be used to
   mark untyped individuals until we know the correct
   missing value */
#ifndef QUNDEF
//#define QUNDEF -999999
#define QUNDEF  (-9.999999999E30)
#endif

// QMISSING is used as the numeric equivalent of the string "NA" for missing quantitative phenotype values.
// It should an enumerated type equal to INT_MIN + 2 (or some other small int).
#ifndef QMISSING
//#define QMISSING -888888
#define QMISSING (-8.888888888E30)
#endif

#ifndef UNDEF
#define UNDEF -1
#endif

#ifndef REC_UNDEF
#define REC_UNDEF "-1"
#endif

#ifndef LINKAGE_UNKNOWN
#define LINKAGE_UNKNOWN "0"
#endif

#define System(str)    show_system_cmd(str, __FILE__, __LINE__)
#define EXIT(no)       Exit(no, __FILE__, __LINE__, #no)


#ifndef MALLOC
#define MALLOC(type)    ((type *) my_malloc(NULL,(size_t) 1,sizeof(type),__FILE__,__LINE__))
#endif

#ifndef MALLOC_STR
#define MALLOC_STR(n,type)    ((type *) my_malloc(NULL,(size_t) 1,n*sizeof(type),__FILE__,__LINE__))
#endif

#ifndef CALLOC
#define CALLOC(nelem,type) ((type *) my_calloc(NULL,nelem,sizeof(type),__FILE__,__LINE__))
#endif

#ifndef CALLOC_PTR
#define CALLOC_PTR(nelem,type) my_calloc(NULL,nelem,sizeof(type),__FILE__,__LINE__)
#endif

#ifndef REALLOC
#define REALLOC(ptr,n,t)     ((t *) my_realloc(ptr,n,sizeof(t),__FILE__,__LINE__))
#endif

/*
  #ifndef FREE
  #define FREE(p) { free(p); p = NULL }
  #endif

  #ifndef FLOAT
  #define FLOAT(x) ((float) (x))
  #endif
*/

#ifndef DOUBLE
#define DOUBLE(x) ((double) (x))
#endif

/*
  #ifndef INT
  #define INT(x) ((int) (x))
  #endif
*/

#ifndef LONG
#define LONG(x) ((long) (x))
#endif

#ifndef ODD
#define ODD(x)  ((x % 2)? 1 : 0)
#endif

#ifndef EVEN
#define EVEN(x) (!ODD(x));
#endif

#ifndef HALF
#define HALF(x) ((ODD(x))? (x+1)/2 : x/2);
#endif

#ifndef CR
#define CR  0x0D
#endif

#ifndef LF
#define LF  0x0A
#endif

#ifdef false
#undef false
#endif

#ifdef true
#undef true
#endif

#ifdef __MWERKS__
#define CODEWARRIOR
#define F_OK    0       /* Test for existence of File */
#endif

#ifdef _WIN
#define _ALLOW_KEYWORD_MACROS	// TRANSITION
#endif
#define false 0
#define true 1
#define EPSILON 0.000001    /* for floating point comparisons */
#define LARGE  999999999.0   /* To set upper limits */


#define INIT   0
#define FILE_NAMES_READ   1
#define ANALYSIS_NAME_READ   2
#define INPUT_FILES_READ   3
#define INSIDE_RECODE   4
#define DONE_RECODE   5
#define LOCI_REORDERED   6
#define TRAIT_SELECTED   7
#define INSIDE_ANALYSIS   8
#define TERM_MEGA2  9

#ifndef NOBOOLEAN
typedef int boolean;
#endif

typedef enum {
  UNKNOWN_SO = 0, HUMAN_SO = 1, HORSE_SO = 2, SHEEP_SO = 3, DOG_SO = 4, MOUSE_SO = 5, COW_SO = 6
} SubjectOrganismType;

typedef enum {
    POSTMAKEPED_PFT = 0, PREMAKEPED_PFT = 1
} PedigreeFileType;

typedef enum {
    PEDIGREE = 0, LOCUS = 1, MAP = 2, OMIT = 3, FREQ = 4, PEN = 5, BED = 6, PHEfl = 7, PMAP = 8
} InputFileType;

typedef enum {
  SEX_AVERAGED_GDMT = 0, SEX_SPECIFIC_GDMT = 1, FEMALE_GDMT = 2, NONE_AVAILABLE_GDMT = -2, UNKNOWN_GDMT = -1
} genetic_distance_map_type;

typedef enum {
  SEX_AVERAGED_MAP = 0, MALE_SEX_MAP = 1, FEMALE_SEX_MAP = 2
} sex_map_types;

typedef enum {
  BATCH_FILE_INPUTMODE = 0, INTERACTIVE_INPUTMODE = 1
} InputModeType;


// Used by write_ghfiles.c:write_recomb_fracs to determine the theta type to process...
typedef enum {
  MALE_THETA, FEMALE_THETA, SEX_AVERAGED_THETA
} theta_type;

typedef enum
{
    UNKNOWN, SL, ML, MULT, LINKAGE, MENDEL, ASPEX, CREATESUMMARY, LIABILITY,
    GENEHUNTER, NUKE, FREQUENCY, SLINK, SPLINK, GENEHUNTERPLUS, SIMULATE, SAGE, TDTMAX,
    SOLAR, HWETEST, VITESSE, PREMAKEPED, NAMES, ANNOTATED
} file_format;

typedef enum
{
    Premakeped, Postmakeped, Raw_premake, Raw_postmake, Annotated
} record_type;

typedef enum
{
    NOSEX_CHROMO_SLT = 0, SEX_NOAUTO_CHROMO_SLT = 1, SEX_AUTO_CHROMO_SLT = 2
} sex_linked_type;


typedef enum
{
    SAVE, FREE
} convert_type;

typedef enum
{
    FATAL, NON_FATAL
} error_type;

/* This will be used to define output file extensions where
   multiple chromosomes are being analyzed and multiple sets
   of files are produced.

   CHROMONUM : extensions are based on chromosome number,
   01 .. -9, 10 .. 23 which is how files are now named,
   and will remain the default.

   NUMBER : numbered 1 through N.
   EXTNUM : If single set of files, then user-provided extension
   is used, otherwise, we number extensions from 1 through N.
*/

/*  typedef enum */
/*  { */
/*    CHROMONUM, NUMBER, EXTNUM; */
/*  } outputfile_ext_type; */

#define INPUT  0
#define OUTPUT 1
#define MaxLinesPerScreen 15
/* MAXCHOICE is a global variable defining the number of output
 * options available within Mega2.
 */
#define MAXENTRY 500
#define MAXALLELES 200
#define MAX_NAMELEN 255
#define NAMELEN 25
#define FILENAME_LENGTH MAX_NAMELEN
#define NUM_OUTFILES 17
#define ALL_LEN 5
#define LIABILITY_LEN 50
#define MAX_CHROMO_NUM 100
#define MAXMARRIAGES 100
#define MAXLOOPMEMBERS 100
/* #define MAXPEDCNT 200 */
/*  #define MAXLOCUS MAXLOCI */
/*  #define FILENAME_LENGTH 100 */
#define NEWLNQ if (!quiet) putchar('\n')
#define LTOP Top->LocusTop
#define newline  printf("\n")
#define warn_unknown(opt)  printf("Unknown option %s.\n", opt)
#define NUMGENOS(n) n*(n+1)/2
#define TOGGLE(v) ((v)? 0 : 1)
#define LIABILITYMULTIPLIER 100000

/*Global vars*/

/* int             quiet; */
/* global variables that control output behaviour*/
extern int             MARKER_SCHEME;
extern const char     *genetic_distance_map_type_string[3];
extern int             LoopOverTrait; /* flag for looping over traits */
extern int             LoopOverChrm; /* flag for looping over chromosomes */
extern int             ManualReorder; /* flag denoting whether option 2 of reorder menu was used */
extern char            ProgName[100]; /* descriptive option name */
extern char            Mega2Log[FILENAME_LENGTH], Mega2LogRun[FILENAME_LENGTH]; /* Log file name */
extern char            Mega2Err[FILENAME_LENGTH], Mega2ErrRun[FILENAME_LENGTH]; /* Error file name */
extern char            Mega2Sim[FILENAME_LENGTH], Mega2SimRun[FILENAME_LENGTH]; /* Error sim log file name */
extern char            Mega2Batch[FILENAME_LENGTH]; /* Batch file name */
extern char            Mega2Keys[FILENAME_LENGTH], Mega2KeysRun[FILENAME_LENGTH]; /* ID-ped-per map */
extern char            Mega2Recode[FILENAME_LENGTH], Mega2RecodeRun[FILENAME_LENGTH]; /* Recoded markers */
extern char            Mega2Reset[FILENAME_LENGTH], Mega2ResetRun[FILENAME_LENGTH]; /* Reset genotypes */
extern const char      *NOPATH; /* what to print for file name iff in HIDEPATH mode */
extern const char      *NOcTIME;
extern char            Mega2Version[20];
extern int             Mega2Ver, Mega2Rev, Mega2Patch; /* Version */
extern char            Mega2WebVersion[50]; /* web-site file name */
extern char            err_msg[2*FILENAME_LENGTH]; /* array for error, warning and log messages */
extern char            **output_paths; /* trait directories prepended with output_dir */
extern char            **trait_paths;  /* list of directory names for each trait */
extern char            InputPath[FILENAME_LENGTH];
extern char            yorn[2][4]; /* two element list "yes", "no" */
extern char            nory[2][4]; /* two element list "no ", "yes" */
extern int             *Labels; /* affected labels (status-liability pairs) */
extern int             NumLabels; /* Number of Labels */
extern int             UntypedPedOpt; /* Which peds to omit based on typing */
extern int             liability_multiplier;
extern int             ErrorSimOpt; /* Flag for whether to introduce errors */

// Determine how to display the pedigree or individual using a 'linkage_ped_rec'.
// In some instances it is a name string, in others it is a number string (e.g., '1-5'),
// in the last instance it is the position at which it was encountered in the input file.
// OrigIds[1] pedigree...
// OrigIds[1]=={2|4} Ped[i].Name (string)
// OrigIds[1]=={3} i (int)
// OrigIds[1]=={5} Ped[i].Num (string)
// OrigIds[0] individual...
// OrigIds[0]=={1|2} OrigID (string)
// OrigIds[0]=={3|4} UniqueID (string)
// OrigIds[0]==?; ID (int)
extern int             OrigIds[2];

extern char           *Outfile_Names[NUM_OUTFILES];
extern int            CreateRunFolder;  /* create time-stamped folder for each run? */
extern int            check_web_ver;
extern int            BatchFileCreated; /* 1 if batchfile is created during run */
extern char           *REC_UNKNOWN;
extern int            MissingIsNa;

/*  outputfile_ext_type OutFileExtType; */
/*  char           **OutFileExts; */

/* globals to store progress */
extern int             TimeStampWritten[2];
extern int             Mega2Status;
extern char           *Mega2OutputPath;
extern int             FirstIterMenu;
/* globals that describe input data */
extern InputModeType   InputMode; /* Not batch mode? */
extern file_format     InputFileFormat; /* Annotated or linkage */
extern char            mega2_path[256]; /* path to mega2 executable */
#define NUMBER_OF_MEGA2_INPUT_FILES         9
extern char            *mega2_input_files[NUMBER_OF_MEGA2_INPUT_FILES];
extern char            mega2_input_file_type[NUMBER_OF_MEGA2_INPUT_FILES][21];
extern int             pedfile_type; /* whether input-file is pre-makeped or not */
extern int             HasLoops;
extern int             HasAff;
extern int             HasMarkers;
extern int             HasQuant;
extern int             LinkageLocusFile; /* added Dec 2007 */

extern int             MaxChromo; /* Largest chrom number */
extern int             NumChromo; /* Number of chromosomes in locus data */
extern int             *PedTreeMito; /* mitochondrial marker indexes when creating pedtree */
/* variables for storing markers without chromosome numbers */
extern int             NumUnmapped, *unmapped_markers, UnmappedSelected;
extern int             AllowUnmapped;

/* reordered_marker_loci = list of marker loci in user-specified or map order
   after locus reordering,
   num_reordered = length of reordered_marker_loci

   global_chromo_entries = chromosome numbers
   chromo_loci_count = number of loci on chromosome (cumulative)
   main_chromocnt = #chromosomes =
   length of global_chromo_entries and chromo_loci_count
*/
extern int             main_chromocnt, *global_chromo_entries, *chromo_loci_count;
extern int             num_reordered, *reordered_marker_loci;

/*  global _trait_entries = list of trait loci indexes
    may contain an entry -1 to denote position of markers
    num_traits = number of traits selected = length of global_trait_entries
    ditto num_covariates and covariates
*/
extern int             num_traits, num_affection, num_quantitative;
extern int             *global_trait_entries;
extern int             num_covariates, *covariates;


extern int             seed1, seed2, seed3; /* random seeds */
extern int             FirstTime; /* First time pedigrees are constructed */
extern double          MissingQuant;
extern double          MissingOutQuant, QuantOutLow, QuantOutHi;
extern int             MissingOutQuantSet;
extern int             NukedMultiplier;

extern int             SelectIndividuals;
extern int             HalfTypedReset; /* 1 if user resets half-typed genotypes */
extern int             NonMendelianReset; /* 1 if user resets inconsistent genotypes */
extern char            RunDate[FILENAME_LENGTH];
extern char            bindir[40];
extern double          FreqMismatchThreshold;

extern int             HasMapFileBeenRead;
extern int             HasFreqFileBeenRead;
extern int             HasPenFileBeenRead;
extern int             HasPedFileBeenRead;
extern int             HasLocFileBeenRead;
/* Global variables that control output data */

extern int             *UntypedPeds; /* vector of flags to indicate whether to leave out a pedigree */
extern int             NumTypedPeds; /* just a count of positive flags */
extern int             *ChrLoci; /* Contains selected marker loci on selected
                                    chromosomes
                                    For LoopOverTrait=0, list may have trait_loci as well */
extern int             NumChrLoci; /* Number of ChrLoci */
extern int             NumChrSite; /* Number of markers */
extern int             human_x, human_xy, human_y, human_mt, human_unknown, human_auto;
extern int             default_output_filenames;
/* Global flags. */
extern int             abortflag;
extern int             SetMarkerPosToSpecial;
extern int             force_numeric_alleles;
/* have to break up ASPEX into 4 options */

extern int lastautosome;

// The SEX_*, MALE_*, PSEUDO_*, MOTO_* variables are defined in utils.cpp
// along with the functions that set them are defined in mega2.cpp
extern int SEX_CHROMOSOME;
extern int MALE_CHROMOSOME;
extern int PSEUDO_X;
extern int MITO_CHROMOSOME;

extern char SEX_CHROMOSOME_STR[4];
extern char MALE_CHROMOSOME_STR[4];
extern char PSEUDO_X_STR[4];
extern char MITO_CHROMOSOME_STR[4];

extern SubjectOrganismType subject_organism;

extern const char write_binary[];
extern const char read_binary[];


#define UNKNOWN_CHROMO 999
#define MISSING_CHROMO 888
#define UNKNOWN_CHROMO_STR    "999"

/* global that control processing of the data */

extern int debug; /* for all the little things */
extern int found_path; /* for the makeped part */
extern int expire_year, expire_date, expire_mo; /* validation */
extern char months[12][10];

extern char awk_str[5];
#define   NLOOP       if (LoopOverTrait == 0 || num_affec < 2)  nloop=1; \
    else nloop = num_affec

#define skip_line(fp, lch)    lch=' '; fcmap(fp, "%=\n", lch)

#define SIMWALK2(analysis) ((analysis)->simwalk2())
#define QTL_DISALLOW ((analysis)->qtl_disallow())

/* This list is not the negation of the above, since
   some of the options do not take trait loci at all */

#define QTL_ALLOW ((analysis)->qtl_allow())
#define NUKE_OPTS ((analysis)->nuke_opts())

#ifndef advance_trp
#define advance_trp(i)                          \
    if (*i == -1) {                             \
        i++;                                    \
        continue;                               \
    }

#endif

#ifndef SKIP_TRI
#define SKIP_TRI(i)                                             \
    if (i>= 0) {if (global_trait_entries[i] == -1) continue; }
#endif

#ifndef SKIP_TRII
#define SKIP_TRII(i)                            \
    if (i == -1) continue
#endif

// (up to 10 displayed)

#ifndef SECTION_ERR
extern int Display_Errors, Display_Messages;

#define MAX_PED_ERRORS 9

#define SECTION_ERR(errors)  	        if (_##errors##_++ == 0) {         \
        fflush(stdout);                                                    \
        if (Display_##errors##_ != 0)                                      \
            errf("\n===== Errors/warnings of type \"" #errors "\": "); \
        Display_Errors = Display_##errors##_ = 1;                          \
        fflush(stdout);                                                    \
    } else if (_##errors##_ == (MAX_PED_ERRORS+2)) {                       \
        fflush(stdout);                                                    \
        errf("===== Too many \"" #errors "\" records, display is temporarily suspended ..", "", 0); \
        fflush(stdout);                                                    \
        Display_Errors = Display_##errors##_ = 0;                          \
    } else {                                                               \
        Display_Errors = Display_##errors##_;                              \
    }
#endif

#ifndef SECTION_ERR_INIT
#define SECTION_ERR_INIT(errors)     int Display_##errors##_ = 1, _##errors##_ = 0;
#endif

#ifndef SECTION_ERR_HEADER
#define SECTION_ERR_HEADER(errors)   Display_##errors##_ = _##errors##_ = 0;
#endif

#ifndef SECTION_ERR_EXTERN
#define SECTION_ERR_EXTERN(errors)   extern int Display_##errors##_, _##errors##_;
#endif

#ifndef SECTION_ERR_FINI 
#define SECTION_ERR_FINI(errors)	                                \
    {                                                                   \
        Display_Errors = Display_##errors##_ = 1;                       \
        fflush(stdout);                                                 \
        if (_##errors##_) errvf("===== %d total records of type \"" #errors "\" are in the ERR log.\n\n", _##errors##_); \
        fflush(stdout);                                                 \
        _##errors##_ = 0;                                               \
    }
#endif

#ifndef SECTION_ERR_COUNT
#define SECTION_ERR_COUNT(errors)  (_##errors##_)
#endif

#ifndef SECTION_ERR_FORCE
#define SECTION_ERR_FORCE(errors)  Display_Errors = 1;
#endif

#ifndef SECTION_LOG
#define SECTION_LOG(errors)  	        if (_##errors##_++ == 0) {         \
        fflush(stdout);                                                    \
        if (Display_##errors##_ != 0)                                      \
            mssgf("\n===== Messages of type \"" #errors "\": "); \
        Display_Messages = Display_##errors##_ = 1;                        \
        fflush(stdout);                                                    \
    } else if (_##errors##_ == (MAX_PED_ERRORS+2)) {                       \
        fflush(stdout);                                                    \
        mssgf("===== Too many \"" #errors "\" records, display is temporarily suspended ..", 0); \
        fflush(stdout);                                                    \
        Display_Messages = Display_##errors##_ = 0;                        \
    } else {                                                               \
        Display_Messages = Display_##errors##_;                            \
    }
#endif

#ifndef SECTION_LOG_INIT
#define SECTION_LOG_INIT(errors)     int Display_##errors##_ = 1, _##errors##_ = 0;
#endif

#ifndef SECTION_LOG_HEADER
#define SECTION_LOG_HEADER(errors)   Display_##errors##_ = _##errors##_ = 0;
#endif

#ifndef SECTION_LOG_EXTERN
#define SECTION_LOG_EXTERN(errors)   extern int Display_##errors##_, _##errors##_;
#endif

#ifndef SECTION_LOG_COUNT
#define SECTION_LOG_COUNT(errors)  (_##errors##_)
#endif

#ifndef SECTION_LOG_FINI 
#define SECTION_LOG_FINI(errors)	                                \
    {                                                                   \
        Display_Messages = Display_##errors##_ = 1;                     \
        fflush(stdout);                                                 \
        if (_##errors##_) msgvf("===== %d total records of type \"" #errors "\" are in the log.\n\n", _##errors##_); \
        fflush(stdout);                                                 \
        _##errors##_ = 0;                                               \
    }
#endif

#define ALLOW_NO_CHR(analysis)  ((analysis)->allow_no_chr())
#define ALLOW_NO_MAP(analysis)  ((analysis)->allow_no_map())
#define ALLOW_SEX_MAP(analysis)  ((analysis)->allow_sex_map())

#define Ignore_Unmapped(m) (((m) == UNKNOWN_CHROMO) && (NumUnmapped < 1 || !AllowUnmapped))

#endif
