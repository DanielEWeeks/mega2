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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdint.h>
#include <fcntl.h>
#include <errno.h>
#ifdef  USE_GETOPT
#include <getopt.h>
#endif

#include "common.h"
#include "typedefs.h"

#include "batch_input_ext.h"
#include "compress_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "net_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "version.h"
/*
     error_messages_ext.h:  close_logs errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
                net_ext.h:  
         user_input_ext.h:  untyped_ped_messg
              utils_ext.h:  EXIT chmod_X_file print_mega2_help print_mega2_version
*/

#ifdef _WIN
#define popen(cmd,mode) _popen(cmd,mode)
#define pclose(fd) _pclose(fd)
#define read(fd,buf,len) _read(fd,buf,len)
#define write(fd,buf,len) _write(fd,buf,len)
#define open _open
#define close _close

#define O_CREAT _O_CREAT
#define O_TRUNC _O_TRUNC
#define O_BINARY _O_BINARY
#define O_RDONLY _O_RDONLY
#define O_WRONLY _O_WRONLY

#elif !(defined(__CYGWIN__)) && !(defined(_WIN32))
#define O_BINARY 0
#endif


#ifdef TEST
#undef EXPIRE
#define SETSEED
#endif

#ifdef HIDEDATE
#define THEDATE "1997-8-29-02-14"
#endif

#define LOG2HTML          "mega2log2html.pl"

/* utils.c

Contents: general utility functions and file utility functions */

/* prototypes */

void           invalid_value_field(int batch_item);
void           unknown_prog(char *prog_name);
void           TCL_help_information(void);
double         factorial(int n);
double         choose(int n, int k);
void           draw_line(void);
void           log_line(void (*log_func)(const char *messg));
int            press_return(void);
void           missing_arg(char *arg);
void           missing_param(char *arg);
void           string_to_upper(char *strg);
float          safe_divide(int a, int b);
void           hello(FILE *fp);
const char*    mklogdir(void);
void           summary_time_stamp(char **input_files, FILE *fp,
				  const char *count_option);
void           script_time_stamp(FILE *script_fp);
double         randomnum( void );
void           seed_random(void);
void           write_time(FILE *fil);
int            imax(int i1, int i2);
double         dmax(double d1, double d2);
int            icomp(int i, int j);
void           exclaim(void);
void           delete_file(const char *flname);
void           copy_file(char *flname1, char *flname2);
void           move_file(char *flname1, char *flname2);
void           makedir(char *dirname);
void           move_logs(char *sumdir) ;
void           goodbye(int exit);
char *         strtail(const char *name, const int n);
const char    *perl_pgm(const char *pl);
static void    append_to_fd(const char *file_name, const char *str, FILE *FDo);
#ifdef EXPIRE
void           check_expiration(void);
#endif

#define MULTI_USE_BUFFER 4096
static         char multi_use_buffer[MULTI_USE_BUFFER];


/*   end prototypes */

/*
 This is the compare function used by the hash table that contains the names.
 The arguments actually passed to the compare function when it is called are
 addresses of the element pointers added to the set. E.g.: If you add char *
 pointers to the set, the compare function will be called with char ** pointers
 as arguments.
 */
int hashtable_strcmp(const void *prev, const void *next) {
    char *s1 = *(char **)prev;
    char *s2 = *(char **)next;
    int val = strcmp(s1, s2);
    return val;                                                                                             
}                                                                                                   

void invalid_value_field(int batch_item)
{
    errorvf("Keyword %s has invalid value in batch file.\n",
            Mega2BatchItems[batch_item].keyword);
    EXIT(BATCH_FILE_ITEM_ERROR);
}

void unknown_prog(char *prog_name)
{
    errorvf("Analysis option %s does not match a known option.\n", prog_name);
    invalid_value_field(5);
}

/* divides an integer with another, returns 0 if either of them is 0 */

float safe_divide(int a, int b)
{
    float result;
    if (b==0)
        result=0;
    else
        result=(float)a/(float)b;
    return result;
}

/* TCL_help_information will be a help database for
   information regarding the TCl parameter file */

void            TCL_help_information(void)
{
    /* just a stub now */
}

double factorial(int n)
/* Based on Numerical Recipes */
{
    static double table[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
    int             k;
    static int      l = 4;
    if (n > 32) {
        errorvf("factorial function only works for n < 32.\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }
    if (n < 0) {
        errorvf("called factorial function with a negative argument.\n");
        EXIT(OUTPUT_FORMAT_ERROR);
    }
    while (l < n) {
        k = l++;
        table[l] = table[k] * l;
    }
    return table[n];
}   /* end of factorial function */

double          choose(int n, int k)
{
/* double           factorial(int n); */

    if (n < k)   return 0;
    return floor(0.5 + factorial(n) / (factorial(k) * factorial(n - k)));
}


void            draw_line(void)
{
    fflush(stdout);
    printf("==========================================================\n");
    return;
}

void            log_line(void (*log_func)(const char *messg))
{
    fflush(stdout);
    (log_func)("===========================================================");
    return;
}

int            press_return(void)
{
    char            nextpage;
    printf(" > press 'RETURN' key for next page, 'q' to end display < ");

    nextpage = (char)getc(stdin);

    if (nextpage == 'Q' || nextpage == 'q') {
        fcmap(stdin, "%=\n", &nextpage);
        return 1;
    } else if (nextpage == '\n') {
        return 0;
    } else {
        ungetc(nextpage, stdin);
        return -1;
    }
}
/*---------------------------------------------------------------+
  | Inform the user of a missing required argument.               |
  | Request to abort.                                             |
  +---------------------------------------------------------------*/
void            missing_arg(char *arg)

{
    fprintf(stderr, "ERROR: Missing required argument %s.\n", arg);
    abortflag = 1;
}


/*---------------------------------------------------------------+
  | Inform the user that a program argument requires a parameter  |
  | and one was not supplied. Request to abort.                   |
  +---------------------------------------------------------------*/
void            missing_param( char           *arg)
{
    fprintf(stderr, "ERROR: Missing required parameter of argument %s.\n", arg);
    abortflag = 1;
}


/*---------------------------------------------------------------+
  | Convert a string to all upper case.                           |
  +---------------------------------------------------------------*/
void            string_to_upper(char           *strg)
{
    for (; *strg; strg++)
      *strg = (char)toupper((unsigned char)*strg);
}



#ifdef EXPIRE
/* License expiration checking code */

void check_expiration(void)
{

    /* Check  date against the __DATE__ string */
    /* will have a compilation tag that says beta */
    /* extern   time_t  time(); */
    time_t time_val;
    char v_date[14], v_month[14], v_year[5];
    char tt[50];
    int j, curr_year,  curr_date,  curr_mo, expired=0;

    expire_year = EXPIRE_YEAR;
    expire_mo   = EXPIRE_MO; /* June */
    expire_date = EXPIRE_DATE;

    time_val = time((time_t *) NULL);
    strcpy(tt, ctime(&time_val));

    for (j=0;j<3;j++)  {
        v_month[j]=tt[j+4];
    }
    v_month[j]='\0';
    for (j=0;j<2;j++)
        {
            v_date[j]=tt[j+8];
        }
    v_date[j]='\0';
    curr_date = atoi(v_date);


    for (j=0;j<4;j++)  {
        v_year[j]=tt[j+20];
    }
    v_year[j]='\0';
    curr_year = atoi(v_year);


    /* find the dist month [0 , 11} to match with curr_mon */
    curr_mo=-1;
    for (j=0; j < 12; j++) {
        if (v_month[0] == months[j][0] &&
           v_month[1] == months[j][1] &&
           v_month[2] == months[j][2]) {
            curr_mo = j;
            break;
        }
    }


    if ((curr_year < 2003) || (curr_year == 2003 && curr_mo < 4)) {
        errorvf("Date appears to be %s-%d, please check system clock!\n",
                months[curr_mo], curr_year);
        EXIT(SYSTEM_ERROR);
    }

    if (curr_year < expire_year) {
        expired=0;
    } else {
        if (curr_year == expire_year) {
            if (curr_mo < expire_mo) {
                expired = 0;
            } else {
                if (curr_mo == expire_mo) {
                    if (curr_date <= expire_date) {
                        expired=0;
                    } else {
                        expired=1;
                    }
                } else {
                    expired=1;
                }
            }
        } else {
            expired=1;
        }
    }

    if (expired) {
        errorvf("MEGA2 %s has expired, please obtain an up-to-date version.\n",
		Mega2Version);
        EXIT(SYSTEM_ERROR);
    }
}
#endif

void chomp(char *str)

{
    char c;

    if (!strcmp(str, "")) {
        return;
    }
    c = str[strlen(str) - 1];
    /* remove newline characters from end of str */
    if (c == CR || c == LF) {
        str[strlen(str) - 1] = ' ';
    }
}

/*---------------------------------------------------------------+
  | Just print out some program identification too.               |
  +---------------------------------------------------------------*/
char sumdir[FILENAME_LENGTH+30];

void            hello(FILE *fp)
{
#ifdef EXPIRE
    extern const char *gcc__VERSION__;
    extern const char *gcc_compile__DATE__;
    extern const char *gcc_compile__TIME__;
    char valid_str[20];
#endif

    fprintf(fp, "==========================================================\n");
    fprintf(fp, "                          MEGA2 %s\n", Mega2Version);
    fprintf(fp, "     Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,\n");
    fprintf(fp, "     Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,\n");
    fprintf(fp, "     Daniel E. Weeks, and University of Pittsburgh\n\n");
#ifndef HIDEFILE
#ifdef EXPIRE
    fprintf(fp, "     Last updated: %s, %s ", gcc_compile__DATE__,  gcc_compile__TIME__);
    sprintf(valid_str, "%s %d, %d", months[expire_mo], expire_date, expire_year);
    fprintf(fp, ", valid until %s.\n", valid_str);
    fprintf(fp, "     Compiled with gcc version %s\n\n", gcc__VERSION__);
#else
    fprintf(fp, "\n");
#endif
#endif
    fprintf(fp, "     Mega2 comes with ABSOLUTELY NO WARRANTY.\n");
    fprintf(fp, "     See LICENSE.txt for terms of copying, modifying & redistributing Mega2.\n");
    fprintf(fp, "==========================================================\n");
    fprintf(fp, "NOTE: If you have previously used explicit numbers for sex chromosomes, BEWARE!\n");
    fprintf(fp, "We have changed the numbers to be compatible with PLINK. 23 still codes for X,\n");
    fprintf(fp, "but 24 codes for Y and 25 Codes for XY.\n\n");

}

const char *mklogdir(void)
{
  // CreateRunFolder is the switch for creating a time-stamped folder for each run.
  // This is the default. Turned off with the -nosave option
    if (CreateRunFolder == 1) {
        strcpy(sumdir, RunDate);
        makedir(sumdir);
    } else {
        strcpy(sumdir, ".");
    }
    return (const char *)sumdir;
}

#ifdef _WIN
int snprintf(char *buf, int cnt, const char *fmt, ...)
{
    va_list ap;
    int ret;

    va_start(ap, fmt);
    ret = vsnprintf(buf, cnt, fmt, ap);
    va_end(ap);

    return ret;
}
#endif

void goodbye(int exit)
{
    FILE *XX, *err;
#ifndef VALGRIND
#ifndef _WIN
    FILE *pfd;
#endif
#endif
    char syscmd[200];
    char *fl_name;
#ifndef HIDESTATUS
    int exit_status;
#endif
/*
    char Mega2BatchRun[FILENAME_LENGTH+30];

    sprintf(Mega2BatchRun, "%s.%s", Mega2Batch, RunDate);
    copy_file(Mega2Batch, Mega2BatchRun);
 */

    if (CreateRunFolder == 1) {
        if (BatchFileCreated == 1) {
            sprintf(err_msg,
                    "Run parameters stored in batch file %s/MEGA2.BATCH",
#ifdef HIDEPATH
                    ""
#else
                    sumdir
#endif
                );
            mssgf(err_msg);
        }

        sprintf(err_msg, "See run summaries in directory %s ",
#ifdef HIDEPATH
                ""
#else
                sumdir
#endif
            );
        mssgf(err_msg);

    } else {
        strcpy(sumdir, ".");
        if (BatchFileCreated == 1) {
            mssgf("Run parameters stored in batch file MEGA2.BATCH");
        }

        sprintf(err_msg, "See run summaries in current directory .");
        mssgf(err_msg);
    }

    sprintf(err_msg, "   %s", Mega2Log);
    if (access(Mega2RecodeRun, F_OK) == 0) {
        strcat(err_msg, ", "); strcat(err_msg, Mega2Recode);
    }
    if (access(Mega2ResetRun, F_OK) == 0) {
        strcat(err_msg, ", "); strcat(err_msg, Mega2Reset);
    }
    if (access(Mega2ErrRun, F_OK) == 0)   {
        strcat(err_msg, ", "); strcat(err_msg, Mega2Err);
    }
    if (access(Mega2KeysRun, F_OK) == 0)  {
        strcat(err_msg, ", "); strcat(err_msg, Mega2Keys);
    }
    if (ErrorSimOpt == 1) {
        strcat(err_msg, ", "); strcat(err_msg, Mega2Sim);
    }

    /*  printf("\n"); */
    mssgf(err_msg);

    fl_name=CALLOC(50+strlen(sumdir), char);

    sprintf(fl_name, "%s/MEGA2.BATCH.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2.ERR.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2.KEYS.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2.LOG.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2file.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2links.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2outputfiles.html", sumdir);
    delete_file(fl_name);
    sprintf(fl_name, "%s/MEGA2run.html", sumdir);
    delete_file(fl_name);

    delete_file("__tmp__");
    sprintf(syscmd, "%s %s %s %s > __tmp__\n",
            perl_pgm(LOG2HTML),
            Mega2LogRun,
            (output_paths == NULL || output_paths[0] == NULL || output_paths[0] == 0) ? "." : output_paths[0],
            sumdir);
//    printf("syscmd: %s\n", syscmd);
//    exit(0);
    // returns the exit status of the shell as returned by waitpid(2), or
    // -1 if an error occurred when invoking fork(2) or waitpid(2)
    move_logs(sumdir);
    // BIG NOTE: this script assumes that MEGA2.BATCH is found in 'sumdir'.
#ifndef HIDESTATUS
    exit_status = System((const char *)syscmd);
#if defined(_WIN) || defined(MINGW) || (! defined(WIFEXITED))
#define WIFEXITED(exit_status) (exit_status & 0x7f)
#define WEXITSTATUS(exit_status) ((exit_status>>8) & 0xff)
#endif
    if (WIFEXITED(exit_status)) {
        if (WEXITSTATUS(exit_status) == 127) {
            // 127 means the execution of the shell failed
            warnvf("Execution of the shell failed for '%s'.\n", LOG2HTML);
            warnvf("No 'MEGA2*.html' files will have been created.\n");
        } else if (WEXITSTATUS(exit_status) == 0) {
            // 0 is normal (default) exit status for a perl script...
            mssgvf("The script '%s' exited normally.\n", LOG2HTML);
        } else {
            warnvf("The script '%s' exited with a status of %d.\n", LOG2HTML, WEXITSTATUS(exit_status));
            warnvf("This indicates that the script encountered a problem.\n");
        }
    } else {
        warnvf("An error occurred while attempting to run '%s'.\n", LOG2HTML);
    }
#else
    (void) System((const char *)syscmd);
#endif
    err_or_warn(&XX, &err);
    sprintf(syscmd, "Can not find '%s' run log.", LOG2HTML);
    append_to_fd("__tmp__", syscmd, err);
/*
    sprintf(syscmd, "cat __tmp__ >> %s", Mega2ErrRun);
    System(syscmd);
*/
    delete_file("__tmp__");
    if (access(fl_name, F_OK) == 0) {
#ifndef HIDESTATUS
      msgvf("To view the HTML-formatted run summaries, open\n%s/%s\nin a web browser.\n", InputPath, fl_name);
#endif
    } else {
#ifndef HIDESTATUS
        warnf("Unable to create html-formatted run summary.");
#endif
#ifndef VALGRIND
#ifndef _WIN
        sprintf(syscmd, "which %s", LOG2HTML);
        pfd = popen(syscmd, "r");
        syscmd[0] = '0';
        (void)fgets(syscmd, sizeof (syscmd), pfd);
        pclose(pfd);
        if (syscmd[0] == '0')
#endif
#endif
        warnvf("The program '%s' may not be installed on your system.\n", LOG2HTML);
	warnvf("Run 'install.sh scripts' to install it.\n");
    }
    free(fl_name);

    log_line(mssgf);
    if (exit == 0) {
        mssgf("If you use Mega2 as part of a published work, please reference ");
        /*
          mssgf(" Mukhopadhyay N, Almasy L, Schroeder M, Mulvihill WP, Weeks DE (2005)");
          mssgf(" Mega2: data-handling for facilitating genetic linkage and association analyses.");
          mssgf(" Bioinformatics. 2005 May 15;21(10):2556-7, PMID: 15746282");
        */
        mssgf(" Baron RV, Kollar C, Mukhopadhyay N, Weeks DE");
        mssgf(" Mega2: validated data-reformatting for linkage and association analyses");
        mssgf(" Source Code for Biology and Medicine.2014, 9:26");
        mssgf(" DOI: 10.1186/s13029-014-0026-y");

        sprintf(err_msg, "as well as the version used, which is currently Version %s",
                Mega2Version);
        mssgf(err_msg);

        log_line(mssgf);
    }
    close_logs();
/*
    delete_file(Mega2BatchRun);
 */
}

static void append_to_fd(const char *file_name, const char *str, FILE *FDo)
{
    size_t lenr, lenw;
    FILE *FDi = fopen(file_name, "r");
    if (FDi == NULL) {
        warnvf("%s (%s)\n", str, file_name);
        return;
    }
    while (1) {
        lenr =  fread(multi_use_buffer, sizeof (multi_use_buffer), 1, FDi);
        if (lenr != sizeof (multi_use_buffer) && ferror(FDi)) {
            switch(errno) {
            default:
                warnvf("append_to_fd: fread(%s, , %d, 1) = %d failed with errno %d (\"%s\")\n",
		       file_name, ((int)sizeof (multi_use_buffer)), lenr,
		       errno, strerror(errno));
                break;
            }
        }
        lenw = fwrite(multi_use_buffer, lenr, 1, FDo);
        if (lenr != lenw && ferror(FDo)) {
            switch(errno) {
            default:
                warnvf("append_to_fd: fwrite(%s, , %d) = %d failed with errno %d (\"%s\")\n",
		       "ErrLog", lenr, lenw, errno, strerror(errno));
                break;
            }
        }
        /*
        warnvf("append_to_fd: fread(%s) %d/fwrite %d\n",
                file_name, lenr, lenw);
        */
        if (feof(FDi)) break;
    }
    fclose(FDi);
}

void   summary_time_stamp(char **input_files, FILE *fp, const char *message)
{

    char mssg[100];
    int i;

    fprintf(fp, "-----------------------------------------------------\n");
    fprintf(fp, "        Mega2 version %s\n", Mega2Version);
#ifdef HIDEDATE
    fprintf(fp, "Run date:                  %s\n", THEDATE);
#else
    fprintf(fp, "Run date:                  %s\n", RunDate);
#endif
    fprintf(fp, "This file created on       ");
    write_time(fp);
    if (input_files != NULL) {
        fprintf(fp, "Input file names\n");

        for (i=0; i< 3; i++) {
            if (input_files[i] != NULL) {
                // DO NOT CHANGE: This line is parsed by '#define LOG2HTML'...
                fprintf(fp, "# %19s:               %s\n",
                        mega2_input_file_type[i],
#ifdef HIDEPATH
                        NOPATH
#else
                        input_files[i]
#endif
);
            }
#ifdef HIDEPATH
            else {
                // DO NOT CHANGE: This line is parsed by '#define LOG2HTML'...
                fprintf(fp, "# %19s:               %s\n",
                        mega2_input_file_type[i], NOPATH);
            }
#endif
        }
#ifndef HIDEPATH
        for (i=3; i< NUMBER_OF_MEGA2_INPUT_FILES; i++) {
            if (input_files[i] != NULL) {
                // DO NOT CHANGE: This line is parsed by '#define LOG2HTML'...
                fprintf(fp, "# %19s:               %s\n",
                        mega2_input_file_type[i],
                        input_files[i]
);
            }
        }
#endif

    }
    if (UntypedPedOpt != -1)
        fprintf(fp, "Untyped pedigree option: %s\n",
                untyped_ped_messg(UntypedPedOpt, mssg));
/*
    if (NonMendelianReset == 1) {
        fprintf(fp,
                "Mendelianly-inconsistent genotypes set to unknowns.\n");
    } else if (HalfTypedReset == 0) {
        fprintf(fp, "Mendelianly-inconsistent genotypes included in output.\n");
    }
    if (HalfTypedReset == 1) {
        fprintf(fp, "Half-typed individuals set to unknowns.\n");
    } else if (HalfTypedReset == 0) {
        fprintf(fp, "Half-typed individuals' genotypes included in output.\n");
    }
*/
    if (strcmp(message, ""))
        fprintf(fp, "%s\n", message);
    fprintf(fp, "---------------------------------------------\n\n");
}

void script_time_stamp(FILE *script_fp)

{
#ifndef HIDESTATUS
    char mssg[100];
#endif
    int i;

    /* mega2_input_files is a global */
    fprintf(script_fp, "#----------------------------------------------\n");
    fprintf(script_fp, "#   Mega2 version %s\n", Mega2Version);
#ifdef HIDEDATE
    fprintf(script_fp, "#   Run date:                %s\n", THEDATE);
#else
    fprintf(script_fp, "#   Run date:                %s\n", RunDate);
#endif
    fprintf(script_fp, "#   This script created on   ");
    write_time(script_fp);
    fprintf(script_fp, "#   Input file names:\n");

    for (i=0; i< 3; i++) {
        if (mega2_input_files[i] != NULL) {
            fprintf(script_fp, "# %19s:              %s\n",
                    mega2_input_file_type[i],
#ifdef HIDEPATH
                    NOPATH
#else
                    mega2_input_files[i]
#endif
);
        }
#ifdef HIDEPATH
        else {
            fprintf(script_fp, "# %19s:              %s\n",
                    mega2_input_file_type[i], NOPATH);
        }
#endif
    }
#ifndef HIDEPATH
    for (i=3; i< NUMBER_OF_MEGA2_INPUT_FILES; i++) {
        if (mega2_input_files[i] != NULL) {
            fprintf(script_fp, "# %19s:              %s\n",
                    mega2_input_file_type[i],
                    mega2_input_files[i]
);
        }
    }
#endif
    /*   fprintf(script_fp, "#    Pedigree file:           %s\n",  */
    /* 	  mega2_input_files[0]); */
    /*   fprintf(script_fp, "#    Map file:                %s\n", mega2_input_files[2]); */


    /*     fprintf(script_fp, "#    Omit file:               %s \n", mega2_input_files[3]); */
    /*   if (mega2_input_files[4] != NULL)  */
    /*     fprintf(script_fp, "#    Frequency file:          %s \n", mega2_input_files[4]); */

#ifndef HIDESTATUS
    fprintf(script_fp, "#    Untyped pedigree option  %s\n",
            untyped_ped_messg(UntypedPedOpt, mssg));
    fprintf(script_fp, "#----------------------------------------------\n");
#endif
    return;

}

void fcopy(const char *file1, const char *file2, const char *write_mode)

{
    /* copies the contents of file1 into file2. Write_mode can be
       "a" for append, or "w" for overwrite */

    FILE *f1, *f2;
    char line[100];

    f1=fopen(file1, "r");
    if (f1 == NULL) {
        warnvf("fcopy: file \"%s\" could not be opened for reading\n", file1);
        return;
    }

    f2=fopen(file2, write_mode);
    if (f2 == NULL) {
        warnvf("fcopy: file \"%s\" could not be opened for writing in mode %s\n", file2, write_mode);
        fclose(f1);
        return;
    }

    (void)fgets(line, 100, f1);
    while (1) {
        fputs(line, f2);
        if (fgets(line, 100, f1)==NULL) break;
    }

    fclose(f1); fclose(f2);
}


/**
 Ask the user if they want chromosome specific files.
 Return values:
  0 == One set of output files per chromosome.
  1 == Combine chromosomes into one output file.
 */
int global_ou_files(int default_opt)
{
    int reply=-1, global_opt=0;
    char stars[3];

//SL
    if (batchANALYSIS && Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].item_read) {
        return (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
    }

    if (DEFAULT_OPTIONS) {
        return (default_opt - 1);
    }

    ((default_opt == 2)?
     strcpy(stars, " *"):   strcpy(stars, "* "));

    while (reply != 0) {
        draw_line();
        printf("Output option menu for multiple chromosomes\n");
        printf("0) Done with menu - please proceed.\n");
        printf("%c 1) One set of output files per chromosome\n",
               stars[0]);
        printf("%c 2) Combine chromosomes into one output file\n",
               stars[1]);
        printf("Enter 0 - 2 > ");
        fcmap(stdin, "%d", &reply); printf("\n");
        if (reply != 1 && reply != 2 && reply != 0)
            printf("Unknown option %d.\n", reply);
        else if (reply == 1) {
            strcpy(stars, "* ");
            global_opt=0;
        } else if (reply == 2) {
            strcpy(stars, " *");
            global_opt=1;
        }
    }

    Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt = global_opt == 0 ? 'y': 'n';

    if (InputMode == INTERACTIVE_INPUTMODE)
        batchf(Loop_Over_Chromosomes);

    return global_opt;
}

void            enter_number(int *num)
{

    char            number[10];
    do
        {
            printf("\nEnter the number > ");
            fcmap(stdin, "%s", number);
            printf("\n");
        }
    while (*number < '0' || *number > '9');

    *num = atoi(number);

    return;
}

double randomnum( void )

{
    /*Global variable used:  seed1*/
    double r;

    /*From FORTRAN code:
      THIS FUNCTION GENERATES INDEPENDENT UNIFORM DEVIATES ON
      THE INTERVAL (0.0,1.0).  SEE THE REFERENCE:  WICHMAN B.A.
      AND HILL I.D.(1982). ALGORITHM 183: AN EFFICIENT AND PORTABLE
      PSEUDO-RANDOM NUMBER GENERATOR. APPLIED STATISTICS 31;188-190.
      SEED1, SEED2, AND SEED3 SHOULD BE SET TO INTEGER VALUES
      BETWEEN 1 AND 30000 BEFORE THE FIRST ENTRY.  INTEGER
      ARITHMETIC UP TO 30323 IS NEEDED ON YOUR COMPUTER.

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER SEED1,SEED2,SEED3
      SAVE SEED2,SEED3
      DATA SEED2,SEED3/2321,18777/
      Note that all three seeds - SEED1, SEED2, and SEED3 - must
      be initialized at the beginning of the main program.
    */

    seed1 = seed1 % 177 * 171 - seed1 / 177 * 2;
    seed2 = seed2 % 176 * 172 - seed2 / 176 * 35;
    seed3 = seed3 % 178 * 170 - seed3 / 178 * 63;

    if (seed1 < 0)
        seed1 += 30269;
    if (seed2 < 0)
        seed2 += 30307;
    if (seed3 < 0)
        seed3 += 30323;
    r = seed1 / 30269.00 + seed2 / 30307.00 + seed3 / 30323.00;
    return (r - (long)r);
    /*  RANDOM:=MOD(R,1.000)  Note: a mod b = a - ((a div b)*b  */
}

int randomindex(int N)

{
    int n=N+1;
    while(n >= N) {
        n=((int)((double)N * randomnum()));
    }
    return n;

}
#ifdef SETSEED
/*---------------------------------------------------+
  |  Seed the random number generator using a constant
  |  product of primes between 1 and 13 + 1
  +---------------------------------------------------*/
void seed_random(void) {
    seed1 = 30031;
}
#else
/*---------------------------------------------------+
  |  Seed the random number generator by using the    |
  |  current time.                                    |
  +---------------------------------------------------*/
void seed_random(void) {
#ifdef __INTEL__
    long secs;
#else
    time_t secs;
#endif


    time(&secs);
    /* printf("random seed: %d\n\n", (int) tv.tv_usec); */
    seed1 = (int) (secs % 30000);
}

#endif

void Exit(int arg, const char *file, const int line, const char *err)
{

#ifdef __INTEL__
    printf("Press any key to exit program.\n");
    fgetc(stdin);
    close_logs();
    exit(arg);
#else
    log_line(mssgf);
    goodbye(arg);
    printf("%s:%d Mega2 terminated. Error \"%s\" (#%d).\n", file, line, err, arg);
    fflush(stdout);
    exit(arg);
#endif

}


void print_only (const char *messg)

{
    printf("%s\n", messg);
    return;

}

time_t NOTIMEval;

void write_time(FILE *fil)

{
/*  extern   time_t  time();
    extern   struct  tm  *localtime(); */
    time_t   time_val;
    struct   tm     *time_stt;
/*    char     *month[12] = */
/*    { */
/*      "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec" */
/*    }; */

#ifdef HIDEDATE
    // There are two additional fields in 'struct tm' that are BSD/GNU extensions
    // tm_gmtoff, and tm_zone that are not used here but lacking these may generate a warning.
    struct tm tm = { 0, 14, 2, 29, 8-1, 1997-1900, 0, 0, 1};
    time_val = mktime(&tm);
    NOTIMEval = time_val;
#else
    time(&time_val);
#endif
    time_stt = localtime (&time_val);

    /*
      (void) fprintf (fil,"%2d-%s-%d  %02d:%02d:%02d\n",
      time_stt->tm_mday,month[time_stt->tm_mon],time_stt->tm_year,
      time_stt->tm_hour,time_stt->tm_min,time_stt->tm_sec);
    */
    fprintf(fil, "%s", asctime(time_stt));
}


int imax(int i1, int i2)
{
    if (i1 > i2) return i1;
    else return i2;
}

int imin(int i1, int i2)
{
    if (i1 < i2) return i1;
    else return i2;
}


double dmax(double d1, double d2)
{
    if (d1 > d2) return d1;
    else return d2;
}

int icomp(int i, int j)

{
    if (i > j) return 1;
    else if (i == j) return 0;
    else return -1;
}

void           exclaim(void)

{
    printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    return;
}



void chmod_X_file(const char *flname)
{
#if defined(MINGW) || defined(_WIN)
/*
    char *syscmd = CALLOC(strlen(flname) + 40, char);
    if(access(flname, F_OK)==0) {
        sprintf(syscmd, "chmod +x %s", flname);
        system(syscmd);
    }
    free(syscmd);
*/
#else
    struct stat stbuf;
    int err;
    err = stat(flname, &stbuf);
    if (err < 0) {
        switch(errno) {
        default:
            warnvf("chmod_X_file: stat(%s, buf) failed with errno %d (\"%s\")\n",
		   flname, errno, strerror(errno));
        }
        return;
    }
    if (S_ISDIR(stbuf.st_mode)) {
        warnvf("chmod_X_file: stat(%s, +x) file is a directory.  Will not change mode to +X\n",
	       flname);
        return;
    }
    stbuf.st_mode |= (S_IXUSR|S_IXGRP|S_IXOTH);
    err = chmod(flname, stbuf.st_mode);
    if (err < 0) {
        switch(errno) {
        default:
            warnvf("chmod_X_file: chmod(%s, +x) failed with errno %d (\"%s\")\n",
		   flname, errno, strerror(errno));
        }
        return;
    }
#endif
}

void delete_file(const char *flname)
{
    int err;
    if (access(flname, F_OK)==0) {
        err = unlink(flname);
        if (err < 0) {
            switch(errno) {
            case EACCES:
            case ENOENT:
            default:
                warnvf("delete_file: unlink(%s) failed with errno %d (\"%s\")\n",
                        flname, errno, strerror(errno));
            }
        }
    }
    return;
}

void copy_file(char *flname1, char *flname2)
{
    int fdi, fdo;
    long lenr, lenw;
    int err;
    if (access(flname1, F_OK)==0) {
        fdi = open(flname1, O_RDONLY | O_BINARY);
        if (fdi < 0) {
            switch(errno) {
            default:
                warnvf("copy_file: open(%s, \"r\") failed with errno %d (\"%s\")\n",
		       flname1, errno, strerror(errno));
            }
            return;
        }
        fdo = open(flname2, O_CREAT | O_TRUNC | O_WRONLY | O_BINARY, 0664);
        if (fdo < 0) {
            switch(errno) {
            case EACCES:
            case EEXIST:
            case EINVAL:
            case EMFILE:
            case ENOENT:
            default:
                warnvf("copy_file: open(%s, O_CREAT | O_TRUNC, 0664) failed with errno %d (\"%s\")\n",
		       flname2, errno, strerror(errno));
            }
            close(fdi);
            return;
        }
        while (1) {
            lenr =  read(fdi, multi_use_buffer, sizeof (multi_use_buffer));
            if (lenr < 0) {
                switch(errno) {
                default:
                    warnvf("copy_file: read(%s, , %d) = %d failed with errno %d (\"%s\")\n",
			   flname1, ((int)sizeof (multi_use_buffer)), lenr,
			   errno, strerror(errno));
                    break;
                }
            }
            lenw = write(fdo, multi_use_buffer, lenr);
            if (lenr != lenw) {
                switch(errno) {
                default:
                    warnvf("copy_file: write(%s, , %d) = %d failed with errno %d (\"%s\")\n",
			   flname2, lenr, lenw,
			   errno, strerror(errno));
                    break;
                }
            }
            /*
            warnvf("copy_file: read(%s) %d/write(%s) %d\n",
                    flname1, lenr, flname2, lenw);
            */
            if (lenr < MULTI_USE_BUFFER) break;
        }
        err = close(fdo);
        if (err < 0) {
            switch(errno) {
            default:
                warnvf("copy_file: close(%s) failed with errno %d (\"%s\")\n",
		       flname2, errno, strerror(errno));
            }
        }
        err = close(fdi);
        if (err < 0) {
            switch(errno) {
            default:
                warnvf("copy_file: close(%s) failed with errno %d (\"%s\")\n",
		       flname1, errno, strerror(errno));
            }
        }
    }
}

void init_file(char *file_name)
{
    int fd = open(file_name, O_CREAT | O_TRUNC, 0664);
    if (fd < 0) {
        switch(errno) {
        case EACCES:
        case EEXIST:
        case EINVAL:
        case EMFILE:
        case ENOENT:
        default:
            warnvf("init_file: open(%s) failed with errno %d (\"%s\")\n",
		   file_name, errno, strerror(errno));
        }
    }
    if (fd >= 0) {
        fd = close(fd);
        if (fd < 0) {
            switch(errno) {
            default:
                warnvf("init_file: close(%s) failed with errno %d (\"%s\")\n",
		       file_name, errno, strerror(errno));
            }
        }
    }
}

void move_file(char *flname1, char *flname2)
{
    int err;
    if (access(flname1, F_OK)==0) {
        if (access(flname2, F_OK) == 0) {
            err = unlink(flname2);
            if (err < 0) {
                switch(errno) {
                case EACCES:
                case ENOENT:
                default:
                    warnvf("move_file: unlink(%s) failed with errno %d (\"%s\")\n",
			   flname2, errno, strerror(errno));
                }
            }            
        }
        err = rename(flname1, flname2);
        if (err < 0) {
            switch(errno) {
            case EACCES:
            case ENOENT:
            case EINVAL:
            default:
                warnvf("move_file: rename(%s,%s) failed with errno %d (\"%s\")\n",
		       flname1, flname2, errno, strerror(errno));
            }
        }
    }
}

void backup_file(char *flname)
{
    char           *backupname;

    if (access(flname, F_OK) == 0) {
        if ((backupname = CALLOC((strlen(flname) + 5), char)) == NULL) {
            errorf("Unable to allocate space.");
            EXIT(MEMORY_ALLOC_ERROR);
        }
        sprintf(backupname, "%s.old", flname);
        move_file(flname, backupname);
        printf("Moved existing %s to %s\n", flname, backupname);
        free(backupname);
    }
}

void makedir(char *dirname)
{
    int err;
    if (access(dirname, F_OK)) {
#ifdef _WIN
        err = _mkdir(dirname);
#elif defined(MINGW)
        err = mkdir(dirname);
#else
        err = mkdir(dirname, 0775);
#endif
        if (err < 0) {
            switch(errno) {
            case EEXIST:
            case ENOENT:
            default:
                warnvf("makedir: mkdir(%s) failed with errno %d (\"%s\")\n",
		       dirname, errno, strerror(errno));
            }
        }
    }
}


int is_dir(char *dirname)
{
    struct stat stbuf;
    int err = stat(dirname, &stbuf);
    if (err < 0) {
        switch(errno) {
        default:
            warnvf("can not stat: stat(%s, buf) failed with errno %d (\"%s\")\n",
                   dirname, errno, strerror(errno));
        }
        return 0;
    }
    if (! (S_IFDIR & stbuf.st_mode)) {
        /*
           warnvf("access_dir: stat(%s, +x) path is NOT a directory.\n", dirname);
        */
        return 0;
    } else 
        return 1;
}

void getRunDate(void)
{
/* Should be defined in time.h
   extern time_t time();
   extern struct tm *localtime(); */
    time_t time_val;
    struct tm *time_struct;
    char *time_str = &(RunDate[0]);
    char day[3], hour[3], min[3];

    time(&time_val);
    time_struct = localtime(&time_val);
    if (time_struct->tm_mday < 10) {
        sprintf(day, "0%d", time_struct->tm_mday);
    } else {
        sprintf(day, "%2d", time_struct->tm_mday);
    }
    if (time_struct->tm_hour < 10) {
        sprintf(hour, "0%d", time_struct->tm_hour);
    } else {
        sprintf(hour, "%2d", time_struct->tm_hour);
    }
    if (time_struct->tm_min < 10) {
        sprintf(min, "0%d", time_struct->tm_min);
    } else {
        sprintf(min, "%2d", time_struct->tm_min);
    }

    sprintf(time_str, "%4d-%d-%s-%s-%s",
            time_struct->tm_year+1900, time_struct->tm_mon+1,
            day, hour, min);
    return;

}

void LogFileNames(void)
{
    sprintf(Mega2Log, "MEGA2.LOG");
    sprintf(Mega2LogRun, "%s/%s", sumdir, Mega2Log);

    sprintf(Mega2Err, "MEGA2.ERR");
    sprintf(Mega2ErrRun, "%s/%s", sumdir, Mega2Err);

    sprintf(Mega2Sim, "MEGA2.SIM");
    sprintf(Mega2SimRun, "%s/%s", sumdir, Mega2Sim);

    sprintf(Mega2Keys, "MEGA2.KEYS");
    sprintf(Mega2KeysRun, "%s/%s", sumdir, Mega2Keys);

    sprintf(Mega2Recode, "MEGA2.RECODE");
    sprintf(Mega2RecodeRun, "%s/%s", sumdir, Mega2Recode);

    sprintf(Mega2Reset, "MEGA2.RESET");
    sprintf(Mega2ResetRun, "%s/%s", sumdir, Mega2Reset);
    return;
}

void move_logs(char *new_location)
{
    char new_file[FILENAME_LENGTH+45];
/*      
    if (access(Mega2LogRun, F_OK) == 0){
        sprintf(new_file, "%s/%s", new_location, Mega2Log);
        move_file(Mega2LogRun,  new_file);
    }

    if (access(Mega2RecodeRun, F_OK) == 0){
        sprintf(new_file, "%s/%s", new_location, Mega2Recode);
        move_file(Mega2RecodeRun,  new_file);
    }

    if (access(Mega2ResetRun, F_OK) == 0){
        sprintf(new_file, "%s/%s", new_location, Mega2Reset);
        move_file(Mega2ResetRun,  new_file);
    }

    if (access(Mega2ErrRun, F_OK) == 0) {
        sprintf(new_file, "%s/%s",  new_location, Mega2Err);
        move_file(Mega2ErrRun,  new_file);
    }

    if (access(Mega2KeysRun, F_OK) == 0){
        sprintf(new_file, "%s/%s", new_location, Mega2Keys);
        move_file(Mega2KeysRun, new_file);
    }

    if (ErrorSimOpt == 1) {
        sprintf(new_file, "%s/%s", new_location, Mega2Sim);
        move_file(Mega2SimRun, new_file);
    }
*/
    if (strcmp(Mega2Batch, "none")) {
        char *p = strrchr(Mega2Batch, '/');
        if (p != NULL) {
            p++;
            sprintf(new_file, "%s/%s", new_location, p);
            copy_file(Mega2Batch, new_file);
        } else {
            if (strcmp(new_location, ".")) {
                sprintf(new_file, "%s/%s", new_location, Mega2Batch);
                copy_file(Mega2Batch, new_file);
            }
            /* else do not copy */
        }

    }
}

void get_line(FILE *fp, char *retline)
{
    int c;
    if ((c=fgetc(fp)) == COMMENT_CHAR) {
        skip_line(fp, c);
        strcpy(retline, "");
    } else {
        ungetc(c, fp);
        strcpy(retline, "");
        (void)fgets(retline, FILENAME_LENGTH-1, fp);
        if (!strcmp(retline, "\n")) {
            strcpy(retline, "");
        }
    }
}

/* void skip_line(FILE *fp)  */

/* { */

/*   int lch; */
/*   fcmap(fp, "%=\n", lch); */
/*   return; */
/* } */

/* Utility function to figure out number of markers by simply
   counting lines */

int linecountO(FILE *fp)
{
    int lch, cnt = 0;

    while(!feof(fp)) {
        lch = fgetc(fp);
        if (lch != '\n') {
            if (lch != '#') {
                cnt++;
            }
            fcmap(fp, "%=\n", lch);  /* Skip to end of line */
        }
    }
    rewind(fp);
/*    printf("linecountO: %d\n", cnt); */
    return cnt;

}

int linecount(FILE *fp)
{
    char buffer[4096], c = '\0';
    char *eof;
    int cnt = 0;
    size_t len;

    /* check pedigree file */
    while(!feof(fp))  {
        *buffer = 0;
        eof = fgets(buffer, sizeof(buffer), fp);
        if (eof == NULL) break;
	// The first character in the line...
        if (c == '\0') c = *buffer;
        len = strlen(buffer);
        if (len == sizeof(buffer) - 1 && buffer[len-1] != '\n') {
            continue;  /* really long line */
        }
	// skip blank lines and comment lines...
        if (c != '\n' && c != '#') cnt++;
	// we will be processing a new line...
        c = '\0';
    }
    rewind(fp);
/*    printf("linecount: %d\n", cnt); */
    return cnt;
}

int columncount(FILE *fp)
{
    int lch, cnt = 0;

    lch = fgetc(fp);
    while(!feof(fp)) {
        if (lch == '\n' && cnt) break;

        while (isspace(lch) && lch != '\n' && !feof(fp)) {
            lch = fgetc(fp);
        }
        if (feof(fp))
            break;
        else if (lch == '#' && !cnt) {
            fcmap(fp, "%=\n", lch);  /* Skip to end of line */
            lch = fgetc(fp);
            continue;
        } else if (lch == '\n') {
            if (cnt) break;
            continue;
        }

        cnt++;
        while (!isspace(lch) && lch != '\n' && !feof(fp)) {
            lch = fgetc(fp);
        }
    }
    rewind(fp);
    return cnt;
}

static int wget_failures(FILE *verf, int *ver, int *rev, int*patch, int *has_bugs)
{
#ifndef NET_UNIX
    if (verf == NULL) return CONNECTION_REFUSED;
    return SUCCESS;
#else
    char web_version[FILENAME_LENGTH];

    *rev=*ver=*patch=*has_bugs=0;
    if (verf == NULL) return CONNECTION_REFUSED;
    (void)fgets(web_version, FILENAME_LENGTH-1, verf);
    if (strstr(web_version, "Version") != NULL) {
        sscanf(web_version, "Version %d", ver);

        strcpy(web_version,"");
        (void)fgets(web_version, FILENAME_LENGTH - 1, verf);
        if (strcmp(web_version, "")) {
            sscanf(web_version, "Revision %d", rev);

            strcpy(web_version,"");
            (void)fgets(web_version, FILENAME_LENGTH - 1, verf);
            if (strcmp(web_version, "")) {
                sscanf(web_version, "Patch %d", patch);
            }
            strcpy(web_version,"");
            (void)fgets(web_version, FILENAME_LENGTH - 1, verf);
            if (strcmp(web_version, "")) {
                sscanf(web_version, "Bugs %d\n", has_bugs);
            } else {
                *has_bugs = 0;
            }
        }
    }
    return SUCCESS;
#endif
}

#ifndef NET_UNIX
static void show_bug(int cnt, char *line, int len, int type, void *ap)
{
    if (line == NULL) msgvf("error reading bug file; try again later\n");
    if (type > 0) msgvf("%s\n", line);
}

static void show_version(int cnt, char *line, int len, int type, void *ap)
{
    if (line == NULL) return;
    if (type > 0) msgvf("%s\n", line);
}
#endif

static void mega2_bug_report_check(const char *bug_report_file, int num_bugs)
{
#ifndef NET_UNIX
    printf("==============================================================\n");
    printf("  WARNING: %d %s been discovered in the released version:\n", num_bugs,
           ((num_bugs > 1)? "bugs have" : "bug has"));
    printf("\n");
    http_response("GET /pub/mega2/.bugs", "watson.hgen.pitt.edu", 80, show_bug, NULL);
    printf("\n");
    printf("==============================================================\n");

    return;
#else
    FILE *bugf;
    char *bug_report_line = multi_use_buffer;
    printf("==============================================================\n");
    printf("  WARNING: %d %s been discovered in the released version:\n", num_bugs,
           ((num_bugs > 1)? "bugs have" : "bug has"));
    printf("\n");
    bugf = http_response_body("GET /pub/mega2/.bugs", "watson.hgen.pitt.edu", 80);
    if (bugf == NULL) {
        printf("error reading bug file; try again later\n");
        return;
    }
    while(!feof(bugf)) {
        strcpy(bug_report_line, "");
        (void)fgets(bug_report_line, sizeof multi_use_buffer, bugf);
        chomp(bug_report_line);
        mssgf(bug_report_line);
    }
    socket_close_fd(bugf);
    printf("\n");
/*   printf("Please see the online Mega2 documentation for more information at\n"); */
/*   printf("    http://watson.hgen.pitt.edu/docs/mega2.html\n"); */
    printf("==============================================================\n");
#endif
}

static void mega2_version_report_check(const char *name)
{
    char *request = multi_use_buffer;
    sprintf(request, "GET /pub/mega2/%s", name);

#ifndef NET_UNIX
    http_response(request, "watson.hgen.pitt.edu", 80, show_version, NULL);
#else
    FILE *bugf;
    char *bug_report_line = multi_use_buffer;
    bugf = http_response_body(name, "watson.hgen.pitt.edu", 80);
    if (bugf == NULL) {
        return;
    }
    while(!feof(bugf)) {
        strcpy(bug_report_line, "");
        (void)fgets(bug_report_line, sizeof multi_use_buffer, bugf);
        chomp(bug_report_line);
        mssgf(bug_report_line);
    }
    socket_close_fd(bugf);
#endif

}

/* function to check the web-version using wget */

struct http_version {
    int ver;
    int rev;
    int patch;
    int has_bugs;
} hv = { -1, -1, -1, 0};

#ifndef NET_UNIX
static void get_mega2_version(int cnt, char *line, int len, int type, void *ap)
{
    struct http_version *hvp = (struct http_version *)ap;

    if (type == 0) return;
    if (strstr(line, "Version") != NULL) {
        sscanf(line, "Version %d", &hvp->ver);

    } else if (strstr(line, "Revision") != NULL) {
        sscanf(line, "Revision %d", &hvp->rev);

    } else if (strstr(line, "Patch") != NULL) {
        sscanf(line, "Patch %d", &hvp->patch);

    } else if (strstr(line, "Bugs") != NULL) {
        sscanf(line, "Bugs %d", &hvp->has_bugs);
    }
    if (debug) msgvf("#%d: %s (L%d) %d %d.%d.%dx%d\n", cnt, line, len, type, hvp->ver, hvp->rev, hvp->patch, hvp->has_bugs);

}
#endif

void mega2_version_check(void)
{
    FILE *verf;
    int newver=0, newrev=0, newpatch=0, verrev, newverrev, has_bugs;
    char ans;
    char version_name[20];

    printf("Checking the latest Mega2 release version (-w or --noweb to skip) \n");

#ifndef NET_UNIX
    http_response("GET /pub/mega2/.version", "watson.hgen.pitt.edu", 80, get_mega2_version, &hv);
    if (hv.ver == -1) {
        verf = NULL;
    } else {
        verf     = stdout; /* some random pointer */
        newver   = hv.ver;
        newrev   = hv.rev;
        newpatch = hv.patch;
        has_bugs = hv.has_bugs;
    }
#else
    verf = http_response_body("GET /pub/mega2/.version", "watson.hgen.pitt.edu", 80);
#endif

    switch(wget_failures(verf, &newver, &newrev, &newpatch, &has_bugs)) {
    case COMMAND_NOT_FOUND:
        printf("wget not installed on your system.\n");
        break;

    case CONNECTION_REFUSED:
    case SERVER_NOT_FOUND:
        printf("watson may be down, try again later for latest version check.\n");
        break;

    case SUCCESS:
        newverrev = 100*newver + 10*newrev;
        verrev = 100*Mega2Ver + 10*Mega2Rev;

#ifdef HIDESTATUS
        printf("OLD %d, NEW %d\n", verrev, newverrev);
#endif

        if (newverrev == 0) {
            printf("Unable to check web version.\n");
        } else if (newverrev == verrev) {

            /* (ver == Mega2Ver && rev == Mega2Rev && patch == Mega2Patch) { */
            if (newpatch == Mega2Patch) {
                printf("Mega2 version is up to date.\n");
            } else if (newpatch < Mega2Patch) {
                printf("Your Mega2 version is newer than the official distribution.\n");
                printf("    Latest version: %d.%d.%d\n", newver, newrev, newpatch);
                printf("    Your   version: %d.%d.%d\n", Mega2Ver, Mega2Rev, Mega2Patch);
            } else {
                printf("There is a new minor revision of Mega2; You might want to get it.\n");
                printf("    Latest version: %d.%d.%d\n", newver, newrev, newpatch);
                printf("    Your   version: %d.%d.%d\n", Mega2Ver, Mega2Rev, Mega2Patch);

                sprintf(version_name, ".bug.%d.%d.%d", Mega2Ver, Mega2Rev, Mega2Patch);
                mega2_version_report_check(version_name); // ... new name
                sprintf(version_name, ".release.%d.%d.%d", Mega2Ver, Mega2Rev, Mega2Patch);
                mega2_version_report_check(version_name);
                draw_line();
            }
        } else if (newverrev > verrev) {
/*       (ver > Mega2Ver || rev > Mega2Rev || patch > Mega2Patch) { */
            printf("Upgrade needed: \n");
            printf("        Release version: %d.%d.%d\n", newver, newrev, newpatch);
            printf("        Your    version: %d.%d.%d\n", Mega2Ver, Mega2Rev, Mega2Patch);
            printf("Exit Mega2 now, while you update [y/n] ? ");

            (void)scanf("%c", &ans); /* newline; */
            if (ans == 'y' || ans == 'Y') {
                exit(0); // without writing to the log...
            } else {
                printf("Please update Mega2 as soon as possible.\n");
            }
        } else if (newverrev < verrev) {
            /* (ver < Mega2Ver || rev < Mega2Rev || patch < Mega2Patch) { */
            printf("You appear to be running a more recent version of Mega2.\n");
            printf("        Release version: %d.%d.%d\n", newver, newrev, newpatch);
            printf("        Your    version: %d.%d.%d\n", Mega2Ver, Mega2Rev, Mega2Patch);
        }

        if (has_bugs) {
            mega2_bug_report_check(".bugs", has_bugs);
        }

        break;

    default:
        printf("Could not check web-version, please try later.\n");
        break;
    }

#ifdef NET_UNIX
    socket_close_fd(verf);
#endif

    return;
}

#ifdef USE_GETOPT
void mega2_opts(int argc, char **argv)
{
    int c;

    static struct option long_options[] =	{
        /* These options don't set a flag.
           We distinguish them by their indices. */
        {"noweb", no_argument, 0, 'w'},
        {"nosave", no_argument, 0, 'x'},
        {"help",     no_argument,       0, 'h'},
        {"version",  no_argument,       0, 'v'},
        {0, 0, 0, 0}
    };

    /* getopt_long stores the option index here. */
    int option_index;

    for (c=1; c < argc; c++) {
        if (!strcmp(argv[c], "-nosave")) {
            printf("Mega2 no longer supports the -nosave option.\n");
            print_mega2_help();
            exit(-1);
        }
    }
    while (1) {
        option_index=0;
        c = getopt_long (argc, argv, "wxvh",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c) {
        case 'w':
            check_web_ver = 0;
            break;

        case 'x':
            CreateRunFolder=0;
            break;

        case 'h':
            print_mega2_help();
            break;

        case 'v':
            print_mega2_version();
            break;

        case '?':
            /* getopt_long already printed an error message. */
            print_mega2_help();
            exit(-1);
            break;

        default:
            exit(-1);
        }
    }

/*   if (!check_web_ver) */
/*     puts ("check_web_ver flag is off"); */

/*   if (!CreateRunFolder) { */
/*     puts("Saving logs is off"); */
/*   } */

    /* There should be exactly one command line argument that is not an option. */
    if (optind < argc)  {
        if ((argc - optind) >= 2) {
            printf("Invalid arguments to Mega2.\n");
            EXIT(INPUT_DATA_ERROR);
        } else {
            strcpy(Mega2Batch, argv[optind]);
        }
    }
    return;
}

#else /*  USE_GETOPT */

extern void set_chromosomes(const int autosome, const int xy, const int mito);
extern void set_chromosomes_horse(); /* { set_chromosomes(31,  0,  0); } */
extern void set_chromosomes_sheep(); /* { set_chromosomes(26,  0,  0); } */
extern void set_chromosomes_dog();   /* { set_chromosomes(38, 41,  0); } */
extern void set_chromosomes_mouse(); /* { set_chromosomes(19,  0,  0); } */
extern void set_chromosomes_cow();   /* { set_chromosomes(29,  0,  0); } */
extern void set_chromosomes_human(); /* { set_chromosomes(22, 25, 26); } */
extern int lastautosome, pseudoautosome, mitoautosome;
extern int missingv_flags;
extern char *quant_in, *quant_out, *affect_in, *affect_out;

void mega2_opts(int argc, char **argv)
{
    int   c;
    char *as;

    argv++;
    while (--argc > 0) {
        as = *argv;
        if ((*as) == '-') {
            if (strcasecmp(as, "-nosave") == 0) {
                printf("Mega2 no longer supports the -nosave option.\n");
                print_mega2_help();
                EXIT(INPUT_DATA_ERROR);
            } else if (*(as+1) == '-') {
                as += 2;
                if ( (strcasecmp(as, "autosome") == 0) || (strcasecmp(as, "pseudo") == 0) ||
                     (strcasecmp(as, "mito") == 0)) {
                    int tmp;
                    if (--argc > 0) {
                        argv++;
                        tmp  = atoi(*argv);
                        if (tmp <= 0) {
                            print_mega2_help();
                            EXIT(INPUT_DATA_ERROR);
                        } else {
                            if (strcasecmp(as, "autosome") == 0)
                                lastautosome = tmp;
                            else if (strcasecmp(as, "pseudo") == 0)
                                pseudoautosome = tmp;
                            else if (strcasecmp(as, "mito") == 0)
                                mitoautosome = tmp;
                        }
                    } else {
                        print_mega2_help();
                        EXIT(INPUT_DATA_ERROR);
                    }
                } else if (strcasecmp(as, "horse") == 0)
                    set_chromosomes_horse();
                 else if (strcasecmp(as, "sheep") == 0)
                     set_chromosomes_sheep();
                 else if (strcasecmp(as, "dog") == 0)
                     set_chromosomes_dog();
                 else if (strcasecmp(as, "mouse") == 0)
                     set_chromosomes_mouse();
                 else if (strcasecmp(as, "cow") == 0)
                     set_chromosomes_cow();
                 else if (strcasecmp(as, "human") == 0)
                     set_chromosomes_human();
                else if (strcasecmp(as, "noweb") == 0)
                        check_web_ver = 0;
                else if (strcasecmp(as, "nosave") == 0)
                    CreateRunFolder=0;
                else if (strcasecmp(as, "bed") == 0)
                    Input_Format = in_format_binary_PED;
                else if (strcasecmp(as, "ped") == 0)
                    Input_Format = in_format_PED;
                else if (strcasecmp(as, "bcf") == 0)
                    Input_Format = in_format_binary_VCF;
                else if (strcasecmp(as, "vcf.gz") == 0)
                    Input_Format = in_format_compressed_VCF;
                else if (strcasecmp(as, "vcf") == 0)
                    Input_Format = in_format_VCF;
                else if (strcasecmp(as, "mega2") == 0)
                    Input_Format = in_format_mega2;
                else if (strcasecmp(as, "linkage") == 0)
                    Input_Format = in_format_linkage;
                else if (strcasecmp(as, "extended_linkage") == 0)
                    Input_Format = in_format_extended_linkage;
                else if (strcasecmp(as, "input") == 0)
                    Input_Format = in_format_traditional;
                else if (strcasecmp(as, "mega2") == 0)
                    Input_Format = in_format_mega2;
                else if (strcasecmp(as, "linkage") == 0)
                    Input_Format = in_format_linkage;
                else if (strcasecmp(as, "hybrid") == 0)
                    Input_Format = in_format_extended_linkage;
                else if (strcasecmp(as, "force_numeric_alleles") == 0)
                    force_numeric_alleles = 1;
                else if (strcasecmp(as, "quant_in") == 0) {
		    argv++; --argc;
		    quant_in = *argv;
		    missingv_flags |= 1;
		} else if (strcasecmp(as, "quant_out") == 0) {
		    argv++; --argc;
		    quant_out = *argv;
		    missingv_flags |= 2;
                } else if (strcasecmp(as, "affect_in") == 0) {
		    argv++; --argc;
		    affect_in = *argv;
		    missingv_flags |= 4;
                } else if (strcasecmp(as, "affect_out") == 0) {
		    argv++; --argc;
		    affect_out = *argv;
		    missingv_flags |= 8;
		} else if (strcasecmp(as, "help") == 0)
                    print_mega2_help();
                else if (strcasecmp(as, "version") == 0)
                    print_mega2_version();
                else {
                    print_mega2_help();
                    EXIT(INPUT_DATA_ERROR);
                }
            } else {
                extern int MARKER_SCHEME;
                while ((c = *++as)) {
                    switch (c) {
                    case 'w': case 'W':
                        check_web_ver = 0;
                        break;
                    case 'x': case 'X':
                        CreateRunFolder=0;
                        break;
                    case 'h': case 'H':
                        print_mega2_help();
                        break;
                    case 'd': case 'D':
                        debug++;
                        break;
#ifdef TEST
                    case 'f': case 'F':
                        test_socket_fd();
                        break;
                    case 's': case 'S':
                        test_socket_s();
                        break;
#endif
                    case '1':
                        MARKER_SCHEME = MARKER_SCHEME_BITS;
                        break;
                    case '2':
                        MARKER_SCHEME = MARKER_SCHEME_BYTE;
                        break;
                    case '3':
                        MARKER_SCHEME = MARKER_SCHEME_PTR;
                        break;
                    case 'v': case 'V':
                        print_mega2_version();
                        break;
                    default:
                        print_mega2_help();
                        EXIT(INPUT_DATA_ERROR);
                        break;
                    }
                }
            }
        } else {
            break;
        }
        argv++;
    }

    if (lastautosome) {
        subject_organism = UNKNOWN_SO;
        set_chromosomes(lastautosome, pseudoautosome, mitoautosome);
    } else {
       set_chromosomes_human(); /* { set_chromosomes(22, 25, 26); } */
    }

    /* There should be exactly one command line argument that is not an option. */
    if (argc == 0) ;
    else if (argc == 1)  {
        strcpy(Mega2Batch, *argv);
    } else {
        printf("Invalid arguments to Mega2.\n");
        print_mega2_help();
        EXIT(INPUT_DATA_ERROR);
    }
    return;
}
#endif /*  USE_GETOPT */

void print_mega2_help(void)

{

    printf("Usage: mega2 [options] [batch-file-name]\n");
    printf("  acceptable options:\n");
    printf("             --cow\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for cow.\n");
    printf("             --dog\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for dog.\n");
    printf("             --horse\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for horse.\n");
    printf("             --human\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for human.\n");
    printf("             --mouse\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for mouse.\n");
    printf("             --sheep\n");
    printf("                set last autosome, pseudo autosome, and mitocondria chromosome numbers for sheep.\n");
    printf("             --autosome <value>\n");
    printf("                set last autosome chromosome number to <value>.\n");
    printf("             --pseudo <value>\n");
    printf("                set pseudo autosome xy chromosome to <value>.\n");
    printf("             --mito <value>\n");
    printf("                set mitocondria chromosome to <value>.\n");
    printf("             --bed\n");
    printf("                input files are in PLINK binary Ped format (bed).\n");
    printf("             --ped\n");
    printf("                input files are in PLINK Ped format (ped).\n");
    printf("             --bcf\n");
    printf("                input files are in binary Variant Call File (VCF) format (bcf).\n");
    printf("             --vcf.gz\n");
    printf("                input files are in compressed Variant Call File (VCF) format (vcf.gz).\n");
    printf("             --vcf\n");
    printf("                input files are in Variant Call File (VCF) format (vcf).\n");
    printf("             --mega2\n");
    printf("                input files are in Mega2 format (tabular files with header line).\n");
    printf("             --linkage\n");
    printf("                input files are in Linkage format.\n");
    printf("             --extended_linkage\n");
    printf("                input files are in Linkage format with a Mega2 names file vs a linkage locus file.\n");

    printf("             --force_numeric_alleles\n");
    printf("                recode alleles as numbers even though analysis can accept letter alleles.\n");

    printf("             -x, --nosave\n");
    printf("                    Do not create a new run-folder.\n");
    printf("             -w, --noweb\n");
    printf("                    Do not check for latest on-line version.\n");
    printf("             -h, --help\n");
    printf("                    Print this message.\n");
    printf("             -v, --version\n");
    printf("                    Print mega2 version number.\n");
    exit(0);

}

void print_mega2_version(void)
{
    printf("Mega2 version %s\n", Mega2Version);
    exit(0);
}

int show_system_cmd(const char *str, const char* file, const int line)
{
#ifdef HIDESTATUS
    printf("%s:%d %s\n", file, line, str); fflush(stdout);
#endif
#ifndef VALGRIND
    return system(str);
#else
    return 0;
#endif
}

#ifdef MINGW
static char perl_pgm_buf[FILENAME_LENGTH+6];
const char *perl_pgm(const char *pl)
{
    FILE *wh;
    char *cp;

    sprintf(perl_pgm_buf, "where %s", pl);
    wh = popen(perl_pgm_buf, "r");
    strcpy(perl_pgm_buf, "perl ");
    (void)fgets(perl_pgm_buf+5, FILENAME_LENGTH, wh);
    fclose(wh);
    perl_pgm_buf[strlen(perl_pgm_buf)-1] = 0; /* nl */
    for (cp = perl_pgm_buf; *cp; cp++) if (*cp == '\\') *cp = '/';
    return perl_pgm_buf;
}
#else
const char *perl_pgm(const char *pl)
{
    return pl;
}
#endif

/*  item_type get_item(FILE *fp, item_value_type item_type, value_struct *value) */

/*  { */

/*    int c; */
/*    char str[FILENAME_LENGTH]; */
/*    char endptr[1]; */
/*    if ((c=fgetc(fp)) == COMMENT_CHAR) { */
/*      skip_line(fp); */
/*      return(COMMENT); */
/*    } */
/*    fscanf(fp, "%s", str); */

/*    switch(item_type) { */
/*    case INT: */
/*      value->Int = strtol(str, &endptr, 10); */
/*      if (*endptr != "NULL") { */
/*        input_file_error(); */
/*      } */
/*      break; */
/*    case CHAR: */
/*      if (strlen(str) != 1) { */
/*        input_file_error(); */
/*      } */
/*      else { */
/*        value->Char = str[0]; */
/*      } */
/*      break; */
/*    case STRING: */
/*      break; */
/*    case FLOAT: */
/*      value->Float = strtod(str, &endptr); */
/*      if (*endptr != "NULL") { */
/*        input_file_error(); */
/*      } */
/*      break; */

/*    default: */
/*      input_file_error(); */
/*      break; */
/*    } */

/*    return(value); */

/*  } */

/*

   From: http://www.cora.nwra.com/~gourlay/software/VortexSystem/my_malloc.c

   GPL'd (Version 2)

   Usage:
   #define CALLOC(nelem,type) my_calloc(NULL,nelem,sizeof(type),__FILE__,__LINE__)

   Old definition:
   #define CALLOC(n, t) ((n > 0)? (t *) calloc((n), sizeof(t)) : NULL)


// ARGUMENTS
//   ptr (in/out): address of the user portion of previously allocated memory
//
//   nelem (in): number of elements to allocate
//
//   elsize (in): size of the element, in bytes
//
//   file (in): name of the file from which mjg_realloc was called
//
//   line (in): line number of the file from which mjg_realloc was called
//
//
*/
void * my_calloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line)
{
    void *mem = ptr;

    /* Check the validity of the input parameters */
    /* Note: nelem and elsize are size_t, which can not be less than zero */
    if (nelem * elsize == 0) {
        sprintf(err_msg, "WARNING: calloc allocating no memory at line %i of file %s\n",line,file);
        mssgf(err_msg);
        sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                "WARNING: my_calloc called with nelem=%Iu elsize=%Iu\n",
#else
                "WARNING: my_calloc called with nelem=%zu elsize=%zu\n",
#endif
                nelem, elsize);
        mssgf(err_msg);
        return NULL;
    }

    mem= (void *) calloc(nelem, elsize);

    if (mem==NULL)
        {
            sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                    "ERROR: ran out of memory: allocation of %Iu bytes failed on line %i of file %s\n",
#else
                    "ERROR: ran out of memory: allocation of %zu bytes failed on line %i of file %s\n",
#endif
                    nelem * elsize, line, file);
            log_line(mssgf);
            mssgf(err_msg);
            EXIT(MEMORY_ALLOC_ERROR);
        }
    return mem;
}

/* Usage:
   #define MALLOC(type)    ((type *) my_malloc(NULL,1,sizeof(type),__FILE__,__LINE__))
   Standard definition:     void * malloc(size_t size);
   Old define:  #define MALLOC(t) ((t *) malloc(sizeof(t)))
*/
void * my_malloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line)
{
    void *mem = ptr;

    /* Check the validity of the input parameters */
    /* Note: nelem and elsize are size_t, which can not be less than zero */
    if (elsize == 0) {
        sprintf(err_msg, "WARNING: malloc allocating no memory at line %i of file %s\n",line,file);
        mssgf(err_msg);
        sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                "WARNING: my_malloc called with elsize=%Iu\n",
#else
                "WARNING: my_malloc called with elsize=%zu\n",
#endif
                elsize);
        mssgf(err_msg);
        return NULL;
    }

    mem= (void *) malloc(elsize);

    if (mem==NULL)
        {
            sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                    "ERROR: ran out of memory: allocation of %Iu bytes failed on line %i of file %s\n",
#else
                    "ERROR: ran out of memory: allocation of %zu bytes failed on line %i of file %s\n",
#endif
                    elsize, line, file);
            log_line(mssgf);
            mssgf(err_msg);
            EXIT(MEMORY_ALLOC_ERROR);
        }
    return mem;
}

/*
  Usage:
  #define REALLOC(ptr,n,t)     my_realloc(ptr,n,sizeof(t),__FILE__,__LINE__)

  Old defintion:
  #define REALLOC(p, n, t) ((t *) realloc((p), (n) * sizeof(t)))

*/
void * my_realloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line)
{
    void *mem = ptr;

    /* Check the validity of the input parameters */
    /* Note: nelem and elsize are size_t, which can not be less than zero */
    if (nelem*elsize == 0) {
        sprintf(err_msg, "WARNING: realloc allocating no memory at line %i of file %s\n",line,file);
        mssgf(err_msg);
        sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                "WARNING: my_realloc called with nelem=%Iu elsize=%Iu\n",
#else
                "WARNING: my_realloc called with nelem=%zu elsize=%zu\n",
#endif
                nelem, elsize);
        mssgf(err_msg);
        return NULL;
    }

    mem= (void *) realloc(ptr,nelem*elsize);

    if (mem==NULL)
        {
            sprintf(err_msg,
#if defined(_WIN) || defined(MINGW)
                    "ERROR: ran out of memory: reallocation of %Iu bytes failed on line %i of file %s\n",
#else
                    "ERROR: ran out of memory: reallocation of %zu bytes failed on line %i of file %s\n",
#endif
                    nelem*elsize, line, file);
            log_line(mssgf);
            mssgf(err_msg);
            EXIT(MEMORY_REALLOC_ERROR);
        }
    return mem;
}

// Return a pointer to the tail of the string which consists of at least 'n' characters
// If the string is shorter than 'n' characters, return a pointer to the entire string
char *strtail(const char *name, const size_t n)
{
    if (strlen(name) > n) {
        // The last n characters of the name...
      return (char *)&(name[strlen(name) - n]);
    }
    return (char *)name;
}

void executable_in_path_does_not_exist_csh(FILE *fp, const char *prog, const char *msg)
{
  // 'which' searches your path for for an executable....
  //fprintf(fp, "if ( \"`/usr/bin/which %s`\" == \"\") then\n", prog);
  fprintf(fp, "if ( \"`which %s`\" == \"\") then\n", prog);
  fprintf(fp, "  echo \"%s\"\n", msg); // doesn't look for " in message
  fprintf(fp, "endif\n");
}

//
// Used to output a check for assignment to an environment variable and setting
// it if it has not been assigned. For example...
// sprintf_checkset_csh(sh_fp, "_MERLIN", "merlin");
// will produce the following output:
// if ( ! $?_MERLIN) then
//  set _MERLIN='merlin'
// endif
void fprintf_env_checkset_csh (FILE *fp, const char *var, const char *default_value)
{
  fprintf(fp, "if ( ! $?%s) then\n", var);
  fprintf(fp, "  set %s='%s'\n", var, default_value);
  fprintf(fp, "endif\n");
}

//
// Used to check the status of a progrem's execution. It will copy the exit status
// into an environment variable and a file (with the same name), and then test for
// a failure exit status (e.g. != 0) and if found will issue an error message and exit.
//
// NOTE: The '$status' (csh) / '$?' (bash) is only valid immediately after the program has
// finished execution so this routine must be called then for it to be able to effctively
// insert it's text into the output stream.
//
// fprintf_status_check_csh(filep, "merlin");
// should produce the following code...
//
// set merlin_status=$status
// echo $merlin_status > merlin_status
// if ($merlin_status != 0) then
//   echo "Run of merlin failed with status code: $merlin_status"
//   exit $merlin_status
// endif
void fprintf_status_check_csh (FILE *fp, const char *prog, const int exit_on_error)
{
  //fprintf(fp, "set %s_status=$?\n", prog);
  fprintf(fp, "set %s_status=$status\n", prog);
  fprintf(fp, "echo $%s_status > %s_status\n", prog, prog);
  fprintf(fp, "if ($%s_status != 0) then\n", prog);
  fprintf(fp, "  echo \"ERROR: Run of '%s' failed with status code: $%s_status\"\n", prog, prog);
  if (exit_on_error != 0) fprintf(fp, "  exit $%s_status\n", prog);
  fprintf(fp, "endif\n");
}


//
// A boolean function that determines if at least one trait is available.
// Generally used for those analysis that require a trait.
int have_trait_b (linkage_locus_top *LTop, bool affect, bool quant) {
    int tr, trait, *trp = global_trait_entries;
    
    // Look for at least one trait because Eigenstrat requires one...
    for (tr=0; tr < num_traits; tr++) {
        linkage_locus_rec *loc;
        trait = *trp++;
        if (trait == -1) continue;
        if (trait == -99) break;
        loc = &(LTop->Locus[trait]);
//      if (loc->Type != AFFECTION && loc->Type != QUANT) continue;
//      return 1;

        if (affect && loc->Type == AFFECTION) return 1;
        if (quant  && loc->Type == QUANT) return 1;
    }
    return 0;
}
