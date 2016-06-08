/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"

#include "utils_ext.h"
/*
              utils_ext.h:  EXIT hello summary_time_stamp
*/


void           open_logs(void);
void           time_stamp_logs(void);
void           close_logs(void);
void           err_or_warn(FILE **errfp, FILE **logfp);
void           errorf(const char *messg);
void           warnf(const char *messg);
void           mssgf(const char *messg);
void           errorvf(const char *fmt, ...);
void           warnvf(const char *fmt, ...);
void           mssgvf(const char *fmt, ...);
void           msgvf(const char *fmt, ...);
void           errsimf(const char *messg);
/* void           input_file_error(); */
void           empty_file(const char *filename, error_type etype);
/* void           memory_alloc_error(); */

static FILE *Mega2logf = NULL;
static FILE *Mega2errf = NULL;

void open_logs(void)
{

    Mega2logf = fopen(Mega2LogRun, "w");
    if (Mega2logf == NULL) {
        fprintf(stderr,
                "Could not open Mega2 log file (not enough access rights?)!\n");
        EXIT(FILE_WRITE_ERROR);

    }
    Mega2errf = fopen(Mega2ErrRun, "w");
    if (Mega2errf == NULL) {
        fprintf(stderr,
                "Could not open Mega2 error file (not enough access rights?)!\n");
        EXIT(FILE_WRITE_ERROR);
    }

    hello(Mega2logf);
    hello(Mega2errf);
}

void time_stamp_logs(void)
{
    summary_time_stamp(mega2_input_files, Mega2logf, "");
    summary_time_stamp(mega2_input_files, Mega2errf, "");
}

void close_logs(void)
{
    fclose(Mega2logf);
    fclose(Mega2errf);

    Mega2logf = NULL;
    Mega2errf = NULL;
}

void err_or_warn(FILE **errfp, FILE **logfp)
{
    *logfp = Mega2logf;
    *errfp = Mega2errf;
}

/* NM- could not implement as variable arguments, because
   message is a char pointer, implemented as a
   global variable */

void errf(const char *messg, const char *shout, int showerrf)
{
    if (showerrf && Mega2errf != NULL) {
        fprintf(Mega2errf, "%s%s\n", shout, messg);
        fflush(Mega2errf);
    }

    if (Display_Errors==1 || Mega2logf == NULL || Mega2errf == NULL) {
        if (Mega2logf != NULL) {
            fprintf(Mega2logf, "%s%s\n", shout, messg);
            fflush(Mega2logf);
        }
        fprintf(stderr, "%s%s\n", shout, messg);
        fflush(stderr);
    }
}

void warnf(const char *messg)
{
    errf(messg, "WARNING: ", 1);
}

/* write in ERROR file as well as LOG file */
void errorf(const char *messg)
{
    errf(messg, "ERROR: ", 1);
}

/* logf is called only  */
void mssgf(const char *messg, int show)
{
    if (show && Mega2logf != NULL) {
        fprintf(Mega2logf, "%s\n", messg);
        fflush(Mega2logf);
    }

    if (Display_Messages == 1 || Mega2logf == NULL) {
        fprintf(stdout, "%s\n", messg);
        fflush(stdout);
    }
}
void mssgf(const char *messg)
{
    mssgf(messg, 1);
}

void errvf(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    if (Mega2errf != NULL) {
        vfprintf(Mega2errf, fmt, ap);
        fflush(Mega2errf);
    }
    va_end(ap);

    va_start(ap, fmt);
    if (Display_Errors==1 || Mega2logf == NULL || Mega2errf == NULL) {
        if (Mega2logf != NULL) {
            vfprintf(Mega2logf, fmt, ap);
            fflush(Mega2logf);
            va_end(ap);
            va_start(ap, fmt);
        }
        vfprintf(stderr, fmt, ap);
        fflush(stderr);
    }
    va_end(ap);
}

void warnvf(const char *fmt, ...)
{
    va_list ap;
    const char *shout = "WARNING: ";

    va_start(ap, fmt);
    if (Mega2errf != NULL) {
        fputs(shout, Mega2errf);
        vfprintf(Mega2errf, fmt, ap);
        fflush(Mega2errf);
    }
    va_end(ap);

    va_start(ap, fmt);
    if (Display_Errors==1 || Mega2logf == NULL || Mega2errf == NULL) {
        if (Mega2logf != NULL) {
            fputs(shout, Mega2logf);
            vfprintf(Mega2logf, fmt, ap);
            fflush(Mega2logf);
            va_end(ap);
            va_start(ap, fmt);
        }
        fputs(shout, stderr);
        vfprintf(stderr, fmt, ap);
        fflush(stderr);
    }
    va_end(ap);
}

void errorvf(const char *fmt, ...)
{
    va_list ap;
    const char *shout = "ERROR: ";

    va_start(ap, fmt);
    if (Mega2errf != NULL) {
        fputs(shout, Mega2errf);
        vfprintf(Mega2errf, fmt, ap);
        fflush(Mega2errf);
    }
    va_end(ap);

    va_start(ap, fmt);
    if (Display_Errors==1 || Mega2logf == NULL || Mega2errf == NULL) {
        if (Mega2logf != NULL) {
            fputs(shout, Mega2logf);
            vfprintf(Mega2logf, fmt, ap);
            fflush(Mega2logf);
            va_end(ap);
            va_start(ap, fmt);
        }
        fputs(shout, stderr);
        vfprintf(stderr, fmt, ap);
        fflush(stderr);
    }
    va_end(ap);
}

/* logf is called only  */
void mssgvf(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    if (Mega2logf != NULL) {
        vfprintf(Mega2logf, fmt, ap);
        fflush(Mega2logf);
    }
    va_end(ap);

    va_start(ap, fmt);
    if (Display_Messages == 1 || Mega2logf == NULL) {
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
    }
    va_end(ap);
}

void msgvf(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    if (Mega2logf != NULL) {
        vfprintf(Mega2logf, fmt, ap);
        fflush(Mega2logf);
    }
    va_end(ap);

    va_start(ap, fmt);
    if (Display_Messages == 1 || Mega2logf == NULL) {
        vfprintf(stdout, fmt, ap);
        fflush(stdout);
    }
    va_end(ap);
}

void errsimf(const char *messg)

{
    FILE *logf;
    int  first_time;

    first_time  = ((access(Mega2SimRun, F_OK) != 0)? 1  : 0);

    logf = fopen(Mega2SimRun, "a");
    if (logf == NULL) {
        fprintf(stderr,
                "Could not open Mega2 log file (not enough access rights?)!\n");
        EXIT(FILE_WRITE_ERROR);
    }

    if (first_time) {
        /* Mega2ErrSim does not exist, first call to logf */
        hello(logf);
        /* By this time, we should have read all the input file names,
           this is only reached when executing the simulate_errors
           routine.
        */
        summary_time_stamp(mega2_input_files, logf, "");
    }

    fprintf(logf, "%s\n", messg);

    fprintf(stdout, "%s\n", messg);
    fclose(logf);
}

void empty_file(const char *filename, error_type etype)
{
    sprintf(err_msg, "Input file %s is empty.", filename);
    errorf(err_msg);
    if (etype == FATAL) {
        EXIT(FILE_READ_ERROR);
    }
}
