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
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  errorvf
     utils_ext.h:           EXIT
 */

#ifndef ANSI
#define ANSI
#endif

int NARGS, STAT;

static void skip_past_char(FILE *filep, const char c);
static int get_sign(FILE *filep, int mode);
static int read_int(FILE *filep);
static double read_double(FILE *filep);
static void read_string(FILE *filep, char *sval, char tch, int mode);


/* make these global so they needn't be passed around */
static int LCH, TERM;

/* skip_past_char()
 *
 * Scan the given stream until the given character
 * (or EOF or TERM) is encountered. Leave LCH undefined
 * if the character was found.
 */

static void skip_past_char(FILE *filep, const char c)

{
    if (!LCHDEF) GETLCH;
    if (LCH_IS_TERM) { NARGS--; DEFTERM; return; }
    while ((LCH != EOF) && (LCH != c) && (LCH_NOT_TERM)) {
        GETLCH;
    }
    if (LCH == c) UNDEFLCH;
}


/* get_sign()
 *
 * Get a sign (from '+''s and '-''s) and return it.
 *
 * If mode == 1 then watch out for leading '.''s.
 * Otherwise mode should be 0.
 */
static int get_sign(FILE *filep, int mode)

{
    register int sign = 1;

    if (!LCHDEF) GETLCH;
    while ((LCH != EOF) && (LCH_NOT_TERM) && (!isdigit(LCH))
           && ((mode == 0) || (LCH != '.'))) {
        if (LCH == '-') sign *= -1;
        GETLCH;
    }
    return sign;
}


/* read_int()
 *
 * Read an integer from the stream.
 */
static int read_int(FILE *filep)

{
    int sign;
    register int dval;

    dval = 0;
    sign = get_sign(filep, 0);
    if ((LCH_IS_TERM) || (LCH == EOF)) { NARGS--; DEFTERM; return 0; }
    while ((LCH_NOT_TERM) && (isdigit(LCH))) {
        dval *= 10;
        dval += LCH - '0';
        GETLCH;
    }
    dval *= sign;
    if ((LCH == EOF) || (LCH_IS_TERM)) DEFTERM;
    return dval;
}


/* read_double()
 *
 * Read a double from the stream. Set LCH to -TERM
 * and return 0.0 if TERM is read before any digit.
 */
static double read_double(FILE *filep)

{
    int i;
    double pow = 1;
    int sign;
    double gval;

    gval = 0.0;
    sign = get_sign(filep, 1);
    if ((LCH_IS_TERM) || (LCH == EOF)) { DEFTERM; NARGS--; return 0; }
    while ((LCH_NOT_TERM) && (isdigit(LCH))) {
        gval *= 10.0;
        gval += LCH - '0';
        GETLCH;
    }
    if ((LCH_IS_TERM) || (LCH == EOF)) { DEFTERM; return gval*sign; }
    if (LCH == '.') {
        GETLCH;
        while ((LCH_NOT_TERM) && (isdigit(LCH))) {
            pow /= 10;
            gval += (LCH - '0') * pow;
            GETLCH;
        }
    }
    gval *= sign;
    if ((LCH == EOF) || (LCH_IS_TERM)) { DEFTERM; return gval; }

    if ((LCH == 'e') || (LCH == 'E')) {
        GETLCH;
        i = read_int(filep);
        if (i > 0) {
            for (pow = 1; i > 0; i--)
                pow *= 10;
        } else {
            for (pow = 1; i < 0; i++)
                pow /= 10;
        }
        gval *= pow;;
    }
    if ((LCH == EOF) || (LCH_IS_TERM)) DEFTERM;
    return gval;
}


/* read_string()
 *
 * Read a string from the stream, skipping leading
 * spaces or not according to mode. Terminate if
 * tch is read (excluding this from the string)
 * or if TERM is read.
 *
 * If tch is '\0', the terminating character is any
 * white space character.
 *
 * Set LCH to -TERM and return if TERM is read while
 * skipping leading spaces.
 *
 * modes: 0 <-> skip whitespace at beginning, 1 <-> don't
 */

static void read_string(FILE *filep, char *sval, char tch, int mode)

{
    if (filep == stdin && feof(filep)) {
        errorvf("input data stream closed.  But I need more data.\n");
        EXIT(FILE_READ_ERROR);
    }
    if (!LCHDEF) GETLCH;
    if (mode == 0) {
        while ((LCH_NOT_TERM) && (LCH != tch) && (isspace(LCH))) GETLCH;
        if (LCH_IS_TERM) { NARGS--; DEFTERM; return; }
    }
    while ((LCH != EOF) && (LCH_NOT_TERM) && (LCH != 13)
           && ((tch == '\0') ? (!isspace(LCH)) : (LCH != tch))) {
        *sval++ = (char)LCH;
        GETLCH;
    }
    if (LCH == tch) UNDEFLCH;
    else if ((LCH == EOF) || (LCH_IS_TERM)) DEFTERM;
    *sval = '\0';
}


/* fcmap()
 *
 * Process formatted input from the stream.
 * (Like fscanf().)
 *
 * Return the last character read, or -TERM if the
 * termination character was read before all data
 * was read, or LCH_UNDEF if a character was last read.
 */

#ifdef ANSI
int fcmap(FILE *filep, const char *fmt, ... )
#else
    int fcmap(filep, fmt, va_alist)
    FILE *filep;
    char *fmt;
    va_dcl
#endif
{
    va_list ap;
    register char *chp, tch;
    void *ptr;
    int mode = 0;

    if (filep == stdin) fflush(stdout);

    UNDEFLCH;
    UNDEFTCH;
    NARGS = 0;
    STAT = 0;

#ifdef ANSI
    va_start(ap, fmt);
#else
    va_start(ap);
#endif
    /* NULL is defined on Darwin as ((void *)0)
     * "a constant zero may be assigned to a pointer, and a pointer may be
     * compared with the constant zero.  The symbolic constant NULL is often
     * used in place of zero..."
     * 'NULL is sometimes defined as ((void *)0), which is a void pointer to the number 0.
     * Yes, it resolves to 0, but it's generally meant in that case to be used to compare
     * with pointers, not chars since a char is not a pointer.'
     * 'The C compiler can define NULL as (void*)0 or simple 0.
     * Moral: NULL is not a zero alias in C, use NULL only for pointer values.'
     * */
    while (*fmt != '\0') {
        NARGS++;
        if (TERMDEF)
            break;
        if (*fmt == '%') {
            switch (*(++fmt)) {
            case 'd':
                ptr = (void *) va_arg(ap, int *);
                *(int *) ptr = read_int(filep);
                break;
            case 'f':
                ptr = (void *) va_arg(ap, float *);
                *(float *) ptr = (float) read_double(filep);
                break;
            case 'g':
                ptr = (void *) va_arg(ap, double *);
                *(double *) ptr = read_double(filep);
                break;
            case 'c':
                chp = va_arg(ap, char *);
                if (!LCHDEF)
                    GETLCH
                        ;
                *chp = (char)LCH;
                UNDEFLCH;
                break;
            case 'S':
                mode = 1;
            case 's':
                tch = *(++fmt);
                if (tch == '\0')
                    fmt--;
                /* try to be smart about this */
                else if (tch == '%') {
                    if (*(++fmt) != '%') {
                        tch = '\0';
                        fmt--;
                    }
                } else if (tch == ' ')
                    tch = '\0';
                read_string(filep, va_arg(ap, char *), tch, mode);
                mode = 0;
                break;
            case '=':
                LCH = va_arg(ap, int);
                /* The call to va_arg above returns an 'int', which is stored in LCH */
                if (LCH == '\0')
                    UNDEFLCH;
                else DEFLCH;
                break;
            case 'T':
                TERM = va_arg(ap, int);
                if (TERM == '\0')
                    UNDEFTCH;
                else DEFTCH;
                break;
            default:
                skip_past_char(filep, *fmt);
                break;
            }
        } else if (*fmt != ' ')
            skip_past_char(filep, *fmt);
        fmt++;
    }
    va_end(ap);
    if (!LCHDEF)
        return 0;
    return LCH;
}
