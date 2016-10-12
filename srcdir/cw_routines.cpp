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
/* functions that have to be defined for codewarrior */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#ifndef DMALLOC
#if (defined(SOLARIS) || defined(DEBUG_CW) || defined(DARWIN_OS) || defined(MINGW))
void cfree(void *ptr);
#endif
#endif

#if defined(__MWERKS__) || defined(DEBUG_CW)

#ifndef CW_DEFS_H
#define CW_DEFS_H

#define R_OK    4       /* Test for Read permission */
#define W_OK    2       /* Test for Write permission */
#define X_OK    1       /* Test for eXecute permission */
#define F_OK    0       /* Test for existence of File */
#define EOSTR   '\0'
#endif

static int safe_tolower(int c)
{
    if (isupper(c))
        return tolower(c);
    else
        return c;
}

#ifndef __INTEL__
int access(const char*path, int amode)

{
    FILE *fp;
    int retval;

    switch(amode) {
    case F_OK:
        retval=((fp=fopen(path, "r+"))? 0 : -1);
        break;
    case R_OK:
        retval=((fp=fopen(path, "r"))? 0 : -1);
        break;
    case W_OK:
        retval=((fp=fopen(path, "a"))? 0 : -1);
        break;
    case X_OK:
        retval=((fp=fopen(path, "r+"))? 0 : -1);
        break;
    default:
        fprintf(stderr, "ERROR: unknown mode!\n");
        EXIT(FILE_WRITE_ERROR);
    }
    if (retval==0) fclose(fp);
    return retval;
}
#endif

int strcasecmp(char *s1, char *s2)

{

    char *sp1, *sp2;

    for (sp1=&s1[0], sp2=&s2[0]; *sp1 != EOSTR && *sp2 != EOSTR; sp1++, sp2++) {

        if (safe_tolower(*sp1) > safe_tolower(*sp2)) return 1;
        if (safe_tolower(*sp1) < safe_tolower(*sp2)) return -1;
    }

    if (*sp1!=EOSTR && *sp2 == EOSTR) return 1;
    if (*sp1==EOSTR && *sp2 != EOSTR) return -1;
    if (*sp1==EOSTR && *sp2==EOSTR) return 0;
}

int strncasecmp(char *s1, char*s2, size_t n)

{
    char *sp1, *sp2;
    int i;
    for (sp1=s1, sp2=s2, i=0; i<n; sp1++, sp2++, i++) {

        if (safe_tolower(*sp1) > safe_tolower(*sp2)) return 1;
        if (safe_tolower(*sp1) < safe_tolower(*sp2)) return -1;
    }
    return 0;
}

char *strdup(char *s1)

{
    char *s2;
    size_t n;

    n=strlen(s1);
    s2=MALLOC_STR((n+1), char);
    strncpy(s2, s1,n);
    s2[n]=EOSTR;
    return s2;
}
#endif

#ifndef DMALLOC
#if (defined(SOLARIS) || defined(DEBUG_CW) || defined(DARWIN_OS) || defined(MINGW))
void cfree(void *ptr)
{
    free(ptr);
}

#endif
#endif

#if (defined(SOLARIS) || defined(MINGW) || defined(_WIN))

/*
 * Get next token from string *stringp, where tokens are possibly-empty
 * strings separated by characters from delim.
 *
 * Writes NULs into the string at *stringp to end tokens.
 * delim need not remain constant from call to call.
 * On return, *stringp points past the last NUL written (if there might
 * be further tokens), or is NULL (if there are definitely no more tokens).
 *
 * If *stringp is NULL, strsep returns NULL.
 */
char *strsep(char **stringp, const char *delim)
{
    char *s;
    const char *spanp;
    int c, sc;
    char *tok;

    if ((s = *stringp) == NULL)
        return (NULL);
    for (tok = s;;) {
        c = *s++;
        spanp = delim;
        do {
            if ((sc = *spanp++) == c) {
                if (c == 0)
                    s = NULL;
                else
                    s[-1] = 0;
                *stringp = s;
                return (tok);
            }
        } while (sc != 0);
    }
    /* NOTREACHED */
}
#endif

#ifdef DEBUG_CW

main(int argc, char *argv[])
{
    /*  printf("%d\n", access("mega2.c_copy", F_OK)); */
    printf("%d %s\n", strcasecmp(argv[1], argv[2]), argv[1]);
    printf("%d\n", strncasecmp(argv[1], argv[2], atoi(argv[3])));
}

#endif
