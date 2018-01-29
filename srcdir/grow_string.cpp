/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

#ifdef MAIN

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "common.h"
#define MAX_NAMELEN 400
#define EXIT exit
#define SYSTEM_ERROR -14
#else

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "utils_ext.h"
/*
              utils_ext.h:  EXIT
*/
#endif


/* This is a function to append a formatted string to
   and existing one (modeled after fcmap).
   IMPORTANT: the variable arg-list should not include
   the existing string that is being appended to.
   Modifies original string.
   Now only handles format strings of the type %{-|.}[0-9]*{s|c|i|f}.
   Will NOT HANDLE %%, as it expects % to be follwoed by -, . or a number, then
   followed by s,c,i,f.

*/


void grow(char *str, const char *format, ...)

{
    char result[MAX_NAMELEN]="", temp[MAX_NAMELEN], format_str[20];

    char *fmp, *s;
    int i;
    double f;
    char c;
    va_list ap;

    va_start(ap, format);
    format_str[0]='\0';
    fmp = &(format_str[0]);

    while(*format) {
        strcpy(temp, "");
        if (*format == '%') {
            *fmp=*format;
            format++;
            while(isdigit((unsigned char)*format) || *format == '.' || *format == '-') {
                *(++fmp) = *format; format++;
            }

            if (*format == 's' || *format == 'c' || *format == 'f' || *format == 'd') {
                *(++fmp)=*format;
                *(++fmp) = '\0';

                switch(*format) {
                case 's':
                    s= va_arg(ap, char*);
                    sprintf(temp, format_str, s);
                    break;
                case 'd':
                    i= va_arg(ap, int);
                    sprintf(temp, format_str, i);
                    break;
                case 'c':
		    c = (char)va_arg(ap, int);
                    sprintf(temp, format_str, (char) c);
                    break;
                case 'f':
                    f=va_arg(ap, double);
                    sprintf(temp, format_str, f);
                    break;
                }
                strcat(result, temp);
                format_str[0] = '\0';
                fmp = &(format_str[0]);
            } else {
                printf("Unrecognized format, aborting!\n");
                EXIT(SYSTEM_ERROR);
            }
        } else {
            sprintf(temp, "%c", *format);
            strcat(result, temp);
        }
        format++;
    }
    strcat(str, result);
    va_end(ap);
    return;

}



#ifdef MAIN
int main()

{
    char  mssg[600];
    int a=10, b =20;
    char c='*';
    int i, id[10];

    strcpy(mssg, "Result of adding ");
    grow(mssg, "%-15s %d + %d = %d %c, remainder of division a/b = %10.5f", "numbers",
         a, b, a+b, c, (double)a/b);
    printf("%s\n", mssg);
    strcpy(mssg, " ");
    for (i=1; i <= 10; i++) {
        id[i-1] = i;
    }
    for (i=0; i < 10; i++) {
        grow(mssg, " %d", id[i]);
    }
    printf("%s\n", mssg);
    return(1);
}

#endif
