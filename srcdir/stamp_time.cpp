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

/***********************************************************************

  program name: stamp_time.c

  function to stamp the date and time the program is executed.


  Status:

  Date: 12/10/89

************************************************************************/

#include <stdio.h>
#include <time.h>

#ifdef HIDEDATE
extern const char *NOcTIME;
#endif

#ifdef __INTEL__
void stamp_time(long t11, FILE *outfile)
{
    long t2, now;
#else
    void stamp_time(time_t t11, FILE *outfile)
{
    time_t t2, now;
#endif
    time(&t2);
    t2 -= t11;
    time(&now);

#ifdef HIDEDATE
    fprintf(outfile, "\nTotal elapsed time: %d''\n", (int)0);
    fprintf(outfile, "Date and time: %s\n", NOcTIME);
#else
    fprintf(outfile, "\nTotal elapsed time: %d''\n", (int)t2);
    fprintf(outfile, "Date and time: %s\n", ctime(&now));
#endif
}
