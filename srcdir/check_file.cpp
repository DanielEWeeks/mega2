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
/* check_file.c */

/************************************************************************

  function to check execution command and file

  status:

  date: 11/14/89

  Modified - Nandita
  Date     - May 1999
           - added an argument for mode of opening outfile
*************************************************************************/

#include <stdio.h>


int check_file (int argc, const char *argv[], FILE **infile, FILE **outfile, const char *write_mode)

{

    int exit_value = 0;

    /* file manipulation */

    if ( argc != 3 ) {
        fprintf (stderr, "Bad command.\nCorrect usage: hwe infile outfile.\n");
        exit_value = 1;
    }

    if ( ( *infile = fopen (argv[1], "r")) == (FILE *) NULL ) {
        fprintf (stderr, "Can't read %s\n", argv[1]);
        exit_value = 2;
    }

    if ( ( *outfile = fopen (argv[2], write_mode)) == (FILE *) NULL ) {
        fprintf (stderr, "Can't write %s\n", argv[2]);
        exit_value = 3;
    }

    return (exit_value);

}
