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

/************************************************************************

  program name: read_data.c

  funtion to read data

  status: in progress

  date: 12/5/89

*************************************************************************/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"


int read_data (int *a, int *no_allele, int *total,
               struct randomization *sample, FILE **infile)

{
    register int i, j, l, err = 1;

    *total = 0;

    if ( fscanf(*infile, "%d", no_allele) != 1) {
        fprintf(stderr, "Please supply number of alleles\n");
        return ( err );
    }

    if ( *no_allele < 3 ) {
        fprintf(stderr, "***Error! Number of alleles less than 3. \n");
        return ( err );
    }

    for ( i = 0; i < *no_allele; ++i ) {
        for ( j = 0; j <= i; ++j ) {
            l = LL(i, j);
            (void)fscanf (*infile, "%d ", &a[l]);
            *total += a[l];
        }
    }

    if ( fscanf(*infile, "%d %d %d \n", &sample->step,
               &sample->group, &sample->size) != 3 ) {
        fprintf( stderr, " Please supply parameters.\n" );
        return ( err );
    }

    if ( sample->step < 1 || sample->group <= 1 ) {
        fprintf( stderr, "***Error in parameter specification.\n" );
        return ( err );
    }

    return ( 0 );

}
