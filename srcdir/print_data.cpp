/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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

/*************************************************************************/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"


void print_data(int *a, int no_allele, struct randomization sample,
		FILE **outfile)

{

    register int i, j, k, l;
    char line[FILENAME_LENGTH];

    line[0] = '-';

    k = 1;

    fprintf (*outfile, "Observed genotype frequencies: \n\n");

    for ( i = 0; i < no_allele; ++i ) {

        for ( j = k; j < k + 5; ++j )
            line[j] = '-';

        line[j] = STR_END;
        k = j;

        fprintf (*outfile, "%s\n", line);

        fprintf (*outfile, "|");

        for ( j = 0; j <= i; ++j ) {
            l = LL(i, j);
            fprintf(*outfile, "%4d|", a[l]);
        }

        fprintf (*outfile, "\n");
    }

    fprintf (*outfile, "%s\n\n", line);
    fprintf (*outfile, "Total number of alleles: %2d\n\n", no_allele);

    fprintf(*outfile, "Number of initial steps: %d\n", sample.step);
    fprintf(*outfile, "Number of chunks: %d\n", sample.group);
    fprintf(*outfile, "Size of each chunk: %d\n\n", sample.size);

}
