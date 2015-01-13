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

/* random_choose.c */

/***********************************************************************

  function to randomly choose two integer numbers, k1 and k2, between 0
  and k - 1.  ( 0 <= k1 < k2 < k )

  status: OK

  date: 9/21/89

************************************************************************/

#include <stdio.h>

#include "common.h"
#include "typedefs.h"

#include "utils_ext.h"
/*
              utils_ext.h:  randomnum
*/


void random_choose (int *k1, int *k2,int k)
{
    register int temp, i, not_find;
    int work[MAX_ALLELE];

    for ( i = 0; i < k; ++i )
        work[i] = i;

    *k1 = (int) (randomnum() * (double) k);

    --k;

    for ( i = *k1; i < k; ++i )
        work[i] = i + 1;

    not_find = 1;

    while ( not_find ) {
        i = (int) (randomnum() * (double) k);
        *k2 = work[i];
        not_find = 0;
    }

    if ( *k1 > *k2 ) {
        temp = *k1;
        *k1 = *k2;
        *k2 = temp;
    }

}
