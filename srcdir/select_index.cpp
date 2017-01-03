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

  program name: select_index.c

  function to randomly choose three integers i1, i2, i3 with
  0 <= i1 < i2 <= no_allele, 0 <= i3 <= no_allele

  status: OK

  date: 1/14/90

************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "random_choose_ext.h"
/*
      random_choose_ext.h:  random_choose
*/


void select_index (Index *index, int no_allele)

/* Index *index;
   int no_allele; */

{


    int i1, i2, j1, j2;
    int k = 0;
    int l = 0;

/* generate row indices */

    random_choose ( &i1, &i2, no_allele );

    index->i1 = i1;
    index->i2 = i2;

/* generate column indices */

    random_choose ( &j1, &j2, no_allele );

    index->j1 = j1;
    index->j2 = j2;

/* calculate Delta = d(i1,j1) + d(i1,j2) + d(i2,j1) + d(i2,j2) */

    if ( i1 == j1 )
        ++k;

    if ( i1 == j2 )
        ++k;

    if ( i2 == j1 )
        ++k;

    if ( i2 == j2 )
        ++k;

    index->type = k;

    if ( ( i1 == j1 ) || ( i2 == j2 ) )
        ++l;

    index->cst = ( l == 1 ) ? pow(2.0, (double) k) : pow(2.0, - (double) k);
}
