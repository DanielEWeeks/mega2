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

/*************************************************************************

  program name: ln_p_value.c

  function to compute the log p value for given genetype frequencies

  status:

  date: 12/6/89

*************************************************************************/

#include <stdio.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "hwe_ext.h"
/*
                hwe_ext.h:  log_factorial
*/


double ln_p_value (int a[LENGTH], int no_allele, double constant)

/*int a[LENGTH];
  int no_allele;
  double constant;*/

{
    register int i, j, l, temp;
    register double ln_prob;
/*  double log_factorial(int k); */

    ln_prob = constant;
    temp = 0;

    for ( i = 0; i < no_allele; ++i ) {
        for ( j = 0; j < i; ++j ) {
            l = LL(i, j);
            temp += a[l];
            ln_prob -= log_factorial ( a[l] );
        }

        l = LL(i, i);
        ln_prob -= log_factorial ( a[l] );
    }

    ln_prob += temp * log ( 2.0 );

    return ( ln_prob );

}
