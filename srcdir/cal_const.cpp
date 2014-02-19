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

/**************************************************************************

  program name: cal_const.c

  function to compute the constant part of the probability function for
  testing the H-W equilibrium

  constant = log N! - log (2N)! + sum(i) log n(i)!

  status: OK

  date: 12/3/89

**************************************************************************/


#include <stdio.h>

#include "common.h"
#include "typedefs.h"

#include "hwe_ext.h"
/*
                hwe_ext.h:  log_factorial
*/


double cal_const (int no_allele, int n[MAX_ALLELE], int total)

/* int no_allele;
   int n[MAX_ALLELE];
   int total; */

{
    double constant;
    register int i;
/*  double log_factorial(int k); */

    constant = log_factorial ( total ) - log_factorial ( 2*total );

    for ( i = 0; i < no_allele; ++i )
        constant += log_factorial ( n[i] );

    return ( constant );

}
