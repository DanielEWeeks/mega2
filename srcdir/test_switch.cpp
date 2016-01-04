/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

/***************************************************************************

  program name: test_switch.c

  function to determine the switchibility

  switch_ind = 0 if non-switchable, 1 if partially-switchable, 2 if switchable.

  And it returns the switch type if switchable.

  switch_type = 0 if D-switchable and 1 if R-switchable;

  But switch_dir = 0 if switch_ind = 2.

  In addition, it returns the probability ratio if appropriete.

  status:

  date: 1/14/90

****************************************************************************/

#include <stdio.h>

#include "common.h"
#include "typedefs.h"


void test_switch(int a[LENGTH], Index index, int *switch_ind, int *switch_type, /* switchability and type of switch */
                 double *p1_rt, double *p2_rt) /* probability ratio */
{
    register int k11, k22, k12, k21;

    *switch_ind = 0;

    k11 = L(index.i1, index.j1);
    k22 = L(index.i2, index.j2);
    k12 = L(index.i1, index.j2);
    k21 = L(index.i2, index.j1);

    if ( index.type <= 1 ) { /* type = 0, 1 */
        if ( a[k11] > 0 && a[k22] > 0 ) {
            *switch_ind = 1;
            *switch_type = 0; /* D-switchable */
            *p1_rt = RATIO(a[k11], a[k12]) *  RATIO(a[k22], a[k21]) * index.cst;
        }
        if ( a[k12] > 0 && a[k21] > 0 ) {
            *switch_ind += 1;
            *switch_type = 1; /* R-switchable */
            *p2_rt = RATIO(a[k12], a[k11]) *  RATIO(a[k21], a[k22]) / index.cst;
        }

    } else {                  /* type = 2 */
        if ( a[k11] > 0 && a[k22] > 0 ) {
            *switch_ind = 1;
            *switch_type = 0; /* D-switchable */
            *p1_rt = RATIO(a[k11],a[k12] + 1.0)*RATIO(a[k22],a[k12]) * index.cst;
        }
        if ( a[k12] > 1 ) {
            *switch_ind += 1;
            *switch_type = 1; /* R-switchable */
            *p2_rt = RATIO(a[k12],a[k11]) * RATIO(a[k12] - 1,a[k22]) / index.cst;
        }

    }

}
