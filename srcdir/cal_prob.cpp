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

/***************************************************************************

  program name: cal_prob.c

  function to calculate

  status:

  date: 12/7/89

***************************************************************************/

#include <stdio.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "do_switch_ext.h"
#include "test_switch_ext.h"
#include "utils_ext.h"
/*
          do_switch_ext.h:  do_switch
        test_switch_ext.h:  test_switch
              utils_ext.h:  randomnum
*/


double cal_prob (int a[LENGTH], Index index,double ln_p_old, int *actual_switch )

/* int a[LENGTH];
   Index index;
   double ln_p_old;
   int *actual_switch; */

{

    double p1_ratio, p2_ratio;
    register double ln_p_new;
    double rand_num;
    int switch_ind, type;

    *actual_switch = 0;

/* determine the switchability and direction of switch for given face */

    test_switch( a, index, &switch_ind, &type, &p1_ratio, &p2_ratio);
/*  printf("%d %d %d swhind=%d\n",index.i1, index.i2,index.i3,switch_ind);*/

    switch (switch_ind)
        {
        case 0:           /* non-switchable */

            ln_p_new = ln_p_old;  /* retain the pattern, probability unchanged */
            break;

        case 1:           /* partially-switchable */

            if ( type == 1 )
                p1_ratio = p2_ratio;
            rand_num = randomnum();

            if ( rand_num < TRANS( p1_ratio ) ) {  /* switch w/ transition P TRANS */
                do_switch ( a, index, type );
                ln_p_new = ln_p_old + log (p1_ratio);  /* ln P_after-switch */
                *actual_switch = 1;
            } else                   /* remain the same w/ P = 1 - TRANS */
                ln_p_new = ln_p_old;           /* probability unchanged */
            break;

        default:          /* fully switchable */
            rand_num = randomnum();

            if ( rand_num <= TRANS(p1_ratio)) {
                do_switch ( a, index, 0 ); /* D-switch */
                ln_p_new = ln_p_old + log (p1_ratio);  /* ln P_after-switch */
                *actual_switch = 2;
            } else if ( rand_num <= TRANS(p1_ratio) + TRANS(p2_ratio) ) {
                do_switch ( a, index, 1 ); /* R-switch */
                ln_p_new = ln_p_old + log (p2_ratio);
                *actual_switch = 2;
            } else
                ln_p_new = ln_p_old;
            break;
        }

    return (ln_p_new);
}
