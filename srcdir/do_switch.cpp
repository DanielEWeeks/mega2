/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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
/*********************************************************************

  program name: do_switch.c

  function to make switch according to given switchability and switch
  type.

  status:  in progress

  date: 1/14/90

**********************************************************************/

#include "common.h"
#include "typedefs.h"


void do_switch ( int a[LENGTH], Index index, int type );

void do_switch ( int a[LENGTH], Index index, int type )

{
    int k11, k22, k12, k21;

    k11 = L(index.i1, index.j1);
    k12 = L(index.i1, index.j2);
    k21 = L(index.i2, index.j1);
    k22 = L(index.i2, index.j2);


    if ( type == 0 ) {  /* D-switch */
        --a[k11];
        --a[k22];
        ++a[k12];
        ++a[k21];
    } else {     /* R-switch */
        ++a[k11];
        ++a[k22];
        --a[k12];
        --a[k21];
    }
}
