/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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
/***************************************************************************

  program name: hwe.h

  header file for hwe.c

  status:

  date: 1/14/90

***************************************************************************/

#ifndef HWE_H
#define HWE_H

#define STR_END  ('\0')
#define  MAX_ALLELE (50)
#define  LENGTH     (MAX_ALLELE * ( MAX_ALLELE + 1 ) / 2)

#define  MIN(x, y)     ((x) < (y)) ? (x) : (y)
#define  RATIO(u, v)   ( (double) (u) ) / ( 1.0 + (double) (v) )
#define  TRANS(x)     (MIN(1.0, x))/2.0  /* transition probability */

#define  LL(a, b)      a * ( a + 1 ) / 2  + b
#define  L(a, b)       ( a < b ) ? b*(b + 1)/2 + a : a*(a+1)/2 + b

#define  EXPECT(a,b,c) ((double) a) / ((double) c) * ((double) b) / 2.0


typedef struct _Index
{
    int i1;
    int i2;
    int j1;
    int j2;
    int type;
    double cst;
} Index;

struct outcome
{
    double p_value;  /* mean p-value */
    double se;       /* standard error of the p-value */
    int swch_count[3];  /* switch counts for partial and full switch */
};

struct randomization
{
    int group; /* total number of chunks */
    int size;  /* size of a chunk */
    int step;  /* number of steps to de-memerization */
};

#endif
