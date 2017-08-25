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

#ifndef HWE_EXT_H
#define HWE_EXT_H

extern double cal_prob(int a[], Index indx,double ln_p_old, int *actual_switch);


extern void hwe_test(char *marker_name, const char *infile_name, const char *outfile_name,
		     double homo_o, double homo_e, int num_genos);

extern double log_factorial(int k);

extern void print_data(int *a, int no_allele, struct randomization sample,
                       FILE **outfile);

extern void select_index(Index *indx, int no_allele);

extern void stamp_time(time_t t11, FILE *outfile);


/*
extern double cal_const(int no_allele, int n[], int total);

extern void cal_n(int no_allele,int a[], int n[]);

extern void get_interval(void);

extern double ln_p_value(int a[], int no_allele, double constant);


extern double rans(void);

*/

#endif
