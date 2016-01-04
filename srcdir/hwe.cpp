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

/***********************************************************************

  program name: hwe.c

  main program

  Status:

  Date: 1/14/90


  Modified by - Nandita
  Date        - May 1999
              -Integrated with mega2
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include "hwe.h"

#include "check_file_ext.h"
#include "hwe_ext.h"
#include "read_data_ext.h"
/*
         check_file_ext.h:  check_file
                hwe_ext.h:  cal_prob print_data select_index stamp_time
          read_data_ext.h:  read_data
*/


#ifdef STANDALONE

main(argc,argv)  /* correct execution of this program:
                    hwe infile outfile               */

int argc;
char *argv[];
#else

void hwe_test(char *marker_name, const char *infile_name, const char *outfile_name,
	      double homo_o, double homo_e, int num_genos)

#endif

{

    int a[LENGTH];
    double ln_p_observed, ln_p_simulated, p_mean, p_square;
    double p_simulated, total_step;
    int no_allele, total, counter, actual_switch;
    Index index;
    struct randomization sample;
    struct outcome result;
    FILE *infile, *outfile;

    const char *append_mode = "a";

    register int i, j;
#ifdef __INTEL__
    long t1;
#else
    time_t t1;
#endif

#ifdef STANDALONE
    if ( check_file ( argc, argv, &infile, &outfile, write_mode ) )
        exit ( 1 );
#else
    int num_args = 3;
    const char *filenames[3];

    filenames[0] = NULL;
    filenames[1] = &(infile_name[0]);
    filenames[2] = &(outfile_name[0]);
    if ( check_file ( num_args, filenames, &infile, &outfile, append_mode) )
        exit ( 1 );
    fprintf(outfile, "=============\n%s\n============\n",
            marker_name);
    fprintf(outfile, "Observed homozygosity: %5.4f\n", homo_o);
    fprintf(outfile, "Expected homozygosity: %5.4f\n", homo_e);
    fprintf(outfile, "Genotypes counted:     %d\n\n", num_genos);

#endif


    if ( read_data ( a, &no_allele, &total, &sample, &infile ) )
        exit (2);

    print_data ( a, no_allele, sample, &outfile );

    ln_p_observed = 0.0;

    ln_p_simulated = ln_p_observed;

    p_mean = p_square = (double) 0.0;

    result.p_value = result.se = (double) 0.0;  /* initialization */

    result.swch_count[0] = result.swch_count[1] = result.swch_count[2] = 0;

    for ( i = 0; i < sample.step; ++i ) {  /* de-memorization for given steps */

        select_index ( &index, no_allele );

        ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);

        ++result.swch_count[actual_switch];
    }

    for ( i = 0; i < sample.group; ++i ) {

        counter = 0;

        for ( j = 0; j < sample.size; ++j ) {

            select_index ( &index, no_allele );

            ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);

            if ( ln_p_simulated <= ln_p_observed )
                ++counter;

            ++result.swch_count[actual_switch];

        }
        p_simulated = (double) counter  / sample.size;
        p_mean += p_simulated;
        p_square += p_simulated * p_simulated;

    }

    p_mean /= sample.group;
    result.p_value = p_mean;
    result.se = p_square / ((double) sample.group)/(sample.group - 1.0)
        - p_mean / ( sample.group - 1.0 ) * p_mean;
    result.se = sqrt ( result.se );

    total_step = sample.step + sample.group * sample.size;

    fprintf(outfile, "Randomization test P-value: %7.4g  (%7.4g) \n",
            result.p_value, result.se);
    fprintf(outfile, "Percentage of partial switches: %6.2f \n",
            result.swch_count[1] / total_step * 100);
    fprintf(outfile, "Percentage of full switches: %6.2f \n",
            result.swch_count[2] / total_step * 100);
    fprintf(outfile, "Percentage of all switches: %6.2f \n",
            (result.swch_count[1] + result.swch_count[2]) / total_step * 100 );

    time(&t1);
    stamp_time ( t1, outfile );

    fclose(infile);   /* added by nandita */

    fclose(outfile);
}
