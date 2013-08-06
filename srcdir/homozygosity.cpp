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

/***********************************************************/
/*                homozygosity.c                           */
/*          computes observed and expected homozygosities  */
/*          Author; Nandita Mukhopadhyay, June 1999        */
/***********************************************************/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "genetic_utils_ext.h"
/*
      genetic_utils_ext.h:  genoindx
*/


double expected_homozygosity(linkage_ped_top *LPedTreeTop, int locus_id)

{
    int ii, allele_cnt;
    linkage_locus_top *LocusTop;
    double allele_freq, homozygosity;

    LocusTop=LPedTreeTop->LocusTop;
    allele_cnt=LocusTop->Locus[locus_id].AlleleCnt;
    homozygosity=0.0;

    for (ii=0; ii<allele_cnt; ii++) {
        allele_freq=LocusTop->Locus[locus_id].Allele[ii].Frequency;
        homozygosity += allele_freq*allele_freq;
    }
    return homozygosity;
}

/* NM- Nov 2007, only count females if xlinked */

double observed_homozygosity(int *geno_count, int num_alleles)

{
    /* works off the geno_count data structure created inside hwe_user_input
       for one locus */
    int ii, num_homozygotes, total_indivs, num_genos;
    double homozygosity;

    /* total individuals */
    total_indivs=0;
    num_genos=NUMGENOS(num_alleles);
    for (ii=0; ii<num_genos; ii++)
        total_indivs += geno_count[ii];

    /* frequency of homozygotes */
    num_homozygotes=0;
    for (ii=0; ii<num_alleles; ii++)
        num_homozygotes += geno_count[genoindx(ii+1,ii+1)];

    /* percentage homozygotes */
    if (total_indivs == 0)
        homozygosity= 0.0;
    else
        homozygosity=(double) num_homozygotes/(double) total_indivs;

    return homozygosity;

}
