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

#ifndef REORDER_LOCI_EXT_H
#define REORDER_LOCI_EXT_H

extern void  Free_reorder_loci(void);

extern linkage_ped_top *ReOrderLoci(linkage_ped_top *Top, int *numchr,
				    analysis_type *analysis);

extern int   ReOrderMappedLoci(linkage_ped_top *Top, int *numchr);


extern int   STR_CHR(const char *dummy);


extern int   display_selections(linkage_ped_top *LocusTop, int *entries,
                                int num_entries, char *asterisks,
                                char *scanned_str, int traits_only, int map_num);

extern int   get_chromosome_list(linkage_locus_top *LTop, int *lst,
                                 int *counts,  int exit_no_chr,
				 analysis_type analysis);

extern void  get_trait_list(linkage_locus_top *LTop, const int check_traits_combine);

extern int   parse_string(char *entrydummy, int *number,
                          int *num_entries, int num_total_loci, int require_e);

extern int   x_linked_check(int numchr, int *chromosomes, analysis_type analysis);

/*
   extern int   *advance_trp(int *current); */

extern int  *check_for_duplicates(int *entrycount, int *number);


/* extern int   *init_trp(void);

extern int   linkage_set_locus_order(linkage_ped_top * Top, int *order);

*/

#endif
