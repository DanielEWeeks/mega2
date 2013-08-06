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

#ifndef SELECT_INDIVIDUALS_EXT_H
#define SELECT_INDIVIDUALS_EXT_H

extern int founders_or_everyone(linkage_ped_top *Top, int locus_id,
				record_type rec, int ped,
				int *member_ids,
				int count_option, int inc_ht, int xlinked);

extern int get_count_option(int halftyped_item, int *include_halftyped, const char *messg);

extern int random_ped_member(linkage_ped_top *LPedTreeTop,
			     int locus_id, int ped, record_type rec,
			     int inc_ht, int *num_alleles, int xlinked);


extern void select_individuals(linkage_ped_top *LPedTreeTop, int locus_id,
			       int option, int ped, int *member_ids,
			       int inc_ht, int xlinked);


#endif
