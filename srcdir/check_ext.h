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

#ifndef CHECK_EXT_H
#define CHECK_EXT_H

extern int check_half_type(ped_tree *PedTree, ped_status *PedStatus,
                           locus_top *LTop1, int ped_num, int locus,
                           int *display_error, int uniqueids,
                           FILE **reset_fp, bool *first);

extern int check_invalid_fam(ped_tree *PedTree, ped_status *PedStatus,
                             locus_top *LTop1, int ped_num, int locus1,
                             int *display_error, int uniqueids,
                             ped_rec **Sibs);

extern int check_locus(locus_rec *Locus, int *disp_err,
		       analysis_type analysis,
		       int *plink_locus_num);

extern int check_out_of_bounds(ped_tree *PedTree, ped_status *PedStatus,
                               locus_top *LTop1, int ped_num, int locus,
                               int *display_error, int uniqueids,
                               FILE **reset_fp, bool *first);

extern int check_ped_relations(ped_tree *PedTree, ped_status *PedStatus);

extern int check_unique_ids(linkage_ped_top *Top);

#endif
