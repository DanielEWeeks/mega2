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

#ifndef PEDTREE_EXT_H
#define PEDTREE_EXT_H

extern ped_top *convert_to_pedtree(linkage_ped_top *Top,
				   int assign_affs);

extern void free_all_including_ped_top(ped_top *Top,
				       void (*PTmpFreeFun)(void *TmpData),
				       void (*ETmpFreeFun)(void *TmpData));

extern int is_typed(ped_rec *Entry, linkage_locus_top *LTop, int mode);

extern void pedtree_markers_check(linkage_ped_top *LPedTreeTop, analysis_type analysis);

extern void reassign_affecteds(ped_top *Top);

extern void remove_untyped_affecteds(ped_tree *Ped, linkage_locus_top *LTop, int mode,
				     int *order, int n);

extern int renumber_ped(ped_tree *Ped);


/*
extern void add_to_offspring(ped_rec *Parent, ped_rec *Child);

extern void clear_ped_tree(ped_tree *Ped);

extern void free_all_from_locus_top(locus_top *LTop);

extern void free_all_from_ped_tree(ped_tree *Ped,
				   void (*PTmpFreeFun)(void *TmpData),
				   void (*ETmpFreeFun)(void *TmpData));

extern ped_top *new_ped_top(void);

*/

#endif
