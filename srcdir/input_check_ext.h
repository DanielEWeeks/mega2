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

#ifndef INPUT_CHECK_EXT_H
#define INPUT_CHECK_EXT_H

extern void            check_renumber_ped(ped_tree *Ped);

extern void            clear_ped_status(ped_status *Stat);


extern void            full_check(ped_top *Top, linkage_ped_top *LPedTop,
				  analysis_type analysis);

extern void             mito_transmission_report(ped_top *PTop,
						 int *num_hetero, int *num_non_maternal);


/*
extern int             UnMappedLociCheck(linkage_ped_top *Top);

*/

#endif
