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

#ifndef MAKENUCS_EXT_H
#define MAKENUCS_EXT_H


extern void  create_nuclear_families(linkage_ped_top **Top,
				     char *outfl_name,
				     char *loutfl_name,
				     analysis_type * analysis,
				     file_format *infl_type,
				     int untyped_ped_opt,
				     linkage_ped_top **NukeTop);

extern void  create_nuked_fids(linkage_ped_top * Topp);

extern void makenucs1(linkage_ped_top *Top, ped_top *PedTop,
		      linkage_ped_top *NewTop);


#endif
