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

#ifndef WRITE_APM_EXT_H
#define WRITE_APM_EXT_H



extern void    create_APMULT(char *pedfl_name, char *locusfl_name, int disease_locus,
			     const char *affdata_strg, linkage_ped_top **LPedTreeTop,
			     ped_top * PedTreeTop,
			     file_format * infl_type, file_format * outfl_type,
			     int *numchr, char *mapfl_name,
			     analysis_type * analysis, char *fl_names[],
			     int untyped_ped_opt);


extern void  create_APM_file(char *pedfl_name, char *locusfl_name, int disease_locus,
			     const char *affdata_strg, linkage_ped_top * LPedTreeTop,
			     ped_top * PedTreeTop,
			     file_format * infl_type, file_format * outfl_type,
			     int *numchr, char *mapfl_name, char *fl_names[],
			     analysis_type * analysis, char *omitfl_name,
			     int untyped_ped_opt);


#endif
