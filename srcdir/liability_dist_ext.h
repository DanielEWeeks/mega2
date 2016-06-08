/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef LIABILITY_DIST_EXT_H
#define LIABILITY_DIST_EXT_H

extern char **get_liablfl_name(int *numchr, int tab_text);

extern void get_output_options(int *ou_percent, int *empty_cols,
			       int *empty_rows, int *tab_text);


extern void liability_summary(char **input_files,
			      linkage_ped_top *LPedTreeTop,
			      int *numchr, int untyped_ped_opt);


extern void write_liable_dist(liable_allele_dist *allele_dist,
			      char **liablefl_name,
			      linkage_ped_top *LPedTreeTop, int num_AFFECTION,
			      int num_NUMBERED, int ou_percent,
			      int empty_rows, int empty_cols, int tab_text);


#endif
