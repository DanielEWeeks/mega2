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

#ifndef MARKER_LOOKUP_EXT_H
#define MARKER_LOOKUP_EXT_H

void clear_marker();

extern int test_and_add_marker(const char *key, int value);

extern void add_marker(const char *key, int value);

extern void make_marker(linkage_locus_top *LTop);

extern int search_marker(const char *marker, int *indx);

extern void debug_dump_marker_dict();


extern void add_allele(const char *key, char *value);

extern char *search_allele(const char *allele);

extern void debug_dump_allele_dict();

extern void allele_prop_reset();

#endif
