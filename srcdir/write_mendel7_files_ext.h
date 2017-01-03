/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef WRITE_MENDEL7_FILES_EXT_H
#define WRITE_MENDEL7_FILES_EXT_H

extern void csv_mendel7_pen_file(char *fl_name,
                                 linkage_ped_top *Top,
                                 int sex_linked);


extern void csv_save_mendel_peds(char *file_name, linkage_ped_top *Top);


extern void csv_write_mendel_locus_file(char *file_name,
                                        linkage_locus_top *LTop,
                                        int sex_linked);


extern void csv_write_mendel_map(char *mapfile, linkage_locus_top *LTop,
                                 analysis_type analysis, int xlikned);
#endif


// http://www.genetics.ucla.edu/software/mendel

// Page 21...
// Locus names can be up to 16 characters long.
#define MENDEL7_MAX_LOCUS_NAME_LEN                    16
// Page 24...
// The definition file also defines fators and their permitted categories (or levels).
// These are listed after all genetic locus entries.
#define MENDEL7_MAX_FACTOR_NAME_LEN                   16
// Page 25...
// Variable names can be up to 16 characters long.
#define MENDEL7_MAX_VARIABLE_NAME_LEN                 16
