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

#ifndef OUTPUT_FILE_NAMES_EXT_H
#define OUTPUT_FILE_NAMES_EXT_H

extern void           CHR_STR(int numchr, char chr_str[]);

extern void           change_output_chr(char *file_name, const int numchr);

extern void           create_mssg(analysis_type analysis);

extern void           default_outfile_names(const analysis_type  analysis,
                                            int *numchr, char *fl_names[],
                                            const char *logdir);

extern void           Free_output_paths(analysis_type analysis);

extern char           *file_status(char *fl_name, char *status);

extern void           print_outfile_mssg(void);


extern void           set_output_paths(analysis_type analysis,
				       linkage_ped_top *Top,
				       char *Mega2OutputPth);

extern void           shorten_path(char *old_path, char *new_path);


#endif
