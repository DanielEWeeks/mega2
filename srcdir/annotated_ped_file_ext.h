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

#ifndef ANNOTATED_PED_FILE_EXT_H
#define ANNOTATED_PED_FILE_EXT_H


extern int check_annotated_file_format(char *input_files[]);


extern void copy_exmap_locmap(linkage_locus_top *LTop,
			      ext_linkage_locus_top *EXLTop,
			      int map_num);

extern void Free_annotated_files(linkage_ped_top *Top);

extern int plink_args(char *str);

extern linkage_ped_top *read_annotated_files(char *ped_file, char *names_file,
					     char *map_file, char *freq_file,
					     char *pen_file, char *omit_file,
                         char *bed_file, char *phe_file,
					     int untyped_ped_opt,
					     analysis_type analysis,
                         plink_info_type *plink_info);  //temp

extern void write_annotated_aff(FILE *filep, int locusnm, linkage_locus_rec *locus,
				linkage_ped_rec *entry);


extern void write_annotated_numbered(FILE *filep, int locusnm,
				     linkage_ped_rec *entry);


extern void write_annotated_quant(FILE *filep, int locusnm,
				  linkage_ped_rec *entry);


#endif
