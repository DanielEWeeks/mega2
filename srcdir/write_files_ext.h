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

#ifndef WRITE_FILES_EXT_H
#define WRITE_FILES_EXT_H


extern void create_linkage_files(linkage_ped_top **Top, int *numchr,
				 char *fl_names[], int untyped_opt);


extern int save_linkage_peds(char *pedfl, linkage_ped_top *Top,
                             analysis_type analysis, int num_loci, int *loci);

extern void write_affection_data(FILE *filep, int locusnm,
                                 linkage_locus_rec *locus,
                                 linkage_ped_rec *entry);


extern void write_affection_tdtmax(FILE *filep, int locusnm,
                                   linkage_locus_rec *locus,
                                   linkage_ped_rec *entry);

extern void write_binary_data(FILE *filep, int locusnm,
                              linkage_locus_rec *locus,
                              linkage_ped_rec *entry);


extern void write_linkage_locfile_inorder(linkage_locus_top *LTop,
					  char *loutfl_name,
					  int locus_num,
					  int *loci, int sex_linked,
					  analysis_type analysis);


extern void write_locus_stats(linkage_locus_top *LTop,
			      file_format locus_file_type);

extern void write_numbered_data(FILE *filep,
                                const int locusnm,
                                linkage_locus_rec *locus,
                                linkage_ped_rec *entry);


extern void write_ped_stats(linkage_ped_top *Top, int pedfl_type);

extern int write_quant_stats(linkage_ped_top *Top, analysis_type analysis);



extern void write_quantitative_data(FILE *filep, int locusnm,
                                    linkage_locus_rec *locus,
                                    linkage_ped_rec *entry);


#endif
