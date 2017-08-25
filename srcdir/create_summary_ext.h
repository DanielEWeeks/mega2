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

#ifndef CREATE_SUMMARY_EXT_H
#define CREATE_SUMMARY_EXT_H


extern int aff_status_entry(int status, int xclass, linkage_locus_rec *Locus);


extern void count_alleles(linkage_ped_top *Top1, allele_freq_struct *ar,
			  int loc, int inc_ht, int comp_geno);

extern void create_summary_file(linkage_ped_top * LPedTreeTop,
				analysis_type *analysis,
				char *fl_names[], int untyped_ped_opt,
				int *numchr);

extern void log_quant_selections(linkage_ped_top *Top);


extern void marker_typing_summary(FILE *sum_fp, linkage_ped_top *LPedTreeTop,
				  int numchr);


#endif
