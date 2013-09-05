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

#ifndef MRECODE_EXT_H
#define MRECODE_EXT_H

extern void assign_dummy_alleles(marker_type *marker_list,
				 linkage_locus_top *LTop);

extern file_format check_locus_file_format(FILE *filep);

extern void convert_to_freq(marker_type *marker_list, linkage_locus_top *LTop,
			    int count_option, analysis_type analysis);

extern void count_classes(pheno_type *pheno_list, linkage_locus_top *LTop);

extern void count_raw_alleles(marker_type *markers, linkage_locus_top *LTop);

extern linkage_ped_top *create_allele_list(linkage_ped_top *Top,
					   int locus,
					   int count_option,
					   marker_type *marker_listi,
                                           allelecnt **member_ids,
					   int count_ht);

extern linkage_ped_top *create_full_marker_data(FILE *lfp,
						char *omit_flname,
						int *untyped_ped_opt,
						linkage_locus_top *LTop,
                                                int *col2locus,
						analysis_type analysis);

extern void default_trait_penetrances(linkage_ped_top *Top, int locus,
                                      pheno_type *pheno_listi);

extern void free_marker_item(allele_list_type *marker_item);

extern void free_marker_item_class(class_list_type *liability_class);

extern void load_classes(pheno_type *pheno, linkage_locus_top *LTop);

extern void need_recoding(marker_type *marker);

extern linkage_locus_top *read_marker_data(FILE *lfp, int type_col,
					   int name_col);

extern linkage_locus_top *read_marker_only_data(FILE *fp, int cols, char **phe_names, int *phe_types);

extern void recode_liability_class(marker_type *marker_list, linkage_locus_top *LTop, int mlocus);

extern void recode_locus_top(marker_type *marker_list,
                             pheno_type *pheno_list,
                             linkage_locus_top *LTop);

extern void recode_ped_top(marker_type *marker_list, linkage_ped_top *Top, plink_info_type *plink_info);

extern void write_recode_summary(marker_type *markers,
				 linkage_locus_top *LTop, int inc_ht);


#endif
