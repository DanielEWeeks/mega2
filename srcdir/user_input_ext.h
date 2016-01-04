/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

#ifndef USER_INPUT_EXT_H
#define USER_INPUT_EXT_H

#include "input_ops.hh"

extern INPUT_FORMAT_t Input_Format;
extern const char *INPUT_FORMAT_STR[];
extern Input_Base *Input;

extern double Imputed_Info_Metric_Threshold;

extern int ReOrderMenu(int num_chromo, int *chromsomes, int *selection);

extern void  affected_by_status(linkage_ped_top *Top, int locus);

extern int analysis_menu1(analysis_type  *analysis);


extern int check_quant_phenotype_data_has_value(linkage_ped_top *Top, double value);

extern void define_affection_labels(linkage_ped_top *Top,
				    analysis_type analysis);

extern void enter_number(int *num);

extern int individual_id_item(int item_number, analysis_type analysis,
                              int current_opt_val, int align_right_pos,
                              int get_or_disp, int has_orig, int has_uniq);

extern int invalid_analysis(int opt);


extern void menu1(file_format *infl_type,
                            char **pedfl_name, char **locusfl_name,
                            char **mapfl_name, char **pmapfl_name,
                            char **input_path, char **omitfl_name,
                            char **freqfl_name, char **penfl_name,
                            char **bedfl_name, char **phefl_name,
                            int *Untyped_ped_opt, int *err_sim_opt,
                            char **mega2_output_path, double *freq_miscmatch_thresh);


extern void ped_ind_defaults(int unique, analysis_type analysis);

extern int pedigree_id_item(int item_number, analysis_type analysis,
                            int current_opt_val, int align_right_pos,
                            int get_disp_log, int has_orig);

extern void set_missing_quant_input(linkage_ped_top *Top, analysis_type analysis);
extern void set_missing_quant_output(linkage_ped_top *Top, analysis_type analysis);

extern void test_modified(int choice);

extern char *untyped_ped_messg(int opt, char *msg);

extern void get_genetic_distance_index(ext_linkage_locus_top *EXLTop);
extern void get_base_pair_position_index(ext_linkage_locus_top *EXLTop);

extern int genetic_distance_index;
//cpk extern genetic_distance_map_type genetic_distance_sex_type_map; // see common.h
extern int genetic_distance_sex_type_map; // see common.h
extern int base_pair_position_index;

/*

extern int  gh_cov_selection(int num_select,
			     int *trait_order, int *local_covariates,
			     linkage_locus_top *LTop);

extern void  seg_affected_by_status(linkage_ped_top *Top, int locus,
				    char *SegFileName,
				    char **input_files);

*/

#endif
