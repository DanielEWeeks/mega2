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

#ifndef GENETIC_UTILS_EXT_H
#define GENETIC_UTILS_EXT_H

extern void         LogFileNames(void);

extern void           backup_file(char *flname1);

extern void         check_map_positions(linkage_ped_top *Top, int numchr,
                                        int check_or_set);

extern char         *chrom_num_to_name(const int chrom_num, char *chrom_name);

extern int             **create_genos(int num_alleles);

extern void           create_marker_sublist(linkage_ped_top *Top, loci_markers *loci_mrkrs,
                                            int num_markers);

extern void           draw_histogram(FILE *filep, int num_alleles,
				     allele_freq_struct *ar,
				     const char *messg,
				     const char *sex);

extern void         exclaim(void);

extern int             genoindx(int i, int j);

extern void         getRunDate(void);

extern double           haldane_theta(double x);

extern double           haldane_to_kosambi(double f);

extern double           haldane_x(double t);

extern void           init_file(char *flname);

extern double          kinship(ped_tree *Ped, ped_rec *a, ped_rec *b);

extern double           kosambi_theta(double t);

extern double           kosambi_to_haldane(double f);

extern const char      *loc_type_name(linkage_locus_type loc_type);


extern void           log_marker_selections(linkage_ped_top *Top,
					    int *entries,
					    int num_entries, char *asterisks,
					    void (*log_func)(const char *messg),
					    int disp_err);

extern void         makedir(char *dirname);

extern void         mega2_version_check(void) ;

extern void         move_file(char *flname1, char *flname2);

extern int          purge_unneeded_pedigrees(ped_top *Top, int *save);

extern float        safe_divide(int a, int b);

extern void         switch_map(linkage_locus_top *LTop, double *position, sex_map_types sex);

extern char        *write_time();


/*

extern old_new        *add_new_to_old(int oldped, int newped, int oldmem, int newmem, old_new *old_new_list);

extern void           copy_file(char *flname1, char *flname2);

extern void            help_sib_map(void);

extern void         log_only(const char *messg);

extern void         move_logs(char *new_location);

*/

#endif
