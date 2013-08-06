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

#ifndef LINKAGE_EXT_H
#define LINKAGE_EXT_H

extern void clean_reordered_markers(linkage_locus_top *LTop, analysis_type analysis);

extern void clear_laffclass(linkage_affection_class *Class);

extern void clear_lallelerec(linkage_allele_rec *Allele);

extern void clear_llocusdata(linkage_locus_data *Data,
			     linkage_locus_type Type);

extern void clear_llocusrec(linkage_locus_rec *Locus,
			    linkage_locus_type Type);

extern void clear_lpedrec(linkage_ped_rec *Entry);

extern void clear_lpedtree(linkage_ped_tree *Ped);

extern int connect_loops(linkage_ped_tree *Ped, linkage_ped_top *Top);

extern void copy_linkage_locus_rec(linkage_locus_rec *from, linkage_locus_rec *to);


extern void copy_linkage_ped_top(linkage_ped_top * Original,
				 linkage_ped_top * Copy, int malloc_it);

extern void delink_linkage_ped_top(linkage_ped_top * Copy);


extern void copy_lpedrec(linkage_ped_rec *From, linkage_ped_rec *To,
			 linkage_ped_top *Top);

extern void copy_lpedrec_raw_alleles(linkage_ped_rec *From,
				     linkage_ped_rec *To,
				     linkage_ped_top *top);

extern void copy_lpedtree1(linkage_ped_tree *From, linkage_ped_tree *To);

extern void count_lgenotypes(linkage_ped_top *Top, size_t *num_inds,
			     size_t *males, size_t *females, size_t *untyped,
			     size_t *typed, size_t *peds_typed, size_t *males_typed,
			     size_t *females_typed, size_t *half_typed);


extern void free_all_from_llocusrec(linkage_locus_rec *Locus);

extern void free_all_from_lpedtop(linkage_ped_top *PTop);

extern void free_all_from_lpedtree(linkage_ped_tree *Ped);

extern void free_all_including_lpedtop(linkage_ped_top **PTop);

extern linkage_loop_rec *get_llooprec(list **List, int Num);

extern void get_loci_on_chromosome(int chromosome);

extern void get_unmapped_loci(int append);

extern void markers_on_chromosome(int num_chr);

extern int is_typed_lentry(linkage_ped_rec Entry, linkage_locus_top *LTop,
			   int mode, double *percent_typed);

extern void  malloc_lpedtree_entries(linkage_ped_tree *Copy, size_t entrycount,
				     size_t locuscnt);


extern linkage_locus_top  *new_llocustop(void);

extern linkage_ped_rec    *new_lpedrec(void);

extern linkage_ped_top    *new_lpedtop(void);

extern linkage_ped_tree   *new_lpedtree(void);

extern void set_link_IDs(linkage_ped_top *Top);

/*
extern void clear_llocustop(linkage_locus_top *LTop);

extern void clear_llooprec(linkage_loop_rec *Loop);

extern void clear_lpedtop(linkage_ped_top *PTop);

extern void copy_linkage_locus_top(linkage_locus_top *from, linkage_locus_top *to);

extern void free_all_from_llocustop(linkage_locus_top *LTop);

extern void free_all_from_lpedrec(linkage_ped_rec *Entry);

extern void free_all_including_llocustop(linkage_locus_top *LTop);

extern void free_chr_loci(void);

extern void            free_linkage_ped_top(linkage_ped_top * Copy);

extern void free_pedrec_raw_alleles(void *Entry, linkage_locus_top *LTop);


extern void malloc_linkage_locus_top(linkage_ped_top *Original,
                                     linkage_ped_top *Copy,
                                     int alloc_alleles);

extern void malloc_linkage_locus_top1(linkage_ped_top *Original,
                                      linkage_ped_top *Copy,
                                      int alloc_alleles,
                                      size_t locuscnt,
				      int *locus_inds);

extern void            malloc_linkage_ped_top(linkage_ped_top * Original,
					      linkage_ped_top * Copy,
					      int alloc_alleles);

extern linkage_loop_rec   *new_llooprec(void);

*/

#endif
