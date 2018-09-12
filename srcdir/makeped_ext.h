/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

#ifndef MAKEPED_EXT_H
#define MAKEPED_EXT_H


extern int check_pre_makeped(FILE *fp, int *num_lines);

extern int copy_pedrec_data(person_node_type *p1, person_node_type *p2,
			    linkage_locus_top *LTop);

extern void copy_preped_peds(linkage_ped_top *From, linkage_ped_top *To);

extern void count_pgenotypes(linkage_ped_top *Top, size_t *num_inds,
			     size_t *num_males, size_t *num_females,
			     size_t *untyped, size_t *typed, size_t *ped_typed,
			     size_t *males_typed, size_t *females_typed,
			     size_t *half_typed);

extern void free_marriage_graph(marriage_graph_type *m_graph);

extern int makeped(linkage_ped_top *Top, analysis_type analysis);

extern linkage_ped_top *read_pre_makeped(FILE *fp, int pedcount,
					 int linecnt,
					 pre_makeped_record *persons,
					 linkage_locus_top *LTop,
                                         int *col2locus);


extern int break_no_founders_menu(void);

/*
extern void break_loops(int ped_count, marriage_graph_type *mped,
			linkage_locus_top *LTop, int break_loop);

extern void clear_prepedtree(marriage_graph_type *m_graph);

extern int compare_edges(const void *e11, const void *e22);

extern void copy_marriage_graph(marriage_graph_type *From,
				marriage_graph_type *To);

extern void copy_marriage_graph_person(person_node_type *From,
                                       person_node_type *To,
                                       linkage_locus_top *LocusTop);

extern int copy_pedrec_data1(linkage_pedrec_data *d1, linkage_pedrec_data *d2,
			     linkage_locus_top *LTop, int locuscnt,
			     int *locus_inds);

extern marriage_graph_type *make_ped_structure(pre_makeped_record *persons,
					       int pedcount, int linecnt);

extern void  malloc_ppedtree_entries(marriage_graph_type *Copy,
				     int entrycount, int locuscnt);

*/

#endif
