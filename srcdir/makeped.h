/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

/* makedped.h - makeped data structures */
#ifndef MAKEPED_H
#define MAKEPED_H

#define PFOUNDER(p) ((p).father == 0 && (p).mother == 0)
/* reads in a pedigree file, checks if it is in pre-makeped format,
   then converts it to post-makeped format, breaking loops if necessary.
   Reference:
   Automatic Selection of Loop Breakers for Genetic Linkage Analysis,
   Ann Becker, Dan Geiger & Alejandro A. Schaffer,
   October 9, 1997
*/

/* define some data structures to store pedigrees */

typedef struct marriage_edge_ {
    int marriage_node1, marriage_node2;
    int edge_person, weight;
    int in_spanning_tree;
    int edge_id;

} marriage_edge_type;

/* this is the structure that will be used to get the
   min span tree for detecting and breaking loops.
   num_nodes = marriage_graph.num_marriages,
   i.e. one node for each marriage */

typedef struct minimizing_graph_ {
    int num_edges, num_nodes, num_edges_out;
    int *node_ids;
    int *node_selected;
    marriage_edge_type  *edges;

} minimizing_graph_type;

/* added field "linked" to tag whether siblings were linked
   while figuring out first_pa_sibs, first_ma_sibs when
   writing linkage_record */

typedef struct _marriage_node_{
    int node_id;
    int male_member, female_member;
    int num_offspring, *offspring;
    int linked;
} marriage_node_type;


/* A marriage pedigree has two sets of nodes: person nodes and marriage nodes.
   If only single marriages are allowed, then any person node has at most
   2 adjacent nodes, the marriage node that produced the person, and the
   marriage that person participates in.
*/
/* struct _person_node_ {
   node_id, indiv, father, mother, gender;
   linkage_pedrec_data *data;
   from_marriage_node_id, to_marriage_node_id;
   degree_genotyped, int removed, loop_breaker_id;

   } person_node_type;
*/

typedef struct _marriage_graph_ {
    int ped, max_id;
    int total_data_cols, num_selected_loci;
    int num_persons, num_marriages;
    person_node_type *persons;
    marriage_node_type *marriages;
    char Name[FILENAME_LENGTH];
} marriage_graph_type;

class   CLASS_ANALYSIS;
typedef CLASS_ANALYSIS *analysis_type;

typedef struct _linkage_ped_top {
    int PedCnt, IndivCnt;
    int OrigIds;
    int UniqueIds;
    int pedfile_type;
    linkage_ped_tree *Ped;          /* will be Ped[] */
    marriage_graph_type *PTop;     /* will be PTop[] */
    linkage_locus_top *LocusTop;
    ext_linkage_locus_top *EXLTop;
    analysis_type analysis;
} linkage_ped_top;


typedef union _entry {
    person_node_type *PEntry;
    linkage_ped_rec *LEntry;
} entry_type;

#endif
