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

/* Changes to make when adding a new output option:

   This is where pre-makeped pedigrees are converted to postmakeped, and
   loops are broken. If your analysis option requires this conversion,
   add the KEYWORD to the second case statement in function
   makeped(), which deals with pre-makeped format pedigree files.

   If, on the other hand, your option requires broken loops to be
   maintained when post-makeped pedigrees are provided, then add the
   KEYWORD to the first case statement. Otherwise, post-makeped
   pedigrees will have their loops identified and reconnected.

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "makeped.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "read_files_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  Mega2LogF errorf mssgf my_calloc my_malloc my_realloc
              fcmap_ext.h:  fcmap
        grow_string_ext.h:  grow
            linkage_ext.h:  clear_lpedrec connect_loops free_all_from_lpedtop new_lpedtop
         read_files_ext.h:  read_numbered_data read_premakeped_affec read_premakeped_bin read_premakeped_quant
              utils_ext.h:  EXIT draw_line log_line
*/


/*--------------- static functions --------------*/
static int read_peds(int line_count, FILE *pfilep,
                     pre_makeped_record *persons,
                     linkage_locus_top *LTop,
                     int *col2locus);
static void free_min_graph(minimizing_graph_type *m_graph);

static void degree_genotyped(person_node_type *entry, linkage_locus_top *LTop);
static int make_marriage_graph(marriage_graph_type *mped);
static void prune_leaf_members(marriage_graph_type *mped, int pedc);
static int make_minimizing_graph(marriage_graph_type mped,
                                 minimizing_graph_type *m_graph);
static void min_span_tree(minimizing_graph_type *m_graph,
			  int no_founder,
			  marriage_graph_type *mped);
static int **find_split_nodes(minimizing_graph_type *m_graph);
static void copy_node_remove_parents(int proband, person_node_type *p1,
				     person_node_type *p2, int new_id,
				     int new_indiv, linkage_locus_top *LTop);
static void reassign_to_marriage(marriage_graph_type *m_graph,
				 person_node_type *p1, person_node_type *p2,
				 int to_marriage_node_id);
static void make_linkage_record(int pid, marriage_graph_type mped,
				linkage_locus_top *LTop,
				linkage_ped_rec *lrec);
static int check_and_reassign_parents(marriage_graph_type *mped);

static int check_pped_connected(minimizing_graph_type *mped,
				marriage_graph_type *ped);
static int check_disconnected_inds(marriage_graph_type *mped);

static int founder_more_typed(marriage_graph_type *m_graph,
			      marriage_edge_type *edges,
			      int nedges);
static int force_no_founders(void);

/*----------exported functions----------- */
int check_pre_makeped(FILE *fp, int *num_lines);
void free_marriage_graph(marriage_graph_type *m_graph);
linkage_ped_top *read_pre_makeped(FILE *fp, int pedcount,
				  int linecount,
				  pre_makeped_record *persons,
				  linkage_locus_top *LTop,
                                  int *col2locus);
int compare_edges(const void *e11, const void *e22);
void break_loops(int ped_count, marriage_graph_type *mped,
		 linkage_locus_top *LTop, int break_loop);
int makeped(linkage_ped_top *Top, analysis_type analysis);
void copy_marriage_graph(marriage_graph_type *From, marriage_graph_type *To);
void copy_marriage_graph_person(person_node_type *From,
				person_node_type *To,
				linkage_locus_top *LocusTop);
void copy_preped_peds(linkage_ped_top *From, linkage_ped_top *To);
void clear_prepedtree(marriage_graph_type *m_graph);
void count_pgenotypes(linkage_ped_top *Top, size_t *num_inds,
		      size_t *num_males, size_t *num_females,
		      size_t *untyped, size_t *typed, size_t *ped_typed,
		      size_t *males_typed, size_t *females_typed,
		      size_t *half_typed);
marriage_graph_type *make_ped_structure(pre_makeped_record *persons,
					int pedcount, int linecount);
int copy_pedrec_data(person_node_type *p1, person_node_type *p2,
                     linkage_locus_top *LTop);

/**************/

static int read_peds(int line_count, FILE *pfilep,
                     pre_makeped_record *persons,
                     linkage_locus_top *LTop,
                     int *col2locus)

{
    /* persons have already been allocated */
    /* person data structure contains
       ped, indiv, father, mother, gender, data.
       We already know the line_count */
    char rest[FILENAME_LENGTH], dummy1[199], dummy2[199];
    int i, lch = ' ', num_read, unique=0;
    int col;
    int expected_col=5+LTop->NumPedigreeCols, last_marker;
    SUPPRESS_MSSG_NESTED_INIT(untyped);
    int untyped = 0, totaltyped = 0;
    int a1, a2;
    const char *ar1, *ar2;

    for (i=0; i<line_count; i++) {
        num_read=fscanf(pfilep, "%d %d %d %d %d",
                        &(persons[i].ped), &(persons[i].indiv),
                        &(persons[i].father), &(persons[i].mother),
                        &(persons[i].gender));
        if (num_read <=0) {
            /* last blank line */
            break;
        }
        if (num_read != 5 && num_read > 0) {
            errorvf("File %s, Line %d: Malformed line.\n", mega2_input_files[0], i+1);
            EXIT(INPUT_DATA_ERROR);
        }
        col = 5;
        if ((persons[i].gender) && (persons[i].gender != MALE_ID) &&
            (persons[i].gender != FEMALE_ID)) {
            errorvf("File %s, Line %d, Col 5: Sex has invalid value %d.\n",
                    mega2_input_files[0], i+1, persons[i].gender);
            EXIT(DATA_INCONSISTENCY);
        }
        persons[i].rec_num = i+1;

        persons[i].pheno  = (LTop->PhenoCnt > 0) ? CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;
        // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
        persons[i].marker = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
        int locus;
        for (;col < expected_col; ) {
            if (col2locus[col] == 0) continue;
            locus = col2locus[col] - 1;
            last_marker = 0;
            lch=fgetc(pfilep);
            while(isspace(lch)) {
                if ((lch == LF ||  lch == CR) && !last_marker) {
                    errorvf("File %s, Unexpected end of line %d\n",
                            mega2_input_files[0], i+1);
                    errorvf("Expecting %d cols, read in %d.\n",
			    expected_col, num_read);
                    EXIT(INPUT_DATA_ERROR);
                }
                lch=fgetc(pfilep);
            }
            ungetc(lch, pfilep);
            if (i == 0) {
                LTop->Locus[locus].col_num = col;
            }
            switch (LTop->Locus[locus].Type) {
            case QUANT:
                lch=read_premakeped_quant(pfilep, locus, &LTop->Pheno[locus], &(persons[i]));
		// Here we are not checking for it being undefined, just "invalid".
		// There should be a better way...
                if (persons[i].pheno[locus].Quant <= QUNDEF) {
                    errorvf("File %s, Line %d, Col %d: Invalid quant data at locus %s.\n",
                            mega2_input_files[0], i+1, num_read+lch,
                            LTop->Locus[locus].Name);
                    EXIT(INPUT_DATA_ERROR);
                }
                col++;
                num_read += 1;
                break;
            case AFFECTION:
                lch=read_premakeped_affec(pfilep, locus, &(LTop->Pheno[locus]), &(persons[i]));
                if (persons[i].pheno[locus].Affection.Status == UNDEF ||
                    persons[i].pheno[locus].Affection.Class == UNDEF) {
                    errorvf("File %s, Line %d, Col %d: Invalid affection data at locus %s.\n",
                            mega2_input_files[0], i+1, num_read+lch,
                            LTop->Locus[locus].Name);
                    EXIT(INPUT_DATA_ERROR);
                }
                if (LTop->Pheno[locus].Props.Affection.ClassCnt > 1) {
                    num_read += 2;
                    col += 2;
                } else {
                    col++;
                    num_read++;
                }
                break;
            case BINARY:
                lch=read_premakeped_bin(pfilep, locus, &(LTop->Locus[locus]), &(persons[i]));
                get_2alleles(persons[i].marker, locus, &a1, &a2);
                if (a1 == UNDEF || a2 == UNDEF) {
                    errorvf("File %s, Line %d, Col %d: Invalid binary data at locus %s.\n",
                            mega2_input_files[0], i+1, num_read+lch,
                            LTop->Locus[locus].Name);
                    EXIT(INPUT_DATA_ERROR);
                }
                col += 2;
                num_read += 2;
                break;
            case NUMBERED:
            case XLINKED:
            case YLINKED:
                lch=read_numbered_data(pfilep, locus,
                                       &LTop->Locus[locus],
                                       (void *) (&(persons[i])),
                                       LTop->PedRecDataType, last_marker);
                switch(LTop->PedRecDataType) {
                    case Premakeped:
                        get_2alleles(persons[i].marker, locus, &a1, &a2);
                        if (a1 == UNDEF || a2 == UNDEF) {
                            errorvf("File %s, Line %d, Col %d: Invalid numbered data at locus %s.\n",
                                    mega2_input_files[0], i+1, num_read+lch,
                                    LTop->Locus[locus].Name);
                            EXIT(INPUT_DATA_ERROR);
                        }
                        break;
                    case Raw_premake:
                        get_2Ralleles(persons[i].marker, locus, &ar1, &ar2);
                        if (!strcmp(ar1, REC_UNDEF) ||
                            !strcmp(ar2, REC_UNDEF)) {
                            errorvf("File %s, Line %d, Col %d: Invalid numbered data at locus %s.\n",
                                    mega2_input_files[0], i+1, num_read+lch,
                                    LTop->Locus[locus].Name);
                            EXIT(INPUT_DATA_ERROR);
                        }
                        break;
                    default:
                        break;
                }
                col += 2;
                num_read += 2;
                break;
            default:
                errorf("Unknown type");
                EXIT(DATA_TYPE_ERROR);
            }
        }
        strcpy(dummy1, " ");
        strcpy(dummy2, " ");
        if (lch != LF && lch != CR) {
            (void)fgets(rest, FILENAME_LENGTH-1, pfilep);
            sscanf(rest, "%s %s", dummy1, dummy2);
            if (!strcasecmp(dummy1, "id:")) {
                unique=1;
                strcpy(persons[i].uniqueid, dummy2);
            } else {
                sprintf(persons[i].uniqueid, "%d", persons[i].indiv);
            }
        }
        totaltyped++;
        persons[i].genocnt  = crunch_notype(&persons[i].marker, LTop);
        if (persons[i].genocnt == 0) {
/*
            SUPPRESS_MSSG_NESTED(untyped);
            warnvf("Untyped person %4d linenum %d: person %d/%d, fa %d, ma %d\n",
                   untyped, i+1, persons[i].ped, persons[i].indiv, persons[i].father, persons[i].mother);
*/
            untyped++;
        }
    }
    SUPPRESS_MSSG_NESTED_FORCE(untyped);
    if (untyped > 0)
        warnvf("Individuals untyped: %d out of %d\n", untyped, totaltyped);
    SUPPRESS_MSSG_NESTED_FINI(untyped);

    count_Missing_Quant_consistency();

    return unique;
}

marriage_graph_type *make_ped_structure(pre_makeped_record *persons,
                                        int pedcount, int linecount)

{
    int i,j, k, *indcount;
    int curr_ped, last_line_read=0;
    pre_makeped_record *person;

    marriage_graph_type *pped;
    pped = CALLOC((size_t) pedcount, marriage_graph_type);

    indcount= CALLOC((size_t) pedcount, int);

    /* make a first pass to determine the number of individuals
       in each pedigree */
    curr_ped=persons[0].ped; i=0; j=0;
    while (i < linecount) {
        if (curr_ped != persons[i].ped) j++;
        indcount[j]++;
        curr_ped=persons[i].ped;
        i++;
    }
    /* Now store the pedigree data */
    for (i=0; i<pedcount; i++) {
        pped[i].num_marriages=0;
        pped[i].marriages=NULL;
        pped[i].ped=persons[last_line_read].ped;
        sprintf(pped[i].Name, "%d", pped[i].ped);
        pped[i].max_id=0;
        pped[i].num_persons = indcount[i];
        pped[i].persons = CALLOC((size_t) indcount[i], person_node_type);
        for (j=0; j < indcount[i]; j++) {
            person=&(persons[last_line_read + j]);
            strcpy(pped[i].persons[j].uniqueid, person->uniqueid);
            pped[i].persons[j].indiv=person->indiv;
            pped[i].persons[j].father=person->father;
            pped[i].persons[j].mother=person->mother;
            pped[i].persons[j].gender=person->gender;
            pped[i].persons[j].pheno=person->pheno;
            pped[i].persons[j].marker=person->marker;
            pped[i].persons[j].node_id=j;
            pped[i].persons[j].genocnt=person->genocnt;
            pped[i].persons[j].proband=((j==0)? 1 : 0);
            pped[i].persons[j].from_marriage_node_id=-1;
            pped[i].persons[j].num_to_marriages=0;
            pped[i].persons[j].first_pa_sib=
                pped[i].persons[j].first_ma_sib=
                pped[i].persons[j].first_off=-1;
            for (k=0; k<MAXMARRIAGES; k++)
                pped[i].persons[j].to_marriage_node_id[k] = -1;
            pped[i].persons[j].removed=0;
            pped[i].persons[j].loop_breaker_id=0;
            if (pped[i].max_id < pped[i].persons[j].indiv) {
                pped[i].max_id = pped[i].persons[j].indiv;
            }
        }
        last_line_read += indcount[i];
    }
    free(indcount);
    free(persons);
    return pped;
}

void free_marriage_graph(marriage_graph_type *m_graph)

{
    int i;

    if (m_graph->persons != NULL) {
        for (i=0; i < m_graph->num_persons; i++) {
            if (m_graph->persons[i].pheno != NULL)
                free(m_graph->persons[i].pheno);
            if (m_graph->persons[i].marker != NULL)
                free(m_graph->persons[i].marker);
        }
        free(m_graph->persons);
    }
    if (m_graph->marriages != NULL) {
        for (i=0; i < m_graph->num_marriages; i++) {
            if (m_graph->marriages[i].offspring != NULL)
                free(m_graph->marriages[i].offspring);
        }
        free(m_graph->marriages);
    }
}

static void free_min_graph(minimizing_graph_type *m_graph)

{

    free(m_graph->node_ids);
    free(m_graph->node_selected);
    free(m_graph->edges);
}

int copy_pedrec_data(person_node_type *p1, person_node_type *p2,
		     linkage_locus_top *LTop)

{
    int i, status=0;
    pheno_pedrec_data   *t1 = p1->pheno,  *t2 = p2->pheno;
    void                *s1 = p1->marker, *s2 = p2->marker;

    for (i=0; i < LTop->PhenoCnt; i++) {
        switch (LTop->Locus[i].Type) {
        case QUANT:
            t2[i].Quant = t1[i].Quant;
            break;
        case AFFECTION:
            status=t1[i].Affection.Status;
            t2[i].Affection.Status = t1[i].Affection.Status;
            t2[i].Affection.Class  = t1[i].Affection.Class;
            break;
        default:
            break;
        }
    }

    if (s1 == NOTYPED_ALLELES) {
        p2->marker = s1;
        return status;
    }

    for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
        switch (LTop->Locus[i].Type) {
        case BINARY:
        case NUMBERED:
        case XLINKED:
        case YLINKED:
            switch(LTop->PedRecDataType) {
            case Raw_premake:
            case Raw_postmake:
                copy_2Ralleles(s2, s1, i);
                break;
            case Premakeped:
            case Postmakeped:
                copy_2alleles(s2, s1, i);
                break;
            default:
                break;
            }
            break;
        default:
            break;
        }
    }
    return status;
}

//
// This routine is used to determine if the pedigree file is in premakeped format
// by checking that there are no lines contianing invalid or unknown sex entries.
// Returns 0 if invalid.
// Returns the number of pedigree records encountered if a valid premakeped formatted file.
//
int check_pre_makeped(FILE *fp, int *num_lines)
{
    int gender, current_ped=-1, pedc=0, c_prev=-100, c;
    int iped,  father, mother;
    size_t si;
    char ped[FILENAME_LENGTH], ind[FILENAME_LENGTH];
    char line[FILENAME_LENGTH];
    int num_valid_sex=0;
    int num_unknown_sex=0;
    int num_invalid_sex=0;
#ifndef HIDEFILE
    int first_unknown;
    int first_invalid = 0; // too hard for compiler
#endif
    *num_lines=0;

    // Here we gather only statistics about the pedigree file from the first five columns:
    // - count of distinct pedigrees
    // - count of valid/invalid/unknown sex entries
    // - on which line in the file the first invalid and unknown sex entry was found
    // NOTE: The actual data read from the file is thrown away...
    while (!feof(fp)) {
        strcpy(line, "");
        (void)fgets(line, FILENAME_LENGTH - 1, fp);
        if (*line == 0) break;
        strcpy(ped, ""); strcpy(ind, "");
        father= 0;
        mother = 0;
        gender = -1;
        /* This should be enough to read in the 1st 5 columns */
        c=sscanf(line, "%s %s %d %d %d", ped, ind, &father, &mother, &gender);
        if (c > 0) {
            if (c < 5) {
                /* This check will allow blank lines, where c <= 0,
                   necessary for trailing blank lines */
                errorvf("File %s, Line %d: Incomplete pedigree record.\n",
                        mega2_input_files[0], *num_lines + 1);
                EXIT(INPUT_DATA_ERROR);
            }
            if (c_prev <= 0 && c_prev != -100) {
                /* This means a blank line in the middle of the pedigree file */
                errorvf("File %s, Line %d: Blank line in pedigree file.\n",
                        mega2_input_files[0], *num_lines +1);
                EXIT(INPUT_DATA_ERROR);
            }
            /* first check to see if the ind and ped are integers */
            for (si=0; si  < strlen(ped); si++) {
                if (!isdigit((int)ped[si])) {
                    errorvf("File %s, Line %d: Non-numeric characters in Ped id (%s)\n",
                            mega2_input_files[0], *num_lines+1, ped);
                    errorf("This version of Mega2 requires integer ids.");
                    EXIT(INPUT_DATA_ERROR);
                }
            }
            for (si=0; si  < strlen(ind); si++) {
                if (!isdigit((int)ind[si])) {
                    errorvf("File %s, Line %d: Non-numeric characters in Person id (%s).\n",
                            mega2_input_files[0], *num_lines+1, ind);
                    errorf("This version of Mega2 requires integer ids.");
                    EXIT(INPUT_DATA_ERROR);
                }
            }

            // set counts and flags according to the current gender...
            if (gender == 0) {
                num_unknown_sex++;
                if (num_unknown_sex == 1) {
#ifndef HIDEFILE
                    first_unknown = *num_lines+1;
#endif
                }
            } else if (gender == 1 || gender == 2) {
                num_valid_sex++;
            } else {
                num_invalid_sex++;
                if (num_invalid_sex == 1) {
#ifndef HIDEFILE
                    first_invalid = *num_lines+1;
#endif
                }
            }

            iped = atoi(ped);
            // Set the count of pedigrees...
            if (current_ped != iped) {
                pedc++; current_ped=iped;
            }

            // Set the number of lines in the file
            (*num_lines)++;
        }

        c_prev = c;
        // Throw away the remainder of the line...
        if (strlen(line) >= 1) {
            /* skip to the end of line */
            if (line[strlen(line)-1] != '\n') {
                /* read until the end */
                fcmap(fp, "%=\n", c);
            }
        }
        if (feof(fp)) {
            break;
        }
    }

    /* Not a premakeped format */
    if (num_unknown_sex > 0 || num_invalid_sex > 0) {
#ifndef HIDEFILE
        mssgvf("Pedigree file %s will be read in as post-makeped format.\n",
               mega2_input_files[0]);
        mssgf("Reason:");
        if (num_unknown_sex > 0) {
            mssgf("Observed 0s in gender column (5th col).");
            mssgvf("First occurrence is in line %d.\n", first_unknown+1);
        }
        if (num_invalid_sex > 0) {
            mssgf("Observed values > 2 in gender column (5th col).");
            mssgvf("First occurrence is in line %d.\n", first_invalid+1);
        }
        log_line(mssgf);
#endif
        return 0;

    } else {
#ifndef HIDEFILE
        mssgvf("Pedigree file %s is in pre-makeped format.\n", mega2_input_files[0]);
#endif
        draw_line();
        return pedc;
    }
}

/*      if (strlen(line) >= (filename_length - 2)) { */
/*        if ((c=fgetc(fp)) != '\n') { */
/*  	fputc(c, fp); */
/*        } */


linkage_ped_top *read_pre_makeped(FILE *fp, int pedcount,
				  int linecount,
				  pre_makeped_record *persons,
				  linkage_locus_top *LTop,
                                  int *col2locus)

{
    /*
      Read in a pre_makeped file and assign ped, individual ids, father and mother
      This function should be invoked after the locus file is read in, and the loci
      selected. Note: have to separate the selection process from the locus_reordering.

      Check if pre-makeped format
      5th column will have at least one 0 if post-makeped (link to first sibling)
      will have only 1s and 2s if pre-makeped (gender)

      for only numbered loci
    */

    int i, unique;
    marriage_graph_type *pped;
    linkage_ped_top *Top;

    if (persons == NULL) {
        if ((persons=CALLOC((size_t) linecount, pre_makeped_record))==NULL) {
            errorf("Could not allocate enough memory, exiting.");
            EXIT(MEMORY_ALLOC_ERROR);
        }
        // NOTE: This never sets LTop->Locus[i].Allele[j].name for entry.cpp::pr_marker_alleles()
        unique=read_peds(linecount, fp, persons, LTop, col2locus);
        if (!unique) {
            for (i=0; i<linecount; i++) {
                sprintf(persons[i].uniqueid, "%d_%d",
                        persons[i].ped, persons[i].indiv);
            }
        }
    } else {
        unique=1;
    }
    pped=make_ped_structure(persons, pedcount, linecount);

    Top=new_lpedtop();
    Top->LocusTop = LTop;

    Top->PedCnt=pedcount;
    Top->pedfile_type=PREMAKEPED_PFT;
/*  Top->UniqueIds=unique; */
    Top->UniqueIds=1;
    Top->OrigIds=0;         /* will be set = 1 iff makeped is run */
    Top->PTop=pped;
#ifndef HIDESTATUS
    msgvf("Data read in from pedigree file:\n");
#endif
    /*  log_line(mssgf); */
    return Top;
}

static void degree_genotyped(person_node_type *entry, linkage_locus_top *LTop)

{

    int m, allele1, allele2, deg, j;
    /* now compute the degree_genotyped based on selected loci */
    /* This is after reordering, so both the loci list and a person's
       phenotype has only the selected loci in the proper order */

    deg=0;
    /*  for (j=0; j < LTop->LocusCnt; j++) { */
    for (j=0; j < num_reordered; j++) {
        m = reordered_marker_loci[j];
        if (LTop->Locus[m].Type == NUMBERED ||
            LTop->Locus[m].Type == BINARY ||
            LTop->Locus[m].Type == XLINKED ||
            LTop->Locus[m].Type == YLINKED) {
            get_2alleles(entry->marker, m, &allele1, &allele2);
            if (allele1 > 0 && allele2 > 0) deg++;
        }
        entry->degree_genotyped=deg;
    }
}

static int make_marriage_graph(marriage_graph_type *mped)

{
    /* for a single pedigree */
    /* assume that all the person nodes have been filled in */

    int found, i, j, off_ind, fm_id;
    int current_marriage=0, num_marriages=0;
    int father=0, mother=0, *offspring_c;
    int ped_disconnected=0, *connected;

    /* first find out the number of marriages */
    if (mped->num_persons == 1) {
        mped->num_marriages = 0;
        mped->marriages = NULL;
        return 0;
    }
    connected = CALLOC((size_t) mped->num_persons, int);

    for (i=0; i<mped->num_persons; i++) {
        connected[i]= 0;
        if (mped->persons[i].father > 0) {
            connected[i]= 1;
            if (mped->persons[i].father != father ||
                mped->persons[i].mother != mother) {
                num_marriages++;
                father=mped->persons[i].father;
                mother=mped->persons[i].mother;
            }
        }
    }
    mped->num_marriages=num_marriages;
    mped->marriages=CALLOC((size_t) num_marriages, marriage_node_type);
    offspring_c=CALLOC((size_t) num_marriages, int);

    /* now find out the number of offspring for each marriage */
    father=0; mother=0;
    num_marriages=0; current_marriage=0;
    for (i=0; i<mped->num_persons; i++) {
        if (mped->persons[i].father > 0) {
            if (mped->persons[i].father != father ||
                mped->persons[i].mother != mother) {
                /* look for the marriage this belongs to */
                found=0;
                for (j=0; j < num_marriages; j++) {
                    if ((mped->marriages[j].male_member == mped->persons[i].father) &&
                        (mped->marriages[j].female_member == mped->persons[i].mother)) {
                        /* marriage found */
                        found=1; current_marriage = j;
                        break;
                    }
                }
                if (!found)  {
                    current_marriage = num_marriages;
                    mped->marriages[current_marriage].node_id=current_marriage;
                    mped->marriages[current_marriage].male_member=
                        mped->persons[i].father;
                    mped->marriages[current_marriage].female_member =
                        mped->persons[i].mother;
                    num_marriages++;
                }
            }
            mped->persons[i].from_marriage_node_id=current_marriage;
            offspring_c[current_marriage]++;
            father = mped->persons[i].father;
            mother = mped->persons[i].mother;
        }
    }

    mped->num_marriages = num_marriages;
    /* now fill out the marriages */
    for (i=0; i<num_marriages; i++) {
        if (offspring_c[i] > 0) {
            mped->marriages[i].offspring=CALLOC((size_t) offspring_c[i], int);
        }
        mped->marriages[i].num_offspring=0;
    }

    /* fill out the offspring */
    for (i=0; i< mped->num_persons; i++) {
        if ((fm_id = mped->persons[i].from_marriage_node_id) >= 0) {
            off_ind=mped->marriages[fm_id].num_offspring;
            mped->marriages[fm_id].offspring[off_ind]=i;
            mped->marriages[fm_id].num_offspring++;
        }
    }

    for (i=0; i < mped->num_persons; i++) {
        mother=0; father=0;
        for (j=0; j < num_marriages; j++) {
            if (mped->persons[i].indiv == mped->marriages[j].male_member) {
                connected[i] = 1;
                /* find all the different female spouses of  this male */
                if (mped->marriages[j].female_member != mother) {
                    mped->persons[i].to_marriage_node_id[mped->persons[i].num_to_marriages]=j;
                    mped->persons[i].num_to_marriages++;
                    mother=mped->marriages[j].female_member;
                }
            }
            if (mped->persons[i].indiv == mped->marriages[j].female_member) {
                connected[i] = 1;
                /* find all the different male spouses of this female */
                if (mped->marriages[j].male_member != father) {
                    mped->persons[i].to_marriage_node_id[mped->persons[i].num_to_marriages]=j;
                    mped->persons[i].num_to_marriages++;
                    father=mped->marriages[j].male_member;
                }
            }
            if (mped->persons[i].num_to_marriages >= MAXMARRIAGES) {
                errorvf("Too many multiple marriages in pedigree %d (only %d allowed).\n",
                        mped->ped, MAXMARRIAGES);
                EXIT(OUTOF_BOUNDS_ERROR);
            }
        }
    }

    for (i=0; i < mped->num_persons; i++) {
        if (connected[i] == 0) {
            errorvf(" Pedigree %s has disconnected individuals: %s\n",
                    mped->Name, mped->persons[i].uniqueid);
            ped_disconnected=1;
//          break;
        }
    }
    free(connected);
    free(offspring_c);
    return ped_disconnected;
}


static void prune_leaf_members(marriage_graph_type *mped, int pedc)
{
    /* take out all the persons (set removed to 1), who do not participate
       in a marriage, or are not the offspring of a marriage */
    int i, j;
    marriage_graph_type curr_ped;

    for (i=0; i< pedc; i++) {
        curr_ped=mped[i];
        if (curr_ped.num_persons == 0) {
            curr_ped.persons[0].removed = 1;
        } else {
            for (j=0; j < curr_ped.num_persons; j++) {
                if (curr_ped.persons[j].num_to_marriages == 0)
                    curr_ped.persons[j].removed = 1;
                /* keep any person who participates in multiple marriages */
                if ((curr_ped.persons[j].from_marriage_node_id < 0) &&
                    (curr_ped.persons[j].num_to_marriages <= 1))
                    curr_ped.persons[j].removed = 1;
            }
        }
    }
}

/* Create a graph whose edges connect two marriage nodes through a person
   who has not been removed in prune_leaf_members for a single pedigree */

static int make_minimizing_graph(marriage_graph_type mped,
                                 minimizing_graph_type *m_graph)

{
    int i, j, k;
    int max_edges = mped.num_marriages *(mped.num_marriages + 1)/2;

    /* number of marriage edges are equal to the number of persons that have not
       been removed from the pedigree graph (N choose 2) where
       total marriage pairs for each person with a multiple marriage
    */
    /* check if this is a nuclear family, then no loops */
    if (mped.num_marriages <= 1) {
        m_graph->num_edges = 0;
        return 0;
    }

    /* at most as many edges as (marriages*marriages-1)/2 */
    m_graph->edges=CALLOC((size_t) max_edges, marriage_edge_type);
    m_graph->node_ids=CALLOC((size_t) mped.num_marriages, int);
    m_graph->node_selected=CALLOC((size_t) mped.num_marriages, int);
    m_graph->num_nodes=mped.num_marriages;
    m_graph->num_edges_out=0;
    for (i=0; i < mped.num_marriages; i++) {
        m_graph->node_ids[i]=mped.marriages[i].node_id;
        m_graph->node_selected[i]=0;
    }

    j=0;
    for (i=0; i < mped.num_persons; i++) {
        if (mped.persons[i].removed == 0) {
            if (mped.persons[i].from_marriage_node_id >=0) {
                /* create edges from "from_marriage" to each "to_marriage" */
                for (k=0; k < mped.persons[i].num_to_marriages; k++) {
                    m_graph->edges[j].marriage_node1=
                        mped.persons[i].from_marriage_node_id;
                    m_graph->edges[j].marriage_node2=
                        mped.persons[i].to_marriage_node_id[k];
                    m_graph->edges[j].edge_person=mped.persons[i].node_id;
                    m_graph->edges[j].weight=mped.persons[i].degree_genotyped;
                    m_graph->edges[j].edge_id=j;
                    j++;
                }
            } else {
                /* if a person has multiple marriages, create n-1 edges
                   between n to_marriages */
                if (mped.persons[i].num_to_marriages > 1) {
                    for (k=0; k < mped.persons[i].num_to_marriages-1; k++) {
                        m_graph->edges[j].marriage_node1 =
                            mped.persons[i].to_marriage_node_id[k];
                        m_graph->edges[j].marriage_node2 =
                            mped.persons[i].to_marriage_node_id[k+1];
                        m_graph->edges[j].edge_person=mped.persons[i].node_id;
                        m_graph->edges[j].weight=mped.persons[i].degree_genotyped;
                        m_graph->edges[j].edge_id=j;
                        j++;
                    }
                }
            }
        }
    }
    m_graph->num_edges = j;
    /* check if this pedigree is connected */
    if (check_pped_connected(m_graph, &mped) > 0) {
        return 1;
    }

    /* find out if there is a loop straight away, by comparing
       number of edges with number of marriage_nodes */
    if (m_graph->num_edges < mped.num_marriages) {
        free_min_graph(m_graph);
        m_graph->num_edges=0;
    }

    return 0;
}


int compare_edges(const void *e11, const void *e22)
{
    const marriage_edge_type *e1, *e2;

    e1 = (const marriage_edge_type *) e11;
    e2 = (const marriage_edge_type *) e22;

    if (e1->weight > e2->weight)
        return 1;
    else {
        if (e1->weight < e2->weight) return -1;

        if (e1->edge_person > e2->edge_person)
            return 1;
        else if (e1->edge_person < e2->edge_person)
            return -1;
        else return 0;

    }
}

/* create a minimum spanning tree of the marriage_graph-
   store the number of edges not chosen.
   1. Sort the edges
   2. Initialize the selected edge set with the very first edge
   3. Do the following until all nodes have been marked selected
   a. Select the first edge in the sorted set, such that
   one of its nodes is already selected, and the other
   one is not.
   b. Mark both nodes selected, and the edge selected
*/

static void min_span_tree(minimizing_graph_type *m_graph,
			  int no_founder,
			  marriage_graph_type *mped)

{
    int i, is_more_typed, all_selected=0;
    int num_in, edge_node1, edge_node2;
    int (*comp_edges)(const void *, const void *);

    comp_edges=compare_edges;
    /* assume that thi will only be invoked on cyclic graphs,
       number of nodes is less than or equal to the number of edges */

    /* have to sort the edges by weight in descending order */
    qsort(((void *) m_graph->edges), (size_t) m_graph->num_edges,
          sizeof(marriage_edge_type), comp_edges);

    is_more_typed=founder_more_typed(mped, m_graph->edges,
                                     m_graph->num_edges);
    if (is_more_typed && no_founder) {
        /* give non_founders a positive weight */
        for(i=0; i < m_graph->num_edges; i++) {
            if (!PFOUNDER(mped->persons[m_graph->edges[i].edge_person])) {
                m_graph->edges[i].weight += 999;
            }
        }
        /* then sort edges again */
        qsort(((void *) m_graph->edges), (size_t) m_graph->num_edges,
              sizeof(marriage_edge_type), comp_edges);
    }

    /* now make a spanning tree */
    num_in=0;
    /* start off with the first node of the first edge
       in the sorted set */
    edge_node1 = m_graph->edges[0].marriage_node1;
    m_graph->node_selected[edge_node1] = 1;
/*    for (i=0; i<m_graph->num_edges; i++) { */
/*      m_graph->edges[i].in_spanning_tree = 0; */
/*      printf("%d edge_node1 %d edge_node2 %d weight\n",  */
/*  	   m_graph->edges[i].marriage_node1, */
/*  	   m_graph->edges[i].marriage_node2, */
/*  	   m_graph->edges[i].weight); */

/*    } */
/*    sleep(3); */

    /* initialization */
    for (i=0; i < m_graph->num_edges; i++)
        m_graph->edges[i].in_spanning_tree = 0;

    while (!all_selected) {
        for (i=0; i < m_graph->num_edges; i++) {
            edge_node1 = m_graph->edges[i].marriage_node1;
            edge_node2 = m_graph->edges[i].marriage_node2;
            /* undirected graph, so consider both directions */
            if (m_graph->node_selected[edge_node1] == 1 &&
                m_graph->node_selected[edge_node2] == 0) {
                m_graph->edges[i].in_spanning_tree=1;
                m_graph->node_selected[edge_node2]=1;
                num_in++;
                break;
            } else {
                if (m_graph->node_selected[edge_node1] == 0 &&
                    m_graph->node_selected[edge_node2] == 1) {
                    m_graph->edges[i].in_spanning_tree=1;
                    m_graph->node_selected[edge_node1]=1;
                    num_in++;
                    break;
                }
            }
        }
        all_selected=1;
        for (i=0; i< m_graph->num_nodes; i++) {
            if (m_graph->node_selected[i] == 0) {
                all_selected=0; break;
            }
        }
        if (all_selected==1) break;
    }
    m_graph->num_edges_out=m_graph->num_edges-num_in;
}

static int **find_split_nodes(minimizing_graph_type *m_graph)

{
    int i;
    int **return_nodes, num_returned_nodes=0;

    return_nodes= CALLOC((size_t) m_graph->num_edges_out, int *);
    for (i=0; i<m_graph->num_edges_out; i++) {
        return_nodes[i] = CALLOC((size_t) 2, int);
    }
    num_returned_nodes=0;

    for (i=0; i<m_graph->num_edges; i++) {
        if (m_graph->edges[i].in_spanning_tree == 0) {
            return_nodes[num_returned_nodes][0]=m_graph->edges[i].edge_person;
            return_nodes[num_returned_nodes][1]=m_graph->edges[i].marriage_node2;
            num_returned_nodes++;
        }
    }

    return return_nodes;
}

static void copy_node_remove_parents(int proband, person_node_type *p1,
				     person_node_type *p2,
				     int new_id, int new_indiv,
				     linkage_locus_top *LTop)
{

    /* Makeped makes the assumption that if a loop person is also the proband,
     * then that person must be in the first loop.
     */
    if (p1->proband == 1) {
        p2->proband=2;
        if (proband !=2) {
            errorf("A loop person who is also a proband must be in the first loop (e.g., loop number 2).");
            errorvf("Loop person %d (unique ID %s) has proband number %d, but is in loop %d.\n",
		    p1->indiv, p1->uniqueid, p1->proband, proband);
            errorf("Edit your pedigree file to change who is the proband");
            errorf("or try a different option in the loop-breaker selection menu.");
            EXIT(LOOP_BREAKING_ERROR);
        }
    } else {
        p2->proband=proband;
        p1->proband=proband;
    }

    p2->node_id=new_id;
    strcpy(p2->uniqueid, p1->uniqueid);
    p2->indiv=new_indiv;
    p2->father=0; p2->mother=0;
    p2->from_marriage_node_id=-1;
    p2->degree_genotyped=p1->degree_genotyped;
    p2->pheno  = (LTop->PhenoCnt > 0) ? CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;
    // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
    p2->marker = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
    copy_pedrec_data(p1, p2, LTop);
    /*
    if (Mega2Status == INSIDE_RECODE) {
        copy_pedrec_data_raw_alleles(p1->data, p2->data, LTop);
    }
    */
    p2->gender=p1->gender;
    p2->loop_breaker_id=p1->indiv;
    p2->genocnt = p1->genocnt;
}

static void reassign_to_marriage(marriage_graph_type *m_graph,
				 person_node_type *p1, person_node_type *p2,
				 int m_id)
{
    /* p1->to_marriage should be made null
       p2->to_marriage should be set to p1->to_marriage
       p1->to_marriage.parents should be updated properly
       parents of the offspring of that marriage should be updated
       accordingly from p1 to p2*/

    int off, i,j;
    marriage_node_type *m;

    m=&(m_graph->marriages[m_id]);
    p2->num_to_marriages=1;
    p2->to_marriage_node_id[0]=m_id;
    if (m->male_member  == p1->node_id) {
        m->male_member = p2->node_id;
        for (i=0; i<m->num_offspring; i++) {
            /* offspring are indexes */
            off=m->offspring[i];
            m_graph->persons[off].father = m_graph->persons[p2->node_id].indiv;
        }
    } else {
        if (m->female_member == p1->node_id) {
            m->female_member = p2->node_id;
            for (i=0; i<m->num_offspring; i++) {
                off=m->offspring[i];
                m_graph->persons[off].mother = m_graph->persons[p2->node_id].indiv;
            }
        } else {
            errorvf("person %d is not in marriage %d (ped %d)!\n",
                    m_graph->persons[p1->node_id].indiv, m_id,
                    m_graph->ped);
            EXIT(INPUT_DATA_ERROR);
        }
    }
    for (i=0; i<p1->num_to_marriages; i++) {
        if (p1->to_marriage_node_id[i] == m_id) {
            for (j=i; j < p1->num_to_marriages - 1; j++)
                p1->to_marriage_node_id[j]=p1->to_marriage_node_id[j+1];
            p1->to_marriage_node_id[p1->num_to_marriages-1]=-1;
            break;
        }
    }
    for (i=p1->num_to_marriages; i<10; i++) {
        p2->to_marriage_node_id[i]=-1;
    }
    p1->num_to_marriages--;
}

void break_loops(int ped_count, marriage_graph_type *mped,
		 linkage_locus_top *LTop, int break_loop)

{
    /* called after the pre-makeped pedigree file has been read in,
       forms the marriage graph, prunes the graph, then finds loop-breakers
       by forming a minimal-spanning tree of the marriage graph
       (minimal because we are using the degree to which a person is
       genotypes, so we want loop-breakers, who are not in the spanning tree
       to have the highest degree of genotyping */

    int i, j, num_persons, no_founder=-1;
    int **ind_ids, proband;
    int ped_err, fatal_err=0;

    minimizing_graph_type *min_graph;
    person_node_type *p1, *p2;

    for (i=0; i<ped_count; i++) {
        for (j=0; j < mped->num_persons; j++)
            degree_genotyped(&(mped->persons[j]), LTop);
        make_marriage_graph(&(mped[i]));
        fatal_err += check_and_reassign_parents(&(mped[i]));
    }

    /* The disconnected pedigrees are a problem only if we need to
       break loops, so is this necessary? */

    if (break_loop==0) return;

    if (fatal_err) {
        errorvf("Found fatal errors in pedigree file, aborting Mega2.\n");
        EXIT(INPUT_DATA_ERROR);
    }

    prune_leaf_members(mped, ped_count);

    for (i=0; i<ped_count; i++) {
        min_graph= MALLOC(minimizing_graph_type);
        ped_err = make_minimizing_graph(mped[i], min_graph);
        if (break_loop == 0) {
            /* checking process for connected pefs frees the min_graph if check fails */
            continue;
        }
        if (ped_err) {
            errorvf("Ped %d: Unable to break loops.\n", mped[i].ped);
            fatal_err += ped_err;
            continue;
        }

        /* else go on to breaking loops */
        if (min_graph->num_edges > 0) {
            if (no_founder == -1)
                no_founder = force_no_founders();
            min_span_tree(min_graph, no_founder, &(mped[i]));
            ind_ids=find_split_nodes(min_graph);
#ifdef DEBUG_MAKEPED1
            for(j=0; j < min_graph->num_edges_out; j++) {
                if (PFOUNDER(mped[i].persons[ind_ids[j][0]]) && no_founder) {
                    printf("Tch tch : broke loop at a founder !!!\n");
                    break;
                }
            }
#endif
            num_persons=mped[i].num_persons;
            mped[i].persons=
                REALLOC(mped[i].persons,
                        (size_t) (num_persons+min_graph->num_edges_out), person_node_type);
            for (j=0; j < min_graph->num_edges_out; j++) {
                if (j==0) {
                    mssgvf("Found loops in pedigree %d.\n", mped[i].ped);
                }
                proband=j+2;
                p1=&(mped[i].persons[ind_ids[j][0]]);
                p2=&(mped[i].persons[num_persons+j]);
                mssgvf("Broke loops in Pedigree %d at Person %d (original id %d)\n",
                       mped[i].ped, ind_ids[j][0]+1, p1->indiv);
                copy_node_remove_parents(proband, p1, p2, (num_persons+j),
                                         (mped[i].max_id+1+j), LTop);

                reassign_to_marriage(&(mped[i]), p1, p2, ind_ids[j][1]);
            }
            mped[i].num_persons=num_persons+min_graph->num_edges_out;
            free_min_graph(min_graph);
        }
        free(min_graph);
    }

    if (fatal_err) {
        errorvf("Found fatal errors in pedigree file, aborting Mega2.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    return;
}

/* Set the first_maternal_sibs if the parent is the mom, or the
   first paternal sibs if the parent is the pa */

static void make_pointers(marriage_graph_type mped, int parent)

{
    int m;
    marriage_node_type mrec_to_next, mrec_to;
    int j;

    int parent_sex = mped.persons[parent].gender;

    for (m=0; m < mped.persons[parent].num_to_marriages; m++) {
        mrec_to=mped.marriages[mped.persons[parent].to_marriage_node_id[m]];
        if (m==0) {
            mped.persons[parent].first_off = mrec_to.offspring[0];
        }
        if (mrec_to.num_offspring > 1) {
            for (j =0; j < (mrec_to.num_offspring-1); j++) {
                if (parent_sex == MALE_ID) {
                    /* set the first_pa_sib */
                    mped.persons[mrec_to.offspring[j]].first_pa_sib = mrec_to.offspring[j+1];
                } else {
                    mped.persons[mrec_to.offspring[j]].first_ma_sib = mrec_to.offspring[j+1];
                }
            }
        }

        /* Now create the links between the last offspring of the current marriage, and
           first offspring of the next marriage */
        if (m < (mped.persons[parent].num_to_marriages-1)) {
            mrec_to_next = mped.marriages[mped.persons[parent].to_marriage_node_id[m+1]];
            if (parent_sex == MALE_ID) {
                /* set the first_pa_sib of the last offspring to the first offspring of
                   the next marriage. */
                mped.persons[mrec_to.offspring[mrec_to.num_offspring-1]].first_pa_sib =
                    mrec_to_next.offspring[0];
            } else {
                mped.persons[mrec_to.offspring[mrec_to.num_offspring-1]].first_ma_sib =
                    mrec_to_next.offspring[0];
            }
        }
    }
}

static void make_linkage_record(int pid, marriage_graph_type mped,
				linkage_locus_top *LTop,
				linkage_ped_rec *lrec)

{
    /* create a linkage (post-makeped) record from a pre-makeped record
       pointed to by the person pid. The sibling and offspring pointers are
       already set.
    */

    int first_pa_sib, first_ma_sib, first_off;
    int father=-1, mother=-1;

    person_node_type prec;
    int mrec_from /*, mrec_to*/ ;


    prec=mped.persons[pid];
    mrec_from=prec.from_marriage_node_id;
    /* mrec_to=prec.to_marriage_node_id[0]; */

    /* find out first maternal(paternal) sibling should be the offspring
       of the from_marriage occuring right after this person */

    first_pa_sib = prec.first_pa_sib;
    first_ma_sib = prec.first_ma_sib;

    /* printf("%d %d %d\n", prec.node_id + 1, first_pa_sib+1, first_ma_sib+1); */

    first_off = prec.first_off;
    if (mrec_from >= 0) {
        father = mped.marriages[mrec_from].male_member;
        mother = mped.marriages[mrec_from].female_member;
    } else {
        father = mother = -1;
    }

    /* NOTE: this breaks the "OLD" value of OrigID */
    /* now create the linkage_record */
    if (prec.loop_breaker_id > 0){
        sprintf(lrec->OrigID, "%d", prec.loop_breaker_id);
    } else {
        sprintf(lrec->OrigID, "%d", prec.indiv);
    }

    strcpy(lrec->UniqueID, prec.uniqueid);
    lrec->ID=pid+1;
    lrec->Sex=prec.gender;
    lrec->OrigProband=prec.proband;
    lrec->Ngeno=prec.degree_genotyped;
    lrec->loopbreakers = NULL;
    lrec->ext_ped_num=UNDEF;
    lrec->ext_per_num=UNDEF;

    /* fill in the links */
    if (father >= 0) lrec->Father=mped.persons[father].node_id+1;
    else lrec->Father=0;

    if (mother >= 0) lrec->Mother=mped.persons[mother].node_id+1;
    else lrec->Mother=0;

    if (first_off >= 0)
        lrec->First_Offspring=mped.persons[first_off].node_id+1;
    else
        lrec->First_Offspring=0;

    if (first_pa_sib >= 0)
        lrec->Next_PA_Sib=mped.persons[first_pa_sib].node_id+1;
    else
        lrec->Next_PA_Sib=0;

    if (first_ma_sib >= 0)
        lrec->Next_MA_Sib=mped.persons[first_ma_sib].node_id+1;
    else
        lrec->Next_MA_Sib=0;

    /* fill in the data */
/*    
    lrec->Data= CALLOC((size_t) LTop->LocusCnt, linkage_pedrec_data);
    lrec->Orig_status = copy_pedrec_data(prec.data, lrec->Data, LTop);
*/
    lrec->Pheno  = prec.pheno;
    lrec->Marker = prec.marker;
    lrec->genocnt   = prec.genocnt;
    prec.pheno  = NULL;
    prec.marker = NULL;
/*  but prec.data is a copy of ... so zero it too */
    mped.persons[pid].pheno = NULL;
    mped.persons[pid].marker = NULL;

    return;
}

#ifdef DEBUG_MAKEPED

/* the overall function to make a post-makeped file from a pre-makeped file:
   1) Check that it is a pre_makeped file
   2) Read in pedigree data into a marriage graph
   3) Break loops (if they exist)
   4) Go through the pedigrees one by one and output the linkage file
   5) Add original ped names and person ids at the end of each record
   returns 1, if a new file was created
*/

int makeped(char *pre_makeped_file, char *output_file, int num_loc,
	    int num_select, int *loci_select)

{
    int ped_count, i, j, k, num_lines;
    marriage_graph_type *big_graph;
    linkage_ped_rec lrec;

    FILE *fout;

    printf("%s\n", pre_makeped_file);
    ped_count=check_pre_makeped(pre_makeped_file, &num_lines);
    if (ped_count ==0) {
        printf("Post makeped format\n");
        return 0;
    } else {
        big_graph= CALLOC((size_t) ped_count, marriage_graph_type);
        for (i=0; i<ped_count; i++) {
            big_graph[i].num_persons=0;
            big_graph[i].ped=0;
            big_graph[i].total_data_cols=num_loc;
            big_graph[i].num_selected_loci=num_select;
            big_graph[i].num_marriages=0;

        }
        ped_count=read_pre_makeped(big_graph, pre_makeped_file, ped_count,
                                   num_loc, num_select, loci_select);

        break_loops(ped_count, big_graph);
        fout=fopen(output_file, "w");

        for (i=0; i< ped_count; i++) {
            for (j=0; j < big_graph[i].num_persons; j++) {
                lrec=make_linkage_record(i, j,
                                         big_graph[i], num_select, loci_select);
                fprintf(fout, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t",
                        lrec.ped_id, lrec.indiv, lrec.father, lrec.mother,
                        lrec.first_off, lrec.first_pat_sib, lrec.first_mat_sib,
                        lrec.gender);
                for (k=0; k < num_select; k++)
                    fprintf(fout, "%d  %d\t", lrec.data[2*k], lrec.data[2*k+1]);
                fprintf(fout, "Ped: %d\tPerson: %d\n",
                        big_graph[i].ped,
                        ((big_graph[i].persons[j].loop_breaker_id > 0)?
                         big_graph[i].persons[j].loop_breaker_id :
                         big_graph[i].persons[j].indiv));

            }
        }
        fclose(fout);
    }
    return 1;
}

int main(int argc, char *argv[])

{
    char infile[18];
    char outfile[18];
    int num_loc, num_select, i;
    int *loci_select;

    /* just to test out */

    printf("Input number of loci...");
    scanf("%d", &num_select); printf("\n");
    loci_select= CALLOC((size_t) num_select, int);
    for (i=0; i < num_select; i++) {
        printf("Input next locus > ");
        scanf("%d", &loci_select[i]); printf("\n");
    }
    printf("Input total number of loci > ");
    scanf("%d", &num_loc); printf("\n");
    makeped("test.in", "test.out",  num_loc, num_select, loci_select);
    exit(0);
}
#endif


static void check_any(linkage_ped_top *Top)
{
/*
    int i, j;

    for (i=0; i < Top->PedCnt; i++) {
        for (j=0; j <Top->Ped[i].EntryCnt; j++) {

            warnvf("genocnt: %s #%d\n",
                   Top->Ped[i].Entry[j].UniqueID,
                   Top->Ped[i].Entry[j].genocnt);
        }
    }
*/
}

int makeped(linkage_ped_top *Top, analysis_type analysis)

{
    /* int quiet; */
    int ped, i, j;
    linkage_ped_rec *Entry;
    linkage_ped_tree *LPed;

/*    FILE *fp; */
    /* Options not requiring broken loops:

       Mendel, Simwalk2, nuclear families,
       APMMult, APM, Aspex,
       solar, sage, genehunter,
       summary files, hwe

       Options requiring loops broken:

       Linkage, Homogeneity, TDT, SLINK, SPLINK,
       Vitesse, SIMULATE
    */

    if (Top->pedfile_type == POSTMAKEPED_PFT) {
	/*
	 * If your option requires broken loops to be
         * maintained when post-makeped pedigrees are provided, then add the
         * KEYWORD to the first case statement. Otherwise, post-makeped
         * pedigrees will have their loops identified and reconnected.
	 *
	 */
        if (analysis->maintain_broken_loops()) {
            Top->IndivCnt = 0;
            for (ped = 0; ped < Top->PedCnt; ped++) {
                LPed = &(Top->Ped[ped]);
                Top->IndivCnt += LPed->EntryCnt;
            }
        } else {
            /* else connect loops and count the total number of individuals */
            Top->IndivCnt = 0;
            for (ped = 0; ped < Top->PedCnt; ped++) {
                LPed = &(Top->Ped[ped]);
                if (LPed->Loops != NULL) {
                    if (connect_loops(LPed, Top) < 0)   {
                        errorf("Fatal error connecting loops. Aborting.");
                        EXIT(LOOP_CONNECTION_ERROR);
                    }
                }
                Top->IndivCnt += LPed->EntryCnt;
            }
        }
        check_any(Top);
        return 1;
    }

    /* else pre-makeped pedigree file */
    /*
     * This is where pre-makeped pedigrees are converted to post-makeped, and
     * loops are broken. If your analysis option requires this conversion,
     * add the KEYWORD to the second case statement in function
     * makeped(), which deals with pre-makeped format pedigree files.
     */
    if (analysis->break_loops())
        break_loops(Top->PedCnt, Top->PTop, Top->LocusTop, 1);
    else 
        break_loops(Top->PedCnt, Top->PTop, Top->LocusTop, 0);
    /* copy from one to the other */
    Top->Ped = CALLOC((size_t) Top->PedCnt, linkage_ped_tree);
    Top->IndivCnt = 0;

    for (i=0; i < Top->PedCnt; i++) {
        Top->Ped[i].EntryCnt=Top->PTop[i].num_persons;
        Top->Ped[i].Num=i+1;
        strcpy(Top->Ped[i].Name, Top->PTop[i].Name);
        Top->Ped[i].Loops=NULL;
        Top->Ped[i].Entry = CALLOC((size_t) Top->Ped[i].EntryCnt, linkage_ped_rec);
        /*    Top->Ped[i].Name = CALLOC(100, char); */
        /* first initialized the linked fields */
        for (j=0; j < Top->PTop[i].num_marriages; j++) {
            Top->PTop[i].marriages[j].linked=0;
        }
        for (j=0; j < Top->PTop[i].num_persons; j++) {
            Top->PTop[i].persons[j].first_off =
                Top->PTop[i].persons[j].first_pa_sib =
                Top->PTop[i].persons[j].first_ma_sib = -1;
        }
        for (j=0; j < Top->PTop[i].num_persons; j++) {
            make_pointers(Top->PTop[i], j);
        }
        for (j=0; j < Top->PTop[i].num_persons; j++) {
            if (Top->PTop[i].persons[j].proband == 1) {
                Top->Ped[i].Proband = Top->PTop[i].persons[j].node_id + 1;
            }
            Entry=&(Top->Ped[i].Entry[j]);
            clear_lpedrec(Entry);
            make_linkage_record(j, Top->PTop[i], Top->LocusTop, Entry);
        }

        Top->IndivCnt += Top->Ped[i].EntryCnt;
    }
    Top->OrigIds=1;
    printf("Created linkage ped tree\n");
    /*    fp=fopen("tmp_ped", "w"); */
    /*    save_linkage_peds(fp, Top, TO_MAKEPED, 0, NULL);  */
    /*    fclose(fp); */
    free_all_from_lpedtop(Top);
    Top->pedfile_type=POSTMAKEPED_PFT;

    /*  free_all_including_lpedtop(Top);
        fp=fopen("tmp_ped", "r");
        Top=NULL;
        quiet=1;
        Top=read_linkage_ped_file(fp, LTop);
        fclose(fp); */
    /*  system("/bin/rm -f tmp_ped"); */
    check_any(Top);
    return 1;
}

void clear_prepedtree(marriage_graph_type *m_graph)
{

    m_graph->num_persons=0;
    m_graph->num_marriages=0;
    m_graph->marriages=NULL;
    m_graph->persons=NULL;
    m_graph->ped=-1;
    m_graph->total_data_cols=0;
    m_graph->num_selected_loci=0;
}

void copy_marriage_graph_person(person_node_type *From,
				person_node_type *To,
				linkage_locus_top *LocusTop)

{
    To->father = From->father;
    To->mother = From->mother;
    To->gender = From->gender;
    To->indiv = From->indiv;
    To->node_id = From->node_id;
    To->proband = From->proband;
    To->from_marriage_node_id =  From->from_marriage_node_id;
    To->num_to_marriages =  From->num_to_marriages;

    To->pheno = From->pheno;
    From->pheno = NULL;
    To->marker = From->marker;
    From->marker = NULL;

    To->genocnt   = From->genocnt;
/*
    int k;

    for (k = 0; k < LocusTop->LocusCnt; k++)    {
        switch (LocusTop->Locus[k].Type)     {
        case NUMBERED:
            To->data[k].Alleles.Allele_1 = From->data[k].Alleles.Allele_1;
            To->data[k].Alleles.Allele_2 = From->data[k].Alleles.Allele_2;
            break;
        case AFFECTION:
            To->data[k].Affection.Status = From->data[k].Affection.Status;
            To->data[k].Affection.Class = From->data[k].Affection.Class;
            break;
        case QUANT:
            To->data[k].Quant = From->data[k].Quant;
            break;
        default:
            break;
        }
    }
*/
    return;
}

void  malloc_ppedtree_entries(marriage_graph_type *Copy, int entrycount,
			      int locuscnt)

{
    /* This will allocate enough space for a premake-ped
       marriage-graph, given the number of entries and information
       about locus data */
    Copy->persons = CALLOC((size_t) entrycount, person_node_type);
   /*
    int j;
    for (j = 0; j < entrycount; j++)
        Copy->persons[j].data = CALLOC((size_t) locuscnt, linkage_pedrec_data);
   */
    return;
}

void copy_marriage_graph(marriage_graph_type *From,
                         marriage_graph_type *To)

{
    /* pedigree-specific info */
    To->ped = From->ped;
    To->num_persons = From->num_persons;
    To->max_id = From->max_id;
}

/* code for copying pre-makeped pedigrees,
   malloc_linkage_ped_top has code for allocating this type
   of pedigree structure */


void copy_preped_peds(linkage_ped_top *From, linkage_ped_top *To)

{
    int i,j;

    if ((From != NULL) && (To != NULL)) {
        for (i = 0; i < From->PedCnt; i++)  {
            /* pedigree-specific info */
            copy_marriage_graph(&(From->PTop[i]),  &(To->PTop[i]));
            /* for each entry */
            for (j = 0; j < From->PTop[i].num_persons; j++)   {
                copy_marriage_graph_person(&(From->PTop[i].persons[j]),
                                           &(To->PTop[i].persons[j]), From->LocusTop);
            }
        }
    }
}

void count_pgenotypes(linkage_ped_top *Top, size_t *num_inds,
		      size_t *males, size_t *females, size_t *num_untyped,
		      size_t *num_typed, size_t *peds_typed, size_t *males_typed,
		      size_t *females_typed, size_t *half_typed)
{

    int i,j,k, ind_count=0, typed=0, halftyped=0, untyped=0;
    int male_count=0, female_count=0;
    int this_person_typed, this_ped_typed, this_male_typed, this_female_typed;
    int numloc, l, first_time=1;
    FILE *logfp;

    *peds_typed = *males_typed = *females_typed = 0;
    *half_typed = 0;

    strcpy(err_msg, "");
    // For each Pedigree...
    for(i=0; i < Top->PedCnt; i++) {
        this_ped_typed = 0;

        ind_count += Top->PTop[i].num_persons;
	// For each person (individual) in the pedigree...
        for (j=0; j < Top->PTop[i].num_persons; j++) {
	  // Count the number of males/females....
            person_node_type *pp = &(Top->PTop[i].persons[j]);
            int gender = pp->gender;
            if (gender == MALE_ID) {
                male_count++;
            } else if (gender == FEMALE_ID) {
                female_count++;
            }
            this_male_typed=0; this_female_typed=0;
            if (Mega2Status < LOCI_REORDERED) {
                numloc = Top->LocusTop->LocusCnt;
            } else {
	      numloc = num_reordered; // A global defined in common.h
            }
	    // For each loci associated with this person...
            for (l=0; l < numloc; l++) {
	      // Account for reordered loci...
                if (Mega2Status < LOCI_REORDERED) {
                    k=l;
                } else {
                    k = reordered_marker_loci[l];
                }
		// !(AFFECTION || QUANT) so it must be: TYPE_UNSET, BINARY, NUMBERED, XINKED, YLINKED
		// Determine if the alleles are known
                if (Top->LocusTop->Locus[k].Type != AFFECTION &&
                    Top->LocusTop->Locus[k].Type != QUANT)    {
                    this_person_typed=0;
                    if (Mega2Status <= INSIDE_RECODE && Top->LocusTop->PedRecDataType == Raw_premake) {
                        this_person_typed = num_typed_2Ralleles(pp->marker, k);
                    } else {
                        this_person_typed = num_typed_2alleles(pp->marker, k);
                    }

                    if (this_person_typed) {
                        this_ped_typed = 1;
                        if (gender == MALE_ID) {
                            this_male_typed=1;
                        } else if (gender == FEMALE_ID) {
                            this_female_typed=1;
                        }
                        if (this_person_typed == 2) {
			  // Both alleles are known...
                            typed++;
                        } else {
			  // Just one of the alleles is known...
                            halftyped++;
                        }
                    }
                    else
                        untyped++;
                }
            }
            if (gender == MALE_ID) {
                *males_typed += this_male_typed;
            } else if (gender == FEMALE_ID) {
                *females_typed += this_female_typed;
            }
        }

        if (! this_ped_typed) {
            /* output the pedigree number into log file */
            if (strlen(err_msg) > 70) {
                logfp = Mega2LogF();
                if (first_time) {
                    fprintf(logfp,
                            "------------------------------------------------------------\n");
                    fprintf(logfp, "Completely untyped pedigrees:\n");
                    first_time=0;
                }
                fprintf(logfp, "%s\n", err_msg);
                strcpy(err_msg, "");
                fflush(logfp);
            }
            strcat(err_msg, Top->PTop[i].Name); strcat(err_msg, " ");
        }
        *peds_typed += this_ped_typed;
    }

    if (strlen(err_msg) > 0) {
        logfp = Mega2LogF();
        if (first_time) {
            fprintf(logfp,
                    "------------------------------------------------------------\n");
            fprintf(logfp, "Completely untyped pedigrees:\n");
            first_time=0;
        }
        fprintf(logfp, "%s\n", err_msg);
        strcpy(err_msg, "");
        fflush(logfp);
    }

    *num_inds=ind_count;
    *males = male_count;
    *females = female_count;
    *num_untyped = untyped;
    *num_typed = typed;
    *half_typed = halftyped;
}

static int check_and_reassign_parents(marriage_graph_type *mped)
{
    int def_undef_par_err=0, no_par_err=0, par_sex_err=0;
    int father=0, mother=0 /*, first_off */;
    int i, j;

    // Check to see if one of the parents is missing...
    if (!AnalysisOpt->allow_missing_parent_in_linkage_input()) {
        for (i=0; i < mped->num_persons; i++) {
            if ((mped->persons[i].father == 0 && mped->persons[i].mother > 0) ||
                (mped->persons[i].father > 0 && mped->persons[i].mother == 0)) {
                errorvf("Ped %s: one parent of entry %s is defined, the other is undefined\n",
                        mped->Name, mped->persons[i].origid);
                def_undef_par_err=1;
            }
        }
    } else {
        warnf("Relaxing missing parent in Linkage input rule for this Analysis type.");
    }

    /* check that all persons defined in a pedigree structure actually
       exist, i.e. father, mother ids match up */
    for (i=0; i < mped->num_persons; i++) {
        father = mother = 0;
        // Check is only made if both parents exist...
        if (mped->persons[i].father > 0 && mped->persons[i].mother > 0) {
            for (j=0; j < mped->num_persons; j++) {
                if (mped->persons[j].indiv == mped->persons[i].father) father = 1;
                if (mped->persons[j].indiv == mped->persons[i].mother) mother = 1;
            }
            if (father == 0) {
                errorvf("Ped %d: father %d of entry %d is not in pedigree\n",
                        mped->ped,
                        mped->persons[i].father, mped->persons[i].indiv);
                no_par_err=1;
            }
            if (mother == 0) {
                errorvf("Ped %d: mother %d of entry %d is not in pedigree\n",
                        mped->ped,
                        mped->persons[i].mother, mped->persons[i].indiv);
                no_par_err=1;
            }
        }
    }

    /* change the father and mother ids of the marriages to indexes */
    for (i=0; i < mped->num_marriages; i++) {
        father=mped->marriages[i].male_member;
        mother=mped->marriages[i].female_member;
	/*first_off = mped->persons[mped->marriages[i].offspring[0]].indiv; */
        for (j=0; j < mped->num_persons; j++) {
            if (father == mped->persons[j].indiv) {
                /* check to see if father's sex is 1 */
/* 	if (mped->persons[j].gender != 1) { */
/* 	  sprintf(err_msg, "Ped %d: father %d of entry %d is not male", */
/* 			  mped->ped, father, first_off); */
/* 	  errorf(err_msg); */
/* 	  par_sex_err=1; */
/* 	} */
/* 	else { */
                mped->marriages[i].male_member = mped->persons[j].node_id;
/* 	} */
            }
            if (mother == mped->persons[j].indiv) {
/* 	if (mped->persons[j].gender != 2) { */
/* 	  sprintf(err_msg, "Ped %d: mother %d of entry %d is not female", */
/* 			  mped->ped, mother, first_off); */
/* 	  errorf(err_msg); */
/* 	  par_sex_err=1; */
/* 	} */
/* 	else { */
                mped->marriages[i].female_member = mped->persons[j].node_id;
/* 	} */
            }
        }
    }

    if (no_par_err || def_undef_par_err || par_sex_err) return 1;
    return 0;
}

static int check_disconnected_inds(marriage_graph_type *mped)

{
    int i, discon=0, *disconnected;

    if (mped->num_marriages < 1) {
        return 0;
    }
    disconnected = CALLOC((size_t) mped->num_persons, int);

    /* first do a simple check for disconnected individuals */
    for (i=0; i< mped->num_persons; i++) {
        if ((mped->persons[i].from_marriage_node_id < 0) &&
            (mped->persons[i].num_to_marriages < 1)) {
            /* person is not a parent or an offspring */
            disconnected[i]=1; discon++;
        }
    }
    if (discon > 1) {
        errorvf("Ped %d: Found %d disconnected individuals:\n",
                mped->ped, discon);
        for (i=0; i<mped->num_persons; i++) {
            if (disconnected[i]==1) {
                errorvf("Person %d\n", mped->persons[i].indiv);
            }
        }
    } else if (discon == 1) {
        for (i=0; i<mped->num_persons; i++) {
            if (disconnected[i]==1) {
                errorvf("Ped %d: Found 1 disconnected individual %d\n",
                        mped->ped, mped->persons[i].indiv);
                break;
            }
        }
        free(disconnected);
        return 1;
    }

    free(disconnected);
    return 0;
}



static int check_pped_connected(minimizing_graph_type *mped,
				marriage_graph_type *ped)
{

    int o, i,num_graphs=1;

    int *graph_num, *person_graph_num;
    int new_marked, all_marked=0;
    char names[100];


    /* now check to see that the pedigree does not contain disconnected
       subgraphs */
    if (check_disconnected_inds(ped)) {
        return 1;
    }

    graph_num = CALLOC((size_t) mped->num_nodes, int);
    person_graph_num = CALLOC((size_t) ped->num_persons, int);
    /*  roots=CALLOC(mped->num_nodes, int); */
    graph_num[mped->edges[0].marriage_node1]=num_graphs;
    /*  roots[0] = ped->marriages[mped->edges[0].marriage_node1].male_member; */
    new_marked=1;

    while (!all_marked) {
        while(new_marked) {
            /* go through the edges as many times as we can find a new node to mark
               mark as many marriage nodes as we can */
            new_marked=0;
            for (i=0; i < mped->num_edges; i++) {
                if ((graph_num[mped->edges[i].marriage_node1] > 0) &&
                    (graph_num[mped->edges[i].marriage_node2] == 0)) {
                    graph_num[mped->edges[i].marriage_node2] =
                        graph_num[mped->edges[i].marriage_node1];
                    new_marked=1;
                } else {
                    if ((graph_num[mped->edges[i].marriage_node2] > 0) &&
                        (graph_num[mped->edges[i].marriage_node1] == 0)) {
                        graph_num[mped->edges[i].marriage_node1] =
                            graph_num[mped->edges[i].marriage_node2];
                        new_marked=1;
                    }
                }
            }
        }
        /* now check if any marriages are not marked */
        all_marked=1;
        for (i=0; i<mped->num_nodes; i++) {
            if (graph_num[i] < 1) {
                all_marked=0;
                /*	roots[num_graphs]=ped->marriages[mped->node_ids[i]].male_member; */
                num_graphs++;
                graph_num[i]=num_graphs;
                break;
            }
        }
    }

    if (num_graphs > 1) {
        for (i=0; i<ped->num_marriages; i++) {
            /* graph num corresponds to ped->marriages
               Mark the members of each marriage_node with the graph number,
               members are ped->marriages[i].male_member,
               female_member, offspring etc.
            */
            person_graph_num[ped->marriages[i].male_member] = graph_num[i];
            person_graph_num[ped->marriages[i].female_member] = graph_num[i];
            for(o=0; o < ped->marriages[i].num_offspring; o++) {
                person_graph_num[ped->marriages[i].offspring[o]] = graph_num[i];
            }
        }

        errorvf("Ped %s: found %d disconnected sub-pedigrees:\n",
                ped->Name, num_graphs);

        for (i=0; i < num_graphs; i++) {
            sprintf(names, "        ");
            o=0;
            while(strlen(names) < 70 && o < ped->num_persons) {
                if (person_graph_num[o] == (i+1)) {
                    grow(names, "  %s",
                         ped->persons[o].uniqueid);
                }
                o++;
            }
            errorvf("Sub-pedigree %d including (but not limited to) individuals\n%s\n",
                    i+1, names);
        }
        free(graph_num); free(person_graph_num);
        return 1;
    }

    free(graph_num); free(person_graph_num);
    return 0;
}

static int founder_more_typed(marriage_graph_type *m_graph,
			      marriage_edge_type *edges,
			      int nedges)

{
    int nonfounder=0;
    int e;
    marriage_edge_type *edgep = &(edges[0]);

    for(e=0; e < nedges ; e++, edgep++) {
        if (nonfounder &&
            PFOUNDER(m_graph->persons[edgep->edge_person])) {
            return 1;
        } else {
            if (!PFOUNDER(m_graph->persons[edgep->edge_person])) {
                nonfounder=1;
            }
        }
    }
    return 0;

}

/* user input on forcing selection of non-founders only */

static int force_no_founders(void)

{
    char yesorno[4];
    int item=-1, force;

    force=0;

    if (DEFAULT_OPTIONS) {
        mssgf("Select ANY individual as loop-breaker (DEFAULT)");
        mssgf("as requested in the batch file.");
        return(force);
    }

    draw_line();
    printf("Mega2 does not restrict loop-breakers ");
    printf("to be non-founders\nby default. ");
    printf("However, some linkage programs may not support\n");
    printf("loop-breakers who are founders.\n");

    while(item != 0) {
        draw_line();
        printf("Loop-breaker selection menu:\n");
        printf("0) Done with this menu - please proceed.\n");
        printf("%c1) Select only non-founders as loop-breakers\n",
               (force? '*': ' '));
        printf("%c2) Select any individual as a loop-breaker\n",
               (force? ' ': '*'));
        printf("Enter 1 or 2 to select, 0 to exit > ");
        fflush(stdout);
        (void)fgets(yesorno, 3, stdin); newline;
        sscanf(yesorno, "%d", &item);
        if (item==1) {
            force = 1;
        } else if (item == 2) {
            force = 0;
        }
    }

    return force;
}
