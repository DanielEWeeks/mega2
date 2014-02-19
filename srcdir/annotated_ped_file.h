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

#ifndef ANNOTATED_PED_FILE_H
#define ANNOTATED_PED_FILE_H

/* number of reserved keywords for column headers */
#define ANN_UNKNOWN "NA"
#define NUM_PEDCOL_NAMES 16
#define NUM_REQ_PEDCOL_NAMES 5 /* the first 5 keywords have to be present */
#define NUM_MAPCOL_NAMES 3
#define NUM_FREQCOL_NAMES 3
#define NUM_PENCOL_NAMES 6
#define NUM_NAMESCOL_NAMES 2
#define NUM_OMITCOL_NAMES 3
#define COLNAMEWIDTH 50       /* maximum width of a column header */
#define MAX_GROUPS 20 /* How many population groups */

/* How many characters to read in at a time */
#define READ_CHUNK 400
#define INIT_COLNAME(colnames, i, type, name)   strcpy(colnames[i].ColName, name); \
    colnames[i].value_type=type;                                        \
    colnames[i].input_col = -1;                                         \
    colnames[i].locus_number = -1;                                      \
    colnames[i].map_number = -1;                                        \
    colnames[i].sex_map_number = -1;

typedef enum _col_hdr_err_type {
    CORRECT, UNKNOWN_HDR, NO_MAIN_PART, NO_EXT, LONG_EXT,
    UNKNOWN_EXT, MISSING_REQ_HDR
} col_hdr_err;

typedef enum _col_value_type {
    INT_AN, FLOAT_AN, CHAR_AN, STRING_AN, IGNORE
} col_value_type;


typedef struct _col_hdr_type_ {
    char ColName[COLNAMEWIDTH]; /* resreved keyword name, Ped, ID, Sex etc. */
    col_value_type value_type; /* what type of data */
    int input_col; /* which column in input file, set after parsing header */
    int output_col; /* which column to output after matching various files*/
    int locus_number; /* relevant only for pedigree file */
    int map_number, sex_map_number; /* relevant only for map file */
    int output_length;
    char output_format[10];
} col_hdr_type;

typedef struct _annotated_file_desc {
    int num_names_cols;
    col_hdr_type *names_file_columns;
    int num_freq_cols;
    col_hdr_type *freq_file_columns;
    int num_pen_cols;
    col_hdr_type *pen_file_columns;
    int num_ped_cols;
    col_hdr_type *ped_file_columns;
    int num_map_cols;
    col_hdr_type *map_file_columns;
    int num_omit_cols;
    col_hdr_type *omit_file_columns;
} annotated_file_desc;

typedef struct _annotated_ped_rec {
    int rec_num, ped_index, per_index;
    /* record number, pedigree number and membership in pedigree */
    char ID[MAX_NAMELEN];
    char Pedigree[MAX_NAMELEN], Father[MAX_NAMELEN], Mother[MAX_NAMELEN];
    char PedID[MAX_NAMELEN], PerID[MAX_NAMELEN];
    char FirstOff[MAX_NAMELEN], NextMatSib[MAX_NAMELEN], NextPatSib[MAX_NAMELEN];
    char Sex;
    int MZTwin, DZTwin, Proband, Group;
    int genocnt;
    int LinkPedID, LinkPerID;
    pheno_pedrec_data   *pheno;
    void                *marker;
} annotated_ped_rec;

extern col_hdr_type    ReservedColnames[NUM_PEDCOL_NAMES];

extern annotated_file_desc AnnotatedFileInfo;

#endif
