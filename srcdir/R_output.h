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


#ifndef R_OUTPUT_H
#define R_OUTPUT_H

#include "typedefs.h"

#define SW2OUT(chrstr, of)  sprintf(of, "STATS-%s.ALL", chrstr)
#define MERLINOUT(chrstr, of) sprintf(of, "merlin_out.%s", chrstr);
#define MERLINTABLE(chrstr, of) sprintf(of, "merlin_table.%s", chrstr);
#define MERLINPDF(chrstr, of) sprintf(of, "merlin_out.%s.pdf", chrstr);

extern char input_stem[4]; /* name for the generated  R_data file */
extern char output_stem[FILENAME_LENGTH];
extern char Rperl_file[20];
extern char Rdatafile[20];
extern char Rplot_function[14];
extern char Rshell_file[20];
extern char Rscript_file[20];
extern char Rmap_file[20];

typedef struct _out_pars {
    int pcomb_traits;
    int fcomb_traits;
    int fcomb_chromo;
    int ncols, nrows;
    double thresh;
    double Yrange[2];
    char Routputfile[FILENAME_LENGTH];
    char YLabel[15];
    char mode;
    int yfixed;
    double legend_cex, axis_cex, tcl;
    linkage_locus_top *LocusTop;
} R_plot_params_type;

typedef struct merlin_opt_ {
    int npl;
    int qtl;
    int pairs;
    int vc;
    int parametric;
    int markernames;
    int tabulate;
    int prefix;
    int pdf;
} merlin_opt_type;


#endif
/* ifndef R_OUTPUT_h */
