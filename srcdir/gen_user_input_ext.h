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

#ifndef GEN_USER_INPUT_EXT_H
#define GEN_USER_INPUT_EXT_H

extern void execute_gen(char *gen_executable_name, char *locus_name, const char *gen_outfile,
			double homo_o, double homo_e, int num_genos);

extern char *gen_input_file(int *geno_count, int num_alleles, char *locus_name,
                            int ninditer, int nmciter, int ntests);

extern void gen_outfile_name(int numchr, char *outfile1, char *outfile2);

extern void gen_run_options(int *niniditer, int *nmciter, int *ntests);


extern void init_gen_outfile(char **input_files, char *outfile_name, int count_option);


#endif
