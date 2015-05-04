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

#ifndef BATCH_INPUT_EXT_H
#define BATCH_INPUT_EXT_H

#include "batch_input.h"

// called in mega2.cpp only...
extern void batchfile_init_Mega2BatchItems(void);
extern void batchfile_process(char *batch_file_name, analysis_type *analysis);
extern void check_batch_items(void);


// Called in many places were the batch item is to be updated
extern void batchf(int batch_item);
extern void batchf(const std::string& keyword);
extern void batchf(batch_item_type *bi);

extern void batch_file_doc(FILE *batchfp);
extern void create_batchfile(void);

extern void Free_batch_items(void);

#endif
