/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"
#include "error_messages_ext.h"
#include "utils_ext.h"
#include "input_ops.hh"

#include "dblite.hh"
#include "dbbatch.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;
extern INPUT_FORMAT_t Input_Format;

static
const char *INPUT_FORMAT_DATA[] = {
     "Mega2 format data with header",
     "Linkage format data",
     "Linkage data with Mega2 names file",
     "PLINK binary PED format data",
     "PLINK PED format data",
     "BCF format data",
     "VCF compressed format data",
     "VCF format data",
     "IMPUTE2 GEN format data",
     "IMPUTE2 BGEN format data",
};

void File_table::stat() {
    int xx = sizeof (INPUT_FORMAT_DATA) / sizeof (char *);
    msgvf("This database was created from %s using the following files:\n",
          Input_Format >= xx ? "just" : INPUT_FORMAT_DATA[Input_Format]);

#ifndef HIDEPATH
    for (int i = 0; i < NUMBER_OF_MEGA2_INPUT_FILES; i++) {
        if (mega2_input_files[i])
            msgvf("\t%20s   %s\n",
                  mega2_input_file_type[i], mega2_input_files[i]);
    }
#endif
}

int File_table::db_getall() {
    int ret = select_stmt && select_stmt->abort();
    int k = 0;
    char *f = (char *)0, *t = (char *)0;

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(k, f, t);

            if (*f != 0) {
                if (mega2_input_files[k]) {
                    free(mega2_input_files[k]);
                }
                mega2_input_files[k] = strdup(f);
            } else
                mega2_input_files[k] = 0;
            strcpy(mega2_input_file_type[k], t);

        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return ret;
}

void dbbatch_file_export(linkage_ped_top *Top) {
    int i;

    MasterDB.begin();

    for (i = 0; i < NUMBER_OF_MEGA2_INPUT_FILES; i++) {
        file_table.insert(i, mega2_input_files[i], mega2_input_file_type[i]);
    }

    MasterDB.commit();
}


void dbbatch_file_import(linkage_ped_top *Top) {

    MasterDB.begin();

    file_table.db_getall();

    MasterDB.commit();
}
