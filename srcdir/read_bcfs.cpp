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


/*
 *  read lines reseting stringstream/vector each line
 */

#include <stdlib.h>

#include <iostream>
#include <fstream>

#include <string>
#include <vector>
#include <set>
#include <algorithm>

#include "common.h"
#include "tod.hh"
#include "error_messages_ext.h"
extern void           Exit(int arg, const char *file, const int line, const char *err);
#include "typedefs.h"
#include "fcmap_ext.h"
#include "batch_input.h"
#include "batch_input_ext.h"
#include "mrecode.h"
#include "mrecode_ext.h"
#include "reorder_loci_ext.h"
#include "plink_ext.h"
#include "annotated_ped_file_ext.h"
#include "annotated_ped_file.h"
#include "read_files_ext.h"
#include "vcftools/mega2_vcftools_interface.h"
#include "utils_ext.h"

#include "str_utils.hh"
#include "input.hh"
#include "input_ops.hh"

#include "read_bcfs.hh"

#ifdef _WIN
#define R_OK 4
#define access(str,type) _access(str,type)
#endif

void ReadBCFs::do_menu_display(int &idx, int line_len, int choiceA[])
{
    printf("%2d) %-*s%s\n", idx, line_len,
           "BCF File Directory:",
           ((!strcmp(*bcfs_path, "."))?"[ Current directory ]" : *bcfs_path));
    choiceA[idx] = site_bcfs_dir_i;
    idx++;

    printf("%2d) %-*s%s\n", idx, line_len,
           "BCF File Template:", *bcfs_template);
    choiceA[idx] = site_bcfs_template_i;
    idx++;
}
int ReadBCFs::do_menu_parse(int choice_)
{
    int ret = 0;
    if(choice_ == site_bcfs_dir_i) {
        while (1) {
            draw_line();
            printf("Please enter BCF directory name > ");
            fcmap(stdin, "%s", *bcfs_path);
            newline;

            if (access(bcfs_path, F_OK)) {
                printf("WARNING: Could not find directory %s\n", *bcfs_path);
                continue;
            } else if (!is_dir(*bcfs_path)) {
                printf("WARNING: %s is not a directory.\n", *bcfs_path);
                printf("Please specify a new or valid directory.\n");
                strcpy(*bcfs_path, ".");
                continue;
            } else if (access(*bcfs_path, W_OK)) {
                printf("WARNING: %s is not a writable directory.\n", *bcfs_path);
                printf("Please specify a new or valid directory.\n");
                continue;
            }
            else
                break;
        }
        ret = 1;
    } else if(choice_ == site_bcfs_template_i) {
        while (1) {
            draw_line();
            printf("To enter a template please enter a value of the form:\n");
            printf("[data?.bcf]\nWhere the wildecard '?' will replace the CHR number for all chromosomes.\n");
            printf("Please enter BCF file template format > ");
            fcmap(stdin, "%s", *bcfs_template);
            newline;

            Vecs bcfsplit;
            split(bcfsplit, *bcfs_template, "?");

            if (bcfsplit.size() != 2) {
                printf("Please include one and only one ? in the template name\n");
                continue;
            } else
                break;
        }
        ret = 1;
    }
    return ret;
}