/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
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

*/

#ifndef WRITE_PSEQ_FILES_EXT_H
#define WRITE_PSEQ_FILES_EXT_H

#include "write_plink_ext.h"


class CLASS_PSEQ: public CLASS_PLINK {
public:
    CLASS_PSEQ() : CLASS_PLINK() {
        _name = "PSEQ";
	strcpy(file_name_stem, "pseq");
    }
   ~CLASS_PSEQ() {}

    bool has_sub_options()    { return false; }

/*  this is necessary because we are inheriting from CLASS_PLINK not CLASS_ANALYSIS */
    virtual const char* output_quant_default_value() { return "NA"; }

    virtual const char* output_affect_default_value() { return "0"; }
    virtual bool output_affect_must_be_numeric() { return true; }

    virtual void create_output_file(linkage_ped_top *LPedTreeTop,
			    analysis_type *analysis,
			    char *file_names[],
			    int untyped_ped_opt,
			    int *numchr,
                            linkage_ped_top **Top2);

    void save_pheno_file(linkage_ped_top *Top,
			 const int pwid,
			 const int fwid);
    void create_sh_file(linkage_ped_top *Top,
			char *file_names[],
			const int numchr);
};

extern CLASS_PSEQ               *TO_PSEQ;

#endif /* WRITE_PSEQ_FILES_EXT_H */
