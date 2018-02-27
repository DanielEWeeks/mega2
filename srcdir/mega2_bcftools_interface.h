/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2018 Robert Baron, Charles P. Kollar,
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
  Mega2 interface to BCFTOOLS.

*/

#ifndef MEGA2_BCFTOOLS_INTERFACE_H
#define MEGA2_BCFTOOLS_INTERFACE_H

#include <string>
#include <vector>
#include "lib/bcftools-1.6/filter.h"
#include <lib/bcftools-1.6/htslib-1.6/htslib/synced_bcf_reader.h>

extern "C" {
    #include "lib/bcftools-1.6/filter.h"
    #include "lib/bcftools-1.6/vcfview.h"
};

#ifndef INTERNAL_MEGA2_BCFTOOLS_INTERFACE

class MEGA2_BCFTOOLS_INTERFACE {

public:
    MEGA2_BCFTOOLS_INTERFACE( ) { }

    ~MEGA2_BCFTOOLS_INTERFACE() { }

    int mega2_main_vcfview(int argc, char *argv[]);

    args_t * get_args(int argc, char *argv[]);
    int test_args(int argc, char *argv[]);

};

extern MEGA2_BCFTOOLS_INTERFACE *mega2_bcftools_interface;


#endif /* INTERNAL_MEGA2_BCFTOOLS_INTERFACE */


#endif /* MEGA2_BCFTOOLS_INTERFACE_H */

