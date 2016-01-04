/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

/* NOTE:
 *
 * To avoid loops in the header files.  NO header file except typdefs.h should include headers that are
 * part of Mega2.  typedefs.h includes those headers that define types that are generally needed.  The heaers
 * are enumerated in a dependency ordering.
 *
 */

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#define USE_PROTOS prototypes
#define ANSI ansi

//In current gcc/glib you can not ignore the value returned by some functions.
//#define IgnoreValue(x) ((void)x)
//should be good enough
template<typename v>
void inline IgnoreValue(v x) {}

#ifdef ANSI
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#ifdef ANSI
#include <stdlib.h>
#endif

#if (!defined(CODEWARRIOR))  && (!defined(_WIN))
#include <unistd.h>
#endif


#ifdef DMALLOC
#include "dmalloc.h"
#endif

#include "errorno.h"
#include "allutil.h"
#include "linkage.h"
#include "linkage_external.h"
#include "makeped.h"
#include "pedtree.h"

extern analysis_type   AnalysisOpt; /* analysis option */

#include "batch_input.h"
#include "input_check.h"
#include "hwe.h"

#include "compress_ext.h"

#endif
