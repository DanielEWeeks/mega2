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

#ifndef COMPRESS_H
#define COMPRESS_H

#define MARKER_SCHEME_BITS   1
#define MARKER_SCHEME_BYTE   2
#define MARKER_SCHEME_PTR    3

extern void *NOTYPED_ALLELES;

extern int  marker_size(int size);

extern void *marker_alloc(size_t size, int offset);

extern void marker_free(void *marker, int offset);



extern void get_2Ralleles(void *mp, int marker, const char **all1, const char **all2);

extern void set_2Ralleles(void *mp, int marker, linkage_locus_rec *locus, const char *all1, const char *all2);

extern int crunch_Rnotype(void **p, linkage_locus_top *LTop);

extern void copy_2Ralleles(void *to, void *from, int marker);

extern void order_heterozygous_allele_raw(linkage_ped_top *Top);



extern void get_2alleles(void *mp, int marker, int *all1, int *all2);

extern void set_2alleles(void *mp, int marker, linkage_locus_rec *locus, int all1, int all2);

extern int  crunch_notype(void **p, linkage_locus_top *LTop);

extern void copy_2alleles(void *to, void *from, int marker);

extern void copy_2alleles(void *to, void *from, int tomarker, int frommarker);

extern void order_heterozygous_allele(linkage_ped_top *Top);

#endif /* COMPRESS_H */
