/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

#include "compress_ext.h"
#include "linkage.h"

#include "dblite.hh"
#include "dbmisc.hh"
#include "dballele.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;

extern int allele_count;
extern int ALLELE_ARRAY;
extern allele_prop **Allele_Array;
extern allele_prop *canonical_allele_internal(const char *v);

extern int  MARKER_SCHEME;
extern Alleles_int *MARKER_SCHEME3_alleles;

void CanonicalAllele_table::db_getall(linkage_ped_top *Top) {
    int ret = select_stmt && select_stmt->abort();
    int k = 0;
    const char *v = (char *)0;

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {

            ret = select(k, v);
            allele_prop *Ara = canonical_allele_internal(v);
            Ara->idx = k;
            Allele_Array[k] = Ara;

        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
}

void MarkerScheme_table::db_getall(linkage_ped_top *Top) {
    int ret = select_stmt && select_stmt->abort();
    int k = 0, allele1 = 0, allele2 = 0;

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {

            ret = select(k, allele1, allele2);
            MARKER_SCHEME3_alleles[k].Allele_1 = allele1;
            MARKER_SCHEME3_alleles[k].Allele_2 = allele2;

        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
}

void dballele_export(linkage_ped_top *Top, bp_order *bp) {
    int i;
    bp_order *b;

    Tod pedexp("export canonical allele/ marker scheme");

    MasterDB.begin();

    int_table.insert("allele_count", allele_count);
    int_table.insert("ALLELE_ARRAY", ALLELE_ARRAY);

    for (i = 0; i < allele_count; i++) {
        canonicalallele_table.insert(Allele_Array[i]->idx,
                                     allele_prop_allele(Allele_Array[i]));
    }

    if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int size = Top->LocusTop->MarkerCnt;
        int offset = Top->LocusTop->PhenoCnt;
        for (i = offset, b = bp_sort + offset; i < offset + size; i++, b++) {
            markerscheme_table.insert(i,
                                      MARKER_SCHEME3_alleles[b->i].Allele_1,
                                      MARKER_SCHEME3_alleles[b->i].Allele_2);
        }
    }

    MasterDB.commit();
    pedexp();

}



void dballele_import(linkage_ped_top *Top) {

    int_table.get("allele_count", allele_count);
    int_table.get("ALLELE_ARRAY", ALLELE_ARRAY);

    if (Allele_Array) free(Allele_Array);
    Allele_Array = CALLOC((size_t) ALLELE_ARRAY, allele_prop *);

    if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int size   = Top->LocusTop->MarkerCnt;
        int offset = Top->LocusTop->PhenoCnt;

        MARKER_SCHEME3_alleles = CALLOC(size + offset, Alleles_int);
    }

    MasterDB.begin();

//xxcanonicalallele_table.db_getall(Top);

    markerscheme_table.db_getall(Top);

    MasterDB.commit();

}
