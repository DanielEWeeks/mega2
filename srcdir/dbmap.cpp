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

#include "common.h"
#include "typedefs.h"
#include "error_messages_ext.h"
#include "utils_ext.h"

#include "compress_ext.h"
#include "linkage.h"

#include "dblite.hh"
#include "dbmisc.hh"
#include "dbmap.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;

extern map<const char *, int,             charsless> Int_hash;
extern map<const char *, double,          charsless> Double_hash;
extern map<const char *, char *,          charsless> Charstar_hash;
extern map<const char *, unsigned char *, charsless> Stuff_hash;

void Map_table::db_getall(linkage_ped_top *Top) {
    int ret = select_stmt && select_stmt->abort();
    int i = 0, j = 0;
    double position = 0, pos_female = 0, pos_male = 0;

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {

            ret = select(i, j, position, pos_female, pos_male);
            Top->EXLTop->EXLocus[i].positions[j]  = position;
            Top->EXLTop->EXLocus[i].pos_female[j] = pos_female;
            Top->EXLTop->EXLocus[i].pos_male[j]   = pos_male;

        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
}

void MapNames_table::db_getall(linkage_ped_top *Top) {
    int ret = select_stmt && select_stmt->abort();
    int j = 0;
    int sex_averaged_map = 0, male_sex_map = 0, female_sex_map = 0;
    char *name = (char *)0;

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {

            ret = select(j, sex_averaged_map, male_sex_map, female_sex_map, name);

            Top->EXLTop->SexMaps[j] = CALLOC((size_t) 3, int);
            Top->EXLTop->SexMaps[j][SEX_AVERAGED_MAP] = sex_averaged_map;
            Top->EXLTop->SexMaps[j][MALE_SEX_MAP]     = male_sex_map;
            Top->EXLTop->SexMaps[j][FEMALE_SEX_MAP]   = female_sex_map;
            Top->EXLTop->MapNames[j] = strdup(name);

        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
}

void dbmap_export(linkage_ped_top *Top, bp_order *bp) {
    int i, j;
    bp_order *b;

    Tod pedexp("export canonical allele/ marker scheme");

    ext_linkage_locus_top *EXLTop = Top->EXLTop;
    ext_linkage_locus_rec *EXLocus = EXLTop->EXLocus;
    MasterDB.begin();

    int_table.insert("MapCnt", EXLTop->MapCnt);
    int_table.insert("GroupCnt", EXLTop->GroupCnt);
//  charstar_table.insert("map_functions", (unsigned char *)EXLTop->map_functions);
    charstar_table.insert("map_functions", EXLTop->map_functions);

    for (j = 0; j < EXLTop->MapCnt; j++) {
        mapnames_table.insert(j, EXLTop->SexMaps[j][SEX_AVERAGED_MAP],
                              EXLTop->SexMaps[j][MALE_SEX_MAP],
                              EXLTop->SexMaps[j][FEMALE_SEX_MAP],
                              EXLTop->MapNames[j]);
    }

    for (i = Top->LocusTop->PhenoCnt, b = bp_sort + Top->LocusTop->PhenoCnt;
         i < Top->LocusTop->LocusCnt; i++, b++) {
        ext_linkage_locus_rec *EXL = EXLocus + b->i;
        for (j = 0; j < EXLTop->MapCnt; j++) {
            map_table.insert(i, j, EXL->positions[j], EXL->pos_female[j], EXL->pos_male[j]);
        }
    }

    MasterDB.commit();
    pedexp();
}

void dbmap_import(linkage_ped_top *Top) {
    int i;
    extern ext_linkage_locus_top *new_EXLTop(linkage_locus_top *LTop);
    ext_linkage_locus_top *EXLTop = new_EXLTop(Top->LocusTop);
    Top->EXLTop = EXLTop;

    linkage_locus_top *LTop = Top->LocusTop;
    EXLTop->LocusCnt = LTop->LocusCnt;
    ext_linkage_locus_rec *EXLocus = (LTop->MarkerCnt > 0) ? 
        (CALLOC((size_t) LTop->MarkerCnt, ext_linkage_locus_rec) - LTop->PhenoCnt) : 0;
    EXLTop->EXLocus = EXLocus;

/*
    const char *a;
    a = "MapCnt";
    if (!map_get(Int_hash, a, EXLTop->MapCnt)) {
        printf("Int read failed for %s\n", a);
    }
    a = "GroupCnt";
    if (!map_get(Int_hash, a, EXLTop->GroupCnt)) {
        printf("Int read failed for %s\n", a);
    }
    a = "map_functions";
    char *aa = (char *)0;
    if (!map_get(Stuff_hash, a, aa)) {
        printf("Charstar read failed for %s\n", a);
    }
    EXLTop->map_functions = strdup(aa);
*/
    int_table.get("MapCnt", EXLTop->MapCnt);
    int mapcnt = EXLTop->MapCnt;
    EXLTop->SexMaps  = CALLOC((size_t)mapcnt, int *);
    EXLTop->MapNames = CALLOC((size_t)mapcnt, char *);
    int_table.get("GroupCnt", EXLTop->GroupCnt);
    charstar_table.get("map_functions", EXLTop->map_functions);

    for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
        EXLocus[i].positions  = CALLOC((size_t)mapcnt, double);
        EXLocus[i].pos_male   = CALLOC((size_t)mapcnt, double);
        EXLocus[i].pos_female = CALLOC((size_t)mapcnt, double);
    }

    MasterDB.begin();

    map_table.db_getall(Top);

    mapnames_table.db_getall(Top);

    MasterDB.commit();

}
