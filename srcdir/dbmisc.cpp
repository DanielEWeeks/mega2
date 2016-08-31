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
#include "input_check.h"
#include "input_ops.hh"

#include "dblite.hh"
#include "dbmisc.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;
map<const char *, int,             charsless> Int_hash;
map<const char *, double,          charsless> Double_hash;
map<const char *, char *,          charsless> Charstar_hash;
map<const char *, unsigned char *, charsless> Stuff_hash;

int Int_table::db_getall() {
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            const char *k = (const char *)0;
            int v = 0;
            ret = select(k, v);
            Int_hash[strdup(k)] = v;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return ret;
}

int Double_table::db_getall() {
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            const char *k = (const char* )0;
            double v = 0;
            ret = select(k, v);
            Double_hash[strdup(k)] = v;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return ret;
}


int Charstar_table::db_getall() {
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            const char *k = (const char *)0;
            char *v = (char *)0;
            ret = select(k, v);
            Charstar_hash[strdup(k)] = strdup(v);
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return ret;
}


int Stuff_table::db_getall() {
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            const char *k = (const char *)0;
            const unsigned char *v = (const unsigned char *)0;
            unsigned char *u = (unsigned char *)0;;
            int bytes = 0;
            ret = select(k, v, bytes);
            if (bytes) {
                u = CALLOC((size_t) bytes, unsigned char);
                memcpy(u, v, bytes);
            }
            Stuff_hash[strdup(k)] = u;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return ret;
}


void dbmisc_export(linkage_ped_top *Top) {

    Tod pedexp("export genotype/phenotype");

    MasterDB.begin();

    extern char *SQLversion;
    charstar_table.insert("DBCreateTime", RunDate);
    charstar_table.insert("DBMega2Version", Mega2Version);
    charstar_table.insert("DBVersion", SQLversion);
    
    extern int seed1, seed2, seed3;
    int_table.insert("Seed1", seed1);
    int_table.insert("Seed2", seed2);
    int_table.insert("Seed3", seed3);

    extern INPUT_FORMAT_t Input_Format;
    int_table.insert("Input_Format", Input_Format);

    // for each genotype (consulting the linkage_ped_top structure)...
    linkage_locus_top *ELTop = Top->LocusTop;

    int_table.insert("PedCnt", Top->PedCnt);
    int_table.insert("UniqueIds", Top->UniqueIds);
    int_table.insert("OrigIds", Top->OrigIds);
    int_table.insert("LocusCnt", ELTop->LocusCnt);
    int_table.insert("PhenoCnt", ELTop->PhenoCnt);
    int_table.insert("MarkerCnt", ELTop->MarkerCnt);
    int_table.insert("Program", ELTop->Program);
    int_table.insert("SexDiff", ELTop->SexDiff);
    int_table.insert("SexLinked", ELTop->SexLinked);
    int_table.insert("MapDistanceType", ELTop->map_distance_type);

//  calculated elsewhere
//  int_table.insert("AllelesCnt", ...)
//  int_table.insert("AffectClassCnt", ...);

    extern int MARKER_SCHEME, MARKER_SCHEME3_offset;
    int_table.insert("MARKER_SCHEME", MARKER_SCHEME);
    int_table.insert("MARKER_SCHEME3_offset", MARKER_SCHEME3_offset);
    int_table.insert("HasMarkers", HasMarkers);

    extern int MaxChromo, NumUnmapped;
    int_table.insert("MaxChromo", MaxChromo);
    int_table.insert("NumUnmapped", NumUnmapped);
    int_table.insert("basefile_type", basefile_type);

    extern int human_unknown, human_x, human_xy, human_y, human_mt;
    int_table.insert("human_unknown", human_unknown);
    int_table.insert("human_x",  human_x);
    int_table.insert("human_xy", human_xy);
    int_table.insert("human_y",  human_y);
    int_table.insert("human_mt", human_mt);

    extern ped_status PedStat;
    int_table.insert("reset_genotype_invalid", PedStat.genotype_invalid);
    int_table.insert("reset_halftyped",        PedStat.halftyped);
    int_table.insert("reset_exceed_allcnt",    PedStat.exceed_allcnt);
    int_table.insert("HalfTypedReset",         HalfTypedReset);

    extern int genetic_distance_index;
    extern int genetic_distance_sex_type_map;
    extern int base_pair_position_index;
    int_table.insert("genetic_distance_index", genetic_distance_index);
    int_table.insert("genetic_distance_sex_type_map", genetic_distance_sex_type_map);
    int_table.insert("base_pair_position_index", base_pair_position_index);

    stuff_table.insert("MaleRecomb", (unsigned char *)(ELTop->MaleRecomb), sizeof(double) * (ELTop->LocusCnt - 1));
    if (ELTop->Program == LINKMAP)
        stuff_table.insert("recomb_frac", (unsigned char *)(ELTop->Run.linkmap.recomb_frac), sizeof(double) * (ELTop->LocusCnt - 1));
    else
        stuff_table.insert("recomb_frac", 0, 0);
    if (ELTop->SexDiff == VARIABLE_SEX_DIFF)
        stuff_table.insert("FemaleRecomb", (unsigned char *)(ELTop->FemaleRecomb), sizeof(double) * (ELTop->LocusCnt - 1));
    else
        stuff_table.insert("FemaleRecomb", 0, 0);

    MasterDB.commit();

    pedexp();
}


linkage_ped_top SQLop;
linkage_locus_top SLop;

void dbmisc_import(linkage_ped_top *Top) {

    Top->LocusTop = &SLop;
    linkage_locus_top *ELTop = &SLop;

    MasterDB.begin();

    int_table.db_getall();
    double_table.db_getall();
    charstar_table.db_getall();
    stuff_table.db_getall();

    MasterDB.commit();

    extern char *DBCreateTime, *DBMega2Version, *DBversion;
    charstar_table.get("DBCreateTime", DBCreateTime);
    charstar_table.get("DBMega2Version", DBMega2Version);
    charstar_table.get("DBVersion", DBversion);


#ifdef TEST
#define SETSEED
#endif
#ifdef SETSEED
    extern int seed1, seed2, seed3;
    int_table.get("Seed1", seed1);
    int_table.get("Seed2", seed2);
    int_table.get("Seed3", seed3);
#endif

    extern INPUT_FORMAT_t Input_Format;
    int Input_tmp = 0;
    int_table.get("Input_Format", Input_tmp);
    Input_Format = (INPUT_FORMAT_t)Input_tmp;

    int_table.get("PedCnt", Top->PedCnt);
    int_table.get("UniqueIds", Top->UniqueIds);
    int_table.get("OrigIds", Top->OrigIds);
    int_table.get("LocusCnt", ELTop->LocusCnt);
    int_table.get("PhenoCnt", ELTop->PhenoCnt);
    int_table.get("MarkerCnt", ELTop->MarkerCnt);

// defined elsewhere
    extern int AllelesCnt, AffectClassCnt;
    int_table.get("AllelesCnt", AllelesCnt);
    int_table.get("AffectClassCnt", AffectClassCnt);


//computed
    extern int MARKER_SCHEME, MARKER_SCHEME3_offset;
    int_table.get("MARKER_SCHEME", MARKER_SCHEME);
    int_table.get("MARKER_SCHEME3_offset", MARKER_SCHEME3_offset);
    int_table.get("HasMarkers", HasMarkers);

    extern int MaxChromo, NumUnmapped;
    int_table.get("MaxChromo", MaxChromo);
    int_table.get("NumUnmapped", NumUnmapped);
    int_table.get("basefile_type", basefile_type);

    extern int human_unknown, human_x, human_xy, human_y, human_mt;
    int_table.get("human_unknown", human_unknown);
    int_table.get("human_x",  human_x);
    int_table.get("human_xy", human_xy);
    int_table.get("human_y",  human_y);
    int_table.get("human_mt", human_mt);

    extern ped_status PedStat;
    int_table.get("reset_genotype_invalid", PedStat.genotype_invalid);
    int_table.get("reset_halftyped",        PedStat.halftyped);
    int_table.get("reset_exceed_allcnt",    PedStat.exceed_allcnt);
    int_table.get("HalfTypedReset",         HalfTypedReset);

//linkage hdr
    int_table.get("Program", ELTop->Program);
    int_table.get("SexDiff", ELTop->SexDiff);
    int_table.get("SexLinked", ELTop->SexLinked);

    int c = 0;
    int_table.get("MapDistanceType", c);
    Top->LocusTop->map_distance_type = (unsigned char)c;

//asked
    extern int genetic_distance_index;
    extern int genetic_distance_sex_type_map;
    extern int base_pair_position_index;
    int_table.get("genetic_distance_index", genetic_distance_index);
    int_table.get("genetic_distance_sex_type_map", genetic_distance_sex_type_map);
    int_table.get("base_pair_position_index", base_pair_position_index);

//needed

}
