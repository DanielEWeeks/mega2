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

#ifndef DBMAP_HH
#define DBMAP_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Map_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Map_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS map_table (pId INTEGER PRIMARY KEY,"
            " marker INTEGER, map INTEGER, position DOUBLE, pos_female DOUBLE, pos_male DOUBLE"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO map_table(marker, map, position, pos_female, pos_male) VALUES(?, ?, ?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT marker, map, position, pos_female, pos_male FROM map_table;");
    }
    int insert(int marker, int map, double position, double pos_female, double pos_male) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, marker, map)
            && insert_stmt->rowbind(idx, position, pos_female, pos_male)

            && insert_stmt->step();
    }
    int select(int &marker, int &map, double &position, double &pos_female, double &pos_male) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, marker, map)
            && select_stmt->row(idx, position, pos_female, pos_male);
    }
    void print(int marker, int map, double position, double pos_female, double pos_male) {
        printf("M %d %d %f %f %f\n", marker, map, position, pos_female, pos_male);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS map_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_map_table IF NOT EXISTS on map_table (key);");
    }

    void db_getall(linkage_ped_top *Top);

};

extern Map_table map_table;

class MapNames_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    MapNames_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS mapnames_table (pId INTEGER PRIMARY KEY,"
            " map INTEGER, sex_averaged_map INTEGER, male_sex_map INTEGER, female_sex_map INTEGER, name TEXT"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO mapnames_table( map, sex_averaged_map, male_sex_map, female_sex_map, name) VALUES(?, ?, ?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT map, sex_averaged_map, male_sex_map, female_sex_map, name FROM mapnames_table;");
    }
    int insert(int map, int sex_averaged_map, int male_sex_map, int female_sex_map, char *name) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, map)
            && insert_stmt->rowbind(idx, sex_averaged_map, male_sex_map, female_sex_map)
            && insert_stmt->rowbind(idx, name)

            && insert_stmt->step();
    }
    int select(int &map, int &sex_averaged_map, int &male_sex_map, int &female_sex_map, char * &name) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, map)
            && select_stmt->row(idx, sex_averaged_map, male_sex_map, female_sex_map)
            && select_stmt->row(idx, name);
    }
    void print(int map, int sex_averaged_map, int male_sex_map, int female_sex_map, char *name) {
        printf("M %d %s %d %d %d\n", map, name, sex_averaged_map, male_sex_map, female_sex_map);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS map_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_map_table IF NOT EXISTS on map_table (key);");
    }

    void db_getall(linkage_ped_top *Top);

};

extern MapNames_table mapnames_table;

extern void dbmap_export(linkage_ped_top *Top);
extern void dbmap_import(linkage_ped_top *Top);

#endif
