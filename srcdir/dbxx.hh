/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

#ifndef DBXX_HH
#define DBXX_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

using namespace std;

extern DBlite MasterDB;

class xx1_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    xx1_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS marker_table (pId INTEGER PRIMARY KEY,"
                        ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO marker_table("
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT "
            "   FROM marker_table WHERE ID = ?;");
    }
    int insert(linkage_ped_rec *p) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx)

            && insert_stmt->step();
    }
    int select(linkage_ped_rec *p) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx )
            ;
    }
    void print(linkage_ped_rec *p) {
        printf("P %2d %3d ", 0, 0);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS xx1_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_xx1_table on xx1_table (Name);");
    }

    int db_getall(linkage_ped_tree *t);

};

extern xx1_table xx1_table;

class xx2_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    xx2_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS xx2_table (pId INTEGER PRIMARY KEY,"
                        ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO xx2_table("
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT "
            "   FROM xx2_table WHERE ID = ?;");
    }
    int insert(linkage_ped_rec *p) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx)

            && insert_stmt->step();
    }
    int select(linkage_ped_rec *p) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx )
            ;
    }
    void print(linkage_ped_rec *p) {
        printf("P %2d %3d ", 0, 0);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS xx2_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_xx2_table on xx2_table (UniqueID);");
    }

    int db_getall(linkage_ped_rec *p);
};

extern xx2_table xx2_table;

#endif
