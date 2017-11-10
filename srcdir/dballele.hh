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

#ifndef DBALLELE_HH
#define DBALLELE_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class CanonicalAllele_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    CanonicalAllele_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS canonicalallele_table (pId INTEGER PRIMARY KEY,"
            " key INTEGER, value TEXT"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO canonicalallele_table(key, value) VALUES(?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, value FROM canonicalallele_table;");
    }
    int insert(int key, char *value) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, value)

            && insert_stmt->step();
    }
    int select(int &key, const char * &value) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, value );
    }
    void print(int key, char *value) {
        printf("V %d %s\n", key, value);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS canonicalallele_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_canonicalallele_table on canonicalallele_table (key);");
    }

    void db_getall(linkage_ped_top *Top);

};

extern CanonicalAllele_table canonicalallele_table;

class MarkerScheme_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    MarkerScheme_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS markerscheme_table (pId INTEGER PRIMARY KEY,"
            " key INTEGER, allele1 INTEGER, allele2 INTEGER"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO markerscheme_table(key, allele1, allele2) VALUES(?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, allele1, allele2 FROM markerscheme_table;");
    }
    int insert(int key, int allele1, int allele2) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, allele1, allele2)

            && insert_stmt->step();
    }
    int select(int &key, int &allele1, int &allele2) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, allele1, allele2 );
    }
    void print(int key, int allele1, int allele2) {
        printf("M %d %d %d\n", key, allele1, allele2);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS markerscheme_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_markerscheme_table on markerScheme_table (key);");
    }

    void db_getall(linkage_ped_top *Top);

};

extern MarkerScheme_table markerscheme_table;

extern void dballele_export(linkage_ped_top *Top, bp_order *bp);
extern void dballele_import(linkage_ped_top *Top);

#endif
