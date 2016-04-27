/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

#ifndef DBMISC_HH
#define DBMISC_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Int_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Int_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS int_table (pId INTEGER PRIMARY KEY,"
            " key TEXT, value INTEGER"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO int_table(key, value) VALUES(?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, value FROM int_table;");
    }
    int insert(const char* key, int value) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, value)

            && insert_stmt->step();
    }
    int select(const char * &key, int &value) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, value );
    }
    void print(const char* key, int value) {
        printf("V %s %d\n", key, value);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS int_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_int_table IF NOT EXISTS on int_table (key);");
    }

    int db_getall();

    void get(const char *key, int &val) {
        extern std::map<const char *, int, charsless> Int_hash;
        if (!map_get(Int_hash, key, val)) {
            printf("Int read failed for %s\n", key);
        }
    }

};

extern Int_table int_table;

class Double_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Double_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS double_table (pId INTEGER PRIMARY KEY,"
            " key TEXT, value DOUBLE"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO double_table(key, value) VALUES(?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, value FROM double_table;");
    }
    int insert(const char* key, double value) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, value)

            && insert_stmt->step();
    }
    int select(const char * &key, double &value) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, value );
    }
    void print(const char* key, double value) {
        printf("V %s %f\n", key, value);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS int_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_double_table IF NOT EXISTS on double_table (key);");
    }

    int db_getall();

};

extern Double_table double_table;

class Charstar_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Charstar_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS charstar_table (pId INTEGER PRIMARY KEY,"
            " key TEXT, value TEXT"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO charstar_table(key, value) VALUES(?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, value FROM charstar_table;");
    }
    int insert(const char* key, char *value) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, value)

            && insert_stmt->step();
    }
    int select(const char * &key, char * &value) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, value );
    }
    void print(const char* key, char *value) {
        printf("V %s %s\n", key, value);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS charstar_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_charstar_table IF NOT EXISTS on charstar_table (key);");
    }

    int db_getall();

};

extern Charstar_table charstar_table;

class Stuff_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Stuff_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS stuff_table (pId INTEGER PRIMARY KEY,"
            " key TEXT, value BLOB"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO stuff_table(key, value) VALUES(?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, value FROM stuff_table;");
    }
    int insert(const char* key, unsigned char *value, int bytes) {
        int idx = 1;
        return insert_stmt 
            && insert_stmt->rowbind(idx, key)
            && insert_stmt->bind(idx++, value, bytes)

            && insert_stmt->step();
    }
    int select(const char * &key, const unsigned char * &value, int &bytes) {
        int idx = 0;
        bytes = 0;
        const void *v = 0;
        int ret = select_stmt 
            && select_stmt->row(idx, key)
            && select_stmt->column(idx++, v, bytes);
        value = (unsigned char *) v;
        return ret;
    }
    void print(const char* key, unsigned char *value, int bytes) {
        printf("V %s %p L%d\n", key, value, bytes);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS stuff_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_stuff_table IF NOT EXISTS on stuff_table (key);");
    }

    int db_getall();

};

extern Stuff_table stuff_table;

extern void dbmisc_export(linkage_ped_top *Top);
extern void dbmisc_import(linkage_ped_top *Top);

#endif
