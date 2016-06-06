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

#ifndef DBBATCH_HH
#define DBBATCH_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Batch_file_parameters {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Batch_file_parameters()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS batch_parameters (Id INTEGER PRIMARY KEY, Key TEXT, Value TEXT, Read INTEGER);"
            );
    }
    void init () {
        insert_stmt = MasterDB.prep("INSERT INTO batch_parameters(Key, Value, Read) VALUES(?, ?, ?);");
        select_stmt = MasterDB.prep("SELECT Key, Value, Read FROM batch_parameters WHERE Key = ?;");
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS batch_parameters;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_batch_parameters IF NOT EXISTS on batch_parameters (key);");
    }

    int db_get(const char *&key, const char *&value, int &read)  {
        int idx = 1;
        int ret = 0;
        ret = select_stmt && select_stmt->abort() && select_stmt->rowbind(idx, key);
        if (ret) {
            idx = 0;
            ret = select_stmt->step();
            if (ret == SQLITE_ROW) {
//                  select_stmt->row(1, value, read);
                select_stmt->row(idx, key, value, read);
            } else if (ret == SQLITE_DONE) {
                value = "NOT FOUND";
                read = 0;
            } else {
                value = "ERROR";
                read = 0;
            }
        }
        return ret;
    }

    int db_get(const char *key, Vecs &vec, int &read)
    {
        int idx = 1;
        int ret = 0;
        vec.clear();
        ret = select_stmt && select_stmt->abort() && select_stmt->rowbind(idx, key);
        const char *value;
        while (ret) {
            idx = 1;
            ret = select_stmt->step();
            if (ret == SQLITE_ROW) {
                ret = select_stmt->row(idx, value, read);
                vec.push_back(value);
            } else if (ret == SQLITE_DONE) {
                ret = 0;
            } else {
                value = "ERROR";
                read = 0;
                ret = 0;
            }
        }
        return ret;
    }
};

extern Batch_file_parameters batch_file_parameters;

class File_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    File_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS file_table (pId INTEGER PRIMARY KEY,"
            " key INTEGER, file TEXT, type TEXT"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO file_table(key, file, type) VALUES(?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT key, file, type FROM file_table;");
    }
    int insert(int key, char *file, char* type) {
        int idx = 1;
        if (file == 0) file = (char *)"";
        return insert_stmt 
            && insert_stmt->rowbind(idx, key, file, type)

            && insert_stmt->step();
    }
    int select(int &key, char * &file, char * &type) {
        int idx = 0;
        return select_stmt 
            && select_stmt->row(idx, key, file, type );
    }
    void print(int key, char *file, char *type) {
        printf("F %d %s %s\n", key, file, type);
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS file_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_file_table IF NOT EXISTS on file_table (key);");
    }

    int db_getall();

    void stat();
};

extern File_table file_table;

extern void dbbatch_file_export(linkage_ped_top *Top);

extern void dbbatch_file_import(linkage_ped_top *Top);

#endif
