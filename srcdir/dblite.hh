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

#ifndef DBLITE_HH
#define DBLITE_HH

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern void msgvf(const char *fmt, ...);

#include <string>

#define SQLITE_CORE 1
#include "sqlite3ext.h"

/*
//#define MT
*/
#ifdef MT
    static void db_mutex_enter(sqlite3 *d) {sqlite3_mutex_enter(sqlite3_db_mutex(d));}
    static void db_mutex_leave(sqlite3 *d)  {sqlite3_mutex_leave (sqlite3_db_mutex(d));}
#else
    static void db_mutex_enter(sqlite3 *d) {}
    static void db_mutex_leave(sqlite3 *d)  {}
#endif

class DBstmt {
private:
    sqlite3      *db;
    const char   *dbfile;
    const char   *sqlprep;
    sqlite3_stmt *stmt;
    int verbose;
public:
    int columns;

    DBstmt(sqlite3 *db_, const char *dbfile_): db(db_), dbfile(dbfile_), stmt(NULL), verbose(0) {}
   ~DBstmt() {
       if (stmt)
           sqlite3_finalize(stmt);
       stmt = 0;
    }

    int ok(int ret, const char *op) {
        if (ret == SQLITE_OK) return 1;
        msgvf("DB %s \"%s\" (file \"%s\") failed: #%d:%s\n",
              op, sqlprep, dbfile, ret, sqlite3_errmsg(db));
        return 0;
    }
    int ok(int ret, const char *op, int flag) {
        if (ret == flag) return 1;
        msgvf("DB %s \"%s\" (file \"%s\") failed: returned %d %s\n",
              op, sqlprep, dbfile, ret, sqlite3_errmsg(db));
        return 0;
    }
    int okt(int ret, int flag) {
        return (ret == flag) ? 1 : 0;
    }
    int okf(int ret, const char *op) {
        msgvf("DB %s \"%s\" (file \"%s\") failed: returned %d %s\n",
              op, sqlprep, dbfile, ret, sqlite3_errmsg(db));
        return 0;
    }

    int abort() {
        db_mutex_enter(db);
        int ret = ok(sqlite3_reset(stmt), "reset");
        db_mutex_leave(db);
        return ret;
    };

    int step() {
        db_mutex_enter(db);
        int ans = sqlite3_step(stmt);

        if (verbose > 2) printf("\n");
        if (okt(ans, SQLITE_ROW)) {
            columns = sqlite3_column_count(stmt);
            if (verbose > 1) printf("step: %d\n", ans);

        } else if (okt(ans, SQLITE_DONE)) {
            sqlite3_reset(stmt);
            if (verbose > 1) printf("done: %d\n\n", ans);

        } else if (okf(ans, "step")) {
            sqlite3_finalize(stmt);
            if (verbose > 1) printf("err!: %d\n", ans);
        }
        db_mutex_leave(db);
        return ans;
    };

    void fini() {
        delete this;
    }

    int step1() {
        int ret = 0;

        return ret;
    }

    int stepall() {
        int ret = 0;

        return ret;
    }


    template<typename a>
    int rowbind(int &k, a v1) {
        return (bind(k++, v1));
    }
    template<typename a, typename b>
    int rowbind(int &k, a v1, b v2) {
        return (bind(k++, v1) && bind(k++, v2));
    }
    template<typename a, typename b, typename c>
    int rowbind(int &k, a v1, b v2, c v3) {
        return (bind(k++, v1) && bind(k++, v2) && bind(k++, v3));
    }
    template<typename a, typename b, typename c, typename d>
    int rowbind(int &k, a v1, b v2, c v3, d v4) {
        return (bind(k++, v1) && bind(k++, v2) && bind(k++, v3) && bind(k++, v4));
    }

    int bind(int col, const char *cp) {
        if (verbose > 1) printf("bind: %d %s\n", col, cp);
        return ok(sqlite3_bind_text(stmt, col, cp, strlen(cp), SQLITE_STATIC), "bind");
//          return ok(sqlite3_bind_text(stmt, col, cp, strlen(cp), SQLITE_TRANSIENT), "bind");
    }
    int bind(int col, char *cp) {
        if (verbose > 1) printf("bind: %d %s\n", col, cp);
        return ok(sqlite3_bind_text(stmt, col, cp, strlen(cp), SQLITE_STATIC), "bind");
//            return ok(sqlite3_bind_text(stmt, col, cp, strlen(cp), SQLITE_TRANSIENT), "bind");
    }
    int bind(int col, int val) {
        if (verbose > 1) printf("bind: %d %d\n", col, val);
        return ok(sqlite3_bind_int(stmt, col, val), "bind");
    }
    int bind(int col, long long val) {
#ifdef MS_PRINTF
        if (verbose > 1) printf("bind: %d %I64d\n", col, val);
#else
        if (verbose > 1) printf("bind: %d %lld\n", col, val);
#endif
        return ok(sqlite3_bind_int64(stmt, col, val), "bind");
    }
    int bind(int col, double val) {
        if (verbose > 1) printf("bind: %d %f\n", col, val);
        return ok(sqlite3_bind_double(stmt, col, val), "bind");
    }
    int bind(int col, const void *vp, int sz) {
        if (verbose > 1) printf("bind: %d %s\n", col, (char *)vp);
        return ok(sqlite3_bind_blob(stmt, col, vp, sz, SQLITE_STATIC), "bind");
    }


    template<typename a>
    int row(int &k, a &v0) {
        return (column(k++, v0));
    }
    template<typename a, typename b>
    int row(int &k, a &v0, b &v1) {
        return (column(k++, v0) && column(k++, v1));
    }
    template<typename a, typename b, typename c>
    int row(int &k, a &v0, b &v1, c &v2) {
        return (column(k++, v0) && column(k++, v1) && column(k++, v2));
    }
    template<typename a, typename b, typename c, typename d>
    int row(int &k, a &v0, b &v1, c &v2, d &v3) {
        return (column(k++, v0) && column(k++, v1) && column(k++, v2) && column(k++, v3));
    }
    template<typename a, typename b, typename c, typename d, typename e, typename f, typename g>
    int row(int &k, a &v0, b &v1, c &v2, d &v3, e &v4, f &v5, g &v6) {
        return (column(k++, v0) && column(k++, v1) && column(k++, v2) && column(k++, v3) && column(k++, v4) && column(k++, v5) && column(k++, v6));
    }

    int column(int col, const char * &val) {
        val = (const char *)sqlite3_column_text(stmt, col);
        if (verbose > 1) printf("colm: %d %s\n", col, val);
        return 1;
    }
    int column(int col, const unsigned char * &val) {
        val = sqlite3_column_text(stmt, col);
        if (verbose > 1) printf("colm: %d %s\n", col, val);
        return 1;
    }
    int column(int col, char * &val) {
        val = (char *)sqlite3_column_text(stmt, col);
        if (verbose > 1) printf("colm: %d %s\n", col, val);
        return 1;
    }
    typedef char strA[25];
    int column(int col, strA &val) {
        strcpy(val, (const char *)sqlite3_column_text(stmt, col));
        if (verbose > 1) printf("colm: %d %s\n", col, val);
        return 1;
    }
    typedef char strB[255];
    int column(int col, strB &val) {
        strcpy(val, (const char *)sqlite3_column_text(stmt, col));
        if (verbose > 1) printf("colm: %d %s\n", col, val);
        return 1;
    }
    int column(int col, int &val) {
        val = sqlite3_column_int(stmt, col);
        if (verbose > 1) printf("colm: %d %d\n", col, val);
        return 1;
    }
    int column(int col, long long &val) {
        val = sqlite3_column_int64(stmt, col);
#ifdef MS_PRINTF
        if (verbose > 1) printf("colm: %d %I64d\n", col, val);
#else
        if (verbose > 1) printf("colm: %d %lld\n", col, val);
#endif
        return 1;
    }
    int column(int col, double &val) {
        val = sqlite3_column_double(stmt, col);
        if (verbose > 1) printf("colm: %d %f\n", col, val);
        return 1;
    }
    int column(int col, const void * &vp, int &sz) {
        vp = sqlite3_column_text(stmt, col);
        sz = sqlite3_column_bytes(stmt, col);
        if (verbose > 1) printf("colm: %d L%d %s\n", col, sz, (char *)vp);
        return 1;
    }


    int prep(const char *sql, int verbose_) {
        return prep(sql, -1, 0, verbose_);
    }
    int prep(const char *sql, int sqllen, const char * *sqlp, int verbose_) {
        int ret = 0;
        verbose = verbose_;
        sqlprep = sql;
        if (verbose > 2) printf("\n");
        if (verbose) printf("prep: %s\n", sql);

        db_mutex_enter(db);
        if (ok(sqlite3_prepare_v2(db, sql, sqllen, &stmt, sqlp), "prep")) ret = 1;
        else stmt = 0;
        db_mutex_leave(db);
        return ret;
    }
};


class DBlite {
private:
    sqlite3 *db;
    const char *dbfile;
    int autocommit;
    int verbose;

public:
    DBlite(): db(0), verbose(0) {
/*
#ifdef MT
        sqlite3_config(SQLITE_CONFIG_SERIALIZED);  // also serializes db and prep
#else
        sqlite3_config(SQLITE_CONFIG_MULTITHREAD); // ok for step
#endif
*/
    }
    DBlite(int verbose_): db(0), verbose(verbose_) {
        if (verbose_)
            msgvf("DB single_thread = %d\n", ! sqlite3_threadsafe());
    }
   ~DBlite() {
       if (db) close();
    }

    int ok(int ret, const char *op) {
        int myret = 0;
        if (ret == SQLITE_OK) myret = 1;
        else
            msgvf("DB %s (file \"%s\") failed: #%d:%s\n",
                  op, dbfile, ret, sqlite3_errmsg(db));
        return myret;
    }
    int ok(int ret, const char *op, const char *msg) {
        if (! msg) return ok(ret, op);
        int myret = 0;
        if (ret == SQLITE_OK) myret = 1;
        else {
            msgvf("DB %s (file \"%s\") failed: #%d:%s\n",
                  op, dbfile, ret, msg);
            sqlite3_free((void *)msg);
        }
        return myret;
    }

    int open(const char *file) {
        int ret = 0;

        if (file == 0 || *file == 0)
            file = (char*)":memory:";

        dbfile = file;
//      db_mutex_enter(db)
//      if (ok(sqlite3_open_v2(file, &db, SQLITE_OPEN_READWRITE, 0), "open"))
        if (ok(sqlite3_open(file, &db), "open")) {
#ifndef HIDEFILE
//          msgvf("DB open (file \"%s\") succeeded. SQLite3 lib version %s\n",
//                dbfile, sqlite3_libversion());
#endif
            ret = 1;
        }
//      db_mutex_leave(db);
        autocommit = sqlite3_get_autocommit(db);
        return ret;
    }

    int open(const char *file, int v) {
        verbose = v;
        return open(file);
    }

    int close() {
        int ret = 0;
        db_mutex_enter(db);
        if (ok(sqlite3_close(db), "close")) ret = 1;
        db = 0;
        db_mutex_leave(db);
        return ret;
    }

    DBstmt *prep(const char *sql) {
        DBstmt *dbstmt = new DBstmt(db, dbfile);

        if (! dbstmt->prep(sql, verbose)) {
            delete dbstmt;
            dbstmt = 0;
        }
        return dbstmt;
    }

    int exec(const char *sql) {
        char *error_msg = 0;
        if (verbose > 2) printf("\n\n");
        if (verbose) printf("exec: %s\n", sql);
        int ans = sqlite3_exec(db, sql, 0, 0, &error_msg);
        return (ok(ans, "sqlite3_exec", error_msg));
    }

    int begin()    { return (autocommit) ? exec("BEGIN TRANSACTION") : 1; }
    int commit()   { return (autocommit) ? exec("COMMIT") : 1; }
    int rollback() { return (autocommit) ? exec("ROLLBACK") : 1; }
};

#endif
