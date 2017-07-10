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

#ifndef DBREFALLELE_H
#define DBREFALLELE_H

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Reference_Allele_Table {

protected:
    DBstmt *insert_stmt;
    DBstmt *select_stmt;

public:
    Reference_Allele_Table( ) { }

    ~Reference_Allele_Table( ) { }

    int create() {
        return MasterDB.exec("CREATE TABLE IF NOT EXISTS ref_allele_table(chr integer, pos  integer, marker integer, ref text, alt text)");
    }
    int init () {
        insert_stmt = MasterDB.prep("INSERT INTO ref_allele_table(chr, pos, marker, ref, alt) VALUES(?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep("SELECT chr, pos, marker, ref, alt FROM ref_allele_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(int chr, int pos, int marker, char *ref, char *alt) {
        int idx = 1;

        return insert_stmt
               && insert_stmt->rowbind(idx, chr, pos, marker)
               && insert_stmt->rowbind(idx, ref, alt)
               && insert_stmt->step();
    }
    int select(int chr, int pos, int marker, char *ref, char *alt) {
        int idx = 0;
        int ret =
                select_stmt->row(idx, chr, pos, marker)
                && select_stmt->row(idx, ref, alt);
        return ret;
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS ref_allele_table;");
    }

    int db_getall(linkage_ped_tree *t, std::map<int, linkage_ped_tree *> &pedigree_hash);

    void read_ref_allele_file(linkage_ped_top *Top, Str filename, bool use_bp_sort, bp_order *bp = NULL);
    Str get_filename();
    //int get_marker(int pos);
    //void insert_into_table(int chr, int pos, Str ref);

};

extern Reference_Allele_Table *reference_allele_table;


class Reference_Flips_Table {
protected:
    DBstmt *insert_stmt;
    DBstmt *select_stmt;

public:
    Reference_Flips_Table( ) { }

    ~Reference_Flips_Table( ) { }

    int create() {
        return MasterDB.exec("CREATE TABLE IF NOT EXISTS ref_allele_flips(marker integer, strand integer, major_minor integer, dummy, oldref, oldalt)");
    }
    int init () {
        insert_stmt = MasterDB.prep("INSERT INTO ref_allele_flips(marker, strand, major_minor, dummy, oldref, oldalt) VALUES(?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep("SELECT marker, strand, major_minor, dummy, oldref, oldalt FROM ref_allele_flips;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(int marker, int strand, int major_minor, int dummy, const char *oldref, const char *oldalt) {
        int idx = 1;

        return insert_stmt
               && insert_stmt->rowbind(idx, marker)
               && insert_stmt->rowbind(idx, strand, major_minor, dummy)
               && insert_stmt->rowbind(idx, oldref, oldalt)
               && insert_stmt->step();
    }
    int select(int marker, int strand, int major_minor, int dummy, const char *oldref, const char *oldalt) {
        int idx = 0;
        int ret =
                select_stmt->row(idx, marker)
                && select_stmt->row(idx, strand, major_minor, dummy)
                && select_stmt->row(idx, oldref, oldalt);
        return ret;
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS ref_allele_flips;");
    }

    int determine_flips(linkage_ped_top *Top, int locus, const char *data_ref, const char *data_alt, char *ref_ref, char *ref_alt, int chromosome, int position);

    void flip_strands(linkage_ped_top *Top);

    int db_getall(linkage_ped_tree *t, std::map<int, linkage_ped_tree *> &pedigree_hash);
};

extern Reference_Flips_Table *reference_flips_table;
#endif