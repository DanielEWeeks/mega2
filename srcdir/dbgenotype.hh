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

#ifndef DBGENOTYPE_HH
#define DBGENOTYPE_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Phenotype_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
    long long xhash(void *v, int sz) {
        if (v == 0) return 0;
        const char *cp = (const char *)v;
        long long r = 5381;
        int c;
        for (int i = 0; i < sz; i++) {
            c = *cp++;
            r = ((r << 5) + r) + c;
        }
        return r;
    }
public:
    Phenotype_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS phenotype_table (pId INTEGER PRIMARY KEY,"
            " person_link INTEGER, bytes INTEGER, data BLOB"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO phenotype_table("
            " person_link, bytes, data"
            ")  VALUES(?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT "
            " person_link, bytes, data"
            "  FROM phenotype_table;");
    }
    int insert(linkage_ped_rec *p, int cnt) {
        int idx = 1;
        int size = cnt * sizeof(pheno_pedrec_data);
//        long long ll = xhash(p->Pheno, size);
//      printf("PHI: %d %llx\n", p->person_link, ll);
        return insert_stmt 
            && insert_stmt->rowbind(idx, p->person_link, size)
            && insert_stmt->bind(idx++, (const void *)(p->Pheno), size)

            && insert_stmt->step();
    }
    int select(int &link, int &bytes, unsigned char * &data) {
        int idx = 0;
        int sz = 0;
        const void *v = 0;
        int ret = select_stmt 
            && select_stmt->row(idx, link, bytes)
            && select_stmt->column(idx++, v, sz);
        data = (unsigned char *) v;
//        long long ll = xhash(data, bytes);
//      printf("PHS: %d %llx\n", link, ll);
        return ret;
    }
    void print(linkage_ped_rec *p, linkage_ped_top *Top) {
        int cnt = Top->LocusTop->PhenoCnt;
        int bytes = cnt * sizeof(pheno_pedrec_data);
        long long ll = xhash(p->Pheno, bytes);
        printf("PHP %2d %3d ", p->pedigree_link, p->person_link);
#ifdef MS_PRINTF
        printf("%d %d %p %I64d\n", cnt, bytes, p->Pheno, ll);
#else
        printf("%d %d %p %lld\n", cnt, bytes, p->Pheno, ll);
#endif
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS phenotype_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_phenotype_table IF NOT EXISTS on phenotype_table (Name);");
    }

    int db_getall(linkage_locus_top *LTop, pheno_pedrec_data **Phenotypes);

};

extern Phenotype_table phenotype_table;

class Genotype_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
    long long xhash(void *v, int sz) {
        if (v == 0) return 0;
        const char *cp = (const char *)v;
        long long r = 5381;
        int c;
        for (int i = 0; i < sz; i++) {
            c = *cp++;
            r = ((r << 5) + r) + c;
        }
        return r;
    }
public:
    Genotype_table()  {}
    int create() {
	return MasterDB.exec(
	    "CREATE TABLE IF NOT EXISTS genotype_table (pId INTEGER PRIMARY KEY,"
            " person_link INTEGER, bytes INTEGER, data BLOB"
            ");"
	    );
    }
    void init () {
	insert_stmt = MasterDB.prep(
            "INSERT INTO genotype_table("
            " person_link, bytes, data"
            ")  VALUES(?, ?, ?);");
	select_stmt = MasterDB.prep(
            "SELECT "
            " person_link, bytes, data"
            "  FROM genotype_table;");
    }
    int insert(linkage_ped_rec *p, int cnt, int offset) {
        extern int marker_size(int);
        int idx = 1;
        int size = marker_size(cnt);
        const void *v = (const void *)marker_start(p->Marker, offset);
//NOTE: p->Marker can be 0
//        long long ll = xhash(p->Marker, size);
//      printf("GHI: %d %llx\n", p->person_link, ll);
        return insert_stmt 
            && insert_stmt->rowbind(idx, p->person_link, size)
            && insert_stmt->bind(idx++, v, size)

            && insert_stmt->step();
    }
    int select(int &link, int &bytes, unsigned char * &data) {
        int idx = 0;
        int sz = 0;
        const void *v = 0;
        int ret = select_stmt 
            && select_stmt->row(idx, link, bytes)
            && select_stmt->column(idx++, v, sz);
        data = (unsigned char *)v;
//        long long ll = xhash(data, bytes);
//      printf("GHS: %d %llx\n", link, ll);
        return ret;
    }
    void print(linkage_ped_rec *p, linkage_ped_top *Top) {
        int cnt = Top->LocusTop->MarkerCnt;
        int bytes = marker_size(cnt);
        long long ll = xhash(p->Marker, bytes);
        printf("GHP %2d %3d ", p->pedigree_link, p->person_link);
#ifdef MS_PRINTF
        printf("%d %d %p %I64d\n", cnt, bytes, p->Marker, ll);
#else
        printf("%d %d %p %lld\n", cnt, bytes, p->Marker, ll);
#endif
    }

    void close() {
	delete insert_stmt;
	delete select_stmt;
    }
    int drop() {
	return MasterDB.exec("DROP TABLE IF EXISTS genotype_table;");
    }

    int index() {
	return MasterDB.exec("CREATE Index Idx_genotype_table IF NOT EXISTS on genotype_table (UniqueID);");
    }

    int db_getall(linkage_locus_top *LTop, void **Genotypes);
};

extern Genotype_table genotype_table;

extern void dbgenotype_export(linkage_ped_top *Top, bp_order *bp);
extern void dbgenotype_import(linkage_ped_top *Top);

#endif
