/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

#ifndef DBPEDIGREE_HH
#define DBPEDIGREE_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;

class Pedigree_table {
protected:
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
//    linkage_ped_rec *Entry;         /* will be Entry[] */
//    listhandle *Loops; /* will be NULL if no loops */
public:
    Pedigree_table()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS pedigree_table (pId INTEGER PRIMARY KEY,"
                " Num Integer, EntryCnt INTEGER, isTyped INTEGER, Name TEXT, PedPre TEXT,"
                " OriginalID INTEGER, origped INTEGER, Proband INTEGER, Loops BLOB, "
                " pedigree_link INTEGER);"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO pedigree_table("
            " Num, EntryCnt, isTyped, Name, PedPre, OriginalID,"
            " origped, Proband, Loops, pedigree_link"
            ") VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " Num, EntryCnt, isTyped, Name, PedPre, "
            " OriginalID, origped, Proband, Loops, pedigree_link"
            "   FROM pedigree_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_ped_tree *t) {
        int idx = 1;

        t->Loops = 0;  // no need to save; BUT it would be a listhandle of arrays; YUCK
                       // i.e. blobs would not work.

        return insert_stmt
            && insert_stmt->rowbind(idx, t->Num, t->EntryCnt, t->IsTyped)
            && insert_stmt->rowbind(idx, t->Name, t->PedPre)
            && insert_stmt->rowbind(idx, t->OriginalID, t->origped, t->Proband)
            && insert_stmt->bind(idx++,  t->Loops, 0)
            && insert_stmt->rowbind(idx, t->pedigree_link)

            && insert_stmt->step();
    }
    int select(linkage_ped_tree *t) {
        int idx = 0;
        int sz = 0;
        const void *v = 0;
        int ret =
            select_stmt->row(idx, t->Num, t->EntryCnt, t->IsTyped)
            && select_stmt->row(idx, t->Name, t->PedPre)
            && select_stmt->row(idx, t->OriginalID, t->origped, t->Proband)
            && select_stmt->column(idx++, v, sz)
            && select_stmt->row(idx, t->pedigree_link);
            t->Loops = 0;
        return ret;
    }
    void print(linkage_ped_tree *t) {
        printf("\nF %2d %d %d %d, ", t->pedigree_link, t->Num, t->EntryCnt, t->IsTyped);
        printf("%s %s, ", t->Name, t->PedPre);
        printf("%d %d %d\n", t->OriginalID, t->origped, t->Proband);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS pedigree_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_pedigree_table on pedigree_table (Name);");
    }

    int db_getall(linkage_ped_tree *t, std::map<int, linkage_ped_tree *> &pedigree_hash);

};

class Pedigree_brkloop_table: public Pedigree_table {
public:
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS pedigree_brkloop_table (pId INTEGER PRIMARY KEY,"
                " Num Integer, EntryCnt INTEGER, isTyped INTEGER, Name TEXT, PedPre TEXT,"
                " OriginalID INTEGER, origped INTEGER, Proband INTEGER, Loops BLOB, "
                " pedigree_link INTEGER);"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO pedigree_brkloop_table("
            " Num, EntryCnt, isTyped, Name, PedPre, OriginalID,"
            " origped, Proband, Loops, pedigree_link"
            ") VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " Num, EntryCnt, isTyped, Name, PedPre, "
            " OriginalID, origped, Proband, Loops, pedigree_link"
            "   FROM pedigree_brkloop_table;");
        return insert_stmt && select_stmt && 1;
    }    
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS pedigree_brkloop_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_pedigree_brkloop_table on pedigree_brkloop_table (Name);");
    }
};

extern Pedigree_table pedigree_table;
extern Pedigree_brkloop_table pedigree_brkloop_table;

class Person_table 
{
protected:
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
/*
    pheno_pedrec_data   *Pheno;
    void                *Marker;
    void *TmpData;
    int *loopbreakers;
*/
    enum { CREATE=1, FIELDS, BIND, COLUMN } ;
public:
    Person_table() {}

    int create() {
        return
            MasterDB.exec(
                "CREATE TABLE IF NOT EXISTS person_table (pId INTEGER PRIMARY KEY,"
                       " UniqueID TEXT, OrigID TEXT, FamName TEXT, PerPre TEXT,"
                       " ID INTEGER, Father INTEGER, Mother INTEGER,"
                       " First_Offspring INTEGER, Next_PA_Sib INTEGER, Next_MA_Sib INTEGER,"
                       " Sex INTEGER, OrigProband INTEGER, genocnt INTEGER,"
                       " Orig_status INTEGER, Ngeno INTEGER, IsTyped INTEGER,"
                       " PercentTyped DOUBLE,"
                       " ext_ped_num INTEGER, ext_per_num INTEGER,"
                       " MZTwin INTEGER, DZTwin INTEGER, GroupX INTEGER, rec_num INTEGER,"
                       " loopbreakers BLOB,"
                       " pedigree_link INTEGER, person_link INTEGER);"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO person_table("
            " UniqueID, OrigID, FamName, PerPre,"
            " ID, Father, Mother, First_Offspring, Next_PA_Sib, Next_MA_Sib,"
            " Sex, OrigProband, genocnt, Orig_status, Ngeno, IsTyped,"
            " PercentTyped, ext_ped_num, ext_per_num,"
            " MZTwin, DZTwin, GroupX, rec_num, loopbreakers, pedigree_link, person_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " UniqueID, OrigID, FamName, PerPre,"
            " ID, Father, Mother, First_Offspring, Next_PA_Sib, Next_MA_Sib,"
            " Sex, OrigProband, genocnt, Orig_status, Ngeno, IsTyped,"
            " PercentTyped, ext_ped_num, ext_per_num,"
            " MZTwin, DZTwin, GroupX, rec_num, loopbreakers,"
            " pedigree_link, person_link" 
            "   FROM person_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_ped_rec *p) {
        int sz = 0, idx = 1;
        if (p->loopbreakers != NULL) {
            sz =(MAXLOOPMEMBERS + 1) * sizeof(int);
        }
        return insert_stmt 
            && insert_stmt->rowbind(idx, p->UniqueID, p->OrigID, p->FamName, p->PerPre)
            && insert_stmt->rowbind(idx, p->ID, p->Father, p->Mother)
            && insert_stmt->rowbind(idx, p->First_Offspring, p->Next_PA_Sib, p->Next_MA_Sib)
            && insert_stmt->rowbind(idx, p->Sex, p->OrigProband, p->genocnt)
            && insert_stmt->rowbind(idx, p->Orig_status, p->Ngeno, p->IsTyped)
            && insert_stmt->rowbind(idx, p->PercentTyped)
            && insert_stmt->rowbind(idx, p->ext_ped_num, p->ext_per_num)
            && insert_stmt->rowbind(idx, p->MZTwin, p->DZTwin, p->Group, p->rec_num)
            && insert_stmt->bind(idx++,  p->loopbreakers, sz)
            && insert_stmt->rowbind(idx, p->pedigree_link, p->person_link)

            && insert_stmt->step();
    }
    int select(linkage_ped_rec *p) {
        int idx = 0;
        int sz = 0;
        const void *v = 0;
        int ret = select_stmt 
            && select_stmt->row(idx, p->UniqueID, p->OrigID, p->FamName, p->PerPre)
            && select_stmt->row(idx, p->ID, p->Father, p->Mother)
            && select_stmt->row(idx, p->First_Offspring, p->Next_PA_Sib, p->Next_MA_Sib)
            && select_stmt->row(idx, p->Sex, p->OrigProband, p->genocnt)
            && select_stmt->row(idx, p->Orig_status, p->Ngeno, p->IsTyped)
            && select_stmt->row(idx, p->PercentTyped)
            && select_stmt->row(idx, p->ext_ped_num, p->ext_per_num)
            && select_stmt->row(idx, p->MZTwin, p->DZTwin, p->Group, p->rec_num)
            && select_stmt->column(idx++, v, sz)
            && select_stmt->row(idx, p->pedigree_link, p->person_link);
        p->loopbreakers = (int *)v;  // might be 0
        return ret;
    }
    void print(linkage_ped_rec *p) {
        printf("P %2d %3d ", p->pedigree_link, p->person_link);
        printf("%4s %4s %4s %4s, ", p->UniqueID, p->OrigID, p->FamName, p->PerPre);
        printf("%2d %2d %2d, ", p->ID, p->Father, p->Mother);
        printf("%2d %2d %2d, ", p->First_Offspring, p->Next_PA_Sib, p->Next_MA_Sib);
        printf("%2d %2d %2d, ", p->Sex, p->OrigProband, p->genocnt);
        printf("%2d %2d %2d, ", p->Orig_status, p->Ngeno, p->IsTyped);
        printf("%.6f, ", p->PercentTyped);
        printf("%2d %2d, ", p->ext_ped_num, p->ext_per_num);
        printf("%2d %2d %2d %2d\n", p->MZTwin, p->DZTwin, p->Group, p->rec_num);
    }

    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS person_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_person_table on person_table (UniqueID);");
    }

    int db_getall(linkage_ped_rec *p);
};

class Person_brkloop_table: public Person_table {
public:
    int create() {
        return
            MasterDB.exec(
                "CREATE TABLE IF NOT EXISTS person_brkloop_table (pId INTEGER PRIMARY KEY,"
                       " UniqueID TEXT, OrigID TEXT, FamName TEXT, PerPre TEXT,"
                       " ID INTEGER, Father INTEGER, Mother INTEGER,"
                       " First_Offspring INTEGER, Next_PA_Sib INTEGER, Next_MA_Sib INTEGER,"
                       " Sex INTEGER, OrigProband INTEGER, genocnt INTEGER,"
                       " Orig_status INTEGER, Ngeno INTEGER, IsTyped INTEGER,"
                       " PercentTyped DOUBLE,"
                       " ext_ped_num INTEGER, ext_per_num INTEGER,"
                       " MZTwin INTEGER, DZTwin INTEGER, GroupX INTEGER, rec_num INTEGER,"
                       " loopbreakers BLOB,"
                       " pedigree_link INTEGER, person_link INTEGER);"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO person_brkloop_table("
            " UniqueID, OrigID, FamName, PerPre,"
            " ID, Father, Mother, First_Offspring, Next_PA_Sib, Next_MA_Sib,"
            " Sex, OrigProband, genocnt, Orig_status, Ngeno, IsTyped,"
            " PercentTyped, ext_ped_num, ext_per_num,"
            " MZTwin, DZTwin, GroupX, rec_num, loopbreakers, pedigree_link, person_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,  ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " UniqueID, OrigID, FamName, PerPre,"
            " ID, Father, Mother, First_Offspring, Next_PA_Sib, Next_MA_Sib,"
            " Sex, OrigProband, genocnt, Orig_status, Ngeno, IsTyped,"
            " PercentTyped, ext_ped_num, ext_per_num,"
            " MZTwin, DZTwin, GroupX, rec_num, loopbreakers,"
            " pedigree_link, person_link" 
            "   FROM person_brkloop_table;");
        return insert_stmt && select_stmt && 1;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS person_brkloop_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index IF NOT EXISTS Idx_person_brkloop_table on person_brkloop_table (Name);");
    }
};

extern Person_table person_table;
extern Person_brkloop_table person_brkloop_table;

extern void dbpedigree_export(linkage_ped_top *Top);
extern void dbpedigree_import(linkage_ped_top *Top);

#endif
