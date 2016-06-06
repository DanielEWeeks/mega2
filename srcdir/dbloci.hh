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

#ifndef DBLOCI_HH
#define DBLOCI_HH

#include "types.hh"
#include "dblite.hh"

#include "common.h"
#include "typedefs.h"

extern DBlite MasterDB;
#if 1
class raw_table {
protected:
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
    const char *table;
    int idx;
public:
    raw_table(const char *tbl) : table(tbl) {}
    virtual const char *create_sql() = 0 ;
    virtual const char *init_insert_sql() = 0 ;
    virtual const char *init_select_sql() = 0 ;
    virtual int index() = 0;
    int create() {
        return MasterDB.exec(create_sql());
    }
    int init () {
        insert_stmt = MasterDB.prep(init_insert_sql());
        select_stmt = MasterDB.prep(init_select_sql());
        return insert_stmt && select_stmt && 1;
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
//      return MasterDB.exec(C("DROP TABLE IF EXISTS " + table + " ;"));
        return MasterDB.exec(C("DROP TABLE IF EXISTS locus_table ;"));
    }
};

class Locus_table : public raw_table {
public:
    Locus_table() : raw_table("locus_table") {};
    const char *create_sql() {
        return
            "CREATE TABLE IF NOT EXISTS locus_table (pId INTEGER PRIMARY KEY,"
            " LocusName TEXT, Type INTEGER, Class INTEGER, number INTEGER,"
            " col_num INTEGER, AlleleCnt INTEGER, locus_link INTEGER"
            ");";
    }
    const char *init_insert_sql () {
        return
            "INSERT INTO locus_table("
            " LocusName, Type, Class, number, col_num, AlleleCnt, locus_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?);";
    }
    const char *init_select_sql () {
        return
            "SELECT "
            " LocusName, Type, Class, number, col_num, AlleleCnt, locus_link"
            "   FROM locus_table;";
    }
    
/*    int io(linkage_locus_rec *p, DBstmt *stmt, int (DBstmt::*callme)(a v1)) {
        return io(p, insert_stmt, &DBstmt::rowbind);
*/
    int insert(linkage_locus_rec *p) {
        idx = 1;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->LocusName, p->Type, p->Class)
            && insert_stmt->rowbind(idx, p->number, p->col_num)
            && insert_stmt->rowbind(idx, p->AlleleCnt, p->locus_link)

            && insert_stmt->step();
    }
    int select(linkage_locus_rec *p) {
        const char *ccp;
        int xt;
        int xc;
        idx = 0;
        int ret = select_stmt->row(idx, ccp, xt, xc)
            && select_stmt->row(idx, p->number, p->col_num)
            && select_stmt->row(idx, p->AlleleCnt, p->locus_link);

        p->Type = (linkage_locus_type) xt;
        p->Class = (linkage_locus_class) xc;
        p->LocusName = strdup(ccp);
        return ret;
    }
    void print(linkage_locus_rec *p) {
        printf("\nL %d %s %d %d, ", p->locus_link, p->LocusName, p->Type, p->Class);
        printf("%d %d, ", p->number, p->col_num);
        printf("%d\n", p->AlleleCnt);
    }
    int index() {
        return MasterDB.exec("CREATE Index Idx_locus_table IF NOT EXISTS on locus_table (UniqueID);");
    }

    int db_getall(linkage_locus_rec *t);
};
#else

class Locus_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Locus_table()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS locus_table (pId INTEGER PRIMARY KEY,"
            " LocusName TEXT, Type INTEGER, Class INTEGER, number INTEGER,"
            " col_num INTEGER, AlleleCnt INTEGER, locus_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO locus_table("
            " LocusName, Type, Class, number, col_num, AlleleCnt, locus_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " LocusName, Type, Class, number, col_num, AlleleCnt, locus_link"
            "   FROM locus_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_locus_rec *p) {
        int idx = 1;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->LocusName, p->Type, p->Class)
            && insert_stmt->rowbind(idx, p->number, p->col_num)
            && insert_stmt->rowbind(idx, p->AlleleCnt, p->locus_link)

            && insert_stmt->step();
    }
    int select(linkage_locus_rec *p) {
        int idx = 0;
        const char *ccp;
        int xt = p->Type;
        int xc = p->Class;

        int ret = select_stmt->row(idx, ccp, xt, xc)
            && select_stmt->row(idx, p->number, p->col_num)
            && select_stmt->row(idx, p->AlleleCnt, p->locus_link);

        p->Type = (linkage_locus_type) xt;
        p->Class = (linkage_locus_class) xc;
        p->LocusName = strdup(ccp);
        return ret;
    }
    void print(linkage_locus_rec *p) {
        printf("\nL %d %s %d %d, ", p->locus_link, p->LocusName, p->Type, p->Class);
        printf("%d %d, ", p->number, p->col_num);
        printf("%d\n", p->AlleleCnt);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS locus_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_locus_table IF NOT EXISTS on locus_table (UniqueID);");
    }

    int db_getall(linkage_locus_rec *t);
};
#endif

extern Locus_table locus_table;

class Allele_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Allele_table()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS allele_table (pId INTEGER PRIMARY KEY,"
            " AlleleName TEXT, Frequency DOUBLE, indexX INTEGER, locus_link INTEGER"
                        ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO allele_table("
            " AlleleName, Frequency, indexX, locus_link"
            ")   VALUES(?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " AlleleName, Frequency, indexX, locus_link"
            "   FROM allele_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_allele_rec *p) {
        int idx = 1;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->AlleleName ? p->AlleleName : "")
            && insert_stmt->rowbind(idx, p->Frequency, p->index)
            && insert_stmt->rowbind(idx, p->locus_link)

            && insert_stmt->step();
    }
    int select(linkage_allele_rec *p) {
        int idx = 0;
        int ret = select_stmt->row(idx, p->AlleleName)
            && select_stmt->row(idx, p->Frequency, p->index)
            && select_stmt->row(idx, p->locus_link);
        if (*(p->AlleleName) == 0) p->AlleleName = 0;
        return ret;
    }
    void print(linkage_allele_rec *p) {
        printf("A %d %s %f %d\n", p->locus_link, p->AlleleName, p->Frequency, p->index);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS allele_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_allele_table IF NOT EXISTS on allele_table (UniqueID);");
    }

    int db_getall(linkage_allele_rec *t);
};

extern Allele_table allele_table;

class Marker_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    Marker_table()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS marker_table (pId INTEGER PRIMARY KEY,"
            " MarkerName TEXT,"
            " Recoded INTEGER, NumAlleles  INTEGER, SelectOpt INTEGER, EstFreq INTEGER,"
            " pos_avg DOUBLE, pos_male DOUBLE, pos_female DOUBLE, error_prob DOUBLE,"
            " chromosome INTEGER, col_num INTEGER, locus_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO marker_table("
            " MarkerName,"
            " Recoded, NumAlleles, SelectOpt, EstFreq,"
            " pos_avg, pos_male, pos_female, error_prob,"
            " chromosome, col_num, locus_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?,   ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " MarkerName,"
            " Recoded, NumAlleles, SelectOpt, EstFreq,"
            " pos_avg, pos_male, pos_female, error_prob,"
            " chromosome, col_num, locus_link"
            "   FROM marker_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(marker_rec *p) {
        int idx = 1;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->MarkerName)
            && insert_stmt->rowbind(idx, p->Props.Numbered.Recoded, p->Props.Numbered.NumAlleles)
            && insert_stmt->rowbind(idx, p->Props.Numbered.SelectOpt, p->Props.Numbered.EstimateFrequencies)
            && insert_stmt->rowbind(idx, p->pos_avg, p->pos_male, p->pos_female, p->error_prob)
            && insert_stmt->rowbind(idx, p->chromosome, p->col_num, p->locus_link)

            && insert_stmt->step();
    }
    int select(marker_rec *p) {
        int idx = 0;
        return select_stmt->row(idx, p->MarkerName)
            && select_stmt->row(idx, p->Props.Numbered.Recoded, p->Props.Numbered.NumAlleles)
            && select_stmt->row(idx, p->Props.Numbered.SelectOpt, p->Props.Numbered.EstimateFrequencies)
            && select_stmt->row(idx, p->pos_avg, p->pos_male, p->pos_female, p->error_prob)
            && select_stmt->row(idx, p->chromosome, p->col_num, p->locus_link);
    }
    void print(marker_rec *p) {
        printf("M %d %s, ", p->locus_link, p->MarkerName);
        printf("%d %d %d %d,",
               p->Props.Numbered.Recoded, p->Props.Numbered.NumAlleles,
               p->Props.Numbered.SelectOpt, p->Props.Numbered.EstimateFrequencies);
        printf("%f %f %f %f, ",  p->pos_avg, p->pos_male, p->pos_female, p->error_prob);
        printf("%d %d\n", p->chromosome, p->col_num);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS marker_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_marker_table IF NOT EXISTS on marker_table (UniqueID);");
    }

    int db_getall(marker_rec *p);
};

extern Marker_table marker_table;


class TraitAff_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    TraitAff_table()  {}
    int create() {
/* int *Labels
typedef struct _linkage_affection_class {
    unsigned char  MaleDef,  FemaleDef,  AutoDef;
    double *MalePen, *FemalePen, *AutoPen;
} linkage_affection_class;
 */
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS traitaff_table (pId INTEGER PRIMARY KEY,"
            " ClassCnt INTEGER, PenCnt INTEGER, NumLabels INTEGER,"
            " Labels BLOB, locus_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO traitaff_table("
            " ClassCnt, PenCnt, NumLabels,"
            " Labels, locus_link"
            ")   VALUES(?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " ClassCnt, PenCnt, NumLabels,"
            " Labels, locus_link"
            "   FROM traitaff_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_affection_data *p) {
        int idx = 1;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->ClassCnt, p->PenCnt, p->NumLabels)
            && insert_stmt->bind(idx++, p->Labels, p->NumLabels * sizeof(int))
            && insert_stmt->rowbind(idx, p->locus_link)

            && insert_stmt->step();
    }
    int select(linkage_affection_data *p) {
        int idx = 0;
        int sz = 0;
        const void *v = 0;
        int ret = select_stmt->row(idx, p->ClassCnt, p->PenCnt, p->NumLabels)
            && select_stmt->column(idx++, v, sz)
            && select_stmt->row(idx, p->locus_link);
        p->Labels = (int *)v;
        return ret;
    }
    void print(linkage_affection_data *p) {
        printf("T %d %d %d %d\n", p->locus_link, p->ClassCnt, p->PenCnt, p->NumLabels);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS traitaff_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_trait_table IF NOT EXISTS on traitaff_table (Name);");
    }

    int db_getall(linkage_affection_data *p);

};

extern TraitAff_table traitaff_table;

class AffectClass_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    AffectClass_table()  {}
    int create() {
/* int *Labels
typedef struct _linkage_affection_class {
    unsigned char  MaleDef,  FemaleDef,  AutoDef;
    double *MalePen, *FemalePen, *AutoPen;
} linkage_affection_class;
 */
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS affectclass_table (pId INTEGER PRIMARY KEY,"
            " MaleDef INTEGER, FemaleDef INTEGER, AutoDef INTEGER,"
            " MalePen BLOB, FemalePen BLOB, AutoPen BLOB,"
            " locus_link INTEGER, class_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO affectclass_table("
            " MaleDef, FemaleDef, AutoDef,"
            " MalePen, FemalePen, AutoPen,"
            " locus_link, class_link"
            ")   VALUES(?, ?, ?, ?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " MaleDef, FemaleDef, AutoDef,"
            " MalePen, FemalePen, AutoPen,"
            " locus_link, class_link"
            "   FROM affectclass_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_affection_class *p, int cnt, int i, int j) {
        int idx = 1;
        int bytes = cnt * sizeof(double);
        void *vm = p->MalePen;
        void *vf = p->FemalePen;
        void *va = p->AutoPen;
        return insert_stmt
            && insert_stmt->rowbind(idx, p->MaleDef, p->FemaleDef, p->AutoDef)

            && insert_stmt->bind(idx++, vm, bytes-sizeof(double))
            && insert_stmt->bind(idx++, vf, bytes)
            && insert_stmt->bind(idx++, va, bytes)

            && insert_stmt->rowbind(idx, i, j)

            && insert_stmt->step();
    }
    int select(linkage_affection_class *p, int &i, int &j, int &cm, int &cf, int &ca) {
        int idx = 0;
        int dm = 0, df = 0, da = 0;
        const void *vm = 0;
        const void *vf = 0;
        const void *va = 0;
        int ret = select_stmt->row(idx, dm, df, da)
            && select_stmt->column(idx++, vm, cm)
            && select_stmt->column(idx++, vf, cf)
            && select_stmt->column(idx++, va, ca)
            && select_stmt->row(idx, i, j);
        p->MaleDef = (unsigned char)dm;
        p->MalePen = (double *)vm;
        p->FemaleDef = (unsigned char)df;
        p->FemalePen = (double *)vf;
        p->AutoDef = (unsigned char)da;
        p->AutoPen = (double *)va;
        return ret;
    }
    void print(linkage_affection_class *p, int cnt, int i, int j) {
        printf("AC %2d %2d %d %d %d\n", i, j, p->MaleDef, p->FemaleDef, p->AutoDef);
        printf("P %d %d %f %f %f\n", i, cnt, p->MalePen[0], p->FemalePen[0], p->AutoPen[0]);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS affectclass_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_trait_table IF NOT EXISTS on affectclass_table (Name);");
    }

    int db_getall(linkage_affection_class *p);
};

extern AffectClass_table affectclass_table;

class ClassPen_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
public:
    ClassPen_table()  {}
    int create() {
/* int *Labels
typedef struct _linkage_penection_class {
    unsigned char  MaleDef,  FemaleDef,  AutoDef;
    double *MalePen, *FemalePen, *AutoPen;
} linkage_penection_class;
 */
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS classpen_table (pId INTEGER PRIMARY KEY,"
            " MalePen BLOB, FemalePen BLOB, AutoPen BLOB,"
            " locus_link INTEGER, class_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO classpen_table("
            " MalePen, FemalePen, AutoPen,"
            " locus_link, class_link"
            ")   VALUES(?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " MalePen, FemalePen, AutoPen,"
            " locus_link, class_link"
            "   FROM classpen_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_affection_class *p, int i, int j, int cnt) {
        int idx = 1;
        int bytes = cnt * sizeof(double);
        void *vm = p->MalePen;
        void *vf = p->FemalePen;
        void *va = p->AutoPen;
        return insert_stmt
            && insert_stmt->bind(idx++, vm, bytes-sizeof(double))
            && insert_stmt->bind(idx++, vf, bytes)
            && insert_stmt->bind(idx++, va, bytes)
            && insert_stmt->rowbind(idx, i, j)

            && insert_stmt->step();
    }
    int select(linkage_affection_class *p, int &i, int &cnt) {
/*
        int idx = 0;
        return select_stmt->column(idx++, p->MalePen, cnt)
            && select_stmt->column(idx++, p->FemalePen, cnt)
            && select_stmt->column(idx++, p->AutoPen, cnt)
            && select_stmt->row(idx, p->locus_link, i);
*/
        return 0;
    }
    void print(linkage_affection_class *p, int i, int cnt) {
        printf("P %d %d %f %f %f\n", i, cnt, p->MalePen[0], p->FemalePen[0], p->AutoPen[0]);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS classpen_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_trait_table IF NOT EXISTS on classpen_table (Name);");
    }

    int db_getall();
};

extern ClassPen_table classpen_table;

class TraitQuant_table {
    DBstmt *insert_stmt;
    DBstmt *select_stmt;
// Mean, Variance, Multiplier,
public:
    TraitQuant_table()  {}
    int create() {
        return MasterDB.exec(
            "CREATE TABLE IF NOT EXISTS traitquant_table (pId INTEGER PRIMARY KEY,"
            " ClassCnt INTEGER, Mean Blob, Variance DOUBLE, Multiplier DOUBLE,"
            " locus_link INTEGER"
            ");"
            );
    }
    int init () {
        insert_stmt = MasterDB.prep(
            "INSERT INTO traitquant_table("
            " ClassCnt, Mean, Variance, Multiplier,"
            " locus_link"
            ")   VALUES(?, ?, ?, ?, ?);");
        select_stmt = MasterDB.prep(
            "SELECT "
            " ClassCnt, Mean, Variance, Multiplier,"
            " locus_link"
            "   FROM traitquant_table;");
        return insert_stmt && select_stmt && 1;
    }
    int insert(linkage_quant_data *p); /* {
        int idx = 1;

        int sz = p->ClassCnt * 3 * sizeof(double);
        unsigned char *m = CALLOC((size_t) sz, unsigned char);

        int mi = 0;
        for (int i = 0; i < p->ClassCnt; i++)
            for (int j = 0; j < 3; j++)
                ((double *)m)[i++] = p->Mean[i][j];

        return insert_stmt
            && insert_stmt->rowbind(idx, p->ClassCnt)
            && insert_stmt->bind(idx++, m, sz)
            && insert_stmt->rowbind(idx, p->Variance[0][0], p->Multiplier)
            && insert_stmt->rowbind(idx, p->locus_link)

            && insert_stmt->step();
    }
*/
    int select(linkage_quant_data *p); /* {
        asm("int $3");
        int idx = 0;
        int sz = 0;
        double v = 0;
        unsigned char *m = 0;
        return select_stmt->row(idx, p->ClassCnt)
            && select_stmt->column(idx++, m, sz)
            && select_stmt->row(idx, v, p->Multiplier)
            && select_stmt->row(idx, p->locus_link);

        int mi = 0;
        p->Mean =   CALLOC((size_t) p->ClassCnt, double *);
        for (int i = 0; i < p->ClassCnt; i++) {
            p->Mean[i] = CALLOC((size_t) 3, double);
            for (int j = 0; j < 3; j++)
                p->Mean[i][j] = ((double *)m)[i++];

        }

        p->Variance    = CALLOC((size_t) 1, double *);
        p->Variance[0] = CALLOC((size_t) 1, double);
        p->Variance[0][0] = v;
        
    }
*/
    void print(linkage_quant_data *p) {
        printf("Q %d %d, ", p->locus_link, p->ClassCnt);
//        printf("%f %f %f\n", p->Mean, p->Variance, p->Multiplier);
    }
    void close() {
        delete insert_stmt;
        delete select_stmt;
    }
    int drop() {
        return MasterDB.exec("DROP TABLE IF EXISTS traitquant_table;");
    }

    int index() {
        return MasterDB.exec("CREATE Index Idx_trait_table IF NOT EXISTS on traitquant_table (Name);");
    }

    int db_getall(linkage_quant_data *p);

};

extern TraitQuant_table traitquant_table;



extern void dblocus_export(linkage_locus_top *LTop);

extern void dblocus_import(linkage_locus_top *LTop);

#endif
