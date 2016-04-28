/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

#include "common.h"
#include "typedefs.h"
#include "error_messages_ext.h"
#include "utils_ext.h"

#include "dblite.hh"
#include "dbgenotype.hh"

#include "Tod.hh"

using namespace std;

extern DBlite MasterDB;
extern map<int, linkage_ped_rec *> Person_hash;

int Phenotype_table::db_getall(linkage_locus_top *LTop, pheno_pedrec_data **Phenotypes) {
    int cnt = 0;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
//            ret = select(p);
            int link = 0, cnt = 0, bytes = 0;
            unsigned char *data = (unsigned char *)0;
            ret = select(link, cnt, bytes, data);
            Phenotypes[link] = CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data);
            memcpy(Phenotypes[link], data, bytes);
            cnt++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return cnt;
}

int Genotype_table::db_getall(linkage_locus_top *LTop, void **Genotypes) {
    int cnt = 0;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            int link = 0, cnt = 0, bytes = 0;
            unsigned char *data = (unsigned char*)0;
            ret = select(link, cnt, bytes, data);
            if (data) {
//              Genotypes[link] = marker_alloc((size_t) LTop->MarkerCnt, 0);
                Genotypes[link] = (void *)CALLOC((size_t) bytes, unsigned char);
                memcpy(Genotypes[link], data, bytes);
            } else {
                Genotypes[link] = 0;
            }
            cnt++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return cnt;
}

void dbgenotype_export(linkage_ped_top *Top) {
    int ped, per;
    linkage_ped_tree *tpedtreep;
    linkage_ped_rec  *tpersonp;

    Tod pedexp("export genotype/phenotype");

    MasterDB.begin();

    // for each genotype (consulting the linkage_ped_top structure)...
    int pers = 0;
    for (ped=0; ped < Top->PedCnt; ped++) {
        tpedtreep = &(Top->Ped[ped]);

        for (per = 0; per < Top->Ped[ped].EntryCnt; per++, pers++) {
            // record for the individual phenotype and genotype
            tpersonp = &(tpedtreep->Entry[per]);

            phenotype_table.insert(tpersonp, Top->LocusTop->PhenoCnt);
            genotype_table.insert(tpersonp, Top->LocusTop->MarkerCnt, Top->LocusTop->PhenoCnt);
        }
    }

    MasterDB.commit();

    pedexp();

}

pheno_pedrec_data **Phenotypes;
void **Genotypes;

void dbgenotype_import(linkage_ped_top *Top) {
    int ped, per, cnt = 0, Entries = 0;


    for (ped = 0; ped < Top->PedCnt; ped ++)
        Entries += (Top->Ped + ped)->EntryCnt;

    Phenotypes = new pheno_pedrec_data * [Entries];
    Genotypes  = new void *[Entries];

    MasterDB.begin();

    phenotype_table.db_getall(Top->LocusTop, Phenotypes);
    genotype_table.db_getall(Top->LocusTop, Genotypes);

    MasterDB.commit();

    for (ped = 0; ped < Top->PedCnt; ped ++) {
        for (per = 0; per < (Top->Ped + ped)->EntryCnt; per++) {
            Top->Ped[ped].Entry[per].Pheno = Phenotypes[cnt];
            Top->Ped[ped].Entry[per].Marker = marker_start(Genotypes[cnt], -Top->LocusTop->PhenoCnt);
            cnt++;
        }
    }

    return;
/*
    for (ped = 0; ped < Top->PedCnt; ped++) {
        int EntryCnt = (Top->Ped + ped)->EntryCnt;
        for (per = 0; per < EntryCnt; per++) {
            phenotype_table.print((Top->Ped + ped)->Entry + per, Top);
            genotype_table.print((Top->Ped + ped)->Entry + per, Top);
        }
    }
 */
}
