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

#include "common.h"
#include "typedefs.h"
#include "error_messages_ext.h"
#include "utils_ext.h"

#include "dblite.hh"
#include "dbpedigree.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;
map<int, linkage_ped_tree *> Pedigree_hash;
map<int, linkage_ped_tree *> Pedigree_brkloop_hash;
//map<int, linkage_ped_rec *> Person_hash;
//map<int, linkage_ped_rec *> Person_brkloop_hash;

int Pedigree_table::db_getall(linkage_ped_tree *t, map<int, linkage_ped_tree *> &pedigree_hash) {
    linkage_ped_tree *ot = t;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(t);
            pedigree_hash[t->pedigree_link] = t;
            t++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return t-ot;
}

int Person_table::db_getall(linkage_ped_rec *p) {
    linkage_ped_rec *op = p;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(p);
//          Person_hash[p->person_link] = p;

            if (p->loopbreakers) {
                printf("LB sel %p\n", p->loopbreakers);
                int *v = p->loopbreakers;
                p->loopbreakers = CALLOC((size_t) MAXLOOPMEMBERS+1, int);
                memcpy(p->loopbreakers, v, (MAXLOOPMEMBERS+1) * sizeof (int));
            }
            p++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return p-op;
}

HMapvi Geno2idx;

void dbpedigree_export(linkage_ped_top *Top) {
    int ped, per;
    linkage_ped_tree *tpedtreep;
    linkage_ped_rec  *tpersonp;

    Tod pedexp("export pedigree/person");

    MasterDB.begin();

    // for each pedigree (consulting the linkage_ped_top structure)...
    int pers = 0, perl;
    for (ped=0; ped < Top->PedCnt; ped++) {

        tpedtreep = &(Top->PedBroken[ped]);
        tpedtreep->pedigree_link  = ped;
        pedigree_brkloop_table.insert(tpedtreep);

        for (per = 0; per < Top->PedBroken[ped].EntryCnt; per++, pers++) {
            // record for the individual (a linkage_ped_rec)...
            tpersonp = &(tpedtreep->Entry[per]);
            tpersonp->pedigree_link = ped;
            tpersonp->person_link   = pers;
            tpersonp->IsTyped       = 0;
            tpersonp->Ngeno         = 0;

            Geno2idx[tpersonp->Pheno] = pers;
            person_brkloop_table.insert(tpersonp);
        }

        tpedtreep = &(Top->PedRaw[ped]);
        tpedtreep->pedigree_link  = ped;
        pedigree_table.insert(tpedtreep);

        for (per = 0; per < Top->PedRaw[ped].EntryCnt; per++) {
            // record for the individual (a linkage_ped_rec)...
            tpersonp = &(tpedtreep->Entry[per]);
            tpersonp->pedigree_link = ped;

            perl = Geno2idx[tpersonp->Pheno];
            tpersonp->person_link = perl;
            tpersonp->IsTyped     = 0;
            tpersonp->Ngeno       = 0;
            person_table.insert(tpersonp);
        }
    }

    MasterDB.commit();

    pedexp();

}


linkage_ped_rec *Persons;
linkage_ped_rec *PersonsBroken;

static int linkage_ped_rec_sort(const void *a, const void *b) {
    linkage_ped_rec *la = (linkage_ped_rec *) a;
    linkage_ped_rec *lb = (linkage_ped_rec *) b;

    if (la->pedigree_link == lb->pedigree_link) {
        return (la->ID - lb->ID);
    } else 
        return (la->pedigree_link - lb->pedigree_link);
}

void dbpedigree_import(linkage_ped_top *Top) {
    int ped, per, Entries = 0, BrkEntries = 0;

    Top->PedRaw    = new linkage_ped_tree [Top->PedCnt];
    Top->PedBroken = new linkage_ped_tree [Top->PedCnt];

    MasterDB.begin();
    pedigree_table.db_getall(Top->PedRaw, Pedigree_hash);
    pedigree_brkloop_table.db_getall(Top->PedBroken, Pedigree_brkloop_hash);

    for (ped = 0; ped < Top->PedCnt; ped ++)
        BrkEntries += (Top->PedBroken + ped)->EntryCnt;  // vs Top->PedRaw (smaller)
    PersonsBroken = new linkage_ped_rec [BrkEntries];
    person_brkloop_table.db_getall(PersonsBroken);

    for (ped = 0; ped < Top->PedCnt; ped ++)
        Entries += (Top->PedRaw + ped)->EntryCnt;  // vs Top->PedRaw (smaller)
    Persons = new linkage_ped_rec [Entries];
    person_table.db_getall(Persons);

    MasterDB.commit();

    Top->IndivCnt = Entries;

    int pedigree_link;
    qsort(PersonsBroken, BrkEntries, sizeof (linkage_ped_rec), linkage_ped_rec_sort);
    pedigree_link = -1;
    for (per = 0; per < BrkEntries; per++) {
        linkage_ped_rec  *personp = &(PersonsBroken[per]);
        if (personp->pedigree_link != pedigree_link) {
            linkage_ped_tree *pedtreep;
            if (map_get(Pedigree_brkloop_hash, personp->pedigree_link, pedtreep))
                pedtreep->Entry = personp;
            else {
                errorvf("pedigree_link from persons (%d) not found in pedigree tree\n",
                        personp->pedigree_link);
                EXIT(0);
            }
            pedigree_link = personp->pedigree_link;
        }
    }

    qsort(Persons, Entries, sizeof (linkage_ped_rec), linkage_ped_rec_sort);
    pedigree_link = -1;
    for (per = 0; per < Entries; per++) {
        linkage_ped_rec  *personp = &(Persons[per]);
        if (personp->pedigree_link != pedigree_link) {
            linkage_ped_tree *pedtreep;
            if (map_get(Pedigree_hash, personp->pedigree_link, pedtreep))
                pedtreep->Entry = personp;
            else {
                errorvf("pedigree_link from persons (%d) not found in pedigree tree\n",
                        personp->pedigree_link);
                EXIT(0);
            }
            pedigree_link = personp->pedigree_link;
        }
    }

    return;
/*
    Entries = 0;
    int EntryCnt;
    for (ped = 0; ped < Top->PedCnt; ped++) {
        pedigree_table.print(Top->Ped + ped);
        EntryCnt = (Top->Ped + ped)->EntryCnt;
        for (per = 0; per < EntryCnt; per++)
            person_table.print(Persons + Entries + per);
        Entries += EntryCnt;
    }
 */
}
