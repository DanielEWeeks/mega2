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
#include "read_files_ext.h"
#include "error_messages_ext.h"
#include "utils_ext.h"

#include "dblite.hh"
#include "dbloci.hh"
#include "dbmisc.hh"

#include "Tod.hh"

using namespace std;

extern DBlite MasterDB;
map<int, linkage_locus_rec *> Locus_hash;
map<pair<int,int>, linkage_affection_class *>AffectClass_hash;

int Locus_table::db_getall(linkage_locus_rec *t) {
    linkage_locus_rec *ot = t;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(t);
            Locus_hash[t->locus_link] = t;
            t++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return t-ot;
}

int Allele_table::db_getall(linkage_allele_rec *t) {
    linkage_allele_rec *ot = t;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(t);
            if (t->AlleleName)
                t->AlleleName = canonical_allele(t->AlleleName);
            t++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return t-ot;
}

int Marker_table::db_getall(marker_rec *p) {
    marker_rec *op = p;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(p);
            p->MarkerName = strdup(p->MarkerName);
            p++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return p-op;
}

int TraitAff_table::db_getall(linkage_affection_data *p) {
    linkage_affection_data *op = p;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(p);

            if (p->Labels) {
                printf("LBL sel %p\n", p->Labels);
                int *v = p->Labels;
                p->Labels = CALLOC((size_t) p->NumLabels, int);
                memcpy(p->Labels, v, p->NumLabels * sizeof (int));

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

int  AffectClass_table::db_getall(linkage_affection_class *p) {
    linkage_affection_class *op = p;
    int l_link, c_link, cm, cf, ca;
    double *ds;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(p, l_link, c_link, cm, cf, ca);

            ds = p->MalePen;
            p->MalePen = (double *)CALLOC((size_t)cm, unsigned char);
            memcpy(p->MalePen, ds, cm);

            ds = p->FemalePen;
            p->FemalePen = (double *)CALLOC((size_t)cf, unsigned char);
            memcpy(p->FemalePen, ds, cf);

            ds = p->AutoPen;
            p->AutoPen = (double *)CALLOC((size_t)ca, unsigned char);
            memcpy(p->AutoPen, ds, ca);

            AffectClass_hash[make_pair(l_link, c_link)] = p;

            p++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }

    return p - op;
}

int TraitQuant_table::insert(linkage_quant_data *p) {
    int idx = 1;

    int sz = p->ClassCnt * 3 * sizeof(double);
    unsigned char *m = CALLOC((size_t) sz, unsigned char);

    int mi = 0;
    for (int i = 0; i < p->ClassCnt; i++)
        for (int j = 0; j < 3; j++)
            ((double *)m)[mi++] = p->Mean[i][j];

    return insert_stmt
        && insert_stmt->rowbind(idx, p->ClassCnt)
        && insert_stmt->bind(idx++, m, sz)
        && insert_stmt->rowbind(idx, p->Variance[0][0], p->Multiplier)
        && insert_stmt->rowbind(idx, p->locus_link)

        && insert_stmt->step();
}

int TraitQuant_table::select(linkage_quant_data *p) {
    int idx = 0;

    int sz = 0;
    double v = 0;
    const void *m = 0;
    int ret = select_stmt->row(idx, p->ClassCnt)
        && select_stmt->column(idx++, m, sz)
        && select_stmt->row(idx, v, p->Multiplier)
        && select_stmt->row(idx, p->locus_link);

    int mi = 0;
    p->Mean =   CALLOC((size_t) p->ClassCnt, double *);
    for (int i = 0; i < p->ClassCnt; i++) {
        p->Mean[i] = CALLOC((size_t) 3, double);
        for (int j = 0; j < 3; j++)
            p->Mean[i][j] = ((double *)m)[mi++];

    }

    p->Variance    = CALLOC((size_t) 1, double *);
    p->Variance[0] = CALLOC((size_t) 1, double);
    p->Variance[0][0] = v;

    return ret;
}

int TraitQuant_table::db_getall(linkage_quant_data *p) {
    linkage_quant_data *op = p;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            ret = select(p++);
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return p-op;
}

linkage_allele_rec *Allele_rec;
linkage_affection_data *Affect_rec;
linkage_affection_class *AffectClass_rec;
linkage_quant_data *Quant_rec;
//tmp
int AllelesCnt = 0;
int AffectClassCnt = 0;

void dblocus_export(linkage_locus_top *Top) {
    Tod pedexp("export loci/marker");

    int LocusCnt  = Top->LocusCnt;

    MasterDB.begin();

    linkage_locus_rec *p = Top->Locus;
    for (int i = 0; i < LocusCnt; i++, p++) {
        p->locus_link = i;
        locus_table.insert(p);
        linkage_allele_rec *a = p->Allele;
        for (int j = 0; j < p->AlleleCnt; j++) {
            a->locus_link = i;
            allele_table.insert(a++);
            AllelesCnt++;
        }
        if (p->Type == NUMBERED || p->Type == XLINKED || p->Type == YLINKED) {
            p->Marker->locus_link = i;
            marker_table.insert(p->Marker);
        } else if (p->Type == AFFECTION) {
            p->Pheno->Props.Affection.locus_link = i;
            traitaff_table.insert(&(p->Pheno->Props.Affection));
            for (int j = 0; j < p->Pheno->Props.Affection.ClassCnt; j++) {
                affectclass_table.insert((p->Pheno->Props.Affection.Class + j),
                                         p->Pheno->Props.Affection.PenCnt,
                                         i, j);
                AffectClassCnt++;
//                classpen_table.insert(&(p->Pheno->Props.Affection.Class + j), i, j)
            }
        } else if (p->Type == QUANT) {
            p->Pheno->Props.Quant.locus_link = i;
            traitquant_table.insert(&(p->Pheno->Props.Quant));
        }
    }

    int_table.insert("AllelesCnt", AllelesCnt);
    int_table.insert("AffectClassCnt", AffectClassCnt);

    MasterDB.commit();

    pedexp();
}


static int linkage_allele_rec_sort(const void *a, const void *b) {
    linkage_allele_rec *la = (linkage_allele_rec *) a;
    linkage_allele_rec *lb = (linkage_allele_rec *) b;

    if (la->locus_link == lb->locus_link) {
        return (la->index - lb->index);
    } else 
        return (la->locus_link - lb->locus_link);
}

static int marker_rec_sort(const void *a, const void *b) {
    marker_rec *ma = (marker_rec *) a;
    marker_rec *mb = (marker_rec *) b;

    return (ma->locus_link - mb->locus_link);
}

static int trait_affect_sort(const void *a, const void *b) {
    linkage_affection_data *ta = (linkage_affection_data *) a;
    linkage_affection_data *tb = (linkage_affection_data *) b;

    return (ta->locus_link - tb->locus_link);
}

static int trait_quant_sort(const void *a, const void *b) {
    linkage_quant_data *qa = (linkage_quant_data *) a;
    linkage_quant_data *qb = (linkage_quant_data *) b;

    return (qa->locus_link - qb->locus_link);
}

void dblocus_import(linkage_locus_top *LTop) {
    int a, m;
    int locus_link = -1;
//  int LocusCnt  = LTop->LocusCnt;
    int MarkerCnt = LTop->MarkerCnt;
    int PhenoCnt  = LTop->PhenoCnt;
    int PhenoAffectCnt, PhenoQuantCnt;

    LTop->Locus = new linkage_locus_rec [LTop->LocusCnt];
    Allele_rec = new linkage_allele_rec [AllelesCnt];
    LTop->Marker = new marker_rec [MarkerCnt];
    LTop->Pheno = new pheno_rec [PhenoCnt];
    Affect_rec = new linkage_affection_data [PhenoCnt];
    AffectClass_rec = new linkage_affection_class [AffectClassCnt];
    Quant_rec = new linkage_quant_data [PhenoCnt];

// LTop parameters
//  asm("int $3");

    extern map<const char *, unsigned char *, charsless> Stuff_hash;

    const char *k;
    unsigned char *d;
    k = "MaleRecomb";
    if (!map_get(Stuff_hash, k, d)) {
        printf("Stuff read failed for %s\n", k);
    }
    LTop->MaleRecomb = (double *)d;

    k = "FemaleRecomb";
    if (!map_get(Stuff_hash, k, d)) {
        printf("Stuff read failed for %s\n", k);
    }
    LTop->FemaleRecomb = (double *)d;

    k = "recomb_frac";
    if (!map_get(Stuff_hash, k, d)) {
        printf("Stuff read failed for %s\n", k);
    }
    LTop->Run.linkmap.recomb_frac = (double *)d;


    MasterDB.begin();

    locus_table.db_getall(LTop->Locus);
    allele_table.db_getall(Allele_rec);
    marker_table.db_getall(LTop->Marker);
    PhenoAffectCnt = traitaff_table.db_getall(Affect_rec);
    PhenoQuantCnt = traitquant_table.db_getall(Quant_rec);

    MasterDB.commit();

    qsort(Allele_rec, AllelesCnt, sizeof (linkage_allele_rec), linkage_allele_rec_sort);
    locus_link = -1;
    for (a = 0; a < AllelesCnt; a++) {
        linkage_allele_rec  *ap = &(Allele_rec[a]);
        if (ap->locus_link != locus_link) {
            linkage_locus_rec *rp;
            if (map_get(Locus_hash, ap->locus_link, rp))
                rp->Allele = ap;
            else {
                errorvf("locus_link from allele (%d) not found in locus\n",
                        ap->locus_link);
                EXIT(0);
            }
            locus_link = ap->locus_link;
        }
    }

    qsort(LTop->Marker, MarkerCnt, sizeof (marker_rec), marker_rec_sort);
    for (m = 0; m < MarkerCnt; m++) {
        marker_rec  *mp = &(LTop->Marker[m]);
        linkage_locus_rec *lp;
        if (map_get(Locus_hash, mp->locus_link, lp))
            lp->Marker = mp;
        else {
            errorvf("locus_link from marker_rec (%d) not found in locus\n",
                    mp->locus_link);
            EXIT(0);
        }
    }
    LTop->Marker -= PhenoCnt;

    affectclass_table.db_getall(AffectClass_rec);

    qsort(Affect_rec, PhenoAffectCnt, sizeof (linkage_affection_data), trait_affect_sort);
    pair<int,int> par;
    for (m = 0; m < PhenoAffectCnt; m++) {
        linkage_affection_data  *pp = &(Affect_rec[m]);
        linkage_affection_class *acp = 0;
        linkage_locus_rec *lp;
        par.first = pp->locus_link;
        if (map_get(Locus_hash, pp->locus_link, lp)) {
            lp->Pheno = &LTop->Pheno[pp->locus_link];  // by construction pp->locus_link < PhenoCnt
            lp->Pheno->TraitName = lp->LocusName;
            lp->Pheno->Props.Affection = *pp;
            lp->Pheno->Props.Affection.Class = CALLOC((size_t) lp->Pheno->Props.Affection.ClassCnt,
                                                      linkage_affection_class);
            for (int j = 0; j < lp->Pheno->Props.Affection.ClassCnt; j++) {
                par.second = j;
                if (map_get(AffectClass_hash, par, acp)) {
                    lp->Pheno->Props.Affection.Class[j] = *acp;
                } else {
                    printf("acp: %d %d .. none\n", pp->locus_link, j);
                }
            }
        } else {
            errorvf("locus_link from affect trait (%d) not found in locus\n",
                    pp->locus_link);
            EXIT(0);
        }
    }
    qsort(Quant_rec, PhenoQuantCnt, sizeof (linkage_quant_data), trait_quant_sort);
    for (m = 0; m < PhenoQuantCnt; m++) {
        linkage_quant_data  *pp = &(Quant_rec[m]);
        linkage_locus_rec *lp;
        if (map_get(Locus_hash, pp->locus_link, lp)) {
            lp->Pheno = &LTop->Pheno[pp->locus_link];  // by construction pp->locus_link < PhenoCnt
            lp->Pheno->TraitName = lp->LocusName;
            lp->Pheno->Props.Quant = *pp;
        } else {
            errorvf("locus_link from marker_rec (%d) not found in locus\n",
                    pp->locus_link);
            EXIT(0);
        }
    }

    delete [] Affect_rec;
    delete [] Quant_rec;

    return;
/*
    for (int l = 0; l < LocusCnt; l++) {
        linkage_locus_rec *p = LTop->Locus + l;
        locus_table.print(p);
        for (a = 0; a < p->AlleleCnt; a++) {
            allele_table.print(p->Allele + a);
        }
        if (p->Type == NUMBERED || p->Type == XLINKED || p->Type == YLINKED) {
            marker_table.print(p->Marker);
        } else if (p->Type == AFFECTION) {
            traitaff_table.print(&(p->Pheno->Props.Affection));
        } else if (p->Type == QUANT) {
            traitquant_table.print(&(p->Pheno->Props.Quant));
        }
    }
 */
}
