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
#include "errorno.h"
#include "utils_ext.h"

#include "dblite.hh"
#include "dbgenotype.hh"

#include "tod.hh"

using namespace std;

extern DBlite MasterDB;
extern map<int, linkage_ped_rec *> Person_hash;


typedef std::list<Pairll>                      List2pl;
typedef std::list<Pairll>::const_iterator      List2plp;
typedef std::map<int,List2pl*>                 Marker_filter;
typedef std::map<int,List2pl*>::const_iterator Mfp;

Marker_filter marker_filter;

//xx
int SHOW = 0;
void marker_filter_add(int bchr, long long bmin, long long bmax)
{
    List2pl *blist;

    if (not map_get(marker_filter, bchr, blist)) {
        blist = new List2pl;
        marker_filter[bchr] = blist;
    }
    blist->push_back(Pairll(bmin, bmax));
}

void mk_marker_filter(linkage_ped_top *Top)
{
    extern int base_pair_position_index;
    extern int genetic_distance_index;

    linkage_locus_top *LTop = Top->LocusTop;
    ext_linkage_locus_top *EXLTop = Top->EXLTop;
    int offset = LTop->PhenoCnt;
    int bchr = -1;
    int i, io = offset;
    long long bmin, bmax, bp;
//    asm("int $3");
    if (base_pair_position_index < 0)
        bmin = bmax = bp = EXLTop->EXLocus[offset].positions[genetic_distance_index];
    else
        bmin = bmax = bp = EXLTop->EXLocus[offset].positions[base_pair_position_index];

    for (i = offset; i < LTop->LocusCnt; i++) {
        if (LTop->Locus[i].Marker->chromosome != bchr && bchr != -1) {

//          marker_filter_add(bchr, bmin, bmax);
            marker_filter_add(bchr, io, i - 1);

            io = i;
            bmin = bmax = EXLTop->EXLocus[i].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];
        }
        bchr = LTop->Locus[i].Marker->chromosome;
        bp = EXLTop->EXLocus[i].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];
        if (bp < bmin) bmin = bp;
        if (bp > bmax) bmax = bp;
    }   
//  marker_filter_add(bchr, bmin, bmax);
    marker_filter_add(bchr, io, i - 1);
}

void use_marker_filter(void *mk, linkage_ped_top *Top, int link, int bchr, int bytes, unsigned char *data)
{
    linkage_locus_top *LTop = Top->LocusTop;
    int offset = LTop->PhenoCnt;

    long long bmin, bmax;
    List2pl  *blist;
    List2plp lp;
    Pairll   pp;

    SHOW = 0;
    if (not map_get(marker_filter, bchr, blist)) {
        if (bchr == -1)
            printf("Marker_filter: person_link %d no data\n", link);
        else
            printf("Marker_filter: person_link %d ignoring chromosome %d\n", link, bchr);
        return;
    }
//xx
    if (SHOW && bchr == 5) {
        printf("bchr %d, j %d\n", bchr, bytes);
        for (int ii = 0; ii < bytes; ii++) {
            unsigned char uc = *(((unsigned char *)data) + ii);
            printf("%x%x ", (uc>>4)&0xf, uc&0xf);
        }
        printf("\n");
        asm("int $3");
    }
//xx
    for (lp = blist->begin(); lp != blist->end(); lp++) {
        pp = *lp;
        bmin = pp.first;
        bmax = pp.second;

        for (long long i = bmin, dataj = offset; i <= bmax; i++, dataj++)
            copy_2alleles_2staging(mk, marker_start(data, -offset), i, dataj);
    }

//xx
    if (SHOW && bchr == 5) {
        printf("bchr %d, j %d\n", bchr, bytes);
        for (int ii = 2*114, j = 0; j < bytes; j++) {
            unsigned char uc = *(((unsigned char *)mk) + ii + j);
            printf("%x%x ", (uc>>4)&0xf, uc&0xf);
        }
        printf("\n");
        asm("int $3");
    }
//xx
    SHOW = 0;
}


int Phenotype_table::db_getall(linkage_locus_top *LTop, pheno_pedrec_data **Phenotypes) {
    int knt = 0;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
//            ret = select(p);
            int link = 0, bytes = 0;
            unsigned char *data = (unsigned char *)0;
            ret = select(link, bytes, data);
            Phenotypes[link] = CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data);
            memcpy(Phenotypes[link], data, bytes);
            knt++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return knt;
}

int Genotype_table::db_getall(linkage_ped_top *Top, void **Genotypes) {
    int knt = 0;
    linkage_locus_top *LTop = Top->LocusTop;

    int ret = select_stmt && select_stmt->abort();
    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
            int link = 0, bytes = 0, chr = -1;
            unsigned char *data = (unsigned char*)0;
//          asm("int $3");
            ret = select(link, chr, bytes, data);
            if (! Genotypes[link]) {
                if (data) {
                    Genotypes[link] = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
//                  Genotypes[link] = (void *)CALLOC((size_t) bytes, unsigned char);
//                  memcpy(Genotypes[link], data, bytes);
                }
            }
            use_marker_filter(Genotypes[link], Top, link, chr, bytes, data);

            knt++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    return knt;
}

void dbgenotype_export(linkage_ped_top *Top, bp_order *bp) {
    linkage_locus_top *LTop = Top->LocusTop;
    int ped, per;
    linkage_ped_tree *tpedtreep;
    linkage_ped_rec  *tpersonp;

    Tod pedexp("export genotype/phenotype");

//  asm("int $3");
    // for each genotype (consulting the linkage_ped_top structure)...
    void *mk = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
    void *mks = marker_start(mk, LTop->PhenoCnt);
    int   sz = marker_size(LTop->MarkerCnt);
    void *sv;
    bp_order *b;
    int bchr;

    MasterDB.begin();

    int pers = 0;
    int i, j = 0;
    int offset = LTop->PhenoCnt;
    for (ped=0; ped < Top->PedCnt; ped++) {
        tpedtreep = &(Top->PedBroken[ped]);  // vs Top->PedRaw

        for (per = 0; per < Top->PedBroken[ped].EntryCnt; per++, pers++) {
            // record for the individual phenotype and genotype
            tpersonp = &(tpedtreep->Entry[per]);

//          if (tpersonp->person_link == 259) asm("int $3");
            phenotype_table.insert(tpersonp, LTop->PhenoCnt);
            sv = tpersonp->Marker;
            memset(mks, 0, sz);
            tpersonp->Marker = mk;
            bchr = -1;
            for (i = offset, j = 0,
                     b = bp + offset; i < LTop->LocusCnt; i++, b++) {

                if (b->chr != bchr && bchr != -1) {  // time to write a chr
//                  asm("int $3"); 
                    if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
                        for (; (j % 4); j++) 
                            set_2alleles(mk, offset + j, NULL, 0, 0); // NULL is not used iff alleles are 0/0 
                    }
                    genotype_table.insert(mks, tpersonp->person_link, bchr, j);
//xx
                    if (SHOW && bchr == 5) {
                        printf("bchr %d, j %d\n", bchr, j);
                        for (int ii = 0; ii < 2*j; ii++) {
                            unsigned char uc = *(((unsigned char *)mks) + ii);
                            printf("%x%x ", (uc>>4)&0xf, uc&0xf);
                        }
                        printf("\n");
                        asm("int $3");
                    }
//xx
                    j = 0;
                }


                if (copy_2alleles_2staging(mk, sv, offset + j, b->i)) {
                    tpersonp->Marker = sv;
                    break;
                }
                bchr = b->chr;
                j++;
            }
//?? if sv is NULL; the above loop is aborted and we end up here
//?? PS we don't know length of each chr run;

            for (; (j % 4); j++) 
                set_2alleles(mk, offset + j, NULL, 0, 0);
            if (tpersonp->Marker == sv)
                genotype_table.insert(NULL, tpersonp->person_link, bchr, j);
            else
                genotype_table.insert(mks, tpersonp->person_link, bchr, j);

            tpersonp->Marker = sv;
        }
    }

    MasterDB.commit();
    marker_free(mk, LTop->PhenoCnt);
    pedexp();

}

pheno_pedrec_data **Phenotypes;
void **Genotypes;

void dbgenotype_import_phenotype(linkage_ped_top *Top) {
    int ped, per;

    MasterDB.begin();
    phenotype_table.db_getall(Top->LocusTop, Phenotypes);
    MasterDB.commit();

    for (ped = 0; ped < Top->PedCnt; ped ++) {
        for (per = 0; per < (Top->PedBroken + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedBroken[ped].Entry[per].person_link;
            Top->PedBroken[ped].Entry[per].Pheno = Phenotypes[person_link];
        }
        for (per = 0; per < (Top->PedRaw + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedRaw[ped].Entry[per].person_link;
            Top->PedRaw[ped].Entry[per].Pheno = Phenotypes[person_link];
        }
    }
}

void dbgenotype_import_genotype(linkage_ped_top *Top) {
    extern void db_open_db();

    char *buf, *bufp;
    int ped, per;

    buf  = CALLOC((size_t) 4 * main_chromocnt + 9 /* chr in () */ + 1 /* 999 */+ 1 /* \0 */, char);
    bufp = buf;
    sprintf(bufp, "chr in (%d", global_chromo_entries[0]);
    bufp += strlen(bufp);
    for (int i = 1; i < main_chromocnt; i++) {
        sprintf(bufp, ", %d", global_chromo_entries[i]);
        bufp += strlen(bufp);
    }
    sprintf(bufp, ")");
    warnvf("Select * from genotype_table where %s\n", buf);
                 
    db_open_db();
    genotype_table.init(buf);
    free(buf);

    MasterDB.begin();
    genotype_table.db_getall(Top, Genotypes);
    MasterDB.commit();

    for (ped = 0; ped < Top->PedCnt; ped ++) {
        for (per = 0; per < (Top->PedBroken + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedBroken[ped].Entry[per].person_link;
//?x
            Top->PedBroken[ped].Entry[per].Marker = marker_start(Genotypes[person_link],
                                                            Top->LocusTop->PhenoCnt
                                                           -Top->LocusTop->PhenoCnt);
        }
        for (per = 0; per < (Top->PedRaw + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedRaw[ped].Entry[per].person_link;
//?x
            Top->PedRaw[ped].Entry[per].Marker = marker_start(Genotypes[person_link],
                                                            Top->LocusTop->PhenoCnt
                                                           -Top->LocusTop->PhenoCnt);
        }
    }

    genotype_table.close();
    MasterDB.close();
}

void dbgenotype_import(linkage_ped_top *Top) {
    int ped, Entries = 0;

    for (ped = 0; ped < Top->PedCnt; ped ++)
        Entries += (Top->PedBroken + ped)->EntryCnt;  // vs Top->PedRaw

    Phenotypes = new pheno_pedrec_data * [Entries];
    Genotypes  = new void *[Entries];

    for (int i = 0; i < Entries; i++) {
        Phenotypes[i] = NULL;
        Genotypes[i] = NULL;
    }

    dbgenotype_import_phenotype(Top);
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
