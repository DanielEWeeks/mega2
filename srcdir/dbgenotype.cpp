/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

extern int base_pair_position_index;
extern int genetic_distance_index;
extern map<int,Pairii> Chr2Locus;


typedef std::list<Pairll>                      List2ll;
typedef std::list<Pairll>::const_iterator      List2llp;
typedef std::map<int,List2ll*>                 Marker_filter;
typedef std::map<int,List2ll*>::const_iterator Mfp;

typedef std::pair<int,Pairll>                   Pairill;
typedef std::list<Pairill>                      List3ill;
typedef std::list<Pairill>::const_iterator      List3illp;
typedef std::map<int,List3ill*>                 Locus_filter;
typedef std::map<int,List3ill*>::const_iterator Lfp;

Marker_filter marker_filter;
Locus_filter  locus_filter;

void marker_filter_add(int bchr, long long bmin, long long bmax)
{
    List2ll *bplist;

    if (bchr != 1 && bchr != 3) return;

    if (! map_get(marker_filter, bchr, bplist)) {
        bplist = new List2ll;
        marker_filter[bchr] = bplist;
    }
//  bplist->push_back(Pairll(bmin, bmax));
    long long del = bmax-bmin;
    del /=5;
    bplist->push_back(Pairll(bmin+3*del, bmin+4*del));
    bplist->push_back(Pairll(bmin+del, bmin+2*del));
    bplist->push_back(Pairll(bmax-del, bmax));
}

void mk_marker_filter(linkage_ped_top *Top)
{
    linkage_locus_top *LTop = Top->LocusTop;
    ext_linkage_locus_top *EXLTop = Top->EXLTop;
    int offset = LTop->PhenoCnt;
    int bchr = -1;
    int i;
    long long bmin, bmax, bp;

    if (1 /*ALL*/) return;

    if (base_pair_position_index < 0)
        bmin = bmax = bp = EXLTop->EXLocus[offset].positions[genetic_distance_index];
    else
        bmin = bmax = bp = EXLTop->EXLocus[offset].positions[base_pair_position_index];

    for (i = offset; i < LTop->LocusCnt; i++) {
        if (LTop->Locus[i].Marker->chromosome != bchr && bchr != -1) {

            marker_filter_add(bchr, bmin, bmax);

            bmin = bmax = EXLTop->EXLocus[i].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];
        }
        bchr = LTop->Locus[i].Marker->chromosome;
        bp = EXLTop->EXLocus[i].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];
        if (bp < bmin) bmin = bp;
        if (bp > bmax) bmax = bp;
    }
    marker_filter_add(bchr, bmin, bmax);
}

SECTION_LOG_INIT(ignore_chromosome);

void use_marker_filter(void *mk, linkage_ped_top *Top, int link, int bchr, int bytes, unsigned char *data)
{
    SECTION_LOG_EXTERN(ignore_chromosome);

    linkage_locus_top *LTop = Top->LocusTop;
    int offset = LTop->PhenoCnt;

    long long bmin, bmax;
    List2ll  *bplist;
    List2llp lp;
    Pairll   pp;

    if (! map_get(marker_filter, bchr, bplist)) {
        SECTION_LOG(ignore_chromosome);
        if (bchr == -1)
            dbgvf("Marker_filter: person_link %d no data\n", link);
        else
            dbgvf("Marker_filter: person_link %d ignoring chromosome %d\n", link, bchr);
        return;
    }
    for (lp = bplist->begin(); lp != bplist->end(); lp++) {
        pp = *lp;
        bmin = pp.first;
        bmax = pp.second;

        for (long long i = bmin, dataj = offset; i <= bmax; i++, dataj++)
            copy_2alleles_2staging(mk, marker_start(data, -offset), i, dataj);
    }
}

void use_locus_filter(void *mk, linkage_ped_top *Top, int link, int bchr, int bytes, unsigned char *data)
{
    SECTION_LOG_EXTERN(ignore_chromosome);

    linkage_locus_top *LTop = Top->LocusTop;
    int offset = LTop->PhenoCnt;
    void *ds = marker_start(data, -offset);

    long long lmin, lmax, k;
    List3ill  *locuslist;
    List3illp  lp;
    Pairill    ill;

    if (! map_get(locus_filter, bchr, locuslist)) {
        SECTION_LOG(ignore_chromosome);
        if (bchr == -1)
            dbgvf("Locus_filter: person_link %d no data\n", link);
        else
            dbgvf("Locus_filter: person_link %d ignoring chromosome %d\n", link, bchr);
        return;
    }

    for (lp = locuslist->begin(); lp != locuslist->end(); lp++) {
        ill = *lp;
        k    = ill.first;
        lmin = ill.second.first;
        lmax = ill.second.second;
//      printf("lnk %d, chr %d, k %lld, [%lld, %lld]\n", link, bchr, k, lmin, lmax);

        for (long long i = lmin; i < lmax; i++, k++)
            copy_2alleles_2staging(mk, ds, k, i + offset);
    }
}

int Phenotype_table::db_getall(linkage_locus_top *LTop, pheno_pedrec_data **Phenotypes) {
    int knt = 0;
    int ret = select_stmt && select_stmt->abort();

    while (ret) {
        ret = select_stmt->step();
        if (ret == SQLITE_ROW) {
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

            ret = select(link, chr, bytes, data);
            if (! Genotypes[link]) {
                if (data) {
                    Genotypes[link] = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
                }
            }
            use_locus_filter(Genotypes[link], Top, link, chr, bytes, data);
            knt++;
        } else if (ret == SQLITE_DONE) {
            ret = 0;
        } else {
            ret = 0;
        }
    }
    SECTION_LOG_FINI(ignore_chromosome);
    return knt;
}

void dbgenotype_export(linkage_ped_top *Top, bp_order *bp) {
    linkage_locus_top *LTop = Top->LocusTop;
    int ped, per;
    linkage_ped_tree *tpedtreep;
    linkage_ped_rec  *tpersonp;

    Tod pedexp("export genotype/phenotype");

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

            phenotype_table.insert(tpersonp, LTop->PhenoCnt);
            sv = tpersonp->Marker;
            memset(mks, 0, sz);
            tpersonp->Marker = mk;
            bchr = -1;
            for (i = offset, j = 0,
                     b = bp + offset; i < LTop->LocusCnt; i++, b++) {

                if (b->chr != bchr && bchr != -1) {  // time to write a chr
                    if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
                        for (; (j % 4); j++) 
                            set_2alleles(mk, offset + j, NULL, 0, 0); // NULL is not used iff alleles are 0/0 
                    }
                    genotype_table.insert(mks, tpersonp->person_link, bchr, j);
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

            if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
                for (; (j % 4); j++) 
                    set_2alleles(mk, offset + j, NULL, 0, 0);
            }
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

    Tod index_genotype("index_genotype");
    genotype_table.index();
    index_genotype();
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

static int Lesslocus = 0;

static int cmp(const void *a, const void *b) { return ( (*(const int *)a) - (*(const int *)b) ); }

static int *convert_marker_filter2locus_filter(linkage_ped_top *Top)
{
    ext_linkage_locus_top *EXLTop = Top->EXLTop;
    int offset = Top->LocusTop->PhenoCnt;
    int bchr, i, chrcnt, iii;
    long long bmin, bmax, bp;

    long cmin, cmax;

    List2ll  *bplist;
    List2llp lp;
    Pairll   pp;

    List3ill *locuslist;

    int cnt, *chrs;

    bool AllChr = false; // seems to work OK; select/use chr's that are required
//  bool AllChr = true;  // seems to work OK; select all chrs, use what is required
    if (AllChr && /*ALL*/ marker_filter.size() == 0) {
        int start, len;
        List3ill *locuslist;
        Lesslocus = 0;
        long long cmin, cmax;
        for (map<int,Pairii>::const_iterator chrp = Chr2Locus.begin();
             chrp != Chr2Locus.end();
             chrp++) {
            bchr  = chrp->first;
            start = chrp->second.first;
            len   = chrp->second.second;

            Lesslocus = start;
            cmin = 0;
            cmax = len;

            if (not map_get(locus_filter, bchr, locuslist)) {
                locuslist = new List3ill;
                locus_filter[bchr] = locuslist;
            }
            locuslist->push_back(make_pair(Lesslocus, Pairll(cmin, cmax)));
        }
        Lesslocus += cmax;
        return 0;
    }
// main_chromocnt of global_chromo_entries
    if (! AllChr) {
        int len;
        List3ill *locuslist;
        Lesslocus = offset;
        chrcnt = 0;
        long long cmin, cmax;
        chrs = CALLOC((size_t) main_chromocnt, int);

        for (i = 0; i < main_chromocnt; i++)
            chrs[i] = global_chromo_entries[i];
        qsort(chrs, main_chromocnt, sizeof (int), cmp);

        for (i = 0; i < main_chromocnt; i++) {
            bchr = chrs[i];
            Pairii &chrp = Chr2Locus[bchr];
//          start = chrp.first;
            len   = chrp.second;

            cmin = 0;
            cmax = len;

            if (not map_get(locus_filter, bchr, locuslist)) {
                locuslist = new List3ill;
                locus_filter[bchr] = locuslist;
            }
            // use compressed Lesslocus vs start
            locuslist->push_back(make_pair(Lesslocus, Pairll(cmin, cmax)));
            Lesslocus += cmax;
        }
        return chrs;
    }

//  cnt = main_chromocnt;
    cnt = marker_filter.size();

    chrs = CALLOC((size_t) cnt, int);
    chrcnt = 0;
    for (Mfp mfp = marker_filter.begin(); mfp != marker_filter.end(); mfp++) {
        chrs[chrcnt++] = mfp->first;
        // convert filters for chosen chromsomes to offsets
        (*(mfp->second)).sort();
    }
    // sort marker_filters for chosen chromosomes by ascending request positions
    qsort(chrs, chrcnt, sizeof (int), cmp);

    for (i = 0; i < chrcnt; i++) {
        bchr = chrs[i];
        Pairii ii = Chr2Locus[bchr];
        printf("chr %d, st %d, # %d\n", bchr, ii.first, ii.second);
        bplist = marker_filter[chrs[i]];
        if (not map_get(locus_filter, bchr, locuslist)) {
            locuslist = new List3ill;
            locus_filter[bchr] = locuslist;
        }
        iii = 0;
        for (lp = bplist->begin(); lp != bplist->end(); lp++) {
            pp = *lp;
            bmin = pp.first;
            bmax = pp.second;
            printf("bp [%lld, %lld] ", bmin, bmax);

            for (; iii < ii.second; iii++) {
                bp = EXLTop->EXLocus[iii + ii.first].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];            
                if (bmin <= bp) {
                   printf("[%d, %lld] %lld; ", iii+ii.first, bp, bmin);
                    cmin = iii;
                    break;
                }
            }
            for (iii++; iii < ii.second; iii++) {
                bp = EXLTop->EXLocus[iii + ii.first].positions[base_pair_position_index < 0 ? genetic_distance_index : base_pair_position_index];            
                if (bmax == bp) {
                    printf("[%d, %lld] %lld; ", iii+ii.first, bp, bmax);
                    cmax = iii;
                    break;
                } else if (bmax < bp) {
                    printf("[%d, %lld] %lld; ", iii+ii.first, bp, bmax);
                    cmax = iii - 1;
                    break;
                }
            }
            printf("\nk %d, [%ld %ld] del %ld\n\n", Lesslocus, cmin, cmax, cmax-cmin);
            locuslist->push_back(make_pair(Lesslocus, Pairll(cmin, cmax)));
            Lesslocus += cmax - cmin;
        }
    }

    printf("Cnt = %d\n", Lesslocus);

    return chrs;
}

int *loci_reorder_w_locus_filter = NULL;

void compress_loci_w_locus_filter(linkage_ped_top *Top)
{
    linkage_locus_top *LTop = Top->LocusTop;
    ext_linkage_locus_top *EXLTop = Top->EXLTop;
    int offset = LTop->PhenoCnt;

    int bchr, i;
    long long lmin, lmax, oldLoc;
    List3ill  *locuslist;
    Pairill    ill;

    Alleles_int *alleles = NULL;

    int cnt = locus_filter.size();

    int *chrs = CALLOC((size_t) cnt, int);
    int j = 0;
    for (Lfp lfp = locus_filter.begin(); lfp != locus_filter.end(); lfp++) {
        chrs[j++] = lfp->first;
    }
    // sort marker_filters for chosen chromosomes by ascending request positions
    qsort(chrs, j, sizeof (int), cmp);

    linkage_locus_rec *NLocus  = new linkage_locus_rec [ Lesslocus ];
    marker_rec        *NMarker = new marker_rec [ Lesslocus - offset ];
    NMarker -= offset;

    ext_linkage_locus_rec *NEXLocus  = new ext_linkage_locus_rec [ Lesslocus - offset ];
    NEXLocus -= offset;

    if (MARKER_SCHEME == MARKER_SCHEME_BITS)
        alleles = CALLOC(Lesslocus, Alleles_int);

    loci_reorder_w_locus_filter = CALLOC((size_t) LTop->LocusCnt, int);
    memset(loci_reorder_w_locus_filter, -1, (size_t) LTop->LocusCnt * sizeof (int));

    for (i = 0; i < offset; i++) {
        NLocus[i] = LTop->Locus[i];
        loci_reorder_w_locus_filter[i] = i;
    }

    int base;
    Pairii ii;
    List3illp lp;
    long long newLoc = offset;

    for (i = 0; i < cnt; i++) {
        bchr = chrs[i];
        map_get(locus_filter, bchr, locuslist);
        ii = Chr2Locus[bchr];
        base = ii.first;

        for (lp = locuslist->begin(); lp != locuslist->end(); lp++) {
            ill = *lp;
//          k = ill.first;
            lmin = ill.second.first;
            lmax = ill.second.second;
//          printf("gceidx %d, chr %d, base %d, k %lld, [%lld, %lld]\n", i, bchr, base, k, lmin, lmax);
            // newLoc is analog of compressed Locusless vs using start
            for (oldLoc = lmin+base; oldLoc < lmax+base; oldLoc++, newLoc++) {
//              if (oldLoc-lmin-base < 2 || lmax+base-oldLoc < 2)
//                  printf("gceidx %lld, chr %d, newLoc %lld, [%lld, %lld]\n", oldLoc, bchr, newLoc, lmin, lmax);

                NLocus[newLoc]   = LTop->Locus[oldLoc];
                NMarker[newLoc]  = LTop->Marker[oldLoc];
                NLocus[newLoc].Marker = &NMarker[newLoc];
                NEXLocus[newLoc] = EXLTop->EXLocus[oldLoc]; 

                if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
                    extern Alleles_int *MARKER_SCHEME3_alleles;
                    alleles[newLoc] = MARKER_SCHEME3_alleles[oldLoc];
                }

                loci_reorder_w_locus_filter[oldLoc] = newLoc;
            }
        }
    }

    delete [] LTop->Locus;
    LTop->Locus     = NLocus;

    LTop->Marker += offset;
    delete [] LTop->Marker;
    LTop->Marker    = NMarker;

//?? free old position, male, female
    EXLTop->EXLocus += offset;
    free(EXLTop->EXLocus);
    EXLTop->EXLocus = NEXLocus;

    if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        extern Alleles_int *MARKER_SCHEME3_alleles;
        free(MARKER_SCHEME3_alleles);
        MARKER_SCHEME3_alleles = alleles;
    }

    LTop->LocusCnt   = Lesslocus;
    LTop->MarkerCnt  = Lesslocus - offset;
    EXLTop->LocusCnt = Lesslocus;

    for (int i = 0; i < num_reordered; i++) {
        reordered_marker_loci[i] = loci_reorder_w_locus_filter[ reordered_marker_loci[i] ];
    }
}

void dbgenotype_import_genotype(linkage_ped_top *Top) {
    extern void db_open_db();

    char *buf, *bufp;
    int ped, per;

    int *chrs = convert_marker_filter2locus_filter(Top);
    int  cnt  = locus_filter.size(); 

    db_open_db();

    if (chrs != 0 && cnt <= 3) {
//      buf  = CALLOC((size_t) 4 * cnt + 9 /* chr in () */ + 1 /* 999 */+ 1 /* \0 */, char);
        buf  = CALLOC((size_t) 15 * cnt + 9, char);
        bufp = buf;
//      sprintf(bufp, "chr in (%d", chrs[0]);
//      sprintf(bufp, "+chr = %d", chrs[0]);
        sprintf(bufp, "chr = %d", chrs[0]);
        bufp += strlen(bufp);
        for (int i = 1; i < cnt; i++) {
//          sprintf(bufp, ", %d", chrs[i]);
            sprintf(bufp, " OR chr = %d", chrs[i]);
            bufp += strlen(bufp);
        }
//      sprintf(bufp, ")");
        dbgvf("Select * from genotype_table where %s\n", buf);
        genotype_table.init(buf);
        free(buf);
        free(chrs);
    } else {
        dbgvf("Select * from genotype_table\n");
        genotype_table.init();
    }


    Tod imgenraw("import raw genotypes");
    MasterDB.begin();
    genotype_table.db_getall(Top, Genotypes);
    MasterDB.commit();
    imgenraw();

    genotype_table.close();
    MasterDB.close();

    for (ped = 0; ped < Top->PedCnt; ped ++) {
        for (per = 0; per < (Top->PedBroken + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedBroken[ped].Entry[per].person_link;
            Top->PedBroken[ped].Entry[per].Marker = Genotypes[person_link];
        }
        for (per = 0; per < (Top->PedRaw + ped)->EntryCnt; per++) {  // vs Top->PedRaw
            int person_link = Top->PedRaw[ped].Entry[per].person_link;
            Top->PedRaw[ped].Entry[per].Marker = Genotypes[person_link];
        }
    }

    compress_loci_w_locus_filter(Top);
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
}
