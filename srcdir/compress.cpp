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

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "linkage.h"
#include "compress_ext.h"
#include "error_messages_ext.h"
#include "utils_ext.h"
#include "annotated_ped_file.h"

extern int      MARKER_SCHEME;
int             MARKER_SCHEME3_offset = -1;
int             MARKER_SCHEME3_check = 1;

Alleles_int    *MARKER_SCHEME3_alleles = NULL;
Alleles_str    *MARKER_SCHEME3_Ralleles = NULL;
int             MARKER_SCHEME3_mask[]  = {0x03, 0x0c, 0x30, 0xc0};
int             MARKER_SCHEME3_shift[] = {   0,    2,    4,    6};
//int           MARKER_SCHEME3_mask[]  = {0xc0, 0x30, 0x0c, 0x03};
//int           MARKER_SCHEME3_shift[] = {   6,    4,    2,    0};

void           *NOTYPED_ALLELES = (void *) 0;

extern allele_prop **Allele_Array;

int marker_size(int size)
{
    if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        return size * sizeof(marker_pedrec_data);
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        return (size + 3) >>2;
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        return size * sizeof(marker_pedrec_char);
    }
}

void *marker_start(void *marker, int offset) {
    if (marker == NULL) return marker;
    if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        return ((void *) ((marker_pedrec_data *)marker + offset));

    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        return marker;

    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        return ((void *) ((marker_pedrec_char *)marker + offset));
    }
}

void *marker_alloc(size_t size, int offset) {
    if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        return ((void *) ((CALLOC(size, marker_pedrec_data)) - offset));
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        if (MARKER_SCHEME3_offset == -1) {
            MARKER_SCHEME3_offset = offset;
// lazy
            MARKER_SCHEME3_alleles  = CALLOC(size + offset, Alleles_int);
            MARKER_SCHEME3_Ralleles = CALLOC(size + offset, Alleles_str);
        }
        return ((void *) CALLOC((size+3)>>2, unsigned char));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        return ((void *) ((CALLOC(size, marker_pedrec_char)) - offset));
    }
//  return ((void *) 0);
}

void marker_free(void *marker, int offset) {
    if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        free(((marker_pedrec_data *)marker) + offset);
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        free((marker_pedrec_char *)marker);
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        free(((marker_pedrec_char *)marker) + offset);
    }
}

void get_2Ralleles(void *mp, int marker, const char **all1, const char **all2) {
    if (mp == NOTYPED_ALLELES) 
        *all2 = *all1 = REC_UNKNOWN;
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        *all1 = mpd[marker].RAlleles.Allele_1;
        *all2 = mpd[marker].RAlleles.Allele_2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        Alleles_str *allelep = &MARKER_SCHEME3_Ralleles[marker];
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
        if (the_bits == 0) {
            *all2 = *all1 = allelep->Allele_1;

        } else if (the_bits == 1) { // 0
            *all2 = *all1 = REC_UNKNOWN;

        } else if (the_bits == 2) { // ne
            *all1 = allelep->Allele_1;
            *all2 = allelep->Allele_2;

        } else { //  3:
            *all2 = *all1 = allelep->Allele_2;
        }
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        *all1 = Allele_Array[mpd[marker].Allele_1]->name;
        *all2 = Allele_Array[mpd[marker].Allele_2]->name;
    }
}

//for ped stat {p/l}genotype()
int num_typed_2Ralleles(void *mp, int marker) {
    if (mp == NOTYPED_ALLELES) 
        return 0;
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        return (!! allelecmp(mpd[marker].RAlleles.Allele_1, REC_UNKNOWN) ) + (!! allelecmp(mpd[marker].RAlleles.Allele_2, REC_UNKNOWN) );
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
        if (the_bits == 1) { // 0
            return 0;
        } else {
            return 2;
        }
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        return (!! mpd[marker].Allele_1) + (!! mpd[marker].Allele_2);
    }
}

void set_2Ralleles(void *mp, int marker, linkage_locus_rec *locus, const char *all1, const char *all2) {
#ifdef ORDER_HETEROZYGOTE
    if (all1 != all2 && strcmp(all1, all2) > 0) {
        const char *tmp = all1;
        all1 = all2;
        all2 = tmp;
    }
#endif /* ORDER_HETEROZYGOTE */
    if (mp == NOTYPED_ALLELES) 
        ; // do nothing
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        mpd[marker].RAlleles.Allele_1 = all1;
        mpd[marker].RAlleles.Allele_2 = all2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        Alleles_str *allelep = &MARKER_SCHEME3_Ralleles[marker];
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = the_byte & ~MARKER_SCHEME3_mask[get_bits];
        int the_field;
        if (all1 == REC_UNKNOWN && all2 == REC_UNKNOWN) {
            the_field = 1;
        } else {
            if (allelep->Allele_1 == 0) {
                allelep->Allele_1 = all1;
                if (all1 != all2 && allelep->Allele_2 == 0)
                    allelep->Allele_2 = all2;
            } else if (allelep->Allele_2 == 0) {
                if (allelep->Allele_1 != all1)
                    allelep->Allele_2 = all1;
                else if (all1 != all2 ) // Alllele1 = all1
                    allelep->Allele_2 = all2;
            }
            if (MARKER_SCHEME3_check) {
                if (all1 == REC_UNKNOWN || all2 == REC_UNKNOWN) {
                    errorvf("You set the maximum number of alleles to 2.\nHalf type genotypes are not allowed in 2 allele mode: %s/%s.\nPlease adjust the \"maximum number of alleles per marker\" option in the initial input menu.\n",
                            all1, all2);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
                const char * estr = "While you set the maximum number of alleles to 2, there are more than two alleles in the data:\nMarker %s has the alleles %s, %s; trying to add %s.\nPlease adjust the \"maximum number of alleles per marker\" option in the initial input menu.\n";
                if ( (all1 != allelep->Allele_1) && (all1 != allelep->Allele_2) ) {
                    errorvf(estr, locus->LocusName, allelep->Allele_1, allelep->Allele_2, all1);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
                if ( (all2 != allelep->Allele_1) && (all2 != allelep->Allele_2) ) {
                    errorvf(estr, locus->LocusName, allelep->Allele_1, allelep->Allele_2, all2);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
            }

            if (all1 != all2)
                the_field = 2;
            else if (all1 == allelep->Allele_1)
                the_field = 0;
            else
                the_field = 3;
        }
        mpd[get_byte] = ((unsigned char) (the_bits | (the_field << MARKER_SCHEME3_shift[get_bits])));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        mpd[marker].Allele_1 = ((unsigned char) allele2allele_prop_idx(all1));
        mpd[marker].Allele_2 = ((unsigned char) allele2allele_prop_idx(all2));
    }
}

int crunch_Rnotype(void **p, linkage_locus_top *LTop) {
    void *mp = *p;
    int marker;
    int cnt = 0;

    if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        for (marker = LTop->PhenoCnt, mpd = &mpd[marker];
             marker < LTop->LocusCnt; marker++, mpd++) {
            if (mpd->RAlleles.Allele_1 != REC_UNKNOWN || 
                mpd->RAlleles.Allele_2 != REC_UNKNOWN) cnt++;
        }
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        unsigned char *mpd = (unsigned char *) mp;
        for (marker = LTop->PhenoCnt; marker < LTop->LocusCnt; marker++) {
            int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
            int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
            int the_byte = mpd[get_byte];
            int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
            if (the_bits != 1) cnt++;
        }
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        for (marker = LTop->PhenoCnt, mpd = &mpd[marker];
             marker < LTop->LocusCnt; marker++, mpd++) {
            // We know that "0" is loaded as the 0th allele.
            if (mpd->Allele_1 || mpd->Allele_2) cnt++;
        }
    }
    if (cnt == 0) {
        marker_free(mp, LTop->PhenoCnt);
        *p = NOTYPED_ALLELES;
    }
    return cnt;
}

void copy_2Ralleles(void *to, void *from, int marker) {
    if (from == NOTYPED_ALLELES) 
        ; // do nothing
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *frompd = (marker_pedrec_data *) from;
        marker_pedrec_data *tompd  = (marker_pedrec_data *) to;
        tompd[marker].RAlleles.Allele_1 = frompd[marker].RAlleles.Allele_1;
        tompd[marker].RAlleles.Allele_2 = frompd[marker].RAlleles.Allele_2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *frompd = (unsigned char *) from;
        unsigned char *tompd  = (unsigned char *) to;
        int from_byte = frompd[get_byte];
        int from_bits = from_byte & MARKER_SCHEME3_mask[get_bits];
        int to_byte   = tompd[get_byte];
        int to_bits   = to_byte & ~MARKER_SCHEME3_mask[get_bits];
        tompd[get_byte] = ((unsigned char) (to_bits | from_bits));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *frompd = (marker_pedrec_char *) from;
        marker_pedrec_char *tompd  = (marker_pedrec_char *) to;
        tompd[marker].Allele_1 = frompd[marker].Allele_1;
        tompd[marker].Allele_2 = frompd[marker].Allele_2;
    }
}

#if 1
void order_heterozygous_allele_raw(linkage_ped_top *Top)
{
#ifdef ORDER_HETEROZYGOTE
    int i, ped, entrycount, per, cnt = 0, all = 0;

    if (MARKER_SCHEME != MARKER_SCHEME_BITS) return;

    for (ped=0; ped < Top->PedCnt; ped++) {
        entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt :
                                                         Top->PTop[ped].num_persons;
        for(per=0; per < entrycount; per++) {
            unsigned char *mp;
            mp = (pedfile_type == POSTMAKEPED_PFT) ? 
                   (unsigned char *) Top->Ped[ped].Entry[per].Marker :
                   (unsigned char *) Top->PTop[ped].persons[per].marker;

            if (mp != NOTYPED_ALLELES) {
                for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
                    Alleles_str *allelep = &MARKER_SCHEME3_Ralleles[i];
                    if (allelep->Allele_1 && allelep->Allele_2 &&
                        (allelep->Allele_1 != allelep->Allele_2) &&
                        (strcmp(allelep->Allele_1, allelep->Allele_2) > 0)) {
                            int get_byte = (i - MARKER_SCHEME3_offset) >> 2;
                            int get_bits = (i - MARKER_SCHEME3_offset) & 3;
                            int mask     = MARKER_SCHEME3_mask[get_bits];
                            int shift    = MARKER_SCHEME3_shift[get_bits];
                            int the_byte = mp[get_byte];
                            int the_bits = (the_byte & mask) >> shift;
                            if (the_bits == 0 || the_bits == 3) {
                                mp[get_byte] ^= (3 << shift);
                            }
                    }
                }
            }
        }
    }

    for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
        Alleles_str *allelep = &MARKER_SCHEME3_Ralleles[i];
        all++;
//      msgvf("%d: l %d, %s/%s\n", all, i, allelep->Allele_1, allelep->Allele_2);
        if (allelep->Allele_1 && allelep->Allele_2 &&
            (allelep->Allele_1 != allelep->Allele_2) &&
            (strcmp(allelep->Allele_1, allelep->Allele_2) > 0)) {
            const char *tmp = allelep->Allele_1;
            allelep->Allele_1 = allelep->Allele_2;
            allelep->Allele_2 = tmp;
            cnt++;
//          msgvf("%d: l %d, %s/%s\n", cnt, i, allelep->Allele_1, allelep->Allele_2);
        }
    }

    if (cnt > 0)
        msgvf("Fix raw alleles: %d/%d heterozygotes flipped to 1/2.\n", cnt, all);

#endif /* ORDER_HETEROZYGOTE */
}
#else
void order_heterozygous_allele_raw(linkage_ped_top *Top)
{
#ifdef ORDER_HETEROZYGOTE
    int i, ped, entrycount, per, cnt = 0, all = 0;

    if (MARKER_SCHEME != MARKER_SCHEME_BITS) return;

    for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
        Alleles_str *allelep = &MARKER_SCHEME3_Ralleles[i];
        all++;
//      msgvf("%d: l %d, %s/%s\n", all, i, allelep->Allele_1, allelep->Allele_2);
        if (allelep->Allele_1 && allelep->Allele_2 &&
            (allelep->Allele_1 != allelep->Allele_2) &&
            (strcmp(allelep->Allele_1, allelep->Allele_2) > 0)) {
            const char *tmp = allelep->Allele_1;
            allelep->Allele_1 = allelep->Allele_2;
            allelep->Allele_2 = tmp;

            int get_byte = (i - MARKER_SCHEME3_offset) >> 2;
            int get_bits = (i - MARKER_SCHEME3_offset) & 3;
            int mask     = MARKER_SCHEME3_mask[get_bits];
            int shift    = MARKER_SCHEME3_shift[get_bits];

            cnt++;
//          msgvf("%d: l %d, %s/%s\n", cnt, i, allelep->Allele_1, allelep->Allele_2);
            for (ped=0; ped < Top->PedCnt; ped++) {
                entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt :
                                                                 Top->PTop[ped].num_persons;
                for(per=0; per < entrycount; per++) {
                    unsigned char *mp;
                    mp = (pedfile_type == POSTMAKEPED_PFT) ? 
                           (unsigned char *) Top->Ped[ped].Entry[per].Marker :
                           (unsigned char *) Top->PTop[ped].persons[per].marker;
                    if (mp == NOTYPED_ALLELES) continue;

                    int the_byte = mp[get_byte];
                    int the_bits = (the_byte & mask) >> shift;
                    if (the_bits == 0 || the_bits == 3) {
                        mp[get_byte] ^= (3 << shift);
                    }
                }
            }
        }
    }
    if (cnt > 0)
        msgvf("Fix raw alleles: %d/%d heterozygotes flipped to 1/2.\n", cnt, all);
#endif /* ORDER_HETEROZYGOTE */
}
#endif

inline void decode_compression(int the_bits, Alleles_int *allelep, int *all1, int *all2)
{
    if (the_bits == 0) {
        *all2 = *all1 = allelep->Allele_1;

    } else if (the_bits == 1) { // 0
        *all2 = *all1 = 0;

    } else if (the_bits == 2) { // ne
        *all1 = allelep->Allele_1;
        *all2 = allelep->Allele_2;

    } else { // 3:
        *all2 = *all1 = allelep->Allele_2;
    }
}

void get_2alleles(void *mp, int marker, int *all1, int *all2) {
    if (mp == NOTYPED_ALLELES) 
        *all2 = *all1 = 0;
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        *all1 = mpd[marker].Alleles.Allele_1;
        *all2 = mpd[marker].Alleles.Allele_2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
//      Alleles_int *allelep = &MARKER_SCHEME3_alleles[marker];
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
        decode_compression(the_bits, &MARKER_SCHEME3_alleles[marker], all1, all2);

    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        *all1 = mpd[marker].Allele_1;
        *all2 = mpd[marker].Allele_2;
    }
}

//for ped stat {p/l}genotype()
int num_typed_2alleles(void *mp, int marker) {
    if (mp == NOTYPED_ALLELES) 
        return 0;
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
        return (!! mpd[marker].Alleles.Allele_1) + (!! mpd[marker].Alleles.Allele_2);
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
        if (the_bits == 1) { // 0
            return 0;
        } else {
            return 2;
        }
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        return (!! mpd[marker].Allele_1) + (!! mpd[marker].Allele_2);
    }
}

void set_2alleles(void *mp, int marker, linkage_locus_rec *locus, int all1, int all2) {
#ifdef ORDER_HETEROZYGOTE
    if (all1 > all2) {
        int tmp = all1;
        all1 = all2;
        all2 = tmp;
    }
#endif /* ORDER_HETEROZYGOTE */
    if (mp == NOTYPED_ALLELES) 
        ; // do nothing
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
         mpd[marker].Alleles.Allele_1 = all1;
         mpd[marker].Alleles.Allele_2 = all2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        Alleles_int *allelep = &MARKER_SCHEME3_alleles[marker];
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *mpd = (unsigned char *) mp;
        int the_byte = mpd[get_byte];
        int the_bits = the_byte & ~MARKER_SCHEME3_mask[get_bits];
        int the_field;
        if (all1 == 0 && all2 == 0) {
            the_field = 1;
        } else {
            if (allelep->Allele_1 == 0) {
                allelep->Allele_1 = all1;
                if (all1 != all2 && allelep->Allele_2 == 0)
                    allelep->Allele_2 = all2;
            } else if (allelep->Allele_2 == 0) {
                if (allelep->Allele_1 != all1)
                    allelep->Allele_2 = all1;
                else if (all1 != all2 ) // Alllele1 = all1
                    allelep->Allele_2 = all2;
            } 
            if (MARKER_SCHEME3_check) {
                if (all1 == 0 || all2 == 0) {
                    errorvf("You set the maximum number of alleles to 2.\nHalf type genotypes are not allowed in 2 allele mode: %d/%d.\nPlease adjust the \"maximum number of alleles per marker\" option in the initial input menu.\n",
                            all1, all2);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
                const char * estr = "While you set the maximum number of alleles to 2, there are more than two alleles in the data:\nMarker %s has the alleles %d, %d; trying to add%d.\nPlease adjust the \"maximum number of alleles per marker\" option in the initial input menu.\n";
                if ( (all1 != allelep->Allele_1) && (all1 != allelep->Allele_2) ) {
                    errorvf(estr, locus->LocusName, allelep->Allele_1, allelep->Allele_2, all1);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
                if ( (all2 != allelep->Allele_1) && (all2 != allelep->Allele_2) ) {
                    errorvf(estr, locus->LocusName, allelep->Allele_1, allelep->Allele_2, all2);
                    EXIT(OUTOF_BOUNDS_ERROR);
                }

            }

            if (all1 != all2)
                the_field = 2;
            else if (all1 == allelep->Allele_1)
                the_field = 0;
            else
                the_field = 3;

        }
        mpd[get_byte] = ((unsigned char) (the_bits | (the_field << MARKER_SCHEME3_shift[get_bits])));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        /*
         * Checks to make sure that all1/2 fit int uchar have been made earlier:
         *  in canonical_allele() and in read_numbered_data();
         */
        marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
        mpd[marker].Allele_1 = ((unsigned char) all1);
        mpd[marker].Allele_2 = ((unsigned char) all2);
    }
}

int crunch_notype(void **p, linkage_locus_top *LTop) {
    void *mp = *p;
    int marker;
    int cnt = 0;

    for (marker = LTop->PhenoCnt; marker < LTop->LocusCnt; marker++) {
        if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
            marker_pedrec_data *mpd = (marker_pedrec_data *) mp;
            if (mpd[marker].Alleles.Allele_1 != 0 ||
                mpd[marker].Alleles.Allele_2 != 0) cnt++; 
        } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
            int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
            int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
            unsigned char *mpd = (unsigned char *) mp;
            int the_byte = mpd[get_byte];
            int the_bits = (the_byte & MARKER_SCHEME3_mask[get_bits]) >> MARKER_SCHEME3_shift[get_bits];
            if (the_bits != 1) cnt++;
        } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
            marker_pedrec_char *mpd = (marker_pedrec_char *) mp;
            if (mpd[marker].Allele_1 || mpd[marker].Allele_2) cnt++;
        }
    }
    if (cnt == 0) {
        marker_free(mp, LTop->PhenoCnt);
        *p = NOTYPED_ALLELES;
    }
    return cnt;
}

void copy_2alleles(void *to, void *from, int marker) {
    if (from == NOTYPED_ALLELES) 
        ; // do nothing
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *frompd = (marker_pedrec_data *) from;
        marker_pedrec_data *tompd  = (marker_pedrec_data *) to;
        tompd[marker].Alleles.Allele_1 = frompd[marker].Alleles.Allele_1;
        tompd[marker].Alleles.Allele_2 = frompd[marker].Alleles.Allele_2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        int get_byte = (marker - MARKER_SCHEME3_offset) >> 2;
        int get_bits = (marker - MARKER_SCHEME3_offset) & 3;
        unsigned char *frompd = (unsigned char *) from;
        unsigned char *tompd  = (unsigned char *) to;
        int from_byte = frompd[get_byte];
        int from_bits = from_byte & MARKER_SCHEME3_mask[get_bits];
        int to_byte   = tompd[get_byte];
        int to_bits   = to_byte & ~MARKER_SCHEME3_mask[get_bits];
        tompd[get_byte] = ((unsigned char) (to_bits | from_bits));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *frompd = (marker_pedrec_char *) from;
        marker_pedrec_char *tompd  = (marker_pedrec_char *) to;
        tompd[marker].Allele_1 = frompd[marker].Allele_1;
        tompd[marker].Allele_2 = frompd[marker].Allele_2;
    }
}

void copy_2alleles(void *to, void *from, int tomarker, int frommarker) {
    if (from == NOTYPED_ALLELES) 
        ; // do nothing
    else if (MARKER_SCHEME == MARKER_SCHEME_PTR) {
        marker_pedrec_data *frompd = (marker_pedrec_data *) from;
        marker_pedrec_data *tompd  = (marker_pedrec_data *) to;
        tompd[tomarker].Alleles.Allele_1 = frompd[frommarker].Alleles.Allele_1;
        tompd[tomarker].Alleles.Allele_2 = frompd[frommarker].Alleles.Allele_2;
    } else if (MARKER_SCHEME == MARKER_SCHEME_BITS) {
        unsigned char *frompd = (unsigned char *) from;
        int get_from_byte = (frommarker - MARKER_SCHEME3_offset) >> 2;
        int get_from_bits = (frommarker - MARKER_SCHEME3_offset) & 3;
        unsigned char *tompd  = (unsigned char *) to;
        int get_to_byte = (tomarker - MARKER_SCHEME3_offset) >> 2;
        int get_to_bits = (tomarker - MARKER_SCHEME3_offset) & 3;

        int from_byte = frompd[get_from_byte];
        int from_bits = from_byte & MARKER_SCHEME3_mask[get_from_bits];
        int to_byte   = tompd[get_to_byte];
        int to_bits   = to_byte & ~MARKER_SCHEME3_mask[get_to_bits];
        tompd[get_to_byte] = ((unsigned char) (to_bits | from_bits));
    } else { // if (MARKER_SCHEME == MARKER_SCHEME_BYTE)
        marker_pedrec_char *frompd = (marker_pedrec_char *) from;
        marker_pedrec_char *tompd  = (marker_pedrec_char *) to;
        tompd[tomarker].Allele_1 = frompd[frommarker].Allele_1;
        tompd[tomarker].Allele_2 = frompd[frommarker].Allele_2;
    }
}

#if 1
void order_heterozygous_allele(linkage_ped_top *Top)
{
#ifdef ORDER_HETEROZYGOTE
    int i, ped, entrycount, per, cnt = 0, all = 0;

    if (MARKER_SCHEME != MARKER_SCHEME_BITS) return;

    for (ped=0; ped < Top->PedCnt; ped++) {
        entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt :
                                                         Top->PTop[ped].num_persons;
        for(per=0; per < entrycount; per++) {
            unsigned char *mp;
            mp = (pedfile_type == POSTMAKEPED_PFT) ? 
                   (unsigned char *) Top->Ped[ped].Entry[per].Marker :
                   (unsigned char *) Top->PTop[ped].persons[per].marker;

            if (mp != NOTYPED_ALLELES) {
                for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
                    Alleles_int *allelep = &MARKER_SCHEME3_alleles[i];
                    if (allelep->Allele_1 && allelep->Allele_2 &&
                        (allelep->Allele_1 > allelep->Allele_2)) {
                            int get_byte = (i - MARKER_SCHEME3_offset) >> 2;
                            int get_bits = (i - MARKER_SCHEME3_offset) & 3;
                            int mask     = MARKER_SCHEME3_mask[get_bits];
                            int shift    = MARKER_SCHEME3_shift[get_bits];
                            int the_byte = mp[get_byte];
                            int the_bits = (the_byte & mask) >> shift;
                            if (the_bits == 0 || the_bits == 3) {
                                mp[get_byte] ^= (3 << shift);
                            }
                    }
                }
            }
        }
    }

    for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
        Alleles_int *allelep = &MARKER_SCHEME3_alleles[i];
        all++;
//      msgvf("%d: l %d, %s/%s\n", all, i, allelep->Allele_1, allelep->Allele_2);
        if (allelep->Allele_1 && allelep->Allele_2 &&
            (allelep->Allele_1 > allelep->Allele_2)) {
            int tmp = allelep->Allele_1;
            allelep->Allele_1 = allelep->Allele_2;
            allelep->Allele_2 = tmp;
            cnt++;
//          msgvf("%d: l %d, %d/%d\n", cnt, i, allelep->Allele_1, allelep->Allele_2);
        }
    }

    if (cnt > 0)
        msgvf("Fix alleles: %d/%d heterozygotes flipped to 1/2.\n", cnt, all);
#endif /* ORDER_HETEROZYGOTE */
}
#else
void order_heterozygous_allele(linkage_ped_top *Top)
{
#ifdef ORDER_HETEROZYGOTE
    int i, ped, entrycount, per, cnt = 0, all = 0;

    if (MARKER_SCHEME != MARKER_SCHEME_BITS) return;

    for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++) {
        Alleles_int *allelep = &MARKER_SCHEME3_alleles[i];
        all++;
        if (allelep->Allele_1 && allelep->Allele_2 &&
            (allelep->Allele_1 > allelep->Allele_2)) {
            int tmp = allelep->Allele_1;
            allelep->Allele_1 = allelep->Allele_2;
            allelep->Allele_2 = tmp;

            int get_byte = (i - MARKER_SCHEME3_offset) >> 2;
            int get_bits = (i - MARKER_SCHEME3_offset) & 3;
            int mask     = MARKER_SCHEME3_mask[get_bits];
            int shift    = MARKER_SCHEME3_shift[get_bits];

            cnt++;
            for (ped=0; ped < Top->PedCnt; ped++) {
                entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt :
                                                                 Top->PTop[ped].num_persons;
                for(per=0; per < entrycount; per++) {
                    unsigned char *mp;
                    mp = (pedfile_type == POSTMAKEPED_PFT) ? 
                           (unsigned char *) Top->Ped[ped].Entry[per].Marker :
                           (unsigned char *) Top->PTop[ped].persons[per].marker;
                    if (mp == NOTYPED_ALLELES) continue;

                    int the_byte = mp[get_byte];
                    int the_bits = (the_byte & mask) >> shift;
                    if (the_bits == 0 || the_bits == 3) {
                        mp[get_byte] ^= (3 << shift);
                    }
                }
            }
        }
    }
    if (cnt > 0)
        msgvf("Fix alleles: %d/%d heterozygotes flipped to 1/2.\n", cnt, all);
#endif /* ORDER_HETEROZYGOTE */
}
#endif

