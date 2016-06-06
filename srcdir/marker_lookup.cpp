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

#include "common.h"
#include "typedefs.h"

void clear_marker();
int test_and_add_marker(const char *key, int value);

void add_marker(const char *key, int value);
void make_marker(linkage_locus_top *LTop);
int  search_marker(const char *marker, int *index);
void debug_dump_marker_dict();

void add_allele(const char *key, char *value);
char *search_allele(const char *allele);
void debug_dump_allele_dict();
void allele_prop_reset();

//#define HASH
#ifdef  HASH

#ifdef _WIN
#include <hash_map>
#else
#include <ext/hash_map>
#endif /* _WIN */

using namespace std;

#ifndef _WIN
using namespace __gnu_cxx;
#endif /* _WIN */

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

typedef hash_map<const char *,int,hash<char *>,eqstr> dict_s2i;
typedef hash_map<const char *,int,hash<char *>,eqstr>::iterator dict_s2ip;
dict_s2i  MARKERS;

typedef hash_map<const char *,char *,hash<char *>,eqstr> dict_s2s;
typedef hash_map<const char *,char *,hash<char *>,eqstr>::iterator dict_s2sp;
dict_s2s  ALLELES;
dict_s2sp ALLELES_fini = ALLELES.end();

#else

#include <map>
using namespace std;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef map<const char *,int,ltstr>  dict_s2i;
typedef map<const char *,int,ltstr>::iterator dict_s2ip;
dict_s2i MARKERS;

typedef map<const char *,char *,ltstr>  dict_s2s;
typedef map<const char *,char *,ltstr>::iterator dict_s2sp;
dict_s2s  ALLELES;
dict_s2sp ALLELES_fini = ALLELES.end();

#endif

void clear_marker()
{
  MARKERS.clear();
}

// Return -1 if the key was not found, otherwise the value of the key...
int test_and_add_marker(const char *key, int value)
{
  size_t cnt = MARKERS.count(key);
  if (cnt == 0) {
     MARKERS[key] = value;
     return -1;
  }
  return MARKERS[key];
}

void add_marker(const char *key, int value)
{
    MARKERS[key] = value;
}

void make_marker(linkage_locus_top *LTop)
{
    int m;

    MARKERS.clear();

    for (m=0; m < LTop->LocusCnt; m++) {
        MARKERS[LTop->Locus[m].LocusName] = m;
    }

    return;
}

void debug_dump_marker_dict()
{
    printf("Dumping MARKERS:\n");
    dict_s2ip n = MARKERS.begin(), e = MARKERS.end();
    for (; n != e; n++) {
        printf("key %s, val %d\n", n->first, n->second);
        fflush(stdout);
    }
}

int search_marker(const char *marker, int *index)
{
    int success = 0;

    dict_s2ip ans_new = MARKERS.find(marker);

    if (ans_new != MARKERS.end()) {
        *index  = ans_new->second;
        success = 1;
    } else {
//        sprintf(err_msg, "marker lookup failed for  \"%s\"\n", marker);
//        warnf(err_msg);
        success = 0;
    }

    return success;
}


void add_allele(const char *key, char *value)
{
    ALLELES[key] = value;
    ALLELES_fini = ALLELES.end();
}

void debug_dump_allele_dict()
{
    dict_s2sp n = ALLELES.begin(), e = ALLELES.end();
    int i = 0;
    for (; n != e; n++) {
        printf("%d: key %s(%p), val %s(%p)\n", i++, n->first, n->first, n->second, n->second);
        fflush(stdout);
    }
}

char *search_allele(const char *allele)
{
    dict_s2sp ans_new = ALLELES.find(allele);

    if (ans_new != ALLELES_fini) {
        return ans_new->second;
    } else {
        return NULL;
    }
}

void allele_prop_reset()
{
    dict_s2sp n = ALLELES.begin(), e = ALLELES.end();
    for (; n != e; n++) {
        allele2allele_prop_prop(n->second) = (void *) 0;
        fflush(stdout);
    }
}
