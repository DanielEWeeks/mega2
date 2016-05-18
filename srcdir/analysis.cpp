/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
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

*/

/*
 * You MUST make the appropriate changes as indicated below to add a new analysis:
 * See changes 1. - 4.
 */
#include <stdio.h>
#include "common.h"
#include "typedefs.h"

#include "ctype.h"
#include "error_messages_ext.h"
#include "utils_ext.h"

#include "analysis.h"
#include "class_old.h"

#include "write_fbat_ext.h"
#include "write_plink_ext.h"
#include "write_eigenstrat_ext.h"
#include "write_pangaea_ext.h"
#include "write_beagle_ext.h"
#include "write_structure_ext.h"
#include "write_pseq_ext.h"
#include "write_shapeit_ext.h"
#include "write_impute2_ext.h"
#include "write_roadtrips_ext.h"
#include "write_mach_ext.h"
/*
 * You MUST make the appropriate changes here to define the new header and class.
 * (item 1.)
 */


CLASS_HAPLOTYPE          *HAPLOTYPE = new CLASS_HAPLOTYPE();
CLASS_LOCATION           *LOCATION = new CLASS_LOCATION();
CLASS_NONPARAMETRIC      *NONPARAMETRIC = new CLASS_NONPARAMETRIC();
CLASS_IBD_EST            *IBD_EST = new CLASS_IBD_EST();
CLASS_MISTYPING          *MISTYPING = new CLASS_MISTYPING();
CLASS_MENDEL             *TO_MENDEL = new CLASS_MENDEL();
CLASS_ASPEX              *TO_ASPEX = new CLASS_ASPEX();
CLASS_CREATE_SUMMARY     *CREATE_SUMMARY = new CLASS_CREATE_SUMMARY();
CLASS_GENOTYPING_SUMMARY *GENOTYPING_SUMMARY = new CLASS_GENOTYPING_SUMMARY();
CLASS_GENEHUNTER         *TO_GeneHunter = new CLASS_GENEHUNTER();
CLASS_APM                *TO_APM = new CLASS_APM();
CLASS_APM_MULT           *TO_APM_MULT = new CLASS_APM_MULT();
CLASS_NUKE               *TO_NUKE = new CLASS_NUKE();
CLASS_LIABLE_FREQ        *TO_LIABLE_FREQ = new CLASS_LIABLE_FREQ();
CLASS_ALLELE_FREQ        *TO_ALLELE_FREQ = new CLASS_ALLELE_FREQ();
CLASS_SLINK              *TO_SLINK = new CLASS_SLINK();
CLASS_SPLINK             *TO_SPLINK = new CLASS_SPLINK();
CLASS_LOD2               *TO_LOD2 = new CLASS_LOD2();
CLASS_GENEHUNTERPLUS     *TO_GeneHunterPlus = new CLASS_GENEHUNTERPLUS();
CLASS_SIMULATE           *TO_SIMULATE = new CLASS_SIMULATE();
CLASS_SAGE               *TO_SAGE = new CLASS_SAGE();
CLASS_TDTMAX             *TO_TDTMAX = new CLASS_TDTMAX();
CLASS_SOLAR              *TO_SOLAR = new CLASS_SOLAR();
CLASS_HWETEST            *TO_HWETEST = new CLASS_HWETEST();
CLASS_VITESSE            *TO_VITESSE = new CLASS_VITESSE();
CLASS_LINKAGE            *TO_LINKAGE = new CLASS_LINKAGE();
CLASS_ALLEGRO            *TO_Allegro = new CLASS_ALLEGRO();
CLASS_GHMLB              *TO_GHMLB = new CLASS_GHMLB();
CLASS_SAGE4              *TO_SAGE4 = new CLASS_SAGE4();
CLASS_PREMAKEPED         *TO_PREMAKEPED = new CLASS_PREMAKEPED();
CLASS_MERLIN             *TO_MERLIN = new CLASS_MERLIN();
CLASS_PREST              *TO_PREST = new CLASS_PREST();
CLASS_PAP                *TO_PAP = new CLASS_PAP();
CLASS_MERLINONLY         *TO_MERLINONLY = new CLASS_MERLINONLY();
CLASS_QUANT_SUMMARY      *QUANT_SUMMARY = new CLASS_QUANT_SUMMARY();
CLASS_LOKI               *TO_LOKI = new CLASS_LOKI();
CLASS_MENDEL4            *TO_MENDEL4 = new CLASS_MENDEL4();
CLASS_SUP                *TO_SUP = new CLASS_SUP();
CLASS_PLINK              *TO_PLINK = new CLASS_PLINK();
CLASS_MENDEL7_CSV        *TO_MENDEL7_CSV = new CLASS_MENDEL7_CSV();
CLASS_CRANEFOOT          *CRANEFOOT = new CLASS_CRANEFOOT();
CLASS_MEGA2ANNOT         *MEGA2ANNOT = new CLASS_MEGA2ANNOT();
CLASS_IQLS               *IQLS = new CLASS_IQLS();

CLASS_FBAT               *FBAT = new CLASS_FBAT();

CLASS_SIMWALK2           *TO_SIMWALK2 = new CLASS_SIMWALK2();
CLASS_SUMMARY            *TO_SUMMARY  = new CLASS_SUMMARY();
CLASS_PANGAEA            *PANGAEA = new CLASS_PANGAEA();

CLASS_BEAGLE             *BEAGLE = new CLASS_BEAGLE();
/**
   These classes define the sub-options associated with the different Genotype file
   that Beagle can produce. See the section on 'Beagle Genotype files'.
 */
CLASS_BEAGLE_UNPHASED_UNRELATED *BEAGLE_UNPHASED_UNRELATED = new CLASS_BEAGLE_UNPHASED_UNRELATED();
CLASS_BEAGLE_UNPHASED_TRIO *BEAGLE_UNPHASED_TRIO = new CLASS_BEAGLE_UNPHASED_TRIO();
CLASS_BEAGLE_UNPHASED_PAIR *BEAGLE_UNPHASED_PAIR = new CLASS_BEAGLE_UNPHASED_PAIR();

CLASS_EIGENSTRAT         *EIGENSTRAT = new CLASS_EIGENSTRAT();
CLASS_STRUCTURE          *STRUCTURE = new CLASS_STRUCTURE();
CLASS_PSEQ               *TO_PSEQ = new CLASS_PSEQ();

CLASS_SHAPEIT            *SHAPEIT = new CLASS_SHAPEIT();
CLASS_IMPUTE2            *IMPUTE2 = new CLASS_IMPUTE2();
CLASS_ROADTRIPS          *ROADTRIPS = new CLASS_ROADTRIPS();
CLASS_MACH               *MACH = new CLASS_MACH();
/*
 * You MUST make the appropriate changes here to define the new analysis object
 * (item 2.)
 */


analysis_types analysis_list[] = {
    { "SimWalk2 format",            TO_SIMWALK2 },
    { "Vintage Mendel format",      TO_MENDEL },
    { "ASPEX format",               TO_ASPEX },
    { "GeneHunter-Plus format",     TO_GeneHunterPlus },
    { "GeneHunter format",          TO_GeneHunter },
    { "APM format [DISABLED]",      TO_APM },
    { "APM-MULT format [DISABLED]", TO_APM_MULT },
    { "Create nuclear families",    TO_NUKE },
    { "SLINK format",               TO_SLINK },
    { "SPLINK format",              TO_SPLINK },
    { "Homogeneity analyses",       TO_LOD2 },
    { "SIMULATE format",            TO_SIMULATE },
    { "Create summary files",       TO_SUMMARY },
    { "Old SAGE format",            TO_SAGE },
    { "TDTMax analyses [DISABLED]", TO_TDTMAX },
    { "SOLAR format",               TO_SOLAR },
    { "Vitesse format",             TO_VITESSE },
    { "Linkage format",             TO_LINKAGE },
    { "Test loci for HWE",          TO_HWETEST },
    { "Allegro format",             TO_Allegro },
    { "MLBQTL format",              TO_GHMLB },
    { "SAGE format",                TO_SAGE4 },
    { "Pre-makeped format",         TO_PREMAKEPED },
    { "Merlin/SimWalk2-NPL format", TO_MERLIN },
    { "PREST format",               TO_PREST },
    { "PAP format",                 TO_PAP },
    { "Merlin format",              TO_MERLINONLY },
    { "Loki format",                TO_LOKI },
    { "Mendel format",              TO_MENDEL7_CSV },
    { "SUP format",                 TO_SUP },
    { "PLINK format",               TO_PLINK },
    { "CRANEFOOT format",           CRANEFOOT },
    { "Mega2 format",               MEGA2ANNOT },
    { "IQLS/Idcoefs format",        IQLS },
    { "FBAT format",                FBAT },
    { "PANGAEA MORGAN format",      PANGAEA },
    { "Beagle format",              BEAGLE },
    { "Eigenstrat format",          EIGENSTRAT },
    { "Structure format",           STRUCTURE },
    { "PSEQ format",                TO_PSEQ },
    { "SHAPEIT format",             SHAPEIT},
//  { "IMPUTE2 format",             IMPUTE2}
    { "RoadTrips format",           ROADTRIPS},
    { "MaCH format",                MACH},
/*
 * You MUST make the appropriate changes here to define the mapping from the name
 * to the new analysis object
 * (item 3.)
 */
};

int count_analysis_list = sizeof(analysis_list) / sizeof (analysis_types);

#include "menu_value_missing.h"
/*
 * You MUST make the appropriate changes in the menu_value_missing.h file to show how the
 * new analysis object represents missing quants and affects
 * (item 4.)
 */


//
// The two routines 'prog_name()', and 'prog_name_to_num()' are intertwined.
// 'prog_name()' is used to output a string to the batch file that is understandable to the user.
// 'prog_name_to_num()' does the inverse converting the string to a number that is used internally.
// So clearly, if you change one, you must change the other.
//
// The two routines 'prog_name()', and 'prog_name_to_num()' are intertwined.
// 'prog_name()' is used to output a string to the batch file that is understandable to the user.
// 'prog_name_to_num()' does the inverse converting the string to a number that is used internally.
// So clearly, if you change one, you must change the other.

//typedef UM<const char *, analysis_type *, HH<const char *> charseq> Prog_to_analysis;
typedef std::map<char *, analysis_type, charsless> Prog_to_analysis;
Prog_to_analysis P2A;

void init_analysis()
{
    analysis_types *al = analysis_list;
    char *key;
    for (int i = 0; i < count_analysis_list; i++, al++) {
        key = strdup(al->option_name);
        for (char *cp = key; *cp; cp++) *cp = tolower(*cp);
        P2A[key] = al->analysis;  // key/string from analysis_list item

        key = strdup(al->analysis->_name);
        for (char *cp = key; *cp; cp++) *cp = tolower(*cp);
        P2A[key] = al->analysis;  // key/string from CLASS_ANALYSIS subclass
    }
}

void prog_name_to_num(char *prog_name, analysis_type *analysis)
{
    *analysis = (analysis_type )0;

    char *key = strdup(prog_name);
    for (char *cp = key; *cp; cp++) *cp = tolower(*cp);
    analysis_type a = NULL;
    if (map_get(P2A, key, a)) {
        *analysis = a;
        return;
    }
    errorvf("Analysis: prog_name_to_num() key (%s) lookup failed\n", key);

#if 0
    switch(tolower((unsigned char)prog_name[0])) {
    case 's':
        /* Simwalk2, Summary, Sage 3 and 4, Simulate, splink, slink, solar,
           sup, structure */

        switch(tolower((unsigned char)prog_name[1])) {  // sX...

        case 'h':
            *analysis = SHAPEIT;
            break;

        case 'i':
            /* simwalk2 options  or simulate */
            if (tolower((unsigned char)prog_name[3]) == 'w') {
                *analysis = TO_SIMWALK2; // 1. SimWalk2 
            } else if (tolower((unsigned char)prog_name[3]) == 'u') {
                *analysis = TO_SIMULATE; // 12. Simulate
            } else {
                unknown_prog(prog_name);
            }
            break;

        case 'u':
            /* summary or sup */
            if (tolower((unsigned char)prog_name[2]) == 'm') {
                *analysis = TO_SUMMARY; // 13. Summary
            } else if (tolower((unsigned char)prog_name[2]) == 'p') {
                *analysis = TO_SUP; // 30. SUP
            } else {
                unknown_prog(prog_name);
            }
            break;

        case '.':
            /* sage 3 or sage 4 */
            if (prog_name[strlen(prog_name) - 3] == '3') {
                *analysis = TO_SAGE; // 14. S.A.G.E.3.0
            } else if (prog_name[strlen(prog_name) - 3] == '4') {
                *analysis = TO_SAGE4; // 22. S.A.G.E.4.0
            } else {
                unknown_prog(prog_name);
            }
            break;

        case 'p':
            *analysis = TO_SPLINK; // 10. Splink
            break;

        case 'l':
            *analysis = TO_SLINK; // 9. Slink
            break;

        case 'o':
            *analysis = TO_SOLAR; // 16. SOLAR
            break;

        case 't':
            *analysis = STRUCTURE; // 39. STRUCTURE
            break;

        default :
            unknown_prog(prog_name);
            break;
        }
        break;

    case 'm':
        /* mendel 3 and 9, mlbqtl, merlin, merlin-simwalk2, mega2 */
        switch(tolower((unsigned char)prog_name[strlen(prog_name) - 1])) {  // m...X

        case 'l':
            if (tolower((unsigned char)prog_name[3]) == 'q') {
	      *analysis = TO_GHMLB; // 21. MLBQTL
            } else if (tolower((unsigned char)prog_name[3]) == 'd') {
                /* mendel 3*/
	      *analysis = TO_MENDEL; // 2. Mendel
            } else {
                unknown_prog(prog_name);
            }
            break;

        case '+':
            *analysis = TO_MENDEL7_CSV; // 29. Mendel7+ 
            break;

        case 'n':
            /* merlin only */
            *analysis = TO_MERLINONLY; // 27. Merlin 
            break;

        case '2':
            if (prog_name[strlen(prog_name) - 2] == 'a') {
	      *analysis = MEGA2ANNOT; // 33. Mega2
            } else {
	      *analysis = TO_MERLIN; // 24. Merlin/SimWalk2 
            }
            break;

        case 'h':
            *analysis = MACH; //MACH
            break;


        default:
            unknown_prog(prog_name);
            break;
        }
        break;

    case 'a':  // aX...
        /* aspex, allegro */
        if (tolower((unsigned char)prog_name[1] == 's')) {
	  *analysis = TO_ASPEX; // 3. Aspex 
        } else if (prog_name[2] == 'l') {
	  *analysis = TO_Allegro; // 20. Allegro 
        } else {
            unknown_prog(prog_name);
        }
        break;

    case 'b':  // bX...
        *analysis = BEAGLE; // 37. Beagle
        break;

    case 'g': // g...X
        /* genehunter, genehunter-plus */
        if (prog_name[strlen(prog_name) - 1] == 's') {
            *analysis = TO_GeneHunterPlus; // 4. GeneHunter-Plus 
        } else if (tolower((unsigned char)prog_name[strlen(prog_name) - 1]) == 'r') {
            *analysis = TO_GeneHunter; // 5. GeneHunter
        } else {
            unknown_prog(prog_name);
        }
        break;

    case 'n':
        *analysis = TO_NUKE; // 8. Nuclear families 
        break;

    case 'h':  // hX
        /* homogeneity, hardy-weinberg */
        if (tolower((unsigned char)prog_name[1]) == 'a') {
            *analysis = TO_HWETEST; // 19. Hardy-Weinberg 
        } else if (tolower((unsigned char)prog_name[1]) == 'o') {
            *analysis = TO_LOD2; // 11. Homogeneity 
        } else {
            unknown_prog(prog_name);
        }
        break;

    case 'p': // p...X
        /* premakeped, prest, pap, plink, pangaea */
        switch(tolower((unsigned char)prog_name[strlen(prog_name) - 1])) {

        case 'd':
            *analysis = TO_PREMAKEPED; // 23. Premakeped 
            break;

        case 't':
            *analysis = TO_PREST; // 25. Prest 
            break;

        case 'p':
            *analysis = TO_PAP; // 26. PAP
            break;

        case 'k':
            *analysis = TO_PLINK; // 31. PLINK 
            break;

        case 'q':
            *analysis = TO_PSEQ; // 40. PSEQ
            break;

        case 'a':
            *analysis = PANGAEA; // 36. PANGAEA MORGAN 
            break;

        default:
            unknown_prog(prog_name);
            break;
        }
        break;

    case 'l': // l...X
        /* loki, linkage */
        if (tolower((unsigned char)prog_name[strlen(prog_name) - 1]) == 'i') {
          *analysis = TO_LOKI; // 28. Loki 
        } else if (tolower((unsigned char)prog_name[strlen(prog_name) - 1]) == 'e') {
          *analysis = TO_LINKAGE; // 18. Linkage 
        } else {
            unknown_prog(prog_name);
        }
        break;

    case 'v':
        *analysis = TO_VITESSE; // 17. Vitesse 
        break;

    case 'c' :
        *analysis = CRANEFOOT; // 32. CRANEFOOT 
        break;

    case 'i':
        if (tolower((unsigned char)prog_name[1]) == 'q') {
            *analysis = IQLS; // 34. IQLS/Idcoefs
        } else if (tolower((unsigned char)prog_name[1]) == 'm') {
            *analysis = IMPUTE2;
        } else {
            unknown_prog(prog_name);
        }
        break;

    case 'f':
        *analysis = FBAT; // 35. FBAT
        break;

    case 'e':
        *analysis = EIGENSTRAT; // 38. EIGENSTRAT
        break;

    case 'r':
        *analysis = ROADTRIPS;
        break;

    default:
        unknown_prog(prog_name);
    }

    msgvf("Analysis: %p (%s) == %p (%s)\n",
          *analysis, (*analysis)->_name, tmp, tmp->_name);
    if (*analysis != tmp) {
        msgvf("Analysis: %p (%s) != %p (%s)\n",
              *analysis, (*analysis)->_name, tmp, tmp->_name);
    }
#endif
}

