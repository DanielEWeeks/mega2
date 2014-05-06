/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2014 Robert Baron, Charles P. Kollar,
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

#include <stdio.h>
#include "common.h"
#include "typedefs.h"

#include "ctype.h"
#include "utils_ext.h"


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
CLASS_EIGENSTRAT         *EIGENSTRAT = new CLASS_EIGENSTRAT();
CLASS_STRUCTURE          *STRUCTURE = new CLASS_STRUCTURE();
CLASS_PSEQ               *TO_PSEQ = new CLASS_PSEQ();

/*
extern CLASS_XXXXX              *XXXXX = new CLASS_XXXXX();
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
    { "Mega2 annotated format",     MEGA2ANNOT },
    { "IQLS/Idcoefs format",        IQLS },
    { "FBAT format",                FBAT },
    { "PANGAEA MORGAN format",      PANGAEA },
    { "Beagle format",              BEAGLE },
    { "Eigenstrat format",          EIGENSTRAT },
    { "Structure format",           STRUCTURE },
    { "PSEQ format",                TO_PSEQ }
};

int count_analysis_list = sizeof(analysis_list) / sizeof (analysis_types);

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

void prog_name_to_num(char *prog_name, analysis_type *analysis)
{
    *analysis = (analysis_type )0;
    switch(tolower((unsigned char)prog_name[0])) {
    case 's':
        /* Simwalk2, Summary, Sage 3 and 4, Simulate, splink, slink, solar,
           sup, structure */

        switch(tolower((unsigned char)prog_name[1])) {  // sX...

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
        *analysis = IQLS; // 34. IQLS/Idcoefs
        break;

    case 'f':
        *analysis = FBAT; // 35. FBAT
        break;

    case 'e':
        *analysis = EIGENSTRAT; // 38. EIGENSTRAT
        break;

    default:
        unknown_prog(prog_name);
    }
}

#include <string.h>

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "output_file_names_ext.h"

#include "aspex_ext.h"
#include "cranefoot_ext.h"
#include "create_summary_ext.h"
#include "hwe_user_input_ext.h"
#include "makenucs_ext.h"
#include "mega2annot_ext.h"
#include "scripts_ext.h"
#include "slink_ext.h"
#include "splink_ext.h"
#include "write_IQLS_ext.h"
#include "write_SUP_files_ext.h"
#include "write_apm_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
#include "write_loki_ext.h"
#include "write_mfiles_ext.h"
#include "write_pap_ext.h"
#include "write_premakeped_ext.h"
#include "write_prest_ext.h"
#include "write_sage_files_ext.h"
#include "write_simulate_files_ext.h"
#include "write_solar_files_ext.h"
#include "write_vitesse_ext.h"


void CLASS_ALLEGRO::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
        
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_ALLEGRO::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "al_ped.%s", num);
    sprintf(file_names[1], "al_dat.%s", num);
    sprintf(file_names[2], "al_in.%s", num);
    sprintf(file_names[3], "allegro.%s.sh", num);
}

void CLASS_ALLEGRO::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_ALLELE_FREQ::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_ALLELE_FREQ::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "freq.%s", num);
}

void CLASS_ALLELE_FREQ::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_APM::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    char *pedfl_name = NULL, *locusfl_name = NULL, *omitfl_name = NULL; // assignment to supress warnings
    int disease_locus = 0;
    create_APM_file(pedfl_name, locusfl_name, disease_locus,
                    "2", LPedTreeTop, PedTreeTop, &infl_type, &outfl_type,
                    numchr, mapfl_name, file_names,
                    analysis, omitfl_name, UntypedPedOpt);
}

void CLASS_APM::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "kin_ml.%s", num);
}

void CLASS_APM::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_APM_MULT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    char *pedfl_name = NULL, *locusfl_name = NULL; // assignment to supress warnings
    int disease_locus = 0;

    create_APMULT(pedfl_name, locusfl_name, disease_locus,
                  "2", &LPedTreeTop, PedTreeTop, &infl_type, &outfl_type,
                  numchr, mapfl_name, analysis, file_names,
                  UntypedPedOpt);
}

void CLASS_APM_MULT::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "kin_mult");
    sprintf(file_names[3], "apmmult.%2s.sh", num);
    sprintf(file_names[14], "apmmult_sum.%2s", num);
}

void CLASS_APM_MULT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[14], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_ASPEX::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "Sib-ibd");      break;
    case 2:  strcpy(subprog, "Sib-tdt");      break;
    case 3:  strcpy(subprog, "Sib-phase");    break;
    case 4:  strcpy(subprog, "Sib-map");      break;
    default:                                  break;
    }
}

void CLASS_ASPEX::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (*sub_prog_name == 0) {
        (*analysis)->_suboption = 0;
        return;
    }
    switch(tolower((unsigned char)sub_prog_name[strlen(sub_prog_name) - 1])) {
    case 'd':
        (*analysis)->_suboption = 1; break;
    case 't':
        (*analysis)->_suboption = 2; break;
    case 'e':
        (*analysis)->_suboption = 3; break;
    case 'p':
        (*analysis)->_suboption = 4; break;
    default :
        (*analysis)->_suboption = 0; break;
    }
}

void CLASS_ASPEX::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_aspex_files(&LPedTreeTop,
                      *numchr, file_names, UntypedPedOpt, Top2);
}

void CLASS_ASPEX::file_names(char **file_names, char *num) {
    sprintf(file_names[4], "asp_in.%s", num);
    sprintf(file_names[5], "asp_dat.%s", num);
}

void CLASS_ASPEX::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_CRANEFOOT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_CRANEFOOT_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_CRANEFOOT::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "crnft_ped.%s", num);
    sprintf(file_names[1], "crnft_control.%s", num);
    sprintf(file_names[2], "crnft_shell.%s.sh", num);
}

void CLASS_CRANEFOOT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_CREATE_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_CREATE_SUMMARY::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "seg_sum.%s", num);
    sprintf(file_names[1], "cnt_sum.%s", num);
    sprintf(file_names[2], "sib_sum.%s", num);
}

void CLASS_CREATE_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
}


void CLASS_GENEHUNTER::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GENEHUNTER::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "gh_ped.%s", num);
    sprintf(file_names[1], "gh_dat.%s", num);
    sprintf(file_names[2], "gh_in.%s", num);
    sprintf(file_names[3], "gh.%s.sh", num);
}

void CLASS_GENEHUNTER::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_GENEHUNTERPLUS::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GENEHUNTERPLUS::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "ghp_ped.%s", num);
    sprintf(file_names[1], "ghp_dat.%s", num);
    sprintf(file_names[2], "ghp_in.%s", num);
    sprintf(file_names[3], "ghp.%s.sh", num);
}

void CLASS_GENEHUNTERPLUS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_GENOTYPING_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_GENOTYPING_SUMMARY::file_names(char **file_names, char *num) {
    strcpy(file_names[0], "genotyping_rate");
}

void CLASS_GENOTYPING_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_GHMLB::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    
    create_gh_file(&LPedTreeTop, &infl_type, &outfl_type, numchr,
                   analysis, file_names, UntypedPedOpt, Top2);
}

void CLASS_GHMLB::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "mlb_ped.%s", num);
    sprintf(file_names[1], "mlb_dat.%s", num);
    sprintf(file_names[2], "mlb_in.%s", num);
    sprintf(file_names[3], "mlb.%s.sh", num);
}

void CLASS_GHMLB::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_HAPLOTYPE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_HAPLOTYPE::file_names(char **file_names, char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_haplo.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_HAPLOTYPE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}


void CLASS_HWETEST::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "Gen");      break;
    case 2:  strcpy(subprog, "HWE");      break;
    case 3:  strcpy(subprog, "Chi-sq");   break;
    case 4:  strcpy(subprog, "Exact");    break;
    case 5:  strcpy(subprog, "Mendel");   break;
    default:                              break;
    }
}

void CLASS_HWETEST::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'g':
        (*analysis)->_suboption = 1; break;
    case 'h':
        (*analysis)->_suboption = 2; break;
    case 'c':
        (*analysis)->_suboption = 3; break;
    case 'e':
        (*analysis)->_suboption = 4; break;
    case 'm':
        (*analysis)->_suboption = 5; break;
    default:
        (*analysis)->_suboption = 0; break;
    }
}

void CLASS_HWETEST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    hwe_user_input(LPedTreeTop, numchr, file_names);
}

void CLASS_HWETEST::file_names(char **file_names, char *num) {
}

void CLASS_HWETEST::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_IBD_EST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_IBD_EST::file_names(char **file_names, char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_ibd.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}


void CLASS_IBD_EST::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_IQLS::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    int format_checksum = 1; //XX
    if (format_checksum == 0) {
        errorvf("Pedigree, names and map file appear to be in LINKAGE format.\n");
        mssgvf("IQLS requires the input data be in annotated format.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    create_IQLS_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_IQLS::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "IQLS_pedigree.%s", num);
    sprintf(file_names[1], "IQLS_marker.%s", num);
    sprintf(file_names[2], "IQLS_parameter.%s", num);
    sprintf(file_names[3], "IQLS.%s.sh", num);
    sprintf(file_names[4], "Idcoefs_pedigree.%s", num);
    sprintf(file_names[5], "Idcoefs_study.%s", num);
}

void CLASS_IQLS::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_LIABLE_FREQ::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_LIABLE_FREQ::file_names(char **file_names, char *num) {
    strcpy(file_names[0], "liabililty_freq");
}

void CLASS_LIABLE_FREQ::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
}

void CLASS_LINKAGE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_linkage_files(&LPedTreeTop, numchr, file_names, UntypedPedOpt);
}

void CLASS_LINKAGE::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "Lpedin.%s", num);
    sprintf(file_names[1], "Ldatain.%s", num);
}

void CLASS_LINKAGE::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_LOCATION::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_LOCATION::file_names(char **file_names, char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_loc.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_LOCATION::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}
void CLASS_LOD2::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_LOD2_format_files(&LPedTreeTop, PedTreeTop,
                             analysis, numchr, file_names,
                             UntypedPedOpt);
}

void CLASS_LOD2::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "lod2ped.%s", num);
    sprintf(file_names[1], "lod2data.%s", num);
    sprintf(file_names[3], "lod2.%s.sh", num);
}

void CLASS_LOD2::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_LOKI::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_LOKI_files(&LPedTreeTop, numchr, &infl_type, file_names,
                      UntypedPedOpt, *analysis);
}

void CLASS_LOKI::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "Loki_ped.%s", num);
    sprintf(file_names[1], "Loki_freq.%s", num);
    sprintf(file_names[2], "Loki_map.%s", num);
    sprintf(file_names[3], "Loki_locus.%s", num);
    sprintf(file_names[4], "Loki_link.%s", num);
    sprintf(file_names[5], "Loki_control.%s", num);
    sprintf(file_names[6], "Loki_param.%s", num);
    sprintf(file_names[7], "Loki.%s.sh", num);
}

void CLASS_LOKI::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
}

void CLASS_MEGA2ANNOT::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_mega2annot_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_MEGA2ANNOT::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "pedin.%s.mega2", num);
    sprintf(file_names[1], "names.%s.mega2", num);
    sprintf(file_names[2], "map.%s.mega2", num);
    sprintf(file_names[3], "frequency.%s.mega2", num);
    sprintf(file_names[4], "penetrance.%s.mega2", num);
}

void CLASS_MEGA2ANNOT::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_MENDEL4::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_mega2annot_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_MENDEL4::file_names(char **file_names, char *num) {
    sprintf(file_names[4], "mendel_map.%s", num);
    sprintf(file_names[0], "mendel_ped.%s", num);
    sprintf(file_names[1], "mendel_locus.%s", num);
    sprintf(file_names[6], "mendel_pen.%s", num);
    sprintf(file_names[3], "mendel_control.%s", num);
    sprintf(file_names[7], "mendel_variable.%s", num);
}

void CLASS_MENDEL4::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[7], numchr);
}

void CLASS_MENDEL7_CSV::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MENDEL7_CSV::file_names(char **file_names, char *num) {
    sprintf(file_names[4], "mendel_map.%s", num);
    sprintf(file_names[0], "mendel_ped.%s", num);
    sprintf(file_names[1], "mendel_locus.%s", num);
    sprintf(file_names[6], "mendel_pen.%s", num);
    sprintf(file_names[3], "mendel_control.%s", num);
    sprintf(file_names[7], "mendel_variable.%s", num);
}

void CLASS_MENDEL7_CSV::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[7], numchr);
}

void CLASS_MENDEL::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MENDEL::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "pedm.%s", num);
    sprintf(file_names[1], "locus.%s", num);
    sprintf(file_names[6], "pen.%s", num);
    sprintf(file_names[2], "batch.%s", num);
    sprintf(file_names[3], "m13bat.%s", num);
}

void CLASS_MENDEL::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_MERLIN::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_merlin_files(&LPedTreeTop, numchr, &infl_type,
                        file_names, UntypedPedOpt, *analysis);
    mssgf("------------------------------------");
    create_mendel_file(&LPedTreeTop, mapfl_name,
                       PedTreeTop, &infl_type, &outfl_type,
                       numchr, analysis, file_names,
                       UntypedPedOpt);
}

void CLASS_MERLIN::file_names(char **file_names, char *num) {
    /* set simwalk2 file names */
    sprintf(file_names[0], "PEDIGREE.%s", num);
    sprintf(file_names[1], "LOCUS.%s", num);
    sprintf(file_names[2], "BATCH2.%s", num);
    sprintf(file_names[3], "npl.%s.sh", num);
    sprintf(file_names[6], "PEN.%s", num);
    sprintf(file_names[4], "SW2_MAP.%s", num);
    /* Merlin file names */
    sprintf(file_names[9], "sw2merlin_freq.%s", num);
    sprintf(file_names[10], "sw2merlin_ped.%s", num);
    sprintf(file_names[11], "sw2merlin_data.%s", num);
    sprintf(file_names[12], "sw2merlin_map.%s", num);
    sprintf(file_names[13], "sw2merlin_order.%s", num);
    strcpy(file_names[14], "merlin2sw2.pl");
}

void CLASS_MERLIN::replace_chr_number(char *file_names[], int numchr) {
    /* simwalk2 files */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[4], numchr);
    /* Merlin files */
    change_output_chr(file_names[9], numchr);
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[12], numchr);
    change_output_chr(file_names[13], numchr);
}

void CLASS_MERLINONLY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_merlin_files(&LPedTreeTop, numchr, &infl_type,
                        file_names, UntypedPedOpt, *analysis);
}

void CLASS_MERLINONLY::file_names(char **file_names, char *num) {
    sprintf(file_names[9], "merlin_freq.%s", num);
    sprintf(file_names[10], "merlin_ped.%s", num);
    sprintf(file_names[11], "merlin_data.%s", num);
    sprintf(file_names[12], "merlin_map.%s", num);
    sprintf(file_names[13], "merlin.%s.sh", num);
    sprintf(file_names[14], "merlin_model");
}

void CLASS_MERLINONLY::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[9], numchr);
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[12], numchr);
    change_output_chr(file_names[13], numchr);
}

void CLASS_MISTYPING::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_MISTYPING::file_names(char **file_names, char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_mis.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_MISTYPING::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_NONPARAMETRIC::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    (void)create_mendel_file(&LPedTreeTop, mapfl_name,
                             PedTreeTop, &infl_type, &outfl_type,
                             numchr, analysis, file_names,
                             UntypedPedOpt);
}

void CLASS_NONPARAMETRIC::file_names(char **file_names, char *num) {
       sprintf(file_names[0], "sw2_pedigree.%s", num);
       sprintf(file_names[1], "sw2_locus.%s", num);
       sprintf(file_names[2], "sw2_batch.%s", num);
       sprintf(file_names[3], "sw2_npl.%s.sh", num);
       sprintf(file_names[6], "sw2_pen.%s", num);
       sprintf(file_names[4], "sw2_map.%s", num);
}

void CLASS_NONPARAMETRIC::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[6], numchr);
}

void CLASS_NUKE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_nuclear_families(&LPedTreeTop, file_names[0], file_names[1],
                            analysis, &infl_type, UntypedPedOpt, Top2);
}

void CLASS_NUKE::file_names(char **file_names, char *num) {
    sprintf(file_names[1], "nuke_data.%s", num);
    sprintf(file_names[0], "nuke_ped.%s", num);
}

void CLASS_NUKE::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_PAP::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_pap_files(&LPedTreeTop, file_names, UntypedPedOpt, numchr, &infl_type);
}

void CLASS_PAP::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "trip.%s", num);
    sprintf(file_names[1], "header.%s", num);
    sprintf(file_names[2], "phen.%s", num);
    sprintf(file_names[3], "popln.%s", num);
    sprintf(file_names[4], "pap.%s.sh", num);
}

void CLASS_PAP::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_PREMAKEPED::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_premakeped_files(&LPedTreeTop, numchr, &infl_type,
                            file_names, UntypedPedOpt);
}

void CLASS_PREMAKEPED::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "Ppedin.%s", num);
    sprintf(file_names[1], "Pdatain.%s", num);
}

void CLASS_PREMAKEPED::replace_chr_number(char *file_names[], int numchr) {
    /* file_nums =2, 0 and 1 */
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
}

void CLASS_PREST::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_prest_files(LPedTreeTop, *analysis,
                       file_names, UntypedPedOpt, numchr,
                       &infl_type);
}

void CLASS_PREST::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "prest_ped.%s", num);
    sprintf(file_names[1], "prest_loc.%s", num);
    sprintf(file_names[2], "prest_geno.%s", num);
    sprintf(file_names[3], "prest_chrom.%s", num);
    sprintf(file_names[4], "prest.%s.sh", num);
}

void CLASS_PREST::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_QUANT_SUMMARY::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_summary_file(LPedTreeTop, analysis, file_names,
                        UntypedPedOpt, numchr);
}

void CLASS_QUANT_SUMMARY::file_names(char **file_names, char *num) {
    strcpy(file_names[0], "phenotyping_rate");
}

void CLASS_QUANT_SUMMARY::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_SAGE4::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SAGE4_file(&LPedTreeTop, PedTreeTop, &infl_type,
                      numchr, analysis, file_names, UntypedPedOpt);
}

void CLASS_SAGE4::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "sage4_ped.%s", num);
    sprintf(file_names[1], "sage4_dat.%s", num);
    sprintf(file_names[6], "sage4_par.%s", num);
    sprintf(file_names[7], "sage4_map.%s", num);
}

void CLASS_SAGE4::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
}

void CLASS_SAGE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;
    char *mapfl_name = NULL; // assignment to supress warning
    int trait_locus_first;

    create_SAGE_file(&LPedTreeTop, mapfl_name, &trait_locus_first,
                     PedTreeTop, &infl_type, &outfl_type,
                     numchr, analysis, file_names, UntypedPedOpt);
}

void CLASS_SAGE::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "sage_ped.%s", num);
    sprintf(file_names[1], "sage_loc.%s", num);
    sprintf(file_names[2], "sage_par.%s", num);
    sprintf(file_names[3], "sage.%s.sh", num);
    sprintf(file_names[6], "sage_cnt.%s", num);
    sprintf(file_names[7], "sage_sibpal.%s", num);
    sprintf(file_names[10], "sage_cntpar.%s", num);
    sprintf(file_names[11], "sage_cnt.%s.sh", num);
}

void CLASS_SAGE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[10], numchr);
    change_output_chr(file_names[11], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[0], numchr);
}

void CLASS_SIMULATE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SIMULATE_format_files(&LPedTreeTop, PedTreeTop, analysis,
                                 &infl_type, &outfl_type,
                                 numchr, mapfl_name, file_names,
                                 UntypedPedOpt);
}

void CLASS_SIMULATE::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "simped.%s", num);
    sprintf(file_names[1], "simdata.%s", num);
    sprintf(file_names[2], "problem.%s", num);
    sprintf(file_names[3], "simulate.%s.sh", num);
}

void CLASS_SIMULATE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_SIMWALK2::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int  selection;
    char select[10];

    if (batchANALYSIS) {
        selection = (*analysis)->_suboption;
    } else {
        selection=0;
        while (selection < 1 || selection > 5) {
            draw_line();
            printf("SimWalk2 program options:\n");
            printf("1) Haplotype analysis\n");
            printf("2) Parametric Linkage analysis\n");
            printf("3) Non-Parametric Linkage analysis\n");
            printf("4) IBD analysis\n");
            printf("5) Mistyping analysis\n");
            printf("Enter selection 1-5 > ");
            fcmap(stdin, "%s", select); newline;
            selection = 0;
            sscanf(select, "%d", &selection);
            if (selection < 1 || selection > 5)
                warn_unknown(select);
        }
        draw_line();
    }
    switch (selection) {
    case 1:
        *analysis = HAPLOTYPE;
        break;
    case 2:
        *analysis = LOCATION;
        break;
    case 3:
        *analysis = NONPARAMETRIC;
        break;
    case 4:
        *analysis = IBD_EST;
        break;
    case 5:
        *analysis = MISTYPING;
        break;
    }
}

void CLASS_SIMWALK2::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'h':
        *analysis = HAPLOTYPE;
        break;
    case 'i':
        *analysis = IBD_EST;
        break;
    case 'n':
        *analysis = NONPARAMETRIC;
        break;
    case 'p':
        *analysis = LOCATION;
        break;
    case 'm':
        *analysis = MISTYPING;
        break;
    default:
        *analysis = NULL;
        break;
    }
}

void CLASS_SLINK::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SLINK_format_files(LPedTreeTop, PedTreeTop, analysis,
                              numchr, file_names,
                              UntypedPedOpt);
}

void CLASS_SLINK::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "simped.%s", num);
    sprintf(file_names[1], "simdata.%s", num);
    sprintf(file_names[2], "slinkin.%s", num);
    sprintf(file_names[3], "slink.%s.sh", num);
}

void CLASS_SLINK::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
}

void CLASS_SOLAR::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {
    infl_type = LPedTreeTop->pedfile_type == POSTMAKEPED_PFT ? LINKAGE : PREMAKEPED;

    create_SOLAR_files(&LPedTreeTop, PedTreeTop, analysis,
                       &infl_type, &outfl_type,
                       numchr, mapfl_name, file_names, UntypedPedOpt);
}

void CLASS_SOLAR::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "solar_ped.%s", num);
    sprintf(file_names[1], "solar_freq.%s", num);
    sprintf(file_names[6], "solar_phen.%s", num);
    sprintf(file_names[8], "solar_marker.%s", num);
    sprintf(file_names[9], "solar_map.%s", num);
    sprintf(file_names[10],"solar_load.%s.tcl", num);
}

void CLASS_SOLAR::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[8], numchr);
    change_output_chr(file_names[9], numchr);
    /*    change_output_chr(file_names[10], numchr); */
}

void CLASS_SPLINK::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SPLINK(&LPedTreeTop, analysis, numchr,
                  file_names, UntypedPedOpt);
}

void CLASS_SPLINK::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "splink_ped.%s", num);
    sprintf(file_names[3], "splink.%s.sh", num);
}

void CLASS_SPLINK::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[3], numchr);
}

void CLASS_SUP::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_SUP_files(&LPedTreeTop, file_names, UntypedPedOpt);
}

void CLASS_SUMMARY::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int  selection;
    char select[10];

    if (batchANALYSIS) {
        selection = (*analysis)->_suboption;
    } else {
        selection = 0;
        while(!(selection == 1 ||
                selection == 2 ||
                selection == 3 ||
                selection == 4 ||
                selection == 5)) {
            draw_line();
            printf("Selection Menu: Summary file options\n");
            printf("1) Create segregation and relative count summary files.\n");
            printf("2) Create allele frequency summary table.\n");
            printf("3) Count alleles and genotypes within groups.\n");
            printf("4) Create genotyping success rate summary.\n");
            printf("5) Create quantitative phenotype summary.\n");
            printf("Enter selection: 1 - 5 > ");
            fcmap(stdin,"%s", select); newline;
            selection=0;
            sscanf(select, "%d", &selection);
            if (selection < 1 || selection > 5)
                warn_unknown(select);
        }
    }

    if (selection == 1) {
        *analysis = CREATE_SUMMARY;
    } else if (selection  == 2) {
        *analysis = TO_ALLELE_FREQ;
    } else if (selection == 3 ) {
        *analysis = TO_LIABLE_FREQ;
    } else if (selection == 4 ) {
        *analysis = GENOTYPING_SUMMARY;
    } else if (selection == 5 ) {
        *analysis = QUANT_SUMMARY;
    }
}

void CLASS_SUMMARY::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    switch(tolower((unsigned char)sub_prog_name[0])) {
    case 'a':
        *analysis = TO_ALLELE_FREQ;
        break;
    case 's':
        *analysis = CREATE_SUMMARY;
        break;
    case 'l':
    case 'c':
        *analysis = TO_LIABLE_FREQ;
        break;
    case 'g':
        *analysis = GENOTYPING_SUMMARY;
        break;
    case 'q':
        *analysis = QUANT_SUMMARY;
        break;
    default:
        *analysis = NULL;
	break;
    }
}

void CLASS_SUP::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "sup_simped.%s", num);
    sprintf(file_names[1], "sup_simdata.%s", num);
    sprintf(file_names[2], "sup_locus.%s", num);
    sprintf(file_names[3], "slinkin.%s", num);
    sprintf(file_names[4], "sup.%s.sh", num);
    sprintf(file_names[5], "sup_outpeds.%s", num);
    sprintf(file_names[6], "sup_mega2_locus.%s", num);
    sprintf(file_names[7], "sup_mega2_map.%s", num);
    sprintf(file_names[8], "sup_mega2_batch.%s", num);
}

void CLASS_SUP::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
}

void CLASS_TDTMAX::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_TDTMAX(&LPedTreeTop, analysis, numchr, file_names,
                  UntypedPedOpt, Top2);
}


void CLASS_TDTMAX::file_names(char **file_names, char *num) {
    strcpy(file_names[0], "tdtmax_data");
    sprintf(file_names[3], "tdtmax.%s.sh", num);
}

void CLASS_TDTMAX::replace_chr_number(char *file_names[], int numchr) {
}

void CLASS_VITESSE::sub_prog_name(int sub_opt, char *subprog) {
    switch(sub_opt) {
    case 1:  strcpy(subprog, "LINKMAP");      break;
    case 2:  strcpy(subprog, "MLINK");        break;
    default:                                  break;
    }
}

void CLASS_VITESSE::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    if (tolower((unsigned char)sub_prog_name[0]) == 'l') {
        (*analysis)->_suboption = 1;
    } else if (tolower((unsigned char)sub_prog_name[0]) == 'm') {
        (*analysis)->_suboption = 2;
    }
}

void CLASS_VITESSE::create_output_file(
    linkage_ped_top *LPedTreeTop,
    analysis_type *analysis,
    char *file_names[],
    int untyped_ped_opt,
    int *numchr,
    linkage_ped_top **Top2) {

    create_vitesse_files(&LPedTreeTop, numchr, file_names, UntypedPedOpt);
}

void CLASS_VITESSE::file_names(char **file_names, char *num) {
    sprintf(file_names[0], "vpedin.%s", num);
    sprintf(file_names[1], "vdatain.%s", num);
    sprintf(file_names[3], "vitesse.%s.sh", num);
    sprintf(file_names[4], "vout.%s", num);
    sprintf(file_names[5], "vstrm.%s", num);
}

void CLASS_VITESSE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[4], numchr);
    change_output_chr(file_names[5], numchr);
}
