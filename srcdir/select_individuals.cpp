/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "mrecode.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "utils_ext.h"
/*
        batch_input_ext.h:  batchf
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
            mrecode_ext.h:  allelecnt (type)
              utils_ext.h:  draw_line log_line randomindex
*/


#define DEFAULT_SelectIndividuals 4
#define HWE_DEFAULT_SelectIndividuals 3

int get_count_option(int halftyped_item, int *include_halftyped, const char *messg);
int count_genotypes(linkage_ped_top *LPedTreeTop, int locus_id, int option,
                    int *geno_count);


/* Note : NM 8/18/05 : this function used to work on all
   pedigrees, changing it to work on a specific pedigree.
   Also implementing uncoded alleles.
   Returns the number of people counted
*/

int founders_or_everyone(linkage_ped_top *Top, int locus_id,
			 record_type rec, int ped, allelecnt *member_ids,
			 int count_option, int inc_ht, int xlinked)


{

    int jj, valid, typed1, typed2, has_typed=0;
    int entry_cnt = 0;
    linkage_ped_rec  *pr = NULL, *prp = NULL;
    person_node_type *pn = NULL, *pnp = NULL;

    /* member ids - vector, set to 1-half typed, 2- fully typed */

    /* find the founder in each pedigree */

    if (rec == Raw_postmake || rec == Postmakeped) {
        entry_cnt=Top->Ped[ped].EntryCnt;
        pr = Top->Ped[ped].Entry;
    } else if (rec == Raw_premake || rec == Premakeped || rec == Annotated) {
        entry_cnt=Top->PTop[ped].num_persons;
        pn = Top->PTop[ped].persons;
    }

    for (jj=0; jj<entry_cnt; jj++) {
        typed1 = typed2 =  valid=0;

        if (rec == Raw_postmake || rec == Postmakeped) {
            prp = &pr[jj];
            if (xlinked && prp->Sex == 1) {
                valid = 0;
            } else {
                valid=((count_option == 4) || IS_LFOUNDER(*prp));
            }
        } else if (rec == Raw_premake || rec == Premakeped || rec == Annotated) {
            pnp = &pn[jj];
            if (xlinked && pnp->gender == 1) {
                valid = 0;
            } else {
                valid=((count_option == 4) || PFOUNDER(*pnp));
            }
        }

        if (valid == 0) {
            ; /* why bother */
        } else
/*
        if (rec == Raw_premake) {
            const char *a1, *a2;
            get_2Ralleles(pnp->marker, locus_id, &a1, &a2);
            typed1=allelecmp(a1, REC_UNKNOWN);
            typed2=allelecmp(a2, REC_UNKNOWN);
        } else if (rec == Raw_postmake) {
            const char *a1, *a2;
            get_2Ralleles(prp->Marker, locus_id, &a1, &a2);
            typed1=allelecmp(a1, REC_UNKNOWN);
            typed2=allelecmp(a2, REC_UNKNOWN);
*/
        if ((rec == Raw_premake) || (rec == Raw_postmake)) {
            typed1=allelecmp(member_ids[jj].all1, REC_UNKNOWN);
            typed2=allelecmp(member_ids[jj].all2, REC_UNKNOWN);
        } else if (rec == Postmakeped) {
            int a1, a2;
            get_2alleles(prp->Marker, locus_id, &a1, &a2);
            typed1=(a1)? 1 : 0;
            typed2=(a2)? 1 : 0;
        } else if (rec == Premakeped || rec == Annotated) {
            int a1, a2;
            get_2alleles(pnp->marker, locus_id, &a1, &a2);
            typed1=(a1)? 1 : 0;
            typed2=(a2)? 1 : 0;
        }

        /* find if this person is genotyped at the locus */
        if (inc_ht) {
            if (valid && (typed1 || typed2)) {
                if (typed1 && typed2) {
                    member_ids[jj].num=2;
                } else {
                    member_ids[jj].num=1;
                }
                has_typed++;
            }
        } else {
            if (valid && typed1 && typed2) {
                member_ids[jj].num=2;
                has_typed++;
            }
        }
    }
    return has_typed;
}

/* option 2: select a random member from peds not already selected */
/* NM: 8-21-08: Since females' genotypes are set to NULL at Y-linked
   loci, we also have to check for NULL */
int random_ped_member(linkage_ped_top *LPedTreeTop,
		      int locus_id, int ped, record_type rec,
		      int inc_ht, int *num_alleles, int xlinked)

{
    int jj, sex;
    int entry_cnt, geno_entry_cnt;
    int typed1=0, typed2=0;
    int *genotyped_members, *num_all; /* num_all to store the number of alleles */
    int homo_test=1;
    int ylinked;

    ylinked = ((LPedTreeTop->LocusTop->Locus[locus_id].Type == YLINKED)? 1: 0);
    /* member ids is a vector of 0/1 - selected/not selected */

    entry_cnt=((rec == Raw_postmake || rec == Postmakeped )?
               LPedTreeTop->Ped[ped].EntryCnt :
               LPedTreeTop->PTop[ped].num_persons);

    genotyped_members = CALLOC((size_t) entry_cnt, int);
    num_all=CALLOC((size_t) entry_cnt, int);

    geno_entry_cnt=0;
    for (jj=0; jj < entry_cnt; jj++) {
        if (xlinked) {
            if (rec == Raw_premake || rec == Premakeped) {
                sex = LPedTreeTop->PTop[ped].persons[jj].gender;
            } else {
                sex = LPedTreeTop->Ped[ped].Entry[jj].Sex;
            }
            if (sex == 1) continue;
        }
        if (rec == Raw_premake) {
            const char *a1, *a2;
            get_2Ralleles(LPedTreeTop->PTop[ped].persons[jj].marker, locus_id, &a1, &a2);
            if (a1 == NULL) {
                typed1 = typed2 = 0;
            } else {
                typed1=(allelecmp(a1, REC_UNKNOWN)? 1: 0);
                typed2=(allelecmp(a2, REC_UNKNOWN)? 1: 0);
                if (ylinked && typed1 && typed2) {
                    homo_test =
                        !allelecmp(a1, a2);
                }
            }
        } else if (rec == Raw_postmake) {
            const char *a1, *a2;
            get_2Ralleles(LPedTreeTop->Ped[ped].Entry[jj].Marker, locus_id, &a1, &a2);
            if (a1 == NULL) {
                typed1=typed2=0;
            } else {
                typed1=(allelecmp(a1, REC_UNKNOWN)? 1: 0);
                typed2=(allelecmp(a2, REC_UNKNOWN)? 1: 0);
                if (ylinked && typed1 && typed2) {
                    homo_test =
                        !allelecmp(a1, a2);
                }
            }
        } else if (rec == Postmakeped) {
            int a1, a2;
            get_2alleles(LPedTreeTop->Ped[ped].Entry[jj].Marker, locus_id, &a1, &a2);
            typed1= (a1)? 1 : 0;
            typed2= (a2)? 1 : 0;
            if (ylinked && typed1 && typed2) {
                homo_test = a1 == a2;
            }
        } else {
            int a1, a2;
            get_2alleles(LPedTreeTop->PTop[ped].persons[jj].marker, locus_id, &a1, &a2);
            typed1=(a1)? 1 : 0;
            typed2=(a2)? 1 : 0;
            if (ylinked && typed1 && typed2) {
                homo_test = a1 == a2;
            }
        }

        if (typed1 && typed2 && homo_test) {
            genotyped_members[geno_entry_cnt]=jj;
            num_all[geno_entry_cnt]=2;
            geno_entry_cnt++;
        } else if (inc_ht && (typed1 || typed2) && homo_test) {
            genotyped_members[geno_entry_cnt]=jj;
            num_all[geno_entry_cnt]=typed1+typed2;
            geno_entry_cnt++;
        }

    }

    jj=-1;
    if (geno_entry_cnt > 1) {
        jj=randomindex(geno_entry_cnt);
        *num_alleles=num_all[jj];
        jj=genotyped_members[jj];
    } else if (geno_entry_cnt == 1) {
        jj=genotyped_members[0];
        *num_alleles=num_all[0];
    }

    /*  printf("%d ped %d per\n", ped, jj, num_all[); */

    free(genotyped_members); free(num_all);

    return jj;
}

/* have a menu option to select which individuals to count */
/* Added hwe select option */

int get_count_option(int halftyped_item, int *include_halftyped, const char *messg)
{
    char choice[10], stars[5];
    int i, menu_select = 0, count_type=0;

    menu_select = InputMode == INTERACTIVE_INPUTMODE;
// allow multiple calls ... only first counts
//  if (InputMode == BATCH_FILE_INPUTMODE) {
        // Default Select Individuals...
        if ( (Mega2Status == INSIDE_ANALYSIS && Mega2BatchItems[Count_HWE_genotypes].items_read) ||
              Mega2BatchItems[Count_Genotypes].items_read) {
            menu_select = 0;
            if (Mega2Status == INSIDE_ANALYSIS) {
                SelectIndividuals = Mega2BatchItems[/* 38 */ Count_HWE_genotypes].value.option;
            } else {
                SelectIndividuals = Mega2BatchItems[/* 34 */ Count_Genotypes].value.option;
            }

            if (halftyped_item && Mega2BatchItems[/* 36 */ Count_Halftyped].items_read) {
                *include_halftyped =
                    ((tolower((unsigned char)Mega2BatchItems[/* 36 */ Count_Halftyped].value.copt) == 'y')?
                     1 : 0);
                menu_select = 0;
            } else {
                *include_halftyped = 0;
                menu_select = 0;
            }
        } else {
            menu_select = 1;
        }
//  }
    if (! menu_select) {

    } else {
        if (Mega2Status == INSIDE_ANALYSIS) {
            SelectIndividuals = HWE_DEFAULT_SelectIndividuals;
        } else {
            SelectIndividuals = DEFAULT_SelectIndividuals;
        }
        count_type = SelectIndividuals;
        *include_halftyped=0;
        while (count_type){
            strcpy(stars, "    ");
            stars[SelectIndividuals-1]='*';
            draw_line();
            printf("%s\n", messg);
            draw_line();
            printf("0) Done with this menu, please proceed.\n");
            printf("%c1) Genotyped founders only\n", stars[0]);
            printf("%c2) Genotyped founders + a randomly chosen genotyped person \n",
                   stars[1]);
            printf("    from pedigrees without genotyped founders.\n");
            i=3;
// next time: choose[i++] = 1:5
// then later switch(choose[count_type]) { 1: ...; 2: ... ; 3: ...; 4: ...; 5: ...; }
            if (Mega2Status == INSIDE_RECODE) {
                printf("%c%d) Genotyped founders + genotyped individuals with unique alleles\n",
                       stars[i-1], i);
                i++;
            }
            printf("%c%d) All genotyped individuals\n", 
                   (Mega2Status == INSIDE_RECODE) ? stars[i-1] : stars[i], i);
            if (halftyped_item) {
                i++;
                printf(" %d) Count half-typed individuals' alleles. [%s]\n",
                       i, yorn[*include_halftyped]);
            }
            if (halftyped_item) {
                printf("Select from options 0 - %d or %d to toggle > ", i-1, i);
            } else {
                printf("Select from options 0 - %d > ", i);
            }
            fcmap(stdin, "%s", choice); newline;
            sscanf(choice, "%d", &count_type);
            switch(count_type) {
            case 0:
                break;
            case 1:
                SelectIndividuals = 1;
                break;
            case 2:
                SelectIndividuals = 2;
                break;
            case 3:
//              SelectIndividuals = 3;
// I believe the above is wrong if INSIDE_ANALYSIS
                if (Mega2Status == INSIDE_RECODE) {
                    SelectIndividuals = 3;
                } else {
                    SelectIndividuals = 4;
                }
                break;
            case 4:
                if (Mega2Status == INSIDE_RECODE) {
                    SelectIndividuals = 4;
                    break;
                }
            case 5:
                if (halftyped_item) {
                    *include_halftyped = TOGGLE(*include_halftyped);
                    break;
                }
            default:
                warn_unknown(choice);
            }
        }
    }
    /* log selection */
    switch(SelectIndividuals) {
    case 1:
        sprintf(err_msg, "Count option: founder alleles");
        break;
    case 2:
        sprintf(err_msg, "Count option: founder + random alleles");
        break;
    case 3:
        if (Mega2Status == INSIDE_RECODE || Mega2Status == ANALYSIS_NAME_READ) {
            sprintf(err_msg, "Count option: founder + all new alleles");
        } else {
            sprintf(err_msg, "Count option: all individuals");
        }
        break;
    case 4:
        if (Mega2Status == INSIDE_RECODE || Mega2Status == ANALYSIS_NAME_READ) {
            sprintf(err_msg, "Count option: all individuals");
        }
        break;
    default:
        break;
    }
    mssgf(err_msg);
    if (halftyped_item) {
        sprintf(err_msg, "Count half-typed individuals' alleles : %s",
                yorn[*include_halftyped]);
        mssgf(err_msg);
    }
    log_line(mssgf);

    if (InputMode == INTERACTIVE_INPUTMODE) {
        if (Mega2Status == INSIDE_RECODE || Mega2Status == ANALYSIS_NAME_READ) {
            BatchValueSet(SelectIndividuals, "Count_Genotypes");
            batchf(Count_Genotypes);
        } else {
            BatchValueSet(SelectIndividuals, "Count_HWE_Genotypes");
            batchf(Count_HWE_genotypes);
        }
        if (halftyped_item) {
            char c = yorn[*include_halftyped][0];
            BatchValueSet(c, "Count_Halftyped");
            batchf(Count_Halftyped);
        }
    }

    return SelectIndividuals;
}

void select_individuals(linkage_ped_top *LPedTreeTop, int locus_id,
			int count_option, int ped, allelecnt *member_ids,
			int inc_ht, int xlinked)
{

    int found,  ind;
    int num_alleles;
    found=founders_or_everyone(LPedTreeTop, locus_id,
                               Postmakeped, ped, member_ids,
                               count_option, inc_ht, xlinked);

    if (!found && count_option == 2) {
        ind=random_ped_member(LPedTreeTop, locus_id, ped, Postmakeped,
                              inc_ht, &num_alleles, xlinked);
        if (ind > -1) {
            member_ids[ind].num=num_alleles;
        }
    }

    return;
}
