/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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


/* This is a try at creating nuclear families
   by looking at the pedtree structure.. In this
   version of mega2, LRec->ID is mapped to PRec->ID,
   and the assumption is that numbering of entries is
   consecutive.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "pedtree_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
#include "append_locus_array_ext.h"
/*
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  connect_loops get_loci_on_chromosome get_unmapped_loci new_lpedtop switch_penetrances
           makenucs_ext.h:  create_nuked_fids
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg print_outfile_mssg
    output_routines_ext.h:  ped_id_width ped_name_width person_id_width
            pedtree_ext.h:  convert_to_pedtree copy_le_static_data free_all_including_ped_top
         user_input_ext.h:  pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line global_ou_files
        write_files_ext.h:  write_affection_data write_binary_data write_numbered_data write_quantitative_data
      write_ghfiles_ext.h:  write_gh_locus_file
*/


static void     set_nuked_pednames(linkage_ped_top *Topp);

typedef struct _marriage_list_item {
    int father, mother, num_off;
    struct _marriage_list_item *next;
} marriage_list_item;

marriage_list_item *paste_marriage_lists(marriage_list_item *list1, marriage_list_item *list2);
marriage_list_item *append_to_marriage_list(marriage_list_item *marriages, int father, int mother);


/*--------------------------------------------------------+
  | Write the pedigree data file for nuclear families      |
  +--------------------------------------------------------*/
static void  write_nuclear_peds(char *outfl, linkage_ped_top *Top2,
				int numchr, int linkage_format)

{
    int             ped, entry1, locus;
    linkage_ped_rec *Entry;
    int             *trp, tr, nloop, num_affec = num_traits;
    char            outfl1[2*FILENAME_LENGTH], fformat[6];
    FILE            *filep;
    int             pwid,  fwid;

    NLOOP;
    trp = &(global_trait_entries[0]);

    if (OrigIds[1] == 4 || OrigIds[1] == 5) {
        create_nuked_fids(Top2);
        if (OrigIds[1] == 4) {
            /* extensions, Ped->Name */
            ped_name_width(Top2, &fwid);
        } else {
            /* multipliers, Ped->Num */
            ped_id_width(Top2, &fwid);
        }
    } else {
        /* Has to be done manually because Ped fields are not used */
        fwid = (int) floor(log10((double)Top2->PedCnt)) + 1;
    }

    /* The person id is the ID field */
    person_id_width(Top2, &pwid);

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;
        sprintf(outfl1, "%s/%s", output_paths[tr], outfl);
        if ((filep = fopen(outfl1, "w")) == NULL) {
            errorvf("Can't write nuclear pedigree file %s\n", outfl1);
            EXIT(FILE_WRITE_ERROR);
        }

        for (ped = 0; ped < Top2->PedCnt; ped++) {
            if (UntypedPeds != NULL) {
                if (UntypedPeds[ped]) {
                    continue;
                }
            }
            for (entry1 = 0; entry1 < Top2->Ped[ped].EntryCnt; entry1++)  {
                Entry = &(Top2->Ped[ped].Entry[entry1]);
                if (OrigIds[1] == 3) {
                    Top2->Ped[ped].Num = ped+1;
                    sprintf(fformat, "%%%dd ", fwid);
                    fprintf(filep, fformat, ped+1);
                } else if (linkage_format || OrigIds[1] == 5) {
                        sprintf(fformat, "%%%dd ", fwid);
                        fprintf(filep, fformat, Top2->Ped[ped].Num);
                } else if (OrigIds[1] == 2) {
                        sprintf(fformat, "%%%ds ", fwid);
                        fprintf(filep, fformat, Top2->Ped[ped].Name);
                } else if (OrigIds[1] == 4) {
                        sprintf(fformat, "%%%ds ", fwid);
                        fprintf(filep, fformat, Top2->Ped[ped].Name);
                } else if (OrigIds[1] == 6) {
                        sprintf(fformat, "%%%ds ", fwid);
                        fprintf(filep, fformat, Top2->Ped[ped].PedPre);
                }

                sprintf(fformat, "%%%dd ", pwid);
                fprintf(filep, fformat, Entry->ID);
                /* write the parents */
                fprintf(filep, fformat, Entry->Father);
                fprintf(filep, fformat, Entry->Mother);
                if (linkage_format) {
                    /* write the first sib */
                    fprintf(filep, fformat, Entry->First_Offspring);

                    /* write the next sibs */
                    if (Entry->Father == 0)  fprintf(filep, fformat, 0);
                    else fprintf(filep, fformat, Entry->Next_PA_Sib);
                    if (Entry->Mother == 0)  fprintf(filep, fformat, 0);
                    else fprintf(filep, fformat, Entry->Next_MA_Sib);
                }
                /* write the sex */
                if (Entry->Sex == MALE_ID)
                    fprintf(filep, " %1d", MALE_ID);
                else if (Entry->Sex == FEMALE_ID)
                    fprintf(filep, " %1d", FEMALE_ID);
                else
                    fprintf(filep, " 0");
                /* write the proband */
                if (linkage_format) {
                    fprintf(filep, " %1d ", Entry->OrigProband);
                }
                if (LoopOverTrait == 1) {
                    /* write the trait locus first */
                    locus = *trp;
                    switch (Top2->LocusTop->Locus[locus].Type)    {
                    case QUANT:
                        // This routine will write the quantitative data or the appropriate
                        // missing quantitative data symbol or value if appropriate...
                        write_quantitative_data(filep, locus,
                                                &(Top2->LocusTop->Locus[locus]),
                                                Entry);
                        break;
                    case AFFECTION:
                        write_affection_data(filep, locus,
                                             &(Top2->LocusTop->Locus[locus]),
					     Entry);
                        break;
                    default:
                        /* should not happen */
                        break;
                    }
                }
                /* now write the marker data */
                for (locus = 0; locus < NumChrLoci; locus++)   {
                    switch (Top2->LocusTop->Locus[ChrLoci[locus]].Type)    {
                    case QUANT:
                        if (LoopOverTrait == 0) {
                            // This routine will write the quantitative data or the appropriate
                            // missing quantitative data symbol or value if appropriate...
                            write_quantitative_data(filep, ChrLoci[locus],
                                                    &(Top2->LocusTop->Locus[ChrLoci[locus]]),
                                                    Entry);
                        }
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            write_affection_data(filep, ChrLoci[locus],
                                                 &(Top2->LocusTop->Locus[ChrLoci[locus]]),
                                                 Entry);
                        }
                        break;
                    case NUMBERED:
                        if (numchr == 0 ||
                            Top2->LocusTop->Marker[ChrLoci[locus]].chromosome != numchr)
                            write_numbered_data(filep, ChrLoci[locus],
                                                &(Top2->LocusTop->Locus[ChrLoci[locus]]),
                                                Entry);
                        break;
                    case BINARY:
                        if (numchr == 0 ||
                            Top2->LocusTop->Marker[ChrLoci[locus]].chromosome != numchr)
                            write_binary_data(filep, ChrLoci[locus],
                                              &(Top2->LocusTop->Locus[ChrLoci[locus]]),
                                              Entry);
                        break;
                    default:
                        errorvf("unknown locus type (locus %d) in nuclear peds.\n",
                                locus + 1);
                        EXIT(DATA_TYPE_ERROR);
                    }
                }
                if (linkage_format) {
                    fprintf(filep, "  Ped: %s Per: %s",
                            Top2->Ped[ped].Name, Entry->OrigID);
                }
                if (Top2->UniqueIds) {
                    fprintf(filep, " ID: %s", Entry->UniqueID);
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);   trp++;
        if (nloop == 1) break;
    }
}

static void nuclear_family_files(char *outfl_name, char *loutfl_name,
				 int single_output, int *linkage_format)
{
    int reply = -1;
    char creply[10];

    *linkage_format=1;
    if (DEFAULT_OUTFILES) {
        draw_line();
        mssgf("Using post-makeped format");
        mssgf("Pedigrees will be numbered consecutively.");
        draw_line();
        return;
    }

    printf("Note: Renaming of pedigrees can be either\n");
    printf("      consecutive renumbering (1 - #nuclear-families)\n");
    printf("      or original names with extensions e.g.\n");
    printf("          100_1, 100_2, 100_3 \n");
    printf("      or numeric with multipliers e.g\n");
    printf("          100001, 100002, 100003 \n");
    printf("      where original pedigree 100 is broken into\n");
    printf("      3 nuclear pedigrees.\n");
    while (reply != 0) {
        draw_line();
        print_outfile_mssg();
        printf("Nuclear-pedigrees file name menu:\n");
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt == 1) {
            printf(" 1) Locus file name:            %-15s\t[%s]\n",
                   loutfl_name,
                   (access(loutfl_name, F_OK) == 0) ? "overwrite" : "new");
            printf(" 2) Pedigree file name:         %-15s\t[%s]\n", outfl_name,
                   (access(outfl_name, F_OK) == 0) ? "overwrite" : "new");
        } else {
            if (single_output == 1) {
                change_output_chr(loutfl_name, 0);
                change_output_chr(outfl_name, 0);
                printf(" 1) Locus file name:            %-15s\t[%s]\n", loutfl_name,
                       (access(loutfl_name, F_OK) == 0) ? "overwrite" : "new");
                printf(" 2) Pedigree file name:         %-15s\t[%s]\n", outfl_name,
                       (access(outfl_name, F_OK) == 0) ? "overwrite" : "new");
            } else {
                change_output_chr(loutfl_name, -9);
                change_output_chr(outfl_name, -9);
                printf(" 1) Locus file name stem:       %-15s\n", loutfl_name);
                printf(" 2) Pedigree file name stem:    %-15s\n", outfl_name);
            }
        }
        pedigree_id_item(3, TO_NUKE, OrigIds[1], 33, 2, 0);
        printf(" 4) Post-makeped pedigrees:     [%s]\n", yorn[*linkage_format]);
        printf("Enter selection: 0 - 4 (3 and 4 to toggle) > ");
        fcmap(stdin, "%s", creply); newline;
        reply=-1;
        sscanf(creply, "%d", &reply);
        test_modified(reply);

        switch (reply) {
        case 1:
            printf("Enter new locus datafile: ");
            fcmap(stdin, "%s", loutfl_name); newline;
            break;
        case 2:
            printf("Enter new pedigree datafile: ");
            fcmap(stdin, "%s", outfl_name); newline;
            break;
        case 0:
            break;
        case 3:
            OrigIds[1] = pedigree_id_item(0, TO_NUKE, OrigIds[1], 0, 1, 0);
            break;
        case 4:
            *linkage_format = TOGGLE(*linkage_format);
            break;
        default:
            warn_unknown(creply);
            break;
        }
    }
    draw_line();
    ((*linkage_format==1)?
     mssgf("Using post-makeped format") :
     mssgf("Using pre-makeped format"));
    pedigree_id_item(0, TO_NUKE, OrigIds[1], 0, 3, 0);
    change_output_chr(outfl_name, global_chromo_entries[0]);
    change_output_chr(loutfl_name, global_chromo_entries[0]);
    draw_line();
}

static void free_marriage_list(marriage_list_item *list)

{
    if (list == NULL)
        return;
    else if (list->next == NULL) {
        return;
    } else {
        free_marriage_list(list->next);
        list->next = NULL;
        free(list);
        list=NULL;
    }
}

/* check parents - check if a parent pair already exists within list,
   if so, then increment that marriage's offspring count by 1 */

static int check_parents(int father, int mother,
			 marriage_list_item *list)

{
    if (list == NULL) return 1;

    while (list != NULL) {
        if (father == 0 || mother == 0) {
            list = list->next;
        } else if (list->father == father &&
                list->mother == mother) {
            (list->num_off)++;
            return 0;
        } else {
            list = list->next;
        }
    }
    return 1;
}

/* paste list2 to list 1, and return a pointer to the combined
   list */

marriage_list_item *paste_marriage_lists(marriage_list_item *list1,
                                         marriage_list_item *list2)

{
    marriage_list_item *itemp = list1;
    if (itemp == NULL) {
        return(list2);
    } else {
        while(itemp->next != NULL) {
            itemp = itemp->next;
        }
        itemp->next = list2;
    }
    return list1;
}

/* add a new marriage to the end of the marriage list,
   given father and mother's IDs, set offspring count to 1 */

marriage_list_item *append_to_marriage_list(marriage_list_item *marriages,
                                            int father, int mother)

{
    marriage_list_item *next_pair = marriages;
    marriage_list_item *new_pair;

    new_pair = CALLOC((size_t) 1, marriage_list_item);
    new_pair->father=father;
    new_pair->mother=mother;
    new_pair->next = NULL;

    if (father > 0 && mother > 0) {
        new_pair->num_off = 1;
    } else {
        new_pair->num_off = 0;
    }
    if (marriages == NULL) {
        marriages = new_pair;
    } else {
        while(next_pair->next != NULL) {
            next_pair = next_pair->next;
        }
        next_pair->next = new_pair;
    }
    return marriages;
}

void makenucs1(linkage_ped_top *Top, ped_top *PedTop,
	       linkage_ped_top *NewTop)

{

    int ped, entry, pr, newentry, father, mother;
    int *num_marriages;
    int created = 0, nuke_pednum=0;
    marriage_list_item *marriages, *all_marriages, *marriagep;
    linkage_ped_tree *To;


    if (PedTop==NULL) {
        PedTop=convert_to_pedtree(Top, 0);
        created=1;
    }

    /* first copy over the locus data */
    NewTop->LocusTop = Top->LocusTop;
    NewTop->UniqueIds = 0;
    NewTop->OrigIds = 0;
    NewTop->pedfile_type = POSTMAKEPED_PFT;

    num_marriages=CALLOC((size_t) Top->PedCnt, int);
    /* Now form the nuclear families */

    all_marriages=NULL;
    NewTop->PedCnt=0;
    for(ped=0; ped<PedTop->PedCnt; ped++) {
        marriages=NULL;
        num_marriages[ped]=0;

/*     if (UntypedPeds != NULL) { */
/*       if (UntypedPeds[ped]) { */
/* 	continue; */
/*       } */
/*     } */

        if (PedTop->PedTree[ped].EntryCnt < 3) {
            sprintf(err_msg, "Family %s has fewer than 3 members, skipping.",
                    PedTop->PedTree[ped].Name);
            warnf(err_msg);
            continue;
        }
        if (PedTop->PedTree[ped].EntryCnt > 1) {
            for(entry=0; entry < PedTop->PedTree[ped].EntryCnt; entry++) {
                if (PedTop->PedTree[ped].Entry[entry].Father != NULL) {
                    father =
                        PedTop->PedTree[ped].Entry[entry].Father->ID;
                    mother =
                        PedTop->PedTree[ped].Entry[entry].Mother->ID;
                    if (check_parents(father, mother, marriages)) {
                        marriages = append_to_marriage_list(marriages, father, mother);
                        num_marriages[ped]++;
                    }
                }
            }
            NewTop->PedCnt += num_marriages[ped];
        } else {
            if (PedTop->PedTree[ped].Entry[0].Sex == 1) {
                marriages=
                    append_to_marriage_list(marriages,
                                            PedTop->PedTree[ped].Entry[0].ID,
                                            0);
            } else {
                marriages=
                    append_to_marriage_list(marriages, 0,
                                            PedTop->PedTree[ped].Entry[0].ID);
            }
            (NewTop->PedCnt)++;
        }
        all_marriages = paste_marriage_lists(all_marriages, marriages);

    }

    /* Now allocate space for the new linkage pedigree structure */

    NewTop->Ped = CALLOC((size_t) NewTop->PedCnt, linkage_ped_tree);
    marriagep = all_marriages;
    for (pr=0; pr < NewTop->PedCnt; pr++) {
        if (marriagep->num_off > 0) {
            NewTop->Ped[pr].EntryCnt = marriagep->num_off + 2;
            if ((NewTop->Ped[pr].Entry =
                CALLOC((size_t) (marriagep->num_off + 2),linkage_ped_rec))
               == NULL) {
                errorf("Unable to allocate memory.\n");
                EXIT(MEMORY_ALLOC_ERROR);
            }
        } else {
            NewTop->Ped[pr].EntryCnt = 1;
            if ((NewTop->Ped[pr].Entry = CALLOC((size_t) 1,linkage_ped_rec)) == NULL) {
                errorf("Unable to allocate memory.\n");
                EXIT(MEMORY_ALLOC_ERROR);
            }
        }
        marriagep = marriagep->next;
    }

    marriagep = all_marriages;
    nuke_pednum = 0;

    NewTop->OrigIds = Top->OrigIds;
    NewTop->UniqueIds = Top->UniqueIds;
    for(ped = 0; ped < Top->PedCnt; ped++) {
        if (PedTop->PedTree[ped].EntryCnt < 3) {
            continue;
        }

        if (num_marriages[ped] == 0) {
            To = &(NewTop->Ped[nuke_pednum]);
            To->Num = nuke_pednum+1;
            To->Loops=NULL;
            strcpy(To->Name, Top->Ped[ped].Name);
            strcpy(To->PedPre, Top->Ped[ped].PedPre);
            To->Proband = 1;
            To->OriginalID = Top->Ped[ped].Num;
            To->origped = ped;
            To->IsTyped = 0;
            copy_le_static_data(&(Top->Ped[ped].Entry[0]), &(To->Entry[0]));
            To->Entry[0].Pheno  = Top->Ped[ped].Entry[0].Pheno;
            To->Entry[0].Marker = Top->Ped[ped].Entry[0].Marker;
            /* adjust some fields for mum */
            To->Entry[0].ID=1;
            To->Entry[0].Father=0;
            To->Entry[0].Mother=0;
            To->Entry[0].Next_PA_Sib=0;
            To->Entry[0].Next_MA_Sib=0;
            To->Entry[0].OrigProband=1;
            To->Entry[0].First_Offspring=0;
            To->Entry[0].ext_per_num = Top->Ped[ped].Entry[0].ID;
            To->Entry[0].ext_ped_num = ped+1;
            nuke_pednum++;
            marriagep = marriagep->next;
        } else {
            for(pr = 0; pr < num_marriages[ped]; pr++) {
                /* now copy the entries */
                To = &(NewTop->Ped[nuke_pednum]);
                /* copy static info, extracted from copy_lpedtree1 */
                To->Num = nuke_pednum + 1;
                To->Loops=NULL;
                strcpy(To->Name, Top->Ped[ped].Name);
                strcpy(To->PedPre, Top->Ped[ped].PedPre);
                To->Proband = 1;
                To->OriginalID = Top->Ped[ped].Num;
                To->origped = ped;
                To->IsTyped = 0;

                mother = marriagep->mother;
                father = marriagep->father;

                copy_le_static_data(PedTop->PedTree[ped].Entry[mother-1].LEntry,
                                    &(To->Entry[0]));
                To->Entry[0].Pheno  = (PedTop->PedTree[ped].Entry[mother-1].LEntry)->Pheno;
                To->Entry[0].Marker = (PedTop->PedTree[ped].Entry[mother-1].LEntry)->Marker;

                /* adjust some fields for mum */
                To->Entry[0].ID=1;
                strcpy(To->Entry[0].PerPre, (PedTop->PedTree[ped].Entry[mother-1].LEntry)->PerPre);
                To->Entry[0].Father=0;
                To->Entry[0].Mother=0;
                To->Entry[0].Next_PA_Sib=0;
                To->Entry[0].Next_MA_Sib=0;
                To->Entry[0].OrigProband=1;
                To->Entry[0].First_Offspring=3;
                To->Entry[0].ext_per_num = mother;
                To->Entry[0].ext_ped_num = ped+1;

                copy_le_static_data(PedTop->PedTree[ped].Entry[father-1].LEntry,
                                    &(To->Entry[1]));
                To->Entry[1].Pheno  = (PedTop->PedTree[ped].Entry[father-1].LEntry)->Pheno;
                To->Entry[1].Marker = (PedTop->PedTree[ped].Entry[father-1].LEntry)->Marker;


                /* adjust some fields for papa */
                To->Entry[1].ID=2;
                strcpy(To->Entry[1].PerPre, (PedTop->PedTree[ped].Entry[father-1].LEntry)->PerPre);
                To->Entry[1].Father=0;
                To->Entry[1].Mother=0;
                To->Entry[1].Next_PA_Sib=0;
                To->Entry[1].Next_MA_Sib=0;
                To->Entry[1].OrigProband=0;
                To->Entry[1].First_Offspring=3;
                To->Entry[1].ext_per_num = father;
                To->Entry[1].ext_ped_num = ped+1;
                newentry=2;

                for(entry=0; entry < PedTop->PedTree[ped].EntryCnt; entry++) {
                    if (PedTop->PedTree[ped].Entry[entry].Father != NULL) {
                        if ((father == PedTop->PedTree[ped].Entry[entry].Father->ID) &&
                           (mother == PedTop->PedTree[ped].Entry[entry].Mother->ID)) {
                            /* 	      copy_lpedrec(&(Top->Ped[ped].Entry[entry]),  */
                            /* 			   &(To->Entry[newentry]),  Top); */
                            copy_le_static_data(PedTop->PedTree[ped].Entry[entry].LEntry,
                                                &(To->Entry[newentry]));
                            To->Entry[newentry].Pheno = (PedTop->PedTree[ped].Entry[entry].LEntry)->Pheno;
                            To->Entry[newentry].Marker = (PedTop->PedTree[ped].Entry[entry].LEntry)->Marker;

                            /* adjust some fields for offspring */
                            To->Entry[newentry].ID=newentry+1;
                            sprintf(To->Entry[newentry].PerPre, "%d", newentry+1);
                            To->Entry[newentry].Father=2;
                            To->Entry[newentry].Mother=1;
                            To->Entry[newentry].Next_PA_Sib=
                                (((newentry+2) > To->EntryCnt)? 0 : newentry+2);
                            To->Entry[newentry].Next_MA_Sib=To->Entry[newentry].Next_PA_Sib;
                            To->Entry[newentry].OrigProband=0;
                            To->Entry[newentry].First_Offspring=0;
                            To->Entry[newentry].ext_per_num =
                                PedTop->PedTree[ped].Entry[entry].ID;
                            To->Entry[newentry].ext_ped_num = ped+1;

                            newentry++;
                        }
                    }
                }
                nuke_pednum++;
                marriagep = marriagep->next;
            }
        }
    }
    free_marriage_list(all_marriages);
    free(num_marriages);

    if (created) {
        free_all_including_ped_top(PedTop, NULL, NULL);
        PedTop=NULL;
    }
}

/* here the main function and subroutines of nuke.c (written by Sean Davis)
   has been modified to become the following option create_nuclear_families ...*/

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void  create_nuclear_families(linkage_ped_top **LPedTop,
			      char *outfl_name,
			      char *loutfl_name,
			      analysis_type *analysis,
			      file_format *infl_type,
			      int untyped_ped_opt,
			      linkage_ped_top **NukeTop)
{

    char           chr_str[3];
    int             single_output=0, sex_linked;
    int             numchr, i, linkage_format;
    linkage_ped_top *Top2;

    linkage_ped_tree *LPed;

    if (*infl_type == LINKAGE) {
        for (i = 0; i < Top->PedCnt; i++)  {
            LPed = &(Top->Ped[i]);
            if (LPed->Loops != NULL)   {
                sprintf(err_msg, "Connecting loops in pedigree %d...", LPed->Num);
                mssgf(err_msg);
                if (connect_loops(LPed, Top) < 0)    {
                    errorf("Fatal error! Aborting.");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    if (main_chromocnt > 1) {
        single_output=global_ou_files(1);
    }
//SL
    if (Top->LocusTop->SexLinked == 2 && single_output == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }
    nuclear_family_files(outfl_name, loutfl_name, single_output, &linkage_format);

    Top2=new_lpedtop();
    makenucs1(Top, NULL, Top2);

    if (UntypedPeds != NULL) {
        free(UntypedPeds);
    }

    UntypedPeds = CALLOC((size_t) Top2->PedCnt, int);

    if (main_chromocnt == 1 || single_output==1) {
        if (main_chromocnt == 1 && global_chromo_entries[0] == UNKNOWN_CHROMO) {
            get_unmapped_loci(2);
        } else {
            /* main_chromocnt > 1 or known chromo */
            get_loci_on_chromosome(0);
            if (main_chromocnt > 0) {
                for (i=0; i < main_chromocnt; i++) {
                    if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
                        get_unmapped_loci(1);
                        break;
                    }
                }
            }
        }

        if (main_chromocnt > 1) {
            change_output_chr(outfl_name, 0);
            change_output_chr(loutfl_name, 0);
        } else {
            change_output_chr(outfl_name, global_chromo_entries[0]);
            change_output_chr(loutfl_name, global_chromo_entries[0]);
        }

        if (linkage_format) {
            Top2->pedfile_type = POSTMAKEPED_PFT;
        } else {
            Top2->pedfile_type = PREMAKEPED_PFT;
        }

        if (main_chromocnt == 1 && global_chromo_entries[0] == SEX_CHROMOSOME) {
            sex_linked = 1;
        } else {
            sex_linked = ((Top2->LocusTop->SexLinked == 1)? 1 : 0);
        }
        write_gh_locus_file(loutfl_name, Top2->LocusTop, *analysis, sex_linked);

        if (OrigIds[1] == 3) {
            /* This is consecutive numbering option,
               set pedigree names to <origped_nucpednum> as is done in
               create_nuked_fids below */
            set_nuked_pednames(Top2);
        }

        omit_peds(untyped_ped_opt, Top2);
        write_nuclear_peds(outfl_name, Top2, 0, linkage_format);
        draw_line();
        create_mssg(TO_NUKE);
        sprintf(err_msg, "        Locus file: %s", loutfl_name);
        mssgf(err_msg);
        sprintf(err_msg, "        Pedigree file: %s", outfl_name);
        mssgf(err_msg);
        draw_line();
        *NukeTop = Top2;
        /*    free_all_including_lpedtop(&Top2); */
        return;
    }

/*   Top_copy=new_lpedtop(); */
/*   copy_linkage_ped_top(Top, Top_copy, 1); */

    draw_line();
    create_mssg(TO_NUKE);
    for (i=0; i< main_chromocnt; i++) {
        numchr=global_chromo_entries[i];
        if (numchr == UNKNOWN_CHROMO) {
            get_unmapped_loci(2);
        } else {
            get_loci_on_chromosome(numchr);
        }
        sex_linked =
            (((numchr == SEX_CHROMOSOME && Top->LocusTop->SexLinked == 2) ||
              (Top->LocusTop->SexLinked == 1))? 1 : 0);


        CHR_STR(numchr, chr_str);

        /* create the chromosome specific lpedtop */
/*     free_all_including_lpedtop(&Top); Top=new_lpedtop(); */
/*     copy_linkage_ped_top(Top_copy, Top, 1); */
/*     ReOrderMappedLoci(Top, &numchr); */

        omit_peds(untyped_ped_opt, Top2);
        if (linkage_format) {
            Top2->pedfile_type = pedfile_type = POSTMAKEPED_PFT;
        } else {
            Top2->pedfile_type = pedfile_type = PREMAKEPED_PFT;
        }
        change_output_chr(outfl_name, numchr);
        change_output_chr(loutfl_name, numchr);
        write_gh_locus_file(loutfl_name, Top->LocusTop, *analysis, sex_linked);
        write_nuclear_peds(outfl_name, Top2, 0, linkage_format);

        sprintf(err_msg, "        Locus file: %s", loutfl_name);
        mssgf(err_msg);
        sprintf(err_msg, "        Pedigree file: %s", outfl_name);
        mssgf(err_msg);
        /* if not the last iteration */

    }
    draw_line();
    *NukeTop = Top2;
}

#undef Top
/* Sets the names of the pedigrees */

static void     set_nuked_pednames(linkage_ped_top *Topp)
{

    int j, i, *count = CALLOC((size_t) Topp->PedCnt, int);
    int ext_ped_num = -1, next_ped=0;
    char tmpName[MAX_NAMELEN];


    /* First count how many nuclear pedigrees each pedigree has */

    for (i = 0; i < Topp->PedCnt; i++) {
        if (i == 0){
            ext_ped_num=0;
            count[ext_ped_num] = 1;
        } else if (Topp->Ped[i].OriginalID == Topp->Ped[i - 1].OriginalID) {
            count[ext_ped_num]++;
        } else {
            ext_ped_num++;
            count[ext_ped_num] = 1;
        }
    }

    for(j=0; j <= ext_ped_num; j++) {
        if (count[j] > 1) {
            if (Topp->OrigIds) {
                for(i=0; i < count[j]; i++) {
                    sprintf(Topp->Ped[i+next_ped].Name, "%d-%d",
                            Topp->Ped[i+next_ped].OriginalID, i+1);
                }
            } else {
                for(i=0; i < count[j]; i++) {
                    sprintf(tmpName, "%s-%d",
                            Topp->Ped[i+next_ped].Name, i+1);
                    sprintf(Topp->Ped[i+next_ped].Name, "%s", tmpName);
                }
            }
        }
        next_ped += count[j];
    }

    free(count);
}
/*------------- following creates new fids for the nuked families -----------*/

void            create_nuked_fids(linkage_ped_top *Topp)
{

    int j, i, *count = CALLOC((size_t) Topp->PedCnt, int);
    int *names = CALLOC((size_t) Topp->PedCnt, int);
    int ext_ped_num = -1, next_ped=0;
    int clash=0, multiplier=NukedMultiplier;
    char tmpName[MAX_NAMELEN];

    /* First count how many nuclear pedigrees each pedigree has */
    for (i = 0; i < Topp->PedCnt; i++) {
        if (i == 0){
            ext_ped_num=0;
            count[ext_ped_num] = 1;
        } else if (Topp->Ped[i].OriginalID == Topp->Ped[i - 1].OriginalID) {
            count[ext_ped_num]++;
        } else {
            ext_ped_num++;
            count[ext_ped_num] = 1;
        }
    }

    for(j=0; j <= ext_ped_num; j++) {
        if (count[j] > 1) {
            if (Topp->OrigIds) {
                for(i=0; i < count[j]; i++) {
                    sprintf(Topp->Ped[i+next_ped].Name, "%d_%d",
                            Topp->Ped[i+next_ped].OriginalID, i+1);
                }
            } else {
                for(i=0; i < count[j]; i++) {
                    sprintf(tmpName, "%s_%d",
                            Topp->Ped[i+next_ped].Name, i+1);
                    sprintf(Topp->Ped[i+next_ped].Name,"%s",tmpName);
                }
            }
        }
        next_ped += count[j];
    }

    while(1) {
        next_ped = 0;
        for(j=0; j <= ext_ped_num; j++) {
            if (count[j] > 1) {
                for(i=0; i < count[j]; i++) {
                    if (Topp->OrigIds) {
                        names[i+next_ped]=
                            Topp->Ped[i+next_ped].OriginalID*multiplier + i + 1;
                    } else {
                        names[i+next_ped]=
                            atoi(Topp->Ped[i+next_ped].Name)*multiplier + i + 1;
                    }
                }
            } else {
                names[next_ped] = ((Topp->OrigIds)?
                                   Topp->Ped[next_ped].OriginalID :
                                   atoi(Topp->Ped[next_ped].Name));

            }
            next_ped += count[j];
        }

        /* Now make sure that the extended names do not clash with any
           others */
        clash=0;
        for (i=0; i < (Topp->PedCnt-1); i++) {
            for (j=(i+1); j < Topp->PedCnt; j++) {
                if (names[i] == names[j]) {
                    clash=1;
                    break;
                }
            }
            if (clash == 1) {
                printf("In creating ped ids for nuclear peds\n");
                printf("Trying to generate unique ped ids\n");
                break;
            }
        }
        if (clash) {
            /* Multiply multipler by 10 */
            multiplier *= 10;
        } else {
            break;
        }
    }

    for (i=0; i < Topp->PedCnt; i++) {
        Topp->Ped[i].Num = names[i];
    }
    free(names);
    free(count);
}
