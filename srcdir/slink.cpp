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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "pedtree_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_ghfiles_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file haldane_theta kosambi_theta
            linkage_ext.h:  get_loci_on_chromosome is_typed_lentry switch_penetrances
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
            pedtree_ext.h:  convert_to_pedtree reassign_affecteds
         user_input_ext.h:  affected_by_status test_modified
              utils_ext.h:  EXIT dmax draw_line script_time_stamp
      write_ghfiles_ext.h:  write_gh_locus_file
*/


#define DEFAULT_sim_pheno 0
#define DEFAULT_min_affecteds 0
#define DEFAULT_all_or_any  0
#define DEFAULT_before_or_after 0


#define  DEFAULT_seed 21341
#define  DEFAULT_replicates 100
#define  DEFAULT_proportion  0.00

/* all routines related to slink_format files */

/*--------------static------------------*/
static void SLINK_file_names(char *outfl_name, char *loutfl_name,
			     char *slinkin_name,
			     char *slink_shell_name);
static void SLINK_simgeno_opt(int *sim_pheno, int *min_affecteds,
			      int *all_or_any, int *before_or_after,
			      int *seed, int *replicates,
			      double *proportion,
			      int has_aff);
static int find_entry_code(linkage_ped_rec *Entry, linkage_locus_top *LTop,
			   int all_any,
			   int before_after, int sim_geno);

static void write_slink_pedigree(ped_top *PTop, linkage_ped_top *Top,
				 char *fl_name,  int disease_loc,
				 int sim_pheno, int all_or_any,
				 int before_after, int min_affecteds,
				 int show);

static void write_slink_shell(char *shell_name, char *loutfl_name,
			      char *outfl_name, char *slinkin_name);

static int select_disease_locus(linkage_locus_top *LTopp);

/*------------exported-------------------*/

void create_SLINK_format_files(linkage_ped_top *LPedTreeTop, ped_top *PedTreeTop,
			       analysis_type * analysis,
			       int *numchr,
			       char *fle_names[], int untyped_ped_opt);

int write_slinkin_file(char *fl_name,
		       int disease_locus_index,
		       int seed, int replicates, double proportion);

static int slink_thetas(linkage_locus_top *LTopp,
			int disease_locus_index);

/*--------------end of prototypes---------------*/


static void SLINK_simgeno_opt(int *sim_pheno, int *min_affecteds,
			      int *all_or_any, int *before_or_after,
			      int *seed, int *replicates,
			      double *proportion,
			      int has_aff)

{

    char copt[10], all_any[3]={"* "}, before_after[3] = {"* "};
    int item=0, i, opt1 = -1, opt2;
    char all_any_str[2][15] = {"Any marker",
                               "All markers"};
    char before_after_str[2][20] = {"Before reordering",
                                    "After reordering"};
    char min_aff_messg[2][100] =
        {"Include families with or without affecteds",
         "Include families with minimum affected count: "};

    double prop;

    *sim_pheno= DEFAULT_sim_pheno;
    *min_affecteds = DEFAULT_min_affecteds;
    *all_or_any = DEFAULT_all_or_any;
    *before_or_after = DEFAULT_before_or_after;
    *seed = DEFAULT_seed;
    *replicates = DEFAULT_replicates;
    *proportion = DEFAULT_proportion;


    if (DEFAULT_OPTIONS) {
        mssgf("SLINK pedigree and phenotype options set to defaults");
        mssgf("as requested in the batch file.");
        return;
    }

    while (opt1 != 0) {
        item=1;
        draw_line();
        printf("SLINK pedigree and phenotype options:\n");
        printf("0) Done with this menu - please proceed\n");
        if (has_aff) {
            printf(" %d) %s", item,
                   ((*min_affecteds)? min_aff_messg[1] : min_aff_messg[0]));
            if (*min_affecteds) { printf("%d\n", *min_affecteds); }
            else { printf("\n"); }
            item++;
        } else {
            *min_affecteds = 0;
        }

        printf(" %d) Set trait phenotypes to unknown: %s\n", item,
               yorn[*sim_pheno]);
        item++;
        printf(" %d) Individual available if genotyped at %s %s\n", item,
               all_any_str[*all_or_any], before_after_str[*before_or_after]);
        item++;
        printf(" %d) Random seed: %d\n", item, *seed);
        item++;
        printf(" %d) Number of replicates: %d\n", item, *replicates);
        item++;
        printf(" %d) Proportion of UNLINKED families: %f\n", item, *proportion);

        printf("Select from options 0 - %d > ", item);
        fcmap(stdin, "%s", copt); newline;
        sscanf(copt, "%d", &opt1);
        if (!has_aff && opt1) {
            opt1++;
        }
        switch(opt1) {
        case 0:
            break;
        case 1:
            printf("Enter minimum number of affecteds > ");
            fcmap(stdin, "%s", copt); newline;
            i=sscanf(copt, "%d", min_affecteds);
            if (i != 1 || *min_affecteds < 0) {
                printf("Value entered is not 0 or a positive integer.\n");
                *min_affecteds=0;
            }
            break;
        case 2:
            *sim_pheno = ((*sim_pheno) ? 0 : 1);
            break;
        case 3:
            opt2 = -1;
            while(opt2 != 0) {
                printf("Individual available if genotyped at: \n");
                printf("0) Done with this menu - please proceed\n");
                printf(" 1) %c %s/ %c %s\n", all_any[0], all_any_str[0],
                       all_any[1], all_any_str[1]);
                printf(" 2) %c %s/ %c %s\n",
                       before_after[0], before_after_str[0],
                       before_after[1], before_after_str[1]);
                printf("Enter 1 or 2 to toggle selection, 0 to exit > ");
                fcmap(stdin, "%s", copt);     newline;
                sscanf(copt, "%d", &opt2);
                switch(opt2) {
                case 0:
                    break;
                case 1:
                    if (*all_or_any == 0) {
                        *all_or_any=1;	strcpy(all_any," *");
                    } else {
                        *all_or_any = 0;	strcpy(all_any, "* ");
                    }
                    break;
                case 2:
                    if (*before_or_after == 0) {
                        *before_or_after=1;	strcpy(before_after, " *");
                    } else {
                        *before_or_after = 0;	strcpy(before_after, "* ");
                    }
                    break;
                default:
                    warn_unknown(copt);
                    break;
                }
            }
            break;
        case 4:
            printf("Enter value of random seed ( 25000 -> 30000 ) > ");
            fcmap(stdin, "%s", copt); newline;
            sscanf(copt, "%d", &opt2);
            if ((opt2 < 25000) || (opt2 > 30000)) {
                warn_unknown(copt);
            } else {
                *seed = opt2;
            }
            break;
        case 5:
            printf("Enter desired number of replicates (minimum 1) > ");
            fcmap(stdin, "%s", copt); newline;
            sscanf(copt, "%d", &opt2);
            if (opt2 < 1) {
                warn_unknown(copt);
            } else {
                *replicates = opt2;
            }
            break;
        case 6:
            printf("Enter proportion of UNLINKED families between 0 and 1\n");
            printf("e.g. enter '0.1' to have 90%% of your families be linked > ");
            fcmap(stdin, "%s", copt); newline;
            sscanf(copt, "%lf", &prop);
            if (prop < 0.0 || prop > 1.0) {
                warn_unknown(copt);
            } else {
                *proportion = prop;
            }
            break;
        default:
            warn_unknown(copt);
            break;
        }
    }

    return;

}

static void SLINK_file_names(char *outfl_name, char *loutfl_name,
			     char *slinkin_name,
			     char *slink_shell_name)

{
    char reply[10], fl_stat[12];
    int user_reply = -1;
    char stem[5];

    if (DEFAULT_OUTFILES) {
        return;
    }

    if (main_chromocnt > 1) {
        change_output_chr(slinkin_name, 0);
        change_output_chr(slink_shell_name, 0);
        change_output_chr(outfl_name, -9);
        change_output_chr(loutfl_name, -9);
        strcpy(stem, "stem");
    } else {
        strcpy(stem, "");
    }

    while (user_reply != 0) {
        draw_line();
        print_outfile_mssg();
        printf("SLINK output file name menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt > 1) {
            printf(" 1) Pedigree file name stem:   %-15s\n", outfl_name);
            printf(" 2) Locus file name stem:      %-15s\n", loutfl_name);
        } else {
            printf(" 1) Pedigree file name:        %-15s\t%s\n", outfl_name,
                   file_status(outfl_name, fl_stat));
            printf(" 2) Locus file name:           %-15s\t%s\n", loutfl_name,
                   file_status(loutfl_name, fl_stat));
        }
        printf(" 3) Parameter file name:       %-15s\t%s\n", slinkin_name,
               file_status(slinkin_name, fl_stat));
        printf(" 4) C-shell file  name:        %-15s\t%s\n", slink_shell_name,
               file_status(slink_shell_name, fl_stat));

        printf("Enter options 1-4 to change file names > ");
        fcmap(stdin, "%s", reply); newline;
        user_reply=-1;
        sscanf(reply, "%d", &user_reply);
        test_modified(user_reply);

        switch (user_reply)  {
        case 1:
            printf("Enter NEW SLINK pedigree file name %s: > ", stem);
            fcmap(stdin, "%s", outfl_name); newline;
            break;
        case 2:
            printf("Enter NEW SLINK locus file name %s: > ", stem);
            fcmap(stdin, "%s", loutfl_name); newline;
            break;
        case 3:
            printf("Enter NEW SLINK parameter file name: > ");
            fcmap(stdin, "%s", slinkin_name); newline;
            break;
        case 4:
            printf("Enter NEW SLINK C-shell file name: > ");
            fcmap(stdin, "%s", slink_shell_name); newline;
            break;
        case 0:
            break;
        default:
            warn_unknown(reply);
            break;
        }
    }

    if (main_chromocnt > 1) {
        change_output_chr(outfl_name, global_chromo_entries[0]);
        change_output_chr(loutfl_name, global_chromo_entries[0]);
    }
    return;

}

/*------------------------------------------------------------*/
static int find_entry_code(linkage_ped_rec *Entry, linkage_locus_top *LTop,
			   int all_any,
			   int before_after, int sim_pheno)

{
    /*                    Marker=0           Marker=1
                          Simchoice =  0     00                 10
                          Simchoice =  1     11                 01
    */

    int ret = -1, marker_avail;

    if (before_after == 0) {
        marker_avail = ((all_any == 0 && Entry->IsTyped > 0) ||
                        (all_any == 1 && Entry->IsTyped==2));
    } else {
        marker_avail = is_typed_lentry(*Entry, LTop, !all_any, NULL);
    }

    if (marker_avail && sim_pheno)   ret = 1;
    if (marker_avail && !sim_pheno)  ret = 2;
    if (!marker_avail && sim_pheno)  ret = 3;
    if (!marker_avail && !sim_pheno) ret = 0;

    return ret;

}
/* NOTE - (Nandita, 5/25/01)
   We were assuming until now that only affection status loci
   can be used as the disease locus, but this is not true. In
   case of QTLs, the user shouldn't be asked who is affected


   k=0;
   for (i=0; i < num_traits; i++) {

   k++;
   }
   }
   aff = CALLOC((size_t) k, int);
   k=0;
   for (i=0; i < num_traits; i++) {
   if (LTopp->Locus[global_trait_entries[i]].Type == AFFECTION) {
   aff[k]=global_trait_entries[i];
   k++;
   }
   }
   disease_locus_index = aff[0];
*/
/* slink options */

static int slink_thetas(linkage_locus_top *LTopp,
			int disease_locus_index)
{
    int i, theta_ct = 0;
    double *theta;
    int num_markers=0, *markers;
    
    if (NumChrLoci > 0) {
        markers = CALLOC((size_t) NumChrLoci, int);
    } else {
        // Yes, and so we just keep exeuting the code below till we segfault?!
        markers = NULL;
    }
    
    /* set the disease_locus number with the index from ChrLoci */
    for(i=0; i < NumChrLoci; i++) {
        if (LoopOverTrait == 0 ||
            LTopp->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTopp->Locus[ChrLoci[i]].Type == BINARY) {
            markers[num_markers]=ChrLoci[i];
            num_markers++;
        }
    }
    
    /* Switch the disease locus to the index if there is one */
    if (LoopOverTrait == 0 && disease_locus_index >= 0) {
        for(i=0; i < NumChrLoci; i++) {
            if (ChrLoci[i] == global_trait_entries[disease_locus_index]) {
                disease_locus_index = i;
                break;
            }
        }
    }
    
    theta = CALLOC((size_t) imax(num_markers, 1), double);
    theta_ct=0;
    if (LoopOverTrait == 1) {
        /* trait is first */
        theta[0]=0.5; theta_ct++;
        for (i = 1; i < num_markers; i++)  {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i]].pos_avg -
                                LTopp->Marker[markers[i - 1]].pos_avg);
            theta[i] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
    } else if (LoopOverTrait == 0 && disease_locus_index < 0) {
        /* How do you get a disease_locus_index < 0 ? */
        for (i = 1; i < num_markers; i++)  {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i]].pos_avg -
                                LTopp->Marker[markers[i - 1]].pos_avg);
            theta[i-1] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
    } else if (LoopOverTrait == 0 && disease_locus_index == 0) {
        theta[0]=0.5; theta_ct++;
        for (i = 2; i < num_markers; i++)  {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i]].pos_avg -
                                LTopp->Marker[markers[i - 1]].pos_avg);
            theta[i-1] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
    } else if (LoopOverTrait == 0 && disease_locus_index == num_markers-1) {
        for (i = 0; i < num_markers-2; i++)  {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i + 1]].pos_avg -
                                LTopp->Marker[markers[i]].pos_avg);
            theta[i] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
        theta[num_markers-2]=0.5;
        theta_ct++;
    } else  {
        for (i=1; i < disease_locus_index; i++) {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i]].pos_avg -
                                LTopp->Marker[markers[i - 1]].pos_avg);
            theta[i-1] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
        theta[disease_locus_index-1]=0.5;
        theta[disease_locus_index] =
        dmax(0.0,
             (LTopp->Marker[markers[disease_locus_index+1]].pos_avg -
              LTopp->Marker[markers[disease_locus_index-1]].pos_avg));
        theta[disease_locus_index] = ((LTopp->map_distance_type == 'h')?
                                      haldane_theta(theta[disease_locus_index]):
                                      kosambi_theta(theta[disease_locus_index]));
        
        theta_ct += 2;
        for (i=disease_locus_index+1; i < num_markers-1; i++) {
            double delta = dmax(0.0,
                                LTopp->Marker[markers[i + 1]].pos_avg -
                                LTopp->Marker[markers[i]].pos_avg);
            theta[i] = (LTopp->map_distance_type == 'h') ? haldane_theta(delta) : kosambi_theta(delta);
            theta_ct++;
        }
    }
    
    /* changed this part- Nandita, sets the MaleRecomb
     instead of writing into the slinkin file, storing the number
     inside the first element */
    if (LTopp->MaleRecomb != NULL) {
        free(LTopp->MaleRecomb);
    }
    LTopp->MaleRecomb = (double *) CALLOC((size_t) theta_ct + 1, double);
    LTopp->MaleRecomb[0]=(double) theta_ct;
    for (i = 1; i <= theta_ct; i++) {
        LTopp->MaleRecomb[i] = theta[i-1];
    }
    free(theta); free(markers);
    
    return disease_locus_index;
}

/*------------------------------------------------------------*/
int write_slinkin_file(char *slinkin_name,
		       int disease_locus_index,
		       int seed,
		       int replicates,
		       double proportion)
{
    int            tr, nloop, num_affec = num_traits;
    char           slinkin[2*FILENAME_LENGTH];
    FILE           *fp;

    /* now write the slinkin file */
    NLOOP;

    for(tr=0; tr<=nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        sprintf(slinkin, "%s/%s", output_paths[tr], slinkin_name);
        if ((fp = fopen(slinkin, "w")) == NULL) {
            errorvf("Could not open '%s' for writing.\n", slinkin);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(fp, "%d\n%d\n%d\n%f\n", seed, replicates,
                disease_locus_index+1, proportion);
        fclose(fp);
        if (nloop == 1) break;
    }
    return disease_locus_index;
}


/*------------------------------------------------------------*/
static void ped_count_message(int ped_cnt, linkage_locus_top *LocusTop,
			      int disease_loc)

{

    if (disease_loc >= 0) {
        sprintf(err_msg, "Saving %d pedigrees for trait %s.", ped_cnt,
                LocusTop->Locus[disease_loc].LocusName);
    } else {
        sprintf(err_msg, "Saving %d pedigrees (no trait).", ped_cnt);
    }

    mssgf(err_msg);
    return;
}

static void   write_slink_pedigree(ped_top *PTop, linkage_ped_top *Top,
				   char *pedfl_name, int disease_loc,
				   int sim_pheno, int all_or_any,
				   int before_or_after, int min_affecteds,
				   int show_numsaved)
{

    int    nloop, tr, *tr_pos, num_affec=num_traits;
    int    ped, entry_code, i, k, is_aff;
    int    *loc_order, num_loc, ped_cnt;
    linkage_ped_rec *tpe;
    char ofl[2*FILENAME_LENGTH];
    FILE *fp;
/* analysis_type   analysis=TO_SLINK; */

    NLOOP;

//  if (num_affec > 0)
//      tr_pos = &global_trait_entries[0];
    tr_pos = global_trait_entries;

    loc_order = CALLOC((size_t) NumChrLoci+1, int);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;

        sprintf(ofl, "%s/%s", output_paths[tr], pedfl_name);
        if ((fp = fopen(ofl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", ofl);
            EXIT(FILE_WRITE_ERROR);
        }
        num_loc=0;
        if (LoopOverTrait == 1) {
            /* copy the disease locus first */
            loc_order[0] = *tr_pos;
            num_loc++;
        }

        for (i=0; i < NumChrLoci; i++) {
            if (LoopOverTrait == 0 ||
                Top->LocusTop->Locus[ChrLoci[i]].Type == NUMBERED ||
                Top->LocusTop->Locus[ChrLoci[i]].Type == BINARY) {
                loc_order[num_loc] = ChrLoci[i]; num_loc++;
            }
        }

        /* take out this check later */
/*     if (k != num_loc) {  */
/*       printf("Locus numbers do not match!\n"); */
/*       exit(-1); */
/*     } */

        /* first find out the affected status wrt current trait */
        if (LoopOverTrait == 1) {
            if (Top->LocusTop->Locus[*tr_pos].Type == AFFECTION) {
                is_aff=1;
                affected_by_status(Top, *tr_pos);
            } else {
                is_aff=0;
            }
            disease_loc=*tr_pos;
        } else {
            if (num_traits > 1) {
                if (Top->LocusTop->Locus[global_trait_entries[disease_loc]].Type == AFFECTION) {
                    is_aff=1;
                    affected_by_status(Top, global_trait_entries[disease_loc]);
                } else {
                    is_aff=0;
                }
            } else {
                is_aff=0;
            }
        }

        ped_cnt=0;
        for (ped=0; ped < PTop->PedCnt; ped++) {
            if (is_aff) {
                /* assign affecteds in pedtree */
                reassign_affecteds(PTop);

                /* find out if this pedigree has minimum number of affecteds */
                if (PTop->PedTree[ped].AffectedCnt < min_affecteds) {
                    sprintf(err_msg,
                            "Not saving pedigree %s (%d affected, less than %d).\n",
                            PTop->PedTree[ped].Name,
                            PTop->PedTree[ped].AffectedCnt,
                            min_affecteds);
                    mssgf(err_msg);
                    continue;
                    /* don't output this pedigree */
                }
            }
            /* write pedigree if disease locus is a QTL or there are at least
               <min_affecteds> affecteds */
            ped_cnt++;
            for (i = 0; i < Top->Ped[ped].EntryCnt; i++) {
                tpe = &(Top->Ped[ped].Entry[i]);
                fprintf(fp, "%-4d %-4d %-4d %-4d %-4d %-4d %-4d %1d %1d ",
                        Top->Ped[ped].Num,
                        tpe->ID,
                        tpe->Father,
                        tpe->Mother,
                        tpe->First_Offspring,
                        tpe->Next_PA_Sib,
                        tpe->Next_MA_Sib,
                        tpe->Sex,
                        tpe->OrigProband);

                entry_code=find_entry_code(tpe,
                                           Top->LocusTop,
                                           all_or_any, before_or_after,
                                           sim_pheno);
                /* now use the loc order to write the loci correctly */
                /* NM- If phenotype is to be simulated as well, write 0 for the status,
                   and output liability class as is, otherwise, output status as well */
                for (k = 0; k < num_loc; k++) {
                    switch(Top->LocusTop->Locus[loc_order[k]].Type) {
                    case AFFECTION:
                        if (sim_pheno && loc_order[k] == disease_loc) {
                            fprintf(fp, "  0");
                        } else {
                            fprintf(fp, " %2d", tpe->Pheno[loc_order[k]].Affection.Status);
                        }
                        if (Top->LocusTop->Pheno[loc_order[k]].Props.Affection.ClassCnt > 1) {
                            fprintf(fp, " %-2d",
                                    tpe->Pheno[loc_order[k]].Affection.Class);
                        }
                        break;
                    case QUANT:
                        if (fabs(tpe->Pheno[loc_order[k]].Quant - MissingQuant) <= EPSILON ||
                            (sim_pheno && loc_order[k] == disease_loc)) {
                            fprintf(fp, "    0.0   ");
                        } else {
                            fprintf(fp, "%10.5f", tpe->Pheno[loc_order[k]].Quant);
                        }
                        break;
                    case NUMBERED:
                    case BINARY:
                        fprintf(fp, "  0  0");
                        break;
                    default:
                        break;
                    }
                }

                fprintf(fp, "%2d\n", entry_code);
                /* here is the 'code' column ' */
            }
        }
        fclose(fp);
        if (show_numsaved) {
            if (LoopOverTrait == 1) {
                ped_count_message(ped_cnt, Top->LocusTop, *tr_pos);
            } else {
                ((disease_loc >= 0)?
                 ped_count_message(ped_cnt, Top->LocusTop, global_trait_entries[disease_loc]) :
                 ped_count_message(ped_cnt, Top->LocusTop, disease_loc));
            }
        }
        if (nloop == 1) break;
        tr_pos++;
    }
    free(loc_order);
    return;

}


void create_SLINK_format_files(linkage_ped_top *LPedTreeTop,
			       ped_top * PedTreeTop,
			       analysis_type *analysis,  int *numchr,
			       char *file_names[], int untyped_ped_opt)
{

    char            *outfl_name, *loutfl_name, *slinkin_name, *slink_shell_name;
    int             i, disease_locus, has_aff;
    int             sim_pheno, min_affecteds, all_or_any, before_or_after;
    int             sex_linked;
    int             seed, replicates;
    double          proportion;
    int             disease_locus_index = 0;

    linkage_locus_top *LTopp = LPedTreeTop->LocusTop;

    outfl_name=file_names[0];
    loutfl_name=file_names[1];
    slinkin_name = file_names[2];
    slink_shell_name=file_names[3];

    SLINK_file_names(outfl_name, loutfl_name, slinkin_name,
                     slink_shell_name);
    /*set_missing_quant_input(LPedTreeTop); */

    disease_locus = select_disease_locus(LTopp);

    has_aff = 0;

    if (ManualReorder) {
        if (disease_locus >= 0) {
            /* find the disease locus */
            for (i = 0; i < NumChrLoci; i++) {
                if (ChrLoci[i] == global_trait_entries[disease_locus]) {
                    disease_locus_index = i;
                    if (LTopp->Locus[ChrLoci[i]].Type == AFFECTION) {
                        has_aff = 1;
                    }
                    break;
                }
            }
        }
    } else {
        if (LoopOverTrait == 0) {
            if (disease_locus >= 0) {
                /* disease_locus = global_trait_entries[disease_locus]; */
                if (LTopp->Locus[global_trait_entries[disease_locus]].Type == AFFECTION) {
                    has_aff = 1;
                }
            }
        } else {
            has_aff = 0;
            for (i = 0; i < num_traits; i++) {
                if (LTopp->Locus[global_trait_entries[i]].Type == AFFECTION) {
                    has_aff = 1;
                    break;
                }
            }
        }
    }

    /* analysis options */
    SLINK_simgeno_opt(&sim_pheno, &min_affecteds, &all_or_any,
                      &before_or_after, &seed, &replicates, &proportion, has_aff);

    /* determine the number of affecteds */
    /* convert to pedtree, but don't assign affecteds */

    PedTreeTop = convert_to_pedtree(LPedTreeTop, 0);

    for (i = 0; i < main_chromocnt; i++) {
        *numchr = global_chromo_entries[i];
        change_output_chr(outfl_name, *numchr);
        change_output_chr(loutfl_name, *numchr);
        change_output_chr(slinkin_name, *numchr);
        get_loci_on_chromosome(*numchr);
        omit_peds(untyped_ped_opt, LPedTreeTop);
	// SINK should (optionally) support sex-specific maps by outputting a simdata.* locus
	// file with sex-specific recombination fractions (this locus file follows the LINKAGE
	// formatting rules). 
        disease_locus_index = slink_thetas(LPedTreeTop->LocusTop, disease_locus);
        write_slinkin_file(slinkin_name, disease_locus_index, seed, replicates,
                           proportion);

        sex_linked = (((*numchr == SEX_CHROMOSOME && LTopp->SexLinked == 2)
                       || (LTopp->SexLinked == 1)) ? 1 : 0);

        write_slink_pedigree(PedTreeTop, LPedTreeTop, outfl_name,
                             disease_locus, sim_pheno, all_or_any, before_or_after,
                             min_affecteds, ((i == 0) ? 1 : 0));

        write_gh_locus_file(loutfl_name, LPedTreeTop->LocusTop, *analysis, sex_linked);

        if (i == 0) create_mssg(*analysis);

        sprintf(err_msg, "   Pedigree file:      %s", outfl_name);
        mssgf(err_msg);
        sprintf(err_msg, "   Locus file:         %s", loutfl_name);
        mssgf(err_msg);

    }

    sprintf(err_msg, "   Parameter file:     %s", slinkin_name);
    mssgf(err_msg);

    write_slink_shell(slink_shell_name, loutfl_name, outfl_name, slinkin_name);
    sprintf(err_msg, "   C-shell file:       %s", slink_shell_name);
    mssgf(err_msg);
    draw_line();
} /* End of create_SLINK_format_files  */


static void write_slink_shell(char *shell_name, char *loutfl_name,
			      char *outfl_name, char *slinkin_name)

{
    char shfl[2*FILENAME_LENGTH], syscmd[100];
    FILE *slink_shell_p;
    int nloop, tr, num_affec=num_traits;
    int i;

    NLOOP;
    for (tr=0; tr <=nloop; tr++) {
        if (tr==0 && nloop > 1) continue;

        sprintf(shfl, "%s/%s", output_paths[tr], shell_name);
        delete_file(shfl);
        if ((slink_shell_p = fopen(shfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", shfl);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(slink_shell_p, "#!/bin/csh -f\n");
        script_time_stamp(slink_shell_p);
        fprintf(slink_shell_p, "unalias rm\n");
        fprintf(slink_shell_p, "rm -f simdata.dat\n");
        fprintf(slink_shell_p, "rm -f simped.dat\n");
        fprintf(slink_shell_p, "rm -f slinkin.dat\n");
        fprintf(slink_shell_p, "cp %s slinkin.dat\n", slinkin_name);
        for (i=0; i < main_chromocnt; i++) {
            change_output_chr(loutfl_name, global_chromo_entries[i]);
            change_output_chr(outfl_name, global_chromo_entries[i]);
            fprintf(slink_shell_p, "cp %s simdata.dat\n", loutfl_name);
            fprintf(slink_shell_p, "cp %s simped.dat\n", outfl_name);
            fprintf(slink_shell_p, "Echo Running slink on chromosome %d\n",
                    global_chromo_entries[i]);
            fprintf(slink_shell_p, "slink\n");
        }

        fclose(slink_shell_p);
        sprintf(syscmd, "chmod +x %s", shfl);
        System(syscmd);
        if (nloop == 1) break;
    }
    return;

}


static int select_disease_locus(linkage_locus_top *LTopp)

{

    int i, j, m_ind = INT_MAX, disease_locus_index = -1, choice_ = -1;
    char choice[10];

    if (LoopOverTrait == 1) {
        disease_locus_index = 0;
    } else if (num_traits == 2) {
        /* Assumes one 'trait' is the markers with global_trait_entries = -1, and
         * the other is the 'trait'.  This is NOT true if only two traits are
         * selected, and no markers are selected.
         */
        disease_locus_index = ((global_trait_entries[0] < 0) ? 1 : 0);
    } else if (num_traits == 1) {
        disease_locus_index = -1;
    } else {
        for (i = 0; i < num_traits; i++) {
            if (global_trait_entries[i] < 0) continue;
            else if (LTopp->Locus[global_trait_entries[i]].Class == TRAIT) {
                disease_locus_index = i;
                break;
            }
        }
        if (DEFAULT_OPTIONS) {
            mssgf("Selected first affection locus in loci list (DEFAULT)");
            mssgf("as requested in the batch file.");
            return (disease_locus_index);
        }
        /* else */
        while (1) {
            draw_line();
            printf("Disease locus selection menu\n");
            printf("0) Done with this menu - please proceed.\n");
            j = 0;
            for (i = 0; i < num_traits; i++) {
                if (global_trait_entries[i] < 0) {
                    m_ind = i;
                    continue;
                }
                if (i == disease_locus_index) {
                    printf("*");
                } else {
                    printf(" ");
                }
                printf("%d) %s\n", j + 1,
                       LTopp->Pheno[global_trait_entries[i]].TraitName);
                j++;
            }
            printf(
                "Select the disease locus by entering a number from 1 - %d > ",
                ((m_ind < INT_MAX) ? num_traits - 1 : num_traits));
            fcmap(stdin, "%s", choice);
            sscanf(choice, "%d", &choice_);
            if (choice_ == 0)
                break;
            if (choice_ >= 1 && choice_ <= j + 1) {
                if (choice_ <= m_ind) {
                    disease_locus_index = choice_ - 1;
                } else {
                    disease_locus_index = choice_;
                }
            } else {
                warn_unknown(choice);
            }
        }
    }

    return disease_locus_index;

}
