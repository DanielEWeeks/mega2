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

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "input_check_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "pedtree_ext.h"
#include "reorder_loci_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc my_malloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  create_marker_sublist delete_file haldane_theta kosambi_theta purge_unneeded_pedigrees
        input_check_ext.h:  check_renumber_ped
            linkage_ext.h:  connect_loops copy_linkage_ped_top free_all_including_lpedtop new_lpedtop
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
            pedtree_ext.h:  convert_to_pedtree free_all_including_ped_top reassign_affecteds remove_untyped_affecteds
       reorder_loci_ext.h:  ReOrderMappedLoci
         user_input_ext.h:  affected_by_status
              utils_ext.h:  EXIT draw_line log_line script_time_stamp
*/



#ifdef USE_PROTOS


void  create_APMULT(char *pedfl_name, char *locusfl_name, int disease_locus,
		    char *affdata_strg, linkage_ped_top **LPedTreeTop,
		    ped_top * PedTreeTop,
		    file_format * infl_type, file_format * outfl_type,
		    int *numchr, char *mapfl_name,
		    analysis_type * analysis, char *file_names[],
		    int untyped_ped_opt);
void  create_APM_file(char *pedfl_name, char *locusfl_name,
		      int disease_locus,
		      const char *affdata_strg, linkage_ped_top * LPedTreeTop,
		      ped_top * PedTreeTop,
		      file_format * infl_type, file_format * outfl_type,
		      int *numchr, char *mapfl_name, char *file_names[],
		      analysis_type * analysis, char *omitfl_name,
		      int untyped_ped_opt);



#endif


static int save_apm_peds_ML(FILE *filep, ped_top *Top, int *save)
{
    int ped, locus, allele, loc, lloc;
    int saved, entry, charcnt, typed;
    ped_tree *Ped;

    /* save the header */
    /* check that ped is valid */
    /* write number of pedigrees, number of loci, and sex-linked keyword */
    saved=0;
    for (ped=0; ped < Top->PedCnt; ped++) {
        if (save[ped]==1) saved++;
    }
    fprintf(filep, "   %d  %d  %s\n", saved, Top->LocusTop->LocusCnt,
            (Top->LocusTop->SexLinked)? X_LINKED_KEYWORD : AUTOSOMAL_KEYWORD);
    for (loc = 0; loc < Top->LocusTop->LocusCnt; loc++) {
        /* write number of alleles and locus name for each locus */
        fprintf(filep, "  %2d   %s\n",
                Top->LocusTop->Locus[loc].AlleleCnt,
                Top->LocusTop->Locus[loc].LocusName);
        /* write allele frequencies for this locus */
        charcnt = 0;
        for (allele = 0; allele < Top->LocusTop->Locus[loc].AlleleCnt; allele++) {
            fprintf(filep, "%1.4g  ",
                    Top->LocusTop->Locus[loc].Allele[allele].Frequency);
            charcnt += 7;
            if (charcnt >= 70) { fputc('\n', filep); charcnt = 0; }
        }
        if (charcnt != 0) fputc('\n', filep);
    }
    /* now save each pedigree */
    for (ped = 0; ped < Top->PedCnt; ped++) {
        if (save[ped]==1) {
            Ped=&(Top->PedTree[ped]);
            fprintf(filep, "Ped %s\n", Ped->Name);
            fprintf(filep, "   %d   %d", Ped->EntryCnt, Ped->AffectedCnt);
            typed=0;
            for (loc = 0; loc < Top->LocusTop->LocusCnt; loc++) {
                lloc = Top->LocusTop->Locus[loc].linkage_loc_num;
                for (entry = 0; entry < Ped->AffectedCnt; entry++) {
                    int a1, a2;
                    get_2alleles(Ped->Affected[entry]->LEntry->Marker, lloc, &a1, &a2);
                    if (a1 != 0 && a2 != 0)
                        { typed++; break; }
                }
            }
            fprintf(filep, "   %d\n", typed);

            /* write mother vector */
            charcnt=0;
            for (entry = 0; entry < Ped->EntryCnt; entry++) {
                if (Ped->Entry[entry].Mother == NULL) fprintf(filep, "   0");
                else fprintf(filep, " %3d", Ped->Entry[entry].Mother->ID);
                charcnt +=4;
                if (charcnt >= 70) {
                    fprintf(filep, "\n"); charcnt=0;
                }
            }
            if (charcnt > 0) fputc('\n', filep);

            /* write father vector */
            charcnt=0;
            for (entry = 0; entry < Ped->EntryCnt; entry++) {
                if (Ped->Entry[entry].Father == NULL) fprintf(filep, "   0");
                else fprintf(filep, " %3d", Ped->Entry[entry].Father->ID);
                charcnt +=4;
                if (charcnt >= 70) {
                    fprintf(filep, "\n"); charcnt=0;
                }
            }
            if (charcnt > 0) fputc('\n', filep);

            /* write sex vector if appropriate */
            charcnt=0;
            if (Top->LocusTop->SexLinked) {
                for (entry = 0; entry < Ped->EntryCnt; entry++) {
                    if (IS_MALE(Ped->Entry[entry])) fprintf(filep, "   1");
                    else if (IS_FEMALE(Ped->Entry[entry])) fprintf(filep, "   2");
                    charcnt +=4;
                    if (charcnt >= 70) {
                        fprintf(filep, "\n"); charcnt=0;
                    }
                }
            }
            if (charcnt > 0) fputc('\n', filep);

            /* write list of all affecteds */

            charcnt = 0;
            for (entry = 0; entry < Ped->AffectedCnt; entry++) {
                fprintf(filep, "   %3d", Ped->Affected[entry]->ID);
                if ((charcnt += 6) >= 70) { fputc('\n', filep); charcnt = 0; }
            }

            if (charcnt > 0) fputc('\n', filep);
            /* write locus number and list of genotypes */
            for (locus = 0; locus < Top->LocusTop->LocusCnt; locus++) {
                /* First check if we are actually typed at this locus. */
                lloc = Top->LocusTop->Locus[loc].linkage_loc_num;
                typed = 0;
                for (entry = 0; (typed == 0) && (entry < Ped->AffectedCnt); entry++) {
                    if ((Ped->Affected[entry]->Allele_1[locus] != 0)
                        && (Ped->Affected[entry]->Allele_2[locus] != 0))
                        typed = 1;
                }
                if (typed == 0) continue;
                fprintf(filep, "%3d", locus+1);

                charcnt = 2;

                /* write affected id's and genotypes */
                for (entry = 0; entry < Ped->AffectedCnt; entry++) {
                    if ((Ped->Affected[entry]->Allele_1[locus] == 0)
                        || (Ped->Affected[entry]->Allele_2[locus] == 0))
                        fprintf(filep, "   0  0");
                    else fprintf(filep, "  %2d %2d",
                                 Ped->Affected[entry]->Allele_1[locus],
                                 Ped->Affected[entry]->Allele_2[locus]);
                    if ((charcnt += 6) >= 70) {
                        fputc('\n', filep); charcnt=0;
                        if (entry+1 < Ped->AffectedCnt) fputc(' ', filep);
                    }
                }
                if (charcnt > 0) fprintf(filep, "\n");
            }
        }
    }

    return 0;
}

static int save_apm_peds_MULT(char *outfl_name, ped_top *Top,
			      linkage_ped_top *LTop,
			      int *order, int num_markers,
			      const char *affdata_strg)
{

    int ped, entry, locus, allele, loc;
    int *trp, tr, nloop, num_affec=num_traits, charcnt;
    int peds_total, peds_saved;
    ped_tree *Ped;
    char *apmult_file;
/*  analysis_type analysis=TO_APM_MULT; */
    FILE *filep;
    int  *save = CALLOC((size_t) Top->PedCnt, int);

    NLOOP;
    trp = &(global_trait_entries[0]);
    apmult_file = CALLOC(strlen(outfl_name)+FILENAME_LENGTH, char);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        affected_by_status(LTop, *trp);
        reassign_affecteds(Top);
        for (ped = 0; ped < Top->PedCnt; ped++) {
            check_renumber_ped(&(Top->PedTree[ped]));
            remove_untyped_affecteds(&(Top->PedTree[ped]),
                                     Top->LocusTop, TYPED_AT_ALL_LOCI, order,
                                     num_markers);
        }
        peds_total = Top->PedCnt;
        peds_saved = peds_total - purge_unneeded_pedigrees(Top, save);

        if (peds_saved == 0) {
            warnf("No pedigrees to save. Output file not created.");
            continue;
        }

        sprintf(apmult_file, "%s/%s", output_paths[tr], outfl_name);

        if ((filep=fopen(apmult_file, "w")) == NULL) {
            errorvf("Unable to open file %s\n", apmult_file);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(filep, "   %d  %d  %s\n", peds_saved, num_markers,
                (Top->LocusTop->SexLinked)?
                X_LINKED_KEYWORD : AUTOSOMAL_KEYWORD);
        for (loc = 0; loc < num_markers; loc++) {
            /* write number of alleles and locus name for each locus */
            locus=order[loc];
            fprintf(filep, "  %2d   %s\n",
                    Top->LocusTop->Locus[locus].AlleleCnt,
                    Top->LocusTop->Locus[locus].LocusName);
            /* write allele frequencies for this locus */
            charcnt = 0;
            for (allele=0; allele <Top->LocusTop->Locus[locus].AlleleCnt; allele++) {
                fprintf(filep, "%1.4g  ",
                        Top->LocusTop->Locus[locus].Allele[allele].Frequency);
                charcnt += 7;
                if (charcnt >= 70) { fputc('\n', filep); charcnt = 0; }
            }
            if (charcnt > 0) fputc('\n', filep);
        }
        /* now save each pedigree if not purged */
        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (save[ped]==1) {
                Ped=&(Top->PedTree[ped]);
                fprintf(filep, "Ped %s\n", Ped->Name);
                fprintf(filep, "   %d   %d\n", Ped->EntryCnt, Ped->AffectedCnt);

                /* write mother vector */
                charcnt=0;
                for (entry = 0; entry < Ped->EntryCnt; entry++) {
                    fputc(' ', filep);
                    if (Ped->Entry[entry].Mother == NULL) fprintf(filep, "   0");
                    else fprintf(filep, " %3d", Ped->Entry[entry].Mother->ID);
                    charcnt +=4;
                    if (charcnt >= 70) {
                        fprintf(filep, "\n"); charcnt=0;
                    }
                }
                if (charcnt > 0) fprintf(filep, "\n");
                /* write father vector */
                charcnt=0;
                for (entry = 0; entry < Ped->EntryCnt; entry++) {
                    fputc(' ', filep);
                    if (Ped->Entry[entry].Father == NULL) fprintf(filep, "   0");
                    else fprintf(filep, " %3d", Ped->Entry[entry].Father->ID);
                    charcnt +=4;
                    if (charcnt >= 70) {
                        fprintf(filep, "\n"); charcnt=0;
                    }
                }
                /* write sex vector if appropriate */
                if (charcnt > 0) fprintf(filep, "\n");
                charcnt=0;
                if (Top->LocusTop->SexLinked) {
                    for (entry = 0; entry < Ped->EntryCnt; entry++) {
                        if (IS_MALE(Ped->Entry[entry])) fprintf(filep, " 1");
                        else
                            if (IS_FEMALE(Ped->Entry[entry])) fprintf(filep, " 2");
                        charcnt +=4;
                        if (charcnt >= 70) {
                            fprintf(filep, "\n"); charcnt=0;
                        }
                    }
                }
                if (charcnt > 0) fprintf(filep, "\n");
                /* write list of all affecteds */
                charcnt = 0;
                for (entry = 0; entry < Ped->AffectedCnt; entry++) {
                    fprintf(filep, "%3d ", Ped->Affected[entry]->ID);
                    charcnt = 4;
                    for (loc = 0; loc < num_markers; loc++) {
                        locus=order[loc];
                        if ((Ped->Affected[entry]->Allele_1[locus] == 0)
                            || (Ped->Affected[entry]->Allele_2[locus] == 0))
                            fprintf(filep, "   0  0");
                        else fprintf(filep, " %2d %2d",
                                     Ped->Affected[entry]->Allele_1[locus],
                                     Ped->Affected[entry]->Allele_2[locus]);
                        if ((charcnt += 6) >= 70) {
                            fputc('\n', filep); charcnt = 0;
                        }
                    }
                    fputc('\n', filep);
                }
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;
    }
    free(apmult_file);
    free(save);
    return 0;
}


/*================== here is future create_APM fcn () =========================
  this fcn will create APM files using the original chapm code and a new user menu interface provided here...*/

void create_APM_file(char *pedfl_name, char *locusfl_name,
		     int disease_locus,
		     const char *affdata_strg, linkage_ped_top * LPedTreeTop,
		     ped_top * PedTreeTop,
		     file_format * infl_type, file_format * outfl_type,
		     int *numchr, char *mapfl_name, char *file_names[],
		     analysis_type * analysis, char *omitfl_name,
		     int untyped_ped_opt)

{

    int     tr, *trp, i, nloop, num_affec= num_traits;
    char    fl_stat[12], outfl_name[2*FILENAME_LENGTH];
    int      peds_total, peds_saved;
    int      reply=-1;
    FILE *filep;
    int            *save;

    /*
     * if affected_by_status is not called, a correct value of AffectedCnt is
     * not determined ...
     */

    if (main_chromocnt > 1)
        change_output_chr(file_names[0], 0);

    if (!DEFAULT_OUTFILES) {
        /* ask user for output file name */
        while(reply != 0) {
            draw_line();
            print_outfile_mssg();
            printf("APM output file name selection menu:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) Output file name:    %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf("Enter option 0 or 1 > ");
            fcmap(stdin, "%d", &reply); newline;
            if (reply==1) {
                printf("Enter new output file name > ");
                fcmap(stdin, "%s", file_names[0]); newline;
            }
        }
    }

    /* omit untyped peds */
    omit_peds(untyped_ped_opt, LPedTreeTop);

    PedTreeTop = convert_to_pedtree(LPedTreeTop, 0);

    save = CALLOC((size_t) PedTreeTop->PedCnt, int);
    NLOOP;
    trp = &(global_trait_entries[0]);
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(outfl_name, "%s/%s", output_paths[tr], file_names[0]);
        affected_by_status(LPedTreeTop, *trp);
        reassign_affecteds(PedTreeTop);
        for (i = 0; i < PedTreeTop->PedCnt; i++) {
            check_renumber_ped(&(PedTreeTop->PedTree[i]));
        }
        for (i = 0; i < PedTreeTop->PedCnt; i++)
            remove_untyped_affecteds(&(PedTreeTop->PedTree[i]),
                                     PedTreeTop->LocusTop,
                                     TYPED_AT_ANY_LOCUS, NULL, 0);
        peds_total = PedTreeTop->PedCnt;
        peds_saved = peds_total - purge_unneeded_pedigrees(PedTreeTop, save);
        if (peds_saved == 0) {
            printf("WARNING: No pedigrees to save for disease locus %s.\n",
                   LPedTreeTop->LocusTop->Pheno[*trp].TraitName);
        } else {
            if ((filep = fopen(outfl_name, "w")) == NULL) {
                errorvf("Could not open file %s.\n", outfl_name);
                EXIT(FILE_WRITE_ERROR);
            }
            save_apm_peds_ML(filep, PedTreeTop, save);
            printf("Saved %d pedigrees out of %d for disease locus %s.\n",
                   peds_saved, peds_total, LPedTreeTop->LocusTop->Pheno[*trp].TraitName);
        }
        if (nloop == 1) break;
        trp++;
    }

    log_line(mssgf);
    create_mssg(*analysis);
    sprintf(err_msg, "      Data file:    %s", file_names[0]);
    mssgf(err_msg);
    log_line(mssgf);

    free(save);
}


static void  create_apmult_script_file(linkage_ped_top *Copy,
				       char *outfl_name,
				       char *apmult_sum, char *script_fl,
				       int first_time,
				       int num_markers, int *markers)
{

    int             i, tr, nloop, num_affec= num_traits;
    int             *numbered, *num_p;
    double          dx, theta;
    FILE            *fp;
    char            script_name[2*FILENAME_LENGTH], syscmd[200];


    numbered=CALLOC((size_t) Copy->LocusTop->LocusCnt, int);
    num_p=&(numbered[0]);

    for(i=0; i< Copy->LocusTop->LocusCnt; i++) {
	// AMP only analyze autosomal data, so XLINKED and YLINKED Loci should not be included...
        if (Copy->LocusTop->Locus[i].Type == NUMBERED ||
           Copy->LocusTop->Locus[i].Type == BINARY) {
            *num_p = i; num_p++;
        }
    }

    NLOOP;

    for (tr=0; tr <=nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(script_name, "%s/%s", output_paths[tr], script_fl);
        if (first_time) {
            delete_file(script_name);
            if ((fp = fopen(script_name, "w")) == NULL)    {
                errorvf("Unable to open file %s\n", script_name);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(fp, "#!/bin/csh -f\n");
            fprintf(fp, "# apmmult script: %s\n", script_fl);
            script_time_stamp(fp);
            fprintf(fp, "unalias rm\n");
            fprintf(fp, "echo \"Starting multilocus APM run - \"\n");
            fprintf(fp,
                    "echo \"    Results will be in apmmult.out and table.out\"\n");
            fprintf(fp, "echo \"    Summary will be in %s\"\n", apmult_sum);
            fprintf(fp, "touch %s\n", apmult_sum);
            fprintf(fp, "date >> %s\n", apmult_sum);
        } else {
            if ((fp = fopen(script_name, "a")) == NULL)    {
                errorvf("Unable to open file %s\n", script_name);
                EXIT(FILE_WRITE_ERROR);
            }
        }
        fprintf(fp, "if (-e %s) then\n", outfl_name);
        fprintf(fp, "rm -f table.out\n");
        fprintf(fp, "echo \"Running apmmult on %s\"\n", outfl_name);
        fprintf(fp, "apmmult << DATA > /dev/null\n");
        fprintf(fp, "%s\n", outfl_name);
        /* print the haldane thetas */
        fprintf(fp, "0\n");
        for (i = 1; i < num_markers; i++) {
	  // what about SEX_SPECIFIC??
            dx = Copy->LocusTop->Marker[numbered[markers[i]]].pos_avg -
                Copy->LocusTop->Marker[numbered[markers[i - 1]]].pos_avg;
            theta = ((Copy->LocusTop->map_distance_type == 'h')?
                     haldane_theta(dx): kosambi_theta(dx));
            fprintf(fp, "%f\n", theta);
        }

        fprintf(fp, "DATA\n");
        fprintf(fp, "grep \"Data file\" table.out >> %s\n", apmult_sum);

        for (i = 1; i < num_markers; i++) {
            dx = Copy->LocusTop->Marker[numbered[markers[i]]].pos_avg -
                Copy->LocusTop->Marker[numbered[markers[i - 1]]].pos_avg;
            /*      theta = ((Copy->LocusTop->map_distance_type == 'h')?
                    haldane_theta(dx): kosambi_theta(dx)); */
            fprintf(fp, "echo \"%10s - %7.5f - %-10s\" >> %-s\n",
                    Copy->LocusTop->Marker[numbered[markers[i-1]]].MarkerName,
                    dx,
                    Copy->LocusTop->Marker[numbered[markers[i]]].MarkerName, apmult_sum);
        }

        fprintf(fp, "grep \"Theta\" table.out >> %s\n", apmult_sum);
        fprintf(fp,
                "%s 'BEGIN { n = 0 } { if ($1 ==\"Overall\") n += 1 ; if ( n>0 ) print $0 }' table.out >> %s\n",
                awk_str, apmult_sum);
        fprintf(fp, "endif\n");
        fclose(fp);
        sprintf(syscmd, "chmod +x %s", script_name);
        System(syscmd);

        if (nloop == 1) break;
    }
    free(numbered);
    return;

}

/*      main routine for APM-MULT */

#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*Top)

void create_APMULT(char *pedfl_name, char *locusfl_name, int disease_locus,
		   const char *affdata_strg, linkage_ped_top **Top,
		   ped_top * PedTreeTop,
		   file_format * infl_type, file_format *outfl_type,
		   int *numchr, char *mapfl_name,
		   analysis_type *analysis, char *file_names[],
		   int untyped_ped_opt)
{

    linkage_ped_tree *LPed;
    char    *apmult_file, *new_name_part;
    int     ch, *order, select=-1, NumberOfMarkers, i, j, first_time;
    int of;
    loci_markers   *markerset;
    char chr[3], stem[5], fl_stat[12];
    linkage_ped_top *Copy, *Copy_marker;

    char **apm_outfile_names, **apm_script_names;
    int *num_files;

    /* apm file name selction */
    change_output_chr(file_names[0], -9);

    if (!DEFAULT_OUTFILES) {
        while (select != 0) {
            draw_line();
            print_outfile_mssg();
            printf("APMult-format file name selection\n");
            printf("0) Done with this menu - please proceed.\n");
            if (main_chromocnt == 1) {
                printf(" 1) Pedigree file name stem:    %-15s\n", file_names[0]);
                strcpy(stem, "");
                printf(" 2) C-shell script file:        %-15s\t%s\n", file_names[3],
                       file_status(file_names[3], fl_stat));
                printf(" 3) Summary file name:          %-15s\t%s\n", file_names[14],
                       file_status(file_names[14], fl_stat));
            } else {
                strcpy(stem, "stem");
                change_output_chr(file_names[3], -9);
                change_output_chr(file_names[14], -9);
                printf(" 1) Pedigree file name stem:        %-15s\n", file_names[0]);
                printf(" 2) C-shell script file name stem:  %s\n", file_names[3]);
                printf(" 3) Summary file name stem:         %s\n", file_names[14]);
            }
            printf("Enter selection: 0 - 3 > ");
            fcmap(stdin, "%d", &select); newline;
            switch (select)  {
            case 0:
                break;
            case 1:
                printf("Enter pedigree file name stem > ");
                fcmap(stdin, "%s", file_names[0]); newline;
                break;
            case 2:
                printf("Enter C-shell script file name %s > ", stem);
                fcmap(stdin, "%s", file_names[3]); newline;
                break;
            case 3:
                printf("Enter summary file name > %s ", stem);
                fcmap(stdin, "%s", file_names[14]); newline;
                break;
            default:
                printf("Unknown option %d\n", select);
                break;
            }
        }
        if (main_chromocnt > 1) {
            /* .sh was stripped, so put a dummy chromosome + .sh */
            sprintf(file_names[3], "%s.01.sh", file_names[3]);
        }
    }

    if (*infl_type == LINKAGE) {
        for (i = 0; i < LPedTreeTop->PedCnt; i++)   {
            LPed = &(LPedTreeTop->Ped[i]);
            if (LPed->Loops != NULL)       {
                sprintf(err_msg, "Connecting loops in pedigree %d...", LPed->Num);
                mssgf(err_msg);
                if (connect_loops(LPed, LPedTreeTop) < 0)     {
                    errorvf("Fatal error connecting loops! Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    NumberOfMarkers = 2; select = -1;
    if (DEFAULT_OPTIONS) {
        mssgf("Number of markers in each sub-order set to 2 (default)");
        mssgf("as requested in batch file.");
    } else {
        while (select != 0) {
            draw_line();
            printf("Select number of markers in each subset:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) Number of markers in suborder: %d \n", NumberOfMarkers);
            printf("Enter 0 or 1 > ");
            fcmap(stdin, "%d", &select); newline;
            if (select == 1) {
                printf(
                    "Enter number of markers in suborder (excluding disease locus) > ");
                fcmap(stdin, "%d", &NumberOfMarkers);      newline;
                if (NumberOfMarkers < 2) {
                    printf("WARNING: number of markers less than 2 \n");
                    NumberOfMarkers=2;
                }
            } else if (select != 0) {
                printf("Unknown option %d.\n", select);
            }
        }
    }

    Copy=new_lpedtop();
    copy_linkage_ped_top(LPedTreeTop, Copy, 1);
    /* allocate some pointers for maximum possible number of files */
    apm_outfile_names = CALLOC((size_t) Copy->LocusTop->LocusCnt,  char *);
    apm_script_names = CALLOC((size_t) main_chromocnt, char *);
    num_files = CALLOC((size_t) main_chromocnt, int);

    of=0;
    for (ch=0; ch < main_chromocnt; ch++) {
        if (main_chromocnt > 1) {
            *numchr=global_chromo_entries[ch];
            abort();
#ifdef DEFUNCT
          NOTE: luckily apm is not currently supported.  and ReOrderMappedLoci has ceased to exist.
                If you need to bring apm back then, Youll have to figure out hot to use
                ReOrderMappedLoci_new() ... like everyone else
            ReOrderMappedLoci(LPedTreeTop, numchr);
#endif
#ifdef MOREDEFUNCT
            to silence abort() above
            change_output_chr(file_names[3], *numchr);
            change_output_chr(file_names[14], *numchr);
#endif
        }
        CHR_STR(*numchr, chr);

        if (NumberOfMarkers > (LPedTreeTop->LocusTop->LocusCnt))
            NumberOfMarkers = LPedTreeTop->LocusTop->LocusCnt;

        markerset = MALLOC(loci_markers);
        create_marker_sublist(LPedTreeTop, markerset, NumberOfMarkers);
        new_name_part = CALLOC((size_t) 8, char);
        order = CALLOC((size_t) NumberOfMarkers, int);

        Copy_marker=new_lpedtop();
        copy_linkage_ped_top(LPedTreeTop, Copy_marker,1);

        apmult_file=CALLOC((size_t) FILENAME_LENGTH+8, char);
        first_time=1;
        for (i = 0; i < markerset->total_files; i++) {
            for (j = 0; j < NumberOfMarkers; j++) {
                order[j]=markerset->MarkerSubsetMatrix[i][j];
            }
            sprintf(new_name_part, "%d_%d", order[0]+1,
                    order[NumberOfMarkers-1]+1);
            sprintf(apmult_file, "%s%s.%s", file_names[0], new_name_part, chr);
            /* omit pedigrees here */
            omit_peds(untyped_ped_opt, LPedTreeTop);
            PedTreeTop = convert_to_pedtree(LPedTreeTop, 0);
            save_apm_peds_MULT(apmult_file, PedTreeTop, LPedTreeTop, order,
                               NumberOfMarkers, affdata_strg);
            free_all_including_ped_top(PedTreeTop,  NULL, NULL);
            create_apmult_script_file(LPedTreeTop, apmult_file, file_names[14],
                                      file_names[3], first_time,
                                      NumberOfMarkers, order);
            first_time=0;
            apm_outfile_names[of] = CALLOC(strlen(apmult_file)+1, char);
            strcpy(apm_outfile_names[of], apmult_file);
            of++;
        }
        num_files[ch] = of;
        apm_script_names[ch] = CALLOC(strlen(file_names[3])+1, char);
        strcpy(apm_script_names[ch], file_names[3]);
 
        free(new_name_part);
        free(order); free(apmult_file);

        for (i = 0; i < markerset->total_files; i++)
            free(markerset->MarkerSubsetMatrix[i]);
        free(markerset);

        if (main_chromocnt > 1) {
            delink_linkage_ped_top(LPedTreeTop);
            free_all_including_lpedtop(&LPedTreeTop);
            LPedTreeTop=new_lpedtop();
            copy_linkage_ped_top(Copy, LPedTreeTop, 1);
        }
        delink_linkage_ped_top(Copy_marker);
        free_all_including_lpedtop(&Copy_marker); Copy_marker=NULL;
    }
    /* display names of output files */
    log_line(mssgf);
    create_mssg(*analysis);
    of=0;
    for (i = 0; i < main_chromocnt; i++) {
        for (j = of; j < num_files[i]; j++) {
            sprintf(err_msg, "    Data file:   %s", apm_outfile_names[j]);
            mssgf(err_msg);
        }
        of = num_files[i];
        sprintf(err_msg, "    Script file: %s", apm_script_names[i]);
        mssgf(err_msg);
        log_line(mssgf);
    }
    for (i=0; i < num_files[main_chromocnt-1]; i++) {
        free(apm_outfile_names[i]);
    }
    free(apm_outfile_names);

    for (i=0; i < main_chromocnt; i++) {
        free(apm_script_names[i]);
    }
    free(apm_script_names);

    delink_linkage_ped_top(Copy);
    free_all_including_lpedtop(&Copy);
    return;

}
#undef LPedTreeTop
