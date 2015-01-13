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

/* write_prest.c
   create files for the option PREST (relationship testing)
   Since the trait locus is not required, let us make this
   simple like allele_frequency computation.

   pedigree file is a pre-makeped format file (requires a dummy
   trait locus)

   locus file is a strippd down linkage-format locus file,
   only has marker loci records and distances between them.

   User has to specify option 1 or 2 (make 2 the default).
*/

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
/*
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  connect_loops get_loci_on_chromosome get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp
        write_files_ext.h:  write_binary_data write_numbered_data
      write_ghfiles_ext.h:  print_recomb_fracs
*/


/* prototypes */

static void save_prest_peds(char *pedfl,
			    linkage_ped_top *Top, const char *genotypes,
			    int chromosome)

{

    int ped, per, loc;
    char poutfl[2*FILENAME_LENGTH];
    FILE *filep;
    char pformat[6], fformat[6];
    int p, pwid, fwid, loc1;
    linkage_ped_rec *Entry;

    field_widths(Top, NULL, &fwid, &pwid, NULL, NULL);
    create_formats(fwid, pwid, fformat, pformat);

    sprintf(poutfl, "%s/%s", output_paths[0], pedfl);

    if ((filep = fopen(poutfl, "w")) == NULL) {
        if (!strcmp(genotypes, "P"))
            errorvf("Unable to open pedigree file '%s'\n", poutfl);
        else
            errorvf("Unable to open genotypes file '%s'\n", poutfl);
        EXIT(FILE_WRITE_ERROR);
    }

    for (ped=0; ped < Top->PedCnt; ped++) {
        if (UntypedPeds != NULL) {
            if (UntypedPeds[ped]) {
                continue;
            }
        }
        for(per=0; per < Top->Ped[ped].EntryCnt; per++) {
            Entry = &(Top->Ped[ped].Entry[per]);
            /* pedigree, person, father, mother, sex and 0 affection status */
            if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                fprintf(filep, fformat, Top->Ped[ped].Name);
            } else if (OrigIds[1] == 3) {
                fprintf(filep, fformat, ped+1);
            } else {
                fprintf(filep, fformat, Top->Ped[ped].Num);
            }

            if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
                fprintf(filep, pformat, Entry->OrigID);
                if (Entry->Father != 0) {
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Father-1].OrigID);
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Mother-1].OrigID);
                }
            } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
                fprintf(filep, pformat, Entry->UniqueID);
                if (Entry->Father != 0) {
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Father-1].UniqueID);
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Mother-1].UniqueID);
                }
            } else {
                fprintf(filep, pformat, Entry->ID);
                if (Entry->Father != 0) {
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Father-1].ID);
                    fprintf(filep, pformat,
                            Top->Ped[ped].Entry[Entry->Mother-1].ID);
                }
            }
            if (Entry->Father == 0){
                /* create the 0 strings */
                for(p=0; p < (pwid-1); p++)    fprintf(filep, " ");
                fprintf(filep, "0 ");
                for(p=0; p< (pwid-1); p++)    fprintf(filep, " ");
                fprintf(filep, "0 ");
            }

            fprintf(filep, "%d 0", Top->Ped[ped].Entry[per].Sex);

            if (chromosome == 0) {
                fprintf(filep,"\n");
                continue;
            }
            /* genotypes file */
            for(loc1=0; loc1 < NumChrLoci; loc1++) {
                loc = ChrLoci[loc1];
                switch(Top->LocusTop->Locus[loc].Type) {
                case AFFECTION:
/*  	  if (!strcmp(genotypes, "P")) { */
/*  	    write_affection_data(filep, loc,  */
/*  				 &(Top->LocusTop->Locus[loc]),  */
/*  				 &(Top->Ped[ped].Entry[per])); */
/*  	  } */

                    break;
                case QUANT:
                    break;

                case BINARY:
                    write_binary_data(filep, loc,
                                      &(Top->LocusTop->Locus[loc]),
                                      &(Top->Ped[ped].Entry[per]));
                    break;

                case NUMBERED:
                    write_numbered_data(filep, loc,
                                        &(Top->LocusTop->Locus[loc]),
                                        &(Top->Ped[ped].Entry[per]));
                    break;

                default:
                    break;
                }
            }
            fprintf(filep,"\n");
        }
    }

    fclose(filep);
    return;

}


//
// PREST (Pedigree RElationship Statistical Test)
// http://fisher.utstat.toronto.edu/sun/Software/Prest/
//
// Detects pedigree errors, cryptic relatedness and relationship misspecification in GWAS or linkage data.
static void prest_locus_file(char *locfl, linkage_ped_top *Top,
			     analysis_type analysis, int chromosome)

{
    int loc, allele, loc1, num_markers;
    char loutfl[2*FILENAME_LENGTH];
    FILE *filep;
    int *markers=CALLOC((size_t) NumChrLoci, int);

    sprintf(loutfl, "%s/%s", output_paths[0], locfl);
    if ((filep = fopen(loutfl, "w")) == NULL) {
        errorvf("Unable to open locus file '%s'\n", loutfl);
        EXIT(FILE_WRITE_ERROR);
    }

    // create an array of only the markers...
    num_markers=0;
    for(loc1=0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        if (Top->LocusTop->Locus[loc].Type == NUMBERED ||
            Top->LocusTop->Locus[loc].Type == BINARY) {
            markers[num_markers] = loc;
            num_markers++;
        }
    }

    fprintf(filep, "%d << no. of markers\n", num_markers);

    for(loc1=0; loc1 < NumChrLoci; loc1++) {
        loc = ChrLoci[loc1];
        switch(Top->LocusTop->Locus[loc].Type) {
        case NUMBERED:
        case BINARY:
	    // For each marker/locu, output the type and allele count...
	    fprintf(filep, "%d %d", (int) Top->LocusTop->Locus[loc].Type - 1,
                    Top->LocusTop->Locus[loc].AlleleCnt);
            if (Top->LocusTop->Locus[loc].Name != NULL)
                fprintf(filep, " # %s", Top->LocusTop->Locus[loc].Name);
            fputc('\n', filep);
	    // On the next line output the frequency of each allele for that marker...
            for (allele = 0; allele < Top->LocusTop->Locus[loc].AlleleCnt; allele++) {
                fprintf(filep, " %.6f",
                        Top->LocusTop->Locus[loc].Allele[allele].Frequency);
            }
            fputc('\n', filep);
            break;

        case QUANT:
        case AFFECTION:
            break;
        default:
            break;
        }
    }

    // [CPK] Prest currently does not model non-autosomal chromosomes, so it makes sense here to use
    // the sex averaged theta (see the Prest documenattion cited above).
    print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, SEX_AVERAGED_THETA);
    free(markers);

    fclose(filep);
    return;

}

static  void prest_chrom_file(char *chromfl, char *poutfl2,
			      char *loutfl, int iter)

{

    FILE *filep;

    char mode[2], chroutfl[2*FILENAME_LENGTH];

    sprintf(chroutfl, "%s/%s", output_paths[0], chromfl);
    ((iter)? strcpy(mode, "a") : strcpy(mode, "w"));

    if ((filep = fopen(chroutfl, mode)) == NULL) {
        errorvf("Unable to open chrom file '%s'\n", chroutfl);
        EXIT(FILE_WRITE_ERROR);
    }

    fprintf(filep, "%s\t%s\n", poutfl2, loutfl);
    fclose(filep);
    return;

}

static void prest_shell_file(char *sh_file, char *pedfile,
			     char *chromfile)

{

    FILE *filep;
    char syscmd[100], shell_file[2*FILENAME_LENGTH];

    sprintf(shell_file, "%s/%s", output_paths[0], sh_file);
    if ((filep = fopen(shell_file, "w")) == NULL) {
        errorvf("Unable to open shell file '%s'\n", shell_file);
        EXIT(FILE_WRITE_ERROR);
    }

    fprintf(filep, "#!/bin/csh -f\n");
    script_time_stamp(filep);
    fprintf(filep, "prest %s %s 1\n", pedfile, chromfile);
    fclose(filep);
    sprintf(syscmd, "chmod +x %s", shell_file);
    System(syscmd);
    return;

}


static void prest_file_names(char *file_names[],
			     analysis_type analysis,
			     int has_orig, int has_uniq)

{

    int item, done=0;
    char cselect[10], fl_stat[12];
    char stem[5];
    int i, first_autosome, num_chromo = main_chromocnt;

    for (i=0; i < main_chromocnt; i++) {
//      commented code will remove SEX_CHROMOSOME from global_chromo_entries ... if desired
//      if (num_chromo != main_chromocnt)
//          global_chromo_entries[i-1] = global_chromo_entries[i];
//      else
        if (global_chromo_entries[i] == SEX_CHROMOSOME) {
            num_chromo=main_chromocnt-1;
//          if (i == main_chromocnt - 1)
            break;
        }
    }

    if (i==0) {
        first_autosome=i+1;
    } else {
        first_autosome = 0;
    }
    if (num_chromo > 1) {
        /* change the extensions on pedigree and chromfile to .all */
        change_output_chr(file_names[0], 0);
        change_output_chr(file_names[3], 0);
        change_output_chr(file_names[4], 0);
        /* remove the extensions on genotypes and locus files */
        change_output_chr(file_names[1], -9);
        change_output_chr(file_names[2], -9);
    }

    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        done = 1;
    }

    while(!done) {
        draw_line();
        print_outfile_mssg();
        printf("PREST file name selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt <= 1) {
            strcpy(stem,"");
            printf(" 1) Pedigree file name:       %-15s\t%s\n",
                   file_names[0], file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:          %-15s\t%s\n",
                   file_names[1], file_status(file_names[1], fl_stat));
            printf(" 3) Genotypes file name:      %-15s\t%s\n",
                   file_names[2], file_status(file_names[2], fl_stat));
            printf(" 4) Chromfiles file name:     %-15s\t%s\n",
                   file_names[3], file_status(file_names[3], fl_stat));
            printf(" 5) Shell file name:          %-15s\t%s\n",
                   file_names[4], file_status(file_names[4], fl_stat));

        } else {
            sprintf(stem, "stem");
            printf(" 1) Pedigree file name:       %-15s\t%s\n",
                   file_names[0], file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name %s:        %-15s\n", stem, file_names[1]);
            printf(" 3) Genotypes file name %s:    %-15s\n", stem, file_names[2]);
            printf(" 4) Chromfiles file name:     %-15s\t%s\n",
                   file_names[3], file_status(file_names[3], fl_stat));
            printf(" 5) Shell file name:          %-15s\t%s\n",
                   file_names[4], file_status(file_names[4], fl_stat));

        }

        individual_id_item(6, analysis, OrigIds[0], 31, 2, 0, 0);
        pedigree_id_item(7, analysis, OrigIds[1], 26, 2, 0);

        printf("Enter option 0-7 > ");

        fflush(stdout);
        (void)fgets(cselect, 9, stdin); newline;
        sscanf(cselect, "%d", &item);
        test_modified(item);

        switch(item) {
        case 0:
            done=1;
            break;

        case 1:
            printf("Enter new pedigree file name > ");
            fcmap(stdin, "%s", file_names[0]); newline;
            break;

        case 2:
            printf("Enter new locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]); newline;
            break;

        case 3:
            printf("Enter new genotypes file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]); newline;
            break;

        case 4:
            printf("Enter new chromfiles file name > ");
            fcmap(stdin, "%s", file_names[3]); newline;
            break;

        case 5:
            printf("Enter new shell file name > ");
            fcmap(stdin, "%s", file_names[4]); newline;
            break;

        case 6:
            OrigIds[0] =
                individual_id_item(0, analysis, OrigIds[0],
                                   35, 1, has_orig, has_uniq);
            break;

        case 7:
            OrigIds[1] = pedigree_id_item(0, analysis, OrigIds[1],
                                          35, 1, has_orig);
            break;

        default:
            warn_unknown(cselect);
            break;
        }
    }

    if (main_chromocnt > 1) {
        change_output_chr(file_names[1], global_chromo_entries[first_autosome]);
        change_output_chr(file_names[2], global_chromo_entries[first_autosome]);
    }

    individual_id_item(0, analysis, OrigIds[0], 0, 3, has_orig, has_uniq);
    pedigree_id_item(0, analysis, OrigIds[1], 0, 3, has_orig);
    return;
}

void create_prest_files(linkage_ped_top *Top,
			analysis_type analysis,
			char *file_names[],
			int untyped_ped_opt,
			int *numchr,
			file_format *infl_type)

{

    int ped, chr;
    int xlinked;

    linkage_ped_tree *LPed;
    linkage_locus_top *LTop = Top->LocusTop;

    if (*infl_type == LINKAGE) {
        for (ped = 0; ped < Top->PedCnt; ped++)  {
            LPed = &(Top->Ped[ped]);
            if (LPed->Loops != NULL)   {
                errorvf("Connecting loops in pedigree %d...\n", LPed->Num);
                if (connect_loops(LPed, Top) < 0)  {
                    errorvf("FATAL: Aborting.\n");
                    EXIT(LOOP_CONNECTION_ERROR);
                }
            }
        }
    }

    prest_file_names(file_names, analysis, Top->OrigIds, Top->UniqueIds);

    /* omit pedigrees based on all selected chromosomes */
    omit_peds(untyped_ped_opt, Top);
    save_prest_peds(file_names[0], Top, "P", 0);

    for(chr=0; chr < main_chromocnt; chr++) {
        *numchr = global_chromo_entries[chr];
        if (main_chromocnt > 1) {
            /* change the extensions on locus and genotype files */
            change_output_chr(file_names[1], *numchr);
            change_output_chr(file_names[2], *numchr);
        }
        xlinked=(((LTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
                  (LTop->SexLinked == 1))? 1 : 0);
        if (xlinked) {
            warnf("Files not created for X-Chromosome.");
            continue;
        }
        if (*numchr == UNKNOWN_CHROMO) {
            get_unmapped_loci(2);
        } else {
            get_loci_on_chromosome(*numchr);
        }
        prest_locus_file(file_names[1], Top, analysis, *numchr);
        save_prest_peds(file_names[2], Top, "G", *numchr);
        prest_chrom_file(file_names[3], file_names[1], file_names[2], chr);
        if (chr == 0) {
            prest_shell_file(file_names[4], file_names[0],
                             file_names[3]);
        }
        if (chr == 0) {
            create_mssg(analysis);
            sprintf(err_msg, "\tPedigree file : %s", file_names[0]);
            mssgf(err_msg);
            sprintf(err_msg, "\tChrom file    : %s", file_names[3]);
            mssgf(err_msg);
        }
        sprintf(err_msg, "\tLocus file    : %s", file_names[1]);
        mssgf(err_msg);
        sprintf(err_msg, "\tGenotype file : %s", file_names[2]);
        mssgf(err_msg);
        if (chr == 0) {
            sprintf(err_msg, "\tShell file    : %s", file_names[4]);
            mssgf(err_msg);
        }
        draw_line();
    }

    return;
}
