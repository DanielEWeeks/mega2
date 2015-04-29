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

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "write_ghfiles_ext.h"
/*
     error_messages_ext.h:  errorf mssgf warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file
            linkage_ext.h:  get_loci_on_chromosome new_lpedtop switch_penetrances
           makenucs_ext.h:  makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  change_output_chr create_mssg file_status print_outfile_mssg
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT chmod_X_file draw_line global_ou_files script_time_stamp
        write_files_ext.h:  save_linkage_peds write_affection_tdtmax
      write_ghfiles_ext.h:  write_gh_locus_file
*/



/* prototypes */

#ifdef USE_PROTOS
void  create_sib_phase_graph_scripts(char *asp_in_name, char *asp_dat_name,
				     double risk,
				     FILE *fp, char *numchr,
				     int count, int xlinked);

void  create_TDTMAX(linkage_ped_top **Top,
		    analysis_type *analysis, int *numchr,
		    char *file_names[], int untyped_ped_opt,
		    linkage_ped_top **NukeTop);

void  create_LOD2_format_files(linkage_ped_top **LPedTreeTop, ped_top *PedTreeTop,
			       analysis_type *analysis, int *numchr,
			       char *file_names[], int untyped_ped_opt);

/* --------------------------------------------------------------------------*/

#endif

void create_sib_phase_graph_scripts(char *asp_in_name, char *asp_dat_name,
				    double risk,
				    FILE *fp, char *chr_str,
				    int count, int xlinked)
{
    /* modified by Nandita 7/21/99 */

    char graph_file[2*FILENAME_LENGTH], aspin_name1[2*FILENAME_LENGTH];

    /* graph_file name is either gph_lod.<chr> or gph.<risk>.<chr> */
    if (count ==0) {
        sprintf(graph_file, "gph_lod.%s", chr_str);
        strcpy(aspin_name1, asp_in_name);
    } else {
        sprintf(graph_file, "gph.%4.2f.%s", risk, chr_str);
        sprintf(aspin_name1, "%s.%4.2f.%s", asp_in_name, risk, chr_str);
    }

    /* common lines */
    fprintf(fp, "/bin/rm -f %s\n", graph_file);

    /* Need the nawk thingy at the beginning */
    fprintf(fp, "echo \"Running sib_phase -v -f %s %s to create %s\"\n",
	    aspin_name1, asp_dat_name, graph_file);
    if (count==0) {
        if (xlinked)
            fprintf(fp,
                    "sib_phase -v -f %s %s | %s '{ print $1,\" \", $5}' > %s\n",
                    aspin_name1, asp_dat_name, awk_str, graph_file);
        else
            fprintf(fp,
                    "sib_phase -v -f %s %s | %s '{ print $1,\" \", $6}' > %s\n",
                    aspin_name1, asp_dat_name, awk_str, graph_file);
    } else {
        fprintf(fp, "sib_phase -v -f %s %s > %s\n", aspin_name1, asp_dat_name,
                graph_file);
    }

    fprintf(fp, "xmgr_map %s >> %s\n", aspin_name1, graph_file);

    /* more common stuff */
    fprintf(fp, "echo \"@    legend on\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend loctype view\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend layout 0\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend vgap 2\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend hgap 1\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend length 4\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box on\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box fill on\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box fill with color\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box fill color 0\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box fill pattern 1\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box color 1\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box linewidth 2\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend box linestyle 1\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend x1 0.87\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend y1 0.8\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend font 4\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend char size 0.760000\" >>  %s\n", graph_file);
    fprintf(fp, "echo \"@    legend linestyle 1\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend linewidth 1\" >> %s\n", graph_file);
    fprintf(fp, "echo \"@    legend color 1\" >> %s\n", graph_file);

    if (count ==0)
        fprintf(fp, "printf \"@    legend string %d  \\042 MLS \\042\\n\" >> %s\n",
                count, graph_file);
    else
        fprintf(fp, "printf \"@    legend string %d  \\042 %5.2f\\042\\n\" >> %s\n",
                count, risk, graph_file);
    return;
}

#ifdef Top
#undef Top
#endif
#define Top (*LPTop)

/*--------------------------------------------------------------------*/
void create_TDTMAX(linkage_ped_top **LPTop,
		   analysis_type *analysis,
		   int *numchr, char *file_names[],
		   int untyped_ped_opt,
		   linkage_ped_top **NukeTop)
{

    char            cchoice[10], allone[2][4]={"all", "one"};
    int             opt, i, j, j1, m, tr, choice_;
    linkage_locus_top *LTop;
    linkage_ped_rec *tpe;
    FILE           *fp, *cfp;
    int             nperm;
    char            pedfl[2*FILENAME_LENGTH], shfl[2*FILENAME_LENGTH];
    char            fl_stat[12], infl[FILENAME_LENGTH];
    linkage_ped_top *Top2;
    int             nloop, num_affec=num_traits, *trp;
    int             al1, al2;

    LTop = Top->LocusTop;
    if (main_chromocnt > 1) change_output_chr(file_names[3], 0);
    if (DEFAULT_OUTFILES) {
        choice_ = 0;
    } else {
        choice_ = -1;
    }
    while (choice_ != 0) {
        draw_line();
        print_outfile_mssg();
        printf("TDTMAX-format file name selection menu:\n");
        draw_line();
        printf("0) Done with menu - proceed program.\n");
        printf(" 1) Pedigree file name stem: %s\n", file_names[0]);
        printf(" 2) C-shell file name :      %s\t%s\n", file_names[3],
               file_status(file_names[3], fl_stat));
        printf("Enter selection:  0, 1 or 2 > ");
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice_);
        switch(choice_) {
        case 1:
            printf("Enter NEW TDTMAX pedigree file name stem: ");
            fcmap(stdin, "%s", file_names[0]); newline;
            break;
        case 2:
            printf("Enter NEW TDTMAX C-shell file name: ");
            fcmap(stdin, "%s", file_names[3]); newline;
            break;
        case 0:
            break;
        default:
            warn_unknown(cchoice);
            break;
        }
    }

    opt=-1;
    choice_ = 0;
    nperm=200;

    if (DEFAULT_OPTIONS) {
        mssgf("TDTMAX analysis parameters set to default values");
        mssgf("as requested in the batch file.");
    } else {
        while (opt != 0) {
            printf("TDTMAX analysis options menu:\n");
            printf("0) Done with this menu - please proceed.\n");
            printf(" 1) All affected sibs/one affected sib per family [%s]\n",
                   allone[choice_]);
            printf(" 2) Number of permutations : %d \n", nperm);
            printf("Enter option 1 to toggle, 2 to change value > ");

            fcmap(stdin, "%s", cchoice); newline;
            sscanf(cchoice, "%d", &opt);
            if (opt < 0 || opt > 2) {
                warn_unknown(cchoice);
            } else if (opt == 1) {
                choice_ = TOGGLE(choice_);
            } else if (opt == 2) {
                printf("Enter number of permutations > ");
                (void)scanf("%d", &nperm); newline;
            }
        }
    }



    printf("OPTION: Create TDTMAX-format file.\n");
    /* omit untyped peds */
    Top2=new_lpedtop();
    makenucs1(Top, NULL, Top2);
    omit_peds(untyped_ped_opt, Top2);
    LTop=Top->LocusTop;
    NLOOP;

    create_mssg(TO_TDTMAX);

    trp= &(global_trait_entries[0]);

    for (tr=0;  tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;
        /* assign the pedigree and cshell file names for each trait */
        sprintf(shfl, "%s/%s", output_paths[tr], file_names[3]);
        delete_file(shfl);
        /* write the header into the shell file */
        if ((cfp = fopen(shfl, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s'\n", shfl);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(cfp, "#!/bin/csh -f \n");
        script_time_stamp(cfp);
        fprintf(cfp, "touch tdtmax.lst\n");
        fprintf(cfp, "date >> tdtmax.lst\n");
        fprintf(cfp, "unalias rm\n");
        if (choice_ == 0)
            fprintf(cfp,
                    "echo \"Analysis option: all affected offspring in families\" >> tdtmax.lst\n");
        else
            fprintf(cfp,
                    "echo \"Analysis option:  only one affected sib per family\" >> tdtmax.lst\n");

        /* now write a pedigree file for each marker
           and corresponding lines into shell file */
        for (j1 = 0; j1 < num_reordered; j1++)   {
            j=reordered_marker_loci[j1];
            if (LTop->Locus[j].Type != NUMBERED &&
                LTop->Locus[j].Type != BINARY) continue;
            sprintf(infl, "%s.%s", file_names[0], LTop->Locus[j].Name);
            sprintf(pedfl, "%s/%s", output_paths[tr], infl);

            if ((fp = fopen(pedfl, "w")) == NULL)	   {
                draw_line();
                errorvf("Unable to open file '%s'\n", pedfl);
                EXIT(FILE_WRITE_ERROR);
            }
            for (m = 0; m < Top2->PedCnt; m++)	   {
                for (i = 0; i < Top2->Ped[m].EntryCnt; i++)	  {
                    tpe = &(Top2->Ped[m].Entry[i]);
                    if ((m+1) > 9999)  {
                        errorf("Number of families has to be less than 9999");
                        errorf("for 'convert' which is a part of the tdtmax package");
                        EXIT(OUTPUT_FORMAT_ERROR);
                    }
                    if (tpe->Father != 0)
                        fprintf(fp, "%4d%3d%2d%2d%2d  ",
                                m+1,
                                tpe->ID,
                                tpe->Father,
                                tpe->Mother,
                                tpe->Sex);
                    else
                        fprintf(fp, "%4d%3d %c %c%2d  ",
                                m+1,
                                tpe->ID,
                                'x',
                                'x',
                                tpe->Sex);
                    write_affection_tdtmax(fp, *trp, &(Top2->LocusTop->Locus[*trp]),
                                           tpe);

                    get_2alleles(Top2->Ped[m].Entry[i].Marker, j, &al1, &al2);
                    if (al1 != 0) {
                        linkage_locus_rec *locus = &Top2->LocusTop->Locus[*trp];
                        fprintf(fp, " %2s %2s\n", format_allele(locus, al1), format_allele(locus, al2));
                    } else
                        fprintf(fp,"  x  x\n");
                }
            }
            fclose(fp);
            if (tr <= 1) {
                sprintf(err_msg, "        Pedigree file: %s", infl);
                mssgf(err_msg);
            }

            /* now print analysis lines into the shell file */
            fprintf(cfp, "echo Running TDTMAX on the marker %s\n",
                    Top2->LocusTop->Locus[j].Name);
            fprintf(cfp, "rm -f data.txt permdata permout permute\n");
            fprintf(cfp, "cp %s data.txt\n", infl);
            fprintf(cfp, "convert << DATA > /dev/null\n");
            fprintf(cfp, "3\n");   /* Implies 4 columns for ped id */
            fprintf(cfp, "%d\n",Top2->LocusTop->Locus[j].AlleleCnt); /* Number of alleles */
            fprintf(cfp, "%d\n", (choice_)? 2 : 1);   /* Option: all vs 1 sib */
            fprintf(cfp, "DATA\n");
            fprintf(cfp, "set s = `wc permdata`\n");
            fprintf(cfp, "echo $s[1] parents in the sample\n");
            fprintf(cfp, "tdtmax << DATA2 > /dev/null\n");
            fprintf(cfp, "%d\n",Top2->LocusTop->Locus[j].AlleleCnt); /* Number of alleles */
            fprintf(cfp, "$s[1]\n");
            fprintf(cfp, "%d\n",nperm);   /* Number of permutations */
            fprintf(cfp, "DATA2\n");
            fprintf(cfp, "echo \"Locus: %s  %10.7f  %d\" >> tdtmax.lst\n",Top2->LocusTop->Marker[j].Name, Top2->LocusTop->Marker[j].pos_avg, Top2->LocusTop->Marker[j].chromosome);
            fprintf(cfp, "cat permout >> tdtmax.lst\n");
            fprintf(cfp, "echo \"----------------------------------\" >> tdtmax.lst\n");
        }  /* for j */
        /* print the footer into the shell file */
        fprintf(cfp, "echo TDTMAX analyses are completed.\n");
        fprintf(cfp, "touch tdtmax.sum\n");
        fprintf(cfp, "date >> tdtmax.sum\n");
        fprintf(cfp, "# Create a summary TDTMAX file\n");
        fprintf(cfp, "%s '\\\n", awk_str);
        fprintf(cfp, "BEGIN { print(\"Locus         Position  Chr     P-value\") } \\\n");
        fprintf(cfp, "/Locus/ { printf(\"%%-15s\\t%%s\\t%%s  \", $2,$3,$4) } \\\n");
        fprintf(cfp, "/value/ { printf(\"\\t%%s\\n\",$7 ) }' tdtmax.lst >> tdtmax.sum\n");
        fprintf(cfp, "echo Please see tdtmax.lst for the detailed results,\n");
        fprintf(cfp, "echo and tdtmax.sum for a summary of the results.\n");
        fclose(cfp);
        /* go on to the next trait */
        if (nloop == 1) break;
        trp++;

    }
    chmod_X_file(shfl);
    sprintf(err_msg, "        C-shell script file: %s", file_names[3]);
    mssgf(err_msg);

    draw_line();
    *NukeTop = Top2;
    /* need to free Top2 here */
    /*  free_all_including_lpedtop(&Top2); */
    /* will be freed outside */

    return;

}

#undef Top

static void lod2_file_names(int single_file, char *file_names[],
			    analysis_type *analysis)

{
    int user_reply=-1;
    char creply[10], stem[5], fl_stat[12];

    // Here -9 is a flag to return only the file name portion...
    // Here 0 is a flag that says to replace <extension> with 'all', keeping <extension> and <rest> if they exist...
    (*analysis)->replace_chr_number(file_names, (!single_file && main_chromocnt > 1 ? -9 : 0));

    if (DEFAULT_OUTFILES) {
        user_reply=0;
    }

    while (user_reply != 0) {
        draw_line();
        print_outfile_mssg();
        printf("Lod2-format file names selection menu:\n");
        draw_line();
        printf("0) Done with this menu - please proceed.\n");
        if (!single_file && main_chromocnt > 1) {
            strcpy(stem, "stem");
            printf(" 1) Pedigree file name stem:     %-15s\n", file_names[0]);
            printf(" 2) Locus file name stem:        %-15s\n", file_names[1]);
            printf(" 3) C-shell file name stem:      %-15s\n", file_names[3]);
        } else {
            strcpy(stem, "");
            printf(" 1) Pedigree file name:          %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:             %-15s\t%s\n", file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 3) C-shell file  :              %-15s\t%s\n", file_names[3],
                   file_status(file_names[3], fl_stat));
        }
        printf("Enter selection: 0, 1, 2, or 3 > ");
        fcmap(stdin, "%s", creply); newline;
        sscanf(creply, "%d", &user_reply);
        test_modified(user_reply);

        switch (user_reply)  {
        case 1:
            printf("Enter NEW pedigree file name %s: > ", stem);
            fcmap(stdin, "%s", file_names[0]); newline;
            break;
        case 2:
            printf("Enter NEW locus file name %s: > ", stem);
            fcmap(stdin, "%s", file_names[1]); newline;
            break;
        case 3:
            printf("Enter NEW C-shell file name %s: > ", stem);
            fcmap(stdin, "%s", file_names[3]); newline;
            break;
        case 0:
            break;
        default:
            warn_unknown(creply);
            break;
        }
    }

    if (!single_file && main_chromocnt > 1) {
        strcat(file_names[3], ".01.sh");
    }
    return;
}

/* the locus file contains one affection locus at the begining,
   then markers. This should match the order in the shell file.
*/

static void lod2_cshell_file(char *shellfl, char *outfl_name,
			     char *loutfl_name,
			     linkage_locus_top *LTopp, int *numchr)

{
    int i1, i, tr, nloop, num_affec=num_traits;
    char lod2_shell_name[2*FILENAME_LENGTH];
    FILE *lod2_shell_p;
    int chromosome, marker;
    char chrstr[4];

    NLOOP;

    if (numchr != NULL) {
        chromosome = *numchr;
        sprintf(chrstr, "%d", *numchr);
    } else {
        chromosome = 0;
        strcpy(chrstr, "all");
    }

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        marker=2;
        sprintf(lod2_shell_name, "%s/%s", output_paths[tr], shellfl);
        delete_file(lod2_shell_name);
        if ((lod2_shell_p = fopen(lod2_shell_name, "w")) == NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", lod2_shell_name);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(lod2_shell_p, "#!/bin/csh -f\n");
        fprintf(lod2_shell_p,
                "# Compute lod2 lod scores under heterogeneity \n");
        script_time_stamp(lod2_shell_p);
        fprintf(lod2_shell_p,
                "# Output: homog.final, final.out, homog.%s.sum\n", chrstr);
        if (HasQuant) {
            if (ITEM_READ(Value_Missing_Quant_On_Output))
              fprintf(lod2_shell_p,
                      "echo Unknown quantitative phenotype value is %s\n",
		      Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
            fprintf(lod2_shell_p,
                    "echo \"Okay to run homog with this unknown value?[yes/no] (default no)\"\n");
            fprintf(lod2_shell_p,
                    "set continue = $<\n");
            fprintf(lod2_shell_p,
                    "if ($continue != \"yes\") then\n");
            fprintf(lod2_shell_p, "   exit(0)\n");
            fprintf(lod2_shell_p, "endif\n");
        }

        fprintf(lod2_shell_p, " if !(-e %s) then\n", outfl_name);
        fprintf(lod2_shell_p, "  echo ERROR input file %s not found\n",
                outfl_name);
        fprintf(lod2_shell_p, "  exit(2)\n");
        fprintf(lod2_shell_p, " endif\n");
        fprintf(lod2_shell_p, " unalias rm \n");
        fprintf(lod2_shell_p, " rm -f incode.dat \n");
        fprintf(lod2_shell_p, " touch incode.dat \n");
        fprintf(lod2_shell_p, " if (-e final.out) then\n");
        fprintf(lod2_shell_p, "    mv final.out final.out.bak\n");
        fprintf(lod2_shell_p, " endif\n");
        fprintf(lod2_shell_p, " rm -f homog.final\n");
        fprintf(lod2_shell_p, " date >> homog.final\n");
        fprintf(lod2_shell_p, " pwd >> homog.final\n");
        for (i1 = 0; i1 < NumChrLoci; i1++) {
            i = ChrLoci[i1];
            if (LTopp->Locus[i].Type == NUMBERED ||
                LTopp->Locus[i].Type == BINARY) {
                chromosome=LTopp->Marker[i].chromosome;
                fprintf(lod2_shell_p, " echo Analyzing marker %s\n",
                        LTopp->Locus[i].Name);
                fprintf(lod2_shell_p, " rm -f outfile.dat \n");
                fprintf(lod2_shell_p,
                        " lsp mlink %s %s 2 1   %d 0 0 1 0.02 0.48 0 \n",
                        outfl_name, loutfl_name, marker);
                fprintf(lod2_shell_p, " rm -f locus.dat \n");
                fprintf(lod2_shell_p, " grep Order lsp.log > locus.dat \n");
                fprintf(lod2_shell_p, " locusorder\n");
                fprintf(lod2_shell_p, " cat lsp.log >> final.out \n");
                fprintf(lod2_shell_p, " unknown > /dev/null \n");
                fprintf(lod2_shell_p, " mlink > /dev/null \n");
                fprintf(lod2_shell_p, " rm -f HOMOG.DAT.%d.%d\n", chromosome, marker);
                fprintf(lod2_shell_p, " rm -f homog.out.%d.%d\n", chromosome, marker);
                fprintf(lod2_shell_p, " rm -f HOMOG.DAT \n");
                fprintf(lod2_shell_p, " rm -f homog.out \n");
                fprintf(lod2_shell_p, " rm -f HOMOG.OUT \n");
                fprintf(lod2_shell_p, " homogpre > /dev/null \n");
                fprintf(lod2_shell_p, " homog > /dev/null \n");
                fprintf(lod2_shell_p, " homogpos > /dev/null \n");
                fprintf(lod2_shell_p, " cat lsp.log >> homog.final \n");
                fprintf(lod2_shell_p, " cat HOMOG.OUT >> homog.final \n");
                fprintf(lod2_shell_p, " cat outfile.dat >> final.out \n");
                fprintf(lod2_shell_p, " mv HOMOG.DAT HOMOG.DAT.%d.%d\n", chromosome, marker);
                fprintf(lod2_shell_p, " mv homog.out homog.out.%d.%d\n",
                        chromosome, marker);
                fprintf(lod2_shell_p, " mv HOMOG.OUT HOMOG.OUT.%d.%d\n",
                        chromosome, marker);
                marker++;
            }
        }
        fprintf(lod2_shell_p, "rm -f homog.%s.sum \n", chrstr);
        fprintf(lod2_shell_p, "%s '\\\n", awk_str);
        fprintf(lod2_shell_p, "# Create a summary file from a homog.final file\\\n");

        fprintf(lod2_shell_p, "BEGIN {system(\"date\");\\\n");
        fprintf(lod2_shell_p, "i=0; \\\n");
        fprintf(lod2_shell_p, "system(\"pwd\");\\\n");
        fprintf(lod2_shell_p, "print FILENAME;\\\n");
        fprintf(lod2_shell_p,
                "printf(\"                Homogeneity         Heterogeneity\\n\");\\\n");
        fprintf(lod2_shell_p,
                "printf(\"   Locus     Lod      Theta     Lod      Theta    Alpha\\n\");\\\n");
        fprintf(lod2_shell_p, "}\\\n");
        fprintf(lod2_shell_p, "/Locus Names/ {i +=1; names[i] = $5; ok[i]=0;}\\\n");
        fprintf(lod2_shell_p, "/H2:/ {lod1[i] = $4; theta1[i] = $6; alpha[i]=$5; ok[i]++;} \\\n");
        fprintf(lod2_shell_p, "/H1:/ {lod2[i] = $4; theta2[i] = $6; ok[i]++;} \\\n");
        fprintf(lod2_shell_p, "END {\\\n");
        fprintf(lod2_shell_p, "for (j=1; j <= i; j++) \\\n");
        fprintf(lod2_shell_p, "  if (ok[j] == 2)\\\n");
        fprintf(lod2_shell_p, "   printf(\"%%8s %%10.4f %%7.2f %%10.4f %%7.2f %%7.2f\\n\", \\\n");
        fprintf(lod2_shell_p, "          names[j],lod2[j]/log(10),theta1[j], \\\n");
        fprintf(lod2_shell_p, "          lod1[j]/log(10),theta2[j], alpha[j]);\\\n");
        fprintf(lod2_shell_p, "else\\\n");
        fprintf(lod2_shell_p, "printf(\"%%8s            FAILED\\n\",names[j]);\\\n");
        fprintf(lod2_shell_p, "   }' homog.final > homog.%s.sum \n", chrstr);
        fprintf(lod2_shell_p, "echo See homog.%s.sum for a summary of the results\n", chrstr);
        fclose(lod2_shell_p);
        chmod_X_file(lod2_shell_name);
        if (nloop == 1) break;
    }
/* printf("    C-shell file  : %s\n", shellfl); */
    return;
}

/* This routine will create a reordered lod2data.## and lod2ped.## pair
   of files, along with a C-shell file to carry out the desired lod2
   analyses */

#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*LTreeTop)

void create_LOD2_format_files(linkage_ped_top **LTreeTop,
			      ped_top *PedTreeTop,
			      analysis_type *analysis, int *numchr,
			      char *file_names[], int untyped_ped_opt)
{

    linkage_locus_top *LTopp;
    int             i, single_file;
    int             xlinked;
    LTopp = LPedTreeTop->LocusTop;

    xlinked = ((LTopp->SexLinked == 1) ? 1: 0);

    if (main_chromocnt > 1) {
        draw_line();
        single_file=global_ou_files(1);
    } else {
        single_file=0;
    }
//SL
    if (LTopp->SexLinked == 2 && single_file == 1) {
        if (batchANALYSIS)
            if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Loop_Over_Chromosomes) and try again.\n");
            else
                errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the batch file (Chromosomes_Multiple) and try again.\n");
        else
            errorvf("You may not analyze autosomal chromosomes and sex chromosomes in the same file.\nPlease fix this in the chromosome reorder menu and try again.\n");

        EXIT(DATA_INCONSISTENCY);
    }

    lod2_file_names(single_file, file_names, analysis);
    LTopp=LPedTreeTop->LocusTop;

    if (single_file) {
        /* omit untyped pedigrees */
        get_loci_on_chromosome(0);
        omit_peds(untyped_ped_opt, LPedTreeTop);
        save_linkage_peds(file_names[0], LPedTreeTop, *analysis, 0, NULL);
        write_gh_locus_file(file_names[1], LPedTreeTop->LocusTop, *analysis, xlinked);
        lod2_cshell_file(file_names[3], file_names[0], file_names[1],
                         LTopp, NULL);

        create_mssg(*analysis);
        sprintf(err_msg,"     Pedigree file : %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg,"     Locus file    : %s", file_names[1]);
        mssgf(err_msg);
        sprintf(err_msg,"     C-shell file  : %s", file_names[3]);
        mssgf(err_msg);
        draw_line();
        return;
    }

    create_mssg(*analysis);
    for (i=0; i < main_chromocnt; i++) {
        *numchr=global_chromo_entries[i];
        change_output_chr(file_names[0], *numchr);
        change_output_chr(file_names[1], *numchr);
        change_output_chr(file_names[3], *numchr);

        get_loci_on_chromosome(*numchr);
        xlinked =
            (((*numchr == SEX_CHROMOSOME && LTopp->SexLinked == 2) ||
              (LTopp->SexLinked == 1))? 1 : 0);

        /* omit untyped pedigrees */
        omit_peds(untyped_ped_opt, LPedTreeTop);
        save_linkage_peds(file_names[0], LPedTreeTop, *analysis, 0, NULL);
        write_gh_locus_file(file_names[1], LPedTreeTop->LocusTop, *analysis, xlinked);
        draw_line();
        sprintf(err_msg,"    Pedigree file : %s", file_names[0]);
        mssgf(err_msg);
        sprintf(err_msg,"    Locus file    : %s", file_names[1]);
        mssgf(err_msg);
        lod2_cshell_file(file_names[3], file_names[0], file_names[1],
                         LTopp, numchr);
        sprintf(err_msg,"     C-shell file  : %s", file_names[3]);
        mssgf(err_msg);

        draw_line();

    }
    return;
}   /* End of create_lod2_format_files  */

#undef LPedTreeTop
