/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

#include "R_output.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "grow_string_ext.h"
#include "output_file_names_ext.h"
#include "utils_ext.h"
#include "plink_ext.h"
#include "user_input_ext.h"

#include "class_old.h"

/*
        batch_input_ext.h:  batchf
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
        grow_string_ext.h:  grow
  output_file_names_ext.h:  CHR_STR change_output_chr file_status
              utils_ext.h:  EXIT draw_line fcopy imin script_time_stamp
*/


char input_stem[FILENAME_LENGTH]; /* name for the generated  R_data file */
char output_stem[FILENAME_LENGTH];
char Rperl_file[FILENAME_LENGTH];
char Rdatafile[FILENAME_LENGTH];
char Rplot_function[FILENAME_LENGTH];
char Rshell_file[FILENAME_LENGTH];
char Rscript_file[FILENAME_LENGTH];
char Rmap_file[FILENAME_LENGTH];

#ifdef TEST
#define BPPI_p(x) 1
#else
#define BPPI_p(x) (x)
#endif


/* One or more models go into one plot.
   We are producing a standard score file
   for allegro and Simwalk2 results
   Format of RSW2
*/

/* always combine stats if there are more than one */

static void remove_old_file_command(FILE *fp, const char *path_name,
                                    const char *file_name);

static void convert_SW2R(R_plot_params_type Rplot_params,
			 int num_stats, int *output_stats,
			 int numchr, char *input_stat_names[],
			 int num_stat_names,
			 linkage_ped_top *Top);
/* static void SW2R_perl_script(); */
static void draw_graphs(R_plot_params_type Rplot_params, FILE *shfp,
			int numchr, analysis_type analysis,
			linkage_ped_top *Top);
static void stat_selection_menu(char *stat_names[], int num_stats,
				int  *output_stats, int *num_select);
static void init_plot_params(R_plot_params_type *Rplot_params);
static void R_plot_params(R_plot_params_type *Rplot_params, int numchr);

static void convert_allegro2R(R_plot_params_type Rplot_params,
			      int num_stats, int *output_stats,
			      int numchr, char *input_stat_names[],
			      int num_stat_names,
			      linkage_ped_top *Top);
static void convert_merlin2R(R_plot_params_type Rplot_params,
			     int num_stats, int *output_stats,
			     int numchr, char *input_stat_names[],
			     int num_stat_names,
			     linkage_ped_top *Top);
static void nplplot_init_header(R_plot_params_type Rplot_params,
				char *Rdatafile, int trait_num, int numchr);



int allegro_R_setup( int numchr, linkage_ped_top *Top);
int sw2_R_setup(analysis_type analysis, int numchr, linkage_ped_top *Top);
int merlin_R_setup(int numchr, linkage_ped_top *Top,
                   merlin_opt_type merlin_opt);


/*=======================================================*/

/*static void */

static void remove_old_file_command(FILE *fp, const char *path_name,
                                    const char *file_name)

{
    if (strcmp(path_name, ".")) {
        fprintf(fp, "if (( -e %s/%s )) then\n", path_name, file_name);
        fprintf(fp, "   /bin/rm -f %s/%s\n", path_name, file_name);
    } else {
        fprintf(fp, "if (( -e %s )) then\n", file_name);
        fprintf(fp, "   /bin/rm -f %s\n", file_name);
    }
    fprintf(fp, "endif\n");
    return;

}

static void convert_SW2R(R_plot_params_type Rplot_params,
			 int num_stats, int *output_stats,
			 int numchr, char *input_stat_names[],
			 int num_stat_names,
			 linkage_ped_top *Top)
{

    /* This function sets up the call to execute the perl script
       set up by SW2R_perl_script which combines a number
       of SW2 output files into a single R-readable file,
       if traits are to be combined, or into one file per trait

       simwalk2 has one file for all models called STATS-<chr>.ALL.
       if pcomb_stats is set, then put all stats from one file
       into one output file, otherwise, they have to be broken up.
       if pcomb_traits is set then combine traits into one file
    */


    /* output_paths[i] contain the trait specific files */

    int tr, trr, nloop, num_affec=num_traits;
    int st, ch;
    char swoutfile[FILENAME_LENGTH], chr_str[3], *syscmd;
    FILE *shfp;
    char rshfl[20+FILENAME_LENGTH];

    sprintf(rshfl, "%s/%s", output_paths[0], Rshell_file);
    shfp = fopen(rshfl, "w");
    if (shfp == NULL) {
        sprintf(err_msg, "Could not open shell script %s for writing.\n",
    		rshfl);
        mssgf(err_msg);
        EXIT(FILE_WRITE_ERROR);
    }
    fprintf(shfp, "#!/bin/csh -f\n");
    /* print some identification */
    fprintf(shfp, "# R-plot C-shell file name: %s\n", 
#ifdef HIDEPATH
            NOPATH
#else
            Rshell_file
#endif
        );
    script_time_stamp(shfp);
    fprintf(shfp, "#\n");
    /* print the statistic names */
    fprintf(shfp, "# Statistic names and numbers \n");
    for (st = 0; st < num_stat_names; st++) {
        fprintf(shfp, "#  %d) %s\n", st+1, input_stat_names[st]);
    }

    NLOOP;
    for(ch=0; ch < main_chromocnt; ch++) {
        if (main_chromocnt > 1) {
            numchr = global_chromo_entries[ch];
        }
        CHR_STR(numchr, chr_str);

        SW2OUT(chr_str, swoutfile);
        change_output_chr(Rdatafile, numchr);

        if (Rplot_params.pcomb_traits) {
            remove_old_file_command(shfp, ".", Rdatafile);
            fprintf(shfp, "set message = `%s -c%d", Rperl_file, output_stats[0]+1);
            for (st=1; st < num_stats; st++) {
                fprintf(shfp, ",%d", output_stats[st]+1);
            }
            trr=0;
            for (tr=0; tr < num_traits; tr++) {
                if (global_trait_entries[tr] < 0) continue;
                if (trr == 0) {
                    fprintf(shfp, " -t%s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[0]].TraitName);
                } else {
                    fprintf(shfp, ",%s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                }
                trr++;
            }

            fprintf(shfp, " -o%s", Rdatafile);
            for (tr=0; tr <= nloop; tr++) {
                if (nloop > 1 && tr == 0) continue;
                fprintf(shfp, " %s/%s",	trait_paths[tr], swoutfile);
                if (nloop == 1) break;
            }
            fprintf(shfp, "`\n");
            /*      fprintf(shfp, "echo\n"); */
            fprintf(shfp, "foreach m ($message)\n");
            fprintf(shfp, "    if (( $m == ERROR )) then\n");
            fprintf(shfp, "       echo \"Failed to convert SimWalk2 output file.\"\n");
            fprintf(shfp, "       exit(-1)\n");
            fprintf(shfp, "    endif\n");
            fprintf(shfp, "end\n");
            fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
            fprintf(shfp, "   exit(-1)\n");
            fprintf(shfp, "endif\n");

        } else {
            if (num_traits > 1) {
                for (tr=1; tr <= nloop; tr++) {
                    if (global_trait_entries[tr-1] < 0) continue;

                    remove_old_file_command(shfp, trait_paths[tr], Rdatafile);
                    fprintf(shfp, "set message = `%s -c%d", Rperl_file,
                            output_stats[0]+1);
                    for (st=1; st < num_stats; st++) {
                        fprintf(shfp, ",%d",
                                output_stats[st]+1);
                    }
                    fprintf(shfp, " -t%s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);

                    fprintf(shfp, " -o%s/%s", trait_paths[tr], Rdatafile);
                    fprintf(shfp, " %s/%s", trait_paths[tr], swoutfile);
                    fprintf(shfp, "`\n");
                    /*	  fprintf(shfp, "echo \n ");*/
                    fprintf(shfp, "if (( ! -f %s/%s )) then\n", trait_paths[tr],
                            Rdatafile);
                    fprintf(shfp, "   exit(-1)\n");
                    fprintf(shfp, "endif\n");
                }
            } else {
                /* sw2 output file is in current directory */
                remove_old_file_command(shfp, ".", Rdatafile);
                fprintf(shfp, "set message = `%s -c%d", Rperl_file,
                        output_stats[0]+1);
                for (st=1; st < num_stats; st++) {
                    fprintf(shfp, ",%d",
                            output_stats[st]+1);
                }
                fprintf(shfp, " -t%s",
                        Rplot_params.LocusTop->Pheno[global_trait_entries[0]].TraitName);
                fprintf(shfp, " -o%s", Rdatafile);
                fprintf(shfp, " %s", swoutfile);
                fprintf(shfp, "\n");
                /*	fprintf(shfp, "echo\n"); */
                fprintf(shfp, "foreach m ($message)\n");
                fprintf(shfp, "    if (( $m == ERROR )) then\n");
                fprintf(shfp, "       echo \"Failed to convert SimWalk2 output file.\"\n");
                fprintf(shfp, "       exit(-1)\n");
                fprintf(shfp, "    endif\n");
                fprintf(shfp, "end\n");
                fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
                fprintf(shfp, "   exit(-1)\n");
                fprintf(shfp, "endif\n");
            }
        }

    }

    draw_graphs(Rplot_params, shfp, numchr, NONPARAMETRIC, Top);
    fclose(shfp);
    syscmd=CALLOC((strlen(rshfl)+11), char);
    sprintf(syscmd, "chmod +x %s", rshfl);
    System(syscmd);
    free(syscmd);
}

static void ps_file_failure(FILE *shfp,
			    char *outputfile,
			    char *Rerror)
{
    fprintf(shfp,
            "   echo \"Unable to create plot file %s\"\n",
            outputfile);
    fprintf(shfp,
            "   echo Check file %s for error messages\n",
            Rerror);
}

static void draw_graphs(R_plot_params_type Rplot_params,
			FILE *shfp, int numchr, analysis_type analysis,
			linkage_ped_top *Top)
{
    /*
       write command to invoke R on multiple R-readable files
       which will be combined into one postscript file or
       displayes onscreen.
       Command is written into the shell file.

       Ymax is the ylimit of each graph,
       threshold is a lod score where a horizontal line will
       be drawn on each graph.
       If traits are to be drawn separately, then the R-readable
       files are found in each trait subdirectory, and abelled with
       the chromosome number.

       graph_out_file_stem.<stat_name> are files to look for
       in the current directory and its subdirectories
       combine_chr_or_trait =
       'C' -> all chromosomes for each trait
       in one file.
       'T' -> all traits for each chromosome
       in one file.

       Input files will be labeled RSW2.<chr>.<tr>


    */
    int trp, np;
    int tr, chr;
    char output_file[2*FILENAME_LENGTH], rscr[20+FILENAME_LENGTH];
    char Rerror_file[FILENAME_LENGTH+3];
    FILE *Rp;
    char legend_str[FILENAME_LENGTH], rdataf[2*FILENAME_LENGTH+1];
    char mega2_map[20+FILENAME_LENGTH];
    
    if (base_pair_position_index >= 0) {
        sprintf(mega2_map, "%s/%s", output_paths[0], Rmap_file);
        // always generate the map file from internal data...
        // annotated, genetic distance in centiMorgans, no comments (R doesn't deal with them), use a 3-column map...
        write_R_PLINK_map_file(Top, mega2_map);
    } else if (! (BPPI_p(base_pair_position_index >= 0)) ) {
        mssgvf("UCSC Custom Tracks could be produced if you provided a physical map.\n");
    }

    if (analysis == NONPARAMETRIC || analysis == TO_Allegro) {
        strcpy(Rplot_params.YLabel, "-log(P-value)");
    } else {
        strcpy(Rplot_params.YLabel, "LOD Score");
    }
    sprintf(rscr, "%s/%s", output_paths[0], Rscript_file);
    Rp=fopen(rscr, "w");
    if (Rp == NULL) {
        sprintf(err_msg, "Could not open shell script %s for writing.\n",
     		rscr);
        mssgf(err_msg);
        EXIT(FILE_WRITE_ERROR);
    }
    /* commands to save the plots */
    sprintf(Rerror_file, "%sout", Rscript_file);
    fprintf(shfp, "R CMD BATCH %s\n", Rscript_file);

    script_time_stamp(Rp);
    fprintf(Rp, "library(nplplot)\n");

    /* Number of output files can be 1 or equal to the number of
       converted Rdata files
    */

    if (Rplot_params.pcomb_traits) {
        /* can choose to combine chromosomes or not */
        if (Rplot_params.fcomb_chromo || main_chromocnt <= 1) {
            /* simplest case, all in one file,
               output file name is set correctly
            */
            if (main_chromocnt > 1) {
                change_output_chr(Rdatafile, global_chromo_entries[0]);
            } else {
                change_output_chr(Rdatafile, numchr);
            }

            fprintf(Rp, "if (%s(c(\"%s\"", Rplot_function, Rdatafile);
            for (chr=1; chr < main_chromocnt; chr++) {
                change_output_chr(Rdatafile, global_chromo_entries[chr]);
                fprintf(Rp, ", \"%s\"", Rdatafile);
            }
            /* Wrong -  Don't need to change the output file */
            if (main_chromocnt > 1) {
                change_output_chr(Rplot_params.Routputfile, 0);
            }
            fprintf(Rp, "), col=%d, row=%d, mode=\"%c\", output=\"%s\", ",
                    Rplot_params.ncols, Rplot_params.nrows,
                    Rplot_params.mode, Rplot_params.Routputfile);

            if (main_chromocnt == 1) {
                nplplot_init_header(Rplot_params, Rdatafile, 0, numchr);
                fprintf(Rp, "headerfiles=\"%s.hdr\"", Rdatafile);
            } else {
                fprintf(Rp, "headerfiles=c(");
                for(chr=0; chr < main_chromocnt; chr++) {
                    change_output_chr(Rdatafile, global_chromo_entries[chr]);
                    nplplot_init_header(Rplot_params, Rdatafile, 0, global_chromo_entries[chr]);
                    if (chr == (main_chromocnt-1)) {
                        fprintf(Rp, "\"%s.hdr\"", Rdatafile);
                    } else {
                        fprintf(Rp, "\"%s.hdr\", ", Rdatafile);
                    }
                }
                fprintf(Rp, ")");
            }
            /* Single call to nplplot.multi */
	    if (BPPI_p(base_pair_position_index >= 0)) {
                // Needed to produce UCSC Custom Tracks
                fprintf(Rp, ", lgnd=\"page\", customtracks=TRUE, mega2mapfile=\"%s\")==F) {\n",
                        Rmap_file);
            } else {
                fprintf(Rp, ", lgnd=\"page\")==F) {\n");
            }
            fprintf(Rp, "   system(\"/bin/rm %s\")\n}\n",
                    Rplot_params.Routputfile);
            fprintf(shfp,"if ( -e %s ) then\n",
                    Rplot_params.Routputfile);
            fprintf(shfp, "  echo \"Created plot file %s\" \n",
                    Rplot_params.Routputfile);
            fprintf(shfp, "else \n");
            ps_file_failure(shfp, Rplot_params.Routputfile,
                            Rerror_file);
            fprintf(shfp,
                    "  echo \"Unable to create plot file %s\"\n",
                    Rplot_params.Routputfile);

            fprintf(shfp, "endif\n");
        } else {
            /* one file per chromosome */
            for (chr=0; chr < main_chromocnt; chr++) {
                if (main_chromocnt == 1) {
                    change_output_chr(Rdatafile, numchr);
                    change_output_chr(Rplot_params.Routputfile,
                                      numchr);
                } else {
                    change_output_chr(Rdatafile, global_chromo_entries[chr]);
                    change_output_chr(Rplot_params.Routputfile,
                                      global_chromo_entries[chr]);
                }
                /* 	strcat(Rplot_params.Routputfile, ".ps"); */
                /* Not necessary since we put in the extension
                   in the menu */
                fprintf(Rp, "if (");
                fprintf(Rp, "%s(\"%s\",", Rplot_function, Rdatafile);
                fprintf(Rp, "col=%d, row=%d, mode=\"%c\", output=\"%s\", ",
                        Rplot_params.ncols, Rplot_params.nrows,
                        Rplot_params.mode, Rplot_params.Routputfile);

                nplplot_init_header(Rplot_params, Rdatafile, 0, global_chromo_entries[chr]);
                /* multiple calls to nplplot.multi, one per chromosome */
	        if (chr == 0 && BPPI_p(base_pair_position_index >= 0)) {
                    // Needed to produce UCSC Custom Tracks
                    fprintf(Rp,
                            "headerfiles=\"%s.hdr\", lgnd=\"page\", customtracks=TRUE, mega2mapfile=\"%s\")==F) {\n",
                            Rdatafile, Rmap_file);
                } else {
                    fprintf(Rp, "headerfiles=\"%s.hdr\", lgnd=\"page\")==F) {\n", Rdatafile);
                }
                fprintf(Rp, "  system(\"/bin/rm %s\")\n}\n",
                        Rplot_params.Routputfile);
                fprintf(shfp,"if ( -e %s ) then\n",
                        Rplot_params.Routputfile);
                fprintf(shfp, "  echo \"Created plot file %s\" \n",
                        Rplot_params.Routputfile);
                fprintf(shfp, "else \n");
                ps_file_failure(shfp, Rplot_params.Routputfile,
                                Rerror_file);
                fprintf(shfp, "endif\n");
            }
        }
    } else {
        np=0;
        if (Rplot_params.fcomb_chromo && Rplot_params.fcomb_traits) {
            /* all plots in one file which should have the correct name
               example call to gplot -> to write down in
               gplot(c("tr1/RSW.01", "tr1/RSW.02", "tr2/RSW.01", "tr2/RSW.02"),
	       output = outfile, yline=thresh, ymin=Ymin, ymax=Ymax)
            */
            /*      create the legend string, loop between 0 and num_traits-1 */
            strcpy(legend_str, "1");
            for (tr=1; tr < num_traits; tr++) {
                /* added skipping, looks like it is used below */
                if (global_trait_entries[tr-1] < 0) continue;
                /* sprintf(legend_str, "%s, %d ", legend_str, tr*main_chromocnt+1); */
                grow(legend_str, ", %d ",  tr*main_chromocnt+1);
            }
            /*      list of input file names, loop between 1 and num_traits */
            fprintf(Rp, "if (");
            fprintf(Rp, "%s(c(", Rplot_function);
            for(tr=1; tr <= num_traits; tr++) {
                if (global_trait_entries[tr-1] < 0) continue;
                if (analysis == TO_MERLINONLY && LoopOverTrait == 0 && num_traits > 2) {
                    /* Trait-specific Rdata files are in current folder */
                    change_output_chr(Rdatafile, -9);
                    sprintf(rdataf, "%s_%s", Rdatafile,
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                } else if (LoopOverTrait == 1) {
                    sprintf(rdataf, "%s/%s", trait_paths[tr], Rdatafile);
                } else {
                    strcpy(rdataf, Rdatafile);
                }

                if (main_chromocnt > 1) {
                    change_output_chr(rdataf, global_chromo_entries[0]);
                    change_output_chr(Rplot_params.Routputfile, 0);
                } else {
                    change_output_chr(rdataf, numchr);
                }
                if (tr > 1) {
                    /* not the first file */
                    fprintf(Rp, ", \"%s\"", rdataf);
                } else {
                    /* The first file (first trait and first chromosome */
                    fprintf(Rp, "\"%s\"", rdataf);
                }
                for (chr=1; chr < main_chromocnt; chr++) {
                    if (main_chromocnt > 1) {
                        change_output_chr(rdataf, global_chromo_entries[chr]);
                    }
                    fprintf(Rp, ", \"%s\"", rdataf);
                }
            }
            fprintf(Rp, "), col=%d, row=%d, mode=\"%c\", output=\"%s\", ",
                    Rplot_params.ncols, Rplot_params.nrows,
                    Rplot_params.mode, Rplot_params.Routputfile);

            fprintf(Rp, "headerfiles=c(");

            np=0;
            for(tr=1; tr <= num_traits; tr++) {
                if (global_trait_entries[tr-1] < 0) continue;
                if (analysis == TO_MERLINONLY && LoopOverTrait == 0 && num_traits > 2) {
                    /* Trait-specific Rdata files are in current folder */
                    change_output_chr(Rdatafile, -9);
                    sprintf(rdataf, "%s_%s", Rdatafile,
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                } else if (LoopOverTrait == 1) {
                    sprintf(rdataf,"%s/%s", trait_paths[tr], Rdatafile);
                } else {
                    strcpy(rdataf, Rdatafile);
                }
                for (chr=0; chr < main_chromocnt; chr++) {
                    if (main_chromocnt > 1) {
                        change_output_chr(rdataf, global_chromo_entries[chr]);
                    } else {
                        change_output_chr(rdataf, numchr);
                    }

                    nplplot_init_header(Rplot_params, rdataf, 0, global_chromo_entries[chr]);
                    if (np > 0) {
                        fprintf(Rp, ", \"%s.hdr\"", rdataf);
                    } else {
                        fprintf(Rp, "\"%s.hdr\"", rdataf);
                        np++;
                    }
                }
            }
            /* single call over all chromosomes */
	    if (BPPI_p(base_pair_position_index >= 0)) {
                // Needed to produce UCSC Custom Tracks
                fprintf(Rp,
                        "), lgnd=c(%s), customtracks=TRUE, mega2mapfile=\"/%s\") == F) {\n",
                        legend_str, Rmap_file);
            } else {
                fprintf(Rp, "), lgnd=c(%s)) == F) {\n", legend_str);
            }
            fprintf(Rp, "  system(\"/bin/rm %s\")\n}\n",
                    Rplot_params.Routputfile);
            fprintf(shfp, "if (-e %s) then\n",
                    Rplot_params.Routputfile);
            fprintf(shfp, "   echo \"Created plot file %s\" \n",
                    Rplot_params.Routputfile);
            fprintf(shfp, "else \n");
            ps_file_failure(shfp, Rplot_params.Routputfile,
                            Rerror_file);
            fprintf(shfp, "endif\n");
        }

        else if (!Rplot_params.fcomb_chromo && Rplot_params.fcomb_traits) {
            np=0;
            /* one gplot call per chromosome */
            for (chr=0; chr < main_chromocnt; chr++) {
                /* 	strcat(Rplot_params.Routputfile, ".ps"); */
                np=0;
                fprintf(Rp, "if (");
                for(tr=1; tr <= num_traits; tr++) {

                    if (global_trait_entries[tr-1] < 0) continue;

                    if (analysis == TO_MERLINONLY && LoopOverTrait == 0 && num_traits > 2) {
                        change_output_chr(Rdatafile, -9);
                        sprintf(rdataf, "%s_%s", Rdatafile,
                                Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                    } else if (LoopOverTrait == 1) {
                        sprintf(rdataf, "%s/%s", trait_paths[tr], Rdatafile);
                    } else {
                        strcpy(rdataf, Rdatafile);
                    }

                    if (main_chromocnt > 1) {
                        change_output_chr(rdataf, global_chromo_entries[chr]);
                        change_output_chr(Rplot_params.Routputfile, global_chromo_entries[chr]);
                    } else {
                        change_output_chr(rdataf, numchr);
                        change_output_chr(Rplot_params.Routputfile, numchr);
                    }

                    if (np == 0) {
                        fprintf(Rp, "%s(c(\"%s\"", Rplot_function, rdataf);
                        np++;
                    } else {
                        fprintf(Rp, ", \"%s\"", rdataf);

                    }
                }
                fprintf(Rp, "), col=%d, row=%d, mode=\"%c\", output=\"%s\", ",
                        Rplot_params.ncols, Rplot_params.nrows,
                        Rplot_params.mode, Rplot_params.Routputfile);

                fprintf(Rp, "headerfiles=c(");
                np=0;
                for(tr=1; tr <= num_traits; tr++) {
                    if (global_trait_entries[tr-1] < 0) continue;

                    if (analysis == TO_MERLINONLY && LoopOverTrait == 0 && num_traits > 2) {
                        /* Trait-specific Rdata files are in current folder */
                        change_output_chr(Rdatafile, -9);
                        sprintf(rdataf, "%s_%s", Rdatafile,
                                Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                    } else if (LoopOverTrait == 1) {
                        sprintf(rdataf, "%s/%s", trait_paths[tr], Rdatafile);
                    } else {
                        strcpy(rdataf, Rdatafile);
                    }

                    if (main_chromocnt > 1) {
                        change_output_chr(rdataf, global_chromo_entries[chr]);
                    } else {
                        change_output_chr(rdataf, numchr);
                    }

                    /* all titles refer to a single chromosome */
                    nplplot_init_header(Rplot_params, rdataf, 0, global_chromo_entries[chr]);
                    if (np==0) {
                        fprintf(Rp, "\"%s.hdr\"", rdataf);
                        np++;
                    } else {
                        fprintf(Rp, ", \"%s.hdr\"", rdataf);
                    }
                }
                if (chr == 0 && BPPI_p(base_pair_position_index >= 0)) {
                    // Needed to produce UCSC Custom Tracks
                    fprintf(Rp, "), lgnd=\"page\", customtracks=TRUE, mega2mapfile=\"%s\") == F) {\n",
                            Rmap_file);
                } else {
                    fprintf(Rp, "), lgnd=\"page\") == F) {\n");
                }
                fprintf(Rp, "  system(\"/bin/rm %s\")\n}\n",
                        Rplot_params.Routputfile);
                fprintf(shfp, "if (-e %s) then\n", Rplot_params.Routputfile);
                fprintf(shfp, "   echo \"Created plot file %s\" \n",
                        Rplot_params.Routputfile);
                fprintf(shfp, "else \n");
                ps_file_failure(shfp, Rplot_params.Routputfile,
                                Rerror_file);
                fprintf(shfp, "endif\n");
            }
        }


        else if (Rplot_params.fcomb_chromo && !Rplot_params.fcomb_traits) {
            /* one or more traits, for multiple traits,
               one gplot call per trait */
            for(tr=0; tr <= num_traits; tr++) {
                if (num_traits > 1 && tr==0) continue;

                if (tr > 0) {
                    trp = tr - 1;
                } else {
                    trp = tr;
                }

                if (global_trait_entries[trp] < 0) continue;

                if (main_chromocnt > 1) {
                    change_output_chr(Rdatafile, global_chromo_entries[0]);
                    /*
                      change_output_chr(Rplot_params.Routputfile,0);
                    */
                } else {
                    change_output_chr(Rdatafile, numchr);
                }
                fprintf(Rp, "if (");
                fprintf(Rp, "%s(c(\"%s/%s\"", Rplot_function,
                        trait_paths[tr], Rdatafile);
                for (chr=1; chr < main_chromocnt; chr++) {
                    change_output_chr(Rdatafile, global_chromo_entries[chr]);
                    /* don't need to change output file name */
                    fprintf(Rp, ", \"%s/%s\"", trait_paths[tr], Rdatafile);
                }

                fprintf(Rp, "), col=%d, row=%d, mode=\"%c\", output=\"%s/%s\", ",
                        Rplot_params.ncols, Rplot_params.nrows,
                        Rplot_params.mode,
                        trait_paths[tr], Rplot_params.Routputfile);

                fprintf(Rp, "headerfiles=c(");
                for (chr=0; chr < main_chromocnt; chr++) {
                    change_output_chr(Rdatafile, global_chromo_entries[chr]);
                    nplplot_init_header(Rplot_params, Rdatafile, tr, global_chromo_entries[chr]);
                    if (chr < (main_chromocnt-1)) {
                        fprintf(Rp, "\"%s/%s.hdr\",", trait_paths[tr], Rdatafile);
                    } else {
                        fprintf(Rp, "\"%s/%s.hdr\")", trait_paths[tr], Rdatafile);

                    }
                }
                /* single file for all chromosomes */
                sprintf(output_file, "%s/%s", trait_paths[tr], Rplot_params.Routputfile);
		if (BPPI_p(base_pair_position_index >= 0 )) {
                    // Needed to produce UCSC Custom Tracks
                    fprintf(Rp, ", lgnd=\"page\", customtracks=TRUE, mega2mapfile=\"%s\") == F) {\n",
                            Rmap_file);
                } else {
                    fprintf(Rp, ", lgnd=\"page\") == F) {\n");
                }
                fprintf(Rp, "  system(\"/bin/rm %s\")\n}\n",
                        output_file);
                fprintf(shfp, "if (-e %s) then\n", output_file);
                fprintf(shfp, "   echo \"Created plot file %s\" \n",
                        output_file);

                fprintf(shfp, "else \n");
                ps_file_failure(shfp, output_file, Rerror_file);
                fprintf(shfp, "endif\n");

                if (num_traits == 1) break;
            }
        } else {
            /* separate calls to gplot per chromosome and trait */
            for(tr=0; tr <= num_traits; tr++) {
                if (num_traits > 1 && tr==0) continue;
                if (tr > 0) {
                    trp = tr - 1;
                } else {
                    trp = tr;
                }

                if (global_trait_entries[trp] < 0) continue;

                for (chr=0; chr < main_chromocnt; chr++) {
                    if (main_chromocnt > 1) {
                        change_output_chr(Rplot_params.Routputfile,
                                          global_chromo_entries[chr]);
                        change_output_chr(Rdatafile, global_chromo_entries[chr]);

                    } else {
                        change_output_chr(Rplot_params.Routputfile, numchr);
                        change_output_chr(Rdatafile, numchr);
                    }

                    sprintf(output_file, "%s/%s", trait_paths[tr], Rplot_params.Routputfile);
                    fprintf(Rp, "if (");
                    fprintf(Rp,
                            "%s(\"%s/%s\", col=%d, row=%d, mode=\"%c\", output=\"%s\", ",
                            Rplot_function, trait_paths[tr],
                            Rdatafile, Rplot_params.ncols,
                            Rplot_params.nrows,
                            Rplot_params.mode, output_file);

                    /* separate calls for each chromosome */
                    nplplot_init_header(Rplot_params, Rdatafile, tr, global_chromo_entries[chr]);
                    fprintf(Rp, "headerfiles=\"%s/%s.hdr\", ", trait_paths[tr], Rdatafile);
                    if (chr == 0 && BPPI_p(base_pair_position_index >= 0)) {
                        // Needed to produce UCSC Custom Tracks
                        fprintf(Rp, "lgnd=\"page\", customtracks=TRUE, mega2mapfile=\"%s\")==F) {\n",
                                Rmap_file);
                    } else {
                        fprintf(Rp, "lgnd=\"page\")==F) {\n");
                    }
                    fprintf(Rp, "  system(\"/bin/rm %s\")\n}\n",
                            output_file);

                    fprintf(shfp, "if (-e %s) then\n",  output_file);
                    fprintf(shfp, "   echo \"Created plot file %s\" \n",
                            output_file);
                    fprintf(shfp, "else \n");
                    ps_file_failure(shfp, output_file, Rerror_file);
                    fprintf(shfp, "endif\n");

                }
                if (num_traits==1) break;
            }
        }
    }
    fclose(Rp);
    fprintf(shfp, "date >> %s\n", Rerror_file);
    fprintf(shfp, "echo Transcript of R session in %s.\n",
            Rerror_file);

    return;
}


static void stat_selection_menu(char *input_stat_names[],
				int num_stats,
				int *output_stats, int *num_select)

{

    int i, j, done;
    char str[FILENAME_LENGTH], *str_p, selected[3]; /* 1-999 stats */
    char num[6];
    int batchRPLOT=0;

    if (Mega2BatchItems[/* 35 */ Rplot_Statistics].items_read) {
        strcpy(str, Mega2BatchItems[/* 35 */ Rplot_Statistics].value.name);
        batchRPLOT=1;
        /* parse the string and set the stat values */
        *num_select=0;
        strcpy(num, strtok(str, " "));
        while(strcmp(num, "e")) {
            i=atoi(num);
            if (i > 0 && i <= num_stats) {
                output_stats[*num_select] = i;
                (*num_select)++;
            } else {
                batchRPLOT=0;
                errorf("Invalid statistic number in batch file.");
                errorf("   Use the menu to select statistics.");
                break;
            }
            strcpy(num, strtok(NULL, " "));
        }
    }
    if (batchRPLOT) {
        mssgf("Statistic selections read in from batch file.");
    } else {
        if (num_stats == 1) {
            *num_select=1;
            output_stats[0]=1;
        } else {
            /* go through the selection process */
            strcpy(str, "");
            draw_line();

            printf("R plot statistic selection menu:\n");
            newline;
            printf("NPL statistics will be automatically plotted into a \n");
            printf("postscript file using R after %s \n",
                   ProgName);
            printf("has been run.\n");
            printf("Please select the statistics to be included in this R plot.\n");
            newline;
            printf("This list can be later modified in the shell script %s before\n",
                   Rshell_file);
            printf("running this script.\n");
            newline;

            draw_line();
            for(i=0; i < num_stats; i++) {
                printf(" %d) %s\n", i+1, input_stat_names[i]);
            }
            draw_line();

            *num_select=0; done=0;
            int eof;
            while(!done) {
                printf("Enter string of statistic numbers ('e' to terminate) > ");
                fflush(stdout);
                eof = fgets(str, 199, stdin) == NULL; newline;
                str_p = &(str[0]);
                while(*str_p != '\0' && eof != 1) {
                    if (*num_select >= num_stats) {
                        done = 1; break;
                    }
                    if (*str_p == 'e' || *str_p == 'E') {
                        str_p++;
                        done = 1; break;
                    } else if (*str_p == '\n') {
                        break;
                    } else if (isspace((int) *str_p)) str_p++;
                    else if (isdigit((int) *str_p)) {
                        j=0;
                        strcpy(selected, "");
                        while(isdigit((int) *str_p)) {
                            selected[j] = *str_p;
                            str_p++; j++;
/*	      if (j > ((int)(log10(num_stats))+1))  {
              printf("Statistic number too big!\n");
              *num_select=0;
              fflush(stdin);
	      } */
                        }
                        selected[j] = '\0';
                        output_stats[*num_select] = atoi(selected); (*num_select)++;
                    } else {
                        printf("Invalid characters in input (%c).\n", *str_p);
                        *num_select=0; fflush(stdin);
                        break;
                    }
                }
                if (*num_select > 0) {
                    for(j=0; j<*num_select; j++) {
                        if (output_stats[j] < 1 || output_stats[j] > num_stats) {
                            printf("Invalid statistic number %d\n", output_stats[j]);
                            printf("Please select again.\n");
                            *num_select=0;
                            done=0;
                            strcpy(str, "");
                        }
                    }
                }
                if (eof) done = 1;
            }
        }
    }

    if (*num_select > 0) {
        strcpy(err_msg, "R output will include the following statistics:");
        if (InputMode == INTERACTIVE_INPUTMODE) {
            strcpy(Mega2BatchItems[/* 35 */ Rplot_Statistics].value.name, "");
        }
        for(j=0; j < *num_select; j++) {
            char stat_name[20];
            grow(Mega2BatchItems[/* 35 */ Rplot_Statistics].value.name, " %d", output_stats[j]);
            output_stats[j]--;
            sscanf(input_stat_names[output_stats[j]], "%s", stat_name);
            strcat(err_msg, " "); strcat(err_msg, stat_name);
        }
        mssgf(err_msg);
        draw_line();
        if (InputMode == INTERACTIVE_INPUTMODE) {
            strcat(Mega2BatchItems[/* 35 */ Rplot_Statistics].value.name, " e");
            batchf(/* 35 */ Rplot_Statistics);
        }
    } else {
        warnf("No statistics were selected, not setting up plots.");
    }
    return;
}

/* put all traits and/or all stats into one plot?
   This is a toggle menu
   pcomb_traits automatically override fcomb_traits */

static void init_plot_params(R_plot_params_type *Rplot_params)

{

    Rplot_params->pcomb_traits=((num_traits <= 1 || LoopOverTrait == 0)?
                                1 : 0);
    /* default = separate traits if not LoopOverTraits */
    Rplot_params->fcomb_traits=Rplot_params->pcomb_traits; /* separate if possible */
    Rplot_params->fcomb_chromo=1; /* combine chromosomes */
    Rplot_params->Yrange[0] = 0.0;
    Rplot_params->Yrange[1] = 3.0;
    Rplot_params->thresh = 2.0;
    Rplot_params->ncols=2;
    Rplot_params->nrows=2;
    Rplot_params->mode='l';
    Rplot_params->yfixed=0;
    /* Keep these fixed for now, will add customization later on */
    strcpy(Rplot_params->YLabel, "");
    Rplot_params->legend_cex = 0.7;
    Rplot_params->axis_cex = 0.9;
    Rplot_params->tcl = -0.3;


    return;
}


static void R_plot_params(R_plot_params_type *Rplot_params, int numchr)


/*  int *pcomb_traits, int *fcomb_traits,  */
/*  		   int *fcomb_chromo,  */
/*  		   float *thresh, float *YRange, */
/*  		   char *Routputfile, int ncols, int nrows)  */

{

    int done=0;
    /*
       pctr <- pcomb_trait
       fctr <- fcomb_trait
       fcchr <- fcomb_chromo
       nc <- num_col
       nr <- num_row
       thr <- threshold
       mo <- orientation
    */
    int item, pctr=0, fctr=0, fcchr=0;
    int nc=0, nr=0, yf=0, ymax, ymin, thr, mo, fname;
    char *pos;
    char cselect[10], opfile_type[15];
    char toggle[25], *tog_ptr, change[25], *chg_ptr;
    char stem[5], fl_stat[12];
    int old_loopovertrait=LoopOverTrait;
    int tog_initiated, plots_per_file;
    char yval[30];

    /* default values dictate separate plots per trait and
       separate files per chromosome and combined plot for statistics,
       so multiple traits mean multiple plots per file at this time.
    */

    if (main_chromocnt > 1) {
        if (!Rplot_params->fcomb_chromo) {
            if (!DEFAULT_ROPTS) {
                change_output_chr(Rplot_params->Routputfile, -9);
            }
        } else {
            change_output_chr(Rplot_params->Routputfile, 0);
        }
    } else {
        change_output_chr(Rplot_params->Routputfile, numchr);
    }

    if (DEFAULT_ROPTS) {
        mssgf("Graphical output options set to default values");
        mssgf("as requested in the batch file.");
        return;
    }

    printf("Customization of graphical output:\n");
    newline;

    if ((old_loopovertrait == 0 && num_traits > 2) ||
       (old_loopovertrait == 1 && num_traits > 1)) {

        printf("** Number of plots displayed on a page can be controlled using <row> and <col>.\n");

        if (old_loopovertrait == 0) {
            printf("** In combined-traits mode, traits cannot be plotted into separate files.\n");
        } else {
            printf("** Traits may be combined into a single file\n");
            printf("   or plotted separately, one trait per plot file.\n");
        }
    }

    if (main_chromocnt > 1) {
        printf("** Choose one file per chromosome, or combine all into one file.\n");
    }


    printf("** Y-axis min and max values can be enforced by setting \"enforce\" to yes\n");
    printf("   otherwise, these will be set according to the plot data.\n");

    printf("** Plot orientation can be toggled between \n");
    printf("   portrait or landscape.\n");

    printf("** Format of output file (ps/pdf) is decided by file-name extension (.ps/.pdf).\n");
    printf("   If no extension is provided, a pdf file will be created.\n");

    newline;
    while(!done) {
        item = 1;
        fname = item;

        if (Rplot_params->pcomb_traits) {
            plots_per_file = 1;
        } else {
            plots_per_file = ((Rplot_params->fcomb_traits)? num_traits :  1);
        }

        plots_per_file = plots_per_file *
            ((Rplot_params->fcomb_chromo)? main_chromocnt : 1);

        /* For 1 - 3 plots per file, don't use 2 columns */
        if (plots_per_file < 4) {
            Rplot_params->ncols = 1;
            Rplot_params->nrows = imin(plots_per_file, Rplot_params->nrows);
        }

        if (main_chromocnt > 1 && !Rplot_params->fcomb_chromo) {
            strcpy(stem, "stem");
        } else {
            strcpy(stem, "");
        }

        /* Use LoopOverTrait to indicate if separate R plot files
           are created per trait or not */
        LoopOverTrait = ((Rplot_params->fcomb_traits)?
                         0 : 1);
        draw_line();
        pctr=fctr=fcchr=nc=nr=yf=0;
        tog_initiated=0;

        printf("R plot parameter selection menu:\n");
        printf("0) Done with this menu - please proceed\n");
        if ((main_chromocnt == 1 || Rplot_params->fcomb_chromo)) {
            /* See if the output file has  a .ps or a .pdf extension */
            pos = strrchr(Rplot_params->Routputfile, '.');
            if (pos == NULL) {
                /* No .ps or .pdf extension, add one */
                strcat(Rplot_params->Routputfile, ".pdf");
                strcpy(opfile_type, "PDF");
            } else if (!strcmp(pos, ".ps") || !strcmp(pos, ".pdf")) {
                if (!strcmp(pos, ".pdf")) {
                    strcpy(opfile_type, "PDF");
                } else {
                    strcpy(opfile_type, "Postscript");
                }
            } else {
                /* No .ps or .pdf extension, add one */
                strcat(Rplot_params->Routputfile, ".pdf");
                strcpy(opfile_type, "PDF");
            }
/* 	if (!strcmp(&(Rplot_params->Routputfile[pos]), ".ps") || */
            printf(" %d) %s output file name          %-15s\t%s\n",
                   item, opfile_type, Rplot_params->Routputfile,
                   file_status(Rplot_params->Routputfile, fl_stat));
        } else {
            change_output_chr(Rplot_params->Routputfile, -9);
            printf(" %d) Plot output file name stem     %-15s\n",
                   fname, Rplot_params->Routputfile);
        }
        item++;

        if ((num_traits > 1 && old_loopovertrait == 1) ||
           (num_traits > 2 && old_loopovertrait == 0))     {
            printf(" %d) Plot traits together in same graph   [%s]\n",
                   item, yorn[Rplot_params->pcomb_traits]);
            pctr=item; item++;

            if (old_loopovertrait == 0) {
                Rplot_params->fcomb_traits=1;
            } else {
                printf(" %d) Plot traits in one output file         [%s]\n",
                       item, yorn[Rplot_params->fcomb_traits]);
                fctr=item; item++;
            }
        }
        if (main_chromocnt > 1) {
            printf(" %d) Combine chromosomes into one file       [%s]\n",
                   item, yorn[Rplot_params->fcomb_chromo]);
            fcchr=item; item++;
        }

        sprintf(yval, "%-4.2f", Rplot_params->Yrange[0]);
        printf(" %d) Minimum Y-axis value                    %4s\n",
               item, yval);
        ymin=item; item++;
        sprintf(yval, "%-4.2f", Rplot_params->Yrange[1]);
        printf(" %d) Maximum Y-axis value                    %4s\n",
               item, yval);
        ymax=item; item++;

        printf(" %d) Horizontal cut-off line at              %-3.2f\n",
               item, Rplot_params->thresh);
        thr=item; item++;

        if (plots_per_file >= 2) {
            printf(" %d) Plots per page: number of rows          %d\n",
                   item, Rplot_params->nrows);
            nr=item; item++;
        }
        if (plots_per_file >= 4) {
            printf(" %d) Plots per page: number of columns       %d\n",
                   item, Rplot_params->ncols);
            nc=item;       item++;
        }
        printf(" %d) Plot orientation             [%s]\n",
               item, ((Rplot_params->mode == 'l')?
                      "landscape" : "portrait"));
        mo=item; item++;
        printf(" %d) Enforce Y-axis bounds even if data does not fit [%s]\n",
               item, yorn[Rplot_params->yfixed]);
        yf=item;

        /* toggle items == pctr, fctr, fcchr, mo, yf */
        strcpy(toggle, "");
        tog_ptr=&(toggle[0]);

        if (pctr) {
	    // need to be careful if the value is > 9
	    *tog_ptr++ = (char)(pctr+'0');
            tog_initiated=1;
        }
        if (fctr) {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ = (char)(fctr+48);
            tog_initiated=1;
        }
        if (fcchr) {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ = (char)(fcchr+48);
            tog_initiated=1;
        }
        if (mo > 9) {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ = (char)(48 + (int)(mo/10));
            *tog_ptr++ = (char)(48 + (mo % 10));
            tog_initiated=1;
        } else {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ = (char)(mo+48);
            tog_initiated=1;
        }
        if (yf > 9) {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ =  (char)(48 + (int)(yf/10));
            *tog_ptr++ = (char)(48 + (yf % 10));
        } else {
            if (tog_initiated) {
                *tog_ptr++ = ',';
            }
            *tog_ptr++ = (char)(yf+48);
        }

        *tog_ptr='\0';
        strcat(toggle, " to toggle, ");

        /* value items ymin, ymax, thr, nr, nc */
        strcpy(change, "");
        chg_ptr=&(change[0]);
        *chg_ptr++ = (char)(ymin+48); *chg_ptr++ = ',';
        *chg_ptr++ = (char)(ymax+48); *chg_ptr++ = ',';
        *chg_ptr++ = (char)(thr+48);

        if (nc) {
            *chg_ptr++ = ',';
            *chg_ptr++ = (char)(nc+48);
        }
        if (nr) {
            *chg_ptr++ = ',';
            *chg_ptr++ = (char)(nr+48);
        }
        *chg_ptr='\0';

        printf("Select from options 0-%d (%s%s to change values) > ",
               item, toggle, change);
        fflush(stdout);

        IgnoreValue(fgets(cselect, 9, stdin)); fflush(stdin);
        sscanf(cselect, "%d", &item);
        newline;
        
        // Process user input...
        if (item == 0) {
            done = 1;
        } else if (item == fname) {
            // get a new plot output file name..
            printf("If a .ps or .pdf extension is not provided\n");
            printf("a '.pdf' extension will be added to the plot file-name.\n");
            printf("Enter new plot file name %s > ", stem);
            fcmap(stdin, "%s", Rplot_params->Routputfile);
            newline;
        } else if (item == pctr) {
            // plot traits together in the same graph...
            Rplot_params->pcomb_traits = TOGGLE(Rplot_params->pcomb_traits);
            Rplot_params->fcomb_traits=((Rplot_params->pcomb_traits  == 1)? 1:
                                        Rplot_params->fcomb_traits);
        } else if (item == fctr) {
            // plot traits in out output file...
            Rplot_params->fcomb_traits = TOGGLE(Rplot_params->fcomb_traits);
            if (Rplot_params->pcomb_traits == 1 && Rplot_params->fcomb_traits == 0) {
                printf("Cannot put traits in separate files unless traits are \n");
                printf("plotted into separate plots (i.e. set item %d to no).\n", pctr);
                Rplot_params->fcomb_traits=1;
            }
        } else if (item == fcchr) {
            // combine chromosome into one file...
            if (!Rplot_params->fcomb_chromo) {
                change_output_chr(Rplot_params->Routputfile, 0);
            }
            Rplot_params->fcomb_chromo = TOGGLE(Rplot_params->fcomb_chromo);
        } else if (item == ymin) {
            printf("Enter new value for Y-axis minimum > ");
            fcmap(stdin, "%g", &(Rplot_params->Yrange[0])); newline;
        } else if (item  == ymax) {
            printf("Enter new value for Y-axis maximum > ");
            fcmap(stdin, "%g", &(Rplot_params->Yrange[1])); newline;
        } else if (item == thr) {
            printf("Enter new Y-value for horizontal line > ");
            fcmap(stdin, "%g", &(Rplot_params->thresh)); newline;
        } else if (item == nr) {
            printf("Enter number of plot rows per page > ");
            fcmap(stdin, "%d", &(Rplot_params->nrows)); newline;
            if (Rplot_params->nrows > plots_per_file) {
                printf("Number of rows per per page exceeds number of plots in file.\n");
            }
        } else if (item == nc) {
            printf("Enter number of plot columns per page > ");
            fcmap(stdin, "%d", &(Rplot_params->ncols)); newline;
        } else if (item == mo) {
            // orientation...
            Rplot_params->mode = ((Rplot_params->mode == 'p')? 'l' : 'p');
        } else if (item == yf) {
            // Enforce Y-axis bounds even if data does not fit...
            Rplot_params->yfixed = TOGGLE(Rplot_params->yfixed);
        } else {
            printf("Invalid selection %s", cselect);
            strcpy(cselect, "");
        }
    } // while (!done)
    
    /* See if the output file has  a .ps extension */
    pos = strrchr(Rplot_params->Routputfile, '.');
    if (pos == NULL || ((strcmp(pos, ".pdf") && strcmp(pos, ".ps")))) {
        if (Rplot_params->fcomb_chromo || main_chromocnt <= 1) {
            /* just put in the .pdf extension */
            strcat(Rplot_params->Routputfile, ".pdf");
        } else {
            /* put in the first chromocome number, then
             stick the extension back in. */
            change_output_chr(Rplot_params->Routputfile,
                              global_chromo_entries[0]);
            strcat(Rplot_params->Routputfile, ".pdf");
        }
    }
    
    LoopOverTrait = old_loopovertrait;
}


int sw2_R_setup(analysis_type analysis, int numchr, linkage_ped_top *Top)
{

    char chr[4], *input_stat_names[5];
    int num_stats, *output_stats, num_select=0, i;
    R_plot_params_type Rplot_params;

    /* get selection of statistics */

    strcpy(Rperl_file, "Rsimwalk2.pl");

    if (main_chromocnt > 1) {
        strcpy(Rshell_file,"Rsimwalk2.all.sh");
        strcpy(Rscript_file,"Rsimwalk2.all.R");
        strcpy(Rmap_file,"Rsimwalk2.all.map");
        strcpy(Rdatafile, "RSW2DATA.01");
    } else {
        CHR_STR(numchr, chr);
        sprintf(Rdatafile, "RSW2DATA.%s", chr);
        sprintf(Rshell_file,"Rsimwalk2.%s.sh", chr);
        sprintf(Rmap_file,"Rsimwalk2.%s.map", chr);
        sprintf(Rscript_file,"Rsimwalk2.%s.R", chr);
    }

    if (analysis == TO_MERLIN) {
        num_stats = 2;
    } else {
        num_stats = 5;
    }

    for(i=0; i<num_stats; i++) {
        input_stat_names[i] = CALLOC((size_t) 10, char);
        switch(i) {
        case 0:
            if (analysis == TO_MERLIN) {
                strcpy(input_stat_names[0], "NPL_PAIR");
            } else {
                strcpy(input_stat_names[0], "BLOCKS");
            }
            break;
        case 1:
            if (analysis == TO_MERLIN) {
                strcpy(input_stat_names[1], "NPL_ALL");
            } else {
                strcpy(input_stat_names[1], "MAX-TREE");
            }
            break;
        case 2:
            strcpy(input_stat_names[2], "ENTROPY");
            break;
        case 3:
            strcpy(input_stat_names[3], "NPL_PAIRS");
            break;
        case 4:
            strcpy(input_stat_names[4], "NPL_ALL");
            break;
        }
    }

    output_stats = CALLOC((size_t) num_stats, int);
    stat_selection_menu(input_stat_names, num_stats,
                        output_stats, &num_select);
    if (num_select == 0) {
        return 0;
    }
    if (analysis == TO_MERLIN) {
        if (num_select >= 1) {
            output_stats[0] += 3;
        }
        if (num_select == 2) {
            output_stats[1] += 3;
        }
    }

    /* initialize the R parameters */
    init_plot_params(&Rplot_params);
    Rplot_params.LocusTop = Top->LocusTop;
    strcpy(Rplot_params.Routputfile, "SW2NPL.01.pdf");

    /* get user options for R params */
    R_plot_params(&Rplot_params, numchr);

    /* set up the perl script */
    /*  SW2R_perl_script();   */

    /* write the cshell script */
    convert_SW2R(Rplot_params, num_select, output_stats,
                 numchr, input_stat_names, num_stats, Top);

    draw_line();
    for(i=0; i<num_stats; i++) {
        free(input_stat_names[i]);
    }
    free(output_stats);
    return 1;
}


static void convert_allegro2R(R_plot_params_type Rplot_params,
			      int num_stats, int *output_stats,
			      int numchr, char *input_stat_names[],
			      int num_stat_names,
			      linkage_ped_top *Top)
{

    /* This function sets up the call to execute the perl script
       set up by allegro_perl_script which combines a number
       of allegro output files into a single R-readable file,
       if traits are to be combined, or into one file per trait

       Allegro has one file for each model called allegro_<model>.01
       Put all stats into one output file,
       if pcomb_traits is set then combine traits into one file
    */


    /* output_paths[i] contain the trait specific files */

    int tr, nloop, num_affec=num_traits;
    int st, ch;
    char *syscmd, chr[3];
    char rshfl[20 + FILENAME_LENGTH];
    FILE *shfp;


    sprintf(rshfl, "%s/%s", output_paths[0], Rshell_file);
    shfp = fopen(rshfl, "w");
    if (shfp == NULL) {
        sprintf(err_msg, "Could not open shell script %s for writing.\n",
     		rshfl);
        mssgf(err_msg);
        EXIT(FILE_WRITE_ERROR);
    }
    fprintf(shfp, "#!/bin/csh -f\n");
    /* print some identification */
    fprintf(shfp, "# R-plot C-shell file name: %s\n",
#ifdef HIDEPATH
            NOPATH
#else
            Rshell_file
#endif
        );
    script_time_stamp(shfp);
    fprintf(shfp, "#\n");
    /* print the names of the statistics */
    fprintf(shfp, "# Statistic file-names and numbers:\n");
    for(st=0; st < num_stat_names; st++) {
        fprintf(shfp, "#   %d)  %s\n",
                st+1, input_stat_names[st]);
    }

    NLOOP;
    for(ch=0; ch < main_chromocnt; ch++) {
        if (main_chromocnt > 1) {
            numchr = global_chromo_entries[ch];
        }

        change_output_chr(Rdatafile, numchr);
        CHR_STR(numchr, chr);

        if (Rplot_params.pcomb_traits) {
            remove_old_file_command(shfp, ".", Rdatafile);
            fprintf(shfp, "set message = `%s -s%s -c%d", Rperl_file,
                    chr, output_stats[0]+1);
            for (st=1; st < num_stats; st++) {
                fprintf(shfp, ",%d", output_stats[st]+1);
            }

            if (global_trait_entries[0] >= 0) {
                fprintf(shfp, " -t%s",
                        Rplot_params.LocusTop->Pheno[global_trait_entries[0]].TraitName);
            }
            for (tr=1; tr < num_traits; tr++) {
                if (global_trait_entries[tr] >= 0) {
                    fprintf(shfp, ",%s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                }
            }
            fprintf(shfp, " -o%s", Rdatafile);
            for (tr=0; tr <= nloop; tr++) {
                if (nloop > 1 && tr == 0) continue;
                fprintf(shfp, " %s", trait_paths[tr]);
                if (nloop == 1) break;
            }
            fprintf(shfp, "`\n");
            /*      fprintf(shfp, "echo\n"); */
            fprintf(shfp, "foreach m ($message)\n");
            fprintf(shfp, "    if (( $m == ERROR )) then\n");
            fprintf(shfp, "       echo \"Failed to convert Allegro output file.\"\n");
            fprintf(shfp, "       exit(-1)\n");
            fprintf(shfp, "    endif\n");
            fprintf(shfp, "end\n");
            fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
            fprintf(shfp, "   exit(-1)\n");
            fprintf(shfp, "endif\n");

        } else {
            if ((LoopOverTrait == 1 && num_traits == 2) || (num_traits > 2)) {
                /* only possible in the LoopOverTrait mode */
                for (tr=0; tr <=nloop; tr++) {
                    if (nloop > 1 && tr == 0) continue;
                    remove_old_file_command(shfp, trait_paths[tr], Rdatafile);
                    fprintf(shfp, "set message = `%s -s%s -c%d", Rperl_file,
                            chr, output_stats[0]+1);
                    for (st=1; st < num_stats; st++) {
                        fprintf(shfp, ",%d", output_stats[st]+1);
                    }
                    fprintf(shfp, " -t%s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                    fprintf(shfp, " -o%s/%s", trait_paths[tr], Rdatafile);
                    fprintf(shfp, " %s", trait_paths[tr]);
                    fprintf(shfp, "`\n");
                    /*	  fprintf(shfp, "echo \n "); */
                    fprintf(shfp, "foreach m ($message)\n");
                    fprintf(shfp, "    if (( $m == ERROR )) then\n");
                    fprintf(shfp, "       echo \"Failed to convert Allegro output file.\"\n");
                    fprintf(shfp, "       exit(-1)\n");
                    fprintf(shfp, "    endif\n");
                    fprintf(shfp, "end\n");
                    fprintf(shfp, "if (( ! -f %s/%s )) then\n", trait_paths[tr],
                            Rdatafile);
                    fprintf(shfp, "   exit(-1)\n");
                    fprintf(shfp, "endif\n");
                    if (nloop == 1) break;
                }
            } else {
                /* allegro output file is in current directory */
                remove_old_file_command(shfp, ".", Rdatafile);
                fprintf(shfp, "set message = `%s -s%s -c%d", Rperl_file,
                        chr, output_stats[0]+1);
                for (st=1; st < num_stats; st++) {
                    fprintf(shfp, ",%d",
                            output_stats[st]+1);
                }
                fprintf(shfp, " -t%s",
                        Rplot_params.LocusTop->Pheno[global_trait_entries[0]].TraitName);
                fprintf(shfp, " -o%s", Rdatafile);
                fprintf(shfp, " . ");
                fprintf(shfp, "`\n");
                /*	fprintf(shfp, "echo \n "); */
                fprintf(shfp, "foreach m ($message)\n");
                fprintf(shfp, "    if (( $m == ERROR )) then\n");
                fprintf(shfp, "       echo \"Failed to convert Allegro output file.\"\n");
                fprintf(shfp, "       exit(-1)\n");
                fprintf(shfp, "    endif\n");
                fprintf(shfp, "end\n");
                fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
                fprintf(shfp, "   exit(-1)\n");
                fprintf(shfp, "endif\n");
            }
        }
    }

    draw_graphs(Rplot_params, shfp, numchr, TO_Allegro, Top);
    fclose(shfp);
    syscmd=CALLOC((strlen(rshfl)+11), char);
    sprintf(syscmd, "chmod +x %s", rshfl);
    System(syscmd);
    free(syscmd);

    return;

}

int allegro_R_setup( int numchr,
		     linkage_ped_top *Top)
{

    char chr[3], *input_stat_names[12];
    int num_stats=12, *output_stats, num_select=0, i;
    R_plot_params_type Rplot_params;

    /* get selection of statistics */

    strcpy(Rperl_file, "Rallegro.pl");

    if (main_chromocnt > 1) {
        strcpy(Rdatafile, "RALLEGRO.01");
        strcpy(Rshell_file,"Rallegro.all.sh");
        strcpy(Rscript_file,"Rallegro.all.R");
        strcpy(Rmap_file,"Rallegro.all.map");
    } else {
        CHR_STR(numchr, chr);
        sprintf(Rdatafile, "RALLEGRO.%s", chr);
        sprintf(Rshell_file,"Rallegro.%s.sh", chr);
        sprintf(Rscript_file,"Rallegro.%s.R", chr);
        sprintf(Rmap_file,"Rallegro.%s.map", chr);
    }

    for(i=0; i < num_stats; i++) {
        input_stat_names[i] = CALLOC((size_t) 150, char);
    }

    strcpy(input_stat_names[0],
           "exppairs_mpt    Exponential multi-point pairs");
    strcpy(input_stat_names[1],
           "exppairs_spt    Exponential single-point pairs");
    strcpy(input_stat_names[2],
           "expall_mpt      Exponential multi-point all");
    strcpy(input_stat_names[3],
           "expall_spt      Exponential single-point all");
    strcpy(input_stat_names[4],
           "linpairs_mpt    Linear multi-point pairs");
    strcpy(input_stat_names[5],
           "linpairs_spt    Linear single-point pairs");
    strcpy(input_stat_names[6],
           "linall_mpt      Linear multi-point all");
    strcpy(input_stat_names[7],
           "linall_spt      Linear single-point all");
    strcpy(input_stat_names[8],
           "par_mpt:LOD     Parametric multi-point LOD scores");
    strcpy(input_stat_names[9],
           "par_spt:LOD     Parametric single-point LOD scores");
    strcpy(input_stat_names[10],
           "par_mpt:HLOD    Parametric multi-point heterogeneity LOD scores");
    strcpy(input_stat_names[11],
           "par_spt:HLOD    Parametric single-point heterogeneity LOD scores");

    output_stats = CALLOC((size_t) num_stats, int);
    stat_selection_menu(input_stat_names, num_stats,
                        output_stats, &num_select);

    if (num_select == 0) return 0;
    /* initialize the R parameters */
    init_plot_params(&Rplot_params);
    Rplot_params.LocusTop = Top->LocusTop;
    /* get user options for R params */
    strcpy(Rplot_params.Routputfile, "Allegro.01.pdf");

    R_plot_params(&Rplot_params, numchr);

    /* set up the perl script */
    /*  allegro_perl_script(num_stats, input_stat_names);  */

    /* write the cshell script */
    convert_allegro2R(Rplot_params, num_select, output_stats,
                      numchr, input_stat_names, num_stats, Top);

    draw_line();
    for(i=0; i<num_stats; i++) {
        free(input_stat_names[i]);
    }
    free(output_stats);
    return 1;
}

int merlin_R_setup(int numchr, linkage_ped_top *Top,
                   merlin_opt_type merlin_opt)
{

    char chr[4],  *input_stat_names[2];
    int num_stats, *output_stats, num_select=0;
/* analysis_type analysis=TO_MERLINONLY; */

    R_plot_params_type Rplot_params;

    /* get selection of statistics */

    sprintf(Rperl_file, "Rmerlin.py");

    if (main_chromocnt > 1) {
        strcpy(Rshell_file,"Rmerlin.all.sh");
        strcpy(Rscript_file,"Rmerlin.all.R");
        strcpy(Rmap_file,"Rmerlin.all.map");
        strcpy(Rdatafile, "RMERLINDATA.01");
    } else {
        CHR_STR(numchr, chr);
        sprintf(Rdatafile, "RMERLINDATA.%s", chr);
        sprintf(Rshell_file,"Rmerlin.%s.sh", chr);
        sprintf(Rscript_file,"Rmerlin.%s.R", chr);
        sprintf(Rmap_file,"Rmerlin.%s.map", chr);
    }

    num_stats = 2;

/*   if (merlin_opt.markernames &&  */
/*      (merlin_opt.npl || merlin_opt.pairs)) { */
/*     warnf("Plotting not yet supported for markerNames option."); */
/*     warnf("R-plots will not be generated."); */
/*     return 0; */
/*   } */

    if (merlin_opt.npl && merlin_opt.pairs) {
        num_stats = 2;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        input_stat_names[1] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "ALL");
        strcpy(input_stat_names[1], "Pairs");
    } else if (merlin_opt.npl) {
        num_stats=1;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "ALL");
    } else if (merlin_opt.qtl) {
        num_stats=1;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "QTL");
    } else if (merlin_opt.pairs) {
        num_stats=1;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "Pairs");
    } else if (merlin_opt.vc) {
        num_stats=1;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "VC");
    } else if (merlin_opt.parametric) {
        num_stats = 2;
        input_stat_names[0] = CALLOC((size_t) 6, char);
        input_stat_names[1] = CALLOC((size_t) 6, char);
        strcpy(input_stat_names[0], "LOD");
        strcpy(input_stat_names[1], "HLOD");
    } else {
        warnf("Plotting not set up, --npl, --pairs, --model and --vc were absent. ");
        return 0;
    }

    output_stats = CALLOC((size_t) num_stats, int);
    stat_selection_menu(input_stat_names, num_stats,
                        output_stats, &num_select);
    if (num_select == 0) return 0;

    /* initialize the R parameters */
    init_plot_params(&Rplot_params);
    Rplot_params.LocusTop = Top->LocusTop;
    strcpy(Rplot_params.Routputfile, "RMerlin.01.pdf");

    /* get user options for R params */
    R_plot_params(&Rplot_params, numchr);
    /* write the cshell script */
    convert_merlin2R(Rplot_params, num_select, output_stats,
                     numchr, input_stat_names, num_stats, Top);

    draw_line();
    if (num_stats > 1) {
        free(input_stat_names[0]);
    }
    if (num_stats == 2) {
        free(input_stat_names[1]);
    }
    free(output_stats);
    return 1;
}

/* Modified 7/21/2010 : if parametric analysis, get trait names
   from the table file itself */

static void convert_merlin2R(R_plot_params_type Rplot_params,
			     int num_stats, int *output_stats,
			     int numchr, char *input_stat_names[],
			     int num_stat_names,
			     linkage_ped_top *Top)
{

    /* This function sets up the call to execute the perl script
       Rmerlin.pl script which combines a number
       of Merlin output files into a single R-readable file,
       if traits are to be combined, or into one file per trait

       Merlin has one file for all models called merlin.<chr>.out.
       if pcomb_stats is set, then put all stats from one file
       into one output file, otherwise, they have to be broken up.
       if pcomb_traits is set then combine traits into one file
    */


    /* output_paths[i] contain the trait specific files */

    int tr, trr, nloop, num_affec=num_traits;
    int st, ch;
    char merlinoutfile[FILENAME_LENGTH], *syscmd, chr_str[3];
    FILE *shfp;
    char rshfl[20+FILENAME_LENGTH];
    char rdataf[2*FILENAME_LENGTH+1];

    sprintf(rshfl, "%s/%s", output_paths[0], Rshell_file);
    shfp = fopen(rshfl, "w");
    if (shfp == NULL) {
        sprintf(err_msg, "Could not open shell script %s for writing.\n",
     		rshfl);
        mssgf(err_msg);
        EXIT(FILE_WRITE_ERROR);
    }
    fprintf(shfp, "#!/bin/csh -f\n");
    /* print some identification */
    fprintf(shfp, "# R-plot C-shell file name: %s\n",
#ifdef HIDEPATH
            NOPATH
#else
            Rshell_file
#endif
        );
    script_time_stamp(shfp);
    fprintf(shfp, "#\n");
    /* print the statistic names */
    fprintf(shfp, "# Statistic names and numbers \n");
    for (st = 0; st < num_stat_names; st++) {
        fprintf(shfp, "#  %d) %s\n", st+1, input_stat_names[st]);
    }

    NLOOP;
    for(ch=0; ch < main_chromocnt; ch++) {
        if (main_chromocnt > 1) {
            numchr = global_chromo_entries[ch];
        }
        CHR_STR(numchr, chr_str);
        MERLINTABLE(chr_str, merlinoutfile);
        change_output_chr(Rdatafile, numchr);

        if (LoopOverTrait == 1 && num_traits >= 2) {
            if (Rplot_params.pcomb_traits) {
                remove_old_file_command(shfp, ".", Rdatafile);
                fprintf(shfp, "set message = `%s -c %s", Rperl_file,
                        input_stat_names[output_stats[0]]);
                for (st=1; st < num_stats; st++) {
                    fprintf(shfp, ",%s", input_stat_names[output_stats[st]]);
                }

                if (global_trait_entries[0] >= 0) {
                    fprintf(shfp, " -t %s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[0]].TraitName);
                }

                for (tr=1; tr < num_traits; tr++) {
                    if (global_trait_entries[tr] >= 0) {
                        fprintf(shfp, ",%s",
                                Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                    }
                }
                fprintf(shfp, " -o %s", Rdatafile);
                for (tr=0; tr <= nloop; tr++) {
                    if (nloop > 1 && tr == 0) continue;
                    fprintf(shfp, " %s/%s",	trait_paths[tr], merlinoutfile);
                    if (nloop == 1) break;
                }
                fprintf(shfp, "`\n");
                /*      fprintf(shfp, "echo\n"); */
                fprintf(shfp, "foreach m ($message)\n");
                fprintf(shfp, "    if (( $m == ERROR )) then\n");
                fprintf(shfp, "       echo\n");
                fprintf(shfp, "       echo \"Failed to convert Merlin output file.\"\n");
                fprintf(shfp, "       exit(-1)\n");
                fprintf(shfp, "    endif\n");
                fprintf(shfp, "end\n");
                fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
                fprintf(shfp, "   exit(-1)\n");
                fprintf(shfp, "endif\n");
            } else {
                /* separate traits into separate plots in LoopOverTrait mode */
                for (tr=0; tr <= nloop; tr++) {
                    if (nloop > 1 && tr == 0) continue;
                    remove_old_file_command(shfp, trait_paths[tr], Rdatafile);
                    fprintf(shfp, "set message = `%s -c %s", Rperl_file,
                            input_stat_names[output_stats[0]]);
                    for (st=1; st < num_stats; st++) {
                        fprintf(shfp, ",%s",
                                input_stat_names[output_stats[st]]);
                    }
                    fprintf(shfp, " -t %s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr-1]].TraitName);
                    fprintf(shfp, " -o %s/%s", trait_paths[tr], Rdatafile);
                    fprintf(shfp, " %s/%s", trait_paths[tr], merlinoutfile);
                    fprintf(shfp, "`\n");
                    /*	  fprintf(shfp, "echo\n"); */
                    fprintf(shfp, "foreach m ($message)\n");
                    fprintf(shfp, "    if (( $m == ERROR )) then\n");
                    fprintf(shfp, "       echo \n");
                    fprintf(shfp, "       echo \"Failed to convert Merlin output file.\"\n");
                    fprintf(shfp, "       exit(-1)\n");
                    fprintf(shfp, "    endif\n");
                    fprintf(shfp, "end\n");
                    fprintf(shfp, "if (( ! -f %s/%s )) then\n", trait_paths[tr],
                            Rdatafile);
                    fprintf(shfp, "   exit(-1)\n");
                    fprintf(shfp, "endif\n");
                    if (nloop == 1) break;
                }
            }
        } else if (LoopOverTrait == 0 && num_traits > 2) {
            if (!Rplot_params.pcomb_traits) {
                /* separate out traits into Rdata files in the current directory
                   tagged with trait names, so that they can be plotted separately */
                for (tr=0; tr < num_traits; tr++) {
                    if (global_trait_entries[tr] < 0) continue;
                    change_output_chr(Rdatafile, -9);
                    sprintf(rdataf, "%s_%s", Rdatafile,
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                    change_output_chr(rdataf, numchr);
                    remove_old_file_command(shfp, ".", rdataf);

                    fprintf(shfp, "set message = %s -c %s", Rperl_file,
                            input_stat_names[output_stats[0]]);
                    for (st=1; st < num_stats; st++) {
                        fprintf(shfp, ",%s",
                                input_stat_names[output_stats[st]]);
                    }
                    fprintf(shfp, " -t %s",
                            Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);

                    fprintf(shfp, " -o %s", rdataf);
                    fprintf(shfp, " %s", merlinoutfile);
                    fprintf(shfp, "\n");
                    /*	  fprintf(shfp, "echo\n"); */
                    fprintf(shfp, "foreach m ($message)\n");
                    fprintf(shfp, "    if (( $m == ERROR )) then\n");
                    fprintf(shfp, "       echo \n");
                    fprintf(shfp, "       echo \"Failed to convert Merlin output file.\"\n");
                    fprintf(shfp, "       exit(-1)\n");
                    fprintf(shfp, "    endif\n");
                    fprintf(shfp, "end\n");
                    fprintf(shfp, "if (( ! -f %s )) then\n", rdataf);
                    fprintf(shfp, "   exit(-1)\n");
                    fprintf(shfp, "endif\n");
                }
            } else {
                /* combine traits into a single Rdata file from merlin output
                   in current directory */
                remove_old_file_command(shfp, ".", Rdatafile);
                fprintf(shfp, "set message = `%s -c %s", Rperl_file,
                        input_stat_names[output_stats[0]]);
                for (st=1; st < num_stats; st++) {
                    fprintf(shfp, ",%s",
                            input_stat_names[output_stats[st]]);
                }

                trr=0;
                for (tr=0; tr < num_traits; tr++) {
                    if (global_trait_entries[tr] < 0)  continue;
                    if (trr == 0) {
                        fprintf(shfp, " -t %s",
                                Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                    } else {
                        fprintf(shfp, ",%s",
                                Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
                    }
                    trr++;
                }
                fprintf(shfp, " -o %s", Rdatafile);
                fprintf(shfp, " %s", merlinoutfile);
                fprintf(shfp, "`\n");
                /*	fprintf(shfp, "echo\n"); */
                fprintf(shfp, "foreach m ($message)\n");
                fprintf(shfp, "    if (( $m == ERROR )) then\n");
                fprintf(shfp, "       echo \n");
                fprintf(shfp, "       echo \"Failed to convert Merlin output file.\"\n");
                fprintf(shfp, "       exit(-1)\n");
                fprintf(shfp, "    endif\n");
                fprintf(shfp, "end\n");
                fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
                fprintf(shfp, "   exit(-1)\n");
                fprintf(shfp, "endif\n");
            }
        } else {
        /* else num_traits <= 2,  LoopOverTrait = 0 or num_traits = 1  */
            /* combine traits into a single Rdata file from merlin output
               in current directory */
            remove_old_file_command(shfp, ".", Rdatafile);
            fprintf(shfp, "set message = `%s -c %s", Rperl_file,
                    input_stat_names[output_stats[0]]);
            for (st=1; st < num_stats; st++) {
                fprintf(shfp, ",%s",
                        input_stat_names[output_stats[st]]);
            }

            if (num_traits == 2) {
                tr = (global_trait_entries[0] < 0) ? 1 : 0;
            } else {
                tr=0;
            }
            fprintf(shfp, " -t %s",
                    Rplot_params.LocusTop->Pheno[global_trait_entries[tr]].TraitName);
            fprintf(shfp, " -o %s", Rdatafile);
            fprintf(shfp, " %s", merlinoutfile);
            fprintf(shfp, "`\n");
            /*	fprintf(shfp, "echo\n"); */
            fprintf(shfp, "foreach m ($message)\n");
            fprintf(shfp, "    if (( $m == ERROR )) then\n");
            fprintf(shfp, "       echo \n");
            fprintf(shfp, "       echo \"Failed to convert Merlin output file.\"\n");
            fprintf(shfp, "       exit(-1)\n");
            fprintf(shfp, "    endif\n");
            fprintf(shfp, "end\n");
            fprintf(shfp, "if (( ! -f %s )) then\n", Rdatafile);
            fprintf(shfp, "   exit(-1)\n");
            fprintf(shfp, "endif\n");
        }
    }

    draw_graphs(Rplot_params, shfp, numchr, TO_MERLINONLY, Top);
    fclose(shfp);
    syscmd=CALLOC((strlen(rshfl)+11), char);
    sprintf(syscmd, "chmod +x %s", rshfl);
    System(syscmd);
    free(syscmd);

    return;

}

static void nplplot_init_header(R_plot_params_type Rplot_params,
				char *RdataFile, int trait, int chromosome)


{

    char *Rheaderfile;
    FILE *fp;

    Rheaderfile = CALLOC(strlen(RdataFile) + strlen(output_paths[trait]) + 7, char);
    sprintf(Rheaderfile, "%s/%s.hdr", output_paths[trait], RdataFile);

    fp=fopen(Rheaderfile, "w");
    if (fp == NULL) {
        sprintf(err_msg, "Could not open R header file %s for writing.\n",
     		Rheaderfile);
        mssgf(err_msg);
        EXIT(FILE_WRITE_ERROR);
    }
    fprintf(fp, "# reference line, y-axis minimum and y-axis maximum\n");
    fprintf(fp, "yline <- %4.2f\n ymin <- %4.2f\n ymax <- %4.2f\n ",
            Rplot_params.thresh,
            Rplot_params.Yrange[0], Rplot_params.Yrange[1]);
    fprintf(fp, "# Enforce y-axis bounds?\n");
    ((Rplot_params.yfixed == 1)?
     fprintf(fp, "yfix <- TRUE\n ") : fprintf(fp, "yfix=FALSE\n"));
    fprintf(fp, "# Plot subtitle\n");
    fprintf(fp, "title <- \"Chromosome %d\"\n", chromosome);
    fprintf(fp, "# Score units\n");
    fprintf(fp, "ylabel <- \"%s\"\n", Rplot_params.YLabel);
    fprintf(fp, "# Automatically adjusted by nplplot.multi.R\n");
    fprintf(fp, "# draw.lgnd <- FALSE\n");
    fprintf(fp, "# xlabl <- \"\"\n");
    fprintf(fp, "# lgndx <- NULL\n");
    fprintf(fp, "# lgndy <- NULL\n");
    fprintf(fp, "# Font scaling for legend\n");
    fprintf(fp, "cex.legend <- %5.4f\n", Rplot_params.legend_cex);
    fprintf(fp, "# Font scaling for axis\n");
    fprintf(fp, "cex.axis <- %5.4f\n", Rplot_params.axis_cex);
    fprintf(fp, "# Tick length for marker labels at the top of plots\n");
    fprintf(fp, "tcl <- %5.4f \n", Rplot_params.tcl);
    fprintf(fp, "# Use colors\n");
    fprintf(fp, "bw <- FALSE \n");
    fprintf(fp, "# Remove NAs before plotting each curve\n");
    fprintf(fp, "na.rm <- TRUE \n");
    fprintf(fp, "#---inserted by Mega2----- \n");
    fclose(fp);
}


void append_R_commands(char *shell_file, char *R_shell_file)

{

    char shfl[20+FILENAME_LENGTH];
    FILE *filep;

    sprintf(shfl, "%s/%s", output_paths[0], shell_file);
    filep = fopen(shfl, "a");

    if (filep == NULL) {
        fprintf(stderr, "Could not open shell script %s for writing.\n",
                shell_file);
        EXIT(FILE_WRITE_ERROR);
    }
    fprintf(filep,
            "echo Executing R shell script %s\n", R_shell_file);
    fprintf(filep, "./%s\n", R_shell_file);
    return;
}
