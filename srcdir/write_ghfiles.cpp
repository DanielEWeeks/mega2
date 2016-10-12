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
#include <math.h>

#include "common.h"
#include "typedefs.h"
#include "R_output.h"

#include "R_output_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "genetic_utils_ext.h"
#include "linkage_ext.h"
#include "makenucs_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"

#include "class_old.h"

/*
           R_output_ext.h:  allegro_R_setup append_R_commands
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
      genetic_utils_ext.h:  delete_file haldane_theta haldane_x kosambi_theta
            linkage_ext.h:  get_loci_on_chromosome new_lpedtop switch_penetrances
           makenucs_ext.h:  create_nuked_fids makenucs1
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  ped_id_width person_id_width
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line log_line script_time_stamp
        write_files_ext.h:  write_affection_data write_binary_data write_numbered_data
*/

#ifdef LTOP
#undef LTOP
#define LTOP Top->LocusTop
#endif


/*=============internal ============*/
static void    ghfile_name_select(int global_cshell, char *file_names[],
				  int *r_setup);
static int    global_cshell_option(void);
static void save_gh_peds(char *outfl_name,linkage_ped_top *Top,
		         double missing_quant, analysis_type analysis);
static void save_mlb_peds(char *outfl_name,linkage_ped_top *Top,
			  double missing_quant, analysis_type analysis);
static void write_al_input_cmd_file(char *gh_in, int *numchr,
				    char *gh_dat, char *gh_ped);
static void write_gh_input_cmd_file(char *gh_in, int *numchr,
				    char *gh_dat, char *gh_ped,
				    analysis_type analysis,
				    linkage_locus_top *LTop,
				    double unknown_quant);
static void write_al_cshell_file(char *gh_script, int *numchr, char *gh_in,
				 int global);
static void write_gh_cshell_file(char *gh_script, int *numchr, char *gh_in,
				 int first_time, analysis_type analysis,
				 linkage_locus_top *LTop);
static void write_gh_quant(FILE *filep, int locusnm,
			   linkage_ped_rec *entry,
			   double missing_quant,
			   analysis_type analysis);
static int gh_locus_file1(char *loutfl_name,
			  linkage_ped_top *Top, analysis_type analysis,
			  int sex_linked);
/*============exports================*/
void            create_gh_file(linkage_ped_top **LPedTreeTop,
			       file_format *infl_type, file_format *outfl_type,
			       int *numchr, analysis_type *analysis,
			       char *file_names[], int untyped_ped_opt,
			       linkage_ped_top **NukeTop);

void print_recomb_fracs(FILE *filep,
			int num_markers, int *markers,
			linkage_locus_top *LTop,
			analysis_type analysis, theta_type tt);

void write_gh_locus_file(char *loutfl_name, linkage_locus_top *LTop,
			 analysis_type analysis, int sex_linked);

/*===================================*/
/*
 * Save all the pedigrees to a file. Return the number of members saved (total).
 */
static void write_gh_quant(FILE *filep, int locusnm,
			   linkage_ped_rec *entry,
			   double missing_quant,
			   analysis_type analysis)
{
    fprintf(filep, "  ");
    if (fabs(entry->Pheno[locusnm].Quant - missing_quant) <= EPSILON) {
        if (analysis == TO_GeneHunter) {
            // http://linkage.rockefeller.edu/soft/gh/
            // Manual Section2. SCAN PEDIGREES
            //
            // A '-' in the phenotype/covariate data indicates missing data - NB:
            // 0 is a real value that a phenotype may take on and DOES NOT represent
            // missing phenotype data.
            fprintf(filep, "    -     ");
        } else  if (ITEM_READ(Value_Missing_Quant_On_Output)) {
            fprintf(filep, "%s", Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        } else {
            fprintf(filep, "%10.5f", missing_quant);
        }
    } else {
        fprintf(filep, "%10.5f", entry->Pheno[locusnm].Quant);
    }
}

static int gh_locus_file1(char *loutfl_name, linkage_ped_top *Top,
			  analysis_type analysis, int sex_linked)
{

    int m, tr, *trp, locus, allele, tmpi, tmpi2 = 0, need_dummy=1;
    int num_markers, *markers=NULL, nloop, num_affec=num_traits;
    linkage_locus_rec *Locus;
    pheno_rec *Pheno;
    FILE *filep;
    char lfl_name[2*FILENAME_LENGTH];


    NLOOP;

    if (num_traits > 0) {
        trp = &(global_trait_entries[0]);
        if (LoopOverTrait == 0) {
            if (num_traits == 1) {
                /* no trait */
                need_dummy = 0;
            } else {
                for (tr = 0; tr < num_traits; tr++) {
                    SKIP_TRI(tr);
                    if (Top->LocusTop->Locus[global_trait_entries[tr]].Type == AFFECTION) {
                        need_dummy = 0;
                        break;
                    }
                }
            }
        }
    } else {
        trp = NULL;
        need_dummy = 0;
    }
    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        if (LoopOverTrait == 1) {
            if (Top->LocusTop->Locus[*trp].Type == AFFECTION)
                need_dummy=0;
            else
                need_dummy = 1;
        }
        sprintf(lfl_name, "%s/%s", output_paths[tr], loutfl_name);

        if ((filep = fopen(lfl_name, "w")) == NULL) {
            errorvf("Unable to open locus file '%s'\n", lfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        /* set num_markers to write at the top of the locus file */
        if (trp == NULL) {
            num_markers = NumChrLoci;
        } else if (LoopOverTrait == 0) {
            if (analysis == TO_GHMLB) {
                num_markers = NumChrLoci - HasQuant + need_dummy;
            } else {
                num_markers = NumChrLoci + need_dummy;
            }
        } else  {
            if (Top->LocusTop->Locus[*trp].Type == QUANT &&
               analysis == TO_GeneHunter) {
                /* dummy affec and quant */
                num_markers = NumChrLoci - num_affec + 2;
            }
            else
                num_markers = NumChrLoci - num_affec + 1; /* only dummy */
        }
        /* write the headers of the locus file */
        fprintf(filep, "%d %d %d %d\n%d %.3f %.3f %d\n",
                num_markers, Top->LocusTop->RiskLocus,
                sex_linked,
                Top->LocusTop->Program, Top->LocusTop->MutLocus,
                Top->LocusTop->MutMale,
                Top->LocusTop->MutFemale, Top->LocusTop->Haplotype);

        for(tmpi=0; tmpi<num_markers; tmpi++)
            fprintf(filep, "%d ", tmpi+1);
        fprintf(filep, "\n");

        if (LoopOverTrait == 1) {
            Locus = &(Top->LocusTop->Locus[*trp]);
            Pheno = &(Top->LocusTop->Pheno[*trp]);
            switch(Locus->Type) {
            case AFFECTION:
                fprintf(filep, "%d %d", (int) Locus->Type - 1,
                        Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s", Locus->LocusName);
                fputc('\n', filep);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
                }
                fputc('\n', filep);
                fprintf(filep, "%d\n", Pheno->Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < Pheno->Props.Affection.ClassCnt; tmpi++) {
                    if (sex_linked) {
                        for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                        fputc('\n', filep);
                        for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                        fputc('\n', filep);
                    } else {
                        for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                        fputc('\n', filep);
                    }
                }
                break;
            case QUANT:
                /* dummy affection */
                fprintf(filep, "1  2 #dummy\n"); /* type, num_alleles, name */
                fprintf(filep, "0.5  0.5\n"); /* freq */
                fprintf(filep, " 1\n"); /* classcnt */
                fprintf(filep, "0.0  0.0  1.0\n"); /* penetrances */
                break;
            default:
                break;
            }
        } else {
            if (need_dummy) {
                /* dummy affection */
                fprintf(filep, "1  2 #dummy\n"); /* type, num_alleles, name */
                fprintf(filep, "0.5  0.5\n"); /* freq */
                fprintf(filep, " 1\n"); /* classcnt */
                fprintf(filep, "0.0  0.0  1.0\n"); /* penetrances */
            }
        }

        for(locus=0; locus < NumChrLoci; locus++) {
            Locus = &(Top->LocusTop->Locus[ChrLoci[locus]]);
            Pheno = &(Top->LocusTop->Pheno[ChrLoci[locus]]);
            switch(Locus->Type) {
            case AFFECTION:
                if (LoopOverTrait == 0) {
                    fprintf(filep, "%d %d", (int) Locus->Type - 1,
                            Locus->AlleleCnt);
                    if (Locus->LocusName != NULL)
                        fprintf(filep, " # %s", Locus->LocusName);
                    fputc('\n', filep);
                    for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                        fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
                    }
                    fputc('\n', filep);
                    fprintf(filep, "%d\n", Pheno->Props.Affection.ClassCnt);
                    for (tmpi = 0; tmpi < Pheno->Props.Affection.ClassCnt; tmpi++) {
                        if (sex_linked) {
                            for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                            fputc('\n', filep);
                            for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                            fputc('\n', filep);
                        } else {
                            for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                            fputc('\n', filep);
                        }
                    }
                }
                break;
            case QUANT:
                break;
            case NUMBERED:
            case BINARY:
                fprintf(filep, "%d %d", (int) Locus->Type - 1,
                        Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s", Locus->LocusName);
                fputc('\n', filep);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %.6f", Locus->Allele[allele].Frequency);
                }
                fputc('\n', filep);
                break;
            default:
                warnvf("Invalid locus type %d\n", (int) Locus->Type - 1);
                break;
            }
        }
        /* Now print the quants and covariates */
        if (analysis == TO_GeneHunter) {
            if (LoopOverTrait == 1) {
                /* print only the trait being looped over
                   (if it is a qtl) */
                if (Top->LocusTop->Locus[*trp].Type == QUANT) {
                    fprintf(filep, "0 2 ");
                    if (Top->LocusTop->Pheno[*trp].TraitName != NULL)
                        fprintf(filep, " # %s",
                                Top->LocusTop->Pheno[*trp].TraitName);
                    fprintf(filep, "\n\n\n\n\n\n");
                    /* 5 blank lines */
                }
            } else {
                /* print all traits that are qtls
                   not marked as covariates */
                if (num_affec > 0) {
                    for (tmpi = 0; tmpi < num_affec; tmpi++) {
                        SKIP_TRI(tmpi);
                        Locus
                            = &(Top->LocusTop->Locus[global_trait_entries[tmpi]]);
                        switch (Locus->Type) {
                        case AFFECTION:
                            break;
                        case QUANT:
                            fprintf(filep, "0 2 ");
                            if (Locus->LocusName != NULL)
                                fprintf(filep, " # %s", Locus->LocusName);
                            fprintf(filep, "\n\n\n\n\n\n");
                            /* 5 blank lines */
                            break;
                        default:
                            break;
                        }
                    }
                }
            }
            /* Now print the covariates */
            if (num_covariates > 0) {
                for(tmpi=0; tmpi < num_covariates; tmpi++) {
                    if (Top->LocusTop->Locus[covariates[tmpi]].LocusName != NULL)
                        fprintf(filep, "4 0  # %s\n",
                                Top->LocusTop->Locus[covariates[tmpi]].LocusName);
                    else
                        fprintf(filep, "4 0 \n");
                }
            }
        }

        /* recombination info */
        fprintf(filep, "%d %d\n", Top->LocusTop->SexDiff,
                Top->LocusTop->Interference);

        /* Now store the markers for the recombination */
        if (LoopOverTrait == 1) {
            /* first figure out the number of loci */
            if (analysis == TO_GHMLB) {
                /* 0 or 1 traits followed by the markers */
                num_markers=
                    NumChrLoci - num_affec +  ((trp != NULL)? 1 : 0);
            } else {
                if (trp != NULL) {
                    if (Top->LocusTop->Locus[*trp].Type == QUANT)
                        num_markers=NumChrLoci - num_affec + 2; /* dummy + quant */
                    else
                        num_markers=NumChrLoci - num_affec + 1; /* affec only */
                } else {
                    num_markers=NumChrLoci;
                }
            }

            markers=CALLOC((size_t) num_markers, int);
            tmpi=0;
            if (trp != NULL) {
                /* for both MLB and GH, add trait in the beginning */
                markers[tmpi++] = *trp;
            }
            /* now the numbered loci */
            for (m=0; m < NumChrLoci; m++) {
                if (Top->LocusTop->Locus[ChrLoci[m]].Type == NUMBERED ||
                   Top->LocusTop->Locus[ChrLoci[m]].Type == BINARY) {
                    markers[tmpi++] = ChrLoci[m];
                }
            }
            if (trp != NULL && analysis == TO_GeneHunter) {
                    if (Top->LocusTop->Locus[*trp].Type == QUANT)
                        markers[tmpi] = *trp; /* add the quant locus again at the end */
	    }
        } else {
            /* all loci */
            num_markers = NumChrLoci;
            markers=CALLOC((size_t) num_markers, int);

            if (analysis == TO_GeneHunter) {
                tmpi=0;
                for (m=0; m < NumChrLoci; m++) {
                    if (Top->LocusTop->Locus[ChrLoci[m]].Type != QUANT) {
                        markers[tmpi++]=ChrLoci[m];
                    } else {
                        /* store quant locus number */
                        tmpi2 = ChrLoci[m];
                    }
                }
            } else {
                if (HasQuant && need_dummy==0) {
                    /* affection as well as quant, skip the quant */
                    tmpi=0;
                    for (m=0; m < NumChrLoci; m++) {
                        if (Top->LocusTop->Locus[ChrLoci[m]].Type != QUANT) {
                            markers[tmpi++]=ChrLoci[m];
                        } else {
                            num_markers--;
                        }
                    }
                } else {
                    for (m=0; m < NumChrLoci; m++) {
                        markers[m]=ChrLoci[m];
                    }
                }
            }

            if (analysis == TO_GeneHunter) {
                if (HasQuant) {
                    for (m=tmpi; m < num_markers; m++) {
                        markers[m]= tmpi2;
                        /* making an assumption here that quants do
                           not have map-distances, so put first quant #quant times
                           at the end */
                    }
                }
                if (need_dummy) fprintf(filep, " 0.50000"); /* for the dummy */
            }
        }

#ifdef USEOLDMAPCODE
	if (sex_linked) {
            print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, MALE_THETA);
            print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, FEMALE_THETA);
	} else {
            print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, SEX_AVERAGED_THETA);
	}
#else /* USEOLDMAPCODE */
	// see user_input.c: get_genetic_distance_index()
	if (genetic_distance_index >= 0) {
            if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, SEX_AVERAGED_THETA);
            } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
                // With a SEX_SPECIFIC_GDMT when LTop->SexDiff == NO_SEX_DIFF, we should NOT print out
                // the male recombination fractions to a LINKAGE-format output file, since on the
                // X-chromosome, LTop->SexDiff must be set to NO_SEX_DIFF...
                if (Top->LocusTop->SexDiff != NO_SEX_DIFF) {
                    print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, MALE_THETA);
                }
                print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, FEMALE_THETA);
            } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
                print_recomb_fracs(filep, num_markers, markers, Top->LocusTop, analysis, FEMALE_THETA);
            }
	} else {
            // give an error if no genteic difference map was provided (at this moment one must be provided)
	}
#endif /* USEOLDMAPCODE */

        if (Top->LocusTop->Interference == NO_MAPPING_FUN)
            fprintf(filep, " %.5f\n", Top->LocusTop->FlankRecomb);

        fprintf(filep,"1 0.10000 0.45000\n");
        fclose(filep);
        free(markers);

        if (nloop == 1) break;
        trp++;
    }

    mssgvf("        Locus file:          %s\n", loutfl_name);
    return need_dummy;
}


//
// The MLB 3.0 software package and examples directory is located here...
// http://www.hgid.net/hgid/site/site.php?rubr=9

/* this has to be separate because the nuclear pedigrees
   are created offspring first
*/
static void  save_mlb_peds(char *outfl_name, linkage_ped_top *Top,
			   double missing_quant, analysis_type analysis)

{
    int tr, *trp, nloop, num_affec=num_traits, need_dummy=1;
    int i, tmpi, *markers=NULL, ped, entry, *locus, locus1, num_markers;
    linkage_ped_rec *Entry1;
    FILE *filep;
    char outfl[2*FILENAME_LENGTH], pformat[6], fformat[6];

    NLOOP;

    trp = &(global_trait_entries[0]);

    markers = CALLOC((size_t) NumChrLoci + 1, int);
    for (i=0, tmpi=0; i<NumChrLoci; i++) {
        if (LoopOverTrait == 0 ||
            Top->LocusTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            Top->LocusTop->Locus[ChrLoci[i]].Type == BINARY) {
            markers[tmpi++] = ChrLoci[i];
        }
    }
    num_markers = tmpi;

    if (LoopOverTrait == 0 && HasQuant && HasAff) {
        need_dummy =0;
    }

    /* Set the OrigIds values,
       we don't offer the user any choices here */

    /* Set the ID width */
    ped_id_width(Top, &tmpi);
    sprintf(fformat, "%%%dd ", tmpi);
    person_id_width(Top, &tmpi);
    sprintf(pformat, "%%%dd ", tmpi);

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;
        sprintf(outfl, "%s/%s", output_paths[tr], outfl_name);

        if ((filep = fopen(outfl, "w")) == NULL) {
            errorvf("Unable to open pedigree file '%s'\n", outfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL && UntypedPeds[ped]) continue;

            /* find the greatest allele number for each locus */

            for (entry = Top->Ped[ped].EntryCnt - 1; entry >= 0; entry--)  {
                Entry1 = &(Top->Ped[ped].Entry[entry]);
                /* write the pedigree and entry numbers */
                /* lwa , original fid ID  */
                fprintf(filep, fformat, Top->Ped[ped].Num);
                fprintf(filep, pformat, Entry1->ID);
                /* write the parents */
                fprintf(filep, pformat, Entry1->Father);
                fprintf(filep, pformat, Entry1->Mother);
                /* write the sex */
                if (Entry1->Sex == MALE_ID) fprintf(filep, "   %d", MALE_ID);
                else if (Entry1->Sex == FEMALE_ID) fprintf(filep, "   %d", FEMALE_ID);
                else fprintf(filep, "   0");
                if (LoopOverTrait == 1) {
                    /* first write the affection locus data */
                    locus1 = *trp;
                    switch(Top->LocusTop->Locus[locus1].Type) {
                    case QUANT:
                        // MLB_3.0_Source/readme.mlbgh
			//  - Data file format :
			//    - For the analysis of quantitative phenotypes, the difference with the 
			// previous format is that the quantitative phenotype is included between 
			// the sex and the binary phenotype. Therefore, the format is the following: 
			// Family-ID-Father-Mother-Sex-Quantitative phenotype-Binary phenotype-Markers.
			// Note that the binary phenotype information is required. If you have no binary
			// phenotype information, you should create a dummy binary variable.
			//
			// So since this locus is quantitative, we must create a dummy affection status
		        // information after writing the quantitative information...
                        write_gh_quant(filep, locus1, Entry1, missing_quant, analysis);
                        fprintf(filep, "  0");
                        break;
                    case AFFECTION:
                        // MLB_3.0_Source/readme.mlbgh
			//  - Data file format :
		        //- For the analysis of binary phenotypes (phenotype option is qualitative) 
                        // the format is the same as the classical GENEHUNTER program, that is :
	                //Family-ID-Father-Mother-Sex-Binary phenotype-Markers.
                        write_affection_data(filep, locus1,
                                             &(Top->LocusTop->Locus[locus1]), Entry1);
                        break;
                    default:
                        printf("Invalid locus type\n");
                        break;
                    }
                }
                /* now write marker genotype data */
                for (locus1 = 0; locus1 < num_markers; locus1++) {
                    locus = &(markers[locus1]);
                    switch(Top->LocusTop->Locus[*locus].Type) {
                    case QUANT:
                        if (LoopOverTrait == 0 ||
                            Top->LocusTop->Locus[*locus].Class == COVARIATE) {
                            write_gh_quant(filep, *locus, Entry1, missing_quant, analysis);
                            if (need_dummy == 1) {
                                /* write an affected phenotype */
                                fprintf(filep, " 2");
                            }
                        }
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0 ||
                            Top->LocusTop->Locus[*locus].Class == COVARIATE) {
                            write_affection_data(filep, *locus,
                                                 &(Top->LocusTop->Locus[*locus]), Entry1);
                        }
                        break;
                    case NUMBERED:
                        write_numbered_data(filep, *locus,
                                            &(Top->LocusTop->Locus[*locus]), Entry1);
                        break;

                    case BINARY:
                        write_binary_data(filep, *locus,
                                          &(Top->LocusTop->Locus[*locus]), Entry1);
                        break;
                    default:
                        warnvf("unknown locus type (locus %d) while writing output\n", *locus+1);
                    }
                }
                fprintf(filep, "\n");
            }
        }
        fclose(filep);
        if (nloop == 1) break;
        trp++;

    }

    free(markers);
    sprintf(err_msg, "        Pedigree file:       %s", outfl_name);
    mssgf(err_msg);
}

static void write_allegro_affection(FILE *filep, int locusnm,
				    linkage_locus_rec *locus, linkage_ped_rec *entry)

{
    fprintf(filep, "  ");
    if (entry->Pheno[locusnm].Affection.Status != UNDEF)
        fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Status);
    else fprintf(filep, "0");
    if (locus->Pheno->Props.Affection.ClassCnt > 1) {
        fprintf(filep, " ");
        if (entry->Pheno[locusnm].Affection.Class != UNDEF)
            fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Class);
        else
            fprintf(filep, "1");
    }
}

static void  save_gh_peds(char *outfl_name, linkage_ped_top *Top,
			  double missing_quant, analysis_type analysis)
{
    int tr, *trp, nloop, num_affec=num_traits, need_dummy=1;
    int i, tmpi, *markers=NULL, ped, entry, *locus, locus1, num_markers;
    linkage_ped_rec *Entry1;
    FILE *filep;
    char outfl[2*FILENAME_LENGTH], pformat[6], fformat[6];

    NLOOP;

    if (num_traits > 1 || (LoopOverTrait == 1 && num_traits > 0)) {
        trp = &(global_trait_entries[0]);
        if (LoopOverTrait == 0 && HasAff) {
            need_dummy=0;
        }
    } else {
        trp=NULL; need_dummy=0;
    }

    // An array of indicies that correspond only to makers...
    markers = CALLOC((size_t) NumChrLoci + 1, int);
    for (i=0, tmpi=0; i<NumChrLoci; i++) {
        if (LoopOverTrait == 0 ||
            Top->LocusTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            Top->LocusTop->Locus[ChrLoci[i]].Type == BINARY) {
            markers[tmpi++] = ChrLoci[i];
        }
    }
    num_markers = tmpi;

    /* Set the OrigIds values,
       we don't offer the user any choices here */

    /* Set the ID width */
    ped_id_width(Top, &tmpi);
    sprintf(fformat, "%%%dd ", tmpi);
    person_id_width(Top, &tmpi);
    sprintf(pformat, "%%%dd ", tmpi);

    for (tr=0; tr <= nloop; tr++) {
        if (tr == 0 && nloop > 1) continue;

        if (LoopOverTrait==1) {
	  need_dummy = (Top->LocusTop->Locus[*trp].Type == AFFECTION) ? 0 : 1;
        }
        sprintf(outfl, "%s/%s", output_paths[tr], outfl_name);

        if ((filep = fopen(outfl, "w")) == NULL) {
            errorvf("Unable to open pedigree file '%s'\n", outfl_name);
            EXIT(FILE_WRITE_ERROR);
        }

        for (ped = 0; ped < Top->PedCnt; ped++) {
            if (UntypedPeds != NULL && UntypedPeds[ped]) continue;

            /* find the greatest allele number for each locus */
            for (entry = 0; entry < Top->Ped[ped].EntryCnt; entry++) {
                Entry1 = &(Top->Ped[ped].Entry[entry]);

		// http://linkage.rockefeller.edu/soft/gh/
		// Manual Section2. SCAN PEDIGREES
		//
		// Each line of this [pedigree] file must have the following structure:
		//    3  12   8   9   1   2   1   1 2   8 3   0 0   4 6   1 3 ...
		//  (a) (b) (c) (d) (e) (f) (g)  (h ------------------------)
		// (a)  pedigree name
		// (b)  individual ID #
		// (c)  father's ID #
		// (d)  mother's ID #
		// (e)  sex (1=MALE, 2=FEMALE)
		// (f)  affectation status (1=UNAFFECTED, 2=AFFECTED)
		// (g)  liability class (OPTIONAL) - classes specified in marker data file
		// (h)  marker genotypes
		//
		// A '0' in any of the disease phenotype or marker genotype positions indicates
		// missing data.
		//
		// A '-' in the phenotype/covariate data indicates missing data - NB:
		// 0 is a real value that a phenotype may take on and DOES NOT represent
		// missing phenotype data.

		// http://www.helsinki.fi/~tsjuntun/autogscan/pedigreefile.html
		// Pedigree file
		//
		// Pre-makeped format [pedigree file] contains following columns, separated by space and/or tab characters:
		// Column 1: Pedigree identifier           { The identifier can be a number
		// Column 2: Individual's ID               { or a character string
		// Column 3: The individual's father       { If the person is a founder, just put a
		// Column 4: The individual's mother       { 0 in each column.
		// Column 5: Sex                           { 1 = Male, 2 = Female, Unknown sex is not permitted
		// Column 6+: Genetic data                 { Disease and marker phenotypes

                fprintf(filep, fformat, Top->Ped[ped].Num); // (a/1)
                fprintf(filep, pformat, Entry1->ID);  // (b/2)
                fprintf(filep, pformat, Entry1->Father); // (c/3)
                fprintf(filep, pformat, Entry1->Mother); // (d/4)
                if (Entry1->Sex == MALE_ID) fprintf(filep, "   %d", MALE_ID);
                else if (Entry1->Sex == FEMALE_ID) fprintf(filep, "   %d", FEMALE_ID);
                else fprintf(filep, "   0"); // (e/5)

                if (LoopOverTrait == 1) {
                    /* first write the affection locus data from global_trait_entries[] (f) */
                    locus1 = *trp;
                    switch(Top->LocusTop->Locus[locus1].Type) {
                    case QUANT:
                        /* don't write quant for GeneHunter */
                        if (analysis == TO_GeneHunter)
			  // GeneHunter (above) only deals with affection status in this position of the file
			  // write the affection status as missing here...
                            fprintf(filep, "  0");
                        else
                            write_gh_quant(filep, locus1, Entry1, missing_quant, analysis);
                        break;
                    case AFFECTION:
                        if (analysis == TO_Allegro) {
                            write_allegro_affection(filep, locus1,
                                                    &(Top->LocusTop->Locus[locus1]), Entry1);
                        } else {
                            write_affection_data(filep, locus1, &(Top->LocusTop->Locus[locus1]), Entry1);
                        }
                        break;
                    default:
                        printf("Invalid locus type\n");
                        break;
                    }
                }

                if (need_dummy && LoopOverTrait == 0 && analysis==TO_GeneHunter) {
                    fprintf(filep, "  0"); /* dummy affection (f) */
                }

                /* now write marker genotype data (h) */
                for (locus = &(markers[0]), tmpi = 0; tmpi < num_markers; locus++, tmpi++) {
                    switch(Top->LocusTop->Locus[*locus].Type) {
                    case QUANT:
                        if (LoopOverTrait == 0 && analysis != TO_GeneHunter) {
                            write_gh_quant(filep, *locus, Entry1, missing_quant, analysis);
                        }
                        break;
                    case AFFECTION:
                        if (LoopOverTrait == 0) {
                            write_affection_data(filep, *locus,
                                                 &(Top->LocusTop->Locus[*locus]), Entry1);
                        }
                        break;
                    case NUMBERED:
                        write_numbered_data(filep, *locus,
                                            &(Top->LocusTop->Locus[*locus]), Entry1);
                        break;

                    case BINARY:
                        write_binary_data(filep, *locus,
                                          &(Top->LocusTop->Locus[*locus]), Entry1);
                        break;
                    default:
                        warnvf("unknown locus type (locus %d) while writing output\n", *locus+1);
                    }
                }
                if (analysis == TO_GeneHunter) {
                    if (LoopOverTrait == 1) {
                        if (Top->LocusTop->Locus[*trp].Type == QUANT) {
                            /* write the single quant */
                            write_gh_quant(filep, *trp, Entry1, missing_quant, analysis);
                        }
                    } else if (!ManualReorder) {
                        /* write all quants */
                        for (i=0; i<num_traits; i++) {
                            SKIP_TRI(i)
                                tmpi = global_trait_entries[i];
                            if (Top->LocusTop->Locus[tmpi].Type == QUANT) {
                                write_gh_quant(filep, tmpi, Entry1, missing_quant, analysis);
                            }
                        }
                    }
                    /* write all covariates */
                    for (i=0; i<num_covariates; i++) {
                        tmpi = covariates[i];
                        write_gh_quant(filep, tmpi, Entry1, missing_quant, analysis);

                    }
                }
                fprintf(filep, "\n");
            }
        }
        trp++;
        fclose(filep);
        if (nloop == 1) break;
    }

    free(markers);
    sprintf(err_msg, "        Pedigree file:       %s", outfl_name);
    mssgf(err_msg);
}


/*
 * Write the locus data file in GeneHunter-format:
 1) Insert a space between the '#' and the locus name.
 2) Put the intermarker recombination fractions according to the map at the
 bottom of the file.
*/

// For a description of the linkage file format...
// http://linkage.rockefeller.edu/soft/linkage/
// 2.6 Description of Loci (DATAFILE)
//
// This routine is used for the following analysis types:
//  TO_LINKAGE, TO_NUKE (Create Nuclear Families), TO_LOD2, TO_PREMAKEPED, TO_SIMULATE
//
// This routine is the same as write_vitesse.c: write_vitesse_locus_file() which takes a vector of selected loci.
// Note, when using a π
void write_gh_locus_file(char *loutfl_name, linkage_locus_top *LTop,
			 analysis_type analysis, int sex_linked)
{
    int tr, *trp, *locus, allele, tmpi, tmpi2, tmpi3, i, warn_message = 0;
    int num_markers, num_markers_to_use, *markers=NULL, nloop, num_affec=num_traits;
    linkage_locus_rec *Locus;
    pheno_rec *Pheno;
    marker_rec *Marker;
    FILE *filep;
    char lfl_name[2*FILENAME_LENGTH];

    NLOOP;

    // The number of traits selected...
    trp = (num_traits > 0) ? &(global_trait_entries[0]) : NULL;

    /* since this routine is used by various analysis options that vary
       in their trait locus requirements, had to do the following :
       1. analysis = TO_NUKE or TO_LINKAGE
       all loci should be written in their given order, hence markers
       contains all sorts of loci.
       2. others
       only marker loci should be written after the affection status locus
    */

    // An array of indicies that correspond only to makers...
    // worst case size, using -1 as a centinal...
    markers = CALLOC((size_t) NumChrLoci + 1, int);
    for (i=0, tmpi=0; i<NumChrLoci; i++) {
        if (LoopOverTrait == 0 || // flag for looping over traits
            LTop->Locus[ChrLoci[i]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[i]].Type == BINARY) {
            markers[tmpi++] = ChrLoci[i];
        }
    }
    // We will always compute the recomnination map from the map file, and
    // not from any recombination data that is read by MEGA2
    num_markers_to_use = tmpi;
    num_markers = num_markers_to_use + ((LoopOverTrait == 0) ? 0 : 1);

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr == 0) continue;
        sprintf(lfl_name, "%s/%s", output_paths[tr], loutfl_name);

        if ((filep = fopen(lfl_name, "w")) == NULL) {
            errorvf("Unable to open locus file '%s'\n", lfl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        // num_markers: the number of the locus in the file
        // RiskLocus: 0 if risk is not to be calculated, or disease locus number (input order) if
	// risk is to be calculated...
        // sex_link: 0 for autosomal data; 1 for sex-linked data
        // Program: 1 CILINK; 2 CMAP; 3 ILINK; 4 LINKMAP; 5 MLINK; 6 LODSCORE; 7 CLODSCORE
        fprintf(filep, "%d %d %d %d\n",	// << no loci, risk locus, sexlinked (if 1), program code
                num_markers, LTop->RiskLocus, sex_linked, LTop->Program);
        // mutLocus: 0 if mutation rates are zero; mutation locus number (input order) for non-zero mutation rates
        // mutmale:  male mutation rate
        // mutfem: female mutation rate
        // disequil: 0 if loci are assumed to be in linkage equilibrium; 1 if loci are in linkage disequilibrium
        fprintf(filep, "%d %.3f %.3f %d\n", // << mut locus, mut mal, mut fem, hap freq (if 1)
                LTop->MutLocus, LTop->MutMale, LTop->MutFemale, LTop->Haplotype);

        // Write the locus order info on one line...
        for (tmpi = 0; tmpi < num_markers; tmpi++) {
            fprintf(filep, "%d ", tmpi+1);
        }
        fprintf(filep, "\n");
	// The loci are described in the order in which they appear in the PEDFILE
	// The descriptions differ according to the type of locus. A numeric code distinguishes
	// each of the types.
	// For the enums associated with Locus->Type see linkage.h:linkage_locus_type
        if (LoopOverTrait == 1) {
            /* affection locus first */
            Locus = &(LTop->Locus[*trp]);
            Pheno = &(LTop->Pheno[*trp]);
            trp++;
            switch (Locus->Type) {
            case AFFECTION:
            case QUANT:
                fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s\n", Locus->LocusName);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                }
                fprintf(filep, "\n");
                break;
            case BINARY:
            case NUMBERED:
                fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s\n", Locus->LocusName);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                }
                fprintf(filep, "\n");
                break;
            default:
                break;
            }
            switch(Locus->Type) {
            case AFFECTION:
                fprintf(filep, "%d\n", Pheno->Props.Affection.ClassCnt);
                for (tmpi = 0; tmpi < Pheno->Props.Affection.ClassCnt; tmpi++) {
                    if (sex_linked) {
                        for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                        fputc('\n', filep);
                        for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                        fputc('\n', filep);
                    } else {
                        for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                            fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                        fputc('\n', filep);
                    }
                }
                break;
            case QUANT:
                fprintf(filep, "%d\n", Pheno->Props.Quant.ClassCnt);
                for (tmpi = Pheno->Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                    for (tmpi2 = 0; tmpi2 < 3; tmpi2++)
                        fprintf(filep," %7.6f", Pheno->Props.Quant.Mean[tmpi-1][tmpi2]);
                fprintf(filep, " << GENOTYPE MEANS\n");
                /* ASSUMES only one trait per qtl */
                fprintf(filep," %7.6f\n", Pheno->Props.Quant.Variance[0][0]);
                fprintf(filep," %7.6f\n", Pheno->Props.Quant.Multiplier);
                break;
            default:
                warnvf("Invalid locus type\n");
                break;
            }
        }
        /* now the numbered, binary and covariates*/
        for (locus = &(markers[0]), tmpi3 = 0; tmpi3 < num_markers_to_use; locus++, tmpi3++) {
            Locus = &(LTop->Locus[*locus]);
            Pheno = &(LTop->Pheno[*locus]);
            Marker = &(LTop->Marker[*locus]);
            switch(Locus->Type) {
            case AFFECTION:
                if (LoopOverTrait == 0 || Locus->Class == COVARIATE) {
                    fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                    if (Locus->LocusName != NULL)
                        fprintf(filep, " # %s", Locus->LocusName);
                    fputc('\n', filep);
                    for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                        fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                    }
                    fputc('\n', filep);
                    fprintf(filep, "%d\n", Pheno->Props.Affection.ClassCnt);
                    for (tmpi = 0; tmpi < Pheno->Props.Affection.ClassCnt; tmpi++) {
                        if (sex_linked) {
                            for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]);
                            fputc('\n', filep);
                            for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]);
                            fputc('\n', filep);
                        } else {
                            for (tmpi2 = 0; tmpi2 < Pheno->Props.Affection.PenCnt; tmpi2++)
                                fprintf(filep, " %.4f", Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]);
                            fputc('\n', filep);
                        }
                    }
                }
                break;

            case QUANT:
                if (LoopOverTrait == 0 || Locus->Class == COVARIATE) {
                    fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                    if (Locus->LocusName != NULL)
                        fprintf(filep, " # %s", Locus->LocusName);
                    fputc('\n', filep);
                    for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                        fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                    }
                    fputc('\n', filep);
                    fprintf(filep, "%d\n", Pheno->Props.Quant.ClassCnt);
                    for (tmpi = Pheno->Props.Quant.ClassCnt; tmpi > 0; tmpi--)
                        for (tmpi2 = 0; tmpi2 < 3; tmpi2++) {
                            // local for debugging....
                            double f = Pheno->Props.Quant.Mean[tmpi-1][tmpi2];
                            fprintf(filep," %.6f", f);
                        }
                    fprintf(filep, " << GENOTYPE MEANS\n");
                    /* ASSUMES only one trait per qtl */
                    fprintf(filep," %.6f\n", Pheno->Props.Quant.Variance[0][0]);
                    fprintf(filep," %.6f\n", Pheno->Props.Quant.Multiplier);
                }
                break;

            case BINARY:
                fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s", Locus->LocusName);
		// THE COMMENT CAN COME HERE AFTER THE LOCUS NAME>...
                fputc('\n', filep);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                }
                fputc('\n', filep);
                fprintf(filep, "%d\n", Marker->Props.Binary.FactorCnt);
                for (tmpi = 0; tmpi < Marker->Props.Binary.FactorCnt; tmpi++) {
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++)
                        fprintf(filep, " %d", (int) Marker->Props.Binary.Factor[tmpi][tmpi2]);
                    fputc('\n', filep);
                }
                break;
            case NUMBERED:
                fprintf(filep, "%d %d", (int) Locus->Type - 1, Locus->AlleleCnt);
                if (Locus->LocusName != NULL)
                    fprintf(filep, " # %s", Locus->LocusName);
                fputc('\n', filep);
                for (allele = 0; allele < Locus->AlleleCnt; allele++) {
                    fprintf(filep, " %7.6f", Locus->Allele[allele].Frequency);
                }
                if (analysis == TO_LINKAGE && Locus->Marker->Props.Numbered.EstimateFrequencies &&
                    Locus->Marker->Props.Numbered.SelectOpt > 0) {
                    fprintf(filep, " %d ", Locus->Marker->Props.Numbered.NumAlleles);
                    switch(Locus->Marker->Props.Numbered.SelectOpt) {
                    case 1:
                        fprintf(filep, "founder alleles");
                        break;
                    case 2:
                        fprintf(filep, "founder + random alleles");
                        break;
                    case 3:
                        fprintf(filep, "founder + all new alleles");
                        break;
                    case 4:
                        fprintf(filep, "all alleles");
                        break;
                    default:
                        break;
                    }
                }
                fputc('\n', filep);
                break;
            default:
                warnvf("Invalid locus type %d\n", (int) Locus->Type);
                break;
            }
        }

        // Recombination Information:
	//
	// sex diference:
	// 0 = no sex-difference
	// 1 = constant sex-difference (the ratio of female/male genetic distance is the same in all intervals)
	// 2 = variable sex-difference (the female/male distance ratio can be different in each interval)
	//
	// interference:
	// 0 = no interference
	// 1 = interference without a mapping function
	// 2 = user-specified mapping function
	// sexdiff interferene << sex difference, interference
        fprintf(filep, "%d %d\n", LTop->SexDiff, LTop->Interference);

        if (genetic_distance_index >= 0) {
	  // a genetic map was specified by the user...
            if (LTop->Interference == NO_MAPPING_FUN)
                fprintf(filep, " %.5f\n", LTop->FlankRecomb);
            switch (LTop->SexDiff) {
                case NO_SEX_DIFF:
                    if (LoopOverTrait == 1) fprintf(filep, " 0.50000");
                    if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                        print_recomb_fracs(filep, num_markers_to_use, markers, LTop, analysis, SEX_AVERAGED_THETA);
                    } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
			       genetic_distance_sex_type_map == FEMALE_GDMT) {
		      // On the X chromosome, LTop->SexDiff must be set to NO_SEX_DIFF. see...
		      // annotated_ped_file.c:copy_exmap_locmap()
		      // Here we should print out only the female recombination fractions to a
		      // LINKAGE-format output file.
                        print_recomb_fracs(filep, num_markers_to_use, markers, LTop, analysis, FEMALE_THETA);
                    } else if (warn_message == 0) {
		        // These messages should be printed only once....
		        warnvf("A sex-averaged genetic map must be used to compute sex averaged recombination fractions.\n");
			warn_message = 1;
                    }
                    break;
                case CONSTANT_SEX_DIFF:
                    fprintf(filep, "%.3f\n", LTop->FMRatio);
                    break;
                case VARIABLE_SEX_DIFF:
                    if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
			genetic_distance_sex_type_map == FEMALE_GDMT) {
                        if (genetic_distance_sex_type_map != FEMALE_GDMT) {
                            if (LoopOverTrait == 1) fprintf(filep, " 0.50000");
                            print_recomb_fracs(filep, num_markers_to_use, markers, LTop, analysis, MALE_THETA);
                        } else if (warn_message == 0) {
                            warnvf("A sex-specific genetic map must be used to compute male recombination fractions.\n");
                            warn_message = 1;
                        }
                        if (LoopOverTrait == 1) fprintf(filep, " 0.50000");
                        print_recomb_fracs(filep, num_markers_to_use, markers, LTop, analysis, FEMALE_THETA);
                    } else if (warn_message == 0) {
                        warnvf("A sex-specific or a female genetic map must be used to compute recombination fractions.\n");
                        warn_message = 1;
                    }
                    break;
                default:
		    warnvf("Unknown sex difference code %d, ignoring\n", LTop->SexDiff);
		    break;
            }
            fprintf(filep,"1 0.10000 0.45000\n");
            fclose(filep);
            if (nloop == 1) break;
        }
    }
    if (analysis == TO_GeneHunterPlus || analysis == TO_Allegro) {
        sprintf(err_msg, "        Locus file:          %s", loutfl_name);
        mssgf(err_msg);
    }
    free(markers);
}

static void          write_al_input_cmd_file(char *gh_in, int *numchr,
					     char *gh_dat,
					     char *gh_ped)
{
    int           i, nloop, num_affec = num_traits;
    char          chr_str[3], gh_in1[2*FILENAME_LENGTH];
    FILE          *fp;

    CHR_STR(*numchr, chr_str);
    NLOOP;

    for (i=0; i<=nloop; i++) {
        if (nloop > 1 && i==0) continue;
        sprintf(gh_in1, "%s/%s", output_paths[i], gh_in);
        if ((fp = fopen(gh_in1, "w")) == NULL)  {
            errorvf("Unable to open command file '%s'\n", gh_in1);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(fp, "PREFILE %s\n", gh_ped);
        fprintf(fp, "DATFILE %s\n", gh_dat);
        fprintf(fp, "STEPS 2\n");
        fprintf(fp, "MODEL mpt par het allegro_par_mpt.%s\n", chr_str);
        fprintf(fp, "MODEL spt par het allegro_par_spt.%s\n", chr_str);
        fprintf(fp, "MODEL mpt exp pairs equal allegro_exppairs_mpt.%s\n", chr_str);
        fprintf(fp, "MODEL spt exp pairs equal allegro_exppairs_spt.%s\n", chr_str);
        fprintf(fp, "MODEL mpt exp all equal allegro_expall_mpt.%s\n", chr_str);
        fprintf(fp, "MODEL spt exp all equal allegro_expall_spt.%s\n", chr_str);
        fprintf(fp, "MODEL mpt lin pairs equal allegro_linpairs_mpt.%s\n", chr_str);
        fprintf(fp, "MODEL spt lin pairs equal allegro_linpairs_spt.%s\n", chr_str);
        fprintf(fp, "MODEL mpt lin all equal allegro_linall_mpt.%s\n", chr_str);
        fprintf(fp, "MODEL spt lin all equal allegro_linall_spt.%s\n", chr_str);
        fclose(fp);

        if (nloop == 1) break;
    }
    sprintf(err_msg, "        Command file:        %s", gh_in);
    mssgf(err_msg);
}

static void   write_gh_input_cmd_file(char *gh_in, int *numchr,
				      char *gh_dat, char *gh_ped,
				      analysis_type analysis,
				      linkage_locus_top *LTop,
				      double unknown_quant_val)
{

    int           *trp, i, nloop, num_affec = num_traits;
    char          chr_str[3], gh_in1[2*FILENAME_LENGTH];
    char          gh_out[FILENAME_LENGTH];
    FILE          *fp;

    CHR_STR(*numchr, chr_str);
    if (analysis == TO_GeneHunter)
        sprintf(gh_out, "gh.%s.out", chr_str);
    else if (analysis == TO_GeneHunterPlus)
        sprintf(gh_out, "ghp.%s.out", chr_str);
    else if (analysis == TO_GHMLB)
        sprintf(gh_out, "mlb.%s.out", chr_str);

    NLOOP;

    trp= &(global_trait_entries[0]);
    for (i=0; i<= nloop; i++) {
        if (nloop > 1 && i==0) continue;
        sprintf(gh_in1, "%s/%s", output_paths[i], gh_in);
        if ((fp = fopen(gh_in1, "w")) == NULL)  {
            errorvf("Unable to open command file '%s'\n", gh_in1);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(fp, "photo %s\n", gh_out);
        if (analysis == TO_GHMLB && LTop->Locus[*trp].Type == QUANT) {
            fprintf(fp, "phenotype quantitative\n");
            if (ITEM_READ(Value_Missing_Quant_On_Output))
                fprintf(fp, "unknown %s\n", Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
            else
                fprintf(fp, "unknown %6.2f\n", unknown_quant_val);
        }
        fprintf(fp, "load %s\n", gh_dat);
        fprintf(fp, "use\n");
        fprintf(fp, "score\n");
        fprintf(fp, "scan %s\n", gh_ped);
        if (analysis != TO_GHMLB) {
            fprintf(fp, "ps on\n");
            fprintf(fp, "total het\n");
            fprintf(fp, "npl_all.%s.ps\n", chr_str);
            fprintf(fp, "lod.%s.ps\n", chr_str);
            fprintf(fp, "info.%s.ps\n", chr_str);
        }
        fprintf(fp, "quit\n");
        fclose(fp);
        if (nloop == 1) break;
        trp++;
    }
    mssgvf("        Command file:        %s\n", gh_in);
}


static void   write_al_cshell_file(char *gh_script, int *numchr, char *gh_in,
				   int global)
{

    int    i, nloop, num_affec=num_traits;
    char   gh_scr[2*FILENAME_LENGTH], *syscmd;
    char   global_scr[2*FILENAME_LENGTH];
    FILE   *fp;

    NLOOP;

    if (strcmp(output_paths[0], ".")) {
        sprintf(global_scr, "%s/%s",  output_paths[0], gh_script);
    } else {
        strcpy(global_scr, gh_script);
    }
    change_output_chr(global_scr, 0);


    for (i=0; i<= nloop; i++) {
        if (i == 0 && global == 1) {
            /*  create a top level shell script */
            if ((fp = fopen(global_scr, "w")) == NULL) {
                errorvf("Unable to open shell file '%s'\n", global_scr);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(fp, "#!/bin/csh -f\n");
            /* print some identification */
            fprintf(fp, "# C-shell file name: %s\n",
#ifdef HIDEPATH
                    NOPATH
#else
                    global_scr
#endif
                );
            script_time_stamp(fp);
            fclose(fp);
            syscmd = CALLOC((strlen(global_scr) + 10), char);
            sprintf(syscmd, "chmod +x %s", global_scr);
            System(syscmd);
            free(syscmd);
        }

        if (nloop > 1 && i==0) continue;
        sprintf(gh_scr, "%s/%s", output_paths[i], gh_script);
        if ((fp = fopen(gh_scr, "w")) == NULL) {
            errorvf("Unable to open %s for writing.\n", gh_scr);
            EXIT(FILE_WRITE_ERROR);
        }

        fprintf(fp, "#!/bin/csh -f\n");
        script_time_stamp(fp);
        fprintf(fp, "echo Starting Allegro run on chromosome %d\n",*numchr);
        fprintf(fp, "allegro %s\n", gh_in);
        fclose(fp);

        syscmd = CALLOC((strlen(gh_scr) + 10), char);
        sprintf(syscmd, "chmod +x %s", gh_scr);
        System(syscmd);
        free(syscmd);

        if (global) {
            /* combined cshell for all chromosomes */
            if ((fp = fopen(global_scr, "a")) == NULL) {
                errorvf("Unable to open %s for appending.\n", global_scr);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(fp, "echo Running Alllegro shell script %s\n", gh_script);
            if (strcmp(trait_paths[i], ".")) {
                fprintf(fp, "cd %s\n", trait_paths[i]);
            }
            fprintf(fp, "./%s\n", gh_script);
            if (strcmp(trait_paths[i], ".")) {
                fprintf(fp, "cd ..\n");
            }
            fclose(fp);
        }
        if (nloop == 1) break;
    }
}

static void   write_gh_cshell_file(char *gh_script, int *numchr, char *gh_in,
				   int first_time, analysis_type analysis,
				   linkage_locus_top *LTop)
{

    int             i, nloop, num_affec=num_traits, i2;
    char            chr_str[3], gh_scr[2*FILENAME_LENGTH], syscmd[100];
    char            ghp[4];
    FILE            *fp;
    int  *locus, tmpi, tmpi2, *locus1, ntick;
    double gap, slop, tot, *position, diff, prev;
    int *markers=NULL, num_markers;
    linkage_locus_rec *Locus;
    int   xlinked;

    xlinked =(((LTop->SexLinked == 2 && *numchr == SEX_CHROMOSOME) ||
               (LTop->SexLinked == 1))? 1 : 0);

    /*  prog_name(analysis, &(prog[0])); */
    NLOOP;
    CHR_STR(*numchr, chr_str);

    if (analysis == TO_GeneHunterPlus) {
        strcpy(ghp, "ghp");
    }
    else strcpy(ghp, "gh");

    for (i=0; i<=nloop; i++) {
        if (nloop > 1 && i==0) continue;
        sprintf(gh_scr, "%s/%s", output_paths[i], gh_script);

        if (first_time==1)  {
            delete_file(gh_scr);
            fp = fopen(gh_scr, "w");
        }
        else fp = fopen(gh_scr, "a");

        if (fp == NULL) {
            errorvf("Unable to open script file '%s'\n", gh_scr);
            EXIT(FILE_WRITE_ERROR);
        }
        if (first_time) {
            fprintf(fp, "#!/bin/csh -f\n");
            script_time_stamp(fp);
        }
        fprintf(fp, "if (-e %s.%s.out.Z) then\n", ghp, chr_str);
        fprintf(fp, "   echo Moving %s.%s.out.Z to %s.%s.out.old.Z\n",
                ghp, chr_str, ghp, chr_str);
        fprintf(fp, "   /bin/mv %s.%s.out.Z %s.%s.out.old.Z \n",
                ghp, chr_str, ghp, chr_str);
        fprintf(fp, "endif\n");
        fprintf(fp, "echo Starting %s run on chromosome %d\n", ProgName, *numchr);
        fprintf(fp, "unalias rm\n");
        fprintf(fp, "rm -f %s.%s.out npl_all.%s.ps lod.%s.ps info.%s.ps\n",
                ghp, chr_str, chr_str, chr_str, chr_str);
        if (analysis == TO_GeneHunterPlus)
            fprintf(fp, "rm -f nullprobs.dat probs.dat kac.dat\n");
        fprintf(fp, "touch %s.%s.sum\n", ghp, chr_str);
        fprintf(fp, "date >> %s.%s.sum\n", ghp, chr_str);

        fprintf(fp, "echo %s run on chromosome %d >> %s.%s.sum\n",
                ProgName, *numchr, ghp, chr_str);
        if (xlinked) {
            fprintf(fp, "x%s <<DATA > /dev/null \n", ghp);
        } else {
            fprintf(fp, "%s <<DATA > /dev/null \n", ghp);
        }
        fprintf(fp, "run %s\n", gh_in);
        fprintf(fp, "DATA\n");
        fprintf(fp, "grep -i error %s.%s.out >> %s.%s.sum\n", ghp, chr_str,
                ghp, chr_str);
        fprintf(fp, "grep -i warn  %s.%s.out >> %s.%s.sum\n", ghp, chr_str,
                ghp, chr_str);
        fprintf(fp, "grep -i fail  %s.%s.out >> %s.%s.sum\n", ghp, chr_str,
                ghp, chr_str);
        fprintf(fp, "%s  '\\\n", awk_str);
        fprintf(fp, "BEGIN { n = 0 } \\\n");
        fprintf(fp, "/ score/ { n = 0 } \\\n");
        fprintf(fp, "/ps on/ { n = 0 } \\\n");
        fprintf(fp, "/file to store/ { n = 0 } \\\n");
        fprintf(fp, "{ if ( n == 1 ) \\\n");
        fprintf(fp, "  { \\\n");
        fprintf(fp, "  print $0 \\\n");
        fprintf(fp, "  } \\\n");
        fprintf(fp, "} \\\n");
        /* markers */
        fprintf(fp, "/Current/ { print $0; n = 1 } \\\n");
        fprintf(fp, "/TOTAL/ { print $0; n = 1 } \\\n");
        fprintf(fp, "/Totalling/ { n = 1 } ' %s.%s.out >> %s.%s.sum\n",
                ghp, chr_str, ghp, chr_str);
        if (analysis == TO_GeneHunterPlus)  {
            fprintf(fp, "kac 500\n");
            fprintf(fp, "echo GeneHunter-PLUS KAC output >> ghp.%s.sum\n",chr_str);
            fprintf(fp, "echo \"Position NPL_Z-Score   Zlr          GHP_Lodscore    delta\" >> ghp.%s.sum\n",chr_str);
            fprintf(fp, "cat kac.dat >> ghp.%s.sum\n",chr_str);
            /* Make XMGR graphics file: */
            fprintf(fp, "#  A post-processing script to be used on\n");
            fprintf(fp, "#  ghp.#.out and kac.dat that are created by \n");
            fprintf(fp, "#  a run of GeneHunter-Plus.  This will extract\n");
            fprintf(fp, "#  the lod score, the heterogeneity lod score, and\n");
            fprintf(fp, "#  the GeneHunter-Plus lod score into xmgr format,\n");
            fprintf(fp, "#  which can then be graphed right away.\n");
            fprintf(fp, "# Assumptions:\n");
            fprintf(fp, "#  1) That ghp.#.out and kac.dat exist\n");
            fprintf(fp, "#  2) That the lod score under heterogeneity has been computed\n");
            fprintf(fp, "#  3) That it is O.K. to write over the file 'tmp'\n");
            fprintf(fp, "unalias rm\n");
            fprintf(fp, "rm -f tmp ghp_xmgr.%s\n",chr_str);
            fprintf(fp, "touch ghp_xmgr.%s\n",chr_str);
            fprintf(fp, "# Extract the desired table from the end of ghp.#.out.\n");
            fprintf(fp, "%s  '\\\n", awk_str);
            fprintf(fp, "BEGIN { n = 0 } \\\n");
            fprintf(fp, "/file to store/ { n = 0 } \\\n");
            fprintf(fp, "{ if ( n == 1 ) \\\n");
            fprintf(fp, "  { \\\n");
            fprintf(fp, "  print $0 \\\n");
            fprintf(fp, "  } \\\n");
            fprintf(fp, "} \\\n");
            fprintf(fp, "/HLOD/ { n = 1 } ' ghp.%s.out > tmp\n", chr_str);
            fprintf(fp, "# Edit it to remove the parentheses and comma\n");
            fprintf(fp, "ed tmp << DATA >& /dev/null\n");
            fprintf(fp, "g/(/s// /\n");
            fprintf(fp, "g/)/s// /\n");
            fprintf(fp, "g/,/s// /\n");
            fprintf(fp, "w\n");
            fprintf(fp, "q\n");
            fprintf(fp, "DATA\n");
            fprintf(fp, "echo \"@ s10 linestyle 5\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@ s11 linestyle 1\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@ s11 linewidth 2\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@ s12 linestyle 3\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend on \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend loctype view \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend layout 0 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend vgap 2 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend hgap 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend length 4 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box on \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box fill on \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box fill with color \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box fill color 0 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box fill pattern 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box color 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box linewidth 2 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend box linestyle 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend x1 0.672154 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend y1 0.94958 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend font 4 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend char size 0.690000 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend linestyle 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend linewidth 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@    legend color 1 \" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \\@    legend string 10 \\\"Lod score\\\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \\@    legend string 11 \\\"HLOD\\\"  >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \\@    legend string 12 \\\"GHP lodscore\\\"  >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TARGET S10\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TYPE XY\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "# Extract the lod score curve into s0\n");
            fprintf(fp, "%s '\\\n", awk_str);
            fprintf(fp, "{ print $1/100.0, $2 } ' tmp >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"&\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TARGET S11\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TYPE XY\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "# Extract the heterogeneity lod score curve into s1\n");
            fprintf(fp, "%s '\\\n", awk_str);
            fprintf(fp, "{print $1/100.0, $4 } ' tmp >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"&\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TARGET S12\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@TYPE XY\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "# Extract the GeneHunter-Plus lod score into s2\n");
            fprintf(fp, "%s '\\\n", awk_str);
            fprintf(fp, "{print $1/100.0, $4 } ' kac.dat >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"&\" >> ghp_xmgr.%s \n", chr_str);
            /* Now add on the marker names & positions for the xmgr file */
            /*  NOTE: Need to add in treatment of gaps  */
            fprintf(fp, "echo \"@autoscale yaxes\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis bar on\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis ticklabel on\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis ticklabel start type spec\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis ticklabel char size 0.75\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis tick major on\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis tick size 0.5\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@zeroxaxis tick minor size 0.25\" >> ghp_xmgr.%s \n", chr_str);

            position = CALLOC((size_t) LTop->LocusCnt + 1, double);
            position[0] = 0.0000;
            
            markers = CALLOC((size_t) LTop->LocusCnt + 1, int);
            for (i2=0, tmpi=0; i2<NumChrLoci; i2++) {
                if (LTop->Locus[ChrLoci[i2]].Type == NUMBERED ||
                    LTop->Locus[ChrLoci[i2]].Type == BINARY) {
                    markers[tmpi++] = ChrLoci[i2];
                }
            }
	    num_markers = tmpi;

            /* Note: Need to be careful re if the input map is Haldane or Kosambi */
            /* tot is the total map distance from the first marker to the last,
               in Morgans */
            /* This won't work:
               tot = LTop->Marker[markers[tmpi-1]].pos_avg - LTop->Marker[markers[0]].pos_avg;
               because GeneHunter-Plus is using the Haldane position, even if
               the input map is in Kosambi */
	    if (genetic_distance_index < 0 || genetic_distance_sex_type_map != SEX_AVERAGED_GDMT) {
	      errorvf("For computation of the map distance a genetic sex averaged map must be specified.\n");
	      EXIT(INPUT_DATA_ERROR);
	    }
            tot = 0;
            ntick = 0;
            if ((LTop->LocusCnt - num_affec) > 1) {
	      for (locus = &(markers[0]), locus1 = &(markers[1]), tmpi = 0;
                   tmpi < num_markers-1;
                   locus++, locus1++, tmpi++) {
                    if (LTop->Marker[*locus].chromosome != LTop->Marker[*locus1].chromosome)
                        diff = 0.50000;
                    else {
                        diff = LTop->Marker[*locus1].pos_avg - LTop->Marker[*locus].pos_avg;
                        if (diff < 0.0)   {
                            warnvf("loci are not in correct order\n");
                            warnvf("    Locus %d (%s) is at position %5.1f\n",
                                    *locus,LTop->Locus[*locus].LocusName,
                                    LTop->Marker[*locus].pos_avg);
                            warnvf("    Locus %d (%s) is at position %5.1f\n",
                                    *locus1,LTop->Locus[*locus1].LocusName,
                                    LTop->Marker[*locus1].pos_avg);
                        } else {
                            diff = ((LTop->map_distance_type == 'h') ? diff : haldane_x(kosambi_theta(diff)));
                            tot += diff;
                            position[++ntick] = tot;
                        }
                    }
                }
            }

            tot /= 100.0;
            gap = tot/50.0;
            /* If the map length is 0.465, then ticks can be as close as 0.01;
               If the map length is 1.000, then ticks can be as close as 0.02. */
            if (tot > 1.0) {
                fprintf(fp, "echo \"@zeroxaxis tick major 0.5\" >> ghp_xmgr.%s \n", chr_str);
                fprintf(fp, "echo \"@zeroxaxis tick minor 0.25\" >> ghp_xmgr.%s \n", chr_str);
                slop = 0.05;
            } else {
                fprintf(fp, "echo \"@zeroxaxis tick major 0.1\" >> ghp_xmgr.%s \n", chr_str);
                fprintf(fp, "echo \"@zeroxaxis tick minor 0.05\" >> ghp_xmgr.%s \n", chr_str);
                slop = 0.01;
            }
            fprintf(fp, "echo \"@xaxis ticklabel type spec\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@xaxis ticklabel char size 0.5\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@xaxis ticklabel layout spec\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@xaxis ticklabel angle 80\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, "echo \"@xaxis tick type spec\" >> ghp_xmgr.%s \n", chr_str);

            /* now the numbered and binary */
            /* Print here the number of ticks */
            /* Looks like one can set the number of user-defined ticks to
               be larger than the actual number of lables without xmgr
               complaining */
            /* Need a loop across ticks and markers, starting with tick zero */
            prev = -1.0;
            ntick = 0;
            for (locus = &(markers[0]), tmpi2 = 0; tmpi2 < num_markers; locus++, tmpi2++) {
                Locus = &(LTop->Locus[*locus]);
                /* Note: This doesn't work if the input map is Kosambi, as GeneHunter-Plus
                   will have been given the Haldane map */
                /* Skip those ticks that were too close to the previous one */
                if (position[tmpi2]/100.0 - prev > gap) {
                        fprintf(fp, "echo \"@xaxis tick %d, %.3f\" >> ghp_xmgr.%s \n",
				ntick, position[tmpi2]/100.0, chr_str);
                        if (Locus->LocusName != NULL)
                            fprintf(fp, "echo \'@xaxis ticklabel %d, \"%s\"\' >> ghp_xmgr.%s \n",
				    ntick, Locus->LocusName, chr_str);
                        ntick++;
                        prev = position[tmpi2]/100.0;
                    }
                /* example:
                   @xaxis tick 1, 0.172
                   @xaxis ticklabel 1, "D22S345"
                */
            }
            fprintf(fp, "echo \"@xaxis tick spec %d\" >> ghp_xmgr.%s \n", ntick, chr_str);
            free(markers);
            free(position);
            /* fprintf(fp, " echo \"@with g0\" >> ghp_xmgr.%s \n", chr_str); */
            /* Need the 'slop' factor here */
            fprintf(fp, " echo \"@world xmin -%.2f\" >> ghp_xmgr.%s \n", slop, chr_str);
            fprintf(fp, " echo \"@world xmax %.2f\" >> ghp_xmgr.%s \n", tot + slop, chr_str);
            fprintf(fp, " echo \"@with line\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, " echo \"@line on\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, " echo \"@line loctype world\" >> ghp_xmgr.%s \n", chr_str);
            /* fprintf(fp, " echo \"@line g0\" >> ghp_xmgr.%s \n", chr_str); */
            /* Put in the line here at -2 for exclusion */
            fprintf(fp, " echo \"@line -%.2f, -2.0, %.2f, -2.0\" >> ghp_xmgr.%s \n",slop,tot+slop, chr_str);
            fprintf(fp, " echo \"@line linewidth 1\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, " echo \"@line linestyle 1\" >> ghp_xmgr.%s \n", chr_str);
            fprintf(fp, " echo \"@line def\" >> ghp_xmgr.%s \n", chr_str);
        }
        fprintf(fp, "echo \"--------------------------------\" >> %s.%s.sum \n",
                ghp, chr_str);
        fprintf(fp, "compress %s.%s.out\n", ghp, chr_str);
        fprintf(fp,"echo GeneHunter%s run completed - See %s.%s.sum for summary\n",
                ((analysis==TO_GeneHunterPlus) ? "-Plus" : ""), ghp, chr_str);
        fprintf(fp, "echo The full output may be found in %s.%s.out.Z%s\n",
                ghp, chr_str,
                ((analysis == TO_GeneHunterPlus)? " and kac.dat" : ""));
        fprintf(fp, "echo See also npl_all.%s.ps lod.%s.ps info.%s.ps \n",
                chr_str, chr_str, chr_str);
        if (analysis == TO_GeneHunterPlus)
            fprintf(fp, "echo The xmgr graphics file may be found in ghp_xmgr.%s\n",
                    chr_str);

        fclose(fp);
        sprintf(syscmd, "chmod +x %s", gh_scr);
        System(syscmd);
        if (nloop == 1) break;
    }
}

/* adapted from write_gh_cshell_file */
static void write_mlb_cshell_file(char *gh_script, int *numchr,
				  char *gh_in, int first_time)

{
    int   i, nloop, num_affec=num_traits;
    FILE            *fp;
    char  chr_str[3], gh_scr[2*FILENAME_LENGTH], syscmd[100];

    NLOOP;
    CHR_STR(*numchr, chr_str);

    for (i=0; i<=nloop; i++) {
        if (nloop > 1 && i==0) continue;
        sprintf(gh_scr, "%s/%s", output_paths[i], gh_script);

        if (first_time)  {
            delete_file(gh_scr);
            fp = fopen(gh_scr, "w");
        }
        else fp = fopen(gh_scr, "a");

        if (fp == NULL) {
            errorvf("Unable to open script file '%s'\n", gh_scr);
            EXIT(FILE_WRITE_ERROR);
        }
        if (first_time) {
            fprintf(fp, "#!/bin/csh -f\n");
            script_time_stamp(fp);
        }

        fprintf(fp, "echo Starting MLBQTL run on chromosome %d\n", *numchr);
        fprintf(fp, "unalias rm\n");
        fprintf(fp, "if (-e mlb.%s.out.Z) then\n", chr_str);
        fprintf(fp, "   echo Moving mlb.%s.out.Z to mlb.%s.out.old.Z\n",
                chr_str, chr_str);
        fprintf(fp, "   /bin/mv mlb.%s.out.Z mlb.%s.out.old.Z \n",
                chr_str, chr_str);
        fprintf(fp, "endif\n");
        fprintf(fp, "touch mlb.%s.sum\n", chr_str);
        fprintf(fp, "date >> mlb.%s.sum\n", chr_str);
        fprintf(fp, "mlb <<DATA > /dev/null \n");
        fprintf(fp, "run %s\n", gh_in);
        fprintf(fp, "DATA\n");
        fprintf(fp, "grep -i error mlb.%s.out >> mlb.%s.sum\n",  chr_str,
                chr_str);
        fprintf(fp, "grep -i warn  mlb.%s.out >> mlb.%s.sum\n",  chr_str,
                chr_str);
        fprintf(fp, "grep -i fail  mlb.%s.out >> mlb.%s.sum\n",  chr_str,
                chr_str);
        fprintf(fp, "%s  '\\\n", awk_str);
        fprintf(fp, "BEGIN { n = 0 } \\\n");
        fprintf(fp, "/ score/ { n = 0 } \\\n");
        fprintf(fp, "/ quit/  { n = 0 } \\\n");
        fprintf(fp, "{ if (n == 1) \\\n");
        fprintf(fp, "     { \\\n");
        fprintf(fp, "     print $0 \\\n");
        fprintf(fp, "     } \\\n");
        fprintf(fp, "} \\\n");
        fprintf(fp, "/Current/ { print $0; n = 1 } \\\n");
        fprintf(fp, "/MLB-/ { print $0; n = 1 } ' mlb.%s.out >> mlb.%s.sum\n",
                chr_str, chr_str);
        fprintf(fp, "echo \"--------------------------------\" >> mlb.%s.sum \n",
                chr_str);
        fprintf(fp, "compress mlb.%s.out\n", chr_str);
        fprintf(fp,"echo MLBQTL run completed - See mlb.%s.sum for summary\n",
                chr_str);
        fprintf(fp, "echo The full output may be found in mlb.%s.out.Z\n",
                chr_str);

        fclose(fp);
        sprintf(syscmd, "chmod +x %s", gh_scr);
        System(syscmd);
        if (nloop == 1) break;
    }
}


static int      global_cshell_option(void)
{
    char creply[10], stars[3];
    int item=-1, opt=1;

    if (DEFAULT_OPTIONS) {
        mssgf("Separate c-shell scripts will be written (DEFAULT)");
        mssgf("as requested in the batch file.");
        return 0;
    }

    sprintf(stars, "* ");
    printf("Global c-shell script option:\n");
    while (item) {
        printf("0) Done with this menu - please proceed.\n");
        printf("%c 1) Separate C-shell script for each chromosome\n",
               stars[0]);
        printf("%c 2) Single C-shell script for all chromosomes\n",
               stars[1]);
        printf("Enter option 0 - 2 > ");
        fcmap(stdin, "%s", creply); printf("\n");
        sscanf(creply, "%d", &item);
        switch(item) {
        case 0:
            break;
        case 1:
            sprintf(stars, "* ");
            opt=1;
            break;
        case 2:
            sprintf(stars, " *");
            opt = 2;
            break;
        default:
            warn_unknown(creply);
            break;
        }
    }
    return opt - 1;
}


static void    ghfile_name_select(int global_cshell, char *file_names[],
				  int *r_setup)

{
    int choice = -1;
    char cchoice[10], stem[5], fl_stat[12];
    int ask_r_setup = *r_setup;

    if (DEFAULT_OUTFILES) return;

    if (main_chromocnt > 1) {
        change_output_chr(file_names[0], -9);
        change_output_chr(file_names[1], -9);
        change_output_chr(file_names[2], -9);
        if (global_cshell)
            change_output_chr(file_names[3], 0);
        else
            change_output_chr(file_names[3], -9);
        strcpy(stem, "stem");
    } else {
        strcpy(stem, "");
    }

    while (choice != 0) {
        printf("0) Done with this menu - please proceed.\n");
        if (main_chromocnt == 1) {
            printf(" 1) Pedigree file name:       %-15s\t%s\n", file_names[0],
                   file_status(file_names[0], fl_stat));
            printf(" 2) Locus file name:          %-15s\t%s\n", file_names[1],
                   file_status(file_names[1], fl_stat));
            printf(" 3) Command file name:        %-15s\t%s\n", file_names[2],
                   file_status(file_names[2], fl_stat));
        } else {
            printf(" 1) Pedigree file name %s:  %-15s\n", stem, file_names[0]);
            printf(" 2) Locus file name %s:     %-15s\n", stem, file_names[1]);
            printf(" 3) Command file name %s:   %-15s\n", stem, file_names[2]);
        }
        if (global_cshell || main_chromocnt == 1) {
            printf(" 4) C-shell file name:        %-15s\t%s\n", file_names[3],
                   file_status(file_names[3], fl_stat));
        }
        else
            printf(" 4) C-shell file name %s:   %-15s\n", stem, file_names[3]);

        if (ask_r_setup) {
            printf(" 5) Generate R-plots for Allegro:  [%s]\n",
                   yorn[*r_setup]);
            printf("Enter option 0 - 5 > ");
        } else {
            printf("Enter option: 0 - 4 > ");
        }

        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice)  {
        case 0:
            break;
        case 1:
            printf("Enter NEW pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]); newline;
            break;
        case 2:
            printf("Enter NEW locus file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]); newline;
            break;
        case 3:
            printf("Enter NEW command file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]); newline;
            break;
        case 4:
            printf("Enter NEW C-Shell script file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]);
            break;
        case 5:
            if (ask_r_setup) {
                *r_setup = TOGGLE(*r_setup);
                break;
            }
        default:
            warn_unknown(cchoice);
            break;
        }
    }

    if (main_chromocnt > 1 && !global_cshell) {
        /* put back the .sh extension on the script file */
        change_output_chr(file_names[3], global_chromo_entries[0]);
        strcat(file_names[3], ".sh");
    }
}

static void init_cshell_file(char *cshell_file, int global)

{
    int tr, nloop, num_affec = num_traits;
    char fl_name[2*FILENAME_LENGTH];
    FILE *gh_script_p;

    NLOOP;

    for (tr=0; tr <= nloop; tr++) {
        if (nloop > 1 && tr==0) continue;
        sprintf(fl_name, "%s/%s", output_paths[tr], cshell_file);
        delete_file(fl_name);
        if ((gh_script_p = fopen(fl_name, "w")) == NULL) {
            errorvf("Unable to open script file '%s'\n", fl_name);
            EXIT(FILE_WRITE_ERROR);
        }
        fprintf(gh_script_p,"#!/bin/csh -f\n");
        script_time_stamp(gh_script_p);
        fclose(gh_script_p);
        if (global) {
            sprintf(fl_name, "%s/%s", output_paths[0], cshell_file);
            if ((gh_script_p = fopen(fl_name, "a")) == NULL) {
                errorvf("Could not open %s for writing.\n", cshell_file);
                EXIT(FILE_WRITE_ERROR);
            }
            fprintf(gh_script_p,
                    "echo Running analysis for trait #%d\n", tr);
            if (nloop > 1) {
                fprintf(gh_script_p, "cd %s\n", trait_paths[tr]);
                fprintf(gh_script_p, "./%s\n", cshell_file);
                fprintf(gh_script_p, "cd ..\n");
            }
            fclose(gh_script_p);
        }
        if (nloop == 1) break;
    }
}

#ifdef LPedTreeTop
#undef LPedTreeTop
#endif
#define LPedTreeTop (*Top)

void            create_gh_file(linkage_ped_top **Top,
			       file_format *infl_type, file_format *outfl_type,
			       int *numchr, analysis_type *analysis,
			       char *file_names[], int untyped_ped_opt,
			       linkage_ped_top **NukeTop)
{

    int             l, created_dummy=0;
    char            *outfl_name, *loutfl_name, *gh_in, *gh_script;
    int             first_time, global_cshell;
    double          missing_quant = QUNDEF;
    int             sex_linked;
    int             *traits=&(global_trait_entries[0]);
    int             R_setup=0, global=0;

    linkage_ped_top *NewTop = NULL;

    /* ask the user if he wants a global shell script */

    if (main_chromocnt > 1 && *analysis != TO_Allegro) {
        global_cshell=global_cshell_option();
        /* strip the chromosome number from the file names */
    } else {
        global_cshell=0;
    }
    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        log_line(mssgf);
        if (*analysis == TO_Allegro) R_setup=1;
    } else {
        print_outfile_mssg();
        if (*analysis == TO_GeneHunter) {
            printf("GeneHunter output file selection menu:\n");
        } else if (*analysis == TO_GeneHunterPlus) {
            printf("GeneHunter-plus output file selection menu:\n");
        } else if (*analysis == TO_Allegro) {
            printf("Allegro output file selection menu:\n");
            R_setup = 1;
        } else if (*analysis == TO_GHMLB) {
            printf("GeneHunter-MLB output file selection menu:\n");
        } else {
            errorvf("INTERNAL: GH analysis selection should not reach this statement.\n");
        }

        ghfile_name_select(global_cshell, file_names, &R_setup);
    }

    // missing_quant would retain it's default value of QUNDEF if there were no
    // quantitative traits, or if this were TO_Allegro, or TO_GeneHunterPlus analysls.
    if (*analysis == TO_GHMLB || *analysis == TO_GeneHunter) {
        // If there are any traits create a 'missing_quant'...
        for (l=0; l < num_traits; l++) {
	  if (global_trait_entries[l] == -1) continue;
	  if (LPedTreeTop->LocusTop->Locus[traits[l]].Type == QUANT) {
	    missing_quant = MissingQuant;
	    break;
	  }
        }
    }

    if (num_covariates > 0)
        missing_quant = MissingQuant;

/*   if (*analysis== TO_GeneHunter) { */
/*     covariates=CALLOC((size_t) num_traits, int); */
/*     for(l=0;  l<num_traits; l++) */
/*       covariates[l]=0; */
/*     gh_cov_selection(num_traits, global_trait_entries, */
/*  		    covariates, LPedTreeTop->LocusTop); */
/*   } */

/* Removed the check for LoopOverTrait==1, always 1 for Allegro */

    if (*analysis == TO_Allegro) {
        global = ((main_chromocnt > 1 || (num_traits > 1 && LoopOverTrait == 1))? 1: 0);
    } else {
        if (main_chromocnt ==  1 || global_cshell==1) {
            global = 0;
            (global_cshell)?
                change_output_chr(file_names[3], 0) :
                change_output_chr(file_names[3], *numchr);
            init_cshell_file(file_names[3], global);
            first_time=0;
        }
    }

    if (R_setup) {
        R_setup = allegro_R_setup(*numchr, *Top);
    }

    if (*analysis == TO_GHMLB) {
        NewTop = new_lpedtop();
        makenucs1(LPedTreeTop, NULL, NewTop);
        create_nuked_fids(NewTop);
    }

    for (l=0; l< main_chromocnt; l++) {
        if (main_chromocnt > 1) {
            if (!global_cshell) first_time=1;
            *numchr=global_chromo_entries[l];
            change_output_chr(file_names[0], *numchr);
            change_output_chr(file_names[1], *numchr);
            change_output_chr(file_names[2], *numchr);
        }

        (global_cshell)?
            change_output_chr(file_names[3], 0) :
            change_output_chr(file_names[3], *numchr);
        get_loci_on_chromosome(*numchr);
        sex_linked =
            (((*numchr == SEX_CHROMOSOME &&
	       LPedTreeTop->LocusTop->SexLinked == 2) ||
              (LPedTreeTop->LocusTop->SexLinked == 1)
	      ) ? 1 : 0);
        if (sex_linked) {
            if (*analysis == TO_GeneHunter) {
                warnf("Files not created for X-Chromosome.");
                continue;
            }
        }

        outfl_name = file_names[0];
        loutfl_name = file_names[1];
        gh_in = file_names[2];
        gh_script = file_names[3];

        if (l==0) {
            create_mssg(*analysis);
        }

        if (*analysis == TO_Allegro)
            write_al_input_cmd_file(gh_in, numchr, loutfl_name, outfl_name);
        else
            write_gh_input_cmd_file(gh_in, numchr, loutfl_name, outfl_name,
                                    *analysis, LPedTreeTop->LocusTop,
                                    missing_quant);
        if (*analysis == TO_GeneHunter || *analysis == TO_GHMLB) {
            created_dummy=gh_locus_file1(loutfl_name, LPedTreeTop, *analysis, sex_linked);
        } else {
            write_gh_locus_file(loutfl_name, LPedTreeTop->LocusTop, *analysis, sex_linked);
        }

        if (*analysis == TO_GHMLB) {
            if (UntypedPeds != NULL) {
                free(UntypedPeds);
            }

            UntypedPeds = CALLOC((size_t) NewTop->PedCnt, int);
            omit_peds(untyped_ped_opt, NewTop);
        } else {
            omit_peds(untyped_ped_opt, LPedTreeTop);
        }

        if (*analysis == TO_GHMLB) {
            save_mlb_peds(outfl_name, NewTop, missing_quant, *analysis);
        } else {
            save_gh_peds(outfl_name, LPedTreeTop, missing_quant,
                         *analysis);
        }
        if (*analysis == TO_Allegro) {
            write_al_cshell_file(gh_script, numchr, gh_in, global);
            if (global) global++;
        } else if (*analysis == TO_GHMLB) {
            write_mlb_cshell_file(gh_script, numchr, gh_in, first_time);
        } else {
            write_gh_cshell_file(gh_script, numchr, gh_in, first_time,
                                 *analysis, LPedTreeTop->LocusTop);
        }

        if (*analysis == TO_Allegro) {
            sprintf(err_msg, "        C-shell script:        %s", gh_script);
            mssgf(err_msg);
        } else {
            if (!global_cshell || (l == (main_chromocnt-1))) {
                sprintf(err_msg, "        C-shell script:      %s", gh_script);
                mssgf(err_msg);
            }
        }

    }


    if (*analysis == TO_Allegro) {
        if (R_setup) {
            mssgf("Additional files:");
            if (global) {
                change_output_chr(file_names[3], 0);
            }
            append_R_commands(file_names[3], Rshell_file);
            sprintf(err_msg, "        R-script file:       %s", Rscript_file);
            mssgf(err_msg);
            sprintf(err_msg, "        Shell file to run R: %s", Rshell_file);
            mssgf(err_msg);
            if (base_pair_position_index >= 0) {
                 // Needed to produce UCSC Custom Tracks which is only done if a physical map file exists.
                 sprintf(err_msg, "        Physical map for R:  %s", Rmap_file);
                 mssgf(err_msg);
            }
            if (global) {
                sprintf(err_msg,
                        "        Combined shell-script: %s", file_names[3]);
                mssgf(err_msg);
            }
        }

        if (R_setup) {
            log_line(mssgf);
            sprintf(err_msg,
                    "Run %s in %s for all analysis and creation of R-plots",
                    file_names[3],
#ifdef HIDEPATH
                    NOPATH
#else
                    ((strcmp(output_paths[0], ".")?
                      output_paths[0] : "current directory"))
#endif
                );
            mssgf(err_msg);
        }
    }

    if (created_dummy==1) {
        draw_line();
        sprintf(err_msg, "No affection status locus was present.");
        mssgf(err_msg);
        sprintf(err_msg, "Added a dummy affection status locus.");
        mssgf(err_msg);
    }
    /* if (*analysis == TO_GeneHunter) free(covariates); */
    if (*analysis == TO_GHMLB && HasQuant) {
        /*   prog_name(*analysis, progname); */
        warnvf("Quantitative phenotypes should be normalized for %s.\n", ProgName);
    }

    if (*analysis == TO_GHMLB) {
        *NukeTop = NewTop;
    }

    /* free_all_including_lpedtop(&LPedTreeTop_copy); */
}

#undef LPedTreeTop

/* We should not need to recreate the list of markers here,
   Checked for input with and without traits, and
   combined chromosomes.
*/

/* void print_recomb_fracs(FILE *filep, int num_markers,  */
/* 			int *markers,  */
/* 			linkage_locus_top *LTop, */
/* 			analysis_type analysis) */


/* { */

/*   int *numbered; */
/*   int i, num_numbered, tmpi; */
/*   int i; */
/*   double *thetas; */

/*   switch (num_markers) { */
/*   case 0: */
/*   case 1: */
/*     break; */
/*   default: */
/*     numbered = CALLOC((size_t) num_markers, int); */
/*     thetas=CALLOC((size_t) num_markers-1, double); */
/*     for (i=0; i < num_markers-1; i++) { */
/*       thetas[i] = -99.0; */
/*     } */
/*     num_numbered=0; */
/*     for (i=0; i < num_markers; i++) { */
/*       if (LTop->Marker[markers[i]].pos_avg >= 0.0) { */
/* 	numbered[num_numbered]=i; num_numbered++; */
/*       } */
/*     } */
/*     if (num_numbered > 1) { */
/*       for (i=0; i < num_numbered-1; i++) { */

/* 	if (LTop->Locus[markers[numbered[i]]].Class == TRAIT || */
/* 	    LTop->Locus[markers[numbered[i+1]]].Class == TRAIT) {	   */
/* 	  thetas[numbered[i]] = 0.5; */
/* 	  continue; */
/* 	}   */

/* 	if (LTop->Marker[markers[numbered[i]]].chromosome == */
/* 	   LTop->Marker[markers[numbered[i+1]]].chromosome) { */
/* 	  thetas[numbered[i]]= */
/* 	    (double) */
/* 	    (LTop->Marker[markers[numbered[i+1]]].pos_avg - */
/* 	     LTop->Marker[markers[numbered[i]]].pos_avg); */
/* 	  if (thetas[numbered[i]] >= 0.0) */
/* 	    thetas[numbered[i]] =  */
/* 	      ((LTop->map_distance_type=='h')? */
/* 	       haldane_theta(thetas[numbered[i]]) :   */
/* 	       kosambi_theta(thetas[numbered[i]])); */
/* 	  else { */
/* 	    if (LTop->Marker[markers[numbered[i]]].pos_avg >= 0.0 && */
/* 		LTop->Marker[markers[numbered[i+1]]].pos_avg >= 0.0) { */
/* 	      sprintf(err_msg,  */
/* 		      "%s and %s are not in order of increasing map distance!", */
/* 		      LTop->Marker[markers[numbered[i]]].MarkerName, */
/* 		      LTop->Marker[markers[numbered[i+1]]].MarkerName); */
/* 	      warnf(err_msg); */
/* 	      sprintf(err_msg, */
/* 		      "Their respective map positions are %7.4g and %7.4g respectively.",  */
/* 		      LTop->Marker[markers[numbered[i]]].pos_avg,  */
/* 		      LTop->Marker[markers[numbered[i+1]]].pos_avg); */
/* 	      warnf(err_msg); */
/* 	    } */
/* 	  } */
/* 	} */
/* 	else */
/* 	  thetas[numbered[i]] = 0.5; */
/*       } */
/*     } */
/*     tmpi=0; */
/*     if (thetas[0] < 0.0) { */
/*       thetas[0] = 0.5; /\* okay to keep the first one as 0.5 *\/ */
/*       fprintf(filep, " %7.5f", thetas[0]); */
/*       tmpi++; */
/*     } */
/*     for (i=tmpi; i < num_markers-1; i++) { */
/*       if (analysis == TO_GeneHunter) { */
/* 	thetas[i]=((thetas[i] < 0.0)? 0.0 : thetas[i]); */
/*       } */
/*       else { */
/* 	thetas[i]=((thetas[i] < 0.0)? 0.5 : thetas[i]); */
/*       } */
/*       fprintf(filep, " %7.5f", thetas[i]); */
/*     } */

/*     free(numbered); free(thetas); */

/*     if (!(analysis == TO_MERLIN || analysis == TO_MERLINONLY || */
/* 	  analysis == TO_PREST)) { */
/*       if (LTop->map_distance_type == 'h') */
/* 	fprintf(filep," Haldane"); */
/*       else */
/* 	fprintf(filep," Kosambi"); */
/*     } */
/*     break; */
/*   } */
/*   fprintf(filep, "\n");   */
/*   return; */
/* } */

/**
 * Print the recombination fractions of all of the markers on all of the chromosomes.
 *
 * @param filep the file descriptor to use as output.
 * @param num_markers the number of markers.
 * @param marders a vector of markers in the order that they are to be output.
 * @param LTop the base of the tree of structures in which the marker data is stored.
 * @param analysis the analysis (output format) being used.
 * @param tt the theta type which is either sex average, male or female.
 *
 * @return void
 */
void print_recomb_fracs(FILE *filep,
			int num_markers, int *markers,
			linkage_locus_top *LTop,
			analysis_type analysis, theta_type tt) {
    int i;
    double *thetas;

    // Need at least 2 markers to compute a delta...
    if (num_markers > 1) {
        thetas=CALLOC((size_t) num_markers-1, double);
        for (i=0; i < num_markers-1; i++) {
            if (LTop->Locus[markers[i]].Class == TRAIT ||
                LTop->Locus[markers[i+1]].Class == TRAIT) {
                // When two traits are located together.
                // This works even if they are on different chromosomes.
                thetas[i] = 0.5; // The null hypothesis between the trait and the first marker....
                continue;
            }

            if (LTop->Marker[markers[i]].chromosome ==
                LTop->Marker[markers[i+1]].chromosome) {
                // It's the same chromosome..
                double delta, delta_abs, marker = 0, marker_next = 0; //silly compiler
                // CPK: compute the delta between the appropriate markers...
                switch (tt) {
                case SEX_AVERAGED_THETA:
                    marker_next = LTop->Marker[markers[i+1]].pos_avg;
                    marker = LTop->Marker[markers[i]].pos_avg;
                    break;
                case MALE_THETA:
                    marker_next = LTop->Marker[markers[i+1]].pos_male;
                    marker = LTop->Marker[markers[i]].pos_male;
                    break;
                case FEMALE_THETA:
                    marker_next = LTop->Marker[markers[i+1]].pos_female;
                    marker = LTop->Marker[markers[i]].pos_female;
                    break;
                }
                delta = marker_next - marker;
                // CPK: use the absolute value of the delta...
                delta_abs = (delta >= 0.0) ? delta : -delta;
                // If the markers are real, are they in the right order?
                if (delta < 0.0 && marker >= 0.0 && marker_next >= 0.0) {
                    warnvf("%s (%.5g) and %s (%.5g) are not ordered by increasing map distance; using abs(distance)!\n",
                           LTop->Marker[markers[i]].MarkerName, marker,
                           LTop->Marker[markers[i+1]].MarkerName, marker_next);
                    //CPK: warnf("Setting the distance between these to zero cM.");
                }
	        // Both these conversion functions return 0 thetas if difference is negative.
                // CPK: However, now that they are getting the absolute value they will not...
                thetas[i] =
                    ((LTop->map_distance_type=='h')?
                     haldane_theta(delta_abs) :  kosambi_theta(delta_abs));
            } else {
	        // It's a new chromosome...
                thetas[i] = 0.5;
            }
        }

	// Print out the thetas after having coputed them.
	// Looks like we could have just printed them out in the first place without the vector...
        for (i=0; i < num_markers-1; i++) {
            fprintf(filep, " %7.5f", thetas[i]);
        }
        free(thetas);

        if (!(analysis == TO_MERLIN || analysis == TO_MERLINONLY || analysis == TO_PREST)) {
            if (LTop->map_distance_type == 'h')
                fprintf(filep," Haldane");
            else
                fprintf(filep," Kosambi");
        }
    }
    switch (tt) {
    case SEX_AVERAGED_THETA:
        break;
    case MALE_THETA:
        fprintf(filep," male"); 
        break;
    case FEMALE_THETA:
        fprintf(filep," female"); 
        break;
    }
    fprintf(filep, "\n");
}
