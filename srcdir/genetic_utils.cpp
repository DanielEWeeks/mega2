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

#include "error_messages_ext.h"
#include "grow_string_ext.h"
#include "utils_ext.h"
#include "user_input_ext.h"
/*
     error_messages_ext.h:  err_or_warn errorf mssgf my_calloc my_malloc warnf
        grow_string_ext.h:  grow
              utils_ext.h:  EXIT draw_line
*/



/* prototypes */
#ifdef USE_PROTOS

double           kinship(ped_tree *Ped, ped_rec *a, ped_rec *b);
void            create_marker_sublist(linkage_ped_top *Top, loci_markers *loci_markers,
				      int NumberOfMarkers);
double           haldane_theta(double x);
double           haldane_x(double t);
double           kosambi_theta(double x);
double           kosambi_x(double t);
double           kosambi_to_haldane(double f);
double           haldane_to_kosambi(double f);
void             switch_map(linkage_locus_top *LTop, double *position, sex_map_types sex);
old_new        *add_new_to_old(int oldped, int newped, int oldmem, int newmem, old_new *old_new_list);
void            help_sib_map(void);
int             purge_unneeded_pedigrees(ped_top *Top, int *save);
int             genoindx(int i, int j);
int             **create_genos(int num_alleles);
void            draw_histogram( FILE *filep, int num_alleles,
				allele_freq_struct *ar, const char *messg,
				const char *sex);

void             check_map_positions(linkage_ped_top *LTop, int numchr,
				     int check_or_set);
const char      *loc_type_name(linkage_locus_type Type);

#endif

/******************* Generic Genetics Routines *******************/

void log_marker_selections(linkage_ped_top *Top, int *entries,
			   int num_entries, char *asterisks,
			   void (*log_func)(const char *messg), int disp_err)
{
    int m, k, entry_count, *number1;
    linkage_locus_top *LocusTop = Top->LocusTop;

    if (entries == NULL) {
        /* ignore num_entries, just display all loci */
        number1 = CALLOC((size_t) LocusTop->LocusCnt, int);
        entry_count = LocusTop->LocusCnt;
        for (k=0; k < entry_count; k++) number1[k]=k;
    } else {
        number1 = CALLOC((size_t) num_entries, int);
        entry_count = num_entries;
        for (k=0; k < entry_count; k++) number1[k]=entries[k];
    }

    if (disp_err ==0) {
        strcpy(err_msg, "Number  Locus name      Locus Type  ");

    } else {
        strcpy(err_msg, "Number  Locus name      Locus Type   Error");
    }

    if (Top->EXLTop != NULL) {
#ifdef USEOLDMAPCODE
        for (m=0; m < Top->EXLTop->MapCnt; m++) {
            strcat(err_msg, Top->EXLTop->MapNames[m]);
        }
    } else {
#endif /* USEOLDMAPCODE */
        strcat(err_msg, " Position");
    }

    (*log_func)(err_msg);

    for (k = 0; k < entry_count; k++)   {
/*     if (LocusTop->Marker[number1[k]].pos_avg >= 0.0 ||  */
/* 	LocusTop->Locus[number1[k]].Type == QUANT || */
/* 	LocusTop->Locus[number1[k]].Type == AFFECTION) { */
        strcpy(err_msg, "");
        if (asterisks != NULL) {
            sprintf(err_msg, "%c%-6d %-15s %-9s",
                    asterisks[k], k+1,
                    LocusTop->Locus[number1[k]].LocusName,
                    loc_type_name(LocusTop->Locus[number1[k]].Type));
        } else {
            sprintf(err_msg, " %-6d %-15s %-9s",
                    k+1,
                    LocusTop->Locus[number1[k]].LocusName,
                    loc_type_name(LocusTop->Locus[number1[k]].Type));
        }

        if (disp_err) {
            if (LocusTop->Locus[number1[k]].Type == NUMBERED &&
                LocusTop->Marker[number1[k]].error_prob > 0.0) {
                grow(err_msg, " %7.5f",
                     LocusTop->Marker[number1[k]].error_prob);
            } else {
                strcat(err_msg, "    -    ");
            }
        }

        /* Now print the positions */
        if (Top->EXLTop != NULL) {
            double pos;
            for (m=0; m < Top->EXLTop->MapCnt; m++) {
#ifdef USEOLDMAPCODE
                if (LocusTop->Locus[number1[k]].Class != MARKER ||
                    LocusTop->Marker[number1[k]].chromosome == MALE_CHROMOSOME) {
                    // If this is not a marker, or it's the Y-chromosome then it's missing...
                    pos = UNKNOWN_POSITION;
                } else if (LocusTop->Marker[number1[k]].chromosome == SEX_CHROMOSOME) {
                    // If no female position, choose the average on the X-chromosome...
                    pos = ((Top->EXLTop->EXLocus[number1[k]].pos_female[m] >= 0.0)?
                           Top->EXLTop->EXLocus[number1[k]].pos_female[m] :
                           Top->EXLTop->EXLocus[number1[k]].positions[m]);
                } else {
                    // Here it's an autosome so we use the sex-averaged position...
                    pos = Top->EXLTop->EXLocus[number1[k]].positions[m];
                }
#else /* USEOLDMAPCODE */
                pos = UNKNOWN_POSITION;
                if (LocusTop->Locus[number1[k]].Class == MARKER) {
                    // If it's a marker...
                    if (LocusTop->Marker[number1[k]].chromosome == SEX_CHROMOSOME &&
                        Top->EXLTop->EXLocus[number1[k]].pos_female[m] >= 0.0 &&
                        (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
                         genetic_distance_sex_type_map == FEMALE_GDMT)) {
                        // If it's a female chromosome, a female map position exists, and the user has specified a usable map...
                        pos = Top->EXLTop->EXLocus[number1[k]].pos_female[m];
                    } else if (LocusTop->Marker[number1[k]].chromosome != MALE_CHROMOSOME &&
                               Top->EXLTop->EXLocus[number1[k]].positions[m] >= 0.0 &&
                               genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                        // If it's not male (e.g., autosome, mitocondrial), a position exists, and the user has specified a sex-averaged map...
                        pos = Top->EXLTop->EXLocus[number1[k]].positions[m];
                    }
                }
#endif /* USEOLDMAPCODE */
                if (pos >= 0.0) {
                    grow(err_msg, " %10.7f", pos);
                } else {
                    strcat(err_msg, "     -     ");
                }
            }
#ifdef USEOLDMAPCODE
	// In the new map code we will never get here, because there will always be an EXLTop after reading a map file...
        } else {
            if (LocusTop->Marker[number1[k]].pos_avg < 0.0 ||
                LocusTop->Marker[number1[k]].chromosome == UNKNOWN_CHROMO) {
                strcat(err_msg, "     -     ");
            } else {
                grow(err_msg, " %10.7f",
                     LocusTop->Marker[number1[k]].pos_avg);
            }
#endif /* USEOLDMAPCODE */
        }
        (*log_func)(err_msg);
        /* } */
    }

    (*log_func)(" ");
    free(number1);
}

/* index into the list of genotypes, which is a flat list as long as
   the number of genotypes */

int genoindx(int i, int j)

{
    int indx;
    if (i>=j) indx=(i*(i+1))/2 -(i-j);
    else
        indx=(j*(j+1))/2 -(j-i);

    return(indx-1);
}

/* create genotypes list from num_alleles,
   assume alleles are numbered from 1 to num_alleles */

int **create_genos(int num_alleles)

{
    int **geno, num_genos, indx;
    int i, j, i1, j1;

    num_genos = NUMGENOS(num_alleles); geno=NULL;
    /* put these two lines in a separate piece of code -allocmatrix */
    geno = (int **) CALLOC_PTR((size_t) num_genos, int *);
    for(i=0; i<num_genos; i++)
        geno[i] = CALLOC((size_t) 2, int);
    /*------------*/
    for (i=0; i<num_alleles; i++){
        for (j=i; j<num_alleles; j++){
            i1=i+1; j1=j+1;
            indx=genoindx(i1, j1);
            geno[indx][0]=i1; geno[indx][1]=j1;
        }
    }
    return geno;
}

/* Compute the kinship coefficient between two individuals a and b */

double      kinship(ped_tree *Ped, ped_rec *a, ped_rec *b)
{
    double      k;
    if (a == NULL || b == NULL)
        return 0.0;
    else if (a->ID == b->ID)
        k = (1.0 + kinship(Ped, a->Father, a->Mother))/2.0;
    else  {
        if (a->ID > b->ID) {
            k = (kinship(Ped, a->Father, b) + kinship(Ped, a->Mother, b))/2.0;
        } else {
            k = (kinship(Ped, b->Father, a) + kinship(Ped, b->Mother, a))/2.0;
        }
    }
    return k;
}

/* End of kinship */

void            create_marker_sublist(linkage_ped_top *Top, loci_markers *markerset,
				      int NumberOfMarkers)
{
    int            i, j;
    int            locus_cnt;

    /* marker sublist matrix only aplies to numbered or binary loci,
       and this should apply to pedtree NOT linkage_ped_tree,
       thus it is okay to simply have consecutive marker positions
       starting from zero */

    locus_cnt=0;
    markerset->NumberInSubOrder = NumberOfMarkers;
    for (i=0; i < Top->LocusTop->LocusCnt; i++)
        if (Top->LocusTop->Locus[i].Type == NUMBERED ||
            Top->LocusTop->Locus[i].Type == BINARY)
            locus_cnt++;

/*    for (i=0,j=0; i < Top->LocusTop->LocusCnt; i++) { */
/*      if (Top->LocusTop->Locus[i].Type == NUMBERED) { */
/*        marker_loci[j]=i; j++; */
/*      } */
/*    }  */

    markerset->total_files=locus_cnt - (NumberOfMarkers-1);
    markerset->MarkerSubsetMatrix = CALLOC((size_t) markerset->total_files, int*);

    for (i = 0; i < markerset->total_files; i++)
        markerset->MarkerSubsetMatrix[i] = CALLOC((size_t) NumberOfMarkers, int);

    for (i = 0; i < markerset->total_files; i++)
        for (j = 0; j < NumberOfMarkers; j++)
            markerset->MarkerSubsetMatrix[i][j] = i+j;
}

/* convert map distance(cM) to recombination fraction by haldane map function */
double           haldane_theta(double x)
{
    if (x <= 0)  return 0.0;
    x = x / 100.0; // convert cM to M
    return(0.5 * (1.0 - exp(-2.0*x)));
}

/* convert recombination fraction to map distance(cM) by haldane map function */
double          haldane_x(double t)
{
    return(100.0*(-0.5 * log(1.0 - 2.0*t)));
}

/* convert map distance(cM) to recombination fraction by kosambi map function */
double           kosambi_theta(double x)
{
    if (x <= 0)  return 0.0;
    x=x/100.0;
    return(0.5*((exp(4.0*x) - 1.0)/(exp(4.0*x) + 1.0)));
}

/* convert recombination fraction to map distance(cM) by kosambi map function */
double           kosambi_x(double t)
{
    return(100*0.25 * (log((1.0+2.0*t)/(1.0 - 2.0*t))));
}

/* convert a haldane map distance to kosambi map distance */
double           haldane_to_kosambi(double f)
{
    double f1;
    if (f <= 0.0) {
        f1=0.0;
    } else {
        f1=haldane_theta(f);
        f1=kosambi_x(f1);
    }
    return(f1);
}

/* convert a kosambi map distance to haldane map distance */
double           kosambi_to_haldane(double f)
{
    double f1;
    if (f <= 0.0) {
        f1=0.0;
    } else {
        f1=kosambi_theta(f);
        f1=haldane_x(f1);
    }
    return(f1);
}

/*
  Convert all the map distances from 'LTop' on ChrLoci for chromosome numchr from Kosambi to Haldane
  or vice versa and return this in the parameter 'position' for the type of sex map specified in 'sex'.
  Uses the global NumChrLoci.
*/

void            switch_map(linkage_locus_top *LTop, double *position, sex_map_types sex)
{
    int l, pl;

    /* lastpos is the last cumulative map position */
    double lastpos = 0.0;
    /* ppos is the previous marker's map position */
    double ppos = 0.0;

    for (l=0, pl=-1; l < NumChrLoci; l++) {

        // Only if this is a marker...
        if (LTop->Locus[ChrLoci[l]].Type == NUMBERED ||
            LTop->Locus[ChrLoci[l]].Type == BINARY
	    //            LTop->Locus[ChrLoci[l]].Type == XLINKED ||
	    //            LTop->Locus[ChrLoci[l]].Type == YLINKED
	    ) {

            switch(sex) {

            case SEX_AVERAGED_MAP :
	        // Don't make this check on the first time through...
	        // Warn the user if the order is not increasing...
                if (((LTop->Marker[ChrLoci[l]].pos_avg - lastpos) < 0.0) &&
                    (pl >= 0)) {
                    warnvf("%s (%.5g) and %s (%.5g) are not ordered by increasing average map distance!\n",
                            LTop->Marker[ChrLoci[pl]].MarkerName, LTop->Marker[ChrLoci[pl]].pos_avg,
                            LTop->Marker[ChrLoci[l]].MarkerName, LTop->Marker[ChrLoci[l]].pos_avg);
                    warnf("Assigning a 0 inter-marker distance, map conversion may be inaccurate.");
                }
                position[l] = ppos +
                    ((LTop->map_distance_type == 'h')?
                     haldane_to_kosambi(LTop->Marker[ChrLoci[l]].pos_avg - lastpos) :
                     kosambi_to_haldane(LTop->Marker[ChrLoci[l]].pos_avg - lastpos));
                lastpos = LTop->Marker[ChrLoci[l]].pos_avg;
                break;

            case MALE_SEX_MAP :
                if (((LTop->Marker[ChrLoci[l]].pos_male - lastpos) < 0.0) &&
                    (pl >= 0)) {
                    warnvf("%s (%.5g) and %s (%.5g) are not ordered by increasing male map distance!\n",
                            LTop->Marker[ChrLoci[pl]].MarkerName, LTop->Marker[ChrLoci[pl]].pos_male,
                            LTop->Marker[ChrLoci[l]].MarkerName, LTop->Marker[ChrLoci[l]].pos_male);
                    warnf("Assigning a 0 inter-marker distance, male map conversion may be inaccurate.");
                }
                position[l] = ppos +
                    ((LTop->map_distance_type == 'h')?
                     haldane_to_kosambi(LTop->Marker[ChrLoci[l]].pos_male - lastpos) :
                     kosambi_to_haldane(LTop->Marker[ChrLoci[l]].pos_male - lastpos));
                lastpos = LTop->Marker[ChrLoci[l]].pos_male;
                break;

            case FEMALE_SEX_MAP :
                if (((LTop->Marker[ChrLoci[l]].pos_female - lastpos) < 0.0) &&
                    (pl >= 0)) {
                    warnvf("%s (%.5g) and %s (%.5g) are not ordered by increasing female map distance!\n",
                            LTop->Marker[ChrLoci[pl]].MarkerName, LTop->Marker[ChrLoci[pl]].pos_female,
                            LTop->Marker[ChrLoci[l]].MarkerName, LTop->Marker[ChrLoci[l]].pos_female);
                    warnf("Assigning a 0 inter-marker distance, female map conversion may be inaccurate.");
                }
                position[l] = ppos +
                    ((LTop->map_distance_type == 'h')?
                     haldane_to_kosambi(LTop->Marker[ChrLoci[l]].pos_female - lastpos) :
                     kosambi_to_haldane(LTop->Marker[ChrLoci[l]].pos_female - lastpos));
                lastpos = LTop->Marker[ChrLoci[l]].pos_female;
                break;

            } // switch(sex) {
            ppos = position[l];
            pl = l;

        } else {
	    // If this is not a marker...
            position[l] = UNKNOWN_POSITION;
        }
    } // for (l=0; l < NumChrLoci; l++) {
}

old_new        *add_new_to_old(int oldped, int newped, int oldmem, int newmem, old_new *old_new_list)

{
    old_new        *newone;
    newone = MALLOC(old_new);

    newone->oldped = oldped;
    newone->newped = newped;
    newone->oldmem = oldmem;
    newone->newmem = newmem;
    newone->next = old_new_list;
    return newone;
}

/*  old_new        *find_old(int newped, int newmem, old_new *old_new_list) */

/*  {    */
/*     old_new        *temp;    */
/*     temp = old_new_list;  */
/*     while (temp != NULL) */
/*     { */
/*     if ((newmem == temp->newmem) && (newped == temp->newped)) */
/*     { */
/*     return temp; */
/*     } */
/*     temp = temp->next; */
/*     } */
/*     return NULL; */
/*  } */


void            help_sib_map(void)
{
    int             i;

    const char           *(tcl_map_help_support[]) =
        {
            " ",
            "'support'-specifies the LOD score cutoff for support",
            "      intervals. The default is 1.0 log units, which",
            "      means that support intervals will cover the range",
            "      of distances that give LOD scores within 1.0 of the",
            "      LOD score of the most likely distance.",
            NULL
        };

    const char           *(tcl_map_help_epsilon[]) =
        {
            " ",
            "'EPSILON' - specifies the convergence criterion for the",
            "            iterative distance refinement routines.",
            "            The default is 0.00001 Morgans.",
            NULL
        };

    const char           *(tcl_map_help_do_shuffle[]) =
        {
            " ",
            "'do_shuffle ' - a boolean value: specifies if sib_map should",
            "                perform the map shuffling function, to verify",
            "                map order and/or place new markers on an existing",
            "                map. The default is FALSE. ",
            NULL
        };


    draw_line();
    printf("HELP: sib_map parameters...\n");
    printf("sib_map [-v] [-c cmd] [-f file][marker data...] [>map_data]\n\n");
    printf("the following TCL commands can specified using the \n");
    printf("-c or -f mechanisms...\n");
    for (i = 0; tcl_map_help_support[i] != NULL; i++)
        printf("%s\n", tcl_map_help_support[i]);
    printf("\n");
    for (i = 0; tcl_map_help_epsilon[i] != NULL; i++)
        printf("%s\n", tcl_map_help_epsilon[i]);
    printf("\n");
    for (i = 0; tcl_map_help_do_shuffle[i] != NULL; i++)
        printf("%s\n", tcl_map_help_do_shuffle[i]);
    draw_line();
    return;
}


/*---------------------------------------------------------------+
  | Check which pedigrees have at least two affecteds, and mark   |
  | them to be saved. Allocate a new array and copy the valid     |
  | pedigrees over. Free the old array, assign the new one, and   |
  | update PedCnt.                                                |
  |                                                               |
  | Return the number of pedigrees deleted.                       |
  +---------------------------------------------------------------*/
int             purge_unneeded_pedigrees(ped_top *Top, int *save)
{

    int             c=0, ped, newpedcnt, peds_killed;
    int             *savep = &save[0];
    /* Check the number of valid pedigrees and record them. */
    newpedcnt = 0;
    peds_killed = 0;

    for (ped = 0; ped < Top->PedCnt; ped++) {
        if (Top->PedTree[ped].AffectedCnt < 2)  {
            *savep = 0;
            peds_killed++;
        } else  {
            *savep = 1;
            newpedcnt++;
        }
        savep++;
    }

    if (peds_killed > 0) {
        SECTION_ERR_INIT(purge_ped);
        SECTION_ERR(purge_ped);
        warnvf("Not saving untyped pedigrees: \n");
        SECTION_ERR_HEADER(purge_ped);
        sprintf(err_msg, " ");
        for(ped=0; ped < Top->PedCnt; ped++) {
            if (save[ped] == 0) {
                strcat(err_msg, Top->PedTree[ped].Name); strcat(err_msg, " ");
                c += (int) strlen(Top->PedTree[ped].Name) + 1;
                if (c >= 70) {
                    SECTION_ERR(purge_ped);
                    errf(err_msg, "", 1);
                    sprintf(err_msg, " ");
                    c=0;
                }
            }
        }
        if (c > 0) {
            SECTION_ERR(purge_ped);
            errf(err_msg, "", 1);
        }
        SECTION_ERR_FINI(purge_ped);
    }
    return peds_killed;

}


void draw_histogram( FILE *filep, int num_alleles, allele_freq_struct *ar,
                     const char *messg, const char *sex)

{
    int i,j, num_asterix;
    double *freq;
    int *counts;


    if (!(strcmp(sex, "male"))) {
        freq = &(ar->male_allele_freq[0]);
        counts = &(ar->Male_Allele_Bin[0]);
    } else if (!(strcmp(sex, "female"))) {
        freq = &(ar->female_allele_freq[0]);
        counts = &(ar->Female_Allele_Bin[0]);
    } else {
        freq = &(ar->allele_freq[0]);
        counts = &(ar->Allele_Bin[0]);
    }
    fprintf(filep, "\n%s:\n", messg);

    fprintf(filep, "Allele#    Input Freq.   Obs. Freq.  Counts  Histogram\n");

    fprintf(filep, "-------------------------------------------------------\n");

    for (i=0; i<num_alleles; i++) {
        fprintf(filep, "%6d       %7.6f     %7.6f    %5d  ", i+1,
                ar->input_freq[i], freq[i],	counts[i]);

        num_asterix=(int)(50.0 * freq[i]);
        for (j=0; j<num_asterix; j++)
            fprintf(filep, "*");
        fprintf(filep, "\n");
    }
    fprintf(filep, "-------------------------------------------------------\n");
    return;
}

/* check if genetic distances in the selected map order is the
   correct. If check_or_set  = 1, then only check,  otherwise
   warn about setting */

void   check_map_positions(linkage_ped_top *Top, int numchr,
			   int check_or_set)

{
    char cselect[5];
    int i, found_neg=0, done=0, opt=0;
    double diff;

    int loc_cnt, *markers;
    linkage_locus_top  *LTop= Top->LocusTop;

    markers=CALLOC((size_t) LTop->LocusCnt, int);
    loc_cnt=0;
    for (i = 0; i < LTop->LocusCnt - 1; i++)  {
        if (LTop->Locus[i].Type == NUMBERED) {
            markers[loc_cnt]=i;
            loc_cnt++;
        }
    }

    for(i=0; i <  loc_cnt-1; i++) {
        if (LTop->Marker[markers[i + 1]].chromosome != LTop->Marker[markers[i]].chromosome)
            continue;
        if ((diff =
             LTop->Marker[markers[i + 1]].pos_avg - LTop->Marker[markers[i]].pos_avg) < 0.0) {
            if (check_or_set == 1) {
                sprintf(err_msg,
                        "Selected marker loci on chromosome %d are not in map order.",
                        numchr);
                warnf(err_msg);
                sprintf(err_msg, "   Locus %d (%s) is at position %10.7f",
                        markers[i]+1, LTop->Marker[markers[i]].MarkerName,
                        LTop->Marker[markers[i]].pos_avg);
                warnf(err_msg);
                sprintf(err_msg, "   Locus %d (%s) is at position %10.7f",
                        markers[i+1]+1,
                        LTop->Marker[markers[i + 1]].MarkerName,
                        LTop->Marker[markers[i + 1]].pos_avg);
                warnf(err_msg);
                found_neg=1;
            } else {
                Display_Errors=1;
                sprintf(err_msg,
                        "Recombination fraction between %s and %s set to 0.0999.",
                        LTop->Marker[markers[i]].MarkerName,
                        LTop->Marker[markers[i+1]].MarkerName);
                warnf(err_msg);
            }
        }
    }

    if (check_or_set != 1)  return;

    if (found_neg) {
        while(!done) {
            printf("Selected markers not in map order - how shall we proceed?\n");
            printf("0) Exit Mega2\n");
            printf(
                " 1) Set inconsistent genetic distances to 0.0999\n");
            printf("Enter option 0 or 1 > ");
            fflush(stdout);
            IgnoreValue(fgets(cselect, 4, stdin)); newline;
            sscanf(cselect, "%d", &opt);
            if (opt==1) {
                done =1;
                mssgf("Inconsistent genetic distances will be set to 0.0999.");
                printf("See MEGA2.LOG for details.\n");
                SetMarkerPosToSpecial=1;
            } else if (opt==0) {
                errorf("Terminating Mega2 because markers are not in map order.");
                EXIT(OUTPUT_FORMAT_ERROR);
            } else {
                warn_unknown(cselect);
                done=0;
            }
        }
    } else {
        SetMarkerPosToSpecial=0;
    }
}


char *chrom_num_to_name(const int chrom_num, char *chrom_name)

{
    if (chrom_num ==  UNKNOWN_CHROMO) {
        strcpy(chrom_name, "U");

    } else if (chrom_num ==  SEX_CHROMOSOME) {
        if (human_x) {
            strcpy(chrom_name, "X");
        }

    } else if (chrom_num == PSEUDO_X) {
        if (human_xy) {
            strcpy(chrom_name, "XY");
        }
  
    } else if (chrom_num ==  MALE_CHROMOSOME) {
        if (human_y) {
            strcpy(chrom_name, "Y");
        }

    } else if (chrom_num == MITO_CHROMOSOME) {
        if (human_mt) {
            strcpy(chrom_name, "MT");
        }

    } else {
        sprintf(chrom_name, "%d", chrom_num);
    }

    return &(chrom_name[0]);
}

const char *loc_type_name(linkage_locus_type Type)

{

    switch(Type) {
    case NUMBERED:
        return("Numbered");
        break;
    case BINARY:
        return("Binary");
        break;
    case XLINKED:
        return("Xlinked");
        break;
    case YLINKED:
        return("Ylinked");
        break;
    case AFFECTION:
        return("Affection");
        break;
    case QUANT:
        return("Quant");
        break;
    default:
        return("Unknown");
        break;
    }
}

/*************** End Of Generic Genetics Routines ****************/
