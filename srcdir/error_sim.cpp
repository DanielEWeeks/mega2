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

/* error_sim.c
   routines to generate errors in genotypes according to a
   specified model
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "error_sim.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "genetic_utils_ext.h"
#include "grow_string_ext.h"
#include "list_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "reorder_loci_ext.h"
#include "utils_ext.h"
/*
        batch_input_ext.h:  batchf invalid_value_field
     error_messages_ext.h:  errorf errsimf mssgf my_calloc warnf
      genetic_utils_ext.h:  log_marker_selections
        grow_string_ext.h:  grow
               list_ext.h:  append_to_list_tail new_list pop_first_list_entry
  output_file_names_ext.h:  change_output_chr file_status
    output_routines_ext.h:  field_widths
       reorder_loci_ext.h:  display_selections
              utils_ext.h:  EXIT draw_line imax randomnum
*/


int ***GenotypeIndexes;
/* for indexing and reverse-indexing genotypes
   for all markers */

double UniformErrProb;
double SW2ErrProb[5];
int MarkerErrorProb; /* flag to specify if map file contains errors */

/* prototype definitions */

static void error_at_locus(int locnum, int loc_index, linkage_ped_top *Top,
			   genotype (*err_model)(genotype gen, int number,
                                                 linkage_locus_rec Loc,
                                                 sw2_error_class *err_type), listhandle **err_list,
			   int *typed, int *numhomo, int *numhetero);
genotype uniform_error_model(genotype gen, int number,
			     linkage_locus_rec Loc, sw2_error_class *err);
genotype marker_error_model(genotype gen, int number,
			    linkage_locus_rec Loc, sw2_error_class *err);
genotype sw2_error_model(genotype gen, int number,
			 linkage_locus_rec Loc, sw2_error_class *err);
static genotype select_this_geno(int index, int num_indices,
                                 int num_total_genos,
                                 genotype true_geno, int loc,
                                 sw2_error_class error_class);
static genotype match_index_to_genotype(int index,
					int loc,
					int num_genos);
static void gen_all_indexes(linkage_locus_top *LTop,
			    int num_selected,
			    int *selected);
static void gen_marker_indexes(int num_alleles,
			       int **indexes);
static int gen_index(genotype gen, int num_alleles);
static int ***alloc_gen_indexes(int num_selected,
				int *selected_loci,
				linkage_locus_top *LTop);
static void free_gen_indexes(linkage_locus_top *LTop,
			     int num_loci,
			     int *selected_loci,
			     int ****indexes_ptr);
static void error_sim_menu(linkage_ped_top *Top,
			   char *err_model, char *fl_names[],
			   int numchr);

static void init_error_probs(void);

static void output_err_files(int *loci, int num_loci, char *error_model,
			     listhandle **error_recs,
			     linkage_ped_top *Top, char *fl_names[],
			     int *NumTyped, int *NumHomo, int *NumHetero);

static void errsim_file_names(const char *type, char *fl_names[], int numchr);
/* export */
void simulate_errors(linkage_ped_top *Top, int numchr, char *fl_names[]);

/* generate errors at a given locus for all pedigrees */

static void error_at_locus(int locnum, int loc_index, linkage_ped_top *Top,
			   genotype (*err_model)(genotype gen, int number,
                                                 linkage_locus_rec Loc,
                                                 sw2_error_class *err_type),
			   listhandle **err_list, int *typed, int *homo,
			   int *hetero)

{

    int ped, entry;
    genotype gen, new_gen;
    int all1, all2, num_errors;
    err_sim_out *err_out_item;
    sw2_error_class err_class;

    errsimf(" ");
    sprintf(err_msg, "Introducing errors at locus %s",
            Top->LocusTop->Locus[locnum].Name);
    errsimf(err_msg);

    *err_list = new_list();

    *typed = 0;
    if (homo != NULL) *homo = 0;
    if (hetero != NULL) *hetero= 0;

    num_errors = 0;
    for(ped=0; ped < Top->PedCnt; ped++) {
        for(entry=0; entry < Top->Ped[ped].EntryCnt; entry++) {
            /* Then generate error at this locus */
            get_2alleles(Top->Ped[ped].Entry[entry].Marker, locnum, &all1, &all2);
            if (!all1 || !all2) {
                continue;
            }
            (*typed)++;
            if (homo != NULL) {
                if (all1 == all2) (*homo)++;
            }
            if (hetero != NULL) {
                if (all1 != all2) (*hetero)++;
            }

            /* always place the smaller allele first */
            if (all1 <= all2) {
                gen.all1 = all1; gen.all2 = all2;
            } else {
                gen.all2 = all1; gen.all1 = all2;
            }
            new_gen = (*err_model)(gen, loc_index,
                                   Top->LocusTop->Locus[locnum], &(err_class));
            if ((new_gen.all1 != gen.all1) ||
               (new_gen.all2 != gen.all2)) {
                /* store the error */
                num_errors++;
                err_out_item = CALLOC((size_t) 1, err_sim_out);
                err_out_item->ped = ped;
                err_out_item->entry = entry;
                err_out_item->original.all1 = gen.all1;
                err_out_item->original.all2 = gen.all2;
                err_out_item->new_genotype.all1 = new_gen.all1;
                err_out_item->new_genotype.all2 = new_gen.all2;
                err_out_item->err_type = err_class;
                append_to_list_tail(*err_list, (list_data)err_out_item);
		// The alleles will print out as numeric...
                sprintf(err_msg,
                        "(Ped %s, Person %s): %d, %d -> %d, %d",
                        Top->Ped[ped].Name,
                        Top->Ped[ped].Entry[entry].OrigID,
                        gen.all1, gen.all2, new_gen.all1, new_gen.all2);
                errsimf(err_msg);
                /* easier to test with a diff command */
                set_2alleles(Top->Ped[ped].Entry[entry].Marker, locnum, 
                             &Top->LocusTop->Locus[locnum], new_gen.all1, new_gen.all2);
            }
        }
    }

    sprintf(err_msg, "Total number of genotypes = %d", *typed);
    errsimf(err_msg);
    sprintf(err_msg, "Number of errors introduced = %d",
            num_errors);
    errsimf(err_msg);
    sprintf(err_msg, "Percentage error rate = %4.3f",
            (double)num_errors/(double)(*typed));
    errsimf(err_msg);

    return;
}

/* Functions to introduce error with a specific probability
   at a genotype */
/* We will use global parameters for the error models */

/* This uses a single parameter E = probability of error at this
   genotype, the specific genotype is determined by allele
   frequencies.
   First generate a single random number (0, 1.0)
   If this number is less than the error_probability value,
   then select one out of n-1 wrong genotypes.
   So, generate a number between 0 and n-1. This number is the
   index of the new genotype. If it is the same or larger than
   the correct genotype's index, we have to increase it by 1
   to index correctly into the sorted list of genotypes.
*/

genotype uniform_error_model(genotype gen, int number,
			     linkage_locus_rec Loc,
			     sw2_error_class *err_type)

{

    genotype new_gen;
    int new_ind, slct_index;
    double slct, select1;

    slct = randomnum();
    slct_index = gen_index(gen, Loc.AlleleCnt);

    if (slct < UniformErrProb) {
        *err_type = Uniform;
        /* Need to change genotype */
        select1 = randomnum();

        new_ind = (int)((double)(NUMGENOS(Loc.AlleleCnt) - 1) *
                        select1);
        if (new_ind >= slct_index) new_ind++;
        new_gen = match_index_to_genotype(new_ind, number,
                                          NUMGENOS(Loc.AlleleCnt));
        return new_gen;
    } else {
        return gen;
    }
}

genotype marker_error_model(genotype gen, int number,
			    linkage_locus_rec Loc,
			    sw2_error_class *err_type)

{

    genotype new_gen;
    double slct;
    int new_ind, slct_index = gen_index(gen, Loc.AlleleCnt);

    slct = randomnum(); (void) randomnum();
    if (slct < Loc.Marker->error_prob) {
        *err_type = Marker;
        /* Need to change genotype */
        new_ind = (int) ((double)((NUMGENOS(Loc.AlleleCnt))-1) *
                         randomnum());
        if (new_ind >= slct_index) new_ind++;

        new_gen =
            match_index_to_genotype(new_ind, number,
                                    (NUMGENOS(Loc.AlleleCnt)));
        return new_gen;
    } else {
        return gen;
    }
}

/* This is a complicated model involving 5 classes of errors
   associated with 5 error parameters.
   E1 <- heterozygote is read in as homozygote with one correct
   allele
   E2 <- heterozygote is read in as heterozygote with one
   wrong allele

   E3 <- homozygote is read in a homozygote with both alleles
   wrong
   or
   heterozygote is read in as heterozygote with both alleles
   incorrect
   E4 <- homozygote is read in as heterozygote with onw allele
   incorrect
   E5 <- heterozygote is read in as homozogote with both alleles
   inocrrect
   or
   homozygote is read in as heterozygote with both alleles
   incorrect

   Total error probability for homozygotes = E3+E4+E5
   Total error probability for heterozygotes = E1+E2+E3+E5
   Steps for generating errors:

   1. For every genotype, depending on whether it is homozygous or
   heterozygous, decide whether to introduce error by looking at the
   total error probabilities.

   2. Divide the [0,1] interval into sub-itervals by looking at the
   ratios of the component errors e.g. for homozygotes the
   interval will be divided into e3:e4:e5. Then select a random
   number in the [0,1] interval to decide which component to select.

   3. Sort the "wrong" genotypes belonging to each component, then
   choose a random integer between [1, K] where K is the number of
   wrong genotypes in each component.

*/



genotype sw2_error_model(genotype gen, int number,
			 linkage_locus_rec Loc,
			 sw2_error_class *err_class)

{

    /* Step 1: homozygote or heterozygote? */

    double slct;
    double select1;
    double select2;
    double total_prob, homo_limits[2], hetero_limits[3];
    int index, num_error_genos = 0, total_genos = NUMGENOS(Loc.AlleleCnt);

    genotype new_gen;


    slct = randomnum();
    select1 = randomnum();
    select2 = randomnum();

    if (gen.all1 == gen.all2) {
        /* homozygote */
        total_prob = SW2ErrProb[2]+SW2ErrProb[3]+
            SW2ErrProb[4];
        homo_limits[0] = SW2ErrProb[2]/total_prob;
        homo_limits[1] = homo_limits[0] + SW2ErrProb[3] / total_prob;

        if (slct < total_prob) {
            if (select1 < homo_limits[0]) {
                /* error of type E3 -
                   number of  possible genotypes = n-1 */

                num_error_genos = Loc.AlleleCnt - 1;
                *err_class = E3_homo;
            } else if (select1 >= homo_limits[0] && select1 < homo_limits[1]) {
                /* Error of type E4 -
                   number of possible genotypes = n-1 */

                num_error_genos = Loc.AlleleCnt - 1;
                *err_class = E4;
            } else if (select1 > homo_limits[1]) {
                /* Error of type E5 -
                   number of possible genotypes = (n-1)(n-2)/2 */

                num_error_genos = (Loc.AlleleCnt - 1)*(Loc.AlleleCnt - 2)/2;
                *err_class = E5_homo;
            }

            index = (int)(select2 * (double) num_error_genos);
            printf("I %d, NE %d, NT %d, G(%d,%d), C %d\n",
                   index, num_error_genos, total_genos,
                   gen.all1, gen.all2, *err_class);
            new_gen = select_this_geno(index, num_error_genos,
                                       total_genos, gen, number, *err_class);

            return new_gen;
        } else {
            return gen;
        }
    } else {

        /* heterozygote */
        total_prob = SW2ErrProb[0]+SW2ErrProb[1]+
            SW2ErrProb[2]+SW2ErrProb[4];

        hetero_limits[0] = SW2ErrProb[0] / total_prob;
        hetero_limits[1] = hetero_limits[0] + SW2ErrProb[1] / total_prob;
        hetero_limits[2] = hetero_limits[1] + SW2ErrProb[2] / total_prob;

        if (slct < total_prob) {
            if (select1 < hetero_limits[0]) {
                /* Error of type E1 -
                   number of possible genotypes = 2 */
                num_error_genos = 2;
                *err_class = E1;
            } else if (select1 >= hetero_limits[0] && select1 < hetero_limits[1]) {
                /* Error of type E2 -
                   number of possible genotypes = 2(n-2) */
                num_error_genos = 2*(Loc.AlleleCnt -  2);
                *err_class = E2;

            } else if (select1 >= hetero_limits[1] && select1 < hetero_limits[2]) {
                /* Error of type E3 -
                   number of possible genotypes = (n - 2)(n-3)/2 */

                num_error_genos = (Loc.AlleleCnt - 2)*(Loc.AlleleCnt - 3)/2;
                *err_class = E3_hetero;
            } else if (select1 > hetero_limits[2]) {
                /* Error of type E5 -

                   number of possible genotypes = n - 2 */
                num_error_genos = Loc.AlleleCnt - 2;
                *err_class = E5_hetero;
            }
            index = (int) (select2 * (double) num_error_genos);
            printf("I %d, NE %d, NT %d, G(%d,%d), C %d\n",
                   index, num_error_genos, total_genos,
                   gen.all1, gen.all2, *err_class);
            new_gen = select_this_geno(index, num_error_genos,
                                       total_genos, gen, number, *err_class);
            return new_gen;

        } else {
            return gen;
        }

    }
}


static genotype select_this_geno(int index, int num_indices,
				 int num_total_genos,
				 genotype true_geno, int loc,
				 sw2_error_class error_class)

{

    int i=0, num;
    genotype obs_geno;
    boolean matches_error;
    obs_geno.all1 = obs_geno.all2 = 0;

    for (num=0; num < num_total_genos; num++) {
        /* Go through the genotypes that satisfy
           the conditions until we have counted up to index
           -  condition is that genotype is homozygotic
           and not equal to true_genotype */

        switch(error_class) {
        case E1:
            if ((GenotypeIndexes[loc][num][1] == GenotypeIndexes[loc][num][2]) &&
               ((GenotypeIndexes[loc][num][1] == true_geno.all1) ||
                (GenotypeIndexes[loc][num][1] == true_geno.all2))) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E2:
            if ((GenotypeIndexes[loc][num][1] != GenotypeIndexes[loc][num][2]) &&
               ((GenotypeIndexes[loc][num][1] == true_geno.all1)  ||
                (GenotypeIndexes[loc][num][2] == true_geno.all1)  ||
                (GenotypeIndexes[loc][num][1] == true_geno.all2)  ||
                (GenotypeIndexes[loc][num][2] == true_geno.all2))) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E3_homo:
            if ((GenotypeIndexes[loc][num][1] == GenotypeIndexes[loc][num][2]) &&
               (GenotypeIndexes[loc][num][1] != true_geno.all1)) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E3_hetero:
            if ((GenotypeIndexes[loc][num][1] != GenotypeIndexes[loc][num][2]) &&
               ((GenotypeIndexes[loc][num][1] != true_geno.all1) &&
                (GenotypeIndexes[loc][num][1] != true_geno.all2) &&
                (GenotypeIndexes[loc][num][2] != true_geno.all1) &&
                (GenotypeIndexes[loc][num][2] != true_geno.all2))) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E4:
            if ((GenotypeIndexes[loc][num][1] != GenotypeIndexes[loc][num][2]) &&
               ((GenotypeIndexes[loc][num][1] == true_geno.all1) ||
                (GenotypeIndexes[loc][num][2] == true_geno.all1))) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E5_homo:
            if ((GenotypeIndexes[loc][num][1] != GenotypeIndexes[loc][num][2]) &&
               ((GenotypeIndexes[loc][num][1] != true_geno.all1) &&
                (GenotypeIndexes[loc][num][2] != true_geno.all1))) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        case E5_hetero:
            if ((GenotypeIndexes[loc][num][1] == GenotypeIndexes[loc][num][2]) &&
               (GenotypeIndexes[loc][num][1] != true_geno.all1)) {
                matches_error = true;
            } else {
                matches_error = false;
            }
            break;

        default:
            matches_error = false;
            break;
        }
        printf("loc %d, num %d\n", loc, num);
        if (matches_error == true) {
            if (i == index) {
                obs_geno.all1 = GenotypeIndexes[loc][num][1];
                obs_geno.all2 = GenotypeIndexes[loc][num][2];
                break;
            } else {
                i++;
                if (i > num_indices) {
                    errorf("Number of valid genotypes exceeded!\n");
                    EXIT(OUTPUT_FORMAT_ERROR);
                }
            }
        }
    }

    return obs_geno;
}


static genotype match_index_to_genotype(int index, int loc,
					int num_genos)

{
    int geno;
    genotype genot;

    for(geno=0; geno < num_genos; geno++) {
        if (index == GenotypeIndexes[loc][geno][0]) {
            genot.all1 = GenotypeIndexes[loc][geno][1];
            genot.all2 = GenotypeIndexes[loc][geno][2];
            break;
        }
    }
    return genot;
}

/* generate indexes for all genotypes */

static void gen_all_indexes(linkage_locus_top *LTop, int num_selected,
			    int *selected)

{

    int loc;

    GenotypeIndexes = alloc_gen_indexes(num_selected,
                                        selected, LTop);

    for (loc = 0; loc < num_selected; loc++) {
        gen_marker_indexes(LTop->Locus[selected[loc]].AlleleCnt,
                           GenotypeIndexes[loc]);
    }
    return;
}

static void gen_marker_indexes(int num_alleles, int **indexes)

{
    int a1, a2;
    int i;
    genotype gen1;

    /* create reverse indexes for each index */
    for (a1 = 1; a1 <= num_alleles; a1++) {
        for (a2 = a1; a2 <= num_alleles; a2++) {
            gen1.all1 = a1;  gen1.all2 = a2;
            i = gen_index(gen1, num_alleles);
            /* store the index */
            indexes[i][0]=i;
            indexes[i][1]=gen1.all1;
            indexes[i][2]=gen1.all2;
        }
    }
    return;
}


static int gen_index(genotype gen1, int num_alleles)

{
    int row, col, ind;

    genotype gen; /* order the alleles */

    /* index into an nxn matrix, numgenos = nxn */

    gen.all1 = ((gen1.all1 <= gen1.all2)? gen1.all1: gen1.all2);
    gen.all2 = ((gen1.all1 > gen1.all2)? gen1.all1: gen1.all2);

    row = gen.all1; /* [1 .. num_alleles] */
    col = gen.all2; /* [1 .. num_alleles] */

    /* To uniquely index into an upper diagonal matrix,
       where row r has (num_alleles - r+1) elements,
       ind = (n-1)*(r-1) - (r-1)*(r-2)/2 + c,
       numgenos = n(n+1)/2
    */

    ind = (num_alleles-1)*(row-1) - (row-1)*(row-2)/2 + col - 1;
    return ind;

}


/* selected_loci are indexes into the
   LTop->Locus[] array. We only store as many allele_index matrices
   as necessary.
*/

static int ***alloc_gen_indexes(int num_selected, int *selected_loci,
				linkage_locus_top *LTop)

{


    int loc, geno;
    int ***indexes;

    /* allocate pointers for all markers */

    indexes = CALLOC((size_t) num_selected, int **);

    /* Now allocate space for the selected loci */
    for (loc = 0; loc < num_selected; loc++) {
        indexes[loc] =
            CALLOC((size_t) NUMGENOS(LTop->Locus[selected_loci[loc]].AlleleCnt), int *);
        for (geno = 0;
             geno < NUMGENOS(LTop->Locus[selected_loci[loc]].AlleleCnt);
             geno++) {
            indexes[loc][geno] = CALLOC((size_t) 3, int);
        }
    }
    return indexes;
}

static void free_gen_indexes(linkage_locus_top *LTop,
			     int num_loci, int *selected_loci,
			     int ****indexes_ptr)
{

    int loc, geno;
    int ***indexes = *indexes_ptr;

    for (loc = 0; loc < num_loci; loc++) {

        for (geno = 0;
             geno < NUMGENOS(LTop->Locus[selected_loci[loc]].AlleleCnt);
             geno++) {

            free(indexes[loc][geno]);
        }
        free(indexes[loc]);
    }

    free(indexes);
    indexes = NULL;

    return;
}

/* Check that the marker numbers are valid, and that
   only marker loci have been chosen.
   isindex indiciates whether selections are
   locus item numbers from the interactive selection.
*/

static int check_loc_selection(linkage_locus_top *LTop,
			       int num_select,
			       int *loc_select, int is_index)

{
    int i, j;
    int found;

    if (is_index) {
        for (i=0; i < num_select; i++) {
            if (loc_select[i] < 1 || loc_select[i] > num_reordered) {
                printf("ERROR: Invalid selection %d.\n", loc_select[i]);
                return 0;
            }
        }
    } else {
        for (i=0; i < num_select; i++) {
            /* check if this locus has been selected in reordering */
            found=0;
            for (j=0; j < num_reordered; j++) {
                if (reordered_marker_loci[j] == loc_select[i]-1) {
                    found=1;
                    loc_select[i] = j+1;
                    break;
                }
            }
            if (found == 0) {
                printf("ERROR: Locus number %d (%s) not in reordered loci.\n",
                       loc_select[i], LTop->Locus[loc_select[i]-1].Name);
                return 0;
            }
        }
    }

    /* No errors */
    return num_select;

}

static void init_error_probs(void)

{

    UniformErrProb = 0.05;
    SW2ErrProb[0]=SW2ErrProb[1]=SW2ErrProb[2]=SW2ErrProb[3]=SW2ErrProb[4]=0.01;
    return;
}

static void errsim_file_names(const char *type, char *file_names[], int numchr)

{

    char filename_str[FILENAME_LENGTH];
    char filename[FILENAME_LENGTH];

    strcpy(filename, "");

    if (!strcmp(type, "init")) {
        change_output_chr(file_names[15], numchr);
        change_output_chr(file_names[16], numchr);
        return;
    }

    if (!strcmp(type, "log")) {
        printf("Enter new genotypes file name > ");
        fflush(stdout);
        IgnoreValue(fgets(filename_str, FILENAME_LENGTH - 1, stdin));
        newline;
        if ((sscanf(filename_str, "%s", filename)) == 1) {
            strcpy(file_names[15], filename);
        }
        return;
    }

    if (!strcmp(type, "sum")) {
        printf("Enter new summary file name > ");
        fflush(stdout);
        IgnoreValue(fgets(filename_str, FILENAME_LENGTH - 1, stdin));
        newline;
        if ((sscanf(filename_str, "%s", filename)) == 1) {
            strcpy(file_names[16], filename);
        }
        return;
    }
}

/* Input error model */
static void error_sim_menu(linkage_ped_top *Top,
			   char *err_model,
			   char *file_names[], int numchr)

{

    /* options
       1) No errors.
       2) Apply error model to selected loci.
       3) Apply error to all except selected loci.
       4) Apply error model to selected loci.
       5) Select error model [default = uniform]

       Models
       a) uniform <- single probability of error -  one parameter
       b) Specified propability  for each marker - <num_loc> params
       c) SW2's error model for all markers  -  5 params

    */

    int i, j, loc, not_done=1, done_loci_select;
    char copt[5], str[FILENAME_LENGTH], *str_p, slct[7];
    int num_markers = num_reordered;
    char loci_item[3] = "  ", model_item[4] = "*  ";
    char model_name[3][20] = {"Uniform",
                              "Marker-specific",
                              "SimWalk2"};

    int locus_opt, model_opt, num_select, *loc_select=NULL, *loc_select_i=NULL; // loc_select relatively safe
    char uerr_str[10];                                                          // reset before use
    float err_prob;
    int item;
    char fl_stat[12];
    char serr_str[51], prob_str[FILENAME_LENGTH];
    float serr[5];
    char *Select;
    int model_it, loc_it, ex_loc_it, prob_it, geno_it, sum_it;

    linkage_locus_top *LTop = Top->LocusTop;

    /* can be Uniform
       Marker-specific
       SW2
    */

    locus_opt = 0;
    if (MarkerErrorProb) {
        model_opt = 2; strcpy(model_item, " * ");
    } else {
        model_opt = 1;
    }
    init_error_probs();

    num_select=0;

    /*  printf("%d batch error\n", batchERROR); */
    if (batchERROR) {
        /* first read the error model option */
        if (Mega2BatchItems[/* 22 */ Error_Model].items_read) {
            switch ((int)Mega2BatchItems[/* 22 */ Error_Model].value.copt) {
            case 'U':
                model_opt = 1;
                break;

            case 'M' :
                model_opt = 2;
                break;

            case 'S' :
                model_opt = 3;
                break;

            default:
                invalid_value_field(22);
                break;
            }
        }

        if (Mega2BatchItems[/* 23 */ Error_Probabilities].items_read) {
            switch(model_opt) {
            case 1:
                UniformErrProb = Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[0];
                break;

            case 2:
                if (MarkerErrorProb) {
                    mssgf("Marker-specific error probabilities read in from map file.");
                } else {
                    warnf("No error probabilities specified in map file");
                    warnf("Errors will not be simulated.");
                }
                break;
            case 3:
                for(i=0; i < 5; i++) {
                    SW2ErrProb[i] = Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[i];
                }
                break;

            default:
                break;
            }

        }

        /* Next read the number of selected loci */
        num_select = Mega2BatchItems[/* 21 */ Error_Loci_Num].value.option;
        /* The value has already been verified to be positive */
        loc_select = CALLOC((size_t) num_select, int);
        loc_select_i = CALLOC((size_t) num_select, int);

        if (Mega2BatchItems[/* 19 */ Error_Loci].items_read) {
            /* Selections to apply error model to */
            locus_opt = 2;
            for(loc=0; loc < num_select; loc++) {
                loc_select_i[loc] =
                    loc_select[loc] = Mega2BatchItems[/* 19 */ Error_Loci].value.mult_opts[loc];
                loc_select_i[loc]--;
            }
        } else if (Mega2BatchItems[/* 20 */ Error_Except_Loci].items_read) {
            /* selections to except from errors */
            locus_opt = 3;
            for(loc=0; loc < num_select; loc++) {
                loc_select_i[loc] =
                    loc_select[loc] = Mega2BatchItems[/* 20 */ Error_Except_Loci].value.mult_opts[loc];
                loc_select_i[loc]--;
            }
        }

        num_select = check_loc_selection(LTop, num_select, loc_select, 0);
        if (num_select == 0) {
            /* Error loci read in from batch file has problems */
            warnf("Locus numbers read in from batch file are invalid.");
            warnf("Continuing in Interactive mode.");
            locus_opt = 0;
            free(loc_select);
        }
    }

    if (num_select == 0) {
        Select = CALLOC((size_t) num_markers, char);
        loc_select = CALLOC((size_t) num_markers, int);
        locus_opt = 2;

        for(loc = 0; loc < num_markers; loc++) {
            Select[loc] = ' ';
        }

        while(not_done) {
            switch(model_opt) {
            case 1:
                /* Uniform */
                sprintf(prob_str, "%4.3f", UniformErrProb);
                break;

            case 2:
                break;

            case 3:
                strcpy(prob_str, "");
                for (i=0; i < 5; i++) {
                    grow(prob_str, "%4.3f ", SW2ErrProb[i]);
                }
                break;

            default:
                break;
            }
            printf("Error model and loci selection menu\n");
            printf("0) Done with this menu - please proceed.\n");
            item=1;
            if (model_opt == 2) {
                printf(" %d) Error model : use error probabilities from map file.\n",
                       item);
                model_it=item;
                item++;
            } else {
                printf(" %d)  Select error model              %-15s\n",
                       item, model_name[model_opt-1]);
                model_it=item;
                item++;

                printf("%c%d) Apply error model to selected loci.\n",loci_item[0], item);
                loc_it=item;
                item++;

                printf("%c%d) Apply error to all except selected loci.\n",loci_item[1], item);
                ex_loc_it = item;
                item++;

                printf(" %d) Change error probability       [%s]\n", item, prob_str);
                prob_it=item;
                item++;
            }
            printf(" %d) Mistyping genotypes file name  %-15s %s\n",
                   item, file_names[15], file_status(file_names[15], fl_stat));
            geno_it=item;
            item++;
            printf(" %d) Mistyping summary file name    %-15s %s\n",
                   item, file_names[16], file_status(file_names[16], fl_stat));
            sum_it=item;
            printf(
                "Enter option 1-%d to input new parameters and 0 to proceed > ",
                item);

            fflush(stdout);
            IgnoreValue(fgets(copt, 4, stdin)); newline;

            sscanf(copt, "%d", &not_done);

            if ((model_opt == 2) && (not_done >= 4)) {
                not_done++;
            }

            switch(not_done) {
            case 0:
                break;

            default:
                if (not_done == loc_it || not_done == ex_loc_it) {

                    printf("Select loci by entering their numbers.\n");
                    printf("Separate each number by a space.\n");
                    printf("Numbers are assumed to refer to loci selected by reordering.\n");
                    printf("However, only numbered loci can be selected.\n");
                    printf("  Enter 'v' to view list of ALL loci.\n");
                    printf("  Enter 'm' to select all marker loci.\n");
                    printf("    Loci marked with '*' are current selections.\n");
                    printf("  Enter 'e' to terminate the selection process.\n");

                    locus_opt = not_done;
                    done_loci_select=0;
                    num_select=0;
                    for(loc=0; loc < num_markers; loc++) {
                        loc_select[loc] = 0;
                    }
                    while(!done_loci_select) {
                        printf("Enter a set of loci numbers ('e' to terminate) > ");
                        fflush(stdout);
                        IgnoreValue(fgets(str, FILENAME_LENGTH - 1, stdin));  newline;
                        str_p = &(str[0]);
                        while(str_p != NULL) {
                            if (num_select >= num_markers) {
                                done_loci_select = 1; break;
                            }

                            if (*str_p == 'v' || *str_p == 'V') {
                                for (i=0; i< num_markers; i++)
                                    Select[i] = ' ';

                                if (num_select > 0) {
                                    for (i=0; i< num_select; i++)
                                        if (loc_select[i] <= num_markers &&
                                            loc_select[i] > 0) {
                                            Select[loc_select[i]-1]='*';
                                        }
                                }
                                printf("Original list of loci\n");
                                display_selections(Top, reordered_marker_loci,
                                                   num_markers, Select, NULL, 0, 0);
                                break;
                            } else if (*str_p == 'm' || *str_p == 'M') {
                                printf("Selecting all markers.\n");
                                done_loci_select = 1;
                                i=0;
                                for (loc=0; loc < num_markers; loc++) {
                                    if (LTop->Locus[reordered_marker_loci[loc]].Class ==
                                       MARKER) {
                                        loc_select[i]=loc + 1;
                                        i++;
                                    }
                                }
                                num_select = num_markers;
                                break;
                            } else if (*str_p == 'e' || *str_p == 'E') {
                                str_p++;
                                done_loci_select = 1; break;
                            } else if (*str_p == '\n') {
                                str_p++;
                                /* check the current selections */
                                num_select = check_loc_selection(LTop,num_select,
                                                                 loc_select, 1);
                                if (num_select == 0) {
                                    strcpy(loci_item, "  ");
                                    locus_opt = 0;
                                }
                                break;
                            } else if (isspace((int) *str_p)) str_p++;
                            else if (isdigit((int) *str_p)) {
                                j=0;
                                strcpy(slct, "");
                                while(isdigit((int) *str_p)) {
                                    slct[j] = *str_p;
                                    str_p++; j++;
                                    if (j > 6) {
				      errorvf("Locus number too big in line... %s\n", str);
				      EXIT(INPUT_DATA_ERROR);
				    }
                                }
                                slct[j] = '\0';
                                loc_select[num_select] = atoi(slct); num_select++;
                            } else {
                                printf("Invalid characters in input (%c).\n", *str_p);
                                num_select=0; locus_opt  = 0;
                                strcpy(loci_item, "  ");
                                fflush(stdin);
                                break;
                            }
                        }
                    }

                    /* verify */
                    num_select = check_loc_selection(LTop, num_select,loc_select,1);

                    if (num_select == 0) {
                        strcpy(loci_item, "  "); locus_opt=0;
                        done_loci_select=1;
                    } else {
                        /* display the new selections in order */
                        if (not_done ==  loc_it) {
                            printf("Selections to apply error model to:\n");
                            strcpy(loci_item, "* ");
                        } else if (not_done == ex_loc_it) {
                            printf("Selections to except from applying error model: \n");
                            strcpy(loci_item, " *");
                        }
                        if (loc_select_i != NULL) {
                            free(loc_select_i);
                        }
                        loc_select_i = CALLOC((size_t) num_select, int);
                        for(i=0; i < num_select; i++) {
                            loc_select_i[i] = reordered_marker_loci[loc_select[i]-1];
                        }
                        display_selections(Top, loc_select_i, num_select, NULL, NULL, 0, 0);
                        newline;
                    }
                } else if (not_done == model_it) {
                    /* Change the error model */
                    printf("Error model selection menu \n");
                    printf("%c1) Uniform\n", model_item[0]);
                    printf("%c2) Marker-specific error probabilities from map file\n",
                           model_item[1]);
                    printf("%c3) Simwalk2 error model\n", model_item[2]);
                    printf("\n Select from options 1 - 3 > ");
                    fflush(stdout);
                    IgnoreValue(fgets(copt, 4, stdin)); newline;
                    sscanf(copt, "%d", &model_opt);
                    strcpy(model_item, "   ");
                    model_item[model_opt-1] = '*';

                } else if (not_done == prob_it) {
                    if (model_opt == 1) {
                        printf("Current uniform error probability = %4.3f\n", UniformErrProb);
                        printf("Enter new value (< 9 char, blank to keep old value) > ");
                        fflush(stdout);
                        IgnoreValue(fgets(uerr_str, 9, stdin)); newline;
                        if (sscanf(uerr_str, "%g", &err_prob) >= 1 ) {
                            if (err_prob >= 0.0){
                                UniformErrProb = err_prob;
                            }
                        }
                    } else if (model_opt == 3) {
                        printf("Current SW2 error probabilities :\n");
                        printf("   %4.3f\n   %4.3f\n   %4.3f\n",
                               SW2ErrProb[0], SW2ErrProb[1], SW2ErrProb[2]);
                        printf("   %4.3f\n   %4.3f\n",
                               SW2ErrProb[3], SW2ErrProb[4]);
                        printf("Enter 5 new error probabilities(<= 50 char) > ");
                        fflush(stdout);
                        IgnoreValue(fgets(serr_str, 50, stdin)); newline;
                        if (sscanf(serr_str, "%g %g %g %g %g",
                                  &serr[0],  &serr[1],  &serr[2],
                                  &serr[3],  &serr[4]) == 5) {
                            SW2ErrProb[0] = ((serr[0] >= 0.0)? serr[0] : SW2ErrProb[0]);
                            SW2ErrProb[1] = ((serr[1] >= 0.0)? serr[1] : SW2ErrProb[1]);
                            SW2ErrProb[2] = ((serr[2] >= 0.0)? serr[2] : SW2ErrProb[2]);
                            SW2ErrProb[3] = ((serr[3] >= 0.0)? serr[3] : SW2ErrProb[3]);
                            SW2ErrProb[4] = ((serr[4] >= 0.0)? serr[4] : SW2ErrProb[4]);
                        }
                    }
                } else if (not_done == geno_it) {
                    errsim_file_names("log", file_names, numchr);
                } else if (not_done == sum_it) {
                    errsim_file_names("sum", file_names, numchr);
                } else {
                    warn_unknown(copt);
                }
            }
        }
    }

    if (model_opt == 2) {
        for(loc=0; loc < num_reordered;  loc++) {
            err_model[loc] = 'M';
        }
    }
    /* At this point loc_select_i contains actual locus numbers */
    else if (num_select > 0) {

        /* log the new selections in order */
        sprintf(err_msg, "Error model = %s", model_name[model_opt - 1]);
        errsimf(err_msg);

        if (locus_opt == 2) {
            for(loc=0; loc < num_reordered;  loc++) {
                err_model[loc] = 'X';
            }

            for(loc=0; loc < num_select;  loc++) {
                err_model[loc_select[loc]-1] = model_name[model_opt-1][0];
            }

            errsimf("Selections to apply error model to:\n");

        } else if (locus_opt == 3) {
            for(loc=0; loc < num_reordered;  loc++) {
                err_model[loc] = model_name[model_opt - 1][0];
            }

            for(loc=0; loc < num_select;  loc++) {
                err_model[loc_select[loc]-1] = 'X';
            }

            errsimf("Selections to except from applying error model: \n");

        }

        /* Now switch loc_select to actual locus numbers */
        if (loc_select_i != NULL) {
            for(loc=0; loc < num_select;  loc++) {
                loc_select[loc] = loc_select_i[loc];
            }
            free(loc_select_i);
        }

        log_marker_selections(Top, loc_select, num_select, NULL, (void (*)(const char *)) errsimf,
                              ((model_opt==2)? 1:0));
        switch(model_opt) {
        case 1:
            sprintf(err_msg, "Input error probability = %4.3f", UniformErrProb);
            errsimf(err_msg);
            break;

        case 2:
            break;

        case 3:
            errsimf("Input error probabilities:");
            for (i=0; i < 5; i++) {
                sprintf(err_msg, "   E%d = %4.3f", (i+1), SW2ErrProb[i]);
                errsimf(err_msg);
            }
            break;
        }
        if (InputMode == INTERACTIVE_INPUTMODE) {
            /* Store the batch items */
            Mega2BatchItems[/* 21 */ Error_Loci_Num].value.option = num_select;
            batchf(Error_Loci_Num);
            if (locus_opt == 2) {
                Mega2BatchItems[/* 19 */ Error_Loci].value.mult_opts = CALLOC((size_t) (num_select + 1), int);
                Mega2BatchItems[/* 19 */ Error_Loci].value.mult_opts[num_select] = -99;
                for (loc=0; loc < num_select; loc++) {
                    Mega2BatchItems[/* 19 */ Error_Loci].value.mult_opts[loc] = loc_select[loc] + 1;
                }
                batchf(Error_Loci);
            } else if (locus_opt == 3) {
                Mega2BatchItems[/* 20 */ Error_Except_Loci].value.mult_opts =
                    CALLOC((size_t) (num_select + 1), int);
                Mega2BatchItems[/* 20 */ Error_Except_Loci].value.mult_opts[num_select] = -99;
                for (loc=0; loc < num_select; loc++) {
                    Mega2BatchItems[/* 20 */ Error_Except_Loci].value.mult_opts[loc] = loc_select[loc] + 1;
                }
                batchf(Error_Except_Loci);
            }
            Mega2BatchItems[/* 22 */ Error_Model].value.copt = model_name[model_opt - 1][0];
            batchf(Error_Model);
            if (Mega2BatchItems[/* 22 */ Error_Model].value.copt == 'U') {
                Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues = CALLOC((size_t) 2, double);
                Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[0] = UniformErrProb;
                Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[1] = -9999.0;
                batchf(Error_Probabilities);
            } else if (Mega2BatchItems[/* 22 */ Error_Model].value.copt == 'S') {
                Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues = CALLOC((size_t) 6, double);
                for(i=0; i< 5; i++) {
                    Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[i] = SW2ErrProb[i];
                }
                Mega2BatchItems[/* 23 */ Error_Probabilities].value.mult_fvalues[5] = -9999.0;
                batchf(Error_Probabilities);
            }
            /* finished writing batch file */
        }
    } else {/* No selections */
        for(loc=0; loc < num_reordered;  loc++) {
            err_model[loc] = 'X';
        }
        warnf("No loci selected for mistyping simulation");
        warnf("Error simulation will not be carried out.");
    }

    draw_line();

    return;
}

/* routine for creation of log files */

static void output_err_files(int *loci, int num_loci, char *err_model,
			     listhandle **error_recs,
			     linkage_ped_top *Top, char *file_names[],
			     int *NumTyped, int *NumHetero, int *NumHomo)

{
    /* file_names[15] = errors
       file_names[16] = summary
    */

    /* d1s228     106      4      1 1       1 4
       d1s228     109      4      1 2       2 2
    */


    FILE *fp1, *fp2;
    err_sim_out *err_sim_item;
    int i, loc, num_errors;
    double frac_serrors[5];
    int loc_wid, fid_wid, pid_wid, all_wid;
    char pformat[6], lformat[6], fformat[6], aformats[6], aformatd[6];
    char err_type_name[16];

    field_widths(Top, Top->LocusTop, &fid_wid, &pid_wid,
                 &all_wid, &loc_wid);
    fid_wid = imax((int) strlen("Pedigree"), fid_wid);
    pid_wid = imax((int) strlen("Person"), pid_wid);
    loc_wid = imax((int) strlen("Locus"), loc_wid);
    all_wid = imax((int) strlen("Orig1"), all_wid);

    sprintf(pformat, "%%%ds ", pid_wid);
    sprintf(fformat, "%%%ds ", fid_wid);
    sprintf(aformats, "%%%ds ", all_wid);
    sprintf(aformatd, "%%%dd ", all_wid);
    sprintf(lformat, "%%-%ds ", loc_wid);

    fp1 = fopen(file_names[15], "w");

    fprintf(fp1, lformat, "Locus");
    fprintf(fp1, fformat, "Pedigree");
    fprintf(fp1, pformat, "Person");
    fprintf(fp1, aformats, "Orig1");
    fprintf(fp1, aformats, "Orig2");
    fprintf(fp1, aformats, " Mis1");
    fprintf(fp1, aformats, " Mis2");
    fprintf(fp1, "   Error type  \n");

    fp2 = fopen(file_names[16], "w");
    fprintf(fp2, lformat, "Locus");
    fprintf(fp2, " Genotypes  Errors  Overall_rate");
    /*  printf("%s \n", strchr(err_model, 'S')); */
    if (strchr(err_model, 'S') != NULL) {
        /* AT least one marker has a simwalk2 model */
        fprintf(fp2, " Obs_E1 Obs_E2 Obs_E3 Obs_E4 Obs_E5\n");
        for(i=0; i < 5; i++) {
            frac_serrors[i]=0.0;
        }
    } else {
        fprintf(fp2, "\n");
    }

    for(loc = 0; loc < num_loci; loc++) {
        num_errors=0;
        while (err_sim_item = (err_sim_out *) pop_first_list_entry(error_recs[loc]),
               err_sim_item != NULL) {
            num_errors++;

            fprintf(fp1, lformat, Top->LocusTop->Locus[loci[loc]].Name);
            fprintf(fp1, fformat, Top->Ped[err_sim_item->ped].Name);
            fprintf(fp1, pformat,
                    Top->Ped[err_sim_item->ped].Entry[err_sim_item->entry].OrigID);
            fprintf(fp1, aformatd, err_sim_item->original.all1);
            fprintf(fp1, aformatd, err_sim_item->original.all2);
            fprintf(fp1, " ");
            fprintf(fp1, aformatd, err_sim_item->new_genotype.all1);
            fprintf(fp1, aformatd, err_sim_item->new_genotype.all2);

            switch (err_sim_item->err_type) {
            case E1:
                strcpy(err_type_name, "E1");
                frac_serrors[0] += ((NumHetero[loc] > 0)?
                                    (1.0/(double)(NumHetero[loc])) : 0.0);
                break;

            case E2:
                strcpy(err_type_name, "E2");
                frac_serrors[1] += ((NumHetero[loc] > 0)?
                                    (1.0/(double)(NumHetero[loc])) : 0.0);
                break;

            case E3_homo:
                strcpy(err_type_name, "E3 Homozygote");
                frac_serrors[2] += ((NumHomo[loc] > 0)?
                                    (1.0/(double)(NumHomo[loc])) : 0.0);
                break;

            case E3_hetero:
                strcpy(err_type_name, "E3 Heterozygote");
                frac_serrors[3] += ((NumHetero[loc] > 0)?
                                    (1.0/(double)(NumHetero[loc])) : 0.0);
                break;

            case E4:
                strcpy(err_type_name, "E4");
                frac_serrors[3] += ((NumHomo[loc] > 0)?
                                    (1.0/(double)(NumHomo[loc])) : 0.0);
                break;

            case E5_homo:
                strcpy(err_type_name, "E5 Homozygote");
                frac_serrors[4] += ((NumHomo[loc] > 0)?
                                    (1.0/(double)(NumHomo[loc])) : 0.0);
                break;

            case E5_hetero:
                strcpy(err_type_name, "E5 Heterozygote");
                frac_serrors[3] += ((NumHetero[loc] > 0)?
                                    (1.0/(double)(NumHetero[loc])) : 0.0);
                break;

            case Uniform:
                strcpy(err_type_name, "Uniform");
                break;

            case Marker:
                strcpy(err_type_name, "Marker Specific");
                break;

            default:
                printf("UNKNOWN ERROR TYPE\n");
                break;
            }

            fprintf(fp1, "%15s\n", err_type_name);
        }

        fprintf(fp2, lformat, Top->LocusTop->Locus[loci[loc]].Name);

        fprintf(fp2, "%10d  %6d         %4.3f",
                NumTyped[loc], num_errors,
                ((NumTyped[loc] > 0)? ((double)num_errors/(double)(NumTyped[loc])) :
                 0.0));
        if (strchr(err_model, 'S') != NULL) {
            for(i=0; i < 5; i++) {
                fprintf(fp2, "  %4.3f", frac_serrors[i]);
            }
        }

        fprintf(fp2, "\n");
    }

    fclose(fp1);
    fclose(fp2);
    mssgf("Mega2 created the following files for mistyping simulation:");
    sprintf(err_msg, "  Genotypes file: %s", file_names[15]);
    mssgf(err_msg);
    sprintf(err_msg, "  Summary file: %s", file_names[16]);
    mssgf(err_msg);
    draw_line();
    return;
}

/* Main routine to process all peigrees */


void simulate_errors(linkage_ped_top *Top, int numchr, char *file_names[])

{

    char *error_model, *err_error_model;
    int loc, loc1, num_error_loci, *selected_loci;
    genotype (*err_model_function)(genotype gen, int number,
                                   linkage_locus_rec Loc,
                                   sw2_error_class *err_type);
    listhandle **error_recs;
    int *NumTyped, *NumHomo, *NumHetero;

    if (main_chromocnt > 1) {
        errsim_file_names("init", file_names, 0);
    } else if (main_chromocnt == 1) {
        errsim_file_names("init", file_names, numchr);
    } else {
        warnf("No mapped loci, cannot simulate errors.");
        warnf("Skipping mistyping simulation step.");
        draw_line();
        return;
    }

    NumTyped = CALLOC((size_t) num_reordered, int);
    NumHomo = CALLOC((size_t) num_reordered, int);
    NumHetero = CALLOC((size_t) num_reordered, int);
    error_model = CALLOC((size_t) num_reordered+1, char);
    error_model[num_reordered] = '\0';
    error_sim_menu(Top, error_model, file_names, numchr);

    /* create the list of loci that will have errors */

    num_error_loci = 0;
    for (loc=0; loc < num_reordered; loc++) {
        if (error_model[loc] == 'U') {
            if (UniformErrProb > 0.0) {
		num_error_loci++;
            } else {
                warnf("Uniform mistyping probability value is 0 and");
                warnvf("Locus %s requires a non-zero uniform mistyping probability.\n",
                        Top->LocusTop->Locus[reordered_marker_loci[loc]].Name);
            }
        } else if (error_model[loc] == 'S') {
            if ((SW2ErrProb[0]+SW2ErrProb[1]+SW2ErrProb[2]+
                SW2ErrProb[3]+SW2ErrProb[4]) > 0.0) {
                num_error_loci++;
            } else {
                warnf("All SimWalk2 mistyping probabilities are 0 and");
                warnvf("Locus %s requires non-zero SimWalk2 mistyping probabilities.\n",
                        Top->LocusTop->Locus[reordered_marker_loci[loc]].Name);
            }
        } else if (error_model[loc] == 'M') {
            if (Top->LocusTop->Marker[reordered_marker_loci[loc]].error_prob >
               0.0) {
                num_error_loci++;
            } else {
                warnvf("Marker-specific mistyping probability for locus %s is 0.\n",
                        Top->LocusTop->Locus[reordered_marker_loci[loc]].Name);
                error_model[loc] = 'X';
            }
        }
    }

    if (num_error_loci == 0 && error_model[0] == 'M') {
        warnf("No locus -zero input mistyping probability,");
        warnf("Proceeding without random error simulation.");
        free(error_model); free(NumTyped);
        free(NumHomo); free(NumHetero);
        return;
    }

    selected_loci = CALLOC((size_t) num_error_loci, int);
    err_error_model = CALLOC((size_t) num_error_loci, char);
    error_recs = CALLOC((size_t) num_error_loci, list_rec*);

    loc1=0;
    for (loc=0; loc < num_reordered; loc++) {
        if (error_model[loc] != 'X') {
            selected_loci[loc1]=reordered_marker_loci[loc];
            err_error_model[loc1] = error_model[loc];
            loc1++;
        }
    }

    /* now create the genotype matrices,
       GenotypeIndex is of dimension LTop->LocusCnt X NUMGENOS X 3
    */

    gen_all_indexes(Top->LocusTop, num_error_loci, selected_loci);

    /* Now run through the pedigrees and change the genotypes according
       to selected model */

    for (loc=0; loc < num_error_loci; loc++) {
        NumTyped[loc] = 0;
        NumHomo[loc] = 0;
        NumHetero[loc] = 0;
        switch(err_error_model[loc]) {
        case 'X' :
            break;

        case 'U':
            err_model_function=uniform_error_model;
            error_at_locus(selected_loci[loc], loc, Top, err_model_function,
                           &(error_recs[loc]), &(NumTyped[loc]),
                           NULL, NULL);
            break;

        case 'M':
            err_model_function = marker_error_model;
            error_at_locus(selected_loci[loc], loc, Top, err_model_function,
                           &(error_recs[loc]), &(NumTyped[loc]),
                           NULL, NULL);
            break;

        case 'S':
            err_model_function = sw2_error_model;
            error_at_locus(selected_loci[loc], loc, Top, err_model_function,
                           &(error_recs[loc]), &(NumTyped[loc]),
                           &(NumHomo[loc]), &(NumHetero[loc]));
            break;

        default:
            break;

        }
    }

    mssgf("Done introducing random errors in selected loci.");
    free_gen_indexes(Top->LocusTop, num_error_loci,
                     selected_loci,
                     &GenotypeIndexes);

    /* Create the log and summary files */

    output_err_files(selected_loci, num_error_loci, error_model,
                     error_recs, Top,
                     file_names, NumTyped, NumHomo, NumHetero);
    free(NumTyped);
}
