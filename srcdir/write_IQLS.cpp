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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA .

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

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "pedtree_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

#include "loop.h"

/*
     create_summary_ext.h:  aff_status_entry
     error_messages_ext.h:  errorf mssgf my_calloc
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_loci_on_chromosomes get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
            pedtree_ext.h:  convert_to_pedtree free_all_including_ped_top
         user_input_ext.h:  test_modified
              utils_ext.h:  EXIT draw_line log_line script_time_stamp
*/

void create_IQLS_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt);
static void IQLS_file_names(char *file_names[], int *glob_files,
                            int has_orig, int has_uniq);
/* static void write_IQLS_pedigree(char *outfl_name, linkage_ped_top *Top,
   int pedfile_type, int pwid, int fwid);
   static void write_IQLS_marker(linkage_ped_top *Top, char *outfl_name); */

static void write_IQLS_parameter(linkage_ped_top *Top, int numchr, char *file_names[]);

/* static void write_Idcoefs_pedigree(linkage_ped_top *Top, char *outfl_name); */
static void write_Idcoefs_study(linkage_ped_top *Top, char *outfl_name,
                                int pedfile_type, int pwid, int fwid);

static void write_IQLS_shell_script(linkage_ped_top *Top,
                                    int numchr, char *file_names[],
                                    int first_time);
int index_renumber_ped(ped_tree *Ped, int *index);

/*
  IQLS option:

  The idcoefs program requires that "The pedigree must be ordered in
  such a way that parents always appear earlier in the file than
  children."  From testing and looking at the idcoefs code, it looks
  like this is a statement about ordering, and not numbering.  Still,
  people need to be reordered from top to bottom.  It does not suffice
  to write out all the founders first, as some internal nodes might be
  non-founders, yet need to be written out before their children.

  The IQLS program requires that the individual IDs must be unique.  Looks
  like IQLS requires that allele labels be coded as letters (A,C,T,G).

  The IQLS program requires IBD coefficients "for every pair of eligible
  individuals within each family (including an individual with
  himself/herself ), where an individual is eligible if he or she has
  either (1) known affection status or (2) non-missing genotype for at
  least one marker. (E.g. an individual with unknown phenotype is still
  eligible if he or she has any non-missing genotype information.)"

  Tasks:

  1) Order each family from 'top' to 'bottom'.
  (Instead of actually re-ordering, could create an index array).

  2) Check size of maximum family.

  3) Based on size of maximum family, generate unique integer IDs.

  4) Figure out who are the 'eligible' individuals - those with known
  affection status or non-missing genotypes or both.

  5) Using the new unique integer IDs, write out a pedigree file for
  idcoefs in the new order, so parents are before children.

  6) Write a study file listing the eligible pairs for idcoefs.

  7) Write IQLS files:
  a) Use original allele labels where possible.
  b) Must have physical map information.

*/

static void write_IQLS_pedigree(char *outfl_name, linkage_ped_top *Top,
                                int pedfile_type, int pwid, int fwid)
//pedfile_type is not used
{
/* IQLS pedigree/phenotype file:
   This file contains the pedigree and phenotype information. Individuals who are not
   listed in this file will not be included in the analysis. The columns in the file should
   be organized as follows:
   1 1 7 6 1 1
   1 2 7 6 2 2
   1 3 7 6 1 2
   1 7 0 0 1 1
   1 6 0 0 2 0
   2 11 18 19 2 2
   2 12 18 19 1 1
   2 18 0 0 1 0
   2 19 15 16 2 1
   2 15 0 0 1 0
   2 16 0 0 2 0
   (1) (2) (3) (4) (5) (6)

   Column (1) family ID (positive integer)
   Column (2) individual ID (positive integer; must be unique)
   Column (3) father's ID (0 if the individual is a founder)
   Column (4) mother's ID (0 if the individual is a founder)
   Column (5) sex (1=male, 2=female)
   Column (6) affection status (0=unknown, 1=unaffected, 2=affected)

   Sampled individuals who are unrelated to anyone else in the sample should be included
   by giving each such person their own unique family ID (as well as unique individual
   ID) and setting both parent IDs to 0. There is no limit on the number of individuals
   nor on the number of families. Each individual should be entered only once. The
   individual ID is required to be unique (e.g. it cannot be reused in a different family).
   Individuals from the same family should appear in a single cluster, though there is
   no requirement on the order of individuals within a family nor on the order in which
   different families are listed.

*/
    struct IQLS_pedigree: public loop::outer, loop::ped_per {

        IQLS_pedigree(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::ped_per(Top) { }
        void make_file() {
            msgvf("        IQLS pedigree file:        %s/%s\n", *_opath, file_names[0]);
            run_loop(file_names[0]);
        }
        void inner()  {
            pr_id();
            pr_parent();
            pr_sex();
            pr_aff();
            pr_nl();
        }
    } *xp = new IQLS_pedigree(Top);

    xp->load_formats(fwid, pwid, -1);
    xp->iterate();

    delete xp;
}

static void write_IQLS_marker(linkage_ped_top *Top, char *outfl_name, int pwid, int fwid, int mwid)
{
/*
  This file contains the marker data. All markers should be on the same chromosome.
  (To analyze more than one chromosome, a separate run must be performed for each
  chromosome, with each chromosome having its own marker data file.) The columns in
  the file should be organized as follows
  marker chromosome position orientation allele1 allele2 1 2 3 7 6 11 12 18 19 15 16
  rs7909677 10 101955 + A G AG AA AA AA AG AG GG GG AG AG AG
  rs9419560 10 142201 + A G AA GG AG AG AG AG NN GG AG AA GG
  rs9419419 10 153707 - T C TC TC CC TC TC TT TC TT TC TT CC
  (1) (2) (3) (4) (5) (6) (7) (8) (9) (10) (11) (12) (13) (14) (15) (16) (17)
  Column (1) marker rs number
  Column (2) chromosome
  Column (3) physical position (in nucleotides)
  Column (4) strand orientation (+) same strand as HapMap, (-) opposite strand
  from HapMap)
  Column (5) nucleotide for allele 0
  Column (6) nucleotide for allele 1
  Columns (7)... marker genotypes (NN for missing genotype)
  The first row of the file must contain the column headings. The headings for the first
  6 columns can be arbitrary, but should not contain any space characters. Columns
  7 and beyond contain marker genotype data for the sampled individuals, and each of
  these columns must have the corresponding individual ID number as the heading.
  The order of the individuals is not required to be the same as the order in the pedigree
  file. The column headings must specify the order. There is no limit on the number
  of markers. However, all markers should be on the same chromosome and should be
  listed in increasing order of their position on the chromosome. All individuals in the
  marker data file should also appear in the phenotype data file, otherwise, they will not
  be included in the analysis.

  The number of columns should be the same for every marker: Use NN for missing
  genotype.

  NOTE:
  If a marker was originally input with text allele labels in the input annotated
  pedigree file, then those original allele labels are used on output into the IQLS
  marker file, and the strand orientation is set to '+'.

  If a marker was originally input with numeric allele labels, then it is output using
  the dummy alleles 'A' and 'G', and the strand orientation is set to '-'.

  As far as I can tell, it appears that the orientation information is not used by the
  IQLS program.

*/

    struct IQLS_marker: public loop::outer, loop::loci_ped_per {

        IQLS_marker(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::loci_ped_per(Top) { }
        void make_file() {
            msgvf("        IQLS marker file:          %s/%s\n", *_opath, file_names[1]);
            run_loop(file_names[1]);
        }
        void file_header() {
            int ped, per;
            linkage_ped_rec *tpe;
            /* Write header line */
            pr_printf("marker    chrom. pos'n  orient. all1 all2");
            for (ped = 0; ped < _Top->PedCnt; ped++) {
                if (UntypedPeds != NULL) {
                    if (UntypedPeds[ped]) {
                        continue;
                    }
                }
                for (per = 0; per < _Top->Ped[ped].EntryCnt; per++) {
                    tpe = &(_Top->Ped[ped].Entry[per]);
                    pr_per(tpe);
                }
            }
            pr_nl();

            /* Check that the selected map is a physical map */
            if (base_pair_position_index < 0) {
                log_line(mssgf);
                errorf("IQLS requires that you select a physical map.");
                EXIT(SYSTEM_ERROR);
            }
        }
        void loci_start() {
            if (_tle->Marker->chromosome != UNKNOWN_CHROMO) {
                pr_printf("%-15s %2d %d ",
                          _tle->Name,
                          _tle->Marker->chromosome,
                          (int) _EXLTop->EXLocus[_locus].positions[base_pair_position_index]);
            }
            if (_tle->Marker->Props.Numbered.Recoded) {
                pr_printf("   +  %4s %4s ",
                          _tle->Allele[0].name,
                          _tle->Allele[_tle->AlleleCnt == 1 ? 0 : 1].name);
            } else {
                pr_printf("   -  %4s %4s ","A","G"); /* Use dummy labels A/G for markers input with numbered alleles */
            }
        }
        void inner() {
            if (_tle->Marker->Props.Numbered.Recoded) {
                pr_printf(" %1s%1s ",
                          (_allele1== 0) ? "N" : _tle->Allele[_allele1 - 1].name,
                          (_allele2== 0) ? "N" : _tle->Allele[_allele2 - 1].name);
            } else {
                char all1[2], all2[2];
                /* If input allele labels were numeric, then use dummy alleles 'A' and 'G' and
                 * set strand orientation to '-'.
                 */
                switch(_allele1) {
                case 0:
                    strcpy(all1, "N");
                    break;
                case 1:
                    strcpy(all1, "A");
                    break;
                case 2:
                    strcpy(all1, "G");
                    break;
                default:
                    strcpy(all1, "X");
                    break;
                }
                switch(_allele2) {
                case 0:
                    strcpy(all2, "N");
                    break;
                case 1:
                    strcpy(all2, "A");
                    break;
                case 2:
                    strcpy(all2, "G");
                    break;
                default:
                    strcpy(all2, "X");
                    break;
                }
                pr_printf(" %1s%1s ", all1, all2);
            }
        }
        void loci_end() {
            pr_nl();
        }
    } *xp = new IQLS_marker(Top);

    xp->load_formats(fwid, pwid, mwid);
    xp->iterate();

    delete xp;
}

static void write_IQLS_parameter(linkage_ped_top *Top, int numchr, char *files[]) {
    struct IQLS_parameter: public loop::chr, loop::null {

        IQLS_parameter(linkage_ped_top *Top) : person_locus_entry(Top), loop::chr(Top), loop::null(Top) {}
        void make_file() {
            msgvf("        IQLS parameter file:       %s/%s\n", *_opath, file_names[2]);
            run_loop(file_names[2]);
        }
        void inner() {
            pr_printf("0.05\n");
            pr_printf("2\n");
        }
    } *xp = new IQLS_parameter(Top);

    xp->iterate();

    delete xp;
}

static void write_IQLS_shell_script(linkage_ped_top *Top, int numchr, char *file_names[],
                                    int first_time)
{
    struct IQLS_shell_script: public loop::outer, sh_util {
        typedef char *str;
        str *file_names;

        IQLS_shell_script(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), sh_util(Top) { }
        void make_file() {
            msgvf("        IQLS shell script file:    %s/%s\n", *_opath, file_names[3]);
            run_loop(file_names[3]);
        }
        void file_header() {
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
	    // This handles the environment variable setup to allow the checking
	    // functions in 'batch_run' to work correctly...
	    fprintf_env_checkset_csh(_filep, "_IDCOEFS", "idcoefs");
	    fprintf_env_checkset_csh(_filep, "_IQLS", "IQLS");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            char cmd[2*FILENAME_LENGTH];
            char out_fl[2*FILENAME_LENGTH];

            if (_numchr > 0)
                pr_printf("echo Running Idcoefs and IQLS on chromosome %d\n", _numchr);

            sprintf(cmd, "$_IDCOEFS -p %s -s %s -o IBDfile\n", file_names[4], file_names[5]);
            sh_rm("IBDfile");
            sh_run("IDCOEFS", cmd);
            fprintf_status_check_csh(_filep, "IDCOEFS", 1);
            sh_status("idcoefs", "IBDfile");

            sprintf(out_fl, "IQLStest_out.%02d", _numchr);
            sprintf(cmd, "$_IQLS -pheno %s -geno %s -r %s -ibd IBDfile\n", file_names[0], file_names[1], file_names[2]);
            sh_rm("IQLStest.out");
            sh_run("IQLS", cmd);
            //fprintf_status_check_csh(_filep, "IQLS", 1); // exits with status == 1
            pr_printf("# IQLS will output the results of all haplotype tests to a file called “IQLStest.out”.\n");
            pr_printf("# See: http://www.stat.uchicago.edu/~mcpeek/software/IQLS/IQLS_Documentation.pdf\n");
            pr_printf("# Section 2.6 Running IQLS\n");
            sh_save_output("IQLS", "IQLStest.out", out_fl);
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new IQLS_shell_script(Top);

    xp->file_names = file_names;

    xp->iterate();

    delete xp;
}

static void write_Idcoefs_pedigree(linkage_ped_top *Top, char *outfl_name,
                                   int pedfile_type, int pwid, int fwid) {
/*
  Pedigree file:
  Each line in the pedigree file should correspond to one individual in
  the pedigree. Each line should be separated into at least three
  columns with white space separating the columns. The first column is
  the ID of the person. The next two columns are the IDs of the parents
  (whether the mother or father comes first doesn't matter). Any
  additional columns beyond these first three are ignored. The IDs must
  all be integers. The ID 0 is reserved and means "unknown." If a person
  is a founder, both parents must have ID 0. All non-founders must have
  both parents specified. The pedigree must be ordered in such a way
  that parents always appear earlier in the file than children.

  NOTE the requirement regarding ORDERING of the pedigree!

*/
    struct Idcoefs_pedigree: public loop::outer, loop::ped_per {
        ped_top *PedTreeTop;
        int     *index;

        Idcoefs_pedigree(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::ped_per(Top) { }
        void make_file() {
            msgvf("        Idcoefs pedigree file:     %s/%s\n", *_opath, file_names[4]);
            run_loop(file_names[4]);
        }
        void trait_start() {
            PedTreeTop = convert_to_pedtree(_Top, 0);
        }
        void ped_start() {
            index = CALLOC((size_t) PedTreeTop->PedTree[_ped].EntryCnt, int);
            index_renumber_ped(&(PedTreeTop->PedTree[_ped]), index);
        }
        void inner() {
            _tpe = &(_Top->Ped[_ped].Entry[index[_per]]);
            pr_per();
            pr_parent();
            pr_printf(" ");
            pr_fam();
            pr_printf(_pformat, _tpe->OrigID);
            pr_nl();
        }
        void ped_end() { free(index); }
        void trait_end() { free_all_including_ped_top(PedTreeTop, NULL, NULL); }
    } *xp = new Idcoefs_pedigree(Top);

    xp->load_formats(fwid, pwid, -1);
    xp->iterate();

    delete xp;
}

int index_renumber_ped(ped_tree *Ped, int *index) {
    int entry, flag;
    register ped_rec *Entry;
    int entry_ID;
    int *place = CALLOC((size_t) Ped->EntryCnt, int);
    int *mother = CALLOC((size_t) Ped->EntryCnt, int);
    int *father = CALLOC((size_t) Ped->EntryCnt, int);

    /* place[] maintains a map of old IDs to new IDs */
    /* initialize place[] */
    for (entry = 0; entry < Ped->EntryCnt; entry++) {
        place[entry] = 0;
        Entry = &(Ped->Entry[entry]);
        mother[entry] = (Entry->Mother == NULL ? 0 : Entry->Mother->ID);
        father[entry] = (Entry->Father == NULL ? 0 : Entry->Father->ID);
    }

    /* Now set new IDs until we get a permissible, complete pedigree */
    /* This is the algorithm Dr. Weeks suggested. */
    entry_ID = 1;
    do {
        flag = 0;
        for (entry = 0; entry < Ped->EntryCnt; entry++) {
            Entry = &(Ped->Entry[entry]);
            if (place[entry] != 0)
                continue;
            if (((mother[entry] == 0) || (place[mother[entry] - 1] != 0))
                && ((father[entry] == 0) || (place[father[entry] - 1] != 0))) {
                place[entry] = entry_ID;
                index[entry_ID-1] = entry;
                entry_ID++;
                flag = 1;
            }
        }
    } while (flag);
    free(mother);
    free(father);

    return (entry_ID);
}

static void write_Idcoefs_study(linkage_ped_top *Top, char *outfl_name,
                                int pedfile_type, int pwid, int fwid)
{
/*
  Studyfile:
  The study file may have either of two formats. The first is that each
  line contains a pair of IDs separated by white space. The program will
  then compute the condensed identity coefficients for each of these
  pairs. The other format is to have one ID per line, in which case
  condensed identity coefficients will be computed for every possible
  pair of IDs. */

/* For IQLS:  need to create a list of eligible individuals.

   This file contains condensed identity coefficients for every pair of eligible individuals
   within each family (including an individual with himself/herself ), where an individual
   is eligible if he or she has either (1) known affection status or (2) non-missing genotype
   for at least one marker. (E.g. an individual with unknown phenotype is still eligible if
   he or she has any non-missing genotype information.)
   8
   IBD coefficients should be included for every pair of eligible individuals who have the
   same family ID (including each individual with himself/herself ). A sampled individual
   who does not share a family ID with anyone else in the sample, would be represented
   in the markid file by a single line that gives the IBD coefficients for the person with
   himself/herself.

*/
    struct IDcoefs_study : public loop::outer, loop::ped_per {
        int *eligible;

        IDcoefs_study(linkage_ped_top *Top) : person_locus_entry(Top), loop::outer(Top), loop::ped_per(Top) { }
        void make_file() {
            msgvf("        Idcoefs study file:        %s/%s\n", *_opath, file_names[5]);
            run_loop(file_names[5]);
        }
        void ped_start () {
            eligible = CALLOC((size_t) (_Top->Ped[_ped].EntryCnt), int);
        }
        void inner() {
            int aff;
            eligible[_per] = 0;
            if (_tpe->IsTyped > 0)
                eligible[_per] = 1;

            /* trait locus */
            switch(_tte->Type) {
            case AFFECTION:
                // linkage.h    :linkage_pedrec_data is a union (Affection(2xint),
                // Quant(2xint), Alleles         (2xint), RAlleles (2xchar*)
                if (_tte->Pheno->Props.Affection.ClassCnt == 1)
                    aff = _tpe->Pheno[_trait].Affection.Status;
                else           
                    aff = aff_status_entry(_tpe->Pheno[_trait].Affection.Status, _tpe->Pheno[_trait].Affection.Class, _tte);
                if (aff > 0)
                    eligible[_per]=1;
                break;
            default:
                break;
            }
        }
        void ped_end() {
            int per, per2;
            linkage_ped_rec *tpe, *tpe2;
            for (per = 0; per < _Top->Ped[_ped].EntryCnt; per++) {
                tpe = &(_Top->Ped[_ped].Entry[per]);
                /* if (eligible[per]) */
                for (per2 = per; per2 < _Top->Ped[_ped].EntryCnt; per2++) {
                    tpe2 = &(_Top->Ped[_ped].Entry[per2]);
                    /* if (eligible[per2]) { */
                    pr_printf("%5s %5s\n", tpe->UniqueID, tpe2->UniqueID);
                    /* } */
                }
            }
            free(eligible);
        }
    } *xp = new IDcoefs_study(Top);

    xp->load_formats(fwid, pwid, -1);
    xp->iterate();

    delete xp;
}

/* This allows the
   pointer to the LPedTop to be renamed for convenience */

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

void create_IQLS_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt)

{
    int combine_chromo;
    int numchr=0, first_time=1;
    int fwid, pwid, mwid;
    analysis_type analysis = IQLS;
    linkage_locus_top *LTop = Top->LocusTop;


    /* get the output file names, also use this function to set output
       related parameters, like whether an overall shell script should
       be created for all chromosome etc. */

    IQLS_file_names(file_names, &combine_chromo, Top->UniqueIds, Top->UniqueIds);
    LoopOverChrm = ! combine_chromo;

    /* This displays a line on the screen stating that the following
       files will be created for IQLS */
    create_mssg(analysis);

    /* This function sets the formatting for pedigree and individual
       ids, as well as locus Ids */

    field_widths(Top, LTop, &fwid, &pwid, NULL, &mwid);

    /* Add in code for setting other run-related parameters that are
       specific to this option, and if these are going to remain
       constant across all chromosomes */

    /* Now write the main loop for creating chromosome-specific files:
     */

    omit_peds(untyped_ped_opt, Top);

    write_IQLS_pedigree(file_names[0], Top, pedfile_type, fwid, pwid);

    write_IQLS_marker(Top, file_names[1], fwid, pwid, mwid);

    write_IQLS_parameter(Top, numchr, file_names);

    write_IQLS_shell_script(Top, numchr, file_names, first_time);

    write_Idcoefs_pedigree(Top, file_names[4], pedfile_type, fwid, pwid);

    write_Idcoefs_study(Top, file_names[5], pedfile_type, fwid, pwid);

}

/* other functions that are usually relevant */

static void IQLS_file_names(char *file_names[], int *combine_chromo,
                            int has_orig, int has_uniq)
{
    int item, choice = -1;
    char cchoice[10], stem[15];
    char fl_stat[12];
    analysis_type analysis = IQLS;

    *combine_chromo=0;
    if (main_chromocnt > 1) {
// This is the batch file item that controls whether you wish to comnine the
// chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
// This is a derective from the user which will override the default...
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
            *combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
    }

    if (DEFAULT_OUTFILES)  {
        mssgf("Output file names set to defaults.");
        return;
    }

    ((num_traits > 1 && LoopOverTrait)?
     strcpy(stem, "in each subdir") : strcpy(stem, ""));

    if (main_chromocnt > 1) {
        if (!(*combine_chromo)) {
            strcpy(stem, "stem");
            // Here -9 is a flag to return only the file name portion...
            analysis->replace_chr_number(file_names, -9);
        } else {
            // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
            analysis->replace_chr_number(file_names, 0);
        }
    }

    while (choice != 0) {
        draw_line();
        printf("  IQLS/Idcoefs file name menu\n");
        print_outfile_mssg();
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        if (main_chromocnt == 1 || *combine_chromo) {
            printf( " 1) IQLS pedigree filename:              %-15s\t%s\n",
                    file_names[0], file_status(file_names[0], fl_stat));
            printf(
                " 2) IQLS marker file name:               %-15s\t%s\n",
                file_names[1], file_status(file_names[1], fl_stat));
            printf(
                " 3) IQLS parameter file name:            %-15s\t%s\n",
                file_names[2], file_status(file_names[2], fl_stat));
            printf(
                " 4) IQLS/Idcoefs shell file name:        %-15s\t%s\n",
                file_names[3], file_status(file_names[3], fl_stat));
            printf(
                " 5) Idcoefs pedigree file name:          %-15s\t%s\n",
                file_names[4], file_status(file_names[4], fl_stat));
            printf(
                " 6) Idcoefs study file name:             %-15s\t%s\n",
                file_names[5], file_status(file_names[5], fl_stat));
            item = 6;
        } else {
            printf(" 1) IQLS pedigree file name %s:          %-15s\n", stem,
                   file_names[0]);
            printf(" 2) IQLS marker file name %s:            %-15s\n", stem,
                   file_names[1]);
            printf(" 3) IQLS parameter file name %s:         %-15s\n",
                   stem, file_names[2]);
            printf(" 4) IQLS/Idcoefs shell file name %s:     %-15s\t\n", stem,
                   file_names[3]);
            printf(" 5) Idcoefs pedigree name %s:            %-15s\t\n", stem,
                   file_names[4]);
            printf(" 6) Idcoefs study name %s:               %-15s\t\n", stem,
                   file_names[5]);
            item = 6;
        }

        printf("Select options 0-%d > ", item);

        choice = -1;
        fcmap(stdin, "%s", cchoice); newline;
        sscanf(cchoice, "%d", &choice);
        test_modified(choice);

        switch (choice) {
        case 0:
            break;
        case 1:
            printf("New IQLS pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[0]);
            printf("\n");
            break;
        case 2:
            printf("New IQLS marker file name %s > ", stem);
            fcmap(stdin, "%s", file_names[1]);
            printf("\n");
            break;
        case 3:
            printf("New IQLS parameter file name %s > ", stem);
            fcmap(stdin, "%s", file_names[2]);
            printf("\n");
            break;
        case 4:
            printf("New IQLS/Idcoefs shell file name %s > ", stem);
            fcmap(stdin, "%s", file_names[3]);
            printf("\n");
            break;
        case 5:
            printf("New Idcoefs pedigree file name %s > ", stem);
            fcmap(stdin, "%s", file_names[4]);
            printf("\n");
            break;
        case 6:
            printf("New Idcoefs study file name %s > ", stem);
            fcmap(stdin, "%s", file_names[5]);
            printf("\n");
            break;
        default:
            warn_unknown(cchoice);
            break;
        }
    }
}
