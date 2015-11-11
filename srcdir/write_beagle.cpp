/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
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
#include <ctype.h>

#include "common.h"
#include "typedefs.h"

#include "loop.h"
#include "sh_util.h"
#include "loop_templates.h"

#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "genetic_utils_ext.h"
#include "write_files_ext.h"

#include "write_beagle_ext.h"

/**
   These classes define the sub-options associated with the different Genotype file
   that Beagle can produce. See the section on 'Beagle Genotype files'.
 */
CLASS_BEAGLE_UNPHASED_UNRELATED *BEAGLE_UNPHASED_UNRELATED = new CLASS_BEAGLE_UNPHASED_UNRELATED();
CLASS_BEAGLE_UNPHASED_TRIO *BEAGLE_UNPHASED_TRIO = new CLASS_BEAGLE_UNPHASED_TRIO();
CLASS_BEAGLE_UNPHASED_PAIR *BEAGLE_UNPHASED_PAIR = new CLASS_BEAGLE_UNPHASED_PAIR();

//
// The emtry.cpp code uses 0 for missing in pr_pheno(), pr_father(), pr_mother()...
#define MISSING_ID_CODE                    "0"
#define MISSING_ALLELE_CODE                MISSING_ID_CODE


/**

   BEAGLE_INCLUDE_PHENOTYPE

   Beagle only needs a phenotype for association testing. Furthermore, for association testing
   it will give an error if the phenotype is 0 (e.g., missing). So, it would seem that there are
   two ways of proceeding:
   1) never include phenotypes, and do not do association testing, or
   2) only include person data that has a phenotype (e.g., it's not missing) informaiton and do association testing.

   Implemented here is #1.
 */


/*
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_loci_on_chromosomes get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp summary_time_stamp

 http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

 2 BEAGLE file formats
 BEAGLE assumes that any file that has a name ending in “.gz” is compressed with gzip.

 2.2 Genotype likelihoods file format (not currently handled)
 This file contains genotype likelihoods for unphased, unrelated data when performing genotype
 phasing and imputation. It is used for genotype phasing and imputation. These likelyhood
 files are individual specific, unlike the frequency files available in Mega2.
*/

static void inner_file_names(char **file_names, const char *num, const char *stem = "beagle");

//
// Comparison functions for qsort which is used in @see write_BEAGLE_marker_file
//
// The comparison function must return an integer less than, equal to, or greater
// than zero if the first argument is considered to be respectively less than,
// equal to, or greater than the second.
//
// Since the comparison is based on the value of something other than the data,
// a pointer to it must be saved by write_BEAGLE_marker_file().
//
// Beagle does not use sex maps, only sex-averaged maps. We generate an error if the
// user sets up for a sex-specific or a female map.

static linkage_ped_top *Top_for_qsort;
static int cmp_marker_basepair_positions (const void *a, const void *b) {
    const int *a_i = (const int *)a;
    const int *b_i = (const int *)b;
    const double a_p = Top_for_qsort->EXLTop->EXLocus[*a_i].positions[base_pair_position_index];
    const double b_p = Top_for_qsort->EXLTop->EXLocus[*b_i].positions[base_pair_position_index];

    return  (a_p < b_p) ? -1 : ((a_p > b_p) ? 1 : 0);
}
static int cmp_marker_genetic_positions (const void *a, const void *b) {
    const int *a_i = (const int *)a;
    const int *b_i = (const int *)b;
    const double a_g = Top_for_qsort->EXLTop->EXLocus[*a_i].positions[genetic_distance_index];
    const double b_g = Top_for_qsort->EXLTop->EXLocus[*b_i].positions[genetic_distance_index];

    int ret = (a_g < b_g) ? -1 : ((a_g > b_g) ? 1 : 0);
    
    // We should add a test against the physical position if there is a tie in the genetic position....
    if (ret == 0 && base_pair_position_index >= 0) ret = cmp_marker_basepair_positions (a, b);

    return  ret;
}

// 'markers_on_chromosome()' called from 'dataloop::loci::data_loop(*_opath, )' lies, it will return all loci
// (including traits), but the traits are returned first. This is a hack to find the first
// marker after those traits. At some point 'markers_on_chromosome()' should be refactored
// to just get the markers. 'dataloop::loci::data_loop(*_opath, )' "deals with this" by throwing away
// anything that is not a marker, but this does not help qsort.
extern int ChrLociFirstMarkerIndex, ChrLociLastMarkerIndex;

static void sort_basepair_positions() {
    // Naturally this sort will fail (EXC_BAD_ACCESS) if there is just one marker on the chromosome.
    // It's not a very interesting case, but it happens in the test code.
    if (ChrLociLastMarkerIndex-ChrLociFirstMarkerIndex > 0)
        qsort((void *)&ChrLoci[ChrLociFirstMarkerIndex],
              (size_t)(ChrLociLastMarkerIndex-ChrLociFirstMarkerIndex+1),
              sizeof(int), cmp_marker_basepair_positions);
}

static void sort_genetic_positions() {
    if (ChrLociLastMarkerIndex-ChrLociFirstMarkerIndex > 0)
        qsort((void *)&ChrLoci[ChrLociFirstMarkerIndex],
              (size_t)(ChrLociLastMarkerIndex-ChrLociFirstMarkerIndex+1),
              sizeof(int), cmp_marker_genetic_positions);
}

//HERE
struct LCLchr: public FLOOPchr {
    LCLchr(linkage_ped_top *Top, const char *f_name, const char *f_mode) : FLOOPchr(Top, f_name, f_mode) {}
    void make_file() {
        if (_fnumchr > lastautosome) return; // because we have already issued a warning.
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct LCLboth: public FLOOPboth {
    LCLboth(linkage_ped_top *Top, const char *f_name, const char *f_mode) : FLOOPboth(Top, f_name, f_mode) {}
    void make_file() {
        if (_fnumchr > lastautosome) {
            if (*file_type) warnf("Beagle analysis can only be run on autosomes.");
            return;
        }
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

/**
 @brief Output the Beagle Marker file
 
 The format of this type of file is described in the Beagle documentat and section listed below...
 
 http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf
 2.4 Marker file format
 
 A <marker file> is used to reconstruct marker order. It lists all markers in a chromosome in order.
 The marker file contains the marker name, followed by the posiiton, followed by the alleles.
 The position can be phsical or genetic (cM) depending on the type of analysis required by Beagle.
 
 The markers must be given in increasing chromosomal order, so they must be sorted. Comparison
 functions for qsort are defined above. The decision to place the sorting here was made for several
 reasons: 1) Beagle seeems to be the only analysis type that requires sorted markers; 2) qsort
 requires access to a static comparison function which requies the use of non-instance variables
 (e.g., an instance method will just not work). A natural place to put this sorting wound be
 loop.cpp:dataloop::loci::data_loop. However, a static method with access to a global 'linkage_ped_top *'
 would be required to be set.
 
 If you are performing HBD/IBD detection, the position of each marker must be given using the centiMorgan
 scale (Mega2 internal default). If you are not performing HBD/IBD detection, the marker positions on a
 chromosome can be base positions, genetic positions, or sequential integers.
 
 We output all maps that are available. However, since genetic map is [currently] required, it will be
 used in the batch file if both physical and genetic maps exist.

 NOTES:
 For association testing, pairs and trios are assumed to have unaffected parents and affected offspring.
 It cannot handle other affection status configurations.

 Beagle will ignore the affection status line (A ...) unless the "trait=" argument is specified.
 The "trait=" argument is used only for association testing.  So for phasing, imputation, and IBD
 detection, Beagle will ignore the affection status and won't complain if the affection status line
 is missing or if individuals have a missing affection status ('0'). When association testing,
 the phased parental haplotypes are coded as 2/1 depending on whether or not the haplotype is transmitted.
 Children must be omitted from the analysis (either by filtering the input data or by using the
 excludesample= argument). Beagle is limited in this way because it only performs Fisher exact tests on 2x2 tables.
 In association testing a 0 (unknown) phenotype will generate an error.
 
 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] file_names List of file neames used by this analysis mode.
 @param[in] pwid       Width of (max characers taken by) person id in the output file.
 @param[in] fwid       Width of (max characers taken by) pedigree/family id in the output file.
 @return void
 */
static void write_BEAGLE_marker_file(linkage_ped_top *Top, char *file_names[],
                                     const int pwid, const int fwid)
{
    // Stash Top away for the qsort cmp functions above...
    Top_for_qsort = Top;

    if (base_pair_position_index >= 0) {
        // If we have the data create the file...

        LCLchr *floopM = new LCLchr(Top, file_names[1], "w");
        floopM->file_type = "        BEAGLE base position marker file:    ";

        struct save_base_pair_position_markers: public dataloop::loci {
            save_base_pair_position_markers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) { }

            // The name 'file_header()' is a bit of a misnomer. It's really the place to do some work
            // just before processing on the data to build the output file begins. You could output a header,
            // or in this case sort the data that will be used to build the output file (see method 'inner()').
            void file_header() {
                sort_basepair_positions();
            }
            void inner() {
                pr_printf("%15s ", _tlocusp->LocusName);
                //pr_marker_name();
                (void)pr_physical_distance(NULL); // for the marker...
                pr_marker_alleles(); // Print the alleles associated with the marker...
                pr_nl();
            }
        } *sbppM = new save_base_pair_position_markers(Top, floopM);
        sbppM->load_formats(fwid, pwid, -1);

        floopM->iterate();

        delete sbppM;
        delete floopM;
    }

    if (genetic_distance_index >= 0) {
        // If we have the data create the file...

        LCLchr *floopgdM = new LCLchr(Top, file_names[2], "w");
        floopgdM->file_type = "        BEAGLE genetic distance marker file: ";

        struct save_genetic_distance_markers: public dataloop::loci {
            save_genetic_distance_markers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci(Top, fl) { }

            void file_header() {
                sort_genetic_positions();
            }
            void inner() {
                int warnp;
                double genetic_distance = get_genetic_distance(&warnp);
                pr_genetic_distance_warning(warnp);
                pr_printf("%15s ", _tlocusp->LocusName);
                // Haldane cM was intended for use in Beagle...
#ifdef DONT_CONVERT_TO_HALDANE_FOR_BEAGLE
                warnf("*****************************************************************");
                warnf("Not converting genetic distance to Haldane for debugging purposes!");
                warnf("*****************************************************************");
#else /* DONT_CONVERT_TO_HALDANE_FOR_BEAGLE */
                if (_LTop->map_distance_type == 'k') genetic_distance = kosambi_to_haldane(genetic_distance);
#endif /* DONT_CONVERT_TO_HALDANE_FOR_BEAGLE */
#if defined(_WIN) || defined(MINGW)
                pr_printf("%.6f ", genetic_distance);
#else
                pr_printf("%.6lf ", genetic_distance);
#endif
                pr_marker_alleles(); // Print the alleles associated with the marker...
                pr_nl();
            }
        } *sgdM = new save_genetic_distance_markers(Top, floopgdM);
        
        sgdM->load_formats(fwid, pwid, -1);

        floopgdM->iterate();

        delete sgdM;
        delete floopgdM;
    }
}

/**
 Beagle Genotype files

 The format of this type of file is described in the Beagle documentat and section listed below...

 http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf
 2.1 Genotypes file format

 The Beagle genotype files have a format where rows are variables and columns are individuals.
 The first column describes the data that will appear on the line (as per below).
 The second column contains the data variable name which must be unique.

 A <Beagle file> contain two sections: header lines and marker lines. Header lines (optional) are
 all of the lines preceeding the marker lines. This imposes the restruction that marker lines are
 requied to appear at the end of the genotype file.
 
 Marker lines (M) require two columns for diploid data. An identifier and an affection status
 are given for each allele (column), and will typically be the same for both alleles.
 The markers in the Beagle file must be in chromosomal order. This implies that a file can
 contain markers for only one chromosome. It is also possible to divide chromosome data among
 multiple files. In this case a <marker file> must be provided because the data file is not
 requied to contain all markers.

 Beagle can process phased and unphased data files. Mega2 will only produce unphased data files
 based on the Beagle analysis sub-option (defined above and in the .h file). These types are:
 unphased unrelated, unphased trio, and unphased pair data.
*/

/**
 @brief Output the Beagle Genotype file; Unphased unrelated data format

 From Section "2.1 Genotypes file format" of the Beagle manual...
 http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

 In Unphased unrelated data each pair of columns (following the first two) will give the
 genotype for each individual.

 Processing of these files requires the 'unphased' command line argument to Beagle.

 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] file_names List of file neames used by this analysis mode.
 @param[in] pwid       Width of (max characers taken by) person id in the output file.
 @param[in] fwid       Width of (max characers taken by) pedigree/family id in the output file.
 @return void
*/
static void write_BEAGLE_genotype_unphased_unrelated_file(linkage_ped_top *Top, char *file_names[],
                                                          const int pwid, const int fwid)
{
    // Pedigree...
    // 'fileloop::both' loops over chromosome, and then trait...

    LCLboth *floopP = new LCLboth(Top, file_names[0], "w");
    floopP->file_type = "        BEAGLE genotype file:                ";

    struct save_fams: public dataloop::ped_per {
        save_fams(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("P Pedigree ");
        }
        void inner() {
            // FUTURE: Also, only write data for individuals where some marker information is known.
            pr_fam(); pr_fam(); // Must be printed out twice because there are two allele columns...
        }
        void file_trailer() {
            pr_nl();
        }
    } *sP = new save_fams(Top, floopP);
    
    sP->load_formats(fwid, pwid, -1);

    floopP->iterate();

    delete sP;
    delete floopP;
    
    // Person...
    // The person or identifier line is not currently required, but will be in future versions of Beagle.

    LCLboth *floopI = new LCLboth(Top, file_names[0], "a");
    floopI->file_type = "";

    struct save_pers: public dataloop::ped_per {
        save_pers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("I ID ");
        }
        void inner() {
            // For diploid data, the identifier and affection status will typically be the same for both alleles.
            // Must be printed out twice because there are two allele columns...
            // The sampling ID must be unique throughout....

            // Here we use UniqueID explicitly...
            pr_printf("%s %s ", _tpersonp->UniqueID, _tpersonp->UniqueID);
            //pr_per(); pr_per();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sI = new save_pers(Top, floopI);
    
    sI->load_formats(fwid, pwid, -1);

    floopI->iterate();

    delete sI;
    delete floopI;
    
    // Father...

    LCLboth *floopPID = new LCLboth(Top, file_names[0], "a");
    floopPID->file_type = "";

    struct save_fathers: public dataloop::ped_per {
        save_fathers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("PID Father ");
        }
        void inner() {
            // The indexes for IDs in linkage_ped_rec (_tpersonp) are 1-based with 0 used for NA.
            // Since something is output for the child/individual then we need to
            // output something here so that the data in the subsequent columns is consistent.
            const char *fid = (_tpersonp->Father != 0) ? _tpedtreep->Entry[_tpersonp->Father-1].UniqueID : MISSING_ID_CODE;
            pr_printf("%s %s ", fid, fid);
            //pr_father(); pr_father();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sPID = new save_fathers(Top, floopPID);
    
    sPID->load_formats(fwid, pwid, -1);

    floopPID->iterate();

    delete sPID;
    delete floopPID;
    
    // Mother...

    LCLboth *floopMID = new LCLboth(Top, file_names[0], "a");
    floopMID->file_type = "";

    struct save_mothers: public dataloop::ped_per {
        save_mothers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("MID Mother ");
        }
        void inner() {
            const char *mid = (_tpersonp->Mother != 0) ? _tpedtreep->Entry[_tpersonp->Mother-1].UniqueID : MISSING_ID_CODE;
            pr_printf("%s %s ", mid, mid);
            //pr_mother(); pr_mother();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sMID = new save_mothers(Top, floopMID);
    
    sMID->load_formats(fwid, pwid, -1);

    floopMID->iterate();

    delete sMID;
    delete floopMID;
    
    // Sex...

    LCLboth *floopC = new LCLboth(Top, file_names[0], "a");
    floopC->file_type = "";

    struct save_sexs: public dataloop::ped_per {
        save_sexs(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }
        
        void file_header() {
            pr_printf("C Sex ");
        }
        void inner() {
            pr_sex(); pr_sex();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sC = new save_sexs(Top, floopC);
    
    sC->load_formats(fwid, pwid, -1);

    floopC->iterate();

    delete sC;
    delete floopC;

#ifdef BEAGLE_INCLUDE_PHENOTYPE
    // Phenotypes...

    LCLboth *floopAT = new LCLboth(Top, file_names[0], "a");
    floopAT->file_type = "";

    struct save_phenotypes: public dataloop::ped_per {
        save_phenotypes(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            // The affection status 'A' is of the form (1 == unaffected, 2 == affected) for each person.
            // It is necessary only when performing association testing.
            // Quantitative traits are marked with a 'T', and Covariates with a 'C'.
            // NOTE: Quantitative and Covariate traits are not currently used by Beagle.

            if (_ttraitp != (linkage_locus_rec *)NULL) pr_printf("%s %s ", (_ttraitp->Type == AFFECTION ? "A" : "T"), _ttraitp->LocusName);
            else warnf("The genotype file contains no phenotypes.");
        }
        void inner() {
            pr_pheno(); pr_pheno();
        }
        void file_trailer() {
            if (_ttraitp != (linkage_locus_rec *)NULL) pr_nl();
        }
    } *sAT = new save_phenotypes(Top, floopAT);

    sAT->_trait_affect = true;
    sAT->load_formats(fwid, pwid, -1);

    floopAT->iterate();

    delete sAT;
    delete floopAT;

#endif /* BEAGLE_INCLUDE_PHENOTYPE */
    
    // Markers...
    // Markers must be listed in chromosomial order, and must all appear at the end of the file.
    // So, they must be written in the same order as written in the marker file @see write_BEAGLE_marker_file
    //
    // NOTE: Mega2 annotated input file format (7.2) uses {M|X|Y} for marker names rather then just 'M' here...

    LCLboth *floopM = new LCLboth(Top, file_names[0], "a");
    floopM->file_type = "";

    struct save_markers: public dataloop::loci_ped_per {
        save_markers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci_ped_per(Top, fl) { }

        void file_header() {
            if (genetic_distance_index >= 0) {
                // Beagle only uses sex-averaged maps, so if we get to here that's what we have
                sort_genetic_positions();
            } else if (base_pair_position_index >= 0) {
                sort_basepair_positions();
            }
        }
        void loci_start() {
            pr_printf("M %s ", _tlocusp->LocusName);
        }
        void inner() {
            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                pr_printf("%s %s ",
                          recode_name(_allele1, MISSING_ALLELE_CODE, "?"),
                          recode_name(_allele2, MISSING_ALLELE_CODE, "?"));
            } else {
                linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                pr_printf("%s %s ", format_allele(Locus, _allele1), format_allele(Locus, _allele2));
            }
        }
        void loci_end() {
            pr_nl();
        }
    } *sM = new save_markers(Top, floopM);
    
    sM->load_formats(fwid, pwid, -1);

    floopM->_trait_affect = true;
    floopM->iterate();

    delete sM;
    delete floopM;
}

/**
 @brief Output the Beagle Genotype file; Unphased [related] trio data format

 From Section "2.1 Genotypes file format" of the Beagle manual...
 http://faculty.washington.edu/browning/beagle/beagle_3.3.2_31Oct11.pdf

 In Unphased [related] trio data each 6-tuple of columns (following the first two) will give the
 genotype of one parent-offspring trio <1st_parent, 2nd_parent, child>. If both Father and Mother
 are not present, then nothing is written for the child.

 Processing of these files requires the 'trios' command line argument to Beagle.

 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] file_names List of file neames used by this analysis mode.
 @param[in] pwid       Width of (max characers taken by) person id in the output file.
 @param[in] fwid       Width of (max characers taken by) pedigree/family id in the output file.
 @return void
*/
static void write_BEAGLE_genotype_unphased_trio_file(linkage_ped_top *Top, char *file_names[],
                                                     const int pwid, const int fwid)
{
    // Pedigree...

    LCLboth *floopP = new LCLboth(Top, file_names[0], "w");
    floopP->file_type = "        BEAGLE genotype file:                ";

    struct save_fams: public dataloop::ped_per {
        save_fams(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("P Pedigree ");
        }
        void inner() {
            pr_fam(); pr_fam();
            pr_fam(); pr_fam();
            pr_fam(); pr_fam();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sP = new save_fams(Top, floopP);
    
    sP->load_formats(fwid, pwid, -1);

    floopP->iterate();

    delete sP;
    delete floopP;
    
    // Father, Mother, Child (6-tuple)...

    LCLboth *floopI = new LCLboth(Top, file_names[0], "a");
    floopI->file_type = "";

    struct save_pers: public dataloop::ped_per {
        save_pers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("I ID ");
        }
        void file_trailer() {
            pr_nl();
        }
        void inner() {
            const char *fid = (_tpersonp->Father != 0) ? _tpedtreep->Entry[_tpersonp->Father-1].UniqueID : MISSING_ID_CODE;
            const char *mid = (_tpersonp->Mother != 0) ? _tpedtreep->Entry[_tpersonp->Mother-1].UniqueID : MISSING_ID_CODE;
            pr_printf("%s %s %s %s %s %s ", fid, fid, mid, mid, _tpersonp->UniqueID, _tpersonp->UniqueID);
            //pr_father(); pr_father();
            //pr_mother(); pr_mother();
            //pr_per(); pr_per();
        }
    } *sI = new save_pers(Top, floopI);
    
    sI->load_formats(fwid, pwid, -1);

    floopI->iterate();

    delete sI;
    delete floopI;
    
    // Sex...

    LCLboth *floopC = new LCLboth(Top, file_names[0], "a");
    floopC->file_type = "";

    struct save_sexs: public dataloop::ped_per {
        save_sexs(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("C Sex ");
        }
        void file_trailer() {
            pr_nl();
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                pr_sex(&_tpedtreep->Entry[_tpersonp->Father-1]); pr_sex(&_tpedtreep->Entry[_tpersonp->Father-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            if (_tpersonp->Mother != 0) {
                pr_sex(&_tpedtreep->Entry[_tpersonp->Mother-1]); pr_sex(&_tpedtreep->Entry[_tpersonp->Mother-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_sex(); pr_sex();
        }
    } *sC = new save_sexs(Top, floopC);
    
    sC->load_formats(fwid, pwid, -1);

    floopC->iterate();

    delete sC;
    delete floopC;
    
#ifdef BEAGLE_INCLUDE_PHENOTYPE
    // Phenotypes...

    LCLboth *floopAT = new LCLboth(Top, file_names[0], "a");
    floopAT->file_type = "";

    struct save_phenotypes: public dataloop::ped_per {
        save_phenotypes(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            // The affection status 'A' is of the form (1 == unaffected, 2 == affected) for each person.
            // It is necessary only when performing association testing.
            // Quantitative traits are marked with a 'T', and Covariates with a 'C'.
            // NOTE: Quantitative and Covariate traits are not currently used by Beagle.
            if (_ttraitp != (linkage_locus_rec *)NULL) pr_printf("%s %s ", (_ttraitp->Type == AFFECTION ? "A" : "T"), _ttraitp->LocusName);
            else warnf("The genotype file contains no phenotypes.");
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                pr_pheno(&_tpedtreep->Entry[_tpersonp->Father-1]); pr_pheno(&_tpedtreep->Entry[_tpersonp->Father-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            if (_tpersonp->Mother != 0) {
                pr_pheno(&_tpedtreep->Entry[_tpersonp->Mother-1]); pr_pheno(&_tpedtreep->Entry[_tpersonp->Mother-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_pheno(); pr_pheno();
        }
        void file_trailer() {
            if (_ttraitp != (linkage_locus_rec *)NULL) pr_nl();
        }
    } *sAT = new save_phenotypes(Top, floopAT);
    
    sAT->load_formats(fwid, pwid, -1);

    floopAT->_trait_affect = true;
    floopAT->iterate();

    delete sAT;
    delete floopAT;
#endif /* BEAGLE_INCLUDE_PHENOTYPE */
    
    // Markers...
    // Markers must be listed in chromosomial order, and must all appear at the end of the file.
    // NOTE: Mega2 annotated input file format (7.2) uses {M|X|Y} for marker names rather then just 'M' here...

    LCLboth *floopM = new LCLboth(Top, file_names[0], "a");
    floopM->file_type = "";

    struct save_markers: public dataloop::loci_ped_per {
        save_markers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci_ped_per(Top, fl) { }

        void file_header() {
            if (genetic_distance_index >= 0) {
                // Beagle only uses sex-averaged maps, so if we get to here that's what we have
                sort_genetic_positions();
            } else if (base_pair_position_index >= 0) {
                sort_basepair_positions();
            }
        }
        void loci_start() {
            pr_printf("M %s ", _tlocusp->LocusName);
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                int _allele1_f, _allele2_f;
                linkage_ped_rec  *_tpe_f = &(_Top->Ped[_ped].Entry[_tpersonp->Father-1]);
                get_2alleles(_tpe_f->Marker, _locus, &_allele1_f, &_allele2_f);
                if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                    pr_printf("%s %s ",
                              recode_name(_allele1_f, MISSING_ALLELE_CODE, "?"),
                              recode_name(_allele2_f, MISSING_ALLELE_CODE, "?"));
                } else {
                    linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                    pr_printf("%s %s ",
                              format_allele(Locus, _allele1_f), format_allele(Locus, _allele2_f));
                }
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }

            if (_tpersonp->Mother != 0) {
                int _allele1_m, _allele2_m;
                linkage_ped_rec  *_tpe_m = &(_Top->Ped[_ped].Entry[_tpersonp->Mother-1]);
                get_2alleles(_tpe_m->Marker, _locus, &_allele1_m, &_allele2_m);
                if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                    pr_printf("%s %s ",
                              recode_name(_allele1_m, MISSING_ALLELE_CODE, "?"),
                              recode_name(_allele2_m, MISSING_ALLELE_CODE, "?"));
                } else {
                    linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                    pr_printf("%s %s ",
                              format_allele(Locus, _allele1_m), format_allele(Locus, _allele2_m));
                }
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }

            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                pr_printf("%s %s ",
                          recode_name(_allele1, MISSING_ALLELE_CODE, "?"),
                          recode_name(_allele2, MISSING_ALLELE_CODE, "?"));
            } else {
                linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                pr_printf("%s %s ",
                          format_allele(Locus, _allele1), format_allele(Locus, _allele2));
            }
        }
        void loci_end() {
            pr_nl();
        }
    } *sM = new save_markers(Top, floopM);
    
    sM->load_formats(fwid, pwid, -1);

    floopM->_trait_affect = true;
    floopM->iterate();

    delete sM;
    delete floopM;
}

/**
 @brief Output the Beagle Genotype file; Unphased [related] pair data format

 In Unphased [related] pair data each 4-tuple of columns (following the first two) will give the
 genotype of one parent-offspring pair <father/mother, child>.

 Processing of these files requires the 'pairs' command line argument to Beagle.

 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] file_names List of file neames used by this analysis mode.
 @param[in] pwid       Width of (max characers taken by) person id in the output file.
 @param[in] fwid       Width of (max characers taken by) pedigree/family id in the output file.
 @return void
*/
static void write_BEAGLE_genotype_unphased_pair_file(linkage_ped_top *Top, char *file_names[],
                                                     const int pwid, const int fwid)
{
    // Pedigree...

    LCLboth *floopP = new LCLboth(Top, file_names[0], "w");
    floopP->file_type = "        BEAGLE genotype file:                ";

    struct save_fams: public dataloop::ped_per {
        save_fams(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("P Pedigree ");
        }
        void inner() {
            pr_fam(); pr_fam(); pr_fam(); pr_fam();

            pr_fam(); pr_fam(); pr_fam(); pr_fam();
        }
        void file_trailer() {
            pr_nl();
        }
    } *sP = new save_fams(Top, floopP);
    
    sP->load_formats(fwid, pwid, -1);

    floopP->iterate();

    delete sP;
    delete floopP;
    
    // Each set of four consecutive columns (beginning with columns 3-6) gives the genotype data
    // for one parent-offspring pair. In each set of four columns, the first two columns give the
    // genotypes for the genotyped parent, and the last two columns give the genotypes for the offspring.

    LCLboth *floopI = new LCLboth(Top, file_names[0], "a");
    floopI->file_type = "";

    struct save_pers: public dataloop::ped_per {
        save_pers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("I ID ");
        }
        void file_trailer() {
            pr_nl();
        }
        void inner() {
            const char *fid = (_tpersonp->Father != 0) ? _tpedtreep->Entry[_tpersonp->Father-1].UniqueID : MISSING_ID_CODE;
            const char *mid = (_tpersonp->Mother != 0) ? _tpedtreep->Entry[_tpersonp->Mother-1].UniqueID : MISSING_ID_CODE;
            pr_printf("%s %s %s %s ",fid, fid, _tpersonp->UniqueID, _tpersonp->UniqueID);

            pr_printf("%s %s %s %s ", mid, mid, _tpersonp->UniqueID, _tpersonp->UniqueID);
        }
    } *sI = new save_pers(Top, floopI);
    
    sI->load_formats(fwid, pwid, -1);

    floopI->iterate();

    delete sI;
    delete floopI;
    
    // Sex...

    LCLboth *floopC = new LCLboth(Top, file_names[0], "a");
    floopC->file_type = "";

    struct save_sexs: public dataloop::ped_per {
        save_sexs(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::ped_per(Top, fl) { }

        void file_header() {
            pr_printf("C Sex ");
        }
        void file_trailer() {
            pr_nl();
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                pr_sex(&_tpedtreep->Entry[_tpersonp->Father-1]); pr_sex(&_tpedtreep->Entry[_tpersonp->Father-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_sex(); pr_sex();

            if (_tpersonp->Mother != 0) {
                pr_sex(&_tpedtreep->Entry[_tpersonp->Mother-1]); pr_sex(&_tpedtreep->Entry[_tpersonp->Mother-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_sex(); pr_sex();
        }
    } *sC = new save_sexs(Top, floopC);
    
    sC->load_formats(fwid, pwid, -1);

    floopC->iterate();

    delete sC;
    delete floopC;
    
#ifdef BEAGLE_INCLUDE_PHENOTYPE
    // Phenotypes...

    LCLboth *floopAT = new LCLboth(Top, file_names[0], "a");
    floopAT->file_type = "";

    struct save_phenotypes: public dataloop::trait_ped_per {
        save_phenotypes(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::trait_ped_per(Top, fl) { }

        void file_header() {
            // The affection status 'A' is of the form (1 == unaffected, 2 == affected) for each person.
            // It is necessary only when performing association testing.
            // Quantitative traits are marked with a 'T', and Covariates with a 'C'.
            // NOTE: Quantitative and Covariate traits are not currently used by Beagle.
            if (_ttraitp != (linkage_locus_rec *)NULL) pr_printf("%s %s ", (_ttraitp->Type == AFFECTION ? "A" : "T"), _ttraitp->LocusName);
            else warnf("The genotype file contains no phenotypes.");
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                pr_pheno(&_tpedtreep->Entry[_tpersonp->Father-1]); pr_pheno(&_tpedtreep->Entry[_tpersonp->Father-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_pheno(); pr_pheno();

            if (_tpersonp->Mother != 0) {
                pr_pheno(&_tpedtreep->Entry[_tpersonp->Mother-1]); pr_pheno(&_tpedtreep->Entry[_tpersonp->Mother-1]);
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            pr_pheno(); pr_pheno();
        }
        void file_trailer() {
            if (_ttraitp != (linkage_locus_rec *)NULL) pr_nl();
        }
    } *sAT = new save_phenotypes(Top, floopAT);
    
    sAT->load_formats(fwid, pwid, -1);

    floopAT->_trait_affect = true;
    floopAT->iterate();

    delete sAT;
    delete floopAT;
#endif /* BEAGLE_INCLUDE_PHENOTYPE */
    
    // Markers...
    // Markers must be listed in chromosomial order, and must all appear at the end of the file.
    // NOTE: Mega2 annotated input file format (7.2) uses {M|X|Y} for marker names rather then just 'M' here...

    LCLboth *floopM = new LCLboth(Top, file_names[0], "a");
    floopM->file_type = "";

    struct save_markers: public dataloop::loci_ped_per {
        save_markers(linkage_ped_top *Top, fileloop::fileloop_data *fl) : dataloop::loci_ped_per(Top, fl) { }

        void file_header() {
            if (genetic_distance_index >= 0) {
                // Beagle only uses sex-averaged maps, so if we get to here that's what we have
                sort_genetic_positions();
            } else if (base_pair_position_index >= 0) {
                sort_basepair_positions();
            }
        }
        void loci_start() {
            pr_printf("M %s ", _tlocusp->LocusName);
        }
        void inner() {
            if (_tpersonp->Father != 0) {
                int _allele1_f, _allele2_f;
                linkage_ped_rec  *_tpe_f = &(_Top->Ped[_ped].Entry[_tpersonp->Father-1]);
                get_2alleles(_tpe_f->Marker, _locus, &_allele1_f, &_allele2_f);
                if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                    pr_printf("%s %s ",
                              recode_name(_allele1_f, MISSING_ALLELE_CODE, "?"),
                              recode_name(_allele2_f, MISSING_ALLELE_CODE, "?"));
                } else {
                    linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                    pr_printf("%s %s ",
                              format_allele(Locus, _allele1_f), format_allele(Locus, _allele2_f));
                }
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }

            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                pr_printf("%s %s ",
                          recode_name(_allele1, MISSING_ALLELE_CODE, "?"),
                          recode_name(_allele2, MISSING_ALLELE_CODE, "?"));
            } else {
                linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                pr_printf("%s %s ",
                          format_allele(Locus, _allele1), format_allele(Locus, _allele2));
            }

            if (_tpersonp->Mother != 0) {
                int _allele1_m, _allele2_m;
                linkage_ped_rec  *_tpe_m = &(_Top->Ped[_ped].Entry[_tpersonp->Mother-1]);
                get_2alleles(_tpe_m->Marker, _locus, &_allele1_m, &_allele2_m);
                if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                    pr_printf("%s %s ",
                              recode_name(_allele1_m, MISSING_ALLELE_CODE, "?"),
                              recode_name(_allele2_m, MISSING_ALLELE_CODE, "?"));
                } else {
                    linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                    pr_printf("%s %s ",
                              format_allele(Locus, _allele1_m), format_allele(Locus, _allele2_m));
                }
            } else {
                pr_printf("%s %s ", MISSING_ID_CODE, MISSING_ID_CODE);
            }
            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                pr_printf("%s %s ",
                          recode_name(_allele1, MISSING_ALLELE_CODE, "?"),
                          recode_name(_allele2, MISSING_ALLELE_CODE, "?"));
            } else {
                linkage_locus_rec *Locus = &(_LTop->Locus[_locus]);
                pr_printf("%s %s ",
                          format_allele(Locus, _allele1), format_allele(Locus, _allele2));
            }
        }
        void loci_end() {
            pr_nl();
        }
    } *sM = new save_markers(Top, floopM);
    
    sM->load_formats(fwid, pwid, -1);

    floopM->_trait_affect = true;
    floopM->iterate();

    delete sM;
    delete floopM;
}


/**
 @brief Output an initial shell control file

 This file will run run Beagle in accordance to the analysis sub-type specified (see above
 text on 'Beagle Genotype files').

 @param[in] Top Root   Data structure references to pedigrees, persons, and loci.
 @param[in] analysis   A class instance associated with the current analysis.
 @param[in] file_names List of file neames used by this analysis mode.
 @return void
*/
static void write_BEAGLE_sh(linkage_ped_top *Top,
                            analysis_type *analysis,
                            char *file_names[])
{
    int top_shell = (LoopOverChrm && main_chromocnt > 1) || (LoopOverTrait && num_traits > 1);
    
    DTshell *sh = 0;
    if (top_shell) {
        sh = new DTshell(Top);
        sh->filep_open(output_paths[0], file_names[4], "w");
        sh->sh_main();
    }
    
    LCLboth *floop = new LCLboth(Top, file_names[3], "w");
    floop->file_type = "        BEAGLE shell file:                   ";

    struct BEAGLE_sh_script: public DTshell {
        BEAGLE_sh_script(linkage_ped_top *Top, fileloop::fileloop_data *fl, analysis_type *analysis) : DTshell(Top, fl) {
          this->analysis = analysis;
        }

        typedef char *str;
        str *file_names;
        DTshell *sh;
        analysis_type *analysis;
        
        void file_header() {
            if (sh) sh->sh_sh(this);
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
#ifdef RUNSHELL_SETUP
            pr_printf("# Allows the user to define their own version of Java.\n");
            fprintf_env_checkset_csh(_filep, "_JAVA", "java");
            pr_printf("# Please set the environment variable BEAGLE_JAR to point to the location of\n");
            pr_printf("# the beagle.jar file. The default is the current directory. You may wish to set\n");
            pr_printf("# this environment variable in you shell's rc file (e.g., .bash_profile) so that it\n");
            pr_printf("# only needs to be done once.\n");
            fprintf_env_checkset_csh(_filep, "BEAGLE_JAR", "beagle.jar");
            pr_printf("if ( ! -f $BEAGLE_JAR ) then\n");
            pr_printf("  echo 'The beagle.jar file was not found - please set your BEAGLE_JAR'\n");
            pr_printf("  echo 'environment variable properly so it can be found.'\n");
            pr_printf("  echo 'If using Bash and ksh you would use something like this:'\n");
            pr_printf("  echo 'export BEAGLE_JAR=/usr/local/src/beagle/beagle.jar'\n");
            pr_printf("  echo 'If using csh you would use something like this:'\n");
            pr_printf("  echo 'setenv BEAGLE_JAR /usr/local/src/beagle/beagle.jar'\n");
            pr_printf("  echo 'For further details, please see the 'Beagle' section of the Mega2 documentation'\n");
            pr_printf("  exit 0\n");
            pr_printf("endif\n");
            pr_printf("\n");
#endif /* RUNSHELL_SETUP */
        }
        void inner () {
            pr_printf("# Process %s format %s data...\n", (*analysis)->_name, (*analysis)->_subname);
            pr_printf("# For further information on Beagle command line options please see:\n");
            pr_printf("# http://faculty.washington.edu/browning/beagle/beagle.html\n");
            pr_printf("# \n");
            pr_printf("echo ' NOTE: If you encounter a 'Java heap space' exception while running Beagle,'\n");
            pr_printf("echo ' please consult the 'Memory Management' section of the Beagle documentation'\n");
            pr_printf("echo ' for further assistance.'\n");
            pr_printf("echo\n");
            // For a complete list of Beagle command line arguments,
            // please see section 3.2 of the Beagle manual.
            // Also, see section 5.1 Memory Management to give Java more memory for it's computations.
            pr_printf("$_JAVA -jar $BEAGLE_JAR");
            
            if ((*analysis)->_suboption == BEAGLE_UNPHASED_UNRELATED_SUBOPTION)
                pr_printf(" unphased=%s", file_names[0]);
            else if ((*analysis)->_suboption == BEAGLE_UNPHASED_TRIO_SUBOPTION)
                // not permitted to have any Mendelian inconsistencies.
                // Genotypes causing Mendelian inconsistencies for a trio must be replaced with missing genotypes.
                pr_printf(" trios=%s", file_names[0]);
            else if ((*analysis)->_suboption == BEAGLE_UNPHASED_PAIR_SUBOPTION)
                // not permitted to have any Mendelian inconsistencies.
                // Genotypes causing Mendelian inconsistencies for a trio must be replaced with missing genotypes.
                pr_printf(" pairs=%s", file_names[0]);
            
            // From the discussion in section 2.4 of the Beagle manual it would seem that the genetic position
            // can be used in all instances, so if have the genetic positions we will give use of them presidence
            // over use of the physical positions.
            if (genetic_distance_index >= 0)
                pr_printf(" markers=%s", file_names[2]);
            else if (base_pair_position_index >= 0)
                pr_printf(" markers=%s", file_names[1]);

            // The missing command line argument is required.
            pr_printf(" missing=%s out=%s", MISSING_ALLELE_CODE, file_names[5]);

            if ((*analysis)->_suboption == BEAGLE_UNPHASED_TRIO_SUBOPTION)
                pr_printf(" redundant=true");

            pr_printf("\n");
            pr_printf("if ( ! -f %s.log ) then\n", file_names[5]);
            pr_printf("  echo 'First Beagle run produced no output log.'\n");
            pr_printf("  exit 1\n");
            pr_printf("endif\n");
            pr_printf("\n");
            pr_printf("if ( -z `tail -1 %s.log | grep finished` ) then\n", file_names[5]);
            pr_printf("  echo 'Beagle log file error.'\n");
            pr_printf("  exit 2\n");
            pr_printf("endif\n");

#ifdef BEAGLE_INCLUDE_PHENOTYPE
            pr_printf("\n");
            pr_printf("\necho\n");
            
            // There needs to be a trait for association testing
            if (_ttraitp != (linkage_locus_rec *)NULL) {
#define OUT_NAME                      "assoc"
                pr_printf("echo 'Association Testing...'\n");

                pr_printf("if ( -e %s.log) then\n", OUT_NAME);
                pr_printf("  mv -f %s.log %s.log.old\n", OUT_NAME, OUT_NAME);
                pr_printf("endif\n");

                // Since the source for Beagle is not available, it's not clear what
		// System.exit(N); code that it returns on a failure.
                pr_printf("$_JAVA -jar $BEAGLE_JAR data=%s.%s.phased.gz trait=%s out=%s",
                          file_names[5], file_names[0],
                          _ttraitp->Name, OUT_NAME);
                if ((*analysis)->_suboption == BEAGLE_UNPHASED_PAIR_SUBOPTION)
                    pr_printf(" diplotypes=false");
                pr_printf("\n");
            } else {
                pr_printf("echo 'No Trait available for Association Testing.'\n");
            }
#endif /* BEAGLE_INCLUDE_PHENOTYPE */
        }
        void file_post() {
            chmod_X_file(path_);
        }
    } *xp = new BEAGLE_sh_script(Top, floop, analysis);
    
    xp->file_names = file_names;
    xp->sh         = sh;

    floop->iterate();
    
    delete xp;
    delete floop;
    
    if (top_shell) {
        mssgvf("        BEAGLE top shell file:               %s/%s\n", output_paths[0], file_names[4]);
        mssgvf("               the above shell runs all shells\n");
        sh->filep_close();
        delete sh;
    }
}

void CLASS_BEAGLE::create_output_file(linkage_ped_top *LPedTreeTop,
                                      analysis_type *analysis,
                                      char *file_names[],
                                      int untyped_ped_opt,
                                      int *numchr,
                                      linkage_ped_top **Top2)
{
    linkage_ped_top *Top = LPedTreeTop;
    int pwid, fwid, mwid;
    int combine_chromo = 0;
    char prefix[100];

    // allow_sex_map() { return false; }
    if (genetic_distance_index >= 0 && genetic_distance_sex_type_map != SEX_AVERAGED_GDMT) {
        errorf("Beagle can only use sex-averaged maps if a genetic map is chosen.");
        EXIT(DATA_TYPE_ERROR);
    }
    
    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    combine_chromo = 0;
    
    get_file_names(file_names, prefix, Top->OrigIds, Top->UniqueIds, &combine_chromo);
    combine_chromo = 0;
    LoopOverChrm   = 1;
    
    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, Top->LocusTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);
    
    write_BEAGLE_marker_file(Top, file_names, pwid, fwid);

    // chosen based on suboption...
    if ((*analysis)->_suboption == BEAGLE_UNPHASED_UNRELATED_SUBOPTION)
        write_BEAGLE_genotype_unphased_unrelated_file(Top, file_names, pwid, fwid);
    else if ((*analysis)->_suboption == BEAGLE_UNPHASED_TRIO_SUBOPTION)
        write_BEAGLE_genotype_unphased_trio_file(Top, file_names, pwid, fwid);
    else if ((*analysis)->_suboption == BEAGLE_UNPHASED_PAIR_SUBOPTION)
        write_BEAGLE_genotype_unphased_pair_file(Top, file_names, pwid, fwid);

    write_BEAGLE_sh(Top, analysis, file_names);
}

void CLASS_BEAGLE::get_file_names(char *file_names[], char *prefix,
                                  int has_orig, int has_uniq, int *combine_chromo)
{
    int i, choice, igl, ipre, ibgl, ipmkr, igmkr, ish, ioui, ioup, isum, isumf;
    char fl_stat[12];
    analysis_type analysis = this;
    
    strcpy(prefix, "beagle");
    
    if (DEFAULT_OUTFILES) {
        mssgf("Output file names set to defaults.");
        choice = 0;
    } else {
        choice = -1;
    }

    if (main_chromocnt > 1) {
        // This is the batch file item that controls whether you wish to comnine the
        // chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
        // This is a derective from the user which will override the default...
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read &&
            tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) != 'y')
            warnvf("For %s analysis, batch file item 'Loop_Over_Chromosomes=y' is assumed.\n",
		 analysis->_name);
    }
    
    // For Beagle, each chromosome must get it's own file...
    *combine_chromo=0;
    analysis->replace_chr_number(file_names, global_chromo_entries[0]);
    
    /* output file name menu */

    igl = ipre = ibgl = ipmkr = igmkr = ish = ioui = ioup = isum = isumf = -1;
    // If the default output files are used we do not go here.
    // Otherwise, enter with choice -- -1.
    while (choice != 0) {
        draw_line();
        print_outfile_mssg();
        printf("Output file names menu:\n");
        printf("0) Done with this menu - please proceed\n");
        i=1;
        printf(" %d) Genotypes file name                       %-15s\t%s\n",
               i, file_names[0],
               file_status(file_names[0], fl_stat));
        ibgl=i++;
        
        if (genetic_distance_index >= 0) {
            printf(" %d) Genetic distance marker file name         %-15s\t%s\n",
                   i, file_names[2],
                   file_status(file_names[2], fl_stat));
            igmkr=i++;
        }

        if (base_pair_position_index >= 0) {
            printf(" %d) Base position marker file name            %-15s\t%s\n",
                   i, file_names[1],
                   file_status(file_names[1], fl_stat));
            ipmkr=i++;
        }

        printf(" %d) File name stem:                           %-15s\n", i, prefix);
        ipre=i++;
        
        printf(" %d) Shell file name:                          %-15s\t%s\n",
               i, file_names[3],
               ((main_chromocnt <= 1 || *combine_chromo == 1) ? file_status(file_names[3], fl_stat) : ""));
        ish=i++;

        // Here we determine the length and overlap of the genotype file...
        // See the Beagle manual 

        printf("Enter options 0-%d > ", i);
        fcmap(stdin, "%d", &choice); printf("\n");
        test_modified(choice);
        
        if (choice < 0) {
            printf("Unknown option %d\n", choice);
        } else if (choice == 0) {
            ;
        } else if (choice == igl) {
            *combine_chromo = TOGGLE(*combine_chromo);
            if (main_chromocnt > 1 && *combine_chromo) {
                // replaces <extension> with 'all', keeping <extension> and <rest> if they exist...
                analysis->replace_chr_number(file_names, 0);
            } else {
                analysis->replace_chr_number(file_names, global_chromo_entries[0]);
            }
            
        } else if (choice == ipre) {
            printf("Enter new file name stem > ");
            fcmap(stdin, "%s", prefix);    newline;
            inner_file_names(file_names, "", prefix);
            if (main_chromocnt > 1 && *combine_chromo)
                analysis->replace_chr_number(file_names, 0);
            else
                analysis->replace_chr_number(file_names, global_chromo_entries[0]);
            
        } else if (choice == ibgl) {
            printf("Enter new genotypes file name > ");
            fcmap(stdin, "%s", file_names[0]);    newline;

        } else if (choice == igmkr) {
            printf("Enter new genetic marker file name > ");
            fcmap(stdin, "%s", file_names[2]);    newline;

        } else if (choice == ipmkr) {
            printf("Enter new physical marker file name > ");
            fcmap(stdin, "%s", file_names[1]);    newline;
            
        } else if (choice == ish) {
            printf("Enter new shell script name %s > ", file_names[3]);
            fcmap(stdin, "%s", file_names[3]);    newline;
            
        } else if (choice == ioui) {
            OrigIds[0] = individual_id_item(0, analysis, OrigIds[0], 35, 1, has_orig, has_uniq);
            individual_id_item(0, analysis, OrigIds[0], 0, 3, has_orig, has_uniq);
            
        } else if (choice == ioup) {
            OrigIds[1] = pedigree_id_item(0, analysis, OrigIds[1], 35, 1, has_orig);
            pedigree_id_item(0, analysis, OrigIds[1], 0, 3, has_orig);
            
        } else {
            printf("Unknown option %d\n", choice);
        }
    }
}

static void inner_file_names(char **file_names, const char *num, const char *stem /* = "beagle" */)
{    
    sprintf(file_names[0], "%s.%s.bgl.gz", stem, num);
    if (base_pair_position_index >= 0)
        sprintf(file_names[1], "%s.%s.pmkr.gz", stem, num);
    if (genetic_distance_index >= 0) 
        sprintf(file_names[2], "%s.%s.gmkr.gz", stem, num);

/* w/o zlib
    sprintf(file_names[0], "%s.%s.bgl", stem, num);
    if (base_pair_position_index >= 0)
        sprintf(file_names[1], "%s.%s.pmkr", stem, num);
    if (genetic_distance_index >= 0) 
        sprintf(file_names[2], "%s.%s.gmkr", stem, num);
 */

    sprintf(file_names[3], "%s.%s.sh", stem, num);
    sprintf(file_names[4], "%s.all.sh", stem);
    sprintf(file_names[5], "%s.%s.out", stem, num);
    sprintf(file_names[6], "%s", stem);
}

void CLASS_BEAGLE::file_names(char **file_names, char *num)
{
    inner_file_names(file_names, num);
}

void CLASS_BEAGLE::replace_chr_number(char *file_names[], int numchr) {
    change_output_chr(file_names[0], numchr);
    if (base_pair_position_index >= 0)
        change_output_chr(file_names[1], numchr);
    if (genetic_distance_index >= 0) 
        change_output_chr(file_names[2], numchr);
    change_output_chr(file_names[3], numchr);
    change_output_chr(file_names[5], numchr);
}

static void illegal_suboption_error_exit() {
    errorf("For 'Beagle' analysis, the batch file item 'Analysis_Sub_Option' must be");
    errorf("specified as either:'Unphased_Unrelated', 'Unphased_Trio', 'Unphased_Pair'.");
    EXIT(DATA_TYPE_ERROR);
}

void CLASS_BEAGLE::interactive_sub_prog_name_to_sub_option(analysis_type *analysis)
{
    int  selection;
    char select[10];
    
    if (batchANALYSIS) {
        if ((*analysis)->_suboption < 1 || (*analysis)->_suboption > 3)
            illegal_suboption_error_exit();
        return;
    }
    
    selection = 0;
    while(!(selection == BEAGLE_UNPHASED_UNRELATED_SUBOPTION ||
            selection == BEAGLE_UNPHASED_TRIO_SUBOPTION ||
            selection == BEAGLE_UNPHASED_PAIR_SUBOPTION)) {
        draw_line();
        printf("Selection Menu: Beagle genotype file creation options:\n");
        printf(" NOTE: The structure of your input data should match the\n"); 
        printf(" the structure expected by your chosen BEAGLE option.\n");
        printf(" Please see the Mega2 documentation for details.\n");
        printf("%1d) Unphased Unrelated genotype file.\n", BEAGLE_UNPHASED_UNRELATED_SUBOPTION);
        printf("%1d) Unphased Trio genotype file.\n", BEAGLE_UNPHASED_TRIO_SUBOPTION);
        printf("%1d) Unphased Pair genotype file.\n", BEAGLE_UNPHASED_PAIR_SUBOPTION);
        printf("Enter selection: 1 - 3 > ");
        fcmap(stdin,"%s", select); newline;
        selection=0;
        sscanf(select, "%d", &selection);
        if (selection < 1 || selection > 3) warn_unknown(select);
    }
    
    // These should be the same as the '_suboption' in the class definitions...
    if (selection == BEAGLE_UNPHASED_UNRELATED_SUBOPTION) {
        *analysis = BEAGLE_UNPHASED_UNRELATED;
    } else if (selection  == BEAGLE_UNPHASED_TRIO_SUBOPTION) {
        *analysis = BEAGLE_UNPHASED_TRIO;
    } else if (selection == BEAGLE_UNPHASED_PAIR_SUBOPTION ) {
        *analysis = BEAGLE_UNPHASED_PAIR;
    }
}

void CLASS_BEAGLE::sub_prog_name_to_sub_option(char *sub_prog_name, analysis_type *analysis) {
    // This should match _suboption...
    switch(tolower((unsigned char)sub_prog_name[9])) {
    case 'u':
    case 'U':
        *analysis = BEAGLE_UNPHASED_UNRELATED;
        break;
    case 't':
    case 'T':
        *analysis = BEAGLE_UNPHASED_TRIO;
        break;
    case 'p':
    case 'P':
        *analysis = BEAGLE_UNPHASED_PAIR;
        break;
    default:
        illegal_suboption_error_exit();
        break;
    }
}
