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
#include <ctype.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"
#include "tod.hh"

#include "loop.h"
#include "sh_util.h"

#include "create_summary_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "linkage_ext.h"
#include "omit_ped_ext.h"
#include "output_file_names_ext.h"
#include "output_routines_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

/*
     create_summary_ext.h:  aff_status_entry marker_typing_summary
     error_messages_ext.h:  mssgf my_calloc warnf
              fcmap_ext.h:  fcmap
            linkage_ext.h:  get_loci_on_chromosome get_loci_on_chromosomes get_unmapped_loci
           omit_ped_ext.h:  omit_peds
  output_file_names_ext.h:  CHR_STR change_output_chr create_mssg file_status print_outfile_mssg
    output_routines_ext.h:  create_formats field_widths
         user_input_ext.h:  individual_id_item pedigree_id_item test_modified
              utils_ext.h:  EXIT draw_line script_time_stamp summary_time_stamp
*/

#include "plink_core_ext.h"

void write_R_PLINK_map_file(linkage_ped_top *LPTop, char *map_file_name);

void CLASS_PLINK_CORE::replace_chr_number(char *file_names[], int numchr) {
    /* PLINK: map, lgen, shell, summary */
    /* PLINK_BINARY: bim, bed, shell, summary */
    if (_suboption == PLINK_SUB_OPTION_PED_INT)
        change_output_chr(file_names[0], numchr); //othereise it remains .all. for a .fam file
    change_output_chr(file_names[1], numchr);
    change_output_chr(file_names[3], numchr);

    change_output_chr(file_names[5], numchr);
    change_output_chr(file_names[6], numchr);
    change_output_chr(file_names[7], numchr);
    change_output_chr(file_names[8], numchr);
}

static void save_PLINK_lgen(const char *genofl_name, linkage_ped_top *Top,
                                const int pwid, const int fwid, const int mwid)
{
    vlpCLASS(plink_core_lgen,chr,loci_ped_per) {
     vlpCTOR(plink_core_lgen,chr,loci_ped_per) { }

        void file_loop() {
            mssgvf("        PLINK lgen file:           %s/%s\n", *_opath, Outfile_Names[3]);
            data_loop(*_opath, Outfile_Names[3], "w");
        }
        void inner() {
            if (!_allele1 || !_allele2) return;
            pr_id();
            pr_marker_name();
            if (_LTop->Marker[_locus].Props.Numbered.Recoded) {
                // How we represent the allele in the output...
                // display the allele name according to the value associated with this person...
                pr_printf("%2s %2s",
                          recode_name(_allele1, "0", "?"),
                          recode_name(_allele2, "0", "?"));
            } else {
                // an allele with no name....
                pr_printf("%2d %2d", _allele1, _allele2);
            }
            pr_nl();
        }
    } *sp = new plink_core_lgen(Top);

    sp->load_formats(fwid, pwid, mwid);
    sp->iterate();

    delete sp;
}

static void save_PLINK_pheno(const char *phenofl_name, linkage_ped_top *Top,
                             const int pwid, const int fwid)
{
    vlpCLASS(plink_core_pheno,once,ped_per_trait) {
     vlpCTOR(plink_core_pheno,once,ped_per_trait) { }

        void file_loop() {
            msgvf("        PLINK phenotype file:      %s/%s\n", *_opath, Outfile_Names[2]);
            data_loop(*_opath, Outfile_Names[2], "w");
        }
        void file_header() {
            int tr;
            /* print the file header */
            pr_printf("FID IID ");
            for (tr = 0; tr < num_traits; tr++) {
                // If this is a marker, skip it...
                if (global_trait_entries[tr] == -1) continue;
                pr_printf(" %s",
                          _Top->LocusTop->Pheno[global_trait_entries[tr]].TraitName);
            }
            pr_nl();
        }
        void per_start() { pr_id(); }
        void per_end()   { pr_nl(); }
        void inner()     { pr_pheno(); }
    } *sp = new plink_core_pheno(Top);

    sp->load_formats(fwid, pwid, -1);

    sp->iterate();

    delete sp;
}

void CLASS_PLINK_CORE::save_pheno_file(linkage_ped_top *Top,
				       const int pwid,
				       const int fwid)
{
  save_PLINK_pheno(Outfile_Names[2], Top, pwid, fwid);
}

// MAP_MISSING is used for an unknown or missing map value...
#define MAP_MISSING               -1.0

//
// This function writes the data associated with a PLINK .MAP file when 'generate_bim_file == 0'.
// The map file contins the following column: Chromoome, SNP, Genetic Distance, Physical Distance.
// It is also used to write the PLINK .BIM file when 'generate_bim_file != 0'.
// The .BIM file contians two additional columns: ALLELE1 name, ALLELE2 name.
//
// LPTop                   the linkage ped top
// map_fp                  an open file descriptor to the map file. The routine will start
//                         writing the map data here.
// use_selected_markers_p  a boolean, write selected markers (== 1), write all markers (== 0)
// morgans_p               a boolean; output M/cM (1 == M; 0 == cM)
// three_column_map_p      do not write the genetic information to the plink map file.
// format_str              a format string to use when there is a genetic map
//                         (should not be terminated with a new line)
// generate_bim_file       true (!= 0) if we are writing a .BIM file, false (== 0) for a .MAP file
static void write_PLINK_map_data(linkage_ped_top *LPTop,
                                 FILE *map_fp,
                                 const int use_selected_markers_p,
                                 const int morgans_p,
                                 const int three_column_map_p,
                                 const char *format_str,
                                 const int generate_bim_file) {
    linkage_locus_top *LTop = LPTop->LocusTop;
    ext_linkage_locus_top *EXLTop = LPTop->EXLTop;
    // so we display certain error messages only once...
    int display_error_sagm = 0, display_error_ssgm = 0, display_error_fgm = 0, display_error_mallele = 0;
    int m, i, num_loci;
    
    // Use selected markers, or use all markers...
    num_loci = (use_selected_markers_p) ? NumChrLoci : LTop->LocusCnt;
    for (i=0; i < num_loci; i++) {
        // Use selected markers, or use all markers...
        m = (use_selected_markers_p) ? ChrLoci[i] : i;
        
        // Ignore anything but a marker.
        // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
        if (LTop->Locus[m].Class == MARKER) {
//          int LType = LTop->Locus[m].Type;
            int chr = LTop->Marker[m].chromosome;
            char *snp = LTop->Locus[m].LocusName;
            // These are the defaults (e.g. nothing was specified)...
            double genetic_distance = MAP_MISSING, base_pair_position = MAP_MISSING;
            
            if (use_selected_markers_p == 0 && chr == MISSING_CHROMO) continue;

            // A PLINK map file can only have 2 alleles, so if there are more we cannot use this marker...
            if (LTop->Locus[m].AlleleCnt > 2) {
                if (display_error_mallele == 0) {
                    warnvf("More than two alleles exist beginning at marker %s.\n", snp);
                    display_error_mallele++;
                }
                continue;
            }
            
            //
            // Credit goes to Robert for finding this interesting feature in preceeding
            // logic of the code! In the instance when there is no informaiton in the map
            // file for a marker then the vector positions pointer will be null (it never
            // get's CALLOCed; See annotated_ped_file.c copy_exmap_locmap()). We use this
            // as a flag to determine if we should output anything for the marker...
            if (EXLTop->EXLocus[m].positions == (double*)NULL) continue;
            
            if (genetic_distance_index >= 0) {
                // The routines that build the EXLTop structure should be making sure that
                // sex maps can only be of the form <a>, <m,f>, <f>.
                if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
                    // we were told to only use the average map, or only an average map was specified...
                    genetic_distance = EXLTop->EXLocus[m].positions[genetic_distance_index];
                    // But tell the user about it if used for a sex linked chromosome...
//                  if ((LType == XLINKED || LType == YLINKED) && display_error_sagm == 0)
                    if ((chr == SEX_CHROMOSOME || chr == MALE_CHROMOSOME) && display_error_sagm == 0) {
                        warnvf("Since only a sex-averaged genetic map was used, these positions have been used for markers (%d, %s) on the X and Y chromosomes.\n", chr, snp);
                        display_error_sagm++;
                    }
                } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
                    // Male, Female are available, so pick the right one given the current locus type...
//                  if (LType == XLINKED) // X or SEX_CHROMOSOME
                    if (chr == SEX_CHROMOSOME) // X or SEX_CHROMOSOME
                        genetic_distance = EXLTop->EXLocus[m].pos_female[genetic_distance_index];
//                  else if (LType == YLINKED) // Y or MALE_CHROMOSOME
                    else if (chr == MALE_CHROMOSOME) // Y or MALE_CHROMOSOME
                        genetic_distance = EXLTop->EXLocus[m].pos_male[genetic_distance_index];
                    else if (display_error_ssgm == 0) {
                        // autozomes will be set to missing is a result of the decision made to allow the user
                        // to select a SEX_SPECIFIC_GDMT
                        warnvf("Since a sex specific genetic map was used, only analysis on the X & Y chromosome is possible.\n");
                        display_error_ssgm++;
                    }
                } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
//                  if (LType == XLINKED) // X or SEX_CHROMOSOME
                    if (chr == SEX_CHROMOSOME) // X or SEX_CHROMOSOME
                        genetic_distance = EXLTop->EXLocus[m].pos_female[genetic_distance_index];
                    else if (display_error_fgm == 0) {
                        warnvf("Since only a female genetic map was used, only analysis on the X chromosome is possible.\n");
                        display_error_fgm++;
                    }
                }
            } // if (genetic_distance_index >= 0) {
            
            // If there is a physical map, then get the value from the sex-averaged spot..
            // This is a fundamental hack of the system...
            if (base_pair_position_index >= 0)
                base_pair_position = EXLTop->EXLocus[m].positions[base_pair_position_index];
            
            // The internal representation is in cM...
            if (genetic_distance == MAP_MISSING) genetic_distance = 0.0;
            else if (morgans_p == 1) genetic_distance /= 100.0;
            
            // This was a '1' but was changed so that it looks the same as genetic_distance
            // In a PLINK map file this would be a '1'.
            if (base_pair_position == MAP_MISSING) base_pair_position = 0.0;
            
            // PLINK uses these default values...
            // PLINK allows the chromosome value of '0' to represent "unplaced"
            if (chr == UNKNOWN_CHROMO) chr = 0;
            
	    // NOTE: if generate_bim_file == true then, three_column_map_p must be false...
            if (three_column_map_p) {
                // just print three columns excluding the genetic_distance....
#if defined(_WIN) || defined(MINGW)
                fprintf(map_fp, "%d\t%s\t%.0f\n", chr, snp, base_pair_position);
#else
                fprintf(map_fp, "%d\t%s\t%.0lf\n", chr, snp, base_pair_position);
#endif
            } else {
                // Print what we have for a genetic distance in the format that we were told to use...
                fprintf(map_fp, format_str, chr, snp, genetic_distance, base_pair_position);
                
                // A .BIM file must have six columns (e.g., start with a four column map file).
                // The additional columns are the names of the alleles...
                if (generate_bim_file != 0) {
                    // When PLINK makes the .BIM file, the minor allele is listed before major allele.
                    // We do not list alleles based on frequency...
                    fprintf(map_fp, "\t%s\t%s",
                            LTop->Locus[m].Allele[0].AlleleName, LTop->Locus[m].Allele[1].AlleleName);
                }
                
                fprintf(map_fp, "\n");
            }
            
        } // if (LTop->Locus[m].Class == MARKER) {
    } //  for (m=0; m < LTop->LocusCnt; m++) {
}

//
// Write some number of 4-column PLINK map files that use selected markers,
// will have comments, will NOT be annotated, and outputs genetic distance in Morgans.
// If generate_bim_file is true (!= 0) output a .bim file, otherwise a .map file
static void write_PLINK_map(linkage_ped_top *LPTop,
                            const char *mfl_name,
                            const int generate_bim_file) {
#ifdef PLINK_MAP_FILE_COMMENTS
    ext_linkage_locus_top *EXLTop = LPTop->EXLTop;
#endif /* PLINK_MAP_FILE_COMMENTS */        

    vlpCLASS(plink_core_map,chr,null) {
     vlpCTOR(plink_core_map,chr,null) { }
        int generate_bim_file;

        void file_loop() {
            mssgvf("        PLINK map file:            %s/%s\n", *_opath, Outfile_Names[1]);
            data_loop(*_opath, Outfile_Names[1], "w");
        }
        void file_header() {
            // It seems that PLINK map (.BIM) files cannot handle comments...
#ifdef PLINK_MAP_FILE_COMMENTS
            // Insert comments by default into the PLINK map file telling about the maps that are used to produce it....
            if (genetic_distance_index >= 0) {
                pr_printf("# Genetic Map '%s' is sex-%s; of type: %s; units: Morgans.\n",
                          EXLTop->MapNames[genetic_distance_index],
                          (genetic_distance_sex_type_map == 0 ? "average" : "specific"),
                          (EXLTop->map_functions[genetic_distance_index] == 'k' ? "Kosambi" : "Haldane")
                    );
            } else {
                pr_printf("# No Genetic Map was specified in the input map file. Using values of 0.\n");
            }
            if (base_pair_position_index >= 0) {
                pr_printf("# Physical Map used: %s.%c\n",
                          EXLTop->MapNames[base_pair_position_index],
                          EXLTop->map_functions[base_pair_position_index]
                    );
            } else {
                pr_printf("# No Physical Map was specified in the input map file. Using values of 0.\n");
            }

            if (generate_bim_file != 0) {
                pr_printf("# This is a .BIM file, the last two columns are the names of allele 1, and allele 2.\n");
            }
#endif /* PLINK_MAP_FILE_COMMENTS */        
            // This file is NOT annotated...
            // Use selected markers, output Morgans, and create a 4-column map file...
        }
        void inner() {
            // Loop through the Loci...
            markers_on_chromosome(_numchr);

#if defined(_WIN) || defined(MINGW)
            write_PLINK_map_data(_Top, _filep, 1, 1, 0, "%d\t%s\t%9.7f\t%.0f", generate_bim_file);
#else
            write_PLINK_map_data(_Top, _filep, 1, 1, 0, "%d\t%s\t%9.7lf\t%.0lf", generate_bim_file);
#endif

        }
    } *sp = new plink_core_map(LPTop);

    sp->generate_bim_file = generate_bim_file;

    sp->iterate();

    delete sp;
}

void  create_PLINK_files(linkage_ped_top **LPedTop,
                         char *file_names[],
                         const int untyped_ped_opt,
                         const int output_format,
                         const char *prefix,
                         analysis_type *analysis)
{
    int numchr=0, pwid, fwid, mwid;
    int combine_chromo = 0;
    int create_geno_summary = 0;
    linkage_ped_top *Top = *LPedTop;
    linkage_locus_top *LTop = Top->LocusTop;
    
    extern int base_pair_position_index;

    // if 'combine_chromo == 0' each chromosome gets it's own file.
    // if 'combine_chromo == 1' all informaiton goes into one file with the '.all' suffix.
    // For the .BED format modes (e.g., output_format > 0), the default should be to
    // write everything to one file (unless only one chromosome has been selected).
    // The user can override this behavior in the menu created by 'PLINK_file_names'
    // (if more than one chromosome has been specified) by changing the menu item
    // 'Combine chromosomes?'
    combine_chromo = (output_format > 0 && main_chromocnt > 1) ? 1 : 0;
    
    if (main_chromocnt > 1) {
        // This is the batch file item that controls whether you wish to comnine the
        // chromosomes in the same file or not. If true (y), each chromosome gets it's own file.
        // This is a derective from the user which will override the default...
        // NOTE that for interactive input, reorder_loci.cpp:ReOrderLoci() will ask the user
        // if they want to write chromosomes to one file or one file for each cromosome IF the user selects
        // multiple chromosomes in the preceeding menu. The answer will be stored in "Loop_Over_Chromosomes".
        if (Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].items_read)
            combine_chromo = (tolower((unsigned char)Mega2BatchItems[/* 50 */ Loop_Over_Chromosomes].value.copt) == 'y') ? 0 : 1;
    }

    // NOTE: comnine_chromo is passed in because the user can change it via a menu.
    if (InputMode == INTERACTIVE_INPUTMODE) {
        (*analysis)->user_queries(file_names, &combine_chromo, &create_geno_summary);
        LoopOverChrm = ! combine_chromo;
    } else
        (*analysis)->batch_in();

//  file_stem may have been reset
    if (! DEFAULT_OUTFILES) {
	(*((CLASS_PLINK_CORE **) analysis))->file_names_w_stem(file_names, (char *)"xx", (*analysis)->file_name_stem, output_format+1);
	add_sumdirs(file_names);
    }


    (*analysis)->batch_show();

    // http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map
    // in PLINK a 3 column map is one that is missing the Genetic Distance column.
    // That is it has: chromosome, Name (snp), Base-pair position (Physical Map).
    if (genetic_distance_index < 0) {
        warnf("No Genetic Map was specified.\nGenetic Map will not be written.\n");
    }
    
    // the analysis parameter is not used by this function...
    create_mssg(*analysis);
    
    // determine the maximum field widths necessary to output certain data...
    field_widths(Top, LTop, &fwid, &pwid, NULL, &mwid);
    
    // Omit pedigrees under certain circumstances...
    omit_peds(untyped_ped_opt, Top);

    // either the 1 chosen chromosome # or .all
    if (main_chromocnt > 1)
            change_output_chr(file_names[0], 0);
    else
            change_output_chr(file_names[0], global_chromo_entries[0]);

    if (num_traits > 1)
        ((plink_analysis_type)*analysis)->save_pheno_file(Top, pwid, fwid);

    // not annotated, genetic distance in Morgans, include comments...
    // generate a .bim file if output_format != 0 otherwise a .map file if output_format == 0
    Tod tod_plmap("plink: write map file");
    write_PLINK_map(Top, file_names[1], output_format == 1 || output_format == 2);
    tod_plmap();

    Tod tod_plfmt("plink: write files: lgen/ped/ped6+bed");
    if (output_format == 0) {
        ((plink_analysis_type)*analysis)->save_pedsix_file(Top, pwid, fwid);
        save_PLINK_lgen(file_names[3], Top, pwid, fwid, mwid);
    } else if (output_format == 3) {
        ((plink_analysis_type)*analysis)->save_ped_file(Top, pwid, fwid, mwid);
    } else {
        ((plink_analysis_type)*analysis)->save_pedsix_file(Top, pwid, fwid);
        // the value of output_format determines the .BED mode
        // (SNP major == 1; Individual major == 2).
        ((plink_analysis_type)*analysis)->save_bed_file(file_names[3], Top, output_format);
    }
    tod_plfmt();

    // Now a method...
    //write_PLINK_sh(numchr, combine_chromo, file_names, &(prefix[0]), Top, output_format);

    if (create_geno_summary) {
        FILE *filep;
        if ((filep=fopen(file_names[5], "w")) == NULL) {
            errorvf("Unable to open %s for writing.\n", file_names[5]);
            EXIT(FILE_WRITE_ERROR);
        }
        summary_time_stamp(mega2_input_files, filep,
                           "IDs are from :Ped, :Per fields if present, or pedigree, person ids.");
        marker_typing_summary(filep, Top, numchr);
        mssgvf("   Genotyping summary file:            %s\n", file_names[5]);
    }

    if (base_pair_position_index < 0) {
        warnf("No physical distances given, setting map positions to 0.");
    }
}


//
// Write the plink map file for consumption by R from internal structures.
// LPTop                     the linkage ped top
// map_file_name             file name with complete path
//
// This routine also uses some globals set in user_input.c: get_genetic_distance_index and
// get_base_pair_position_index. The globals are: genetic_distance_index, genetic_distance_sex_type_map,
// and base_pair_position_index. These variables can also be set in the batch file...
//
// NOTE: We can't call 'LPTop' 'Top' as is the custom, because above us
// '#define Top ...' wonderful code isn't it...
void write_R_PLINK_map_file(linkage_ped_top *LPTop, char *map_file_name) {
    ext_linkage_locus_top *EXLTop = LPTop->EXLTop;
    FILE *fp;
    
    // The user chooses the base_pair_position_index in the routine
    // user_input.c:get_base_pair_position_index(ext_linkage_locus_top *EXLTop)
    // which is run in mega2.c:main() before this function is called...
    // IMPORTANT NOTE: DO NOT REMOVE THIS CHECK THE CODE BELOW DEPENDS ON A VALID BASE_PAIR_POSITION_INDEX
    if (base_pair_position_index < 0) {
        errorvf("For the analysis type specified, a Physical Map is required, but none was chosen.\n");
        EXIT(INPUT_DATA_ERROR);
    }
    
    if ((fp = fopen(map_file_name, "w")) == (FILE *)NULL) {
        errorvf("Unable to open %s for writing.\n", map_file_name);
        EXIT(FILE_WRITE_ERROR);
    }
    
    // This is an annotated file, with NO comments...
    fprintf(fp, "Chromosome\tName\t%s.%c\n",
            EXLTop->MapNames[base_pair_position_index],
            EXLTop->map_functions[base_pair_position_index]
            );
    
    // Use all markers, output centi-Morgans, create a 3-column map file, generate a .BIM file...
#if defined(_WIN) || defined(MINGW)
    write_PLINK_map_data(LPTop, fp, 0, 0, 1, "%d\t%s\t%7.5f\t%.0f", 0);
#else
    write_PLINK_map_data(LPTop, fp, 0, 0, 1, "%d\t%s\t%7.5lf\t%.0lf", 0);
#endif
    fclose(fp);
}


// the value of output_format determines how the files are written:
// .LGEN mode (== 0),
// .BED mode: (SNP major == 1; Individual major == 2).
void CLASS_PLINK_CORE::user_queries(char **file_names,
				    int *combine_chromo, int *create_summary)
{
    int i, choice = -1, icc=-1, isum=-1, istem=-1;
    
    while (choice != 0) {
        print_outfile_mssg();
        draw_line();
        printf("0) Done with this menu - please proceed\n");
        i=1;
        
        if (main_chromocnt > 1) {
            printf(" %d) Combine chromosomes?                      %s\n",
                   i, yorn[*combine_chromo]);
            icc=i++;
        }
        
        printf(" %d) Create genotyping summary?                %s\n",
               i, yorn[*create_summary]);
        isum=i++;
        
        printf(" %d) Change file names stem?                   %s\n",
               i, this->file_name_stem);
        istem=i++;
        
        printf("Enter options 0-%d > ", i-1);
        fcmap(stdin, "%d", &choice); printf("\n");
        
        if (choice == 0) {
            ; // OK...
        } else if (choice == icc) {
            *combine_chromo = TOGGLE(*combine_chromo);
        } else if (choice == isum) {
            *create_summary = TOGGLE(*create_summary);
        } else if (choice == istem) {
            printf("Enter new stem for the output file names > ");
            fcmap(stdin, "%s", this->file_name_stem);    newline;
            // It doesn't matter what the parameter 'num' in the method file_names() is. It will get
            // changed to the appropriate thing later in the code. The method should be rewritten
            // globally without num and a place holder inserted instead.
            this->gen_file_names(file_names, (char *)"xx");
        } else {
            printf("Unknown option %d\n", choice);
        }
    }
}
