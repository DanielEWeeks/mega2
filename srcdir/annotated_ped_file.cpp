/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

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
/* provides:
   read_annotated_ped_file()
   write_annotated_ped_file()
   create_names_file()
   read_frequency_file()
   write_frequency_file()
*/

/* read in an annotated pedigree file with headers. If the pedigree file
   has a header line containing some recognized keywords, then recognize
   it as an annotated file.
   Create the names file from the header.
   If a frequency file is provided, then read in the frequencies as well.
*/

/* To do:
   Need to read in the sic columns of penetrance file, then store the
   penetrances properly.

   Have a way of indexing correctly into the trait_list

   Copy the penetrances over into the Linkage locus structure,
   using functions similar to count_raw_alleles, and recode_linkage_locus_top

*/


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include "common.h"
#include "typedefs.h"
#include "tod.hh"
#include "mrecode.h"

#include "annotated_ped_file.h"
#include "annotated_ped_file_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "linkage.h"
#include "list_ext.h"
#include "makeped_ext.h"
#include "marker_lookup_ext.h"
#include "mrecode_ext.h"
#include "omit_ped_ext.h"
#include "phe_lookup_ext.h"
#include "read_files_ext.h"
#include "reorder_loci_ext.h"
#include "select_individuals_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "plink_ext.h"

#include "str_utils.hh"
#include "read_impute.hh"

#include "class_old.h"

/*
     error_messages_ext.h:  errorf mssgf my_calloc my_malloc warnf
              fcmap_ext.h:  fcmap
        grow_string_ext.h:  grow
            linkage_ext.h:  clear_lpedrec clear_lpedtree get_llooprec new_lpedtop
        linkage_ext_ext.h:  new_ex_llocustop
               list_ext.h:  append_to_list_tail new_list pop_first_list_entry
            makeped_ext.h:  copy_pedrec_data
      marker_lookup_ext.h:  make_marker search_marker
            mrecode_ext.h:  assign_dummy_alleles convert_to_freq count_classes count_raw_alleles create_allele_list default_trait_penetrances free_marker_item need_recoding read_marker_data read_marker_only_data recode_locus_top recode_ped_top write_recode_summary read_m2_map_as_names_file
           omit_ped_ext.h:  omit_peds
         read_files_ext.h:  read_aff_phen read_numbered_data read_quant_phen untype_locus
       reorder_loci_ext.h:  get_chromosome_list
 select_individuals_ext.h:  get_count_option
         user_input_ext.h:  set_missing_quant_input
              utils_ext.h:  EXIT columncount log_line
        write_files_ext.h:  write_locus_stats write_ped_stats
*/

col_hdr_type    ReservedColnames[NUM_PEDCOL_NAMES];

std::vector<Vecc> VecAlleles;

Mapci ped2idx;

#if 0
VeccDB PedName;
VeccDB PedIdName;
VeccDB PedLinkIdName;
VeccDB PerName;
VeccDB PerIdName;
VeccDB PerLinkIdName;
VeccDB PerUniqName;
#endif



static linkage_ped_top *read_common_ped_file(FILE *filep, char *pedfile,
                                             plink_info_type *plink_info,
                                             linkage_locus_top *LTop,
                                             annotated_file_desc *file_desc,
                                             int phecols,
                                             int num_groups,
                                             int *groups,
                                             int num_ped_records,
                                             int has_extra_ids);

annotated_file_desc AnnotatedFileInfo;

#ifdef DEFUNCT
static void ann_field_widths(linkage_ped_top *TTop,
			     annotated_file_desc *file_desc);
#endif
static void init_reserved_pedcol_names(void);

static int parse_variable_width_hdr(FILE *file,
				    col_hdr_type *reserved_colnames,
				    listhandle *userdef_ped_cols,
				    const int num_reserved_cols,
				    const int required_cols,
				    const char **valid_header_item_extensions);

static int read_plink_map_as_names_file(char *map_file, linkage_locus_top **LTop,
                                        int cols, char **phe_names, int *phe_types,
                                        annotated_file_desc *file_desc);

static ext_linkage_locus_top *read_plink_map_file(const char *map_file,
                                                  linkage_locus_top *LTop,
                                                  annotated_file_desc *file_desc,
                                                  plink_info_type *plink_info);

static linkage_ped_top *read_plink_ped_file(char *pedfile,
                                            char *bedfile,
                                            plink_info_type *plink_info,
                                            int phecols,
                                            linkage_locus_top *LTop,
                                            ext_linkage_locus_top *EXLTop,
                                            int bp_map,
                                            annotated_file_desc *file_desc,
                                            int num_groups,
                                            int *groups);

static int parse_phe_types(char *phe_file, char ***phe_names, int **phe_types);

/* PLINK parameters */
#include "plink_ext.h"
/* 
    plink, no_fid, no_parents, no_pheno, map3, 
    cM, missing_pheno, geneticMapType, pheno_value, trait
*/
PLINK_t PLINK = { 0, not_plink_format, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"};

/* end PLINK parameters */

// //////////////////////////////////////////////////////////////// //



// full header check routines...
static int is_missing_required_hdr(col_hdr_type *col_names,
                                   int num_elements);
static void error_hdr(char *hdr,
                      int num_known,
                      col_hdr_type *known_hdrs);


// 7.2.2. Annotated pedigree file (column headers)
//
// Phenotype and genotype column headers are defined by the user using the following convention: numbered
// autosomal (.M) and x-linked (.X|.Y) loci genotype column names are named as 'marker-name.marker-type.extension',
// where 'marker-name' one of the locus names within the names file, 'marker-type' is the corresponding locus
// type, and the 'extension' 1 or 2 denotes the allele-number. Currently, the allele-number is of no significance
// with respect to the haplotype, in future, we will consider extending Mega2 to process haplotypes.
//
// Binary traits are given an 'extension' A for the status column header, and L for the liability class header.
// Traits with a single liability class does not need this liability column. Quantitative traits and covariates have
// either a Q or a C as the 'extension'.
// For an affection trait, the headers should be as follows:
// a. If there is a single liability class, then you only need one phenotype column containing the affected status,
//    and the header for this should have the single extension A .
// b. If there are multiple liability classes, then the status column header should have 2 extensions L.1, and the
//    liability class column should have the extensions L.2.

// These are the valid extensions for header items for the particular types of files.
// They are used by the routine 'is_valid_hdr_ext()', and case is not important.
// A null pointer marks the end of the list...
// NOTE: Since the remainder of this code uses .A.1 and .A.2 for .L.1 and .L.2 we
// are going to change the manual to match the code.
static const char *pedigree_valid_header_item_extensions[] = {
    ".T", ".C", ".A", ".M.1", ".M.2", ".X.1", ".X.2", ".Y.1", ".Y.2", ".A.1", ".A.2", ".L.1", ".L.2", ".T.1", ".T.2", NULL
};

static const char *map_valid_header_item_extensions[] = {
  ".K.A", ".K.M", ".K.F", ".H.A", ".H.M", ".H.F", ".P", NULL
};

// header extension check routines...
static int is_valid_hdr_ext(const char *colname,
                            const char **valid_header_item_extensions);
static void error_hdr_ext(const char *colname,
                          const int column_index,
                          const char **valid_header_item_extensions);

// divide the header into it's component parts...
static void parse_hdr(const char *colname,
                      char *name, char *ext1, char *ext2,
                      const char **valid_header_item_extensions);



static void skip_long_line(FILE *fp, int c)
{
    while(c != '\n') {
        c=fgetc(fp);
    }
}

Mapcc canonicalName;
const char *canonicalColName(const char *p)
{
    char *r;
    const char *a;
    if (map_get(canonicalName, p, a)) {
    } else {
        r = new char[strlen(p)+1];
        strcpy(r, p);
        a = r;
        canonicalName[a] = a;
    }
    return a;
}
#define SAME(a,b) (a == b)

struct _TKN { 
    const char *MT;
    const char *Zero;
    const char *NA;

    const char *Pedigree;
    const char *ID;
    const char *Father;
    const char *Mother;
    const char *Sex;
    const char *PedID;
    const char *LinkPedID;
    const char *PerID;
    const char *LinkPerID;
    const char *FirstOff;
    const char *NextMatSib;
    const char *NextPatSib;
    const char *MZTwin;
    const char *DZTwin;
    const char *Proband;
    const char *Group;
} TKN;

static void init_tokens()
{
    TKN.MT   = canonicalColName("");
    TKN.Zero = canonicalColName("0");
    TKN.NA   = canonicalColName("NA");

    TKN.Pedigree = canonicalColName("Pedigree");
    TKN.ID = canonicalColName("ID");
    TKN.Father = canonicalColName("Father");
    TKN.Mother = canonicalColName("Mother");
    TKN.Sex = canonicalColName("Sex");
    TKN.PedID = canonicalColName("PedID");
    TKN.LinkPedID = canonicalColName("LinkPedID");
    TKN.PerID = canonicalColName("PerID");
    TKN.LinkPerID = canonicalColName("LinkPerID");
    TKN.FirstOff = canonicalColName("FirstOff");
    TKN.NextMatSib = canonicalColName("NextMatSib");
    TKN.NextPatSib = canonicalColName("NextPatSib");
    TKN.MZTwin = canonicalColName("MZTwin");
    TKN.DZTwin = canonicalColName("DZTwin");
    TKN.Proband = canonicalColName("Proband");
    TKN.Group = canonicalColName("Group");
}

static void copy_colname(col_hdr_type *from, col_hdr_type *to)
{
//  strcpy(to->ColName, from->ColName);
    to->ColName        = from->ColName;
    to->value_type     = from->value_type;
    to->input_col      = from->input_col;
    to->output_length  = from->output_length;
    to->output_col     = from->output_col;
    to->locus_number   = from->locus_number;
    to->map_number     = from->map_number;
    to->sex_map_number = from->sex_map_number;
}

/* reorder the file_descriptor for pedigree columns , so that the
   colnum headers are in the same order as they appear in the pedigree
   file.
   Since the pedigree file is likely to have a large number of columns,
   this is desirable, so that we don't have to look for the correct entry
   in the pedigree file descriptor list.
   Also check each column header that is not a reserved header against the
   names file, and flag those that are not in the names file. Signal an error
   if this happens.
*/
static void loctype_to_descriptor(linkage_locus_top *LTop, int i, char loctype[])
{
    linkage_locus_rec Locus = LTop->Locus[i];
    linkage_locus_type Type = Locus.Type;

    switch(Type) {
    case QUANT:
        if (Locus.Class == COVARIATE) {
            loctype[0]=loctype[1]='C';
        } else {
            loctype[0]=loctype[1]='T';
        }
        break;
    case AFFECTION:
        if (LTop->Pheno[i].Props.Affection.ClassCnt != 1) {
            loctype[0]='A'; loctype[1]='L';
        } else {
            loctype[0]=loctype[1]='A';
        }
        break;
    case BINARY:
    case NUMBERED:
        loctype[0]=loctype[1]= 'M';
        break;
    case XLINKED:
/*      loctype[0]=loctype[1]= 'X'; */
        loctype[0]= 'X';
        loctype[1]= 'M';
        break;
    case YLINKED:
/*      loctype[0]=loctype[1]= 'Y'; */
        loctype[0]= 'Y';
        loctype[1]= 'M';
        break;
    default:
        loctype[0]=loctype[1]= 'U';
        break;
    }
    return;
}

static col_hdr_type *pedcol_check_reverse_index(col_hdr_type *reserved_colnames,
                                                int num_userdef_cols,
                                                listhandle *userdef_colnames,
                                                linkage_locus_top *LTop,
                                                int *num_reserved_cols,
                                                int *has_extra_ids)
{
    /* Store the column headers in order of their appearance in the pedigree file
       to make the reading in easier */

    int i, phenotype_col=0, output_col=0;
    col_hdr_type *colname_item;
    col_hdr_type *ped_all_colnames;
    int locus_found, inconsistent_type=0, name_not_found=0;
    char main_part[FILENAME_LENGTH], ext1[FILENAME_LENGTH], ext2[FILENAME_LENGTH], loctype[2];

    /* time_t t1, t2; */

    *num_reserved_cols=0;
    *has_extra_ids = 0;
    pedfile_type = PREMAKEPED_PFT;
    basefile_type = pedfile_type;

    ped_all_colnames = CALLOC((size_t)num_userdef_cols, col_hdr_type);
#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: ped_all_colnames: %d MB (%d x %d)\n",
          num_userdef_cols * sizeof(col_hdr_type) / 1024 / 1024,
          num_userdef_cols,  sizeof(col_hdr_type));
#endif

    for (i=0; i < NUM_PEDCOL_NAMES; i++) {
        if (reserved_colnames[i].input_col > -1) {
            colname_item = &(ped_all_colnames[reserved_colnames[i].input_col]);
            copy_colname(&(reserved_colnames[i]), colname_item);
/*       printf("from %d, to %d: Name %s\n", reserved_colnames[i].input_col, */
/* 	     colname_item->input_col, colname_item->ColName); */
            colname_item->output_col = output_col;
            (*num_reserved_cols)++;
            output_col++;
            /*      printf("%s\n", reserved_colnames[i].ColName); */
            if (SAME(reserved_colnames[i].ColName, TKN.FirstOff)) {
                pedfile_type = POSTMAKEPED_PFT;
                basefile_type = pedfile_type;
            }
            if (SAME(reserved_colnames[i].ColName, TKN.LinkPedID)) {
                *has_extra_ids |= 2;
            }
            if (SAME(reserved_colnames[i].ColName, TKN.LinkPerID)) {
                *has_extra_ids |= 1;
            }
            if (SAME(reserved_colnames[i].ColName, TKN.PedID)) {
                *has_extra_ids |= 8;
            }
            if (SAME(reserved_colnames[i].ColName, TKN.PerID)) {
                *has_extra_ids |= 4;
            }
        }
    }

    /* for the user-defined columns, we don't need to handle IGNORE type
       columns because the output column is set according to the locus file
       This was incorrect, we do need to handle ignore columns */

    /* t1 = time((time_t *) NULL); */
    while((colname_item =
           (col_hdr_type *) pop_first_list_entry(userdef_colnames))
          != NULL) {
        if (colname_item->value_type != IGNORE) {

            //locus_found=0;
            parse_hdr(colname_item->ColName, main_part, ext1, ext2, pedigree_valid_header_item_extensions);
            locus_found = search_marker(main_part, &i);
            colname_item->locus_number = i;
            if (!locus_found) {
                errorvf("Column name %s (column %d) is neither a standard keyword nor defined in names file.\n",
                        main_part, colname_item->input_col+1);
                name_not_found++;
            } else {
                colname_item->output_col = output_col;
                output_col++; phenotype_col++;
                loctype_to_descriptor(LTop, i, loctype);
                if (loctype[0] != ext1[0] && loctype[1] != ext1[0]) {
                    errorvf("Invalid header extension %s (column %d) in %s\n",
                            ext1, colname_item->input_col+1, colname_item->ColName);
                    sprintf(err_msg, "Allowed type %c",  loctype[0]);
                    if (loctype[0] != loctype[1]) {
                        grow(err_msg, ", or %c", loctype[1]);
                    }
                    errorf(err_msg);
                    inconsistent_type++;
                }
            }
        }
        copy_colname(colname_item,
                     &(ped_all_colnames[colname_item->input_col]));
        free(colname_item);
    }
    free(userdef_colnames);
    /* t2=time((time_t *) NULL);
       printf("Time taken = %g\n", difftime(t2, t1)); */
    if (name_not_found > 0) {
        errorf("Found mismatch between pedigree file header and names file: ");
        errorf("     one or more locus-names in header missing from names file.");
    } else if (phenotype_col < LTop->NumPedigreeCols) {
        errorf("Found fewer than required pedigree phenotype columns, ");
        errorvf("    need %d columns, found %d columns.\n", LTop->NumPedigreeCols, phenotype_col);
    }

    if (name_not_found || inconsistent_type) {
        errorf("Found invalid extension after locus names inside pedigree file header.");
        EXIT(FILE_HEADER_ERROR);
    }

    if (phenotype_col < LTop->NumPedigreeCols) {
        errorf("Please correct errors and restart Mega2.");
        EXIT(DATA_TYPE_ERROR);
    }
    return ped_all_colnames;
}

/* fill in the output columns within the file_desc, so that when we read in
   marker data, this is ordered according to the names file, i.e. the locus top
   order.
*/

static int parse_pedigree_header(FILE *pedfile,
                                 col_hdr_type *reserved_pedcols,
                                 listhandle *userdef_pedcols)
{
    int colnum;
    /* ped, person, father, mother, sex */

    colnum = parse_variable_width_hdr(pedfile,
				      reserved_pedcols,
                                      userdef_pedcols,
				      NUM_PEDCOL_NAMES, 5,
				      pedigree_valid_header_item_extensions);
/*   for (i=0; i < 5; i++) { */
/*     printf("%s %d\n", reserved_pedcols[i].ColName, */
/* 	   reserved_pedcols[i].input_col); */
/*   } */
    return colnum;
}

/* returns the number of non-reserved columns */
static int parse_variable_width_hdr(FILE *file,
                                    col_hdr_type *reserved_colnames,
                                    listhandle *userdef_ped_cols,
                                    const int num_reserved_cols,
                                    const int required_cols,
                                    const char **valid_header_item_extensions)
{
    int i, col_num=0;
    char buffer[READ_CHUNK+1];
    char *token;//, col_name[40];
    const char *col_name;
    col_hdr_type *user_col;
    int reserved, invalid_extension=0;
    int reading=1;
    char *overflow = &buffer[READ_CHUNK];
    size_t overflow_size, read_size, read_length;

    /* Outer loop to get entire line in one  more chunks */

    while(reading) {
        strcpy(buffer, "");
        if (overflow < &buffer[READ_CHUNK]) {
            strcpy(buffer, overflow);
            overflow_size = &buffer[READ_CHUNK-1] - overflow;
        } else {
            overflow_size = 0;
        }
        read_size = READ_CHUNK - overflow_size;

        IgnoreValue(fgets(buffer + overflow_size, (int)read_size, file));
        read_length = strlen(buffer);
        /* check if line overflows the buffer, and there is a string that
           spans two consecutive reads */
        if (read_length == READ_CHUNK - 1) {
            overflow = buffer + read_length;
            while (*--overflow != ' ' && *overflow != '\t' &&  *overflow != '\n') ;
            *overflow++ = 0;
        }
        if (read_length < (READ_CHUNK-1) || buffer[strlen(buffer)-1] == '\n') {
            reading=0;
        }

        token = strtok(buffer, "\t \n");

        /* Inner loop to get tokens from each chunk */
        while (token != NULL) {
//          strcpy(col_name, token);
            col_name = canonicalColName(token);
	    // 7.2 Annotated input file formats
            // A column is read in and ignored if its name contains “X.” as the first two characters.
	    // Note that the “#” character at the beginning of the column name also works, but the “X.”
	    // convention allows files to be compatible with R. 
            if (col_name[0] == '#' || (col_name[0] == 'X' && col_name[1] == '.')) {
                /* create a column header for the user defined column */
                user_col=CALLOC((size_t)1, col_hdr_type);
//              sscanf(col_name, "%s", user_col->ColName);
                user_col->ColName = col_name;
                user_col->value_type = IGNORE;
                user_col->input_col = col_num;
                user_col->output_length = 0;
                append_to_list_tail(userdef_ped_cols,
                                    (list_data)user_col);
            } else {
                /*	printf("token %s\n", col_name); */
                reserved=0;
                for (i=0; i < num_reserved_cols; i++) {
                    if (SAME(reserved_colnames[i].ColName, col_name)) {
                        reserved_colnames[i].input_col = col_num;
                        reserved = 1;
                        break;
                    }
                }
                if (!reserved) {
		    if (is_valid_hdr_ext(col_name, valid_header_item_extensions) <= 0) {
                        error_hdr_ext(col_name, col_num + 1, valid_header_item_extensions);
                        invalid_extension = 1;
		    }
                    /* create a column header for the user defined column */
                    user_col=CALLOC((size_t) 1, col_hdr_type);
//                  sscanf(col_name, "%s", user_col->ColName);
                    user_col->ColName = col_name;
                    user_col->value_type = STRING_AN;
                    user_col->input_col = col_num;
                    user_col->output_length = 0;
                    user_col->locus_number = -1;
                    append_to_list_tail(userdef_ped_cols, (list_data)user_col);
                }
            }
            col_num++;
            token = strtok(NULL, "\t \n");
        }
    }

/*   sprintf(err_msg, "Read in %d column headers\n", */
/* 	  col_num); */
/*   mssgf(err_msg); */
    if (is_missing_required_hdr(reserved_colnames, required_cols) || invalid_extension) {
        errorf("Found errors in the Mega2 headers, please correct and restart Mega2.\n");
        errorf("Missing required headers.\n");
        EXIT(FILE_HEADER_ERROR);
    }

    return(col_num);
}

// CPK: Here we create a hash table that contians all of the allele
// strings, the idea is to malloc only once. It seems that malloc
// takes up about 2/3's of the time in process_binary_genotype.
// Later in process_binary_genotype we will look for the one malloc'ed
// version of the allele and use it in every locus. We must be careful
// not to free these more than once.
//
static void create_allele_strings_hashtable(plink_info_type *plink_info) {
    int i;
    char allele_string[2];
    allele_string[1] = '\0';
    
    // add all of the allele strings to the hash table...
    for (i=0; i < plink_info->allele_count; i++) {
        allele_string[0] = plink_info->alleles[i];
	(void)canonical_allele(allele_string);
    }
    
    // make sure that this special case is in the hash table...
    allele_string[0] = '0';
    (void)canonical_allele(allele_string);
}

static char **canonical_allele_cache = (char **)NULL;

// CPK: The routine takes advantage of the fact that alleles are only one character long
// in the PLINK world. So, we create a vector to cache the canonical allele string addresses
// returned from 'search_allele' so that we only call it once (when the cache is empty).
// Thank you Robert for suggesting this!
static char *search_allele_cache(const char *allele) {
    char *canonical_allele = canonical_allele_cache[((unsigned char)allele[0])];
    if (canonical_allele == (char *)NULL) {
        // canonical allele not in the cache, so get it, and add it...
        canonical_allele = search_allele((char *)allele);
        canonical_allele_cache[((unsigned char)allele[0])] = canonical_allele;
    }
    return canonical_allele;
}

//
// CPK: Here the Genotype information is reconstituted from the allele[1,2} names
// and the .bed file information. This routine should be called for
// SNP_major mode, or for SNP_minor mode.
static void process_binary_genotype(FILE *bed_filep, annotated_ped_rec *entry,
                                    int mrkindex, int *SNP_count, int *SNP_data,
                                    linkage_locus_rec *locus,
                                    const char *allele1, const char *allele2) {
    char *allele1_from_table, *allele2_from_table;
    int ra_val = 0;
    
    //if (*SNP_count % 8 == 0) *SNP_data = fgetc(bed_filep);
    if ((*SNP_count & 0x7) == 0) *SNP_data = fgetc(bed_filep);
    // each nibble (2bits) is effectively read backwards...
    // NOTE: Since we are only reading bytes, we will not get into trouble over endian issues.
    ra_val |= (*SNP_data & 0x01) != 0 ? 0x02 : 0x00;
    ra_val |= (*SNP_data & 0x02) != 0 ? 0x01 : 0x00;
    
    switch (ra_val) {
        case 0x00:
            // 00  Homozygote "1"/"1"
	    allele1_from_table = search_allele_cache(allele1);
            set_2Ralleles(entry->marker, mrkindex, locus, allele1_from_table, allele1_from_table);
            break;
        case 0x01:
            // 01  Heterozygote
	    allele1_from_table = search_allele_cache(allele1);
	    allele2_from_table = search_allele_cache(allele2);
            set_2Ralleles(entry->marker, mrkindex, locus, allele1_from_table, allele2_from_table);
            break;
        case 0x03:
            // 11  Homozygote "2"/"2"
	    allele2_from_table = search_allele_cache(allele2);
            set_2Ralleles(entry->marker, mrkindex, locus, allele2_from_table, allele2_from_table);
            break;
        case 0x02:
            // 10  Missing genotype
	    allele1_from_table = search_allele_cache("0");
            set_2Ralleles(entry->marker, mrkindex, locus, allele1_from_table, allele1_from_table);
            break;
    }
    
    *SNP_data >>= 2;
    *SNP_count += 2;
}

/* NM: 8-21-08: Changed this routine to read the */

static const char *INIT_COLNAME(col_hdr_type *colnames, int i,
                                col_value_type type, const char *name)
{
//  strcpy(colnames[i].ColName, name);
    name = canonicalColName(name);
    colnames[i].ColName = name;
    colnames[i].value_type=type;
    colnames[i].input_col = -1;
    colnames[i].locus_number = -1;
    colnames[i].map_number = -1;
    colnames[i].sex_map_number = -1;
    colnames[i].output_length = 0;

    return name;
}

static int read_annotated_pedrec(FILE *filep,
                                 char *pedfile,
                                 plink_info_type *plink_info,
                                 int phe_cols,
                                 annotated_ped_rec *entry,
                                 linkage_locus_top *LTop,
                                 annotated_file_desc *file_desc,
                                 char **phe_vals,
                                 int *curr_ped_index, int *curr_per_index,
                                 char **ped_names, int *num_errors)
{

    /* the file descriptor has been reorganized to match the column order
       in the input pedigree file, and locus number have been filled in.
    */
    int i, lch, found=-1;
    char token[FILENAME_LENGTH];
    const char *tokenp = 0;
    int has_pedid = 0;
    int has_link_pedid = 0;
    int has_perid = 0;
    int has_link_perid = 0;
    int has_phenotype = 0;
    char c;
    int phe_cnt = 0; //silly compiler

    col_hdr_type *ped_col_names = file_desc->ped_file_columns;
    *num_errors=0;

    entry->Group = 0;

    /* get rid of leading spaces */
    lch=fscanf(filep, "%c", &(c));
    while (isspace((unsigned char)c) && !feof(filep)) {
        lch=fscanf(filep, "%c", &(c));
    }
    if (c == '\n') {
        /* empty line */
        return -2;
    }

    ungetc(c, filep);

    // CPK: Here we go through the actual number of columns in the file.
    // For a .ped file this is usually 6 plus the number of alleles;
    // for a .bim file this is usually just 6...
//  int REC_NUM = entry->rec_num - 1;

    for (i=0; i < file_desc->num_ped_cols; i++) {

        if (ped_col_names[i].locus_number < 0) {
            lch=fcmap(filep, "%s", token);
            tokenp = token;
            tokenp = canonicalColName(tokenp);
        }

        if (ped_col_names[i].value_type == IGNORE) {
            int mrkindex=ped_col_names[i].locus_number;
            lch=fcmap(filep, "%s", token);
            switch(LTop->Locus[mrkindex].Type) {
                case NUMBERED:
                case XLINKED:
                case YLINKED:
                    lch=fcmap(filep, "%s", token);
                    annot_ignore_numbered_data(entry->rec_num, &(LTop->Locus[mrkindex]),
                                               entry->marker, mrkindex);  /* this will set the pair */
                    
                    i++;
                    break;
                    
                case AFFECTION:
                    plink_annot_string_aff_phen(entry->rec_num, &(LTop->Pheno[mrkindex]),
                                                &entry->pheno[mrkindex], "0");
                    break;
                case QUANT:
                    plink_annot_string_quant_phen(entry->rec_num, &(LTop->Pheno[mrkindex]),
                                                  &entry->pheno[mrkindex], "NA");
                    break;
                default:
                    break;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.Pedigree) ) {
#if 1
//          PedName[REC_NUM] = tokenp;
            entry->Pedigree = tokenp;
           
            if (map_get(ped2idx, tokenp, found)) {
                entry->ped_index = found;
                entry->per_index = ++(*curr_per_index);
            } else {
                entry->ped_index = ++(*curr_ped_index);
                *curr_per_index = 1;
                entry->per_index = *curr_per_index;
                ped2idx[tokenp] = entry->ped_index;
            }
//          PED_NUM = entry->ped_index - 1;
//          PER_NUM = entry->per_index - 1;
            if (PLINK.no_fid) {
//              PerName[REC_NUM] = tokenp;
                entry->ID = entry->Pedigree;
            }
#else
            lch=fcmap(filep, "%s", entry->Pedigree);
            found = -1;
            for (pi=0; pi < *curr_ped_index; pi++) {
                if (!strcmp(ped_names[pi], entry->Pedigree)) {
                    found=pi; break;
                }
            }

            if (found == -1) {
                /* new pedigree and first person */
                entry->ped_index = *curr_ped_index + 1;
                *curr_per_index = 1;
                ped_names[*curr_ped_index] = strdup(entry->Pedigree);
                (*curr_ped_index)++;
                /*	printf("%d, %d\n", *curr_ped_index, *curr_per_index); */
            } else {
                /* pedigree already exists, as this person */
                entry->ped_index = found + 1;
                (*curr_per_index)++;
                /* 	printf("%d, %d\n", *curr_ped_index, *curr_per_index); */
            }
            if (PLINK.no_fid) {
                strcpy(entry->ID, entry->Pedigree);
            }
#endif
        } else if (SAME(ped_col_names[i].ColName, TKN.ID)) {
            PLINK.individuals += 1;
            if (PLINK.no_fid == 0)
                entry->ID = tokenp;
#ifdef IDS
            printf("Ped/per: %s/%s; pedi/peri %d/%d; ",
                   entry->Pedigree, entry->ID, *curr_ped_index, *curr_per_index);
#endif
//
            if (phe_cols > 0)
                phe_cnt = phesearch(entry->Pedigree, entry->ID, phe_vals);
//
//        } else if (SAME(ped_col_names[i].ColName, Father)) {
        } else if (SAME(ped_col_names[i].ColName, TKN.Father)) {
            if (PLINK.no_parents)
                entry->Father = TKN.Zero;
            else
//              lch=fcmap(filep, "%s", entry->Father);
                entry->Father = tokenp;
        } else if (SAME(ped_col_names[i].ColName, TKN.Mother)) {
            if (PLINK.no_parents)
                entry->Mother = TKN.Zero;
            else
                entry->Mother = tokenp;
        } else if (SAME(ped_col_names[i].ColName, TKN.Sex)) {
            entry->Sex = *tokenp++;
            while(isspace((unsigned char)entry->Sex)) {
                entry->Sex = *tokenp++;
            }
            if (tolower((unsigned char)entry->Sex) == 'm' || entry->Sex == '1') {
                PLINK.males += 1;
            } else if (tolower((unsigned char)entry->Sex) == 'f' || entry->Sex == '2') {
                PLINK.females += 1;
            } else {
                PLINK.unspecified_sex += 1;
                entry->Sex = '0' ;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.PedID)) {
            entry->PedID = tokenp;
            if (!SAME(entry->PedID, TKN.NA)) {
                has_pedid=1;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.LinkPedID)) {
//          lch=fcmap(filep, "%d", &(entry->LinkPedID));
            entry->LinkPedID = atoi(tokenp);
            if (entry->LinkPedID > 0) {
                has_link_pedid=1;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.PerID)) {
            entry->PerID = tokenp;
            if (!SAME(entry->PedID, TKN.NA)) {
                has_perid=1;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.LinkPerID)) {
            entry->LinkPerID = atoi(tokenp);
            if (entry->LinkPerID > 0) {
                has_link_perid=1;
            }
        } else if (SAME(ped_col_names[i].ColName, TKN.FirstOff)) {
            entry->FirstOff = tokenp;
        } else if (SAME(ped_col_names[i].ColName, TKN.NextMatSib)) {
            entry->NextMatSib = tokenp;
        } else if (SAME(ped_col_names[i].ColName, TKN.NextPatSib)) {
            entry->NextPatSib = tokenp;

        } else if (SAME(ped_col_names[i].ColName, TKN.MZTwin)) {
            entry->MZTwin = atoi(tokenp);
        } else if (SAME(ped_col_names[i].ColName, TKN.DZTwin)) {
            entry->DZTwin = atoi(tokenp);
        } else if (SAME(ped_col_names[i].ColName, TKN.Proband)) {
            entry->Proband = atoi(tokenp);
        } else if (SAME(ped_col_names[i].ColName, TKN.Group)) {
            entry->Group = atoi(tokenp);

        } else for (; i < file_desc->num_ped_cols; i++) {
            // Here we process the Phenotypes and Genotypes (.ped file)...
            int mrkindex=ped_col_names[i].locus_number;
	    // POSSIBLE BUG: Should the 'IGNORE' case here be a continue rather than a break,
	    // because you still need to get the other items on the line?
            if (mrkindex < 0 || i == file_desc->num_ped_cols ||
                ped_col_names[i].value_type == IGNORE) {
                i--;
                break;
            }
            /*      printf("i, %d, mrkindex, %d, ", i, mrkindex); */
            /* set the input column numbers of the loci only once */
            if (entry->rec_num == 1) {
                LTop->Locus[mrkindex].col_num = ped_col_names[i].input_col;
            }
            switch(LTop->Locus[mrkindex].Type) {
            case NUMBERED:
            case XLINKED:
            case YLINKED:
                if (entry->rec_num == 1)
                    LTop->Marker[mrkindex].col_num = ped_col_names[i].input_col;
                // If this is a .fam file we will never get here because the alleles
                // are not contained in this file. We will add the alleles below...
                lch=read_numbered_data(filep, mrkindex, &LTop->Locus[mrkindex],
                                       (void *) entry, Annotated, 0);
                i++;
                break;
            case AFFECTION:
                if (entry->rec_num == 1)
                    LTop->Pheno[mrkindex].col_num = ped_col_names[i].input_col;

                has_phenotype = 1;
                if (mrkindex < phe_cols) {
//                  phesearch(entry->Pedigree, entry->ID, phe_vals);
                    plink_annot_string_aff_phen(entry->rec_num, &(LTop->Pheno[mrkindex]),
                                                &entry->pheno[mrkindex],
                                                mrkindex < phe_cnt ? phe_vals[mrkindex] : "0");
                } else
                    lch=read_aff_phen(filep, mrkindex, &(LTop->Pheno[mrkindex]),
                                      (void *)entry, Annotated);

                if (entry->pheno[mrkindex].Affection.Status == UNDEF) {
                    errorvf("File %s, Line %d : Invalid status at locus %s\n",
                            pedfile, entry->rec_num, LTop->Marker[mrkindex].MarkerName);
                    (*num_errors)++;
                }

                if (LTop->Pheno[mrkindex].Props.Affection.ClassCnt != 1) {
                    if (entry->pheno[mrkindex].Affection.Class == UNDEF) {
                        errorvf("Ped File \"%s\", Line %d : Invalid liability class at locus %s\n\n",
                                pedfile, entry->rec_num, LTop->Marker[mrkindex].MarkerName);
                        (*num_errors)++;
                    }
                    i++;
                }
                break;
            case QUANT:
                if (entry->rec_num == 1)
                    LTop->Pheno[mrkindex].col_num = ped_col_names[i].input_col;

                if (mrkindex < phe_cols) {
//                    phesearch(entry->Pedigree, entry->ID, phe_vals);
                    plink_annot_string_quant_phen(entry->rec_num, &(LTop->Pheno[mrkindex]),
                                                  &entry->pheno[mrkindex],
                                                  mrkindex < phe_cnt ? phe_vals[mrkindex] : "NA");
                } else
                    lch=read_quant_phen(filep, mrkindex, &(LTop->Pheno[mrkindex]), (void *)entry, Annotated);

		// Here we are not checking for it being undefined, just "invalid".
		// There should be a better way...
                if (entry->pheno[mrkindex].Quant == QUNDEF) {
                    errorvf("Ped File \"%s\", Line %d : Invalid quantitative phenotype at locus %s\n",
                            pedfile, entry->rec_num, LTop->Marker[mrkindex].MarkerName);
                    (*num_errors)++;
                }
                break;
            default:
                break;
            }
            // Check for premature EOR during the MARKER, AFFECTION, QUANT case loop....
            // We will get here for the MARKER data if PLINK.plink == PED_format.
            // Remember that the .fam file in the binary_PED_format has no marker data.
            // Either way (binary or not) we will get here for the AFFECTION.
            if ((lch == '\n' || lch == '\r') && i < (file_desc->num_ped_cols-1)) {
                *curr_ped_index = -1;
            }
        } /* else for */
        // check for premature EOR during the outer loop...
        if ((lch == '\n' || lch == '\r') && i < (file_desc->num_ped_cols-1)) {
            // ??? should we 'break;' here and do this test when out of the loop?
            *curr_ped_index = -1;
        }
    } /* for */

    if (SAME(entry->Mother, TKN.Zero) && SAME(entry->Father, TKN.Zero) == 0) {
        PLINK.founders += 1;
    } else {
        PLINK.non_founders += 1;
    }
    if (has_phenotype != 0) {
        PLINK.phenotyped_individuals += 1;
    }
    
    // CPK: Here we add the alleles if this is a .fam file (we found a .bed file), and that .bed file 
    // is in Person_major mode (where each line contains all of the SNPs for an individual),
    // then we process the SNPs as we find them in the .bed file...
    if   (PLINK.plink == binary_PED_format &&
          plink_info != NULL &&
          plink_info->bed_filep != (FILE *)NULL &&
          plink_info->SNP_major != 0x01) {

        Tod tod_pl_nmj("plink read indiv major bed");
        int allele_i = 0;               // an index for the alleles array
        char allele1[2], allele2[2];    // for converting allele characters to string arrays
        allele1[1] = '\0';
        allele2[1] = '\0';
        // cycle through the locusts searching for allele markers....
        for (i=0; i < file_desc->num_ped_cols + plink_info->allele_count; i++) {
            int mrkindex=ped_col_names[i].locus_number;
            if (mrkindex > 0 && LTop->Locus[mrkindex].Class == MARKER) {
                linkage_locus_type type = LTop->Locus[mrkindex].Type;
                if (type == NUMBERED || type == XLINKED || type == YLINKED){
                    // we found an allele marker...
                    int SNP_count = 0;  // always start on a new byte
                    int SNP_data;       // byte from the binary file
                    allele1[0] = plink_info->alleles[allele_i++];
                    allele2[0] = plink_info->alleles[allele_i++];
                    process_binary_genotype(plink_info->bed_filep, entry, mrkindex, &SNP_count, &SNP_data, 
                                            &LTop->Locus[mrkindex], allele1, allele2);
                    i++;
                }
            }
            
        }
        tod_pl_nmj();
//        @@ plink process indiv major: 0.157030 original design
//        @@ plink process indiv major: 0.036451 at 2 bits at a time
//        @@ plink process indiv major: 0.014141 at byte at a time

    }

/*   printf("%s, %c\n", entry->ID, entry->Sex); */
    if (lch != '\n') {
        /* skip to the end of the file */
        fcmap(filep, "%=\n", lch);
    }

    if (!has_link_perid) {
        entry->LinkPerID = *curr_per_index;
    }
    if (!has_link_pedid) {
        entry->LinkPedID = *curr_ped_index;
    }

    if (!has_pedid) {
/*      if (has_link_pedid) {
            sprintf(entry->PedID, "%d", entry->LinkPedID);
        } else {
*/            entry->PedID = entry->Pedigree;
/*      }*/
    }

    if (!has_perid) {
/*      if (has_link_perid) {
            sprintf(entry->PerID, "%d", entry->LinkPerID);
        } else {
*/            entry->PerID = entry->ID;
/*      }*/
    }

//  entry->per_index = *curr_per_index;
#ifdef IDS
    printf("Ped/per: %s/%s; Lped/Lperi %d/%d; %d; %s\n",
           entry->PedID, entry->PerID, entry->LinkPedID, entry->LinkPerID, *curr_per_index,
           entry->ID);
#endif
    return *curr_ped_index;
}

static int sort_byped(const void *p1, const void *p2)
{
    const annotated_ped_rec *person1, *person2;

    person1 = (const annotated_ped_rec *) p1;
    person2 = (const annotated_ped_rec *) p2;

    if (person1->ped_index < person2->ped_index) {
        return -1;
    } else if (person1->ped_index == person2->ped_index) {
        if (person1->LinkPerID < person2->LinkPerID) {
            return -1;
        } else if (person1->LinkPerID == person2->LinkPerID) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 1;
    }
}

/* Have to reassign the LinkPerIDs, since these are record numbers now.
   If pedigree records are not continuous, this will cause problems.
*/

static int sort_and_check(int num_recs, annotated_ped_rec *persons,
                          int num_groups, int *groups,
                          int has_link_perid, int has_link_pedid)
{
    int gfound, g, p;
    int (*sort_persons)(const void *, const void *) = sort_byped;
    int grp_err =0,  par_err=0, sex_err=0;
    int person_num=1, curr_ped = -1;

    qsort((void *)persons, (size_t) num_recs, sizeof(annotated_ped_rec),
          sort_persons);

    for (p=0; p < num_recs; p++) {
        if (num_groups > 0) {
            /* store the group */
            gfound=-1;
            for (g=0; g < num_groups; g++) {
                if (persons[p].Group == groups[g]) {
                    gfound = g;
                    break;
                }
            }
            if (gfound == -1) {
                errorvf("Invalid group number %d in record %d.\n",
                        persons[p].Group, persons[p].rec_num);
                grp_err++;
            }
        }
        /* check gender field */
        if (!(tolower((unsigned char)persons[p].Sex) == 'm' ||
              tolower((unsigned char)persons[p].Sex) == 'f' ||
              persons[p].Sex == '1' ||
              persons[p].Sex == '2')) {
            if (persons[p].Sex != '0') {
                errorvf("Invalid sex %c in record %d.\n",
                        persons[p].Sex, persons[p].rec_num);
            }
            sex_err++;
        }
        /* check both parents are defined or NA */
        if (SAME(persons[p].Father, TKN.NA) &&
            !SAME(persons[p].Mother, TKN.NA)) {
            errorvf("Father defined & mother undefined in record %d.\n",
                    persons[p].rec_num);
            par_err++;
        } else if (!SAME(persons[p].Father, TKN.NA) &&
                   SAME(persons[p].Mother, TKN.NA)) {
            errorvf("Mother defined & father undefined in record %d.\n",
                    persons[p].rec_num);
            par_err++;
        }
        /* Reassign LinkPerID if this was created by Mega2 */
        if (persons[p].ped_index != curr_ped) {
            person_num = 1;
            curr_ped = persons[p].ped_index;
        }
        if (!has_link_perid) {
            persons[p].LinkPerID = person_num;
            person_num++;
        }
        if (!has_link_pedid) {
            persons[p].LinkPedID = persons[p].ped_index;
        }
    }

    if (grp_err) {
        errorvf("Found %d records with invalid group number.\n", grp_err);
    }

    if (sex_err) {
        errorvf("Found %d records with invalid sex value.\n", sex_err);
    }

    if (par_err) {
        errorvf("Found %d records with errors in parent fields.\n", par_err);
    }

    return(sex_err+par_err+grp_err);
}

/* this is equivalent to the make_ped_structure() function
   in makeped1.c.
   We have to convert the parent names to ids
*/
static int UNIQUEid = 0; /* unless LinkPerID/LinkPedID is set */

void unique_id(char *reslt, const char *ped, const char *per) {
    if (UNIQUEid || PLINK.no_fid)
        strcpy(reslt, per);
    else
        sprintf(reslt, "%s_%s", ped, per);
}

char *unique_id_per(char *reslt, char *ped) {
    if (UNIQUEid || PLINK.no_fid)
        return reslt;
    char *mid = strrchr(reslt, '_');
    if (mid)
        return mid + 1;
    else
        return reslt;
}

static marriage_graph_type *copy_annotated_to_premake(linkage_locus_top *LTop,
                                                      int num_persons,
                                                      annotated_ped_rec *persons,
                                                      int num_peds,
                                                      int unique)
{
    int curr_ped;
    int ped, per, k;
    size_t p, lower_p, upper_p;
    int *fathers, *mothers;

    marriage_graph_type *pped;
    fathers = CALLOC((size_t)num_persons, int);
    mothers = CALLOC((size_t)num_persons, int);

    pped = CALLOC((size_t)num_peds, marriage_graph_type);
    for (ped = 0; ped < num_peds; ped++) {
        pped[ped].num_marriages = 0;
        pped[ped].marriages = NULL;
        pped[ped].max_id = 0;
        pped[ped].num_persons = 0;
    }

    for (per = 0; per < num_persons; per++) {
        ped = persons[per].ped_index - 1;
        if (pped[ped].num_persons == 0) {
            strcpy(pped[ped].Name, persons[per].Pedigree);
//pp        strcpy(pped[ped].Name, persons[per].PedID);
            strcpy(pped[ped].PedPre, persons[per].PedID);
        }
        pped[ped].num_persons++;
    }

    for (ped = 0; ped < num_peds; ped++) {
        pped[ped].persons =
	  CALLOC((size_t) pped[ped].num_persons, person_node_type);
    }

    /* lower_p, and upper_p are indices into the persons array.
       These are the limits of the search for parent IDs
    */

    curr_ped = persons[0].ped_index - 1; /* index starts from 1 */
    lower_p = 0;
    upper_p = pped[curr_ped].num_persons;

    for (per = 0; per < num_persons; per++) {
        if (per == 0) {
            pped[curr_ped].ped = persons[per].LinkPedID;
        }
        if ((persons[per].ped_index - 1) != curr_ped) {
            /* increment lower_p */
            lower_p += pped[curr_ped].num_persons;
            curr_ped = persons[per].ped_index - 1; /* index starts from 1 */
            /* set upper_p */
            upper_p += pped[curr_ped].num_persons;
            pped[curr_ped].ped = persons[per].LinkPedID;
        }
        /*    printf("%d - %d\n", lower_p, upper_p); */
        if (SAME(persons[per].Father, TKN.NA)) {
            fathers[per] = 0;
            mothers[per] = 0;
        } else {
            for (p = lower_p; p < upper_p; p++) {
                if (SAME(persons[p].ID, persons[per].Father)) {
                    fathers[per] = persons[p].LinkPerID;
                }
                if (SAME(persons[p].ID, persons[per].Mother)) {
                    mothers[per] = persons[p].LinkPerID;
                }
            }
        }
    }

    curr_ped = persons[0].ped_index - 1;
    p = 0;
    per = 0;
    while (per < num_persons) {
        while (curr_ped == persons[per].ped_index - 1) {
            pped[curr_ped].persons[p].indiv = persons[per].LinkPerID;
            pped[curr_ped].persons[p].father = fathers[per];
            pped[curr_ped].persons[p].mother = mothers[per];

            unique_id(pped[curr_ped].persons[p].uniqueid, persons[per].Pedigree, persons[per].ID);
            strcpy(pped[curr_ped].persons[p].origid, persons[per].PerID);
            strcpy(pped[curr_ped].persons[p].perpre, persons[per].PerID);
//          strcpy(pped[curr_ped].persons[p].pedpre, persons[per].Pedigree); // do this later
            strcpy(pped[curr_ped].persons[p].famname, persons[per].Pedigree);
//          strcpy(pped[curr_ped].persons[p].pedpre, persons[per].PedID);  // do this later
//pp        strcpy(pped[curr_ped].persons[p].famname, persons[per].PedID); // cause lpedtop does
	    if ((persons[per].Sex == '1') || (tolower((unsigned char)persons[per].Sex) == 'm')) {
                pped[curr_ped].persons[p].gender = 1;
            } else if ((persons[per].Sex == '2') || (tolower((unsigned char)persons[per].Sex) == 'f')) {
                pped[curr_ped].persons[p].gender = 2;
            } else {
                pped[curr_ped].persons[p].gender = 0;
            }

            pped[curr_ped].persons[p].pheno  = persons[per].pheno;
            pped[curr_ped].persons[p].marker = persons[per].marker;
            persons[per].pheno  = NULL;
            persons[per].marker = NULL;
            pped[curr_ped].persons[p].MZTwin  = persons[per].MZTwin;
            pped[curr_ped].persons[p].DZTwin  = persons[per].DZTwin;
            pped[curr_ped].persons[p].Group   = persons[per].Group;
            pped[curr_ped].persons[p].genocnt = persons[per].genocnt;

            /* initialize some fields for makeped */
            pped[curr_ped].persons[p].node_id = (int) p;
            pped[curr_ped].persons[p].proband = ((p == 0) ? 1 : 0);
            pped[curr_ped].persons[p].from_marriage_node_id = -1;
            pped[curr_ped].persons[p].num_to_marriages = 0;
            for (k = 0; k < MAXMARRIAGES; k++)
                pped[curr_ped].persons[p].to_marriage_node_id[k] = -1;
            pped[curr_ped].persons[p].removed = 0;
            pped[curr_ped].persons[p].loop_breaker_id = 0;
            if (pped[curr_ped].max_id < pped[curr_ped].persons[p].indiv) {
                pped[curr_ped].max_id = pped[curr_ped].persons[p].indiv;
            }

            per++;
            p++;
            if (per >= num_persons)
                break;
        }
        if (per >= num_persons) {
            break;
        } else {
            curr_ped = persons[per].ped_index - 1;
            p = 0;
        }
    }

    free(fathers);
    free(mothers);

    return pped;
}


static linkage_ped_tree *copy_annotated_to_lpedtop(linkage_locus_top *LTop,
                                                   int num_persons,
                                                   annotated_ped_rec *persons,
                                                   int num_peds, int unique)
{
    int i, curr_ped;
    int ped, per;
    size_t p, lower_p, upper_p;
    int pa, mom, first_off, next_p_sib, next_m_sib;
    int **links = CALLOC((size_t)num_persons, int *);
    linkage_loop_rec *loop;

    linkage_ped_tree *lped = CALLOC((size_t)num_peds, linkage_ped_tree);
    linkage_ped_rec *entry;


    pa=0;  mom=1;
    first_off=2;   next_p_sib=3;   next_m_sib = 4;

    for (ped=0; ped < num_peds; ped++) {
        clear_lpedtree(&(lped[ped]));
    }

    //ped=0;
    for (per = 0; per < num_persons; per++) {
        links[per] = CALLOC((size_t)5, int);
        ped = persons[per].ped_index-1;
        if (lped[ped].EntryCnt == 0) {
            lped[ped].Num = persons[per].LinkPedID;
            strcpy(lped[ped].Name, persons[per].Pedigree);
            strcpy(lped[ped].PedPre, persons[per].Pedigree);
//pp        strcpy(lped[ped].Name, persons[per].PedID);
//pp        strcpy(lped[ped].PedPre, persons[per].PedID);
        }
        lped[ped].EntryCnt++;
    }

    for (ped=0; ped < num_peds; ped++) {
        lped[ped].Entry=CALLOC((size_t) lped[ped].EntryCnt,
                               linkage_ped_rec);
    }

    /* lower_p, and upper_p are indices into the persons array.
       These are the limits of the search for parent IDs
    */

    curr_ped = persons[0].ped_index - 1; /* index starts from 1 */
    lower_p=0;
    upper_p = lped[curr_ped].EntryCnt;

    for (per=0; per < num_persons; per++) {
        if ((persons[per].ped_index - 1) != curr_ped) {
            /* increment lower_p */
            lower_p += lped[curr_ped].EntryCnt;
            curr_ped = persons[per].ped_index - 1; /* index starts from 1 */
            /* set upper_p */
            upper_p  += lped[curr_ped].EntryCnt;
        }

        for (p=lower_p; p < upper_p; p++) {
            if (SAME(persons[per].Father, TKN.NA)) {
                links[per][pa] = 0;
                links[per][mom] = 0;
            } else {
                if (SAME(persons[p].ID, persons[per].Father)) {
                    links[per][pa]=persons[p].LinkPerID;
                }
                if (SAME(persons[p].ID, persons[per].Mother)) {
                    links[per][mom]=persons[p].LinkPerID;
                }
            }
            if (SAME(persons[per].FirstOff, TKN.NA)) {
                links[per][first_off]=0;
            } else {
                if (SAME(persons[p].ID, persons[per].FirstOff)) {
                    links[per][first_off]=persons[p].LinkPerID;
                }
            }

            if (SAME(persons[per].NextPatSib, TKN.NA)) {
                links[per][next_p_sib]=0;
            } else {
                if (SAME(persons[p].ID, persons[per].NextPatSib)) {
                    links[per][next_p_sib]=persons[p].LinkPerID;
                }
            }

            if (SAME(persons[per].NextMatSib, TKN.NA)) {
                links[per][next_m_sib]=0;
            } else {
                if (SAME(persons[p].ID, persons[per].NextMatSib)) {
                    links[per][next_m_sib]=persons[p].LinkPerID;
                }
            }
        }
    }

    curr_ped = persons[0].ped_index - 1;
    p=0; per = 0;
    while (per < num_persons) {
        while (curr_ped == persons[per].ped_index - 1) {
            entry=&(lped[curr_ped].Entry[p]);
            clear_lpedrec(entry);

            unique_id(entry->UniqueID, persons[per].Pedigree, persons[per].ID);
            strcpy(entry->OrigID, persons[per].PerID);
            strcpy(entry->PerPre, persons[per].PerID);
            strcpy(entry->FamName, persons[per].PedID);
            entry->ID = persons[per].LinkPerID;
            entry->Father = links[per][pa];
            entry->Mother = links[per][mom];
            entry->First_Offspring = links[per][first_off];
            entry->Next_PA_Sib = links[per][next_p_sib];
            entry->Next_MA_Sib = links[per][next_m_sib];
            entry->OrigProband = persons[per].Proband;
            /* Have to store the loop breakers */
            if (persons[per].Proband >= 1) {
                if (persons[per].Proband == 1) {
                    /* could be only the proband or
                       a loop-breaker proband */
                    lped[curr_ped].Proband = entry->ID;
                }
                loop = get_llooprec(&(lped[curr_ped].Loops),
                                    persons[per].Proband - 1);
                /* check that the first n entries are defined */
                for (i=0; i<loop->numEntry; i++) {
                    if (loop->Entry[i]==UNDEF) {
                        errorvf("INTERNAL error setting loop ids.\n");
                        EXIT(SYSTEM_ERROR);
                    }
                }
                if (loop->numEntry >= MAXLOOPMEMBERS) {
                    errorvf("Exceeded loop limit for an individual\n");
                    EXIT(EXCEEDED_FIELD_MAX);
                }

                loop->Entry[loop->numEntry] = persons[per].LinkPerID;
/* 	  printf("%d proband %d ID, %d Father, %d Mother\n",  */
/* 		 persons[per].Proband, entry->ID, entry->Mother, entry->Father); */
                loop->numEntry++;
            }
            if ((persons[per].Sex=='1') ||
                (tolower((int)persons[per].Sex) == 'm')) {
                entry->Sex = 1;
            } else {
                entry->Sex = 2;
            }

            entry->Pheno  = persons[per].pheno;
            entry->Marker = persons[per].marker;
            persons[per].pheno  = NULL;
            persons[per].marker = NULL;
            entry->MZTwin=persons[per].MZTwin;
            entry->DZTwin=persons[per].DZTwin;
            entry->Group=persons[per].Group;
            entry->genocnt = persons[per].genocnt;
            per++; p++;
            if (per >= num_persons) {
                break;
            }
        }
        if (per >= num_persons) {
            break;
        } else {
            curr_ped = persons[per].ped_index-1;
            p=0;
        }
    }

    for (per =0; per < num_persons; per++) {
        free(links[per]);
    }
    free(links);

    return lped;
}

static linkage_ped_top *read_annotated_ped_file(char *pedfile,
                                                linkage_locus_top *LTop,
                                                annotated_file_desc *file_desc,
                                                int num_groups,
                                                int *groups)
{
    int num_userdef_cols, num_reserved_cols, c;
    int has_extra_ids = 0;
    int num_ped_records;
    listhandle *userdef_colnames = new_list();
    col_hdr_type *ped_all_colnames;

    FILE *filep = pedfile ? fopen(pedfile, "r") : NULL;
    if (filep == NULL) {
        errorvf("could not open %s for reading!\n", pedfile);
        EXIT(FILE_READ_ERROR);
    }

    init_reserved_pedcol_names();
#ifndef HIDEFILE
    mssgvf("Reading Mega2 format pedigree file: %s\n", pedfile);
#endif
    num_userdef_cols = parse_pedigree_header(filep,
                                             ReservedColnames,
                                             userdef_colnames);

    ped_all_colnames =
        pedcol_check_reverse_index(ReservedColnames, num_userdef_cols, userdef_colnames,
                                   LTop, &num_reserved_cols, &has_extra_ids);

    file_desc->ped_file_columns = ped_all_colnames;

    file_desc->num_ped_cols = num_userdef_cols;
    c=fgetc(filep);
/*  skip_long_line(filep, c);  wrong as parse_pedigree_header reads entire line. */
    num_ped_records = linecount(filep);
/*   printf("%d\n", num_ped_records); */
    c=fgetc(filep);
    skip_long_line(filep, c);   /* skip header */
    /* set pedigree file type */

/* BUG */
    UNIQUEid = has_extra_ids;
    return read_common_ped_file(filep, pedfile, NULL, LTop, file_desc,
                                0, num_groups, groups,
                                num_ped_records, has_extra_ids);
}

#ifndef UCHAR_MAX
#define UCHAR_MAX            256
#endif

static linkage_ped_top *read_common_ped_file(FILE *filep, char *pedfile,
                                             plink_info_type *plink_info,
                                             linkage_locus_top *LTop,
                                             annotated_file_desc *file_desc,
                                             int phecols,
                                             int num_groups,
                                             int *groups,
                                             int num_ped_records,
                                             int has_extra_ids)
{
    int num_peds=0, num_ped_per=0, p, c;
    int num_pheno_errs, retval, num_err=0;
    char **ped_names, **phe_vals = 0;
    linkage_ped_top *Top;
    annotated_ped_rec *persons;
    int untyped = 0, totaltyped = 0;
    int check_ungenotyped = 0;
    int xcf = Input_Format == in_format_binary_VCF ||
	      Input_Format == in_format_compressed_VCF ||
              Input_Format == in_format_VCF ;

    // allele_count should be == 0 when using a .ped file...

    if (phecols > 0)
        phe_vals = CALLOC((size_t)phecols, char *); /* buffer */

    ped_names = CALLOC((size_t)num_ped_records, char *);
#if 0
    PedName.resize(num_ped_records);
    PedIdName.resize(num_ped_records);
    PedLinkIdName.resize(num_ped_records);
    PerName.resize(num_ped_records);
    PerIdName.resize(num_ped_records);
    PerLinkIdName.resize(num_ped_records);
    PerUniqName.resize(num_ped_records);
#endif
/*
 *  Note: annotated_ped_rec is 2352 bytes because of tall the MAX_NAMELEN allocated strings.
 *  It is freed at the end of this function.
 */
    if ((persons = CALLOC((size_t)num_ped_records, annotated_ped_rec)) == NULL) {
        errorf("Could not allocate enough memory, exiting.");
        EXIT(MEMORY_ALLOC_ERROR);
    }

    // If we are processing a PLINK file allocate the allele chache now...
    if (PLINK.plink) {
        // Allocate an array of character pointers large enough to hold every character that we find.
        // Later we will index this array by the ascii value of the character to find the pointer to
        // the canonical character string.
        // NOTE: This assumes 8-bit ASCII characters.
      if ((canonical_allele_cache = CALLOC((size_t)UCHAR_MAX, char *)) == (char **)NULL) {
        errorf("Could not allocate enough memory, exiting.");
        EXIT(MEMORY_ALLOC_ERROR);
      }
    }

//  SECTION_ERR_INIT(untyped_msg);
    p=0;
    Tod tod_cp_all_eof("read ped file");
    Tod tod_cp_eof(20);
    while (!feof(filep)) {
        c=fgetc(filep);
        if (c == -1) {
            break;
        } else if (c == '#') {
            skip_long_line(filep, c);
            continue;
        } else if (c == '\n' || !isascii(c)) {
            continue;
        } else {
            ungetc(c, filep);
        }
        persons[p].pheno  = (LTop->PhenoCnt > 0) ? CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;
        // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
        persons[p].marker = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;

        persons[p].rec_num = p+1;

        tod_cp_eof.reset();
        retval = read_annotated_pedrec(filep, pedfile, plink_info, phecols,
                                       &(persons[p]), LTop,
                                       file_desc, phe_vals,
                                       &num_peds, &num_ped_per,
                                       ped_names, &num_pheno_errs);

        tod_cp_eof("ped line");
        totaltyped++;

        if ( ((PLINK.plink != binary_PED_format) && (! xcf)) ||
             (PLINK.plink == binary_PED_format && plink_info->SNP_major != 0x01) ) {
//              NB. VCF files look like SNP_major == 1
            persons[p].genocnt = crunch_Rnotype(&persons[p].marker, LTop);
            if (persons[p].genocnt == 0) {
/*
                SECTION_ERR(untyped_msg);
                warnvf("Untyped person %4d linenum %d: person %s, fa %s, ma %s\n",
                       untyped, p+1, persons[p].ID, persons[p].Father, persons[p].Mother);
*/
                untyped++;
            }
        }

        if (retval == -2) {
            /* empty line */
            continue;
        } else if (retval == -1) {
            errorvf("Ped File \"%s\": Line %d does not have the required number of items.\n",
                    pedfile, p+1);
            num_err++;
            p++;
        } else {
            num_peds = retval;
            p++;
        }
        num_err += num_pheno_errs;
    }
//  SECTION_ERR_FINI(untyped_msg);
    tod_cp_all_eof();
    num_ped_records = p;
#ifdef SHOWSTATUS
    int fudge;
    if (MARKER_SCHEME == MARKER_SCHEME_BITS)
        fudge = LTop->LocusCnt * (sizeof (Alleles_int) + sizeof (Alleles_str));
    else
        fudge = 0;
    msgvf("ALLOC SPACE: ALL pedrec: %d MB (%d x %d)\n",
          (p * marker_size(LTop->MarkerCnt) + fudge) / 1024 / 1024,
          p,  marker_size(LTop->MarkerCnt));
#endif
/*   printf("%d\n", num_ped_records); */
    
    Tod tod_cp_smj("plink read snp major bed");
    if   (PLINK.plink == binary_PED_format &&
          plink_info != NULL &&
          plink_info->bed_filep != (FILE *)NULL &&
          plink_info->SNP_major == 0x01) {
        // CPK: We only process the SNPs here if in SNP_major mode
        // where each line is all of the individuals for a SNP
        col_hdr_type *ped_col_names = file_desc->ped_file_columns;
        int allele_i = 0;               // an index for the alleles array
        char allele1[2], allele2[2];    // for converting allele characters to string arrays
        int pp, i;
        allele1[1] = '\0';
        allele2[1] = '\0';
        // cycle through the locusts searching for allele markers....
        // cycle through the locusts searching for allele markers....
        Tod tod_plink(20);
        for (i=0; i < file_desc->num_ped_cols + plink_info->allele_count; i++) {
            int mrkindex=ped_col_names[i].locus_number;
            if (mrkindex > 0 && LTop->Locus[mrkindex].Class == MARKER) {
                linkage_locus_type type = LTop->Locus[mrkindex].Type;
                if (type == NUMBERED || type == XLINKED || type == YLINKED){
                    // we found an allele marker...
                    int SNP_count = 0;  // always start on a new byte
                    int SNP_data;       // byte from the binary file
                    tod_plink.reset();
                    allele1[0] = plink_info->alleles[allele_i++]; // allele characters from the string
                    allele2[0] = plink_info->alleles[allele_i++];
                    for (pp=0; pp < num_ped_records; pp++) {
                        // looping through the individuals...
                        process_binary_genotype(plink_info->bed_filep, &(persons[pp]), mrkindex, &SNP_count, &SNP_data, 
                                                &LTop->Locus[mrkindex], allele1, allele2);
                    }
                    tod_plink("cycle for individs");
                    i++;
                }
            }
        }
        check_ungenotyped = 1;
    }
    tod_cp_smj();

    if (xcf) {
//      if (! getenv("_")) { asm("int $3"); }
        std::string alternative_key = std::string(Mega2BatchItems[/* 57 */ VCF_Marker_Alternative_INFO_Key].value.name);

        Tod vcfpe(" VCFtools_process_entries(persons, num_ped_records, ...)");
        VCFtools_process_entries(persons, (unsigned int)num_ped_records, LTop, alternative_key,
                                 "chr", VecAlleles);
        VCFtools_close();
        check_ungenotyped = 1;
        vcfpe();
    }

    // CPK: At this point we have finished processing the .bed file...
    if   (PLINK.plink == binary_PED_format &&
          plink_info != NULL &&
          plink_info->bed_filep != (FILE *)NULL) {
        fclose(plink_info->bed_filep);
        plink_info->bed_filep = (FILE *)NULL;
    }

    Top = mk_ped_top(persons, num_ped_records, LTop, num_peds, untyped,
                     totaltyped, groups, num_groups, has_extra_ids,
                     num_err, check_ungenotyped);

    if (phecols > 0) free(phe_vals);
    if (ped_names) {
        for (c = 0; c < num_peds; c++)
            free(ped_names[c]);
        free(ped_names);
    }
    for (c = 0; c < num_ped_records; c++) {
        /* annotated_ped_rec has linkage_pedrec_data (for ped locus data) */
        if (persons[c].pheno != NULL) free(persons[c].pheno);
        if (persons[c].marker != NULL) free(persons[c].marker);
    }
    free(persons);

    return Top;
}

linkage_ped_top *mk_ped_top(annotated_ped_rec *persons, int num_ped_records, 
                            linkage_locus_top *LTop, int num_peds,
                            int untyped, int totaltyped,
                            int *groups, int num_groups, int has_extra_ids,
                            int num_err, int check_ungenotyped)
{
    marriage_graph_type *ppeds;
    linkage_ped_tree *lpeds;
    linkage_ped_top *Top;

    Tod tod_cp_cru("check all markers untyped");
    SECTION_ERR_INIT(untyped_msg);
    if (check_ungenotyped) {
        for (int pp=0; pp < num_ped_records; pp++) {
            // looping through the individuals...
            persons[pp].genocnt = crunch_Rnotype(&persons[pp].marker, LTop);
            if (persons[pp].genocnt == 0) {
/*
                SECTION_ERR(untyped_msg);
                warnvf("Untyped person %4d snp major record %d: person %s, fa %s, ma %s\n",
                       untyped, pp, persons[pp].ID, persons[pp].Father, persons[pp].Mother);
*/
                untyped++;
            }
        }
    }
    tod_cp_cru();

    SECTION_ERR_FORCE(untyped_msg);
    if (untyped > 0)
        warnvf("Individuals untyped: %d out of %d\n", untyped, totaltyped);
    SECTION_ERR_FINI(untyped_msg);

    if (num_err == 0) {
        Tod tod_cp_hmm("sort, copy ann to premake/lpedtop");
        sort_and_check(num_ped_records, persons, num_groups, groups,
                       has_extra_ids & 1, has_extra_ids & 2);

        Top = new_lpedtop();
        Top->pedfile_type = pedfile_type;
        if (pedfile_type == PREMAKEPED_PFT) {
            LTop->PedRecDataType = Raw_premake;
            ppeds = copy_annotated_to_premake(LTop, num_ped_records, persons, num_peds,
                                              has_extra_ids & (1|2));
            Top->PTop = ppeds;
        } else {
            LTop->PedRecDataType = Raw_postmake;
            lpeds = copy_annotated_to_lpedtop(LTop, num_ped_records, persons, num_peds,
                                              has_extra_ids & (1|2));
            Top->Ped=lpeds;
        }

/*      Top->UniqueIds = (has_extra_ids & (1|2)) > 0;     1 == LinkPerId; 2 == LinkPedId */
        Top->UniqueIds = 1;
/*      UNIQUEid = (has_extra_ids & (1|2)) > 0; */
        Top->OrigIds =   (has_extra_ids & (4|8)) > 0;  /* 4 ==     PerId; 8 ==     PedId */
        Top->PedCnt = num_peds;
        tod_cp_hmm();

    } else {
        Top=NULL;
    }

    return Top;
}

#define OMITFL_NAME         "ylinked_females.omit"

void clear_YLINKED_females(linkage_ped_top *Top, int raw_allele, int hdr)
{
    linkage_locus_top *LTop = Top->LocusTop;
    int mrkindex;

    int pedi, pedE, peri, perE;
    void *marker;
    int sex;
    int ylinked_females = 0;
    char *omitfl_name = (char *)NULL;
    FILE *omitfp = (FILE *)NULL;

    // CPK: Now set the Y-linked genotypes of females to unknown...
    //
    // We do this here because if we are processing a binary_PED_format file in SNP Major mode
    // we will not have the alleles till now.
    //
    // A potential "issue" is the we are doing this here before the .omit file has been read,
    // so the user may have these entries in their omit file, or they may not. They may have some
    // of them, but not all of them. Probably the best thing to do is to move this check to after
    // the point where the omit file is read, and then check for any YLINKED Female locus that have
    // not been set to unknown. Currently we just list them all, warn the user that these may not all
    // be in their current (if any) .omit file, and leave it at that.

    for (mrkindex=0; mrkindex < LTop->LocusCnt; mrkindex++) {
        if (LTop->Locus[mrkindex].Type == YLINKED) {
            pedE = Top->PedCnt;
            for (pedi = 0; pedi < pedE; pedi++) {
/*              ped = (Top->pedfile_type == PREMAKEPED_PFT) ? &Top->PTop[pedi] : &Top->Ped[pedi]; */
                perE = (Top->pedfile_type == PREMAKEPED_PFT) ? Top->PTop[pedi].num_persons : Top->Ped[pedi].EntryCnt;
                for (peri = 0; peri < perE; peri++) {
                    sex   = (Top->pedfile_type == PREMAKEPED_PFT) ? Top->PTop[pedi].persons[peri].gender : // PRALLELE
                                                                    Top->Ped[pedi].Entry[peri].Sex;
                    marker  = (Top->pedfile_type == PREMAKEPED_PFT) ? Top->PTop[pedi].persons[peri].marker : // PRALLELE
                                                                    Top->Ped[pedi].Entry[peri].Marker;
                    if (sex == 2) {
                        if (raw_allele) {
                            const char  *all1, *all2;
                            get_2Ralleles(marker, mrkindex, &all1, &all2);
                            if (! (allelecmp( all1, REC_UNKNOWN) ||
                                   allelecmp( all2, REC_UNKNOWN) ) )
                                continue;
                        } else {
                            int all1, all2;
                            get_2alleles(marker, mrkindex, &all1, &all2);
                            if ( all1 == 0 && all2 == 0 )
                                continue;
                        }

                        // We found a Y-linked female!
                        if (omitfp == (FILE *)NULL) {
                            // Open the file if it has not been opened before now...
                            omitfl_name=(char *)CALLOC(23+strlen(sumdir), char);
                            sprintf(omitfl_name, "%s/%s", sumdir, OMITFL_NAME);
                            // open the omit file if this is the first time through...
                            if ((omitfp = fopen(omitfl_name, "w")) == (FILE *)NULL && ylinked_females == 0) {
			        errorvf("Unable to open omit file '%s' for writing Y-linked female record(s).\n",
				        omitfl_name);
                                EXIT(FILE_WRITE_ERROR);
                            } else {
                                if (hdr) {
                                    // Write the header, the reading code expects to see it...
                                    fprintf(omitfp, "Pedigree Individual Marker\n");
                                }
                            }
                        }
                        ylinked_females++;
                        if (omitfp != (FILE *)NULL) {
                            // If we were able to open the omit file, note this female with a Y-linked chromosome...
                            if (Top->pedfile_type == PREMAKEPED_PFT) {
                                fprintf(omitfp, "%s %s %s\n",
                                        Top->PTop[pedi].persons[peri].famname,
                                        unique_id_per(Top->PTop[pedi].persons[peri].uniqueid,
                                                      Top->PTop[pedi].persons[peri].famname),
                                        LTop->Marker[mrkindex].MarkerName);
                            } else {
                                fprintf(omitfp, "%s %s %s\n",
                                        Top->Ped[pedi].Entry[peri].FamName,
                                        unique_id_per(Top->Ped[pedi].Entry[peri].UniqueID,
                                                      Top->Ped[pedi].Entry[peri].FamName),
                                        LTop->Marker[mrkindex].MarkerName);
                            }
                        }
                    }
                } /* peri */
            } /* pedi */
        }
    } /* mrkindex */
    // If we needed the omit file, close it now and tell the user of it's existance...
    // It would not have been opened if we didn't need to use it above...
    if (omitfp != (FILE *) NULL) {
        fclose(omitfp);
        warnvf("An omit file '%s' containing %d Y-linked female genotype(s) to be ignored has been generated!\n",
               omitfl_name, ylinked_females);
        warnf("You may wish to check this against your current omit list.");
        free(omitfl_name);
    }
}

#ifdef DEFUNCT
/*----------------------------------------------------------------------------
  Function to compute necessary field widths for writing the
  annotated pedigree file.
  Modified from field_widths() in output_routines.c
  Returns 1 more than the width in order to leave one column gap.
  ---------------------------------------------------------------------------- */

static void allele_name_width(linkage_locus_top *LTop,
                              ext_linkage_locus_top *EXLTop,
                              size_t *AllWidth)
{
    int m, a;

    *AllWidth=0;
    for (m=0; m < LTop->LocusCnt; m++) {
        if (LTop->Locus[m].Type == NUMBERED) {
            *AllWidth=2;
            for (a=0; a < LTop->Locus[m].AlleleCnt; a++) {
                *AllWidth =
                    ((strlen(EXLTop->EXLocus[m].RAlleles[a]) > *AllWidth)?
                     strlen(EXLTop->EXLocus[m].RAlleles[a]):
                     *AllWidth);
            }
        }
    }
}

static void ann_field_widths(linkage_ped_top *TTop,
                             annotated_file_desc *file_desc)
{

    int  FidWidth, FnameWidth, PnameWidth;
    int  PidWidth;
    size_t AllnameWidth;
    int ProWidth, TwinWidth;
    int col;
    linkage_locus_top *LTop = TTop->LocusTop;
    ext_linkage_locus_top *EXLTop = TTop->EXLTop;

    ped_id_width(TTop, &FidWidth);
    ped_name_width(TTop, &FnameWidth);
    person_name_width(TTop, &PnameWidth);
    person_id_width(TTop, &PidWidth);
    allele_name_width(LTop, EXLTop, &AllnameWidth);
/*   proband_width(TTop, &ProWidth); */
/*   twin_width(TTop, &TwinWidth); */

    /* pedigree file */
    const char *ID = canonicalColName("ID");
    const char *Father = canonicalColName("Father");
    const char *Mother = canonicalColName("Mother");
    const char *Pedigree = canonicalColName("Pedigree");
    const char *PerID = canonicalColName("PerID");
    const char *PedID = canonicalColName("PedID");
    const char *DZTwin = canonicalColName("DZTwin");
    const char *MZTwin = canonicalColName("MZTwin");
    const char *Proband = canonicalColName("Proband");

    for (col=0; col < file_desc->num_ped_cols; col++) {
        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.ID) ||
            SAME(file_desc->ped_file_columns[col].ColName, TKN.Father) ||
            SAME(file_desc->ped_file_columns[col].ColName, TKN.Mother)) {
            file_desc->ped_file_columns[col].output_length = PidWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", PidWidth);
            continue;
        }

        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.Pedigree)) {
            file_desc->ped_file_columns[col].output_length = FidWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", FidWidth);
            continue;
        }

        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.PerID)) {
            file_desc->ped_file_columns[col].output_length = PnameWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", PnameWidth);
            continue;
        }
        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.PedID)) {
            file_desc->ped_file_columns[col].output_length = FnameWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", FnameWidth);
            continue;
        }

        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.DZTwin) ||
            SAME(file_desc->ped_file_columns[col].ColName, TKN.MZTwin)) {
            file_desc->ped_file_columns[col].output_length = TwinWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", TwinWidth);
            continue;
        }

        if (SAME(file_desc->ped_file_columns[col].ColName, TKN.Proband)) {
            file_desc->ped_file_columns[col].output_length = ProWidth;
            sprintf(file_desc->ped_file_columns[col].output_format,
                    "%%%ds", ProWidth);
            continue;
        }
        /* else everything is marker data */
        /* Possible error/problem: casting 'AllnameWidth' to an int from size_t */
        file_desc->ped_file_columns[col].output_length = (int) AllnameWidth;
        sprintf(file_desc->ped_file_columns[col].output_format,
                "%%%ds", (int) AllnameWidth);
    }
}
#endif

static int parse_names_file_header(FILE *filep,
                                   annotated_file_desc *file_desc)
{
    /* Names file has two columns always: Name, Type in any order */
    /* extra information is allowed on each line after the two columns */

    char *nextstr, col_name[READ_CHUNK+1];
    int num_read=0;
    char buffer[READ_CHUNK + 1];

    /* initialize the value type and column names of the two elements */
    file_desc->num_names_cols=2;
    file_desc->names_file_columns = CALLOC((size_t) 2, col_hdr_type);

    INIT_COLNAME(file_desc->names_file_columns, 0, CHAR_AN, "Type");
    INIT_COLNAME(file_desc->names_file_columns, 1, STRING_AN, "Name");

    strcpy(buffer, "");
    IgnoreValue(fgets(buffer, READ_CHUNK, filep));
/*  fcmap(filep, "%=\n", lch); can't mix fcmap & fgets; latter already read entire line */
    nextstr = strtok(buffer, "\t \n");
    while(1) {
        sprintf(col_name, "%s", nextstr);
        num_read++;
        if (!strcasecmp(col_name, "Name")) {
            file_desc->names_file_columns[1].input_col = num_read;

        } else if (!strcasecmp(col_name, "Type")) {
            file_desc->names_file_columns[0].input_col = num_read;
        }
        nextstr = strtok(NULL, "\t \n");
        if (nextstr == NULL) {
            break;
        }
    }

    if (file_desc->names_file_columns[1].input_col == -1 ||
        file_desc->names_file_columns[0].input_col == -1) {
        /* not in names_file format */
        return 0;
    }

    /*  free(buffer); */
    if (is_missing_required_hdr(file_desc->names_file_columns, 2)) {
        errorf("Found errors in names file header, please correct and restart Mega2.");
        EXIT(DATA_INCONSISTENCY);
    }
    return num_read;
}

/* calls read_marker_data() in mrecode.c
   Changes the read_marker_data() function to take in type and
   name columns */

static int read_annotated_names_file(char *names_file,
                                     linkage_locus_top **LTop,
                                     annotated_file_desc *file_desc)
{
    int annotated_format;

    FILE *fp = names_file ? fopen(names_file, "r") : NULL;
    if (fp == NULL) {
        errorvf("could not open %s for reading!\n", names_file);
        EXIT(FILE_READ_ERROR);
    }

#ifndef HIDEFILE
    mssgvf("Reading Mega2 format names file: %s\n", names_file);
#endif
    annotated_format = parse_names_file_header(fp, file_desc);
    if (annotated_format > 0) {
        *LTop = read_marker_data(fp,
                                 file_desc->names_file_columns[0].input_col - 1,
                                 file_desc->names_file_columns[1].input_col - 1);
    } else {
        errorvf("Invalid names file header: %s\n", names_file);
    }
    return annotated_format;
}


/* This is again a variable length header, we can have several maps */

static int parse_map_file_header(FILE *mapfp,
                                 col_hdr_type *reserved_mapcols,
                                 listhandle *userdef_mapcols)
{
    int colnum =
        parse_variable_width_hdr(mapfp,
				 reserved_mapcols,
				 userdef_mapcols,
				 NUM_MAPCOL_NAMES, 2,
				 map_valid_header_item_extensions);
    return colnum;
}


static int get_map_names(int num_cols,
                         col_hdr_type *map_colnames,
                         col_hdr_type *reserved_colnames,
                         const int allocate_additional_maps,
                         ext_linkage_locus_top *EXLTop)
{
    int i, j, map_ind,  mapnum=0, has_err=0;
    char MapName[FILENAME_LENGTH], ext1[FILENAME_LENGTH], ext2[FILENAME_LENGTH];
    char **map_names = CALLOC((size_t)num_cols, char *);
//cpk    sex_map_types **Sex = CALLOC((size_t)num_cols, sex_map_types *);
    int **Sex = CALLOC((size_t)num_cols, int *);
    char *Func = CALLOC((size_t)num_cols, char);
    char Yorn[4];
    
    // First store the map names, map functions, and whether Male or female position, ignoring the reserved column names...
    for (i=0; i < num_cols; i++) {
        // If this is a reserved column, then skip it...
        if (map_colnames[i].value_type == IGNORE ||
            !strcasecmp(map_colnames[i].ColName, reserved_colnames[0].ColName) ||
            !strcasecmp(map_colnames[i].ColName, reserved_colnames[1].ColName) ||
            !strcasecmp(map_colnames[i].ColName, reserved_colnames[2].ColName)) {
            continue;
        }
        // parse_hdr will split the column name in the map file in to three parts...
        // <MapName, ext1, ext2> == <'Map', [k | h | p], [m | f| a | <none>]>
        parse_hdr(map_colnames[i].ColName, MapName, ext1, ext2, map_valid_header_item_extensions);
        // Look to see if we have recorded this MapName in the map_names vector...
        // The user does not have to place the maps for any give name next to each other in the map file...
        map_ind=-1;
        if (mapnum > 0) {
            // look for the map index...
            for (j=0;  j < mapnum; j++) {
                if (strcmp(map_names[j], MapName) == 0) {
                    map_ind=j; break; // Found it...
                }
            }
        }
        
        if (map_ind == -1) {
            // We have not recorded the MapName in the map_names vector, so do it now...
            map_names[mapnum]=strdup(MapName);
            //Func[mapnum]=(char)tolower((unsigned char)ext1[0]);
            // Func[mapnum] should be { 'h' | 'k'} according to ext1[0]...
            switch((unsigned char)ext1[0]) {
                case 'h': // Haldane
                case 'H':
                    Func[mapnum]='h';
                    break;
                case 'k': // Kosambi
                case 'K':
                    Func[mapnum]='k';
                    break;
                case 'p': // Physical Map
                case 'P':
                    Func[mapnum]='p';
                    break;
                case '\0':
                    // This is what we get back from parse_hdr if the map is not of the form 'Mapname.a.b'...
                    Func[mapnum]='\0';
                    break;
                default:
                    // Should this be an error? What about old map names which didn't have to be in the form 'a.b.c'?
                     errorvf("Function '%s' associated with map '%s' is not a genetic map ('h', 'k') or a physical map 'p'.\n",
			     ext1, MapName);
		     EXIT(DATA_TYPE_ERROR);
                    break;
            }
            // Sex[mapnum] is a vector <'a','m','f'> which contains a '1' in the slot according to ext2[0]
            // If this is a physical position 'p', it seems that the map data goes into the SEX_AVERAGED_MAP position...
//cpk            Sex[mapnum]=CALLOC((size_t) 3, sex_map_types);
            Sex[mapnum]=CALLOC((size_t) 3, int);
            //Sex[mapnum][0]=Sex[mapnum][1]=Sex[mapnum][2]=0; // CALLOC already initializes to zero
            map_ind=mapnum;
            mapnum++;
            
        } else if (Func[map_ind] == 'p' && ext1[0] == 'p') {
            // We saw this map before, but you can only have one physical map for a given Mapname.
            // so that means that this map name was used at least twice for a physical map.
            errorvf("Duplicate Physical Map entry for %s\n", map_names[map_ind]);
            EXIT(DATA_TYPE_ERROR);
        }
        
        // Fill in the 'Sex' vector based on the content of 'ext2[0]'...
        // for the map designations see... common.h: sex_map_types
        switch((unsigned char)ext2[0]) {
            case 'm': // Male
            case 'M':
                if (Func[map_ind] == 'h' || Func[map_ind] == 'k') {
                    if (Sex[map_ind][MALE_SEX_MAP] == 1) {
		        errorvf("Duplicate Male sex map entry for genetic map %s.%c\n", 
                                map_names[map_ind], Func[map_ind]);
			EXIT(DATA_TYPE_ERROR);
                    }
                    Sex[map_ind][MALE_SEX_MAP] = 1;
                } else {
		    errorvf("Male sex specified without a genetic map ('h', or 'k') for %s.%c.\n",
                            map_names[map_ind], Func[map_ind]);
		    EXIT(DATA_TYPE_ERROR);
                }
                break;
            case 'f': // Female
            case 'F':
                if (Func[map_ind] == 'h' || Func[map_ind] == 'k') {
                    if (Sex[map_ind][FEMALE_SEX_MAP] == 1) {
		        errorvf("Duplicate Female sex map entry for genetic map %s.%c\n", 
                                map_names[map_ind], Func[map_ind]);
			EXIT(DATA_TYPE_ERROR);
                    }
                    Sex[map_ind][FEMALE_SEX_MAP]= 1;
                } else {
		    errorvf("Female sex specified without a genetic map ('h', or 'k') for %s.%c.\n",
                            map_names[map_ind], Func[map_ind]);
		    EXIT(DATA_TYPE_ERROR);
                }
                break;
            case 'a': // Average
            case 'A':
                if (Func[map_ind] == 'h' || Func[map_ind] == 'k') {
                    if (Sex[map_ind][SEX_AVERAGED_MAP] == 1) {
                        errorvf("Duplicate Sex Average map entry for genetic map %s.%c\n", 
                                map_names[map_ind], Func[map_ind]);
			EXIT(DATA_TYPE_ERROR);
                    }
                    Sex[map_ind][SEX_AVERAGED_MAP] = 1;
                } else {
                    errorvf("Sex Average specified without a genetic map ('h', or 'k') for %s.%c.\n",
                            map_names[map_ind], Func[map_ind]);
		    EXIT(DATA_TYPE_ERROR);
                }
                break;
            case '\0':
                // This is what you get back from parse_hdr if the map is not of the form 'Mapname.a.b'.
                // It will be the canse when Func[mapnum] == 'p'
                break;
            default:
                errorvf("Sex '%s' associated with %s.%s.x not one of ('m','f','a').\n",
			ext2, MapName, ext1);
		EXIT(DATA_TYPE_ERROR);
                break;
        }
    }
    
    // Here we are finished parsing the map_colnames...
    if (mapnum == 0) {
        mssgvf("Map file does not contain any map positions,\n");
        mssgvf("No map information will be stored.\n");
#if 0
	// We really shouldn't be trying to talk with the user if this is batch mode...
        if (InputMode == BATCH_INPUTMODE) {
            warnf("Continuing without map information, analysis options will be limited.");
            return 0;
        }
#endif
        printf("Proceed without a map ? (y/n) ");
        IgnoreValue(scanf("%c", &(Yorn[0]))); fflush(stdin);
        newline;
        if (Yorn[0] == 'y' || Yorn[0] == 'Y') {
            warnf("Continuing without map information, analysis options will be limited.");
            return 0;
        } else {
            errorvf("Exiting mega2 because no marker map provided.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    }
    
    // Copy the map Name (functions and sex if they exist) to the EXLTop structure...
    EXLTop->MapCnt = mapnum;
// +1 null at end
    EXLTop->map_functions = CALLOC((size_t)(mapnum + allocate_additional_maps + 1), char);
    EXLTop->MapNames = CALLOC((size_t)(mapnum + allocate_additional_maps), char*);
//cpk    EXLTop->SexMaps = CALLOC((size_t)mapnum, sex_map_types*);
    EXLTop->SexMaps = CALLOC((size_t)(mapnum+ allocate_additional_maps), int*);
    
    for (i=0; i < mapnum; i++) {
        // cpk: These checks can only be made after all of the header information is read in
        // because MEGA2 allows map infomation to be specified in any column order within the map file...
        if (Func[i] == 'h' || Func[i] == 'k') {
            // The only Genetic map that we don't process is a male only map...
            // See comment header for function user_input.c:get_genetic_distance_index()
            if (Sex[i][SEX_AVERAGED_MAP] == 0 && Sex[i][FEMALE_SEX_MAP] == 0 && Sex[i][MALE_SEX_MAP] == 1) {
                warnvf("The Genetic Distance map %s.%c will be ignored because only a male map was supplied.\n",
                       map_names[i], Func[i]);
            }
        }
        /* store the map name, and function */
        EXLTop->MapNames[i]=strdup(map_names[i]);
        EXLTop->map_functions[i]=Func[i];
	// For a physical map, the data is stored in the SEX_AVERAGED_MAP of SexMaps...
//cpk        EXLTop->SexMaps[i]=CALLOC((size_t) 3, sex_map_types);
        EXLTop->SexMaps[i] = (int *)CALLOC((size_t) 3, int);
        EXLTop->SexMaps[i][SEX_AVERAGED_MAP] = Sex[i][SEX_AVERAGED_MAP];
        EXLTop->SexMaps[i][MALE_SEX_MAP] = Sex[i][MALE_SEX_MAP];
        EXLTop->SexMaps[i][FEMALE_SEX_MAP] = Sex[i][FEMALE_SEX_MAP];
        free(map_names[i]);
        free(Sex[i]);
    }
    free(map_names); free(Func); free(Sex);
    
    // Now parse the header from the map file all over again to fill in the other
    // places in the 'map_colnames' structure.
    for (i=0; i < num_cols; i++) {
        // If this is a reserved colum then skip it...
        if (map_colnames[i].value_type == IGNORE ||
	    // These should be the SAME THING that was copied to 'reserved_colnames' in read_annotated_mam_file() and used above.
	    // but here the actual values are used here?!
            !strcasecmp(map_colnames[i].ColName, "Chromosome") ||
            !strcasecmp(map_colnames[i].ColName, "Name") ||
            !strcasecmp(map_colnames[i].ColName, "Error")) {
            if (!strcasecmp(map_colnames[i].ColName, "Error")) {
                has_err = 1;
            }
            continue;
        }
	// Why are we doing this again in the same routine????
        parse_hdr(map_colnames[i].ColName, MapName, ext1, ext2, map_valid_header_item_extensions);
        // Look to see if we have recorded this MapName in the map_names vector...
        for (j=0; j < mapnum; j++) {
            if (strcmp(EXLTop->MapNames[j], MapName) == 0) {
                break; // Found it...
            }
        }
        // Set the map_number in map_colnames according to how the name was first found in the map file...
        map_colnames[i].map_number=j;
        // Also record the index used in this map for the sex map (data)...
        if (ext2[0] == '\0') {
            // This should be a physical map. Again the data for a physical map is stored
            //in the SEX_AVERAGED_MAP slot...
           map_colnames[i].sex_map_number = SEX_AVERAGED_MAP; // 0
           msgvf("Input Map name: %s, type: physical map\n", MapName);
        } else {
	    // for the map designations see... common.h: sex_map_types
            const char *sexx = "";
            switch(tolower((unsigned char)ext2[0])) {
                case 'm':
                    map_colnames[i].sex_map_number=MALE_SEX_MAP;
                    sexx = "male";
                    break;
                case 'f':
                    map_colnames[i].sex_map_number=FEMALE_SEX_MAP;
                    sexx = "female";
                    break;
	        case 'a':
                default:
                    map_colnames[i].sex_map_number=SEX_AVERAGED_MAP;
                    sexx = "average";
                    break;
            }
            const char *units = "centiMorgans";
            if (Input_Format == in_format_binary_PED || Input_Format == in_format_PED) {
                units = "Morgans";
                if (strcasecmp("Map", MapName) == 0)
                    if (PLINK.cM)
                        units = "centiMorgans";
            }
            msgvf("Input Map name: %s, type: %s genetic map, units: %s %s\n",
                  MapName, sexx,
                  ((tolower((unsigned char)ext1[0]) == 'h') ? "haldane" : "kosambi"),
                  units);
        }
        /* write all the sex-maps for the same map name together */
        map_colnames[i].output_col=
        // I don't understand this because 'map_number' for a genetic of physical map could be '0'?
        2 + has_err +
        // this gives that order that allows the maps to be written out together,
        // if there are multiple columns for a map.
        3*map_colnames[i].map_number + map_colnames[i].sex_map_number;
    }
    
    return mapnum;
}

void copy_exmap_locmap(linkage_locus_top *LTop,
                       ext_linkage_locus_top *EXLTop,
                       int map_num)
{
    int m;

    /* have to insert a selection menu here */
/*   if (EXLTop->MapCnt > 1) { */
/*     map_num = select_output_map(EXLTop); */
/*   } */

    LTop->map_distance_type = EXLTop->map_functions[map_num];
    
    if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
      // this will output two vectors of male and female in the recombination fractions...
      if (LTop->SexLinked == 2) {
      /* On the X chromosome, the LINKAGE sex difference parameter
         must be set to NO_SEX_DIFF. */
        LTop->SexDiff = NO_SEX_DIFF;
      } else {
        LTop->SexDiff = VARIABLE_SEX_DIFF;
      }
    } else {
      // this will output just one vector...
      LTop->SexDiff = NO_SEX_DIFF;
    }
    
#ifdef USEOLDMAPCODE
    for (m = LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
// hmm EXLTop -> Marker
        if (EXLTop->EXLocus[m].positions == NULL) {
            // Since there is no sex-averaged position set sex-specific positions as well
            // to unknown...
            LTop->Marker[m].position =
                LTop->Marker[m].pos_male =
                LTop->Marker[m].pos_female = UNKNOWN_POSITION;
        } else if (LTop->map_distance_type == 'p') {
            /* simply copy the average map position, without testing
               for chromosome number */
            LTop->Marker[m].position = EXLTop->EXLocus[m].positions[map_num];
        } else {
            if (LTop->Marker[m].chromosome == SEX_CHROMOSOME) {
                /* set male position to NA */
                if (EXLTop->SexMaps[map_num][MALE_SEX_MAP] == 1) {
                    LTop->Marker[m].pos_male = UNKNOWN_POSITION;
                }
                /* set the average and female positions */
                if (EXLTop->SexMaps[map_num][FEMALE_SEX_MAP] == 1 &&
                    EXLTop->EXLocus[m].pos_female[map_num] >= 0.0) {
                    LTop->Marker[m].pos_female =
                        LTop->Marker[m].position = EXLTop->EXLocus[m].pos_female[map_num];
                } else {
                    LTop->Marker[m].pos_female =
                        LTop->Marker[m].position = EXLTop->EXLocus[m].positions[map_num];
                }
            } else if (LTop->Marker[m].chromosome == MALE_CHROMOSOME) {
                if (EXLTop->SexMaps[map_num][MALE_SEX_MAP] == 1) {
                    LTop->Marker[m].pos_male =
                        LTop->Marker[m].position = EXLTop->EXLocus[m].pos_male[map_num];
                } else {
                    LTop->Marker[m].pos_male =
                        LTop->Marker[m].position = EXLTop->EXLocus[m].positions[map_num];
                }

                LTop->Marker[m].pos_female = UNKNOWN_POSITION;
            } else {
                /* Copy over all the maps that are defined */
                LTop->Marker[m].position = EXLTop->EXLocus[m].positions[map_num];
                if (EXLTop->SexMaps[map_num][MALE_SEX_MAP] == 1) {
                    LTop->Marker[m].pos_male = EXLTop->EXLocus[m].pos_male[map_num];
                } else {
                    LTop->Marker[m].pos_male = EXLTop->EXLocus[m].positions[map_num];
                }
                if (EXLTop->SexMaps[map_num][FEMALE_SEX_MAP] == 1) {
                    LTop->Marker[m].pos_female = EXLTop->EXLocus[m].pos_female[map_num];
                } else {
                    LTop->Marker[m].pos_female = EXLTop->EXLocus[m].positions[map_num];
                }
            }
            /* For X-linked loci, if there is no average map, copy the Sex-map */
            if (LTop->Marker[m].chromosome == SEX_CHROMOSOME && LTop->Marker[m].position < 0 &&
                LTop->Marker[m].pos_female >= 0.0) {
                LTop->Marker[m].position = LTop->Marker[m].pos_female;
            }
        }
    }
    if (EXLTop->SexMaps[map_num][MALE_SEX_MAP] == 1 && EXLTop->SexMaps[map_num][FEMALE_SEX_MAP] == 1) {
        HasSexMaps=1;
    } else {
        HasSexMaps=0;
    }

    if (EXLTop->SexMaps[map_num][SEX_AVERAGED_MAP] == 1) {
        HasAvgMaps = 1;
    } else {
        HasAvgMaps = 0;
    }
#else /* USEOLDMAPCODE */
    
    for (m = LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
        if (LTop->Locus[m].Type == AFFECTION || 
            LTop->Locus[m].Type == QUANT) {
            // This will be the case when this Locus is a trait... which can never be
            LTop->Marker[m].pos_avg =
            LTop->Marker[m].pos_male =
            LTop->Marker[m].pos_female = UNKNOWN_POSITION;
        } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT ||
                   genetic_distance_sex_type_map == FEMALE_GDMT) {
            // Here the user chose a sex specific map...
            LTop->Marker[m].pos_avg = UNKNOWN_POSITION;
            LTop->Marker[m].pos_female = EXLTop->EXLocus[m].pos_female[map_num];
            if (genetic_distance_sex_type_map == FEMALE_GDMT) {
                LTop->Marker[m].pos_male = UNKNOWN_POSITION;
            } else {
                LTop->Marker[m].pos_male = EXLTop->EXLocus[m].pos_male[map_num];
            }

	    // NOTICE: This is do so there will be a sex-averated position
	    // if we missed something in the code. We need to have a review to
	    // catch these...
	    if (EXLTop->SexMaps[map_num][SEX_AVERAGED_MAP] == 1)
	      LTop->Marker[m].pos_avg = EXLTop->EXLocus[m].positions[map_num];

        } else if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            LTop->Marker[m].pos_avg = EXLTop->EXLocus[m].positions[map_num];
            // In this case, there are no sex specific positions...
            LTop->Marker[m].pos_male =
            LTop->Marker[m].pos_female = UNKNOWN_POSITION;
        } else {
            // We will not come into this function with a physical map....
            errorvf("Internal, no genetic distance sex type map is available.\n");
            EXIT(SYSTEM_ERROR);
        }
    }
#endif /* USEOLDMAPCODE */
}

extern int exceeded_max_morgan_value(const double position);

static void init_ext_linkage_locus_top (ext_linkage_locus_top *EXLTop,
                                        linkage_locus_top *LTop)
{
    // This is actually a bit bigger than it needs to be because it also contains phenotype locus...
    EXLTop->LocusCnt = LTop->LocusCnt;
    // Note that CALLOC initializes the storage to '0' (e.g., clears them)...
    EXLTop->EXLocus = (LTop->MarkerCnt > 0) ? (CALLOC((size_t) LTop->MarkerCnt, ext_linkage_locus_rec) - LTop->PhenoCnt) : 0;
#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: EXLocus: %d MB (%d x %d)\n",
          LTop->MarkerCnt * sizeof(ext_linkage_locus_rec) / 1024 / 1024,
          LTop->MarkerCnt,  sizeof(ext_linkage_locus_rec));
#endif
}

/**
   @brief Create a new Extended Linkage Locus Top on the heap optionally initializing it

   @param LTop       If 'LTop' is not NULL then initialize as well.
   @return ext_linkage_locus_top *
 */
ext_linkage_locus_top *new_EXLTop(linkage_locus_top *LTop)
{
    ext_linkage_locus_top *EXLTop = MALLOC(ext_linkage_locus_top);
    
    EXLTop->GroupCnt=1; // This is never read
    /*   EXLTop->Locus = NULL; */
    EXLTop->map_functions = NULL;
    EXLTop->MapNames = NULL;
    EXLTop->MapCnt = 0; // This is filled from the results of get_map_names()
    EXLTop->SexMaps = NULL;
    
    if (LTop) init_ext_linkage_locus_top(EXLTop, LTop);
    else {
        EXLTop->LocusCnt = 0;
        EXLTop->EXLocus = NULL;
    }
    
    return EXLTop;
}

/**
   @brief Create and initialize the External Linkage Locus Top Positions for one marker

   @param EXLTop     A pointer to the External Linkage Locus Top structure.
   @param mrk_num    An index into the EXLocus vector.
   @param num_maps   The number of maps to create entries for (int > 0).
   @return void
 */
static void new_EXLTop_positions(ext_linkage_locus_top *EXLTop,
                                 const int mrk_num,
                                 const int num_maps)
{
    int l;

    if (num_maps <= 0) {
      errorvf("INTERNAL: Cannot create map with zero entries.\n");
      EXIT(DATA_INCONSISTENCY);
    }
    
    // Build the vectors....
    EXLTop->EXLocus[mrk_num].positions = CALLOC((size_t)num_maps, double);
    EXLTop->EXLocus[mrk_num].pos_male = CALLOC((size_t)num_maps, double);
    EXLTop->EXLocus[mrk_num].pos_female = CALLOC((size_t)num_maps, double);

    // Initialize the vectors...
    for (l = 0; l < num_maps; l++) {
        EXLTop->EXLocus[mrk_num].positions[l] =
        EXLTop->EXLocus[mrk_num].pos_female[l] =
        EXLTop->EXLocus[mrk_num].pos_male[l] = UNKNOWN_POSITION;
    }
}

/**
   @brief Create entries for markers without positions
 
   Markers are processed before the map file is, and our internal loci structures contain
   chromosome position information. It is valid to exclude markers from the map file
   that have been specified earlier (e.g., in the ped file). However, since the chromosome position
   information structures are built in they will not have been built in those cases. Other
   parts of the code assume their existance. So, after the map file has been read, we need to
   loop over all of the markers and check to see if it is associated with a valid chromosome.
   If not, we patch in (as unknown) chromosome position information.

   @param LTop       If 'LTop' is not NULL then initialize as well.
   @param EXLTop     A pointer to the External Linkage Locus Top structure.
   @param num_maps   The number of maps to create entries for (int > 0).
   @return           The count of markers missing from the map file
*/
static int create_entries_for_markers_without_positions(linkage_locus_top *LTop,
                                                        ext_linkage_locus_top *EXLTop,
                                                        const int num_maps)
{
    int i, mrk_missing_from_map = 0;
    SECTION_ERR_INIT(locus_dropped);
    
    for (i = 0; i < LTop->LocusCnt; i++) {
        /* Check if maps are provided for all numbered loci */
        if (LTop->Locus[i].Class == MARKER) {
            if (LTop->Locus[i].number < 0) {
                sprintf(err_msg, "Locus %s is not in map file; flushed.", LTop->Locus[i].LocusName);
                SECTION_ERR(locus_dropped);
                warnf(err_msg);
                LTop->Marker[i].chromosome = MISSING_CHROMO;
                mrk_missing_from_map++;
                new_EXLTop_positions(EXLTop, i, num_maps);
            }
        }
    }
#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: EXLocus maps: %d MB (%d x %d)\n",
          LTop->MarkerCnt * num_maps * 3 * 8 / 1024 / 1024,
          LTop->MarkerCnt,  num_maps * 3 * 8);
#endif

    SECTION_ERR_FINI(locus_dropped);

    return mrk_missing_from_map;
}


static ext_linkage_locus_top *read_common_map_file(FILE *mapfp,
                                                   const char *map_file,
                                                   linkage_locus_top *LTop,
                                                   const int allocate_additional_maps,
                                                   col_hdr_type reserved_colnames[NUM_MAPCOL_NAMES],
                                                   annotated_file_desc *file_desc,
                                                   char **alleles)
{
    int ii = 1, line = 0, j, l, m;
    col_hdr_type *map_all_colnames = &(file_desc->map_file_columns[0]);
    int *colnums; /* to store the column number for each position */
    char dummy[FILENAME_LENGTH], dname[FILENAME_LENGTH];
    int mrk_num, chr=0;
    double error_prob, *positions;
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    int **valid_map_p;
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
    int num_maps;
    int num_positions, num_read;
    int lch, has_error; /* simulated genotyping error column? */
    SECTION_ERR_INIT(chrm_errors);
    SECTION_ERR_INIT(max_morgans);
    SECTION_ERR_INIT(locus_extra);
    SECTION_ERR_INIT(dup_locus);
    SECTION_ERR_INIT(bad_position);
    
    /* These are all error flags */
    int mrk_missing_from_map, mrk_missing_from_names=0, duplicate_mrk=0;
    int missing_columns=0, bad_chromo=0, bad_allele=0;
    int display_Y_message=1;
    char yesorno[4];
    int alleles_i = 0;
    int total_maps_to_allocate;

    ext_linkage_locus_top *EXLTop = new_EXLTop(LTop);
    
    // Get the names of maps from the headers.
    // These will be Genetic and Physical maps, and will not include trait markers!
    EXLTop->MapCnt = num_maps =
        get_map_names(PLINK.plink == binary_PED_format ? file_desc->num_map_cols - 2 : file_desc->num_map_cols,
                      map_all_colnames,
                      reserved_colnames,
                      allocate_additional_maps,
                      EXLTop);
    
    // Will this cause problems later in the code since 'create_entries_for_markers_without_positions()' is never called???
    if (num_maps == 0) {
#ifndef HIDESTATUS
        mssgf("Done reading, no valid maps found!");
#endif
        return EXLTop;
    }

    // Maps from the Mega2 map file and maps to be filled in at a later time...

    total_maps_to_allocate = num_maps + allocate_additional_maps;
    
#ifndef HIDESTATUS
    mssgvf("Found %d possible maps in the %s file.\n", num_maps, map_file);
#endif
    printf("Now checking each record in map file %s ...\n", map_file);
    
    human_x = human_y = human_xy = human_unknown = human_mt = human_auto = 0;
    
    colnums = CALLOC((size_t)file_desc->num_map_cols, int);
    positions = CALLOC((size_t)file_desc->num_map_cols, double);
  
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    // Here we keep track of whether the map contains any valid map data (e.g., not all '0')...
    // As of Bug 338 - 'Valid map': a genetic map with positions of all zeros is valid
    // This constraint is relaxed
    valid_map_p = CALLOC((size_t)total_maps_to_allocate, int *);
    for (i = 0; i < total_maps_to_allocate; i++) 
        valid_map_p[i] = CALLOC((size_t)sizeof(sex_map_types), int);
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */

    
    /* Now read in the rest of the file, the part that contains the data */
    /* NOTE: not worth "canonicalName-ifying these tokens (dummy) because they are
       typically different fp numbers.  (unless they are NA or 0)
     */
    MaxChromo = NumUnmapped = 0;
    int a1 = 0, a2 = 0;
    if (PLINK.plink == binary_PED_format && alleles != (char **)NULL) {
        if (PLINK.map3) {
            a1 = 3; a2 = 4;
        } else {
            a1 = 4; a2 = 5;
        }
    }
    while (!feof(mapfp)) {
        line++;
        error_prob = 0.0;
        for (m=0; m < file_desc->num_map_cols; m++) {
            positions[m] = UNKNOWN_POSITION;
        }
        // Throw away comment lines...
        if ((chr = fgetc(mapfp)) == -1) {
            break;
        } else if (chr == '#') {
            skip_line(mapfp, chr);
            continue;
        } else if (chr == '\n') {
            continue;
        } else {
            ungetc(chr, mapfp);
        }
        
        m=0;
        num_read = 0;
        for (j=0; j < file_desc->num_map_cols; j++) {
            /* read in the column if not ignore */
            if (map_all_colnames[j].value_type == IGNORE) {
                lch=fcmap(mapfp, "%s", dummy);
            } else if (!strcasecmp(map_all_colnames[j].ColName, "Name")) {
                lch=fcmap(mapfp, "%s", dname);
            } else if (!strcasecmp(map_all_colnames[j].ColName, "Chromosome")) {
                lch = fcmap(mapfp, "%s", dummy);
                if ((chr = STR_CHR(dummy)) == -1) {
                    sprintf(err_msg, "Chromosome number %s on line %d is invalid.",
                            dummy, line);
                    SECTION_ERR(chrm_errors);
                    errorf(err_msg);
                    bad_chromo++;
                } else if (chr == 0) {
                    // are we handeling this right for PLINK format where '0' represents "unplaced"?
                    chr=UNKNOWN_CHROMO;
                }
                
                MaxChromo = ((chr != UNKNOWN_CHROMO && chr > MaxChromo)? chr : MaxChromo);
            } else if (!strcasecmp(map_all_colnames[j].ColName, "Error")) {
                lch=fcmap(mapfp, "%s", dummy); error_prob = atof(dummy);
            } else {
                lch=fcmap(mapfp, "%s", dummy);
                // This is for processing the allele columns in a PLINK .bim file...
                if (PLINK.plink == binary_PED_format && alleles != (char **)NULL && j>=a1 && j<=a2) {
                    if (strlen(dummy) == 1) {
                        (*alleles)[alleles_i++] = dummy[0];
                    } else {
                        (*alleles)[alleles_i++] = '*';
                        if (!strcmp(dummy, "dummy") || !strcmp(dummy, "dummy1") || !strcmp(dummy, "dummy2") ) {
                            // leave it alone
                        } else {
                            errorvf("For a PLINK .bim file; each allele (\"%s\") must be only one character!\n", dummy);
                            bad_allele++;
                        }
                    }
                } else {
                    if (!strcmp(dummy, "NA")) {
                        positions[m] = UNKNOWN_POSITION;
                    } else {
                        positions[m] = atof(dummy);
                        // This check should only be done for the 3rd column of a .map or .bim
                        // file where the Genetic Distance is found (e.g. Morgans or Centimorgans)
                        // PLINK expects the 3rd column the MAP/BIM file to contain genetic distances
                        // in Morgan units. If the genetic distances are in cM instead of Morgans,
                        // add the --cm flag (PLINK.cM == 1).
                        // Internally MEGA2 uses cM so convert from Morgans here if applicable...
                        // This is not done if this is a 3 column map (e.g., --map3)...
                        if (PLINK.plink != not_plink_format && j==2 && PLINK.cM != 1 && PLINK.map3 != 1) {
                            if (exceeded_max_morgan_value(positions[m])) {
                                SECTION_ERR(max_morgans);
                                warnvf("Line %d of the map file contains a genetic distance value of %lf in Morgans.\n",
                                       line, positions[m]);
                            }
                            positions[m] *= 100.0; // Convert M to cM. There are 100 cM in 1 M.
                        }
                        // To exclude a SNP from analysis, set the 4th column (physical base-pair position)
                        // to any negative value (only work for MAP files, not for binary BIM files).
                        if (PLINK.plink && j==3) {
                            if (positions[m] < 0.0) {
                                if (PLINK.plink == binary_PED_format) {
                                    // This should be a fatal error...
                                    errorvf("Line %d; For a PLINK .bim file the base-pair position cannot be excluded from analysis (e.g. made negative).\n", line);
                                    EXIT(FILE_READ_ERROR);
                                }
                            } else {
                                PLINK.markers_included += 1;
                            }
                        }
                    }
                    // holds the index of columns that are either traits or positions (genetic/physical)...
                    colnums[m] = j;
                    m++;
                }
            }
            num_read++;
            if ((lch == '\n' || lch == '\r')  &&
                num_read < (file_desc->num_map_cols-1)) {
                warnvf( "Line %d has too few columns (expected %d, read %d).\n",
                        line, file_desc->num_map_cols, num_read);
                missing_columns++;
            }
        } // for (j=0; j < file_desc->num_map_cols; j++) {

        // The number of columns in the file that contain map information...
        num_positions = m;
        
        //
        // HERE WE HAVE JUST READ ONE LINE/ROW OF DATA FROM THE FILE...
        //
        
        if (chr == UNKNOWN_CHROMO) {
            for (l=0; l < num_positions; l++) positions[l] = UNKNOWN_POSITION;
        }
        
        /* Now check the chromo and marker name, store the positions */
        mrk_num=-1;
        // This is a hash lookup of the names (markers and phenotypes) in LTop
        // Here we look up the marker that we have just found from the .map file.
        // This tells us if we have an intersection with LTop and the map file.
        search_marker(dname, &mrk_num);
        if (mrk_num == -1) {
            // No intersection...
            // Here it's "extra". It's in the map file but not in the data used to make LTop.
            // NOTE: that if the 'dname' is not in LTop then we never get here...
            /* Locus is not in locus file, warn and skip */
            SECTION_ERR(locus_extra);
            if (Input->xcf) 
                warnvf("Locus %s (line %d in %s) has been filtered from the VCF file; ignoring this locus.\n",
                       dname, line, map_file);
            else if (Input->impute)
                warnvf("Locus %s (line %d in %s) extra locus data in this file; ignoring this locus.\n",
                       dname, line, map_file);
            else
                warnvf("Locus %s (line %d in %s) is not in the locus file; ignoring this locus.\n",
                       dname, line, map_file);
            mrk_missing_from_names++;
        } else {
            // Intersection...
            if (LTop->Locus[mrk_num].number > 0) {
                /* already encountered */
                sprintf(err_msg,
                        "Duplicate marker name in %s on line %d.",
                        dname, line);
                SECTION_ERR(dup_locus);
                errorf(err_msg);
                duplicate_mrk++;
            } else {
                /* first encounter */
                LTop->Locus[mrk_num].number = ii++;
                if (chr == SEX_CHROMOSOME) {
                    human_x++;
                    if (LTop->Locus[mrk_num].Type != XLINKED) {
                        LTop->Locus[mrk_num].Type = XLINKED;
                    }
                } else if (chr == MALE_CHROMOSOME) {
                    human_y++;
                    if (LTop->Locus[mrk_num].Type != YLINKED) {
                        LTop->Locus[mrk_num].Type = YLINKED;
                    }
                } else if (chr == PSEUDO_X) {
                    human_xy++;
                } else if (chr == MITO_CHROMOSOME) {
                    human_mt++;
                    // PLINK allows the value of '0' to represent "unplaced"
                } else if (chr == UNKNOWN_CHROMO) {
                    human_unknown++;
                    NumUnmapped++;
                } else {
                    human_auto++;
                }
                
                LTop->Marker[mrk_num].chromosome = chr;
                LTop->Marker[mrk_num].error_prob = error_prob;

                // Create the data structure so that it can be populated with map data (below from this file)...
                // Create an extra slot at the end for the VCF map (to be filled in later from data found in
                // the VCF file) if we are processing a VCF file.
                new_EXLTop_positions(EXLTop, mrk_num, total_maps_to_allocate);
                
                // num_posiions is the number of maps in the input record...
                for (l=0; l < num_positions; l++) {
                    if (chr != UNKNOWN_CHROMO) {
                        // The index into colnums[], and posiitons[] are relative to the maps as they appear in the input record.
                        // So, if the first map is in position 2 (third column of the file), then colnums[0] == 2.
                        // this is used with map_all_colnames[] because it is relative to the beginning of the input record.
                        // So, map_all_colnames[colnums[0]] get's data for the first map in the input record without regard to what
                        // column it appeared in in the input record.
                        // However, positions[0] get's the value of the first map (e.g., Map.h.f, or Map.p).
                        
                        // p (map_number) is relative to the map 'names' as they are encountered in the input
                        // record (e.g., MapA.h.m,MapA.h.f,MapA.h.a == 0, MapB.?.? == 1) where MapA and MapB have unique map numbers.
                        int p = map_all_colnames[colnums[l]].map_number;
                        double position = positions[l];
                        
                        j = map_all_colnames[colnums[l]].sex_map_number;
                        
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
                        // mark this map as containing valid data...
                        if (position != 0.0) valid_map_p[p][j] = 1;
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */

                        // a map (one column from the map file) can only be of one type.
                        // though a [genetic] MapName my have as many as three columns if it is a genetic map.
                        
                        // Chromosome Name Map.k.f Map.k.m MapB.p Map.k.a MapA.h.a Bp1.p Bp2.p Bp3.p
                        switch(j) {
                            case SEX_AVERAGED_MAP:
                                // IMPORTANT NOTE:
                                // The data for a physical map (e.g., Map.p) will go into the SEX_AVERAGED_MAP
                                // column (e.g, sex_map_number == 0)
                                EXLTop->EXLocus[mrk_num].positions[p] = position;
                                break;
                            case MALE_SEX_MAP:
                                EXLTop->EXLocus[mrk_num].pos_male[p] = position;
                                break;
                            case FEMALE_SEX_MAP:
                                EXLTop->EXLocus[mrk_num].pos_female[p] = position;
                                break;
                            default:
                                /* should not happen */
                                break;
                        }
                    }
                }

                /* Next check that the correct positions have been specified for each map */
                for (l=0; l < num_maps; l++) {
                    has_error=0;
                    strcpy(err_msg, "");
                    if (EXLTop->map_functions[l] == 'p') {
                        if (EXLTop->EXLocus[mrk_num].positions[l] > 0.0) {
                            ;
                        } else {
                            sprintf(err_msg,
                                    "Marker %s (map %s line %d): Missing physical map position.",
                                    dname, EXLTop->MapNames[l], line);
                            has_error=1;
                        }
                        continue;
                    }
                    
//#if 0 /* no chr check */
                    // This is no longer done, because the user is forced to pick a
                    // sex-averaged, or sex-specific map later in the program.

                    /* else condition on chomosome number rather than marker type
                     to check the genetic maps. */
                    if (chr == SEX_CHROMOSOME) {
                            if (EXLTop->EXLocus[mrk_num].pos_female[l] >= 0.0) {
                                /* case 1: female map is provided */
                                /* other map positions should not be specified */
                                if (EXLTop->EXLocus[mrk_num].pos_male[l] >= 0.0 ||
                                    EXLTop->EXLocus[mrk_num].positions[l] >= 0.0) {
                                    //has_error=1;
                                    //sprintf(err_msg,
                                    //      "XLINKED %s (map %s line %d): using female map, ignoring rest.",
                                    //      dname, EXLTop->MapNames[l], line);
                                }
                            } else {
                                /* case 2: female map is negative or NA */
                                has_error=(PLINK.plink) ? 0: 1; // only avg gen position
                                if (fabs(EXLTop->EXLocus[mrk_num].pos_female[l] - UNKNOWN_POSITION) <= EPSILON) {
                                    // MAP MESSAGES ONLY ONCE!!!
                                    sprintf(err_msg,
                                            "XLINKED %s (map %s, line %d): Missing female position",
                                            dname, EXLTop->MapNames[l], line);
                                } else if (EXLTop->EXLocus[mrk_num].pos_female[l] < 0.0) {
                                    sprintf(err_msg,
                                            "XLINKED %s (map %s, line %d): Negative female position",
                                            dname, EXLTop->MapNames[l], line);
                                }
                                /* Can use the average map, if provided */
                                if (EXLTop->EXLocus[mrk_num].positions[l] < 0.0) {
                                    EXLTop->EXLocus[mrk_num].pos_female[l] = UNKNOWN_POSITION;
                                    strcat(err_msg, " with no average position given.");
                                    errorf(err_msg);
                                    EXIT(DATA_INCONSISTENCY);
                                } else {
                                    EXLTop->EXLocus[mrk_num].pos_female[l] =
                                    EXLTop->EXLocus[mrk_num].positions[l];
                                    strcat(err_msg, ", using average position.");
                                }
                            }
                    } else if (chr == MALE_CHROMOSOME) {
                            if (EXLTop->EXLocus[mrk_num].pos_male[l] >= 0.0) {
                                if (display_Y_message) {
                                    warnf("Male positions for Y will be used only for ordering loci.");
                                    display_Y_message = 0;
                                }
                                /* case 1: male map is provided */
                                /* other map positions should not be specified */
                                if (EXLTop->EXLocus[mrk_num].pos_female[l] >= 0.0 ||
                                    EXLTop->EXLocus[mrk_num].positions[l] >= 0.0) {
                                    has_error=1;
                                    sprintf(err_msg,
                                            "YLINKED %s (map %s line %d): using male map, ignoring rest.",
                                            dname, EXLTop->MapNames[l], line);
                                }
                            } else {
                                /* case 2: male map is negative or NA */
                                has_error=(PLINK.plink) ? 0: 1; // only avg gen position
                                if (fabs(EXLTop->EXLocus[mrk_num].pos_male[l] - UNKNOWN_POSITION) <= EPSILON) {
                                    sprintf(err_msg,
                                            "YLINKED %s (map %s, line %d): Missing male position",
                                            dname, EXLTop->MapNames[l], line);
                                } else if (EXLTop->EXLocus[mrk_num].pos_male[l] < 0.0) {
                                    sprintf(err_msg,
                                            "YLINKED %s (map %s, line %d): Negative male position",
                                            dname, EXLTop->MapNames[l], line);
                                }
                                /* Can use the average map, if provided */
                                if (EXLTop->EXLocus[mrk_num].positions[l] < 0.0) {
                                    EXLTop->EXLocus[mrk_num].pos_male[l] = UNKNOWN_POSITION;
                                    strcat(err_msg, " with no average position given.");
                                } else {
                                    if (display_Y_message) {
                                        warnf("Average positions for Y will be used only for ordering loci.");
                                        display_Y_message = 0;
                                    }
                                    EXLTop->EXLocus[mrk_num].pos_male[l] =
                                    EXLTop->EXLocus[mrk_num].positions[l];
                                    strcat(err_msg, ", using average position.");
                                }
                            }
                    } else {
                            /* 	case NUMBERED : */
                            /* all positions can be non-zero */
                            if (EXLTop->EXLocus[mrk_num].positions[l] < 0.0) {
                                if (fabs(EXLTop->EXLocus[mrk_num].positions[l] - UNKNOWN_POSITION) <= EPSILON) {
                                    sprintf(err_msg,
                                            "AUTO %s (map %s, line %d): Missing average position ",
                                            dname, EXLTop->MapNames[l], line);
                                } else if (EXLTop->EXLocus[mrk_num].positions[l] < 0.0) {
                                    sprintf(err_msg,
                                            "AUTO %s (map %s, line %d): Negative average position.",
                                            dname, EXLTop->MapNames[l], line);
                                }
                                if (EXLTop->EXLocus[mrk_num].pos_male[l] >= 0.0 &&
                                    EXLTop->EXLocus[mrk_num].pos_female[l] >= 0.0) {
                                    /* don't flag as an error */
                                    EXLTop->EXLocus[mrk_num].positions[l] =
                                    (EXLTop->EXLocus[mrk_num].pos_male[l] + EXLTop->EXLocus[mrk_num].pos_female[l])/2.0;
                                    
                                } else {
                                    has_error=1;
                                    EXLTop->EXLocus[mrk_num].positions[l] =
                                    EXLTop->EXLocus[mrk_num].pos_female[l] =
                                    EXLTop->EXLocus[mrk_num].pos_male[l] = UNKNOWN_POSITION;
                                    strcat(err_msg, " with no sex-specific positions given.");
                                }
                            }
                    }
                    
                    if (has_error) {
                        SECTION_ERR(bad_position);
                        warnf(err_msg);
                    }
//#endif 0 /* no chr check */
                }
            } /* new locus, i.e. not duplicated */
        } /* valid locus present in names file */
    } // while (!feof(mapfp)) {

    //
    // THIS IS THE PLACE to compact EXLTop and LTop because you have ???????? in question
    // just processed the map file and so you know the final set of markers that
    // you will be using for this run of Mega2.

    SECTION_ERR_FINI(chrm_errors);
    SECTION_ERR_FINI(max_morgans);
    SECTION_ERR_FINI(locus_extra);
    SECTION_ERR_FINI(dup_locus);
    SECTION_ERR_FINI(bad_position);
    Display_Errors = 1;

#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    // valid_map_p is a predicate that tells if the map is valid...
    EXLTop->valid_map_p = valid_map_p;
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
    
    LTop->SexLinked = (human_x + human_y) ? ((human_auto + human_xy) ? 2 : 1) : 0;
    /*
      2 => autosomes and x &/or y read
      1 => only x &/or y
      0 => autosomes
    */
    fclose(mapfp);
    free(positions);
    free(colnums);

    // Assign position data to those markers not listed in the map file....
    mrk_missing_from_map = create_entries_for_markers_without_positions(LTop, EXLTop, total_maps_to_allocate);

    if (bad_chromo) {
        warnvf("Invalid chromosome numbers in map file, see %s for details.\n",
              Mega2Err);
    }
    
    if (missing_columns) {
        warnvf("Found lines with less than required number of fields, see %s for details.\n",
              Mega2Err);
    }
    
    if (mrk_missing_from_map) {
        warnvf("Some marker loci in the %s file are missing from map file, see %s for details.\n",
              (Input->xcf ? "VCF" : (Input->impute ? "Impute2" : "names")),
              Mega2Err);
    }
    
    if (mrk_missing_from_names) {
        warnvf("Some marker loci in map file %s file, see %s for details.\n",
               (Input->xcf ? "have been filtered from the VCF" : 
                (Input->impute ? "are missing from Imputed" :
                 "are missing from the names")),
               Mega2Err);
    }
    
    /* terminate or continue */
    if (missing_columns || duplicate_mrk || bad_allele) {
        EXIT(INPUT_DATA_ERROR);
    } else if (bad_chromo || mrk_missing_from_map || mrk_missing_from_names) {
        if (InputMode == INTERACTIVE_INPUTMODE || Mega2BatchItems[/* 28 */ Default_Ignore_Nonfatal].items_read == 0) {
            printf("Do you wish to continue with Mega2? (y/n) > ");
            fflush(stdin); fcmap(stdin, "%s", yesorno); newline;
        } else {
            yesorno[0] = Mega2BatchItems[/* 28 */ Default_Ignore_Nonfatal].value.copt;
        }
        if (yesorno[0] == 'n' || yesorno[0] == 'N') { EXIT(EARLY_TERMINATION); }
    } else {
#ifndef HIDEFILE
        sprintf(err_msg, "Done reading map file: %s\n", map_file);
        mssgf(err_msg);
#endif
    }
    
    /*   copy_exmap_locmap(LTop, EXLTop, 0); */

    /* Add a one in case there are any unmapped markers with an Unknown chromosome */
    global_chromo_entries=CALLOC((size_t)MaxChromo + 1, int);
    chromo_loci_count = CALLOC((size_t)MaxChromo + 1, int);

    // 'NumUnmapped' is the count of markers in the map file where the chromosome
    // is specified as unknown. 'mrk_missing_from_map' is the count of markers specified
    // elsewhere, but with no entry in the map file. 
    if (NumUnmapped > 0) {
        unmapped_markers = CALLOC((size_t)NumUnmapped, int);
    }
    
    return EXLTop;
}   /* end of read_common_map_file() */


static ext_linkage_locus_top *read_annotated_map_file(const char *map_file,
                                                      linkage_locus_top *LTop,
                                                      std::vector<m2_map> additional_maps,
                                                      annotated_file_desc *file_desc)
{
    listhandle *userdef_colnames = new_list();
    col_hdr_type reserved_colnames[NUM_MAPCOL_NAMES];
    col_hdr_type *colname_item, *map_all_colnames;
    int num_userdef_cols;
    //int chr;

    if (map_file == 0 || *map_file == 0) return (ext_linkage_locus_top *)NULL; // for _win

    FILE *mapfp = fopen(map_file, "r");

    if (mapfp == NULL) return (ext_linkage_locus_top *)NULL;

#ifndef HIDEFILE
    printf("Reading map file %s ...", map_file);
#endif /*  HIDEFILE */

    INIT_COLNAME(reserved_colnames, 0, STRING_AN, "Chromosome");
    INIT_COLNAME(reserved_colnames, 1, STRING_AN, "Name");
    INIT_COLNAME(reserved_colnames, 2, FLOAT_AN, "Error");
    num_userdef_cols = parse_map_file_header(mapfp, reserved_colnames, userdef_colnames);

    /* Check if there is an Error column */
    file_desc->num_map_cols = num_userdef_cols;
    file_desc->map_file_columns =
        CALLOC((size_t)file_desc->num_map_cols, col_hdr_type);
    map_all_colnames = &(file_desc->map_file_columns[0]);

    /* reverse index as for pedigree file, so that the order of column names
       is in order of input file
    */
    /* chromosome */
    colname_item =  &(map_all_colnames[reserved_colnames[0].input_col]);
    copy_colname(&(reserved_colnames[0]), colname_item);
    colname_item->output_col=0;

    /* marker name */
    colname_item =  &(map_all_colnames[reserved_colnames[1].input_col]);
    copy_colname(&(reserved_colnames[1]), colname_item);
    colname_item->output_col=1;

    if (reserved_colnames[2].input_col > -1) {
        colname_item =  &(map_all_colnames[reserved_colnames[2].input_col]);
        copy_colname(&(reserved_colnames[2]), colname_item);
        colname_item->output_col=2;
    }

    /* userdefined columns i.e. Maps */
    while((colname_item =
           (col_hdr_type *) pop_first_list_entry(userdef_colnames))
          != NULL) {
        if (additional_maps.size()) {
            for (size_t i = 0; i < additional_maps.size(); i++) {
                m2_map map = additional_maps[i];
                const char *illegal_annotated_map_name = map.get_name().c_str();
                // Don't use any of the map names that are fixed given the type of input...
                if (strcasecmp(colname_item->ColName, illegal_annotated_map_name) == 0) {
                    errorvf("The map file %s may not contain a map named %s\n", map_file, illegal_annotated_map_name);
                    errorf("when a VCF file is specified as input.");
                    EXIT(DATA_INCONSISTENCY);
                }
            }
        }
        copy_colname(colname_item, &(map_all_colnames[colname_item->input_col]));
        free(colname_item);
    }
    free(userdef_colnames);

    return read_common_map_file(mapfp, map_file, LTop, (int)additional_maps.size(), reserved_colnames, file_desc, (char **)NULL);
}


/* file format =
   Marker-name allele-name pop-name frequency
   Flag error if frequency not a floating-point number
   Comment lines start with a #
*/

static int parse_frequency_file_header(FILE *filep,
                                       annotated_file_desc *file_desc)
{
    char buffer[READ_CHUNK + 1];
    char *col_name;
    int num_read=0, lch = 0;

    file_desc->freq_file_columns = CALLOC((size_t)NUM_FREQCOL_NAMES, col_hdr_type);
    INIT_COLNAME(file_desc->freq_file_columns,
                 0, STRING_AN, "Name");
    INIT_COLNAME(file_desc->freq_file_columns,
                 1, STRING_AN, "Allele");
/*   INIT_COLNAME(file_desc->freq_file_columns, */
/* 	       2, STRING_AN, "Group"); */
    INIT_COLNAME(file_desc->freq_file_columns,
                 2, FLOAT_AN, "Frequency");

    IgnoreValue(fgets(buffer, READ_CHUNK, filep));
    if (strlen(buffer) >= READ_CHUNK && buffer[strlen(buffer)-1] != '\n') {
        fcmap(filep, "%=\n", lch);
    }
    col_name = strtok(buffer, "\t \n ");

    while(num_read < NUM_FREQCOL_NAMES) {
        if (col_name == NULL) break;

        if (!strcasecmp(col_name, "Name")) {
            file_desc->freq_file_columns[0].input_col = num_read;
        } else if (!strcasecmp(col_name, "Allele")) {
            file_desc->freq_file_columns[1].input_col = num_read;
        }
/*     else if (!strcasecmp(col_name, "Group")) { */
/*       file_desc->freq_file_columns[2].input_col = num_read; */
/*     } */
        else if (!strcasecmp(col_name, "Frequency")) {
            file_desc->freq_file_columns[2].input_col = num_read;
        } else {
            error_hdr(col_name, NUM_FREQCOL_NAMES, file_desc->freq_file_columns);
            return 0;
        }
        num_read++;
        col_name = strtok(NULL, "\t \n");
    }
    file_desc->num_freq_cols = num_read;
    return num_read;
}

static void freq_insert_into_allele_list(allele_list_type **marker_item,
                                         char *allele_name,
                                         const double frequency_val,
                                         const char *marker_name)
{
    allele_list_type *marker_itemp = *marker_item;
    allele_list_type *marker_itemp1 = *marker_item;
    allele_list_type *new_entry;
    
    new_entry = CALLOC((size_t)1, allele_list_type);
    /* all_val is a random string in the file */
    new_entry->allele_freq.AlleleName = canonical_allele(allele_name);
    new_entry->allele_freq.freq = frequency_val;
    
    new_entry->allele_freq.count =
    new_entry->allele_freq.founder_count =
    new_entry->allele_freq.random_count =
    new_entry->allele_freq.unique_count =
    new_entry->allele_freq.everyone_count = 0;
    
    if (marker_itemp == NULL) {
        new_entry->next=NULL;
        marker_itemp = new_entry;
        marker_itemp->allele_freq.index = 1;
        *marker_item = marker_itemp;
    } else {
        while((marker_itemp != NULL) &&
              // POSSIBLE BUG: SHouldn't this be '== 0'???
              (strcmp(marker_itemp->allele_freq.AlleleName, allele_name) <= 0)) {
            marker_itemp1 = marker_itemp;
            marker_itemp = marker_itemp->next;
        }
        
        if (allelecmp(marker_itemp1->allele_freq.AlleleName,allele_name)==0) {
            /* warn the user */
            warnvf("Duplicate allele record for %s, marker %s\n", allele_name, marker_name);
            warnf("Ignoring duplicate entries.");
            free(new_entry);
        } else {
            if (marker_itemp == NULL) {
                /* Need to add a new at the end */
                new_entry->next=NULL;
                new_entry->allele_freq.index =
                marker_itemp1->allele_freq.index+1;
                marker_itemp1->next = new_entry;
            } else if (marker_itemp == marker_itemp1) {
                /* Add a new at the beginning */
                new_entry->next=marker_itemp1;
                new_entry->allele_freq.index = 1;
                *marker_item = new_entry;
            } else {
                /* Need to insert between itemp and item1 */
                new_entry->next = marker_itemp1->next;
                new_entry->allele_freq.index =
                marker_itemp1->allele_freq.index + 1;
                marker_itemp1->next = new_entry;
            }
            
            while(marker_itemp != NULL) {
                marker_itemp->allele_freq.index++;
                marker_itemp = marker_itemp->next;
            }
        }
    }
}

#ifndef num_alleles
#define num_alleles data.num_alleles
#endif

#ifndef num_classes
#define num_classes data.num_classes
#endif

/* Read the frequency information into the marker_list data structure used
   in recode. This function simply replaces the convert_to_freq function
   used within recode.
   There is extra code in here for defining populations or groups that has
   been commented out for now. Ext_locus_top is also not used. Number of
   groups is set to 1.
*/

//
// Read and process an annotated frequency file.
//
// Some errors are fatel, and with them the program will exit.
//
// returns:
//  1  if the file was processes without errors,
//  2  if some markers were absent fom the file,
//  0  if no frequency file was found
static int read_annotated_freq_file(char *freq_file_name,
                                    linkage_locus_top *LTop,
                                    annotated_file_desc *file_desc,
                                    marker_type *marker_list,
                                    pheno_type *pheno_list,
                                    int *num_groups, int **groups)
{
    char frequency_line[READ_CHUNK+1];
    int j, mrk_index, num_read, line_num;
    FILE *filep;
    double freq;
    char strs[3][FILENAME_LENGTH];
    char *name, *allname;
    
    /*   int *allele_count; /\* Allelecounts per group for the marker *\/ */
    /*   int *line_indexes; /\* unique marker number for each line *\/ */
    /*   int *freq_to_ltop_map;  */
    /*   int has_grp=0, num_freq_cols; */
    /*   int grp, grp_err;  */
    /*   int num_pops = 0; */
    /*   ext_linkage_locus_top *EXLTop = new_ex_llocustop(); */
    
    *num_groups = 1;
    *groups = CALLOC((size_t)1, int);
    (*groups)[0] = 0;
    
    if (freq_file_name == NULL) {
        warnf("No frequency file provided.");
        warnf("Allele frequencies for these will be estimated from data.");
        return 0;
    }
    
#ifndef HIDEFILE
    sprintf(err_msg, "Reading Mega2 format frequency file: %s", freq_file_name);
    mssgf(err_msg);
#endif
    
    if ((filep=freq_file_name ? fopen(freq_file_name, "r") : NULL) == NULL) {
        errorvf("could not open %s for reading!\n", freq_file_name);
        EXIT(FILE_READ_ERROR);
    }
    
    // Make sure that the frequency file header has the appropriate columns...
    if ((file_desc->num_freq_cols =
         parse_frequency_file_header(filep, file_desc)) < NUM_FREQCOL_NAMES) {
        errorvf("Frequency file %s is not in correct format.\n", freq_file_name);
	EXIT(DATA_INCONSISTENCY);
    }
    line_num=1;
    
    while (!feof(filep)) {
        strcpy(frequency_line, "");
        IgnoreValue(fgets(frequency_line, READ_CHUNK, filep));
        if (*frequency_line == 0) break;
        line_num++;
        // don't process empty or comment lines...
        if (frequency_line[0] == '#' || frequency_line[0] == '\n') continue;
        
        // Name, Allele, Frequency (though they can be in any order)...
        num_read = sscanf(frequency_line, "%s %s %s", strs[0], strs[1], strs[2]);
        if (num_read < file_desc->num_freq_cols) {
            errorvf("Line %d: incorrect number of columns in freq file %s.\n",
		    line_num, freq_file_name);
            EXIT(DATA_INCONSISTENCY);
        }
        
        // The indirection is because the 'Name' can appear in any column and
        // parse_frequency_file_header sorts out what column it is actually in.
        name = strs[file_desc->freq_file_columns[0].input_col];
        int ret;
        ret = search_marker(name, &mrk_index);
        if (ret != 1) {
            warnvf("Line %d of frequency file %s references locus %s.\n",
                   line_num, freq_file_name, name);
            continue;
        }
        
        // This should be the 'Allele'...
        allname = strs[file_desc->freq_file_columns[1].input_col];
        
        /* group */
        /* 	if (num_freq_cols == 4) { */
        /* 	  grp=atoi(strs[file_desc->freq_file_columns[2].input_col]); */
        /* 	  /\* set the number of populations as we go *\/ */
        /* 	  num_pops = ((num_pops < grp)? grp : num_pops); */
        /* 	} */
        /* 	else { */
        /* 	  grp = 0; */
        /* 	} */
        
        // This should be the 'Frequency'...
        freq=atof(strs[file_desc->freq_file_columns[2].input_col]);
        
        if (mrk_index < LTop->PhenoCnt) {
            freq_insert_into_allele_list(&(pheno_list[mrk_index].first_allele),
                                         allname, freq, name);
            pheno_list[mrk_index].estimate_frequencies = 0;
        } else {
            freq_insert_into_allele_list(&(marker_list[mrk_index].first_allele),
                                         allname, freq, name);
            marker_list[mrk_index].estimate_frequencies = 0;
        }
    }
    fclose(filep);

    mrk_index=0;
    for (j=LTop->PhenoCnt; j < LTop->LocusCnt; j++) {
        if (marker_list[j].estimate_frequencies == 0) {
            mrk_index++;
        }
    }
    /* Possible error: if MarkerCnt is out of bounds for the 'int' type */
    if (mrk_index < (int) LTop->MarkerCnt) {
        /* estimate frequencies */
        warnf("Some markers are absent from frequency file.");
        warnf("Allele frequencies for these markers will be estimated from data.");
	return 2;
    }
    return 1;
}

static int parse_penetrance_file_header(FILE *filep, char *pen_file_name,
                                        annotated_file_desc *file_desc)
{

    char buffer[READ_CHUNK + 1];
    char *col_name;
    int num_read=0;

    file_desc->pen_file_columns = CALLOC((size_t)NUM_PENCOL_NAMES, col_hdr_type);
    INIT_COLNAME(file_desc->pen_file_columns,
                 0, STRING_AN, "Name");
    INIT_COLNAME(file_desc->pen_file_columns,
                 1, STRING_AN, "Class");
    INIT_COLNAME(file_desc->pen_file_columns,
                 2, STRING_AN, "Sex");
    INIT_COLNAME(file_desc->pen_file_columns,
                 3, FLOAT_AN, "Pen.11");
    INIT_COLNAME(file_desc->pen_file_columns,
                 4, FLOAT_AN, "Pen.12");
    INIT_COLNAME(file_desc->pen_file_columns,
                 5, FLOAT_AN, "Pen.22");

    IgnoreValue(fgets(buffer, READ_CHUNK, filep));
/* Not used:
 * read_length=strlen(buffer); */

/*   if (read_length >= READ_CHUNK &&  buffer[read_length-1] != '\n') { */
/*     fcmap(filep, "%=\n", lch); */
/*   } */

    col_name = strtok(buffer, "\t \n ");

    while(num_read < NUM_PENCOL_NAMES) {
        if (col_name == NULL) {
            break;
        }
        if (!strcasecmp(col_name, "Name")) {
            file_desc->pen_file_columns[0].input_col = num_read;
        } else if (!strcasecmp(col_name, "Class")) {
            file_desc->pen_file_columns[1].input_col = num_read;
        } else if (!strcasecmp(col_name, "Sex")) {
            file_desc->pen_file_columns[2].input_col = num_read;
            msgvf("Penetrance File \"%s\", \"Sex\" is not the right heading term.\nConsider using \"Type\" instead.\n",
                  pen_file_name);
        } else if (!strcasecmp(col_name, "Type")) {
            file_desc->pen_file_columns[2].input_col = num_read;
        } else if (!strcasecmp(col_name, "Pen.11")) {
            file_desc->pen_file_columns[3].input_col = num_read;
        } else if (!strcasecmp(col_name, "Pen.12")) {
            file_desc->pen_file_columns[4].input_col = num_read;
        } else if (!strcasecmp(col_name, "Pen.22")) {
            file_desc->pen_file_columns[5].input_col = num_read;
        } else {
            error_hdr(col_name, NUM_PENCOL_NAMES, file_desc->pen_file_columns);
            return 0;
        }
        num_read++;
        col_name = strtok(NULL, "\t \n");
    }
    file_desc->num_pen_cols = num_read;
    return num_read;
}

static class_list_type *new_class_itemp(int class_num)
{
    class_list_type *new_entry = CALLOC((size_t)1, class_list_type);

    new_entry->l_class.num =  class_num;
    (new_entry->l_class).autosomal_pen =
        (new_entry->l_class).male_pen = (new_entry->l_class).female_pen = NULL;

    new_entry->next = NULL;
    return new_entry;
}

static void fill_penetrances(char *trait, liability_class *lclass,
                             char *sex, double *pen)
{
    switch(sex[0]) {
    case 'a':
    case  'A':
        if (lclass->autosomal_pen != NULL) {
            sprintf(err_msg, "Ignoring duplicate penetrances for %s class %d",
                    trait, lclass->num);
            warnf(err_msg);
        } else {
            lclass->autosomal_pen = CALLOC((size_t)3, double);
            lclass->autosomal_pen[0] = pen[0];
            lclass->autosomal_pen[1] = pen[1];
            lclass->autosomal_pen[2] = pen[2];
        }
    break;

    case 'm':
    case 'M' :
        if (lclass->male_pen != NULL) {
            sprintf(err_msg, "Ignoring duplicate male penetrances for %s class %d",
                    trait, lclass->num);
            warnf(err_msg);
        } else {
            lclass->male_pen = CALLOC((size_t)2, double);
            lclass->male_pen[0] = pen[0];
            lclass->male_pen[1] = pen[1];
        }
    break;
    case 'f':
    case 'F':
        if (lclass->female_pen != NULL) {
            sprintf(err_msg, "Ignoring duplicate female penetrances for %s class %d",
                    trait, lclass->num);
            warnf(err_msg);
        } else {
            lclass->female_pen = CALLOC((size_t)3, double);
            lclass->female_pen[0] = pen[0];
            lclass->female_pen[1] = pen[1];
            lclass->female_pen[2] = pen[2];
        }
    break;
    }
}


/* This is similar to reading in the frequency file,
   we discover new liability classes as we read the penetrance
   file. Assume that trait_list has already been allocated after
   reading in the names file.
   Returns 0 if succesful in inerting
   1 if this class_name is NA and previously encountered class-name
   not NA
   2 class_names match but dupliacte penetrances for the same phenotype
*/
static int pen_insert_into_class_list(char *trait, class_list_type **class_list,
                                      char *class_name, int class_num, char *sex,
                                      double *pen)
{
    class_list_type *class_itemp = *class_list;
    class_list_type *class_itemp1 = *class_list;
    class_list_type *new_entry;

    if (class_num == 0) {
        /* this trait has only one class */
        if (class_itemp == NULL) {
            /* insert the penetrances */
            class_itemp = new_class_itemp(class_num);
            *class_list = class_itemp;
        }
        fill_penetrances(trait, &((*class_list)->l_class), sex, pen);
    } else {
        /* find the spot to insert in */
        if (class_itemp == NULL) {
            /* list is empty */
            class_itemp = new_class_itemp(class_num);
            *class_list = class_itemp;
            fill_penetrances(trait, &(class_itemp->l_class), sex, pen);
        } else {
            while((class_itemp != NULL) && (class_itemp->l_class.num < class_num) )  {
                class_itemp1 = class_itemp;
                class_itemp = class_itemp->next;
            }

            if (class_itemp != NULL && (class_itemp->l_class.num == class_num) ) {
                fill_penetrances(trait, &(class_itemp->l_class), sex, pen);
                return 0;
            }
            /* else a new item has to be created */
            new_entry = new_class_itemp(class_num);
            fill_penetrances(trait, &(new_entry->l_class), sex, pen);

            if (class_itemp == NULL) {
                /* Need to add the new at the end */
                class_itemp1->next = new_entry;
            } else {
                /* itemp->name > name */
                if (class_itemp == class_itemp1) {
                    /* insert at the beginning of the list */
                    new_entry->next=class_itemp1;
                    *class_list=new_entry;
//                  class_itemp=new_entry;
                } else {
                    /* need to insert new between itemp and itemp1 */
                    new_entry->next = class_itemp1->next;
                    class_itemp1->next = new_entry;
                }
            }
        }
    }

    return 0;
}


static int read_annotated_pen_file(char *pen_file_name,
                                   linkage_locus_top *LTop,
                                   annotated_file_desc *file_desc,
                                   pheno_type *pheno_list)
{

    char penetrance_line[READ_CHUNK+1];
    int num_read, line_num;
    int mrk_index;
    FILE *filep;
    double pen[3];
    char strs[6][70];
    char sex[5], *name, *class_name;
    int class_num = -1;
    int num_err=0;

#ifndef HIDEFILE
    sprintf(err_msg, "Reading Mega2 format penetrance file: %s", pen_file_name);
    mssgf(err_msg);
#endif

    if ((filep=pen_file_name ? fopen(pen_file_name, "r") : NULL) == NULL) {
        errorvf("could not open %s for reading!\n", pen_file_name);
        EXIT(FILE_READ_ERROR);
    }

    if (parse_penetrance_file_header(filep, pen_file_name, file_desc) < 5) {
        warnvf("Input file %s not in correct format.\n", pen_file_name);
        fclose(filep);
        return -1;
    }
    line_num=1;
    while(!feof(filep)) {
        strcpy(penetrance_line, "");
        IgnoreValue(fgets(penetrance_line, READ_CHUNK, filep));
        if (*penetrance_line == 0) break;
        line_num++;
        if (penetrance_line[0] == '#') {
            /* comment line */
/*          skip_line(filep, ch);  NO! fgets read the line! */
        } else if (penetrance_line[0] == '\n') {
        } else {
            if (file_desc->pen_file_columns[2].input_col == -1) {
                num_read = sscanf(penetrance_line, "%s %s %s %s %s",
                                  strs[0], strs[1], strs[2],
                                  strs[3], strs[4]);
            } else {
                num_read = sscanf(penetrance_line, "%s %s %s %s %s %s",
                                  strs[0], strs[1], strs[2],
                                  strs[3], strs[4], strs[5]);
            }
            if (num_read < file_desc->num_pen_cols) {
                if (num_read < 1) {
                    continue;
                } else {
                    num_err++;
                    errorvf("Line %d: too few columns in penetrance file.\n", line_num);
                    continue;
                }
            } else { /* required number of columns */
                /* name */
                name = strs[file_desc->pen_file_columns[0].input_col];
                int ret;
                ret = search_marker(name, &mrk_index);
                if (ret != 1) {
                    errorvf("Line %d: trait %s in penetrance file not in the names file.\n",
                            line_num, name);
                    num_err++;
                    continue;
                }
                /* class name */
                class_name = strs[file_desc->pen_file_columns[1].input_col];
                /* sex */
                if (file_desc->pen_file_columns[2].input_col == -1) {
                    strcpy(sex, "auto");
                } else if (!strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "autosomal") || 
                           !strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "autosome")) {
                    strcpy(sex, "auto");
                } else if (!strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "average")) {
                    strcpy(sex, "auto");
                    msgvf("Penetrance File \"%s\", \"average\" is not the right term for penetrance of the trait %s on line %d.\nConsider using \"autosomal\" instead.\n",
                            pen_file_name, name, line_num);
                } else if (!strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "male")) {
                    strcpy(sex, "male");
                } else if (!strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "female")) {
                    strcpy(sex, "fem");
                } else {
                    errorvf("Penetrance File \"%s\", Invalid sex for trait %s on line %d\n",
                            pen_file_name, name, line_num);
                    num_err++;
                    continue;
                }

                /* penetrances */
                pen[0]=atof(strs[file_desc->pen_file_columns[3].input_col]);
                if (!strcasecmp(strs[file_desc->pen_file_columns[2].input_col], "male")) {
                    pen[1]=atof(strs[file_desc->pen_file_columns[5].input_col]);
                    pen[2]=0;
                } else {
                    pen[1]=atof(strs[file_desc->pen_file_columns[4].input_col]);
                    pen[2]=atof(strs[file_desc->pen_file_columns[5].input_col]);
                }

                if (pen[0] < 0.0 || pen[0] > 1.0 ||
                    pen[1] < 0.0 || pen[1] > 1.0 ||
                    pen[2] < 0.0 || pen[2] > 1.0) {
                    errorvf("Penetrance File \"%s\", Invalid penetrance values for trait %s on line %d\n",
                            pen_file_name, name, line_num);
                    num_err++;
                    continue;
                }

                char *endp;
                class_num = (int)strtol(class_name, &endp, 10);
                if (class_num <= 0 || *endp) {
                    num_err++;
                    errorvf("Penetrance File \"%s\", Line %d: Trait %s should have an integer class label not \"%s\".\n", 
                            pen_file_name, line_num, name, class_name);
                } else if (pen_insert_into_class_list(name,
                                                      &(pheno_list[mrk_index].class_list),
                                                      class_name, class_num, sex, &(pen[0])) == 1) {
                    num_err++;
                    errorvf("Penetrance File \"%s\", Line %d: Trait %s previously encountered with known class label.\n", 
                            pen_file_name, line_num, name);
                    errorf("Traits with unknown class (NA) can have only one class.");
                }
            }
        }
    }

    fclose(filep);
    return num_err;
}

/*--------------------end of penetrance file routines----------------------*/


static void init_reserved_pedcol_names(void)
{
    /* Initialize keyword names and types,
       then also fill in the column values with -1 to indicate absent
    */

    /* The first 5 are required headers */
//??
//    INIT_COLNAME(ReservedColnames, 1, STRING_AN, "Pedigree");
//    INIT_COLNAME(ReservedColnames, 0, STRING_AN, "ID");

    INIT_COLNAME(ReservedColnames, 0, STRING_AN, "Pedigree");
    INIT_COLNAME(ReservedColnames, 1, STRING_AN, "ID");
    INIT_COLNAME(ReservedColnames, 2, STRING_AN, "Father");
    INIT_COLNAME(ReservedColnames, 3, STRING_AN, "Mother");
    INIT_COLNAME(ReservedColnames, 4, CHAR_AN,   "Sex");
    /* The rest are optional headers */
    INIT_COLNAME(ReservedColnames, 5, STRING_AN, "PerID");
    INIT_COLNAME(ReservedColnames, 6, STRING_AN, "PedID");
    INIT_COLNAME(ReservedColnames, 7, INT_AN, "MZTwin");
    INIT_COLNAME(ReservedColnames, 8, INT_AN, "DZTwin");
    INIT_COLNAME(ReservedColnames, 9, INT_AN, "Proband");
    INIT_COLNAME(ReservedColnames, 10, INT_AN, "Group");
    INIT_COLNAME(ReservedColnames, 11, STRING_AN, "FirstOff");
    INIT_COLNAME(ReservedColnames, 12, STRING_AN, "NextMatSib");
    INIT_COLNAME(ReservedColnames, 13, STRING_AN, "NextPatSib");
    INIT_COLNAME(ReservedColnames, 14, STRING_AN, "LinkPerID");
    INIT_COLNAME(ReservedColnames, 15, STRING_AN, "LinkPedID");
    return;
}

static int is_missing_required_hdr(col_hdr_type *col_names, int num_elements)
{
    int i, error=0;
    for (i=0; i < num_elements; i++) {
        if (col_names[i].input_col == -1) {
            errorvf("Missing required column %s\n", col_names[i].ColName);
            error=1;
        }
    }
    return error;
}

static void error_hdr(char *hdr,
                      int num_known,
                      col_hdr_type *known_hdrs)
{
    int i;
    errorvf("Unknown header %s\n", hdr);
    errorf("Allowed headers are:");
    sprintf(err_msg, "%s", known_hdrs[0].ColName);
    for (i=1; i < num_known; i++) {
        grow(err_msg, ", %s", known_hdrs[i].ColName);
    }
    errorf(err_msg);
}

//
// Because a header name may now have the '.' character in it, it is nolonger possible
// to identify an extension in a colname (header name) by looking form the beginning.
// Some extensions have one field (e.g., '.P') and some have two (e.g., '.H.A').
// So we must compare the right handed portion of a name with a list of valid extensions.
// Return -1 if no valid extension was found, or an index into the colname (header_name)
// where the first valid header extension was found.
// NOTE: Remember, the first valid extension is found so ordering is important in certain cases
static int is_valid_hdr_ext(const char *colname,
                            const char **valid_header_item_extensions)
{
    const char **p = valid_header_item_extensions;
    size_t len_colname = strlen(colname);

    // The valid_header_item_extension list is assumed to be NULL terminated...
    while (*p != NULL) {
        // this is a pointer to one of the valid extensions...
        const char *ext = *p++;
        size_t len_ext = strlen(ext);
        size_t ext_index_in_colname = len_colname - len_ext;
        // Case should not matter on the comparison...
        if (ext_index_in_colname > 0 && strcasecmp(ext, &colname[ext_index_in_colname]) == 0)
            return (int)ext_index_in_colname;
    }
    
    return -1; // not found
}

// Create an error message stating that the header 'colname' was found not to have
// a valid extension.
static void error_hdr_ext(const char *colname,
                          const int column_index,
                          const char **valid_header_item_extensions)
{
    char mssg[1000];
    const char **p = valid_header_item_extensions;
    errorvf("The suffix of the column '%s' in column #%d is invalid!\n", colname, column_index);
    strcpy(mssg, "Valid suffixes are:");
    while (*p != NULL) {
        const char *ext = *p++;
        sprintf(&mssg[strlen(mssg)], ((*p != NULL) ? " %s," :  " %s"), ext);
    }
    errorf(mssg);
}

//
// Search for the extensions from right to left, so that ext2 is the string
// after the last period in colname, and ext1 is the string between the last
// and the second to last period in colname (see comments for 'is_valid_hdr_ext()').
// If only one prefix is found, then ext1 will be filled and ext2 set to the empty
// string. If no prefix is found then ext1 and ext1 will be set to the empty string.
static void parse_hdr(const char *colname,
                      char *name, char *ext1, char *ext2,
                      const char **valid_header_item_extensions)
{
    long i, name_i, ext_i, ext1_len, colname_len = (long)strlen(colname);
    
    strcpy(name, "");
    *ext1=*ext2='\0';
    
    // did we find a valid extension in this header???
    if ((name_i = is_valid_hdr_ext(colname, valid_header_item_extensions)) <= 0) return;
    
    // name_i is the index into colname where the extension begins...
    memcpy(name, colname, (size_t)name_i);
    name[name_i] = '\0';
    
    for (ext_i = -1, i = name_i + 1; i < colname_len; i++) {
        // look for the next '.' in the string...
        if (colname[i] == '.') {
            ext_i = i;
            break;
        }
    }
    if (ext_i == -1) {
        // we only found one of the two extensions...
        strcpy(ext1, &colname[name_i+1]);
        return;
    }
    
    // we found two extensions....
    ext1_len = colname_len - ext_i - 1;
    memcpy(ext1, &colname[name_i+1], ext1_len);
    ext1[ext1_len] = '\0';
    
    strcpy(ext2, &colname[ext_i+1]);
}

/**
 @brief Omit file data processing support for annotated files.
 
 The annotated file preprocessing code needs to read the file header to
 determine the order of <ped, per, loci> before calling the data processing routine.
 @see omit_file_data_processing
 
 @param[in,out] Top         The pedigree information @see http://watson.hgen.pitt.edu/doxygen/struct__linkage__ped__top.html
 @param[in]     omitfl_name Is the name of the omit file to read from.
 @param[in,out] file_desc   Structure that describes column semantics.
 @return void
 */
static void annotated_omit_file(linkage_ped_top *Top,
                                const char *omitfl_name,
                                annotated_file_desc *file_desc)
{
    char buffer[READ_CHUNK + 1], *token, col_name[MAX_NAMELEN];
    int num_read=0, lch=0;
    FILE *omitfp;
    
#ifndef HIDEFILE
    msgvf("Reading Mega2 format omit file: %s\n", omitfl_name);
#endif
    if ((omitfp=omitfl_name ? fopen(omitfl_name, "r") : NULL) == NULL) {
        errorvf("could not open %s for reading!\n", omitfl_name);
        EXIT(FILE_READ_ERROR);
    }
    
    /* read the omit file header */
    file_desc->num_omit_cols = NUM_OMITCOL_NAMES;
    file_desc->omit_file_columns = CALLOC((size_t)NUM_OMITCOL_NAMES, col_hdr_type);
    INIT_COLNAME(file_desc->omit_file_columns,
                 0, STRING_AN, "Pedigree");
    INIT_COLNAME(file_desc->omit_file_columns,
                 1, STRING_AN, "Individual");
    INIT_COLNAME(file_desc->omit_file_columns,
                 2, STRING_AN, "Marker");
    
    IgnoreValue(fgets(buffer, READ_CHUNK, omitfp));
    if (buffer[strlen(buffer)-1] != '\n') {
        fcmap(omitfp, "%=\n", lch);
    }
    token = strtok(buffer, "\t \n");
    
    // Search for the required header items...
    while (num_read < 3) {
        if (token != NULL) strcpy(col_name, token);
        if (!strcasecmp(col_name, "Pedigree")) {
            file_desc->omit_file_columns[0].input_col = num_read;
        } else if (!strcasecmp(col_name, "Individual")) {
            file_desc->omit_file_columns[1].input_col = num_read;
        } else if (!strcasecmp(col_name, "Marker")) {
            file_desc->omit_file_columns[2].input_col = num_read;
        } else {
            error_hdr(col_name, NUM_OMITCOL_NAMES, file_desc->omit_file_columns);
        }
        num_read++;
        token = strtok(NULL, "\t \n");
    }
    
    if (is_missing_required_hdr(file_desc->omit_file_columns, 3)) {
        errorvf("In omit file (%s). Three white space separated\n", omitfl_name);
        errorf("header items are required: Pedigree, Individual, and Marker.");
        EXIT(FILE_HEADER_ERROR);
    }
    
    // Process the data in the file...
    Display_Messages = Display_Errors = 1;
    omit_file_data_processing(Top, omitfl_name, omitfp,
                              file_desc->omit_file_columns[0].input_col,
                              file_desc->omit_file_columns[1].input_col,
                              file_desc->omit_file_columns[2].input_col,
                              1);
    Display_Messages = Display_Errors = 1;
    
    fclose(omitfp);
}

/* static void copy_group0_freq(linkage_locus_top *LocusTop,  */
/* 			     ext_linkage_locus_top *EXLTop) */

/* { */
/*   int loc, allele; */
/*   int# group = 0; */

/*   for (loc=0; loc < LocusTop->LocusCnt; loc++) { */
/*   } */
/*   return; */
/* } */

#ifdef _WIN
#define snprintf _snprintf
#endif
#ifndef HIDESTATUS
static void output_stats() {
  char msg[1000];
  // Markers read from the map file...
  snprintf(msg, 1000, "Reading pedigree information from %s", mega2_input_files[0]);
  mssgf(msg);
  snprintf(msg, 1000, "%d individuals read from %s", PLINK.individuals, mega2_input_files[0]);
  mssgf(msg);
  snprintf(msg, 1000, "%d individuals with nonmissing phenotypes", PLINK.phenotyped_individuals);
  mssgf(msg);
  snprintf(msg, 1000, "%d cases, %d controls, %d missing", PLINK.cases, PLINK.controls, PLINK.missing);
  mssgf(msg);
  snprintf(msg, 1000, "%d males, %d females, %d of unspecified sex", PLINK.males, PLINK.females, PLINK.unspecified_sex);
  mssgf(msg);
  snprintf(msg, 1000, "%d founders, %d non-founders found", PLINK.founders, PLINK.non_founders);
  mssgf(msg);
}
#endif /* HIDESTATUS */

//
// If a Mega2 map file was specified, the routine 'read_common_map_file()' will create an extra
// '.posiitons' entry (at the end) that it ignores when processing a VCF (/PLINK) file.
// If no Meta2 map file was specified, then a '.positions' in EXLTop needs to be created with
// just one position entry when map_i == 0. 'map_i' points to the zero based map entry.
// This routine will insert a SEX_AVERAGED_MAP or physical map into the map slot 'map_i'.
static void insert_m2_map_into_EXLTop(ext_linkage_locus_top *EXLTop,
                                      linkage_locus_top *LTop,
                                      const int map_i,
                                      m2_map& map)
{
    int i, j;
    
    EXLTop->map_functions[map_i] = map.get_function();
    EXLTop->MapNames[map_i] = strdup(map.get_name().c_str());
    EXLTop->MapNames[map_i][map.get_name().size()-2] = 0;
    EXLTop->SexMaps[map_i] = CALLOC((size_t) 3, int);
    // For a physical map, the data is stored as a SEX_AVERAGED_MAP...
    EXLTop->SexMaps[map_i][SEX_AVERAGED_MAP] = 1;
    // read_common_map_file() only accounts for it's maps, so we need to increase the count here...
    EXLTop->MapCnt++; // see annotated_ped_file.cpp::get_map_names()
    
    if (map_i == 0) MaxChromo = NumUnmapped = 0;
    
    // For each locus that is a MARKER fill it with 'map' data...
    // There is other book keeping as well that can only be done if the user actually selects this map.
    // See: 'create_entries_for_markers_without_positions()'
    for (i=0, j=0; i < LTop->LocusCnt; i++) {
        // We can process entries in order because 'read_common_marker_data()' loads the markers
        // from 'names[]' which would be in the order that they are read from the VCF file...
        // .Class is set when the Names file 'type' is processed in 'read_common_marker_data()'
        if (LTop->Locus[i].Class == MARKER) {
            // Make the map entry for one map (this one) in EXLTop if there are no others...
            if (map_i == 0) new_EXLTop_positions(EXLTop, i, 1);
            
            // Since the order of entries in the 'names[]' vector should match that of the id[] vector,
            // we can just insert the map positions in order (of 'j')...
            //
            // The chromosome is not known when the names file if read, so it is filled in when processing
            // the map file. It could still be 'UNKNOWN_CHROMO', and we don't give an error if it is...
            // We also don't keep track of the count of male, female, XY, and MT as read_common_map_file() does...
            m2_map_entry map_entry = map.get_entry(j++);
            
            double position = map_entry.get_POS();
            // Must set .number > 0 phenotype array
            EXLTop->EXLocus[i].positions[map_i] = position;
            // Validate the map entry for this one map in EXLTop because there are no others...
            if (map_i == 0) {
                // We can set .chromosome and .number now because there is only one map to select.
                int chr = map_entry.get_chr();
                LTop->Marker[i].chromosome = chr;
                MaxChromo = ((chr != UNKNOWN_CHROMO && chr > MaxChromo)? chr : MaxChromo);
                if (chr == UNKNOWN_CHROMO) NumUnmapped++;
                
                //LTop->Marker[mrk_num].error_prob // Since Marker[] is CALLOCed, this is '0'...
                // .number is initialized to -1 in 'read_common_marker_data()'; should run from i - n...
                LTop->Locus[i].number = j;
                //LTop->SexLinked = 0;
            }
        }
    }
    
    // Verify an assumption...
    if (j != (int)map.size()) {
        errorvf("INTERNAL Names and Map entries dissagree in count.\n");
        EXIT(SYSTEM_ERROR);
    }
    
    // This is the special case where there was no EXLTop created by 'read_annotated_map_file()'
    // and 'read_common_map_file()'.
    if (map_i == 0) {
//      EXLTop->LocusCnt = (int)map.size();
        EXLTop->LocusCnt = LTop->LocusCnt;
        /* Add a one in case there are any unmapped markers with an Unknown chromosome */
        global_chromo_entries=CALLOC((size_t)MaxChromo + 1, int);
        chromo_loci_count = CALLOC((size_t)MaxChromo + 1, int);
        // 'NumUnmapped' is the count of markers in the map file where the chromosome
        // is specified as unknown. 'mrk_missing_from_map' is the count of markers specified
        // elsewhere, but with no entry in the map file.
        if (NumUnmapped > 0) {
            unmapped_markers = CALLOC((size_t)NumUnmapped, int);
        }
    }
}

//
// Look for a valid genetic map in EXLTop, and return a boolian result
static int valid_genetic_map_exists_in_EXLTop(ext_linkage_locus_top *EXLTop)
{
    for (int i=0; i<EXLTop->MapCnt; i++)
        if ((EXLTop->map_functions[i] == 'h' || EXLTop->map_functions[i] == 'k') &&
            (EXLTop->SexMaps[i][SEX_AVERAGED_MAP] != 0 || EXLTop->SexMaps[i][FEMALE_SEX_MAP] != 0))
            return 1;
    return 0;
}

//
// This is a temporary fix to get around the fact that all of the analysis currently
// need a genetic map. Here we create a sex averaged genetic map.
static void insert_zero_sex_average_genetic_map_in_EXLTop(ext_linkage_locus_top *EXLTop,
                                                          linkage_locus_top *LTop)
{
    // Create an additional map 'at the end'...
    const int map_i = EXLTop->MapCnt++;
    // Allocate an additional slot in EXLTop for the map...
    // NOTE: When extending a region allocated with calloc(3), realloc(3)
    // does not guarantee that the additional memory is also zero-filled.
    EXLTop->map_functions = (char *)REALLOC(EXLTop->map_functions, EXLTop->MapCnt+1, char);
    EXLTop->MapNames = (char **)REALLOC(EXLTop->MapNames, EXLTop->MapCnt, char *);
    EXLTop->SexMaps = (int **)REALLOC(EXLTop->SexMaps, EXLTop->MapCnt, int *);
    
    EXLTop->map_functions[map_i] = 'h';
    EXLTop->MapNames[map_i] = strdup("DummyMap.h.a");
    EXLTop->SexMaps[map_i] = (int *)CALLOC((size_t) 3, int);
    EXLTop->SexMaps[map_i][SEX_AVERAGED_MAP] = 1;
    
    for (int i=0; i < LTop->LocusCnt; i++)
        if (LTop->Locus[i].Class == MARKER) {
            // [re]Allocate space for an additional .pos* slot for the map structure...
            EXLTop->EXLocus[i].positions = (double *)REALLOC(EXLTop->EXLocus[i].positions,
                                                             EXLTop->MapCnt,
                                                             double);
            EXLTop->EXLocus[i].pos_male = (double *)REALLOC(EXLTop->EXLocus[i].pos_male,
                                                            EXLTop->MapCnt,
                                                            double);
            EXLTop->EXLocus[i].pos_female = (double *)REALLOC(EXLTop->EXLocus[i].pos_female,
                                                              EXLTop->MapCnt,
                                                              double);
            // Fill the slot that we just made room for...
            EXLTop->EXLocus[i].positions[map_i] = 0; // sex averaged data goes here
            EXLTop->EXLocus[i].pos_female[map_i] = UNKNOWN_POSITION;
            EXLTop->EXLocus[i].pos_male[map_i] = UNKNOWN_POSITION;
        }
}

m2_map save_vcf_map;

SECTION_ERR_INIT(y_female);
linkage_ped_top *read_annotated_files(char *ped_file, char *names_file,
                                      char *map_file, char *freq_file,
                                      char *pen_file, char *omit_file,
                                      char *bed_file, char *phe_file,
                                      int untyped_ped_opt,
                                      analysis_type analysis,
                                      plink_info_type *plink_info)
{
    Tod tod_anfiles("read annotated files");
    int num_groups, count_option, count_halftyped=0;
    int *groups, j;
    int ann_files, i;
    linkage_ped_top *Top;
    linkage_locus_top *LTop;
    ext_linkage_locus_top *EXLTop;
    marker_type *marker_list;
    pheno_type *pheno_list;
    class_list_type *lclass;
    int   phe_cols = 0;
    int   bp_map   = 0;
    int   duplicate_mrk=0;
    const char *names_fn = (PLINK.plink == binary_PED_format) ? mega2_input_files[6] : mega2_input_files[0];

    m2_map vcf_map;
    m2_map impute_map;

    int  xcf = Input_Format == in_format_binary_VCF ||
               Input_Format == in_format_compressed_VCF ||
               Input_Format == in_format_VCF;

    init_tokens();

    Input->GetOps()->do_init(Input);

    if (Input->GetOps()->use_getops()) {
        LTop = Input->GetOps()->do_names(names_fn);
	ann_files = 1;
    } else if (PLINK.plink || xcf) {
        char **phe_names = NULL;
        int *phe_types   = NULL;
        int tot_cols = phe_cols = parse_phe_types(phe_file, &phe_names, &phe_types);

        if (PLINK.no_pheno == 0) {  /* there is an extra row ... if no_pheno == 0 */
            phe_names[tot_cols] = CALLOC(strlen(PLINK.trait)+1, char);
            strcpy(phe_names[tot_cols], PLINK.trait);
            phe_types[tot_cols] = PLINK.traitType;
            tot_cols++;
        }
        if (xcf) {
            std::string alternative_key = std::string(Mega2BatchItems[/* 57 */ VCF_Marker_Alternative_INFO_Key].value.name);

            // This is a physical map that is 'attached to the side' of the VCF file...
            // It will later be coppied to EXLTop, and the map object will be deleted by C++
            // when it goes out of scope.
            Tod vcfgm("VCF get map");
            vcf_map = VCFtools_get_map(alternative_key, "chr");
            vcfgm();

            // Store the map so that we can pull out the VCF reference alleles and drop then into a
            // file to be read by PLINK using '--reference-allele fn'. It is unclear at this point
            // whether the user will choose this map file or not.
            save_vcf_map = vcf_map;
            
            // Process the genotype (from VCF file 'vcf_map') and phenotype 'phe_*' marker data,
            // loading it into LTop...
            Tod vcfmn("VCF map as names");
            read_m2_map_as_names_file(vcf_map, &LTop, tot_cols, phe_names, phe_types);
            ann_files = 1;
            vcfmn();
        } else {
            Tod tod_pmap("read plink map as names");
            ann_files = read_plink_map_as_names_file(names_file ? names_file : map_file,
                                                     &LTop, tot_cols, phe_names, phe_types,
                                                     &AnnotatedFileInfo);
            names_fn = names_file ? names_file : map_file;
            tod_pmap();
        }
        for (i = 0; i < tot_cols; i++) {
            free(phe_names[i]);
        }
        free(phe_names);
        free(phe_types);
    } else {
        ann_files = read_annotated_names_file(names_file, &LTop, &AnnotatedFileInfo);
        names_fn = names_file;
    }
    free(AnnotatedFileInfo.names_file_columns);

    // CPK: Build the MARKER set and in doing so search for duplicate markers....
    clear_marker();
    Display_Errors = 1;
    Tod tod_init("hash markers names");
    SECTION_ERR_INIT(duplicate_markers);
    for (i=0; i < LTop->LocusCnt; i++) {
        char *name = LTop->Locus[i].LocusName;
        int i_old = test_and_add_marker(name, i);
        if (i_old >= 0) {
            SECTION_ERR(duplicate_markers);
            errorvf("Duplicate marker name '%s' in '%s'\n", name, names_fn);
            duplicate_mrk++;
        }
    }
    SECTION_ERR_FINI(duplicate_markers);
    tod_init();

    if (duplicate_mrk) {
        printf("Duplicate locus names, see %s for details.\n",  Mega2Err);
        errorf("Found fatal errors in locus file.");
        errorf("Please fix errors and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

    if (!ann_files) {
        return NULL;
    }

    if (LTop->MarkerCnt > 0) {
        Tod tod_mlist("init marker_list");
        marker_list = (CALLOC((size_t) LTop->MarkerCnt, marker_type)) - LTop->PhenoCnt;
#ifdef SHOWSTATUS
        msgvf("ALLOC SPACE: marker_list: %d MB (%d x %d)\n",
              LTop->MarkerCnt * sizeof(marker_type) / 1024 / 1024,
              LTop->MarkerCnt,  sizeof(marker_type));
#endif
        for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
            marker_list[i].first_allele = NULL;
            /* set estimate frequencies by default */
            marker_list[i].estimate_frequencies = 1;
            marker_list[i].recode_alleles = 0;
            marker_list[i].num_alleles = 0;
        }
        tod_mlist();
    } else
        marker_list = 0;

    if (LTop->PhenoCnt > 0) {
        Tod tod_plist("init pheno_list");
        pheno_list=CALLOC((size_t) LTop->PhenoCnt, pheno_type);
        for (i=0; i < LTop->PhenoCnt; i++) {
            pheno_list[i].first_allele = NULL;
            /* set estimate frequencies by default */
            pheno_list[i].estimate_frequencies = 1;
            pheno_list[i].num_alleles = 0;
            if (LTop->Locus[i].Type == AFFECTION) {
                pheno_list[i].num_classes = (int) LTop->Pheno[i].Props.Affection.ClassCnt;
            } else if (LTop->Locus[i].Type == QUANT) {
                pheno_list[i].num_classes = (int) LTop->Pheno[i].Props.Quant.ClassCnt;
            }
        }
        tod_plist();
    } else
        pheno_list = 0;

/* (rvb) 11/10/11 
    Since read_X_map_file() gets NULL (below), it returns a fresh EXLTop
    line below is unnecessary and leaks memory

    EXLTop = new_ex_llocustop();
*/
    std::vector<m2_map> additional_maps;
    if (Input->GetOps()->use_getops()) {
        Input->GetOps()->do_map(additional_maps);
        EXLTop = NULL; // but additional_maps.size() > 0 so see below; this makes compiler happy
    } else if (PLINK.plink) {
        // if this is PLINK format (double negative)
        // CPK: If the .map file is a .bim file we gather the alleles into the
        // vector that is described earlier in this routine...
        Tod tod_pmf("read plink map file");
        EXLTop = read_plink_map_file(names_file ? names_file : map_file,
                                     LTop,
                                     &AnnotatedFileInfo,
                                     plink_info);
        tod_pmf();
    } else if (xcf) {
        // Create an extra map slot in EXLTop for the map from the VCF file...
        additional_maps.push_back(vcf_map);
        EXLTop = NULL; // but additional_maps.size() > 0 so see below; this makes compiler happy
    } else {
        // only mega2 is left (and linkage ;-) No additional maps...
        EXLTop = read_annotated_map_file(map_file, LTop, additional_maps, &AnnotatedFileInfo);
        if (EXLTop == (ext_linkage_locus_top *)NULL) {
            // The map file was empty and there were no additional maps.
            // We error here now because this is not VCF or PLINK, so it's not OK to have no Mega2 map file...
            errorvf("could not open %s for reading!\n", map_file);
            EXIT(FILE_READ_ERROR);
        }
    }
    if (additional_maps.size() > 0) {
        EXLTop = read_annotated_map_file(map_file, LTop, additional_maps, &AnnotatedFileInfo);
        
        if (EXLTop == (ext_linkage_locus_top *)NULL) {
            // This will be the case if the map file was empty, or it contained no maps.
            // So we need to create the EXLTop data structure...
            EXLTop = new_EXLTop(LTop);
            // SEE: annotated_ped_file.cpp::get_map_names()
            // Continue to set up EXLTop for the '1' map in 'vcf_map'...
            EXLTop->map_functions = CALLOC((size_t)2, char); // e.g., 'h', 'k', or 'p'
            EXLTop->MapNames = CALLOC((size_t)1, char*);
            EXLTop->SexMaps = CALLOC((size_t)1, int *);
        }

        if (EXLTop->MapCnt == 0) {
            insert_m2_map_into_EXLTop(EXLTop, LTop, 0, additional_maps[0]);
            insert_zero_sex_average_genetic_map_in_EXLTop(EXLTop, LTop);
            // In this case, the batch file items must be set since the user does not specify these maps...

            //Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option = 0;
            base_pair_position_index = 0;
            Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].items_read=1;

            //Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option = 1;
            genetic_distance_index = 1;
            Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].items_read=1;

            //Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = SEX_AVERAGED_GDMT;
            genetic_distance_sex_type_map = SEX_AVERAGED_GDMT;
            Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].items_read=1;
        } else {
            // 'read_common_map_file()' will only initialize map positions (EXLTop .pos*)
            // if it finds an entry in the Mega2 map file. Because the VCF map file may contain a super set
            // of the marker from a mega2 map, position data for these markers will not have been created.
            // This is OK if any of the maps from the mega2 map file are used because the mega2 map
            // is used to 'subset' the markers found in the VCF file.
            for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++)
                if (EXLTop->EXLocus[i].positions == NULL)
                    new_EXLTop_positions(EXLTop, i, EXLTop->MapCnt+1);
            // When processing a vcf file, we have asked 'read_common_map_file()' above to create
            // an additional map slot, It will fill the first N-1 map slots with maps from the Mega2
            // annotated map file. We will fill the slot at the end (e.g., indexed by EXLTop->MapCnt)
            // which 'read_common_map_file()' creates but ignores.
            insert_m2_map_into_EXLTop(EXLTop, LTop, EXLTop->MapCnt, additional_maps[0]);
            // Since a genetic map is currently required for all analysis types (will fix this in the future)...
            if (!valid_genetic_map_exists_in_EXLTop(EXLTop)) {
                insert_zero_sex_average_genetic_map_in_EXLTop(EXLTop, LTop);

                // Since there was no genetic map in the Mega2 Map File, we added one at the end...
                genetic_distance_index = EXLTop->MapCnt - 1;
                Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].items_read=1;

                genetic_distance_sex_type_map = SEX_AVERAGED_GDMT;
                Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].items_read=1;
            }
        }
    }

    free(AnnotatedFileInfo.map_file_columns);

    if (EXLTop->MapCnt > 0) {
        HasMapFileBeenRead=1;
    } else {
        HasMapFileBeenRead=0;
        errorf("Fatal errors in map file.");
        errorf("Please correct and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

    Tod tod_cl("get_chromosome_list");

    NumChromo=get_chromosome_list(LTop, global_chromo_entries,
                                  chromo_loci_count, 0, analysis);
    tod_cl();

    Tod tod_lstat("write_locus_stats");
    write_locus_stats(LTop, ANNOTATED);
    tod_lstat();
    
    Tod tod_fp("read more files: freq & penetrance");
    HasFreqFileBeenRead = 
        read_annotated_freq_file(freq_file, LTop, &AnnotatedFileInfo,
			       marker_list, pheno_list, &num_groups, &groups);
    free(AnnotatedFileInfo.freq_file_columns);

    if (pen_file != NULL) {
        if (read_annotated_pen_file(pen_file, LTop,
                                    &AnnotatedFileInfo,
                                    pheno_list)) {
            HasPenFileBeenRead = -1;
        } else {
            HasPenFileBeenRead=1;
        }
        free(AnnotatedFileInfo.pen_file_columns);
    } else {
        HasPenFileBeenRead = 0;
    }
    load_classes(pheno_list, LTop);

    /* do it early so read_quant_phen has it available */
    set_missing_quant_input((linkage_ped_top *)0, analysis);

    for (i = 0; i < EXLTop->MapCnt; i++) {
        if (EXLTop->map_functions[i] == 'p') {
            bp_map = i;
            break;
        }
    }
    tod_fp();

    // CPK: Create the hash set that is used when processing plink binary data.
    // In this manner we malloc only one copy of the allele name.
    // Later we need to make sure that we only free these allele names only once!
    if (PLINK.plink == binary_PED_format && plink_info != NULL && plink_info->alleles != (char *)NULL) {
        Tod tod_hash("create_allele_strings_hashtable");
        create_allele_strings_hashtable(plink_info);
        tod_hash();
    }

    SECTION_ERR_EXTERN(FLOAT_AFFECT);

    if (Input->GetOps()->use_getops()) {
        pedfile_type = PREMAKEPED_PFT;
        basefile_type = pedfile_type;
        Top = Input->GetOps()->do_ped(LTop);
    } else if (PLINK.plink ||
	Input_Format == in_format_binary_VCF ||
        Input_Format == in_format_compressed_VCF ||
	Input_Format == in_format_VCF) {
        pedfile_type = PREMAKEPED_PFT;
        basefile_type = pedfile_type;
        // CPK: If the .ped file is a .fam file, then we process the alleles as per the .bed file..
        Tod tod_ppf("read plink ped file");
        Top = read_plink_ped_file(ped_file, bed_file, plink_info,
                                  phe_cols, LTop, EXLTop, bp_map, &AnnotatedFileInfo,
                                  num_groups, groups);
        tod_ppf();
#ifndef HIDESTATUS
	//mssgf("");
	if (PLINK.plink) {
        mssgvf("%d (of %d) markers to be included from %s\n",
               PLINK.markers_total ,PLINK.markers_included, map_file);
	}
	output_stats();
#endif /* HIDESTATUS */
    } else
        Top = read_annotated_ped_file(ped_file, LTop, &AnnotatedFileInfo,
                                      num_groups, groups);
    if (Top->PTop == NULL) {
        Top->PedRaw    = Top->Ped;
        Top->PedBroken = Top->Ped;
    }

    if (AnnotatedFileInfo.ped_file_columns)
        free(AnnotatedFileInfo.ped_file_columns);

    SECTION_ERR_FORCE(FLOAT_AFFECT);
    if (SECTION_ERR_COUNT(FLOAT_AFFECT) > 0) {
        warnvf("There were %d instances of decimal numbers read where affection status were expected.\n",
               SECTION_ERR_COUNT(FLOAT_AFFECT) );
        if (PLINK.plink) {
            warnvf("Perhaps you should rerun Mega2 and specify that the fam/ped trait is quantitative.\n");
        }
    }
    SECTION_ERR_FINI(FLOAT_AFFECT);

    count_Missing_Quant_consistency();

    if (Top == NULL) {
        HasPedFileBeenRead=0;
        errorf("Fatal errors in pedigree file.");
        errorf("Please correct and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

    HasPedFileBeenRead=1;
    Top->LocusTop = LTop;
    //    set_missing_quant_output(Top, analysis);
    
    if (EXLTop != NULL) {
        Top->EXLTop = EXLTop;
    }

    if (HasPenFileBeenRead == -1) {
        errorf("Fatal errors in penetrance file.");
        errorf("Please correct and restart.");
        EXIT(INPUT_DATA_ERROR);
    }

    log_line(mssgf);
    mssgf("Input pedigree data contains:");

    Tod tod_hraw("order_heterozygous_allele_raw ");
    order_heterozygous_allele_raw(Top);
    tod_hraw();

    Tod tod_pstat("write_ped_stats");
    write_ped_stats(Top);
    tod_pstat();

    Tod tod_omit("read annotated omit file");
    if (omit_file != NULL)	 {
        annotated_omit_file(Top, omit_file, &AnnotatedFileInfo);
        free(AnnotatedFileInfo.omit_file_columns);
    }
    if (UntypedPeds == NULL) {
        UntypedPeds = CALLOC((size_t) Top->PedCnt, int);
    }

    omit_peds(untyped_ped_opt, Top);

    clear_YLINKED_females(Top, 1, 1);
    tod_omit();

    if (database_dump || ! database_read) {
        extern void pedtree_markers_check(linkage_ped_top *LPedTreeTop, analysis_type analysis);

        Tod tod_makeped("makeped0");
        makeped(Top, analysis);  // if --db, might connect loops based on analysis
        tod_makeped();

        pedtree_markers_check(Top, analysis);
    }

    /* recode and compute frequencies if necessary */
    Mega2Status = INSIDE_RECODE;
    if (analysis != TO_PLINK || database_dump) {
        if (HasFreqFileBeenRead == 0 || HasFreqFileBeenRead == 2) {
            count_option =
                get_count_option(1, &count_halftyped,
                                 "Select individuals to compute allele frequencies\n       for recoded marker loci:");

        } else 
            count_option=4;
    } else {
        count_option=4;
    }

    {
        int ped, entrycount;
        allelecnt **member_ids;
        Tod tod_cal_init("init member_ids (for allele/freq counting)");
        member_ids=CALLOC((size_t) Top->PedCnt, allelecnt *);
        for (ped=0; ped < Top->PedCnt; ped++) {
            entrycount = (pedfile_type == POSTMAKEPED_PFT) ? Top->Ped[ped].EntryCnt:
                                               Top->PTop[ped].num_persons;
            member_ids[ped]=CALLOC((size_t) entrycount, allelecnt);
        }
        tod_cal_init(); // 606.971166
        Tod tod_cal_phen("init pheno_list: alleles, num_alleles");
        for (i=0; i < LTop->PhenoCnt; i++) {  // this code should be in mrecode.cpp:*_top()
            if (pheno_list[i].estimate_frequencies) {
                pheno_type *pheno_listi = &(pheno_list[i]);
                pheno_listi->first_allele = CALLOC((size_t) 1, allele_list_type);
                if (count_option > 0)
                    pheno_listi->first_allele->allele_freq.freq = 0.5;
                pheno_listi->first_allele->next = CALLOC((size_t) 1, allele_list_type);
                if (count_option > 0)
                    pheno_listi->first_allele->next->allele_freq.freq = 0.5;
                pheno_listi->first_allele->next->next = NULL;
                pheno_listi->num_alleles = 2;
            }
        }
        tod_cal_phen();
        Tod tod_cal("create allele list");
        Tod tod_fr(20); // 0.000,688
        Tod tod_fr_x4(20);
        SECTION_ERR_EXTERN(y_female);
        for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
            tod_fr.reset();
            tod_fr_x4.reset();
            if (marker_list[i].estimate_frequencies) {
                create_allele_list(Top, (int) i, &(marker_list[i]),
                                   member_ids, count_halftyped);
                tod_fr("create allele list 1 freq");

                if ( (count_option == 1 || count_option == 2 || count_option == 3
                      /* || analysis == TO_HWETEST || analysis == TO_SIMULATE */) &&
                     LTop->Locus[i].number != -1 ) {
                    count_allele_list(Top, (int) i, count_option, &(marker_list[i]),
                                      member_ids, count_halftyped);
                    tod_fr_x4("create allele list 4 freq");
                }
            }
        }
        for(ped=0; ped < Top->PedCnt; ped++) {
            free(member_ids[ped]);
        }
        free(member_ids);
        SECTION_ERR_FINI(y_female);
        tod_cal();
    }
    Tod tod_cra("count_raw_alleles/check recoding");
    count_raw_alleles(marker_list, Top->LocusTop);

    for (i=0; i < LTop->LocusCnt; i++) {
        if (LTop->Locus[i].Type == NUMBERED ||
            LTop->Locus[i].Type == XLINKED ||
            LTop->Locus[i].Type == YLINKED) {
            need_recoding(&(marker_list[i]));
        }

/*     if (LTop->Locus[i].Type == XLINKED || */
/* 	LTop->Locus[i].Type == YLINKED) { */
/*       LTop->Locus[i].Type = NUMBERED; */
/*     } */
    }
    tod_cra();

    Tod tod_cvt("convert_to_freq");
    if (HasFreqFileBeenRead == 0 || HasFreqFileBeenRead == 2) {
        convert_to_freq(marker_list, Top->LocusTop, count_option, analysis);

        char *bimalleles = NULL;
        if (PLINK.plink == binary_PED_format && 
            plink_info != NULL && plink_info->alleles != (char *)NULL)
            bimalleles = plink_info->alleles;

//  if (Input_Format == in_format_binary_VCF ||
//      Input_Format == in_format_compressed_VCF ||
//      Input_Format == in_format_VCF)

//        assign_dummy_alleles(marker_list, Top->LocusTop, 0);
        assign_dummy_alleles(marker_list, Top->LocusTop, bimalleles, VecAlleles);
    }
    tod_cvt();

    Tod tod_rsum("write_recode_summary");
    write_recode_summary(marker_list, Top->LocusTop, count_halftyped);
    tod_rsum();

    Tod tod_dtp("default_trait_penetrance");
    if (HasPenFileBeenRead == 0) {
        for (j=0; j < num_traits; j++) {
            default_trait_penetrances(Top, global_trait_entries[j],
                                      &(pheno_list[global_trait_entries[j]]));
        }
    }
    tod_dtp();

    Tod tod_rped("recode_ped_top"); // 199.470821
    Top->Ped = Top->PedBroken;  /* count doubleganger */
    recode_ped_top(marker_list, Top, plink_info);
    Top->Ped = Top->PedRaw;     /* don't count doubleganger */
    Mega2Status = DONE_RECODE;
    tod_rped();

    Tod tod_rloc("recode locus top");
    recode_locus_top(marker_list, pheno_list, Top->LocusTop);
    tod_rloc();

    Mega2Status = DONE_RECODE;

    log_line(mssgf);
#ifndef HIDESTATUS
    mssgf("Pedigree data summary after recoding:");
#endif
    Tod tod_pstat2("write_ped_stats 2");
    write_ped_stats(Top);
    tod_pstat2();


    Tod tod_free("free markers_item, etc");
    for (i = LTop->PhenoCnt; i < LTop->LocusCnt; i++) {
       free_marker_item(marker_list[i].first_allele);
    }

    for (i=0; i < LTop->PhenoCnt; i++) {
       lclass = pheno_list[i].class_list;
       while(lclass != NULL) {
           free(lclass->l_class.autosomal_pen);
           free(lclass->l_class.male_pen);
           free(lclass->l_class.female_pen);
           lclass = lclass->next;
       }
       free_marker_item_class(pheno_list[i].class_list);
    }

/* (rvb) 11/10/11
    As I read a group array is created in read_annotated_freq_file() and sort_and_check()
    verifies that groups on the people (if used) are in the group array.  I do not believe
    I've seen any code that puts things in the group array.
*/
    free(groups);
    if (marker_list) {
        free(marker_list + Top->LocusTop->PhenoCnt);
        marker_list=NULL;
    }
    free(pheno_list);  pheno_list=NULL;
    InputFileFormat = ANNOTATED;
    tod_free();

    tod_anfiles();
    return Top;
}

/*
static void Free_AnnotatedFileInfo(void) {
    free(AnnotatedFileInfo.names_file_columns);
    free(AnnotatedFileInfo.freq_file_columns);
    free(AnnotatedFileInfo.pen_file_columns);
    free(AnnotatedFileInfo.ped_file_columns);
    free(AnnotatedFileInfo.map_file_columns);
    free(AnnotatedFileInfo.omit_file_columns);
    free(UntypedPeds);
}
*/

static void Free_map_names(ext_linkage_locus_top *EXLTop)
{
    int i;
    for (i = 0; i < EXLTop->MapCnt; i++)
        free(EXLTop->SexMaps[i]);
    free(EXLTop->map_functions);
    for (i = 0; i < EXLTop->MapCnt; i++)
        free(EXLTop->MapNames[i]);
    free(EXLTop->MapNames);
    free(EXLTop->SexMaps);
}

static void Free_ped(linkage_ped_top *PTop)
{
    linkage_locus_top *LTop;
    linkage_locus_rec *Locus;

    ext_linkage_locus_top *EXLTop;
    ext_linkage_locus_rec *EXLocus;  

    linkage_ped_tree  *LPedT;
    linkage_ped_rec   *LPed;
    linkage_ped_rec   *Entry;

    marriage_graph_type *PPeds;
    person_node_type    *Per = NULL; //silly compiler

    int p, l, i;

    if (PTop->pedfile_type == PREMAKEPED_PFT) {  /* Raw_premake: PTop marriage_graph_type */
        /* makeped should have changed everything over to a linkage_ped_tree */
        msgvf("Free_ped: Can not possibly get here\n");

        PPeds = PTop->PTop;
        msgvf("PPeds %x # %d\n", PPeds, PTop->PedCnt);
        for (p = 0; p < PTop->PedCnt; p++) {
            for (i = 0; i < PPeds[p].num_persons; i++) {
                Per = &(PPeds[p].persons[i]);
                free(Per->pheno);
                free(Per->marker);
            }
            free(Per);
        }
        free(PPeds);
    } else {                  /* Raw_postmake: Ped linkage_ped_tree EntryCnt Entry[i] linkage_ped_tree*/
        LPedT = PTop->PedBroken; /* UniqueID[], OrigID[], FamName[] sib-pointers */
        for (p = 0; p < PTop->PedCnt; p++) {
            LPed = LPedT[p].Entry;
            for (i = 0; i < LPedT[p].EntryCnt; i++) {
                Entry = &LPed[i];
//                if (Entry->Pheno != NULL)
//                    free(Entry->Pheno);
                if (Entry->Marker != NOTYPED_ALLELES)
                    marker_free(Entry->Marker, PTop->LocusTop->PhenoCnt);
                if (Entry->loopbreakers != NULL)
                    free(Entry->loopbreakers);
            }
            if (! database_read) free(LPed);
        }
        free(LPedT);
    }

    LTop = PTop->LocusTop;
    pheno_rec *Pheno;
    if (! database_read) {
        for (l = 0; l < LTop->LocusCnt; l++) {
            Locus = &(LTop->Locus[l]);
            Pheno = &(LTop->Pheno[l]);
            switch (Locus->Type) {
            case NUMBERED:
            case XLINKED:
            case YLINKED:
                if (Locus->AlleleCnt) free(Locus->Allele);
                break;
            case AFFECTION:
    //            if (Locus->AlleleCnt) free(Allele);
    /*          wierd; old code sets Penetrance;  New (mrecode.c) code sets AvgPen then assigns Penetrance = AvgPen ... So */
                for (i = 0; i < Pheno->Props.Affection.ClassCnt; i++) {
                    if (Pheno->Props.Affection.Class[i].AutoPen != NULL)
                        free(Pheno->Props.Affection.Class[i].AutoPen);
                    if (Pheno->Props.Affection.Class[i].FemalePen != NULL)
                        free(Pheno->Props.Affection.Class[i].FemalePen);
                    if (Pheno->Props.Affection.Class[i].MalePen != NULL)
                        free(Pheno->Props.Affection.Class[i].MalePen);
                }
                free(Pheno->Props.Affection.Class);                  /* read_annotated_names_files */
                if (Locus->AlleleCnt) free(Locus->Allele);
                break;
            case QUANT:
                if (Locus->AlleleCnt) free(Locus->Allele);
                free(Pheno->Props.Quant.Mean[0]);
                free(Pheno->Props.Quant.Mean);
                free(Pheno->Props.Quant.Variance[0]);
                free(Pheno->Props.Quant.Variance);
                break;
            case TYPE_UNSET:
            case BINARY:
                if (Locus->AlleleCnt) free(Locus->Allele);
                break;
            }
        }
        free(LTop->Locus);
        free(LTop->Pheno);
        if (LTop->Marker != 0)
            free( ((marker_rec *)LTop->Marker) + LTop->PhenoCnt );
    } else {
        extern linkage_allele_rec *Allele_rec;
        extern linkage_affection_class *AffectClass_rec;

        delete [] Allele_rec;
        delete [] AffectClass_rec;
        delete [] (((marker_rec *)LTop->Marker) + LTop->PhenoCnt);
        delete [] LTop->Pheno;
        delete [] LTop->Locus;
    }
    free(LTop->MaleRecomb);

    EXLTop = PTop->EXLTop;                                /* read_annotated_map_files */
    if (EXLTop != NULL) {
        for (l = LTop->PhenoCnt; l < EXLTop->LocusCnt; l++) {
            EXLocus = &(EXLTop->EXLocus[l]);
            free(EXLocus->positions);
            free(EXLocus->pos_male);
            free(EXLocus->pos_female);                       /* read_annotated_map_files */
        }
        if (EXLTop->EXLocus)
            free( EXLTop->EXLocus + LTop->PhenoCnt );
    }
    if (! database_read) free(LTop);

/* Free_map_names
    free(EXLTop->map_functions);
    free(EXLTop->MapNames);
    for (l = 0; l < EXLTop->MapCnt; l++)
        free(EXLTop->SexMaps[l]);
    free(EXLTop->SexMaps);
*/
    free(EXLTop);

    if (! database_read) free(PTop);

    free(global_chromo_entries);
    free(chromo_loci_count);
    if (NumUnmapped > 0) {
        free(unmapped_markers);
    }
}

void Free_annotated_files(linkage_ped_top *Top)
{
//  Free_AnnotatedFileInfo(); //Now doing this immediately after op
    if (Top->EXLTop != NULL)
        Free_map_names(Top->EXLTop);
    Free_ped(Top);
}

/* function to check if pedigree, names and map file are in annotated
   format. Read the first line from each file, and check for some
   keywords (code borrowed from parse_variable_width_header).

*/

int check_annotated_file_format(char *input_files[])
{
    FILE *fp;
    char buffer[READ_CHUNK+1];
    int ifl, flags[2];
#ifndef HIDEFILE
    int format_flags[4];
#endif
    int read_line, checksum=0;
    char *token;
    char *overflow;
    size_t overflow_size, read_size, read_length;

#ifndef HIDESTATUS
    mssgf("Checking format of input files ....");
#endif
    for (ifl=0; ifl < 4; ifl++) {
#ifndef HIDEFILE
        format_flags[ifl]=0;
#endif
        if (input_files[ifl] == NULL) continue;
	if ((fp = fopen(input_files[ifl], "r")) == NULL) {
            errorvf("could not open %s for reading!\n", input_files[ifl]);
            EXIT(FILE_READ_ERROR);
	}
        flags[0]=flags[1]=0;
        read_line=1;
        overflow = &buffer[READ_CHUNK];
        while(read_line) {
            strcpy(buffer, "");
            if (overflow < &buffer[READ_CHUNK]) {
                strcpy(buffer, overflow);
                overflow_size = &buffer[READ_CHUNK-1] - overflow;
            } else {
                overflow_size = 0;
            }
            read_size = READ_CHUNK - overflow_size;

            IgnoreValue(fgets(buffer + overflow_size, (int)read_size, fp));
            read_length = strlen(buffer);
            /* check if line overflows the buffer, and there is a string that
               spans two consecutive reads */
            if (read_length == READ_CHUNK - 1) {
                overflow = buffer + read_length;
                while (*--overflow != ' ' && *overflow != '\t' &&  *overflow != '\n') ;
                *overflow++ = 0;
            }
            if ((read_length < (READ_CHUNK-1)) || (buffer[strlen(buffer)-1] == '\n')) {
                /* reached end of line */
                read_line=0;
            }

            token = strtok(buffer, "\t \n");
            while (token != NULL && (flags[0] == 0 || flags[1] == 0)) {
                switch(ifl) {
               case 0:
                    /* pedigree file */
                    /* Check if header has Pedigree and ID */
                    if (!strcasecmp(token, "Pedigree")) {
                        flags[0] = 1;
                    } else if (!strcasecmp(token, "ID")) {
                        flags[1] = 1;
                    }
                    break;

                case 1:
                    /* names file */
                    if (!strcasecmp(token, "Name")) {
                        flags[0] = 1;
                    } else if (!strcasecmp(token, "Type")) {
                        flags[1] = 1;
                    }
                    break;

                case 2:
                    /* map file, check the first header that is not a keyword */
                    if (strcasecmp(token, "Chromosome") &&
                        strcasecmp(token, "Name") &&
                        strcasecmp(token, "Error")) {
                        /* Not any of these keywords, check the extensions */
                        if (is_valid_hdr_ext(token, map_valid_header_item_extensions) > 0) {
                            flags[0] = flags[1] = 1;
                        }
                    }
                    break;
                case 3:
                    /* omit file, just check to see if the first token matches a header */
                    if (!strcasecmp(token, "Pedigree") ||
                        !strcasecmp(token, "Individual") ||
                        !strcasecmp(token, "Marker")) {
                        flags[0]=flags[1] = 1;
                    }
                    break;
                }
                token = strtok(NULL, "\t \n");
            }

            if (flags[0] && flags[1]) {
                /* no need to read further */
#ifndef HIDEFILE
                format_flags[ifl] = 1;
#endif
                checksum++; break;
            }
        }
#ifndef HIDEFILE
        if (format_flags[ifl] == 1) {
            mssgvf("\t%s is in MEGA2 format.\n", input_files[ifl]);
        } else {
            mssgvf("\t%s is in LINKAGE format.\n", input_files[ifl]);
        }
#endif
        fclose(fp);
    }

    strcpy(err_msg, "");
    /* end of 1st line */

    /* Now see if all files are in annotated format */
    return checksum;
}

void write_annotated_aff(FILE *filep, int locusnm, linkage_locus_rec *locus,
                         linkage_ped_rec *entry)
{
    fprintf(filep, "  ");
    if (entry->Pheno[locusnm].Affection.Status == UNDEF ||
        entry->Pheno[locusnm].Affection.Status == 0) {
        fprintf(filep, "NA");
    } else {
        fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Status);
    }

    if (locus->Pheno->Props.Affection.ClassCnt > 1) {
        fprintf(filep, " ");
        if (entry->Pheno[locusnm].Affection.Class == UNDEF ||
            entry->Pheno[locusnm].Affection.Class == 0) {
            fprintf(filep, "NA");
        } else {
            fprintf(filep, "%d", entry->Pheno[locusnm].Affection.Class);
        }
    }
}

void write_annotated_quant(FILE *filep, int locusnm,
                           linkage_ped_rec *entry)
{
    fprintf(filep, "  ");
    if (fabs(entry->Pheno[locusnm].Quant - MissingQuant) <= EPSILON) {
        fprintf(filep, "    NA    ");
    } else {
        fprintf(filep, "%10.5f", entry->Pheno[locusnm].Quant);
    }
    return;
}


void write_annotated_numbered(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int all1, all2;
    get_2alleles(entry->Marker, locusnm, &all1, &all2);

    fprintf(filep, " ");

    if (all1 == 0)
        fputs("NA ", filep);
    else
        fprintf(filep, "%2s ", format_allele(locus, all1));

    if (all2 == 0)
        fputs("NA", filep);
    else
        fprintf(filep, "%2s", format_allele(locus, all2));
}



#ifdef num_alleles
#undef num_alleles
#endif

#ifdef num_classes
#undef num_classes
#endif


// //////////////////////////////////////////////////////////////// //

void
PLINK_clr(enum PLINK_FORMAT plink)
{
    memset(&PLINK, 0, sizeof (PLINK_t));
    PLINK.plink = plink;
//  if (!plink) return;
    strcpy(PLINK.trait, "default");
    PLINK.missing_pheno = 1;
    PLINK.pheno_value = -9;
}

int PLINK_args(char *str, int xcf)
{
    size_t lstr = strlen(str)+1;
    char *par = CALLOC(lstr, char);
    char *tok;
    int err = 0;

    strcpy(par, str);
    tok = strtok(par, "\t ");

    while(1) {
        if (tok == NULL) break;
        if (!strcasecmp(tok, "clear")) {
            memset(&PLINK, 0, sizeof(PLINK_t));
        } else {
            if (! xcf) {
                if (!strcasecmp(tok, "--bfile")) {
                    if (InputMode == INTERACTIVE_INPUTMODE) {
                        printf("Usage: Bad parameter for interactive input  %s\n", tok);
                        err++;
                    } else
                        PLINK.plink = binary_PED_format;
                } else if (!strcasecmp(tok, "--file")) {
                    if (InputMode == INTERACTIVE_INPUTMODE) {
                        printf("Usage: Bad parameter for interactive input  %s\n", tok);
                        err++;
                    } else
                        PLINK.plink = PED_format;
                } else if (!strcasecmp(tok, "--no-fid")) {
                    PLINK.no_fid = 1;
                } else if (!strcasecmp(tok, "--no-parents")) {
                    PLINK.no_parents = 1;
                } else if (!strcasecmp(tok, "--no-pheno")) {
                    PLINK.no_pheno = 1;
                } else if (!strcasecmp(tok, "--map3")) {
                    PLINK.map3 = 1;
                } else if (!strcasecmp(tok, "--cM")) {
                    PLINK.cM = 1;
                } else if (!strcasecmp(tok, "--missing-phenotype")) {
                    tok = strtok(NULL, "\t ");
                    PLINK.missing_pheno = 1;
                    PLINK.pheno_value = atof(tok);
                } else if (!strcasecmp(tok, "--trait")) {
                    tok = strtok(NULL, "\t ");
                    strcpy(PLINK.trait, tok);
                } else if (!strcasecmp(tok, "--affectionstatus")) {
                    PLINK.traitType = 0;
                } else if (!strcasecmp(tok, "--quantitative")) {
                    PLINK.traitType = 1;
                } else if (!strcasecmp(tok, "--kosambi")) {
                    PLINK.geneticMapType = 0;
                } else if (!strcasecmp(tok, "--haldane")) {
                    PLINK.geneticMapType = 1;
                } else if (!strcasecmp(tok, "--1")) {
                    PLINK.phenotypeCoding = 1;
                } else {
                    printf("Usage: Bad parameter %s\n", tok);
                    err++;
                }
            } else {
                if (!strcasecmp(tok, "--missing-phenotype")) {
                    tok = strtok(NULL, "\t ");
                    PLINK.missing_pheno = 1;
                    PLINK.pheno_value = atof(tok);
                } else if (!strcasecmp(tok, "--trait")) {
                    tok = strtok(NULL, "\t ");
                    strcpy(PLINK.trait, tok);
                } else if (!strcasecmp(tok, "--affectionstatus")) {
                    PLINK.traitType = 0;
                } else if (!strcasecmp(tok, "--quantitative")) {
                    PLINK.traitType = 1;
                } else {
                    printf("Usage: Bad parameter %s\n", tok);
                    err++;
                }
            }
        }
        tok = strtok(NULL, "\t ");
    }
    PLINK.xcf = xcf;


    if (Input_Format == in_format_binary_VCF ||
        Input_Format == in_format_compressed_VCF ||
	Input_Format == in_format_VCF) {
        // nada
    } else if (PLINK.plink != binary_PED_format && PLINK.plink != PED_format) {
//for NEW batch files w/o --bfile/file
        if (Input_Format == in_format_binary_PED)
            PLINK.plink = binary_PED_format;
        else if (Input_Format == in_format_PED)
            PLINK.plink = PED_format;
        else {
            printf("Usage: You must specify PLINK file type: --file or --bfile\n");
            err++;
        }
    }
    if (PLINK.no_pheno == 0 && PLINK.trait[0] == 0) {
        printf("Usage: You must specify a trait name for the pedigree file phenotype\n");
        err++;
    }
    if (err > 0) {
        PLINK_usage(xcf);
    }

    PLINK_str(str, FILENAME_LENGTH);

    free(par);
    return err? 0 : 1;
}

void
PLINK_usage(int xcf)
{
    if (! xcf) {
        printf("\nPLINK options: --no-fid --no-parents --no-pheno --1 --map3 ...\n");
        printf("    --missing-phenotype <value> --cM --kosambi/--haldane ...\n");
        printf("    --trait <value> --affectionstatus/--quantitative\n\n");
    } else {
        printf("\nPLINK options: --missing-phenotype <value> --trait <value> --affectionstatus/--quantitative\n\n");
    }
}

void
PLINK_str(char *ans, int len)
{
    ans[0] = 0;
//  if (PLINK.plink == binary_PED_format)
//      strcat(ans, "--bfile");
//  else if (PLINK.plink == PED_format)
//      strcat(ans, "--file");
    if (PLINK.no_fid)
        strcat(ans, " --no-fid");
    if (PLINK.no_parents)
        strcat(ans, " --no-parents");
    if (PLINK.no_pheno)
        strcat(ans, " --no-pheno");
    if (PLINK.map3)
        strcat(ans, " --map3");
    if (PLINK.geneticMapType != 0)
        strcat(ans, " --haldane");
    if (PLINK.cM)
        strcat(ans, " --cM");
    if (PLINK.phenotypeCoding)
        strcat(ans, " --1");
    if (PLINK.traitType)
        strcat(ans, " --quantitative");

    if (PLINK.missing_pheno != 0) {
        if ( ((int)PLINK.pheno_value) == PLINK.pheno_value) {
            sprintf(&ans[strlen(ans)], " --missing-phenotype %d", (int)PLINK.pheno_value);
        } else {
            sprintf(&ans[strlen(ans)], " --missing-phenotype %.4f", PLINK.pheno_value);
        }
    }
    if (PLINK.trait[0] != 0) {
        strcat(ans, " --trait ");
        strcat(ans, PLINK.trait);
    }
}

static int parse_phe_types(char *phe_file, char ***phe_names, int **phe_types)
{
    FILE  *fp = phe_file ? fopen(phe_file, "r") : NULL;
    int    cols, Xcols = 1;

    if (fp == NULL) {
        *phe_names = (char **)CALLOC((size_t) Xcols, char **);
        *phe_types = (int *)CALLOC((size_t) Xcols, int);
        return 0;
    }
    fclose(fp);

    cols = phemake(phe_file);
    if (cols < 0) {
        errorvf("phe file read failure: %s\n", phe_file);
        EXIT(FILE_READ_ERROR);
    }

    if (PLINK.no_pheno == 0) /* reserve extra row for ped file trait */
        Xcols = cols + 1;    /* it will be there but not counted */
    else
        Xcols = cols;

#ifndef HIDEFILE
    msgvf("reading phenotype file %s ... (%d columns)\n", phe_file, cols);
#endif

    *phe_names = (char **)CALLOC((size_t) Xcols, char **);
    *phe_types = (int *)CALLOC((size_t) Xcols, int);

    phehdr(*phe_names, *phe_types);

    if (cols && PLINK.map3 == 0) {
        int i;
        for (i = 0; i < cols; i++) {
            if ((*phe_types)[i] == 0 && strcmp((*phe_names)[i], PLINK.trait) == 0) {
                errorvf("Pedigree file trait name (\"%s\") matches a name in the Phenotype file.\n", PLINK.trait);
                errorvf("This is not allowed; one of the two names must be changed.\n");
                EXIT(DATA_INCONSISTENCY);
            }
        }
    }

    return cols;
}


static int parse_plink_map_as_names_file_header(FILE *filep,
                                                annotated_file_desc *file_desc)
{
    /* Names file has two columns always: Name, Type in any order */
    /* extra information is allowed on each line after the two columns */

    /* initialize the value type and column names of the two elements */
    file_desc->num_names_cols=2;
    file_desc->names_file_columns = CALLOC((size_t) 2, col_hdr_type);

    INIT_COLNAME(file_desc->names_file_columns, 0, CHAR_AN, "Type");
    INIT_COLNAME(file_desc->names_file_columns, 1, STRING_AN, "Name");

// Type
    file_desc->names_file_columns[0].input_col = 0;
// Name
    file_desc->names_file_columns[1].input_col = 2;

    if (file_desc->names_file_columns[1].input_col == -1 ||
        file_desc->names_file_columns[0].input_col == -1) 	{
        /* not in names_file format */
        return 0;
    }

    /*  free(buffer); */
    if (is_missing_required_hdr(file_desc->names_file_columns, 2)) {
        errorf("Found errors in names file header, please correct and restart Mega2.");
        EXIT(DATA_INCONSISTENCY);
    }
    return 2;
}

static int read_plink_map_as_names_file(char *map_file, linkage_locus_top **LTop,
					int cols, char **phe_names, int *phe_types,
					annotated_file_desc *file_desc)
{
    int annotated_format;

    FILE *fp = map_file ? fopen(map_file, "r") : NULL;
    if (fp == NULL) {
        errorvf("could not open %s for reading!\n", map_file);
        EXIT(FILE_READ_ERROR);
    }

#ifndef HIDEFILE
    sprintf(err_msg, "Reading PLINK map file for names: %s", map_file);
    mssgf(err_msg);
#endif
    annotated_format = parse_plink_map_as_names_file_header(fp, file_desc);
    if (annotated_format > 0) {
//xx    if (PLINK.no_pheno == 0) cols++;
        *LTop = read_marker_only_data(fp, cols, phe_names, phe_types);
    }
    return annotated_format;
}

static ext_linkage_locus_top *read_plink_map_file(const char *map_file,
						  linkage_locus_top *LTop,
						  annotated_file_desc *file_desc,
						  plink_info_type *plink_info)
{
    int col;
    int num_map_recs;
    col_hdr_type tmp[1];
    col_hdr_type reserved_colnames[NUM_MAPCOL_NAMES];
    col_hdr_type *map_all_colnames;
    int num_userdef_cols;
//
    FILE *mapfp = map_file ? fopen(map_file, "r") : NULL;
    if (mapfp == NULL) {
        errorvf("could not open %s for reading!\n", map_file);
        EXIT(FILE_READ_ERROR);
    }

//  skip_line(mapfp, chr);
    num_map_recs = linecount(mapfp);
    
    // NOTE: To exclude a SNP from analysis, set the 4th column (physical base-pair position)
    // to any negative value (this will only work for MAP files, not for binary BIM files).
    // So if the pase-pair position is non-negative PLINK.markers_included will be incremented.
    PLINK.markers_total = num_map_recs;
    
    if (PLINK.plink == binary_PED_format) {
        // http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        // If the map. file is actually a .bim file, then we need to save the alleles information
        plink_info->alleles = (char *)CALLOC((size_t) 2*num_map_recs, char);
        plink_info->allele_count = 2*num_map_recs;
    }
    
// decide 3 or 4 columns, by reading file
    num_userdef_cols = columncount(mapfp);
#ifndef HIDEFILE
    msgvf("Reading map file %s ... (%d columns)\n", map_file, num_userdef_cols);
#endif
    file_desc->num_map_cols = num_userdef_cols;
    file_desc->map_file_columns =
        CALLOC((size_t)file_desc->num_map_cols, col_hdr_type);
    map_all_colnames = &(file_desc->map_file_columns[0]);

    INIT_COLNAME(reserved_colnames, 0, STRING_AN, "Chromosome");
    INIT_COLNAME(reserved_colnames, 1, STRING_AN, "Name");
    INIT_COLNAME(reserved_colnames, 2, FLOAT_AN, "Error");

    reserved_colnames[0].input_col = 0;
    copy_colname(reserved_colnames + 0, map_all_colnames + 0);
    map_all_colnames[0].output_col = 0;

    reserved_colnames[1].input_col = 1;
    copy_colname(reserved_colnames + 1, map_all_colnames + 1);
    map_all_colnames[1].output_col = 1;

    col = 2;
    // Exactly 6 columns for a .bim file...
    if ( (file_desc->num_map_cols == 4 && (PLINK.plink != binary_PED_format)) || 
         (file_desc->num_map_cols == 6 && (PLINK.plink == binary_PED_format))) {
        if (PLINK.map3 == 1) {
            errorvf("You have indicated the map file has 3 columns (via --map3), but it has 4.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    } else if ( (file_desc->num_map_cols == 3 && (PLINK.plink != binary_PED_format)) || 
                (file_desc->num_map_cols == 5 && (PLINK.plink == binary_PED_format))) {
        if (PLINK.map3 != 1) {
            errorvf("The map file only has 3 columns.  You must specify the PLINK parameter --map3.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    }

    if (PLINK.map3 == 0) {
        if (PLINK.geneticMapType == 0) {
            INIT_COLNAME(tmp, 0, FLOAT_AN, "Map.k.a"); //kosambi
        } else {
            INIT_COLNAME(tmp, 0, FLOAT_AN, "Map.h.a"); //haldane
        }
        tmp[0].input_col = col;
        copy_colname(tmp, map_all_colnames + col);
        map_all_colnames[col].output_col = col;
        col++;
    }

    INIT_COLNAME(tmp, 0, INT_AN, "BP.p");
    tmp[0].input_col = col;
    copy_colname(tmp, map_all_colnames + col);
    map_all_colnames[col].output_col = col;
    col++;

//11
    if (PLINK.plink == binary_PED_format) {
        map_all_colnames[col].ColName = TKN.MT;
        map_all_colnames[col].output_col = col;
        col++;
        map_all_colnames[col].ColName = TKN.MT;
        map_all_colnames[col].output_col = col;

    }

    return read_common_map_file(mapfp, map_file, LTop, 0, reserved_colnames, file_desc, PLINK.plink == binary_PED_format ? &plink_info->alleles : NULL);
}

static linkage_ped_top *read_plink_ped_file(char *pedfile,
                                            char *bedfile,
                                            plink_info_type *plink_info,
                                            int phecols,
                                            linkage_locus_top *LTop,
                                            ext_linkage_locus_top *EXLTop,
                                            int bp_map,
                                            annotated_file_desc *file_desc,
                                            int num_groups,
                                            int *groups)
{
    int  num_userdef_cols, i, j;
    int  has_extra_ids = 0;
    int  output_col = 0;
    int  num_ped_records;
    char loctype[2];
    linkage_locus_rec llr;
    ext_linkage_locus_rec llx;
    double pos;

    col_hdr_type *ped_all_colnames;
    col_hdr_type *colname_item, tmp[1];

    int  SNP_major = 0;
    FILE *bed_filep = (FILE *)NULL;

    Tod pf_hdr("read plink ped file");
    FILE *filep = pedfile ? fopen(pedfile, "r") : NULL;
    if (filep == NULL) {
        errorvf("could not open %s for reading!\n", pedfile);
        EXIT(FILE_READ_ERROR);
    }

    init_reserved_pedcol_names();

    num_userdef_cols = columncount(filep); /* may include 1 trait; if so then was put in LocusCnt */

    // http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
    // CPK: If the .ped file is actually a .fam file, then open the .bed file to get the genotype inforamtion.
    // We also do some checks on the .bed file for the proper version number, and read the SNP/Individual flag...
    if (PLINK.plink == binary_PED_format) {
#ifndef HIDEFILE
        msgvf("Checking PLINK Binary format file: %s\n", bedfile);
#endif
        if (plink_info->allele_count == 0) {
            msgvf("WARNING... No alleles were discovered when the PLINK .bim file was read.\n");
        }
        bed_filep = bedfile ? fopen(bedfile, read_binary) : NULL;
        if (bed_filep != (FILE *)NULL) {
            char magic[2];
            magic[0] = (char)fgetc(bed_filep);
            magic[1] = (char)fgetc(bed_filep);
            // look for the magic numbers for a PLINK V1.00 .bed file
            if (magic[0] == 0x6c && magic[1] == 0x1b) {
                SNP_major = (char)fgetc(bed_filep);
                /*
                 a value of 00000001 indicates SNP-major (i.e. list all individuals for first SNP, all individuals
                 for second SNP, etc) whereas a value of 00000000 indicates individual-major (i.e. list all SNPs
                 for the first individual, list all SNPs for the second individual, etc). The default is SNP-major mode.
                 */
#ifndef HIDEFILE
                if (SNP_major == 0x01) {
                    msgvf("Detected PLINK .bed file v1.00 SNP-major mode.\n");
                } else {
                    msgvf("Detected PLINK .bed file v1.00 Individual-major mode.\n");
                }
#endif
            } else {
                msgvf("The .bed file found is not v1.00 format.\n");
            }
        } else {
            msgvf("A PLINK .fam file was specified, but no corresponding .bed file was found.\n");
        }
    } else {
#ifndef HIDEFILE
        msgvf("Reading PLINK .ped file: %s (%d columns).\n", pedfile, num_userdef_cols);
#endif
    }
    if (PLINK.plink == binary_PED_format) {
        extern int MARKER_SCHEME3_check;
        plink_info->SNP_major = SNP_major;
        plink_info->bed_filep = bed_filep;
        MARKER_SCHEME3_check = 0;
    }
    if (phecols > 0)                    /* but not in ped file */
        num_userdef_cols += phecols;
    if (PLINK.no_fid)
        num_userdef_cols += 1;          /* in "virtual file" but not in actual file */
    if (PLINK.no_parents)
        num_userdef_cols += 2;          /* in "virtual file" but not in actual file */

#ifndef HIDEFILE
//z msgvf("Reading PLINK format ped file: %s (%d columns)\n", pedfile, num_userdef_cols);
#endif

    // Need to add in the alleles that were gathered from reading the .bim file,
    // since this information will not be found in the .map file.
    // Note that allele_count = 0 if we have read a .map file...
//z i = num_userdef_cols + (PLINK.plink == binary_PED_format ? plink_info->allele_count : 0);
    i = num_userdef_cols + (PLINK.plink == PED_format ? 0 : 2 * LTop->MarkerCnt);
    ped_all_colnames = CALLOC((size_t)(i), col_hdr_type);
#ifdef SHOWSTATUS
    msgvf("ALLOC SPACE: ped_all_colnames: %d MB (%d x %d)\n",
          i * sizeof(col_hdr_type) / 1024 / 1024,
          i,  sizeof(col_hdr_type));
#endif
    for (i=0; i < 5; i++) {
        ReservedColnames[i].input_col = output_col;
        colname_item = &(ped_all_colnames[ReservedColnames[i].input_col]);
        copy_colname(&(ReservedColnames[i]), colname_item);
        colname_item->output_col = output_col;
        output_col++;
    }
//11
    for (; i < num_userdef_cols ; i++) {
        ped_all_colnames[i].ColName = TKN.MT;
    }
//z  for (i=0; i < LTop->LocusCnt && i < num_userdef_cols ; i++)
    for (i=0; i < LTop->LocusCnt; i++) {
        llr = LTop->Locus[i];
        llx = EXLTop->EXLocus[i];
        loctype_to_descriptor(LTop, i, loctype);
        if ( (loctype[0] == 'M' && loctype[1] == 'M') ||
             (loctype[0] == 'X' && loctype[1] == 'X') ||
             (loctype[0] == 'X' && loctype[1] == 'M') ||
             (loctype[0] == 'Y' && loctype[1] == 'Y') ||
             (loctype[0] == 'Y' && loctype[1] == 'M') ||
             (loctype[0] == 'U' && loctype[1] == 'U')) {
            for (j = 0; j < /*llr.AlleleCnt*/2; j++) {
                INIT_COLNAME(tmp, 0, STRING_AN, llr.LocusName);
                tmp[0].input_col = output_col;
                copy_colname(tmp, ped_all_colnames + output_col);
                ped_all_colnames[output_col].output_col = output_col;
                ped_all_colnames[output_col].locus_number = i /* locus */;
                pos = llx.positions[bp_map];
                if (pos < 0 && pos != UNKNOWN_POSITION) ped_all_colnames[output_col].value_type = IGNORE;
                output_col++;
            }
        } else if ( (loctype[0] == 'A' && loctype[1] == 'A')) {
            INIT_COLNAME(tmp, 0, STRING_AN, llr.LocusName);
            tmp[0].input_col = output_col;
            copy_colname(tmp, ped_all_colnames + output_col);
            ped_all_colnames[output_col].output_col = output_col;
            ped_all_colnames[output_col].locus_number = i /* locus */;
            output_col++;
        } else {
            INIT_COLNAME(tmp, 0, STRING_AN, llr.LocusName);
            tmp[0].input_col = output_col;
            copy_colname(tmp, ped_all_colnames + output_col);
            ped_all_colnames[output_col].output_col = output_col;
            ped_all_colnames[output_col].locus_number = i /* locus */;
            output_col++;
        }
    }

    pf_hdr();
    file_desc->ped_file_columns = ped_all_colnames;
    file_desc->num_ped_cols = num_userdef_cols;

    num_ped_records = linecount(filep);

    UNIQUEid = has_extra_ids;
    return read_common_ped_file(filep, pedfile, plink_info, LTop, file_desc,
                                phecols, num_groups, groups,
                                num_ped_records, has_extra_ids);
}

////////////////////////////////////////////////////////////////
// Eventually all map functions should be in their own file NOT here
////////////////////////////////////////////////////////////////

using namespace std;

m2_map::m2_map(const string name, const char function) {
    if (!(function == 'p' || function == 'h' || function == 'k')) {
        errorf("The map function must be one of: p, h, or k.");
        EXIT(SYSTEM_ERROR);
    }
    if (name.length() < 3) {
        errorf("The map name lenth must not be zero.");
        EXIT(SYSTEM_ERROR);
    }
    full_name= name + "." + string(1, function);
};

// see linkage_ext.h
const linkage_locus_type m2_map_entry::get_locus_type() {
    if (chr == SEX_CHROMOSOME) return XLINKED;
    else if (chr == MALE_CHROMOSOME) return YLINKED;
    return NUMBERED; // BINARY???
};

const string m2_map_entry::get_chr_type_string() {
    if (chr == SEX_CHROMOSOME) return "X";
    else if (chr == MALE_CHROMOSOME) return "Y";
    return "M"; // NOTE: XY and MT are lumped in here.
};

void m2_map_entry::set_chr(const string CHROM) {
    if (CHROM == "chrX" || CHROM == "X") {
        chr = SEX_CHROMOSOME;
    } else if (CHROM == "chrY" || CHROM == "Y") {
        chr = MALE_CHROMOSOME;
    } else if (CHROM == "chrXY" || CHROM == "XY") {
        chr = PSEUDO_X;
    } else if (CHROM == "chrMT" || CHROM == "MT") {
        chr = MITO_CHROMOSOME;
    } else {
        int chrm;
        istringstream ss(CHROM);
        // If we can convert it into a number then use it, otherwise designate it as unknown.
        if (!(ss >> chrm)) this->chr = UNKNOWN_CHROMO; // "U"
        else this->chr = chrm;
    }
};

void m2_map::gc()
{
    this->entries.clear();
    vector<m2_map_entry>().swap(this->entries);
}
