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
#include <math.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"
#include "mrecode.h"
#include "annotated_ped_file.h"
#include "error_sim.h"

#include "error_messages_ext.h"
#include "compress_ext.h"
#include "fcmap_ext.h"
#include "grow_string_ext.h"
#include "linkage_ext.h"
#include "list_ext.h"
#include "makeped_ext.h"
#include "marker_lookup_ext.h"
#include "mrecode_ext.h"
#include "omit_ped_ext.h"
#include "plink_ext.h"
#include "read_files_ext.h"
#include "reorder_loci_ext.h"
#include "utils_ext.h"
#include "user_input_ext.h"
#include "batch_input_ext.h"
#include "write_files_ext.h"

/*
     error_messages_ext.h:  empty_file errorf mssgf my_calloc my_malloc warnf
              fcmap_ext.h:  fcmap
        grow_string_ext.h:  grow
            linkage_ext.h:  clear_laffclass clear_lallelerec clear_llocusdata clear_llocusrec copy_lpedrec copy_lpedrec_raw_alleles copy_lpedtree1 free_all_from_lpedtree get_llooprec malloc_lpedtree_entries new_llocustop new_lpedrec new_lpedtop new_lpedtree
               list_ext.h:  append_to_list_tail list_insert_when_told new_list pop_first_list_entry
            makeped_ext.h:  check_pre_makeped read_pre_makeped
      marker_lookup_ext.h:  add_allele search_allele
            mrecode_ext.h:  check_locus_file_format create_full_marker_data read_marker_data
           omit_ped_ext.h:  omit_peds
         read_files_ext.h:  read_aff_phen read_quant_phen
       reorder_loci_ext.h:  get_chromosome_list
              utils_ext.h:  EXIT draw_line log_line
        write_files_ext.h:  write_locus_stats write_ped_stats
*/

#ifdef _WIN
#define isnan(x) _isnan(x)
#endif

/* static functions */
static int read_quantitative_data(FILE *filep, int locusnm,
                                  pheno_rec *locus, linkage_ped_rec *entry);
static int read_affection_data(FILE *filep, int locusnm,
			       pheno_rec *locus,
			       linkage_ped_rec *entry);
static int read_binary_data(FILE *filep, int locusnm,
			    linkage_locus_rec *locus,
			    linkage_ped_rec *entry);
static int read_linkage_record(FILE *filep, linkage_ped_rec *entry,
			       linkage_locus_top *LTop, lrec_data *lrecdata,
			       int *col2locus, int linenum, char *filename,
                               int *origuniq);
static int read_linkage_list_insert_fun(list_data EntryData,
					list_data NewData, int flag);
static linkage_locus_top *read_linkage_locus_file(FILE *filep, int linkagecols, int **col2locus);

/* to be used externally */
int read_numbered_data(FILE *filep, int locusnm,
                       linkage_locus_rec *locus, void *entry,
		       record_type rec, int last);

int read_premakeped_affec(FILE *filep, int locusnm, linkage_locus_rec *locus,
			  pre_makeped_record *entry);
int read_premakeped_quant(FILE *filep, int locusnm,
			  pre_makeped_record *entry);
int read_premakeped_bin(FILE *filep, int locusnm, linkage_locus_rec *locus,
			pre_makeped_record *entry);
int read_premakeped_num(FILE *filep, int locusnm, linkage_locus_rec *locus1,
			pre_makeped_record *entry);
linkage_ped_top *read_linkage2(char *pedfl_name, char *locusfl_name,
			       char *omitfl_name, char *mapfl_name,
			       int untyped_ped_opt,
			       analysis_type analysis);
static ext_linkage_locus_top *read_map_file(char *mapfl_name,  linkage_locus_top *LTop1);
linkage_ped_top *read_linkage_ped_file(FILE *filep, linkage_locus_top *LTop, int *col2loc);

void check_size_dos(void);
int check_dos_file(FILE *fp);
int check_empty(char *file_name);

void premakeped_omit_file(linkage_ped_top *Top,
			  const char *omitfl_name,
                          const int raw_allele);

SUPPRESS_MSSG_NESTED_INIT(BAD_QMISSING_OUT);
SUPPRESS_MSSG_NESTED_INIT(FLOAT_AFFECT);
SUPPRESS_MSSG_NESTED_INIT(illegal_affect);
SUPPRESS_MSSG_NESTED_INIT(illegal_quant);

/**************************************/
static int check_ped_file_cols(int cols, FILE *fp, char *file_name)
{

    /* special function to see if the pedigree file has the correct
       number of columns
    */

    int num_read, linenum, lch;
    char rest[400];

    linenum=1;
    while(!feof(fp)) {
        num_read=0;
        while(num_read < cols) {
            lch=getc(fp);
            if (feof(fp)) {
                rewind(fp);
                return 1;
            }
            while(isspace(lch)) {
                if ((lch == LF ||  lch == CR) && (num_read < cols)) {
                    sprintf(err_msg, "File %s, Unexpected end of line %d",
                            file_name, linenum);
                    errorf(err_msg);
                    sprintf(err_msg, "Expecting %d cols, read in %d.",
                            cols, num_read);
                    errorf(err_msg);
                    return -2;
                }
                lch=fgetc(fp);
            }
            ungetc(lch, fp);
            while(!isspace(lch)) {
                lch=fgetc(fp);
            }
            ungetc(lch, fp);
            num_read++;
        }
        (void)fgets(rest, 399, fp);
        linenum++;
    }

    rewind(fp);
    return 1;
}

int read_premakeped_quant(FILE *filep, int locusnm, 
                          pheno_rec *locus, pre_makeped_record *entry)
{
    int lch = read_quant_phen(filep, locusnm, locus, (void *) entry, Premakeped);
    return lch;
}

static int read_quantitative_data(FILE *filep, int locusnm, pheno_rec *locus, linkage_ped_rec *entry)
{
    int lch = read_quant_phen(filep, locusnm, locus, (void *) entry, Postmakeped);
    return lch;
}

static int read_affection_data(FILE *filep, int locusnm, pheno_rec *locus, linkage_ped_rec *entry)
{
    int lch=read_aff_phen(filep, locusnm, locus, (void *)entry,
                          Postmakeped);
    return lch;
}

static void check_quant_vs_MissingQuants(double quant)
{
    if (MissingOutQuantSet && (quant != QUNDEF && quant != QMISSING)) {
        if (isnan(QuantOutLow) || quant <= QuantOutLow)
            QuantOutLow = quant;
        if (isnan(QuantOutHi)  || quant >= QuantOutHi)
            QuantOutHi  = quant;

        if (!isnan(QuantOutLow) && MissingOutQuant >= QuantOutLow && MissingOutQuant <= QuantOutHi) {
            SUPPRESS_MSSG_NESTED(BAD_QMISSING_OUT);
            errorvf("The missing quantitative output value (%g) is among the values found \n", MissingOutQuant);
            errorvf("in the set of quantitative values [%g..%g].\n", QuantOutLow, QuantOutHi);
        }
    }
}

/**
   This routine will set the value of Quant to that of the numeric representation
   of "NA".
 */
int plink_annot_string_quant_phen(int line, pheno_rec *locus,
				  pheno_pedrec_data *pedrec, const char *quantstr)
{
    char *endptr;
    double quant;
    int err = 0;

    if (strcasecmp(quantstr, "NA") == 0)
        quant = QMISSING;
    else if (inSet(quantstr, Input->MissingCodesSet)) {
        quant = QMISSING;
    } else {
        /* try to read in the values */
        quant = strtod(quantstr, &endptr);
        if (endptr == quantstr || (*endptr)) {
            /* conversion failed */
            quant=QUNDEF;
            SUPPRESS_MSSG_NESTED(illegal_quant);
            errorvf("Entry %d: %s Invalid quantitative phenotype for locus %s\n",
                    line, quantstr, locus->Name);
            err = 1;
        } else if (quant == QMISSING || quant == QUNDEF) {
            errorvf("Line %d, trait %s, value %s:\n", line, locus->Name, quantstr);
            errorvf("Internal Quantitative Missing Value Consistency Error.  Get Help.\n");
            EXIT(OUTOF_BOUNDS_ERROR);
        }
    }

    if (PLINK.missing_pheno) {
//      if (quant == PLINK.pheno_value) quant = QMISSING;
        if (fabs(quant - MissingQuant) < EPSILON) quant = QMISSING;
    }
    if (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].items_read) {
//      if (quant == MissingQuant) quant = QMISSING;
        if (fabs(quant - MissingQuant) < EPSILON) quant = QMISSING;

    }

    check_quant_vs_MissingQuants(quant);

    pedrec->Quant = quant;

    return err;
}

int read_quant_phen(FILE *filep, int locusnm, 
                    pheno_rec *locus, void *ventry,
                    record_type rec)

{
    int lch;
    char quantstr[20], *endptr;
    double quant;
    /* only one of the below is valid */
    linkage_ped_rec *entry = (linkage_ped_rec *) ventry;
    pre_makeped_record *pentry = (pre_makeped_record *) ventry;
    annotated_ped_rec *anentry = (annotated_ped_rec *) ventry;

    lch=fcmap(filep, "%s", quantstr);
    if (rec == Annotated && !strcasecmp(quantstr, "NA")) {
        /* set to missing (initialized by set_missing_quant_input() */
        quant = QMISSING;
    } else {
        /* try to read in the values */
        quant = strtod(quantstr, &endptr);
        if (endptr == quantstr || (*endptr)) {
            /* conversion failed */
            quant=QUNDEF;
        } else if (quant == QMISSING || quant == QUNDEF) {
            errorvf("Line %d, trait %s, value %s:\n", anentry->rec_num, locus->Name, quantstr);
            errorvf("Internal Quantitative Missing Value Consistency Error.  Get Help.\n");
            EXIT(OUTOF_BOUNDS_ERROR);
        }
    }

    if (PLINK.missing_pheno) {
//      if (quant == PLINK.pheno_value) quant = QMISSING;
        if (fabs(quant - MissingQuant) < EPSILON) quant = QMISSING;
    }
    if (Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].items_read) {
//      if (quant == MissingQuant) quant = QMISSING;
        if (fabs(quant - MissingQuant) < EPSILON) quant = QMISSING;
    }

    check_quant_vs_MissingQuants(quant);

    switch(rec) {
    case Annotated:
        anentry->pheno[locusnm].Quant = quant;
        break;
    case Premakeped:
    case Raw_premake:
        pentry->pheno[locusnm].Quant = quant;
        break;
    case Postmakeped:
    case Raw_postmake:
        entry->Pheno[locusnm].Quant = quant;
        break;
    }
    return lch;
}


int plink_annot_string_aff_phen(int line, pheno_rec *locus,
                                pheno_pedrec_data *pedrec, const char *cstatus)
{
    int status;
    int err = 0;
    char *endptr = 0;

    pedrec->Affection.Status = UNDEF;
    pedrec->Affection.Class = UNDEF;

    if (!strcmp(cstatus, "NA")) {
        status = 0;
    } else {
        status = -1;
        status = strtol(cstatus, &endptr, 10);
        if (*endptr) {
            SUPPRESS_MSSG_NESTED(FLOAT_AFFECT);
            if (*endptr == '.')
                errorvf("Line %d, %s: Floating point number used for an affection status %s, setting to unknown.\n",
                        line, locus->Name, cstatus);
            else
                errorvf("Line %d, %s: bad number used for an affection status %s, setting to unknown.\n",
                        line, locus->Name, cstatus);
            status = 0;
            err = 1;
        }
    }

    if (status != 0 && status != 1 && status != 2) {
        if (!strcmp(cstatus, Mega2BatchItems[/* 58 */ Value_Missing_Affect_On_Input].value.name))
            status = 0;
        else if (inSet(cstatus, Input->MissingCodesSet)) 
            status = 0;
/*
        if (PLINK.missing_pheno && status == PLINK.pheno_value) {
            ; / * status = 0; * /
*/
        else {
#if 0
            sprintf(err_msg,
                    "Line %d, %s: column %d: Illegal affection status %s, setting to unknown.",
                    line, locus->Name, locus->col_num, cstatus);
            errorf(err_msg);
#else
            SUPPRESS_MSSG_NESTED(illegal_affect);
            err = 1;
            sprintf(err_msg,
                    "Line %d, %s: Illegal affection status %s, setting to unknown.",
                    line, locus->Name, cstatus);
            errorf(err_msg);
#endif
        }
        status = 0;
    }

    pedrec->Affection.Status = status;
    pedrec->Affection.Class  = 1;
    return err;
}

int read_aff_phen(FILE *filep, int locusnm,
		  pheno_rec *locus, void *ventry,
		  record_type rec)
{
    int lch, status, lclass;
    int *estatus=NULL, *eclass=NULL;
    char *endptr = 0;
    char cstatus[128], cclass[128];


    linkage_ped_rec *entry = (linkage_ped_rec *) ventry;
    pre_makeped_record *pentry = (pre_makeped_record *) ventry;
    annotated_ped_rec *anentry = (annotated_ped_rec *) ventry;

    if (rec == Postmakeped) {
        estatus = &(entry->Pheno[locusnm].Affection.Status);
        eclass = &(entry->Pheno[locusnm].Affection.Class);
    } else if (rec == Premakeped) {
        estatus = &(pentry->pheno[locusnm].Affection.Status);
        eclass = &(pentry->pheno[locusnm].Affection.Class);
    } else if (rec == Annotated) {
        estatus = &(anentry->pheno[locusnm].Affection.Status);
        eclass = &(anentry->pheno[locusnm].Affection.Class);
    }

    *estatus = UNDEF;
    *eclass = UNDEF;

    lch = fcmap(filep, "%s", cstatus);
    
    if (rec == Annotated) {
        if (!strcmp(cstatus, "NA")) {
            status = 0;
        } else {
            status = -1;
            status = strtol(cstatus, &endptr, 10);
        }
    } else {
        status = -1;
        status = strtol(cstatus, &endptr, 10);
    }
    if (endptr && *endptr) {
        SUPPRESS_MSSG_NESTED(FLOAT_AFFECT);
        if (*endptr == '.')
            errorvf("Line %d, %s: Floating point number used for an affection status %s, setting to unknown.\n",
                    anentry->rec_num, locus->Name, cstatus);
        else
            errorvf("Line %d, %s: bad number used for an affection status %s, setting to unknown.\n",
                    anentry->rec_num, locus->Name, cstatus);
        status = 0;
    }

    if (FCMAP_NARGS < 1) {
        errorf("Unexpected END-OF-LINE character reading affection status");
        return lch;
    }

    if (PLINK.plink != 0) {
        // We are processing a PLINK file...
        if (PLINK.phenotypeCoding == 0) {
            // The usual values are: 2, affected; 1, unaffected; 0 or -9 missing
            if (status == 2) PLINK.cases++;
            else if (status == 1) PLINK.controls++;
            else if (status == 0 || status == -9) PLINK.missing++;
            // At this point we assume that this is a quantitative trait
            // I didn't see any place where a decimal point was looked for above?!
        } else {
            if (status == 1) PLINK.cases++;
            else if (status == 0) PLINK.controls++;
            else {
                sprintf(err_msg,
                        "Line %d, column %d: Affection status not 0/1 when --1 PLINK option specified?",
                        anentry->rec_num, locus->col_num);
                errorf(err_msg);
            }
        }
    }
    
    if (status != 0 && status != 1 && status != 2) {
        if (rec == Annotated) {
            if (!strcmp(cstatus, Mega2BatchItems[/* 58 */ Value_Missing_Affect_On_Input].value.name))
                status = 0;
/*
            if (PLINK.missing_pheno && status == PLINK.pheno_value) {
                status = 0;
*/
            else {
                sprintf(err_msg,
                        "Line %d, column %d: Illegal affection status %s, setting to unknown.",
                        anentry->rec_num, locus->col_num, cstatus);
                errorf(err_msg);
            }
        } else {
            sprintf(err_msg,
                    "Line %d, column %d: Illegal affection status %s, setting to unknown.",
                    ((rec == Premakeped)? pentry->rec_num : entry->rec_num),
                    locus->col_num, cstatus);
            errorf(err_msg);
        }
        status = 0;
    }

    *estatus=status;

    if (locus->Props.Affection.ClassCnt >= 2) {
        if ((lch == EOF) || (lch == LF)) {
            errorf("Unexpected END-OF-LINE character reading liability class");
            lclass=UNDEF;
        }

        lch = fcmap(filep, "%s", cclass);

        if (rec == Annotated) {
            if (!strcmp(cclass, "NA")) {
                lclass = 0;
            } else {
                lclass = -1;
                sscanf(cclass, "%d", &lclass);
            }
        } else {
            lclass = -1;
            sscanf(cclass, "%d", &lclass);
        }
        if (FCMAP_NARGS < 1) {
            errorf("Unexpected END-OF-LINE character reading liability class");
            return lch;
        }
        if (rec == Annotated) {
            if (lclass < 1 || (locus->Props.Affection.ClassCnt > 1 && lclass > locus->Props.Affection.ClassCnt)) {
                errorvf("Trait \"%s\" liability class \"%s\" out of range.\n",
                        locus->Name, cclass);
                *eclass = UNDEF;
            } else if (locus->Props.Affection.Class[lclass-1].AutoDef   == 1 &&
                       locus->Props.Affection.Class[lclass-1].FemaleDef == 1 &&
                       locus->Props.Affection.Class[lclass-1].MaleDef   == 1) {
                errorvf("Trait \"%s\" liability class \"%s\" penetrance is undefined.\n",
                        locus->Name, cclass);
                *eclass = UNDEF;
            } else {
                *eclass=lclass;
            }
        } else {
            if (lclass < 1 || (locus->Props.Affection.ClassCnt > 1 && lclass > locus->Props.Affection.ClassCnt)) {
                errorvf("Trait \"%s\" Affection class %d out of range.\n",
                        locus->Name, lclass);
                *eclass = UNDEF;
            } else {
                *eclass=lclass;
            }
        }
    } else {
        *eclass=1;
    }
    return lch;
}

static int read_binary_data(FILE *filep, int locusnm, linkage_locus_rec *locus, linkage_ped_rec *entry)
{
    int bf, allele, a1 = -1, a2 = -1, lch = 2;

    for (allele = 0; (allele < locus->AlleleCnt) && (lch != EOF); allele++) {
        lch = fcmap(filep, "%d", &bf);
        if (bf && (FCMAP_NARGS > 0)) {
            if (a1 >= 0) {
                if (a2 >= 0) warnf("illegal binary factor combination.");
                else a2 = allele;
            } else a1 = allele;
        }
    }

    if (((lch == EOF) || (lch == LF)) && (allele < locus->AlleleCnt - 1)) {
        warnf("Unexpected END-OF-FILE character reading binary factors.");
        set_2alleles(entry->Marker, locusnm, locus, UNDEF, UNDEF);
        return lch;
    }

    a1++;  a2++;
    if ((a1 < 0) || (a1 > locus->AlleleCnt)) {
        warnvf("allele %d is out of range\n", a1);
        a1 = UNDEF;
    }
    if (a2 == 0) {
        set_2alleles(entry->Marker, locusnm, locus, a1, a1);
    } else {
        if ((a2 < 0) || (a2 > locus->AlleleCnt)) {
            warnvf("allele %d is out of range\n", a2);
            a2 = UNDEF;
        }
        set_2alleles(entry->Marker, locusnm, locus, a1, a2);
    }
    return lch;
}

/* This function has been changed a lot.
   1) It is now used for reading all types of pedigree files.
   2) Argument rec tells read_numbered_data what type of allele
   data to expect, premakeped, postmakeped and raw or recoded.
   3) Passing in the person record as a generic (void *) that can
   be cast accordingly).
   4) Read in the genotype data according to recor_type into either
   RAlleles or Alleles.
   5) This function is available globally for use by read_pre_makeped()
   6) 3/12/04 - Also reading in the alleles as strings in all cases
   to check for illegal (non-digit characters for recoded alleles).
*/
void annot_ignore_numbered_data(int line, linkage_locus_rec *locus,
                                void *pedrec, int loc)
{
    set_2Ralleles(pedrec, loc, locus, REC_UNKNOWN, REC_UNKNOWN);
    
}

int allele_count = 0;

allele_prop *Allele_Array[ALLELE_ARRAY];

char *canonical_allele(const char *ra)
{
    extern int MARKER_SCHEME;

    char *cra;
    allele_prop *Ara;
    cra = search_allele(ra);
    if (cra == NULL) {
        int len = strlen(ra);
        if (len <= ALL_LEN)
            Ara = CALLOC((size_t) 1, allele_prop);
        else
            Ara = (allele_prop *)CALLOC(((size_t) sizeof(allele_prop)) + len - ALL_LEN , char);
        cra = allele_prop_allele(Ara);
        strcpy(cra, ra);
        add_allele(cra, cra);  /* need key on heap */
        Allele_Array[allele_count] = Ara;
        Ara->idx = allele_count++;

        /*
         * Test on string (i.e. annotated) values
         */
        if (allele_count >= ALLELE_ARRAY) {
            errorvf("This version of Mega2 only supports %d unique alleles.\n", ALLELE_ARRAY);
            EXIT(OUTOF_BOUNDS_ERROR);
        }

        char *endptr;
        int i1 = strtod(cra, &endptr);
        if (endptr == cra || (*endptr)) {
            // non integer;
        } else if (i1 >= ALLELE_ARRAY && MARKER_SCHEME != MARKER_SCHEME_PTR) {
            errorvf("This version of Mega2 only supports alleles with value less than %d.\n", ALLELE_ARRAY);
            EXIT(OUTOF_BOUNDS_ERROR);
        }

/*      printf("ra: %s %p %p %p\n", ra, Ara, cra, &allele2allele_prop_prop(cra)); */

        allele_prop_reset();
    }
    return cra;
}

int read_numbered_data(FILE *filep, int locusnm,
                       linkage_locus_rec *locus,
		       void *ventry,
		       record_type rec, int last_marker)
{

    int a1=0, a2=0, lch, i, c;
    char ra1[ALL_LEN], ra2[ALL_LEN];
    char *cra1 = NULL, *cra2 = NULL; //silly compiler
    char *rap;
    int undef=0;

    /* only one of the below is valid */
    linkage_ped_rec *entry = (linkage_ped_rec *) ventry;
    pre_makeped_record *pentry = (pre_makeped_record *) ventry;
    annotated_ped_rec *anentry = (annotated_ped_rec *) ventry;

/*    lch = fcmap(filep, "%s", ra1);*/
    rap = ra1;
    while (((c = fgetc(filep)) != -1) && (!isalnum(c)) && c != 10) ;
    if (c == -1 || c == 10) {
        lch = 10;
    } while ((*rap++ = (char)c)) {
      c = (int)fgetc(filep);
        if (c == -1 || c == 10) {
            lch = 10;
            c = 0;
        } else if (!isalnum(c)) {
            lch = 32;
            c = 0;
        }
    }
    ungetc(lch, filep);

    if (rec == Premakeped || rec == Postmakeped) {
        for (i=0; i < (int) strlen(ra1); i++) {
            if (!isdigit((int)ra1[i])) {
                errorvf("Illegal (non-numeric) character in genotype for marker %d.\n", locusnm+1);
                undef=1;
            }
        }

        /* else everything is okay */
        char *endptr;
        a1 = strtod(ra1, &endptr);
        if (endptr == ra1 || (*endptr)) {
            undef = 1; // non integer; can not happen here
        } else if (a1 >= ALLELE_ARRAY && MARKER_SCHEME != MARKER_SCHEME_PTR) {
            errorvf("This version of Mega2 only supports alleles with value less than %d.\n", ALLELE_ARRAY);
            EXIT(OUTOF_BOUNDS_ERROR);
        }
    } else if (rec == Raw_premake || rec == Raw_postmake || rec == Annotated) {
/* some day A == a
        for (i=0; i < (int) strlen(ra1); i++)
            ra1[i] = toupper(ra1[i]);
*/
        cra1 = canonical_allele(ra1);
    }

/*    lch = fcmap(filep, "%s", ra2);*/
    rap = ra2;
    while (((c = fgetc(filep)) != -1) && (!isalnum(c)) && c != 10) ;
    if (c == -1 || c == 10) {
        lch = 10;
    } while ((*rap++ = (char)c)) {
      c = (int)fgetc(filep);
        if (c == -1 || c == 10) {
            lch = 10;
            c = 0;
        } else if (!isalnum(c)) {
            lch = 32;
            c = 0;
        }
    }
    ungetc(lch, filep);

    if (rec == Premakeped || rec == Postmakeped) {
        /* check for all digits */
        for (i=0; i < (int) strlen(ra2); i++) {
            if (!isdigit((int)ra2[i])) {
                errorvf("Illegal (non-numeric) character in genotype for marker %d.\n", locusnm+1);
                undef=1;
            }
        }

        char *endptr;
        a2 = strtod(ra2, &endptr);
        if (endptr == ra2 || (*endptr)) {
            undef = 1; // non integer; can not happen here
        } else if (a2 >= ALLELE_ARRAY && MARKER_SCHEME != MARKER_SCHEME_PTR) {
            errorvf("This version of Mega2 only supports alleles with value less than %d.\n", ALLELE_ARRAY);
            EXIT(OUTOF_BOUNDS_ERROR);
        }
    } else if (rec == Raw_premake || rec == Raw_postmake || rec == Annotated) {
/* some day A == a
        for (i=0; i < (int) strlen(ra2); i++)
            ra2[i] = toupper(ra2[i]);
*/
        cra2 = canonical_allele(ra2);
    }

    /*
      if ((((lch == EOF) || (lch == LF)) && (!last_marker)) ||
      (FCMAP_NARGS < 1)) {
      sprintf(err_msg,
      "End of line reading numbered allele 2 for marker %d!",
      locusnm+1);
      errorf(err_msg);
      undef=1;
      }
    */

    if (undef) {
        switch(rec) {
        case Postmakeped:
            set_2alleles(entry->Marker, locusnm, locus, UNDEF, UNDEF);
            break;
        case Premakeped:
            set_2alleles(pentry->marker, locusnm, locus, UNDEF, UNDEF);
            break;
        case Raw_premake:
            set_2Ralleles(pentry->marker, locusnm, locus, REC_UNDEF, REC_UNDEF);
            break;
        case Raw_postmake:
            set_2Ralleles(entry->Marker, locusnm, locus, REC_UNDEF, REC_UNDEF);
            break;
        default:
            break;
        }
    } else {
        switch(rec) {
        case Postmakeped:
            set_2alleles(entry->Marker, locusnm, locus, a1, a2);
            break;
        case Premakeped:
            set_2alleles(pentry->marker, locusnm, locus, a1, a2);
            break;
        case Raw_premake:
            set_2Ralleles(pentry->marker, locusnm, locus, cra1, cra2);
            break;
        case Raw_postmake:
            set_2Ralleles(entry->Marker, locusnm, locus, cra1, cra2);
            break;
        case Annotated:
            set_2Ralleles(anentry->marker, locusnm, locus, cra1, cra2);
            break;

        }
    }
    return lch;
}



int read_premakeped_affec(FILE *filep, int locusnm,
			  pheno_rec *locus,
			  pre_makeped_record *entry)

{
    int lch = read_aff_phen(filep, locusnm, locus,
                            (void *)entry, Premakeped);
    return lch;
}

int read_premakeped_bin(FILE *filep, int locusnm, linkage_locus_rec *locus,
			pre_makeped_record *entry)

{
    /* since I can't make out what the original code is doing,
       here is my version: the sum of n bits must be less than
       or equal to 2
    */

    int allele, a1=0, a2=0, bf, sum, lch = 0;

    sum=0;
    for (allele = 0; allele < locus->AlleleCnt; allele++) {
        lch = fcmap(filep, "%d", &bf);
        if (lch == EOF) {
            errorf("Unexpected END-OF-FILE character reading binary alleles");
            if (a1 > 0 && a2 > 0) {
                set_2alleles(entry->marker, locusnm, locus, a1, a2);
                return lch;
            } else {
                set_2alleles(entry->marker, locusnm, locus, UNDEF, UNDEF);
                return lch;
            }
        }
        if (bf==1 && (FCMAP_NARGS > 0)) {
            sum++;
            if (sum > 2) {
                errorf("illegal binary factor combination, setting to unknown\n");
                a1=0; a2=0;
            } else {
                if (a1==0)
                    a1=allele+1; /* set the first allele to the position */
                else
                    a2=allele+1; /* set the second allele to the current position */
            }
        }
    }
    set_2alleles(entry->marker, locusnm, locus, a1, a2);
    return lch;
}

int read_premakeped_num(FILE *filep, int locusnm,
			linkage_locus_rec *locus1,
			pre_makeped_record *entry)
{

    int a1, a2, lch;

    lch = fcmap(filep, "%d", &a1);
    if ((lch == EOF) || (lch == LF) || (FCMAP_NARGS < 1)) {
        errorf("Unexpected END-OF-LINE reading numbered allele 1");
        set_2alleles(entry->marker, locusnm, locus1, UNDEF, UNDEF);
        return lch;
    }

    lch = fcmap(filep, "%d", &a2);
    if (((lch == EOF) || (lch == LF)) && (FCMAP_NARGS < 1)) {
        errorf("Unexpected END-OF-LINE reading numbered allele 2");
        set_2alleles(entry->marker, locusnm, locus1, UNDEF, UNDEF);
        return lch;
    }

    if (a1 < 0) {
        errorvf("Numbered allele 1 (%d) out of range at locus %s (# %d).\n",
                a1, locus1->Name, locusnm+1);
        set_2alleles(entry->marker, locusnm, locus1, UNDEF, UNDEF);
        return lch;
    }

    if ((locus1->AlleleCnt > 0) && (a1 > locus1->AlleleCnt)) {
        errorvf("Locus %s is listed as having %d alleles in the locus file,\n",
                locus1->Name, locus1->AlleleCnt);
        errorvf("but person %d in pedigree %d has genotype %d/%d.\n",
                entry->indiv, entry->ped, a1, a2);
        set_2alleles(entry->marker, locusnm, locus1, UNDEF, UNDEF);
        return lch;
    }

    if (a2 < 0) {
        errorvf("Numbered allele 2 (%d) out of range at locus %s (# %d).\n",
                a2, locus1->Name, locusnm+1);
        return lch;
    }
    if ((locus1->AlleleCnt > 0) && (a1 > locus1->AlleleCnt)) {
        errorvf("Locus %s is listed as having %d alleles in the locus file,\n",
                locus1->Name, locus1->AlleleCnt);
        errorvf("but person %d in pedigree %d has genotype %d/%d.\n",
                entry->indiv, entry->ped, a1, a2);
        set_2alleles(entry->marker, locusnm, locus1, a1, UNDEF);
        return lch;
    }

    set_2alleles(entry->marker, locusnm, locus1, a1, a2);
    return lch;
}

/*
 * Read one LINKAGE record, returning the number of loci.
 *
 * Return 0 for no error, 1 for expected end of file,
 * -1 for unexpected (fatal) end of file, -2 for other
 * fatal error.
 *
 * The locus top must be intialized and contain the
 * locus types. The lrecdata must be malloc'd.
 */
static int read_linkage_record(FILE *filep, linkage_ped_rec *entry,
                               linkage_locus_top *LTop, lrec_data *lrecdata,
                               int *col2locus, int linenum, char *pedfile,
                               int *origuniq)
{
    int locus;
    int lch=' ';
    char dummy[6][66], rest[400];
    int i, last_marker, col;
    int expected_col=9+LTop->NumPedigreeCols;
    int a1, a2;
    const char *ar1, *ar2;

    fcmap(filep, "%d %d %d %d %d %d %d %d %d",
          &(lrecdata->PedID),
          &(entry->ID),
          &(entry->Father), &(entry->Mother),
          &(entry->First_Offspring),
          &(entry->Next_PA_Sib), &(entry->Next_MA_Sib),
          &(entry->Sex),
          &(lrecdata->Proband));

    entry->rec_num = linenum;
    col=9;  // zero based or 10 if first column is 1.
    /* Store Proband value in OrigProband */
    entry->OrigProband = lrecdata->Proband;
    /* interpret an EOF here as ok but not to be ignored */
    if (feof(filep)) return 1;

    /* check that the data read makes some sense */
    if ((entry->Sex) && (entry->Sex != MALE_ID) &&
        (entry->Sex != FEMALE_ID)) {
        sprintf(err_msg,
                "File %s, Line %d, Col 8: Invalid value (%d) in Sex field.",
                pedfile, linenum, entry->Sex);
        errorf(err_msg);
        return -2;
    }

    if (lrecdata->Proband < 0) {
        sprintf(err_msg,
                "File %s, Line %d, Col 9: Invalid value (%d) in Proband field.",
                pedfile, linenum, lrecdata->Proband);
        errorf(err_msg);
        return -2;
    }

    /* allocate the space for genotypes */
    entry->Pheno  = (LTop->PhenoCnt > 0) ? CALLOC((size_t) LTop->PhenoCnt, pheno_pedrec_data) : 0;
    // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
    entry->Marker = (LTop->MarkerCnt > 0) ? marker_alloc((size_t) LTop->MarkerCnt, LTop->PhenoCnt) : 0;
    /* Now get genotypes, classes, etc. */
    for (;col < expected_col; ) {
        if (col2locus[col] == 0) continue;
        locus = col2locus[col] - 1;
        last_marker = 0;  // not used any more
        lch = fgetc(filep);
        while(isspace(lch)) {
            if ((lch == LF ||  lch == CR) && !last_marker) {
                sprintf(err_msg,
                        "File %s, Unexpected end of line %d",
                        pedfile, linenum);
                errorf(err_msg);
                sprintf(err_msg, "Expecting %d cols, read in %d.",
                        expected_col, col);
                errorf(err_msg);
                return -2;
            }
            lch=fgetc(filep);
        }
        ungetc(lch, filep);

        switch (LTop->Locus[locus].Type) {
        case QUANT:
            lch = read_quantitative_data(filep, locus, &LTop->Pheno[locus], entry);
            if (entry->Pheno[locus].Quant == QUNDEF) {
                sprintf(err_msg,
                        "File %s, Line %d, Ped %d: Invalid quant data for entry %d at locus %s",
                        pedfile, linenum, lrecdata->PedID, entry->ID, LTop->Locus[locus].Name);
                errorf(err_msg);
                return -2;
            }
            col++;
            break;
        case AFFECTION:
            lch = read_affection_data(filep, locus, &(LTop->Pheno[locus]), entry);
            if (entry->Pheno[locus].Affection.Status < 0 || entry->Pheno[locus].Affection.Status > 2) {
                sprintf(err_msg,
                        "File %s, Line %d, Ped %d: entry %d has unknown status (%d) at locus %s",
                        pedfile, linenum, lrecdata->PedID, entry->ID, entry->Pheno[locus].Affection.Status,
                        LTop->Locus[locus].Name);
                errorf(err_msg);
                return -2;
            }
            if (entry->Pheno[locus].Affection.Class == UNDEF) {
                sprintf(err_msg,
                        "File %s, Line %d, Ped %d: entry %d has unknown liability class at locus %s",
                        pedfile, linenum, lrecdata->PedID, entry->ID, LTop->Locus[locus].Name);
                errorf(err_msg);
                return -2;
            }
            /* lwa ...*/
            /*      entry->Orig_status = entry->Pheno[locus].Affection.Status; */
            if (LTop->Pheno[locus].Props.Affection.ClassCnt > 1) {
                col += 2;
            } else {
                col++;
            }
            break;
        case BINARY:
            lch = read_binary_data(filep, locus, &(LTop->Locus[locus]),
                                   entry);
            get_2alleles(entry->Marker, locus, &a1, &a2);
            if ((a1 == UNDEF) || (a2 == UNDEF)) {
                sprintf(err_msg,
                        "File %s, Line %d, Ped %d: Invalid binary data for entry %d at locus %s",
                        pedfile, linenum, lrecdata->PedID, entry->ID,LTop->Locus[locus].Name);
                errorf(err_msg);
                return -2;
            }
            col += 2;
            break;
        case NUMBERED:
        case XLINKED:
        case YLINKED:
            lch = read_numbered_data(filep, locus,
                                     &(LTop->Locus[locus]), (void *) entry,
                                     LTop->PedRecDataType,
                                     last_marker);
            switch(LTop->PedRecDataType) {
                case Postmakeped:
                    get_2alleles(entry->Marker, locus, &a1, &a2);
                    if (a1 == UNDEF || a2 == UNDEF) {
                        sprintf(err_msg,
                                "File %s, Line %d, Ped %d: Invalid numbered data for entry %d at locus %s",
                                pedfile, linenum, lrecdata->PedID, entry->ID,LTop->Locus[locus].Name);
                        errorf(err_msg);
                        return -2;
                    }
                    break;
                case Raw_postmake:
                    get_2Ralleles(entry->Marker, locus, &ar1, &ar2);
                    if ((strcmp(ar1, REC_UNDEF) == 0) ||
                        (strcmp(ar2, REC_UNDEF) == 0)) {
                        sprintf(err_msg,
                                "File %s, Line %d, Ped %d: Invalid numbered data for entry %d at locus %s",
                                pedfile, linenum, lrecdata->PedID, entry->ID,LTop->Locus[locus].Name);
                        errorf(err_msg);
                        return -2;
                    }
                    break;
                default:
                    break;
            }
            col += 2;
            break;
        default:
            errorf("Undefined locus type! Locus data must be read first!");
            return -2;
        }
    }

    for(i=0; i<6; i++) {
        strcpy(dummy[i], " ");
    }

/* This code will break if the extra stuff is too long! */
    if (lch != '\n') {
        (void)fgets(rest, 399, filep);
        i=sscanf(rest, "%s %s %s %s %s %s", dummy[0], dummy[1],
                 dummy[2], dummy[3], dummy[4], dummy[5]);
    }
    else i=0;
    switch(i) {
    case 2:
        /* only ID */
        if (!strcasecmp(dummy[0], "id:")) {
            *origuniq=1;
            strcpy(entry->UniqueID, dummy[1]);
            sprintf(entry->OrigID, "%d", entry->ID);
            sprintf(entry->FamName, "%d", lrecdata->PedID);
        } else {
            *origuniq=0;
        }
        break;
    case 4:
        /* Only ped: and Per: */
        if (!(strcasecmp(dummy[0], "ped:")) && !(strcasecmp(dummy[2], "per:"))) {
            *origuniq=2;
            strcpy(entry->OrigID, dummy[3]);
            strcpy(entry->FamName, dummy[1]);
        } else {
            *origuniq=0;
            sprintf(entry->OrigID, "%d", entry->ID);
            sprintf(entry->FamName, "%d", lrecdata->PedID);
        }
        sprintf(entry->UniqueID, "%d_%d", lrecdata->PedID, entry->ID);
        break;
    case 6:
        if ((strcasecmp(dummy[0], "ped:")==0) && (strcasecmp(dummy[2], "per:")==0) &&
            (strcasecmp(dummy[4], "id:")==0)) {
            *origuniq=3;
            strcpy(entry->UniqueID, dummy[5]);
            strcpy(entry->OrigID, dummy[3]);
            strcpy(entry->FamName, dummy[1]);
        } else if (!(strcasecmp(dummy[2], "ped:")) && !(strcasecmp(dummy[4], "per:")) &&
                 !(strcasecmp(dummy[0], "id:"))) {
            *origuniq=3;
            strcpy(entry->OrigID, dummy[5]);
            strcpy(entry->FamName, dummy[3]);
            strcpy(entry->UniqueID, dummy[1]);
        } else {
            *origuniq=0;
            sprintf(entry->UniqueID, "%d_%d", lrecdata->PedID, entry->ID);
            sprintf(entry->OrigID, "%d", entry->ID);
            sprintf(entry->FamName, "%d", lrecdata->PedID);
        }
        break;
    default:
        *origuniq=0;
        sprintf(entry->UniqueID, "%d_%d", lrecdata->PedID, entry->ID);
        sprintf(entry->OrigID, "%d", entry->ID);
        sprintf(entry->FamName, "%d", lrecdata->PedID);
        break;
    }
    return 0;
}


static int read_linkage_list_insert_fun(list_data EntryData, list_data NewData, int flag)
{
    if ((EntryData == NULL)
        || (((linkage_ped_rec *)EntryData)->ID
            > ((linkage_ped_rec *)NewData)->ID))
        return LIST_INSERT_BEFORE_MODE;
    return LIST_DONT_INSERT_YET;
}


/*
 * Read the pedigree file. Create the pedigree data structures anew.
 * Return a pointer to a new pedigree top or NULL on error.
 *
 * The locus information must be read first.
 *
 * Assumes that all the entries in one pedigree are among all the
 * others in that pedigree (ie. no entries for ped B among entries of ped A).
 * modified: Aug 7, 2001: err is now the number of fields at the
 end of the  linkage record.

*/
linkage_ped_top *read_linkage_ped_file(FILE *filep,
				       linkage_locus_top *LTop,
                                       int *col2locus)
{
    linkage_ped_top *Top;
    list *PedList, *EntryList;
    register linkage_ped_rec *NewEntry, *NewEntry2;
    register linkage_ped_tree *NewPed=NULL;
    int LastPedID = -9;
    lrec_data lrecdata;
    int EntryCnt, PedCnt, linenum=1;
    int entry, ped;
    int err;
    linkage_loop_rec *loop;
    int i;
    int origuniq = 0;
    SUPPRESS_MSSG_NESTED_INIT(untyped);
    int untyped = 0, totaltyped = 0;

    if (filep == NULL) return NULL;

    Top = new_lpedtop();
    Top->LocusTop = LTop;

    PedList = new_list();
    PedCnt = 0;
    EntryList = new_list();
    EntryCnt = 0;

    err = 0;

    while (err >= 0) {
        NewEntry = new_lpedrec();
        err = read_linkage_record(filep, NewEntry, LTop, &lrecdata,
                                  col2locus,
                                  linenum,
                                  mega2_input_files[0],
                                  &origuniq);
        if (err < 0) break;

        if (err != 1) {
            totaltyped++;
            NewEntry->genocnt = crunch_notype(&NewEntry->Marker, LTop);
            if (NewEntry->genocnt == 0) {
/*
                SUPPRESS_MSSG_NESTED(untyped);
                warnvf("Untyped person %4d linenum %d: person %s/%d, fa %d, ma %d\n",
                       untyped, linenum, NewEntry->FamName, NewEntry->ID, NewEntry->Father, NewEntry->Mother);
*/
                untyped++;
            }
        }

        linenum++;
        if (PedCnt == 0) {
            PedCnt++;
            NewPed = new_lpedtree();
            append_to_list_tail(PedList, (list_data) NewPed);
            LastPedID = lrecdata.PedID;
            NewPed->Num = lrecdata.PedID;
            strcpy(NewPed->Name, NewEntry->FamName);
        }
        if ((LastPedID == lrecdata.PedID) && (err == 0)) {
            if (lrecdata.Proband >= 1) {
                if (lrecdata.Proband == 1)
                    NewPed->Proband = NewEntry->ID;
                loop = get_llooprec(&(NewPed->Loops), lrecdata.Proband - 1);
                /* check that the first n entries are defined */
                for (i=0; i<loop->numEntry; i++) {
                    if (loop->Entry[i]==UNDEF) {
                        errorf("INTERNAL error setting loop ids.");
                        EXIT(SYSTEM_ERROR);
                    }
                }
                if (loop->numEntry >= MAXLOOPMEMBERS) {
                    errorf("Exceeded loop limit for an individual\n");
                    EXIT(OUTOF_BOUNDS_ERROR);
                }
                loop->Entry[loop->numEntry] = NewEntry->ID;
                loop->numEntry++;
            }
            EntryCnt++;
            if (list_insert_when_told(EntryList, (list_data) NewEntry,
                                      read_linkage_list_insert_fun)
                == NULL) {
                warnf("internal error inserting member into ordered list.");
                append_to_list_tail(EntryList, (list_data) NewEntry);
            }
        } else {
            /* New pedigree or end of data */
            NewPed->EntryCnt = EntryCnt;
            NewPed->Entry = CALLOC((size_t) EntryCnt, linkage_ped_rec);
            /* There should be a better way... */
            for (entry = 0; entry < EntryCnt; entry++) {
                NewEntry2 = (linkage_ped_rec *) pop_first_list_entry(EntryList);
                /* NM- don't think this is necessary can't see anyplace this
                   is being used as an index */
                /*
                  if (NewEntry2->ID != entry+1) {
                  fprintf(stderr,
		  "WARNING: newly read linkage entry ID number %d is not its ordinal %d\n",
		  NewEntry2->ID, entry+1);
                  fprintf(stderr, "Changed it. (This may cause a problem with loops.)\n");
                  }
                */
                /* code here similar to makenucs */

                copy_lpedrec(NewEntry2, &(NewPed->Entry[entry]), Top);
                NewEntry2->Pheno = NULL;
                NewEntry2->Marker = NULL;
                free(NewEntry2);
            }
            /* EntryList is now empty */
            if (err == 0) {
                EntryCnt = 1;   /* last entry read belongs to this pedigree */
                append_to_list_tail(EntryList, (list_data) NewEntry);
                PedCnt++;
                NewPed = new_lpedtree();
                LastPedID = lrecdata.PedID;
                NewPed->Num = lrecdata.PedID;
                strcpy(NewPed->Name, NewEntry->FamName);
                if (lrecdata.Proband >= 1) {
                    if (lrecdata.Proband == 1) {
                        /* could be only the proband or
                           a loop-breaker proband */
                        NewPed->Proband = NewEntry->ID;
                    } else {
                        HasLoops=1;
                    }
                    loop = get_llooprec(&(NewPed->Loops), lrecdata.Proband - 1);
                    loop->Entry[loop->numEntry] = NewEntry->ID;
                    loop->numEntry++;
                }
                append_to_list_tail(PedList, (list_data) NewPed);
            } else  {
                err = -3;
                free(NewEntry);
            }
        }
    }

    SUPPRESS_MSSG_NESTED_FORCE(untyped);
    if (untyped > 0)
        warnvf("Individuals untyped: %d out of %d\n", untyped, totaltyped);
    SUPPRESS_MSSG_NESTED_FINI(untyped);

    if ((err == -1) || (err == -2)) {
        errorvf("fatal error reading linkage record!\n");
        EXIT(INPUT_DATA_ERROR);
    }

    count_Missing_Quant_consistency();

    Top->PedCnt = PedCnt;
    Top->Ped = CALLOC((size_t) PedCnt, linkage_ped_tree);
    for (ped = 0; ped < PedCnt; ped++) {
        NewPed = (linkage_ped_tree *) pop_first_list_entry(PedList);
        copy_lpedtree1(NewPed, &(Top->Ped[ped]));
        malloc_lpedtree_entries(&(Top->Ped[ped]), (size_t) NewPed->EntryCnt, (size_t) LTop->LocusCnt);
        for(entry=0; entry < NewPed->EntryCnt; entry++) {
            copy_lpedrec(&(NewPed->Entry[entry]), &(Top->Ped[ped].Entry[entry]), Top);
            NewPed->Entry[entry].Pheno  = NULL;
            NewPed->Entry[entry].Marker = NULL;
          /*
            if (Mega2Status == INSIDE_RECODE)
                copy_lpedrec_raw_alleles(&(NewPed->Entry[entry]),
                                         &(Top->Ped[ped].Entry[entry]), Top);
          */
        }
        free_all_from_lpedtree(NewPed, LTop->PhenoCnt); free(NewPed);
    }

    /* Set the OrigIds and UniqueIds fields */
    switch(origuniq) {
    case 3:
        Top->OrigIds = 1;
        Top->UniqueIds = 1;
        break;
    case 2:
        Top->OrigIds = 1;
        Top->UniqueIds = 1;
        break;
    case 1:
        Top->OrigIds = 0;
        Top->UniqueIds = 1;
        break;
    case 0:
        Top->OrigIds = 0;
        Top->UniqueIds = 1;
        break;
    default:
        break;
    }

    Top->pedfile_type=POSTMAKEPED_PFT;
#ifndef HIDESTATUS
    sprintf(err_msg, "Data read in from pedigree file:");
    mssgf(err_msg);

    if (Top->OrigIds) {
        sprintf(err_msg, "All pedigree records have \"Ped\" and \"Per\" fields.");
        mssgf(err_msg);
    }
    if (Top->UniqueIds) {
        sprintf(err_msg, "All pedigree records have \"ID\" fields.");
        mssgf(err_msg);
    }
#endif

    /* PedList is now empty */
    free_list(PedList);
    free_list(EntryList);

    return Top;
}

void count_Missing_Quant_consistency()
{
    SUPPRESS_MSSG_NESTED_FORCE(BAD_QMISSING_OUT);
    if (SUPPRESS_MSSG_COUNT(BAD_QMISSING_OUT)) {
        errorvf("Output Quantitative Missing Value was found among the quantitative data %d times\n",
                SUPPRESS_MSSG_COUNT(BAD_QMISSING_OUT));
//tmp   EXIT(INPUT_DATA_ERROR);
    }
    SUPPRESS_MSSG_NESTED_FINI(BAD_QMISSING_OUT);
}


static linkage_locus_top *read_linkage_locus_file(FILE *filep, int linkagecols, int **col2locus)
{
    linkage_locus_top *LTop;
    int locus1, allele, tmpi, tmpi2;
    int ltype;
    int *order;
    int *colnm, colnxt;
    int lch = 0;
    /*   int word_length, *order */
    linkage_locus_rec *Locus, *LocusBeg;
    int loc;
    pheno_rec *Pheno;
    marker_rec *Marker;
#define DUMMYSTRG_LENGTH    80
    char dummystrg[DUMMYSTRG_LENGTH], *dummy, marker_name[MAX_NAMELEN];
    char cloc_num[10];
    int c;

    if (filep == NULL) return NULL;

    LinkageLocusFile=1;
    LTop = new_llocustop();
    
    /* 
     This is a completely bogus comment. The first line is read
     but there is never the check described below??? Should there be???
     
     first decide whether this is in a names file format
     by reading the first string, if it is a number
     then we will assume the linkage format, otherwise if it is
     one of M, T, L, Q, C, then it is names-file format (similar
     */
    
    /*
     For a description of the linkage file formet please see...
     http://linkage.rockefeller.edu/soft/linkage/
     2.6 Description of Loci (DATAFILE)
     
     Read the first two lines of the file (throwing away the comments).
     The lines should look something like this...
5 0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM
0 0.0 0.0 0  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)
    */
    fcmap(filep, "%d %d %d %d\n%d %g %g %d\n",
          &(LTop->LocusCnt), &(LTop->RiskLocus), &(LTop->SexLinked), &(LTop->Program),
          &(LTop->MutLocus), &(LTop->MutMale), &(LTop->MutFemale), &(LTop->Haplotype)
          );

    LTop->Locus = CALLOC((size_t) LTop->LocusCnt, linkage_locus_rec);
// over commit for now
    LTop->Pheno  = CALLOC((size_t) LTop->LocusCnt /* LTop->PhenoCnt */, pheno_rec);
    // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
    // LTop->Marker = (CALLOC((size_t) LTop->MarkerCnt, marker_rec)) - LTop->PhenoCnt;
    LTop->Marker = CALLOC((size_t) LTop->LocusCnt /* LTop->MarkerCnt */, marker_rec);

    order = CALLOC((size_t) LTop->LocusCnt, int);
    colnm = CALLOC((size_t) LTop->LocusCnt, int);

    /*
     Read the next line in the file which should be the locus order info which
     consists of the locus numbers separated by white space.
     It should look something linek this...
     1  2 3 4 5
     Check to see that all of the loci are specified, and that they
     are digits. Store the list into the array 'order'
     */
    c=fgetc(filep);
    for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
        clear_llocusrec(&(LTop->Locus[locus1]), TYPE_UNSET);
        sprintf(cloc_num, " ");
        while(isspace(c) || c == LF || c == CR) {
            if (locus1 > 0 &&(c == LF || c == CR)) {
                errorf("Locus order (line number 3) is too short");
                errorvf("Read in %d locus numbers, %d expected.\n",
                        locus1, LTop->LocusCnt);
                EXIT(INPUT_DATA_ERROR);
            }
            c=fgetc(filep);
        }
        if (isdigit(c)) {
            lch=0;
            while(isdigit(c)) {
	        cloc_num[lch++] = (char)c;
                c=fgetc(filep);
            }
            cloc_num[lch++]='\0';
            order[locus1] = atoi(cloc_num);
        }
    }
    ungetc(c, filep);
    if (c != LF && c != CR) {
        /* REPLACE WITH SKIP_LINE */
        fcmap(filep, "%=\n", lch);
    }
    
    /*
     The remaining "blocks" of lines in the file describe the locus
     */
    int pheno1 = 0, marker1 = 0;
    colnxt = linkagecols;
    for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
        Locus = &(LTop->Locus[locus1]);
        Locus->number = 1;
        colnm[locus1] = colnxt;
        Pheno = &(LTop->Pheno[pheno1]);
        Marker = &(LTop->Marker[marker1]);
        /*
         The first digit of the line describes the locus type using a numeric code:
         0  =      Quantitative variable
         1  =      Affection status
         2  =      Binary factors
         3  =      Numbered alleles
         
         The next digit is the total number of Alleles.
         
         The following are some example lines...
1   2   # TRAIT
0   2  # Q1
3   2   # M1
3   2   # M2
3   2   # M3
         */
        lch = fcmap(filep, "%d %d", &ltype, &(Locus->AlleleCnt));
        
        /*
         The remainder of the line gives the name of the locus following a '#' character...
         */
        strcpy(dummystrg, "");
        (void)fgets(dummystrg, DUMMYSTRG_LENGTH - 1, filep);
        // Throw away the white space leading up to the '#' character...
        dummy = &(dummystrg[0]);
        while(isspace((int) *dummy) && dummy != NULL){
            dummy++;
        }
        if (*dummy == '#'){
            dummy++;
            // Throw away the white space between the '#' character and the name...
            for(;(isspace((int) *dummy) && *dummy );dummy++);
        } else {
            errorf("In locus file, trying to read locus name, but no # symbol found.");
            errorvf("Locus number     : %d\n", (locus1 + 1));
            errorvf("Number of alleles: %d\n", Locus->AlleleCnt);
            errorvf("String read in   : %s\n", dummystrg);
            EXIT(INPUT_DATA_ERROR);
        }

        if (dummy == NULL) {
            errorf("In locus file, trying to read locus name, but no # symbol found.");
            errorvf("Locus number     : %d\n", (locus1 + 1));
            errorvf("Number of alleles: %d\n", Locus->AlleleCnt);
            errorvf("String read in   : %s\n", dummystrg);
            EXIT(INPUT_DATA_ERROR);
        }

        if (dummystrg[strlen(dummystrg)-1] != '\n') {
            fcmap(filep, "%=\n");
        }
        
        // Read and store the marker...
        sscanf(dummy, "%s", marker_name);
        Locus->Name = strdup(marker_name);
        
        /*
         Up to this point we were simply reading the first line of the locus
         description block. Next we read specific informaiton which is found on
         separate lines, and it get's more complicated. See..
         http://linkage.rockefeller.edu/soft/linkage/
         In the section "Description of Loci"
         */
        
        /*
         Read the gene frequencies which is always the first record in the block e.g.,
0.990000 0.010000   << GENE FREQUENCIES
         */
        Locus->Allele = CALLOC((size_t) Locus->AlleleCnt, linkage_allele_rec);
        for (allele = 0; allele < Locus->AlleleCnt; allele++) {
            clear_lallelerec(&(Locus->Allele[allele]));
            lch = fcmap(filep, "%g", &(Locus->Allele[allele].Frequency));
            /*        if (Locus->Allele[allele].Frequency > 1.0) { */
            /*  	errorf */
        }
        fcmap(filep, "%=\n", lch);               /* skip any comments     */
        /*
         The remaining lines in the block are dependent on the locus type...
         */
        switch(ltype) {
        case 0:     /* Quantitative Variable */
            Locus->Type = QUANT;
            Locus->Class = TRAIT;
            Locus->Pheno = Pheno;
            Locus->Marker = 0;
            Pheno->Name = Locus->Name;
            pheno1++;
            clear_llocusdata(Locus, QUANT);
            /* Read the number of traits. e.g.,
1 << NO. OF TRAITS
             */
            fcmap(filep, "%d\n", &tmpi);
            Locus->Pheno->Props.Quant.ClassCnt = tmpi;
            if (Locus->Pheno->Props.Quant.ClassCnt != 1) {
                if (Locus->Pheno->Props.Quant.ClassCnt == 0) {
                    sprintf(err_msg, "Number of liability classes should not be zero,");
                    warnf(err_msg);
                    sprintf(err_msg, "Setting number of classes to one.");
                    warnf(err_msg);
                    Locus->Pheno->Props.Quant.ClassCnt = 1;
                } else {
                    errorf("Only one liability class allowed for quantitative loci.");
                    EXIT(INPUT_DATA_ERROR);
                }
            }
            colnxt++;
            /* 
             Read the genotype means. e.g.,
1.000 10.000 20.000 << GENOTYPE MEANS
            */
            Locus->Pheno->Props.Quant.Mean = CALLOC((size_t) tmpi, double *);
            for (; tmpi > 0; tmpi--) {
                double m;
                Locus->Pheno->Props.Quant.Mean[tmpi-1] = CALLOC((size_t) 3, double);
                // read the three values that follow...
                for (tmpi2 = 0; tmpi2 < 2; tmpi2++) {
                    lch = fcmap(filep, "%g", &m);
#ifdef TEST
                    m = 0;
#endif
                    Locus->Pheno->Props.Quant.Mean[tmpi-1][tmpi2] = m;
                }
                lch = fcmap(filep, "%g\n", &m);
#ifdef TEST
                    m = 0;
#endif
                Locus->Pheno->Props.Quant.Mean[tmpi-1][2] = m;
            }
            /* Assume only one trait per quantitative trait locus */
            Locus->Pheno->Props.Quant.Variance = CALLOC((size_t) 1, double *);
            Locus->Pheno->Props.Quant.Variance[0] = CALLOC((size_t) 1, double);
            /*
1.000  << VARIANCE - COVARIANCE MATRIX
             */
            lch = fcmap(filep, "%g\n", &(Locus->Pheno->Props.Quant.Variance[0][0]));
            /*
1.000  << MULTIPLIER FOR VARIANCE IN HETEROZYGOTES
             */
            lch = fcmap(filep, "%g\n", &(Locus->Pheno->Props.Quant.Multiplier));
            LTop->NumPedigreeCols++;
            break;
        case 1:     /* Affection Status */
            Locus->Type = AFFECTION;
            Locus->Class = TRAIT;
            Locus->Pheno = Pheno;
            Locus->Marker = 0;
            Pheno->Name = Locus->Name;
            pheno1++;
            clear_llocusdata(Locus, AFFECTION);
            fcmap(filep, "%d\n", &tmpi);
            if (tmpi == 0) {
                errorvf("Zero liability class found for affection status locus %s\n",
                        Locus->Name);
                EXIT(INPUT_DATA_ERROR);
            }
            Locus->Pheno->Props.Affection.ClassCnt = tmpi;
            Locus->Pheno->Props.Affection.Class = CALLOC((size_t) tmpi, linkage_affection_class);
            Locus->Pheno->Props.Affection.PenCnt = Locus->AlleleCnt * (Locus->AlleleCnt + 1) / 2;

            for (tmpi = 0; tmpi < Locus->Pheno->Props.Affection.ClassCnt; tmpi++) {
                clear_laffclass(&(Locus->Pheno->Props.Affection.Class[tmpi]));

                if (LTop->SexLinked) {
                    Locus->Pheno->Props.Affection.Class[tmpi].FemaleDef = 0;
                    Locus->Pheno->Props.Affection.Class[tmpi].FemalePen =
                        CALLOC((size_t) Locus->Pheno->Props.Affection.PenCnt, double);
                    for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++) {
                        lch = fcmap(filep, "%g",
                                    &(Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[tmpi2]));
                    }
                    fcmap(filep, "%=\n", lch);

                    Locus->Pheno->Props.Affection.Class[tmpi].MaleDef = 0;
                    Locus->Pheno->Props.Affection.Class[tmpi].MalePen =
                        CALLOC((size_t) Locus->AlleleCnt, double);
                    for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++) {
                        lch = fcmap(filep, "%g",
                                    &(Locus->Pheno->Props.Affection.Class[tmpi].MalePen[tmpi2]));
                    }
                    fcmap(filep, "%=\n", lch);

                    // Default Autosome
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoDef = 1;
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoPen =
                        CALLOC((size_t) Locus->Pheno->Props.Affection.PenCnt, double);
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[0] = DOM_G11;
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[1] = DOM_G12;
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[2] = DOM_G22;

                } else {
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoDef = 0;
                    Locus->Pheno->Props.Affection.Class[tmpi].AutoPen =
                        CALLOC((size_t) Locus->Pheno->Props.Affection.PenCnt, double);
                    for (tmpi2 = 0; tmpi2 < Locus->Pheno->Props.Affection.PenCnt; tmpi2++) {
                        lch = fcmap(filep, "%g",
                                    &(Locus->Pheno->Props.Affection.Class[tmpi].AutoPen[tmpi2]));
                    }
                    fcmap(filep, "%=\n", lch);

                    // Default SEX
                    Locus->Pheno->Props.Affection.Class[tmpi].FemaleDef = 1;
                    Locus->Pheno->Props.Affection.Class[tmpi].FemalePen =
                        CALLOC((size_t) Locus->Pheno->Props.Affection.PenCnt, double);
                    Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[0] = DOM_G11;
                    Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[1] = DOM_G12;
                    Locus->Pheno->Props.Affection.Class[tmpi].FemalePen[2] = DOM_G22;

                    Locus->Pheno->Props.Affection.Class[tmpi].MaleDef = 1;
                    Locus->Pheno->Props.Affection.Class[tmpi].MalePen =
                        CALLOC((size_t) Locus->AlleleCnt, double);
                    Locus->Pheno->Props.Affection.Class[tmpi].MalePen[0] = DOM_G11;
                    Locus->Pheno->Props.Affection.Class[tmpi].MalePen[1] = DOM_G22;
                }
            }
            if (Locus->Pheno->Props.Affection.ClassCnt > 1) {
                LTop->NumPedigreeCols += 2;
                colnxt += 2;
            } else {
                LTop->NumPedigreeCols++;
                colnxt++;
            }
            break;
        case 2:     /* Binary Factors */
            Locus->Type = BINARY;
            Locus->Class = MARKER;
            Locus->Marker = Marker;
            Locus->Pheno = 0;
            Marker->Name = Locus->Name;
            marker1++;
            clear_llocusdata(Locus, BINARY);
            fcmap(filep, "%d\n", &tmpi);
            Locus->Marker->Props.Binary.FactorCnt = tmpi;
            Locus->Marker->Props.Binary.Factor = CALLOC((size_t) tmpi, char *);
            for (tmpi = 0; tmpi < Locus->Marker->Props.Binary.FactorCnt; tmpi++) { /* read factors */
                int f;
                Locus->Marker->Props.Binary.Factor[tmpi] = CALLOC((size_t) Locus->AlleleCnt, char);
                for (tmpi2 = 0; tmpi2 < Locus->AlleleCnt; tmpi2++) {
                    lch = fcmap(filep, "%d", &f);
                    Locus->Marker->Props.Binary.Factor[tmpi][tmpi2] = (char) f;
                }
                fcmap(filep, "%=\n", lch);   /* be sure we go past the newline */
            }
            LTop->NumPedigreeCols += 2;
            colnxt += 2;
            break;
        case 3:     /* Numbered Alleles */
            HasMarkers=1;
            Locus->Type = NUMBERED;
            Locus->Class = MARKER;
            Locus->Marker = Marker;
            Locus->Pheno = 0;
            Marker->Name = Locus->Name;
            marker1++;
            clear_llocusdata(Locus, NUMBERED);
            LTop->NumPedigreeCols += 2;
            colnxt += 2;
            break;
        default:
            break;
        }
        if (LTop->RiskLocus == locus1+1)
            fcmap(filep, "%d\n", &(LTop->RiskAllele));
    }

    /* Read the recombination information. see...
     http://linkage.rockefeller.edu/soft/linkage/
     Section "Recombination Information"
     
     Sex Difference:
     0    =    no sex-difference
     1    =    constant sex-difference (the ratio of female/male genetic distance is the same intervals)
     2    =    variable sex-difference (the female/male distance ratio can be different in each interval)
     
     The interference option can take the following values:
     0    =    no interference
     1    =    interference without a mapping function
     2    =    user-specified mapping function
     
     e.g., 
0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)
     */
    fcmap(filep, "%d %d\n",
          &(LTop->SexDiff), &(LTop->Interference));
    /*
     Read the recomninarion values. e.g.,
0.1 0.1 0.1 0.1  << RECOMBINATION VALUES     
     */
    LTop->MaleRecomb = CALLOC((size_t) LTop->LocusCnt - 1, double);
    /*
     Get either the recombination rates between adjacent loci in the two segments, or additionally the recombination rate between the flanking loci.
     */
    for (tmpi = 0; tmpi < LTop->LocusCnt - 1; tmpi++)
        lch = fcmap(filep, "%g", &(LTop->MaleRecomb[tmpi]));
    if (LTop->Interference == NO_MAPPING_FUN)
        lch = fcmap(filep, "%g", &(LTop->FlankRecomb));
    else LTop->FlankRecomb = UNDEF;

    fcmap(filep, "%=\n", lch);

    switch (LTop->SexDiff) {
    case NO_SEX_DIFF:
        /*
         When sex-difference is "0," one recombination rate is given for each of the nlocus-1 segments. e,g,
0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)
0.1 0.1 0.1 0.1  << RECOMBINATION VALUES  
         */
        LTop->FMRatio = UNDEF;
        LTop->FemaleRecomb = NULL;
        break;
    case CONSTANT_SEX_DIFF:
        /*
         When sex-difference is "1," the male recombination rates are given on one line, and the female/male genetic distance is specified on the next line, e.g.,
1  0                << sex difference, interference
0.1  0.2  0.1       << male recombination
2.0                 << female/male ratio of genetic distance
        */
        fcmap(filep, "%g\n", &(LTop->FMRatio));
        LTop->FemaleRecomb = NULL;
        break;
    case VARIABLE_SEX_DIFF:
        /*
         When the sex-difference option is "2", the male recombination rates are followed on the next line by female recombination rates. e.g.,
2  0                << sex difference, interference
0.1  0.2  0.1       << male recombination
0.2  0.1  0.2       << female recombination
        */
        LTop->FemaleRecomb = CALLOC((size_t) LTop->LocusCnt - 1, double);
        for (tmpi = 0; tmpi < LTop->LocusCnt - 1; tmpi++)
            lch = fcmap(filep, "%g", &(LTop->FemaleRecomb[tmpi]));
        fcmap(filep, "%=\n", lch);
        break;
    default:
        warnvf("Unknown sex difference code %d.\nSetting code to NO_SEX_DIFF\n", LTop->SexDiff);
        LTop->SexDiff=NO_SEX_DIFF;
        LTop->FMRatio = UNDEF;
        LTop->FemaleRecomb = NULL;
        break;
    }
    if (LTop->Program == LINKMAP) {
        fcmap(filep, "%d %g %d",
	      &(LTop->Run.linkmap.trait_marker),
              &(LTop->Run.linkmap.stop_theta),
              &(LTop->Run.linkmap.num_evals));
    }

    LTop->Run.linkmap.recomb_frac=NULL; /* leave it to vitesse
                                           to allocate */

    if (LTop->Program == MLINK) {
/* was f f d */
        fcmap(filep, "%g %g %g", &(LTop->Run.mlink.start_theta),
              &(LTop->Run.mlink.increment),
              &(LTop->Run.mlink.stop_theta));
        LTop->Run.mlink.num_evals=0; /* for now assume this will be put in later
                                        inside Vitesse */
    }

    LocusBeg = Locus = CALLOC((size_t) LTop->LocusCnt, linkage_locus_rec);
    *col2locus = CALLOC((size_t) linkagecols + LTop->NumPedigreeCols, int);
    loc = 1; // must subtract 1 to use loc index
    for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
        if (LTop->Locus[locus1].Type == AFFECTION || LTop->Locus[locus1].Type == QUANT) {
            colnxt = colnm[order[locus1]-1];
            (*col2locus)[colnxt] = loc++;
            LTop->Locus[locus1].col_num = colnxt;
            LTop->Locus[locus1].Pheno->col_num = colnxt;
            *Locus++ = LTop->Locus[locus1];
        }
    }
    for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
        if (LTop->Locus[locus1].Type == NUMBERED || LTop->Locus[locus1].Type == BINARY) {
            colnxt = colnm[order[locus1]-1];
            (*col2locus)[colnxt] = loc++;
            LTop->Locus[locus1].col_num = colnxt;
            LTop->Locus[locus1].Marker->col_num = colnxt;
            *Locus++ = LTop->Locus[locus1];
        }
    }
    free(LTop->Locus);
    LTop->Locus = LocusBeg;
// now recommit
    LTop->PhenoCnt = pheno1;
    LTop->Pheno  = REALLOC(LTop->Pheno, (size_t) LTop->PhenoCnt, pheno_rec);
    // NOTE: No space is allocated in marker for entries 0..LTop->PhenoCnt-1
    LTop->MarkerCnt = marker1;
    LTop->Marker = (REALLOC(LTop->Marker, (size_t) LTop->MarkerCnt, marker_rec))  - LTop->PhenoCnt;

    pheno1 = 0;
    marker1 = LTop->PhenoCnt;
    for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
        if (LTop->Locus[locus1].Type == AFFECTION || LTop->Locus[locus1].Type == QUANT) {
            LTop->Locus[locus1].Pheno = &LTop->Pheno[pheno1++];
        } else if (LTop->Locus[locus1].Type == NUMBERED || LTop->Locus[locus1].Type == BINARY) {
            LTop->Locus[locus1].Marker = &LTop->Marker[marker1++];
        }
    }
    free(order);

    printf("Data read in from Locus file:\n");
/*   printf("Total number of loci =  %d, number of marker loci = %d\n", */
/* 	 LTop->LocusCnt, marker_count); */
/*   log_line(mssgf); */
/*   get_trait_list(LTop); */

    return LTop;
}


/*
 * Read the linkage locus, map and pedigree files. Set EXLTop to NULL,
 * since linkage-format files do not allow multiple maps.
 */
extern void clear_YLINKED_females(linkage_ped_top *Top, int raw_allele, int hdr);

linkage_ped_top *read_linkage2(char *pedfl_name, char *locusfl_name,
			       char *omitfl_name, char *mapfl_name,
			       int untyped_ped_opt, analysis_type analysis)

{
    FILE           *pfilep, *lfilep;
    linkage_ped_top *Top;
    int  line_count, ped_count, linkagecols;
    char *pedfile = mega2_input_files[0];
    file_format locus_file_format;
    linkage_locus_top *LTop;
    ext_linkage_locus_top *EXLTop=NULL;
    int *col2locus;

    extern void set_missing_quant_input(linkage_ped_top *Top, analysis_type analysis);

    /* if command line disease locus is present, print it*/
    /* removed by Nandita - assume no command line
       and that files have already been input correctly
       via the input menu */

    if ((pfilep = fopen(pedfl_name, "r")) == NULL) {
        errorvf("Could not open pedigree file %s\n", pedfl_name);
        EXIT(FILE_READ_ERROR);
    }
#ifndef HIDEFILE
    log_line(mssgf); // WHY is this line here????
    sprintf(err_msg,"Checking format of pedigree file %s.", pedfl_name);
    mssgf(err_msg);
#endif
    // Returns the count of distinct pedigrees if this is a validly formatted premakeped
    // file, otherwise 0. Only looks at the first 5 colums of the file...
    ped_count=check_pre_makeped(pfilep, &line_count);
    linkagecols = ped_count > 0 ? 5 : 9;
    // Since check_pre_makeped read the file, we need to rewind to the beginning.
    rewind(pfilep);

    if ((lfilep = fopen(locusfl_name, "r")) == NULL) {
        errorvf("Could not open locus file %s \n", locusfl_name);
        EXIT(FILE_READ_ERROR);
    }
    // check_locus_file_format returns either LINKAGE, or NAMES (annotated)...
    if ((locus_file_format = check_locus_file_format(lfilep)) == NAMES) {
        // Process the annotated file...
        
        LTop = read_marker_data(lfilep, -1,-1);
        int locus1;
        col2locus = CALLOC((size_t) linkagecols + LTop->NumPedigreeCols, int);
        for (locus1 = 0; locus1 < LTop->LocusCnt; locus1++) {
            linkage_locus_rec *Locus = &LTop->Locus[locus1];
            col2locus[Locus->col_num + linkagecols] = locus1 + 1;
        }

        if (mapfl_name != NULL) {
            EXLTop = read_map_file(mapfl_name, LTop);
        }

        NumChromo=get_chromosome_list(LTop, global_chromo_entries,
                                      chromo_loci_count, 0, analysis);
        write_locus_stats(LTop, locus_file_format);

        Top = create_full_marker_data(pfilep, omitfl_name,
                                      &untyped_ped_opt,
                                      LTop, col2locus, analysis);
        Top->EXLTop = EXLTop;
    } else {
        // Process what we hope to be linkage file format...
        
        LTop = read_linkage_locus_file(lfilep, linkagecols, &col2locus);
        if (mapfl_name != NULL) {
            EXLTop = read_map_file(mapfl_name, LTop);
        }

        NumChromo=get_chromosome_list(LTop, global_chromo_entries,
                                      chromo_loci_count, 0, analysis);
        write_locus_stats(LTop, locus_file_format);

        if (ped_count > 0)    {
            pedfile_type=PREMAKEPED_PFT;
            LTop->PedRecDataType = Premakeped;
            if (check_ped_file_cols(5+LTop->NumPedigreeCols, pfilep, pedfile) != 1) {
                errorvf("Data input error while checking pedegree file %s \n", pedfl_name);
                EXIT(INPUT_DATA_ERROR);
            }
            set_missing_quant_input((linkage_ped_top *) NULL, analysis);
            Top =  read_pre_makeped(pfilep, ped_count, line_count, NULL, LTop, col2locus);
            Top->EXLTop = EXLTop;
        } else {
            pedfile_type=POSTMAKEPED_PFT;
            LTop->PedRecDataType = Postmakeped;
            if (check_ped_file_cols(9+LTop->NumPedigreeCols, pfilep, pedfile) != 1) {
	        errorvf("Data input error while checking pedegree file %s \n", pedfl_name);
                EXIT(INPUT_DATA_ERROR);
            }
            set_missing_quant_input((linkage_ped_top *) NULL, analysis);
            Top = read_linkage_ped_file(pfilep, LTop, col2locus);
            SUPPRESS_MSSG_NESTED_FORCE(FLOAT_AFFECT);
            if (SUPPRESS_MSSG_COUNT(FLOAT_AFFECT) > 0) {
                warnvf("There were %d instances of decimal numbers read where affects were expected.\n",
                       SUPPRESS_MSSG_COUNT(FLOAT_AFFECT));
            }
            SUPPRESS_MSSG_NESTED_FINI(FLOAT_AFFECT);
            Top->EXLTop = EXLTop;
        }
        log_line(mssgf);
        mssgf("Input pedigree data contains:");

        order_heterozygous_allele(Top);
        write_ped_stats(Top, pedfile_type);
        if (omitfl_name != NULL) {
            premakeped_omit_file(Top, omitfl_name, 0);
        }
        if (UntypedPeds == NULL) {
            UntypedPeds = CALLOC((size_t) Top->PedCnt, int);
        }
        omit_peds(untyped_ped_opt, Top);

        clear_YLINKED_females(Top, 0, 0);

        if (Top->PedCnt > 0) {
            log_line(mssgf);
#ifndef HIDESTATUS
            if (omitfl_name != NULL) {
                mssgf("After reading in omit file and excluding untyped pedigrees :");
            } else {
                mssgf("After excluding untyped pedigrees :");
            }
#endif
            write_ped_stats(Top, pedfile_type);
        }
    }

    fclose(pfilep); fclose(lfilep);

    if (Top == NULL) {
        errorvf("Reading pedigree file.\n");
        EXIT(DATA_INCONSISTENCY);
    }
    /* removed some code from here, should be in the read_locus_file function */

    // set this so that the omit_peds() calls inside the options files do not have any effect
    untyped_ped_opt=2;
    
    return Top;
}


/**
   Append the "pedigree identification string" to 'err_msg'.
 */
static void other_pedid(const int pedid, linkage_ped_top *Top)
{
    if (Top->pedfile_type == POSTMAKEPED_PFT && Top->OrigIds) {
        grow(err_msg, " (Ped: %s)", Top->Ped[pedid].Name);
    }
}

#ifdef defunct
/**
   Append the "appropriate person identification string" to 'err_msg'.
 */
static void other_perids(const int pedid, const int perid,
                         linkage_ped_top *Top)
{
    if (Top->pedfile_type == POSTMAKEPED_PFT) {
        if (Top->OrigIds && Top->UniqueIds) {
            grow(err_msg, " (Per: %s, ID: %s)",
                 Top->Ped[pedid].Entry[perid].OrigID,
                 Top->Ped[pedid].Entry[perid].UniqueID);
        } else if (Top->OrigIds) {
            grow(err_msg, " (Per: %s)",
                 Top->Ped[pedid].Entry[perid].OrigID);
        } else if (Top->UniqueIds) {
            grow(err_msg, " (ID: %s)",
                 Top->Ped[pedid].Entry[perid].UniqueID);
        }
    } else {
        if (Top->UniqueIds) {
            grow(err_msg, " (ID: %s)",
                 Top->PTop[pedid].persons[perid].uniqueid);
        }
    }
}
#endif

/**
 @brief Untype the locus whether it be a Genotype or Phenotype.
 
 Added by Nandita - 5/20/02 Untyping for trait loci
 */
static void  untype_locus(linkage_locus_rec *locus,
                          pheno_pedrec_data *pdata,
                          void *mdata,
                          int locusnm,
                          const int raw_allele)
{
    linkage_locus_type ltype = locus->Type;
    switch(ltype) {
        case NUMBERED:
        case XLINKED:
        case YLINKED:
        case BINARY:
            // e.g., .Class == MARKER
            if (raw_allele) {
                set_2Ralleles(mdata, locusnm, locus, REC_UNKNOWN, REC_UNKNOWN);
            } else {
                set_2alleles(mdata, locusnm, locus, 0, 0);
            }
            break;
        case AFFECTION:
            pdata->Affection.Status = 0;
            break;
            
        case QUANT:
//          pdata->Quant = QMISSING;
            pdata->Quant = MissingQuant;
            break;
            
        default:
            // This is an internal error...
            break;
    }
}

/**
 @brief Untype Loci for both premakeped and postmaked formats.
 @see uptyped_locus()
 
 There is some inefficiency in not pulling out the 'Top->P...[person_i]' part while
 looping over locus. However, with the many looping situations now required, the code
 reads more cleanly.
 
 @param[in,out] Top        The pedigree information @see http://watson.hgen.pitt.edu/doxygen/struct__linkage__ped__top.html
 @parem[in]     ped_i      Is the index of the pedigree.
 @parem[in]     per_i      Is the index of the person.
 @parem[in]     locus_i    Is the index of the locus.
 @param[in]     raw_allele Is a flag, @see untype_locus()
 @return void
 */
static void untype_locus_for_pedfile_type(linkage_ped_top *Top,
                                          const int ped_i,
                                          const int person_i,
                                          const int locus_i,
                                          const int raw_allele)
{
    if (Top->pedfile_type == PREMAKEPED_PFT) {
        untype_locus(&Top->LocusTop->Locus[locus_i],
                     &(Top->PTop[ped_i].persons[person_i].pheno[locus_i]),
                     Top->PTop[ped_i].persons[person_i].marker,
                     locus_i,
                     raw_allele);
    } else {
        untype_locus(&Top->LocusTop->Locus[locus_i],
                     &(Top->Ped[ped_i].Entry[person_i].Pheno[locus_i]),
                     Top->Ped[ped_i].Entry[person_i].Marker,
                     locus_i,
                     raw_allele);
    }
}

/**
 @brief For a given loci 'loci_i' this routine will iterate over pedigrees and/or persons as appropriate
 depending on the values of omitped, and omitper (see param definitions below).
 
 @param[in,out] Top        The pedigree information @see http://watson.hgen.pitt.edu/doxygen/struct__linkage__ped__top.html
 @param[in]     omitped    If '0' then all pedigrees have been chosen (will iterate over pedigrees).
 @parem[in]     omitped_i  Is the index of the pedigree if omitped is != 0.
 @param[in]     omitper    If '0' then all persons have been chosen (will iterate over persons).
 @parem[in]     omitper_i  Is the index of the persons if omitped and omitper are != 0.
 @parem[in]     loci_i     Is the index of the locus to untype.
 @param[in]     raw_allele Is a flag, see routing untype_locus();
 @param[in,out] displayed_messages A pointer to the count of the number of error messages written.
 @return void
 */
static void untype_locus_for_pedfile_type_ped_per(linkage_ped_top *Top,
                                                  const int omitped,
                                                  const int omitped_i,
                                                  const int omitper,
                                                  const int omitper_i,
                                                  const int loci_i,
                                                  const int raw_allele,
                                                  int *displayed_messages)
{
    if (omitper == 0 && omitped == 0) { // All Ped; All Per
        int ped_i;
        for (ped_i=0; ped_i < Top->PedCnt; ped_i++) {
            // The number of persons is different for each pedigree...
            int num_pers = (Top->pedfile_type == PREMAKEPED_PFT ?
                            Top->PTop[ped_i].num_persons : Top->Ped[ped_i].EntryCnt);
            int per_i;
            for (per_i=0; per_i < num_pers; per_i++)
                untype_locus_for_pedfile_type(Top, ped_i, per_i, loci_i, raw_allele);
        }
        
    } else if (omitper == 0 && omitped != 0) { // Ped; All Per
        int num_pers = (Top->pedfile_type == PREMAKEPED_PFT ?
                        Top->PTop[omitped_i].num_persons : Top->Ped[omitped_i].EntryCnt);
        int per_i;
        for (per_i=0; per_i < num_pers; per_i++)
            untype_locus_for_pedfile_type(Top, omitped_i, per_i, loci_i, raw_allele);
        
    } else if (omitper != 0 && omitped == 0) { // All Ped; Per
        int ped_i;
        for (ped_i=0; ped_i < Top->PedCnt; ped_i++) {
            int num_pers = (Top->pedfile_type == PREMAKEPED_PFT ?
                            Top->PTop[ped_i].num_persons : Top->Ped[ped_i].EntryCnt);
            // look for this person number in the current pedigree....
            int i, per_i=-1;
            for (i=0; i < num_pers; i++) {
                int current_per = (Top->pedfile_type == PREMAKEPED_PFT ?
                                   Top->PTop[ped_i].persons[i].indiv : Top->Ped[ped_i].Entry[i].ID);
                if (current_per == omitper) { per_i=i; break; }
            }
            // See if the person number was found in this pedigree...
            if (per_i != -1)
                untype_locus_for_pedfile_type(Top, ped_i, per_i, loci_i, raw_allele);
            else {
                // Because of the double loop, it would be time consuming to check for this sooner.
                // It would also make this option less useful if there were a few exceptions and we errored
                // out here since there is no method to enumerate the exclusions...
                SUPPRESS_MSSG(*displayed_messages);
                if (*displayed_messages > MAX_PED_ERRORS) Display_Errors = 0;
                sprintf(err_msg, "Could not locate following individual for untyping in");
                other_pedid(ped_i, Top);
                warnf(err_msg);
                (*displayed_messages)++;
            }
        }
        
    } else { // Ped Per
        untype_locus_for_pedfile_type(Top, omitped_i, omitper_i, loci_i, raw_allele);
    }
}

/**
 @brief Omit file data processing support.
 
 This function is called just after the header is read from an annotated file,
 or upon opening a linkage (no header). The annotated preprocessing code needs
 to determine the order of <ped, per, loci>. For linkage the order of
 <ped, per, loci> is assumed to be in that order (e.g., <0, 1, 2>).
 
 The function will process each line in the file in turn untyping the loci
 according to the table below. The code will handel all eight combinations
 of the three variables.
 
 A comment may be introduced with the '#' character. Also, [for historical
 resaons] all text after the three required data itmes will be ignored even
 without the user of the comment character.
 
 @param[in,out] Top         The pedigree information @see http://watson.hgen.pitt.edu/doxygen/struct__linkage__ped__top.html
 @param[in]     omitfl_name Is the name of the omit file to read from.
 @param[in]     omitfp      Is a file pointer opened to the data area. For linkage, this is the beginning of the file.
                            For annotated, this is the line after the header.
 @param[in]     ped_i       Index of the pedigree column in the input file (0-2).
 @param[in]     per_i       Index of the person column in the input file (0-2). 
 @param[in]     loci_i      Index of the loci column in the input file (0-2).
 @param[in]     raw_allele  Is a flag, see routing untype_locus()
 @return void
 
 ============================================================
 All eight possible combinations possible are handled. These are defined as follows:
 
 Ped Per Loci   # * Untyped the individual in the pedigree at the loci.
 
 Ped All Loci   # * Untyped all individuals in the pedigree at the loci.
 
 All Per Loci   # ** Untyped the individual in all pedigrees at the loci.
 
 All All Loci   # Untyped all individuals in all pedigrees at the loci.
 
 Ped Per All    # Untyped the individual in the pedigree at all MARKER loci.
 
 Ped All All    # * Untyped all individuals in the pedigree at all MARKER loci.
 
 All Per All    # * ** Untyped the individual in all pedigrees at all MARKER loci.
 
 All All All    # Untyped all individuals in all pedigrees at all MARKER loci.
                # This will generate a warning that you are deleting all MARKER loci.
 
 @note (Loci == All) applies ONLY to Markers, Traits are not affected.
 
 * These combinations existed prior to 2/4/2013 and the refactoring of this code
 based on the feature enhancement number 354.
 
 ** A warning will be issued if no such individual exists for the pedigree.
 */
void omit_file_data_processing(linkage_ped_top *Top,
                               const char *omitfl_name, FILE *omitfp,
                               const int ped_i, const int per_i, const int loci_i,
                               const int raw_allele)
{
    int displayed_messages1 = 0, displayed_messages = 0, omitfl_lineno=0;
    int checked_for_locus_named_All = 0;
    extern char *unique_id_per(char *reslt, char *ped);

    if (ped_i < 0 || ped_i > 2 || per_i < 0 || per_i > 2 || loci_i < 0 || loci_i > 2) {
        errorvf("INTERNAL: Omit file (%s) line (%d) ped/per/loci index out of bounds.\n",
                omitfl_name, omitfl_lineno);
        EXIT(INPUT_DATA_ERROR);
    }
    
    do {
        char line[3][MAX_NAMELEN], *omitped_str, *omitper_str, *omitloci, *comment_p;
        int omitped, omitped_i, omitper, omitper_i, omitloci_i;
        int num_read, lch = 0;
        
        omitfl_lineno++; // keeping track of for error reporting...
        
        // Step 1: read a line (ped, person, locus) and check for well formed input.
        // 'ped' and 'person' can be either "All" or a number greater than or equal to zero.
        // 'locus' can be either "All" or a locus name string.
        
        // If this is an annotated file, we don't know the order of ped, per, loci...
        num_read=fscanf(omitfp, "%s %s %s", line[0], line[1], line[2]);
        if (feof(omitfp)) break;
        fcmap(omitfp, "%=\n", lch);  // Skip all remaining data on this line...
        
        if (line[0][0] == '#') continue; // Skip comment lines...
        if (num_read != 3) {
            errorvf("Omit file (%s) line (%d) three items must be specified: ped, person, locus\n",
                    omitfl_name, omitfl_lineno);
            EXIT(INPUT_DATA_ERROR);
        }
        // Remove the comment character if it was used and not preceeded by a space at the end of the data...
        if ((comment_p=strchr(line[2], '#')) != NULL) *comment_p = '\0';
        
        // Assign ped, per, loci according to the order given in the formal parameter list...
        omitped_str = line[ped_i];
        omitper_str = line[per_i];
        omitloci = line[loci_i];
       
        if (strcasecmp(omitped_str, "All") == 0)
            omitped = 0; // All
        else
            omitped = 1;
        if ( (strcasecmp(omitper_str, "all") == 0) || (strcmp(omitper_str, "0") == 0) )
            // Note we already support the "all" concept on the person column with "0"
            // This just makes "all" more explicit.
            omitper = 0; // All
        else
            omitper = 1;

        // Step 1a: Check for the cases that don't or might not make sense....
        
        if (omitped == 0 && omitper == 0 && strcasecmp(omitloci, "all") == 0) {
            warnvf("Omit file (%s) line (%d) will cause all genotypes and phenotypes to be excluded.\n",
                   omitfl_name, omitfl_lineno);
        }
        
        // Step 2: Now that the input is known to be well formed, see if the input matches any of the
        // data that it should. Only check if the input specifies a sepcific datum.
        //
        // Variables postpended with *_i give the appropriate valid index of that datum.
        
        omitloci_i=-1; // index of the omit loci if not all...
        if (strcasecmp(omitloci, "All") != 0) { // if "All" then don't check for a specific datum...
            int i;
            // Search for the specific loci...
            for (i=0; i < Top->LocusTop->LocusCnt; i++)
                if (!(strcasecmp(Top->LocusTop->Locus[i].Name, omitloci))) { omitloci_i=i;  break; }
            if (omitloci_i == -1) {
                errorvf("Omit file (%s) line (%d) Locus (%s) not found.\n",
                        omitfl_name, omitfl_lineno, omitloci);
                EXIT(DATA_INCONSISTENCY);
            }
        } else if (checked_for_locus_named_All == 0) {
            int i;
            // Is there a Marker named "All" in the input data?
            for (i = Top->LocusTop->PhenoCnt; i < Top->LocusTop->LocusCnt; i++)
                if (strcasecmp(Top->LocusTop->Marker[i].Name, "All") == 0) {

                    errorvf("Omit file (%s) line (%d) a Locus named 'All' exists.\n",
                            omitfl_name, omitfl_lineno);
                    errorf("In this case 'All' cannot be used in the locus column of the omit file.");
                    EXIT(DATA_INCONSISTENCY);
                }
            checked_for_locus_named_All = 1; // only need to perform this check once...
        }
        
        omitped_i=-1; // index of the omit ped if not all...
        if (omitped != 0) { // if "All" then don't check for a specific datum...
            int i;
            // Search for the specific pedigree...
            for (i=0; i < Top->PedCnt; i++) {
                char *this_ped = (Top->pedfile_type == PREMAKEPED_PFT ? Top->PTop[i].Name : Top->Ped[i].Name);
                if (strcmp(this_ped, omitped_str) == 0) { omitped_i=i; break; }
            }
            if (omitped_i == -1) {
                errorvf("Omit file (%s) line (%d) Pedigree id (%d) not found.\n",
                        omitfl_name, omitfl_lineno, omitped);
                EXIT(DATA_INCONSISTENCY);
            }
        }
        
        // Mega2 documentation referenced below is from http://watson.hgen.pitt.edu/docs/mega2_html/mega2.html
        // Section "7.4.6 Omit file"
        
        // From "Omit file "documentation:
        // If the person number is zero, then all marker genotypes will be set to unknown for the entire pedigree.
        omitper_i=-1; // index of the omit person if not all...
                      // With 'All All' we cannot easily check if the person is missing, and so we do not
                      // at this point. We give a warning when processing later if a person is missing from a pedigree.
                      // Here we only check for the case when we are given a specific person and pedigree.
        if (omitper != 0 && omitped != 0) {
            int num_pers = (Top->pedfile_type == PREMAKEPED_PFT ?
                            Top->PTop[omitped_i].num_persons : Top->Ped[omitped_i].EntryCnt);
            int i;
            for (i=0; i < num_pers; i++) {
                char *this_per = (Top->pedfile_type == PREMAKEPED_PFT ?
                                  Top->PTop[omitped_i].persons[i].uniqueid : Top->Ped[omitped_i].Entry[i].UniqueID);
                if (strcmp(this_per, omitper_str) == 0) { omitper_i=i; break; }
                if (strcmp(unique_id_per(this_per, omitped_str), omitper_str) == 0) { 
                    omitper_i=i; break;
                }
            }
            if (omitper_i == -1) {
                errorvf("Omit file (%s) line (%d) Person id (%s) not found in Pedigree (%s).\n",
                        omitfl_name, omitfl_lineno, omitper_str, omitped_str);
                EXIT(DATA_INCONSISTENCY);
            }
        }
        
        // Step 3: Now that we have valid input, tell the user the persons and pedigrees to be untyped...
        
        if (omitper == 0 && omitped == 0) {
            sprintf(err_msg, "Untyped all individuals in all pedigrees");
        } else if (omitper == 0 && omitped != 0) {
            sprintf(err_msg, "Untyped all individuals in pedigree %s", omitped_str);
//          other_pedid(omitped_i, Top);
        } else if (omitper != 0 && omitped == 0) {
            sprintf(err_msg, "Untyped individual %s in all pedigrees", omitper_str);
        } else {
            sprintf(err_msg, "Untyped individual %s in pedigree %s", omitper_str, omitped_str);
//          other_perids(omitped_i, omitper_i, Top);
//          other_pedid(omitped_i, Top);
        }
        
        // Step 4: Iterate over the locus to be untyped...
        
        if (strcasecmp(omitloci, "All") == 0) {
            int loci_i;
            strcat(err_msg, " at all marker loci.");
            SUPPRESS_MSSG(displayed_messages1);
            if (displayed_messages1 > MAX_PED_ERRORS) Display_Messages = 0;
            mssgf(err_msg);
            displayed_messages1++;
            
            // From "Omit file "documentation:
            // "If All is used, the person or pedigree indicated will be untyped at all marker loci.
            // Trait phenotypes will not be set to unknown."
            // "Please note that when the locus column contains the keyword “All”, it still refers
            // to only marker loci, trait loci are left untouched."
            for (loci_i = Top->LocusTop->PhenoCnt; loci_i< Top->LocusTop->LocusCnt; loci_i++)
                if (Top->LocusTop->Locus[loci_i].Class == MARKER)
                    untype_locus_for_pedfile_type_ped_per(Top, omitped, omitped_i, omitper, omitper_i, loci_i,
                                                          raw_allele, &displayed_messages);
        } else { // Just the given loci...
            grow(err_msg, " at locus %s.", omitloci);
            SUPPRESS_MSSG(displayed_messages1);
            if (displayed_messages1 > MAX_PED_ERRORS) Display_Messages = 0;
            mssgf(err_msg);
            displayed_messages1++;
            
            // From "Omit file "documentation:
            // "The omit file can now be used to set trait phenotypes to unknown."
            untype_locus_for_pedfile_type_ped_per(Top, omitped, omitped_i, omitper, omitper_i, omitloci_i,
                                                  raw_allele, &displayed_messages);
        }
        
    } while (!feof(omitfp));
}

/**
 @brief Omit file data processing support for linkage format files.
 
 For linkage files the order of the data on each line is fixed. The assumed order of
 <ped, per, loci> as <0, 1, 2> will be used

 @param[in,out] Top         The pedigree information @see http://watson.hgen.pitt.edu/doxygen/struct__linkage__ped__top.html
 @param[in]     omitfl_name Is the name of the omit file to read from.
 @param[in]     raw_allele  Is a flag, see routing untype_locus();
 @return void
*/
void premakeped_omit_file(linkage_ped_top *Top, const char *omitfl_name, const int raw_allele)
{
    FILE *omitfp;
    
#ifndef HIDEFILE
    msgvf("Reading omit file %s ...\n", omitfl_name);
#endif
    if ((omitfp=fopen(omitfl_name, "r")) == NULL) {
        errorvf("could not open %s for reading!\n", omitfl_name);
        EXIT(FILE_READ_ERROR);
    }

    Display_Messages = Display_Errors = 1;
    omit_file_data_processing(Top, omitfl_name, omitfp, 0, 1, 2, raw_allele);
    Display_Messages = Display_Errors = 1;
    
    fclose(omitfp);
}


// copy data found in LTop to EXLTop
static  ext_linkage_locus_top *make_EXLTop_from_LTop(linkage_locus_top *LTop, int male_col, int female_col)
{
    ext_linkage_locus_top *EXLTop;
    int m;
    // ext_linkage_locus_top is defined in linkage_external.h
    EXLTop = CALLOC((size_t)1, ext_linkage_locus_top); // malloc and clear
    EXLTop->LocusCnt = LTop->LocusCnt;
    EXLTop->MapCnt = 1; /* For different locus maps */
    
    /* Names of each map taken from Map file header */
    EXLTop->MapNames = CALLOC((size_t)1, char*);
    EXLTop->MapNames[0]=strdup("Map"); // KOSAMBI or HALDANE
    
    /* will be map_function[], one for each map */
    EXLTop->map_functions = CALLOC((size_t)1, char);
    EXLTop->map_functions[0] = LTop->map_distance_type; // 'k' or 'h'
    
//cpk    EXLTop->SexMaps = CALLOC((size_t)1, sex_map_types*);
//cpk    EXLTop->SexMaps[0]=CALLOC((size_t) 3, sex_map_types);
    EXLTop->SexMaps = CALLOC((size_t)1, int*);
    EXLTop->SexMaps[0]=CALLOC((size_t) 3, int);
    EXLTop->SexMaps[0][SEX_AVERAGED_MAP] = 1; // required for linkage map file
    EXLTop->SexMaps[0][MALE_SEX_MAP] = (male_col > 0 ? 1 : 0);
    EXLTop->SexMaps[0][FEMALE_SEX_MAP] = (female_col > 0 ? 1 : 0);
    
#ifdef ALL_ZERO_GENETIC_MAP_INVALID
    // a predicate that tells if the map is valid...
    EXLTop->valid_map_p = CALLOC((size_t)1, int *);
    EXLTop->valid_map_p[0] = CALLOC((size_t)sizeof(sex_map_types), int);
    EXLTop->valid_map_p[0][SEX_AVERAGED_MAP] = 1;
    EXLTop->valid_map_p[0][MALE_SEX_MAP] = (male_col > 0 ? 1 : 0);
    EXLTop->valid_map_p[0][FEMALE_SEX_MAP] = (female_col > 0 ? 1 : 0);
#endif /* ALL_ZERO_GENETIC_MAP_INVALID */
    
    if (LTop->MarkerCnt > 0) {
        EXLTop->EXLocus = (CALLOC((size_t) LTop->MarkerCnt, ext_linkage_locus_rec)) - LTop->PhenoCnt;
        for (m = LTop->PhenoCnt; m < LTop->LocusCnt; m++) {
            // Since LTop is linkage, it only supports one map....
            EXLTop->EXLocus[m].positions = CALLOC((size_t)1, double);
            EXLTop->EXLocus[m].pos_male = CALLOC((size_t)1, double);
            EXLTop->EXLocus[m].pos_female = CALLOC((size_t)1, double);
            EXLTop->EXLocus[m].positions[0] = LTop->Marker[m].pos_avg;
            EXLTop->EXLocus[m].pos_female[0] = LTop->Marker[m].pos_female;
            EXLTop->EXLocus[m].pos_male[0] = LTop->Marker[m].pos_male;
        }
    } else 
        EXLTop->EXLocus = 0;
    Display_Errors = 1;

    //
    // Set the map. If it was read in, make sure that it is accruate. If not, then give a default...

    // If specified, 'None' is the only acceptable value, if not then make it 'None'...
    if (ITEM_READ(Value_Base_Pair_Position_Index)) {
        if (base_pair_position_index >= 0) {
            errorvf("A Value_Base_Pair_Position_Index of %d has been specified.\n", base_pair_position_index);
            errorvf("A LINKAGE map file does not have provisions for a physical map.\n");
            EXIT(INPUT_DATA_ERROR);
        }
        base_pair_position_index = -2;
    } else {
        // There is NEVER a physical map with LINKAGE (choose 'None')...
        base_pair_position_index = -2;
        Mega2BatchItems[/* 47 */ Value_Base_Pair_Position_Index].value.option = base_pair_position_index;
        batchf(Value_Base_Pair_Position_Index);
    }
    
    // If specified, '0' is the only acceptable value, if not then make it '0'...
    if (ITEM_READ(Value_Genetic_Distance_Index)) {
        if (genetic_distance_index != 0) {
            errorvf("A Value_Genetic_Distance_Index of %d has been specified.\n", genetic_distance_index);
            errorvf("However, for a LINKAGE map file there is ALWAYS only one\n");
	    errorvf("(genetic) map and by defnition it is the first.\n");
            EXIT(INPUT_DATA_ERROR);
        }
    } else {
        // The LINKAGE map file there is ALWAYS only one (genetic) map & (by defn) would be the first...
        genetic_distance_index = 0;
        Mega2BatchItems[/* 46 */ Value_Genetic_Distance_Index].value.option = genetic_distance_index;
        batchf(Value_Genetic_Distance_Index);
    }
    
    // If it's specified use it if it is valid, if not then make it sex-averaged...
    if (ITEM_READ(Value_Genetic_Distance_SexTypeMap)) {
        // If a map was specified, make sure that there is one in the input...
        if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT &&
            !(male_col > 0 && female_col > 0)) {
            // Here the user specified a sex specific map in the batch file...
            errorvf("A sex-specific Value_Genetic_Distance_SexTypeMap was specified,\n");
            errorvf("but none was available in the input map file.\n");
            EXIT(INPUT_DATA_ERROR);
        } else if (genetic_distance_sex_type_map == FEMALE_GDMT &&
                   !(female_col > 0)) {
            // Here the user specified a female map in the batch file...
            errorvf("A female Value_Genetic_Distance_SexTypeMap was specified,\n");
            errorvf("but none was available in the input map file.\n");
            EXIT(INPUT_DATA_ERROR);
        } else if (genetic_distance_sex_type_map == NONE_AVAILABLE_GDMT ||
                   genetic_distance_sex_type_map == UNKNOWN_GDMT) {
            mssgvf("Using a sex-averated Value_Genetic_Distance_SexTypeMap.\n");
            genetic_distance_sex_type_map = SEX_AVERAGED_GDMT;
        }
    } else {
        // No type was read, We force them to use the sex-averated map even if they have specified male and female
        // genetic positions in the genetic map. This is "by definition".
        genetic_distance_sex_type_map = SEX_AVERAGED_GDMT;
        Mega2BatchItems[/* 48 */ Value_Genetic_Distance_SexTypeMap].value.option = genetic_distance_sex_type_map;
        batchf(Value_Genetic_Distance_SexTypeMap);
    }
                
    return EXLTop;
}

//
// Linkage map files are of the following form. (see
// http://watson.hgen.pitt.edu/docs/mega2_html/mega2.html
// section 7.3.4. Map file)
// CHROMOSOME KOSAMBI NAME [ERROR]
// CHROMOSOME KOSAMBI NAME MALE FEMALE [ERROR]
// The genetic distances are given in cM. The keyword <kosambi|haldane>
// specifies the sex-averaged map. The 4th column [ERROR] can be added
// to the map file specifying mistyping probabilities for each marker.
// These values are utilized within the Genotyping error simulation option.
// Sex-specific markers may also be added. The documentation shows them
// in the above order, but MALE, FEMALE, and ERROR can come in any order
// according to the way that the file is parsed in the code. It appears
// that CHROMOSOME, KOSAMNI, and NAME are fixed to be the first, second,
// and third columns of the file respectively.
// All input is case insensitive.

/* Modified by Nandita - Mar 8, 2002
   to  include reeading in a 4th column of genotyping error
   probabilities, only if the heading has a 4th keyword called
   error (case-insensitive) */
static ext_linkage_locus_top *read_map_file(char *mapfl_name, linkage_locus_top *LTop1)
{
    /* check the map file for the following :
       ERRORS:
       0) Names are not duplicated
       1) All numbered or binary loci have valid chromosome numbers
       2) All numbered or binary loci have non-negative positions
       WARNINGS:
       3) Positions on the same chromosome are monotonic
       4) Loci in map file are in the locus file
       5) Loci in the locus file are in the map file
       Each condition generates a warning except 0)

    */
    /* Note: Reporting errors across mega2 should be standardized */

    FILE     *filep;
    int      num, i, lch=0, found=0, j;
    int      *num_chr_loc, num_read;
    int      is_fatal, is_non_fatal, fatal=0, non_fatal=0;
    double   *loc_pos;
    char     line[2*FILENAME_LENGTH]="", CM[20]="", yesorno[10];
    char     dum[20]="", rest[3][20];
    int      err_col=0, male_col=-1, female_col=-1, num_extras=0;
    int      linenum=0, bad_line, ui=0;
    char     numeric_str[5][20];
    int      displayed_messages = 0;
    // these variables hold information read from the linkage file...
    char     dname[FILENAME_LENGTH];
    double   position, male_pos, female_pos, extras[3];

    Display_Errors=Display_Messages=1;
    MarkerErrorProb=0;
    /* this is the first time we have a chance to figure out the number
       of chromosomes in the data, also store the chromosome numbers
       inside global_chromo_entries */

#ifndef HIDEFILE
    mssgvf("Reading in map file %s.\n", mapfl_name);
#endif
    MaxChromo=0;
    if ((filep = fopen(mapfl_name, "r")) == NULL) {
        errorvf("Could not open map file %s\n", mapfl_name);
        EXIT(FILE_READ_ERROR);
    }
    strcpy(rest[0], "");   strcpy(rest[1], "");   strcpy(rest[2], "");

    (void)fgets(line, 2*FILENAME_LENGTH-1, filep);
    // CHROMOSOME KOSAMBI NAME [MALE FEMALE ERROR]
    if (sscanf(line, "%s %s %s %s %s %s",
               dum, CM, dum, rest[0], rest[1], rest[2]) < 3) {
        errorvf("The header line of the map file is in the wrong format.\n");
        errorvf("It should contain at least three keywords.\n");
        EXIT(INPUT_DATA_ERROR);
    }

    // The last three entries on the line (if they exist) MALE FEMALE ERROR
    // can apparently be in any order...
    for(i=0; i<3; i++) {
        if (strcasecmp(rest[i], "error") == 0) {
            err_col=i+1;
            num_extras++;
            MarkerErrorProb=1;
        } else if (strcasecmp(rest[i], "male") == 0) {
            male_col=i+1;
            num_extras++;
        } else if (strcasecmp(rest[i], "female") == 0) {
            female_col=i+1;
            num_extras++;
        }
    }

    // Read the entire file once to decide the maximum chromosome value.
    // This will be the first value (column) in each data line (row).
    // Also to decide number of unmapped loci (second field == UNKNOWN_CHROMO.
    while (!feof(filep)) {
        num=0;
        (void)fscanf(filep, "%d", &num);
        if (num == UNKNOWN_CHROMO) {
            ui++;
        } else if (num > MaxChromo) MaxChromo=num;
        fcmap(filep, "%=\n", lch);
    }
    /*  printf("largest chromosome number %d\n", MaxChromo); */
    fclose(filep);
    if (!strncasecmp(CM, "kosambi", (size_t)7)) {
        LTop1->map_distance_type = 'k';
    } else {
        LTop1->map_distance_type = 'h';
    }

    for (i = 0; i < LTop1->LocusCnt; i++) {
        LTop1->Locus[i].number = -1; /* Not in map file to begin with */
        if (LTop1->Locus[i].Class == MARKER)
            LTop1->Marker[i].pos_avg = UNKNOWN_POSITION;
    }

    human_x = human_y = human_xy = human_unknown = human_mt = human_auto = 0;

    /* initialize the num_chr_loc */
    loc_pos = CALLOC((size_t) MaxChromo, double);
    num_chr_loc = CALLOC((size_t) MaxChromo, int);
    for (i=0;  i < MaxChromo; i++) {
        loc_pos[i]= UNKNOWN_POSITION; // missing
        num_chr_loc[i]=0;
    }

    if ((filep = fopen(mapfl_name, "r")) == NULL) {
        errorvf("Could not open map file %s\n", mapfl_name);
        EXIT(FILE_READ_ERROR);
    }

    // This time, throw away the header line...
    fcmap(filep, "%=\n", lch); 
    linenum=1;

    do { // For each data line in the file...
        bad_line=0;
        linenum++;
        strcpy(line, "");
        (void)fgets(line, 2*FILENAME_LENGTH-1, filep);
        for (i=0; i < 5; i++) {
            strcpy(numeric_str[i], "");
        }
        // CHROMOSOME KOSAMBI NAME [MALE FEMALE ERROR]
        num_read = sscanf(line, "%s %s %s %s %s %s",
			  numeric_str[0], numeric_str[1],
			  dname,
                          numeric_str[2], numeric_str[3], numeric_str[4]);

/*     num_read = sscanf(line, "%d %f %s %f %f %f", &num, &position, dname, */
/* 		      &extras[0], &extras[1], &extras[2]);  */

        if (num_read <= 0) continue; // skip blank lines

	// The line must have the three required columns plus MAKE,FEMALE,ERROR if specified...
        if (num_read < (num_extras + 3)) {
            sprintf(err_msg, "Line %d: Premature end of line in map file", linenum);
            SUPPRESS_MSSG(displayed_messages);
            if (displayed_messages > MAX_PED_ERRORS) {
                Display_Errors = 0;
            }
            errorf(err_msg);
            displayed_messages++;
            bad_line++;
        }
	
	// This should be the chromosome number...
        if (atoi(numeric_str[0]) < 0) {
            sprintf(err_msg,
                    "Line %d: invalid chromosome number (setting to unknown).",
                    linenum);

            SUPPRESS_MSSG(displayed_messages);
            if (displayed_messages > MAX_PED_ERRORS) {
                Display_Errors = 0;
            }
            errorf(err_msg);
            displayed_messages++;
            non_fatal++;
            num=UNKNOWN_CHROMO; // default value
        } else {
            num=atoi(numeric_str[0]);
        }

	// This should be the sex-averaged map...
        if (atof(numeric_str[1]) < 0 && num != UNKNOWN_CHROMO) {
            sprintf(err_msg,
                    "Line %d: non-numeric or negative value in the sex-averaged position column.",
                    linenum);
            SUPPRESS_MSSG(displayed_messages);
            if (displayed_messages > MAX_PED_ERRORS) {
                Display_Errors = 0;
            }
            warnf(err_msg);
            displayed_messages++;
            position = UNKNOWN_POSITION; // default value
        } else {
            position = atof(numeric_str[1]);
        }
        male_pos=female_pos=UNKNOWN_POSITION; // in case they are not specified

	// Read columns 4-6 if they exist. These would be: MALE, FEMALE, ERROR
	// not necessarly in that order...
        for (i=2; i<= 4; i++) {
            if (atof(numeric_str[i]) < 0 && num!= UNKNOWN_CHROMO) {
                sprintf(err_msg,
                        "Line %d: non-numeric or negative values in",
                        linenum);
                SUPPRESS_MSSG(displayed_messages);
                if (displayed_messages > MAX_PED_ERRORS) {
                    Display_Errors = 0;
                }
                errorf(err_msg);
                errorf("         error probability/sex-specific position columns.");
                displayed_messages++;
                extras[i-2] = UNKNOWN_POSITION;
            } else {
                extras[i-2] = atof(numeric_str[i]);
            }
        }

        if (bad_line > 0) continue; // throw away this line

	// At this point the line is read, and we have information for all
        // of the columns, Now, enter this information for this marker...
        found=0;
	// Find the SNP associated with this data...
        for (i = 0; i < LTop1->LocusCnt; i++) {
            // Determine if we have seen this marker NAME before...
            if (strcasecmp(LTop1->Locus[i].Name, dname) == 0) {
                found=1; // marker was seen before
                /* Check for duplicate marker names */
                if (LTop1->Locus[i].number > -1) {
                    sprintf(err_msg,
                            "Duplicate marker name in %s in map file.",
                            dname);
                    SUPPRESS_MSSG(displayed_messages);
                    if (displayed_messages > MAX_PED_ERRORS) {
                        Display_Errors = 0;
                    }
                    errorf(err_msg);
                    displayed_messages++;
                    fatal++;
                } else {
	            // marker first encountered
                    is_fatal= is_non_fatal = 0;
                    if (LTop1->Locus[i].Class == MARKER)
                        LTop1->Marker[i].chromosome = num;
                    for(j = 1; j <= num_extras; j++) {
                        /* read in all the extras */
                        if (j == err_col && ErrorSimOpt==1) {
                            if (extras[j-1] < 0.0) {
                                non_fatal++;
                            } else {
                                if (LTop1->Locus[i].Class == MARKER)
                                    LTop1->Marker[i].error_prob = extras[j-1];
                            }
                        }
                        if (j==male_col) {
                            male_pos = extras[j-1];
                        }
                        if (j==female_col) {
                            female_pos = extras[j-1];
                        }
                    }

                    /* Set marker type against chromosome number */
                    if (LTop1->SexLinked == 1) {
                        if (num == SEX_CHROMOSOME) {
                            LTop1->Locus[i].Type = XLINKED;

                        } else if (num == MALE_CHROMOSOME) {
                            LTop1->Locus[i].Type = YLINKED;

                        } else {
			    // Yes, but shouldn't we assign the value NUMBERED to .Type here?
                        }
                    }
                    /* Now check map positions against marker type */
                    switch(LTop1->Locus[i].Type) {

                    case NUMBERED:
		      // shouldn't this test be (num >= 1 && num <= 22)?
                        if (num != SEX_CHROMOSOME || LTop1->SexLinked != 1) {
                            /* autosomal, read in all usable positions */
                            if (position < 0.0) {
                                sprintf(err_msg,
                                        "AUTO %s: Negative average map position.",
                                        LTop1->Locus[i].Name);
                                position = male_pos = female_pos = 0;
                                is_fatal++;
                            }
                            break;
                        }

                    case XLINKED:
                        if (female_pos >= 0.0 && position < 0.0 && male_pos < 0.0) {
                            /* correct, no message */
                            ;
                        } else if (female_pos >= 0.0 && position >= 0.0) {
                            sprintf(err_msg,
                                    "XLINKED %s: Average map position will be ignored for X-linked marker.",
                                    LTop1->Locus[i].Name);
                            is_non_fatal++;
                        } else if (position >= 0.0 && female_pos < 0.0) {
                            sprintf(err_msg,
                                    "XLINKED %s: Negative or missing female position (average map will be used).",
                                    LTop1->Locus[i].Name);
                            is_non_fatal++;
                        } else {
                            /* no usable positions */
                            sprintf(err_msg,
                                    "XLINKED %s: No usable map positions.",
                                    LTop1->Locus[i].Name);
                            is_fatal++;
                        }
                        break;

                    case YLINKED:
                        /* if any position is > 0, warn */
                        if (position > 0.0 || male_pos > 0.0 || female_pos > 0.0) {
                            sprintf(err_msg,
                                    "YLINKED %s: Ignoring non-zero map position.",
                                    LTop1->Locus[i].Name);
                            is_non_fatal++;
                        }
                        position = male_pos = female_pos = UNKNOWN_POSITION;
                        break;

                    default:
                        break;
                    }

		    // Also filled in above:
		    // LTop1->SexLinked; LTop1->LocusCnt; LTop1->map_distance_type
                    // LTop1->Marker[i].chromosome; LTop1->Locus[i].Name; LTop1->Locus[i].Type; LTop1->Locus[i].error_prob
                    LTop1->Marker[i].pos_avg = position;
                    LTop1->Marker[i].pos_female = female_pos;
                    LTop1->Marker[i].pos_male = male_pos;
		    // num_chr_loc is an array containing the number of markers found for each chromosome.
		    // So '.number' is a unique number 1-N for a marker in that '.chromosome'.
                    LTop1->Locus[i].number = num_chr_loc[num-1]++;

                    if (is_non_fatal) {
                        SUPPRESS_MSSG(displayed_messages)
                            if (displayed_messages > MAX_PED_ERRORS) {
                                Display_Errors = 0;
                            }
                        warnf(err_msg);
                        non_fatal++;
                        displayed_messages++;
                    } else if (is_fatal) {
                        SUPPRESS_MSSG(displayed_messages)
                            if (displayed_messages > MAX_PED_ERRORS) {
                                Display_Errors = 0;
                            }
                        errorf(err_msg);
                        fatal++;
                        displayed_messages++;
                    }

                    /* store the index of the locus in the map file */
                    if (position >= loc_pos[num-1]) {
                        loc_pos[num-1]=position;
                    } else {
                        non_fatal++;
                        sprintf(err_msg,
                                "Markers not in map order in the map file (locus %s).",
                                dname);

                        SUPPRESS_MSSG(displayed_messages)
                            if (displayed_messages > MAX_PED_ERRORS) {
                                Display_Errors = 0;
                            }
                        warnf(err_msg);
                        displayed_messages++;
                    }
                }
            }
        } // Find the SNP associated with this data...

        if (found == 0) {
            non_fatal++;
            sprintf(err_msg,
                    "Locus %s in the map file is not in the locus file.",
                    dname);
            SUPPRESS_MSSG(displayed_messages)
                if (displayed_messages > MAX_PED_ERRORS) {
                    Display_Errors = 0;
                }
            warnf(err_msg);
            displayed_messages++;
        }
    }  while (!feof(filep)); // For each data line in the file...
    Display_Errors=1;

    fclose(filep);
    free(loc_pos); free(num_chr_loc);

    // Go through the markers and complain if it was never assigned to a chromosome.
    // This could be because it was not found in the map file, or it was a duplicate marker.
    // This same code occurs in annotated_ped_file.cpp:create_entries_for_markers_without_positions().
    displayed_messages = 0;
    for (i = 0; i < LTop1->LocusCnt; i++) {
        if (LTop1->Locus[i].Class == MARKER) {
            if (LTop1->Locus[i].number < 0) {
                non_fatal++;
                sprintf(err_msg, "Locus %s is not in map file.", LTop1->Locus[i].Name);
                SUPPRESS_MSSG(displayed_messages);
                if (displayed_messages > MAX_PED_ERRORS) {
                    Display_Errors = 0;
                }
                warnf(err_msg);
                displayed_messages++;
                LTop1->Marker[i].chromosome = MISSING_CHROMO;
                // no EXLTop entry is created at this point so
                // annotated_ped_file.cpp:new_ext_linkage_locus_top_positions() is not called.
            }
        }
    }
    Display_Errors=1;

    if (fatal > 0 || non_fatal > 0) {
        errorvf("Found %d errors in map file, %d fatal errors, %d non-fatal errors.\n",
                fatal+non_fatal, fatal, non_fatal);
    } else {
#ifndef HIDESTATUS
        mssgf("Data read in from map file.\n");
#endif
    }

    if (fatal) { EXIT(DATA_INCONSISTENCY);   }

    if (non_fatal) {
        if (InputMode || Mega2BatchItems[/* 28 */ Default_Ignore_Nonfatal].items_read == 0) {
            printf("Do you wish to continue with Mega2? (y/n) > ");
            fflush(stdin); fcmap(stdin, "%s", yesorno); newline;
        } else {
            yesorno[0] = Mega2BatchItems[/* 28 */ Default_Ignore_Nonfatal].value.copt;
        }
        if (yesorno[0] == 'n' || yesorno[0] == 'N') { EXIT(EARLY_TERMINATION); }
    }


    global_chromo_entries=CALLOC((size_t) MaxChromo, int);
    chromo_loci_count = CALLOC((size_t) MaxChromo, int);

    NumUnmapped = ui;
    /* NM - not counting missing from map */
    /* Here ui = # missing from map + # marked U */

    if (ui > 0) {
        unmapped_markers = CALLOC((size_t) NumUnmapped, int);
    }

    if (!(male_col > 0 && female_col > 0)) {
	// BIG NOTE: if there is no sex-specific informaion, copy the sex-averaged info to it...
        // TO DO: need to go through the code and make sure that this is no longer required.
        for (i = LTop1->PhenoCnt; i < LTop1->LocusCnt; i++) {
            LTop1->Marker[i].pos_male =
	      LTop1->Marker[i].pos_female =
	      LTop1->Marker[i].pos_avg;
        }
    }

    return make_EXLTop_from_LTop(LTop1, male_col, female_col);

}   /* end of read_map_file */


void check_size_dos(void)
{
    int i, dos_format[4], empty[4];
    char file_name[9];
    FILE *fp;

    for (i=0; i < 4; i++) {
        dos_format[i]=0; empty[i]=0;
        switch(i) {
        case 0:
            sprintf(file_name, "Pedigree");
            break;
        case 1:
            sprintf(file_name, "Locus");
            break;
        case 2:
            sprintf(file_name, "Map");
            break;
        case 3:
            sprintf(file_name, "Omit");
            break;
        default:
            break;
        }

        if (mega2_input_files[i] != NULL && strcmp(mega2_input_files[i], "-")) {

            empty[i] = check_empty(mega2_input_files[i]);

            if (empty[i]) {
                empty_file(file_name, NON_FATAL);
                errorf("Please re-start mega2 after correcting above errors.");
                continue;
            }

	    if ((fp = fopen(mega2_input_files[i], "r")) == NULL) {
	      errorvf("could not open %s for reading!\n", mega2_input_files[i]);
	      EXIT(FILE_READ_ERROR);
	    }

            dos_format[i] = check_dos_file(fp);
            fclose(fp);

            switch(dos_format[i]) {
            case 0:
#if defined(_WIN) || defined(MINGW)
#ifndef HIDESTATUS
                sprintf(err_msg, "%s file appears to be a Unix format file.",
                        file_name);
                mssgf(err_msg);
#endif
#endif
                break;

            case -1:
                sprintf(err_msg, "%s file appears to be a Mac-format file.",
                        file_name);
                errorf(err_msg);
                sprintf(err_msg, "Please convert to Unix format and restart Mega2.");
                errorf(err_msg);
                break;

            case 1:
                sprintf(err_msg, "%s file appears to be a DOS format file.",
                        file_name);
                errorf(err_msg);
                sprintf(err_msg, "Please convert to Unix format and restart Mega2.");
                errorf(err_msg);
                break;

            default:
                break;
            }
        }
    }

    if (empty[0]==0 && empty[1]==0 &&
       empty[2]==0 && empty[3]==0) {
        ;
    } else {
        EXIT(INPUT_DATA_ERROR);
    }

    if (dos_format[0] == 0 && dos_format[1]==0 &&
       dos_format[2]==0 && dos_format[3]==0) {
        return;
    } else {
        if (dos_format[0] == 1 || dos_format[1] == 1 ||
            dos_format[2] == 1) {
            mssgf("To convert DOS format files, get the dos2unix utility.");
            mssgf("Alternatively, you can use perl as shown below:");
            mssgf("  perl -i -pe \'s/\\r//g\' file_name");
        } else if (dos_format[0] == -1 || dos_format[1] == -1 ||
                 dos_format[2] == -1) {
            mssgf("To convert Mac OS format files, get the dos2unix (or mac2unix) utility.");
            mssgf("Alternatively, you can use perl as shown below:");
            mssgf("  perl -pi -e \'s/\\r/\\n/g\' file_name");
        }

        EXIT(INPUT_DATA_ERROR);
    }
}

#ifdef old
int check_dos_file(FILE *fp)

{
    char c, nextc;

    /* check pedigree file */

    while(!feof(fp))  {
        c=fgetc(fp);
        if (c == CR) {
            if ((nextc=fgetc(fp)) == LF) {
                /*	printf("DOS format file\n"); */
                return 1;
            } else {
                /*	printf("CR without LFS!\n");
                        fclose(fp);
                        dos_format[0] = -1;
                */
                return -1;
            }
        }
    }
    /* UNix format */
    return 0;
}
#else
int check_dos_file(FILE *fp)

{
    char buffer[4096];
    char c, cold, *cp;

    /* check pedigree file */

    cold = 0;
    while(!feof(fp))  {
        *buffer = 0;
        (void)fgets(buffer, sizeof(buffer), fp);
        if (*buffer == 0) break;
        for (cp = buffer; (c = *cp++) != 0; ) {
            if (cold == CR) {
                if (c == LF) {
                /*	printf("DOS format file\n"); */
                    return 1;
                } else {
                /*	printf("CR without LFS!\n");
                        fclose(fp);
                        dos_format[0] = -1;
                */
                    return -1;
                }
            }
            cold = c;
        }
    }
    /* UNix format */
    return 0;
}
#endif

int check_empty(char *file_name)

{
    FILE *fp;

    fp=fopen(file_name, "r");
    (void) fgetc(fp);
    if (feof(fp)) {
        fclose(fp);
        return 1;
    } else {
        fclose(fp); return 0;
    }
}
