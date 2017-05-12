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
#include <ctype.h>
#include <ctime>
#include <cstring>

#include "common.h"
#include "typedefs.h"

#include "dbrefallele.h"
#include "fcmap_ext.h"

#include "zlib-1.2.8/zlib.h"
#include "dblite.hh"
#include "dbmisc.hh"
#include "dbrefallele.h"

#include "user_input_ext.h"
#include "error_messages_ext.h"
#include "read_files_ext.h"


extern DBlite MasterDB;

extern int  db_exists_db();
extern void db_open_db();
extern char DBfile[255];

SECTION_LOG_INIT(ref_mismatch);
SECTION_LOG_INIT(ref_not_available);

void Reference_Allele_Table::read_ref_allele_file(linkage_ped_top *Top, Str filename, bool use_bp_sort, bp_order *bp) {
    mssgvf("Loading reference allele file %s into database\n",filename.c_str());
    //we should only get here if we don't have a table (no need to drop)
    //but drop just to be safe/for testing
    drop();

    //make our new table
    //printf("Creating database table\n");
    create();

    //make our flip table
    Reference_Flips_Table *reference_flips_table = new Reference_Flips_Table();
    reference_flips_table->create();
    reference_flips_table->init();

    //now we typically do this on database creation and use the first page
    //get our filename of our reference file
    //Str filename = get_filename();

    //create our db statements
    init();

    //begin reading our gzipped file
    //we define a buffer
    int length = 0x1000;
    //get a gzipfile and open it
    //printf("Reading database file.\n");
    gzFile file;
    file = gzopen(filename.c_str(),"r");
    //check for errors
    if (! file) {
        fprintf (stderr, "gzopen of '%s' failed\n", filename.c_str());
        exit (EXIT_FAILURE);
    }

    //begin a database transaction (for the first buffer)
    MasterDB.begin();

    int locus = Top->LocusTop->PhenoCnt;
    if (use_bp_sort)
        bp += locus;



    // we want to break here since this will cause an error
    if(base_pair_position_index == 0 ) {
        printf("base pair position index = %d",base_pair_position_index);
    }

    int chr = 0;
    int pos = 0;
    char * dummy = (char *)"*";

    int success = 0;
    int fail = 0;

    printf("Matching position values between dataset and reference, this may take a while especially for larger GWAS datasets. (a minute or more)\n");
    //read our buffer
    if(use_bp_sort) {
        while (1) {
            int err;
            int bytes_read;
            char buffer[length];
            bytes_read = gzread(file, buffer, length - 1);
            buffer[bytes_read] = '\0';

            //Str ref = "";
            //int marker = 0;
            //split our buffer by lines
            char *token = std::strtok(buffer, " \n");
            while (token != NULL) {
                // if chr is not set then we found a chr
                if (chr == 0)
                    chr = atoi(token);
                    //if chr is set but pos isn't we found a pos
                else if (pos == 0)
                    pos = atoi(token);
                    //if both are set we found a whole entry
                else {
                    char *ref = new char[255];
                    char *alt = new char[255];
                    //set ref to the current token
                    if(token != NULL)
                        strcpy(ref, token);
                    else
                        strcpy(ref, dummy);
                    //grab the next token, this will be our alternate allele
                    token = std::strtok(NULL, " \n");
                    if (token != NULL)
                        strcpy(alt, token);
                    else
                        strcpy(alt, dummy);
                    //printf("%s/%s\n",ref,alt);
                    if (locus == Top->LocusTop->LocusCnt - 1)
                        break;
                    int position = bp->pos;
                    int chromosome = bp->chr;
                    //if the chromosome is the same, and position +/- .1% is the same (not many precise matches)
                    //printf("Internal Chromsome = %d, Reference Chromosome = %d\n", Top->LocusTop->Locus[locus].Marker->chromosome, chr);
                    //printf("Internal Position = %d, Reference Position = %d\n",position, pos);
                    if (chromosome == chr && pos == position) {
                        //we insert our values, the position internally is inserted so we can select on it
                        insert(chromosome, position, locus, ref, alt);
                        reference_flips_table->determine_flips(Top, locus, Top->LocusTop->Locus[locus].Allele[0].AlleleName,Top->LocusTop->Locus[locus].Allele[1].AlleleName, ref, alt, chromosome, position);
                        //printf("Ref:%s/%s    Data:%s/%s\n",ref,alt,Top->LocusTop->Locus[locus].Allele[0].AlleleName,Top->LocusTop->Locus[locus].Allele[0].AlleleName);
                        bp++;
                        locus++;
                        success++;
                    } else if (pos > position && chromosome == chr) {
                        //if we find a value too large insert a dummy and increment
                        SECTION_LOG(ref_not_available);
                        mssgvf("chr%d:%d has no reference value. \n",chromosome, position);
                        insert(chromosome, position, locus, dummy, dummy);
                        bp++;
                        locus++;
                        fail++;
                    }
                    else if( chr > chromosome) {
                        SECTION_LOG(ref_not_available);
                        mssgvf("chr%d:%d has no reference value. \n",chromosome, position);
                        insert(chromosome, position, locus, dummy, dummy);
                        bp++;
                        locus++;
                        fail++;
                        continue;
                    }
                    chr = 0;
                    pos = 0;

                }
                token = std::strtok(NULL, " \n");
            }

            //when the buffer is done
            if (bytes_read < length - 1) {
                //we commit transactions by buffer for speed (rather than by row)
                //is file done?
                if (gzeof(file)) {
                    break;
                } else {
                    const char *error_string;
                    error_string = gzerror(file, &err);
                    if (err) {
                        fprintf(stderr, "Error: %s.\n", error_string);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

        //we use this after the database read as bp_sort can't be accessed but the values are sorted.
        //so this is similar but refences top and is only used after db_read
    else {
        while (1) {
            int err;
            int bytes_read;
            char buffer[length];
            bytes_read = gzread(file, buffer, length - 1);
            buffer[bytes_read] = '\0';

            char *token = std::strtok(buffer, " \n");
            while (token != NULL) {
                if (chr == 0)
                    chr = atoi(token);
                else if (pos == 0)
                    pos = atoi(token);
                else {
                    char *ref = new char[255];
                    char *alt = new char[255];
                    if(token != NULL)
                        strcpy(ref, token);
                    else
                        strcpy(ref, dummy);
                    token = std::strtok(NULL, " \n");
                    if (token != NULL)
                        strcpy(alt, token);
                    else
                        strcpy(alt, dummy);

                    if (locus == Top->LocusTop->LocusCnt)
                        break;
                    int position = Top->EXLTop->EXLocus[locus].positions[base_pair_position_index];
                    int chromosome = Top->LocusTop->Locus[locus].Marker->chromosome;
                    if(chromosome == chr && position == pos) {
                        insert(chromosome, position, locus, ref, alt);
                        locus++;
                    } else if (pos > position) {
                        SECTION_LOG(ref_not_available);
                        mssgvf("chr%d:%d has no reference value. \n",chromosome, position);
                        insert(chromosome, position, locus, dummy, dummy);
                        locus++;
                    }
                    chr = 0;
                    pos = 0;
                }
                token = std::strtok(NULL, " \n");
            }
            if (bytes_read < length - 1) {
                if (gzeof(file)) {
                    break;
                } else {
                    const char *error_string;
                    error_string = gzerror(file, &err);
                    if (err) {
                        fprintf(stderr, "Error: %s.\n", error_string);
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }

    SECTION_LOG_FINI(ref_mismatch);
    SECTION_LOG_FINI(ref_not_available);

    mssgvf("Matched %d/%d variants in the dataset to the provided reference panel.\n", success,success+fail);

    //final commit just in case
    MasterDB.commit();
    //delete our insert statement
    close();
    reference_flips_table->close();
    //close db connection
    //MasterDB.close();
    //close gzipped file
    gzclose (file);
}

//user input to get a the filename
Str Reference_Allele_Table::get_filename(){
    char input[255];
    while(1) {
        printf("You can use an external reference panel to get a set of reference alleles.\n");
        printf("This process is described in the section called 'External Reference Allele Panel in the Database'\n");
        printf("n the Mega2 documentation.\n\n");
        printf("Reference panels are 3 column files of CHR POS REF that are then gzipped.\n");
        printf("They can be constructed by hand or using a shell script included with Mega2\n");
        printf("called GetRefAlleles.sh.  Additionally we provide a reference of 1000 genomes\n");
        printf("most recent build at ____________. \n\n");


        printf("Enter filename of reference panel > ");
        fcmap(stdin, "%s", input);
        newline;

        if (FILE *file = fopen(input, "r")) {
            fclose(file);
            printf("File selected %s is valid.\n Opening with zlib and importing data into database.\n", input);
            return input;
        }
        else
            printf("Not a valid file.\n");
    }
}

void Reference_Flips_Table::determine_flips(linkage_ped_top *Top, int locus, const char *data_ref, const char *data_alt, char *ref_ref, char *ref_alt, int chromosome, int position) {
    int strand = 0;
    int major_minor = 0;
    int dummy = 0;
    char *dummycanon;
    char *canonA;
    char *canonC;
    char *canonT;
    char *canonG;

    dummycanon = canonical_allele("dummy");
    canonA = canonical_allele("A");
    canonC = canonical_allele("C");
    canonG = canonical_allele("G");
    canonT = canonical_allele("T");

    char *canondr = canonical_allele(data_ref);
    char *canonda = canonical_allele(data_alt);
    char *canonrr = canonical_allele(ref_ref);
    char *canonra = canonical_allele(ref_alt);

    //biallelic
    if (Top->LocusTop->Locus[locus].AlleleCnt == 2) {
        //printf("%s / %s \t %s / %s\t", canondr, canonda, canonrr, canonra);
        if (canondr == dummycanon)
            strand = 0;
        else {
            //First we check our data
            //AT TA CG and GC are ambiguous strand flips
            //so we check for those values in the data
            if((canondr == canonA && canonda == canonT) || (canondr == canonT && canonda == canonA)
                    || (canondr == canonC && canonda == canonG) || (canondr == canonC && canonda == canonG)) {
                strand = 0;
            }
            //these are split up for readability
            //check for GT -> AC and TG -> CA
            else if (((canondr == canonG && canonda == canonT) || (canondr == canonT && canonda == canonG)) &&
                     ((canonrr == canonA && canonra == canonC) || (canonrr == canonC && canonra == canonA))) {
                strand = 1;
            //check for AC -> GT and CA -> TG
            } else if (((canondr == canonA && canonda == canonC) || (canondr == canonC && canonda == canonA)) &&
                       ((canonrr == canonG && canonra == canonT) || (canonrr == canonT && canonra == canonG))) {
                strand = 1;
            }
            //check for CT -> GA and TC -> AG
            else if (((canondr == canonC && canonda == canonT) || (canondr == canonT && canonda == canonC)) &&
                     ((canonrr == canonG && canonra == canonA) || (canonrr == canonA && canonra == canonG))) {
                strand = 1;
            //check for AG -> TC and GA -> CT
            } else if (((canondr == canonA && canonda == canonG) || (canondr == canonG && canonda == canonA)) &&
                       ((canonrr == canonT && canonra == canonC) || (canonrr == canonC && canonra == canonT))) {
                strand = 1;
            //now we handle dummy values
            //if we have an A/. and the reference contains a T etc., we choose the reference
            //also for the case of A/. and the reference is A/X we add the reference
            }else if (canonda == dummycanon) {
                //flip the strand either way
                //flip major minor only if it's the reference alt value
                //I'm still calling this "strand" as the logic for the change will be the same regardless for the labeling
                if(canondr == canonA) {
                    if (canonrr == canonA) {
                        major_minor = 0;
                        strand = 1;
                    } else if (canonra == canonA) {
                        major_minor = 1;
                        strand = 1;
                    }
                }
                else if(canondr == canonC) {
                    if (canonrr == canonC) {
                        major_minor = 0;
                        strand = 1;
                    } else if (canonra == canonC) {
                        major_minor = 1;
                        strand = 1;
                    }
                }
                else if(canondr == canonG) {
                    if (canonrr == canonG) {
                        major_minor = 0;
                        strand = 1;
                    } else if (canonra == canonG) {
                        major_minor = 1;
                        strand = 1;
                    }
                }
                else if(canondr == canonT) {
                    if (canonrr == canonT) {
                        major_minor = 0;
                        strand = 1;
                    } else if (canonra == canonT) {
                        major_minor = 1;
                        strand = 1;
                    }
                }
                dummy = 1;
            }

            else
                strand = 0;

            if (strand == 0) {
                if (canonrr == canonda)
                    major_minor = 1;
                else
                    major_minor = 0;
            } else {
                if (canonrr == canonA && canondr != canonT && canonda == canonT)
                    major_minor = 1;
                else if (canonrr == canonC && canondr != canonG && canonda == canonG)
                    major_minor = 1;
                else if (canonrr == canonG && canondr != canonC && canonda == canonC)
                    major_minor = 1;
                else if (canonrr == canonT && canondr != canonA && canonda == canonA)
                    major_minor = 1;
                else
                    major_minor = 0;
            }
        }
    }
    else {
        strand = 0;

        for(int i = 0; i < Top->LocusTop->Locus[locus].AlleleCnt; i++){
            if(canonrr == Top->LocusTop->Locus[locus].Allele[i].AlleleName)
                major_minor = 1;
        }
    }

    //printf("%d / %d\n",strand, major_minor);



    if(major_minor == 0 && strand == 0 && canonrr != canondr) {
        SECTION_LOG(ref_mismatch);
        mssgvf("chr%d:%d alleles (%s, %s) not trivially comparable, and do not match reference alleles (%s, %s). \n", chromosome, position,data_ref,data_alt, ref_ref,ref_alt);
    }

    insert(locus, strand, major_minor, dummy);
//do we want this?
//    if(strand == 1){
//        DBstmt *update1;
//        DBstmt *update2;
//        char update_string1[255];
//        char update_string2[255];
//        if(canondr == canonG && canonda == canonT) {
//            sprintf(update_string1,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 1;", canonC, locus);
//            sprintf(update_string2,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 2;", canonA, locus);
//        }
//        else if(canondr == canonT && canonda == canonG) {
//            sprintf(update_string1,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 1;", canonA, locus);
//            sprintf(update_string2,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 2;", canonC, locus);
//
//        }
//        else if(canondr == canonC && canonda == canonA) {
//            sprintf(update_string1,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 1;", canonG, locus);
//            sprintf(update_string2,"UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 2;", canonT, locus);
//        }
//        else if(canondr == canonA && canonda == canonC) {
//            sprintf(update_string1, "UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 1;",canonT, locus);
//            sprintf(update_string2, "UPDATE allele_table SET AlleleName = '%s' WHERE locus_link = %d AND indexX = 2;",canonG, locus);
//        }
//
//        update1 = MasterDB.prep(update_string1);
//        update1->step();
//        delete update1;
//
//        update2 = MasterDB.prep(update_string2);
//        update2->step();
//        delete update2;
//
//    }

}

void Reference_Flips_Table::flip_strands(linkage_ped_top *Top) {
    HMapii strand_flips;
    HMapii major_minor_flips;
    HMapis references;
    HMapis alternates;

    char *dummycanon;
    char *canonAllele;
    char *canonA;
    char *canonC;
    char *canonT;
    char *canonG;


    dummycanon = canonical_allele("dummy");
    canonA = canonical_allele("A");
    canonC = canonical_allele("C");
    canonG = canonical_allele("G");
    canonT = canonical_allele("T");

    db_open_db();
    MasterDB.begin();
    DBstmt *select;
    char select_string[255];
    sprintf(select_string, "SELECT marker, ref, alt FROM ref_allele_table;");
    select = MasterDB.prep(select_string);
    int ret = select && select->abort();

    //get position and reference and put them into our map
    while (ret) {
        int marker = 0;
        char *reference;
        char *alternate;
        ret = select->step();
        if (ret == SQLITE_ROW) {
            select->column(0, marker);
            select->column(1, reference);
            select->column(2, alternate);
            references[marker] = canonical_allele(reference);
            alternates[marker] = canonical_allele(alternate);
        } else
            break;
    }

    MasterDB.commit();
    delete select;

    DBstmt *select2;
    char select_string2[255];
    sprintf(select_string2, "SELECT marker, strand, major_minor FROM ref_allele_flips");
    select2 = MasterDB.prep(select_string2);
    int ret2 = select2 && select2->abort();

    while (ret2) {
        int locus = 0;
        int strand_flip = 0;
        int major_minor = 0;
        ret2 = select2->step();
        if (ret2 == SQLITE_ROW) {
            select2->column(0, locus);
            select2->column(1, strand_flip);
            select2->column(2, major_minor);
            strand_flips[locus] = strand_flip;
            major_minor_flips[locus] = major_minor;
        } else
            break;
    }


    delete select2;



    for(int locus = Top->LocusTop->PhenoCnt; locus < Top->LocusTop->LocusCnt; locus++){
        if (strand_flips[Top->LocusTop->Locus[locus].locus_link] == 1) {
                //first deal with cases like A/.
                if(Top->LocusTop->Locus[locus].Allele[1].AlleleName == dummycanon) {
                    //two cases, The allele in the dataset is in the reference
                    //or the allele is a compliment T/. -> A/something
                    //So we technically set this wrong, this is so the alleles and person level markers will change later
                    //this also got more complicated than it needed to be since I couldn't set Allelename = refere[] etc.
                    if(strand_flips[Top->LocusTop->Locus[locus].locus_link] == 1) {
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonA)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonA;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonC)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonC;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonG)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonG;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonT)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonT;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonA)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonA;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonC)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonC;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonG)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonG;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonT)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonT;
                    }
                    else{
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonA)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonA;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonC)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonC;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonG)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonG;
                        if(references[Top->LocusTop->Locus[locus].locus_link] == canonT)
                            Top->LocusTop->Locus[locus].Allele[0].AlleleName = canonT;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonA)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonA;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonC)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonC;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonG)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonG;
                        if(alternates[Top->LocusTop->Locus[locus].locus_link] == canonT)
                            Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonT;
                        //if only...
                        //Top->LocusTop->Locus[locus].Allele[0].AlleleName = references[Top->LocusTop->Locus[locus].locus_link];
                        //Top->LocusTop->Locus[locus].Allele[1].AlleleName = alternates[Top->LocusTop->Locus[locus].locus_link];
                    }
                }
                //now other cases
                else {
                    for (int i = 0; i <= 1; i++) {
                        canonAllele = canonical_allele(Top->LocusTop->Locus[locus].Allele[i].AlleleName);
                        //printf("Before: %s\n",canonAllele);
                        if (canonAllele == canonA)
                            Top->LocusTop->Locus[locus].Allele[i].AlleleName = canonT;
                        else if (canonAllele == canonC)
                            Top->LocusTop->Locus[locus].Allele[i].AlleleName = canonG;
                        else if (canonAllele == canonG)
                            Top->LocusTop->Locus[locus].Allele[i].AlleleName = canonC;
                        else if (canonAllele == canonT)
                            Top->LocusTop->Locus[locus].Allele[i].AlleleName = canonA;
                        //canonAllele = canonical_allele(Top->LocusTop->Locus[locus].Allele[i].AlleleName);
                        //printf("After: %s\n",canonAllele);
                    }
                }
            }
            if (major_minor_flips[Top->LocusTop->Locus[locus].locus_link] == 1) {
                int refindex, placeholderindex;
                const char *refallelename, *placeholderallelename;
                double reffrequency, placeholderfrequency;
                int reflocuslink, placeholderlocuslink;
                int reference_allele_position = 0;
                for (int allele = 0; allele < Top->LocusTop->Locus[locus].AlleleCnt; allele++) {
                    if ((references[Top->LocusTop->Locus[locus].locus_link]) == Top->LocusTop->Locus[locus].Allele[allele].AlleleName) {
                        reference_allele_position = allele;
                    }
                }
                if (reference_allele_position != 0) {
                    refindex = Top->LocusTop->Locus[locus].Allele[reference_allele_position].index;
                    refallelename = Top->LocusTop->Locus[locus].Allele[reference_allele_position].AlleleName;
                    reflocuslink = Top->LocusTop->Locus[locus].Allele[reference_allele_position].locus_link;
                    reffrequency = Top->LocusTop->Locus[locus].Allele[reference_allele_position].Frequency;

                    placeholderindex = Top->LocusTop->Locus[locus].Allele[0].index;
                    placeholderallelename = Top->LocusTop->Locus[locus].Allele[0].AlleleName;
                    placeholderlocuslink = Top->LocusTop->Locus[locus].Allele[0].locus_link;
                    placeholderfrequency = Top->LocusTop->Locus[locus].Allele[0].Frequency;

                    Top->LocusTop->Locus[locus].Allele[0].index = refindex;
                    Top->LocusTop->Locus[locus].Allele[0].AlleleName = refallelename;
                    Top->LocusTop->Locus[locus].Allele[0].locus_link = reflocuslink;
                    Top->LocusTop->Locus[locus].Allele[0].Frequency = reffrequency;

                    Top->LocusTop->Locus[locus].Allele[reference_allele_position].index = placeholderindex;
                    Top->LocusTop->Locus[locus].Allele[reference_allele_position].AlleleName = placeholderallelename;
                    Top->LocusTop->Locus[locus].Allele[reference_allele_position].locus_link = placeholderlocuslink;
                    Top->LocusTop->Locus[locus].Allele[reference_allele_position].Frequency = placeholderfrequency;


                    for (int ped = 0; ped < Top->PedCnt; ped++) {
                        linkage_ped_tree *tpedtreep = &(Top->PedRaw[ped]);
                        for (int per = 0; per < Top->PedRaw[ped].EntryCnt; per++) {
                            linkage_ped_rec *tpersonp = &(tpedtreep->Entry[per]);
                            //this causes a seg fault.

                            void *mk = tpersonp->Marker;
                            int a1, a2;

                            get_2alleles(mk, locus, &a1, &a2);;
                            if (a1 == 1)
                                a1 = 2;
                            else if (a1 == 2)
                                a1 = 1;
                            if (a2 == 1)
                                a2 = 2;
                            else if (a2 == 2)
                                a2 = 1;

                            set_2alleles(mk, locus, &Top->LocusTop->Locus[locus], a1, a2);
                        }
                    }
                }
        }
    }

}