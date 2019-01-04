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

#include <cstring>

#include "common.h"
#include "typedefs.h"

#include "dbrefallele.h"
#include "fcmap_ext.h"

#include "user_input_ext.h"
#include "error_messages_ext.h"
#include "read_files_ext.h"

#include "zlib-1.2.8/zlib.h"

#include "dbrefallele.h"

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
    create();

    //make our flip table
    Reference_Flips_Table *reference_flips_table = new Reference_Flips_Table();
    reference_flips_table->create();
    reference_flips_table->init();

    //create our db statements
    init();

    //begin reading our gzipped file
    //we define a buffer
#define BUFLENGTH (0x1000)
    //get a gzipfile and open it
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
    int ttl = 0;

    int flipfail =0;

    int refset = 0;
    int altset = 0;

    printf("Matching position values between dataset and reference, this may take a while especially for larger GWAS datasets. (a minute or more)\n");
    //read our buffer
    char buffer[BUFLENGTH];
    char *token;
    char *ref = (char *)"";
    char *alt = (char *)"";

    int off = 0;
//  int bufc = 0; // needed by printf("TK%d ...
    int c;
    char *bufp, *obufp;
    if(use_bp_sort) {
        while (1) {
            int err;
            int bytes_read;
            bytes_read = gzread(file, buffer+off, BUFLENGTH-off - 1);
            bytes_read += off;
            buffer[bytes_read] = '\0';

            bufp = obufp = buffer;
//          while ( (c = *bufp++) && c != '\n') ;
            for (c = *bufp++; c && c != '\n'; c = *bufp++) ;
            bufp[-1] = 0;
            token = std::strtok(obufp, " \t\n");

//        printf("TK%d %s\n", ++bufc, token);

            while (c != 0) {
                //this was minus 1... but that seems to have caused a row not to process.
                if (locus == Top->LocusTop->LocusCnt)
                    break;

                // if chr is not set then we found a chr
                if (chr == 0)
                    chr = atoi(token);
                    //if chr is set but pos isn't we found a pos
                else if (pos == 0)
                    pos = atoi(token);
                    //if both are set we found a whole entry
                else if(refset == 0) {
                    ref = token;
                    refset = 1;
                }
                else if(altset == 0) {
                    alt = token;
                    altset = 1;
                }

                if( chr != 0  && pos != 0 && refset == 1 && altset == 1) {
                    int position = bp->pos;
                    int chromosome = bp->chr;

                    ttl++;
//            printf("LN: %d %d %s %s\n", chr, pos, ref, alt);
                    if (chromosome == chr && pos == position && strlen(ref) != 0 && strlen(alt) != 0) {
                        //we insert our values, the position internally is inserted so we can select on it
                        insert(chromosome, position, locus, ref, alt);
                        flipfail = reference_flips_table->determine_flips(Top, locus, Top->LocusTop->Locus[locus].Allele[0].AlleleName,Top->LocusTop->Locus[locus].Allele[1].AlleleName, ref, alt, chromosome, position);
                        if(flipfail == 1)
                            fail++;
                        success++;
                        bp++;
                        locus++;
//           printf("MCH: %d %d %s %s\n", chr, pos, ref, alt);
                    } else if (pos > position && chromosome == chr) {
                        //if we find a value too large insert a dummy and increment
                        SECTION_LOG(ref_not_available);
                        mssgvf("chr%d:%d %s has no reference value. \n",chromosome, position, Top->LocusTop->Locus[locus].LocusName);
                        insert(chromosome, position, locus, dummy, dummy);
                        bp++;
                        locus++;
                        fail++;
//           printf("ref>: %d %d %s %s\n", chr, pos, ref, alt);
                        if (bp->pos == pos)
                            continue;
                    }
                    else if( chr > chromosome) {
                        SECTION_LOG(ref_not_available);
                        mssgvf("chr%d:%d %s has no reference value. \n",chromosome, position, Top->LocusTop->Locus[locus].LocusName);
                        insert(chromosome, position, locus, dummy, dummy);
                        bp++;
                        locus++;
                        fail++;
//           printf("chr>: %d %d %s %s\n", chr, pos, ref, alt);
                        continue;
                    }
                    chr = 0;
                    pos = 0;
                    refset = 0;
                    altset = 0;
                }
                token = std::strtok(NULL, " \t\n");
                if (token == NULL) {
                    obufp = bufp;
//                  while ( (c = *bufp++) && c != '\n') ;
                    for (c = *bufp++; c && c != '\n'; c = *bufp++) ;
                    bufp[-1] = 0;
//           printf("SS: %s\n", obufp);
                    if (c != 0)
                        token = std::strtok(obufp, " \t\n");
                }
            }
            off = strlen(obufp);
            strcpy(buffer, obufp);

            //when the buffer is done
            if (bytes_read < BUFLENGTH - 1) {
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
    SECTION_LOG_FINI(ref_mismatch);
    SECTION_LOG_FINI(ref_not_available);

    mssgvf("Matched %d/%d variants in the dataset to the provided reference panel.\n", success,success+fail);
    msgvf("Total records read from ref allele file = %d\n", ttl);
    //final commit just in case
    MasterDB.commit();
    //delete our insert statement
    close();
    reference_flips_table->close();
    gzclose (file);
}

//user input to get a the filename
Str Reference_Allele_Table::get_filename(){
    char input[255];
    while(1) {
        printf("You can use an external reference panel to get a set of reference alleles.\n");
        printf("This process is described in the section called\n");
        printf("'External Reference Allele Panel in the Database'\n");
        printf("n the Mega2 documentation.\n");
        printf("Reference panels are 4 column files of CHR POS REF ALT that are then gzipped.\n");
        printf("They can be constructed by hand or using a shell script included with Mega2\n");
        printf("called GetRefAlleles.sh.  Additionally we provide a reference of 1000 genomes\n");
        printf("most recent build at https://watson.hgen.pitt.edu/mega2/refs/ \n\n");


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

int Reference_Flips_Table::determine_flips(linkage_ped_top *Top, int locus, const char *data_ref, const char *data_alt, char *ref_ref, char *ref_alt, int chromosome, int position) {
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

    int failed = 0;

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
                    } else if (canonra == canonA) {
                        major_minor = 1;
                    }
                }
                else if(canondr == canonC) {
                    if (canonrr == canonC) {
                        major_minor = 0;
                    } else if (canonra == canonC) {
                        major_minor = 1;
                    }
                }
                else if(canondr == canonG) {
                    if (canonrr == canonG) {
                        major_minor = 0;
                    } else if (canonra == canonG) {
                        major_minor = 1;
                    }
                }
                else if(canondr == canonT) {
                    if (canonrr == canonT) {
                        major_minor = 0;
                    } else if (canonra == canonT) {
                        major_minor = 1;
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

    if(major_minor == 0 && strand == 0 && canonrr != canondr) {
        SECTION_LOG(ref_mismatch);
        mssgvf("chr%d:%d %s alleles (%s, %s) not resolvable, ref. alleles (%s, %s). \n", chromosome, position,Top->LocusTop->Locus[locus].LocusName,data_ref,data_alt, ref_ref,ref_alt);
        failed = 1;
    }

    insert(locus, strand, major_minor, dummy, data_ref, data_alt);
    return failed;
}

void Reference_Flips_Table::flip_strands(linkage_ped_top *Top) {
    HMapii strand_flips;
    HMapii major_minor_flips;
    HMapis references;
    HMapis alternates;
    HMapii dummies;

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
        int dummy = 0;
        ret = select->step();
        if (ret == SQLITE_ROW) {
            select->column(0, marker);
            select->column(1, reference);
            select->column(2, alternate);
            select->column(3, dummy);
            references[marker] = canonical_allele(reference);
            alternates[marker] = canonical_allele(alternate);
            dummies[marker] = dummy;
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
        //first deal with cases like A/.
        //with ref being A/Known
        if(Top->LocusTop->Locus[locus].Allele[1].AlleleName == dummycanon) {
            if (references[Top->LocusTop->Locus[locus].locus_link] == Top->LocusTop->Locus[locus].Allele[0].AlleleName) {
                if (alternates[Top->LocusTop->Locus[locus].locus_link] == canonA)
                    Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonA;
                if (alternates[Top->LocusTop->Locus[locus].locus_link] == canonC)
                    Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonC;
                if (alternates[Top->LocusTop->Locus[locus].locus_link] == canonG)
                    Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonG;
                if (alternates[Top->LocusTop->Locus[locus].locus_link] == canonT)
                    Top->LocusTop->Locus[locus].Allele[1].AlleleName = canonT;
            }
        }

        if (strand_flips[Top->LocusTop->Locus[locus].locus_link] == 1) {
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
            }
        }



            //case where we have a dummy and nondummy but the known value is the reference allele
        if (major_minor_flips[Top->LocusTop->Locus[locus].locus_link] == 0 && dummies[Top->LocusTop->Locus[locus].locus_link] == 1) {
            for (int ped = 0; ped < Top->PedCnt; ped++) {
                linkage_ped_tree *tpedtreep = &(Top->PedRaw[ped]);
                for (int per = 0; per < Top->PedRaw[ped].EntryCnt; per++) {
                    linkage_ped_rec *tpersonp = &(tpedtreep->Entry[per]);
                    //this causes a seg fault.

                    void *mk = tpersonp->Marker;
                    int a1, a2;

                    //in flipping so value stays the same
                    //unless it's dummy in which case it becomes 2
                    get_2alleles(mk, locus, &a1, &a2);
                    if (a1 == 0)
                        a1 = 2;
                    else if (a1 == 1)
                        a1 = 1;
                    if (a2 == 0)
                        a2 = 2;
                    else if (a2 == 1)
                        a2 = 1;

                    set_2alleles(mk, locus, &Top->LocusTop->Locus[locus], a1, a2);
                }
            }
        }

        //dummy and nondummy but known value is alt allele
        if (major_minor_flips[Top->LocusTop->Locus[locus].locus_link] == 1 && dummies[Top->LocusTop->Locus[locus].locus_link] == 1) {
            for (int ped = 0; ped < Top->PedCnt; ped++) {
                linkage_ped_tree *tpedtreep = &(Top->PedRaw[ped]);
                for (int per = 0; per < Top->PedRaw[ped].EntryCnt; per++) {
                    linkage_ped_rec *tpersonp = &(tpedtreep->Entry[per]);
                    //this causes a seg fault.

                    void *mk = tpersonp->Marker;
                    int a1, a2;

                    get_2alleles(mk, locus, &a1, &a2);
                    if (a1 == 0)
                        a1 = 1;
                    else if (a1 == 1)
                        a1 = 1;
                    if (a2 == 0)
                        a2 = 1;
                    else if (a2 == 1)
                        a2 = 2;

                    set_2alleles(mk, locus, &Top->LocusTop->Locus[locus], a1, a2);
                }
            }
        }

        //major_minor flip with no dummy values
        if (major_minor_flips[Top->LocusTop->Locus[locus].locus_link] == 1 && dummies[Top->LocusTop->Locus[locus].locus_link] == 0) {
                int refindex, placeholderindex;
                const char *refallelename, *placeholderallelename;
                double reffrequency, placeholderfrequency;
                int reflocuslink, placeholderlocuslink;
                int reference_allele_position = 0;
                for (int allele = 0; allele < Top->LocusTop->Locus[locus].AlleleCnt; allele++) {
                    if ((references[Top->LocusTop->Locus[locus].locus_link]) ==
                        Top->LocusTop->Locus[locus].Allele[allele].AlleleName) {
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
