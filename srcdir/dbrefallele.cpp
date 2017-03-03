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


extern DBlite MasterDB;

extern int  db_exists_db();
extern void db_open_db();
extern char DBfile[255];


void Reference_Allele_Table::read_ref_allele_file(linkage_ped_top *Top, Str filename, bool use_bp_sort, bp_order *bp) {
    mssgvf("Loading reference allele file %s into database\n",filename.c_str());
    //we should only get here if we don't have a table (no need to drop)
    //but drop just to be safe/for testing
    drop();

    //make our new table
    //printf("Creating database table\n");
    create();


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
                    if (locus == Top->LocusTop->LocusCnt - 1)
                        break;
                    int position = bp->pos;
                    int chromosome = bp->chr;
                    //if the chromosome is the same, and position +/- .1% is the same (not many precise matches)
                    //printf("Internal Chromsome = %d, Reference Chromosome = %d\n", Top->LocusTop->Locus[locus].Marker->chromosome, chr);
                    //printf("Internal Position = %d, Reference Position = %d\n",position, pos);
                    if (chromosome == chr && pos == position) {
                        //we insert our values, the position internally is inserted sow e can select on it
                        insert(chromosome, position, locus, token);
                        bp++;
                        locus++;
                        success++;
                    } else if (pos > position && chr == chromosome) {
                        //if we find a value too large insert a dummy and increment
                        insert(chromosome, position, locus, dummy);
                        bp++;
                        //locus++;
                        fail++;
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
                    if (locus == Top->LocusTop->LocusCnt)
                        break;
                    int position = Top->EXLTop->EXLocus[locus].positions[base_pair_position_index];
                    int chromosome = Top->LocusTop->Locus[locus].Marker->chromosome;
                    if(chromosome == chr && position == pos) {
                        insert(chromosome, position, locus, token);
                        locus++;
                    } else if (pos > position) {
                        insert(chromosome, position, locus, dummy);
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

    mssgvf("Successfully matched %d/%d variants in the dataset to the provided reference panel.\n", success,success+fail);

    //final commit just in case
    MasterDB.commit();
    //delete our insert statement
    close();
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

//old way probably removing this
//void Reference_Allele_Table::insert_into_table(int chr, int pos, Str ref) {
    //db_open_db();
    //MasterDB.begin();

    //DBstmt *insert;
    //char insert_string[255];
    //sprintf(insert_string,"INSERT INTO ref_allele_table (chr, pos, marker, ref) VALUES (%d,%d,%d,'%s')",chr,pos,0,ref.c_str());
    //insert = MasterDB.prep(insert_string);
    //insert->step();
    //MasterDB.commit();
    //insert->fini();
    //MasterDB.close();
//}

//old way probably removing this
//int Reference_Allele_Table::get_marker(int pos) {
    //db_open_db();
    //MasterDB.begin();
//
//    DBstmt *select;
//    char select_string[255];
//    //I think we want to change this to not just 1
//    sprintf(select_string,"SELECT marker FROM map_table WHERE position = %d.0 AND map = 1",pos);
//    //printf("%s\n",select_string);
//    select = MasterDB.prep(select_string);
//    int ret = select && select->abort();
//    while (ret){
//        ret = select->step();
//        //marker found
//        if (ret == SQLITE_ROW) {
//            select->column(0, ret);
//        } else if (ret == SQLITE_DONE) {
//            ret = 0;
//        } else {
//            ret = 0;
//        }
//    }
//    MasterDB.commit();

    //return ret;

//}