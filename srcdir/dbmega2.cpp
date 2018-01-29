/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

#include "types.hh"

#include "common.h"
#include "typedefs.h"
#include "error_messages_ext.h"
#include "utils_ext.h"
#include "write_files_ext.h"
#include "tod.hh"

#include "dblite.hh"

#include "dbbatch.hh"
#include "dbmisc.hh"
#include "dbpedigree.hh"
#include "dballele.hh"
#include "dbloci.hh"
#include "dbmap.hh"
#include "dbgenotype.hh"
#include "dbrefallele.h"

using namespace std;

char DBfile[FILENAME_LENGTH] = "";
DBlite MasterDB;
int MasterDBreset = 1;

int database_dump = 0;
int database_read = 0;
int database_off  = 0;

char *DBCreateTime;
char *DBMega2Version;
char *DBversion;
char *SQLversion;

//const char *select1(const char *select)
template<typename cs>
void select1(const char *select, cs &sval)
{
    DBstmt *sver = MasterDB.prep(select);
    if (sver) {
        int sstep = sver->step();
        if (sstep == SQLITE_ROW) {
            sver->column(0, sval);
        } else if (sstep == SQLITE_DONE) {
        }
    }
    delete sver;
}

const char *selectrow(const char *select)
{
    const char *sval = "UNKNOWN";
    DBstmt *sver = MasterDB.prep(select);
    if (sver) {
        int sstep = sver->step();
        if (sstep == SQLITE_ROW) {
            sver->column(0, sval);
        } else if (sstep == SQLITE_DONE) {
        }
    }
    delete sver;
    return sval;
}

const char *selectcol(const char *select)
{
    const char *sval = "UNKNOWN";
    DBstmt *sver = MasterDB.prep(select);
    if (sver) {
        int sstep = sver->step();
        if (sstep == SQLITE_ROW) {
            sver->column(0, sval);
        } else if (sstep == SQLITE_DONE) {
        }
    }
    delete sver;
    return sval;
}

bp_order *bp_sort;

static int bp_fwd (const void *p1, const void *p2) {
    const bp_order *tp1 = (const bp_order *)p1;
    const bp_order *tp2 = (const bp_order *)p2;
    if (tp1->chr == tp2->chr) {
        if (tp1->pos == tp2->pos) {
            return tp1->i - tp2->i;
        } else 
            return tp1->pos - tp2->pos; 
    } else
        return tp1->chr - tp2->chr;
}

static void sort_bp_order(linkage_ped_top *Top)
{
    extern int base_pair_position_index;
    extern int genetic_distance_index;

    bp_sort = CALLOC(Top->LocusTop->LocusCnt, bp_order);
    bp_order *bpp;
    int offset = Top->LocusTop->PhenoCnt;
    int i;

    for (i = 0, bpp = bp_sort; i < offset; i++, bpp++) {
        bpp->i  = i;
        bpp->chr = bpp->pos = 0;
    }

    for (i = offset, bpp = bp_sort+offset; i < Top->LocusTop->LocusCnt; i++, bpp++) {
        bpp->i = i;
        bpp->chr = Top->LocusTop->Locus[i].Marker->chromosome;
        if (base_pair_position_index < 0)
            bpp->pos = Top->EXLTop->EXLocus[i].positions[genetic_distance_index];
        else {
            bpp->pos = Top->EXLTop->EXLocus[i].positions[base_pair_position_index];
        }
    }
    qsort((void *)(bp_sort+offset), Top->LocusTop->MarkerCnt, sizeof (bp_order), bp_fwd);
}

Int_table int_table;
Double_table double_table;
Charstar_table charstar_table;
Stuff_table stuff_table;

Batch_file_parameters batch_file_parameters;
File_table file_table;

Pedigree_table pedigree_table;
Pedigree_brkloop_table pedigree_brkloop_table;
Person_table person_table;
Person_brkloop_table person_brkloop_table;

CanonicalAllele_table canonicalallele_table;
MarkerScheme_table markerscheme_table;

Locus_table locus_table;
Allele_table allele_table;
Marker_table marker_table;
TraitAff_table traitaff_table;
AffectClass_table affectclass_table;
ClassPen_table classpen_table;
TraitQuant_table traitquant_table;

MapNames_table mapnames_table; 
Map_table map_table; 

Phenotype_table phenotype_table;
Genotype_table genotype_table;

void db_index_all() {
    int_table.index();
    double_table.index();
    charstar_table.index();
    stuff_table.index();

    batch_file_parameters.index();
    file_table.index();

    pedigree_table.index();
    person_table.index();

    pedigree_brkloop_table.index();
    person_brkloop_table.index();

    canonicalallele_table.index();
    markerscheme_table.index();

    locus_table.index();
    allele_table.index();
    marker_table.index();
    traitaff_table.index();
    affectclass_table.index();
//    classpen_table.index();
    traitquant_table.index();

    mapnames_table.index();
    map_table.index();

    phenotype_table.index();
    genotype_table.index();

}

void dbmega2_export(linkage_ped_top *Top)
{
    msgvf("Dumping SQLite3 DB ");
#ifndef HIDEFILE
    msgvf("to file \"%s\"\n", DBfile);
#else
    msgvf("\n");
#endif

    sort_bp_order(Top);

    dbmisc_export(Top);

    dbbatch_file_export(Top);

    dbpedigree_export(Top);

    dballele_export(Top, bp_sort);  //sort markerscheme_table

    dblocus_export(Top->LocusTop, bp_sort); // sort locus_table, allele_table, marker_table

    dbmap_export(Top, bp_sort); // sort map_table

    dbgenotype_export(Top, bp_sort); // sort row contents of

    if (mega2_input_files[REFfl] != NULL) {
        if(fopen(mega2_input_files[REFfl], "r") != NULL) {
            Reference_Allele_Table *reference_allele_table = new Reference_Allele_Table();
            reference_allele_table->read_ref_allele_file(Top, mega2_input_files[REFfl], true, bp_sort);
        }
        else{
            mssgvf("Can't find Reference Allele File: %s", mega2_input_files[REFfl]);
        }
    }
}

void dbmega2_stat(linkage_ped_top *Top)
{
    extern int genetic_distance_index;
    extern int genetic_distance_sex_type_map;
    extern int base_pair_position_index;
    extern void show_reset_input();

    log_line(mssgf);
#ifndef HIDEFILE
    msgvf("The path to this SQLite3 database is %s.\n", DBfile);
#endif
    msgvf("This database was created using Mega2 version %s.\n", DBMega2Version);
    msgvf("This database was created using   SQLite3 %s on %s.\n", DBversion, DBCreateTime);
    msgvf("This database was processed using SQLite3 %s on %s.\n", SQLversion, RunDate);

    file_table.stat();
    msgvf("This database contains:\n\t%d persons (%d pedigrees)\n",
          Top->IndivCnt, Top->PedCnt);
    msgvf("\t%d markers\n", Top->LocusTop->MarkerCnt);
    if (Top->LocusTop->PhenoCnt != 1)
        msgvf("\t%d traits\n", Top->LocusTop->PhenoCnt);
    else
        msgvf("\t1 trait\n");
    
    if (genetic_distance_index > -1)
/*
        msgvf("\tgenetic distance(map name/type/idx) \"%s\"/%s/%d, Sex map type %s\n",
              Top->EXLTop->MapNames[genetic_distance_index], 
              Top->EXLTop->map_functions[genetic_distance_index] == 'k' ? "kosambi" : "haldane", 
              genetic_distance_index, 
              genetic_distance_sex_type_map == 0 ? "AVERAGED_MAP" : 
              (genetic_distance_sex_type_map == 1 ? "MALE_MAP" : "FEMALE_SEX_MAP"));
*/
        msgvf("\tgenetic distance(map name/type) \"%s\"/%s, Sex map type %s\n",
              Top->EXLTop->MapNames[genetic_distance_index], 
              Top->EXLTop->map_functions[genetic_distance_index] == 'k' ? "kosambi" : "haldane", 
              genetic_distance_sex_type_map == 0 ? "AVERAGED_MAP" : 
              (genetic_distance_sex_type_map == 1 ? "MALE_MAP" : "FEMALE_SEX_MAP"));

    if (base_pair_position_index > -1)
/*
        msgvf("\tbase pair distance (map name/type/idx) \"%s\"/%s/%d\n\n",
              Top->EXLTop->MapNames[base_pair_position_index], 
              Top->EXLTop->map_functions[base_pair_position_index] == 'p' ? "base pair" : "???",
              base_pair_position_index);
*/
        msgvf("\tbase pair distance(map name) \"%s\"\n",
              Top->EXLTop->MapNames[base_pair_position_index]);

    if(mega2_input_files[REFfl]!= NULL){
        string hgbuild;
        BatchValueGet(hgbuild,"human_genome_build");
        msgvf("\treference allele panel: %s\n", mega2_input_files[REFfl]);
        msgvf("\treference allele build: %s\n", hgbuild.c_str());
    }


    show_reset_input();

    draw_line();
//  write_locus_stats(Top->LocusTop, UNKNOWN);
}

void dbmega2_import(linkage_ped_top *Top)
{
    extern void mk_marker_filter(linkage_ped_top *Top);

    msgvf("Reading SQLite3 DB ");
#ifndef HIDEFILE
    msgvf("from file \"%s\"\n", DBfile);
#else
    msgvf("\n");
#endif

    dbmisc_import(Top);

    dbbatch_file_import(Top);

    dbpedigree_import(Top);

    Tod import_allele("import_allele");
    dballele_import(Top);
    import_allele();

    Tod import_locus("import_locus");
    dblocus_import(Top->LocusTop);
    import_locus();


    Tod import_map("import_map");
    dbmap_import(Top);
    mk_marker_filter(Top);
    import_map();

    dbgenotype_import(Top);

    dbmega2_stat(Top);
}

void db_drop_all() {
    int_table.drop();
    double_table.drop();
    charstar_table.drop();
    stuff_table.drop();

    batch_file_parameters.drop();
    file_table.drop();

    pedigree_table.drop();
    person_table.drop();

    pedigree_brkloop_table.drop();
    person_brkloop_table.drop();

    canonicalallele_table.drop();
    markerscheme_table.drop();

    locus_table.drop();
    allele_table.drop();
    marker_table.drop();
    traitaff_table.drop();
    affectclass_table.drop();
//  classpen_table.drop();
    traitquant_table.drop();

    mapnames_table.drop();
    map_table.drop();

    phenotype_table.drop();
    genotype_table.drop();

}

int db_exists_db() {
    char *dbf = DBfile;

    if (DBfile[0] == 0) { // not on commandline
        BatchValueGet(dbf, "DBfile_name");
    }

    if (strchr(DBfile, '/')) {
        /* leave alone*/;
    } else if ( (database_dump && database_read) && Mega2OutputPath) {
        int l = strlen(Mega2OutputPath);
        memmove(DBfile+l+1, DBfile, strlen(DBfile)+1);
        strcpy(DBfile, Mega2OutputPath);
        DBfile[l] = '/';
        BatchValueSet(dbf, "DBfile_name");
    } 

    FILE *f = fopen(DBfile, "r");
    if (f == NULL) {
        return 0;
    }
    fclose(f);
    return 1;
}

//allows check for if a table is present
//currently this is used for the ref_allele_table but should be extensible
int db_table_exists(const char *table){
    int exists = 0;
    MasterDB.open(DBfile);

    if(db_exists_db()){
        DBstmt *select;

        char select_string[255] = "";
        strcat(select_string,"SELECT name FROM sqlite_master WHERE type='table' AND name='");
        strcat(select_string, table);
        strcat(select_string,"';");

        select = MasterDB.prep(select_string);
        int ret = select && select->abort();
        while (ret) {
            ret = select->step();
            if (ret == SQLITE_ROW)
                exists = 1;
            else
                exists = 0;
            break;
        }

        delete select;
    }
    else {
        exists = 0;
        EXIT(FILE_NOT_FOUND);
    }

    return exists;
}

void db_open_db() {
    extern void delete_file(const char *);
    extern void add_sumdir(char *);
    if (!database_dump && !database_read) return;

    if (database_dump) {
        if (db_exists_db()) {
#ifndef HIDEFILE
            msgvf("Database file \"%s\" will be backed up.\n", DBfile);
#endif
            backup_file(DBfile);
        } else {
#ifndef HIDEFILE
            msgvf("Database file \"%s\" will be created.\n", DBfile);
#endif
        }
    } else {
        if (! db_exists_db()) {
            msgvf("\n");
            errorvf("Database file \"%s\" not found.\n", DBfile);
            EXIT(EARLY_TERMINATION);
        }
    }

    if (! MasterDB.open(DBfile, debug)) {
            EXIT(FILE_NOT_FOUND);
    }
    if (!database_dump && database_read) {
        int int_table_cnt = -1;
        select1("select count(*) from int_table;", int_table_cnt);
        if (int_table_cnt == 0) {
            EXIT(FILE_NOT_FOUND);
        }
    }

    select1("SELECT SQLITE_VERSION();", SQLversion);
#ifndef HIDEFILE
//  msgvf("SQLITE3 DB (%s) Version = %s\n", DBfile, SQLversion);
#endif
}

void db_init_all() {
    MasterDB.exec("PRAGMA page_size = 4096;");

    MasterDB.begin();

    if (database_dump) {
        if (MasterDBreset) db_drop_all();

        int_table.create();
        double_table.create();
        charstar_table.create();
        stuff_table.create();

        batch_file_parameters.create();
        file_table.create();

        pedigree_table.create();
        person_table.create();

        pedigree_brkloop_table.create();
        person_brkloop_table.create();

        canonicalallele_table.create();
        markerscheme_table.create();

        locus_table.create();
        allele_table.create();
        marker_table.create();
        traitaff_table.create();
        affectclass_table.create();
//        classpen_table.create();
        traitquant_table.create();

        mapnames_table.create();
        map_table.create();

        phenotype_table.create();
        genotype_table.create();
    }

    int_table.init();
    double_table.init();
    charstar_table.init();
    stuff_table.init();

    batch_file_parameters.init();
    file_table.init();

    pedigree_table.init();
    person_table.init();

    pedigree_brkloop_table.init();
    person_brkloop_table.init();

    canonicalallele_table.init();
    markerscheme_table.init();

    locus_table.init();
    allele_table.init();
    marker_table.init();
    traitaff_table.init();
    affectclass_table.init();
//    classpen_table.init();
    traitquant_table.init();

    mapnames_table.init();
    map_table.init();

    phenotype_table.init();
    genotype_table.init();

    MasterDB.commit();

/*
    batch_file_parameters.db_set("k1", "vv1", 11);
    batch_file_parameters.db_set("k2", "vv2", 22);
    batch_file_parameters.db_set("k3", "vv3", 33);

    const char *v, *k;
    Vecs ve;
    int i;

    k = "k0"; batch_file_parameters.db_get(k, v, i); printf("get: %s, %s, %d\n", "k0", v, i);
    k = "k1"; batch_file_parameters.db_get(k, v, i); printf("get: %s, %s, %d\n", "k1", v, i);
    k = "k2"; batch_file_parameters.db_get(k, v, i); printf("get: %s, %s, %d\n", "k2", v, i);
    k = "k3"; batch_file_parameters.db_get(k, v, i); printf("get: %s, %s, %d\n", "k3", v, i);

    batch_file_parameters.db_get("k2", ve, i);
    printf("get: %s, %s, %d\n", "k2", v, i);
*/
}

/*
 * This was part of the import process but got moved to after is called dbgenotype_import_genotype to make sure
 * the marker structure is attached to the pedigree person structure.
 */
void strand_flip_reference_alleles(linkage_ped_top *Top) {
    if(_strand_flips) {
        Reference_Flips_Table *ref_flips_table = new Reference_Flips_Table();
        ref_flips_table->flip_strands(Top);
    }
}

void db_fini_all() {

    if (!database_dump && !database_read) return;

    int_table.close();
    double_table.close();
    charstar_table.close();
    stuff_table.close();

    batch_file_parameters.close();
    file_table.close();

    pedigree_table.close();
    person_table.close();

    pedigree_brkloop_table.close();
    person_brkloop_table.close();

    canonicalallele_table.close();
    markerscheme_table.close();

    locus_table.close();
    allele_table.close();
    marker_table.close();
    traitaff_table.close();
    affectclass_table.close();
//  classpen_table.close();
    traitquant_table.close();

    mapnames_table.close();
    map_table.close();

    phenotype_table.close();
    genotype_table.close();

    MasterDB.close();
}

void db_dump_all() {
}
