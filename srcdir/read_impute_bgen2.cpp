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


using namespace std;

#include <stdio.h>

#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <zlib.h>

#include "genfile/bgen/bgen.hpp"

#include "common.h"
#include "tod.hh"
#include "error_messages_ext.h"
extern void           Exit(int arg, const char *file, const int line, const char *err);
#include "typedefs.h"
#include "fcmap_ext.h"
#include "batch_input.h"
#include "batch_input_ext.h"
#include "mrecode.h"
#include "mrecode_ext.h"
#include "reorder_loci_ext.h"
#include "plink_ext.h"
#include "read_files_ext.h"
#include "annotated_ped_file_ext.h"
#include "annotated_ped_file.h"

#include "input_ops.hh"
#include "read_impute.hh"
#include "read_impute_bgen2.hh"


SECTION_ERR_INIT(bad_marker_name);
void ReadBgen2::read_input_file()
{
    BgenParser bgen;
    bgen.m_filename = impute_file;
    bgen.rip = this;

    bgen.read_input_file();
}

void ReadBgen2::build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons)
{
    BgenParserGenotypeReadHelper gh;
    gh.m_filename = impute_file;
    gh.rip = this;

    build_internal_genotypes(LTop, persons, gh);
}

void BgenParserGenotypeReadHelper::genotypes_init()
{
//  asm("int $3");
    open();
    summarise( std::cerr ) ;
    get_sample_ids( [this]( std::string const& id ) { } );

/*
    open(C(impute_file));
    pass = 1;
    read_header();
    pass = 2;
    if (layout < 2) {
        block.cdata = new unsigned char[6 * samples];
        block.rdata = new unsigned char[6 * samples];
        if (layout == 0)
            scale = 10000;
        else if (layout == 1)
            scale = 32768;
    } else if (layout == 2) {
    } else {
//      ERROR
    }
*/
}

boolean BgenParserGenotypeReadHelper::genotypes_marker_hdr(int mrk_idx, string& hmm, string& rsid, string& pos,
                                                     const char* &A, const char* &B)
{
    bool ret;
    string chromosome;
    uint32_t position;
    vector<string> alleles;

    ret = read_variant( &chromosome, &position, &rsid, &alleles );

    hmm = "---";
    fix_marker_pos(pos, position);
    A = C(alleles[0]);
    B = C(alleles[1]);

    prob_sample = 0;
    return ret;
/*
    if (layout == 0) {
        process_10(mrk_idx);
    } else if (layout == 1) {
        process_11(mrk_idx);
    } else if (layout == 2) {
        process_12(mrk_idx);
    }

    if ((flags & compressF) > 0)
        zp = block.cdata;
    else
        zp = 0;

    hmm = "---";
    rsid = block.rsid;
    fix_marker_pos(pos);
    A = C(block.allele[0]);
    B = C(block.allele[1]);

*/
}

void BgenParserGenotypeReadHelper::genotypes_skip()
{
    ignore_probs();
}

void BgenParserGenotypeReadHelper::genotypes_sample_prob(Token::d3& nums)
{
    if (prob_sample == 0)
        read_probs(&probs);

    nums[0] = probs[prob_sample][0];
    nums[1] = probs[prob_sample][1];
    nums[2] = probs[prob_sample][2];
    
    prob_sample++;

    if (prob_sample == m_context.number_of_samples) {
        prob_sample = 0;
    }
/*
    nums[0] = ((double)read_ushort(zp)) / scale;
    nums[1] = ((double)read_ushort(zp)) / scale;
    nums[2] = ((double)read_ushort(zp)) / scale;
*/
}

void BgenParserGenotypeReadHelper::genotypes_end()
{
    close();
/*
    if (layout < 2) {
        delete [] block.cdata;
        delete [] block.rdata;
    }
*/
}


void BgenParser::fix_marker_pos(string& ccpos, unsigned long pos)
{
    char cpos[10];
    sprintf(cpos, "%ld", pos);
    ccpos = string(cpos);
}

void BgenParser::validate_marker_name(string& rsid, string& chrm, string& rsid_field, const string& ccpos)
{
    VecsDB fields;
    split(fields, rsid_field, ":", 3);

    if (fields.size() == 1) {
        rsid = fields[0];
        if (rsid == "." || rsid == "NA" || rsid == "na") {
            chrm = rip->oxford_single_chr;
            rsid = "chr" + rip->oxford_single_chr + "_" + ccpos;
        }
    } else {
        if (fields[0].compare(0, 3, "chr") == 0)
            fields[0].erase(0, 3);
        if (inMap(fields[0], Input->G.chrm_set)) {
            chrm = fields[0];
            rsid = "chr" + fields[0] + "_" + ccpos;
        } else if (fields[0].compare(0, 2, "rs") == 0 && rip->oxford_single_chr != "--") {
            chrm = rip->oxford_single_chr;
            rsid = fields[0];
        } else {
            SECTION_ERR(bad_marker_name);
            errorvf("impute2 file: bad marker name %s\n", C(rsid_field));
            badname++;
        }
    }
}

// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
void BgenParser::read_input_file()
{
//  asm("int $3");
    try {
        open( ) ;

        rip->markers_all      = m_context.number_of_variants;
        rip->markers_filtered = m_context.number_of_variants;

        rip->people_all       = m_context.number_of_samples;
        rip->people_filtered  = m_context.number_of_samples;

        summarise( std::cerr ) ;

/*
        // Output header
        std::cout << "##fileformat=VCFv4.2\n"
            << "FORMAT=<ID=GP,Type=Float,Number=G,Description=\"Genotype call probabilities\">\n"
            << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" ;
        get_sample_ids(
            []( std::string const& id ) { std::cout << "\t" << id ; }
        ) ;
        std::cout << "\n" ;
*/        
        get_sample_ids(
            [this]( std::string const& id ) { sample.samples.push_back(id) ; }
            );
        
        // Output variants
        std::string chromosome ;
        uint32_t position ;
        std::string rsid ;
        std::vector< std::string > alleles ;
//        std::vector< std::vector< double > > probs ;

    
        while( read_variant( &chromosome, &position, &rsid, &alleles )) {
/*
              std::cout << chromosome << '\t'
                << position << '\t'
                << rsid << '\t' ;
            for( std::size_t i = 0; i < alleles.size(); ++i ) {
                std::cout << ( i > 0 ? "," : "" ) << alleles[i] ;
            }
            std::cout << "\t.\t.\t.\tGP" ;
*/
            string ccpos;
            fix_marker_pos(ccpos, position);

            string rsid_name, chrm;
            validate_marker_name(rsid_name, chrm, rsid, ccpos);
            rip->markers.push_back(new ImpMarker(rsid_name, chrm, ccpos, alleles[0], alleles[1], rip->read_info));

            ignore_probs() ;
/*
            read_probs( &probs ) ;

            for( std::size_t i = 0; i < probs.size(); ++i ) {
              std::cout << '\t' ;
                for( std::size_t j = 0; j < probs[i].size(); ++j ) {
                    std::cout << ( j > 0 ? "," : "" ) ;
                    if( probs[i][j] == -1 ) {
                        std::cout << "." ;
                    } else {
                        std::cout << probs[i][j] ;
                    }
                }
            }
          std::cout << "\n" ;
*/
        }
//        return 0 ;
    }
    catch( genfile::bgen::BGenError const& e ) {
        std::cerr << "!! Uh-oh, error parsing bgen file: ";
        std::cerr << e.what() << "\n";
	throw;
//        return -1 ;
    }
}
