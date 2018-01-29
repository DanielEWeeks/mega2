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


#include <stdio.h>

#include <string.h>
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
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

using namespace std;


SECTION_ERR_INIT(bad_marker_name);
void ReadBgen2::read_input_file()
{
    BgenParser bgen;
    bgen.m_filename = impute_file;
    bgen.n_rip = this;

    bgen.read_input_file();

    SECTION_ERR_FINI(bad_marker_name);
}

void ReadBgen2::build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons, vector<Vecc> &veca)
{
    BgenParserGenotypeReadHelper gh;
    gh.m_filename = impute_file;
    gh.n_rip = this;

    build_internal_genotypes(LTop, persons, gh, veca);
}

void BgenParserGenotypeReadHelper::genotypes_init()
{
//  asm("int $3");
    open();
    summarise( std::cerr ) ;
    class noop {
    public:
        noop(BgenParserGenotypeReadHelper *)  {}
       ~noop()  {}
        void operator() ( std::string const& id ) { }
    };
    noop xx(this);
    get_sample_ids( xx );
#if 0
    get_sample_ids(
/*
        [this]( std::string const& id ) { }
*/
        noop(this) );
#endif
}

bool BgenParserGenotypeReadHelper::genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid, 
                                                           std::string& pos, vector<string>& alleles)
{
    bool ret;
    uint32_t position;

    ret = read_variant( &chrm, &position, &rsid, &alleles );

    hmm = "---";
    fix_marker_pos(pos, position);

    n_prob_sample = 0;
    return ret;
}

void BgenParserGenotypeReadHelper::genotypes_skip()
{
    ignore_probs();
}

bool BgenParserGenotypeReadHelper::genotypes_sample_prob(ProbQ& Q)
{
    if (n_prob_sample == 0) {
        read_probs(&n_probs);
        if (m_first) {
            msgvf("samples %d, alleles %d, entries %d, ploidy %d, phased %d\n",
                  (int)m_samples, (int)m_alleles,(int) m_entries, (int)m_ploidy,
                  (int)m_phased);
            m_first = 0;
        }
        if (m_ploidy != 2 || m_phased != 0) {
            errorvf("Mega2 only supports bgen files that are unphased with ploidy = 2\n");
            EXIT(INPUT_DATA_ERROR);
        }
    }

    int idx = 0;
/* an alternate ordering ...
    for (int i = 1; i <= (int) m_alleles; i++)
        for (int j = 1; j <= i; j++) 
            Q.push(ProbID(n_probs[n_prob_sample][idx++], j, i));
*/
    if (m_phased == 0) {
        for (int i = 1; i <= (int) m_alleles; i++)
            for (int j = i; j <= (int) m_alleles; j++) 
                Q.push(ProbID(n_probs[n_prob_sample][idx++], i, j));
    } else {
        double p, pi, pj;
        for (int i = 1; i <= (int) m_alleles; i++)
            for (int j = i; j <= (int) m_alleles; j++) {
                pi = n_probs[n_prob_sample][i];
                pj = n_probs[n_prob_sample][m_alleles+j];
                p = pi * pj * ((i == j) ? 1 : 2);
                Q.push(ProbID(p, i, j));
            }
    }
/*
    for (int k = 0; k < idx; k++) {
        ProbID C = Q.top();
        Q.pop();
        printf("%d prob %.6f [%d, %d]\n", k, C.dd, C.ii, C.jj);
    }
*/
    n_prob_sample++;

    if (n_prob_sample == m_context.number_of_samples) {
        n_prob_sample = 0;
    }

    return true;
}

void BgenParserGenotypeReadHelper::genotypes_end()
{
    close();
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
    rsid = "";
    split(fields, rsid_field, n_rip->rsid_sep, 3);

    if (fields[0].compare(0, 3, "chr") == 0)
        fields[0].erase(0, 3);
    if (inMap(fields[0], Input->G.chrm_set)) {
        chrm = fields[0];
        rsid = "chr" + chrm + "_" + ccpos;
    } else if (chrm == "NA" && n_rip->oxford_single_chr != "--")
        chrm = n_rip->oxford_single_chr;

    if (rsid != "") {
    } else {
/*
        if (fields[0].compare(0, 2, "rs") != 0) {
            SECTION_ERR(bad_marker_name);
            errorvf("impute2 file: odd rsid(%s) for marker name; using %s for marker name\n",
                    C(rsid_field), C(fields[0]));
            n_badname++;
        }
*/
        rsid = fields[0];
    }
}

std::ostream& BgenParser::summarise( std::ostream& o ) const
{
     std::string cmprss[] = {"uncompressed", "zlib compressed", "zstd compressed", "??? compressed" };

    msgvf("BgenParser: bgen file (");
    msgvf( m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" );
    msgvf(", ");
//  msgvf( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" );
    msgvf("%s", C(cmprss[m_context.flags & genfile::bgen::e_CompressedSNPBlocks]));
    msgvf( ") with %d samples and %d variants.\n",
           (int)m_context.number_of_samples,
           (int)m_context.number_of_variants);
    return o ;
}

// This example program reads data from a bgen file specified as the first argument
// and outputs it as a VCF file.
void BgenParser::read_input_file()
{
//  asm("int $3");
    try {
        open( ) ;

        n_rip->markers_all      = m_context.number_of_variants;
        n_rip->markers_filtered = m_context.number_of_variants;

        n_rip->people_all       = m_context.number_of_samples;
        n_rip->people_filtered  = m_context.number_of_samples;

        summarise( std::cerr ) ;

        get_sample_ids(
/*
            [this]( std::string const& id ) { n_sample.samples.push_back(id) ; }
*/
            fn0(this)
            );
        
        // Output variants
        std::string chromosome ;
        uint32_t position ;
        std::string rsid ;
        std::vector< std::string > alleles ;
    
        while( read_variant( &chromosome, &position, &rsid, &alleles )) {
            string ccpos;
            fix_marker_pos(ccpos, position);

            string rsid_name;
            validate_marker_name(rsid_name, chromosome, rsid, ccpos);
            n_rip->markers.push_back(new ImpMarker(rsid_name, chromosome, ccpos, alleles, n_rip->read_info));

            ignore_probs() ;
        }
    }
    catch( genfile::bgen::BGenError const& e ) {
        std::cerr << "!! Uh-oh, error parsing bgen file: ";
        std::cerr << e.what() << "\n";
	throw;
    }
}
