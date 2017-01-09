/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef READ_IMPUTE_BGEN2_HH
#define READ_IMPUTE_BGEN2_HH

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <memory>
#include "genfile/types.hpp"
#include "genfile/bgen/bgen.hpp"

#include "read_impute.hh"

class ReadBgen2 : public ReadImputed {
public:
    virtual void read_input_file();
    virtual void build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons);
};

// ProbSetter is a callback object appropriate
// for passing to bgen::read_genotype_data_block().
// See the comment about bgen::read_genotype_data_block() for more info.
// Its purpose is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
class ProbSetter {
public:
    typedef std::vector< std::vector< double > > Data ;
    ProbSetter( Data* result ):
        m_result( result ),
        m_sample_i(0)
    {}
        
    // Called once allowing us to set storage.
    void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
        m_result->clear() ;
        m_result->resize( number_of_samples ) ;

        n_samples = number_of_samples;
        n_alleles = number_of_alleles;
    }
    
    // If present with this signature, called once after initialise()
    // to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
    // This enables us to set up storage for the data ahead of time.
    void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
        for( std::size_t i = 0; i < m_result->size(); ++i ) {
            m_result->at( i ).reserve( max_entries ) ;
        }
    }

    // Called once per sample to determine whether we want data for this sample
    bool set_sample( std::size_t i ) {
        m_sample_i = i ;
        // Yes, here we want info for all samples.
        return true ;
    }
    
    // Called once per sample to set the number of probabilities that are present.
    void set_number_of_entries(
        std::size_t ploidy,
        std::size_t number_of_entries,
        genfile::OrderType order_type,
        genfile::ValueType value_type
    ) {
        assert( value_type == genfile::eProbability ) ;
        m_result->at( m_sample_i ).resize( number_of_entries ) ;
        m_entry_i = 0 ;

        n_entries = number_of_entries;
        n_ploidy  = ploidy ;
        n_phased  = order_type == genfile::ePerOrderedHaplotype ||
                    order_type == genfile::ePerPhasedHaplotypePerAllele;
    }

    void set_value( uint32_t, double value ) {
        m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
    }

    void set_value( uint32_t, genfile::MissingValue value ) {
        // Here we encode missing probabilities with -1
        m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
    }

public:
    std::size_t n_samples;
    std::size_t n_alleles;
    std::size_t n_entries;  // valueCount
    std::size_t n_ploidy;
    std::size_t n_phased;
/*
    uint32_t const valueCount
            = phased
            ? (ploidy * numberOfAlleles)
            : ((numberOfAlleles == 2) ? (ploidy+1) : 
                  impl::n_choose_k( ploidy + numberOfAlleles - 1, numberOfAlleles - 1 )) ;
    uint32_t const storedValueCount = valueCount - ( phased ? ploidy : 1 ) ;

*/

private:
    Data* m_result ;
    std::size_t m_sample_i ;
    std::size_t m_entry_i ;
} ;

// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
class BgenParser;
class BgenParser {
public:    
    BgenParser():
        m_filename( "" ),
        m_state( e_NotOpen ),
        m_have_sample_ids( false ) {}

    class f0 {
    public:
        f0(BgenParser *hiss) : that(hiss) {}
        ~f0() {}

        void operator() ( std::string id )
            { that->m_sample_ids.push_back( id ) ; }
    private:
        BgenParser *that;
    };
    void open( )
    {
        // Open the stream
        m_stream = new std::ifstream( C(m_filename), std::ifstream::binary );
        if( !*m_stream ) {
            throw std::invalid_argument( m_filename ) ;
        }
        m_state = e_Open ;
        m_first   = 1;

        // Read the offset, header, and sample IDs if present.
        genfile::bgen::read_offset( *m_stream, &m_offset ) ;
        genfile::bgen::read_header_block( *m_stream, &m_context ) ;
        if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
            genfile::bgen::read_sample_identifier_block(
                *m_stream, m_context,
/*
                [this]( std::string id ) { m_sample_ids.push_back( id ) ; }
*/
                f0(this)
            ) ;
            m_have_sample_ids = true ;
        }
        
        // Jump to the first variant data block.
        m_stream->seekg( m_offset + 4 ) ;

        // We keep track of state (though it's not really needed for this implementation.)
        m_state = e_ReadyForVariant ;
    }

    void close() {
        m_stream->close();
        delete m_stream;
    }

    std::ostream& summarise( std::ostream& o ) const;
/*
{
        o << "BgenParser: bgen file ("
            << ( m_context.flags & genfile::bgen::e_v12Layout ? "v1.2 layout" : "v1.1 layout" )
            << ", "
            << ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
            << " with " 
            << m_context.number_of_samples << " samples and "
            << m_context.number_of_variants << " variants.\n" ;
        return o ;
}
*/

    std::size_t number_of_samples() const {
        return m_context.number_of_samples ;
    }

    // Report the sample IDs in the file using the given setter object
    // (If there are no sample IDs in the file, we report a dummy identifier).
    template< typename Setter >
    void get_sample_ids( Setter setter ) {
        n_sample.count = m_context.number_of_samples;
        if( m_have_sample_ids ) {
            for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
                setter( m_sample_ids[i] ) ;
            }
        } else {
            char numb[24];
            for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
                sprintf(numb, "%d", (int)(i+1));
//              setter( "(unknown_sample_" + std::to_string( i+1 ) + ")" ) ;
                setter( "(unknown_sample_" + std::string(numb) + ")" ) ;
            }
        }
    }

    class f1 {
    public:
        f1(std::vector< std::string >* &alleles) : f1_alleles(alleles) {}
        ~f1() {}

        void operator()( std::size_t n)
             { f1_alleles->resize( n ) ; }
    private:
        std::vector< std::string >* &f1_alleles;
    };

    class f2 {
    public:
        f2(std::vector< std::string >* &alleles) : f2_alleles(alleles) {}
        ~f2() {}

        std::string const &operator()( std::size_t i, std::string const& allele )
            { f2_alleles->at(i) = allele; return allele;}
    private:
        std::vector< std::string >* &f2_alleles;
    };

    // Attempt to read identifying information about a variant from the bgen file, returning
    // it in the given fields.
    // If this method returns true, data was successfully read, and it should be safe to call read_probs()
    // or ignore_probs().
    // If this method returns false, data was not successfully read indicating the end of the file.
    bool read_variant(
        std::string* chromosome,
        uint32_t* position,
        std::string* rsid,
        std::vector< std::string >* alleles
    ) {
        assert( m_state == e_ReadyForVariant ) ;
        std::string SNPID ; // read but ignored in this toy implementation
        
        if(
            genfile::bgen::read_snp_identifying_data(
                *m_stream, m_context,
                &SNPID, rsid, chromosome, position,
/*
                [&alleles]( std::size_t n ) { alleles->resize( n ) ; },
                [&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
*/
                f1(alleles),
                f2(alleles)
            )
        ) {
            m_state = e_ReadyForProbs ;
            return true ;
        } else {
            return false ;
        }
    }

    // Read genotype probability data for the SNP just read using read_variant()
    // After calling this method it should be safe to call read_variant() to fetch
    // the next variant from the file.
    void read_probs( std::vector< std::vector< double > >* probs ) {
        assert( m_state == e_ReadyForProbs ) ;
        ProbSetter setter( probs ) ;
        genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
            *m_stream,
            m_context,
            setter,
            &m_buffer1,
            &m_buffer2
        ) ;
        m_state = e_ReadyForVariant ;
        m_samples = setter.n_samples;
        m_alleles = setter.n_alleles;
        m_entries = setter.n_entries;
        m_ploidy  = setter.n_ploidy;
        m_phased  = setter.n_phased;
    }

    // Ignore genotype probability data for the SNP just read using read_variant()
    // After calling this method it should be safe to call read_variant()
    // to fetch the next variant from the file.
    void ignore_probs() {
        genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
        m_state = e_ReadyForVariant ;
    }

    void fix_marker_pos(std::string& ccpos, unsigned long pos);
    void validate_marker_name(std::string& rsid, std::string& chrm, std::string& rsid_field, const std::string& ccpos);
    void read_input_file();

    std::string m_filename ;
    ReadImputed *n_rip;
    unsigned long n_badname;

    unsigned int n_prob_sample;
    std::vector< std::vector< double > > n_probs ;

    struct Sampleblock {
        unsigned count;
        VecsDB   samples;
    } n_sample;

private:
//  std::string const m_filename ;  moved above  (public)
    std::ifstream *m_stream ;
//  genfile::bgen::Context m_context ; // moved next (protected)
    uint32_t m_offset ;
protected:
    // bgen::Context object holds information from the header block,
    // including bgen flags
    genfile::bgen::Context m_context ;
    std::size_t m_samples;
    std::size_t m_alleles;
    std::size_t m_entries;
    std::size_t m_ploidy;
    std::size_t m_phased;
    std::size_t m_first;

private:

    // We keep track of our state in the file.
    // Not strictly necessary for this implentation but makes it clear that
    // calls must be read_variant() followed by read_probs() (or ignore_probs())
    // repeatedly.
    enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
    State m_state ;

    // If the BGEN file contains samples ids, they will be read here.
    bool m_have_sample_ids ;
    std::vector< std::string > m_sample_ids ;
    
    // Buffers, these are used as working space by bgen implementation.
    std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
} ;


class BgenParserGenotypeReadHelper : public GenotypeReadHelper, public BgenParser
{
public:
    virtual void genotypes_init();
    virtual bool genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid, 
                                         std::string& pos, std::vector<std::string>& alleles);
    virtual void genotypes_skip();
    virtual bool genotypes_sample_prob(ProbQ& Q);
    virtual void genotypes_end();
};


#endif
