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

#ifndef READ_IMPUTE_BGEN_HH
#define READ_IMPUTE_BGEN_HH

#include <string>
#include "zlib.h"
#include "read_impute.hh"

class ReadBgen : public ReadImputed {
public:
    virtual void read_input_file();
    virtual void build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons);
};

class ReadBgenFile {
public:
    ReadBgenFile() : debug(0) {};
    ReadBgenFile(int i) : debug(i) {};
   ~ReadBgenFile() {};

    void read_input_file();

    void open(const char *path);
    void close();
    void read_bytes(unsigned char *buf, long len);
    void seek_bytes(long len);

    unsigned char  read_uchar();
    unsigned char  read_uchar(unsigned char *& zp);
    unsigned short read_ushort();
    unsigned short read_ushort(unsigned char *& zp);
    unsigned long  read_ulong();
    unsigned long  read_ulong(unsigned char *& zp);

    void read_header();
    void read_samples();
    void read_str(unsigned char *buf, long len);
    void read_embedded_uchar(Str& str, long len, unsigned char& l);
    void read_embedded_ushort(Str& str, unsigned short& len);
    void read_embedded_ulong(Str& str, unsigned long& len);

    void fix_marker_pos(std::string& ccpos);
    void validate_marker_name(std::string& rsid, const std::string& ccpos);

    void read_snpblock_header_10();
    void read_compressed_block_10();
    void read_expanded_block_10(unsigned char *bp);
    void process_10(int i);

    void read_snpblock_header_11();
    void read_compressed_block_11();
    void read_expanded_block_11();
    void process_11(int i);

    void read_snpblock_header_12();
    void read_compressed_block_12();
    void read_expanded_block_12(unsigned long DC, unsigned char *bp);
    void process_12(int i);

public:
    unsigned long offset;
    unsigned long header;
    unsigned long snps;
    unsigned long samples;
    unsigned long  rsvd;
//  char free[header-20];
    unsigned long  flags;
#define compressF  1
#define layoutFld  (0xF<<2)
#define sidF       0x80000000

    unsigned long layout;
    long badname;
    long scale;
    unsigned char *zp;
    int pass;

    struct SNPblock {
        unsigned long  N;
        unsigned short LS;
        unsigned short LR;
        unsigned short LC;
        unsigned long  pos;
        unsigned long  alleles;
        unsigned long  Lcompress;
        unsigned char  uchrm;

        std::string snp;
        std::string rsid;
        std::string chrm;
        VecsDB allele;
        unsigned char *rdata;
        unsigned char *cdata;

    } block;

    struct Sampleblock {
        unsigned length;
        unsigned count;
        VecsDB   samples;
    } sample;

    std::string impute_file;
    ReadImputed *rip;

private:
    gzFile fd;
    unsigned char ibuf[4];
    int debug;
};


class ReadBgenGenotypeReadHelper : public GenotypeReadHelper, public ReadBgenFile
{
public:
    virtual void genotypes_init();
    virtual boolean genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid, 
                                         std::string& pos, std::vector<std::string>& alleles);
    virtual void genotypes_sample_prob(ProbQ& Q);
    virtual void genotypes_end();
};


#endif
