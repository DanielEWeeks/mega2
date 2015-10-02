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
#include "read_impute_bgen.hh"


SECTION_ERR_INIT(bad_name);
void ReadBgen::read_input_file()
{
    ReadBgenFile bgen;
    bgen.impute_file = impute_file;
    bgen.rip = this;

    bgen.read_input_file();

    SECTION_ERR_FINI(bad_name);
}

void ReadBgen::build_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons)
{
    ReadBgenGenotypeReadHelper gh;
    gh.impute_file = impute_file;
    gh.rip = this;

    build_internal_genotypes(LTop, persons, gh);
}

void ReadBgenGenotypeReadHelper::genotypes_init()
{
    open(C(impute_file));
    pass = 2;
    read_header();
    if (layout < 2) {
        block.cdata = new unsigned char[6 * samples];
        block.rdata = new unsigned char[6 * samples];
        if (layout == 0)
            scale = 10000;
        else if (layout == 1)
            scale = 32768;
//  } else if (layout == 2) {
    } else {
        errorvf("read_imputed_bgen: layout not supported %d\n", layout);
        EXIT(DATA_INCONSISTENCY);
    }
}

boolean ReadBgenGenotypeReadHelper::genotypes_marker_hdr(int mrk_idx, std::string& hmm, std::string& chrm, std::string& rsid, 
                                                         std::string& pos, vector<string>& alleles)
{
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
    chrm = block.chrm;
    fix_marker_pos(pos);

    alleles.clear();
    alleles.push_back(block.allele[0]);
    alleles.push_back(block.allele[1]);

    return true;
}

void ReadBgenGenotypeReadHelper::genotypes_sample_prob(ProbQ& Q)
{
    Token::d3 nums;
    nums[0] = ((double)read_ushort(zp)) / scale;
    nums[1] = ((double)read_ushort(zp)) / scale;
    nums[2] = ((double)read_ushort(zp)) / scale;

    Q.push(ProbID(nums[0], 1, 1));
    Q.push(ProbID(nums[1], 1, 2));
    Q.push(ProbID(nums[2], 2, 2));

}

void ReadBgenGenotypeReadHelper::genotypes_end()
{
    if (layout < 2) {
        delete [] block.cdata;
        delete [] block.rdata;
    }
}

void ReadBgenFile::read_input_file()
{
//  asm("int $3");
    unsigned long i;

    open(C(impute_file));
    pass = 1;

    read_header();

    if (layout < 2) {
        block.cdata = new unsigned char[6 * samples];
        block.rdata = new unsigned char[6 * samples];

        if ( layout == 1) {
            scale = 32768;
            for (i = 0; i < snps; i++)
                process_11(i);
        } else if (layout == 0) {
            scale = 10000;
            for (i = 0; i < snps; i++)
                process_10(i);
        }
        delete [] block.rdata;
        delete [] block.cdata;
//  } else if ( layout == 2) {
//      for (i = 0; i < snps; i++)
//          process_12(i);
    } else {
        errorvf("read_imputed_bgen: layout not supported %d\n", layout);
        EXIT(DATA_INCONSISTENCY);
    }

    close();

}


void ReadBgenFile::open(const char *path)
{
    fd = gzopen(path, "rb");
    if (fd == NULL) {
        errorvf("read_imputed_bgen: bgen file can not be opened: \"%s\"\n", path);
        EXIT(FILE_NOT_FOUND);
    }
}

void ReadBgenFile::close()
{
    long gzsz = gzclose_r(fd);
    if (gzsz == Z_OK) return;
    errorvf("read_imputed_bgen libz: gz_close_r bad return value %d\n", gzsz);
    EXIT(DATA_INCONSISTENCY);
}

void ReadBgenFile::read_bytes(unsigned char *buf, long len)
{
    long gzsz = gzread(fd, buf, len);
    if (gzsz == len) return;

    if (gzsz > 0 && gzsz < len && gzeof(fd)) {
        errorvf("read_imputed_bgen libz: gzread short return value %d\n", gzsz);
        EXIT(DATA_INCONSISTENCY);
    } else {
        errorvf("read_imputed_bgen libz: gzread bad return value %d\n", gzsz);
        EXIT(DATA_INCONSISTENCY);
    }
}

void ReadBgenFile::seek_bytes(long len)
{
    long gzsz = gzseek(fd, len, /* SEEK_CUR */ 1);

    if (gzsz < 0) {
        errorvf("read_imputed_bgen libz: gzseek bad return value %d\n", gzsz);
        EXIT(DATA_INCONSISTENCY);
    }
}

unsigned char ReadBgenFile::read_uchar()
{
    unsigned char *bp = NULL;
    return read_uchar(bp);
}

unsigned char ReadBgenFile::read_uchar(unsigned char *& bp)
{
    if (bp == 0)
        read_bytes(ibuf, 1);
    else {
        ibuf[0] = *bp++;
    }
    return ibuf[0];
}

unsigned short ReadBgenFile::read_ushort()
{
    unsigned char *bp = NULL;
    return read_ushort(bp);
}

unsigned short ReadBgenFile::read_ushort(unsigned char *& bp)
{
    unsigned short ret = 0;
    if (bp == 0)
        read_bytes(ibuf, 2);
    else {
        ibuf[0] = *bp++;
        ibuf[1] = *bp++;
    }
    unsigned char *cp = ibuf;
    for (int i = 0; i < 16; i+=8)
        ret |= *cp++ << i;
    return ret;
}

unsigned long ReadBgenFile::read_ulong()
{
    unsigned char *bp = NULL;
    return read_ulong(bp);
}

unsigned long ReadBgenFile::read_ulong(unsigned char *& bp)
{
    unsigned long ret = 0;
    if (bp == 0)
        read_bytes(ibuf, 4);
    else {
        ibuf[0] = *bp++;
        ibuf[1] = *bp++;
        ibuf[2] = *bp++;
        ibuf[3] = *bp++;
    }
    unsigned char *cp = ibuf;
    for (int i = 0; i < 32; i+=8)
        ret |= *cp++ << i;
    return ret;
}

#define BGEN (((((('b' << 8) | 'g') << 8 ) | 'e') << 8 ) | 'n')
void ReadBgenFile::read_header()
{
    offset  = read_ulong();
    header  = read_ulong();
    snps    = read_ulong();
    samples = read_ulong();
    rsvd    = read_ulong();
    if (rsvd != 0 && rsvd != BGEN) {
        errorvf("read_imputed_bgen libz: read_header bad rsvd value %d\n", rsvd);
    }
    if (header != 20)
        gzseek(fd, 4 + header - 4, 0);

    flags   = read_ulong();
    if (flags & sidF)
        read_samples();
    layout = (flags & layoutFld) >> 2;

    if (offset != header)
        gzseek(fd, offset + 4, 0);

    if (pass == 1) {
        rip->markers_all = snps;
        rip->markers_filtered = snps;

        rip->people_all = samples;
        rip->people_filtered = samples;
    }
}

void ReadBgenFile::read_samples()
{
    unsigned short size;

    sample.length = read_ulong();
    sample.count  = read_ulong();
    for (unsigned long i = 0; i < sample.count; i++) {
        string str;
        read_embedded_ushort(str, size);
        sample.samples.push_back(str);
    }
}

void ReadBgenFile::read_str(unsigned char *buf, long len)
{
    read_bytes(buf, len);
}

void ReadBgenFile::read_embedded_uchar(Str& str, long len, unsigned char& l)
{
    read_bytes(ibuf, 1);
    l = ibuf[0];
    if (len == 0) len = l;
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, l);
    delete [] tmp;
}

void ReadBgenFile::read_embedded_ushort(Str& str, unsigned short& len)
{
    len  = read_ushort();
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, len);
    delete [] tmp;
}

void ReadBgenFile::read_embedded_ulong(Str& str, unsigned long& len)
{
    len  = read_ulong();
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, len);
    delete [] tmp;
}

void ReadBgenFile::fix_marker_pos(string& ccpos) {
    char cpos[10];
    sprintf(cpos, "%ld", block.pos);
    ccpos = string(cpos);
}

void ReadBgenFile::validate_marker_name(string& rsid, const string& ccpos) {
    VecsDB fields;
    rsid = "";
    split(fields, block.rsid, rip->rsid_sep, 3);
    if (fields[0].compare(0, 3, "chr") == 0)
        fields[0].erase(0, 3);
    if (inMap(fields[0], Input->G.chrm_set)) {
        block.chrm = fields[0];
        rsid = "chr" + block.chrm + "_" + ccpos;
    } else if (block.chrm == "NA" && rip->oxford_single_chr != "--")
        block.chrm = rip->oxford_single_chr;

    if (rsid != "") {
    } else {
        if (fields[0].compare(0, 2, "rs") != 0) {
            SECTION_ERR(bad_name);
            errorvf("impute2 file: odd rsid(%s) for marker name; using %s for marker name\n",
                    C(block.rsid), C(fields[0]));
            badname++;
        }
        rsid = fields[0];
    }
}

void ReadBgenFile::read_snpblock_header_10()
{
    unsigned char uchar;

    block.allele.clear();

    if (layout < 2)
        block.N   = read_ulong();
    read_bytes(ibuf, 1);
    long fill = ibuf[0];

    read_embedded_uchar(block.snp,  fill, uchar);
    block.LS = uchar;
    read_embedded_uchar(block.rsid, fill, uchar);
    block.LR = uchar;

    read_bytes(ibuf, 1);
    block.uchrm = ibuf[0];
    char uchr[16];
    sprintf(uchr, "%d", block.uchrm);
    block.LC = (unsigned short)strlen(uchr);
    block.chrm = string(uchr);

    block.pos = read_ulong();

    block.alleles = 2;
    Str str, str2;
    read_embedded_uchar(str, 0, uchar);
    block.allele.push_back(str);

    read_embedded_uchar(str2, 0, uchar);
    block.allele.push_back(str2);
    Vecs alleles(2);
    alleles.clear();
    alleles.push_back(str);
    alleles.push_back(str2);

    if (pass == 1) {
        string ccpos;
        fix_marker_pos(ccpos);

        string rsid;
        validate_marker_name(rsid, ccpos);
        rip->markers.push_back(new ImpMarker(rsid, block.chrm, ccpos, alleles, 0));
    }
}

void ReadBgenFile::read_compressed_block_10()
{
    unsigned long dest_len = 6 * block.N;
    unsigned long CB = read_ulong();

    if (pass == 1) {
        seek_bytes(CB);
        return;
    }
//  asm("int $3");

    read_bytes(block.rdata, CB);
    int z_stat = uncompress(block.cdata, &dest_len, block.rdata, CB);
    if (z_stat == Z_MEM_ERROR) {
        errorvf("read_imputed_bgen libz: uncompress mem_error %d\n", z_stat);
    } else if (z_stat != Z_OK) {
        errorvf("read_imputed_bgen libz: uncompress misc error %d\n", z_stat);
    } else if (dest_len != 6 * block.N) {
        errorvf("read_imputed_bgen libz: uncompress short read %d vs %d\n", dest_len, 6 * block.N);
    }

    read_expanded_block_10(block.cdata);

}

void ReadBgenFile::read_expanded_block_10(unsigned char *bp = NULL)
{
    unsigned long dest_len = 6 * block.N;

    if (pass == 1) {
        seek_bytes(dest_len);
        return;
    }

    if (pass == 2)
        return;
/* never gets here
    double num[3];
    for (unsigned long i = 0; i < block.N; i++) {
        num[0] = ((double)read_ushort(bp)) / scale;
        num[1] = ((double)read_ushort(bp)) / scale;
        num[2] = ((double)read_ushort(bp)) / scale;

        if (debug) {
            cout << "sample#" << i << " ";
            cout << num[0] << " ";
            cout << num[1] << " ";
            cout << num[2] << "\n";
        }
    }
*/
}

void ReadBgenFile::process_10(int snpn)
{
    read_snpblock_header_10();

    if (debug) {
        cout << "snp#" << snpn << " ";
        cout << block.rsid << " ";
        cout << block.chrm << ":" << block.pos << " ";
        for (unsigned long i = 0; i < block.alleles; i++)
            cout << block.allele[i] << " ";
        cout << "\n";
    }
    if ( (flags & compressF) > 0)
        read_compressed_block_10();
    else
        read_expanded_block_10();
}


void ReadBgenFile::read_snpblock_header_11()
{
    unsigned long size;

    block.allele.clear();

    if (layout < 2)
        block.N   = read_ulong();

    read_embedded_ushort(block.snp,  block.LS);
    read_embedded_ushort(block.rsid, block.LR);
    read_embedded_ushort(block.chrm, block.LC);

    block.pos = read_ulong();
    if (layout < 2)
        block.alleles = 2;
    else
        block.alleles = read_ushort();

    for (unsigned long i = 0; i < block.alleles; i++) {
        string str;
        read_embedded_ulong(str, size);
        block.allele.push_back(str);
    }

    if (pass == 1) {
        string ccpos;
        fix_marker_pos(ccpos);

        string rsid;
        validate_marker_name(rsid, ccpos);
        rip->markers.push_back(new ImpMarker(rsid, block.chrm, ccpos, block.allele, 0));
    }
}

void ReadBgenFile::read_compressed_block_11()
{
    read_compressed_block_10();
}

void ReadBgenFile::read_expanded_block_11()
{
    read_expanded_block_10();
}

void ReadBgenFile::process_11(int snpn)
{
    read_snpblock_header_11();
    if (debug) {
        cout << "snp#" << snpn << " ";
        cout << block.rsid << " ";
        cout << block.chrm << ":" << block.pos << " ";
        for (unsigned long i = 0; i < block.alleles; i++)
            cout << block.allele[i] << " ";
        cout << "\n";
    }

    if ( (flags & compressF) > 0)
        read_compressed_block_11();
    else
        read_expanded_block_11();
}

void ReadBgenFile::read_compressed_block_12()
{
//  asm("int $3");
    unsigned long CB = read_ulong();
    unsigned long DC = read_ulong();
    unsigned long dest_len = DC;
    block.cdata = new unsigned char[DC];
    CB -= 4;  //subtract out 4 for DC
    block.rdata = new unsigned char[CB];
    read_bytes(block.rdata, CB);
    int z_stat = uncompress(block.cdata, &dest_len, block.rdata, CB);
    if (z_stat == Z_MEM_ERROR) {
        errorvf("read_imputed_bgen libz: uncompress mem_error %d\n", z_stat);
    } else if (z_stat != Z_OK) {
        errorvf("read_imputed_bgen libz: uncompress misc error %d\n", z_stat);
    } else if (dest_len != 6 * block.N) {
        errorvf("read_imputed_bgen libz: uncompress short read %d vs %d\n", dest_len, 6 * block.N);
    }

    read_expanded_block_12(DC, block.cdata);

    delete [] block.rdata;
    delete [] block.cdata;
}

void ReadBgenFile::read_expanded_block_12(unsigned long DC, unsigned char *bp = NULL)
{
#ifdef SOMEDAY
    unsigned long Nx, allelesx, min, max, phase, bits;
    unsigned char *ploidy;

    Nx       = read_ulong(bp);
    allelesx = read_ushort(bp);
    min      = read_uchar(bp);
    max      = read_uchar(bp);
    ploidy   = new unsigned char [ Nx ];
    read_bytes(ploidy, Nx);
    phase    = read_uchar(bp);
    bits     = read_uchar(bp);
#endif
}

void ReadBgenFile::process_12(int snpn)
{
    read_snpblock_header_11();
    if (debug) {
        cout << "snp#" << snpn << " ";
        cout << block.rsid << " ";
        cout << block.chrm << ":" << block.pos << " ";
        for (unsigned long i = 0; i < block.alleles; i++)
            cout << block.allele[i] << " ";
        cout << "\n";
    }
    if ( (flags & compressF) > 0)
        read_compressed_block_12();
    else {
        unsigned long DC = read_ulong();  // not 6 * N anymore ...
        read_expanded_block_12(DC);
    }
}
