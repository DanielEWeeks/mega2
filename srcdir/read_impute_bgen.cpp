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

#include "input.hh"
#include "read_impute.hh"
#include "read_impute_bgen.hh"

SECTION_ERR_INIT(bad_marker_name);
void ReadBgen::do_init(Input_BGEN *inp)
{
    
    asm("int $3");

    Input_Impute *ii = inp;  // up chain
    ReadImputed& ri  = ii->Obj;
    rip  = &ri;

    ri.input = inp;
    ri.files(*inp->input_files.bedfl, *inp->input_files.pedfl);
    bgen_file = *inp->input_files.bedfl;

    ri.show_settings();

    read_bgen_file();
    SECTION_ERR_FINI(bad_marker_name);
    if (badname > 0) {
        errorvf("read_imputed_bgen_file: can not determine chromosomes for %d markers\n",
                badname);
        EXIT(DATA_INCONSISTENCY);
    }

    rip->read_info = false;  // this should be default ... but not for win mvc
    if (!rip->info_file.empty()) {
        ifstream infs;
        infs.open(rip->info_file.c_str());
        if (! infs.is_open() ) {
            errorvf("read_imputed_file: Info file can not be opened: \"%s\"\n", C(rip->info_file));
            EXIT(FILE_NOT_FOUND);
        }
        infs.close();
        rip->read_info = true;
    }

    if (ri.read_info) ri.read_info_file();

    ri.check_indelsNdups();

    ri.read_sample_file();
}

linkage_ped_top *ReadBgen::do_ped(linkage_locus_top *LTop)
{
    int num_peds = 0;

    annotated_ped_rec *persons = rip->build_impute2_ped(LTop, &num_peds);

    build_bgen_genotypes(LTop, persons);

    linkage_ped_top *Top;
    Top = mk_ped_top(persons, rip->people.size(), LTop, num_peds,
                     /*untyped*/ 0, /*totaltyped*/ rip->people_filtered,
                     /*groups*/ NULL, 0, 0, 
                     /*num_err*/0, 1);

    return Top;

}

void ReadBgen::read_bgen_file()
{
    unsigned long i;

    open(bgen_file);
    pass = 1;

    read_header();

    if (flags & sidF)
        read_samples();

    layout = (flags & layoutFld) >> 2;
    if (layout < 2) {
        block.cdata = new unsigned char[6 * samples];
        block.rdata = new unsigned char[6 * samples];

        if ( layout == 1)
            for (i = 0; i < snps; i++)
                process_11(i);
        else if (layout == 0)
            for (i = 0; i < snps; i++)
                process_10(i);

        delete [] block.rdata;
        delete [] block.cdata;
    } else if ( layout == 2) {
        for (i = 0; i < snps; i++)
            process_12(i);
    } else {
//      ERROR
    }

    close();

}


void ReadBgen::open(const char *path)
{
    fd = gzopen(path, "rb");
    if (fd == NULL) {
//        ERROR
    }
}

void ReadBgen::close()
{
    long gzsz = gzclose_r(fd);
    if (gzsz == Z_OK) return;

    {
//      ERROR
    }
}

void ReadBgen::read_bytes(unsigned char *buf, long len)
{
    long gzsz = gzread(fd, buf, len);
    if (gzsz == len) return;

    if (gzsz > 0 && gzsz < len && gzeof(fd)) {
//short ??
    } else {
//      ERROR
    }
}

void ReadBgen::seek_bytes(long len)
{
    long gzsz = gzseek(fd, len, /* SEEK_CUR */ 1);

    if (gzsz > 0 && gzsz < len && gzeof(fd)) {
//short ??
    } else {
//      ERROR
    }
}

unsigned long ReadBgen::read_uchar()
{
    unsigned char *zp = NULL;
    return read_uchar(zp);
}

unsigned long ReadBgen::read_uchar(unsigned char *& zp)
{
    if (zp == 0)
        read_bytes(ibuf, 1);
    else {
        ibuf[0] = *zp++;
    }
    return ibuf[0];
}

unsigned long ReadBgen::read_ushort()
{
    unsigned char *zp = NULL;
    return read_ushort(zp);
}

unsigned long ReadBgen::read_ushort(unsigned char *& zp)
{
    unsigned long ret = 0;
    if (zp == 0)
        read_bytes(ibuf, 2);
    else {
        ibuf[0] = *zp++;
        ibuf[1] = *zp++;
    }
    unsigned char *cp = ibuf;
    for (int i = 0; i < 16; i+=8)
        ret |= *cp++ << i;
    return ret;
}

unsigned long ReadBgen::read_ulong()
{
    unsigned char *zp = NULL;
    return read_ulong(zp);
}

unsigned long ReadBgen::read_ulong(unsigned char *& zp)
{
    unsigned long ret = 0;
    if (zp == 0)
        read_bytes(ibuf, 4);
    else {
        ibuf[0] = *zp++;
        ibuf[1] = *zp++;
        ibuf[2] = *zp++;
        ibuf[3] = *zp++;
    }
    unsigned char *cp = ibuf;
    for (int i = 0; i < 32; i+=8)
        ret |= *cp++ << i;
    return ret;
}

#define BGEN (((((('b' << 8) | 'g') << 8 ) | 'e') << 8 ) | 'n')
void ReadBgen::read_header()
{
    offset  = read_ulong();
    header  = read_ulong();
    snps    = read_ulong();
    samples = read_ulong();
    rsvd    = read_ulong();
    if (rsvd != 0 && rsvd != BGEN) {
// error
    }
    if (header != 20)
        gzseek(fd, 4 + header - 4, 0);

    flags   = read_ulong();

    if (offset != header)
        gzseek(fd, offset + 4, 0);

    if (pass == 1) {
        rip->markers_all = snps;
        rip->markers_filtered = snps;

        rip->people_all = samples;
        rip->people_filtered = samples;
    }
}

void ReadBgen::read_samples()
{
    unsigned short size;

    sample.length = read_ulong();
    sample.count  = read_ulong();
    for (unsigned long i; i < sample.count; i++) {
        string str;
        read_embedded_ushort(str, size);
        sample.samples.push_back(str);
    }
}

void ReadBgen::read_str(unsigned char *buf, long len)
{
    read_bytes(buf, len);
}

void ReadBgen::read_embedded_uchar(Str& str, long len, unsigned char& l)
{
    read_bytes(ibuf, 1);
    l = ibuf[0];
    if (len == 0) len = l;
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, l);
    delete [] tmp;
}

void ReadBgen::read_embedded_ushort(Str& str, unsigned short& len)
{
    len  = read_ushort();
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, len);
    delete [] tmp;
}

void ReadBgen::read_embedded_ulong(Str& str, unsigned long& len)
{
    len  = read_ulong();
    char *tmp = new char [len];
    read_bytes((unsigned char *)tmp, len);
    str  = string(tmp, 0, len);
    delete [] tmp;
}

void ReadBgen::fix_marker_pos(string& ccpos) {
    char cpos[10];
    sprintf(cpos, "%ld", block.pos);
    ccpos = string(cpos);
}

void ReadBgen::validate_marker_name(string& rsid, const string& ccpos) {
    VecsDB fields;
    split(fields, block.rsid, ":", 3);

    if (fields.size() == 1) {
        rsid = fields[0];
        if (rsid == "." || rsid == "NA" || rsid == "na") {
            block.chrm = rip->oxford_single_chr;
            rsid = "chr" + rip->oxford_single_chr + "_" + ccpos;
        }
    } else {
        if (fields[0].compare(0, 3, "chr") == 0)
            fields[0].erase(0, 3);
        if (inMap(fields[0], Input->G.chrm_set)) {
            block.chrm = fields[0];
            rsid = "chr" + fields[0] + "_" + ccpos;
        } else if (fields[0].compare(0, 2, "rs") == 0 && rip->oxford_single_chr != "--") {
            block.chrm = rip->oxford_single_chr;
            rsid = fields[0];
        } else {
            SECTION_ERR(bad_marker_name);
            errorvf("impute2 file: bad marker name %s\n", C(block.rsid));
            badname++;
        }
    }
}

void ReadBgen::read_snpblock_header_10()
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
    block.LC = strlen(uchr);
    block.chrm = string(uchr);

    block.pos = read_ulong();

    block.alleles = 2;
    Str str, str2;
    read_embedded_uchar(str, 0, uchar);
    block.allele.push_back(str);

    read_embedded_uchar(str2, 0, uchar);
    block.allele.push_back(str2);

    if (pass == 1) {
//      for v1.0 & v1.1
//      some checks please ...
        string ccpos;
        fix_marker_pos(ccpos);

        string rsid;
        validate_marker_name(rsid, ccpos);
        rip->markers.push_back(new ImpMarker(rsid, block.chrm, ccpos, str, str2, 0));
    }
}

void ReadBgen::read_compressed_block_10(long scale = 10000)
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
//        NO MEM
    } else if (z_stat != Z_OK) {
//        CORRUPTION
    } else if (dest_len != 6 * block.N) {
//        short fill
    }

    read_expanded_block_10(scale, block.cdata);

}

void ReadBgen::read_expanded_block_10(long scale = 10000, unsigned char *zp = NULL)
{
    unsigned long dest_len = 6 * block.N;

    if (pass == 1) {
        seek_bytes(dest_len);
        return;
    }

    if (pass == 2)
        return;

//  long AA, AB, BB;
    double num[3];
    for (unsigned long i = 0; i < block.N; i++) {
        num[0] = read_ushort(zp);
        num[1] = read_ushort(zp);
        num[2] = read_ushort(zp);

/*
        AA = num[0];
        AB = num[1];
        BB = num[2];
*/
        num[0] /= scale;
        num[1] /= scale;
        num[2] /= scale;

        if (debug) {
            cout << "sample#" << i << " ";
            cout << num[0] << " ";
            cout << num[1] << " ";
            cout << num[2] << "\n";
        }
    }
}

void ReadBgen::process_10(int i)
{
    read_snpblock_header_10();

    if (debug) {
        cout << "snp#" << i << " ";
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


void ReadBgen::read_snpblock_header_11()
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
//      for v1.0 & v1.1
//      some checks please ...
        string ccpos;
        fix_marker_pos(ccpos);

        string rsid;
        validate_marker_name(rsid, ccpos);
        rip->markers.push_back(new ImpMarker(rsid, block.chrm, ccpos, block.allele[0], block.allele[1], 0));
    }
}

void ReadBgen::read_compressed_block_11()
{
    read_compressed_block_10(32768);
}

void ReadBgen::read_expanded_block_11()
{
    read_expanded_block_10(32768);
}

void ReadBgen::process_11(int i)
{
    read_snpblock_header_11();
    if (debug) {
        cout << "snp#" << i << " ";
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

void ReadBgen::read_compressed_block_12()
{
    asm("int $3");
    unsigned long CB = read_ulong();
    unsigned long DC = read_ulong();
    unsigned long dest_len = DC;
    block.cdata = new unsigned char[DC];
    CB -= 4;  //subtract out 4 for DC
    block.rdata = new unsigned char[CB];
    read_bytes(block.rdata, CB);
    int z_stat = uncompress(block.cdata, &dest_len, block.rdata, CB);
    if (z_stat == Z_MEM_ERROR) {
//        NO MEM
    } else if (z_stat != Z_OK) {
//        CORRUPTION
    } else if (dest_len != 6 * block.N) {
//        short fill
    }

    read_expanded_block_12(DC, block.cdata);

    delete [] block.rdata;
    delete [] block.cdata;
}

void ReadBgen::read_expanded_block_12(unsigned long DC, unsigned char *zp = NULL)
{
#ifdef SOMEDAY
    unsigned long Nx, allelesx, min, max, phase, bits;
    unsigned char *ploidy;

    Nx       = read_ulong(zp);
    allelesx = read_ushort(zp);
    min      = read_uchar(zp);
    max      = read_uchar(zp);
    ploidy   = new unsigned char [ Nx ];
    read_bytes(ploidy, Nx);
    phase    = read_uchar(zp);
    bits     = read_uchar(zp);
#endif
}

void ReadBgen::process_12(int i)
{
    read_snpblock_header_11();
    if (debug) {
        cout << "snp#" << i << " ";
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

void ReadBgen::build_bgen_genotypes(linkage_locus_top *LTop, annotated_ped_rec *persons)
{
    asm("int $3");
    int dbg = 0;

//    string line;

    string hmm, rsid, pos;
    const char *A, *B;
//  double nums[3];
    Token::d3 nums;

    open(bgen_file);
    pass = 2;
    read_header();
    int scale;
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

/*
    ifs.open(impute_file);
    if (! ifs.is_open() ) {
        errorvf("read_imputed_genotype_file: Can not open \"%s\" file\n",
                impute_file);
        EXIT(FILE_NOT_FOUND);
    }
*/
//    int line_n = 0;
//    Token token(3);
    ImpMarker   *mp;

    int   mrk_idx = rip->sample_file_hdr2b.size() - 5 -1;
    linkage_locus_rec *locus;

    char *callele1;
    char *callele2;
    char *callele0  = canonical_allele(C("0"));

    Tod tod_gen("impute genotypes");
    Tod tod_per(30);
    Tod tod_im_line(30);
    Tod tod_im_line_cpy(30);
    SECTION_ERR_INIT(genotype_missing_fraction);
    SECTION_ERR(genotype_missing_fraction);
    warnvf("%8s %8s %8s %8s  %s\n         %8s %8s %8s %8s  %s\n         %8s %8s %8s %8s  %s\n",
           "untyped", "uncertain", "good", "geno-", "Marker + chr:pos",
           "marker", "hard", "hard", "typing", "",
           "0/0", "call", "call", "rate", "");
    for (Vecmarkerpp mpp = (rip->markers).begin(); mpp != (rip->markers).end(); mpp++) {
        tod_im_line.reset();
        tod_im_line_cpy.reset();

        if (layout == 0) {
            process_10(mrk_idx);
        } else if (layout == 1) {
            process_11(mrk_idx);
        } else if (layout == 2) {
            process_12(mrk_idx);
        }

        hmm = "---";
        rsid = block.rsid;
        fix_marker_pos(pos);
        tod_im_line("impute read line  ");
/*
        getline(ifs, line);
        tod_im_line("impute read line  ");
        if (ifs.eof()) break;
        token.set(line);
        tod_im_line("impute set line   ");

//	ifs >> hmm; // +A+B+ nums() all read via >> @ 14.62 sec
//                               vector<double> nums  12.26
//                               double nums[3]       11.59
        token.more(hmm);
        token.more(rsid);
        token.more(pos);
*/
//        line_n++;
        mp = *mpp;

        if (pos != mp->pos) {
            errorvf("internal error: impute_bgen_file (\"%s\") second pass does not match first pass at marker #%d %s\n",
                    bgen_file, mrk_idx, C(rsid));
                EXIT(DATA_INCONSISTENCY);
        }
        if (mp->skip) {
            continue;
        }
        mrk_idx++;  // hence the -1 above
        locus = &LTop->Locus[mrk_idx];

        A = C(block.allele[0]);
        B = C(block.allele[1]);
/*
        token.more(A);
        token.more(B);
*/
        callele1    = canonical_allele(A);
        callele2    = canonical_allele(B);

        int i;
        unsigned char *zp;
        int sam = 0;
        int zero = 0, uncertain = 0, good = 0;
        if ((flags & compressF) > 0)
            zp = block.cdata;
        else
            zp = 0;

        tod_per.reset();
        annotated_ped_rec *entry = persons;
	for(int p = 0 ; p < rip->people_filtered; p++, entry++) {
            sam++;

            nums[0] = ((double)read_ushort(zp)) / scale;
            nums[1] = ((double)read_ushort(zp)) / scale;
            nums[2] = ((double)read_ushort(zp)) / scale;

/*
            token.getDC(nums);
*/
            
            if (nums[0] > nums[1]) {
                i = 0;
            } else {
                i = 1;
            }
            if (nums[2] > nums[i]) {
                i = 2;
            }

            if (nums[i] == 0) {
                i = 3;
                zero++;
            } else if (nums[i] <= 1 - rip->hard_call_uncertainty) {
                i = 3;
                uncertain++;
            } else {
                good++;
            }

            if (dbg) {
                cout << "#";
//              cout << line_n << "@";
                cout << sam << ": ";
                cout << mp->name << " ";
                cout << nums[0] << " " << nums[1] << " " << nums[2] << "  " << i << "\n";
            }

            if (i == 0)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele1);
            else if (i == 1)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele1, callele2);
            else if (i == 2)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele2, callele2);
            else if (i == 3)
                set_2Ralleles(entry->marker, mrk_idx, locus, callele0, callele0);
        }
        tod_per("impute person loop");
        if (good < rip->genotype_missing_fraction * rip->people_filtered) {
/* some day
??          mp->skip = true;
            LTop->Marker[mrk_idx].chromosome = MISSING_CHROMO;
            locus->??
            skip_count++;
*/
            SECTION_ERR(genotype_missing_fraction);
            warnvf("%8d %8d %8d %8.3f  %s %s:%s\n",
                   zero, uncertain, good, ((double)good)/rip->people_filtered,
                   C(mp->name), C(mp->chr), C(mp->pos));
        }
    }
    tod_gen();
    SECTION_ERR_FINI(genotype_missing_fraction);

    if (layout < 2) {
        delete [] block.cdata;
        delete [] block.rdata;
    }
}



#ifdef TEST
main(int argc, char **argv)
{

    ReadBgen bgen;
    int i;

    bgen.open(argv[1]);

//  bgen.open("/Users/rbaron/test.bgen");
    asm("int $3");

    bgen.read_header();

    if (bgen.flags & sidF)
        bgen.read_samples();

    bgen.layout = (bgen.flags & layoutFld) >> 2;
    if (bgen.layout < 2) {
        bgen.block.cdata = new unsigned char[6 * bgen.samples];
        bgen.block.rdata = new unsigned char[6 * bgen.samples];

        if ( bgen.layout == 1)
            for (i = 0; i < bgen.snps; i++)
                bgen.process_11(i);
        else if (bgen.layout == 0)
            for (i = 0; i < bgen.snps; i++)
                bgen.process_10(i);

        delete [] bgen.block.rdata;
        delete [] bgen.block.cdata;
    } else if ( bgen.layout == 2) {
        for (i = 0; i < bgen.snps; i++)
            bgen.process_12(i);
    } else {
//      ERROR
    }
}
#endif
