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

#include <string>
#include "Rcpp.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::Matrix<STRSXP> getgenotypes_1(NumericVector locus_arg,
                                    NumericVector hocus_arg,
                                    List          genotype_arg,
                                    List          allele_arg,
                                    List          markerscheme_arg,
                                    NumericVector phenocnt_arg)
{
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    int locus_size = loci.size();
    Rcpp::NumericVector hoci(hocus_arg);
    int hocus_size = hoci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // geno[, 2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::List markerschemes(markerscheme_arg);
    Rcpp::IntegerVector allele1_map(markerschemes[2]);   // markerscheme_table[, 3]
    Rcpp::IntegerVector allele2_map(markerschemes[3]);   // markerscheme_table[, 4]

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, locus_size);

    if (locus_size != hocus_size) {
        Rf_error("First vector arguments should be the same length, but are %d vs %d\\n",
                   locus_size, hocus_size);
        return mtx;
    }

    int locus, hocus, marker, byte, a1map, a2map, t0;
    for (int i = 0; i < locus_size; i++) {
        locus = loci[i];
        hocus = hoci[i];
        marker = hocus - pheno;

        byte = marker / 4;
        marker = marker - 4 * byte;

        if (debug) Rprintf("locus %d, hocus %d, pheno: %d, marker: %d, byte %d, offset %d\\n",
                           locus, hocus, pheno, hocus-pheno, byte, marker);

        a1map = allele1_map(locus - pheno);
        a2map = allele2_map(locus - pheno);
        if (a2map == 0) a2map = a1map;
//
// allele_table is a single data.frame
        Rcpp::CharacterVector aAlleleName(allele[1]);
        std::string allele1(aAlleleName[2 * locus + a1map - 1]);
        std::string allele2(aAlleleName[2 * locus + a2map - 1]);
        Rcpp::IntegerVector aindexX(allele[3]);
        if (debug) Rprintf("indexX: %d %d\\n", aindexX[2 * locus], aindexX[2 * locus+1]);
//
// allele_table is list of data.frames
//      Rcpp::List alleleList(allele[locus]);
//      Rcpp::CharacterVector aAlleleName(alleleList[1]);
//      std::string allele1(aAlleleName[a1map - 1]);
//      std::string allele2(aAlleleName[a2map - 1]);
//      Rcpp::IntegerVector aindexX(alleleList[3]);
//      if (debug) Rprintf("indexX: %d %d\\n", aindexX[a1map - 1], aindexX[a2map - 1]);
//
        if (debug) Rprintf("allele%d/%d: %s%s; ", a1map, a2map,
                              allele1.c_str(), allele2.c_str());

        decode_allele[0] = allele1 + allele1;
        decode_allele[1] = "00";
        decode_allele[2] = allele1 + allele2;
        decode_allele[3] = allele2 + allele2;

        for (int j = 0; j < genotype_sample_size; j++) {

            Rcpp::RawVector rv(genotype_sample[j]);
            t0 = (rv.size() != 0) ? rv[byte]: 0x55;
/*
          if (debug && j <= 3)
              Rprintf("byte %d, t0 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                       byte, t0,
                       rv(byte-6), rv(byte-5), rv(byte-4), rv(byte-3), rv(byte-2), rv(byte-1),
                       rv(byte-0), rv(byte+1), rv(byte+2), rv(byte+3), rv(byte+4), rv(byte+5));
*/
          if (debug && j <= 3) Rprintf("byte %d, t0 %x \\n",  byte, t0);
          if (marker == 0) {
            mtx(j, i) = decode_allele[(t0 & 0x03) >> 0];
          } else if (marker == 1) {
            mtx(j, i) = decode_allele[(t0 & 0x0c) >> 2];
          } else if (marker == 2) {
            mtx(j, i) = decode_allele[(t0 & 0x30) >> 4];
          } else if (marker == 3) {
            mtx(j, i) = decode_allele[(t0 & 0xc0) >> 6];
          }
        }
    }
    return mtx;
}

// [[Rcpp::export]]
Rcpp::Matrix<STRSXP> getgenotypes_2(NumericVector locus_arg,
                                    List          genotype_arg,
                                    List          allele_arg,
                                    NumericVector phenocnt_arg)
{
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    int locus_size = loci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // geno[, 2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<int> decode_allele(4);

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::CharacterMatrix mtx(genotype_sample_size, locus_size);

    int locus, marker;
    for (int i = 0; i < locus_size; i++) {
        locus = loci[i];
        marker = locus - pheno;

        if (debug) Rprintf("locus %d, pheno: %d, marker: %d\\n",
                           locus, pheno, marker);

        Rcpp::List alleleList(allele[locus]);
        Rcpp::CharacterVector aAlleleName(alleleList[1]);

        int allele1, allele2;
        for (int j = 0; j < genotype_sample_size; j++) {

            Rcpp::RawVector rv(genotype_sample[j]);
            if (rv.size() != 0) {
                allele1 = rv[2 * marker];
                allele2 = rv[2 * marker + 1];
            } else {
                allele1 = allele2 = 0;  // indicate empty
            }

/*
            if (debug && j <= 3)
                Rprintf("marker %d, t0 %x t1 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                         marker, allele1, allele2,
                         rv(marker-6), rv(marker-5), rv(marker-4), rv(marker-3), rv(marker-2), rv(marker-1),
                         rv(marker-0), rv(marker+1), rv(marker+2), rv(marker+3), rv(marker+4), rv(marker+5));
*/
            if (debug && j <= 3) Rprintf("marker %d, allele1 %x, allele2 %x \\n",  marker, allele1, allele2);

            std::string allele1s(aAlleleName[allele1 > 0 ? allele1 - 1 : 0]);
            std::string allele2s(aAlleleName[allele2 > 0 ? allele2 - 1 : 0]);

            mtx(j, i) = ((allele1 != 0) ? allele1s : "0") + ((allele2 != 0) ? allele2s : "0");

        }
    }
    return mtx;
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix getgenotypesraw_1(NumericVector locus_arg,
                                      NumericVector hocus_arg,
                                      List          genotype_arg,
                                      List          allele_arg,
                                      List          markerscheme_arg,
                                      NumericVector phenocnt_arg)
{
   int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    int locus_size = loci.size();
    Rcpp::NumericVector hoci(hocus_arg);
    int hocus_size = hoci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // geno[, 2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<int> decode_allele(4);

    Rcpp::List markerschemes(markerscheme_arg);
    Rcpp::IntegerVector allele1_map(markerschemes[2]);   // markerscheme_table[, 3]
    Rcpp::IntegerVector allele2_map(markerschemes[3]);   // markerscheme_table[, 4]

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::IntegerMatrix mtx(genotype_sample_size, locus_size);

    if (locus_size != hocus_size) {
        Rf_error("First vector arguments should be the same length, but are %d vs %d\\n",
                   locus_size, hocus_size);
        return mtx;
    }

    int locus, hocus, marker, byte, a1map, a2map, t0;
    for (int i = 0; i < locus_size; i++) {
        locus = loci[i];
        hocus = hoci[i];
        marker = hocus - pheno;

        byte = marker / 4;
        marker = marker - 4 * byte;

        if (debug) Rprintf("locus %d, hocus %d, pheno: %d, marker: %d, byte %d, offset %d\\n",
                           locus, hocus, pheno, hocus-pheno, byte, marker);

        a1map = allele1_map(locus - pheno);
        a2map = allele2_map(locus - pheno);
        if (a2map == 0) a2map = a1map;
//
// allele_table is a single data.frame
        Rcpp::IntegerVector aAllele(allele[3]);
        int allele1 = (aAllele[2 * locus + a1map - 1]);
        int allele2 = (aAllele[2 * locus + a2map - 1]);
        Rcpp::IntegerVector aindexX(allele[3]);
        if (debug) Rprintf("indexX: %d %d\\n", aindexX[2 * locus], aindexX[2 * locus+1]);
//
// allele_table is list of data.frames
//      Rcpp::List alleleList(allele[locus]);
//      Rcpp::IntegerVector aAlleleIndex(alleleList[3]);
//      int allele1(aAlleleIndex[a1map - 1]);
//      int allele2(aAlleleIndex[a2map - 1]);
//      Rcpp::IntegerVector aindexX(alleleList[3]);
//      if (debug) Rprintf("indexX: %d %d\\n", aindexX[a1map - 1], aindexX[a2map - 1]);
//
        if (debug) Rprintf("allele%d/%d: %d%d; ", a1map, a2map,
                              allele1, allele2);


        decode_allele[0] = (allele1 << 16) | allele1;
        decode_allele[1] = 0;
        decode_allele[2] = (allele1 << 16) | allele2;
        decode_allele[3] = (allele2 << 16) | allele2;

        for (int j = 0; j < genotype_sample_size; j++) {

            Rcpp::RawVector rv(genotype_sample[j]);
            t0 = (rv.size() != 0) ? rv[byte]: 0x55;
/*
          if (debug && j <= 3)
              Rprintf("byte %d, t0 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                       byte, t0,
                       rv(byte-6), rv(byte-5), rv(byte-4), rv(byte-3), rv(byte-2), rv(byte-1),
                       rv(byte-0), rv(byte+1), rv(byte+2), rv(byte+3), rv(byte+4), rv(byte+5));
*/
          if (debug && j <= 3) Rprintf("byte %d, t0 %x \\n",  byte, t0);
          if (marker == 0) {
            mtx(j, i) = decode_allele[(t0 & 0x03) >> 0];
          } else if (marker == 1) {
            mtx(j, i) = decode_allele[(t0 & 0x0c) >> 2];
          } else if (marker == 2) {
            mtx(j, i) = decode_allele[(t0 & 0x30) >> 4];
          } else if (marker == 3) {
            mtx(j, i) = decode_allele[(t0 & 0xc0) >> 6];
          }
        }
    }
    return mtx;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getgenotypesraw_2(NumericVector locus_arg,
                                      List          genotype_arg,
                                      List          allele_arg,
                                      NumericVector phenocnt_arg)
{
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    int locus_size = loci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // geno[, 2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<int> decode_allele(4);

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::IntegerMatrix mtx(genotype_sample_size, locus_size);

    int locus, marker;
    for (int i = 0; i < locus_size; i++) {
        locus = loci[i];
        marker = locus - pheno;

        if (debug) Rprintf("locus %d, pheno: %d, marker: %d\\n",
                           locus, pheno, marker);

        Rcpp::List alleleList(allele[locus]);
        Rcpp::CharacterVector aAlleleIndex(alleleList[3]);
        int allele1, allele2;

        for (int j = 0; j < genotype_sample_size; j++) {

            Rcpp::RawVector rv(genotype_sample[j]);
            if (rv.size() != 0) {
                allele1 = rv[2 * marker];
                allele2 = rv[2 * marker + 1];
            } else
                allele1 = allele2 = 0;  // indicate empty
/*
          if (debug && j <= 3)
              Rprintf("marker %d, t0 %x t1 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                       marker, allele1, allele2,
                       rv(marker-6), rv(marker-5), rv(marker-4), rv(marker-3), rv(marker-2), rv(marker-1),
                       rv(marker-0), rv(marker+1), rv(marker+2), rv(marker+3), rv(marker+4), rv(marker+5));
*/
          if (debug && j <= 3) Rprintf("marker %d, allele1 %x, allele2 %x \\n",  marker, allele1, allele2);

//        index1 = aAlleleIndex[ allele1 - 1 ];
//        index2 = aAlleleIndex[ allele2 - 1 ];

          mtx(j, i) = (allele1 << 16) | allele2;

        }
    }
    return mtx;
}

// [[Rcpp::export]]
Rcpp::Matrix<STRSXP> getgenotypes_Ri(NumericVector locus_arg,
                                     NumericVector hocus_arg,
                                     List          genotype_arg,
                                     List          allele_arg,
                                     List          markerscheme_arg,
                                     NumericVector phenocnt_arg)
{
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    int locus_size = loci.size();
    Rcpp::NumericVector hoci(hocus_arg);
    int hocus_size = hoci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // genotype[, 2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::List markerschemes(markerscheme_arg);
    Rcpp::IntegerVector allele1_map(markerschemes[2]);   // markerscheme_table[, 3]
    Rcpp::IntegerVector allele2_map(markerschemes[3]);   // markerscheme_table[, 4]

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, locus_size);

    if (locus_size != hocus_size) {
        Rf_error("First vector arguments should be the same length, but are %d vs %d\\n",
                 locus_size, hocus_size);
        return mtx;
    }

    int locus, hocus, marker, byte, a1map, a2map, t0;
    for (int j = 0; j < genotype_sample_size; j++) {

        Rcpp::RawVector rv(genotype_sample[j]);

        for (int i = 0; i < loci.size(); i++) {

            locus = loci[i];
            hocus = hoci[i];
            marker = hocus - pheno;

            byte = marker / 4;
            marker = marker - 4 * byte;

            if (debug) Rprintf("locus %d, hocus %d, pheno: %d, marker: %d, byte %d, offset %d\\n",
                               locus, hocus, pheno, hocus-pheno, byte, marker);

            a1map = allele1_map(locus - pheno);
            a2map = allele2_map(locus - pheno);

            Rcpp::CharacterVector aAlleleName(allele[1]);
            std::string allele1(aAlleleName[2 * locus + a1map - 1]);
            std::string allele2(aAlleleName[2 * locus + a2map - 1]);
            if (debug) Rprintf("allele%d/%d: %s%s; ", a1map, a2map,
                                  allele1.c_str(), allele2.c_str());

            Rcpp::IntegerVector aindexX(allele[3]);
            if (debug) Rprintf("indexX: %d %d\\n", aindexX[2 * locus], aindexX[2 * locus+1]);

            decode_allele[0] = allele1 + allele1;
            decode_allele[1] = "00";
            decode_allele[2] = allele1 + allele2;
            decode_allele[3] = allele2 + allele2;

            t0 = (rv.size() != 0) ? rv[byte]: 0x55;
/*
            if (debug && j <= 3)
                Rprintf("byte %d, t0 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                         byte, t0,
                         rv(byte-6), rv(byte-5), rv(byte-4), rv(byte-3), rv(byte-2), rv(byte-1),
                         rv(byte-0), rv(byte+1), rv(byte+2), rv(byte+3), rv(byte+4), rv(byte+5));
*/
            if (debug && j <= 3) Rprintf("byte %d, t0 %x \\n",  byte, t0);
            if (marker == 0) {
              mtx(j, i) = decode_allele[(t0 & 0x03) >> 0];
            } else if (marker == 1) {
              mtx(j, i) = decode_allele[(t0 & 0x0c) >> 2];
            } else if (marker == 2) {
              mtx(j, i) = decode_allele[(t0 & 0x30) >> 4];
            } else if (marker == 3) {
              mtx(j, i) = decode_allele[(t0 & 0xc0) >> 6];
            }
        }
    }
    return mtx;
}

// [[Rcpp::export]]
Rcpp::NumericVector getgenotypes_forperson(RawVector raw_arg)
{
    Rcpp::RawVector rv(raw_arg);

    Rcpp::NumericVector v(rv.size() * 4);

    if (rv.size() == 0)
        return v;

    int j = 0;
    for (int i = 0; i < rv.size(); i++) {
        int t0 = rv[i];
        int t1 = (t0 & 0x03) >> 0;
        int t2 = (t0 & 0x0c) >> 2;
        int t3 = (t0 & 0x30) >> 4;
        int t4 = (t0 & 0xc0) >> 6;
        v[j++] = t1;
        v[j++] = t2;
        v[j++] = t3;
        v[j++] = t4;
    }

    return v;
}

//////////////////////////////////////////////
// library(inline)                          //
// library(Rcpp)                            //
//                                          //
// getgenotypesa1 = inline::cxxfunction(    //
//     signature(locus_arg = "numeric",     //
//               hocus_arg = "numeric",     //
//               genotype_arg = "list",     //
//               allele_arg = "list",       //
//               markerscheme_arg = "list", //
//               phenocnt_arg = "numeric"), //
//     body =                               //
//     '                                    //
// //getgenotypesa1                         //
//     int debug = 0;                       //
//         ...                              //
//     return mtx;                          //
//     ',                                   //
//   verbose = TRUE,                        //
//   plugin = "Rcpp")                       //
//////////////////////////////////////////////
