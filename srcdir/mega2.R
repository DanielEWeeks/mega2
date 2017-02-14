
#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================



library(DBI)
library(RSQLite)

library(inline)
library(Rcpp)

concat = function(..., sep="") { return (paste(..., sep=sep)) }
lhead  = function(obj, ...)    { print(length(obj)); head(obj, ...) }



# c=(merge(allele_table[allele_table$indexX==1,], markerscheme_table[,c(1,2,3)], by.x="locus_link", by.y="key"))
# 
# d=(merge(allele_table[allele_table$indexX==2,], markerscheme_table[,c(1,2,4)], by.x="locus_link", by.y="key"))
# 
# e=(merge(c,d,by="locus_link"))
# f=(e[,c(1,3,4,7,9,10,13)])
# 
# g=(merge(locus_table, f, by="locus_link")
# h=g[,c(-2,-4,-5,-6,-7)]
# 
# head(merge(h, map_table[map_table$map==0,], by.x="locus_link", by.y="marker"))
# head(merge(h, map_table[map_table$map==1,], by.x="locus_link", by.y="marker"))
# 
#   locus_link  LocusName AlleleCnt AlleleName.x Frequency.x allele1 AlleleName.y
# 1          1 rs10458597         2            C  0.99967825       1            T
# 2          2 rs12565286         2            C  0.97610648       1            G
# 3          3 rs12082473         2            C  0.97674792       1            T
# 4          4  rs3094315         2            C  0.03109971       1            T
# 5          5  rs2286139         2            C  0.12110504       1            T
# 6          6 rs11240776         2            A  1.00000000       1        dummy
#    Frequency.y allele2 pId map position pos_female pos_male
# 1 0.0003217503       2   2   1   564621     -99.99   -99.99
# 2 0.0238935215       2   4   1   721290     -99.99   -99.99
# 3 0.0232520847       2   6   1   740857     -99.99   -99.99
# 4 0.9689002886       2   8   1   752566     -99.99   -99.99
# 5 0.8788949566       2  10   1   761732     -99.99   -99.99
# 6 0.0000000000       0  12   1   765269     -99.99   -99.99
# 
# ================
#    >  traitaff_table
#   pId ClassCnt PenCnt NumLabels Labels locus_link
# 1   1        1      3         0                 0
# 
# > affectclass_table
#   pId MaleDef FemaleDef AutoDef
# 1   1       1         1       0
#                                                          MalePen
# 1 9a, 99, 99, 99, 99, 99, a9, 3f, cd, cc, cc, cc, cc, cc, ec, 3f
#                                                                                        FemalePen
# 1 9a, 99, 99, 99, 99, 99, a9, 3f, cd, cc, cc, cc, cc, cc, ec, 3f, cd, cc, cc, cc, cc, cc, ec, 3f
#                                                                                          AutoPen
# 1 9a, 99, 99, 99, 99, 99, a9, 3f, cd, cc, cc, cc, cc, cc, ec, 3f, cd, cc, cc, cc, cc, cc, ec, 3f
#   locus_link class_link
# 1          0          0
# > affectclass_table[1,5]
# [[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f
# 
# > affectclass_table[1,6]
# [[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f cd cc cc cc cc cc ec 3f
# 
# > affectclass_table[1,5][[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f
# > readBin(affectclass_table[1,5][[1]], double(), 2, size=8)
# [1] 0.05 0.90
# > readBin(affectclass_table[1,6][[1]], double(), 3, size=8)
# [1] 0.05 0.90 0.90
# > readBin(affectclass_table[1,7][[1]], double(), 3, size=8)
# [1] 0.05 0.90 0.90
# 
# 
# ================================================================
#   j=(merge(pedigree_table, person_table, by="pedigree_link"))
# 
# k=(j[,c(1,37,3,4,13,6,16,18,19,23)])
# 
# head(k)
#   pedigree_link person_link Num EntryCnt      UniqueID   Name PerPre Father
# 1             0           0   1        1 SG0026_SG0026 SG0026 SG0026      0
# 2             1           1   2        1 SG0001_SG0001 SG0001 SG0001      0
# 3             2           2   3        1 SG0002_SG0002 SG0002 SG0002      0
# 4             3           3   4        1 SG0003_SG0003 SG0003 SG0003      0
# 5             4           4   5        1 SG0004_SG0004 SG0004 SG0004      0
# 6             5           5   6        1 SG0005_SG0005 SG0005 SG0005      0
#   Mother Sex
# 1      0   2
# 2      0   2
# 3      0   2
# 4      0   2
# 5      0   2
# 6      0   1
# 
#  phenotype_table
#       pId person_link cnt bytes                           data
# 1       1           0   1     8 00, 00, 00, 00, 01, 00, 00, 00
# 2       2           1   1     8 00, 00, 00, 00, 01, 00, 00, 00
# 
# x=function(y) {readBin(y, integer(), 2, size=4)}
# lapply(Phenotype_table$data, x)
#                        
  
TBLS = c("int_table",
         "double_table",
         "charstar_table",
         "stuff_table",
         "batch_parameters",
         "file_table",

         "pedigree_table",
         "person_table",
         "pedigree_brkloop_table",
         "person_brkloop_table",

         "canonicalallele_table",
         "locus_table",
         "allele_table",
         "marker_table",
         "markerscheme_table",
         "map_table",
         "mapnames_table",

         "traitaff_table",
         "affectclass_table",
         "traitquant_table",

         "phenotype_table",
         "genotype_table")

#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbgap2.db")

dbmega2_import = function(dbname="/Users/rbaron/mega2/test/mexnly/change_chrom/bcf/dbmega2.db") {
    con = dbConnect(RSQLite::SQLite(), dbname=dbname);

    for (tbl in TBLS) {
        if (dbExistsTable(con, tbl)) {
            assign(tbl, dbReadTable(con, tbl), pos=globalenv());
            cat(tbl, dbListFields(con, tbl), sep="\t", end="\n");
#            print(head(get(tbl, pos=globalenv())))
        }
    }
}

fnsig = signature(a = "raw")

fnsrc = '
  Rcpp::RawVector rv(a);

  Rcpp::NumericVector v(rv.size() * 4);

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
'

fn = inline::cxxfunction(fnsig, fnsrc, plugin = "Rcpp");

# gg=fn(g)

# a1=f$AlleleName.x
# a2=f$AlleleName.y

get_per = function(pid=1) {

  a1 = allele_table[allele_table$indexX==1,2][-1];
  a2 = allele_table[allele_table$indexX==2,2][-1];

  rv = genotype_table[pid, 5][[1]];
  rv4 = fn(rv);

#  0 1|1
#  1 0|0
#  2 1|2
#  3 2|2

  return
    ifelse(rv4==0, concat(a1, a1),
           ifelse(rv4==1, concat("00"),
                  ifelse(rv4==2, concat(a1, a2), concat(a2,a2))
                  )
           )
}

mkr = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  byte = marker - 1;
  byte = as.integer(byte / 4);
  marker = marker - 1 - byte * 4;
  rv  = matrix(unlist(genotype_table[, 5]),nrow=dim(person_table)[1],byrow=TRUE)[,byte+1]
  rv4 = fn(rv)
  rv4= matrix(rv4,ncol=4,byrow=TRUE)
  rv1 = rv4[,marker+1]

  return
    ifelse(rv1==0, concat(a1, a1),
           ifelse(rv1==1, concat("00"),
                  ifelse(rv1==2, concat(a1, a2), concat(a2,a2))
                  )
           )
}

# ################################################################


ftest = inline::cxxfunction(
  signature(a = "List"),

  '
  Rcpp::List ll(a);

  Rcpp::NumericVector v(3119);

  Rprintf("this is a test A %d\\n", ll.size());

  v[0] = ll.size();
  Rprintf("this is a test 1 %d\\n", ll.size());
  v = ll[1];
  Rprintf("this is a test 2 %d\\n", ll.size());

  Rcpp::List lx(ll[4]);
  Rprintf("this is a test 3\\n");


  Rprintf("this is a test B %d\\n", lx.size());
//return lx;

  Rcpp::RawVector rv(lx[1]);

  Rprintf("this is a test C %d\\n", rv.size());

  for (int i = 0; i < lx.size(); i++) {

    Rcpp::RawVector rv(lx[i]);

    Rprintf("%2d: ", i);
    for (int j = 0; j < 16; j++) {
         Rprintf("%02x ", rv[j]);
    }
    Rprintf("\\n");

  }

//  return rv;
  ',

  plugin = "Rcpp")

################################################################

p = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  return
    fp(genotype_table, marker, c(concat(a1, a1), "00", concat(a1, a2), concat(a2, a2)));
 
}

m = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  rv1 = fo(genotype_table, marker);

  return
    ifelse(rv1==0, concat(a1, a1),
           ifelse(rv1==1, concat("00"),
                  ifelse(rv1==2, concat(a1, a2), concat(a2,a2))
                  )
           )
  
}

n = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  rv1 = fo(genotype_table, marker);

  return
    c(concat(a1, a1), "00", concat(a1, a2), concat(a2, a2))[as.integer(rv1)+1]
 
}

fo = inline::cxxfunction(
  signature(genotype = "List", marker = "NumericVector"),

  '
  Rcpp::List ll(genotype);
  Rcpp::NumericVector mv(marker);
  int markerm1 = (mv[0]) - 1;

  int byte = markerm1 / 4;
  markerm1 = markerm1 - 4 * byte;

  Rcpp::List lx(ll[4]);           // genotype[,5]

  Rcpp::RawVector out(lx.size());
  int j = 0;

  for (int i = 0; i < lx.size(); i++) {

    Rcpp::RawVector rv(lx[i]);

    int t0 = rv[byte];
    if (markerm1 == 0) {
      out[j++] = (t0 & 0x03) >> 0;
    } else if (markerm1 == 1) {
      out[j++] = (t0 & 0x0c) >> 2;
    } else if (markerm1 == 2) {
      out[j++] = (t0 & 0x30) >> 4;
    } else if (markerm1 == 3) {
      out[j++] = (t0 & 0xc0) >> 6;
    }

  }

  return out;
  ',

  plugin = "Rcpp")

################################################################

p = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  return
    fp(genotype_table, marker, c(concat(a1, a1), "00", concat(a1, a2), concat(a2, a2)));
 
}

fp = inline::cxxfunction(
  signature(genotype = "List", marker = "NumericVector", geno = "CharacterVector"),

  '
  Rcpp::List ll(genotype);
  Rcpp::NumericVector mv(marker);
  int markerm1 = (mv[0]) - 1;
  Rcpp::CharacterVector g(geno);

  int byte = markerm1 / 4;
  markerm1 = markerm1 - 4 * byte;

  Rcpp::List lx(ll[4]);           // genotype[,5]

  Rcpp::CharacterVector out(lx.size());
  int j = 0;

  for (int i = 0; i < lx.size(); i++) {

    Rcpp::RawVector rv(lx[i]);

    int t0 = rv[byte];
    if (markerm1 == 0) {
      out[j++] = g[(t0 & 0x03) >> 0];
    } else if (markerm1 == 1) {
      out[j++] = g[(t0 & 0x0c) >> 2];
    } else if (markerm1 == 2) {
      out[j++] = g[(t0 & 0x30) >> 4];
    } else if (markerm1 == 3) {
      out[j++] = g[(t0 & 0xc0) >> 6];
    }

  }

  return out;
  ',

  plugin = "Rcpp")

################################################################

q = function(marker=1,pheno=1) {

  a1 = allele_table[allele_table$indexX==1,2][marker+pheno];
  a2 = allele_table[allele_table$indexX==2,2][marker+pheno];

  return
    fq(marker, genotype_table, allele_table, pheno)
 
}

fq = inline::cxxfunction(
  signature(marker = "NumericVector", genotype = "List", allele = "List", phen = "NumericVector"),

  '
  Rcpp::NumericVector mv(marker);
  int markerm1 = (mv[0]) - 1;

  Rcpp::List ll(genotype);

  Rcpp::List al(allele);
  std::vector<std::string> g(4);

  int byte = markerm1 / 4;
  markerm1 = markerm1 - 4 * byte;

  Rcpp::List lx(ll[4]);           // genotype[,5]

  Rcpp::NumericVector ph(phen);
  int pheno = ph[0];
  int off = mv[0] - 1 + pheno;
  Rprintf("pheno: %d, off: %d\\n", pheno, off);

  Rcpp::CharacterVector aAlleleName(al[1]);
  std::string al1(aAlleleName[2*off]);
  std::string al2(aAlleleName[2*off+1]);
  Rcpp::IntegerVector aindexX(al[3]);
  Rprintf("al: %d %d %d %d\\n", aindexX[0], aindexX[1], aindexX[2], aindexX[3]);
  Rprintf("al2: %s%s\\n", al1.c_str(), al2.c_str());

  g[0] = al1 + al1;
  g[1] = "00";
  g[2] = al1 + al2;
  g[3] = al2 + al2;

  Rcpp::CharacterVector out(lx.size());
  int j = 0;

  for (int i = 0; i < lx.size(); i++) {

    Rcpp::RawVector rv(lx[i]);

    int t0 = rv[byte];
    if (markerm1 == 0) {
      out[j++] = g[(t0 & 0x03) >> 0];
    } else if (markerm1 == 1) {
      out[j++] = g[(t0 & 0x0c) >> 2];
    } else if (markerm1 == 2) {
      out[j++] = g[(t0 & 0x30) >> 4];
    } else if (markerm1 == 3) {
      out[j++] = g[(t0 & 0xc0) >> 6];
    }

  }

  return out;
  ',

  plugin = "Rcpp")

################################################################

getmarker_Cbind = function(marker=1, pheno=1) {

  return
    getmarker_Cbindi(marker, genotype_table, allele_table, pheno)
 
}

getmarker_Cbindi = inline::cxxfunction(
    signature(marker_arg = "NumericVector", genotype_arg = "List", allele_arg = "List", pheno_arg = "NumericVector"),

    '
    int debug = 0;

    Rcpp::NumericVector markers(marker_arg);
    int marker_size = markers.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[4]);           // genotype[,5]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::NumericVector phenos(pheno_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, 0);

    int marker;
    for (int i = 0; i < marker_size; i++) {
        marker = markers[i];
        int markerm1 = marker - 1;

        int byte = markerm1 / 4;
        markerm1 = markerm1 - 4 * byte;

        int locus = marker - 1 + pheno;
        if (debug) Rprintf("marker-1 %d, pheno: %d, locus: %d, byte %d, markerm1 %d\\n",
                           marker-1, pheno, locus, byte, markerm1);

        Rcpp::CharacterVector aAlleleName(allele[1]);
        std::string allele1(aAlleleName[2*locus]);
        std::string allele2(aAlleleName[2*locus+1]);
        if (debug) Rprintf("allele1/2: %s%s; ", allele1.c_str(), allele2.c_str());

        Rcpp::IntegerVector aindexX(allele[3]);
        if (debug) Rprintf("indexX: %d %d\\n", aindexX[2*locus], aindexX[2*locus+1]);

        decode_allele[0] = allele1 + allele1;
        decode_allele[1] = "00";
        decode_allele[2] = allele1 + allele2;
        decode_allele[3] = allele2 + allele2;

        Rcpp::CharacterVector oneCol(genotype_sample_size);
        int k = 0;

        for (int j = 0; j < genotype_sample_size; j++) {

          Rcpp::RawVector rv(genotype_sample[j]);

          int t0 = rv[byte];
          if (markerm1 == 0) {
            oneCol[k++] = decode_allele[(t0 & 0x03) >> 0];
          } else if (markerm1 == 1) {
            oneCol[k++] = decode_allele[(t0 & 0x0c) >> 2];
          } else if (markerm1 == 2) {
            oneCol[k++] = decode_allele[(t0 & 0x30) >> 4];
          } else if (markerm1 == 3) {
            oneCol[k++] = decode_allele[(t0 & 0xc0) >> 6];
          }
        }
        mtx = Rcpp::cbind(mtx, oneCol);
    }
    return mtx;
    ',

  plugin = "Rcpp")

################################################################

getmarker_R = function(marker=1, pheno=1) {

  return
    getmarker_Ri(marker, genotype_table, allele_table, pheno)
 
}

getmarker_Ri = inline::cxxfunction(
    signature(marker_arg = "NumericVector", genotype_arg = "List", allele_arg = "List", pheno_arg = "NumericVector"),

    '
    int debug = 0;

    Rcpp::NumericVector markers(marker_arg);
    int marker_size = markers.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[4]);           // genotype[,5]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::NumericVector phenos(pheno_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, marker_size);

    int marker;
    for (int j = 0; j < genotype_sample_size; j++) {

        Rcpp::RawVector rv(genotype_sample[j]);

//      Rcpp::CharacterVector oneRow(marker_size);

        for (int i = 0; i < markers.size(); i++) {

            marker = markers[i];
            int markerm1 = marker - 1;

            int byte = markerm1 / 4;
            markerm1 = markerm1 - 4 * byte;

            int locus = marker - 1 + pheno;
            if (debug) Rprintf("marker-1 %d, pheno: %d, locus: %d, byte %d, markerm1 %d\\n",
                                marker-1, pheno, locus, byte, markerm1);

            Rcpp::CharacterVector aAlleleName(allele[1]);
            std::string allele1(aAlleleName[2*locus]);
            std::string allele2(aAlleleName[2*locus+1]);
            if (debug) Rprintf("allele1/2: %s%s; ", allele1.c_str(), allele2.c_str());

            Rcpp::IntegerVector aindexX(allele[3]);
            if (debug) Rprintf("indexX: %d %d\\n", aindexX[2*locus], aindexX[2*locus+1]);

            decode_allele[0] = allele1 + allele1;
            decode_allele[1] = "00";
            decode_allele[2] = allele1 + allele2;
            decode_allele[3] = allele2 + allele2;

            int t0 = rv[byte];
            if (markerm1 == 0) {
              mtx(j, i) = decode_allele[(t0 & 0x03) >> 0];
            } else if (markerm1 == 1) {
              mtx(j, i) = decode_allele[(t0 & 0x0c) >> 2];
            } else if (markerm1 == 2) {
              mtx(j, i) = decode_allele[(t0 & 0x30) >> 4];
            } else if (markerm1 == 3) {
              mtx(j, i) = decode_allele[(t0 & 0xc0) >> 6];
            }
        }
    }
    return mtx;
    ',

  plugin = "Rcpp")

################################################################
getmarker = getmarker_C = function(marker=1, pheno=1) {

  return
    getmarker_Ci(marker, genotype_table, allele_table, pheno)
 
}

getmarker_Ci = inline::cxxfunction(
    signature(marker_arg = "NumericVector", genotype_arg = "List", allele_arg = "List", pheno_arg = "NumericVector"),

    '
    int debug = 0;

    Rcpp::NumericVector markers(marker_arg);
    int marker_size = markers.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[4]);           // genotype[,5]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::NumericVector phenos(pheno_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, marker_size);

    int marker;
    for (int i = 0; i < marker_size; i++) {
        marker = markers[i];
        int markerm1 = marker - 1;

        int byte = markerm1 / 4;
        markerm1 = markerm1 - 4 * byte;

        int locus = marker - 1 + pheno;
        if (debug) Rprintf("marker-1 %d, pheno: %d, locus: %d, byte %d, markerm1 %d\\n",
                           marker-1, pheno, locus, byte, markerm1);

        Rcpp::CharacterVector aAlleleName(allele[1]);
        std::string allele1(aAlleleName[2*locus]);
        std::string allele2(aAlleleName[2*locus+1]);
        if (debug) Rprintf("allele1/2: %s%s; ", allele1.c_str(), allele2.c_str());

        Rcpp::IntegerVector aindexX(allele[3]);
        if (debug) Rprintf("indexX: %d %d\\n", aindexX[2*locus], aindexX[2*locus+1]);

        decode_allele[0] = allele1 + allele1;
        decode_allele[1] = "00";
        decode_allele[2] = allele1 + allele2;
        decode_allele[3] = allele2 + allele2;

        Rcpp::CharacterVector oneCol(genotype_sample_size);

        for (int j = 0; j < genotype_sample_size; j++) {

          Rcpp::RawVector rv(genotype_sample[j]);

          int t0 = rv[byte];
          if (markerm1 == 0) {
            mtx(j, i) = decode_allele[(t0 & 0x03) >> 0];
          } else if (markerm1 == 1) {
            mtx(j, i) = decode_allele[(t0 & 0x0c) >> 2];
          } else if (markerm1 == 2) {
            mtx(j, i) = decode_allele[(t0 & 0x30) >> 4];
          } else if (markerm1 == 3) {
            mtx(j, i) = decode_allele[(t0 & 0xc0) >> 6];
          }
        }
    }
    return mtx;
    ',

  plugin = "Rcpp")
