
#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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
#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db")
#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/rs2/dbmega2.db")

mk_markers_with_skip= function(mapselect=1) {
    markersPerChr = sapply(split(marker_table$chromosome, marker_table$chromosome), length)
    extra_markers = cumsum(4*floor((markersPerChr+3)/4) - markersPerChr)
    extra_markers = c(0, extra_markers)
    names(extra_markers)=NULL
    marker_table$locus_link_fill = marker_table$locus_link + extra_markers[marker_table$chromosome]
    assign("marker_table", marker_table, pos=globalenv());

    markers = merge(marker_table[ , c("locus_link","locus_link_fill","MarkerName","chromosome")],
                    map_table[ map_table$map == mapselect, c( "marker", "position")],
                    by.x="locus_link", by.y="marker")
    assign("markers", markers, pos=globalenv())
}

mk_unified_genotype_table = function(mapselect=1) {
    samples = split(genotype_table, genotype_table$person_link)
    samplesize = length(samples)
    person_link = unique(genotype_table$person_link)

    df = data.frame(row.names=1:samplesize,
                    person_link=person_link,
                    data=vector("raw", samplesize))

    for (i in 1:samplesize) {

        chrOrder = order(samples[[i]][,3], decreasing=FALSE)
        v = unlist(samples[[i]][chrOrder,5])
        df$data[i] = list(v)
    }
  
    assign("unified_genotype_table", df, pos=globalenv());

    mk_markers_with_skip(mapselect)
}

dbmega2_import = function(dbname="/Users/rbaron/mega2/test/mexnly/change_chrom/bcf/dbmega2.db",
                          mapselect=1) {
    con = dbConnect(RSQLite::SQLite(), dbname=dbname);

    for (tbl in TBLS) {
        if (dbExistsTable(con, tbl)) {
            cat(tbl, dbListFields(con, tbl), sep="\t", end="\n");
            assign(tbl, dbReadTable(con, tbl), pos=globalenv());
            print(dim(get(tbl, pos=globalenv())))
        }
    }
    mk_unified_genotype_table(mapselect)
}

geno_i = inline::cxxfunction(
    signature(a = "raw"), 

    '
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
    ',

  plugin = "Rcpp")

# gg=geno_i(g)

# a1=f$AlleleName.x
# a2=f$AlleleName.y

get_per = function(pid=1) {

  a1 = allele_table[allele_table$indexX==1,2][-1];
  a2 = allele_table[allele_table$indexX==2,2][-1];

  rv = genotype_table[pid, 5][[1]];
  rv4 = geno_i(rv);

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

################################################################

getlocus_R = function(locus=1, hocus=1, pheno=1) {

  return
    getlocus_Ri(locus, hocus, unified_genotype_table, allele_table, markerscheme_table, pheno)
 
}

getlocus_Ri = inline::cxxfunction(
    signature(locus_arg = "NumericVector", hocus_arg = "NumericVector", genotype_arg = "List", allele_arg = "List", markerscheme_arg = "List", phenocnt_arg = "NumericVector"),

    '
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    Rcpp::NumericVector hoci(locus_arg);
    int locus_size = loci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // genotype[,2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::List markerschemes(markerscheme_arg);
    Rcpp::IntegerVector allele1_map(markerschemes[2]);   // markerscheme_table[,3]
    Rcpp::IntegerVector allele2_map(markerschemes[3]);   // markerscheme_table[,4]

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, locus_size);

    if (locus_size != hoci.size()) {
        Rprintf("First vector arguments should be the same length, but are %d vs %d\\n",
                 locus_size, hoci.size());
        return mtx;
    }

    int locus, hocus, marker, byte, a1map, a2map;
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
            std::string allele1(aAlleleName[2*locus + a1map - 1]);
            std::string allele2(aAlleleName[2*locus + a2map - 1]);
            if (debug) Rprintf("allele%d/%d: %s%s; ", a1map, a2map, allele1.c_str(), allele2.c_str());

            Rcpp::IntegerVector aindexX(allele[3]);
            if (debug) Rprintf("indexX: %d %d\\n", aindexX[2*locus], aindexX[2*locus+1]);

            decode_allele[0] = allele1 + allele1;
            decode_allele[1] = "00";
            decode_allele[2] = allele1 + allele2;
            decode_allele[3] = allele2 + allele2;

            int t0 = rv[byte];
            if (0 && debug && j <= 3)
                Rprintf("byte %d, t0 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                    byte, t0,
                    rv(byte-6), rv(byte-5), rv(byte-4), rv(byte-3), rv(byte-2), rv(byte-1),
                    rv(byte-0), rv(byte+1), rv(byte+2), rv(byte+3), rv(byte+4), rv(byte+5));
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
    ',

  plugin = "Rcpp")

################################################################

getlocus = getlocus_C = function(locus=1, hocus=1, pheno=1) {

  return
    getlocus_Ci(locus, hocus, unified_genotype_table, allele_table, markerscheme_table, pheno)
 
}

getlocus_Ci = inline::cxxfunction(
    signature(locus_arg = "NumericVector", hocus_arg = "NumericVector",
              genotype_arg = "List", allele_arg = "List", markerscheme_arg = "List",
              phenocnt_arg = "NumericVector"),

    '
    int debug = 0;

    Rcpp::NumericVector loci(locus_arg);
    Rcpp::NumericVector hoci(hocus_arg);
    int locus_size = loci.size();

    Rcpp::List genotype(genotype_arg);
    Rcpp::List genotype_sample(genotype[1]);           // geno[,2]
    int genotype_sample_size = genotype_sample.size();

    Rcpp::List allele(allele_arg);
    std::vector<std::string> decode_allele(4);

    Rcpp::List markerschemes(markerscheme_arg);
    Rcpp::IntegerVector allele1_map(markerschemes[2]);   // markerscheme_table[,3]
    Rcpp::IntegerVector allele2_map(markerschemes[3]);   // markerscheme_table[,4]

    Rcpp::NumericVector phenos(phenocnt_arg);
    int pheno = phenos[0];

    Rcpp::Matrix<STRSXP> mtx(genotype_sample_size, locus_size);

    if (locus_size != hoci.size()) {
        Rprintf("First vector arguments should be the same length, but are %d vs %d\\n",
                 locus_size, hoci.size());
        return mtx;
    }

    int locus, hocus, marker, byte, a1map, a2map;
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

        Rcpp::CharacterVector aAlleleName(allele[1]);
        std::string allele1(aAlleleName[2*locus + a1map - 1]);
        std::string allele2(aAlleleName[2*locus + a2map - 1]);
        if (debug) Rprintf("allele%d/%d: %s%s; ", a1map, a2map, allele1.c_str(), allele2.c_str());

        Rcpp::IntegerVector aindexX(allele[3]);
        if (debug) Rprintf("indexX: %d %d\\n", aindexX[2*locus], aindexX[2*locus+1]);

        decode_allele[0] = allele1 + allele1;
        decode_allele[1] = "00";
        decode_allele[2] = allele1 + allele2;
        decode_allele[3] = allele2 + allele2;

        for (int j = 0; j < genotype_sample_size; j++) {

          Rcpp::RawVector rv(genotype_sample[j]);

          int t0 = rv[byte];
          if (0&& debug && j <= 3) Rprintf("byte %d, t0 %x %x %x %x %x %x %x %x %x %x %x %x %x\\n",
                  byte, t0,
                  rv(byte-6), rv(byte-5), rv(byte-4), rv(byte-3), rv(byte-2), rv(byte-1),
                  rv(byte-0), rv(byte+1), rv(byte+2), rv(byte+3), rv(byte+4), rv(byte+5));
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
    ',

  plugin = "Rcpp")

################################################################

allelediff = function(aa=aa, bb=bb, n=24) {
    for (i in 1:24) {
        print(sum(
                  ( (substr(aa[,i],1,1)==substr(bb[,i],1,1)) &
                    (substr(aa[,i],2,2)==substr(bb[,i],2,2)) )  |

                  ( (substr(aa[,i],1,1)==substr(bb[,i],2,2)) &
                     (substr(aa[,i],2,2)==substr(bb[,i],1,1)) )
                  )
              )
      }
}

################################################################

tst1 = function() {

    cat("getlocus_C", system.time((cc=getlocus_C(1:1000,1:1000,1))), "\n");
    cat("getlocus_R", system.time((dd=getlocus_R(1:1000,1:1000,1))), "\n");
                
    print(all(cc==dd));

    cat("getlocus_C", system.time((cc=getlocus_C(500000:501000,500000:501000,1))), "\n");
    cat("getlocus_R", system.time((dd=getlocus_R(500000:501000,500000:501000,1))), "\n");
    print(all(cc==dd));

# 25X slower
#   cat("getlocus_Cbind", system.time((ee=getlocus_Cbind(1:1000))), "\n");
#   print(all(cc==ee));
#   cat("getlocus_Cbind", system.time((ee=getlocus_Cbind(500000:501000))), "\n");
#   print(all(cc==ee));
    
    ## getlocus_C 0.845 0.008 0.854 0 0 
    ## getlocus_R 1.234 0.023 1.265 0 0 
    ## getlocus_Cbind 16.954 9.304 26.325 0 0 
    ## [1] TRUE
    ## [1] TRUE
    ## getlocus_C 0.965 0.018 0.984 0 0 
    ## getlocus_R 1.211 0.018 1.234 0 0 
    ## getlocus_Cbind 17.339 9.259 26.669 0 0 
    ## [1] TRUE
    ## [1] TRUE

    ## system.time((a=mkr(500500)))
    ##    user  system elapsed 
    ##  24.128   0.841  25.228 
}
################################################################
