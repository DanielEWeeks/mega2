library(DBI)
library(RSQLite)

library(Rcpp)
library(inline)

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

per = function(pid=1) {

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
