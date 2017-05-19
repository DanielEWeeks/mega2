
#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

library(mega2)

#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbgap2.db")
#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db")
#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/rs2/dbmega2.db")

# c=(merge(allele_table[allele_table$indexX==1, ], markerscheme_table[ , 1:3], by.x="locus_link", by.y="key"))
# 
# d=(merge(allele_table[allele_table$indexX==2, ], markerscheme_table[ ,c(1, 2, 4)], by.x="locus_link", by.y="key"))
# 
# e=(merge(c, d, by="locus_link"))
# f=(e[ , c(1, 3, 4, 7, 9, 10, 13)])
# 
# g=(merge(locus_table, f, by="locus_link")
# h=g[ , c(-2, -4, -5, -6, -7)]
# 
# head(merge(h, map_table[map_table$map==0, ], by.x="locus_link", by.y="marker"))
# head(merge(h, map_table[map_table$map==1, ], by.x="locus_link", by.y="marker"))
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
#  >  traitaff_table
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
# > affectclass_table[1, 5]
# [[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f
# 
# > affectclass_table[1,6]
# [[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f cd cc cc cc cc cc ec 3f
# 
# > affectclass_table[1, 5][[1]]
#  [1] 9a 99 99 99 99 99 a9 3f cd cc cc cc cc cc ec 3f
# > readBin(affectclass_table[1, 5][[1]], double(), 2, size=8)
# [1] 0.05 0.90
# > readBin(affectclass_table[1, 6][[1]], double(), 3, size=8)
# [1] 0.05 0.90 0.90
# > readBin(affectclass_table[1, 7][[1]], double(), 3, size=8)
# [1] 0.05 0.90 0.90
# 
# 
# ================================================================
#   j=(merge(pedigree_table, person_table, by="pedigree_link"))
# 
# k=(j[ , c(1, 37, 3, 4, 13, 6, 16, 18, 19, 23)])
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

################################################################


#' allelediff
#'
#' @param aa ...
#' @param bb ...
#' @param n ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
allelediff = function(aa=aa, bb=bb, n=24) {
    for (i in 1:24) {
        print(sum(
                  ( (substr(aa[, i], 1, 1) == substr(bb[ , i], 1, 1)) &
                    (substr(aa[, i], 2, 2) == substr(bb[ , i], 2, 2)) )  |

                  ( (substr(aa[, i], 1, 1) == substr(bb[ , i], 2, 2)) &
                    (substr(aa[, i], 2, 2) == substr(bb[ , i], 1, 1)) )
                  )
              )
      }
}

################################################################

#' tst1
#'
#' @param db ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
tst1 = function(db = "~/mega2/data/samoaqwas.db") {

    print("## getgenotypes by column and by row; compare two ranges")

    ENV = dbmega2_import(db, verbose = 1)
    
    cat("getgenotypes",   system.time((cc=getgenotypes  (ENV$markers[1:1000,]))), "\n")
    cat("getgenotypes_R", system.time((dd=getgenotypes_R(ENV$markers[1:1000,]))), "\n")
                
    print(all(cc==dd))

    cat("getgenotypes",   system.time((cc=getgenotypes  (ENV$markers[500000:501000,]))), "\n")
    cat("getgenotypes_R", system.time((dd=getgenotypes_R(ENV$markers[500000:501000,]))), "\n")
    print(all(cc==dd))

# 25X slower
#   cat("getgenotypes_Cbind", system.time((ee=getgenotypes_Cbind(1:1000))), "\n")
#   print(all(cc==ee))
#   cat("getgenotypes_Cbind", system.time((ee=getgenotypes_Cbind(500000:501000))), "\n")
#   print(all(cc==ee))
    
    ## getgenotypes_C 0.845 0.008 0.854 0 0 
    ## getgenotypes_R 1.234 0.023 1.265 0 0 
    ## [1] TRUE
    ## getgenotypes_C 0.965 0.018 0.984 0 0 
    ## getgenotypes_R 1.211 0.018 1.234 0 0 
    ## [1] TRUE

}
################################################################
