
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

#' Mega2VCF package
#'
#' @description This package reads a Mega2 SQLite3 database into R dataframes and
#'	generates a VCF file with the same contents
#'
#' @author Robert V Baron
#' @docType package
#' @name Mega2VCF-package
NULL

#' generate a VCF file
#'
#' @description
#'  Generate a VCF file from the specified Mega2 SQLite database.
#'
#' @param db specify SQLite database to load
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @return None
#'
#' @importFrom mega2 dbmega2_import getgenotypesraw mkfam setfam
#' @importFrom utils write.table
#' @export
#'
#' @examples
#'\dontrun{
#' Mega2VCF()
#'}
Mega2VCF = function(db="ped1.db", mapno = 0) {

    ENV$LocusCnt = ENV$int_table[ENV$int_table$key == 'LocusCnt', 3]

    ENV = dbmega2_import("ped1.db")
    setfam(mkfam())

    unlink("foo")

    mkVCFhdr(ENV)

    j = 0
#   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    n1 = ENV$allele_table[ENV$allele_table$indexX==1,][-(ENV$PhenoCnt),]
    n2 = ENV$allele_table[ENV$allele_table$indexX==2,][-(ENV$PhenoCnt),]
## Mega2 "mis-feature"
    n1$AlleleName[n1$AlleleName %in% c("dummy")] = '.'
    n2$AlleleName[n2$AlleleName %in% c("dummy")] = '.'

    M = ENV$LocusCnt - ENV$PhenoCnt;
    while (TRUE) {
#      if (M < 1000) browser()
#      if (M != ENV$LocusCnt - ENV$PhenoCnt) break
#       print(M)
        if (M <= 0) break
        N = ((j*1000+1):(j*1000+1000))
        L = length(N)
        if (M < L) {
            L = M
            N = ((j*1000+1):(j*1000+L))
        }
        M = M - L

        chrm = ENV$markers[N, c("chromosome", "position", "MarkerName")]
        names(chrm) = c("#CHROM", "POS", "ID")
##new
        REF = n1$AlleleName[N]
        RF  = n1$Frequency[N]
        ALT = n2$AlleleName[N]
        AF  = n2$Frequency[N]
        XXF = REF
        XF  = RF
        
        flip = RF < AF
        doFlip = sum(flip)
        if (doFlip) {
            whichFlip = which(flip)
            REF[whichFlip] = ALT[whichFlip]
            ALT[whichFlip] = XXF[whichFlip]

            RF[whichFlip] = AF[whichFlip]
            AF[whichFlip] = XF[whichFlip]
        }
##newer
        GPos = ENV$map_table[ENV$map_table$map==mapno, c("position", "pos_female", "pos_male")][N, ]
        GPosPos = sprintf("%.2f", GPos$position)
        GPosFem = ifelse (GPos$pos_female != -99.99, sprintf("%f", GPos$pos_female), ".")
        GPosMal = ifelse (GPos$pos_male != -99.99,   sprintf("%f", GPos$pos_male),   ".")
##newer
##new
        chrm = cbind(chrm,
                     REF,
                     ALT,
                     QUAL=rep(".",      times=L),
                     FILTER=rep("PASS", times=L),
                     INFO=paste0("CM=", GPosPos, ",", GPosFem, ",", GPosMal,
                                 ";RF=", sprintf("%f", RF),
                                 ";AF=", sprintf("%f", AF), ";"),
                     FORMAT=rep("GT",   times=L)
               )

        cr = getgenotypesraw(ENV$markers[N, ])
        a1 = t(cr)
        a2 = a1
        dm = dim(a1)
        a1 = bitwShiftR(a1, 16)
        attr(a1, "dim") = dm
        x3 = a1
        a2 = bitwAnd(a2, 65535)
        attr(a2, "dim") = dm
        if (doFlip) {
            a1[, whichFlip] = a2[, whichFlip]
            a2[, whichFlip] = x3[, whichFlip]
        }

        a3 = ifelse(a1 != 0, as.character(a1-1), ".")
        a4 = ifelse(a2 != 0, as.character(a2-1), ".")
        a5 = paste0(a3, "/", a4)
        attr(a5, "dim") = dm
        a6 = cbind(chrm, a5)
        
        if (j == 0) {
            cat(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
                paste0(ENV$fam$PedPre, "_", ENV$fam$PerPre), sep="\t", append=TRUE, file="foo")
            cat("\n", append=TRUE, file="foo")
        }

## Mega2 "mis-feature"
        a6$POS = paste0(a6$POS, " ")
## Mega2 "mis-feature"
        ncols = ncol(a6)
        a6[, ncols] = paste0(a6[, ncols], "\t")
##
        write.table(a6, file="foo", quote=FALSE, append=TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
        j = j + 1
        message(j)
    }
}
vcf = Mega2VCF

#' generate required VCF header
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.  Finally, generate
#'   ##contig entries for each chromosome.
#'
#' @param ENV "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkVCFhdr(ENV)
#'}
mkVCFhdr = function (ENV) {
    cat('##fileformat=VCFv4.1\n', file="foo", append=TRUE)
    cat('##filedate=20170407\n', file="foo", append=TRUE)
    cat('##source=MEGA2\n', file="foo", append=TRUE)
    cat('##INFO=<ID=CM,Number=3,Type=Float,Description="Genetic Distance in centimorgans (avg, male, female)">\n', file="foo", append=TRUE)
    cat('##INFO=<ID=RF,Number=1,Type=Float,Description="Allele Frequency of reference allele">\n', file="foo", append=TRUE)
    cat('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency of alternate allele(s)">\n', file="foo", append=TRUE)
    cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file="foo", append=TRUE)
    cat('##FILTER=<ID=PASS,Description="Passed variant FILTERs">\n', file="foo", append=TRUE)

    j = 0
    for (i in sapply(split(ENV$markers, ENV$markers$chromosome), function(x) max(x$position)+1)) {
        j = j + 1
        cat('##contig=<ID=', j, ',length=', i, ',assembly=B37>\n', file="foo", append=TRUE, sep="")
    }
}
