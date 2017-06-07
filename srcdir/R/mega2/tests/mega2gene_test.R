
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

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

#dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db", verbose = 1)

pos  = NULL
geno = NULL


#' tst10
#'
#' @param genes ...
#' @param type ...
#' @param fuzz ...
#'
#' @return NO
#' @importFrom utils head
#' @export
#'
#' @examples
#'\dontrun{
#'}
tst10 = function(ENV, genes = c("ELL2", "CARD15"), type = "TX", fuzz = 0) {

    print("## test latest applyFnToGenes")
    
    mkfam(envir = ENV)
    show = function(g, m, r, e) {
        print(r)
        print(m)
        print(head(g))
    }
    applyFnToGenes(show, genes = genes, type = type, fuzz = fuzz,
              matrix(c(11, 50000000, 50100000,
                       11, 60000000, 60100000), ncol = 3, nrow = 2, byrow = T),
              markers = ENV$markers[! duplicated(ENV$markers$chromosome), 3],
              chrs=c(24, 26),
              envir = ENV)
}

#' mkmarkers
#'
#' @param genes ...
#' @param type ...
#' @param fuzz ...
#' @param ranges ...
#' @param chrs ...
#' @param marks ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
mkmarkers = function (genes = c("ELL2", "CARD15"),
                      type = "TX",
                      fuzz = 0,
                      ranges = matrix(ncol = 3, nrow = 0),
                      chrs = vector("integer", 0),
                      mrkrs = vector("character", 0),
                      envir = ENV) {

    ## older version of applyFnToGenes with globals

    ## dbconn(gene)/dbConn(txdb)
    ## dbReadTable(dbconn(), "tbl")

    marker_table = envir$marker_table
    map_table    = envir$map_table
    markers = merge(
                   marker_table[ , c("locus_link", "locus_link_fill", "MarkerName", "chromosome")],
                   map_table[ map_table$map == 1, c( "marker", "position")],
                   by.x = "locus_link", by.y = "marker")
#    envir$markers = markers
    assign("markers", markers, pos= globalenv())
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    if (type=="TX")
        COLS = c("TXNAME", "TXID", "TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
    else
        COLS = c("EXONNAME", "EXONID", "EXONSTRAND", "EXONCHROM", "EXONSTART", "EXONEND")

    seqlevels(txdb) = paste("chr", c(1:22, "X", "Y", "M"), sep="")
    
    genedb=org.Hs.eg.db

    pa = select(genedb, keys = genes, columns = c("ALIAS", "ENTREZID", "SYMBOL"), keytype = "ALIAS")
    pb = select(txdb, keys = pa[,2], columns = COLS, keytype = "GENEID")
    range = merge(pa, pb, by.x = "ENTREZID", by.y = "GENEID")
    range[ , 6] = as.integer( sub("chr", "", range[ , 6]))

    #chrs
    for (chr in chrs) { range = rbind(range,
#           ENTREZID  ALIAS SYMBOL  TXID     TXNAME TXCHROM TXSTRAND  TXSTART    TXEND
                                    list("-", "-", "-", "-", "-",
                                         chr, "-", 0, 1000000000) )
                      }
    #ranges
    if (length(ranges))
    for (i in 1:dim(ranges)[1]) { range = rbind(range,
#           ENTREZID  ALIAS SYMBOL  TXID     TXNAME TXCHROM TXSTRAND  TXSTART    TXEND
                                    list("-", "-", "-", "-", "-",
                                         ranges[i, 1], "-", ranges[i, 2], ranges[i, 3]) )
                      }

    #marks
    rows = dim(range)[1]
    if (length(mrkrs)) {
        range = rbind(range,
                        list("-", "-", "-", "-", "-",
                             0, "-", 0, 0) )
        pos = vector("list", rows+1)
        geno = vector("list", rows+1)
    } else {
        pos = vector("list", rows)
        geno = vector("list", rows)
    }
    assign("range", range, pos= globalenv())
    for (i in 1:rows) {
        print(i)
        if (is.na(range[i, 6]) || is.na(range[i, 8]) || is.na(range[i, 9]) ) next
        pos[[i]] = markers[ markers$chromosome == range[i, 6] & markers$position <= (range[i, 9] + fuzz) & markers$position >= (range[i, 8] - fuzz), ]
#       pos[[i]]$locus_link_fill = pos[[i]]$locus_link + chr_gap_skip[pos[[i]]$chromosome]
#       rownames(pos[[i]]) = NULL
        
        int_table = envir$int_table
        geno[[i]] = getgenotypes(pos[[i]], envir)
    }

    pos[[rows+1]] = markers[markers$MarkerName %in% mrkrs, ]
#   pos[[rows+1]]$locus_link_fill = pos[[rows+1]]$locus_link + chr_gap_skip[pos[[rows+1]]$chromosome]
#   rownames(pos[[rows+1]]) = NULL
    assign("pos", pos, pos = globalenv())

    geno[[rows+1]] = getgenotypes(pos[[rows+1]], envir)
    assign("geno", geno, pos = globalenv())
}


#' tst11
#'
#' @param genes ...
#' @param type ...
#' @param fuzz ...
#'
#' @return NO
#' @importFrom utils head 
#' @export
#'
#' @examples
#'\dontrun{
#'}
tst11 = function(ENV, genes = c("ELL2", "CARD15"), type = "TX", fuzz = 0) {
    print("## test earlier version applyFnToGenes; same results")
    mkfam(envir = ENV)
    mkmarkers(genes = genes, type = type, fuzz = fuzz,
              matrix(c(11, 50000000, 50100000,
                       11, 60000000, 60100000), ncol = 3, nrow = 2, byrow = T),
              mrkrs = ENV$markers[! duplicated(ENV$markers$chromosome), 3],
              chrs=c(24,26),
              envir = ENV)

    for (i in 1:dim(range)[1]) {
        print(range[i,])
        print(pos[[i]])
        print(head(geno[[i]]))
    }
}

#' tst12
#'
#' @param genes ...
#' @param type ...
#' @param fuzz ...
#'
#' @return NO
#' @importFrom utils head 
#' @export
#'
#' @examples
#'\dontrun{
#'}
tst12 = function(ENV, genes = c("ELL2"), type = "TX", fuzz = 0) {
    print("## test earlier version applyFnToGenes specific markers")
    mkfam(envir = ENV)
    mkmarkers(genes = genes, type = type, fuzz = fuzz,
              mrkrs = c("rs6587762", "rs7521920",
                        "rs10181821", "rs10195681", "rs7594567", "rs4637157"),
              envir = ENV)

    for (i in 1:dim(range)[1]) {
        print(range[i, ])
        print(pos[[i]])
        print(head(geno[[i]]))
    }
}
