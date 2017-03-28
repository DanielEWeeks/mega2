
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

## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/TxDb.Hsapiens.UCSC.hg19.knownGene/extdata/TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")


## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("org.Hs.eg.db")

#' mkfam
#'
#' @param brkloop ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
mkfam = function (brkloop = F) {

    if (brkloop) {
        ped = ENV$pedigree_brkloop_table
        per = ENV$person_brkloop_table
    } else {
        ped = ENV$pedigree_table
        per = ENV$person_table
    }
    dofam = function(per) {
        per$Father=per[match(per$Father, per$OrigID), "PerPre"]
        per[is.na(per$Father), "Father"] = 0
        per$Mother=per[match(per$Mother, per$OrigID), "PerPre"]
        per[is.na(per$Mother), "Mother"] = 0
        per
    }
    per=unsplit(lapply(split(per, per$pedigree_link), dofam), per$pedigree_link)
    perplus = merge(ped[ , c("pedigree_link", "PedPre")],
                    per[ , c("pedigree_link", "person_link", "PerPre", "Father", "Mother", "Sex")],
                    by = c("pedigree_link"))

    trait = ENV$phenotype_table[ , c("person_link", "data")]
    trait$trait = sapply(trait$data, function (x) { readBin(x, integer(), 2, size = 4) })[1, ]
    ENV$fam = merge(perplus, trait[ , c("person_link", "trait")], by = "person_link")
}

###############

#' applyFnToGenes
#'
#' @param op ...
#' @param genes ...
#' @param type ...
#' @param fuzz ...
#' @param ranges ...
#' @param chrs ...
#' @param marks ...
#'
#' @return NO
## @import TxDb.Hsapiens.UCSC.hg19.knownGene
## @import org.Hs.eg.db
## @importClassesFrom GenomicFeatures TxDb
## @importMethodsFrom AnnotationDbi columns keys keytypes seqlevels select
## @importMethodsFrom AnnotationDbi select
#' @export
#'
#' @examples
#'\dontrun{
#'}
applyFnToGenes = function (op = function (geno, mrkrs, rng) {},
                           genes = c("ELL2", "CARD15"),
                           type = "TX",
                           fuzz = 0,
                           ranges = matrix(ncol = 3, nrow = 0),
                           chrs = vector("integer", 0),
                           marks = vector("character", 0)) {

    ## dbconn(gene)/dbConn(txdb)
    ## dbReadTable(dbconn(), "tbl")

    loadNamespace(TxDb.Hsapiens.UCSC.hg19.knownGene)
    loadNamespace(org.Hs.eg.db)

    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    if (type=="TX")
        COLS = c("TXNAME", "TXID", "TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
    else
        COLS = c("EXONNAME", "EXONID", "EXONSTRAND", "EXONCHROM", "EXONSTART", "EXONEND")

    seqlevels(txdb) = paste("chr", c(1:22, "X", "Y", "M"), sep="")
    
    genedb = org.Hs.eg.db

    pa = select(genedb, keys = genes, columns = c("ALIAS", "ENTREZID", "SYMBOL"), keytype = "ALIAS")
    pb = select(txdb, keys = pa[ , 2], columns = COLS, keytype = "GENEID")
    range = merge(pa, pb, by.x = "ENTREZID", by.y = "GENEID")
    range[ ,6] = as.integer( sub("chr", "", range[ ,6]))

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
                                         ranges[i,1], "-", ranges[i,2], ranges[i,3]) )
                      }

    applyFnToRanges(op, range, c(6, 8, 9))

    #marks
    if (length(marks)) {
        applyFnToMarkers(op, ENV$markers[ENV$markers$MarkerName %in% marks, ])
    }
}

################

#' applyFnToRanges
#'
#' @param op ...
#' @param range ...
#' @param indices ...
#' @param fuzz ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
applyFnToRanges = function (op = function (geno, mrkrs, rng) {},
                            range = NULL,
                            indices = NULL,
                            fuzz = 0 ) {

    if (is.null(range)) {
        range   = ENV$refGene
        indices = ENV$refindices
    }
    rows = nrow(range)
    if (rows) {

        start = range[ , indices[2]]
        end   = range[ , indices[3]]
        chrm = as.integer(sub("chr", "", range[ , indices[1]]))

        for (i in 1:rows) {

            if (is.na(chrm[i]) || is.na(start[i]) || is.na(end[i]) ) next

            markersub = ENV$markers[ ENV$markers$chromosome == chrm[i] & ENV$markers$position <= (end[i] + fuzz) & ENV$markers$position >= (start[i] - fuzz), ]

            if (nrow(markersub)) {
                geno = getgenotypes(markersub)
                op(geno, markersub, range[i,])
            } else {
                message("No markers in range:  chr", chrm[i], " between ", start[i], " and ", end[i], "\n")
            }
        }
    }
}

#' applyFnToMarkers
#'
#' @param op ...
#' @param markers ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
applyFnToMarkers = function (op = function (geno, mrkrs, rng) {},
                             markers) {

    geno = getgenotypes(markers)
    op(geno, markers, NULL)

}
