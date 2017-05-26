
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

#' Mega2pedgene package
#'
#' @description
#'  This package runs the \bold{pedgene} program via the \bold{Mega2} BioConductor Framework.
#'  The framework result is stored in a file and can then be matched to a collection of data
#'  previously analyzed by manual methods.  The results can be compared.
#'
#' @author
#'  Robert V Baron, rvb5@pitt.edu
#'  Maintainer: Robert V Baron, <rvb5@pitt.edu>
#' @docType package
#' @name Mega2pedgene-package
#'
#'@references
#'  This optional section can contain literature or other references for
#'  background information.
#'
#'@seealso
#'  Optional links to other man pages
#'
#'@examples
#'  \dontrun{
#'     ## Optional simple examples of the most important functions
#'     ## These can be in \dontrun{} and \donttest{} blocks.
#'  }
NULL

#library(mega2)
#library(pedgene)

#' get data for pedgene run using \bold{Mega2} framework
#'
#' @description
#'  This populates the \bold{Mega2} data.frames from the specified database.  It also
#'  prunes the samples to only include members that have a definite case or control
#'  status.  Undefined samples are ignored; this is necessary for \emph{pedgene}.
#'
#' @param db specifies a \bold{Mega2} SQLite database containing study data.
#'
#' @param verbose default is passed to the \bold{Mega2} framework.  1 indicates that
#'  diagnostic printouts should be enabled.
#'
#' @return None
#' @importFrom mega2 dbmega2_import mkfam setfam setRanges
#' @importFrom utils read.table write.table
#' @export
#'
#' @note
#'  This also records schaidPed and pedPer that are used later in the \emph{Dopedgene} calculation.
#'
#' @examples
#'\dontrun{
#' init_pedgene("ped3.db", verbose = TRUE)
#'}
init_pedgene = function (db = "ped3.db", verbose = FALSE) {

    ENV = dbmega2_import(db, verbose = verbose)

##  NonmissingPheID=get(load("../NonmissingPheID.RData"))
##  This is a rather complicated way to essentiall accomplish the two lines below
    ## non=read.table("ped3.famphe", header=F)
    ## fam = mkfam()
    ## pl=merge(fam[ , c(1, 3, 4)], non[ , 2:3],
    ##          by.x = c("PedPre", "PerPre"), by.y = c("V2", "V3"))

    ## fam = fam[fam[ , 1] %in% pl[ , 3], ]
    ## row.names(fam) = NULL

    fam = mkfam()
    fam = fam[fam$trait != 0, ]
    setfam(fam)  # also updates unified_genotype_table

#   refRange = read.table("../refseq_genes.txt", header= TRUE, stringsAsFactors= FALSE)
#   refRange$cdsStart <-refRange$cdsEnd <- NULL
#   refRange = refRange[nchar(refRange$chrom) <= 5, ]
#   refRange = refRange[refRange$chrom != "chrX" & refRange$chrom != "chrY", ]
#   refRange = refRange[!duplicated(refRange), ]
#   colnames(refRange) = c("XX", "name2", "chrom", "txStart", "txEnd")
#   colnames(refRange) = c("XX", "SYMBOL", "TXCHROM", "TXSTART", "TXEND")
#   row.names(refRange) = NULL
#   write.table(refRange, file=ped3.ref, quote=F, row.names=F)

# if you want your own ranges
#   refRanges = read.table("ped3.ref", header = T)
#   refRanges = refRanges[! duplicated(refRanges$SYMBOL), ]
#   setRanges(refRanges, 3:5)

    ENV$schaidPed = ENV$fam[ , c(-1, -2)]
    colnames(ENV$schaidPed) = c("ped", "person", "father", "mother", "sex", "trait")
    ENV$pedPer = ENV$schaidPed[ , 1:2]
    ENV$mt = matrix(c(11, 12, 21, 22, 0,    0, 1, 1, 2, 0), nrow = 5, ncol = 2)

}

pedgene_results <- data.frame(chr = character(0), gene = character(0), nvariants = numeric(0),
                      start = numeric(0), end = numeric(0),
                      pKernel_BT = numeric(0), pBurden_BT = numeric(0),
                      pKernel_MB = numeric(0), pBurden_MB = numeric(0),
                      pKernel_UW = numeric(0), pBurden_UW = numeric(0),
                      geneID = numeric(0), stringsAsFactors = FALSE)

#' execute the pedgene function on a subset of the gene transcript ranges
#'
#' @param gs a subsequence of the gene transcript ranges to calculate the \emph{Dopedgene} function
#' on.
#'
#' @return None
#' @importFrom mega2 applyFnToRanges
#' @export
#'
#' @examples
#'\dontrun{
#' run()
#'
#' run(1:10)
#'}
run_pedgene = function (gs = 1:100) {

    unlink("k_Schaid_rare.txt")

#   applyFnToRanges(DOpedgene)
    ENV = (environment(dbmega2_import))$ENV
    applyFnToRanges(DOpedgene, ENV$refRanges[gs, ], ENV$refIndices)
}

#' pedgene call back function
#'
#' @description ...
#'
#' @param geno_arg A matrix with one row per \emph{fam} pedigree member and one column for each marker that is selected.  Each datum are the two characters indicating the genotype for the marker/member.
#'
#' @param markers_arg a data.frame with the following 5 variables:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#'
#' @param range_arg one row of a ranges_arg.  The latter is a data.frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @return None
#' @importFrom pedgene pedgene
#' @export
#'
#' @examples
#'\dontrun{
#'    applyFnToRanges(DOpedgene, ENV$refRanges[gs, ], ENV$refIndices)
#'}
DOpedgene = function(geno_arg, markers_arg, range_arg) {

    ENV = (environment(dbmega2_import))$ENV

    markerNames = markers_arg$MarkerName
    gene  <- as.character(range_arg$SYMBOL)

    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    for (k in 1:(di[2])) {
        vec = ENV$mt[match(as.integer(geno_arg[ , k]), ENV$mt), 2]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (ENV$verbose)
            cat(gene, markerNames[k], g0, g1, g2, "\n")
        if (g0 < g2) {
           geno[ , k] = 2 - vec
        } else {
           geno[ , k] =     vec
        }
    }

    geno = matrix(geno, nrow = di[1])
    maf = colMeans(geno)
    pos = markerNames[maf > 0]
    if (length(pos) >= 2) {       # at least 2 non-polymorphic variants #
        geno <- geno[ , maf > 0]     # remove nonpolymorphic variants #
        nsnp    <- ncol(geno)
        weight <- rep(1, ncol(geno))

        pedgeno <- cbind(ENV$pedPer, geno)

        BT <- pedgene(ENV$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= FALSE, method= "kounen")
        pKernel_BT <- BT$pgdf$pval.kernel
        pBurden_BT <- BT$pgdf$pval.burden

        MB <- pedgene(ENV$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= TRUE, method= "kounen")
        pKernel_MB <- MB$pgdf$pval.kernel
        pBurden_MB <- MB$pgdf$pval.burden

        UW <- pedgene(ENV$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= weight, weights.mb= TRUE, method= "kounen", acc.davies=1e-9)
        pKernel_UW <- UW$pgdf$pval.kernel
        pBurden_UW <- UW$pgdf$pval.burden

        ## read out the results ##
        chr   <- as.character(range_arg$TXCHROM)
        start <- range_arg$TXSTART
        end   <- range_arg$TXEND
        zzz = 1
        if (ENV$verbose)
            cat(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz, "\n")

        pedgene_results[1, ] <- c(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                          pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz)

        write.table(pedgene_results, file="k_Schaid_rare.txt", append= TRUE, row.names= FALSE, col.names= FALSE, quote= FALSE)
    }
}
