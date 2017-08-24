
#   Mega2R: Mega2 for R.
#
#   Copyright 2017, University of Pittsburgh. All Rights Reserved.
#
#   Contributors to Mega2R: Robert V. Baron and Daniel E. Weeks.
#
#   This file is part of the Mega2R program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#
#   Mega2R is distributed in the hope that it will be useful, but WITHOUT
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

#library(pedgene)

#' load Mega2 SQLite database and perform initialization for pedgene usage
#'
#' @description
#'  This populates the \bold{R} data frames from the specified \bold{Mega2 R} database.  It then
#'  prunes the samples to only include members that have a definite case or control
#'  status.  Undefined samples are ignored; this is necessary for CRAN \code{pedgene}.
#'
#' @param db specifies the path of a \bold{Mega2} SQLite database containing study data.
#'
#' @param filename filename to store results data frame.  By default "pedgene.txt" is used.
#'
#' @param verbose TRUE indicates that diagnostic printouts should be enabled.
#'  This value is saved in the returned environment.
#'
#' @return "environment" containing data frames from an SQLite database and some computed values.
#'
#' @importFrom utils read.table write.table
#' @export
#'
#' @note
#'  \emph{init_pedgene} calculates schaidPed and pedPer that are used later in the \emph{Dopedgene} calculation.
#'
#'  It also initializes the dataframe \emph{envir$pedgene_results} to zero rows.
#'
#' @examples
#'\dontrun{
#' init_pedgene("ped3.db", verbose = TRUE)
#'}
init_pedgene = function (db = NULL, filename = NULL, verbose = FALSE) {

    if (is.null(db))
        stop("You must specify a database argument!\n", call. = FALSE)

    envir = dbmega2_import(db, verbose = verbose)

    fam = mkfam(envir = envir)
    fam = fam[fam$trait != 0, ]
    setfam(fam, envir = envir)  # also updates unified_genotype_table

    envir$schaidPed = envir$fam[ , c(-1, -2)]
    colnames(envir$schaidPed) = c("ped", "person", "father", "mother", "sex", "trait")
    envir$pedPer = envir$schaidPed[ , 1:2]
    envir$mt = matrix(c(11, 12, 21, 22, 0,    0, 1, 1, 2, 0), nrow = 5, ncol = 2)

    if (! is.null(filename))
        envir$pedgene_filename = filename
    else
        envir$pedgene_filename = "pedgene.txt"

    envir$pedgene_results <- data.frame(chr = character(0), gene = character(0),
                                        nvariants = numeric(0),
                                        start = numeric(0), end = numeric(0),
                                        sKernel_BT = numeric(0), pKernel_BT = numeric(0),
                                        sBurden_BT = numeric(0), pBurden_BT = numeric(0),
#                                       call_BT    = character(0),

                                        sKernel_MB = numeric(0), pKernel_MB = numeric(0),
                                        sBurden_MB = numeric(0), pBurden_MB = numeric(0),
#                                       call_MB    = character(0),

                                        sKernel_UW = numeric(0), pKernel_UW = numeric(0),
                                        sBurden_UW = numeric(0), pBurden_UW = numeric(0),
#                                       call_UW    = character(0),

                                        stringsAsFactors = FALSE)
    return (envir)
}


#' execute the CRAN pedgene function on a subset of the default gene transcript ranges
#'
#' @description
#' Execute the pedgene function on the first \emph{gs} default gene transcript ranges (gs = 1:100).
#'  Update the \emph{envir$pedgene_results} data frame with the results.
#"
#' @param gs a subrange of the default transcript ranges over which to calculate the \emph{Dopedgene} function.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#' @export
#'
#' @note
#'  This code starts by deleting the output file set in \code{init_pedgene} ("pedgene.txt" by default).  Then \code{Dopedgene}
#'  is applied to all the appropriate ranges.  Finally, the data frame of results, \emph{envir$pedgene_results}, is written
#'  to the output file.
#'
#' @examples
#'\dontrun{
#' run_pedgene()
#'
#' run_pedgene(1:10)
#'}
run_pedgene = function (gs = 1:100, envir = ENV) {

    unlink(envir$pedgene_filename)

    applyFnToRanges(DOpedgene, envir$refRanges[gs, ], envir$refIndices, envir = envir)

    write.table(envir$pedgene_results, file=envir$pedgene_filename,
                row.names= FALSE, col.names= TRUE, quote= FALSE)
}

#' pedgene call back function
#'
#' @description
#'  First, ignore call backs that have less than two markers.  Second, convert the genotypes
#'  patterns of 1/1, 1/2 (and 2/1) and 2/2 in the genotype matrix
#'  to the numbers 0, 1, 2 for each marker. (Reverse, the order iff allele "1" has the
#'  minor allele frequency.)  Finally, prepend the pedigree and person columns of the family data
#'  to processed genotype matrix.  Finally, invoke \code{pedgene} with the family data and converted
#'  genotype matrix for several different weights.  Save the kernel and burden, value and p-value for each
#'  measurement in \emph{envir$pedgene_results}.
#'
#' @param geno_arg A character matrix with one row per \emph{fam} pedigree member and one column for each marker in
#'  markers_arg.  Each cell contains the two characters indicating the nucleotides for the marker.
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
#' @param range_arg one row of a ranges_arg.  The latter is a data frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#' @importFrom pedgene pedgene
#' @export
#'
#' @note
#'  This function accumulates output in the data frame, \emph{envir$pedgene_results}.  It will
#'  print out the lines as they are generated if \emph{envir$verbose} is TRUE.  It does not write anything
#'  to a file.  You must save the data frame or the "observations" you need by yourself.
#'
#' @examples
#'\dontrun{
#'    applyFnToRanges(DOpedgene, ENV$refRanges[gs, ], ENV$refIndices, ENV)
#'}
DOpedgene = function(geno_arg, markers_arg, range_arg, envir = ENV) {

    markerNames = markers_arg$MarkerName
    gene  <- as.character(range_arg$SYMBOL)

    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    for (k in 1:(di[2])) {
        vec = envir$mt[match(as.integer(geno_arg[ , k]), envir$mt), 2]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (envir$verbose)
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

        pedgeno <- cbind(envir$pedPer, geno)

        BT <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= FALSE, method= "kounen")
        sKernel_BT <- BT$pgdf$stat.kernel
        pKernel_BT <- BT$pgdf$pval.kernel
        sBurden_BT <- BT$pgdf$stat.burden
        pBurden_BT <- BT$pgdf$pval.burden
#       call_BT    <- BT$call

        MB <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= TRUE, method= "kounen")
        sKernel_MB <- MB$pgdf$stat.kernel
        pKernel_MB <- MB$pgdf$pval.kernel
        sBurden_MB <- MB$pgdf$stat.burden
        pBurden_MB <- MB$pgdf$pval.burden
#       call_MB    <- MB$call

        UW <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= weight, weights.mb= TRUE, method= "kounen", acc.davies=1e-9)
        sKernel_UW <- UW$pgdf$stat.kernel
        pKernel_UW <- UW$pgdf$pval.kernel
        sBurden_UW <- UW$pgdf$stat.burden
        pBurden_UW <- UW$pgdf$pval.burden
#       call_UW    <- UW$call

        ## read out the results ##
        chr   <- as.character(range_arg$TXCHROM)
        start <- range_arg$TXSTART
        end   <- range_arg$TXEND

        result = list(chr, gene, nsnp, start, end,
                   sKernel_BT, pKernel_BT, sBurden_BT, pBurden_BT,
                   sKernel_MB, pKernel_MB, sBurden_MB, pBurden_MB,
                   sKernel_UW, pKernel_UW, sBurden_UW, pBurden_UW)
        lastp1 = nrow(envir$pedgene_results) + 1
        envir$pedgene_results[lastp1,] = result
        if (envir$verbose) {
            print(envir$pedgene_results[lastp1, ])
        }

    } else {
        if (envir$verbose)
            message("Only one markers in range.  Ignored!\n")
    }
}
