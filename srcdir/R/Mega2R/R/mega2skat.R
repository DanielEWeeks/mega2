
#   Mega2R: Mega2 for R.
#
#   Copyright 2017, University of Pittsburgh. All Rights Reserved.
#
#   Contributors to Mega2R: Robert V. Baron and Daniel E. Weeks.
#
#   This file is part of the Mega2R program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later
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

# library(SKAT)

#' load Mega2 SQLite database and perform initialization for SKAT usage
#'
#' @description
#'  This populates the \bold{R} data frames from the specified \bold{Mega2} SQLite database.  It then
#'  prunes the samples to only include members that have a definite case or control
#'  status.  Undefined samples are ignored.
#'
#' @param db specifies the path of a \bold{Mega2} SQLite database containing study data.
#'
#' @param filename filename to store results data frame.  By default "SKAT.txt" is used.
#'
#' @param verbose TRUE indicates that diagnostic printouts should be enabled.
#'  This value is saved in the returned environment.
#'
#' @param allMarkers TRUE means use all markers in a given transcript even if there is no
#'  variation.  FALSE means ignore markers that show no variation; this is the default.
#'
#' @return "environment" containing data frames from an SQLite database and some computed values.
#'
#' @export
#'
#' @note
#'  \emph{init_SKAT} creats a data frame, \emph{envir$phe}, of phenotype observations.
#'  In addition, it initializes a matrix to aid
#'   in translating a genotype allele matrix to a genotype count matrix.
#'
#'  It also initializes the data frame \emph{envir$SKAT_results} to zero rows and
#'  can be give a filename to write out the data frame to.
#'
#' @examples
#'\dontrun{
#' init_SKAT("ped3.db", verbose = TRUE)
#'}
init_SKAT = function (db = NULL, filename = NULL, verbose = FALSE, allMarkers = FALSE) {

    if (is.null(db))
        stop("You must specify a database argument!\n", call. = FALSE)

    envir = dbmega2_import(db, verbose = verbose)

    fam = mkfam(envir = envir)
    fam = fam[fam$trait != 0, ]
    setfam(fam, envir = envir)  # also updates unified_genotype_table

    envir$phe  = mkphenotype(envir)

    envir$SKAT_filename = filename

    envir$SKAT_results = data.frame(chr = character(0), gene = character(0),
                                     nvariants = numeric(0), start = integer(0), end = integer(0),
                                     skat = numeric(0),
                                     stringsAsFactors = FALSE)

    envir$mt = matrix(c(11, 12, 21, 22, 0,    0, 1, 1, 2, 9), nrow = 5, ncol = 2)

    envir$allMarkers = allMarkers

   return (envir)
}

#' execute the CRAN SKAT function on a subset of the gene transcripts
#'
#' @description
#' Execute the SKAT function on the first \emph{gs} default gene transcripts (gs = 1:100).
#'  Update the \emph{envir$SKAT_results} data frame with the results.
#'
#' @param f SKAT_Null_Model formula.  If this is non NULL, envir$obj is initialized by calling
#'  SKAT_Null_Model(f, out_type = ty).  If you need to specify additional arguments to the Model
#'  viz. (data, Adjustment, n.Resampling, type.Resampling)
#'  or need to use a different model viz. SKAT_NULL_emmaX, \cr SKAT_Null_Model_ChrX
#'  set the formula to NULL, then before Mega2SKAT is called, build the model you need and
#'  assign it to ENV$obj.
#'
#' @param ty type of phenotype C/D = Continuous/Binary 5 (internal type 1/2)
#'
#' @param gs a subrange of the default transcripts (refRanges) over which to calculate
#'  the \emph{DOSKAT} function.
#'
#' @param genes a list of genes over which to calculate the \emph{DOSKAT} function.
#'  The value, "*", means use all the transcripts in the selected Bioconductor database.
#'  If genes is NULL, the gs range of the internal \emph{refRanges} will be used.
#'
#' @param skat alternate SKAT function, viz. SKATBinary, SKAT_CommonRare.  If it is also
#'  necessary is to pass additional arguments to the SKAT function, they may be added to the end
#'  of the Mega2SKAT function and will be passed.  See examples
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @param ... extra arguments for skat
#'
#' @return None
#'
#' @importFrom utils read.table write.table
#' @importFrom SKAT SKAT SKATBinary SKAT_CommonRare
#' @importFrom SKAT SKAT_Null_Model SKAT_NULL_emmaX SKAT_Null_Model_ChrX
#' @export
#'
#' @note
#'  This code starts by deleting the output file set in \code{init_SKAT} ("SKAT.txt" by default).
#'  Then \code{SKAT_Null_Model} is called if the formula, f, is not NULL.  A helper function
#'  \code{SKAT4arg} is defined for the 4 argument callback function which in turn calls
#'  \code{DOSKAT} with the appropriate arguments (including those additional to the
#'  \code{Mega2SKAT} function).
#'  Finally, the data frame of results, \emph{envir$SKAT_results}, is written
#'  to the output file.
#'
#' @examples
#'\dontrun{
#' Mega2SKAT(ENV$phe ~ 1, "D")
#'
#' Mega2SKAT(ENV$phe ~ 1, "D", gs=1:10)
#'
#' Mega2SKAT(NULL, NULL, gs=1:1000, kernel = "linear.weighted", weights.beta=c(0.5,0.5))
#'}
Mega2SKAT = function (f, ty, gs = 1:100, genes=NULL, skat = SKAT::SKAT, envir = ENV, ...) {

    if (! is.null(envir$SKAT_filename)) unlink(envir$SKAT_filename)

    if (! is.null(f))
        envir$obj = SKAT_Null_Model(f, out_type = ty)
    envir$skat = skat
    SKAT4arg = function(geno_arg, markers_arg, range_arg, envir) {
#       do.call(DOSKAT, (list(geno_arg, markers_arg, range_arg, envir, ...)) )
        DOSKAT(geno_arg, markers_arg, range_arg, envir, ...)
    }

    if (is.null(genes)) 
        applyFnToRanges(SKAT4arg, envir$refRanges[gs, ], envir$refIndices, envir = envir)
    else
        applyFnToGenes(SKAT4arg, genes, envir = envir)

    if (! is.null(envir$SKAT_filename))
        write.table(envir$SKAT_results, file=envir$SKAT_filename,
                    row.names= FALSE, col.names= TRUE, quote= FALSE)
}

#' SKAT call back function
#'
#' @description
#'  Convert the genotype allele  patterns of 1/1, 1/2 (and 2/1) and 2/2 in the genotype matrix
#'  to the numbers 0, 1, 2 for each marker. (Reverse, the order iff allele "1" has the
#'  minor allele frequency.)  Ignore markers that have no variants (unless allMarkers is TRUE).
#'  Finally, invoke \code{SKAT} with the converted genotype matrix, Null model saved in envir$obj,
#'  and any additionally supplied arguments.
#'  Save information about the range and the p.value calculated by \code{SKAT} 
#'  in \emph{envir$SKAT_results}.
#'
#' @param geno_arg A character matrix with one row per \emph{fam} pedigree member and one column
#' foreach marker in markers_arg.
#' Each cell contains the two characters indicating the nucleotides for the marker.
#'
#' @param markers_arg a data.frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#' }
#'
#' @param range_arg one row of a ranges_arg.  The latter is a data frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @param ... extra arguments for SKAT
#'
#' @return None
#'
#' @export
#'
#' @note
#'  This function accumulates output in the data frame, \emph{envir$SKAT_results}.  It will
#'  print out the lines as they are generated if \emph{envir$verbose} is TRUE.  It does not write
#'  the data frame to a file.  You must save the data frame.
#'  You also must initialize the data frame when necessary.
#'
#' @examples
#'\dontrun{
#'    applyFnToRanges(DOSKAT, ENV$refRanges[gs, ], ENV$refIndices, ENV)
#'}
# SKAT(<formula>, <out_type>, kernel = "linear.weighted", weights.beta=c(0.5,0.5))
#
DOSKAT = function(geno_arg, markers_arg, range_arg, envir, ...) {

    lastp1 = nrow(envir$SKAT_results) + 1

    markerNames = markers_arg$MarkerName
    gene  = as.character(range_arg$SYMBOL)

    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    kk = 0
    for (k in 1:(di[2])) {
        vec = envir$mt[match(as.integer(geno_arg[ , k]), envir$mt), 2]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (envir$verbose)
            cat(gene, markerNames[k], g0, g1, g2, "\n")

        if (g0 < g2) {
           if (g2 < (di[1] - 0) || envir$allMarkers) {
               kk = kk + 1
               geno[ , kk] = 2 - vec
           }
        } else {
           if (g0 < (di[1] - 0) || envir$allMarkers) {
               kk = kk + 1
               geno[ , kk] =     vec
           }
        }
    }

    geno = geno[,1:kk]
    nsnps = kk
    if (kk == 0) return (NULL)
    if (kk == 1) geno = matrix(geno, nrow = di[1], ncol = kk)
    
    skat = XSKAT(envir$skat, geno, envir$obj, ...)
#browser(skipCalls=2)

    chr   <- as.character(range_arg$TXCHROM)
    start <- range_arg$TXSTART
    end   <- range_arg$TXEND

    result = list(chr, gene, nsnps, start, end, skat$p.value
                  )
    envir$SKAT_results[lastp1,] = result
    if (envir$verbose) {
        print(envir$SKAT_results[lastp1, ])
    }
}

XSKAT = function(skat, ...) {
    skat(...)
}
