
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

#library(mega2)
#library(GenABEL)

#' Mega2GenABEL package
#'
#' @description This package reads a Mega2 SQLite3 database into R dataframes and
#'	generates input data for GenABEL from these same frames.
#'
#' @author Robert V Baron
#' @docType package
#' @name Mega2GenABEL-package
NULL

#' Load Mega2 database and initialize family structure
#'
#' @description
#'  \emph{dbmega2_import} the specified database and load the environment, \emph{ENV}, with the
#'  table data.  Also run \emph{mkfam} to initialize the family structure and then \emph{setfam}
#'  to modify the \emph{unified_genotype_table} to match the family.  By default this will remove
#'  samples that were replicated to break loops in the pedigree, see \emph{mkfam} for details.
#'
#' @param db specify SQLite database to load
#'
#' @param ... aditional arguments to pass to \emph{dbmega2_import}
#'
#' @return ENV an environment that contains all the tables created from the SQLite tables.
#'
#' @importFrom mega2 read.Mega2DB
#' @export
#'
#' @note This functions just calls the same named function in the mega2 package
#'
#' @examples
#'\dontrun{
#' read.Mega2DB("database.db")
#'}
read.Mega2DB = function(db, ...) {

    return (mega2::read.Mega2DB(db, ...))
}

#' generate a PLINK TPED file for GenABEL
#'
#' @description
#'  Generate a PLINK TPED file from the specified Mega2 SQLite database.  The file is named "prefix".tped
#'  If the markers arg is.null(), the entire ENV$markers set is used otherwise markers arg MUST
#'  be a subset of the ENV$markers data.frame -- same columns, but pruned rows.  
#'
#' @param prefix prefix for .tped file name
#'
#' @param ENV "environment" containing SQLite database and other globals
#'
#' @param markers markers selected to be in output file
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @return None
#'
#' @importFrom mega2 getgenotypes
#' @importFrom utils write.table
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELtped("foo")
#'
#' mkGenABELtped("foo", ENV$markers[ENV$markers$chromosome >= 20,])
#'}
mkGenABELtped = function(prefix, ENV, markers=NULL, mapno = 0) {
    file = paste0(prefix, ".tped")
    
    unlink(file)

    if (is.null(markers)) markers = ENV$markers

    allele_table = ENV$allele_table[ENV$allele_table$locus_link %in% markers$locus_link,]
    map_table = ENV$map_table[ENV$map_table$marker %in% markers$locus_link,]
    M = nrow(markers)
    C = 1000

    block = data.frame(matrix(0, nrow = C, ncol = nrow(ENV$fam) + 4), stringsAsFactors=FALSE)
    blockcol = ncol(block)

    ppl = paste0(ENV$fam$PedPre, "_", ENV$fam$PerPre)
#   ppl = rbind(ppl, ppl)
    names(block) = c("#CHROM", "ID", "GEN", "POS", ppl)

    SVallele1 = ENV$allele_table$AlleleName[ENV$allele_table$index == 1]
    ENV$allele_table$AlleleName[ENV$allele_table$index == 1] =
        paste0(SVallele1, " ")
    SVallele2 = ENV$allele_table$AlleleName[ENV$allele_table$index == 2]
    ENV$allele_table$AlleleName[ENV$allele_table$index == 2] =
        paste0(SVallele2, " ")

    j = 0
    while (TRUE) {
        if (M <= 0) break
        R = ((j*C+1):(j*C + C))
        L = length(R)
        if (M < L) {
            L = M
            R = ((j*C+1):(j*C+L))
        }
        M = M - L
        BR = 1:L

print(system.time ({        
        GPos = map_table[map_table$map==mapno, c("position")][R]
        GPosPos = sprintf("%.2f", GPos)

        block[BR , 1] = markers[R , 4]
        block[BR , 2] = markers[R , 3]
        block[BR , 3] = GPosPos
        block[BR , 4] = markers[R , 5]

        cr = getgenotypes(markers[R, ])             # 7.17%
        a1 = t(cr)                                  # 0.86%

#       di = dim(a1)
#       line = paste(substr(a1, start=1, stop=1), substr(a1, start=2, stop=2))
#       dim(line) = di
#       block[BR, 5:blockcol] = line

        block[BR, 5:blockcol] = a1

      }))

print(system.time ({        
        write.table(block[BR, ], file=file, sep="\t", quote=FALSE,     # 48.49%
                    append=TRUE, row.names=FALSE, col.names=FALSE)
 }))
        j = j + 1
        message(j)
    }

    ENV$allele_table$AlleleName[ENV$allele_table$index == 1] = SVallele1
    ENV$allele_table$AlleleName[ENV$allele_table$index == 2] = SVallele2
}

#' generate required VCF header
#'
#' @description
#'  Call functions to: create .tped file, .tfam file and .phe file.
#'  Call the GenABEL functions to process these files.
#'
#' @param prefix prefix for vcf file name
#'
#' @param markers data.frame of markers being processed
#'
#' @importFrom GenABEL convert.snp.tped load.gwaa.data
#' @importFrom mega2 getENV
#' @return gwaa.class-object of previously read(.Mega2DB) database
#'
#' @export
#'
#' @examples
#'\dontrun{
#' ENV <- read.Mega2DB("my.db")
#'
#' Mega2GenABEL(prefix, NULL)
#'}
Mega2GenABEL = function (prefix, markers = NULL) {

    ENV = getENV()

    if (is.null(markers)) markers = ENV$markers

    mkGenABELtped(prefix, ENV, markers)

    mkGenABELtfam(prefix, ENV, markers)

    convert.snp.tped(tpedfile=paste0(prefix,".tped"),
                     tfamfile=paste0(prefix,".tfam"),
                     outfile=paste0(prefix, "tped.raw"),
                     strand="u")

    mkGenABELphe(prefix, ENV, markers)

    return (load.gwaa.data(phenofile=paste0(prefix,".phe"),
                           genofile=paste0(prefix, "tped.raw"),
                           force = TRUE)
            )
}

#' generate required fam family for PLINK TPED (.tfam) file
#'
#' @description
#'  Generate the six column .tfam file used with the .tped file.  Note: Only the person id
#'  column appears to be used by GenABEL.
#'
#' @param prefix prefix for vcf file name
#'
#' @param ENV "environment" containing SQLite database and other globals
#'
#' @param markers data.frame of markers being processed
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELtfam(prefix, ENV, NULL)
#'}
mkGenABELtfam = function (prefix, ENV, markers) {
    file = paste0(prefix, ".tfam")

#   -9 vs 0 for case/control

# mega2 mis
#   ENV$fam[ENV$fam[ , 8] ==0, 8] = -9
    fam = ENV$fam
    fam[ , "PerPre"] = paste(fam[ , "PedPre"], fam[ , "PerPre"], sep="_")
    
    write.table(fam[, -(1:2)], file=file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required PLINK (.phe) file
#'
#' @description
#'  Generate the .phe (phenotype) file for PLINK which is used by GenAbel.  The person
#'  must match that specified in the .tfam file
#'
#' @param prefix prefix for vcf file name
#'
#' @param ENV "environment" containing SQLite database and other globals
#'
#' @param markers data.frame of markers being processed
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELphe(prefix, ENV, NULL)
#'}
mkGenABELphe = function (prefix, ENV, markers) {
    file = paste0(prefix, ".phe")

    unlink(file)

    fam = ENV$fam
    fam$id = paste(fam[ , "PedPre"], fam[ , "PerPre"], sep="_")
   
# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = data.frame(fam$id, stringsAsFactors=FALSE)
    names(out) = "id"
# This (NA / 1 male / 0 female) is different from linkage: NA / 1 male / 2 female
    out$sex = c(NA, 1, 0)[fam$Sex + 1]
    hdr = c("id", "sex")
   
    phenotype_table = ENV$phenotype_table

    raw = unlist(ENV$phenotype_table[,4])
    raw = matrix(raw, ncol=8, byrow=T)
    nrows     = nrow(raw)
    nrowpheno = nrow(out)

    for (i in 1:ENV$PhenoCnt) {
        hdr = c(hdr, ENV$locus_table[i, 2]) # 2 == LocusName

# phenotype_table contains a blob which is a list of entries.  An entry is either an 8 byte
#  double for quant, or two 4 byte ints for affect
        if (ENV$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            col = vector("integer", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[ENV$PhenoCnt*(j-1)+i, 1:4], integer(), n=1, size=4)
            }
#           col[col==0] = NA
            out$col = col
            names(out) = hdr
        } else if (ENV$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            col = vector("numeric", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[ENV$PhenoCnt*(j-1)+i, 1:8], numeric(), n=1, size=8)
            }
            col[col==-99] = NA
            out$col = col
            names(out) = hdr
        }
    }

    cat(hdr,  file=file, sep="\t")
    cat("\n", file=file, append=TRUE)

    write.table(out, file=file, sep="\t", quote=FALSE, append=TRUE,
                row.names=FALSE, col.names=FALSE)
}
