
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

#' generate gwaa.data-class object
#'
#' @description
#'  Call functions to: create .tped file, .tfam file and .phe file.
#'  Call the GenABEL functions to process these files which will additionally
#'  create a GenABEL .raw file.
#'
#' @param prefix prefix for generated file names
#'
#' @param markers data frame of markers being processed
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return gwaa.class-object of previously read.Mega2DB database
#'
#' @importFrom GenABEL convert.snp.tped load.gwaa.data
#' @export
#'
#' @examples
#'\dontrun{
#' ENV <- read.Mega2DB("my.db")
#'
#' gwaa = Mega2GenABEL(prefix, NULL)
#' str(gwaa)
#' head(summary(gwaa))
#'}
Mega2GenABEL = function (prefix, markers = NULL, mapno = 0, envir = ENV) {

    if (is.null(markers)) markers = envir$markers

    mkGenABELtped(prefix, markers, mapno = mapno, envir)

    mkGenABELtfam(prefix, envir)

    convert.snp.tped(tpedfile=paste0(prefix,".tped"),
                     tfamfile=paste0(prefix,".tfam"),
                     outfile=paste0(prefix, "tped.raw"),
                     strand="u")

    mkGenABELphe(prefix, envir)

    return (load.gwaa.data(phenofile=paste0(prefix,".phe"),
                           genofile=paste0(prefix, "tped.raw"),
                           force = TRUE)
            )
}

#' generate a PLINK TPED file for GenABEL
#'
#' @description
#'  Generate a PLINK TPED file from the specified Mega2 SQLite database.  The file is named "prefix".tped
#'  If the markers arg is.null(), the entire envir$markers set is used otherwise markers arg MUST
#'  be a subset of the envir$markers data.frame -- same columns, but pruned rows.  
#'
#' @param prefix prefix for .tped file name
#'
#' @param markers markers selected to be in output file
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @param envir "environment" containing SQLite database and other globals
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
mkGenABELtped = function(prefix, markers=NULL, mapno = 0, envir) {
    file = paste0(prefix, ".tped")
    
    unlink(file)

    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    map_table = envir$map_table[envir$map_table$marker %in% markers$locus_link,]
    M = nrow(markers)
    C = 1000

    block = data.frame(matrix(0, nrow = C, ncol = nrow(envir$fam) + 4), stringsAsFactors=FALSE)
    blockcol = ncol(block)

    ppl = paste0(envir$fam$PedPre, "_", envir$fam$PerPre)
#   ppl = rbind(ppl, ppl)
    names(block) = c("#CHROM", "ID", "GEN", "POS", ppl)

    SVallele1 = envir$allele_table$AlleleName[envir$allele_table$index == 1]
    envir$allele_table$AlleleName[envir$allele_table$index == 1] =
        paste0(SVallele1, " ")
    SVallele2 = envir$allele_table$AlleleName[envir$allele_table$index == 2]
    envir$allele_table$AlleleName[envir$allele_table$index == 2] =
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

## print(system.time ({        
        GPos = map_table[map_table$map==mapno, c("position")][R]
        GPosPos = sprintf("%.2f", GPos)

        block[BR , 1] = markers[R , 4]
        block[BR , 2] = markers[R , 3]
        block[BR , 3] = GPosPos
        block[BR , 4] = markers[R , 5]

        cr = getgenotypes(markers[R, ], envir = envir ) # 7.17%
        a1 = t(cr)                                      # 0.86%

#       di = dim(a1)
#       line = paste(substr(a1, start=1, stop=1), substr(a1, start=2, stop=2))
#       dim(line) = di
#       block[BR, 5:blockcol] = line

        block[BR, 5:blockcol] = a1

##      }))

## print(system.time ({        
        write.table(block[BR, ], file=file, sep="\t", quote=FALSE,     # 48.49%
                    append=TRUE, row.names=FALSE, col.names=FALSE)
## }))
        j = j + 1
        if (envir$verbose) message(".", appendLF = FALSE)
    }

    envir$allele_table$AlleleName[envir$allele_table$index == 1] = SVallele1
    envir$allele_table$AlleleName[envir$allele_table$index == 2] = SVallele2
}

#' generate required fam family for PLINK TPED (.tfam) file
#'
#' @description
#'  Generate the six column .tfam file used with the .tped file.  Note: Only the person id
#'  column appears to be used by GenABEL.
#'
#' @param prefix prefix for generated file name
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELtfam(prefix, envir)
#'}
mkGenABELtfam = function (prefix, envir) {
    file = paste0(prefix, ".tfam")

#   -9 vs 0 for case/control

# mega2 mis
#   envir$fam[envir$fam[ , 8] ==0, 8] = -9
    fam = envir$fam
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
#' @param prefix prefix for generated file name
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELphe(prefix, envir)
#'}
mkGenABELphe = function (prefix, envir) {
    file = paste0(prefix, ".phe")

    unlink(file)

    fam = envir$fam
    fam$id = paste(fam[ , "PedPre"], fam[ , "PerPre"], sep="_")
   
# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = data.frame(fam$id, stringsAsFactors=FALSE)
    names(out) = "id"
# This (NA / 1 male / 0 female) is different from linkage: NA / 1 male / 2 female
    out$sex = c(NA, 1, 0)[fam$Sex + 1]
    hdr = c("id", "sex")
   
    phenotype_table = envir$phenotype_table

    raw = unlist(envir$phenotype_table[,4])
    raw = matrix(raw, ncol=8, byrow=T)
    nrows     = nrow(raw)
    nrowpheno = nrow(out)

    for (i in 1:envir$PhenoCnt) {
        hdr = c(hdr, envir$locus_table[i, 2]) # 2 == LocusName

# phenotype_table contains a blob which is a list of entries.  An entry is either an 8 byte
#  double for quant, or two 4 byte ints for affect
        if (envir$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            col = vector("integer", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:4], integer(), n=1, size=4)
            }
#           col[col==0] = NA
            out$col = col
            names(out) = hdr
        } else if (envir$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            col = vector("numeric", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:8], numeric(), n=1, size=8)
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
