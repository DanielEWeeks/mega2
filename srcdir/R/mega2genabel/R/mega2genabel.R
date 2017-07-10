
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

#library(mega2r)
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

    file = paste0(prefix, ".phe")
    unlink(file)
    out = mkGenABELphe(envir)
    write.table(out, file=file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=TRUE)


#x  return (load.gwaa.data(phenofile=paste0(prefix,".phe"),
    return (gwaaO(phenofile=paste0(prefix,".phe"),
                           genofile=paste0(prefix, "tped.raw"),
                           force = TRUE,
                           envir = envir)
            )
}

#' generate gwaa.data-class object
#'
#' @description
#'  Directly write a gwaa.data-class object from ENV tables
#'
#' @param markers data frame of markers being processed
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return gwaa.class-object of previously read.Mega2DB database
#'
#' @export
#'
#' @examples
#'\dontrun{
#' ENV <- read.Mega2DB("my.db")
#'
#' gwaa = Mega2ENVGenABEL(NULL)
#' str(gwaa)
#' head(summary(gwaa))
#'}
Mega2ENVGenABEL = function (markers = NULL, force = TRUE, makemap = FALSE,
                         sort = TRUE, mapno = 0, envir = ENV) {
#browser()
    if (is.null(markers)) markers = envir$markers

    gwaa(markers = markers, force = force, 
                   makemap = makemap, sort = sort, 
                   envir = envir)

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
#' @importFrom mega2r getgenotypes
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

        cr = getgenotypes(markers[R, ], sepstr = " ", envir = envir ) # 7.17%
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
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' mkGenABELphe(envir)
#'}
mkGenABELphe = function (envir) {

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

#    cat(hdr,  file=file, sep="\t")
#    cat("\n", file=file, append=TRUE)

#   write.table(out, file=file, sep="\t", quote=FALSE, append=TRUE,
#               row.names=FALSE, col.names=FALSE)
    out
}

#' generate GenABEL coding vector 
#'
#' @description
#'  Each element of the coding vector is a lookup in the GenABEL:::alleleID.codes() array.
#'  The alleles are ordered so that with the greater frequency appears first.
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' Mega2GenABELcoding(envir)
#'}
Mega2GenABELcoding = function(markers = NULL, envir = ENV) {
    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    mm = merge(x=allele_table[allele_table$indexX == 1,],
               y=allele_table[allele_table$indexX == 2,],
               by="locus_link")
    nn=ifelse(mm$Frequency.x > mm$Frequency.y,
              paste0(mm$AlleleName.x, mm$AlleleName.y), 
              paste0(mm$AlleleName.y, mm$AlleleName.x))
    envir$xGTy = mm$Frequency.x > mm$Frequency.y
    envir$Frx = mm$Frequency.x
    envir$Fry = mm$Frequency.y

    nn[mm$Frequency.x == 0 & mm$Frequency.y == 0] = '12'

    fx = mm$Frequency.x == mm$Frequency.y
#   nn[fx] = paste0(mm[fx, "AlleleName.x"], mm[fx, "AlleleName.y"])

    fy = mm$Frequency.x == 0 & mm$Frequency.y == 1
    nn[fy] = paste0(mm[fy, "AlleleName.y"], mm[fy, "AlleleName.y"])

    fz = mm$Frequency.x == 1 & mm$Frequency.y == 0
    nn[fz] = paste0(mm[fz, "AlleleName.x"], mm[fz, "AlleleName.x"])

#this is what genabel does for 0/0
    nn[nn=="00"] = "12"
    cc = alleleID.codes()
    ar = match(nn, cc, nomatch=NA)
    if (any(is.na(ar))) {
        warning("Problems mapping genotypes to GenABEL encoding ", mm[is.na(ar)], "/n")
        ar[is.na(ar)] = which(cc == '12')
    }
    as.raw(ar)
}

#' generate GenABEL compressed genotype matrix
#'
#' @description
#'  The matrix is (# of samples / 4 ) x (# of markers).  (Round samples to multiple of 4.
#'  Each byte stores data for 4 samples; intfn() generates the 4 - 2 bit encodings.
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @importFrom mega2r getgenotypesraw
#' @importFrom stats aggregate
#' @export
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' Mega2GenABELconvert(envir)
#'}
Mega2GenABELconvert = function(markers = NULL, envir = ENV) {
# browser("convert")
    if (is.null(markers)) markers = envir$markers
    nmarkers = nrow(markers)

    y=rep(1:(ceiling(nrow(envir$fam)/4)), each=4)
    zz = getgenotypesraw(markers, envir=envir)
    fil= rep(0, length(y)-nrow(zz))
    rag = matrix(raw(0), nrow = length(y)/4, ncol = nmarkers)

 print (system.time ({        
#   for (m in 1:nrow(markers)) 
    for (ms in seq(1, nmarkers, 100)) {
#
      cat(ms, "  ", sep= " ");  print (system.time ({        
        for (m in seq(ms, ms+100-1, 1)) {
          if (m > nmarkers) break
            intf = intfn(m, envir = envir)
            z = c(zz[ , m], fil)
            gg = aggregate(z, by=list(y), intf)
            rag[ , m] = as.raw(gg$x)
        }
#
      }))
    }
 }))
    rag
}

intfn = function(m, envir = ENV) {
    if (envir$Frx[m] > envir$Fry[m]) {
        a = 3; b = 1
#    } else if (envir$Frx[m] == envir$Fry[m]) {
#        a = 3; b = 1
    } else {
        a = 1; b = 3
    }
    fun = function (v) {
        B2 = 0
#       for (el in v) {
        {
            el = v[1]
            if (el == 131074) #22
                b2 = a
            else if (el == 65537) #11
                b2 = b
            else if (el == 131073 || el == 65538) #21 #12
                b2 = 2
            else if (el == 0)
                b2 = 0
            B2 = bitwOr(B2, bitwShiftL(b2, 6))
        }
#       for (el in v) {
        {
            el = v[2]
            if (el == 131074) #22
                b2 = a
            else if (el == 65537) #11
                b2 = b
            else if (el == 131073 || el == 65538) #21 #12
                b2 = 2
            else if (el == 0)
                b2 = 0
            B2 = bitwOr(B2, bitwShiftL(b2, 4))
        }
#       for (el in v) {
        {
            el = v[3]
            if (el == 131074) #22
                b2 = a
            else if (el == 65537) #11
                b2 = b
            else if (el == 131073 || el == 65538) #21 #12
                b2 = 2
            else if (el == 0)
                b2 = 0
            B2 = bitwOr(B2, bitwShiftL(b2, 2))
        }
#       for (el in v) {
        {
            el = v[4]
            if (el == 131074) #22
                b2 = a
            else if (el == 65537) #11
                b2 = b
            else if (el == 131073 || el == 65538) #21 #12
                b2 = 2
            else if (el == 0)
                b2 = 0
            B2 = bitwOr(B2, b2)
        }
    B2
    }
}

#' @importFrom GenABEL snp.data
gwaa = function (markers = NULL, force = TRUE, 
    makemap = FALSE, sort = TRUE, id = "id", envir = ENV)
{
    if (is.null(markers)) markers = envir$markers

    dta = mkGenABELphe(envir = envir)
    dta = gwaaCheckPhe(dta, id)
    ids = paste(envir$fam$PedPre, envir$fam$PerPre, sep="_")
    nids <- length(ids)
    nbytes <- ceiling(nids/4)
    cat("ids loaded...\n")

    mnams = markers$MarkerName
    nsnps <- length(mnams)
    cat("marker names loaded...\n")

    chrom = as.character(markers$chromosome)
    chrom <- as.factor(chrom)
    gc(verbose = FALSE)
    cat("chromosome data loaded...\n")

    pos = markers$position
    cat("map data loaded...\n")
  ver=1 # ?? ??
    if (ver == 0) {
        coding <- new("snp.coding", as.raw(rep(1, length(pos))))
        strand <- new("snp.strand", as.raw(rep(0, length(pos))))
    }
#    else
    {
        coding = Mega2GenABELcoding(markers = markers, envir=envir)
        class(coding) <- "snp.coding"
        cat("allele coding data loaded...\n")

        strand  = raw(nsnps)
        class(strand) <- "snp.strand"
        cat("strand data loaded...\n")
    }

    rdta = Mega2GenABELconvert(markers = markers, envir=envir)
    dim(rdta) <- c(nbytes, nsnps)
 #?? ?? generate compress on person NOT marker
    rdta <- new("snp.mx", rdta)

    gc(verbose = FALSE)
    dta = gwaaCheckPersons(dta, ids)

    gc(verbose = FALSE)
    a <- snp.data(nids = nids, rawdata = rdta, idnames = ids, 
                  snpnames = mnams, chromosome = chrom, map = pos, coding = coding, 
                  strand = strand, male = dta$sex)
    cat("snp.data object created...\n")

    rm(rdta, ids, mnams, chrom, pos, coding, strand)
    gc(verbose = FALSE)

    gwaaEpilog(a, dta, force, makemap, sort)
}
