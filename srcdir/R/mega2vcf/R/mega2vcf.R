
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
#'	generates a VCF file from these same frames.
#'
#' @author Robert V Baron
#' @docType package
#' @name Mega2VCF-package
NULL

#' generate a VCF file
#'
#' @description
#'  Generate a VCF file from the specified Mega2 SQLite database.  The file is named "prefix".vcf
#'  If the markers arg is.null(), the entire ENV$markers set is used otherwise markers arg MUST
#'  be a subset of the ENV$markers data.frame -- same columns, but pruned rows.  
#'
#' @param prefix prefix for vcf file name
#'
#' @param markers markers selected to be in output file
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @param allowFlip REF/ALT of higher frequency occurs first
#'
#' @return None
#'
#' @importFrom mega2 getENV getgenotypesraw 
#' @importFrom utils write.table
#' @export
#'
#' @examples
#'\dontrun{
#' read.Mega2DB("my.db")
#'
#' Mega2VCF("foo")
#'
#' Mega2VCF("foo", ENV$markers[ENV$markers$chromosome >= 20,])
#'}
Mega2VCF = function(prefix, markers=NULL, mapno = 0, allowFlip = FALSE) {

    unlink(paste0(prefix, ".vcf"))

    ENV = getENV()

    if (is.null(markers)) markers = ENV$markers

    mkVCFhdr(prefix, ENV, markers)

    allele_table = ENV$allele_table[ENV$allele_table$locus_link %in% markers$locus_link,]
    map_table = ENV$map_table[ENV$map_table$marker %in% markers$locus_link,]
    M = nrow(markers)
    C = 1000

    QUAL   = rep(".",    times=C)
    FILTER = rep("PASS", times=C)
    FORMAT = rep("GT",   times=C)

    n1 = allele_table[allele_table$indexX==1,]
    n2 = allele_table[allele_table$indexX==2,]
## Mega2 "mis-feature"
    n1$AlleleName[n1$AlleleName %in% c("dummy")] = '.'
    n2$AlleleName[n2$AlleleName %in% c("dummy")] = '.'

    j = 0
    while (TRUE) {
#      if (M != ENV$LocusCnt - ENV$PhenoCnt) break
        if (M <= 0) break
        N = ((j*C+1):(j*C + C))
        L = length(N)
        if (M < L) {
            L = M
            N = ((j*C+1):(j*C+L))
        }
        M = M - L

print(system.time ({        
        chrm = markers[N, c("chromosome", "position", "MarkerName")]
        names(chrm) = c("#CHROM", "POS", "ID")

        REF = n1$AlleleName[N]
        RF  = n1$Frequency[N]
        ALT = n2$AlleleName[N]
        AF  = n2$Frequency[N]
        XXF = REF
        XF  = RF
        
        flip = RF < AF
        doFlip = sum(flip)
        if (allowFlip && doFlip) {
            whichFlip = which(flip)
            REF[whichFlip] = ALT[whichFlip]
            ALT[whichFlip] = XXF[whichFlip]

            RF[whichFlip] = AF[whichFlip]
            AF[whichFlip] = XF[whichFlip]
        }

        GPos = map_table[map_table$map==mapno, c("position", "pos_female", "pos_male")][N, ]
        GPosPos = sprintf("%.2f", GPos$position)
        GPosFem = rep(".", L)
        GPosFem[GPos$pos_female != -99.99] = sprintf("%f", GPos$pos_female)
        GPosMal = rep(".", L)
        GPosMal[GPos$pos_male   != -99.99] = sprintf("%f", GPos$pos_male)

#   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        chrm = cbind(chrm,
                     REF,
                     ALT,
                     QUAL[1:L],
                     FILTER[1:L],
                     INFO=paste0("CM=", GPosPos, ",", GPosFem, ",", GPosMal,
                                 ";RF=", sprintf("%f", RF),
                                 ";AF=", sprintf("%f", AF),
## Mega2 "mis-feature"
##                               ifelse(flip, ",", ""),
                                 ";"),
                     FORMAT[1:L]
               )

        cr = getgenotypesraw(markers[N, ])
        a1 = t(cr)
        a2 = a1
        dm = dim(a1)
        a1 = bitwShiftR(a1, 16)
        attr(a1, "dim") = dm
        x3 = a1
        a2 = bitwAnd(a2, 65535)
        attr(a2, "dim") = dm
        if (allowFlip && doFlip) {
            a1[whichFlip, ] = match(a1[whichFlip, ], c(2, 1), nomatch=0)
            a2[whichFlip, ] = match(a2[whichFlip, ], c(2, 1), nomatch=0)
        } 

        a3 = as.character(a1-1)
        a3[a1 == 0] = "."
        a4 = as.character(a2-1)
        a4[a2 == 0] = "."

        a5 = paste0(a3, "/", a4)
        attr(a5, "dim") = dm
        a6 = cbind(chrm, a5)
        
        if (j == 0) {
            cat(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"),
                paste0(ENV$fam$PedPre, "_", ENV$fam$PerPre),
                sep="\t", append=TRUE, file=paste0(prefix, ".vcf"))
            cat("\n", append=TRUE, file=paste0(prefix, ".vcf"))
        }

## Mega2 "mis-feature"
        a6$POS = paste0(a6$POS, " ")
## Mega2 "mis-feature"
        ncols = ncol(a6)
        a6[, ncols] = paste0(a6[, ncols], "\t")
 }))

print(system.time ({        
        write.table(a6, file=paste0(prefix, ".vcf"), sep="\t", quote=FALSE,
                    append=TRUE, row.names=FALSE, col.names=FALSE)
 }))
        j = j + 1
        message(j)
    }
}
vcf = function(fil = "foo", ...) system.time(Mega2VCF(fil, ...))

#' generate required VCF header
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFhdr(prefix, ENV, NULL)
#'}
mkVCFhdr = function (prefix, ENV, markers) {
    file = paste0(prefix, ".vcf")

    if (is.null(markers)) markers = ENV$markers

    mkVCFfam(prefix, ENV, markers)
    mkVCFfreq(prefix, ENV, markers)
    mkVCFmap(prefix, ENV, markers)
    mkVCFpen(prefix, ENV, markers)
    mkVCFphe(prefix, ENV, markers)

    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)
    cat('##filedate=20170407\n', file=file, append=TRUE)
    cat('##source=MEGA2\n', file=file, append=TRUE)
    cat('##INFO=<ID=CM,Number=3,Type=Float,Description="Genetic Distance in centimorgans (avg, male, female)">\n', file=file, append=TRUE)
    cat('##INFO=<ID=RF,Number=1,Type=Float,Description="Allele Frequency of reference allele">\n', file=file, append=TRUE)
    cat('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency of alternate allele(s)">\n', file=file, append=TRUE)
    cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file=file, append=TRUE)
    cat('##FILTER=<ID=PASS,Description="Passed variant FILTERs">\n', file=file, append=TRUE)

    maxes = sapply(split(markers, markers$chromosome), function(x) max(x$position)+1)
    for (i in 1:length(maxes)) {
        cat('##contig=<ID=', names(maxes)[i], ',length=', maxes[i], ',assembly=B37>\n',
            file=file, append=TRUE, sep="")
    }
}

#' generate required VCF family (.fam) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFfam(prefix, ENV, NULL)
#'}
mkVCFfam = function (prefix, ENV, markers) {
    file = paste0(prefix, ".fam")

#   -9 vs 0 for case/control

# mega2 mis
    ENV$fam[ENV$fam[ , 8] ==0, 8] = -9
    write.table(ENV$fam[, -(1:2)], file=file, sep="\t", quote=FALSE,
                    row.names=FALSE, col.names=FALSE)
}

#' generate required VCF frequency (.freq) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFfreq(prefix, ENV, NULL)
#'}
mkVCFfreq = function (prefix, ENV, markers) {
    file = paste0(prefix, ".freq")

#    unlink(file)

    cat('Name\tAllele\tFrequency\n', file=file)
    allele_pheno = merge(ENV$locus_table, ENV$allele_table[1:2,], by="locus_link")
    allele_pheno = allele_pheno[, c("LocusName", "indexX", "Frequency")]
#std
    allele_pheno$Frequency = sprintf("%.4f", allele_pheno$Frequency)
    write.table(allele_pheno,
                file=file, sep="\t", append=TRUE, quote=FALSE,
                row.names=FALSE, col.names=FALSE)

    allele_table = ENV$allele_table[ENV$allele_table$locus_link %in% markers$locus_link,
                                    c("locus_link", "indexX", "Frequency")]
    alleles = merge(markers[, c("locus_link", "MarkerName")],
                    allele_table[, c("locus_link", "indexX", "Frequency")], by="locus_link")
#std
    alleles = alleles[alleles$Frequency != 0, ]
    alleles$Freq4 = sprintf("%.4f", alleles$Frequency)
    write.table(alleles[ , c(-1, -4)],
                file=file, sep="\t", append=TRUE, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required Mega2 map (.map) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFmap(prefix, ENV, NULL)
#'}
mkVCFmap = function (prefix, ENV, markers) {
    file = paste0(prefix, ".map")

    unlink(file)

    map_table = ENV$map_table[ENV$map_table$marker %in% markers$locus_link,]
    mapnames_table = ENV$mapnames_table
    TBL = markers[, c("chromosome", "MarkerName")]
    hdr = paste("Chromosome", "Name", sep="\t")
    
    for (m in mapnames_table$map) {
        if ( (mapnames_table[m+1, "male_sex_map"] == 0) &
             (mapnames_table[m+1, "female_sex_map"] == 0) ) {
            POS  = map_table[map_table$map == m, "position"]

            if (mapnames_table[m+1, "sex_averaged_map"] == 0) {
                TBL  = cbind(TBL, POS)
                hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.p')
            } else {
#std
                POSS = sprintf("%.6f", POS)
                TBL  = cbind(TBL, POSS)
                hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.a')
            }
          }

        if ( mapnames_table[m+1, "female_sex_map"] != 0) {
            POSF  = map_table[map_table$map == m, "pos_female"]
#std
            POSFS = sprintf("%.6f", POSF)
            TBL  = cbind(TBL, POSFS)
            hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.f')
        }

        if ( mapnames_table[m+1, "male_sex_map"] != 0) {
            POSM  = map_table[map_table$map == m, "pos_male"]
#std
            POSMS = sprintf("%.6f", POSM)
            TBL  = cbind(TBL, POSMS)
            hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.m')
        }
    }
                                 
    cat(hdr, "\n", file=file, sep="")
    
    write.table(TBL, file=file, sep="\t", quote=FALSE, append=TRUE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required Mega2 penetrance (.pen) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFpen(prefix, ENV, NULL)
#'}
mkVCFpen = function (prefix, ENV, markers) {
    file = paste0(prefix, ".pen")

    unlink(file)

    cat('Name\tClass\tPen.11\tPen.12\tPen.22\tType\n', file=file, append=TRUE)
    
    all = merge(merge(ENV$locus_table[1:ENV$PhenoCnt,], ENV$traitaff_table, by="locus_link"),
                ENV$affectclass_table, by="locus_link")
    ord = order(all$locus_link, all$class_link)
    for (i in ord) {
#       malepen   = all[i, ]$MalePen[[1]]
#       mpen  = readBin(malepen, numeric(), n = length(malepen)/8, size = 8, endian = .Platform$endian)
#       mpens = sprintf("%.4f", mpen)
#       cat(all[i, ]$LocusName, all[i, ]$class_link+1, mpens, file=file, append=TRUE, sep="\t")
#       cat("\t\tmale\n", file=file, append=TRUE)

#       femalepen = all[i, ]$FemalePen[[1]]
#       fpen  = readBin(femalepen, numeric(), n = length(femalepen)/8, size = 8, endian = .Platform$endian)
#       fpens = sprintf("%.4f", fpen)
#       cat(all[i, ]$LocusName, all[i, ]$class_link+1, fpens, file=file, append=TRUE, sep="\t")
#       cat("\tfemale\n", file=file, append=TRUE)

        autopen   = all[i, ]$AutoPen[[1]]
        apen  = readBin(autopen, numeric(), n = length(autopen)/8, size = 8, endian = .Platform$endian)
        apens = sprintf("%.4f", apen)
        cat(all[i, ]$LocusName, all[i, ]$class_link+1, apens, file=file, append=TRUE, sep="\t")
        cat("\tautosomal\n", file=file, append=TRUE)
    }
}

#' generate required PLINK (.phe) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
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
#' mkVCFphe(prefix, ENV, NULL)
#'}
mkVCFphe = function (prefix, ENV, markers) {
    file = paste0(prefix, ".phe")

    unlink(file)

    phenotype_table = ENV$phenotype_table
    hdr = 'FID\tIID'
    
# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = ENV$fam[3:4]

    for (i in 1:ENV$PhenoCnt) {
        hdr = paste0(hdr, "\t", ENV$locus_table[i, 2]) # 2 == LocusName

        off = (i-1) * 8

        raw = unlist(ENV$phenotype_table[,4])
        raw = matrix(raw, ncol=8, byrow=T)

        nrows = nrow(raw)
        if (ENV$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            col = vector("integer", nrows)
            for (j in 1:nrows) {
                col[j] = readBin(raw[j, (off+1):(off+4)], integer(), n=1, size=4)
            }
            col[col==0] = -9
            out$col = col
        } else if (ENV$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            col = vector("numeric", nrows)
            for (j in 1:nrows) {
                col[j] = readBin(raw[j, (off+1):(off+8)], numeric(), n=1, size=8)
            }
            out$col = col
        }
    }
    out$SAMPLEID = paste(ENV$fam[,3], ENV$fam[,4], sep="_")
    
    cat(hdr, "\tSAMPLEID\n", file=file, append=TRUE)

    write.table(out, file=file, sep="\t", quote=FALSE, append=TRUE,
                row.names=FALSE, col.names=FALSE)
}
