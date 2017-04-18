
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
Mega2VCF = function(prefix, markers=NULL, mapno = 0) {

    unlink(paste0(prefix, ".vcf"))

    ENV = getENV()

    if (is.null(markers)) markers = ENV$markers

    mkVCFhdr(prefix, ENV, markers)

    allele_table = ENV$allele_table[ENV$allele_table$locus_link %in% markers$locus_link,]
    map_table = ENV$map_table[ENV$map_table$marker %in% markers$locus_link,]
#   M = ENV$LocusCnt - ENV$PhenoCnt;
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
        GPos = map_table[map_table$map==mapno, c("position", "pos_female", "pos_male")][N, ]
        GPosPos = sprintf("%.2f", GPos$position)
        GPosFem = rep(".", L)
        GPosFem[GPos$pos_female != -99.99] = sprintf("%f", GPos$pos_female)
        GPosMal = rep(".", L)
        GPosMal[GPos$pos_male   != -99.99] = sprintf("%f", GPos$pos_male)

##newer
##new
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
                                 ifelse(flip, ",", ""), ";"),
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
        if (doFlip) {
#             for (xx in whichFlip) {
#                 a1[xx,] = match(a1[xx, ], c(2, 1), nomatch=0)
#                 a2[xx,] = match(a2[xx, ], c(2, 1), nomatch=0)
#             }
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
##
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

browser()
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

    j = 0
    for (i in sapply(split(markers, markers$chromosome), function(x) max(x$position)+1)) {
        j = j + 1
        cat('##contig=<ID=', j, ',length=', i, ',assembly=B37>\n', file=file, append=TRUE, sep="")
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

   #######################################
   ## 100    4524      50      51 2 1   ##
   ## 100    4525      50      51 2 1   ##
   ## 100    4526      50      51 2 -9  ##
   ## 100    4527      50      51 2 2   ##
   ## 100    4528      50      51 1 2   ##
   ## 100      50      52      53 1 -9  ##
   ## 100      51       0       0 2 -9  ##
   ## 100      52       0       0 1 -9  ##
   ## 100      53       0       0 2 -9  ##
   ## 100      54      52      53 1 -9  ##
   ## 100      55       0       0 2 -9  ##
   ## 100    7648      54      55 2 2   ##
   #######################################

    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)

    if (is.null(markers)) markers = ENV$markers
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

    ##############################################
##     Name	Allele	Frequency           ##
## default	1	0.5000              ##
## default	2	0.5000              ##
## exm1517553	1	0.9996              ##
## exm1517553	2	0.0004              ##
## exm1517555	1	0.9982              ##
## exm1517555	2	0.0018              ##
## exm1517564	1	0.9996              ##
## exm1517564	2	0.0004              ##
##############################################


    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)

    if (is.null(markers)) markers = ENV$markers
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

#############################################################
## Chromosome	Name	Map.k.a	BP.p	                   ##
## 20	exm1517553	  0.000000	68363	           ##
## 20	exm1517555	  0.000000	68396	           ##
## 20	exm1517564	  0.000000	76771	           ##
## 20	exm1517584	  0.000000	126149	           ##
## 20	exm1517590	  0.000000	126214	           ##
#############################################################

    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)

    if (is.null(markers)) markers = ENV$markers
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

##############################################
## Name	Class	Pen.11	Pen.12	Pen.22	Type
## default	1	0.0500	0.9000	0.9000	autosomal
##############################################


    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)

    if (is.null(markers)) markers = ENV$markers
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

##############################################
## FID	IID	default	SAMPLEID           ##
## 100	4524	1 	100_4524           ##
## 100	4525	1 	100_4525           ##
## 100	4526	-9 	100_4526           ##
## 100	4527	2 	100_4527           ##
## 100	4528	2 	100_4528           ##
## 100	50	-9 	100_50             ##
## 100	51	-9 	100_51             ##
## 100	52	-9 	100_52             ##
## 100	53	-9 	100_53             ##
## 100	54	-9 	100_54             ##
## 100	55	-9 	100_55             ##
## 100	7648	2 	100_7648           ##
#############################################


    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)

    if (is.null(markers)) markers = ENV$markers
}
