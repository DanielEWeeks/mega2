
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

## source("../mega2.R")
## dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db")
## tst1()


## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/TxDb.Hsapiens.UCSC.hg19.knownGene/extdata/TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

library("TxDb.Hsapiens.UCSC.hg19.knownGene")

## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("org.Hs.eg.db")

library("org.Hs.eg.db")

pedgene.ped = NULL

library(pedgene)

init = function () {

#    source("/Users/rbaron/mega2/bb/srcdir/mega2.R")
#    source("/Users/rbaron/mega2/bb/srcdir/R/pedgene.R")
    dbmega2_import("ped3.db")
#    rm("genotype_table", globalenv())

#   NonmissingPheID=get(load("../NonmissingPheID.RData"))
    non=read.table("ped3.famphe", header=F)
    mkped(T)
    pl=merge(ped.Y[,c(1,3,4)], non[, 2:3], by.x=c("PedPre", "PerPre"), by.y=c("V2", "V3"))

    ped.Y = ped.Y[ped.Y[,1] %in% pl[,3], ]
    unified_genotype_table = unified_genotype_table[unified_genotype_table$person_link %in% pl[,3],]
    row.names(ped.Y) = NULL
    row.names(unified_genotype_table) = NULL
    assign("ped.Y", ped.Y, globalenv())
    assign("unified_genotype_table", unified_genotype_table, globalenv())

#   refGene = read.table("../refseq_genes.txt", header= TRUE, stringsAsFactors= FALSE)
#   refGene$cdsStart <-refGene$cdsEnd <- NULL
#   refGene = refGene[nchar(refGene$chrom) <= 5, ]
#   refGene = refGene[refGene$chrom != "chrX" & refGene$chrom != "chrY", ]
#   refGene = refGene[!duplicated(refGene), ]
#   row.names(refGene) = NULL
#   colnames(refGene) = c("XX", "name2", "chrom", "txStart", "txEnd")

    refGene=read.table("ped3.ref",header=F)
    colnames(refGene) = c("XX", "name2", "chrom", "txStart", "txEnd")
    refGene = refGene[! duplicated(refGene$name2), ]
    assign("refGene", refGene, pos=globalenv())

    markers = merge(marker_table[ , c("locus_link","locus_link_fill","MarkerName","chromosome")],
                    map_table[ map_table$map == 1, c( "marker", "position")],
                    by.x="locus_link", by.y="marker")
    assign("markers", markers, pos=globalenv())
}

run = function (gs=1:100) {
    zzz = 0
    results <- data.frame(chr= character(0), gene= character(0), nvariants= numeric(0), 
                          start= numeric(0), end= numeric(0), 
                          pKernel_BT= numeric(0), pBurden_BT= numeric(0), 
                          pKernel_MB= numeric(0), pBurden_MB= numeric(0),
                          pKernel_UW= numeric(0), pBurden_UW= numeric(0), 
                          geneID= numeric(0), stringsAsFactors= FALSE)
    assign("zzz", zzz, pos=globalenv())
    assign("results", results, pos=globalenv())
    unlink("k_Schaid_rare.txt")

    applyFnToRanges(dopedgene, refGene[gs,], indices=3:5)
}

mkped = function (brkloop=T) {

    if (brkloop) {
        ped = pedigree_brkloop_table
        per = person_brkloop_table
    } else {
        ped = pedigree_table
        per = person_table
    }
    ped.X = merge(ped[,c("pedigree_link","PedPre")],
                  per[,c("pedigree_link","person_link","PerPre","Father","Mother","Sex")],
                  by=c("pedigree_link"))

    trait = phenotype_table[ , c("person_link", "data")]
    trait$trait = sapply(trait$data, function (x) { readBin(x, integer(), 2, size=4) })[1,]
    ped.Y = merge(ped.X, trait[ , c("person_link", "trait")], by="person_link")
    assign("ped.Y", ped.Y, pos=globalenv())
    
    pedgene.ped = ped.Y[ , c(-1, -2)]
}

mkmarkers = function (genes=c("ELL2", "CARD15"),
                      type="TX",
                      fuzz=0,
                      ranges=matrix(ncol=3,nrow=0),
                      chrs=vector("integer", 0),
                      marks=vector("character", 0)) {

    ## dbconn(gene)/dbConn(txdb)
    ## dbReadTable(dbconn(), "tbl")

      markers = merge(marker_table[ , c("locus_link","locus_link_fill","MarkerName","chromosome")],
                    map_table[ map_table$map == 1, c( "marker", "position")],
                    by.x="locus_link", by.y="marker")
    assign("markers", markers, pos=globalenv())

    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    if (type=="TX")
        COLS = c("TXNAME", "TXID", "TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
    else
        COLS = c("EXONNAME", "EXONID", "EXONSTRAND", "EXONCHROM", "EXONSTART", "EXONEND")

    seqlevels(txdb) = concat("chr", c(1:22, "X", "Y", "M"))
    
    genedb=org.Hs.eg.db

    pa = select(genedb, keys=genes, columns=c("ALIAS", "ENTREZID", "SYMBOL"), keytype="ALIAS")
    pb = select(txdb, keys = pa[,2], columns=COLS, keytype="GENEID")
    range = merge(pa, pb, by.x="ENTREZID", by.y="GENEID")
    range[,6] = as.integer( sub("chr", "", range[,6]))

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

    #marks
    rows = dim(range)[1]
    if (length(marks)) {
        range = rbind(range,
                        list("-", "-", "-", "-", "-",
                             0, "-", 0, 0) )
        pos = vector("list", rows+1)
        geno = vector("list", rows+1)
    } else {
        pos = vector("list", rows)
        geno = vector("list", rows)
    }
    assign("range", range, pos=globalenv())

    for (i in 1:rows) {
        print(i)
        if (is.na(range[i,6]) || is.na(range[i,8]) || is.na(range[i,9]) ) next
        pos[[i]] = markers[ markers$chromosome == range[i,6] & markers$position <= (range[i,9] + fuzz) & markers$position >= (range[i,8] - fuzz), ]
#       pos[[i]]$locus_link_fill = pos[[i]]$locus_link + chr_gap_skip[pos[[i]]$chromosome]
#       rownames(pos[[i]]) = NULL
        
        geno[[i]] = getlocus(pos[[i]]$locus_link, pos[[i]]$locus_link_fill,
                              int_table[int_table$key == 'PhenoCnt',][1,3])
    }

    pos[[rows+1]] = markers[markers$MarkerName %in% marks, ]
#   pos[[rows+1]]$locus_link_fill = pos[[rows+1]]$locus_link + chr_gap_skip[pos[[rows+1]]$chromosome]
#   rownames(pos[[rows+1]]) = NULL
    assign("pos", pos, pos=globalenv())

    geno[[rows+1]] = getlocus(pos[[rows+1]]$locus_link,
                              pos[[i+1]]$locus_link_fill,
                              int_table[int_table$key == 'PhenoCnt',][1,3])
    assign("geno", geno, pos=globalenv())

}

################

tst2 = function(genes=c("ELL2", "CARD15"), type="TX", fuzz=0) {
    mkped()
    mkmarkers(genes=genes, type=type, fuzz=fuzz)

    for (i in 1:dim(range)[1]) {
        print(range[i,])
        print(pos[[i]])
        print(head(geno[[i]]))
    }
}

tst3 = function(genes=c("ELL2", "CARD15"), type="TX", fuzz=0) {
    mkped()
    mkmarkers(genes=genes, type=type, fuzz=fuzz,
              matrix(c(11,50000000,50100000,11,60000000,60100000),ncol=3,nrow=2,byrow=T),
              marks=markers[! duplicated(markers$chromosome), 3],
              chrs=c(24,26))

    for (i in 1:dim(range)[1]) {
        print(range[i,])
        print(pos[[i]])
        print(head(geno[[i]]))
    }
}

tst31 = function(genes=c("ELL2", "CARD15"), type="TX", fuzz=0) {
    mkped()
    mkmarkers1(genes=genes, type=type, fuzz=fuzz,
              matrix(c(11,50000000,50100000,11,60000000,60100000),ncol=3,nrow=2,byrow=T),
              marks=markers[! duplicated(markers$chromosome), 3],
              chrs=c(24,26))
}

tst4 = function(genes=c("ELL2"), type="TX", fuzz=0) {
    mkped()
    mkmarkers(genes=genes, type=type, fuzz=fuzz,
              marks=c("rs6587762", "rs7521920",
                      "rs10181821", "rs10195681", "rs7594567", "rs4637157") )

    for (i in 1:dim(range)[1]) {
        print(range[i,])
        print(pos[[i]])
        print(head(geno[[i]]))
    }
}

###############

applyFnToGenes = function (genes=c("ELL2", "CARD15"),
                           type="TX",
                           fuzz=0,
                           ranges=matrix(ncol=3,nrow=0),
                           chrs=vector("integer", 0),
                           marks=vector("character", 0)) {

    show = function(g, m, r) {
        print(r)
        print(m)
        print(head(g))
    }
    
    ## dbconn(gene)/dbConn(txdb)
    ## dbReadTable(dbconn(), "tbl")

    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    if (type=="TX")
        COLS = c("TXNAME", "TXID", "TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
    else
        COLS = c("EXONNAME", "EXONID", "EXONSTRAND", "EXONCHROM", "EXONSTART", "EXONEND")

    seqlevels(txdb) = concat("chr", c(1:22, "X", "Y", "M"))
    
    genedb=org.Hs.eg.db

    pa = select(genedb, keys=genes, columns=c("ALIAS", "ENTREZID", "SYMBOL"), keytype="ALIAS")
    pb = select(txdb, keys = pa[,2], columns=COLS, keytype="GENEID")
    range = merge(pa, pb, by.x="ENTREZID", by.y="GENEID")
    range[,6] = as.integer( sub("chr", "", range[,6]))

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

    applyFnToRanges(show, range, c(6,8,9))

    #marks
    if (length(marks)) {
        positions = markers[markers$MarkerName %in% marks, ]
        applyFnToRanges(show, pos=positions)
    }

}

################

zzz = 0

applyFnToRanges = function (op=function (arg, mks) {},
                       range = matrix(ncol=3,nrow=0),
                       indices=1:3,
                       pos=NULL,
                       fuzz=0  ) {
    pheno = int_table[int_table$key == 'PhenoCnt',][1,3]
    rows = nrow(range)
    if (rows) {

        start = range[ , indices[2]]
        end   = range[ , indices[3]]
        chrm = as.integer(sub("chr", "", range[ , indices[1]]))

        for (i in 1:rows) {
            print(i)
            if (is.na(chrm[i]) || is.na(start[i]) || is.na(end[i]) ) next

            posx = markers[ markers$chromosome == chrm[i] & markers$position <= (end[i] + fuzz) & markers$position >= (start[i] - fuzz), ]

            if (nrow(posx)) {
                geno = getlocus(posx$locus_link, posx$locus_link_fill, pheno)
                op(geno, posx, range[i,])
            } else {
                cat("No markers in range:  chr", chrm[i], " between ", start[i], " and ", end[i], "\n")
            }
        }
    }

    if (! is.null(pos)) {
        geno = getlocus(pos$locus_link, pos$locus_link_fill, pheno)
        op(geno, pos, NULL)
    }

}

###############

# pedp = function(g, m) {cat(dim(g), m, "\n")}

dopedgene = function(genoChar, markerFrame, rng) {

    markerList = markerFrame$MarkerName
    schaidPed = ped.Y[ , c(-1, -2)]
    colnames(schaidPed) = c("ped", "person", "father", "mother", "sex", "trait")
    pedPer = schaidPed[, 1:2]

    mt = matrix(c(11, 12, 21, 22, 0, 1, 1, 2), nrow=4,ncol=2)
    di = dim(genoChar)
    genoInt = matrix(0, nrow=(di[1]), ncol=di[2])
    for (k in 1:(di[2])) {
        vec = mt[match(as.integer(genoChar[,k]), mt), 2]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        cat(g0, g1, g2, "\n")
        if (g0 < g2) {
           genoInt[, k] = 2 - vec
        } else {
           genoInt[, k] =     vec
        }
    }

    genoInt = matrix(genoInt, nrow=di[1])
    maf = colMeans(genoInt)
    pos = markerList[maf > 0]
    if (length(pos) >= 2) {       # at least 2 non-polymorphic variants #    
        genoInt <- genoInt[ ,maf > 0]     # remove nonpolymorphic variants #
        nsnp    <- ncol(genoInt)
        weight <- rep(1, ncol(genoInt))

        pedgeno <- cbind(pedPer, genoInt)

        BT <- pedgene(schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= FALSE, method= "kounen") 
        pKernel_BT <- BT$pgdf$pval.kernel
        pBurden_BT <- BT$pgdf$pval.burden

        MB <- pedgene(schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= TRUE, method= "kounen") 
        pKernel_MB <- MB$pgdf$pval.kernel
        pBurden_MB <- MB$pgdf$pval.burden

        UW <- pedgene(schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= weight, weights.mb= TRUE, method= "kounen", acc.davies=1e-9) 
        pKernel_UW <- UW$pgdf$pval.kernel
        pBurden_UW <- UW$pgdf$pval.burden 

        ## read out the results ##
        chr   <- as.character(rng$chrom)
        gene  <- as.character(rng$name2)
        start <- rng$txStart
        end   <- rng$txEnd

        zzz = zzz + 1
        results = get("results")
        results[1, ] <- c(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                          pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz)
        cat(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                          pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz, "\n")
        write.table(results, file="k_Schaid_rare.txt", append= TRUE, row.names= FALSE, col.names= FALSE, quote= FALSE)
    }
}
