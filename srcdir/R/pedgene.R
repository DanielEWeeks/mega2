
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

library(pedgene)

init_pedgene = function (db = "ped3.db") {

    dbmega2_import(db)

#   NonmissingPheID=get(load("../NonmissingPheID.RData"))
    non=read.table("ped3.famphe", header=F)
    mkpedigree()
    pl=merge(ped.Y[ , c(1, 3, 4)], non[ , 2:3],
             by.x = c("PedPre", "PerPre"), by.y = c("V2", "V3"))

    ped.Y = ped.Y[ped.Y[ , 1] %in% pl[ , 3], ]
    unified_genotype_table = unified_genotype_table[unified_genotype_table$person_link %in% pl[ , 3], ]
    row.names(ped.Y) = NULL
    row.names(unified_genotype_table) = NULL
    assign("ped.Y", ped.Y, globalenv())
    assign("unified_genotype_table", unified_genotype_table, globalenv())

#   refGene = read.table("../refseq_genes.txt", header= TRUE, stringsAsFactors= FALSE)
#   refGene$cdsStart <-refGene$cdsEnd <- NULL
#   refGene = refGene[nchar(refGene$chrom) <= 5, ]
#   refGene = refGene[refGene$chrom != "chrX" & refGene$chrom != "chrY", ]
#   refGene = refGene[!duplicated(refGene), ]
#   colnames(refGene) = c("XX", "name2", "chrom", "txStart", "txEnd")
#   colnames(refGene) = c("XX", "SYMBOL", "TXCHROM", "TXSTART", "TXEND")
#   row.names(refGene) = NULL
#   write.table(refGene, file=ped3.ref, quote=F, row.names=F)
    refGene=read.table("ped3.ref", header = T)
    refGene = refGene[! duplicated(refGene$SYMBOL), ]
    assign("refGene", refGene, pos = globalenv())
}

run = function (gs = 1:100) {
    zzz = 0
    results <- data.frame(chr = character(0), gene = character(0), nvariants = numeric(0), 
                          start = numeric(0), end = numeric(0), 
                          pKernel_BT = numeric(0), pBurden_BT = numeric(0), 
                          pKernel_MB = numeric(0), pBurden_MB = numeric(0),
                          pKernel_UW = numeric(0), pBurden_UW = numeric(0), 
                          geneID = numeric(0), stringsAsFactors = FALSE)
    assign("zzz", zzz, pos = globalenv())
    assign("results", results, pos = globalenv())
    unlink("k_Schaid_rare.txt")

    applyFnToRanges(DOpedgene, refGene[gs, ], indices = 3:5)
}

DOpedgene = function(genoChar, markerFrame, rng) {

    markerList = markerFrame$MarkerName
    schaidPed = ped.Y[ , c(-1, -2)]
    colnames(schaidPed) = c("ped", "person", "father", "mother", "sex", "trait")
    pedPer = schaidPed[ , 1:2]

    gene  <- as.character(rng$SYMBOL)

    mt = matrix(c(11, 12, 21, 22, 0, 1, 1, 2), nrow = 4, ncol = 2)
    di = dim(genoChar)
    genoInt = matrix(0, nrow = (di[1]), ncol = di[2])
    for (k in 1:(di[2])) {
        vec = mt[match(as.integer(genoChar[ , k]), mt), 2]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        cat(gene, markerList[k], g0, g1, g2, "\n")
        if (g0 < g2) {
           genoInt[ , k] = 2 - vec
        } else {
           genoInt[ , k] =     vec
        }
    }

    genoInt = matrix(genoInt, nrow = di[1])
    maf = colMeans(genoInt)
    pos = markerList[maf > 0]
    if (length(pos) >= 2) {       # at least 2 non-polymorphic variants #    
        genoInt <- genoInt[ , maf > 0]     # remove nonpolymorphic variants #
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
        chr   <- as.character(rng$TXCHROM)
        start <- rng$TXSTART
        end   <- rng$TXEND
        zzz = zzz + 1
        cat(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                          pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz, "\n")
        results = get("results")
        results[1, ] <- c(chr, gene, nsnp, start, end, pKernel_BT, pBurden_BT,
                          pKernel_MB, pBurden_MB, pKernel_UW, pBurden_UW, zzz)
        write.table(results, file="k_Schaid_rare.txt", append= TRUE, row.names= FALSE, col.names= FALSE, quote= FALSE)
    }
}
