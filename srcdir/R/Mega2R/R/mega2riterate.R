
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

## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/TxDb.Hsapiens.UCSC.hg19.knownGene/extdata/TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")


## sqlite3 /Library/Frameworks/R.framework/Versions/3.2/Resources/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite
## source("https://bioconductor.org/biocLite.R")
## biocLite("org.Hs.eg.db")


#' assemble pedigree information into a data frame
#'
#' @description
#' Generate a data frame with a row for each person.  The observations are:
#' \tabular{ll}{
#' \strong{pedigree} \tab family pedigree name\cr
#' \strong{person} \tab person name\cr
#' \strong{father} \tab father of person\cr
#' \strong{mother} \tab mother of person\cr
#' \strong{sex} \tab sex of person\cr
#' \strong{trait} \tab value of case/control phenotype for person
#' }
#'
#' @usage
#' mkfam(brkloop = FALSE, traitname = "default", envir = ENV)
#'
#' @param brkloop I haven't needed to set this TRUE yet.  Maybe never will.
#'  If loops are broken, a person will be replaced by a dopple ganger in the same family with a
#'  different father/mother.  The number of persons per family will be different when there
#'  are broken loops.  Also, the person_link numbers will be different for all the persons after
#'  the first loop is broken.
#'
#' @param traitname name for trait to use as case/control value; by default, "default"
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @export
#'
#' @return data frame
#' that is described above
#'
## @details
## The observations are: \describe{
## \item{pedigree}{family pedigree name}
## \item{person}{person name}
## \item{father}{father of person}
## \item{mother}{mother of person}
## \item{sex}{sex of person}
## \item{trait}{value of case/control trait for person}
## }
#'
#' @note
#' The columns of this data frame come by selecting the values after merging the data frames: \emph{pedigree_table}, \emph{person_table}, and \emph{trait_table}.
#'
#' @note
#' Also, the father and mother columns from \emph{person_table} are translated from the row index in the
#' \emph{person_table} to the corresponding name.
#'
#' @note
#'  This function stores the data frame in the 'environment' and also returns it.
#'  The function \code{setfam()} stores the data frame
#'  into the 'environment' and adjusts the \emph{genotype_table} and the \emph{phenotype_table}.
#'
#' @examples
#'\dontrun{
#' fam = mkfam()
#'}
mkfam = function (brkloop = FALSE, traitname = "default", envir = ENV) {

    if (brkloop) {
        ped = envir$pedigree_brkloop_table
        per = envir$person_brkloop_table
    } else {
        ped = envir$pedigree_table
        per = envir$person_table
    }
    dofam = function(per) {
        per$Father=per[match(per$Father, per$ID), "PerPre"]
        per[is.na(per$Father), "Father"] = 0
        per$Mother=per[match(per$Mother, per$ID), "PerPre"]
        per[is.na(per$Mother), "Mother"] = 0
        per
    }
    per=unsplit(lapply(split(per, per$pedigree_link), dofam), per$pedigree_link)
    perplus = merge(ped[ , c("pedigree_link", "PedPre")],
                    per[ , c("pedigree_link", "person_link", "PerPre", "Father", "Mother", "Sex")],
                    by = c("pedigree_link"))


    trait = envir$phenotype_table[ , c("person_link", "data")]
    cc = envir$locus_table[envir$locus_table$LocusName == traitname, "locus_link"]
    if (length(cc) == 0) cc = 0
    trait$trait = sapply(trait$data,
                         function (x) { readBin(x[(8*cc+1):8*(cc+1)], integer(), 2, size = 4) })[1, ]

    envir$fam = merge(perplus, trait[ , c("person_link", "trait")], by = "person_link")
}


#' replace the pedigree data frame
#'
#' @description
#' You should first modify the \emph{fam} data frame to filter the members you need to remove.
#'  (For example, you might  want to delete members that have an unknown case/control status.)
#'  This function takes a new data frame of pedigree information and replaces the \emph{fam}
#'  data frame in the 'environment' with it.  Additionally,
#'  changing \emph{fam} data frame will filter the genotypes data frame to only contain persons
#'  matching those in the \emph{fam} data frame.  \code{setfam} also filters for the phenotype data
#'  records.
#'
#' @param fam data frame of family information filtered from \emph{fam}
#' data frame (generated by \code{mkfam}).
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#'\dontrun{
#' setfam(fam)
#'}
setfam = function (fam, envir = ENV) {
    envir$fam = fam

    envir$unified_genotype_table = envir$unified_genotype_table[envir$unified_genotype_table$person_link %in% fam[ , 1], ]
    row.names(envir$unified_genotype_table) = NULL

    envir$phenotype_table = envir$phenotype_table[(envir$phenotype_table$person_link %in% fam[ ,1]), ]
    row.names(envir$phenotype_table) = NULL
    
}

#' load Mega2 database and initialize family data frame and markers data frame
#'
#' @description
#'  Call \code{dbmega2_import()} with the specified database and create an 'environment', with the
#'  SQLite table data loaded into data frames.
#'  Also run \code{mkfam()} to create the pedigree data frame \emph{fam} and then store it with \code{setfam()}.
#'  \code{setfam()} modifies the \emph{unified_genotype_table} (and \emph{phenotype_table}) to match the family members
#'  that remain.
#' @note
#' By default, \code{mkfam} will remove one of each
#'  person that was replicated to break loops in the pedigree, see \code{mkfam} for details.
#'  If you want to leave loops broken, the code is available, but you will have to write your own
#'  version of read.Mega2DB with a different invokation of \code{mkfam()}.
#'
#' @param db specify SQLite database to load
#'
#' @param ... additional arguments to pass to \code{dbmega2_import}
#'
#' @return an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @export
#'
#' @examples
#'\dontrun{
#' ENV = read.Mega2DB("database.db")
#'}
read.Mega2DB = function(db, ...) {

    envir = dbmega2_import(db, ...)

    setfam(mkfam(envir = envir), envir = envir)

    return (envir)
}

#' apply a function to the genotypes (markers) in each gene transcript and/or base pair range
#'
#' @description
#'  This function generates base pair ranges from its input arguments.
#'  Each range specifies a chromosome, a start
#'  base pair and end base pair.  Typically, a range could be a gene transcript, though
#'  it could be a whole chromosome, or a run of base pairs on a chromosome.  Once the
#'  ranges are generated, \code{applyFnToRanges} is called to find all the
#'  rows (i.e. markers) from the \emph{markers} data frame that fall in each range.  For these
#'  markers, a matrix of the genotypes is generated.  Finally, the \code{op} function is called for
#'  each range with the genotypes matrix, markers, range, and 'environment'.
#'
#' @usage
#' applyFnToGenes(op           = function (geno, markers, range, envir) {},
#'                genes_arg    = c("*"),
#'                ranges_arg   = matrix(ncol = 3, nrow = 0),
#'                chrs_arg     = vector("integer", 0),
#'                markers_arg  = vector("character", 0),
#'                type_arg     = "TX",
#'                fuzz_arg     = 0,
#'                envir        = ENV)
#'
#' @param op Is a function of four arguments.  It will be called repeatedly by
#' \code{applyFnToGenes} in a try/catch context.  The arguments are:
#' \describe{
#' \item{geno}{A matrix with one row per \emph{fam} pedigree member and one column for each marker that is selected.  Each datum is two characters indicating the genotype for the marker/member.}
#' \item{markers}{Marker data for each marker selected.  A marker is a data frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding marker genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#' }
#' \item{range}{An indicator of which range argument these markers correspond to.}
#' \item{envir}{An 'environment' holding Mega2R data frames and state data.}
#' }
#'
#' @param genes_arg a character vector of gene names.
#'  All the transcripts identified with the specified gene in BioConductor Annotation,\cr
#'  \bold{TxDb.Hsapiens.UCSC.hg19.knownGene}, are selected.  This produces multiple "range"
#'  elements containing chromosome, start base pair, end base pair.  (If the gene name is "*",
#'  all the transcript will be selected.) Note: BioCoductor Annotation
#'  \bold{org.Hs.eg.db} is used to convert from gene name to ENTREZ gene id.
#'
#' @param ranges_arg an integer matrix of three columns.  The columns define a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param chrs_arg an integer vector of chromosome numbers.  All of the base pairs on each
#'  chromosomes will be selected as a single range.
#'
#' @param markers_arg a data frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding marker genotype data in the\cr
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#' }
#'
#' @param type_arg a character vector of length 1 that contains \bold{"TX"} or does not.  If it is
#'  \bold{"TX"}, which is the default, the \bold{TX} fields of BioConductor Annotation,\cr
#'  \bold{TxDb.Hsapiens.UCSC.hg19.knownGene} are used to define the base pair ranges and chromosome.
#'  Otherwise, the \bold{CDS} fields are used.
#'
#' @param fuzz_arg is an integer vector of length one or two.  The first argument is used to reduce
#'  the start base pair selected from each transcript and the second to increase the end base pair
#'  position.  (If only one value is present, it is used for both adjustments.)  Note: The values
#'  can be positive or negative.
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#'
#' @note
#'  If you want subsequent calls to \code{op} to share information, data can be placed in
#'  a data frame that is added to the 'environment'.
#'
#' @importMethodsFrom GenomeInfoDb 'seqlevels<-'
#' @importFrom AnnotationDbi select keys
#' @export
#'
#' @examples
#'\dontrun{
#' #    show = function(g, m, r, e) {
#' #        print(r)
#' #        print(m)
#' #        print(head(g))
#' #    }
#' #    ENV = read.Mega2DB()
#'
#'    # apply function "show" to all transcripts on genes ELL2 and CARD15
#'    applyFnToGenes(show, genes_arg = c("ELL2", "CARD15"))
#'
#'    # apply function "show" to all genotypes on chromosomes 11 for two base
#'    # pair ranges
#'    applyFnToGenes(show, ranges_arg = matrix(c(11, 50000000, 50100000,
#'                       11, 60000000, 60100000), ncol = 3, nrow = 2, byrow = T))
#'
#'    # apply function "show" to all genotypes for first marker in each chromosome
#'    applyFnToGenes(show, markers_arg = ENV$markers[! duplicated(ENV$markers$chromosome), 3])
#'
#'    # apply function "show" to all genotypes on chromosomes 24 and 26
#'    applyFnToGenes(show, chrs_arg=c(24, 26))
#'
#'}
applyFnToGenes = function (op = function (geno, markers, range, envir) {},
                           genes_arg = c("*"),
                           ranges_arg  = matrix(ncol = 3, nrow = 0),
                           chrs_arg  = vector("integer", 0),
                           markers_arg  = vector("character", 0),
                           type_arg = "TX",
                           fuzz_arg = 0,
                           envir = ENV) {

    env1=loadNamespace(envir$txdb)
    txdb = get(envir$txdb, env1)
    if (type_arg == "TX")
        COLS = c("GENEID", "TXNAME", "TXID", "TXSTRAND", "TXCHROM", "TXSTART", "TXEND")
    else
        COLS = c("EXONNAME", "EXONID", "EXONSTRAND", "EXONCHROM", "EXONSTART", "EXONEND")

    seqlevels(txdb) = paste("chr", c(1:22, "X", "Y", "M"), sep="")

    env2=loadNamespace(envir$entrezGene)
    genedb = get(envir$entrezGene, env2)

    if (genes_arg[1] == "*") {
        entrez = select(genedb, keys = keys(genedb, keytype="ENTREZID"),
                        columns = c("ALIAS", "ENTREZID", "SYMBOL"), keytype = "ENTREZID")
        entrez = entrez[ !duplicated(entrez[ , 1]), ]

        pb = select(txdb, keys = keys(txdb, keytype="TXID"),
                    columns = COLS, keytype = "TXID")
        pb = pb[!duplicated(pb[ , c(4,6,7)]), ]
    } else {
        entrez = select(genedb, keys = genes_arg,
                        columns = c("ALIAS", "ENTREZID", "SYMBOL"), keytype = "ALIAS" )
        entrez = entrez[ !duplicated(entrez[ , 2]), ]

        pb = select(txdb, keys = entrez[ , 2], columns = COLS, keytype = "GENEID")
        pb = pb[!duplicated(pb[ , c(4,6,7)]), ]
    }

    chr2int = data.frame(chr = c(1:26, 23:26))
    chr2int$txchrom = paste("chr", chr2int$chr, sep = "")
    chr2int[27:30,2] = c("chrX", "chrY", "chrXY", "chrM")
##  deleted chrUn
##  allow stuff after chr#_
##  xx = function(l) { l[1] }
##  pb$TXCHROM = chr2int$chr[match(sapply(strsplit(pb$TXCHROM, "_"), xx), chr2int$txchrom)]
##  else don't
    pb$TXCHROM = chr2int$chr[match(pb$TXCHROM, chr2int$txchrom)]
    pb = pb[! is.na(pb$TXCHROM), ]

    ranges = merge(entrez, pb, by.x = "ENTREZID", by.y = "GENEID", all.y = TRUE)
    if (genes_arg[1] == "*") {
        ranges[is.na(ranges[,2]),2] = ranges[is.na(ranges[,2]),5]
        ranges[is.na(ranges[,3]),3] = ranges[is.na(ranges[,3]),5]
    }

    if (length(fuzz_arg) == 2) {
        if (fuzz_arg[1] != 0 | fuzz_arg[2] != 0) {
            ranges[ , 8] = ranges[ , 8] - fuzz_arg[1]
            ranges[ , 9] = ranges[ , 9] + fuzz_arg[2]
        }
    } else if (fuzz_arg[1] != 0) {
        ranges[ , 8] = ranges[ , 8] - fuzz_arg[1]
        ranges[ , 9] = ranges[ , 9] + fuzz_arg[1]
    }

    #chrs
    for (chr in chrs_arg) { ranges = rbind(ranges,
#           ENTREZID  ALIAS SYMBOL  TXID     TXNAME TXCHROM TXSTRAND  TXSTART    TXEND
                                    list("-", "-", "-", "-", "-",
                                         chr, "-", 0, 1000000000) )
                      }
    #ranges
    if (length(ranges_arg ))
    for (i in 1:dim(ranges_arg)[1]) { ranges = rbind(ranges,
#           ENTREZID  ALIAS SYMBOL  TXID     TXNAME TXCHROM TXSTRAND  TXSTART    TXEND
                                    list("-", "-", "-", "-", "-",
                                         ranges_arg[i, 1], "-", ranges_arg[i, 2], ranges_arg[i, 3]) )
                      }

    applyFnToRanges(op, ranges, c(6, 8, 9), envir = envir)

    #marks
    if (length(markers_arg)) {
        applyFnToMarkers(op, envir$markers[envir$markers$MarkerName %in% markers_arg, ], envir = envir)
    }
}

#' set default range data: chromosome and start/end base pair
#'
#' @description
#'  This function sets the default list of ranges used by \code{applyFnToRanges}.  \code{applyFnToRanges}
#'  examines each range and the set of markers that fall within the range will be
#'  processed.
#'
#' @param ranges a data frame that contains at least 4 observations: a name, a chromosome, a start
#' base pair position and an end base pair position.
#'
#' @param indices a vector of 3 integers that specify the location of chromosome column, start base pair
#'  column and end base pair column of range data frame.
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' setRanges(range, c(3, 4, 5))
#'}
setRanges = function (ranges, indices, envir = ENV) {
    envir$refRanges  = ranges
    envir$refIndices = indices
}

#' set default name of transcription database and name of database mapping gene name to entrez gene id
#'
#' @description
#'  This function takes two string parameters:
#'   one to specify entrez gene ids to transcripts,
#'   the other to map gene names to entrez gene id's.
#'
#' @param txdb name of Bioconductor transcription database.
#'
#' @param entrezGene name of Bioconductor mapping of gene name or gene alias to entrez gene id

#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#'
#' @export
#'
#' @note
#'  \bold{Mega2R} will take care to load the necessary databases, but you will have to
#'  install them from Bioconductor.  This is explained at length in the package Vignette.
#'
#' @examples
#'\dontrun{
#' setAnnotations("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db")
#'}
setAnnotations = function (txdb, entrezGene, envir = ENV) {
    envir$txdb = txdb
    envir$entrezGene = entrezGene
}


#' apply a function to all the genotypes of markers in each of several specified ranges
#'
#' For each set of
#'  markers, a matrix of the genotypes is generated.  Finally, the \code{op} function is called for
#'  each range with the genotypes matrix, markers, range, and 'environment'.
#'
#' @usage
#' applyFnToRanges(op          = function (geno, markers, range, envir) {},
#'                 ranges_arg  = NULL,
#'                 indices_arg = NULL,
#'                 fuzz_arg    = 0,
#'                 envir       = ENV)
#' 
#' @param op Is a function of four arguments.  It will be called repeatedly by
#' \code{applyFnToRanges} in a try/catch context.  The arguments are:
#' \describe{
#' \item{geno}{A matrix with one row per \emph{fam} pedigree member and one column for each marker that is selected.  Each datum is two characters indicating the genotype for the marker/member.}
#' \item{markers}{Marker data for each marker in \strong{geno}.  A marker is a data frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding marker genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#' }
#' \item{range}{An indicator of which range argument of \code{applyFnToRanges} these markers correspond to.}
#' \item{envir}{An 'environment' holding Mega2R data frames and state data.}
#' }
#'
#' @param ranges_arg is a data frame that contains at least 4 observations: a name, a chromosome, a 
#'  start base pair position and an end base pair position.
#'
#' @param indices_arg is a vector of 3 integers that specify the location of chromosome, start base
#'  pair column and end base pair column of the ranges_arg data frame.
#'
#' @param fuzz_arg is an integer vector of length one or two.  The first argument is used to reduce
#'  the start base pair selected from each range and the second to increase the end base pair
#'  position.  (If only one value is present, it is used for both changes.) Note: The values
#'  can be positive or negative.
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @note
#'  If the \emph{ranges_arg} and \emph{indices_arg} are NULL or missing, then the default ranges that have been set by \code{setRanges}
#'  are used.  If \code{setRanges} has not been called, a default set of the ranges is used.
#'
#' @export
#'
#' @examples
#'\dontrun{
#' #    show = function(g, m, r, e) {
#' #        print(r)
#' #        print(m)
#' #        print(head(g))
#' #    }
#' #    ENV = read.Mega2DB()
#'
#'    # apply function "show" to all genotypes on chromosomes 11 for two base
#'    # pair ranges
#'    applyFnToRanges(show, 
#'                    ranges_arg =
#'                    matrix(c(11, 50000000, 50100000,
#'                             11, 60000000, 60100000),
#'                            ncol = 3, nrow = 2, byrow = T),
#'                    indices_arg = 1:3)
#'
#'}
applyFnToRanges = function (op          = function (geno, markers, range, envir) {},
                            ranges_arg  = NULL,
                            indices_arg = NULL,
                            fuzz_arg    = 0,
                            envir = ENV) {

    if (is.null(ranges_arg)) {
        ranges  = envir$refRanges
        indices = envir$refIndices
    } else {
        ranges  = ranges_arg
        indices = indices_arg
    }

    rows = nrow(ranges)
    if (rows) {
        start = ranges[ , indices[2]]
        end   = ranges[ , indices[3]]
        if (length(fuzz_arg) == 2) {
            if (fuzz_arg[1] != 0 | fuzz_arg[2] != 0) {
                start = start - fuzz_arg[1]
                end   = end   + fuzz_arg[2]
            }
        } else if (fuzz_arg[1] != 0) {
            start = start - fuzz_arg[1]
            end   = end   + fuzz_arg[1]
        }
        chrm = as.integer(sub("chr", "", ranges[ , indices[1]]))

        Uranges    = data.frame(start=integer(rows), end=integer(rows), chrm=integer(rows), i=integer(rows))
        Umarkersub = vector("list", rows)

        for (i in 1:rows) {

            if (is.na(chrm[i]) || is.na(start[i]) || is.na(end[i]) ) next

            markersub = envir$markers[ envir$markers$chromosome == chrm[i] & envir$markers$position <= (end[i]) & envir$markers$position >= (start[i]), ]

            Umarkersub[[i]] = markersub

            if (nrow(markersub)) {
                markersubpos     = sort(markersub$position)
                Uranges$start[i] = markersubpos[1]
                Uranges$end[i]   = markersubpos[length(markersubpos)]
                Uranges$chrm[i]  = chrm[i]
                Uranges$i[i]     = i
            } else {
                if (envir$verbose)
                    message("tryFn() No markers in range: ", paste(ranges[i, ], collapse=", "))
                next
            }
        }

        Uranges =  Uranges[ ! duplicated(Uranges[,1:3]), ]

        for (i in Uranges$i)
        {
            if (i == 0) next

             markersub = Umarkersub[[i]]
             geno = getgenotypes(markersub, envir = envir)
             tryFn(op, geno, markersub, ranges[i, ], envir)
        }
    }
}

# warning = function( ...) {
#   ...
#    withRestarts({
#        .Internal(.signalCondition(cond, message, call))
#        .Internal(.dfltWarn(message, call))
#    }, muffleWarning = function() NULL)

tryFn1 = function(op, geno, markersub, ranges, envir) {
    withRestarts( #aka    tryCatch
        withCallingHandlers (
            { op(geno, markersub, ranges, envir) },
              error   = function(e) {
                  message("tryFn() <simpleError:: ", conditionMessage(e), ">")
                  invokeRestart("abort")
              },
              warning = function(w) {
                  message("tryFn() <simpleWarning:: ", conditionMessage(w), ">")
                  invokeRestart("muffleWarning")
             }
        ),
      abort = function () {
          message("tryFn() aborting transcript: ", paste(ranges, collapse=", "))
      }
    )
}

tryFn = function(op, geno, markersub, ranges, envir) {
    tryCatch(
        withCallingHandlers (
            { op(geno, markersub, ranges, envir) },
              warning = function(w) {
                  message("tryFn() <simpleWarning:: ", conditionMessage(w), ">")
                  invokeRestart("muffleWarning")
             }
        ),
        error = function(e) {
                  message("tryFn() <simpleError:: ", conditionMessage(e), ">")
                  message("tryFn() aborting transcript: ", paste(ranges, collapse=", "))
              }
    )
}

#' apply a function to the genotypes from a set of markers
#'
#' @description
#'  A matrix of the genotypes for all the specified markers is generated.  Then, the call back function, \code{op},
#'   is called with the genotypes, markers, NULL (for the range), and the 'environment'.
#'
#' @usage
#' applyFnToMarkers(op      = function (geno, markers, range, envir) {},
#'                 markers_arg,
#'                 envir = ENV)
#'
#' @param op Is a function of four arguments.  It will be called once by
#' \code{applyFnToMarkers} in a try/catch context.  The arguments are:
#' \describe{
#' \item{geno}{A matrix with one row per \emph{fam} pedigree member and one column for each marker that is selected.  Each datum is two characters indicating the genotype for the marker/member.}
#' \item{markers}{Marker data for each marker in \strong{geno}.  A marker is a data frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding marker genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#' }
#' \item{range}{NULL: to indicate no explicit range was specified.}
#' \item{envir}{An 'environment' holding Mega2R data frames and state data.}
#' }
#'
#' @param markers_arg a data frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding marker genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#'
#' @param envir an 'environment' that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#'\dontrun{
#' #    show = function(g, m, r, e) {
#' #        print(r)
#' #        print(m)
#' #        print(head(g))
#' #    }
#' #    ENV = read.Mega2DB()
#'
#'    # apply function "show" to all genotypes in chromosome 20, 21, 22, and 23
#'    applyFnToMarkers(show, ENV$markers[ENV$markers$chromosome %IN% 20:23),])
#'
#'}
applyFnToMarkers = function (op = function (geno, markers, range, envir) {},
                             markers_arg,
                             envir = ENV) {

    geno = getgenotypes(markers_arg, envir = envir)
    tryFn(op, geno, markers_arg, NULL, envir)
}
