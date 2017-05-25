
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

#' Mega2 package
#'
#' @description This package reads a Mega2 SQLite3 database into R dataframes and
#'	makes the contained genotypes/phenotypes/linkage data available for analysis.
#'
#' @author Robert V Baron
#' @docType package
#' @name Mega2-package
NULL

#library(DBI)
#library(RSQLite)

ENV = new.env(parent = emptyenv())

#' Mega2 SQLite3 tables
#'
#' @description This indicates which SQLite3 tables and possibly which subset of fields
#'	to fetch.
#'
#' @author Robert V Baron
#' @docType data
#' @name Mega2-TBLS
TBLS = c("int_table",
#        "double_table",
#        "charstar_table",
#        "stuff_table",
#        "batch_parameters",
#        "file_table",

         "pedigree_table",
         "person_table",
         "pedigree_brkloop_table",
         "person_brkloop_table",

#        "canonicalallele_table",
         "locus_table",
         "allele_table",
         "marker_table",
         "markerscheme_table",
         "map_table",
         "mapnames_table",

         "traitaff_table",
         "affectclass_table",
#        "ptraitquant_table",

         "phenotype_table",
         "genotype_table"

## Derived tables
##        unified_genotype_table
##        markers
  )

#' Mega2 SQLite3 table filter
#'
#' @description For the for mentioned tables in some cases we will only extract a subset of columns
#'
#' @author Robert V Baron
#' @docType data
#' @name Mega2-TBLSFilter
#' @note For the data.frames below, only the specified fields are loaded from the SQLite tables
TBLSFilter = list(
          locus_table    = "pId, LocusName, Type, AlleleCnt, locus_link",

          marker_table   = "pId, MarkerName, pos_avg, pos_female, pos_male, chromosome, locus_link",
          pedigree_table = "pId, Num, EntryCnt, Name, PedPre, OriginalID, origped, pedigree_link",
          person_table   = "pId, UniqueID, OrigID, FamName, PerPre, ID, Father, Mother, Sex, pedigree_link, person_link",
          pedigree_brkloop_table = "pId, Num, EntryCnt, Name, PedPre, OriginalID, origped, pedigree_link",
          person_brkloop_table   = "pId, UniqueID, OrigID, FamName, PerPre, ID, Father, Mother, Sex, pedigree_link, person_link",

          traitaff_table = "pId, ClassCnt, PenCnt, locus_link"
  )

#' make markers data frame
#'
#' @description Create the markers data frame as a subset of the markers_table.  It contains 5
#'	observations: \describe{
#'	  \item{locus_link:}{locus offset of this marker}
#'	  \item{locus_link_fill:}{locus offset plus an accumulating fudge factor that jumps
#'          with each new chromosome so that the count of markers per chromosome is force to be
#'          a multiple of 4.  (This number corresponds to the offset of this marker in the
#'          \code{unified_genotype_table}.)}
#'	  \item{MarkerName:}{name of this marker}
#'	  \item{chromosome:}{chromosome number of this marker}
#'	  \item{position:}{base pair position of this marker (selected by bpPosMap[below])}
#'       }
#'
#' @param bpPosMap An integer that indicates the map index to use when selecting the
#'	chromosome/position fields from the map_table to merge with the marker_table.
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mk_markers_with_skip(1)
#'}
mk_markers_with_skip = function(bpPosMap = 1) {
    if (ENV$MARKER_SCHEME == 1) {
        markersPerChr = sapply(split(ENV$marker_table$chromosome, ENV$marker_table$chromosome), length)
        idx = as.integer(names(markersPerChr))
        m = rep(0, max(idx))
        m[idx] = markersPerChr
        markersPerChr = m
        extra_markers = cumsum(4 * floor((markersPerChr + 3) / 4) - markersPerChr)
        extra_markers = c(0, extra_markers)
        names(extra_markers) = NULL
    } else if (ENV$MARKER_SCHEME == 2) {
        extra_markers = vector("integer", length(unique(ENV$marker_table$chromosome))+1)
    }
    ENV$marker_table$locus_link_fill = ENV$marker_table$locus_link + extra_markers[ENV$marker_table$chromosome]
    if (any(is.na(ENV$marker_table$locus_link_fill))) {
        stop("Internal Error in mk_markers_with_skip.  Bad fill values.", call. = FALSE)
    }

    map_table = ENV$map_table[ ENV$map_table$map == bpPosMap, c( "marker", "position")]
    if (nrow(map_table) == 0) {
        message("No entry for map == ", bpPosMap, " in map_table.  Using map == 0 instead.")
        map_table = ENV$map_table[ ENV$map_table$map == 0, c( "marker", "position")]
    }
    ENV$markers = merge(ENV$marker_table[ , c("locus_link", "locus_link_fill", "MarkerName", "chromosome")],
                        map_table,
                        by.x = "locus_link", by.y = "marker")
}

#' mk_unified_genotype_table aggregate separate genotype vectors for each chromosome to one extended vector containing all the chromosomes.
#'
#' @description The genotype_table contains for each person a separate record for each chromosome.
#'  We need a single vector for each person.
#'
#' @return None
#'
#' @keywords internal
#'
mk_unified_genotype_table = function() {
    samples = split(ENV$genotype_table, ENV$genotype_table$person_link)
    samplesize = length(samples)
    person_link = unique(ENV$genotype_table$person_link)

    df = data.frame(row.names = 1:samplesize,
                    person_link = person_link,
                    data = vector("raw", samplesize))

    for (i in 1:samplesize) {

        chrOrder = order(samples[[i]][ , 3], decreasing = FALSE)
        v = unlist(samples[[i]][chrOrder, 5])
        df$data[i] = list(v)
    }
  
    ENV$unified_genotype_table = df
    ENV$genotype_table = NULL
    rm(list = "genotype_table", envir = ENV)
}

#' dbmega2_import read Mega2 SQLite tables into R
#'
#' @description read fields of tables necessary later tables processing.
#'
#' @usage
#' dbmega2_import(dbname,
#'                bpPosMap = 1,
#'                verbose = FALSE)
#'
#' @param dbname file path to SQLite database.
#'
#' @param bpPosMap specified which map in the map_table should be used for marker chromosome/position.
#'
#' @param verbose For this function, print out statistics on the name/size of each table read and show column headers.
#'
#' @return ENV an environment that contains all the tables created from the SQLite tables.
#'
#' @importFrom RSQLite dbConnect dbExistsTable dbReadTable dbListFields SQLITE_RO
#' @export
#'
#' @examples
#'\dontrun{
#' dbmega2_import(verbose = TRUE)
#'
#' dbmega2_import("foo.db", verbose = TRUE)
#'}
dbmega2_import = function(dbname,
                          bpPosMap = 1,
                          verbose = FALSE) {
    con = tryCatch(dbConnect(RSQLite::SQLite(), dbname = dbname, flags = SQLITE_RO),
                   error = function(xx) { stop("DB open failed: ", dbname, call. = FALSE) })

    gc(verbose = FALSE)

    ENV$verbose = verbose

    for (tbl in TBLS) {
        if (dbExistsTable(con, tbl)) {
            filter = TBLSFilter[tbl][[1]]
            if (is.null(filter))
                assign(tbl, dbReadTable(con, tbl), pos = ENV)
            else
                assign(tbl, dbReadTable(con, tbl, select.cols = filter), pos = ENV)
            if (ENV$verbose) {
                cat(tbl, dim(get(tbl, pos=ENV)), sep = "\t", end = "\n")
                if (is.null(filter))
                    cat(tbl, dbListFields(con, tbl), sep = "\t", end = "\n")
                else
                  cat(tbl, strsplit(filter, ", ")[[1]], sep = "\t", end = "\n")
                cat(end = "\n")
            }
        }
    }

    ENV$PhenoCnt      = ENV$int_table[ENV$int_table$key == 'PhenoCnt', 3]
    ENV$LocusCnt      = ENV$int_table[ENV$int_table$key == 'LocusCnt', 3]
    ENV$MARKER_SCHEME = ENV$int_table[ENV$int_table$key == 'MARKER_SCHEME', 3]
    if (ENV$MARKER_SCHEME > 2) {
        stop("Only compressions levels of 1 or 2 are allowed. (",
             ENV$MARKER_SCHEME, ")", call. = FALSE)
    }
    
    mk_unified_genotype_table()
    mk_markers_with_skip(bpPosMap)
    if (ENV$MARKER_SCHEME == 2) {
        message("Partitioninging allele_table by locus_link\n");
        ENV$locus_allele_table = split(ENV$allele_table, ENV$allele_table$locus_link)
    } else {
        ENV$locus_allele_table = NULL
    }

    gc(verbose = FALSE)
    return (ENV)
}

#' return ENV environment from package
#'
#' Mega2 uses an environment, \emph{ENV} to store all the tables it reads in.  In addition,
#'  some package globals are stored there, viz, constants: LocusCnt, PhenoCnt, MARKER_SCEME and
#'  tables: markers and unified_genotype_table (and possibly locus_allele_table)
#'
#' This same \emph{ENV} is returned by \emph{dbmega2_import}
#'
#' @return ENV an environment that contains all the tables created from the SQLite tables.
#' @export
#'
#' @examples
#'\dontrun{
#' ENV = getMega2ENV()
#'}
getMega2ENV = function () {
    return (ENV)
}

#' reset ENV environment; this flushes all tables/state.
#'
#' Mega2 uses an environment, \emph{ENV} to store all the tables it reads in.  This is reset
#'  to an empty env() and then the gc is run.
#'
#' @return None
#' @export
#'
#' @examples
#'\dontrun{
#' resetENV()
#'}
resetMega2ENV = function () {

    rm(list = ls(ENV), envir = ENV)

    unlockBinding("ENV", environment(resetMega2ENV))
    assign("ENV", new.env(parent = emptyenv()), environment(resetMega2ENV))
      lockBinding("ENV", environment(resetMega2ENV))

    gc(verbose = FALSE)

    ENV$refRanges  = refRanges
    ENV$refIndices = refIndices

    ENV$txdb       = "TxDb.Hsapiens.UCSC.hg19.knownGene"
    ENV$entrezGene = "org.Hs.eg.db"

    invisible()
}

#' show the association between mapno and mapname
#'
#' Mega2 allows several different physical and genetic maps to be stored and used to select
#'  distances.  This function show the association between mapno and mapname
#'
#' @return None
#' @export
#'
#' @examples
#'\dontrun{
#' showMapNames()
#'}
showMapNames = function () {
    ENV$mapnames_table[ , c(6, 2)]
}

## geno_i = inline::cxxfunction(
##     signature(a = "raw"), 

##     '
##     Rcpp::RawVector rv(a);

##     Rcpp::NumericVector v(rv.size() * 4);

##     int j = 0;
##     for (int i = 0; i < rv.size(); i++) {
##       int t0 = rv[i];
##       int t1 = (t0 & 0x03) >> 0;
##       int t2 = (t0 & 0x0c) >> 2;
##       int t3 = (t0 & 0x30) >> 4;
##       int t4 = (t0 & 0xc0) >> 6;
##       v[j++] = t1;
##       v[j++] = t2;
##       v[j++] = t3;
##       v[j++] = t4;
##     }
##     return v;
##     ',

##   plugin = "Rcpp")

# gg=geno_i(g)

# a1=f$AlleleName.x
# a2=f$AlleleName.y

#' getgenotype_person
#'
#' @description ...
#'
#' @param pid ...
#'
#' @return None
#' @export
#'
#' @examples
#'\dontrun{
#' # genotypes for all markers for n'th person in genotype table
#' getgenotype_person(1)
#'
#' # genotypes for all markers for range of persons in genotype table
#' getgenotype_person(m:n)
#'}
getgenotype_person = function(pid = 1) {

  a1 = ENV$allele_table[ENV$allele_table$indexX==1, 2][-1]
  a2 = ENV$allele_table[ENV$allele_table$indexX==2, 2][-1]

  rv = ENV$genotype_table[pid, 5][[1]]
  rv4 = getgenotypes_forperson(rv)

#  0 1|1
#  1 0|0
#  2 1|2
#  3 2|2

  return
    ifelse(rv4==0, paste0(a1, a1),
           ifelse(rv4==1, paste0("00"),
                  ifelse(rv4==2, paste0(a1, a2), paste0(a2, a2))
                  )
           )
}

################################################################

#' fetch genotype matrix for specified markers (assemble by rows)
#'
#' @description
#'  This function calls the C++ function that does all the heavy lifting.  It passes the
#'  locus_index and the locus_offset in the \emph{unified_genotype_table} from the
#'  \emph{markers_arg} argument.  It also gathers other data.frames that are in the "global"
#'  \bold{ENV} environment. One frame contains a bit vector of compressed genotype information,
#'  another contains the alleles for each marker, and finally there are some bookkeeping related
#'  data.
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
#' @return a matrix of genotypes represented as a nucleotide pair.  There is one column for each
#'  marker in \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @export
#' @useDynLib mega2
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector
#'  there are two bits for each genotype.  This function creates an output matrix by selecting
#'  from each row all the needed markers and then repeating for each person.
#'
#' @examples
#'\dontrun{
#' # genotypes for all persons in markers data.frame argument
#' getgenotypes_R(ENV$markers)
#'
#' # genotypes for all persons in chromosome n
#' getgenotypes_R(ENV$markers[ENV$markers$chromosome == n,])
#'}
getgenotypes_R = function(markers_arg) {

  return
    if (ENV$MARKER_SCHEME == 1) {
        getgenotypes_Ri(markers_arg$locus_link, markers_arg$locus_link_fill,
                        ENV$unified_genotype_table, ENV$allele_table, ENV$markerscheme_table,
                        ENV$PhenoCnt)
    } else {  # must be == 2

    }
}


#' fetch genotype matrix for specified markers
#'
#' @description
#'  This function calls the C++ function that does all the heavy lifting.  It passes the
#'  locus_index and the locus_offset in the \emph{unified_genotype_table} from the
#'  \emph{markers_arg} argument.  It also gathers other data.frames that are in the "global"
#'  \bold{ENV} environment. One frame contains a bit vector of compressed genotype information,
#'  another contains the alleles for each marker, and finally there are some bookkeeping related
#'  data.
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
#' @return a matrix of genotypes represented as two allele pairs.  There is one column for each
#'  marker in \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @export
#' @useDynLib mega2
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector
#'  there are two bits for each genotype.  This function creates an output matrix by fixing
#'  the marker and collecting genotype information for each person and then repeating for
#'  all the needed markers.  Currently, this appears slightly faster than the other scan in
#'  \emph{genotype_Ri}
#'
#' @examples
#'\dontrun{
#' # genotypes for all persons in markers data.frame argument
#' getgenotypes(ENV$markers)
#'
#' # genotypes for all persons in chromosome n
#' getgenotypes(ENV$markers[ENV$markers$chromosome == n,])
#'}
getgenotypes = function(markers_arg) {

  return
    if (ENV$MARKER_SCHEME == 1) {
        getgenotypes_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                       ENV$unified_genotype_table, ENV$allele_table, ENV$markerscheme_table,
                       ENV$PhenoCnt)
    } else {  # must be == 2
        getgenotypes_2(markers_arg$locus_link,
                       ENV$unified_genotype_table, ENV$locus_allele_table,
                       ENV$PhenoCnt)
    }
}

getgenotypes_C = getgenotypes

#' fetch genotype matrix for specified markers
#'
#' @description
#'  This function calls the C++ function that does all the heavy lifting.  It passes the
#'  locus_index and the locus_offset in the \emph{unified_genotype_table} from the
#'  \emph{markers_arg} argument.  It also gathers other data.frames that are in the "global"
#'  \bold{ENV} environment. One frame contains a bit vector of compressed genotype information,
#'  another contains the alleles for each marker, and finally there are some bookkeeping related
#'  data.
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
#' @return a matrix of genotypes represented as two allele pairs.  There is one column for each
#'  marker in \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @export
#' @useDynLib mega2
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector
#'  there are two bits for each genotype.  This function creates an output matrix by fixing
#'  the marker and collecting genotype information for each person and then repeating for
#'  all the needed markers.  Currently, this appears slightly faster than the other scan in
#'  \emph{genotype_Ri}
#'
#' @examples
#'\dontrun{
#' # two ints in upper/lower half integer representing allele for all persons in
#' # markers data.frame argument
#' getgenotypesraw(ENV$markers)
#'
#' # two ints in upper/lower half integer representing allele # for all persons in chromosome n
#' getgenotypesraw(ENV$markers[ENV$markers$chromosome == n,])
#'}
getgenotypesraw = function(markers_arg) {

  return
    if (ENV$MARKER_SCHEME == 1) {
        getgenotypesraw_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                          ENV$unified_genotype_table, ENV$allele_table, ENV$markerscheme_table,
                          ENV$PhenoCnt)
    } else {  # must be == 2
        getgenotypesraw_2(markers_arg$locus_link,
                          ENV$unified_genotype_table, ENV$locus_allele_table,
                          ENV$PhenoCnt)
    }
}

Rcpp::sourceCpp("src/getgenotypes.cpp")
