
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



#library(DBI)
#library(RSQLite)

ENV = new.env(parent = emptyenv())

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
#        "locus_table",
         "allele_table",
         "marker_table",
         "markerscheme_table",
         "map_table",
         "mapnames_table",

#        "traitaff_table",
#        "affectclass_table",
#        "ptraitquant_table",

         "phenotype_table",
         "genotype_table"

## Derived tables
##        unified_genotype_table
##        markers
  )

#' mk_markers_with_skip
#'
#' @param mapselect ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
mk_markers_with_skip = function(mapselect = 1) {
    markersPerChr = sapply(split(ENV$marker_table$chromosome, ENV$marker_table$chromosome), length)
    extra_markers = cumsum(4 * floor((markersPerChr + 3) / 4) - markersPerChr)
    extra_markers = c(0, extra_markers)
    names(extra_markers) = NULL
    ENV$marker_table$locus_link_fill = ENV$marker_table$locus_link + extra_markers[ENV$marker_table$chromosome]

    ENV$markers = merge(ENV$marker_table[ , c("locus_link", "locus_link_fill", "MarkerName", "chromosome")],
                        ENV$map_table[ ENV$map_table$map == mapselect, c( "marker", "position")],
                        by.x = "locus_link", by.y = "marker")
}

#' mk_unified_genotype_table
#'
#' @param mapselect ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
mk_unified_genotype_table = function(mapselect = 1) {
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

    mk_markers_with_skip(mapselect)
}

#' dbmega2_import
#'
#' @param dbname ...
#' @param mapselect ...
#'
#' @return ENV
#' @importFrom RSQLite dbConnect dbExistsTable dbReadTable dbListFields SQLITE_RO
#' @export
#'
#' @examples
#'\dontrun{
#'}
dbmega2_import = function(dbname = "/Users/rbaron/mega2/test/mexnly/change_chrom/bcf/dbmega2.db",
                          mapselect = 1) {
    con = tryCatch(dbConnect(RSQLite::SQLite(), dbname = dbname, flags = SQLITE_RO),
                   error = function(xx) { stop("DB open failed: ", dbname, call. = FALSE) })

    for (tbl in TBLS) {
        if (dbExistsTable(con, tbl)) {
            assign(tbl, dbReadTable(con, tbl), pos = ENV)
            cat(tbl, dim(get(tbl, pos=ENV)), sep = "\t", end = "\n")
            cat(tbl, dbListFields(con, tbl), sep = "\t", end = "\n")
            cat(end = "\n")
        }
    }
    mk_unified_genotype_table(mapselect)
    return (ENV)
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
#' @param pid ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
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

#' getgenotypes_R
#'
#' @param markers ...
#'
#' @return NO
#' @export
#' @useDynLib mega2
#'
#' @examples
#'\dontrun{
#'}
getgenotypes_R = function(markers) {

  return
    getgenotypes_Ri(markers$locus_link, markers$locus_link_fill,
                    ENV$unified_genotype_table, ENV$allele_table, ENV$markerscheme_table,
                    ENV$int_table[ENV$int_table$key == 'PhenoCnt', ][1, 3])
 
}


#' getgenotypes
#'
#' @param markers ... 
#'
#' @return NO
#' @export
#' @useDynLib mega2
#'
#' @examples
#'\dontrun{
#'}
getgenotypes = getgenotypes_C = function(markers) {

  return
    getgenotypes_Ci(markers$locus_link, markers$locus_link_fill,
                    ENV$unified_genotype_table, ENV$allele_table, ENV$markerscheme_table,
                    ENV$int_table[ENV$int_table$key == 'PhenoCnt', ][1, 3])
 
}

Rcpp::sourceCpp("src/getgenotypes.cpp")
