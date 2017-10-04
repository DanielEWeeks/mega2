
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


FILES = c("MEGA2.BATCH.seqsimr", "MEGA2.BATCH.srdta", "MEGA2.BATCH.vcf",
          "Mega2r.map", "Mega2r.ped", "seqsimr.db", "srdta.db",
          "Mega2r.map.gz", "Mega2r.ped.gz", "seqsimr.db.gz", "srdta.db.gz")

GENED = c("SEQ.phe", "SEQ.tfam", "SEQ.tped", "SEQtped.raw",
          "SRD.phe", "SRD.tfam", "SRD.tped", "SRDtped.raw",
          "srdta.db.old", "srdta.map", "srdta.ped", "srdta.phe",
          "pedgene.txt")

FILES.gz = c("Mega2r.map", "Mega2r.ped", "seqsimr.db", "srdta.db")

#' dump tutorial data
#'
#' @description This function retrieves data stored in the Mega2rtutorial (inst/exdata).  It
#'	dumps them in the current (or specified) directory.
#'
#' @param dir The directory to store the tutorial data to.  By default, this is ".".
#'
#' @export
#' @return None
#'
#' @examples
#'\dontrun{
#' dump_mega2rtutorial_data()
#'}
dump_mega2rtutorial_data = function(dir = ".") {
    for (file in FILES) {
        from = system.file("exdata", file, package="Mega2R")
        to   = paste(dir, file, sep="/")
        file.copy(from, to, copy.mode = TRUE, copy.date = TRUE)
    }
    for (file in FILES.gz) {
#      R.utils::gunzip(file)
       in.gz = gzfile(paste0(file, ".gz"), "rb")
       out   = file(file, "wb")
       while (TRUE) {
           rv = readBin(in.gz, "raw", n = 4096, size = 1)
           ln = length(rv)
           writeBin(rv, out)
           if (ln != 4096) break;
       }
       close(in.gz)
       close(out)
    }
}

#' remove tutorial data
#'
#' @description This function removes the Mega2rtutorial (inst/exdata) data that was
#'	dumped in the specified directory.
#'
#' @param dir The directory to remove the tutorial data from.  By default, this is ".".
#'
#' @export
#' @return None
#'
#' @examples
#'\dontrun{
#' clean_mega2rtutorial_data()
#'}
clean_mega2rtutorial_data = function(dir = ".") {
    for (file in c(FILES, GENED)) {
        to = paste(dir, file, sep="/")
        unlink(to)
    }
    unlink("vcfr", recursive = TRUE)
}

