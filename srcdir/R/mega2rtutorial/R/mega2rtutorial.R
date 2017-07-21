
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

#' mega2rtutorial package
#'
#' @description This package retrieves data stored in the Mega2rtutorial and 
#'	dumps them in the current directory.
#'
#' @author Robert V Baron
#' @docType package
#' @name mega2rtutorial-package
NULL

#' dump tutorial data
#'
#' @description
#'
#' @export
#' @return None
#'
#' @examples
#'\dontrun{
#' dump_mega2rtutorial_data()
#'}
dump_mega2rtutorial_data = function(dir = ".") {
    for (file in c("MEGA2.BATCH.seqsimr", "MEGA2.BATCH.srdta",
                   "MEGA2.BATCH.vcf",
                   "Mega2r.map", "Mega2r.ped", "seqsimr.db")) {
        from = system.file("exdata", file, package="mega2rtutorial")
        to   = paste(dir, file, sep="/")
        file.copy(from, to, copy.mode = TRUE, copy.date = TRUE)
    }
}
