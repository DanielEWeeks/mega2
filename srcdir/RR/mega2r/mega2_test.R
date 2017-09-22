
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

library(Mega2R)

#' allelediff
#'
#' @param aa ...
#' @param bb ...
#' @param n ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
allelediff = function(aa=aa, bb=bb, n=24) {
    for (i in 1:24) {
        print(sum(
                  ( (substr(aa[, i], 1, 1) == substr(bb[ , i], 1, 1)) &
                    (substr(aa[, i], 2, 2) == substr(bb[ , i], 2, 2)) )  |

                  ( (substr(aa[, i], 1, 1) == substr(bb[ , i], 2, 2)) &
                    (substr(aa[, i], 2, 2) == substr(bb[ , i], 1, 1)) )
                  )
              )
      }
}

################################################################

#' tst1
#'
#' @param db ...
#'
#' @return NO
#' @export
#'
#' @examples
#'\dontrun{
#'}
tst1 = function(ENV) {

    print("## getgenotypes by column and by row; compare two ranges")

    cat("getgenotypes",   system.time((cc=getgenotypes  (ENV$markers[1:1000,], "", ENV) )), "\n")
    cat("getgenotypes_R", system.time((dd=Mega2R:::getgenotypes_R(ENV$markers[1:1000,], "", ENV) )), "\n")
                
    print(all(cc==dd))

    cat("getgenotypes",   system.time((cc=getgenotypes  (ENV$markers[500000:501000,], envir=ENV) )), "\n")
    cat("getgenotypes_R", system.time((dd=Mega2R:::getgenotypes_R(ENV$markers[500000:501000,], envir=ENV) )), "\n")
    print(all(cc==dd))

# 25X slower
#   cat("getgenotypes_Cbind", system.time((ee=getgenotypes_Cbind(1:1000))), "\n")
#   print(all(cc==ee))
#   cat("getgenotypes_Cbind", system.time((ee=getgenotypes_Cbind(500000:501000))), "\n")
#   print(all(cc==ee))
    
    ## getgenotypes_C 0.845 0.008 0.854 0 0 
    ## getgenotypes_R 1.234 0.023 1.265 0 0 
    ## [1] TRUE
    ## getgenotypes_C 0.965 0.018 0.984 0 0 
    ## getgenotypes_R 1.211 0.018 1.234 0 0 
    ## [1] TRUE

}
################################################################
