
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

goo = function(rng = 10000:11000) {

    ENV=dbmega2_import("ped1.db")

    cc = getgenotypes(ENV$markers[rng,])
    c1r = getgenotypesraw(ENV$markers[rng,])

    c1x = array("0", dim=dim(cc))
    c1x [c1r == 131074] = "22"
    c1x [c1r == 65538] = "12"
    c1x [c1r == 65537] = "11"
    cat("all(cc == c1x) ")
    print(all(cc == c1x))

    ENV=dbmega2_import("ped2.db")
    c2a = getgenotypes(ENV$markers[rng,])
    cat("all(ccr==c2a) ")
    print(all(cc == c2a))
    c2r = getgenotypesraw(ENV$markers[rng,])
    cat("all(c1r==c2r) ")
    print(all(c1r == c2r))

    c2x = array("0", dim=dim(cc))
    c2x [c2r == 131074] = "22"
    c2x [c2r == 65538] = "12"
    c2x [c2r == 65537] = "11"
    cat("all(cc==c2x) ")
    print(all(cc==c2x))
}
