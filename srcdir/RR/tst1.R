
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

## test compressions 1 & 2 and raw vs neucleotide

goo = function(rng = 10000:11000, file="~/mega2/data/yj1.db") {

    print("test compressions 1 & 2 and raw vs neucleotide alleles")
    
    ENV=dbmega2_import(file)

# get letter and raw compression 1
    cc = getgenotypes(ENV$markers[rng,], envir=ENV)
    c1raw = getgenotypesraw(ENV$markers[rng,], envir=ENV)

    c1cnv = array("0", dim=dim(cc))
    c1cnv [c1raw == 131074] = "22"
    c1cnv [c1raw == 65538] = "12"
    c1cnv [c1raw == 65537] = "11"
    c1cnv [c1raw == 0]     = "00"
    cat("all(cc == c1cnv) ")
    print(all(cc == c1cnv))

# get letter and raw compression 2
    ENV=dbmega2_import("~/mega2/data/yj2.db")
    c2a = getgenotypes(ENV$markers[rng,], envir=ENV)
#   cmp compression 1 and 2
    cat("all(cc==c2a) ")
    print(all(cc == c2a))
    c2raw = getgenotypesraw(ENV$markers[rng,], envir=ENV)
    cat("all(c1raw==c2raw) ")
    print(all(c1raw == c2raw))

    c2cnv = array("0", dim=dim(cc))
    c2cnv [c2raw == 131074] = "22"
    c2cnv [c2raw == 65538] = "12"
    c2cnv [c2raw == 65537] = "11"
    c2cnv [c2raw == 0]     = "00"
    cat("all(cc==c2cnv) ")
    print(all(cc==c2cnv))
}

goo()

rm(goo)
