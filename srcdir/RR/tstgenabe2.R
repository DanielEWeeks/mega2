#   Mega2: Manipulation Environment for Genetic Analysis
#  
#   Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.
#  
#   Contributors to Mega2: Robert Baron and Daniel E. Weeks.
#  
#   This file is part of the Mega2 program, which is free software you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later
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
require("GenABEL")

#     1. GenABEL:::alleleID.codes()
#         [1] "12" "AB" "AT" "AG" "AC" "A-" "TA" "TG" "TC" "T-" "GA" "GT" "GC" "G-" "CA"
#        [16] "CT" "CG" "C-" "-A" "-T" "-G" "-C" "21" "BA" "ID" "DI" "11" "22" "BB" "II"
#        [31] "DD" "AA" "TT" "GG" "CC" "--"
#     2. GenABEL/src/GAlib/convert_snp_tped.cpp
#         1. bit order 1|2|3|4 on person;                 vs    4|3|2|1 on marker
#         2. 0/0->0, 1/1->1 (maj), 1/2->2, 2/2->3(min)    vs    0/0->1, 1/1->0, 1/2->2, 2/2->3

print("################################################################")
print("mex")
ENV=read.Mega2DB("~/mega2/data/mexnly.db", verbose = FALSE)
mega2=Mega2GenABEL(envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("msat")
ENV=read.Mega2DB("~/mega2/data/msatnly.db", verbose = FALSE)
# NOT biallelic

print("################################################################")
print("mome")
ENV=read.Mega2DB("~/mega2/data/momenly.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("tst")
ENV=read.Mega2DB("~/mega2/data/tstnly.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("mec")
ENV=read.Mega2DB("~/mega2/data/mecnly.db", verbose = FALSE)
mega2=Mega2GenABEL(envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("m15k")
ENV=read.Mega2DB("~/mega2/data/m15knly.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("samoaqwas")
#ENV=read.Mega2DB("~/mega2/data/samoaqwas.db", verbose = FALSE)
#mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome %in% c(1, 11, 21),], envir = ENV)
#mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome %in% c(1, 11, 21),], envir = ENV)
#Mega2GenABELtst(mega,mega2)

print("################################################################")
print("seqsimp")
ENV=read.Mega2DB("~/mega2/data/seqsimp.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[1:20000,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[1:20000,], envir = ENV)
Mega2GenABELtst(mega,mega2)
print("################################################################")
print("seqsimr")
ENV=read.Mega2DB("~/mega2/data/seqsimr.db", verbose = FALSE)
mega2=Mega2GenABEL(envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("yj1")
ENV=read.Mega2DB("~/mega2/data/yj1.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
Mega2GenABELtst(mega,mega2,FALSE)

print("yj2")
ENV=read.Mega2DB("~/mega2/data/yj2.db", verbose = FALSE)
mega2=Mega2GenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome != 999,], envir = ENV)
Mega2GenABELtst(mega,mega2,FALSE)
print("################################################################")
