
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

library(mega2r)

print("################################################################")
print("mex")
ENV=read.Mega2DB("~/mega2/data/mexnly.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/mex", envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("msat")
ENV=read.Mega2DB("~/mega2/data/msatnly.db", verbose = FALSE)
# NOT biallelic
#mega2=Mega2GenABEL("tstGNA/msat/vcf/vcf.16", ENV$markers[ENV$markers$chromosome == 16,], alleleOrder = "default", envir = ENV)
#mega2=Mega2GenABEL("tstGNA/msat/vcfmin/vcf.16", ENV$markers[ENV$markers$chromosome == 16,], alleleOrder = "minor", envir = ENV)
#mega2=Mega2GenABEL("tstGNA/msat/vcfmaj/vcf.16", ENV$markers[ENV$markers$chromosome == 16,], alleleOrder = "major", envir = ENV)

print("################################################################")
print("mome")
ENV=read.Mega2DB("~/mega2/data/momenly.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/mome/mome", ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("tst")
ENV=read.Mega2DB("~/mega2/data/tstnly.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/tst/tst", ENV$markers[ENV$markers$chromosome == 22,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome == 22,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("mec")
ENV=read.Mega2DB("~/mega2/data/mecnly.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/mec/mec", envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("m15k")
ENV=read.Mega2DB("~/mega2/data/m15knly.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/m15k/m15k", ENV$markers[ENV$markers$chromosome == 16,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome == 16,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("samoaqwas")
#ENV=read.Mega2DB("~/mega2/data/samoaqwas.db", verbose = FALSE)
#mega2=Mega2GenABEL("tstGNA/samoaqwas/samoaqwas", ENV$markers[ENV$markers$chromosome %in% c(1, 11, 21),], envir = ENV)
#mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome %in% c(1, 11, 21),], envir = ENV)
#Mega2GenABELtst(mega,mega2)

print("################################################################")
print("seqsimp")
ENV=read.Mega2DB("~/mega2/data/seqsimp.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/seqsimp/seqsimp", ENV$markers[1:20000,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[1:20000,], envir = ENV)
Mega2GenABELtst(mega,mega2)
print("################################################################")
print("seqsimr")
ENV=read.Mega2DB("~/mega2/data/seqsimr.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/seqsimr/seqsimr", envir = ENV)
mega =Mega2ENVGenABEL(envir = ENV)
Mega2GenABELtst(mega,mega2)

print("################################################################")
print("yj1")
ENV=read.Mega2DB("~/mega2/data/yj1.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/yj/yj1", ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
Mega2GenABELtst(mega,mega2)

print("yj2")
ENV=read.Mega2DB("~/mega2/data/yj2.db", verbose = FALSE)
mega2=Mega2GenABEL("tstGNA/yj/yj2", ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
mega =Mega2ENVGenABEL(ENV$markers[ENV$markers$chromosome == 1,], envir = ENV)
Mega2GenABELtst(mega,mega2)
print("################################################################")
