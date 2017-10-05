
#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
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

go = function() {

    print("## run 200 ranges of pedgene to compare with YJiang's results; and profile")
    print("## result environment returned in ENV")


    library(Mega2R)

    envir = init_SKAT("~/mega2/data/yj1.db", verbose = TRUE)

    Mega2SKAT(envir$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), envir = envir)

    Mega2SKAT(envir$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), envir = envir, gene=c("ELL2", "CARD2", "ARMS2", "CFH"))

    Mega2SKAT(envir$phe[, 3] - 1 ~ 1, "D", kernel = "linear.weighted", weights.beta=c(0.5,0.5), envir = envir, gene=c("PLEKHA1", "MTHFD1L"))

#    10 PLEKHA1 4 124145585 124191871 0.000112462098709182
#    10 PLEKHA1 4 124151819 124191871 0.000112462098709182
#    10 PLEKHA1 4 124134094 124191871 0.000112462098709182

#    6 MTHFD1L 9 151187469 151423023 0.000597049939828251
#    6 MTHFD1L 9 151186815 151423023 0.000597049939828251
#    6 MTHFD1L 9 151196379 151423023 0.000597049939828251
#    6 MTHFD1L 9 151186815 151413723 0.000597049939828251


#    MTHFD1L exm586499 1274 1 0 
#    MTHFD1L exm2266269 369 657 249 
#    MTHFD1L rs11754661 1132 143 0 
#       chr    gene nvariants     start       end       skat
#    79   6 MTHFD1L         3 151186815 151220136 0.09650606

#    MTHFD1L exm586499 1274 1 0 
#    MTHFD1L exm2266269 369 657 249 
#    MTHFD1L rs11754661 1132 143 0 
#    MTHFD1L exm586522 1275 0 0 
#    MTHFD1L exm2257596 1275 0 0 
#    MTHFD1L exm586534 1273 2 0 
#    MTHFD1L rs6922269 667 510 98 
#    MTHFD1L exm586543 1273 2 0 
#    MTHFD1L exm586578 1275 0 0 
#    MTHFD1L exm586580 1230 45 0 
#    MTHFD1L exm586582 1274 1 0 
#    MTHFD1L exm586591 1275 0 0 
#    MTHFD1L exm2262071 607 527 141 
#       chr    gene nvariants     start       end         skat
#    80   6 MTHFD1L         9 151186815 151413723 0.0005015113

#    PLEKHA1 rs6585827 474 593 208 
#    PLEKHA1 exm861540 1272 3 0 
#    PLEKHA1 exm861542 637 502 136 
#    PLEKHA1 exm861558 1274 1 0 
#       chr    gene nvariants     start       end         skat
#    81  10 PLEKHA1         4 124134094 124191871 0.0001276224


        envir

}

ENV=go()

