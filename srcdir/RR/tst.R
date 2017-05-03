
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

go = function(db = 0) {
#    if (db) browser()

    library(mega2)

    if (TRUE) {
##  ENV = dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db", verbose = 1)
        source("../R/mega2/tests/mega2_test.R")
        tst1()

        source("../R/mega2/tests/mega2gene_test.R")
        tst10()
        tst11()
        tst12()
    }

}

go(1)

