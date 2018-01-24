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

# devtools::document("Mega2R")
# devtools::install("Mega2R")
# devtools::check("Mega2R")

# devtools::run_examples("Mega2R")
# R CMD Rd2pdf Mega2R

# devtools::build_vignettes("Mega2R")


##../RR/cran_run



source("../RR/tst.R")

source("../RR/tst1.R")

source("../RR/tstseq.R") 

source("../RR/tstgena.R")

source("../RR/tstgenabel.R")

source("../RR/tstgenabe2.R")

source("../RR/tstpedgene.R") 

source("../RR/tstskat.R") 

source("../RR/tstvcf.R") 

system("cd ../RR/tstVCF; ../tstvcfdiff")
