#!/bin/csh -f
# C-shell file name: ../example_output_post/crnft_shell.all.sh
#----------------------------------------------
#   Mega2 version 5.0.2
#   Run date:                2019-5-19-12-08
#   This script created on   Sun May 19 12:08:23 2019
#   Input file names:
#       Pedigree file:              pedin.ex
#          Locus file:              datain.ex
#            Map file:              map.ex
#           Omit file:              omit.ex
#    Untyped pedigree option  Include all pedigrees whether typed or not
#----------------------------------------------
echo Running CRANEFOOT on chromosome 5
 if (-e crnft_chr5.ps) then
    /bin/rm crnft_chr5.ps
 endif
 if (-e crnft_chr5.topology.txt) then
    /bin/rm crnft_chr5.topology.txt
 endif

 cranefoot crnft_control.05
 echo Running CRANEFOOT on chromosome 6
 if (-e crnft_chr6.ps) then
    /bin/rm crnft_chr6.ps
 endif
 if (-e crnft_chr6.topology.txt) then
    /bin/rm crnft_chr6.topology.txt
 endif

 cranefoot crnft_control.06
 