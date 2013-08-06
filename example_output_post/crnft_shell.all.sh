#!/bin/csh -f
# C-shell file name: ../example_output_post/crnft_shell.all.sh
#----------------------------------------------
#   Mega2 version 4.5.5
#   Run date:                2012-6-21-10-30
#   This script created on   Thu Jun 21 10:30:37 2012
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
 