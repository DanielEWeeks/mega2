#!/bin/nawk -f
# Rearrange the map file so that the locus name is in the
# first column.
BEGIN { print("locus   chromosome      cm") } 
{ if (FNR > 1) printf("%-10s %5d %10.2f\n", $3,$1,$2) }
