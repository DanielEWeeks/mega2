#!/bin/nawk -f
##### GENOTYPE ERROR: Pedigree 100   Locus  102.   Name D3S2432. #####
/GENOTYPE ERROR/ { printf("%-d\t0\t%s\n",$5,$9) }
