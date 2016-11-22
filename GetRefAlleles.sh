#!/bin/sh

#  GetRefAlleles.sh
#  Allows the construction of a simple reference allele file
#  The output will be of the format:
#  CHR POS REF
#  This exists for the purposes of creating
#  reference alleles from external datasets
#
#  $1 PATH to BCF-tools head directory
#  This is used to get /bcftools
#
#  $2, $3 the first and second half respectively
#  for the vcf.gz filenames (before and after chr#)
#  Generally data will be split by chromosome 1 2 3...
#  This is designed to loop over files to hit all the vcf files
#  For instance if your files were called:
#  ALL.chr2.phase3_genotypes.vcf.gz
#  So your call to this script would be:
#  bash GetRefAlleles.sh ./bcftools-1.3.1 ALL.chr .phase3_genotypes.vcf.gz
#
#  This format is based on 1000 Genomes data that can be found:
#  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/
#  But should be usable for your own data if organized correctly
#  If your data is in one directory named in a onsistent way
#  And split by chromosome into vcf/vcf.gz files you can choose
#  Reference Alleles using this script.

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    if [ -f $2$i$3 ]
    then
        echo "$2$i$3 was found"
        echo "Using BCF-Tools to extract Reference Alleles for Chromosome $i"
        $1/bcftools query -f '%CHR %POS %REF\n' $2$i$3 >> ref_alleles.txt
    else
        echo "$2$i$3 not found."
    fi
done

#What to do for X Y XY MT?  Is there a consistent easy way to do this?
