#!/bin/csh -f

#  GetRefAlleles.sh
#  Allows the construction of a simple reference allele file
#  The output will be of the format:
#  CHR POS REF
#  This exists for the purposes of creating
#  reference alleles from external datasets
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

if ( $?BCFTOOLS ) then
    set BCFTOOLS_def=1
else
    set BCFTOOLS_def=0
endif
echo
if ( "`type -t bcftools`" == "file" ) then
    echo set bcftools_program=`type -p bcftools`
    set bcftools_program=`type -p bcftools`
else if ( $BCFTOOLS_def && -x "BCFTOOLS/bcftools" ) then
    echo set bcftools_program=BCFTOOLS/bcftools
    set bcftools_program=BCFTOOLS/bcftools
else
    echo The BCFTOOLS/bcftools executable was not found -
    echo please set your BCFTOOLS environment variable properly so bcftools can be found.
    echo
        if (! $BCFTOOLS_def) then
            echo BCFTOOLS is not defined.
        else
            echo BCFTOOLS is set to "$BCFTOOLS".
        endif
    echo
    echo If using Bash and ksh you would use something like this:
    echo export BCFTOOLS=dir_to_bcftools
    echo
    echo If using csh you would use something like this:
    echo setenv BCFTOOLS dir_to_bcftools
    echo
    exit 0
endif

echo
set i = 1
set inner = '1-22'
while ($i <= 22 )
    if ( -f $1$i$2 ) then
        echo "$1$i$2 was found"
        echo "Using BCF-Tools to extract Reference Alleles for Chromosome $i"
$bcftools_program query -f '%CHROM %POS %REF\n' $1$i$2 >> reference_alleles_$1$inner$2
    else
        echo "$1$i$2 not found."
    endif
    @ i++
end
