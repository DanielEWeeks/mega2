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

if( $1 == "") then
    echo
    echo "This script enables the construction of a simple reference allele file from vcf.gz reference panels for Mega2 to use."
    echo
    echo "First, this script uses BCFTools to parse vcf.gz files quickly, it can be obtained at http://www.htslib.org/download/."
    echo
    echo "Secondly, you must first obtain a reference panel, this script looks in the current directory for those files."
    echo
    echo "A reference panel can be constructed from any data available, or from distributed public data panels such as 1000 genomes, shapeit etc.  As long as those files are in vcf.gz format organized by chromosome, with standardized naming."
    echo
    echo "This script requires two command line arguements based on the vcf.gz files you choose to use."
    echo
    echo "The two command line arguemnts are 1) everything before the chromsome number, 2) everything after the chromosome number."
    echo "If data is provided in this form this script should loop over all the chromosomes and create a file of reference alles for chromosomes 1-22. "
    echo
    echo "    An exmaple: suppose I downloaded from 1000 genomes ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/, I would get files that look like: ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz."
    echo
    echo "In this case the script would be run with the command 'GetRefAlleles.sh ALL.chr .phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' which should output: 'RefAlleles.b37.ALL.chr1-22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.txt.gz' Which is a 3 collumn file of CHR POS and REF, which is gzipped to save size.  The build number is included for Mega2's purposes and is obtained from the contig field (if possible)."
    echo
    echo "This new file can be used in Mega2 under the menu option Reference Allele file."
    exit 0
endif

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
    echo set bcftools_program=$BCFTOOLS/bcftools
    set bcftools_program=$BCFTOOLS/bcftools
else
    echo The $BCFTOOLS/bcftools executable was not found -
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
set build = '.'
set noextension = $2:r:r
if ( -f $1$i$2 ) then
    set build=`$bcftools_program view -h $1$i$2 | grep -m 1 contig | awk -F, '{ print $2 }' | awk -F'[=&]' '{print $2}'`
endif

while ($i <= 22 )
    if ( -f $1$i$2 ) then
        echo "$1$i$2 was found"
        echo "Using BCF-Tools to extract Reference Alleles for Chromosome $i"
        $bcftools_program query -f '%CHROM %POS %REF\n' $1$i$2 >> RefAlleles.$build.$1$inner$noextension.txt
    else
        echo "$1$i$2 not found."
    endif
    @ i++
end

echo "Gzipping the result..."
gzip RefAlleles.$build.$1$inner$noextension.txt

echo "Done creating: RefAlleles.$1$inner$2.$noextension.txt.gz"
echo "This file is ready for use with Mega2 as a Reference Allele File."
