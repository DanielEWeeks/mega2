#!/bin/csh -f

#   Mega2: Manipulation Environment for Genetic Analysis.
#  
#   Copyright 2017-2019, University of Pittsburgh. All Rights Reserved.
#  
#   Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   and Daniel E. Weeks.
#  
#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

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
    echo "This script enables the construction of a simple reference allele file"
    echo "from vcf.gz reference panels for Mega2 to use."
    echo
    echo "First, this script uses BCFTools to parse vcf.gz files quickly, it can"
    echo " be obtained at http://www.htslib.org/download/."
    echo
    echo "Secondly, you must obtain a reference panel, this script looks"
    echo "in the current directory for those files."
    echo
    echo "A reference panel can be constructed from any data available, or from"
    echo "distributed public data panels such as 1000 genomes, Minimac3, etc."
    echo "As long as those files are in vcf.gz format organized by chromosome,"
    echo "with standardized naming."
    echo
    echo "This script requires two command line arguments based on the vcf.gz"
    echo "files you choose to use."
    echo
    echo "The arguments are 1) everything before the chromosome number,"
    echo "2) everything after the chromosome number."
    echo "If data is provided in this form this script should loop over all the "
    echo "chromosomes and create reference allele file for chromosomes 1-22."
    echo
    echo "An example: suppose I downloaded from 1000 genomes "
    echo "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, I would 
    echo "get files that look like: "
    echo "ALL.chr1.phase3.20130502.genotypes.vcf.gz."
    echo
    echo "In this case the script would be run with the command "
    echo "'GetRefAlleles.sh ALL.chr .phase3.20130502.genotypes.vcf.gz' "
    echo "the script will output: "
    echo "'RefAlleles.b37.ALL.chr1-22.phase3.20130502.genotypes.vcf.gz.txt.gz'"
    echo "Which is a 3 column file of CHR POS and REF, that is Gzipped"
    echo "to save size.  The build number is included for Mega2's purposes"
    echo "and is obtained from the contig field (if possible)."
    echo
    echo "This new file can be used in Mega2 under the menu option called"
    echo "Reference Allele file."
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
else
    if( $BCFTOOLS_def && -x "$BCFTOOLS/bcftools" ) then
        echo set bcftools_program=$BCFTOOLS/bcftools
        set bcftools_program=$BCFTOOLS/bcftools
    else
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
        $bcftools_program query -f '%CHROM %POS %REF %ALT{0}\n' $1$i$2 >> RefAlleles.$build.$1$inner$noextension.txt
    else
        echo "$1$i$2 not found."
    endif
    @ i++
end

echo "Gzipping the result..."
gzip RefAlleles.$build.$1$inner$noextension.txt

echo "Done creating: RefAlleles.$build.$1$inner$noextension.txt"
echo "This file is ready for use with Mega2 as a Reference Allele File."
