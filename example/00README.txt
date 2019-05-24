This directory contains sample data in various formats:
     annotated bcf bcf2 bed impute ped post pre preannotated vcf
(Note: post is mega2 format postmakeped; pre and preannotated are mega2
 format premakeped)

For each input format, there are two MAKE.BATCH_xxx files: one that reads that
format's files, for example MAKE_BATCH_bed.  It stores the Mega2 Metadata in
../example_db_xxx (viz. ../example_db_bed) and stores the database file in
../example_db/dbXxx.db (viz. ../example_db/dbBed.db).  The second BATCH file
process the created database and prepares for Mendel analysis. The BATCH file
is named MEGA2.BATCH_xxx2mendel (MEGA2.BATCH_bed2mendel).  The Mega2 metadata
and the files needed for Mendel are in the ../example_output_xxx (viz.
../example_output_bed).
(Note: for format post, we also have a BATCH file that generates input for Cranefoot.
It is not run.)

Two lines in the go.sh script run the two BATCH files.  These lines are followed
by several comment lines that indicate which actual input files are read to
build the database.

If go.sh is given the argument "save", all the directories ../example_db_xxx and
../example_output_xxx are renamed to ../example_db_xxx.save and
../example_output_xxx.save respectively.
(Any original .save files are deleted.)

If go.sh is given the argument "diff", each directory is diff'ed with the older
directory.save version.  Expect to see different timestamps in the results but
that should be it.

NOTE: A long time ago, Mega2 ordered HETEROZYGOUS genotypes lexically.  This 
"feature" is not used any more if the input specifies an allele ordering, viz.
PLINK bed, BCF, BCF2, VCF, BCF, and IMPUTE.  If you prefer the original ordering for 
compatibility reasons with old data, compile Mega2 with FLAGS=-DSORT_HETEROZYGOUS

