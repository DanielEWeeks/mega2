#
# NOTE:
# Mega2 should be compiled with make FLAGS=-DORDER_HETEROZYGOTE

if [[ "$1" == "diff" ]] ; then
    for i in annotated bcf bcf2 bed impute ped post pre preannotated vcf
    do
        echo diffing ../example_db_$i with ../example_db_$i.save
        diff -rsb ../example_db_$i ../example_db_$i.save

        echo diffing ../example_output_$i with ../example_output_$i.save
        diff -rsb ../example_db_$i ../example_db_$i.save
    done
    exit 0
fi

if [[ "$1" == "save" ]] ; then
    for i in annotated bcf bcf2 bed impute ped post pre preannotated vcf
    do
        echo moving ../example_db_$i to ../example_db_$i.save
        mv ../example_db_$i ../example_db_$i.save

        echo moving ../example_output_$i to ../example_output_$i.save
        mv ../example_db_$i ../example_db_$i.save
    done
fi

if [[ ! -d ../example_db ]] ; then
    mkdir ../example_db
fi

echo Mega2 should be compiled with make FLAGS=-DORDER_HETEROZYGOTE

mega2 --nosave MEGA2.BATCH_annotated
mega2 --nosave MEGA2.BATCH_annotated2mendel
#pedin.annotated
#names.annotated
#map.annotated
#frequency.annotated
#penetrance.annotated
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_annotated

mega2 --nosave MEGA2.BATCH_bed
mega2 --nosave MEGA2.BATCH_bed2mendel 
#bed.fam
#bed.bim
#bed.bed
#ped.phe
#ped.frequency
#ped.penetrance
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_bed

mega2 --nosave MEGA2.BATCH_ped
mega2 --nosave MEGA2.BATCH_ped2mendel 
#ped.ped
#ped.map
#ped.phe
#ped.frequency
#ped.penetrance
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_ped

mega2 --nosave MEGA2.BATCH_post
mega2 --nosave MEGA2.BATCH_post2cranefoot
#datain.ex
#pedin.ex
#omit.ex
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_post

mega2 --nosave MEGA2.BATCH_pre
mega2 --nosave MEGA2.BATCH_pre2mendel
#pedin.pre.05
#map.05
#datain.05
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_pre

mega2 --nosave MEGA2.BATCH_preannotated
mega2 --nosave MEGA2.BATCH_preannotated2mendel
#pedin.preannotated
#names.preannotated
#map.preannotated
#frequency.annotated
#penetrance.annotated
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_preannotated

mega2 --nosave MEGA2.BATCH_bcf
mega2 --nosave MEGA2.BATCH_bcf2mendel
#bed.fam
#map.preannotated
#study.bcf
#study.phe
#study.freq
#study.pen
#study.bcf.bcfidx
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_bcf

mega2 --nosave MEGA2.BATCH_vcf
mega2 --nosave MEGA2.BATCH_vcf2mendel
#bed.fam
#map.preannotated
#study.vcf
#study.freq
#study.pen
#study.phe
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_vcf

mega2 --nosave MEGA2.BATCH_impute
mega2 --nosave MEGA2.BATCH_impute2mendel
#impute.impute
#map.preannotated
#impute.sample
#study.freq
#study.pen
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_impute

mega2 --nosave MEGA2.BATCH_bcf2
mega2 --nosave MEGA2.BATCH_bcf22mendel
#bed.fam
#map.preannotated
#studyBCF2.bcf
#studyBCF2.bcf.csi
#study.phe
#study.freq
#study.pen
#mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_bcf


exit 0


map.ex
pedin.05
study.fam
study.map
study.vcf.vcfidx
study1.6.chr06.bcf
study1.6.map
