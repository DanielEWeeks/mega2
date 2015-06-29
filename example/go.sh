#
# NOTE:
# Mega2 should be compiled with make FLAGS=-DORDER_HETEROZYGOTE

Mega2 --nosave MEGA2.BATCH.annotated
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_annotated

Mega2 --nosave MEGA2.BATCH.bed
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_bed

Mega2 --nosave MEGA2.BATCH.ped
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_ped

Mega2 --nosave MEGA2.BATCH.post
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_post

Mega2 --nosave MEGA2.BATCH.pre
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_pre

Mega2 --nosave MEGA2.BATCH.preannotated
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_preannotated

Mega2 --nosave MEGA2.BATCH.bcf
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_bcf

Mega2 --nosave MEGA2.BATCH.vcf
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_vcf

Mega2 --nosave MEGA2.BATCH.impute
mv MEGA2.{ERR,KEYS,LOG} *.html ../example_output_impute

