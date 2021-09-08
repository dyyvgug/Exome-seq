#!/bin/bash
# Author: Yingying Dong. 2021-9. Analyzing  SRP161770.
mkdir mutation

numbers=( 'SRR7829660' 'SRR7829678' 'SRR7829664' 'SRR7829675' 'SRR7829663' 'SRR7829670' 'SRR7829668' 'SRR7829673' 'SRR7829659' 'SRR7829677' 'SRR7829665' 'SRR7829676' 'SRR7829667' 'SRR7829672' 'SRR7829666' 'SRR7829669' 'SRR7829662' 'SRR7829671')

i=0
j=0
until [ ! $i -lt 8 ]
do

echo "patient $i is ${numbers[$j]} "
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" Mutect2 \
	-R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
	-I ${numbers[$j]}_bqsr.bam \
	-I ${numbers[`expr $j + 1`]}_bqsr.bam \
	-tumor  ${numbers[$j]}\
	-normal ${numbers[`expr $j + 1`]} \
	--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
	-O ./mutation/${numbers[$j]}_mut.vcf \
	2>> mutect2_log.txt
i=`expr $i + 1`
j=`expr $j + 2`

done

