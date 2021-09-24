#!/bin/bash

mkdir mutation
for item in $(ls *.sam)
do
	ref=/proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
	# use HaplotypeCaller to identify variance
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
		-ERC GVCF \
		--native-pair-hmm-threads 60 \
		-R $ref \
		-I ${item%.*}_bqsr.bam \
		-O ./mutation/${item%.*}.g.vcf.gz \
		2>>hc_log.txt
	
done

j=0
for i in $(ls ./mutation/*.gz)
do
	gvcf_list[j]=$i
	j=`expr $j + 1`
done

gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" CombineGVCFs \
	-R $ref
	-V ./mutation/${gvcf_list[0]}
	-V ./mutation/${gvcf_list[1]}
	-V ./mutation/${gvcf_list[2]}
	-O ./mutation/cohort.g.vcf.gz
	2>>./mutation/combine_log.txt
	
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
	-R $ref
	-V ./mutation/cohort.g.vcf.gz
	-O ./mutation/cohort.vcf.gz
	2>>./mutation/genotype_log.txt
