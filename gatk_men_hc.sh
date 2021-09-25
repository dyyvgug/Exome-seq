#!/bin/bash

mkdir mutation
ref=/proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
for item in $(ls *.sam)
do
	# use HaplotypeCaller to identify variance
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
		-ERC GVCF \
		--native-pair-hmm-threads 60 \
		-R $ref \
		-I ${item%.*}_bqsr.bam \
		-O ./mutation/${item%.*}.g.vcf.gz \
		2>>hc_log.txt
	
done

cd mutation
j=0
for i in $(ls *.gz)
do
	gvcf_list[j]=$i
	j=`expr $j + 1`
done

gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" CombineGVCFs \
	-R $ref \
	-V ./${gvcf_list[0]} \
	-V ./${gvcf_list[1]} \
	-V ./${gvcf_list[2]} \
	-O ./cohort.g.vcf.gz \
	2>>./combine_log.txt
	
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
	-R $ref \
	-V ./cohort.g.vcf.gz \
	-O ./cohort.vcf.gz \
	2>>./genotype_log.txt
