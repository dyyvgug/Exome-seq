#!/bin/bash

mkdir mutation
for item in $(ls *.sam)
do

	# use HaplotypeCaller to identify variance
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
		-ERC GVCF \
		--native-pair-hmm-threads 60 \
		-R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
		-I ${item%.*}_bqsr.bam \
		-O ./mutation/${item%.*}.g.vcf.gz \
		2>>hc_log.txt
	
done
