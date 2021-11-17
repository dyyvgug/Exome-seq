#!/bin/bash

for item in $(ls *.sam)
do
	echo "coverage_${item%.*}"
		gatk \
			DepthOfCoverage \
			-R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
			-O ${item%.*} \
			-I ${item%.*}.bam \
			-L chr1
done
