#!/bin/bash

for item in $(ls *.sam)
do
	echo "coverage_${item%.*}"
			
			bedtools genomecov -d -ibam ${item%.*}.bam -g /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta > ${item%.*}_cov.txt
			qualimap bamqc -bam ${item%.*}.bam -outfile ${item%.*}.pdf -outformat PDF
done
