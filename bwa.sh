#!/bin/bash

mkdir aligned

for item in $(ls *.sra)
do
	echo "bwa_${item%.*}"
	bwa mem -t 52 -M -R "@RG\tID:$sample\tSM:$sample\tLB:ES\tPL:Illumina" /DYY/Homo_sapiens/hg38/BWAIndex/genome.fa ./clipper_fastq/${item%.*}_1.fq ./clipper_fastq/${item%.*}_2.fq  > ./aligned/${item%.*}.sam 2>> ./aligned/mapping_repo.txt
	
done
