#!/bin/bash

for item in $(ls *.sra)
do
	echo "sra_${item%.*}"

   fastq-dump --gzip --split-3 -O ./fastq/ -A ${item%.*}.sra
done

# QC
cd ./fastq/
fastqc ./*
multiqc ./*


