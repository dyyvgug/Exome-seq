#!/bin/bash

mkdir stat
for item in $(ls *.sam)
do
	echo "sam_${item%.*}"
	samtools sort -@ 30 -o ${item%.*}.bam ${item%.*}.sam
	samtools flagstat ${item%.*}.bam > ${item%.*}.bam.flag
done
