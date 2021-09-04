#!/bin/bash

for item in $(ls *.sam)
do
		AddOrReplaceReadGroups \
			I=${item%.*}_dup_fixed.bam \
			O=${item%.*}_dup_fixed_add.bam \
			RGID=${item%.*} \
			RGLB=WES \
			RGPL=ILLUMINA \
			RGPU=${item%.*} \
			RGSM=${item%.*}

done
