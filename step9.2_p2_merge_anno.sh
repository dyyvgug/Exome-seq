#!/bin/bash

for item in $(ls *.avinput)
do		
		#rare variants
		echo "annovar filter rare ${item%.*}"
		annotate_variation.pl \
		-thread 52 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_rare \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.01 -reverse \
		2>>anno_fil_rare.log
		
		#frequent variants
		echo "annovar filter frequent ${item%.*}"
		annotate_variation.pl \
		-thread 52 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_fre \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.05 \
		2>>anno_fil_fre.log

done

