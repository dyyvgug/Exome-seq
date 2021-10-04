#!/bin/bash

if [ ! -d ./rare ]
then mkdir ./rare
fi

for item in $(ls *.avinput)
do		
		#rare variants
		echo "annovar filter ${item%.*}"
		annotate_variation.pl \
		-thread 52 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./rare/${item%.*} \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.01 -reverse \
		2>>anno_fil.log

done

