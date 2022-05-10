#!/bin/bash

if [ ! -d ./anno ]
then mkdir ./anno
fi

for item in $(ls *.vcf)
do
        echo "annovar convert ${item%.*}"
        convert2annovar.pl -format vcf4 ./${item%.*}.vcf > ./${item%.*}.avinput
		
		#frequent variants
		echo "annovar annotate ${item%.*}"
        table_annovar.pl \
        -thread 52 \
        ./${item%.*}.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*} \
        -remove \
        -protocol refGene,clinvar_20200316 \
        -operation g,g \
        -nastring . \
        -vcfinput
done

