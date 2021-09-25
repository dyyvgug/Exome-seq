#!/bin/bash

for item in $(ls *.vcf)
do

        echo "annovar convert ${item%.*}"
        # Strict filtering
        convert2annovar.pl -format vcf4 ./${item%.*}_filter.vcf > ./${item%.*}_filter.avinput
		# General filtering
		convert2annovar.pl -format vcf4 ./${item%.*}_somatic.vcf > ./${item%.*}_somatic.avinput
		
		# Strict filtering
		echo "annovar annotate ${item%.*}"
        table_annovar.pl \
        ./filter/${item%.*}_filter.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*} \
        -remove \
        -protocol refGene,clinvar_20200316,cytoBand \
        -operation g,f,r \
        -nastring . \
        -vcfinput
        # General filtering
		echo "annovar ${item%.*}"
        table_annovar.pl \
        ./filter/${item%.*}_somatic.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*}_ge \
        -remove \
        -protocol refGene,clinvar_20200316,cytoBand \
        -operation g,f,r \
        -nastring . \
        -vcfinput
		
done
mkdir for_maf_ge
mkdir for_maf
mv *ge.*txt for_maf_ge/
mv *.txt for_maf/

