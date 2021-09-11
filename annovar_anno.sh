#!/bin/bash
for item in $(ls *.vcf)
do
        echo "annovar SNP ${item%.*}"
        table_annovar.pl \
        ./filter/${item%.*}_filter.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*}_snp \
        -remove \
        -protocol refGene,clinvar_20200316,cytoBand \
        -operation g,f,r \
        -nastring . \
        -vcfinput

done