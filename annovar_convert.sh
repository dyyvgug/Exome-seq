#!/bin/bash
for item in $(ls *.vcf)
do

        echo "annovar convert ${item%.*}"
        
        convert2annovar.pl -format vcf4 ./${item%.*}_filter.vcf > ./${item%.*}_filter.avinput
 
done


