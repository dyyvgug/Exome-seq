#!/bin/bash

for i in *.hg38_multianno.txt;   
do       
	sample=`echo $i|awk -F '.' '{print $1}'`; 
      	cut -f '1-10' $i|sed '1d'|sed "s/$/\t${sample}/">>merge.txt;   
done

sed -i '1s/^/Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tTumor_Sample_Barcode\n/' merge.txt
