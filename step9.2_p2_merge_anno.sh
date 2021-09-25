#!/bin/bash


cut -f '1-10' *.txt|sed '1d'>>for_maf.txt   

sed -i '1s/^/Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tTumor_Sample_Barcode\n/' for_maf.txt
