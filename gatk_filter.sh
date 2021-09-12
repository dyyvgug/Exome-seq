#!/bin/bash
mkdir filter
for item in $(ls *.vcf)
do

        echo "gatk filter ${item%.*}"
        gatk FilterMutectCalls \
                -R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
                -V ${item%.*}.vcf \
                -O ./filter/${item%.*}_somatic.vcf
		# Strict filtering, only take PASS lable(site contains at least one allele that passes filters)
                cat ./filter/${item%.*}_somatic.vcf \
                | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' \
                > ./filter/${item%.*}_filter.vcf
done
