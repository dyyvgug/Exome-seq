#!/bin/bash
mkdir only_pass
cd only_pass
cat ./trio_snp_indel_VQSR.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' \
                > ./only_pass/trio_only_pass.vcf
