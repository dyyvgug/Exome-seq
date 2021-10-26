#!usr/bin/bash
 
gatk=gatk
 
#references
ref=/proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
gatk_bundle=/proj/y.dong/GATK/hg38

dbsnp=$gatk_bundle/dbsnp_146.hg38.vcf.gz
indel=$gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
G1000=$gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
hapmap=$gatk_bundle/hapmap_3.3.hg38.vcf.gz
omni=$gatk_bundle/1000G_omni2.5.hg38.vcf.gz
mills=$gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

if [ ! -d ./filter ]
then mkdir ./filter
fi

cd filter
# SNP mode
time $gatk VariantRecalibrator \
    -R $ref \
    -V ../cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    --resource:omini,known=false,training=true,truth=false,prior=12.0 $omni \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 $G1000 \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
    -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.0 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
    --rscript-file trio_snps.plots.R \
    --tranches-file trio_snps.tranches \
    -O trio_snps.recal \
    2>>snp.log
 
time $gatk ApplyVQSR \
    -R $ref \
    -V ../cohort.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file trio_snps.tranches \
    --recal-file trio_snps.recal \
    -mode SNP \
    -O trio_snps_VQSR.vcf.gz \
    2>>snp_apply.log
 
# Indel mode
time $gatk VariantRecalibrator \
    -R $ref \
    -V trio_snps_VQSR.vcf.gz \
    --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=15.0 $mills \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    --rscript-file trio_indels.plots.R \
    --tranches-file trio_indels.tranches \
    -O trio_indels.recal \
    2>>indel.log
 
time $gatk ApplyVQSR \
    -R $ref \
    -V trio_snps_VQSR.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file trio_indels.tranches \
    --recal-file trio_indels.recal \
    -mode INDEL \
    -O trio_snp_indel_VQSR.vcf.gz \
    2>>indel_apply.log
