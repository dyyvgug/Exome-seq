#!usr/bin/bash
 
gatk=gatk
 
#references
ref=/proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
gatk_bundle=/proj/y.dong/GATK/hg38

dbsnp=$gatk_bundle/dbsnp_146.hg38.vcf.gz
indel=$gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
G1000=$gatk_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
hapmap=$gatk_bundle/hapmap_3.3.hg38.vcf.gz
omini=$gatk_bundle/1000G_omni2.5.hg38.vcf.gz
mills=$gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

if [ ! -d ./filter ]
then mkdir ./filter
fi
 
# SNP mode
time $gatk VariantRecalibrator \
    -R $ref \
    -V cohort.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
    --resource:omini,known=false,training=true,truth=false,prior=12.0 $omini \
    --resource:1000G,known=false,training=true,truth=true,prior=10.0 $G1000 \
    --resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp \
    -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 95.0 -tranche 92.0 -tranche 90.0 -tranche 80.0 \
    --rscript-file ./filter/trio_snps.plots.R \
    --tranches-file ./filter/trio_snps.tranches \
    -O ./filter/trio_snps.recal \
    2>>./filter/snp_log.txt
 
time $gatk ApplyVQSR \
    -R $ref \
    -V cohort.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ./filter/trio_snps.tranches \
    --recal-file ./filter/trio_snps.recal \
    -mode SNP \
    -O ./filter/trio_snps_VQSR.vcf.gz \
    2>>./filter/snp_apply_log.txt
 
# Indel mode
time $gatk VariantRecalibrator \
    -R $ref \
    -V ./filter/trio_snps_VQSR.vcf.gz \
    --max-gaussians 4 \
    -resource:mills,known=false,training=true,truth=true,prior=15.0 $mills \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode INDEL \
    --rscript-file ./filter/trio_indels.plots.R \
    --tranches-file ./filter/trio_indels.tranches \
    -O ./filter/trio_indels.recal \
    2>>./filter/indel_log.txt
 
time $gatk ApplyVQSR \
    -R $ref \
    -V ./filter/trio_snps_VQSR.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ./filter/trio_indels.tranches \
    --recal-file ./filter/trio_indels.recal \
    -mode INDEL \
    -O ./filter/trio_snp_indel_VQSR.vcf.gz \
    2>>./filter/indel_apply_log.txt
