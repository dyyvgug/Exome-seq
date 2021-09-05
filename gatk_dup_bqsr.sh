#!/bin/bash
# author: Yingying Dong. 2021-8-25.

for item in $(ls *.sam)
do
	echo "GATK_${item%.*}"
	ref = /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
	# Duplicates Marking, chip-seq,dap-seq,call snp need this step.
	gatk --java-options "-Xmx10G -Djava.io.tmpdir=./" MarkDuplicates \
		-I ${item%.*}.bam \
		-O ${item%.*}_dup.bam \
		-M ${item%.*}.metrics \
		2>> dup_log.txt
	 # Fix mate information of the pair-end.
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation \
		-I ${item%.*}_dup.bam \
		-O ${item%.*}_dup_fixed.bam \
		-SO coordinate \
		2>>duFix_log.txt
	# samtools index 
	samtools index ${item%.*}_dup_fixed.bam
	# base quality score recalibration
    	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator -R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta -I ${item%.*}_dup_fixed_add.bam --known-sites /proj/y.dong/GATK/hg38/dbsnp_146.hg38.vcf.gz --known-sites /proj/y.dong/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ${item%.*}_recal.table 2>>bqsr_recal_log.txt
    
    	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR -R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta -I ${item%.*}_dup_fixed_add.bam -bqsr ${item%.*}_recal.table -O ${item%.*}_bqsr.bam 2>>apply_bqsr_log.txt
	# check BAM files for error 
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ValidateSamFile -I ${item%.*}_bqsr.bam 2>> check_bam_log.txt
done
