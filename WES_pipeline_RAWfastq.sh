#!/bin/bash
#SBATCH -J ES
#SBATCH -N 2
#SBATCH --ntasks-per-node=30
#SBATCH -o wes.out
#SBATCH -e wes.err
#SBATCH --mem=20G

mkdir ./fastq/
mv *.gz ./fastq/
cd ./fastq/
fastqc ./*
multiqc ./*

for item in $(ls *.fastq.gz)
do
        echo "generate ${item%.fastq.*}"
        touch ${item%.fastq.*}.sra
done

mkdir clipper_fastq

ls *_R1.fastq.gz >1
ls *_R2.fastq.gz >2

paste 1 2 > config
cat config |while read id
do
	echo "clipper_${id}"
	arr=(${id})
    	fq1=${arr[0]}
    	fq2=${arr[1]}
	trim_galore -q 33 --phred33 --length 20 -e 0.1 --stringency 3 --paired $fq1 $fq2 -o ./clipper_fastq/
done

cd ./clipper_fastq/
rename 's/_trimmed//' *
rename 's/_val_1//' *
rename 's/_val_2//' *

cd ..
mkdir aligned

for item in $(ls *.sra)
do
	echo "bwa_${item%.*}"
	bwa mem -t 60 -M -R "@RG\tID:${item%.*}\tSM:${item%.*}\tLB:WES\tPL:Illumina" /proj/y.dong/genome/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.6.0/genome.fa ./clipper_fastq/${item%.*}_R1.fq.gz ./clipper_fastq/${item%.*}_R2.fq.gz  > ./aligned/${item%.*}.sam 2>> ./aligned/mapping_repo.txt
	
done

cd ./aligned/
#mkdir stati
for item in $(ls *.sam)
do
	echo "sam_${item%.*}"

	samtools sort -@ 60 -o ${item%.*}.bam ${item%.*}.sam
	samtools flagstat ${item%.*}.bam > ./stati/${item%.*}.bam.flag

# need lots of time
#	echo "coverage_${item%.*}"
#			
#	bedtools genomecov -d -ibam ${item%.*}.bam -g /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta > ./stati/${item%.*}_cov.txt
#	qualimap bamqc -bam ${item%.*}.bam -outfile ./stati/${item%.*}.pdf -outformat PDF

	echo "gatk_dup_bqsr_${item%.*}"
	echo "GATK_${item%.*}"
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

	# base quality score recalibration
    	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
		-R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
		-I ${item%.*}_dup_fixed.bam \
		--known-sites /proj/y.dong/GATK/hg38/dbsnp_146.hg38.vcf.gz \
		--known-sites /proj/y.dong/GATK/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
		-O ${item%.*}_recal.table \
		2>>bqsr_recal_log.txt
    
    	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
		-R /proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta \
		-I ${item%.*}_dup_fixed.bam \
		-bqsr ${item%.*}_recal.table \
		-O ${item%.*}_bqsr.bam \
		2>>apply_bqsr_log.txt
	# check BAM files for error 
	
	gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ValidateSamFile \
		-I ${item%.*}_bqsr.bam \
		2>> check_bam_log.txt
	
	# samtools index 
	samtools index ${item%.*}_bqsr.bam

done

mkdir mutation
ref=/proj/y.dong/GATK/hg38/Homo_sapiens_assembly38.fasta
for item in $(ls *.sam)
do
	echo "HC_${item%.*}" 
	# use HaplotypeCaller to identify variants
		gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
		-ERC GVCF \
		--native-pair-hmm-threads 60 \
		-R $ref \
		-I ${item%.*}_bqsr.bam \
		-O ./mutation/${item%.*}.g.vcf.gz \
		2>>hc_log.txt
	
done

cd mutation
# combine trio
j=0
for i in $(ls *.gz)
do
	gvcf_list[j]=$i
	j=`expr $j + 1`
done

gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" CombineGVCFs \
	-R $ref \
	-V ./${gvcf_list[0]} \
	-V ./${gvcf_list[1]} \
	-V ./${gvcf_list[2]} \
	-O ./cohort.g.vcf.gz \
	2>>./combine_log.txt
	
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
	-R $ref \
	-V ./cohort.g.vcf.gz \
	-O ./cohort.vcf.gz \
	2>>./genotype_log.txt


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
echo "vcf_filter_${item%.*}"
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

gunzip *.gz

mkdir only_pass
cat ./trio_snp_indel_VQSR.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' \
                > ./only_pass/trio_only_pass.vcf
cd only_pass

if [ ! -d ./anno ]
then mkdir ./anno
fi


for item in $(ls *.vcf)
do
        echo "annovar convert ${item%.*}"
        convert2annovar.pl -format vcf4 ./${item%.*}.vcf > ./${item%.*}.avinput
		
		#frequent variants
	echo "annovar annotate ${item%.*}"
        table_annovar.pl \
        -thread 60 \
        ./${item%.*}.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*} \
        -remove \
        -protocol refGene,clinvar_20200316 \
        -operation g,g \
        -nastring . \
        -vcfinput
done

for item in $(ls *.avinput)
do		
		#rare variants
		echo "annovar filter rare ${item%.*}"
		annotate_variation.pl \
		-thread 60 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_rare \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.01 -reverse \
		2>>anno_fil_rare.log
		
		#frequent variants
		echo "annovar filter frequent ${item%.*}"
		annotate_variation.pl \
		-thread 60 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_fre \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.05 \
		2>>anno_fil_fre.log

done

cd ..
# not only pass
if [ ! -d ./anno ]
then mkdir ./anno
fi


for item in $(ls *.vcf)
do
        echo "annovar convert ${item%.*}"
        convert2annovar.pl -format vcf4 ./${item%.*}.vcf > ./${item%.*}.avinput
		
		#frequent variants
	echo "annovar annotate ${item%.*}"
        table_annovar.pl \
        -thread 60 \
        ./${item%.*}.avinput \
        /proj/y.dong/annovar/humandb \
        -buildver hg38 \
        -out ./anno/${item%.*} \
        -remove \
        -protocol refGene,clinvar_20200316 \
        -operation g,g \
        -nastring . \
        -vcfinput
done

for item in $(ls *.avinput)
do		
		#rare variants
		echo "annovar filter rare ${item%.*}"
		annotate_variation.pl \
		-thread 60 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_rare \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.01 -reverse \
		2>>anno_fil_rare.log
		
		#frequent variants
		echo "annovar filter frequent ${item%.*}"
		annotate_variation.pl \
		-thread 60 \
		./${item%.*}.avinput \
		-filter -dbtype gnomad_exome \
		-buildver hg38 \
		-out ./${item%.*}_fre \
		/proj/y.dong/annovar/humandb \
		-score_threshold 0.05 \
		2>>anno_fil_fre.log

done
