#!/bin/bash
INDIVIDUAL=NA12878
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
GENOME_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
REFERENCE_GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${DATA_PATH}/bam
mkdir -p ${DATA_PATH}/bam_assigned

for SAMPLE in RELA_Zhao_rep1 RELB_Zhao_rep1 RELB_Zhao_rep2 IRF4_MM1S_Loven_2013 BATF_ENCSR000BGT_rep1 BATF_ENCSR000BGT_rep2 IRF4_ENCSR000BGY_rep1 IRF4_ENCSR000BGY_rep2 JUNB_ENCSR897MMC_rep1 JUNB_ENCSR897MMC_rep2 POU2F2_ENCSR000BGP_rep1 POU2F2_ENCSR000BGP_rep2 POU2F2_ENCSR000BGP_rep3 RUNX3_ENCSR000BRI_rep1 RUNX3_ENCSR000BRI_rep2 SPI1_ENCSR000BGQ_rep2 SPI1_ENCSR000BGQ_rep3 Input_ENCSR000BGH_rep1 Input_ENCSR000BGH_rep2 Input_ENCSR000BGH_rep3 Input_ENCSR000BVP_rep1 Input_ENCSR000BVP_rep3 Input_ENCSR000BVP_rep4 Input_ENCSR136WQZ_rep1 Input_ENCSR890ZMI_rep1 
do
	# Mapping reads
	bwa mem -t 16 -M ${REFERENCE_GENOME} \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}.fastq.gz |\
		samtools view -@ 16 -h -S -b -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_ref_Aligned.sortedByCoord.mapq.bam

	bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap1.fa \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}.fastq.gz |\
		samtools view -@ 16 -h -S -b -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_hap1_Aligned.sortedByCoord.mapq.bam

	bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap2.fa \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}.fastq.gz |\
		samtools view -@ 16 -h -S -b -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_hap2_Aligned.sortedByCoord.mapq.bam

	# Remove duplicated reads
	java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
		I=${DATA_PATH}/bam/${SAMPLE}_ref_Aligned.sortedByCoord.mapq.bam \
		O=${DATA_PATH}/bam/${SAMPLE}_ref_Aligned.sortedByCoord.mapq.nodups.bam \
		M=${DATA_PATH}/bam/${SAMPLE}_ref_dedupMetrics.txt

	java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
		I=${DATA_PATH}/bam/${SAMPLE}_hap1_Aligned.sortedByCoord.mapq.bam \
		O=${DATA_PATH}/bam/${SAMPLE}_hap1_Aligned.sortedByCoord.mapq.nodups.bam \
		M=${DATA_PATH}/bam/${SAMPLE}_hap1_dedupMetrics.txt

	java -jar ${PICARD} MarkDuplicates REMOVE_DUPLICATES=true \
		I=${DATA_PATH}/bam/${SAMPLE}_hap2_Aligned.sortedByCoord.mapq.bam \
		O=${DATA_PATH}/bam/${SAMPLE}_hap2_Aligned.sortedByCoord.mapq.nodups.bam \
		M=${DATA_PATH}/bam/${SAMPLE}_hap2_dedupMetrics.txt

	# Remove intermediate files
	rm ${DATA_PATH}/bam/${SAMPLE}_*_Aligned.sortedByCoord.mapq.bam
done
