#!/bin/bash
INDIVIDUAL=NA12878
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
GENOME_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
REFERENCE_GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${DATA_PATH}/bam
mkdir -p ${DATA_PATH}/bam_assigned

for SAMPLE in ELF1_ENCSR841NDX_rep1 ELF1_ENCSR841NDX_rep2 Input_ENCSR956WYO_rep1 Input_ENCSR956WYO_rep2 CTCF_Kasowski_rep1 CTCF_Kasowski_rep2 H3K27ac_Kasowski_rep1 H3K27ac_Kasowski_rep2 H3K27ac_Kasowski_rep3 H3K27ac_Kasowski_rep4 H3K27me3_Kasowski_rep1 H3K27me3_Kasowski_rep2 H3K4me1_Kasowski_rep1 H3K4me1_Kasowski_rep2 H3K4me1_Kasowski_rep3 H3K4me3_Kasowski_rep1 H3K4me3_Kasowski_rep2 H3K4me3_Kasowski_rep3 Input_Kasowski_rep1
do
	# Mapping paired-end reads
	bwa mem -t 16 -M ${REFERENCE_GENOME} \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R2.fastq.gz |\
		samtools view -@ 16 -h -S -b -f 2 -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_ref_Aligned.sortedByCoord.mapq.bam

	bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap1.fa \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R2.fastq.gz |\
		samtools view -@ 16 -h -S -b -f 2 -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_hap1_Aligned.sortedByCoord.mapq.bam

	bwa mem -t 16 -M ${GENOME_PATH}/${INDIVIDUAL}/genome/${INDIVIDUAL}_hap2.fa \
		${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R1.fastq.gz ${DATA_PATH}/fastq/${INDIVIDUAL}/${SAMPLE}_R2.fastq.gz |\
		samtools view -@ 16 -h -S -b -f 2 -q 10 |\
		samtools sort -@ 16 - > ${DATA_PATH}/bam/${SAMPLE}_hap2_Aligned.sortedByCoord.mapq.bam

	# Mark duplicated read pairs
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

