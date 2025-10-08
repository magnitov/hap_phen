#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/

cd ${DATA_PATH}/ChIPseq/bam_assigned

for SAMPLE in ELF1_ENCSR841NDX Input_ENCSR956WYO CTCF_Kasowski H3K27ac_Kasowski H3K27me3_Kasowski H3K36me3_Kasowski H3K4me1_Kasowski H3K4me3_Kasowski Input_Kasowski
do
	samtools merge -f -@ 32 ${SAMPLE}_hap1.bam ${SAMPLE}_rep[0-9]_hap1.bam
	samtools merge -f -@ 32 ${SAMPLE}_hap2.bam ${SAMPLE}_rep[0-9]_hap2.bam
	samtools merge -f -@ 32 ${SAMPLE}_unassigned_hap1.bam ${SAMPLE}_rep[0-9]_unassigned_hap1.bam
	samtools merge -f -@ 32 ${SAMPLE}_unassigned_hap2.bam ${SAMPLE}_rep[0-9]_unassigned_hap2.bam

	samtools sort -@ 32 -n ${SAMPLE}_hap1.bam |\
		bedtools pairtobed -type either -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap1.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_hap1.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_hap2.bam |\
		bedtools pairtobed -type either -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap2.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_hap2.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_unassigned_hap1.bam |\
		bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap1.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_unassigned_hap1.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_unassigned_hap2.bam |\
		bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap2.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_unassigned_hap2.variants.bam
            
	samtools index -@ 32 ${SAMPLE}_hap1.variants.bam
	samtools index -@ 32 ${SAMPLE}_hap2.variants.bam
	samtools index -@ 32 ${SAMPLE}_unassigned_hap1.variants.bam
	samtools index -@ 32 ${SAMPLE}_unassigned_hap2.variants.bam
done
