#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/

cd ${DATA_PATH}/ChIPseq/bam_assigned

for SAMPLE in RELA_Zhao RELB_Zhao BATF_ENCSR000BGT IRF4_ENCSR000BGY JUNB_ENCSR897MMC POU2F2_ENCSR000BGP RUNX3_ENCSR000BRI SPI1_ENCSR000BGQ Input_ENCSR000BGH Input_ENCSR000BVP Input_ENCSR136WQZ Input_ENCSR890ZMI
do
	# Merge replicates together
	samtools merge -f -@ 32 ${SAMPLE}_hap1.bam ${SAMPLE}_rep[0-9]_hap1.bam
	samtools merge -f -@ 32 ${SAMPLE}_hap2.bam ${SAMPLE}_rep[0-9]_hap2.bam
	samtools merge -f -@ 32 ${SAMPLE}_unassigned_hap1.bam ${SAMPLE}_rep[0-9]_unassigned_hap1.bam
	samtools merge -f -@ 32 ${SAMPLE}_unassigned_hap2.bam ${SAMPLE}_rep[0-9]_unassigned_hap2.bam

	# Overlap with variants
	samtools sort -@ 32 -n ${SAMPLE}_hap1.bam |\
		bedtools intersect -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap1.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_hap1.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_hap2.bam |\
		bedtools intersect -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap2.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_hap2.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_unassigned_hap1.bam |\
		bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap1.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_unassigned_hap1.variants.bam
	samtools sort -@ 32 -n ${SAMPLE}_unassigned_hap2.bam |\
		bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/NA12878/genome/NA12878_hap2.vcf |\
		samtools sort -@ 32 - > ${SAMPLE}_unassigned_hap2.variants.bam

	# Index
	samtools index -@ 32 ${SAMPLE}_hap1.variants.bam
	samtools index -@ 32 ${SAMPLE}_hap2.variants.bam
	samtools index -@ 32 ${SAMPLE}_unassigned_hap1.variants.bam
	samtools index -@ 32 ${SAMPLE}_unassigned_hap2.variants.bam
done
