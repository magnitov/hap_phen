#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/

cd ${DATA_PATH}/TTseq/bam_assigned

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2
	do      
		samtools sort -@ 32 -n ${SAMPLE}_${REPLICATE}_hap1.bam |\
			bedtools pairtobed -type either -abam stdin -b ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap1.vcf |\
			samtools sort -@ 32 - > ${SAMPLE}_${REPLICATE}_hap1.variants.bam
		samtools sort -@ 32 -n ${SAMPLE}_${REPLICATE}_hap2.bam |\
			bedtools pairtobed -type either -abam stdin -b ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap2.vcf |\
			samtools sort -@ 32 - > ${SAMPLE}_${REPLICATE}_hap2.variants.bam
		samtools index -@ 16 ${SAMPLE}_${REPLICATE}_hap1.variants.bam
		samtools index -@ 16 ${SAMPLE}_${REPLICATE}_hap2.variants.bam
            
		samtools sort -@ 32 -n ${SAMPLE}_${REPLICATE}_unassigned_hap1.bam |\
			bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap1.vcf |\
			samtools sort -@ 32 - > ${SAMPLE}_${REPLICATE}_unassigned_hap1.variants.bam
		samtools sort -@ 32 -n ${SAMPLE}_${REPLICATE}_unassigned_hap2.bam |\
			bedtools pairtobed -type neither -abam stdin -b ${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_hap2.vcf |\
			samtools sort -@ 32 - > ${SAMPLE}_${REPLICATE}_unassigned_hap2.variants.bam
		samtools index -@ 16 ${SAMPLE}_${REPLICATE}_unassigned_hap1.variants.bam
		samtools index -@ 16 ${SAMPLE}_${REPLICATE}_unassigned_hap2.variants.bam
	done
done
