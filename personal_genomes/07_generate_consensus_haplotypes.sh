#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	mkdir -p ${DATA_PATH}/${SAMPLE}/genome

	bcftools consensus --fasta-ref ${GENOME} -c ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap1.chain --haplotype 1 \
		${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.fa
	bcftools consensus --fasta-ref ${GENOME} -c ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap2.chain --haplotype 2 \
		${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz > ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.fa
done
