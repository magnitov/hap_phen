#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	CrossMap.py vcf --no-comp-alleles \
		${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap1.chain ${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz \
		${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.fa ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1.vcf
	CrossMap.py vcf --no-comp-alleles \
		${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_ref_to_hap2.chain ${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz \
		${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.fa ${DATA_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2.vcf
done
