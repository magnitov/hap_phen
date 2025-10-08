#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

for SAMPLE in NA12878 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486 NA18983 HG01241 HG02601 HG03464
do
	bcftools concat --threads 24 -O v -o ${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf \
		${DATA_PATH}/${SAMPLE}/output/chr*.hap.phased.mvp.VCF
	bgzip ${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf
	tabix -p vcf ${DATA_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz
done

