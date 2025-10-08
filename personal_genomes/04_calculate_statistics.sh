#!/bin/bash
VCF_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

for SAMPLE in NA12878 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486 NA18983 HG01241 HG02601 HG03464
do
	rtg vcfstats ${VCF_PATH}/${SAMPLE}/${SAMPLE}.phased.het.vcf.gz > ${VCF_PATH}/${SAMPLE}/vcfstats_phased_${SAMPLE}.txt
done
