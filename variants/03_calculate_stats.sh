#!/bin/bash

VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	rtg vcfstats ${VCF_PATH}/${SAMPLE}/${SAMPLE}.vcf > ${VCF_PATH}/${SAMPLE}/vcfstats_${SAMPLE}.txt
	rtg vcfstats ${VCF_PATH}/${SAMPLE}/${SAMPLE}.autosomal.vcf > ${VCF_PATH}/${SAMPLE}/vcfstats_${SAMPLE}_autosomal.txt
done
