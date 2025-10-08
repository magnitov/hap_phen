#!/bin/bash
TRUE_SET_VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants/validation
PHASED_VCF_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
CALC_HAP_STATS=/DATA/users/m.magnitov/software/HapCUT2-1.3.3/utilities/calculate_haplotype_statistics.py

mkdir -p ${PHASED_VCF_PATH}/validation

for SAMPLE in NA12878 NA20509 HG03486 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733
do
	mkdir -p ${PHASED_VCF_PATH}/validation/${SAMPLE}
	for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
	do
		python ${CALC_HAP_STATS} -i -v1 ${PHASED_VCF_PATH}/${SAMPLE}/output/${CHROM}.hap.phased.mvp.VCF -v2 ${TRUE_SET_VCF_PATH}/${SAMPLE}/${SAMPLE}.HGSVC.phased.filter.${CHROM}.het.vcf > ${PHASED_VCF_PATH}/validation/${SAMPLE}/hap_stats_${CHROM}.txt
		python ${CALC_HAP_STATS} -i -v1 ${PHASED_VCF_PATH}/${SAMPLE}/output/${CHROM}.hap.phased.mvp.VCF -v2 ${TRUE_SET_VCF_PATH}/${SAMPLE}/${SAMPLE}.HGSVC.phased.SNV.filter.${CHROM}.het.vcf > ${PHASED_VCF_PATH}/validation/${SAMPLE}/hap_stats_${CHROM}_SNV.txt
		python ${CALC_HAP_STATS} -i -v1 ${PHASED_VCF_PATH}/${SAMPLE}/output/${CHROM}.hap.phased.mvp.VCF -v2 ${TRUE_SET_VCF_PATH}/${SAMPLE}/${SAMPLE}.HGSVC.phased.INDEL.filter.${CHROM}.het.vcf > ${PHASED_VCF_PATH}/validation/${SAMPLE}/hap_stats_${CHROM}_INDEL.txt
	done
done
