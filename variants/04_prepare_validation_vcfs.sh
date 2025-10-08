#!/bin/bash
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants
GENOME=/DATA/users/m.magnitov/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
SNPSIFT=/DATA/users/m.magnitov/software/snpEff_v4_3p/SnpSift.jar

###
### extract chr1 variants from NA12878 Platinum Genomes
###
java -jar ${SNPSIFT} filter "( CHROM = 'chr1' ) & isHet( GEN[0] )" NA12878.platinum.phased.vcf.gz > NA12878.platinum.phased.chr1.het.vcf
mv NA12878.platinum.phased* ${VCF_PATH}/NA12878

###
### extract variants from HGSVC haplotype-resolved assemblies
###
mkdir -p ${VCF_PATH}/validation
for SAMPLE in NA12878 NA20509 HG03486 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733
do
	mkdir -p ${VCF_PATH}/validation/${SAMPLE}
	for VAR_TYPE in SNV INDEL
	do
		bcftools view --threads 64 --trim-alt-alleles -s ${SAMPLE} -O z -o ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.vcf.gz HGSVC_20200814_freeze3.${VAR_TYPE}.phased.vcf.gz
		java -jar ${SNPSIFT} filter -v "(countVariant() > 0 )" ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.vcf.gz > ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.filter.vcf
		for CHROMOSOME in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
		do
			java -jar ${SNPSIFT} filter "( CHROM = '${CHROMOSOME}' ) & isHet( GEN[0] )" ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.filter.vcf > ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.filter.${CHROMOSOME}.het.vcf
		done
		mv ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.filter.chr*.het.vcf ${VCF_PATH}/validation/${SAMPLE}/
		rm ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.filter.vcf ${SAMPLE}.HGSVC.phased.${VAR_TYPE}.vcf.gz
	done
done

###
### merge SNVs and INDELs
###
for SAMPLE in NA12878 NA20509 HG03486 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733
do
	cd  ${VCF_PATH}/validation/${SAMPLE}/
	for CHROMOSOME in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
	do
		bcftools concat -o ${SAMPLE}.HGSVC.phased.filter.${CHROMOSOME}.het.vcf ${SAMPLE}.HGSVC.phased.INDEL.filter.${CHROMOSOME}.het.vcf ${SAMPLE}.HGSVC.phased.SNV.filter.${CHROMOSOME}.het.vcf
	done
done
