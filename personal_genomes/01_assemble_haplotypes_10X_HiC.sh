#!/bin/bash
SAMPLE=$1
TENX_BAM=/DATA/users/m.magnitov/hap_phen/10X/${SAMPLE}/outs/phased_possorted_bam.rmdup.rmbx.bam
HIC_PATH=/DATA/users/m.magnitov/hap_phen/HiC/aligned_hapcut/${SAMPLE}/split
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants/${SAMPLE}/split
HAPLOTYPES_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.fa
HAPCUT=/DATA/users/m.magnitov/software/HapCUT2-1.3.3

mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}
mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}/unlinked_10X_frag
mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}/10X
mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}/HiC
mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}/HiC_10X
mkdir -p ${HAPLOTYPES_PATH}/${SAMPLE}/output

for CHROM in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX
do
	# extract unlinked 10X haplotype informative reads
	${HAPCUT}/build/extractHAIRS --10X 1 \
		--bam ${TENX_BAM} \
		--VCF ${VCF_PATH}/${CHROM}.vcf \
		--indels 1 --region ${CHROM} --ref ${GENOME} > ${HAPLOTYPES_PATH}/${SAMPLE}/unlinked_10X_frag/frags_${CHROM}

	# link haplotype informative reads into 10X molecule haplotype fragments
	python ${HAPCUT}/utilities/LinkFragments.py \
		-f ${HAPLOTYPES_PATH}/${SAMPLE}/unlinked_10X_frag/frags_${CHROM} \
		-v ${VCF_PATH}/${CHROM}.vcf -b ${TENX_BAM} \
		-o ${HAPLOTYPES_PATH}/${SAMPLE}/10X/${CHROM}

	# convert Hi-C bam files to haplotype fragment files
	${HAPCUT}/build/extractHAIRS --HiC 1 \
		--bam ${HIC_PATH}/${SAMPLE}.REF_${CHROM}.bam \
		--VCF ${VCF_PATH}/${CHROM}.vcf \
		--indels 1 --ref ${GENOME} > ${HAPLOTYPES_PATH}/${SAMPLE}/HiC/${CHROM}

	# concatenate Hi-C fragment file and 10X fragment file
	cat ${HAPLOTYPES_PATH}/${SAMPLE}/HiC/${CHROM} ${HAPLOTYPES_PATH}/${SAMPLE}/10X/${CHROM} > ${HAPLOTYPES_PATH}/${SAMPLE}/HiC_10X/${CHROM}

	# assemble haplotypes from combined Hi-C + 10X haplotype fragments
	${HAPCUT}/build/HAPCUT2 --hic 1 \
		--fragments ${HAPLOTYPES_PATH}/${SAMPLE}/HiC_10X/${CHROM} \
		--vcf ${VCF_PATH}/${CHROM}.vcf \
		--output ${HAPLOTYPES_PATH}/${SAMPLE}/output/${CHROM}.hap \
		--htrans_data_outfile ${HAPLOTYPES_PATH}/${SAMPLE}/output/${CHROM}.htrans_model
done
