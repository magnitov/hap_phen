#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/simulations
HAPLOTYPES_PATH=/DATA/users/m.magnitov/hap_phen/simulations/haplotypes
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants/NA12878/split/chr1.vcf
HAPCUT=/DATA/users/m.magnitov/software/HapCUT2-1.3.3
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.fa

mkdir -p ${HAPLOTYPES_PATH}
mkdir -p ${HAPLOTYPES_PATH}/unlinked_lr_frag_10X
mkdir -p ${HAPLOTYPES_PATH}/unlinked_lr_frag_stlfr
mkdir -p ${HAPLOTYPES_PATH}/unlinked_lr_frag_tellseq
mkdir -p ${HAPLOTYPES_PATH}/lr_10X
mkdir -p ${HAPLOTYPES_PATH}/lr_stlfr
mkdir -p ${HAPLOTYPES_PATH}/lr_tellseq
mkdir -p ${HAPLOTYPES_PATH}/hic

for SAMPLING_NUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	for DEPTH in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
	do
		TENX_BAM=${DATA_PATH}/subsample_10X/NA12878_10X_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam
		STLFR_BAM=${DATA_PATH}/subsample_stLFR/NA12878_stLFR_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam
		TELLSEQ_BAM=${DATA_PATH}/subsample_TELLseq/NA12878_TELLseq_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam
		HIC_BAM=${DATA_PATH}/subsample_HiC/NA12878_HiC_chr1_${DEPTH}x_ds${SAMPLING_NUM}.bam

		# extract unlinked linked-reads haplotype informative reads
		${HAPCUT}/build/extractHAIRS --10X 1 --bam ${TENX_BAM} --VCF ${VCF_PATH} \
			--indels 1 --region chr1 --ref ${GENOME} > ${HAPLOTYPES_PATH}/unlinked_lr_frag_10X/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM}
		${HAPCUT}/build/extractHAIRS --10X 1 --bam ${STLFR_BAM} --VCF ${VCF_PATH} \
			--indels 1 --region chr1 --ref ${GENOME} > ${HAPLOTYPES_PATH}/unlinked_lr_frag_stlfr/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM}
		${HAPCUT}/build/extractHAIRS --10X 1 --bam ${TELLSEQ_BAM} --VCF ${VCF_PATH} \
			--indels 1 --region chr1 --ref ${GENOME} > ${HAPLOTYPES_PATH}/unlinked_lr_frag_tellseq/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM}

		# link haplotype informative reads into linked-reads molecule haplotype fragments
		python ${HAPCUT}/utilities/LinkFragments.py -f ${HAPLOTYPES_PATH}/unlinked_lr_frag_10X/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM} \
			-v ${VCF_PATH} -b ${TENX_BAM} -o ${HAPLOTYPES_PATH}/lr_10X/chr1_${DEPTH}x_ds${SAMPLING_NUM} -d 20000 -m 40
		python ${HAPCUT}/utilities/LinkFragments.py -f ${HAPLOTYPES_PATH}/unlinked_lr_frag_stlfr/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM} \
			-v ${VCF_PATH} -b ${STLFR_BAM} -o ${HAPLOTYPES_PATH}/lr_stlfr/chr1_${DEPTH}x_ds${SAMPLING_NUM} -d 20000 -m 40
		python ${HAPCUT}/utilities/LinkFragments.py -f ${HAPLOTYPES_PATH}/unlinked_lr_frag_tellseq/frags_chr1_${DEPTH}x_ds${SAMPLING_NUM} \
			-v ${VCF_PATH} -b ${TELLSEQ_BAM} -o ${HAPLOTYPES_PATH}/lr_tellseq/chr1_${DEPTH}x_ds${SAMPLING_NUM} -d 20000 -m 40

		# convert Hi-C bam files to haplotype fragment files
		${HAPCUT}/build/extractHAIRS --HiC 1 --bam ${HIC_BAM} --VCF ${VCF_PATH} \
			--indels 1 --ref ${GENOME} > ${HAPLOTYPES_PATH}/hic/chr1_${DEPTH}x_ds${SAMPLING_NUM}
	done
done
