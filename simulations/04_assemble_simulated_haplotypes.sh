#!/bin/bash

HAPLOTYPES_PATH=/DATA/users/m.magnitov/hap_phen/simulations/haplotypes
VCF_PATH=/DATA/users/m.magnitov/hap_phen/variants/NA12878/split/chr1.vcf
HAPCUT=/DATA/users/m.magnitov/software/HapCUT2-1.3.3

mkdir -p ${HAPLOTYPES_PATH}
mkdir -p ${HAPLOTYPES_PATH}/hic_lr_10X
mkdir -p ${HAPLOTYPES_PATH}/hic_lr_stlfr
mkdir -p ${HAPLOTYPES_PATH}/hic_lr_tellseq
mkdir -p ${HAPLOTYPES_PATH}/output_hic_10X
mkdir -p ${HAPLOTYPES_PATH}/output_hic_stlfr
mkdir -p ${HAPLOTYPES_PATH}/output_hic_tellseq

for SAMPLING_NUM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
	for DEPTH_HiC in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
	do
		for DEPTH_LR in 1 2 3 4 5 6 7 8 10 12 16 20 24 32
		do
			# concatenate Hi-C fragment file and linked-reads fragment file
			cat ${HAPLOTYPES_PATH}/hic/chr1_${DEPTH_HiC}x_ds${SAMPLING_NUM} \
				${HAPLOTYPES_PATH}/lr_10X/chr1_${DEPTH_LR}x_ds${SAMPLING_NUM} > \
				${HAPLOTYPES_PATH}/hic_lr_10X/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM}
			
			cat ${HAPLOTYPES_PATH}/hic/chr1_${DEPTH_HiC}x_ds${SAMPLING_NUM} \
				${HAPLOTYPES_PATH}/lr_stlfr/chr1_${DEPTH_LR}x_ds${SAMPLING_NUM} > \
				${HAPLOTYPES_PATH}/hic_lr_stlfr/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM}
			
			cat ${HAPLOTYPES_PATH}/hic/chr1_${DEPTH_HiC}x_ds${SAMPLING_NUM} \
				${HAPLOTYPES_PATH}/lr_tellseq/chr1_${DEPTH_LR}x_ds${SAMPLING_NUM} > \
				${HAPLOTYPES_PATH}/hic_lr_tellseq/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM}

			# assemble haplotypes from combined Hi-C + linked-reads haplotype fragments
			${HAPCUT}/build/HAPCUT2 --hic 1 \
				--fragments ${HAPLOTYPES_PATH}/hic_lr_10X/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM} \
				--vcf ${VCF_PATH} \
				--output ${HAPLOTYPES_PATH}/output_hic_10X/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM}.hap \
				--htrans_data_outfile ${HAPLOTYPES_PATH}/output_hic_10X/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_10X_ds${SAMPLING_NUM}.htrans_model

			${HAPCUT}/build/HAPCUT2 --hic 1 \
				--fragments ${HAPLOTYPES_PATH}/hic_lr_stlfr/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM} \
				--vcf ${VCF_PATH} \
				--output ${HAPLOTYPES_PATH}/output_hic_stlfr/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM}.hap \
				--htrans_data_outfile ${HAPLOTYPES_PATH}/output_hic_stlfr/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_stlfr_ds${SAMPLING_NUM}.htrans_model

			${HAPCUT}/build/HAPCUT2 --hic 1 \
				--fragments ${HAPLOTYPES_PATH}/hic_lr_tellseq/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM} \
				--vcf ${VCF_PATH} \
				--output ${HAPLOTYPES_PATH}/output_hic_tellseq/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM}.hap \
				--htrans_data_outfile ${HAPLOTYPES_PATH}/output_hic_tellseq/chr1_${DEPTH_HiC}x_hic_${DEPTH_LR}x_tellseq_ds${SAMPLING_NUM}.htrans_model
		done
	done
done
