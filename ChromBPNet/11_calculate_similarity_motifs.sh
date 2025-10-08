#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/chromBPNet

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	cat ${DATA_PATH}/tf_modisco/motifs_${SAMPLE}_hg38/*.pfm > ${DATA_PATH}/tf_modisco/motifs_${SAMPLE}_hg38/merged.pfm
	uniprobe2meme < ${DATA_PATH}/tf_modisco/motifs_${SAMPLE}_hg38/merged.pfm > ${DATA_PATH}/tf_modisco/motifs_${SAMPLE}_hg38/merged.meme
done

mkdir -p ${DATA_PATH}/tf_modisco_clusters
meme2meme ${DATA_PATH}/tf_modisco/motifs_*_hg38/merged.meme > ${DATA_PATH}/tf_modisco_clusters/all_samples_hg38.tf_modisco.meme
tomtom -dist kullback -motif-pseudo 0.1 -min-overlap 1 -text \
	${DATA_PATH}/tf_modisco_clusters/all_samples_hg38.tf_modisco.meme ${DATA_PATH}/tf_modisco_clusters/all_samples_hg38.tf_modisco.meme \
	> ${DATA_PATH}/tf_modisco_clusters/all_samples_hg38.tomtom.txt
