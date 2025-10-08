#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/chromBPNet
JASPAR2024_MOTIF_DATABASE=/DATA/databases/JASPAR2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

mkdir -p ${DATA_PATH}/tf_modisco

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	modisco motifs -i ${DATA_PATH}/contrib/${SAMPLE}_hg38.profile_scores.h5 -n 1000000 \
		-o ${DATA_PATH}/tf_modisco/${SAMPLE}_hg38.modisco_results.h5
	modisco report -i ${DATA_PATH}/tf_modisco/${SAMPLE}_hg38.modisco_results.h5 \
		-o ${DATA_PATH}/tf_modisco/report_${SAMPLE}_hg38/ -s ${DATA_PATH}/tf_modisco/report_${SAMPLE}_hg38/ \
		-m ${JASPAR2024_MOTIF_DATABASE} -n 5
	modisco meme -i ${DATA_PATH}/tf_modisco/${SAMPLE}_hg38.modisco_results.h5 -t PFM \
		-o ${DATA_PATH}/tf_modisco/${SAMPLE}_hg38.modisco_results.meme
done
