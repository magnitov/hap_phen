#!/bin/bash
PRED_PATH=/DATA/users/m.magnitov/hap_phen/chromBPNet
GENOMES_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	# Lift predicted ATAC-seq signals
	CrossMap.py bigwig ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain \
		${PRED_PATH}/pred/${SAMPLE}_hap1_chrombpnet_nobias.bw \
		${PRED_PATH}/pred/${SAMPLE}_hap1_chrombpnet_nobias.hg38.bw
	CrossMap.py bigwig ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain \
		${PRED_PATH}/pred/${SAMPLE}_hap2_chrombpnet_nobias.bw \
		${PRED_PATH}/pred/${SAMPLE}_hap2_chrombpnet_nobias.hg38.bw
	mv ${PRED_PATH}/pred/${SAMPLE}_hap1_chrombpnet_nobias.hg38.bw.bw ${PRED_PATH}/pred/${SAMPLE}_hap1_chrombpnet_nobias.hg38.bw
	mv ${PRED_PATH}/pred/${SAMPLE}_hap2_chrombpnet_nobias.hg38.bw.bw ${PRED_PATH}/pred/${SAMPLE}_hap2_chrombpnet_nobias.hg38.bw
	rm ${PRED_PATH}/pred/${SAMPLE}_hap*_chrombpnet_nobias.hg38.bw.sorted.bgr

	# Lift predicted contribution scores
	CrossMap.py bigwig ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain \
		${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.bw \
		${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.hg38.bw
	CrossMap.py bigwig ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain \
		${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.bw \
		${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.hg38.bw
	mv ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.hg38.bw.bw ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.hg38.bw
	mv ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.hg38.bw.bw ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.hg38.bw
	rm ${PRED_PATH}/contrib/${SAMPLE}_hap*.profile_scores.hg38.bw.sorted.bgr

	/DATA/users/m.magnitov/software/bigWigToBedGraph ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.bw ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.bed
	/DATA/users/m.magnitov/software/bigWigToBedGraph ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.bw ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.bed
	/DATA/users/m.magnitov/software/bigWigToBedGraph ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.hg38.bw ${PRED_PATH}/contrib/${SAMPLE}_hap1.profile_scores.hg38.bed
	/DATA/users/m.magnitov/software/bigWigToBedGraph ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.hg38.bw ${PRED_PATH}/contrib/${SAMPLE}_hap2.profile_scores.hg38.bed
done
