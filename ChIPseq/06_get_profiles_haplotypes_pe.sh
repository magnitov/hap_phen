# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

#
# Input normalisation
#
for SAMPLE in CTCF_Kasowski H3K27ac_Kasowski H3K27me3_Kasowski H3K4me1_Kasowski H3K4me3_Kasowski
do
	for HAPLOTYPE in hap1 hap2
	do
		bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 50 --smoothLength 100 --operation ratio \
			-b1 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-b2 ${DATA_PATH}/bigwig_haplotypes/Input_Kasowski_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.input.bw
	done
done
 
for SAMPLE in ELF1_ENCSR841NDX
do
	for HAPLOTYPE in hap1 hap2
	do
		bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 50 --smoothLength 100 --operation ratio \
			-b1 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-b2 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR956WYO_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.input.bw
	done
done