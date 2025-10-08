# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

#
# RPGC normalisation
#
for SAMPLE in RELA_Zhao RELB_Zhao
do
	for HAPLOTYPE in hap1 hap2
	do
		bamCoverage -p 32 --blackListFileName ${BLACKLIST} --binSize 100 --smoothLength 200 \
			--effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
			-b ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.rpgc.bw
	done
done

#
# Input normalisation
#
for SAMPLE in RUNX3_ENCSR000BRI
do
	for HAPLOTYPE in hap1 hap2
	do
		bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 100 --smoothLength 200 --operation ratio \
			-b1 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-b2 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR000BVP_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.input.bw
	done
done

for SAMPLE in SPI1_ENCSR000BGQ
do
	for HAPLOTYPE in hap1 hap2
	do
		bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 100 --smoothLength 200 --operation ratio \
			-b1 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-b2 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR000BGH_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.input.bw
	done
done

for SAMPLE in JUNB_ENCSR897MMC
do
	for HAPLOTYPE in hap1 hap2
	do
		samtools merge -@ 32 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR136WQZ_ENCSR890ZMI_${HAPLOTYPE}.hg38.sorted.bam \
			${DATA_PATH}/bigwig_haplotypes/Input_ENCSR136WQZ_${HAPLOTYPE}.hg38.sorted.bam \
			${DATA_PATH}/bigwig_haplotypes/Input_ENCSR890ZMI_${HAPLOTYPE}.hg38.sorted.bam
		samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR136WQZ_ENCSR890ZMI_${HAPLOTYPE}.hg38.sorted.bam
		bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 100 --smoothLength 200 --operation ratio \
			-b1 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			-b2 ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR136WQZ_ENCSR890ZMI_${HAPLOTYPE}.hg38.sorted.bam \
			-o ${DATA_PATH}/bigwig_haplotypes/ChIPseq_${SAMPLE}_${HAPLOTYPE}.hg38.input.bw
		rm ${DATA_PATH}/bigwig_haplotypes/Input_ENCSR136WQZ_ENCSR890ZMI_${HAPLOTYPE}.hg38.sorted.bam*
	done
done
