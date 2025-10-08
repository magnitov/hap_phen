# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

mkdir -p ${DATA_PATH}/bigwig

#
# Merge replicates
#
for SAMPLE in POU2F2_ENCSR000BGP IRF4_MM1S_Loven_2013 IRF4_ENCSR000BGY CTCF_Kasowski
do
	samtools merge -f -@ 32 ${DATA_PATH}/bigwig/${SAMPLE}.hg38.bam ${DATA_PATH}/bam/${SAMPLE}_rep[0-9]_ref_Aligned.sortedByCoord.mapq.nodups.bam
	samtools index -@ 32 ${DATA_PATH}/bigwig/${SAMPLE}.hg38.bam
done

#
# RPGC-normalisation
#
for SAMPLE in POU2F2_ENCSR000BGP IRF4_MM1S_Loven_2013 IRF4_ENCSR000BGY
do
	bamCoverage -p 32 -b ${DATA_PATH}/bigwig/${SAMPLE}.hg38.bam \
		--blackListFileName ${BLACKLIST} --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig/ChIPseq_${SAMPLE}.hg38.rpgc.bw
done

#
# Input normalisation
#
for SAMPLE in CTCF_Kasowski
do
	bamCompare -p 32 --blackListFileName ${BLACKLIST} --binSize 50 --smoothLength 100 --operation ratio \
		-b1 ${DATA_PATH}/bigwig/${SAMPLE}.hg38.bam \
		-b2 ${DATA_PATH}/bigwig/Input_Kasowski.hg38.bam \
		-o ${DATA_PATH}/bigwig/ChIPseq_${SAMPLE}.hg38.input.bw
done
