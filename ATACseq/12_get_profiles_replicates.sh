DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
GENOME_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

mkdir -p ${DATA_PATH}/bigwig

# Bigwigs for each replicate mapped to reference
for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2 rep3
	do
		samtools index -@ 32 ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.nodups.bam
		bamCoverage -p 32 -b ${DATA_PATH}/bam/${SAMPLE}_${REPLICATE}_ref_Aligned.sortedByCoord.mapq.nodups.bam \
			--blackListFileName ${BLACKLIST} --extendReads --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
			-o ${DATA_PATH}/bigwig/${SAMPLE}_${REPLICATE}.bw 
	done

	samtools merge -@ 32 ${DATA_PATH}/bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam \
		${DATA_PATH}/bam/${SAMPLE}_rep*_ref_Aligned.sortedByCoord.mapq.nodups.bam
	samtools index -@ 32 ${DATA_PATH}/bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam
	bamCoverage -p 32 -b ${DATA_PATH}/bam/${SAMPLE}_merged_ref_Aligned.sortedByCoord.mapq.nodups.bam \
		--blackListFileName ${BLACKLIST} --extendReads --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig/ATACseq_${SAMPLE}.bw 
done

