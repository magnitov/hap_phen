# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
GENOMES_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENOME_CHROMSIZES=/DATA/users/m.magnitov/genomes/hg38.chrom.sizes
BLACKLIST=/DATA/users/m.magnitov/genomes/hg38-blacklist.v2.bed

mkdir -p ${DATA_PATH}/bigwig_haplotypes

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2 rep3
	do
		for HAPLOTYPE in unassigned_hap1 unassigned_hap2 #hap1 hap2
		do
			# Sort BAM files by name
			samtools sort -@ 32 -n -o ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam \
				${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam 

			# Convert BAM to BED
			bamToBed -i ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam | awk -F"/" '{ print $1 }' |\
				bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max | awk '{ print $2"\t"$3"\t"$4"\t"$1 }' \
				> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.txt

			# Remove intermediate files
			rm ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam
		done
	done

	cat ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_rep*_hap1.txt | awk '{ if ($3-$2 <= 2000) print $0 }' |\
		sort -k1,1 -k2,2n > ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1.bed
	cat ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_rep*_hap2.txt | awk '{ if ($3-$2 <= 2000) print $0 }' |\
		sort -k1,1 -k2,2n > ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2.bed
	cat ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_rep*_unassigned_hap1.txt | awk '{ if ($3-$2 <= 2000) print $0 }' |\
		sort -k1,1 -k2,2n > ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1.bed
	cat ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_rep*_unassigned_hap2.txt | awk '{ if ($3-$2 <= 2000) print $0 }' |\
		sort -k1,1 -k2,2n > ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2.bed

	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1.bed \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted.bed
	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2.bed \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted.bed

	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap1_to_ref.chain \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1.bed \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted.bed
	CrossMap.py bed ${GENOMES_PATH}/${SAMPLE}/genome/${SAMPLE}_hap2_to_ref.chain \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2.bed \
		${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted.bed

	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2500) print $0 }' | sort -k1,1 -k2,2n \
		> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted_tmp.bed
	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2500) print $0 }' | sort -k1,1 -k2,2n \
		> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted_tmp.bed
	mv ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted_tmp.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted.bed
	mv ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted_tmp.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted.bed

	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2500) print $0 }' | sort -k1,1 -k2,2n \
		> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted_tmp.bed
	bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted.bed |\
		awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2500) print $0 }' | sort -k1,1 -k2,2n \
		> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted_tmp.bed
	mv ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted_tmp.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted.bed
	mv ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted_tmp.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted.bed

	bedToBam -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1_lifted.bed -g ${GENOME_CHROMSIZES} >\
		${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.bam
	bedToBam -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2_lifted.bed -g ${GENOME_CHROMSIZES} >\
		${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.bam
	samtools sort -@ 32 -o ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.sorted.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.bam
	samtools sort -@ 32 -o ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.sorted.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.bam        
	samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.sorted.bam
	samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.sorted.bam

	bedToBam -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1_lifted.bed -g ${GENOME_CHROMSIZES} >\
		${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.bam
	bedToBam -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2_lifted.bed -g ${GENOME_CHROMSIZES} >\
		${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.bam
	samtools sort -@ 32 -o ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.sorted.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.bam
	samtools sort -@ 32 -o ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.sorted.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.bam
	samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.sorted.bam
	samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.sorted.bam

	bamCoverage -p 32 -b ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.sorted.bam \
		--blackListFileName ${BLACKLIST} --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig_haplotypes/ATACseq_${SAMPLE}_hap1.hg38.bw
	bamCoverage -p 32 -b ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.sorted.bam \
		--blackListFileName ${BLACKLIST} --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig_haplotypes/ATACseq_${SAMPLE}_hap2.hg38.bw

	bamCoverage -p 32 -b ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.sorted.bam \
		--blackListFileName ${BLACKLIST} --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig_haplotypes/ATACseq_${SAMPLE}_unassigned_hap1.hg38.bw
	bamCoverage -p 32 -b ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.sorted.bam \
		--blackListFileName ${BLACKLIST} --binSize 50 --effectiveGenomeSize 2913022398 --normalizeUsing RPGC \
		-o ${DATA_PATH}/bigwig_haplotypes/ATACseq_${SAMPLE}_unassigned_hap2.hg38.bw

	rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap*_lifted.bed.unmap
	rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap1.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap2.bed
	rm ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap1.hg38.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_hap2.hg38.bam
	rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap*_lifted.bed.unmap
	rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap1.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_unassigned_hap2.bed
	rm ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap1.hg38.bam ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_unassigned_hap2.hg38.bam
done
