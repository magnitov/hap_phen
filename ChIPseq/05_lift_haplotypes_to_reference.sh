# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ChIPseq
GENOMES_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes
GENOME_CHROMSIZES=/DATA/users/m.magnitov/genomes/hg38.chrom.sizes

mkdir -p ${DATA_PATH}/bigwig_haplotypes

for SAMPLE in ELF1_ENCSR841NDX Input_ENCSR956WYO CTCF_Kasowski H3K27ac_Kasowski H3K27me3_Kasowski H3K36me3_Kasowski H3K4me1_Kasowski H3K4me3_Kasowski Input_Kasowski RELA_Zhao RELB_Zhao BATF_ENCSR000BGT IRF4_ENCSR000BGY JUNB_ENCSR897MMC POU2F2_ENCSR000BGP RUNX3_ENCSR000BRI SPI1_ENCSR000BGQ Input_ENCSR000BGH Input_ENCSR000BVP Input_ENCSR136WQZ Input_ENCSR890ZMI
do
	for HAPLOTYPE in hap1 hap2
	do
		# Sort BAM files by name
		echo "Preparing files for lifting"
		samtools sort -@ 32 -n -o ${DATA_PATH}/bam_assigned/${SAMPLE}_${HAPLOTYPE}.variants.name.bam \
			${DATA_PATH}/bam_assigned/${SAMPLE}_${HAPLOTYPE}.variants.bam 

		# Convert BAM to BED
		bamToBed -i ${DATA_PATH}/bam_assigned/${SAMPLE}_${HAPLOTYPE}.variants.name.bam | awk -F"/" '{ print $1 }' |\
			bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max | awk '{ print $2"\t"$3"\t"$4"\t"$1 }' |\
			awk '{ if ($3-$2 <= 2000) print $0 }' | sort -k1,1 -k2,2n \
			> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}.bed

		# Lift from haplotype to reference
		echo "Lifting files"
		CrossMap.py bed ${GENOMES_PATH}/NA12878/genome/NA12878_${HAPLOTYPE}_to_ref.chain \
			${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}.bed \
			${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted.bed
		bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted.bed |\
			awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2000) print $0 }' | sort -k1,1 -k2,2n \
			> ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted_tmp.bed
		mv ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted_tmp.bed ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted.bed

		# Convert BED to BAM
		echo "Processing lifted files"
		bedToBam -i ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted.bed -g ${GENOME_CHROMSIZES} >\
			${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.bam
		samtools sort -@ 32 -o ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam \
			${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.bam
		samtools index -@ 32 ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.sorted.bam

		# Remove intermediate files
		rm ${DATA_PATH}/bam_assigned/${SAMPLE}_${HAPLOTYPE}.variants.bam
		rm ${DATA_PATH}/bam_assigned/${SAMPLE}_${HAPLOTYPE}.variants.name.bam
		rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_hap*_lifted.bed.unmap
		rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}.bed
		rm ${DATA_PATH}/bigwig_haplotypes/fragments_${SAMPLE}_${HAPLOTYPE}_lifted.bed
		rm ${DATA_PATH}/bigwig_haplotypes/${SAMPLE}_${HAPLOTYPE}.hg38.bam
	done
done
