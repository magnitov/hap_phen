# !/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/ATACseq
VARIANTS_PATH=/DATA/users/m.magnitov/hap_phen/personal_genomes

mkdir -p ${DATA_PATH}/counts_variants

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2 rep3
	do
		for HAPLOTYPE in hap1 hap2
		do
			# Sort BAM files by name
			samtools sort -@ 32 -n -o ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam \
				${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam 

			# Convert BAM to BED (paired-end reads), convert to fragments, then overlap peaks with fragments
			# Individual peaks
			bedtools bamtobed -i ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam | awk -F"/" '{ print $1 }' | \
				bedtools groupby -g 4 -c 1,2,3 -o distinct,min,max | awk '{ print $2"\t"$3"\t"$4"\t"$1 }' | awk '{ if ($3-$2 <= 2000) print $0 }' |\
				bedtools intersect -a ${VARIANTS_PATH}/${SAMPLE}/genome/${SAMPLE}_${HAPLOTYPE//unassigned_/}.vcf -b stdin -c |\
				cut -f1,2,3,11 > ${DATA_PATH}/counts_variants/counts_variants_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.txt

			# Remove temporary files
			rm ${DATA_PATH}/bam_assigned/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.name.bam
		done
	done
	paste ${DATA_PATH}/counts_variants/counts_variants_${SAMPLE}_* | cut -f1-4,8,12,16,20,24 > ${DATA_PATH}/counts_variants/allelic_counts_per_variant_${SAMPLE}.txt
done
