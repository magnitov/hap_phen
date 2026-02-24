#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/

#grep -v '^#' NA12878_hap1.vcf | awk '{ print $1"\t1kg\tvariant\t"$2"\t"$2"\t.\t+\t.\tvariant_id="$1"_"$2"_"$3 }' > NA12878_hap1.gtf
#grep -v '^#' NA12878_hap2.vcf | awk '{ print $1"\t1kg\tvariant\t"$2"\t"$2"\t.\t+\t.\tvariant_id="$1"_"$2"_"$3 }' > NA12878_hap2.gtf

mkdir -p ${DATA_PATH}/TTseq/counts_variants

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2
	do
		for HAPLOTYPE in hap1 hap2
		do
			# Haplotype-specific reads
			htseq-count -n 32 --stranded no --nonunique all --order pos --type variant -i variant_id \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts_variants/counts_variants_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.txt
			htseq-count -n 32 --stranded no --nonunique all --order pos --type variant -i variant_id \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/${SAMPLE}_${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts_variants/counts_variants_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.txt
		done
	done
done
