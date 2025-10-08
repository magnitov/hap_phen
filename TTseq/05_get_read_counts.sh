#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/

mkdir -p ${DATA_PATH}/TTseq/counts

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2
	do
		for HAPLOTYPE in hap1 hap2
		do
			# Haplotype-specific reads
			htseq-count -n 32 --stranded no --nonunique all --order pos --type gene \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts/counts_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.txt
			htseq-count -n 32 --stranded no --nonunique all --order pos --type gene \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts/counts_${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.txt
			
			# Unassigned reads
			htseq-count -n 32 --stranded no --nonunique all --order pos --type gene \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_unassigned_${HAPLOTYPE}_forward.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts/counts_${SAMPLE}_${REPLICATE}_unassigned_${HAPLOTYPE}_forward.txt
			htseq-count -n 32 --stranded no --nonunique all --order pos --type gene \
				${DATA_PATH}/TTseq/bam_assigned_split/${SAMPLE}_${REPLICATE}_unassigned_${HAPLOTYPE}_reverse.bam \
				${DATA_PATH}/personal_genomes/${SAMPLE}/genome/gencode.v42.basic.${SAMPLE}.${HAPLOTYPE}.gtf > \
				${DATA_PATH}/TTseq/counts/counts_${SAMPLE}_${REPLICATE}_unassigned_${HAPLOTYPE}_reverse.txt
		done
	done
done
