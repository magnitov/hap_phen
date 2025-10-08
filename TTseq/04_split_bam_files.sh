#!/bin/bash
DATA_PATH=/DATA/users/m.magnitov/hap_phen/TTseq

cd ${DATA_PATH}/bam_assigned
mkdir -p ${DATA_PATH}/bam_assigned_split

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464
do
	for REPLICATE in rep1 rep2
	do
		for HAPLOTYPE in hap1 hap2 unassigned_hap1 unassigned_hap2
		do
			samtools view -@ 32 -b -f 128 -F 16 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam > ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for1.bam
			samtools view -@ 32 -b -f 80 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam > ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for2.bam
			samtools view -@ 32 -b -f 144 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam > ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev1.bam
			samtools view -@ 32 -b -f 64 -F 16 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}.variants.bam > ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev2.bam
		done   
        
		for HAPLOTYPE in hap1 hap2 unassigned_hap1 unassigned_hap2
		do
			samtools merge -@ 32 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward_unsorted.bam \
				${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for1.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for2.bam
			samtools sort -@ 32 -o ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward_unsorted.bam
			rm ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for1.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_for2.bam
			rm ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward_unsorted.bam

			samtools merge -@ 32 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse_unsorted.bam \
				${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev1.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev2.bam
			samtools sort -@ 32 -o ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse_unsorted.bam
			rm ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev1.bam ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_rev2.bam
			rm ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse_unsorted.bam

			samtools index -@ 16 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_forward.bam
			samtools index -@ 16 ${SAMPLE}_${REPLICATE}_${HAPLOTYPE}_reverse.bam
		done
        
	mv ${SAMPLE}_${REPLICATE}_*_forward.bam* ${DATA_PATH}/bam_assigned_split/
	mv ${SAMPLE}_${REPLICATE}_*_reverse.bam* ${DATA_PATH}/bam_assigned_split/
	done
done

