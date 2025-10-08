#!/bin/bash
INPUT_PATH=/DATA/users/m.magnitov/hap_phen/HiC/fastq
OUTPUT_PATH=/DATA/users/m.magnitov/hap_phen/HiC/aligned_hapcut
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.fa
BAMTOOLS=/DATA/users/m.magnitov/software/bamtools-2.5.2/build/bin/bamtools
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${OUTPUT_PATH}

for SAMPLE in NA12878 NA18983 HG01241 HG02601 HG03464 NA19238 NA19239 NA19240 HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 NA20509 HG03486
do
	mkdir -p ${OUTPUT_PATH}/${SAMPLE}
	mkdir -p ${OUTPUT_PATH}/${SAMPLE}/split

	echo "Mapping Hi-C reads"
	bwa mem -t 32 -5SPM ${GENOME} ${INPUT_PATH}/${SAMPLE}/${SAMPLE}_R1.fastq.gz ${INPUT_PATH}/${SAMPLE}/${SAMPLE}_R2.fastq.gz | \
		samtools view -@ 32 -h -S -b | \
		samtools sort -@ 32 - > ${OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.bam

	echo "Removing duplicates"
	java -jar ${PICARD} MarkDuplicates \
		INPUT=${OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.bam \
		OUTPUT=${OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.rmdup.bam \
		METRICS_FILE=${OUTPUT_PATH}/${SAMPLE}/dedup.metrics.txt \
		REMOVE_DUPLICATES=true ASSUME_SORTED=true

	echo "Splitting BAM file"
	${BAMTOOLS} split -in ${OUTPUT_PATH}/${SAMPLE}/${SAMPLE}.sorted.rmdup.bam \
		-reference -stub ${OUTPUT_PATH}/${SAMPLE}/split/${SAMPLE}
done
