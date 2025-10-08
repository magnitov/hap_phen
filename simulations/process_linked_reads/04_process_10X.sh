LONGRANGER=/DATA/users/m.magnitov/software/longranger-2.2.2/longranger
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0
INPUT_PATH=/DATA/users/m.magnitov/hap_phen/10X/fastq
OUTPUT_PATH=/DATA/users/m.magnitov/hap_phen/simulations/linked_reads
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar
SAMTOOLS=/DATA/users/m.magnitov/software/samtools-1.15/samtools

# Align 10X linked-reads data
${LONGRANGER} align --sample bamtofastq --id 10X --fastqs ${INPUT_PATH}/NA12878 --reference ${GENOME} --localcores 32 --localmem 256

# Remove duplicates
java -jar ${PICARD} MarkDuplicates \
	INPUT=${OUTPUT_PATH}/10X/outs/possorted_bam.bam \
	OUTPUT=${OUTPUT_PATH}/10X/outs/possorted_bam.rmdup.bam \
	METRICS_FILE=${OUTPUT_PATH}/10X/outs/picard.dedup.metrics.txt \
	REMOVE_DUPLICATES=true ASSUME_SORTED=true

# Remove reads without barcode
${SAMTOOLS} view -@ 64 -h -S -b -d BX \
	-o ${OUTPUT_PATH}/10X/outs/possorted_bam.rmdup.rmbx.bam \
	${OUTPUT_PATH}/10X/outs/possorted_bam.rmdup.bam
${SAMTOOLS} index -@ 64 ${OUTPUT_PATH}/10X/outs/possorted_bam.rmdup.rmbx.bam
rm ${OUTPUT_PATH}/10X/outs/possorted_bam.rmdup.bam

