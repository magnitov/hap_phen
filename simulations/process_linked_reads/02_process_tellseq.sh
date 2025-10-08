# Based on the script from tellsort Docker image by Universal Sequencing and Chen et al., Genome Research (2020)

INPUT_PATH=/DATA/users/m.magnitov/hap_phen/simulations/linked_reads/fastq
OUTPUT_PATH=/DATA/users/m.magnitov/hap_phen/simulations/linked_reads/TELLseq
GENOME=/DATA/users/m.magnitov/genomes/refdata-hg38-2.1.0/fasta/genome.fa
ADD_BX=/DATA/users/m.magnitov/software/TELLseq-0.1.3/add_bx
PICARD=/DATA/users/m.magnitov/software/picard_2.27.4.jar

mkdir -p ${OUTPUT_PATH}

### adding barcodes to read name
${ADD_BX} ${INPUT_PATH}/TELLseq_NA12878_R1.fastq.gz ${INPUT_PATH}/TELLseq_NA12878_I1.fastq.gz > ${OUTPUT_PATH}/TELLseq_NA12878_R1.BX.fastq
${ADD_BX} ${INPUT_PATH}/TELLseq_NA12878_R2.fastq.gz ${INPUT_PATH}/TELLseq_NA12878_I1.fastq.gz > ${OUTPUT_PATH}/TELLseq_NA12878_R2.BX.fastq

### mapping with bwa
bwa mem -C -t 32 ${GENOME} ${OUTPUT_PATH}/TELLseq_NA12878_R1.BX.fastq ${OUTPUT_PATH}/TELLseq_NA12878_R2.BX.fastq | \
	samtools view -@ 16 -h -S -b | \
	samtools sort -@ 16 -m 5000000000 - > ${OUTPUT_PATH}/TELLseq_NA12878.sorted.bam 

#### removing duplicates in sorted bam file with picard 
java -jar ${PICARD} MarkDuplicates \
	I=${OUTPUT_PATH}/TELLseq_NA12878.sorted.bam \
	O=${OUTPUT_PATH}/TELLseq_NA12878.sorted.rmdup.bam \
	M=${OUTPUT_PATH}/dedup_metrics.txt \
	REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT TMP_DIR=./

#### index final bam file
samtools index ${OUTPUT_PATH}/TELLseq_NA12878.sorted.rmdup.bam

#### remove intermediate files
rm ${OUTPUT_PATH}/TELLseq_NA12878_R1.BX.fastq ${OUTPUT_PATH}/TELLseq_NA12878_R2.BX.fastq ${OUTPUT_PATH}/TELLseq_NA12878.sorted.bam

